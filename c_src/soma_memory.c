/* Copyright (C) 2016-2021 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016-2017 Marcel Langenberg
   Copyright (C) 2016 Fabien Leonforte
   Copyright (C) 2016 Juan Orozco
   Copyright (C) 2016 Yongzhi Ren

 This file is part of SOMA.

 SOMA is free software: you can redistribute it and/or modify
 it under the terms of the GNU Lesser General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SOMA is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Lesser General Public License for more details.

 You should have received a copy of the GNU Lesser General Public License
 along with SOMA.  If not, see <http://www.gnu.org/licenses/>.
*/

/*! \file rng_alternative.h
  \brief Definition of alternative pseudo random numbers generators for SOMA.
  PCG32 is the default RNG generation engine, if a different one is used, the definitions for the alternative are in this file.
*/

#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <assert.h>

#include "soma_memory.h"
#include "soma_util.h"

int init_soma_memory(struct SomaMemory *state, const uint64_t length, const size_t typelength)
{
    state->used = 0;
    state->ptr = NULL;
    state->device_present = false;
    state->typelength = typelength;

    state->length = 0;
    if (length > 0)
        {
            // Trick to force reallocate_soma_memory() to do the initial initialization
            // used to 0, because that is the memcpy length
            // length to target to get this length of states
            return reallocate_soma_memory(state, length);
        }

    return 0;
}

int free_soma_memory(struct SomaMemory *state)
{
    int error = 0;
    if (state->device_present)
        error |= copyout_soma_memory(state);

    free(state->ptr);
    state->ptr = NULL;
    state->length = 0;
    state->used = 0;
    state->typelength = 0;
    return error;
}

int reallocate_soma_memory(struct SomaMemory *state, const uint64_t min_increase)
{
    assert((int)min_increase > 0);
    bool needs_copyin = false;
    int error = 0;
    if (state->device_present)
        {
            needs_copyin = true;
            error |= copyout_soma_memory(state);
        }

    assert(state->typelength > 0);

    const uint64_t optionA = state->length * 1.05 + min_increase;
    const uint64_t optionB = state->length + min_increase * 10;

    const uint64_t new_length = optionA < optionB ? optionA : optionB;
    void *tmp = realloc(state->ptr, new_length * state->typelength);
    if (tmp == NULL)
        {
            printf("ERROR: %s:%d %lu %lu %lu %d\n", __FILE__, __LINE__, state->length, min_increase, new_length,
                   (int)state->typelength);
            return -1;
        }
    state->ptr = tmp;
    state->length = new_length;
    if (needs_copyin)
        error |= copyin_soma_memory(state);

    return error;
}

int copyin_soma_memory(struct SomaMemory *state)
{
    if (state->length > 0)
        {
            if (state->device_present)
                return update_device_soma_memory(state);

            assert(state->ptr != NULL);
            assert(state->typelength > 0);
#pragma acc enter data copyin(state->ptr[0:state->length*state->typelength])
            state->device_present = true;
        }
    return 0;
}

int copyout_soma_memory(struct SomaMemory *state)
{
    if (state->length > 0)
        {
            assert(state->device_present);
            assert(state->ptr != NULL);
            assert(state->typelength > 0);
#pragma acc exit data copyout(state->ptr[0:state->length*state->typelength])
            state->device_present = false;
        }
    return 0;
}

int update_self_soma_memory(struct SomaMemory *state)
{
    if (state->length > 0)
        {
            assert(state->device_present);
            assert(state->ptr != NULL);
            assert(state->typelength > 0);
#pragma acc update self(state->ptr[0:state->used*state->typelength])
        }
    return 0;
}

int update_device_soma_memory(struct SomaMemory *state)
{
    if (state->length > 0)
        {
            assert(state->device_present);
            assert(state->ptr != NULL);
            assert(state->typelength > 0);
#pragma acc update device(state->ptr[0:state->used*state->typelength])
        }
    return 0;                   //Silence compiler warning
}

uint64_t get_new_soma_memory_offset(struct SomaMemory *state, const uint64_t n)
{
    assert(n > 0);
    if (state->used + n >= state->length)
        if (reallocate_soma_memory(state, n + 1) != 0)
            {
                fprintf(stderr, "ERROR: invalid memory extension %s:%d, %lu %lu %lu.\n", __FILE__, __LINE__,
                        state->used, n, state->length);
                return UINT64_MAX;
            }
    state->used += n;
    return state->used - n;
}
