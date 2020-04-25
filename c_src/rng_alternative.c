/* Copyright (C) 2016-2020 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016 Marcel Langenberg
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

//! \file rng_alternative.c
//! \brief Implementation of rng_alternative.h

/*Random number generator Mersenne-Twister*/

#include <limits.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include "phase.h"
#include "rng_alternative.h"
#include "rng.h"

uint64_t get_new_alternative_rng_offset(struct Phase *p)
{
    switch (p->args.pseudo_random_number_generator_arg)
        {
        case pseudo_random_number_generator__NULL:
            break;
        case pseudo_random_number_generator_arg_PCG32:
            break;
        case pseudo_random_number_generator_arg_MT:
            /*intentionally falls through */
        case pseudo_random_number_generator_arg_TT800:
            if (p->rh.used - 1 >= p->rh.length)
                reallocate_rng_heavy(p);
            p->rh.used += 1;
            return p->rh.used - 1;
        }
    return UINT64_MAX;
}

int copyin_rng_heavy(struct Phase *p)
{
    if (p->rh.allocated_device)
        return update_device_rng_heavy(p);

    switch (p->args.pseudo_random_number_generator_arg)
        {
        case pseudo_random_number_generator__NULL:
            break;
        case pseudo_random_number_generator_arg_PCG32:
            break;
        case pseudo_random_number_generator_arg_MT:
            assert(p->rh.mt_state != NULL);
#pragma acc enter data copyin(p->rh.mt_state[0:p->rh.length])
        case pseudo_random_number_generator_arg_TT800:
            assert(p->rh.tt800_state != NULL);
#pragma acc enter data copyin(p->rh.tt800_state[0:p->rh.length])
            break;
        }
    p->rh.allocated_device = true;
    return 0;
}

int copyout_rng_heavy(struct Phase *p)
{
    assert(p->rh.allocated_device);
    switch (p->args.pseudo_random_number_generator_arg)
        {
        case pseudo_random_number_generator__NULL:
            break;
        case pseudo_random_number_generator_arg_PCG32:
            break;
        case pseudo_random_number_generator_arg_MT:
            assert(p->rh.mt_state != NULL);
#pragma acc exit data copyout(p->rh.mt_state[0:p->rh.length])
        case pseudo_random_number_generator_arg_TT800:
            assert(p->rh.tt800_state != NULL);
#pragma acc exit data copyout(p->rh.tt800_state[0:p->rh.length])
            break;
        }
    p->rh.allocated_device = false;
    return 0;
}

int update_self_rng_heavy(struct Phase *p)
{
    assert(p->rh.allocated_device);
    switch (p->args.pseudo_random_number_generator_arg)
        {
        case pseudo_random_number_generator__NULL:
            break;
        case pseudo_random_number_generator_arg_PCG32:
            break;
        case pseudo_random_number_generator_arg_MT:
            assert(p->rh.mt_state != NULL);
#pragma acc update self(p->rh.mt_state[0:p->rh.used])
        case pseudo_random_number_generator_arg_TT800:
            assert(p->rh.tt800_state != NULL);
#pragma acc update self(p->rh.tt800_state[0:p->rh.length])
            break;
        }
    return 0;
}

int update_device_rng_heavy(struct Phase *p)
{
    assert(p->rh.allocated_device);
    switch (p->args.pseudo_random_number_generator_arg)
        {
        case pseudo_random_number_generator__NULL:
            break;
        case pseudo_random_number_generator_arg_PCG32:
            break;
        case pseudo_random_number_generator_arg_MT:
            assert(p->rh.mt_state != NULL);
#pragma acc update device(p->rh.mt_state[0:p->rh.used])
        case pseudo_random_number_generator_arg_TT800:
            assert(p->rh.tt800_state != NULL);
#pragma acc update device(p->rh.tt800_state[0:p->rh.length])
            break;
        }
    return 0;
}

int reallocate_rng_heavy(struct Phase *p)
{
    if (p->rh.allocated_device)
        copyout_rng_heavy(p);
    const uint64_t new_length = p->rh.length * 1.05 + 1;
    switch (p->args.pseudo_random_number_generator_arg)
        {
        case pseudo_random_number_generator__NULL:
            break;
        case pseudo_random_number_generator_arg_PCG32:
            break;
        case pseudo_random_number_generator_arg_MT:
            ;
            MERSENNE_TWISTER_STATE *tmp_mt =
                (MERSENNE_TWISTER_STATE *) malloc(new_length * sizeof(MERSENNE_TWISTER_STATE));
            MALLOC_ERROR_CHECK(tmp_mt, new_length * sizeof(MERSENNE_TWISTER_STATE));
            memcpy(p->rh.mt_state, tmp_mt, p->rh.used * sizeof(MERSENNE_TWISTER_STATE));
            free(p->rh.mt_state);
            p->rh.mt_state = tmp_mt;
            p->rh.length = new_length;
            break;
        case pseudo_random_number_generator_arg_TT800:
            ;
            MTTSTATE *tmp_tt800 = (MTTSTATE *) malloc(new_length * sizeof(MTTSTATE));
            MALLOC_ERROR_CHECK(tmp_tt800, new_length * sizeof(MTTSTATE));
            memcpy(p->rh.tt800_state, tmp_tt800, p->rh.used * sizeof(MTTSTATE));
            free(p->rh.tt800_state);
            p->rh.tt800_state = tmp_tt800;
            p->rh.length = new_length;
            break;
        }
    copyin_rng_heavy(p);
    return 0;
}

int init_rng_heavy(struct Phase *p, const uint64_t target_length)
{
    // Indepent of RNG set to NULL to be free resistant
    p->rh.mt_state = NULL;
    p->rh.tt800_state = NULL;
    p->rh.allocated_device = false;

    // Trick to force reallocate_rng_heavy() to do the initial initialization
    // used to 0, because that is the memcpy length
    p->rh.used = 0;
    // length to target to get this length of states
    p->rh.length = target_length;
    return reallocate_rng_heavy(p);
}

int free_rng_heavy(struct Phase *p)
{
    free(p->rh.mt_state);
    p->rh.mt_state = NULL;
    free(p->rh.tt800_state);
    p->rh.tt800_state = NULL;
    p->rh.used = 0;
    p->rh.length = 0;
    return 0;
}

int soma_seed_rng_mt(PCG_STATE * rng, MERSENNE_TWISTER_STATE * mt_rng)
{
    mt_rng->internal_index = MTMAX_num_int_state + 1;
    mt_rng->A[0] = 0;
    mt_rng->A[1] = 0x9908b0df;
    for (int i = 0; i < MTMAX_num_int_state; i++)
        {
            mt_rng->internal_state[i] = pcg32_random(rng);
        }
    return 0;
}

#pragma acc routine(soma_mersenne_twister) seq
unsigned int soma_mersenne_twister(MERSENNE_TWISTER_STATE * mt_rng)
{

    unsigned int M = 397;
    uint32_t HI = 0x80000000;
    uint32_t LO = 0x7fffffff;
    uint32_t e;

    /* The Mersenne-Twister stete of 624 is seeded with soma_seed_rng_mt()
     * which is called in the init.c by function init_values()
     */
    if (M > MTMAX_num_int_state)
        M = MTMAX_num_int_state / 2;

    if (mt_rng->internal_index >= MTMAX_num_int_state)
        {
            /* Berechne neuen Zustandsvektor */
            uint32_t h;

            for (unsigned int i = 0; i < MTMAX_num_int_state - M; ++i)
                {
                    h = (mt_rng->internal_state[i] & HI) | (mt_rng->internal_state[i + 1] & LO);        // Crashes HERE!!!
                    mt_rng->internal_state[i] = mt_rng->internal_state[i + M] ^ (h >> 1) ^ mt_rng->A[h & 1];
                }
            for (unsigned int i = MTMAX_num_int_state - M; i < MTMAX_num_int_state - 1; ++i)
                {
                    h = (mt_rng->internal_state[i] & HI) | (mt_rng->internal_state[i + 1] & LO);
                    mt_rng->internal_state[i] =
                        mt_rng->internal_state[i + (M - MTMAX_num_int_state)] ^ (h >> 1) ^ mt_rng->A[h & 1];
                }

            h = (mt_rng->internal_state[MTMAX_num_int_state - 1] & HI) | (mt_rng->internal_state[0] & LO);
            mt_rng->internal_state[MTMAX_num_int_state - 1] =
                mt_rng->internal_state[M - 1] ^ (h >> 1) ^ mt_rng->A[h & 1];
            mt_rng->internal_index = 0;
        }

    e = mt_rng->internal_state[mt_rng->internal_index++];
    /* Tempering */
    e ^= (e >> 11);
    e ^= (e << 7) & 0x9d2c5680;
    e ^= (e << 15) & 0xefc60000;
    e ^= (e >> 18);

    return e;
}

unsigned int soma_rng_uint_max_mt()
{
    return 0x80000000;
}

int soma_seed_rng_tt800(PCG_STATE * rng, MTTSTATE * tt800_rng)
{

    tt800_rng->internal_index = MTMAX_num_int_state + 1;
    tt800_rng->A[0] = 0;
    tt800_rng->A[1] = 0x8ebfd028;

    for (int k = 0; k < TT_num_int_state; k++)
        {
            tt800_rng->internal_state[k] = (uint32_t) pcg32_random(rng);
        }
    return 0;
}

#pragma acc routine(soma_tt800) seq
unsigned int soma_rng_tt800(MTTSTATE * itt800_rng)
{

    uint32_t M = 7;
    uint32_t e;
    uint32_t k = 0;
    if (itt800_rng->internal_index >= TT_num_int_state)
        {

            if (itt800_rng->internal_index > TT_num_int_state)
                {
                    uint32_t r = 9;
                    uint32_t s = 3402;
                    for (k = 0; k < TT_num_int_state; ++k)
                        {
                            r = 509845221 * r + 3;
                            s *= s + 1;
                            itt800_rng->internal_state[k] = s + (r >> 10);
                        }
                }
            for (k = 0; k < TT_num_int_state - M; ++k)
                {
                    const uint32_t tmp =
                        (itt800_rng->internal_state[k] >> 1) ^ itt800_rng->A[itt800_rng->internal_state[k] & 1];
                    itt800_rng->internal_state[k] = itt800_rng->internal_state[k + M] ^ tmp;
                }
            for (k = TT_num_int_state - M; k < TT_num_int_state; ++k)
                {
                    const uint32_t tmp =
                        (itt800_rng->internal_state[k] >> 1) ^ itt800_rng->A[itt800_rng->internal_state[k] & 1];
                    itt800_rng->internal_state[k] = itt800_rng->internal_state[k + (M - TT_num_int_state)] ^ tmp;
                }
            itt800_rng->internal_index = 0;
        }

    e = itt800_rng->internal_state[itt800_rng->internal_index++];
    /* Tempering */
    e ^= (e << 7) & 0x2b5b2500;
    e ^= (e << 15) & 0xdb8b0000;
    e ^= (e >> 16);
    return e;
}
