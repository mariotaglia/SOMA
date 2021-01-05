/* Copyright (C) 2016-2021 Ludwig Schneider
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

//! \file rng.c
//! \brief Implementation of rng.h

#include "rng.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "phase.h"
#include "rng_alternative.h"

#pragma acc routine(pcg32_random) seq
uint32_t pcg32_random(PCG_STATE * rng)
{
    const uint64_t old = rng->state;
    // Advance internal state
    rng->state = ((uint64_t) rng->state) * 0X5851F42D4C957F2DULL;
    rng->state += (rng->inc | 1);
    const uint32_t xorshifted = ((old >> 18u) ^ old) >> 27u;
    const uint32_t rot = old >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

int soma_seed_rng(PCG_STATE * rng, uint64_t seed, uint64_t stream)
{
    rng->inc = stream * 2 + 1;
    rng->state = 0;
    pcg32_random(rng);
    rng->state += seed;
    pcg32_random(rng);
    //Improve quality of first random numbers
    pcg32_random(rng);
    return 0;
}

unsigned int soma_rng_uint_max()
{
    return 4294967295U;
}

unsigned int soma_rng_uint(RNG_STATE * state, const Phase * const p)
{
    switch (p->args.pseudo_random_number_generator_arg)
        {
        case pseudo_random_number_generator__NULL:
            return (unsigned int)pcg32_random(&(state->default_state));
            break;
        case pseudo_random_number_generator_arg_PCG32:
            return (unsigned int)pcg32_random(&(state->default_state));
            break;
        case pseudo_random_number_generator_arg_MT:
            ;
            MERSENNE_TWISTER_STATE *mt_state = p->rh.mt.ptr;
            mt_state += state->alternative_rng_offset;
            return (unsigned int)soma_mersenne_twister(mt_state);
            break;
        case pseudo_random_number_generator_arg_TT800:
            ;
            TT800STATE *tt800_state = p->rh.tt800.ptr;
            tt800_state += state->alternative_rng_offset;
            return (unsigned int)soma_rng_tt800(tt800_state);
            break;
        }
    return -1;
}

soma_scalar_t soma_rng_soma_scalar(RNG_STATE * rng, const struct Phase *const p)
{
    return ldexp(soma_rng_uint(rng, p), -32);
}

/*! generate 3D vector, 2 times Box-Mueller Transform, discards one value
*/
void soma_normal_vector(RNG_STATE * rng, const struct Phase *const p, soma_scalar_t * x,
                        soma_scalar_t * y, soma_scalar_t * z)
{
    soma_scalar_t u1, u2, u3, u4, r1, r2;

    u1 = 2 * soma_rng_soma_scalar(rng, p) - 1.;
    u2 = 2 * soma_rng_soma_scalar(rng, p) - 1.;
    u3 = 2 * soma_rng_soma_scalar(rng, p) - 1.;
    u4 = 2 * soma_rng_soma_scalar(rng, p) - 1.;

    r1 = u1 * u1 + u2 * u2;
    r2 = u3 * u3 + u4 * u4;

    while (r1 > 1)
        {
            u1 = 2 * soma_rng_soma_scalar(rng, p) - 1.;
            u2 = 2 * soma_rng_soma_scalar(rng, p) - 1.;
            r1 = u1 * u1 + u2 * u2;
        }

    while (r2 > 1)
        {
            u3 = 2 * soma_rng_soma_scalar(rng, p) - 1.;
            u4 = 2 * soma_rng_soma_scalar(rng, p) - 1.;
            r2 = u3 * u3 + u4 * u4;
        }

    const soma_scalar_t root1 = sqrt(-2.0 * log(r1) / r1);
    const soma_scalar_t root2 = sqrt(-2.0 * log(r2) / r2);

    *x = root1 * u1;
    *y = root1 * u2;
    *z = root2 * u3;
}

unsigned int rng_state_serial_length(const struct Phase *const p)
{
    unsigned int length = 0;
    length += sizeof(PCG_STATE);
    if (p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_MT)
        length += sizeof(MERSENNE_TWISTER_STATE);
    if (p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_TT800)
        length += sizeof(TT800STATE);
    return length;
}

int serialize_rng_state(struct Phase *const p, const RNG_STATE * const state, unsigned char *const buffer)
{
    unsigned int position = 0;
    //default state
    memcpy(buffer + position, &(state->default_state), sizeof(PCG_STATE));
    position += sizeof(PCG_STATE);

    if (p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_MT)
        {
            memcpy(buffer + position, ((MERSENNE_TWISTER_STATE *) p->rh.mt.ptr) + state->alternative_rng_offset,
                   sizeof(MERSENNE_TWISTER_STATE));
            position += sizeof(MERSENNE_TWISTER_STATE);
        }

    if (p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_TT800)
        {
            memcpy(buffer + position, ((TT800STATE *) p->rh.tt800.ptr) + state->alternative_rng_offset,
                   sizeof(TT800STATE));
            position += sizeof(TT800STATE);
        }
    return position;
}

int deserialize_rng_state(struct Phase *const p, RNG_STATE * const state, const unsigned char *const buffer)
{
    unsigned int position = 0;
    //default state
    memcpy(&(state->default_state), buffer + position, sizeof(PCG_STATE));
    position += sizeof(PCG_STATE);

    state->alternative_rng_offset = UINT64_MAX;

    if (p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_MT)
        {
            state->alternative_rng_offset = get_new_soma_memory_offset(&(p->rh.mt), 1);
            if (state->alternative_rng_offset == UINT64_MAX)
                {
                    fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__, __LINE__,
                            p->info_MPI.world_rank, p->n_polymers);
                    return -1;
                }
            memcpy(((MERSENNE_TWISTER_STATE *) p->rh.mt.ptr) + state->alternative_rng_offset, buffer + position,
                   sizeof(MERSENNE_TWISTER_STATE));
            position += sizeof(MERSENNE_TWISTER_STATE);
        }

    if (p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_TT800)
        {
            state->alternative_rng_offset = get_new_soma_memory_offset(&(p->rh.tt800), 1);
            if (state->alternative_rng_offset == UINT64_MAX)
                {
                    fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__, __LINE__,
                            p->info_MPI.world_rank, p->n_polymers);
                    return -1;
                }
            memcpy(((TT800STATE *) p->rh.tt800.ptr) + state->alternative_rng_offset, buffer + position,
                   sizeof(TT800STATE));
            position += sizeof(TT800STATE);
        }

    return position;
}

int seed_rng_state(struct RNG_STATE *const state, const unsigned int seed, const unsigned int stream,
                   const struct Phase *const p)
{
    soma_seed_rng(&(state->default_state), seed, stream);
    switch (p->args.pseudo_random_number_generator_arg)
        {
        case pseudo_random_number_generator__NULL:
            break;
        case pseudo_random_number_generator_arg_PCG32:
            break;
        case pseudo_random_number_generator_arg_MT:;
            MERSENNE_TWISTER_STATE *mt_state = p->rh.mt.ptr;
            soma_seed_rng_mt(&(state->default_state), mt_state + state->alternative_rng_offset);
            break;
        case pseudo_random_number_generator_arg_TT800:;
            TT800STATE *tt800_state = p->rh.tt800.ptr;
            soma_seed_rng_tt800(&(state->default_state), tt800_state + state->alternative_rng_offset);
            break;
        }                       //end switch
    return 0;
}
