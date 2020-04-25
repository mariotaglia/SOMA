/* Copyright (C) 2016-2019 Ludwig Schneider
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
/*! Random number generator PCG32
  \param rng State for PRNG
  \return PRN
*/
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
            return (unsigned int)soma_mersenne_twister(p->rh.mt_state + state->alternative_rng_offset);
            break;
        case pseudo_random_number_generator_arg_TT800:
            return (unsigned int)soma_rng_tt800(p->rh.tt800_state + state->alternative_rng_offset);
            break;
        }
    return -1;
}

soma_scalar_t soma_rng_soma_scalar(RNG_STATE * rng, const struct Phase *const p)
{
    return ldexp(soma_rng_uint(rng, p), -32);
}

/*! generate 3D vector, 2 times Box-Mueller Transform, discards one value
  \param rng RNG State
  \param rng_type Type of the PRNG
  \param x result for X
  \param y result for Y
  \param z result for Z
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

/*! generate 3D vector, with a distribution that just has the 2nd and 4th moment of a gaussian
  \param rng RNG State
  \param rng_type Type of the PRNG
  \param x result for X
  \param y result for Y
  \param z result for Z
*/
void soma_normal_vector2(RNG_STATE * rng, const struct Phase *const p, soma_scalar_t * x,
                         soma_scalar_t * y, soma_scalar_t * z)
{
    // the two factors are connected to ensure a 2nd moment of 1:
    // sfactor = (3-3*qfactor**2/7.0)**0.5

    soma_scalar_t u1 = 2.0 * soma_rng_soma_scalar(rng, p) - 1.;
    soma_scalar_t u2 = 2.0 * soma_rng_soma_scalar(rng, p) - 1.;
    soma_scalar_t u3 = 2.0 * soma_rng_soma_scalar(rng, p) - 1.;
    soma_scalar_t u4 = 2.0 * soma_rng_soma_scalar(rng, p) - 1.;
    soma_scalar_t u5 = 2.0 * soma_rng_soma_scalar(rng, p) - 1.;
    soma_scalar_t u6 = 2.0 * soma_rng_soma_scalar(rng, p) - 1.;
    *x = 1.97212 * u1 * u1 * u1 + 1.1553052583624814 * u2;
    *y = 1.97212 * u3 * u3 * u3 + 1.1553052583624814 * u4;
    *z = 1.97212 * u5 * u5 * u5 + 1.1553052583624814 * u6;
}

unsigned int rng_state_serial_length(const struct Phase *const p)
{
    unsigned int length = 0;
    length += sizeof(PCG_STATE);
    if (p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_MT)
        length += sizeof(MERSENNE_TWISTER_STATE);
    if (p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_TT800)
        length += sizeof(MTTSTATE);
    return length;
}

int serialize_rng_state(const struct Phase *const p, const RNG_STATE * const state, unsigned char *const buffer)
{
    unsigned int position = 0;
    //default state
    memcpy(buffer + position, &(state->default_state), sizeof(PCG_STATE));
    position += sizeof(PCG_STATE);

    if (p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_MT)
        {
            memcpy(buffer + position, state->mt_state, sizeof(MERSENNE_TWISTER_STATE));
            position += sizeof(MERSENNE_TWISTER_STATE);
        }

    if (p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_TT800)
        {
            memcpy(buffer + position, state->tt800_state, sizeof(MTTSTATE));
            position += sizeof(MTTSTATE);
        }
    return position;
}

int deserialize_rng_state(const struct Phase *const p, RNG_STATE * const state, const unsigned char *const buffer)
{
    unsigned int position = 0;
    //default state
    memcpy(&(state->default_state), buffer + position, sizeof(PCG_STATE));
    position += sizeof(PCG_STATE);

    state->mt_state = NULL;
    state->tt800_state = NULL;

    if (p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_MT)
        {
            state->mt_state = (MERSENNE_TWISTER_STATE *) malloc(sizeof(MERSENNE_TWISTER_STATE));
            MALLOC_ERROR_CHECK(state->mt_state, sizeof(MERSENNE_TWISTER_STATE));

            memcpy(state->mt_state, buffer + position, sizeof(MERSENNE_TWISTER_STATE));
            position += sizeof(MERSENNE_TWISTER_STATE);
        }

    if (p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_TT800)
        {
            state->tt800_state = (MTTSTATE *) malloc(sizeof(MTTSTATE));
            MALLOC_ERROR_CHECK(state->tt800_state, sizeof(MTTSTATE));

            memcpy(state->tt800_state, buffer + position, sizeof(MTTSTATE));
            position += sizeof(MTTSTATE);
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
            soma_seed_rng_mt(&(state->default_state), p->rh.mt_state + state->alternative_rng_offset);
            break;
        case pseudo_random_number_generator_arg_TT800:;
            soma_seed_rng_tt800(&(state->default_state), p->rh.tt800_state + state->alternative_rng_offset);
            break;
        }                       //end switch
    return 0;
}
