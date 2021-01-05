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

//! \file rng_alternative.c
//! \brief Implementation of rng_alternative.h

/*Random number generator Mersenne-Twister*/

#include <stdlib.h>
#include <assert.h>
#include "phase.h"
#include "rng_alternative.h"
#include "rng.h"

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

int soma_seed_rng_tt800(PCG_STATE * rng, TT800STATE * tt800_rng)
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
unsigned int soma_rng_tt800(TT800STATE * itt800_rng)
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
