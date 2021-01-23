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

#ifndef SOMA_ALTERNATIVE_RNG_H
#define SOMA_ALTERNATIVE_RNG_H

struct Phase;
#include <stdint.h>
#include <stdbool.h>
#include "soma_config.h"
#include "soma_util.h"
#include "rng.h"
#include "soma_memory.h"

//! \brief Number of internal states of the Mersenne-Twister
#define MTMAX_num_int_state 624
//! \brief State of the random number generator (Mersenne-Twister)
typedef struct MERSENNE_TWISTER_STATE {
    uint32_t A[2];              //!<\brief "Static" mask for bitwise operations
    uint32_t internal_index;    //!<\brief Internal index of the Mersenne-Twister
    uint32_t internal_state[MTMAX_num_int_state];       //!<\brief Internal state of the Mersenne-Twister
} MERSENNE_TWISTER_STATE;

//! \brief Number of internal states of the TT800
#define TT_num_int_state 27
//! \brief State of the random number generator (TT800)
typedef struct TT800STATE {
    uint32_t A[2];              //!<\brief "Static" mask for bitwise operations
    uint32_t internal_index;    //!<\brief Internal index of the TT800
    uint32_t internal_state[TT_num_int_state];  //!<\brief Internal state of the TT800
} TT800STATE;

//! Struct to hold the global arrays of heany data for alternative PRNGs
//!
//! Each Phase has one of these structs. And in case an alternative RNG is used, this struct contains the memory for the states.
//! If the memory is not used the pointers are initialized to NULL
typedef struct RNG_HEAVY {
    SomaMemory mt;              //!< Array and meta data for MersenneTwister states
    SomaMemory tt800;           //!< Array abd meta data for TT800 states
} RNG_HEAVY;

//!\brief Set the seed of Mersenne-Twister with the PCG32
//!
//!\param rng
//!\param mt_rng
//!\return uint32
int soma_seed_rng_mt(PCG_STATE * rng, MERSENNE_TWISTER_STATE * mt_rng);
// !Mersenne Twister with state of 624 integers
// \return  as uint in range [0:soma_rng_uint_max_mt)

//!\brief Mersenne-Twister
//!
//!\param mt_rng is the struct which contains the internal state of the random number generator
//!\return uint32
#pragma acc routine(soma_mersenne_twister)
unsigned int soma_mersenne_twister(MERSENNE_TWISTER_STATE * mt_rng);
//! Status function to get the max random number.
//!
//! \return Maximum generated rng by soma_mersenne_twister()
#pragma acc routine(soma_rng_uint_mt)
unsigned int soma_rng_uint_max_mt();

//! Function initializes the internal state of thr reduced Mersenne-Twister TT800 with the PCG32
//! \param rng  struct which contains all information about PCG32
//! \param mt_rng  is the struct which contains the internal state of the random number generator
//! \return int
int soma_seed_rng_tt800(PCG_STATE * rng, TT800STATE * mt_rng);

//!\brief Function which uses the reduced Mersenne-Twister TT800
//!\param mt_rng is the struct which contains the internal state of the random number generator
//!\return uint32
#pragma acc routine(soma_rng_tt800) seq
unsigned int soma_rng_tt800(TT800STATE * mt_rng);

#endif                          //SOMA_RNG_ALTERNATIVE_H
