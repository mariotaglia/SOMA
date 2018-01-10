//The code for the PCG random number generation is derived work from
//the original PCG software "http://www.pcg-random.org/" the license
//is Apache version 2. A license text is found in the file
// Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>
//"PCG_LICENSE.txt"

/* Copyright (C) 2016-2018 Ludwig Schneider
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
#ifndef SOMA_RNG_H
#define SOMA_RNG_H
#include <stdint.h>
#include "cmdline.h"
struct Phase;
struct RNG_STATE;

/*! \file rng.h
  \brief Definition of pseudo random number generation wrappers for soma.
*/



//! \brief State of the random number generator (PCG)
typedef struct{
    uint64_t state; //!<\brief internal state of the PCG generator
    uint64_t inc;   //!<\brief stream of the PCG generator
    }PCG_STATE;

//! \brief Number of internal states of the Mersenne-Twister
#define MTMAX_num_int_state 624
//! \brief State of the random number generator (Mersenne-Twister)
typedef struct MERSENNE_TWISTER_STATE{
    uint32_t  A[2];                                //!<\brief "Static" mask for bitwise operations
    uint32_t  internal_index;                      //!<\brief Internal index of the Mersenne-Twister
    uint32_t  internal_state[MTMAX_num_int_state]; //!<\brief Internal state of the Mersenne-Twister
    }MERSENNE_TWISTER_STATE;

//! \brief Number of internal states of the TT800
#define TT_num_int_state 27
//! \brief State of the random number generator (TT800)
typedef struct MTTSTATE{
    uint32_t  A[2];                             //!<\brief "Static" mask for bitwise operations
    uint32_t  internal_index;                   //!<\brief Internal index of the TT800
    uint32_t  internal_state[TT_num_int_state]; //!<\brief Internal state of the TT800
    }MTTSTATE;

//! \brief Struct which contains the random number generators.
//!
//!Option to select the pseudo random number  enerator  (possible values="PCG32", "MT" default=`PCG32')
//!if you prefer a different rng add it here. Modify the init.*, rng.*, and add an enum option to soma.ggo
typedef struct RNG_STATE {
    PCG_STATE              default_state;   //!<\brief PCG32
    MERSENNE_TWISTER_STATE *mt_state;       //!<\brief Mersenne-Twister
    MTTSTATE                 *tt800_state;    //!<\brief Reduced Mersenne-Twister TT800
    }RNG_STATE;

//! Wrapper for seeding the global random number generation.
//!
//! \param stream Stream of the RNG.
//! \param rng PCG_STATE to
//! \param seed Seed for the global rng.
//! \return Error code.
int soma_seed_rng(PCG_STATE * rng, uint64_t seed, uint64_t stream);

//! Wrapper for any prng we use for soma.
//!
//! \pre rng has been seeded.
//! \param state RNG_STATE to use and modify for PRNG
//! \param rng_type Type of the used PRNG
//! \return prng as uint in range [0:soma_rng_uint_max)
#pragma acc routine(soma_rng_uint) seq
unsigned int soma_rng_uint(RNG_STATE * state, enum enum_pseudo_random_number_generator rng_type);

//! Status function to get the max random number.
//!
//! \return Maximum generated rng by soma_rng_uint
unsigned int soma_rng_uint_max(void);

//!\brief Set the seed of Mersenne-Twister with the PCG32
//!
//!\param rng
//!\param mt_rng
//!\return uint32
#pragma acc routine(soma_seed_rng_mt) seq
int soma_seed_rng_mt(PCG_STATE *rng, MERSENNE_TWISTER_STATE *mt_rng);
// !Mersenne Twister with state of 624 integers
// \return  as uint in range [0:soma_rng_uint_max_mt)


//!\brief Mersenne-Twister
//!
//!\param mt_rng is the struct which contains the internal state of the random number generator
//!\return uint32
#pragma acc routine(soma_mersenne_twister)
unsigned int soma_mersenne_twister(MERSENNE_TWISTER_STATE *mt_rng);
//! Status function to get the max random number.
//!
//! \return Maximum generated rng by soma_mersenne_twister()
#pragma acc routine(soma_rng_uint_mt)
unsigned int soma_rng_uint_max_mt();

//! Wrapper function for float random numbers.
//! \param rng_type enum which carries information about the selected random number generator
//! \param rng struct which contains all information about the internal states of the rngs
//! \pre rng has been seeded.
//! \return prng in range [0,1)
#pragma acc routine(soma_rng_soma_scalar) seq
soma_scalar_t  soma_rng_soma_scalar(RNG_STATE * rng, enum enum_pseudo_random_number_generator rng_type);

//! Function that adds a 3D gaussian vector to the vector (x,y,z)
//! \param rng struct which contains all information about the internal states of the rngs
//! \param rng_type enum which carries information about the selected random number generator
//! \param x coordinate of the vector
//! \param y coordinate of the vector
//! \param z coordinate of the vector
//! \pre rng has been seeded
#pragma acc routine(soma_normal_vector) seq
void soma_normal_vector(RNG_STATE * rng, enum enum_pseudo_random_number_generator rng_type, soma_scalar_t *x, soma_scalar_t *y, soma_scalar_t *z);


//! Function that generates 3D vector (x,y,z), with a distribution that just has the 2nd and 4th moment of a gaussian
//! \param rng struct which contains all information about the internal states of the rngs
//! \param rng_type enum which carries information about the selected random number generator
//! \param x coordinate of the vector
//! \param y coordinate of the vector
//! \param z coordinate of the vector
//! \pre rng has been seeded
void soma_normal_vector2(RNG_STATE * rng,enum enum_pseudo_random_number_generator rng_type, soma_scalar_t *x, soma_scalar_t *y, soma_scalar_t *z);

//! Function initializes the internal state of thr reduced Mersenne-Twister TT800 with the PCG32
//! \param rng  struct which contains all information about PCG32
//! \param mt_rng  is the struct which contains the internal state of the random number generator
//! \return int
#pragma acc routine(soma_seed_rng_tt800) seq
int soma_seed_rng_tt800(PCG_STATE *rng, MTTSTATE *mt_rng);


//!\brief Function which uses the reduced Mersenne-Twister TT800
//!\param mt_rng is the struct which contains the internal state of the random number generator
//!\return uint32
#pragma acc routine(soma_rng_tt800) seq
unsigned int soma_rng_tt800(MTTSTATE *mt_rng);

//! Obtain the number of bytes, which are necessary to serialize a RNG_STATE.
//!
//! the current system configuration might influence the result. Especially,
//! special RNGs. (Deep copy included.)
//! \param p System configuration.
//! \return Number of bytes.
unsigned int rng_state_serial_length(const struct Phase*const p);

//! Serialize an RNG_STATE to a raw memory buffer.
//!
//! \param p System.
//! \param state State to serialize.
//! \param buffer Preallocated buffer to store the outcome.
//! \pre Allocation of buffer with return value of rng_state_serial_length() minimum.
//! \note Ownership and allocation status is unchanged.
//! \return Number of written bytes. If < 0 Errorcode.
int serialize_rng_state(const struct Phase*const p,const RNG_STATE*const state, unsigned char*const buffer);

//! Deserialize an RNG_STATE from a raw memory buffer.
//!
//! \param p System.
//! \param state State to initialize by memory buffer.
//! \param buffer Initialized memory buffer to read.
//! \pre You are owner of \a state. And there is no deep
//! copy data allocated. Otherwise, you create memory leaks.
//! \post You are owner of the state including deep copy data, because deep copy data is allocated.
//! \return Number of written bytes. If < 0 Errorcode.
int deserialize_rng_state(const struct Phase*const p,RNG_STATE*const state, const unsigned char*const buffer);

//! Function to initialize an RNG_STATE
//!
//! \param state State to initialize
//! \param seed Seed used for initialization
//! \param stream Stream used for initialization
//! \param rng_type Type of the RNG to use
//! \return Errorcode
int init_rng_state(struct RNG_STATE*const state,const unsigned int seed, const unsigned int stream,const enum enum_pseudo_random_number_generator rng_type);

//! Function to allocate memory for the RNG_STATE if necessary, and enter data for device
//!
//! \param state RNG_STATE to init
//! \param rng_type Type of the PRNG in use
//! \return Errorcode
int allocate_rng_state(struct RNG_STATE*const state,const enum enum_pseudo_random_number_generator rng_type);

//! Function to enter copyin memory for the RNG_STATE
//!
//! \param state RNG_STATE to copyin
//! \param rng_type Type of the PRNG in use
//! \return Errorcode
int copyin_rng_state(struct RNG_STATE*const state,const enum enum_pseudo_random_number_generator rng_type);

//! Free a memory of the RNG_STATE, and exit data for device
//!
//! \param state RNG_STATE to free
//! \param rng_type Type of the PRNG in use
//! \return Errorcode
int deallocate_rng_state(struct RNG_STATE*const state,const enum enum_pseudo_random_number_generator rng_type);

//! Function to exit copyout memory for the RNG_STATE
//!
//! \param state RNG_STATE to copyout
//! \param rng_type Type of the PRNG in use
//! \return Errorcode
int copyout_rng_state(struct RNG_STATE*const state,const enum enum_pseudo_random_number_generator rng_type);

//! Seed the RNG state
//!
//! \param seed New seed of the RNG
//! \param stream Stream for the PCG32 generator
//! \param state RNG_STATE to seed
//! \param rng_type Type of the PRNG in use
//! \return Errorcode
int seed_rng_state(struct RNG_STATE*const state,const unsigned int seed, const unsigned int stream,const enum enum_pseudo_random_number_generator rng_type);

//! Update the RNG state to the Device
//!
//! \param state RNG_STATE to update
//! \param rng_type Type of the PRNG in use
//! \return Errorcode
int update_device_rng_state(struct RNG_STATE*const state,const enum enum_pseudo_random_number_generator rng_type);

//! Update the RNG state to the Self
//!
//! \param state RNG_STATE to update
//! \param rng_type Type of the PRNG in use
//! \return Errorcode
int update_self_rng_state(struct RNG_STATE*const state,const enum enum_pseudo_random_number_generator rng_type);

#endif//SOMA_RNG_H
