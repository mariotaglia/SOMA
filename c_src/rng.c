/* Copyright (C) 2016-2018 Ludwig Schneider
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
    rng->inc = stream*2 + 1;
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

/*Random number generator Mersenne-Twister*/
int soma_seed_rng_mt(PCG_STATE *rng, MERSENNE_TWISTER_STATE *mt_rng)
    {
    mt_rng->internal_index = MTMAX_num_int_state +1;
    mt_rng->A[0] = 0;
    mt_rng->A[1] = 0x9908b0df;
    for (int i=0;i< MTMAX_num_int_state;i++){
	mt_rng->internal_state[i] = pcg32_random(rng);
	}
    return 0;
    }

#pragma acc routine(soma_mersenne_twister) seq
unsigned int soma_mersenne_twister(MERSENNE_TWISTER_STATE *mt_rng)
    {

    unsigned int M = 397;
    uint32_t HI = 0x80000000;
    uint32_t LO = 0x7fffffff;
    uint32_t  e;

    /* The Mersenne-Twister stete of 624 is seeded with soma_seed_rng_mt()
     * which is called in the init.c by function init_values()
     */
    if (M>MTMAX_num_int_state) M= MTMAX_num_int_state/2;

    if (mt_rng->internal_index >= MTMAX_num_int_state) {
	/* Berechne neuen Zustandsvektor */
	uint32_t h;

	for (unsigned int i=0; i<MTMAX_num_int_state-M; ++i) {
	    h = (mt_rng->internal_state[i] & HI) | (mt_rng->internal_state[i+1] & LO); // Crashes HERE!!!
	    mt_rng->internal_state[i] = mt_rng->internal_state[i+M] ^ (h >> 1) ^ mt_rng->A[h & 1];
	    }
	for (unsigned int i=MTMAX_num_int_state-M; i<MTMAX_num_int_state-1; ++i) {
	    h = (mt_rng->internal_state[i] & HI) | (mt_rng->internal_state[i+1] & LO);
	    mt_rng->internal_state[i] = mt_rng->internal_state[ i + (M-MTMAX_num_int_state) ] ^ (h >> 1) ^ mt_rng->A[h & 1];
	    }

	h = (mt_rng->internal_state[MTMAX_num_int_state-1] & HI) | (mt_rng->internal_state[0] & LO);
	mt_rng->internal_state[MTMAX_num_int_state-1] = mt_rng->internal_state[M-1] ^ (h >> 1) ^ mt_rng->A[h & 1];
	mt_rng->internal_index = 0;
	}

    e = mt_rng->internal_state[mt_rng->internal_index++];
    /* Tempering */
    e ^= (e >> 11);
    e ^= (e <<  7) & 0x9d2c5680;
    e ^= (e << 15) & 0xefc60000;
    e ^= (e >> 18);

    return e;
    }

unsigned int soma_rng_uint_max_mt()
    {
    return 0x80000000;
    }

unsigned int soma_rng_uint(RNG_STATE * state, enum enum_pseudo_random_number_generator rng_type)
    {
    switch(rng_type){
    case pseudo_random_number_generator__NULL     : return (unsigned int) pcg32_random(&(state->default_state));  break;
    case pseudo_random_number_generator_arg_PCG32 : return (unsigned int) pcg32_random(&(state->default_state));  break;
    case pseudo_random_number_generator_arg_MT    : return (unsigned int) soma_mersenne_twister(state->mt_state); break;
    case pseudo_random_number_generator_arg_TT800 : return (unsigned int) soma_rng_tt800(state->tt800_state); break;
        }
    return -1;
    }


soma_scalar_t soma_rng_soma_scalar(RNG_STATE * rng, enum enum_pseudo_random_number_generator rng_type)
    {
    return ldexp(soma_rng_uint(rng, rng_type), -32);
    }

/*! generate 3D vector, 2 times Box-Mueller Transform, discards one value
  \param rng RNG State
  \param rng_type Type of the PRNG
  \param x result for X
  \param y result for Y
  \param z result for Z
*/
void soma_normal_vector(RNG_STATE * rng, enum enum_pseudo_random_number_generator rng_type,  soma_scalar_t *x, soma_scalar_t *y, soma_scalar_t *z)
    {
    soma_scalar_t u1, u2, u3, u4, r1, r2;

    u1 = 2*soma_rng_soma_scalar(rng, rng_type )-1.;
    u2 = 2*soma_rng_soma_scalar(rng, rng_type )-1.;
    u3 = 2*soma_rng_soma_scalar(rng, rng_type)-1.;
    u4 = 2*soma_rng_soma_scalar(rng, rng_type)-1.;

    r1  = u1*u1+u2*u2;
    r2  = u3*u3+u4*u4;

    while(r1>1){
	u1 = 2*soma_rng_soma_scalar(rng, rng_type)-1.;
	u2 = 2*soma_rng_soma_scalar(rng, rng_type)-1.;
	r1  = u1*u1+u2*u2;
	}

    while(r2>1){
	u3 = 2*soma_rng_soma_scalar(rng, rng_type)-1.;
	u4 = 2*soma_rng_soma_scalar(rng, rng_type)-1.;
	r2  = u3*u3+u4*u4;
	}

    const soma_scalar_t root1  = sqrt( -2.0 * log(r1) / r1 );
    const soma_scalar_t root2  = sqrt( -2.0 * log(r2) / r2 );

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
void soma_normal_vector2(RNG_STATE * rng, enum enum_pseudo_random_number_generator rng_type,   soma_scalar_t *x, soma_scalar_t *y, soma_scalar_t *z)
    {
    // the two factors are connected to ensure a 2nd moment of 1:
    // sfactor = (3-3*qfactor**2/7.0)**0.5

    soma_scalar_t u1 = 2.0*soma_rng_soma_scalar(rng, rng_type )-1.;
    soma_scalar_t u2 = 2.0*soma_rng_soma_scalar(rng, rng_type )-1.;
    soma_scalar_t u3 = 2.0*soma_rng_soma_scalar(rng, rng_type )-1.;
    soma_scalar_t u4 = 2.0*soma_rng_soma_scalar(rng, rng_type )-1.;
    soma_scalar_t u5 = 2.0*soma_rng_soma_scalar(rng, rng_type)-1.;
    soma_scalar_t u6 = 2.0*soma_rng_soma_scalar(rng, rng_type )-1.;
    *x = 1.97212*u1*u1*u1 + 1.1553052583624814 * u2;
    *y = 1.97212*u3*u3*u3 + 1.1553052583624814 * u4;
    *z = 1.97212*u5*u5*u5 + 1.1553052583624814 * u6;
    }

int soma_seed_rng_tt800(PCG_STATE *rng, MTTSTATE *tt800_rng){

    tt800_rng->internal_index = MTMAX_num_int_state +1;
    tt800_rng->A[0] = 0;
    tt800_rng->A[1] = 0x8ebfd028;

    for (int k=0;k<TT_num_int_state;k++){
	tt800_rng->internal_state[k] = (uint32_t) pcg32_random(rng);
	}
    return 0;
    }
#pragma acc routine(soma_tt800) seq
unsigned int soma_rng_tt800(MTTSTATE *itt800_rng){

    uint32_t M=7;
    uint32_t e;
    uint32_t k=0;
    if (itt800_rng->internal_index >= TT_num_int_state) {

	if (itt800_rng->internal_index > TT_num_int_state) {
	    uint32_t r = 9;
	    uint32_t s = 3402;
	    for (k=0; k<TT_num_int_state; ++k) {
		r = 509845221 * r + 3;
		s *= s + 1;
		itt800_rng->internal_state[k] = s + (r >> 10);
		}
	    }
	for (k=0; k<TT_num_int_state-M; ++k){
	    itt800_rng->internal_state[k] = itt800_rng->internal_state[k+M] ^ (itt800_rng->internal_state[k] >> 1) ^ itt800_rng->A[itt800_rng->internal_state[k] & 1];
	    }
	for (k=TT_num_int_state-M; k<TT_num_int_state; ++k){
	    itt800_rng->internal_state[k] = itt800_rng->internal_state[k+(M-TT_num_int_state)] ^ (itt800_rng->internal_state[k] >> 1) ^ itt800_rng->A[itt800_rng->internal_state[k] & 1];
	    }
	itt800_rng->internal_index = 0;
	}

    e = itt800_rng->internal_state[itt800_rng->internal_index++];
    /* Tempering */
    e ^= (e <<  7) & 0x2b5b2500;
    e ^= (e << 15) & 0xdb8b0000;
    e ^= (e >> 16);
    return e;
    }

unsigned int rng_state_serial_length(const struct Phase*const p)
    {
    unsigned int length = 0;
    length += sizeof(PCG_STATE);
    if( p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_MT )
	length += sizeof(MERSENNE_TWISTER_STATE);
    if( p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_TT800 )
	length += sizeof(MTTSTATE);
    return length;
    }

int serialize_rng_state(const struct Phase*const p,const RNG_STATE*const state, unsigned char*const buffer)
    {
    unsigned int position = 0;
    //default state
    memcpy( buffer +position , &(state->default_state), sizeof(PCG_STATE) );
    position += sizeof(PCG_STATE);

    if( p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_MT)
	{
	memcpy( buffer+position , state->mt_state, sizeof(MERSENNE_TWISTER_STATE) );
	position += sizeof(MERSENNE_TWISTER_STATE);
	}

    if(p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_TT800)
	{
	memcpy(buffer + position , state->tt800_state, sizeof(MTTSTATE));
	position += sizeof(MTTSTATE);
	}
    return position;
    }

int deserialize_rng_state(const struct Phase*const p,RNG_STATE*const state, const unsigned char*const buffer)
   {
    unsigned int position = 0;
    //default state
    memcpy( &(state->default_state),buffer +position , sizeof(PCG_STATE) );
    position += sizeof(PCG_STATE);

    state->mt_state = NULL;
    state->tt800_state = NULL;

    if( p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_MT)
	{
	state->mt_state = (MERSENNE_TWISTER_STATE*) malloc(sizeof(MERSENNE_TWISTER_STATE));
	MALLOC_ERROR_CHECK(state->mt_state, sizeof(MERSENNE_TWISTER_STATE));

	memcpy( state->mt_state, buffer+position , sizeof(MERSENNE_TWISTER_STATE) );
	position += sizeof(MERSENNE_TWISTER_STATE);
	}

    if(p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_TT800)
	{
	state->tt800_state = (MTTSTATE*)malloc(sizeof(MTTSTATE));
	MALLOC_ERROR_CHECK(state->tt800_state, sizeof(MTTSTATE));

	memcpy(state->tt800_state,buffer + position , sizeof(MTTSTATE));
	position += sizeof(MTTSTATE);
	}

    return position;
    }

int init_rng_state(struct RNG_STATE*const state,const unsigned int seed, const unsigned int stream,const enum enum_pseudo_random_number_generator rng_type)
    {
    allocate_rng_state(state, rng_type);
    seed_rng_state(state, seed, stream, rng_type);
    copyin_rng_state(state, rng_type);
    return 0;
    }

int allocate_rng_state(struct RNG_STATE*const state,const enum enum_pseudo_random_number_generator rng_type)
    {
    state->mt_state = NULL;
    state->tt800_state = NULL;
    switch( rng_type )
	{
	case pseudo_random_number_generator__NULL: break;
	case pseudo_random_number_generator_arg_PCG32: break;
	case pseudo_random_number_generator_arg_MT: ;
	    MERSENNE_TWISTER_STATE *mt_state_tmp = (MERSENNE_TWISTER_STATE*)malloc(sizeof(MERSENNE_TWISTER_STATE));
	    if(mt_state_tmp == NULL){
	    fprintf(stderr, "ERROR: By malloc MERSENNE_TWISTER_STATE , %s %d ",
		    __FILE__, __LINE__ );
	    return -1;
	    }
	    state->mt_state = mt_state_tmp;
	    break;
	case pseudo_random_number_generator_arg_TT800: ;
	    MTTSTATE *tt800_state_tmp = (MTTSTATE*)malloc(sizeof(MTTSTATE));
	    if(tt800_state_tmp == NULL){
		fprintf(stderr, "ERROR: By malloc TT800 %s %d",__FILE__, __LINE__ );
		return -1;
		}
	    state->tt800_state = tt800_state_tmp;
	    break;
	}//end switch
    return 0;
    }

int deallocate_rng_state(struct RNG_STATE*const state,const enum enum_pseudo_random_number_generator rng_type)
    {
    switch( rng_type )
	{
	case pseudo_random_number_generator__NULL: break;
	case pseudo_random_number_generator_arg_PCG32: break;
	case pseudo_random_number_generator_arg_MT: ;
	    free(state->mt_state);
	    break;
	case pseudo_random_number_generator_arg_TT800: ;
	    free(state->tt800_state);
	    break;
	}//end switch
    return 0;
    }

int copyin_rng_state(struct RNG_STATE*const state,const enum enum_pseudo_random_number_generator rng_type)
    {
    switch( rng_type )
	{
	case pseudo_random_number_generator__NULL: break;
	case pseudo_random_number_generator_arg_PCG32: break;
	case pseudo_random_number_generator_arg_MT: ;
#pragma acc enter data copyin(state->mt_state[0:1])
	    break;
	case pseudo_random_number_generator_arg_TT800: ;
#pragma acc enter data copyin(state->tt800_state[0:1])
	    break;
	}//end switch
    return 0*state->default_state.state;
    }

int copyout_rng_state(struct RNG_STATE*const state,const enum enum_pseudo_random_number_generator rng_type)
    {
    switch( rng_type )
	{
	case pseudo_random_number_generator__NULL: break;
	case pseudo_random_number_generator_arg_PCG32: break;
	case pseudo_random_number_generator_arg_MT: ;
#pragma acc exit data copyout(state->mt_state[0:1])
	    break;
	case pseudo_random_number_generator_arg_TT800: ;
#pragma acc exit data copyout(state->tt800_state[0:1])
	    break;
	}//end switch
    return 0*state->default_state.state;
    }

int seed_rng_state(struct RNG_STATE*const state,const unsigned int seed, const unsigned int stream,const enum enum_pseudo_random_number_generator rng_type)
    {
    soma_seed_rng(&(state->default_state),seed,stream);
    switch( rng_type )
	{
	case pseudo_random_number_generator__NULL: break;
	case pseudo_random_number_generator_arg_PCG32: break;
	case pseudo_random_number_generator_arg_MT: ;
	    soma_seed_rng_mt(&(state->default_state),state->mt_state);
	    break;
	case pseudo_random_number_generator_arg_TT800: ;
	    soma_seed_rng_tt800(&(state->default_state),state->tt800_state);
	    break;
	}//end switch
    return 0;
    }

int update_device_rng_state(struct RNG_STATE*const state,const enum enum_pseudo_random_number_generator rng_type)
    {
#pragma acc update device(state->default_state[0:1])
    switch( rng_type )
	{
	case pseudo_random_number_generator__NULL: break;
	case pseudo_random_number_generator_arg_PCG32: break;
	case pseudo_random_number_generator_arg_MT: ;
#pragma acc update device(state->mt_state[0:1])
	    break;
	case pseudo_random_number_generator_arg_TT800: ;
#pragma acc update device(state->tt800_state)
	    break;
	}//end switch
    return 0*state->default_state.state;
    }

int update_self_rng_state(struct RNG_STATE*const state,const enum enum_pseudo_random_number_generator rng_type)
    {
#pragma acc update self(state->default_state[0:1])
    switch( rng_type )
	{
	case pseudo_random_number_generator__NULL: break;
	case pseudo_random_number_generator_arg_PCG32: break;
	case pseudo_random_number_generator_arg_MT: ;
#pragma acc update self(state->mt_state[0:1])
	    break;
	case pseudo_random_number_generator_arg_TT800: ;
#pragma acc update self(state->tt800_state)
	    break;
	}//end switch
    return 0*state->default_state.state;
    }
