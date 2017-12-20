/* Copyright (C) 2016-2017 Ludwig Schneider

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

//! \file soma_util.c
//! \brief Implementation of soma_util.h

#include "soma_util.h"
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include "cmdline.h"
#include "phase.h"
#include "mpiroutines.h"
#include <math.h>

unsigned int get_bond_type(const uint32_t info)
    {
    return info<<28>>29;
    }

int get_offset(const int32_t info)
    {
    return info>>4;
    }

unsigned int get_end(const uint32_t info)
    {
    return info & 1;
    }

uint32_t get_info(const int offset,const unsigned int bond_type,const unsigned int end)
    {
    int ret = offset;
    ret <<= 3;
#ifndef _OPENACC
    assert( bond_type < 1<<3);
#endif//_OPENACC
    ret |= bond_type;
    ret <<= 1;
#ifndef _OPENACC
    assert( end < 1<<1);
#endif//_OPENACC
    ret |= end;
    return ret;
    }

int get_bondlist_offset(const int32_t info_bl)
    {
    return info_bl>>8;
    }

unsigned int get_particle_type(const uint32_t info_bl)
    {
    return info_bl<<24>>24;
    }

uint32_t get_info_bl(const unsigned int offset_bl,const unsigned int type)
    {
    unsigned int ret = offset_bl;
    ret <<= 8;
#ifndef _OPENACC
    assert( type < 1<<8);
#endif//_OPENACC
    ret |= type;
    return ret;
    }

struct som_args post_process_args(struct som_args*args,const struct Info_MPI*const mpi_info)
    {
    if( mpi_info->current_core == 0)
	{
	cmdline_parser_print_version();
	printf("\n");
	}

    if( args->timesteps_arg < 0)
	{
	if(mpi_info->current_core == 0)
	    fprintf(stderr,"WARNING: negative number of timesteps given. Is set to 0.\n");
	args->timesteps_arg = 0;
	}

    if( ! args->ana_file_given )
	{
	if(mpi_info->current_core == 0)
	    fprintf(stderr,"WARNING: No ana-file specified.\n"
		    "WARNING: This causes that NO observables are going to be analysed during the whole run.\n"
		    "WARNING: This might be useful for timing runs or pure equilibration,\n"
		    "WARNING: but in the latter case we nonetheless highly recommend to specify an analyse file.\n");
	}
#ifdef OPENACC
    if( (!args->gpus_given && !args->only_gpu_given) || (args->gpus_given && args->gpus_arg <= 0) )
	{
	if(mpi_info->current_core == 0)
	    fprintf(stderr,"WARNING: No GPU usage specified, but this version has been compiled with OpenACC support.\n"
		    "WARNING: Are you sure, that you do not want to set a GPU? If not, try the options --gpus or --only-gpu.\n");
	}
    if( args->gpus_given && args->gpus_arg < 0 )
	{
	if(mpi_info->current_core == 0)
	    fprintf(stderr,"WARNING: Negative number of gpus specified. Option will be ignored.\n");
	args->gpus_given = false;
	}
#endif
    if( args->screen_output_interval_arg < 0)
	{
	if(mpi_info->current_core == 0)
	    fprintf(stderr,"WARNING: Negative number of seconds for screen ouput given. Screen ouput is switched off.\n");
	args->screen_output_interval_arg = 0;
	}

    if(args->autotuner_restart_period_arg < 0)
	{
	if( mpi_info->current_core == 0)
	    fprintf(stderr,"WARNING: Negative number for autotuner restart given, is switched off.\n");
	args->autotuner_restart_period_arg = 0;
	}

    if(args->load_balance_arg < 0)
	{
	if( mpi_info->current_core == 0)
	    fprintf(stderr,"WARNING: Negative number for load-balance freq given, is switched off.\n");
	args->load_balance_arg = 0;
	}

    uint32_t fixed_seed;
    if( ! args->rng_seed_given || args->rng_seed_arg < 0)
	fixed_seed = time(NULL);
    else
	fixed_seed = args->rng_seed_arg;

    MPI_Bcast(&fixed_seed,1,MPI_UINT32_T,0,mpi_info->SOMA_MPI_Comm);
    args->rng_seed_arg = fixed_seed;
    if(mpi_info->current_core == 0)
	printf("All %d ranks use fixed seed %u.\n",mpi_info->Ncores,args->rng_seed_arg);

#if (ENABLE_MIC != 1)
    if(args->bond_minimum_image_convention_flag)
	{
	fprintf(stderr, "WARNING: %s:%d compilation without MIC support, but requested via CMDline. Request is ignored.\n", __FILE__, __LINE__);
	}
#endif//ENABLE_MIC

    if(mpi_info->current_core == 0)
	cmdline_parser_dump(stdout, args);
    return *args;
    }

unsigned int get_number_bond_type(const struct Phase*const p,const enum Bondtype btype)
    {
    unsigned int counter = 0;
    for(unsigned int p_type=0; p_type < p->n_poly_type; p_type++)
	{
	const unsigned int N = p->poly_arch[ p->poly_type_offset[ p_type ] ];
	for(unsigned int mono = 0 ; mono < N; mono++)
	    {
	    const int start = get_bondlist_offset(
		p->poly_arch[ p->poly_type_offset[ p_type] + 1 + mono ]);
	    if( start > 0)
		{
		int i=start;
		unsigned int end;
		do{
		    const uint32_t info = p->poly_arch[i++];
		    end = get_end(info);
		    const unsigned int bond_type = get_bond_type(info);
		    if( bond_type == btype )
			counter++;
		    }while(end==0);
		}
	    }
	}
    return counter;
    }

int reseed(struct Phase*const p,const unsigned int seed)
    {
    uint64_t n_polymer_offset;
    MPI_Scan( &(p->n_polymers), &n_polymer_offset, 1,MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_MPI_Comm);
    n_polymer_offset -= p->n_polymers;

    //Reset PRNG to initial state
    for(uint64_t i=0; i < p->n_polymers;i++)
	{
	seed_rng_state(&(p->polymers[i].poly_state), seed,
		       i+n_polymer_offset, p->args.pseudo_random_number_generator_arg);
	update_device_rng_state(&(p->polymers[i].poly_state), p->args.pseudo_random_number_generator_arg);
	}
    return 0;
    }
