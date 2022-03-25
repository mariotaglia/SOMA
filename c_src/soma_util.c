/* Copyright (C) 2016-2021 Ludwig Schneider

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
#include <string.h>

unsigned int get_bond_type(const uint32_t info)
{
    return info << 28 >> 29;
}

int get_offset(const int32_t info)
{
    return info >> 4;
}

unsigned int get_end(const uint32_t info)
{
    return info & 1;
}

uint32_t get_info(const int offset, const unsigned int bond_type, const unsigned int end)
{
    int ret = offset;
    ret <<= 3;
#ifndef _OPENACC
    assert(bond_type < 1 << 3);
#endif                          //_OPENACC
    ret |= bond_type;
    ret <<= 1;
#ifndef _OPENACC
    assert(end < 1 << 1);
#endif                          //_OPENACC
    ret |= end;
    return ret;
}

int get_bondlist_offset(const int32_t info_bl)
{
    return info_bl >> 8;
}

unsigned int get_particle_type(const struct Phase *const p, const uint64_t i, const unsigned int j)
{
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    return (unsigned int)((uint8_t *) p->ph.monomer_types.ptr)[p->polymers[i].monomer_type_offset + j];
#else                           //ENABLE_MONOTYPE_CONVERSIONS
    return get_particle_type_of_poly_arch(p->poly_arch[p->poly_type_offset[p->polymers[i].type] + 1 + j]);
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
}

unsigned int get_particle_type_of_poly_arch(const uint32_t info_bl)
{
    return info_bl << 24 >> 24;
}

uint32_t get_info_bl(const unsigned int offset_bl, const unsigned int type)
{
    unsigned int ret = offset_bl;
    ret <<= 8;
#ifndef _OPENACC
    assert(type < 1 << 8);
#endif                          //_OPENACC
    ret |= type;
    return ret;
}

int post_process_args(struct som_args *args, const unsigned int world_rank)
{
    if (args->timesteps_arg < 0)
        {
            if (world_rank == 0)
                fprintf(stderr, "WARNING: negative number of timesteps given. Is set to 0.\n");
            args->timesteps_arg = 0;
        }

    if (!args->ana_file_given)
        {
            if (world_rank == 0)
                fprintf(stderr, "WARNING: No ana-file specified.\n"
                        "WARNING: This causes that NO observables are going to be analysed during the whole run.\n"
                        "WARNING: This might be useful for timing runs or pure equilibration,\n"
                        "WARNING: but in the latter case we nonetheless highly recommend to specify an analyse file.\n");
        }
#ifdef OPENACC
    if ((!args->gpus_given && !args->only_gpu_given) || (args->gpus_given && args->gpus_arg <= 0))
        {
            if (world_rank == 0)
                fprintf(stderr,
                        "WARNING: No GPU usage specified, but this version has been compiled with OpenACC support.\n"
                        "WARNING: Are you sure, that you do not want to set a GPU? If not, try the options --gpus or --only-gpu.\n");
        }
    if (args->gpus_given && args->gpus_arg < 0)
        {
            if (world_rank == 0)
                fprintf(stderr, "WARNING: Negative number of gpus specified. Option will be ignored.\n");
            args->gpus_given = false;
        }
#endif
    if (args->screen_output_interval_arg < 0)
        {
            if (world_rank == 0)
                fprintf(stderr,
                        "WARNING: Negative number of seconds for screen ouput given. Screen ouput is switched off.\n");
            args->screen_output_interval_arg = 0;
        }

    if (args->autotuner_restart_period_arg < 0)
        {
            if (world_rank == 0)
                fprintf(stderr, "WARNING: Negative number for autotuner restart given, is switched off.\n");
            args->autotuner_restart_period_arg = 0;
        }

    if (args->load_balance_arg < 0)
        {
            if (world_rank == 0)
                fprintf(stderr, "WARNING: Negative number for load-balance freq given, is switched off.\n");
            args->load_balance_arg = 0;
        }
    if (args->accepted_load_inbalance_arg < 0.01)
        {
            if (world_rank == 0)
                fprintf(stderr, "WARNING: accepted load inbalance too low - using 0.01.\n");
            args->accepted_load_inbalance_arg = 0.01;
        }

    if (args->N_domains_arg < 1)
        {
            if (world_rank == 0)
                fprintf(stderr, "WARNING: Non positive number for domain decompostion given. Using 1 domain.\n");
            args->N_domains_arg = 1;
        }
    //Set to no domain decompostion defaults
    if (args->N_domains_arg == 1)
        {
            args->domain_buffer_arg = 0;
            args->rcm_update_arg = 0;
        }

    if (args->domain_buffer_arg < 0)
        {
            if (world_rank == 0)
                fprintf(stderr, "WARNING: Negative number for domain buffer given. Using 0 buffer cells.");
            args->domain_buffer_arg = 0;
        }

    if (args->N_domains_arg > 1 && args->rcm_update_arg < 1)
        {
            if (world_rank == 0)
                fprintf(stderr,
                        "WARNING: Non positive number for update of molecule center of mass given. Using frequency of 1.\n");
            args->rcm_update_arg = 1;
        }

#if (ENABLE_MIC != 1)
    if (args->bond_minimum_image_convention_flag)
        {
            fprintf(stderr,
                    "ERROR: %s:%d compilation without MIC support, but requested via CMDline. Request is ignored.\n",
                    __FILE__, __LINE__);
            return -2;
        }
#endif                          //ENABLE_MIC
#if ( ENABLE_DOMAIN_DECOMPOSITION != 1 || ENABLE_MPI != 1  )
    if (args->N_domains_arg > 1)
        {
            fprintf(stderr,
                    "ERROR: Domain decomposition  support requested via CMDline, but SOMA is compiled without the support /MPI. Recompile SOMA with the appropriate option for domain decomposition.\n");
            return -5;
        }
#endif                          //ENABLE_DOMAIN_DECOMPOSITION

    if (strstr(args->final_file_arg, ".h5") == NULL)
        {
            if (world_rank == 0)
                {
                    fprintf(stderr,
                            "WARNING: '.h5' could not be found in the filename of the final configuration.\nWARNING: Proceed at your own risk. filename=%s",
                            args->final_file_arg);
                }
        }

    return 0;
}

unsigned int get_number_bond_type(const struct Phase *const p, const enum Bondtype btype)
{
    unsigned int counter = 0;
    for (unsigned int p_type = 0; p_type < p->n_poly_type; p_type++)
        {
            const unsigned int N = p->poly_arch[p->poly_type_offset[p_type]];
            for (unsigned int mono = 0; mono < N; mono++)
                {
                    const int start = get_bondlist_offset(p->poly_arch[p->poly_type_offset[p_type] + 1 + mono]);
                    if (start > 0)
                        {
                            unsigned int end = 0;
                            for (int i = start; end == 0; i++)
                                {
                                    const uint32_t info = p->poly_arch[i];
                                    end = get_end(info);
                                    const unsigned int bond_type = get_bond_type(info);
                                    if (bond_type == btype)
                                        counter++;
                                }
                        }
                }
        }
    return counter;
}

int reseed(struct Phase *const p, const unsigned int seed)
{
    uint64_t n_polymer_offset = 0;
#if ( ENABLE_MPI == 1 )
    MPI_Scan(&(p->n_polymers), &n_polymer_offset, 1, MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);
    n_polymer_offset -= p->n_polymers;
#endif                          //ENABLE_MPI

    //Reset PRNG to initial state
    for (uint64_t i = 0; i < p->n_polymers; i++)
        {
            seed_rng_state(&(p->polymers[i].poly_state), seed, i + n_polymer_offset, p);
        }

    return 0;
}
