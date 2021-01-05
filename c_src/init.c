/* Copyright (C) 2016-2021 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016-2017 Marcel Langenberg
   Copyright (C) 2016 Fabien Leonforte
   Copyright (C) 2016 Juan Orozco
   Copyright (C) 2016 Yongzhi Ren
   Copyright (C) 2016 N. Harshavardhan Reddy

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

//! \file init.c
//! \brief Implementation of init.h

#include "init.h"
#include <math.h>
#include <time.h>
#include <assert.h>
#include <memory.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#if ( ENABLE_MPI == 1 )
#include <mpi.h>
#endif                          //ENABLE_MPI
#include <stdio.h>
#ifdef _OPENACC
#include <openacc.h>
#endif                          //_OPENACC
#ifdef _OPENMP
#include <omp.h>
#endif                          //_OPENMP
#include <hdf5.h>
#include "soma_config.h"
#include "phase.h"

int set_openacc_devices(struct Phase *const p)
{
    int ret = 0;
    p->info_MPI.gpu_id = -1;
#ifdef _OPENACC

    bool on_host = false;
    if (p->args.gpus_given && p->args.gpus_arg == 0)
        on_host = true;
    if ((!p->args.gpus_given && !p->args.only_gpu_given))
        on_host = true;
    if (on_host)
        {
            acc_init(acc_device_host);
            acc_set_device(acc_device_host);
            printf("INFO: rank %d runs single core CPU.\n", p->info_MPI.world_rank);
        }
    else
        {
            unsigned int my_gpu_rank;
            if (p->args.only_gpu_given)
                my_gpu_rank = p->args.only_gpu_arg;
            else
                {
                    assert(p->args.gpus_given);
                    assert(p->args.gpus_arg > 0);
                    my_gpu_rank = p->info_MPI.world_rank % p->args.gpus_arg;
                }

            acc_init(acc_device_nvidia);
            const unsigned int num_gpus = acc_get_num_devices(acc_device_nvidia);
            if (my_gpu_rank >= num_gpus)
                {
                    fprintf(stderr, "ERROR: rank %d tried to set gpuId %u, but only %u are available.\n",
                            p->info_MPI.world_rank, my_gpu_rank, num_gpus);
                    return -1;
                }
            acc_set_device_num(my_gpu_rank, acc_device_nvidia);
            const unsigned int check_gpu = acc_get_device_num(acc_device_nvidia);
            if (check_gpu != my_gpu_rank)
                {
                    fprintf(stderr, "WARNING: rank %d tried to run GPU %u, but now it runs GPU %u.\n",
                            p->info_MPI.world_rank, my_gpu_rank, check_gpu);
                }
            printf("INFO: rank %d runs GPU %u.\n", p->info_MPI.world_rank, my_gpu_rank);
            p->info_MPI.gpu_id = my_gpu_rank;
#ifdef ENABLE_NCCL
            ncclUniqueId id;
            if (p->info_MPI.world_rank == 0)
                ncclGetUniqueId(&id);
            MPI_Bcast(&id, sizeof(id), MPI_BYTE, 0, p->info_MPI.SOMA_comm_world);
            ncclCommInitRank(&(p->info_MPI.SOMA_nccl_world), p->info_MPI.world_size, id, p->info_MPI.world_rank);

            if (p->info_MPI.sim_rank == 0)
                ncclGetUniqueId(&id);
            MPI_Bcast(&id, sizeof(id), MPI_BYTE, 0, p->info_MPI.SOMA_comm_sim);
            ncclCommInitRank(&(p->info_MPI.SOMA_nccl_sim), p->info_MPI.sim_size, id, p->info_MPI.sim_rank);

            if (p->info_MPI.domain_rank == 0)
                ncclGetUniqueId(&id);
            MPI_Bcast(&id, sizeof(id), MPI_BYTE, 0, p->info_MPI.SOMA_comm_domain);
            ncclCommInitRank(&(p->info_MPI.SOMA_nccl_domain), p->info_MPI.domain_size, id, p->info_MPI.domain_rank);
#endif                          //ENABLE_MPI_CUDA

        }
    if (p->args.omp_threads_given && p->info_MPI.world_rank == 0)
        {
            fprintf(stderr, "ERROR: You passed an OMP option for an OPENACC compiled program. This"
                    "is not possible. If you want multiple CPU usage for OpenACC builds use"
                    "the multicore option, if the compiler supports it.\n");
            return -1;
        }
#else                           //_OPENACC
    if ((p->args.only_gpu_given || (p->args.gpus_given && p->args.gpus_arg != 0)) && p->info_MPI.world_rank == 0)
        {
            fprintf(stderr, "WARNING: The command line arguments request a GPU use,"
                    " but SOMA has been compiled without OpenACC support.\n"
                    "\t If you want GPU acceleration, recompile SOMA with OpenACC support.\n"
                    "\t This simulation will run on the CPU.\n");
            ret += 1;
        }
#ifdef _OPENMP
    omp_set_dynamic(0);
    //Fallback option set OMP ranks to 1
    unsigned int nthreads = 1;

    if (p->args.omp_threads_given && p->args.omp_threads_arg > 1)
        {
            nthreads = p->args.omp_threads_arg;
            if (p->args.omp_threads_arg > omp_get_num_procs())
                {
                    fprintf(stderr,
                            "WARNING: world rank %d tried to use %d OMP threads, but only %d are available. Resetting to max. number.\n",
                            p->info_MPI.world_rank, p->args.omp_threads_arg, omp_get_num_procs());
                    nthreads = omp_get_num_procs();
                }
        }
    omp_set_num_threads(nthreads);
    printf("INFO: world rank %d runs %u CPU OMP threads.\n", p->info_MPI.world_rank, nthreads);
#else
    if (p->args.omp_threads_given && p->args.omp_threads_arg > 1)
        {
            fprintf(stderr,
                    "WARNING: world rank %d tried to use %d OMP threads, but the binary is compiled without OMP support.\n",
                    p->info_MPI.world_rank, p->args.omp_threads_arg);
        }
#endif                          //_OPENMP

#endif                          //_OPENACC
    return ret;
}
