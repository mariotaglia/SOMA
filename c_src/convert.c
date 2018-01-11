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

/* SOMA  */
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <stdlib.h>
#include "soma_util.h"
#include "phase.h"
#include "io.h"
#include "mpiroutines.h"
#include "init.h"
#include "test.h"

//! \file convert.c
//! \brief Implementation of CONVERT executable.

//! Main function of convert.c
//! \private
//!
//! \param argc Argument counter
//! \param argv Argument vector
//! \return Errorcode
int main(int argc, char *argv[])
    {
    if (MPI_Init(NULL, NULL) != MPI_SUCCESS)
        {
        fprintf(stderr, "MPI_ERROR (start)\n");
        return -1;
        }

    Phase phase;
    Phase *const p = &phase;
   /* initialize MPI */
    const int error = MPI_Comm_dup( MPI_COMM_WORLD, &(p->info_MPI.SOMA_comm_world) );
    if( error != MPI_SUCCESS)
        {
        fprintf(stderr,"MPI Error cannot duplicate MPI_COMM_WORLD %s:%d\n",__FILE__,__LINE__);
        return -1;
        }
    const int error2 = MPI_Comm_dup( p->info_MPI.SOMA_comm_world, &(p->info_MPI.SOMA_comm_sim) );
    if( error2 != MPI_SUCCESS)
        {
        fprintf(stderr,"MPI Error cannot duplicate MPI_COMM_WORLD %s:%d\n",__FILE__,__LINE__);
        return -1;
        }
    p->args.N_domains_arg = 1;
    p->args.domain_buffer_arg = 0;

    init_MPI(p);

    if ( !(argc == 2 || argc ==3) ) {
        if (p->info_MPI.world_rank == 0)
            fprintf(stderr, "Usage: %s filename-to-convert-to-hdf5 <geometry-file>\n",
                    argv[0]);
        finalize_MPI(&(p->info_MPI));
        return 0;               //mpi restrictions -- no errorcode
        }
    if (p->info_MPI.world_size != 1) {
        if (p->info_MPI.world_rank == 0)
            fprintf(stderr, "Use the convert tool only with 1 core.%d \n",
                    p->info_MPI.world_size);
        finalize_MPI(&(p->info_MPI));
        return 0;               //mpi restrictions -- no errorcode
        }

    //Set the args in a way, that the host is used.
    p->args.gpus_given = false;
    p->args.only_gpu_given = false;
    p->args.gpus_arg=0;
    p->args.omp_threads_given = false;
    p->args.omp_threads_arg = 1;
    p->args.coord_file_arg = NULL;
    const int devices = set_openacc_devices(p);
    MPI_ERROR_CHECK(devices,"Cannot set openacc devices");

    /* read in configuration with routine from io */
    const int read = read_old_config(p, argv[1]);
    MPI_ERROR_CHECK(read,"Unsucessfully read a configuration. Try to correct your input file.");

    if(argc > 2)
        {
        const int geometry = read_old_geometry(p, argv[2]);
        MPI_ERROR_CHECK(geometry,"Unsucessfully read the geometry. Try to correct your input file.");
        }

    /* initialize phase */
    const int init_v = init_phase(p);
    MPI_ERROR_CHECK(init_v,"Init values faild.");
    const int ana = init_ana(p, NULL,NULL);
    MPI_ERROR_CHECK(ana,"Init ana failed.");

    const int test_pt = test_particle_types(p);
    MPI_ERROR_CHECK(test_pt,"test_particle_type failed.");

    const int test_51 = test_area51_violation(p);
    MPI_ERROR_CHECK(test_51,"Area51 failed.");

    //String magic to get the hdf5 filename
    char *filename;
    const unsigned int new_len = (strlen(argv[1]) + 4);
    filename = (char *) malloc(new_len * sizeof(char));
    if (filename == NULL) {
        printf("Malloc error\n");
        finalize_MPI(&(p->info_MPI));
        return 0;
        }
    memset(filename, '\0', new_len * sizeof(char));
    strcpy(filename, argv[1]);
    unsigned int first_dot;
    for (first_dot = new_len - 1; first_dot != 0; first_dot--)
        if (filename[first_dot] == '.')
            break;
    filename[first_dot + 1] = 'h';
    filename[first_dot + 2] = '5';
    filename[first_dot + 3] = '\0';

    const int write = write_config_hdf5(p, filename);
    MPI_ERROR_CHECK(write,"Write failed.");
    free(filename);

    /* deallocate all memory */
    free_phase(p);

    finalize_MPI(&(p->info_MPI));
    }
