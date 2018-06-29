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
/* TO do Typedef for curand */
/* Clean up the hdf5 output*/

//! \file soma.c
//! \brief Implementation of the main executable SOMA.

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "io.h"
#include "mc.h"
#include "mesh.h"
#include "mpiroutines.h"
#include "init.h"
#include "test.h"
#include "ana.h"
#include "signal.h"
#include "rng.h"
#include "generate_positions.h"

int wait_for_debugger(void)
    {
    volatile int i=0;
    while( i == 0)
        ;
    return i;
    }

//! Main Function of the Executable SOMA
//! \private
//!
//! \param argc Argument counter
//! \param argv Argument vector
//! \return Errorcode
int main(int argc, char *argv[])
    {

    //wait_for_debugger();

    Phase phase;
    Phase *const p = &phase;

    if (MPI_Init(NULL, NULL) != MPI_SUCCESS) {
      fprintf(stderr, "MPI_ERROR (start)\n");
      return -1;
    }
    /* initialize MPI */

    const int error = MPI_Comm_dup( MPI_COMM_WORLD, &(p->info_MPI.SOMA_comm_world) );
    if( error != MPI_SUCCESS)
        {
        fprintf(stderr,"MPI Error cannot duplicate MPI_COMM_WORLD %s:%d\n",__FILE__,__LINE__);
        return -1;
        }
        {int tmp; MPI_Comm_rank( p->info_MPI.SOMA_comm_world, &tmp); p->info_MPI.world_rank=tmp;}
    const int error2 = MPI_Comm_dup( p->info_MPI.SOMA_comm_world, &(p->info_MPI.SOMA_comm_sim) );
    if( error2 != MPI_SUCCESS)
        {
        fprintf(stderr,"MPI Error cannot duplicate MPI_COMM_WORLD %s:%d\n",__FILE__,__LINE__);
        return -1;
        }

    const int args_success = cmdline_parser(argc,argv,&(p->args));
    if( args_success < 0 )
        {
        fprintf(stderr,"Process %d failed to read the cmdline. Exiting.\n",p->info_MPI.world_rank);
        return -2;
        }
    const int post_args= post_process_args( &(p->args), p->info_MPI.world_rank);
    if( post_args < 0)
        {
        fprintf(stderr,"Post processing the arguments on rank %d failed. Exiting.\n",p->info_MPI.world_rank);
        return -3;
        }

    const int mpi_init = init_MPI(p);
    if( mpi_init != 0)
        {
        fprintf(stderr,"ERROR: Unable to setup MPI %s:%d\n",__FILE__,__LINE__);
        exit(mpi_init);
        }

    const int signal_success = init_soma_signal();
    MPI_ERROR_CHECK(signal_success, "Signal init");

    const int mpi_args = check_status_on_mpi(p,args_success);
    if( mpi_args != 0 )
        {
        finalize_MPI(&(p->info_MPI));
        return 0;
        }
    if(p->args.move_type_arg == move_type_arg_TRIAL)
        {MPI_ERROR_CHECK(1,"ERROR: Trial move type is currently not working.");}

    if(p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_TT800)
        {MPI_ERROR_CHECK(1,"ERROR: TT800 PRNG is currently not working.");}

    const int open_acc = set_openacc_devices(p);
    if( check_status_on_mpi(p,open_acc) != 0){
        if(p->info_MPI.sim_rank == 0)
            fprintf(stderr,"ERROR: cannot set openacc devices.\n");
        finalize_MPI(&(p->info_MPI));
        return 1;
        }

    const unsigned int N_steps = p->args.timesteps_arg;
    /* read in configuration with routine from io */
    const int read = read_config_hdf5(p, p->args.coord_file_arg);
    MPI_ERROR_CHECK(read, "Cannot read coord file.");

    /* initialize phase */
    const int init = init_phase(p);
    MPI_ERROR_CHECK(init, "Cannot init values.");
    if( !p->bead_data_read )
        {

        if(p->info_MPI.sim_rank == 0)
            printf("INFO: Generating new bead initial data. "
                   "If your configuration contains rings "
                   "equilibration might take long.\n");
        const int new_beads = generate_new_beads(p);
        MPI_ERROR_CHECK(new_beads, "Cannot genrate new bead data.");
        //Reset the RNG to initial starting conditions.
        reseed(p, p->args.rng_seed_arg);
        }

    const int init_domain_chains_status = send_domain_chains(p,true);
    MPI_ERROR_CHECK(init_domain_chains_status, "Sending chains for domain decomposition failed.");

    if( ! p->args.skip_tests_flag)
        {
        const int test_p = test_particle_types(p);
        MPI_ERROR_CHECK(test_p, "Partile type test failed.");

        const int test51 = test_area51_violation(p);
        MPI_ERROR_CHECK(test51, "Area51 test failed.");

        const int test51_exact = test_area51_exact(p);
        if( ! p->args.nonexact_area51_flag )
            MPI_ERROR_CHECK(test51_exact, "Area51 exact test failed.");


        const int indepent_sets = test_independet_sets(p);
        MPI_ERROR_CHECK(indepent_sets, "Indepent Set test failed.");

        const int chains_domain = test_chains_in_domain(p);
        MPI_ERROR_CHECK(chains_domain, "Chains in domain test failed");
        }

    int stop_iteration = false;
    for (unsigned int i = 0; i < N_steps; i++) {
        const int mc_error = monte_carlo_propagation(p, 1);
        if( mc_error != 0)
            {
            fprintf(stderr,"ERROR %d in monte_carlo_propagation on rank %d.\n"
                    ,mc_error,p->info_MPI.world_rank);
            exit(mc_error);
            }
        analytics(p);
        screen_output(p,N_steps);
        if(p->args.load_balance_arg > 0 && i % p->args.load_balance_arg  == (unsigned int) p->args.load_balance_arg -1 )
            load_balance_mpi_ranks(p);
        if( p->args.N_domains_arg > 1 && p->args.rcm_update_arg > 0 && i % p->args.rcm_update_arg == (unsigned int) p->args.rcm_update_arg -1)
            {
            const int missed_chains = send_domain_chains(p, false);
            if( missed_chains != 0)
                exit(missed_chains);
            }

        stop_iteration = check_signal_stop();
        if( ! p->args.no_sync_signal_flag)
            {
            //Sync all mpi cores
            MPI_Allreduce(MPI_IN_PLACE,&stop_iteration,1,MPI_INT,MPI_SUM,p->info_MPI.SOMA_comm_world);
            }
        if(stop_iteration)
            {
            if(p->info_MPI.world_rank == 0)
                fprintf(stdout,"Signal to stop iteration at time %d catched by rank %d.\n",p->time,p->info_MPI.world_rank);
            break;
            }
    }
    const int missed_chains = send_domain_chains(p, false);
    if( missed_chains != 0)
        exit(missed_chains);

    const char *filename;
    const char normal[] = "end.h5";
    const char exitfile[] = "exit.h5";
    if( ! stop_iteration)
        filename = normal;
    else
        filename = exitfile;

    const int write = write_config_hdf5(p, filename);
    MPI_ERROR_CHECK(write, "Cannot write final configuration.");
    if( !stop_iteration && ! p->args.skip_tests_flag)
        {
        const int test51 = test_area51_violation(p);
        MPI_ERROR_CHECK(test51, "Area51 test failed.");
        const int test51_exact = test_area51_exact(p);
        if(! p->args.nonexact_area51_flag )
            MPI_ERROR_CHECK(test51_exact, "Area51 exact test failed.");

        const int chains_domain = test_chains_in_domain(p);
        MPI_ERROR_CHECK(chains_domain, "Chains in domain test failed");
        }

    /* deallocate all memory */
    free_phase(p);

    printf("Rank: %d \t polymers %ld\n",p->info_MPI.world_rank, p->n_polymers);

    finalize_MPI(&(p->info_MPI));
    if(p->info_MPI.world_rank == 0)
        printf("SOMA finished execution without errors.\n");
    return 0;
}
