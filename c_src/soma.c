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
#include "rng.h"
#include "walltime.h"
#include "generate_positions.h"
#include "server.h"
#include "phase.h"
#include "err_handling.h"

//! Main Function of the Executable SOMA
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

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    struct som_args args;
    const int args_success = cmdline_parser(argc, argv, &args);
    if (args_success < 0)
        {
            fprintf(stderr, "Process %d failed to read the cmdline. Exiting.\n", world_rank);
            return -2;
        }
    else if (args_success == 1) //Help and version output
        {
            return MPI_Finalize();
        }
    const int post_args = post_process_args(&args, world_rank);
    if (post_args < 0)
        {
            fprintf(stderr, "Post processing the arguments on rank %d failed. Exiting.\n", world_rank);
            return -3;
        }

    // get a role (server or simrank, which domain to cover...)
    struct role_info * ri;
    ri = assign_role(args.N_servers_arg, args.N_domains_arg);

    // get some global constants (like n_polymers_global)
    struct global_consts gc;
    read_consts_from_config(&gc, args.coord_file_arg);

    if (ri->is_server)
        {
            // run this rank as a server
            const struct server_info * serv_inf = ri->detailed_info;

            soma_scalar_t *field_scaling_type = (soma_scalar_t *)malloc(sizeof(soma_scalar_t) * gc.n_types);
            RET_ERR_ON_NULL(field_scaling_type, "Malloc");

            receive_field_scaling_type(field_scaling_type, gc.n_types, serv_inf);

            Ana_Info ai;
            unsigned int * end_mono;

            init_ana(&ai, &end_mono, &gc, field_scaling_type,
                args.ana_file_arg, args.coord_file_arg,
                serv_inf->server_comm.comm);

            //early exit for testing
            {
                printf("worldrank %d (server) reached end of initialization\n", world_rank);
                MPI_Finalize();
                return 0;
            }

            MPI_Finalize();
            free(field_scaling_type);
            free(end_mono);
        }
    else
        {
            // run this rank as a simulation rank
            const struct sim_rank_info *sim_inf = ri->detailed_info;
            Phase phase;
            Phase *const p = &phase;
            p->args = args;

            p->info_MPI.SOMA_comm_world = MPI_COMM_WORLD;
            {
                int tmp;
                MPI_Comm_rank (p->info_MPI.SOMA_comm_world, &tmp);
                p->info_MPI.world_rank = tmp;
            }
            p->info_MPI.SOMA_comm_sim = sim_inf->sim_comm.comm;

            const int mpi_init = init_MPI (p, sim_inf);
            if (mpi_init != 0)
                {
                    fprintf (stderr, "ERROR: Unable to setup MPI %s:%d\n", __FILE__, __LINE__);
                    exit (mpi_init);
                }

            if (p->args.move_type_arg == move_type_arg_TRIAL)
                {
                    MPI_ERROR_CHECK(1, "ERROR: Trial move type is currently not working.");
                }

            if (p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_TT800)
                {
                    MPI_ERROR_CHECK(1, "ERROR: TT800 PRNG is currently not working.");
                }

            const int open_acc = set_openacc_devices (p);
            if (check_status_on_mpi (p, open_acc) != 0)
                {
                    if (p->info_MPI.sim_rank == 0)
                        fprintf (stderr, "ERROR: cannot set openacc devices.\n");
                    finalize_MPI (&(p->info_MPI));
                    return 1;
                }
            const unsigned int N_steps = p->args.timesteps_arg;

            /* read in configuration with routine from io */
            int status;
            status = process_consts_into_phase (p, &gc, sim_inf);
            MPI_ERROR_CHECK(status, "cannot put simulation constants into the phase");
            status = read_poly_and_fields_hdf5 (p, p->args.coord_file_arg, sim_inf);
            MPI_ERROR_CHECK(status, "Cannot read coord file.");

            /* initialize phase */
            const int init = init_phase (p, sim_inf, &gc);
            MPI_ERROR_CHECK(init, "Cannot init values.");
            if (!p->bead_data_read)
                {

                    if (p->info_MPI.sim_rank == 0)
                        printf ("INFO: Generating new bead initial data. "
                                "If your configuration contains rings " "equilibration might take long.\n");
                    const int new_beads = generate_new_beads (p);
                    MPI_ERROR_CHECK(new_beads, "Cannot genrate new bead data.");
                    //Reset the RNG to initial starting conditions.
                    reseed (p, p->args.rng_seed_arg);
                }
#if (ENABLE_MPI == 1)
            const int init_domain_chains_status = send_domain_chains (p, true);
            MPI_ERROR_CHECK(init_domain_chains_status, "Sending chains for domain decomposition failed.");
#endif                          //ENABLE_MPI

            if (!p->args.skip_tests_flag)
                {
                    const int test_p = test_particle_types (p);
                    MPI_ERROR_CHECK(test_p, "Partile type test failed.");

                    const int test51 = test_area51_violation (p);
                    MPI_ERROR_CHECK(test51, "Area51 test failed.");

                    test_area51_exact (p);

                    const int indepent_sets = test_independet_sets (p);
                    MPI_ERROR_CHECK(indepent_sets, "Indepent Set test failed.");

                    const int chains_domain = test_chains_in_domain (p);
                    MPI_ERROR_CHECK(chains_domain, "Chains in domain test failed");
                }

            // early exit for testing
            {
                printf("Worldrank %d (simrank) reached the main loop\n", world_rank);
                MPI_Finalize();
                return 0;
            }


            int stop_iteration = false;
            for (unsigned int i = 0; i < N_steps; i++)
                {
                    analytics(p);
                    const int mc_error = monte_carlo_propagation(p, 1);
                    if (mc_error != 0)
                        {
                            fprintf(stderr, "ERROR %d in monte_carlo_propagation on rank %d.\n", mc_error,
                                    p->info_MPI.world_rank);
                            exit(mc_error);
                        }
                    screen_output(p, N_steps);
#if ( ENABLE_MPI == 1 )
                    if (p->args.load_balance_arg > 0
                        && i % p->args.load_balance_arg == (unsigned int)p->args.load_balance_arg - 1)
                        load_balance_mpi_ranks(p);
                    if (p->args.N_domains_arg > 1 && p->args.rcm_update_arg > 0
                        && i % p->args.rcm_update_arg == (unsigned int)p->args.rcm_update_arg - 1)
                        {
                            const int missed_chains = send_domain_chains(p, false);
                            if (missed_chains != 0)
                                exit(missed_chains);
                        }
#endif                          //ENABLE_MPI

                    stop_iteration = check_walltime_stop();
#if ( ENABLE_MPI == 1 )
                    if (!p->args.no_sync_signal_flag)
                        {
                            //Sync all mpi cores
                            MPI_Allreduce(MPI_IN_PLACE, &stop_iteration, 1, MPI_INT, MPI_SUM, p->info_MPI.SOMA_comm_world);
                        }
#endif                          //ENABLE_MPI

                    if (stop_iteration)
                        {
                            if (p->info_MPI.world_rank == 0)
                                fprintf(stdout, "Environment to stop iteration at time %d catched by rank %d.\n", p->time,
                                        p->info_MPI.world_rank);
                            break;
                        }
                }



#if ( ENABLE_MPI == 1 )
            const int missed_chains = send_domain_chains(p, false);
            if (missed_chains != 0)
                exit(missed_chains);
#endif                          //ENABLE_MPI


            const int write = write_config_hdf5(p, p->args.final_file_arg);
            MPI_ERROR_CHECK(write, "Cannot write final configuration.");


            if (!stop_iteration && !p->args.skip_tests_flag)
                {
                    const int test51 = test_area51_violation(p);
                    MPI_ERROR_CHECK(test51, "Area51 test failed.");

                    test_area51_exact(p);

                    const int chains_domain = test_chains_in_domain(p);
                    MPI_ERROR_CHECK(chains_domain, "Chains in domain test failed");
                }

            free_phase(p);

            printf("Rank: %d \t polymers %ld\n", p->info_MPI.world_rank, p->n_polymers);

            finalize_MPI(&(p->info_MPI));

        }

    free_global_consts(&gc);
    free_role_info(ri);
    if (world_rank == 0)
        {
            printf("SOMA finished execution without errors.\n");
        }
    return 0;
}
