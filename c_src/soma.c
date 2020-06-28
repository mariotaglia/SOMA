#define calc_avg 0
#define save_time_series 1

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
#include "polytype_conversion.h"

//! Main Function of the Executable SOMA
//! \private
//!
//! \param argc Argument counter
//! \param argv Argument vector
//! \return Errorcode
int main(int argc, char *argv[])
{
    Phase phase;
    Phase *const p = &phase;
#if ( ENABLE_MPI == 1 )
    if (MPI_Init(NULL, NULL) != MPI_SUCCESS)
        {
            fprintf(stderr, "MPI_ERROR (start)\n");
            return -1;
        }
    /* initialize MPI */

    const int error = MPI_Comm_dup(MPI_COMM_WORLD, &(p->info_MPI.SOMA_comm_world));
    if (error != MPI_SUCCESS)
        {
            fprintf(stderr, "MPI Error cannot duplicate MPI_COMM_WORLD %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    {
        int tmp;
        MPI_Comm_rank(p->info_MPI.SOMA_comm_world, &tmp);
        p->info_MPI.world_rank = tmp;
    }
    const int error2 = MPI_Comm_dup(p->info_MPI.SOMA_comm_world, &(p->info_MPI.SOMA_comm_sim));
    if (error2 != MPI_SUCCESS)
        {
            fprintf(stderr, "MPI Error cannot duplicate MPI_COMM_WORLD %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
#endif                          //ENABLE_MPI

    const int args_success = cmdline_parser(argc, argv, &(p->args));
    if (args_success < 0)
        {
            fprintf(stderr, "Process %d failed to read the cmdline. Exiting.\n", p->info_MPI.world_rank);
            return -2;
        }
    else if (args_success == 1) //Help and version output
        {
            return MPI_Finalize();
        }
    const int post_args = post_process_args(&(p->args), p->info_MPI.world_rank);
    if (post_args < 0)
        {
            fprintf(stderr, "Post processing the arguments on rank %d failed. Exiting.\n", p->info_MPI.world_rank);
            return -3;
        }

    const int mpi_init = init_MPI(p);
    if (mpi_init != 0)
        {
            fprintf(stderr, "ERROR: Unable to setup MPI %s:%d\n", __FILE__, __LINE__);
            exit(mpi_init);
        }

    if (p->args.move_type_arg == move_type_arg_TRIAL)
        {
            MPI_ERROR_CHECK(1, "ERROR: Trial move type is currently not working.");
        }

    if (p->args.pseudo_random_number_generator_arg == pseudo_random_number_generator_arg_TT800)
        {
            MPI_ERROR_CHECK(1, "ERROR: TT800 PRNG is currently not working.");
        }

    const int open_acc = set_openacc_devices(p);
    if (check_status_on_mpi(p, open_acc) != 0)
        {
            if (p->info_MPI.sim_rank == 0)
                fprintf(stderr, "ERROR: cannot set openacc devices.\n");
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
    init_nanoparticle_rng(p,&p->nanoparticles[0]);
    calc_np_field_total(p);
    if (!p->bead_data_read)
        {
          nanoparticle_area51_switch(p, 1);
            if (p->info_MPI.sim_rank == 0)
                printf("INFO: Generating new bead initial data. "
                       "If your configuration contains rings " "equilibration might take long.\n");
            const int new_beads = generate_new_beads(p);
            MPI_ERROR_CHECK(new_beads, "Cannot genrate new bead data.");
            //Reset the RNG to initial starting conditions.
            reseed(p, p->args.rng_seed_arg);
        }
#if ( ENABLE_MPI == 1 )
    const int init_domain_chains_status = send_domain_chains(p, true);
    MPI_ERROR_CHECK(init_domain_chains_status, "Sending chains for domain decomposition failed.");
#endif                          //ENABLE_MPI
    if (!p->args.skip_tests_flag)
        {
            const int test_p = test_particle_types(p);
            MPI_ERROR_CHECK(test_p, "Partile type test failed.");

            const int test51 = test_area51_violation(p);
            MPI_ERROR_CHECK(test51, "Area51 test failed.");

            test_area51_exact(p);

            const int indepent_sets = test_independet_sets(p);
            MPI_ERROR_CHECK(indepent_sets, "Indepent Set test failed.");

            const int chains_domain = test_chains_in_domain(p);
            MPI_ERROR_CHECK(chains_domain, "Chains in domain test failed");

            const int polytype_conversion = test_poly_conversion(p);
            MPI_ERROR_CHECK(polytype_conversion, "Polytype conversion test failed");
        }
    nanoparticle_area51_switch(p, 0);
    //    test_nanoparticle(p, &p->nanoparticles[0]);
    int stop_iteration = false;
#if(calc_avg==1)
      soma_scalar_t * avg_e_b=calloc(NUMBER_SOMA_BOND_TYPES , sizeof(soma_scalar_t));
      soma_scalar_t * avg_e_nb=calloc(p->n_types, sizeof(soma_scalar_t));
      soma_scalar_t e_b=0;
      soma_scalar_t e_b0=0;	
      soma_scalar_t e_nb=0;
#endif
#if(save_time_series==1)
      int frequency=10;
      int eq_time=N_steps/2;
      int nt=(N_steps-eq_time)/frequency;
      uint16_t * time_series=NULL;
      time_series = (uint16_t *) calloc(p->nx*p->ny*p->nz *nt, sizeof(uint16_t));
      int k=0;
#endif
    update_density_fields(p);
    update_omega_fields(p);
    for (unsigned int i = 0; i < N_steps; i++)
      {
	const int mc_error = monte_carlo_propagation(p, 1);
	if (mc_error != 0)
	  {
	    fprintf(stderr, "ERROR %d in monte_carlo_propagation on rank %d.\n", mc_error,
		    p->info_MPI.world_rank);
	    exit(mc_error);
	  }

	screen_output(p, N_steps);
#if(save_time_series==1)
	  if(i>=eq_time&&i%frequency==0)
	    {
	      memcpy(&time_series[k*p->n_cells], p->fields_unified, p->n_cells*sizeof(uint16_t));
	      k++;
	    }
#endif
#if(calc_avg==1)
	  if(i>=N_steps-N_steps/5){
	    average_field(p, p->umbrella_field, (soma_scalar_t) N_steps-(N_steps-N_steps/5));
	    e_nb+=calc_non_bonded_energy_exact(p)/(N_steps-(N_steps-N_steps/5));
	    calc_non_bonded_energy(p,&e_b0);
	    e_b+=e_b0/(N_steps/5);
	  }
#endif	    
      }
#if(calc_avg==1)
      p->xn[0]=e_b;
      p->k_umbrella[0]=e_nb;
#endif
#if(save_time_series==1)
      FILE *f = fopen("p", "w");
      for(int i=0;i<p->n_cells*nt;i++)
	fprintf(f, "%u\n",time_series[i]);
      fclose(f);
#endif
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
    if (p->info_MPI.world_rank == 0)
        printf("SOMA finished execution without errors.\n");
    return 0;
}
