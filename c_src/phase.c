/* Copyright (C) 2016-2019 Ludwig Schneider

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

//! \file phase.c
//! \brief Implementation of phase.h

#include "phase.h"
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "init.h"
#include "independent_sets.h"
#include "mesh.h"
#include "err_handling.h"
#include "send.h"


void free_global_consts(struct global_consts * gc)
{
    free(gc->xn);
    free(gc->A);
    free(gc->k_umbrella);
    free(gc->poly_type_offset);
    free(gc->poly_arch);
    free(gc->cm_a);
}


int process_consts_into_phase(struct Phase * p, const struct global_consts * gc,
    const struct sim_rank_info * sri)
{

    // n_polymers_global
    p->n_polymers_global = gc->n_polymers_global;
    //Distribute the polymers to different cores.
    uint64_t n_polymers = p->n_polymers_global / sri->sim_comm.size;
    if ((unsigned int)sri->sim_comm.rank < p->n_polymers_global % sri->sim_comm.size)
        n_polymers += 1;
    p->n_polymers = n_polymers;
    p->n_polymers_storage = p->n_polymers;

    // n_beads
    p->reference_Nbeads = gc->reference_Nbeads;

    // n_types
    p->n_types = gc->n_types;

    // xn
    size_t xn_size = p->n_types * p->n_types * sizeof(soma_scalar_t);
    p->xn = (soma_scalar_t *) malloc(p->n_types * p->n_types * sizeof(soma_scalar_t));
    RET_ERR_ON_NULL(p->xn, "Malloc");
    memcpy(p->xn, gc->xn, xn_size);

    // A array (diffusivity of particles)
    size_t A_size = p->n_types * sizeof(soma_scalar_t);
    p->A = (soma_scalar_t *)malloc(A_size);
    RET_ERR_ON_NULL(p->A, "Malloc");
    memcpy(p->A, gc->A, A_size);

    // time
    p->start_time = gc->start_time;
    p->time = p->start_time;

    //hamiltonian
    p->hamiltonian = gc->hamiltonian;

    // k_umbrella
    size_t k_umbrella_size = p->n_types * sizeof(soma_scalar_t);
    p->k_umbrella = (soma_scalar_t * ) malloc(k_umbrella_size);
    RET_ERR_ON_NULL(p->k_umbrella, "Malloc");
    memcpy(p->k_umbrella, gc->k_umbrella, k_umbrella_size);

    // nx, ny, nz
    p->nx = gc->nx;
    p->ny = gc->ny;
    p->nz = gc->nz;

    if (p->nx % sri->n_domains != 0)
        {
            fprintf(stderr, "ERROR: %s:%d\n\t"
                    "The nx %d number is not divible by the number of domains %d\n",
                    __FILE__, __LINE__, p->nx, sri->n_domains);
            return -3;
        }

    // lx, ly, lz
    p->Lx = gc->Lx;
    p->Ly = gc->Ly;
    p->Lz = gc->Lz;

    // n_polymer_type
    p->n_poly_type = gc->n_poly_type;

    // polymer architecture
    size_t poly_type_offset_size = p->n_poly_type * sizeof(int);
    p->poly_type_offset = (int *)malloc(poly_type_offset_size);
    RET_ERR_ON_NULL(p->poly_type_offset, "Malloc");
    memcpy(p->poly_type_offset, gc->poly_type_offset, poly_type_offset_size);

    p->poly_arch_length = gc->poly_arch_length;
    size_t poly_arch_size = p->poly_arch_length * sizeof(uint32_t);
    p->poly_arch = (uint32_t *) malloc(poly_arch_size);
    RET_ERR_ON_NULL(p->poly_arch, "Malloc");
    memcpy(p->poly_arch, gc->poly_arch, poly_arch_size);

    // cm_a
    if (gc->cm_a == NULL)
        {
            p->cm_a = NULL; //null is a valid value meaning deactivated
        }
    else
        {
            size_t cm_a_size = p->n_poly_type * sizeof(soma_scalar_t);
            p->cm_a = (soma_scalar_t *)malloc(cm_a_size);
            RET_ERR_ON_NULL(p->cm_a, "Malloc");
            memcpy(p->cm_a, gc->cm_a, cm_a_size);
        }

    // harmonic_normb
    p->harmonic_normb_variable_scale = gc->harmonic_normb_variable_scale;

    return 0;


}

int init_phase (struct Phase *const p, const struct sim_rank_info *sri, const struct global_consts *gc)
{
    print_version(p->info_MPI.sim_rank);
    p->present_on_device = false;
    p->start_time = p->time;
    p->start_clock = time(NULL);
    p->n_accepts = 0;
    p->n_moves = 0;
    p->n_tries_cm = 0;
    p->n_acc_cm = 0;
    p->end_mono = NULL;
    p->tps_elapsed_time = 1;    //Bla default, bigger 0
    p->tps_elapsed_steps = 1;   //Bla default, bigger 0

    uint64_t n_polymer_offset = 0;
#if ( ENABLE_MPI == 1 )
    MPI_Scan(&(p->n_polymers), &n_polymer_offset, 1, MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);
    n_polymer_offset -= p->n_polymers;
#endif                          //ENABLE_MPI

    for (uint64_t i = 0; i < p->n_polymers; i++)
        {
            p->polymers[i].set_states = NULL;
            p->polymers[i].set_permutation = NULL;

            allocate_rng_state(&(p->polymers[i].poly_state), p->args.pseudo_random_number_generator_arg);
            seed_rng_state(&(p->polymers[i].poly_state), p->args.rng_seed_arg,
                           i + n_polymer_offset, p->args.pseudo_random_number_generator_arg);
        }

    // Max safe move distance
    p->max_safe_jump = p->Lx / p->nx < p->Ly / p->ny ? p->Lx / p->nx : p->Ly / p->ny;
    p->max_safe_jump = p->max_safe_jump < p->Lz / p->nz ? p->max_safe_jump : p->Lz / p->nz;
    p->max_safe_jump *= 0.95;

    // Reference Harmonic Spring Cste
    const soma_scalar_t harmonic_spring_Cste = 1.0 / sqrt(3.0 * (p->reference_Nbeads - 1.0));
    //Reference energy scale for harmonic springs.
    p->harmonic_normb = 1.0 / (2.0 * harmonic_spring_Cste * harmonic_spring_Cste);

#if ( ENABLE_MPI == 1 )
    uint64_t n_polymers_global_sum;
    MPI_Allreduce(&(p->n_polymers), &n_polymers_global_sum, 1, MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);
    assert(p->n_polymers_global == n_polymers_global_sum);
#endif                          //ENABLE_MPI

    p->n_cells = p->nx * p->ny * p->nz;
    unsigned int my_domain;
    if (sri == NULL)
        {
            my_domain = p->info_MPI.sim_rank / p->args.N_domains_arg;
        }
    else
        {
            my_domain = sri->my_domain;
        }
    p->local_nx_low = (p->nx / p->args.N_domains_arg * my_domain) - p->args.domain_buffer_arg;
    p->local_nx_high = p->local_nx_low + p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg;
    p->n_cells_local = (p->local_nx_high - p->local_nx_low) * p->ny * p->nz;

    //Check if it is a valid domain decomposition
    // if this condition doesn't hold, the buffers from domain 0 and 2 would overlap!
    if (p->args.N_domains_arg * 2 * p->args.domain_buffer_arg > (int)p->nx)
        {
            fprintf(stderr, "ERROR: invalid domain decomposition. %s:%d\n", __FILE__, __LINE__);
            fprintf(stderr, "\t N(N_domains)= %d\tb(domain_buffer)= %d\t nx= %d\n", p->args.N_domains_arg,
                    p->args.domain_buffer_arg, p->nx);
            fprintf(stderr, "\t N * 2 * b <= nx \t not fulfilled\n");
            return -2;
        }

    //Allocate Fields
    p->left_tmp_buffer = NULL;
    p->right_tmp_buffer = NULL;

    if (p->args.N_domains_arg > 1 && p->info_MPI.domain_rank == 0)
        {
            p->left_tmp_buffer = (uint16_t *) malloc(p->args.domain_buffer_arg * p->ny * p->nz * sizeof(uint16_t));
            if (p->left_tmp_buffer == NULL)
                {
                    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
                    return -1;
                }
            p->right_tmp_buffer = (uint16_t *) malloc(p->args.domain_buffer_arg * p->ny * p->nz * sizeof(uint16_t));
            if (p->right_tmp_buffer == NULL)
                {
                    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
                    return -1;
                }
        }

    p->fields_unified = (uint16_t *) malloc(p->n_cells_local * p->n_types * sizeof(uint16_t));
    if (p->fields_unified == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    p->old_fields_unified = (uint16_t *) malloc(p->n_cells_local * p->n_types * sizeof(uint16_t));
    if (p->old_fields_unified == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    p->fields_32 = (uint32_t *) malloc(p->n_types * p->n_cells_local * sizeof(uint32_t));
    if (p->fields_32 == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    p->omega_field_unified = (soma_scalar_t *) malloc(p->n_cells_local * p->n_types * sizeof(soma_scalar_t));
    if (p->omega_field_unified == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    p->tempfield = (soma_scalar_t *) malloc(p->n_cells_local * sizeof(soma_scalar_t));
    if (p->tempfield == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    p->num_bead_type = (uint64_t *) malloc(p->n_types * sizeof(uint64_t));
    if (p->num_bead_type == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    p->num_bead_type_local = (uint64_t *) malloc(p->n_types * sizeof(uint64_t));
    if (p->num_bead_type_local == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    p->field_scaling_type = (soma_scalar_t *) malloc(p->n_types * sizeof(soma_scalar_t));
    if (p->field_scaling_type == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    // Set all values to zero

    p->num_all_beads = 0;
    p->num_all_beads_local = 0;
    for (unsigned int i = 0; i < p->n_types; i++)
        p->num_bead_type_local[i] = 0;
    // Determine number of  different bead types
    for (uint64_t j = 0; j < p->n_polymers; j++)
        {                       /*Loop over polymers */
            const unsigned int N = p->poly_arch[p->poly_type_offset[p->polymers[j].type]];
            for (unsigned int k = 0; k < N; k++)
                {               /*Loop over monomers */
                    const unsigned int type =
                        get_particle_type(p->poly_arch[p->poly_type_offset[p->polymers[j].type] + 1 + k]);
                    p->num_bead_type_local[type] += 1;
                    p->num_all_beads_local += 1;
                }
        }

    // Share p->num_all_beads
    MPI_Allreduce(&(p->num_all_beads_local), &(p->num_all_beads), 1, MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);

    // Share p->num_bead_type
    for (unsigned int i = 0; i < p->n_types; i++)
        {
            MPI_Allreduce(&(p->num_bead_type_local[i]), &(p->num_bead_type[i]),
                          1, MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);
        }
    // Check if uint16_t density field is enough
    soma_scalar_t check_short = (double) p->num_all_beads / p->n_cells;

    if (check_short > (((double)USHRT_MAX / 100) * 95))
        {
            fprintf(stderr, "ERROR: Cell-density above 95 Percent of USHRT_MAX\n");
            return -1;
        }
    // setting the geometry field, or if it's already initialized measure the free space
    uint64_t ncells = (p->nx / p->args.N_domains_arg) * p->ny * p->nz;
    if (p->info_MPI.domain_rank == 0)   //Only domain root calculats something else than 0
        {
            if (p->area51 != NULL)
                {
                    // substract the number of non free cells for the correct density scaling
                    for (int x = p->local_nx_low + p->args.domain_buffer_arg;
                         x < p->local_nx_high - p->args.domain_buffer_arg; x++)
                        for (unsigned int y = 0; y < p->ny; y++)
                            for (unsigned int z = 0; z < p->nz; z++)
                                {
                                    uint64_t cell = cell_coordinate_to_index(p, x, y, z);
                                    if (p->area51[cell] > 0)
                                        ncells--;
                                }
                }
        }
    else
        ncells = 0;             //Not domain rank == 0, count only once per domain
#if ( ENABLE_MPI == 1 )
    MPI_Allreduce(MPI_IN_PLACE, &ncells, 1, MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);
#endif                          //ENABLE_MPI

    // Loop to calculate scaling parameter
    for (unsigned int i = 0; i < p->n_types; i++)
        p->field_scaling_type[i] = (ncells / ((soma_scalar_t) p->num_all_beads));
    // Info for Ulrich: programm does take excluded volume into account now!
    p->n_accepts = 0;
    p->n_moves = 0;

    p->R = (soma_scalar_t *) malloc(p->n_types * sizeof(soma_scalar_t));
    for (unsigned int i = 0; i < p->n_types; i++)
        {
            //! \todo kBT required.
            p->R[i] = sqrt(p->A[i] * 2);
        }

    // initialize inverse simulation cell parameters
    p->iLx = 1.0 / p->Lx;
    p->iLy = 1.0 / p->Ly;
    p->iLz = 1.0 / p->Lz;

    p->sets = NULL;             // Default init of the sets
    p->max_set_members = 0;
    if (p->args.iteration_alg_arg == iteration_alg_arg_SET)
        generate_independet_sets(p);

    init_autotuner(&(p->mc_autotuner));
    init_autotuner(&(p->cm_mc_autotuner));

    copyin_phase(p);
    p->num_long_chain = mc_set_init(p);

    // call update_fields routine

    if (p->bead_data_read)
        {
#if ( ENABLE_MPI == 1 )
            const int init_domain_chains_status = send_domain_chains(p, true);
            if (init_domain_chains_status != 0)
                return init_domain_chains_status;
#endif                          //ENABLE_MPI
            update_density_fields(p);
            memcpy(p->old_fields_unified, p->fields_unified, p->n_cells_local * p->n_types * sizeof(uint16_t));
        }

    int ret = 0;
    if (p->args.coord_file_arg != NULL) //Is it a full init Phase?
        {
            //int init_ana(Ana_Info *ana_info, unsigned int * *end_mono, const global_consts *gc, soma_scalar_t *field_scaling_type, const char * const filename, const char * const coord_filename, MPI_Comm writer_comm)
            int status;
            status = send_field_scaling_type(p->field_scaling_type, p->n_types, sri);
            MPI_ERROR_CHECK(status, "sending field-scaling to server failed");
            ret = init_ana(&(p->ana_info), &(p->end_mono), gc, p->field_scaling_type,
                p->args.ana_file_arg, p->args.coord_file_arg, MPI_COMM_NULL);
        }
    else
        {
            ret = init_ana(&p->ana_info, &(p->end_mono), gc, p->field_scaling_type,
                NULL, NULL, MPI_COMM_NULL);
        }

    return ret;
}


int copyin_phase(struct Phase *const p)
{
    if (p->present_on_device)
        {
            fprintf(stderr, "WARNING: %s:%d copyin of the phase, but system seems to be already present on device.\n",
                    __FILE__, __LINE__);
        }

#ifdef _OPENACC
#pragma acc enter data copyin(p[0:1])
#pragma acc enter data copyin(p->xn[0:p->n_types*p->n_types])
#pragma acc enter data copyin(p->polymers[0:p->n_polymers_storage])
#pragma acc enter data copyin(p->fields_unified[0:p->n_types*p->n_cells_local])
#pragma acc enter data copyin(p->old_fields_unified[0:p->n_types*p->n_cells_local])
#pragma acc enter data copyin(p->fields_32[0:p->n_types*p->n_cells_local])
    if (p->area51 != NULL)
        {
#pragma acc enter data copyin(p->area51[0:p->n_cells_local])
        }
#pragma acc enter data copyin(p->omega_field_unified[0:p->n_cells_local*p->n_types])
    if (p->external_field_unified != NULL)
        {
#pragma acc enter data copyin(p->external_field_unified[0:p->n_cells_local*p->n_types])
        }
    if (p->umbrella_field != NULL)
        {
#pragma acc enter data copyin(p->umbrella_field[0:p->n_cells_local*p->n_types])
        }
#pragma acc enter data copyin(p->tempfield[0:p->n_cells_local])
#pragma acc enter data copyin(p->num_bead_type[0:p->n_types])
#pragma acc enter data copyin(p->num_bead_type_local[0:p->n_types])
#pragma acc enter data copyin(p->A[0:p->n_types])
#pragma acc enter data copyin(p->R[0:p->n_types])
#pragma acc enter data copyin(p->field_scaling_type[0:p->n_types])
#pragma acc enter data copyin(p->k_umbrella[0:p->n_types])
#pragma acc enter data copyin(p->poly_type_offset[0:p->n_poly_type])
#pragma acc enter data copyin(p->poly_arch[0:p->poly_arch_length])

    if (p->cm_a != NULL)
        {
#pragma acc enter data copyin(p->cm_a[0:p->n_poly_type])
        }
    if (p->sets != NULL)
        {
#pragma acc enter data copyin(p->sets[0:p->n_poly_type])
            for (unsigned int i = 0; i < p->n_poly_type; i++)
                {
#pragma acc enter data copyin(p->sets[i].set_length[0:p->sets[i].n_sets])
#pragma acc enter data copyin(p->sets[i].sets[0:p->sets[i].n_sets*p->sets[i].max_member])
                }
        }
    for (uint64_t i = 0; i < p->n_polymers; i++)
        {
            Polymer *const poly = &(p->polymers[i]);
            copyin_polymer(p, poly);
        }
#endif                          //_OPENACC

    p->present_on_device = true;
    return p->n_polymers * 0 + 1;
}

int copyout_phase(struct Phase *const p)
{
    if (!p->present_on_device)
        {
            fprintf(stderr, "WARNING: %s:%d copyout of the phase, but system seems to be not present on device.\n",
                    __FILE__, __LINE__);
        }
#ifdef _OPENACC

#pragma acc exit data copyout(p->xn[0:p->n_types*p->n_types])
#pragma acc exit data copyout(p->fields_unified[0:p->n_types*p->n_cells_local])
#pragma acc exit data copyout(p->old_fields_unified[0:p->n_types*p->n_cells_local])
#pragma acc exit data copyout(p->fields_32[0:p->n_types*p->n_cells_local])
    if (p->area51 != NULL)
        {
#pragma acc exit data copyout(p->area51[0:p->n_cells_local])
        }
#pragma acc exit data copyout(p->omega_field_unified[0:p->n_cells_local*p->n_types])
    if (p->external_field_unified != NULL)
        {
#pragma acc exit data copyout(p->external_field_unified[0:p->n_cells_local*p->n_types])
        }
    if (p->umbrella_field != NULL)
        {
#pragma acc exit data copyout(p->umbrella_field[0:p->n_cells_local*p->n_types])
        }
#pragma acc exit data copyout(p->tempfield[0:p->n_cells_local])
#pragma acc exit data copyout(p->num_bead_type[0:p->n_types])
#pragma acc exit data copyout(p->num_bead_type_local[0:p->n_types])
#pragma acc exit data copyout(p->A[0:p->n_types])
#pragma acc exit data copyout(p->R[0:p->n_types])
#pragma acc exit data copyout(p->field_scaling_type[0:p->n_types])
#pragma acc exit data copyout(p->k_umbrella[0:p->n_types])
#pragma acc exit data copyout(p->poly_type_offset[0:p->n_poly_type])
#pragma acc exit data copyout(p->poly_arch[0:p->poly_arch_length])

    if (p->cm_a != NULL)
        {
#pragma acc exit data copyout(p->cm_a[0:p->n_poly_type])
        }
    if (p->sets != NULL)
        {
            for (unsigned int i = 0; i < p->n_poly_type; i++)
                {
#pragma acc exit data copyout(p->sets[i].set_length[0:p->sets[i].n_sets])
#pragma acc exit data copyout(p->sets[i].sets[0:p->sets[i].n_sets*p->sets[i].max_member])
                }
#pragma acc exit data copyout(p->sets[0:p->n_poly_type])
        }
    for (uint64_t i = 0; i < p->n_polymers; i++)
        {
            Polymer *const poly = &(p->polymers[i]);
            copyout_polymer(p, poly);
        }
#pragma acc exit data copyout(p->polymers[0:p->n_polymers_storage])
    //Use here the delete to not overwrite stuff, which only changed on CPU
#pragma acc exit data delete(p[0:1])
#endif                          //_OPENACC

    p->present_on_device = false;
    return p->n_polymers * 0 + 1;
}

int free_phase(struct Phase *const p)
{
    copyout_phase(p);

    /* de-allocate fields */
    free(p->left_tmp_buffer);
    free(p->right_tmp_buffer);
    free(p->omega_field_unified);
    free(p->tempfield);
    free(p->fields_unified);
    free(p->old_fields_unified);
    free(p->fields_32);
    free(p->num_bead_type);
    free(p->num_bead_type_local);
    free(p->field_scaling_type);
    free(p->k_umbrella);
    free(p->A);
    free(p->R);
    free(p->end_mono);
    free(p->cm_a);

    if (p->args.coord_file_arg != NULL) //Required option. So always true, if args used.
        cmdline_parser_free(&(p->args));

    if (p->sets != NULL)
        {
            for (unsigned int i = 0; i < p->n_poly_type; i++)
                {
                    free(p->sets[i].set_length);
                    free(p->sets[i].sets);
                }
            free(p->sets);
        }

    /* free polymers */
    for (uint64_t i = 0; i < p->n_polymers; i++)
        {
            Polymer *const poly = &(p->polymers[i]);
            free_polymer(p, poly);
        }

    free(p->polymers);

    free(p->poly_type_offset);
    free(p->poly_arch);

    free(p->xn);

    if (p->area51 != NULL)
        {
            free(p->area51);
        }

    if (p->external_field_unified != NULL)
        {
            free(p->external_field_unified);
        }

    if (p->umbrella_field != NULL)
        {
            free(p->umbrella_field);
        }

    close_ana(&(p->ana_info));
    return 0;
}

int update_self_phase(const Phase * const p, int rng_update_flag)
{
    static unsigned int last_time_call = 0;
    if (last_time_call == 0 || p->time > last_time_call)
        last_time_call = p->time;
    else                        //Quick exit, because the property has already been calculated for the time step.
        return 1;

    // Not pointer members are expected to not change on device
#pragma acc update self(p->xn[0:p->n_types*p->n_types])

    for (uint64_t i = 0; i < p->n_polymers; i++)
        update_self_polymer(p, p->polymers + i, rng_update_flag);

#pragma acc update self(p->fields_unified[0:p->n_cells_local*p->n_types])
#pragma acc update self(p->old_fields_unified[0:p->n_types*p->n_cells_local])
#pragma acc update self(p->fields_32[0:p->n_types*p->n_cells_local])

    if (p->area51 != NULL)
        {
#pragma acc update self(p->area51[0:p->n_cells_local])
        }
#pragma acc update self(p->omega_field_unified[0:p->n_cells_local*p->n_types])
    if (p->external_field_unified != NULL)
        {
#pragma acc update self(p->external_field_unified[0:p->n_cells_local*p->n_types])
        }
    if (p->umbrella_field != NULL)
        {
#pragma acc update self(p->umbrella_field[0:p->n_cells_local*p->n_types])
        }

#pragma acc update self(p->tempfield[0:p->n_cells_local])
#pragma acc update self(p->num_bead_type[0:p->n_types])
#pragma acc update self(p->num_bead_type_local[0:p->n_types])
#pragma acc update self(p->A[0:p->n_types])
#pragma acc update self(p->R[0:p->n_types])
#pragma acc update self(p->field_scaling_type[0:p->n_types])
#pragma acc update self(p->k_umbrella[0:p->n_types])
#pragma acc update self(p->poly_type_offset[0:p->n_poly_type])
#pragma acc update self(p->poly_arch[0:p->poly_arch_length])
    if (p->cm_a != NULL)
        {
#pragma acc update self(p->cm_a[0:p->n_poly_type])
        }

    //SETS are not updated to host

    return p->n_polymers * 0 + 1;
}

int mc_set_init(Phase * const p)
{
    int num_long_chain = 0;
    unsigned int *poly_order = (unsigned int *)malloc((int)p->n_polymers * sizeof(unsigned int));
    if (poly_order == NULL)
        {
            fprintf(stderr, "ERROR: malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    memset(poly_order, 0, (int)p->n_polymers * sizeof(unsigned int));

    Polymer *const first_poly = &p->polymers[0];
    const unsigned int poly_type = first_poly->type;
    uint32_t length_poly_start = p->poly_arch[p->poly_type_offset[poly_type]];
    if (length_poly_start > p->num_all_beads / p->args.long_chain_threshold_arg)
        num_long_chain++;
    for (uint64_t poly_i = 1; poly_i < p->n_polymers; poly_i++)
        {
            Polymer *const this_poly = &p->polymers[poly_i];
            const unsigned int poly_type = this_poly->type;
            uint32_t length_poly_i = p->poly_arch[p->poly_type_offset[poly_type]];
            if (length_poly_i > p->num_all_beads / p->args.long_chain_threshold_arg)
                num_long_chain++;

            int i = poly_i - 1;
            while (i >= 0 && length_poly_i > p->poly_arch[p->poly_type_offset[(&p->polymers[poly_order[i]])->type]])
                i--;

            i = poly_i - 1 - i;
            if (i != 0)
                memmove(poly_order + poly_i - i, poly_order + poly_i - i + 1, i * sizeof(unsigned int));

            poly_order[poly_i - i] = poly_i;
        }

    for (int index = 0; index < num_long_chain; index++)
        {
            if (poly_order[index] != (unsigned int)index)
                {
                    copyout_phase(p);
                    exchange_polymer(p, poly_order[index], index);
                    copyin_phase(p);
                }
        }
    free(poly_order);
    return num_long_chain;
}
