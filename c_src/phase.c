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

//! \file phase.c
//! \brief Implementation of phase.h

#include "phase.h"
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include "init.h"
#include "independent_sets.h"
#include "mesh.h"
#include "polytype_conversion.h"
#include "monotype_conversion.h"
#include "mobility.h"
#include "self_documentation.h"
#include "poly_heavy.h"

int init_phase(struct Phase *const p)
{
    p->present_on_device = false;
    p->start_time = p->time;
    gettimeofday(&(p->start_clock), NULL);
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

    switch (p->args.pseudo_random_number_generator_arg)
        {
        case pseudo_random_number_generator__NULL:
            /* intentionally falls through */
        case pseudo_random_number_generator_arg_PCG32:
            init_soma_memory(&(p->rh.mt), 0, 0);
            init_soma_memory(&(p->rh.tt800), 0, 0);
            break;
        case pseudo_random_number_generator_arg_MT:
            init_soma_memory(&(p->rh.mt), p->n_polymers, sizeof(MERSENNE_TWISTER_STATE));
            init_soma_memory(&(p->rh.tt800), 0, 0);
            break;
        case pseudo_random_number_generator_arg_TT800:
            init_soma_memory(&(p->rh.mt), 0, 0);
            init_soma_memory(&(p->rh.tt800), p->n_polymers, sizeof(TT800STATE));
            break;
        }
    for (uint64_t i = 0; i < p->n_polymers; i++)
        {
            p->polymers[i].set_states_offset = UINT64_MAX;
            p->polymers[i].set_permutation_offset = UINT64_MAX;
            switch (p->args.pseudo_random_number_generator_arg)
                {
                case pseudo_random_number_generator__NULL:
                    /* intentionally falls through */
                case pseudo_random_number_generator_arg_PCG32:
                    p->polymers[i].poly_state.alternative_rng_offset = UINT64_MAX;
                    break;
                case pseudo_random_number_generator_arg_MT:
                    p->polymers[i].poly_state.alternative_rng_offset = get_new_soma_memory_offset(&(p->rh.mt), 1);
                    break;
                case pseudo_random_number_generator_arg_TT800:
                    p->polymers[i].poly_state.alternative_rng_offset = get_new_soma_memory_offset(&(p->rh.tt800), 1);
                    break;
                }
            seed_rng_state(&(p->polymers[i].poly_state), p->args.rng_seed_arg, i + n_polymer_offset, p);
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
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    p->local_nx_low = (p->nx / p->args.N_domains_arg * my_domain) - p->args.domain_buffer_arg;
    p->local_nx_high = p->local_nx_low + p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg;
    p->n_cells_local = (p->local_nx_high - p->local_nx_low) * p->ny * p->nz;

    //Check if it is a valid domain decomposition
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
            p->left_tmp_buffer = (uint16_t *) malloc((p->args.domain_buffer_arg * p->ny * p->nz) * sizeof(uint16_t));
            if (p->left_tmp_buffer == NULL)
                {
                    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
                    return -1;
                }
            p->right_tmp_buffer = (uint16_t *) malloc((p->args.domain_buffer_arg * p->ny * p->nz) * sizeof(uint16_t));
            if (p->right_tmp_buffer == NULL)
                {
                    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
                    return -1;
                }
        }
    p->fields_unified = (uint16_t *) malloc(p->n_cells_local * p->n_types * sizeof(uint16_t));
    p->fields_unified = (uint16_t *) malloc((p->n_cells_local * p->n_types) * sizeof(uint16_t));
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

    // Set all values to zero
    p->num_all_beads = 0;
    p->num_all_beads_local = 0;

    // Determine number of  different bead types
    for (uint64_t j = 0; j < p->n_polymers; j++)
        {                       /*Loop over polymers */
            const unsigned int N = p->poly_arch[p->poly_type_offset[p->polymers[j].type]];
            p->num_all_beads_local += N;
        }
#if ( ENABLE_MPI == 1 )
    // Share p->num_all_beads
    MPI_Allreduce(&(p->num_all_beads_local), &(p->num_all_beads), 1, MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);
#else
    p->num_all_beads = p->num_all_beads_local;
#endif                          //ENABLE_MPI

    // Check if uint16_t density field is enough
    soma_scalar_t check_short = p->num_all_beads / p->n_cells;

    if (check_short > ((USHRT_MAX / 100) * 95))
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

    p->n_accessible_cells = ncells;
    // Loop to calculate scaling parameter
    // Note the *= the field is initialized with the density weights in read_hdf5_config.
    // default value = 1.
    for (unsigned int i = 0; i < p->n_types; i++)
        p->field_scaling_type[i] *= (ncells / ((soma_scalar_t) p->num_all_beads));

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

    p->sd.Ndoc = 0;
    p->sd.simdoc = NULL;
    p->sd.simdoc_internal = NULL;
    if (p->args.coord_file_arg != NULL) //This is only set to NULL if convert is used.
        {
            init_self_documentation(p, p->args.coord_file_arg, &(p->sd));
        }
    if (p->info_MPI.world_rank == 0)
        print_self_documentation(&(p->sd), stdout);

    int ret = 0;
    if (p->args.coord_file_arg != NULL) //Is it a full init Phase?
        ret = init_ana(p, p->args.ana_file_arg, p->args.coord_file_arg);
    else
        ret = init_ana(p, NULL, NULL);

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
#pragma acc enter data copyin(p->fields_unified[0:(p->n_types*p->n_cells_local)])
#pragma acc enter data copyin(p->old_fields_unified[0:(p->n_types*p->n_cells_local)])
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
#ifdef ENABLE_MPI_CUDA
    //in this case also copy in the buffers:

#pragma acc enter data copyin(p->left_tmp_buffer[0:(p->args.domain_buffer_arg*p->ny*p->nz) ])
#pragma acc enter data copyin(p->right_tmp_buffer[0:(p->args.domain_buffer_arg*p->ny*p->nz)])
#endif                          //ENABLE_MPI_CUDA
#endif                          //_OPENACC

    copyin_poly_conversion(p);
    copyin_mono_conversion(p);
    copyin_mobility(p);
    switch (p->args.pseudo_random_number_generator_arg)
        {
        case pseudo_random_number_generator__NULL:
            break;
        case pseudo_random_number_generator_arg_PCG32:
            break;
        case pseudo_random_number_generator_arg_MT:
            copyin_soma_memory(&(p->rh.mt));
            break;
        case pseudo_random_number_generator_arg_TT800:
            copyin_soma_memory(&(p->rh.tt800));
            break;
        }
    copyin_polymer_heavy(p);

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
#ifdef ENABLE_MPI_CUDA
#pragma acc exit data copyout(p->left_tmp_buffer[0:p->args.domain_buffer_arg*p->ny*p->nz])
#pragma acc exit data copyout(p->right_tmp_buffer[0:p->args.domain_buffer_arg*p->ny*p->nz])
#endif                          //ENABLE_MPI_CUDA

#pragma acc exit data copyout(p->polymers[0:p->n_polymers_storage])
    //Use here the delete to not overwrite stuff, which only changed on CPU
#pragma acc exit data delete(p[0:1])
#endif                          //_OPENACC

    copyout_poly_conversion(p);
    copyout_mono_conversion(p);
    copyout_mobility(p);
    switch (p->args.pseudo_random_number_generator_arg)
        {
        case pseudo_random_number_generator__NULL:
            break;
        case pseudo_random_number_generator_arg_PCG32:
            break;
        case pseudo_random_number_generator_arg_MT:
            copyout_soma_memory(&(p->rh.mt));
            break;
        case pseudo_random_number_generator_arg_TT800:
            copyout_soma_memory(&(p->rh.tt800));
            break;
        }
    copyout_polymer_heavy(p);

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

    free(p->polymers);
    free(p->poly_type_offset);
    free(p->poly_arch);
    free(p->xn);

    if (p->area51 != NULL)
        free(p->area51);

    if (p->external_field_unified != NULL)
        free(p->external_field_unified);

    if (p->umbrella_field != NULL)
        free(p->umbrella_field);

    free_poly_conversion(p);
    free_mono_conversion(p);
    free_mobility(p);

    free_self_documentation(&(p->sd));
    switch (p->args.pseudo_random_number_generator_arg)
        {
        case pseudo_random_number_generator__NULL:
            break;
        case pseudo_random_number_generator_arg_PCG32:
            break;
        case pseudo_random_number_generator_arg_MT:
            free_soma_memory(&(p->rh.mt));
            break;
        case pseudo_random_number_generator_arg_TT800:
            free_soma_memory(&(p->rh.tt800));
            break;
        }
    free_polymer_heavy(p);

    close_ana(&(p->ana_info));
    return 0;
}

int update_self_phase(Phase * const p, int rng_update_flag)
{
    static unsigned int last_time_call = 0;
    if (last_time_call == 0 || p->time > last_time_call)
        last_time_call = p->time;
    else                        //Quick exit, because the property has already been calculated for the time step.
        return 1;

    // Not pointer members are expected to not change on device
#pragma acc update self(p->xn[0:p->n_types*p->n_types])
#pragma acc update self(p->polymers[0:p->n_polymers])

    update_self_polymer_heavy(p, rng_update_flag);
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

    update_self_poly_conversion(p);
    update_self_mono_conversion(p);
    update_self_mobility(p);
    if (rng_update_flag)
        {
            switch (p->args.pseudo_random_number_generator_arg)
                {
                case pseudo_random_number_generator__NULL:
                    break;
                case pseudo_random_number_generator_arg_PCG32:
                    break;
                case pseudo_random_number_generator_arg_MT:
                    update_self_soma_memory(&(p->rh.mt));
                    break;
                case pseudo_random_number_generator_arg_TT800:
                    update_self_soma_memory(&(p->rh.tt800));
                    break;
                }
        }

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
