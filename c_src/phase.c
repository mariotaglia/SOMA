/* Copyright (C) 2016-2018 Ludwig Schneider

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

int init_phase(struct Phase * const p)
    {
    print_version(p->info_MPI.world_rank);
    p->start_time = p->time;
    p->start_clock = time(NULL);
    p->n_accepts = 0;
    p->n_moves =0;
    p->n_tries_cm = 0;
    p->n_acc_cm = 0;
    p->end_mono = NULL;
    p->tps_elapsed_time = 1; //Bla default, bigger 0
    p->tps_elapsed_steps = 1; //Bla default, bigger 0

    uint64_t n_polymer_offset;
    MPI_Scan( &(p->n_polymers), &n_polymer_offset, 1,MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);
    n_polymer_offset -= p->n_polymers;

    for(uint64_t i=0; i < p->n_polymers;i++)
        {
        p->polymers[i].set_states = NULL;
        p->polymers[i].set_permutation = NULL;

        allocate_rng_state(&(p->polymers[i].poly_state), p->args.pseudo_random_number_generator_arg);
        seed_rng_state(&(p->polymers[i].poly_state), p->args.rng_seed_arg,
                       i+n_polymer_offset, p->args.pseudo_random_number_generator_arg);
        }

    // Max safe move distance
    p->max_safe_jump = p->Lx/p->nx < p->Ly / p->ny ? p->Lx/p->nx : p->Ly / p->ny;
    p->max_safe_jump = p->max_safe_jump < p->Lz / p->nz ? p->max_safe_jump : p->Lz/p->nz;
    p->max_safe_jump *= 0.95;

    // Reference Harmonic Spring Cste
    const soma_scalar_t harmonic_spring_Cste =
        1.0 / sqrt(3.0 * (p->reference_Nbeads - 1.0));
    //Reference energy scale for harmonic springs.
    p->harmonic_normb =
        1.0 / (2.0 * harmonic_spring_Cste * harmonic_spring_Cste);

    uint64_t n_polymers_global_sum;
    MPI_Allreduce(&(p->n_polymers), &n_polymers_global_sum, 1,
                  MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);
    assert(p->n_polymers_global == n_polymers_global_sum);

    p->n_cells = p->nx * p->ny * p->nz;
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    p->local_nx_low = (p->nx/p->args.N_domains_arg * my_domain) - p->args.domain_buffer_arg;
    p->local_nx_high = (p->nx/p->args.N_domains_arg * my_domain) + p->args.domain_buffer_arg;
    p->n_cells_local = (p->local_nx_high - p->local_nx_low) * p->ny * p->nz;

    //Check if it is a valid domain decomposition
    if( p->args.N_domains_arg*2*p->args.domain_buffer_arg > (int)p->nx)
        {
        fprintf(stderr,"ERROR: invalid domain decomposition. %s:%d\n",__FILE__,__LINE__);
        fprintf(stderr,"\t N(N_domains)= %d\tb(domain_buffer)= %d\t nx= %d\n",p->args.N_domains_arg,p->args.domain_buffer_arg,p->nx);
        fprintf(stderr,"\t N * 2 * b <= nx \t not fulfilled\n");
        return -2;
        }

    //Allocate Fields

    p->fields_unified =     (uint16_t *) malloc(p->n_cells_local*p->n_types*sizeof(uint16_t));
    if (p->fields_unified == NULL) {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
        }

    p->old_fields_unified =     (uint16_t *) malloc(p->n_cells_local*p->n_types*sizeof(uint16_t));
    if (p->old_fields_unified == NULL) {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
        }
    p->fields_32 = (uint32_t*)malloc(p->n_types * p->n_cells_local * sizeof(uint32_t));
    if (p->fields_32 == NULL) {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
        }
    p->omega_field_unified = (soma_scalar_t*)malloc(p->n_cells_local * p->n_types * sizeof(soma_scalar_t));
    if (p->omega_field_unified == NULL) {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
        }


    p->tempfield =
        (soma_scalar_t *) malloc(p->n_cells_local * sizeof(soma_scalar_t));
    if (p->tempfield == NULL) {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
        }

    p->num_bead_type =
        (uint64_t *) malloc(p->n_types * sizeof(uint64_t));
    if (p->num_bead_type == NULL) {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
        }
    p->num_bead_type_local =
        (uint64_t *) malloc(p->n_types * sizeof(uint64_t));
    if (p->num_bead_type_local == NULL) {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
        }
    p->field_scaling_type = (soma_scalar_t *) malloc(p->n_types * sizeof(soma_scalar_t));
    if (p->field_scaling_type == NULL) {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
        }

    // Set all values to zero

    p->num_all_beads = 0;
    p->num_all_beads_local = 0;
    for (unsigned int i = 0; i < p->n_types; i++)
        p->num_bead_type_local[i] = 0;
    // Determine number of  different bead types
    for (uint64_t j = 0; j < p->n_polymers; j++) {      /*Loop over polymers */
        const unsigned int N = p->poly_arch[ p->poly_type_offset[p->polymers[j].type] ];
        for (unsigned int k = 0; k < N; k++) {  /*Loop over monomers */
            const unsigned int type = get_particle_type(
                p->poly_arch[ p->poly_type_offset[p->polymers[j].type]+1+k]);
            p->num_bead_type_local[type] += 1;
            p->num_all_beads_local += 1;
            }
        }
    // Share p->num_all_beads
    MPI_Allreduce(&(p->num_all_beads_local), &(p->num_all_beads), 1,
                  MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);

    // Share p->num_bead_type
    for (unsigned int i = 0; i < p->n_types; i++) {
        MPI_Allreduce(&(p->num_bead_type_local[i]), &(p->num_bead_type[i]),
                      1, MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);
        }
    // Check if uint16_t density field is enough
    soma_scalar_t check_short = p->num_all_beads/p->n_cells;

    if (check_short > ((USHRT_MAX/100)*95)){
        fprintf(stderr, "ERROR: Cell-density above 95 Percent of USHRT_MAX\n");
        return -1;
        }
    // setting the geometry field, or if it's already initialized measure the free space
    uint64_t ncells = (p->nx/p->args.N_domains_arg)*p->ny*p->nz;
    if( p->info_MPI.domain_rank == 0) //Only domain root calculats something else than 0
        {
        if (p->area51 != NULL) {
            // substract the number of non free cells for the correct density scaling
            for (uint64_t i = 0; i < p->n_cells; i++)
                if ( p->area51[i] > 0 ) ncells--;
            }
        ncells = 0;
        }
    MPI_Allreduce( MPI_IN_PLACE, &ncells, 1, MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);

    // Loop to calculate scaling parameter
    for (unsigned int i = 0; i < p->n_types; i++)
        p->field_scaling_type[i] =
            (ncells / ((soma_scalar_t) p->num_all_beads));
    // Info for Ulrich: programm does take excluded volume into account now!
    p->n_accepts = 0;
    p->n_moves = 0;

    p->R = (soma_scalar_t*) malloc( p->n_types * sizeof(soma_scalar_t));
    for (unsigned int i = 0; i < p->n_types; i++){
        //! \todo kBT required.
        p->R[i] = sqrt( p->A[i] * 2);
        }

    // initialize inverse simulation cell parameters
    p->iLx = 1.0/p->Lx;
    p->iLy = 1.0/p->Ly;
    p->iLz = 1.0/p->Lz;

    p->sets = NULL; // Default init of the sets
    p->max_set_members = 0;
    if( p->args.iteration_alg_arg == iteration_alg_arg_SET)
        generate_independet_sets(p);


    init_autotuner(&(p->mc_autotuner));
    init_autotuner(&(p->cm_mc_autotuner));

    copyin_phase(p);
    // call update_fields routine

    if(p->bead_data_read)
        {
        update_density_fields(p);
        memcpy(p->old_fields_unified, p->fields_unified, p->n_cells_local*p->n_types*sizeof(uint16_t));
        }

    int ret=0;
    if(p->args.coord_file_arg != NULL) //Is it a full init Phase?
        ret = init_ana(p,p->args.ana_file_arg,p->args.coord_file_arg);
    else
        ret = init_ana(p,NULL,NULL);

    return ret;
    }

int copyin_phase(struct Phase*const p)
    {
#ifdef _OPENACC
#pragma acc enter data copyin(p[0:1])
#pragma acc enter data copyin(p->xn[0:p->n_types][0:p->n_types])
#pragma acc enter data copyin(p->polymers[0:p->n_polymers_storage])
#pragma acc enter data copyin(p->fields_unified[0:p->n_types*p->n_cells_local])
#pragma acc enter data copyin(p->old_fields_unified[0:p->n_types*p->n_cells_local])
#pragma acc enter data copyin(p->fields_32[0:p->n_types*p->n_cells_local])
    if (p->area51 != NULL){
#pragma acc enter data copyin(p->area51[0:p->n_cells_local])
        }
#pragma acc enter data copyin(p->omega_field_unified[0:p->n_cells_local*p->n_types])
    if (p->external_field_unified != NULL){
#pragma acc enter data copyin(p->external_field_unified[0:p->n_cells_local*p->n_types])
        }
    if (p->umbrella_field != NULL){
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

    if(p->cm_a != NULL)
        {
#pragma acc enter data copyin(p->cm_a[0:p->n_poly_type])
        }
    if( p->sets != NULL)
        {
#pragma acc enter data copyin(p->sets[0:p->n_poly_type])
        for(unsigned int i=0; i < p->n_poly_type; i++)
            {
#pragma acc enter data copyin(p->sets[i].set_length[0:p->sets[i].n_sets])
#pragma acc enter data copyin(p->sets[i].sets[0:p->sets[i].n_sets*p->sets[i].max_member])
            }
        }
    for(uint64_t i=0; i < p->n_polymers; i++)
        {
        Polymer*const poly = &(p->polymers[i]);
        copyin_polymer(p, poly);
        }
    return 0;
#else
    return p->n_polymers*0+1;
#endif//_OPENACC
    }

int copyout_phase(struct Phase*const p)
    {
#ifdef _OPENACC
#pragma acc exit data delete(p[0:1])
#pragma acc exit data delete(p->xn[0:p->n_types][0:p->n_types])
#pragma acc exit data delete(p->polymers[0:p->n_polymers_storage])
#pragma acc exit data delete(p->fields_unified[0:p->n_types*p->n_cells_local])
#pragma acc exit data delete(p->old_fields_unified[0:p->n_types*p->n_cells_local])
#pragma acc exit data delete(p->fields_32[0:p->n_types*p->n_cells_local])
    if (p->area51 != NULL){
#pragma acc exit data delete(p->area51[0:p->n_cells_local])
        }
#pragma acc exit data delete(p->omega_field_unified[0:p->n_cells_local*p->n_types])
    if (p->external_field_unified != NULL){
#pragma acc exit data delete(p->external_field_unified[0:p->n_cells_local*p->n_types])
        }
    if (p->umbrella_field != NULL){
#pragma acc exit data delete(p->umbrella_field[0:p->n_cells_local*p->n_types])
        }
#pragma acc exit data delete(p->tempfield[0:p->n_cells_local])
#pragma acc exit data delete(p->num_bead_type[0:p->n_types])
#pragma acc exit data delete(p->num_bead_type_local[0:p->n_types])
#pragma acc exit data delete(p->A[0:p->n_types])
#pragma acc exit data delete(p->R[0:p->n_types])
#pragma acc exit data delete(p->field_scaling_type[0:p->n_types])
#pragma acc exit data delete(p->k_umbrella[0:p->n_types])
#pragma acc exit data delete(p->poly_type_offset[0:p->n_poly_type])
#pragma acc exit data delete(p->poly_arch[0:p->poly_arch_length])

    if(p->cm_a != NULL)
        {
#pragma acc exit data delete(p->cm_a[0:p->n_poly_type])
        }
    if( p->sets != NULL)
        {
#pragma acc exit data delete(p->sets[0:p->n_poly_type])
        for(unsigned int i=0; i < p->n_poly_type; i++)
            {
#pragma acc exit data delete(p->sets[i].set_length[0:p->sets[i].n_sets])
#pragma acc exit data delete(p->sets[i].sets[0:p->sets[i].n_sets*p->sets[i].max_member])
            }
        }

    for(uint64_t i=0; i < p->n_polymers; i++)
        {
        Polymer*const poly = &(p->polymers[i]);
        copyout_polymer(p, poly);
        }
    return 0;
#else
    return p->n_polymers*0 +1;
#endif//_OPENACC
    }

int free_phase(struct Phase * const p)
    {
    copyout_phase(p);

    /* de-allocate fields */
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

    if(p->args.coord_file_arg != NULL) //Required option. So always true, if args used.
        cmdline_parser_free(&(p->args));

    if(p->sets != NULL)
        {
        for(unsigned int i=0; i < p->n_poly_type; i++)
            {
            free(p->sets[i].set_length);
            free(p->sets[i].sets);
            }
        free(p->sets);
        }

    /* free polymers */
    for (uint64_t i = 0; i < p->n_polymers; i++)
        {
        Polymer*const poly = &(p->polymers[i]);
        free_polymer(p, poly);
        }

    free(p->polymers);

    free(p->poly_type_offset);
    free(p->poly_arch);

    /* de allocate XN interaction matrix */
    for (unsigned int i = 0; i < p->n_types; i++)
        free(p->xn[i]);
    free(p->xn);

    if (p->area51 != NULL) {
        free(p->area51);
        }

    if (p->external_field_unified != NULL){
        free(p->external_field_unified);
        }

    if (p->umbrella_field != NULL){
        free(p->umbrella_field);
        }

    close_ana(&(p->ana_info));
    return 0;
    }

int update_self_phase(const Phase * const p)
    {
    static unsigned int last_time_call = 0;
    if (last_time_call == 0 || p->time > last_time_call)
        last_time_call = p->time;
    else                        //Quick exit, because the property has already been calculated for the time step.
        return 1;

    // Not pointer members are expected to not change on device

#pragma acc update self(p->xn[0:p->n_types][0:p->n_types])

    for(uint64_t i=0; i< p->n_polymers; i++)
        update_self_polymer(p, p->polymers+i);

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
    if(p->cm_a != NULL)
        {
#pragma acc update self(p->cm_a[0:p->n_poly_type])
        }

    //SETS are not updated to host

    return p->n_polymers*0+1;
    }
