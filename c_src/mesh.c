/* Copyright (C) 2016-2019 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016-2017 Marcel Langenberg
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

//! \file mesh.c
//! \brief Implementation of mesh.h

#include "mesh.h"
#include <stdbool.h>
#include <assert.h>
#include "mpiroutines.h"
#include "soma_config.h"

/* #pragma acc routine seq */
/* inline void increment_16_bit_uint(uint16_t*const ptr) */
/* { */
/*   int x = (((size_t)ptr) & 2)>>1 ; */
/*   uint32_t * ptr32 = (uint32_t *) (ptr-x) ; */
/* #pragma acc atomic update */
/*   *ptr32 += (1 << (x*16)) ; */
/* } */

/*! Function to communicate the density field accross all ranks.
   If a domain decomposition is used, the communication pattern is accordingly adopted.
   \param p Phase describing the system.
*/
void communicate_density_fields(const struct Phase *const p)
{
#if ( ENABLE_MPI == 1 )
    //Const cast!
    mpi_divergence((struct Phase * const)p);

    if (p->info_MPI.sim_size > 1)
        {
            if (p->args.N_domains_arg == 1)
                {
#ifndef ENABLE_MPI_CUDA
#pragma acc update self(p->fields_unified[0:p->n_cells_local*p->n_types])
                    MPI_Allreduce(MPI_IN_PLACE, p->fields_unified, p->n_cells_local * p->n_types, MPI_UINT16_T, MPI_SUM,
                                  p->info_MPI.SOMA_comm_sim);
#pragma acc update device(p->fields_unified[0:p->n_cells_local*p->n_types])
#else                           //ENABLE_MPI_CUDA
                    uint16_t *fields_unified = p->fields_unified;
#pragma acc host_data use_device(fields_unified)
                    {
                        MPI_Allreduce(MPI_IN_PLACE, fields_unified, p->n_cells_local * p->n_types, MPI_UINT16_T,
                                      MPI_SUM, p->info_MPI.SOMA_comm_sim);
                    }
#pragma acc update self(p->fields_unified[0:p->n_cells_local*p->n_types])
#endif                          //ENABLE_MPI_CUDA
                }
            else                //Communication for domain decomposition
                {
                    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;

#pragma acc update self(p->fields_unified[0:p->n_cells_local*p->n_types])
                    //Sum up all values of a single domain to the root domain
                    if (p->info_MPI.domain_rank == 0)
                        MPI_Reduce(MPI_IN_PLACE, p->fields_unified, p->n_cells_local * p->n_types, MPI_UINT16_T,
                                   MPI_SUM, 0, p->info_MPI.SOMA_comm_domain);
                    else
                        MPI_Reduce(p->fields_unified, NULL, p->n_cells_local * p->n_types, MPI_UINT16_T, MPI_SUM, 0,
                                   p->info_MPI.SOMA_comm_domain);
                    if (p->info_MPI.domain_rank == 0)
                        {
                            assert(p->left_tmp_buffer != NULL);
                            assert(p->right_tmp_buffer != NULL);
                            const unsigned int ghost_buffer_size = p->args.domain_buffer_arg * p->ny * p->nz;
                            const int left_rank =
                                (((my_domain - 1) +
                                  p->args.N_domains_arg) % p->args.N_domains_arg) * p->info_MPI.domain_size +
                                p->info_MPI.domain_rank;
                            const int right_rank =
                                (((my_domain + 1) +
                                  p->args.N_domains_arg) % p->args.N_domains_arg) * p->info_MPI.domain_size +
                                p->info_MPI.domain_rank;
                            MPI_Request req[4];
                            MPI_Status stat[4];

                            //Loop over type, because of the memory layout [type][x][y][z] -> x is not slowest moving dimension
                            for (unsigned int type = 0; type < p->n_types; type++)
                                {
                                    //Send buffer to right, recv from left
                                    MPI_Isend(p->fields_unified + (p->n_cells_local - ghost_buffer_size) +
                                              type * p->n_cells_local, ghost_buffer_size, MPI_UINT16_T, right_rank, 0,
                                              p->info_MPI.SOMA_comm_sim, req + 0);
                                    MPI_Irecv(p->left_tmp_buffer, ghost_buffer_size, MPI_UINT16_T, left_rank, 0,
                                              p->info_MPI.SOMA_comm_sim, req + 1);
                                    //Send buffer to left recv from right
                                    MPI_Isend(p->fields_unified + type * p->n_cells_local,
                                              ghost_buffer_size, MPI_UINT16_T, left_rank, 1, p->info_MPI.SOMA_comm_sim,
                                              req + 2);
                                    MPI_Irecv(p->right_tmp_buffer, ghost_buffer_size, MPI_UINT16_T, right_rank, 1,
                                              p->info_MPI.SOMA_comm_sim, req + 3);

                                    MPI_Waitall(4, req, stat);

                                    //Add the recv values to main part
                                    for (unsigned int i = 0; i < ghost_buffer_size; i++)
                                        {
                                            p->fields_unified[ghost_buffer_size + i + type * p->n_cells_local] +=
                                                p->left_tmp_buffer[i];
                                            p->fields_unified[p->n_cells_local - 2 * ghost_buffer_size + i +
                                                              type * p->n_cells_local] += p->right_tmp_buffer[i];
                                        }

                                    //Update the buffers of the neighbors
                                    MPI_Isend(p->fields_unified + (p->n_cells_local - 2 * ghost_buffer_size) +
                                              type * p->n_cells_local, ghost_buffer_size, MPI_UINT16_T, right_rank, 2,
                                              p->info_MPI.SOMA_comm_sim, req + 0);
                                    MPI_Irecv(p->fields_unified + type * p->n_cells_local, ghost_buffer_size,
                                              MPI_UINT16_T, left_rank, 2, p->info_MPI.SOMA_comm_sim, req + 1);

                                    MPI_Isend(p->fields_unified + ghost_buffer_size + type * p->n_cells_local,
                                              ghost_buffer_size, MPI_UINT16_T, left_rank, 3, p->info_MPI.SOMA_comm_sim,
                                              req + 2);
                                    MPI_Irecv(p->fields_unified + (p->n_cells_local - ghost_buffer_size) +
                                              type * p->n_cells_local, ghost_buffer_size, MPI_UINT16_T, right_rank, 3,
                                              p->info_MPI.SOMA_comm_sim, req + 3);
                                    MPI_Waitall(4, req, stat);
                                }
                        }

                    //Update all domain ranks with the results of the root domain rank
                    MPI_Bcast(p->fields_unified, p->n_cells_local * p->n_types, MPI_UINT16_T, 0,
                              p->info_MPI.SOMA_comm_domain);
#pragma acc update device(p->fields_unified[0:p->n_cells_local*p->n_types])

                    //Avoid false loadbalance
                    MPI_Barrier(p->info_MPI.SOMA_comm_domain);
                }
        }
#endif                          //ENABLE_MPI
}

int update_density_fields(const struct Phase *const p)
{
    static unsigned int last_time_call = 0;
    if (last_time_call == 0 || p->time > last_time_call)
        last_time_call = p->time;
    else                        //Quick exit, because the property has already been calculated for the time step.
        return 0;

    int error_flags[1] = { 0 }; //error_flag[0] indicates domain errors
#pragma acc enter data copyin(error_flags[0:1])

    const uint64_t n_indices = p->n_types * p->n_cells_local;

#pragma acc parallel loop independent present(p[0:1])
#pragma omp parallel for
    for (uint64_t index = 0; index < n_indices; index++)        /*Loop over all fields according to monotype */
        p->fields_32[index] = 0;
    const uint64_t n_polymers = p->n_polymers;

#pragma acc parallel loop gang num_gangs(n_polymers) vector_length(128) present(p[0:1])
#pragma omp parallel for
    for (uint64_t i = 0; i < n_polymers; i++)
        {                       /*Loop over polymers */
            const unsigned int N = p->poly_arch[p->poly_type_offset[p->polymers[i].type]];
#pragma acc loop vector
            for (unsigned int j = 0; j < N; j++)
                {               /*Loop over monomers */
                    const unsigned int monotype =
                        get_particle_type(p->poly_arch[p->poly_type_offset[p->polymers[i].type] + 1 + j]);

                    const unsigned int index = coord_to_index_unified(p, p->polymers[i].beads[j].x,
                                                                      p->polymers[i].beads[j].y,
                                                                      p->polymers[i].beads[j].z, monotype);
                    if (index < p->n_cells_local * p->n_types)
                        {
#pragma acc atomic update
#pragma omp atomic
                            p->fields_32[index] += 1;
                        }
                    else
                        {
                            error_flags[0] = i + 1;
                        }
                }
        }

#pragma acc exit data copyout(error_flags[0:1])
    if (error_flags[0] != 0)
        {
            fprintf(stderr, "ERROR: Domain error. %d world-rank %d"
                    " A particle has left the buffer domain."
                    " Restart your simulation with larger buffers. %s:%d\n", error_flags[0], p->info_MPI.world_rank,
                    __FILE__, __LINE__);
            return error_flags[0];
        }

#pragma acc parallel loop independent present(p[0:1])
    for (uint64_t index = 0; index < n_indices; index++)
        p->fields_unified[index] = p->fields_32[index];

    /*Share the densityfields -> needed because Hamiltonian may not only quadratic order */
    communicate_density_fields(p);

    /* Calculate the added up densities */

    /*Use first type to initialize the fields-> saves set zero routine */
    soma_scalar_t rescale_density = p->field_scaling_type[0];
#pragma acc parallel loop present(p[0:1])
#pragma omp parallel for
    for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
        p->tempfield[cell] = rescale_density * p->fields_unified[cell];

    for (unsigned int T_types = 1; T_types < p->n_types; T_types++)
        {
            rescale_density = p->field_scaling_type[T_types];
#pragma acc parallel loop present(p[0:1])
#pragma omp parallel for
            for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
                p->tempfield[cell] += rescale_density * p->fields_unified[T_types * p->n_cells_local + cell];
            /*!\todo p->ncells as a temporary variable */
        }
//#pragma acc update self(p->tempfield[0:p->n_cells_local])
    return 0;
}

void update_omega_fields(const struct Phase *const p)
{
    static unsigned int last_time_call = 0;
    if (last_time_call == 0 || p->time > last_time_call)
        last_time_call = p->time;
    else                        //Quick exit, because the property has already been calculated for the time step.
        return;
    switch (p->hamiltonian)
        {
        case SCMF0:
            update_omega_fields_scmf0(p);
            break;
        case SCMF1:
            update_omega_fields_scmf1(p);
            break;
        case SCMF2;
            update_omega_fields_scmf2(p);
            break;
        default:
            fprintf(stderr, "ERROR: %s:%d Unkown hamiltonian specified %d.\n", __FILE__, __LINE__, p->hamiltonian);
        }
}

//! Set the self interaction terms to the omega fields. This
//! includes the intercation with the external field and the umbrella
//! potential.
//! \private Helper function
//! \param p Phase of the system to init the omega fields
void self_omega_field(const struct Phase *const p)
{
    // soma_scalar_t density=0;
    // for(int c=0;c<(p->n_types * p->n_cells_local);c++)
    //   density+=p->fields_unified[c];
    // density=density/(p->Lx*p->Ly*p->Lz);

    /*Densityfields are shorts and unscaled according to bead_type */
    /*Tempfields save the complete densities and remain as type soma_scalar_t -> used for insothermal Compressibility */

    const soma_scalar_t inverse_refbeads = 1.0 / p->reference_Nbeads;
    //Add cos and sin series with time dependency
    soma_scalar_t external_field_time = 0;

    if (p->time == 0 || p->serie_length == 0)
        {                       //to fix serie_length seg fault
            external_field_time = 1;
        }
    else
        {
            for (unsigned int serie_index = 0; serie_index < p->serie_length; serie_index++)
                {
                    external_field_time +=
                        p->cos_serie[serie_index] * cos(2 * M_PI * serie_index / p->period * p->time);
                    external_field_time +=
                        p->sin_serie[serie_index] * sin(2 * M_PI * serie_index / p->period * p->time);
                }
        }
    // Compressibility part + external fields
    for (unsigned int T_types = 0; T_types < p->n_types; T_types++)     /*Loop over all fields according to monotype */
        {
#pragma acc parallel loop present(p[:1])
#pragma omp parallel for
            for (uint64_t cell = 0; cell < p->n_cells_local; cell++)    /*Loop over all cells, max number of cells is product of nx, ny,nz */
                {
                    switch(p->hamiltonian){
                        case SCMF0://Intentionally falls through
                        case SCMF1:

                            p->omega_field_unified[cell + T_types * p->n_cells_local] =
                                inverse_refbeads * (p->xn[T_types * p->n_types + T_types] * (p->tempfield[cell] - 1.0));
                            break;
                        case SCMF2:;
                            const soma_scalar_t field = p->fields_unified[cell]; 
                            const soma_scalar_t wn = p->wn[(T_types * p->n_types + T_types) * p->n_types + T_types];
                            const soma_scalar_t xn = p->xn[T_types * p->n_types + T_types];
                            p->omega_field_unified[cell + T_types * p->n_cells_local] = inverse_refbeads * ((wn * field + xn) * field);
                            break;
                    //Should be finee................................................................................................................................................................................. 
                    }
                    /* the external field is defined such that the energy of a
                       chain of refbeads in this field is x k_B T, thus the
                       normalization per bead */
                    if (p->external_field_unified != NULL)
                        {
                            p->omega_field_unified[cell + T_types * p->n_cells_local] +=
                                inverse_refbeads * p->external_field_unified[cell +
                                                                             T_types * p->n_cells_local] *
                                external_field_time;
                        }
                    //umbrella part
                    if (p->umbrella_field != NULL)
                        {
                            p->omega_field_unified[cell + T_types * p->n_cells_local] +=
                                -inverse_refbeads * p->field_scaling_type[T_types] * p->k_umbrella[T_types] *
                                (p->umbrella_field[cell + T_types * p->n_cells_local] -
                                 p->fields_unified[cell + T_types * p->n_cells_local]);
                        }
                }
        }
}

//! Add the pair interactions to the omega fields via the SCMF0 hamiltonian.
//! \private Helper function
//! \param p Phase of the system to add the omega fields
void add_pair_omega_fields_scmf0(const struct Phase *const p)
{
    const soma_scalar_t inverse_refbeads = 1.0 / p->reference_Nbeads;

    // XN part

    for (unsigned int T_types = 0; T_types < p->n_types; T_types++)
        {                       /*Loop over all fields according to monotype */
            for (unsigned int S_types = T_types + 1; S_types < p->n_types; S_types++)
                {
                    // precalculate the normalization for this type combination
                    soma_scalar_t dnorm = -0.5 * inverse_refbeads * p->xn[T_types * p->n_types + S_types];
#pragma acc parallel loop present(p[:1])
#pragma omp parallel for
                    for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
                        {
                            soma_scalar_t interaction =
                                dnorm * (p->field_scaling_type[T_types] *
                                         p->fields_unified[cell + T_types * p->n_cells_local] -
                                         p->field_scaling_type[S_types] * p->fields_unified[cell +
                                                                                            S_types *
                                                                                            p->n_cells_local]);
                            p->omega_field_unified[cell + T_types * p->n_cells_local] += interaction;
                            p->omega_field_unified[cell + S_types * p->n_cells_local] -= interaction;
                        }
                }
        }
}

//! Add the pair interactions to the omega fields via the SCMF1 hamiltonian.
//! \private Helper function
//! \param p Phase of the system to add the omega fields
void add_pair_omega_fields_scmf1(const struct Phase *const p)
{
    const soma_scalar_t inverse_refbeads = 1.0 / p->reference_Nbeads;

    // XN part

    for (unsigned int T_types = 0; T_types < p->n_types; T_types++)
        {                       /*Loop over all fields according to monotype */
            for (unsigned int S_types = T_types + 1; S_types < p->n_types; S_types++)
                {
#pragma acc parallel loop present(p[:1])
#pragma omp parallel for
                    for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
                        {
                            const soma_scalar_t normT = inverse_refbeads * p->xn[T_types * p->n_types + S_types];
                            const soma_scalar_t rhoS =
                                p->fields_unified[cell + S_types * p->n_cells_local] * p->field_scaling_type[S_types];
                            const soma_scalar_t normS = inverse_refbeads * p->xn[S_types * p->n_types + T_types];
                            const soma_scalar_t rhoT =
                                p->fields_unified[cell + T_types * p->n_cells_local] * p->field_scaling_type[T_types];
                            p->omega_field_unified[cell + T_types * p->n_cells_local] += normT * rhoS;
                            p->omega_field_unified[cell + S_types * p->n_cells_local] += normS * rhoT;
                        }
                }
        }
}

//! Add the triple interactions to the omega fields via the SCMF2 hamiltonian.
//! \private Helper function
//! \param p Phase of the system to add the omega fields
void add_triple_omega_fields_scmf2(const struct Phase *const p)
{
    const soma_scalar_t inverse_refbeads = 1.0 / p->reference_Nbeads;

    // XN part


    for (unsigned int T_types = 0; T_types < p->n_types; T_types++)
        {                       /*Loop over all fields according to monotype */
            for (unsigned int S_types = T_types; S_types < p->n_types; S_types++)
                {
                    for (unsigned int R_types = S_types + 1; R_types < p->n_types; R_types++)
                        {
#pragma acc parallel loop present(p[:1])
#pragma omp parallel for
                            for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
                                {
                                    const soma_scalar_t norm = inverse_refbeads * p->wn[(T_types * p->n_types + S_types) * p->n_types + R_types];
                                    const soma_scalar_t rhoS =
                                        p->fields_unified[cell + S_types * p->n_cells_local] * p->field_scaling_type[S_types];
                                    const soma_scalar_t rhoT =
                                        p->fields_unified[cell + T_types * p->n_cells_local] * p->field_scaling_type[T_types];
                                    const soma_scalar_t rhoR =
                                        p->fields_unified[cell + R_types * p->n_cells_local] * p->field_scaling_type[R_types];
                                    p->omega_field_unified[cell + S_types * p->n_cells_local] += norm * rhoT * rhoR;
                                    p->omega_field_unified[cell + T_types * p->n_cells_local] += norm * rhoS * rhoR;
                                    p->omega_field_unified[cell + R_types * p->n_cells_local] += norm * rhoS * rhoT;
                                }
                        }
                }
        }
}

void update_omega_fields_scmf0(const struct Phase *const p)
{
    self_omega_field(p);
    add_pair_omega_fields_scmf0(p);
}

void update_omega_fields_scmf1(const struct Phase *const p)
{
    self_omega_field(p);
    add_pair_omega_fields_scmf1(p);
}

void update_omega_fields_scmf2(const struct Phase *const p)
{
    self_omega_field(p);
    add_pair_omega_fields_scmf1(p);
    add_triple_omega_fields_scmf2(p);
}
