/* Copyright (C) 2016-2021 Ludwig Schneider
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


int mod(int a, int b); // modulus

/* #pragma acc routine seq */
/* inline void increment_16_bit_uint(uint16_t*const ptr) */
/* { */
/*   int x = (((size_t)ptr) & 2)>>1 ; */
/*   uint32_t * ptr32 = (uint32_t *) (ptr-x) ; */
/* #pragma acc atomic update */
/*   *ptr32 += (1 << (x*16)) ; */
/* } */

void communicate_simple(const struct Phase *const p)
{
#ifndef ENABLE_MPI_CUDA
#pragma acc update self(p->fields_unified[0:p->n_cells_local*p->n_types])
    MPI_Allreduce(MPI_IN_PLACE, p->fields_unified, p->n_cells_local * p->n_types, MPI_UINT16_T, MPI_SUM,
                  p->info_MPI.SOMA_comm_sim);
#pragma acc update device(p->fields_unified[0:p->n_cells_local*p->n_types])
#else                           //ENABLE_MPI_CUDA
#ifdef ENABLE_NCCL
    if (p->info_MPI.gpu_id >= 0)
        {
            uint32_t *fields_32 = p->fields_32;
#pragma acc host_data use_device(fields_32)
            {
                // NCCL does not support unit16 so far, hence we are using
                ncclAllReduce(fields_32, fields_32, p->n_cells_local * p->n_types, ncclUint32,
                              ncclSum, p->info_MPI.SOMA_nccl_sim, acc_get_cuda_stream(acc_get_default_async()));
            }
#pragma acc wait
            copy_density_32_to_16(p);
        }
    else
        {
            MPI_Allreduce(MPI_IN_PLACE, p->fields_unified, p->n_cells_local * p->n_types, MPI_UINT16_T, MPI_SUM,
                          p->info_MPI.SOMA_comm_sim);
        }

#else                           //ENABLE_NCCL
    uint16_t *fields_unified = p->fields_unified;
#pragma acc host_data use_device(fields_unified)
    {
        MPI_Allreduce(MPI_IN_PLACE, fields_unified, p->n_cells_local * p->n_types, MPI_UINT16_T,
                      MPI_SUM, p->info_MPI.SOMA_comm_sim);
    }
#endif                          //ENABLE_NCCL

    //#pragma acc update self(p->fields_unified[0:p->n_cells_local*p->n_types])
    // No update of data on cpu required since it is only needed for analytics (updated only when needed)
#endif                          //ENABLE_MPI_CUDA
}

void communicate_domain_decomposition(const struct Phase *const p)
{
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;

#ifdef ENABLE_NCCL
    if (p->info_MPI.gpu_id >= 0)
        {
            uint32_t *fields_32 = p->fields_32;
#pragma acc host_data use_device(fields_32)
            {
                ncclReduce(fields_32, fields_32, p->n_cells_local * p->n_types, ncclUint32, ncclSum,
                           0, p->info_MPI.SOMA_nccl_domain, acc_get_cuda_stream(acc_get_default_async()));
            }
#pragma acc wait
            copy_density_32_to_16(p);
        }
    else
        {
            if (p->info_MPI.domain_rank == 0)
                MPI_Reduce(MPI_IN_PLACE, p->fields_unified, p->n_cells_local * p->n_types, MPI_UINT16_T,
                           MPI_SUM, 0, p->info_MPI.SOMA_comm_domain);
            else
                MPI_Reduce(p->fields_unified, NULL, p->n_cells_local * p->n_types, MPI_UINT16_T, MPI_SUM, 0,
                           p->info_MPI.SOMA_comm_domain);
        }
#endif                          //ENABLE_NCCL

    uint16_t *fields_unified = p->fields_unified;
    uint16_t *left_tmp_buffer = p->left_tmp_buffer;
    uint16_t *right_tmp_buffer = p->right_tmp_buffer;
#ifndef ENABLE_MPI_CUDA
#pragma acc update self(p->fields_unified[0:p->n_cells_local*p->n_types])
#else                           //ENABLE_MPI_CUDA
#pragma acc host_data use_device(fields_unified,left_tmp_buffer,right_tmp_buffer)
    {
#endif                          //ENABLE_MPI_CUDA

#ifndef ENABLE_NCCL
        //Sum up all values of a single domain to the root domain
        if (p->info_MPI.domain_rank == 0)
            MPI_Reduce(MPI_IN_PLACE, fields_unified, p->n_cells_local * p->n_types, MPI_UINT16_T,
                       MPI_SUM, 0, p->info_MPI.SOMA_comm_domain);
        else
            MPI_Reduce(fields_unified, NULL, p->n_cells_local * p->n_types, MPI_UINT16_T, MPI_SUM, 0,
                       p->info_MPI.SOMA_comm_domain);
#endif                          //ENABLE_NCCL

        if (p->info_MPI.domain_rank == 0)
            {
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
                        MPI_Isend(fields_unified + (p->n_cells_local - ghost_buffer_size) +
                                  type * p->n_cells_local, ghost_buffer_size, MPI_UINT16_T, right_rank,
                                  0, p->info_MPI.SOMA_comm_sim, req + 0);
                        MPI_Irecv(left_tmp_buffer, ghost_buffer_size, MPI_UINT16_T, left_rank, 0,
                                  p->info_MPI.SOMA_comm_sim, req + 1);
                        //Send buffer to left recv from right
                        MPI_Isend(fields_unified + type * p->n_cells_local,
                                  ghost_buffer_size, MPI_UINT16_T, left_rank, 1, p->info_MPI.SOMA_comm_sim, req + 2);
                        MPI_Irecv(right_tmp_buffer, ghost_buffer_size, MPI_UINT16_T, right_rank, 1,
                                  p->info_MPI.SOMA_comm_sim, req + 3);
                        MPI_Waitall(4, req, stat);

                        //Add the recv values to main part
#ifdef ENABLE_MPI_CUDA
#pragma acc parallel loop present(p[0:1])
#endif                          //ENABLE_MPI_CUDA
                        for (unsigned int i = 0; i < ghost_buffer_size; i++)
                            {
                                p->fields_unified[ghost_buffer_size + i + type * p->n_cells_local] +=
                                    p->left_tmp_buffer[i];
                                p->fields_unified[p->n_cells_local - 2 * ghost_buffer_size + i +
                                                  type * p->n_cells_local] += p->right_tmp_buffer[i];
                            }
                        //Update the buffers of the neighbors
                        MPI_Isend(fields_unified + (p->n_cells_local - 2 * ghost_buffer_size) +
                                  type * p->n_cells_local, ghost_buffer_size, MPI_UINT16_T, right_rank,
                                  2, p->info_MPI.SOMA_comm_sim, req + 0);
                        MPI_Irecv(fields_unified + type * p->n_cells_local, ghost_buffer_size,
                                  MPI_UINT16_T, left_rank, 2, p->info_MPI.SOMA_comm_sim, req + 1);

                        MPI_Isend(fields_unified + ghost_buffer_size + type * p->n_cells_local,
                                  ghost_buffer_size, MPI_UINT16_T, left_rank, 3, p->info_MPI.SOMA_comm_sim, req + 2);
                        MPI_Irecv(fields_unified + (p->n_cells_local - ghost_buffer_size) +
                                  type * p->n_cells_local, ghost_buffer_size, MPI_UINT16_T, right_rank,
                                  3, p->info_MPI.SOMA_comm_sim, req + 3);
                        MPI_Waitall(4, req, stat);
                    }
            }
        //Update all domain ranks with the results of the root domain rank
        MPI_Bcast(fields_unified, p->n_cells_local * p->n_types, MPI_UINT16_T, 0, p->info_MPI.SOMA_comm_domain);

#ifndef ENABLE_MPI_CUDA
#pragma acc update device(p->fields_unified[0:p->n_cells_local*p->n_types])
#else                           //ENABLE_MPI_CUDA
    }                           //Closing OpenACC device pointer region for MPI_CUDA
#endif                          //ENABLE_MPI_CUDA
}

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
                communicate_simple(p);
            else
                communicate_domain_decomposition(p);
        }
    //Avoid false load balancing
    MPI_Barrier(p->info_MPI.SOMA_comm_sim);
#endif                          //ENABLE_MPI
}

int copy_density_32_to_16(const struct Phase *const p)
{
    const uint64_t n_indices = p->n_types * p->n_cells_local;
#pragma acc parallel loop independent present(p[0:1])
    for (uint64_t index = 0; index < n_indices; index++)
        p->fields_unified[index] = p->fields_32[index];
    return 0;
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
            Monomer *beads = p->ph.beads.ptr;
            beads += p->polymers[i].bead_offset;
#pragma acc loop vector
            for (unsigned int j = 0; j < N; j++)
                {               /*Loop over monomers */
                    const unsigned int monotype = get_particle_type(p, i, j);

                    const uint64_t cell = coord_to_index(p, beads[j].x, beads[j].y, beads[j].z);
                    const uint64_t index = cell_to_index_unified(p, cell, monotype);
                    if (cell < p->n_cells_local)        //assuming monotype is correct. Otherwise insert (&& monotype < p->n_types)
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

    copy_density_32_to_16(p);
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

    if (p->efieldsolver != -1) {
        update_invblav(p); // update invblav (inverse of average Bjerrum length)
        update_d_invblav(p); // update dinvblav (derivative of inverse of average Bjerrum length respect to number of segments)
	update_rhoF(p);  // update polymer charge density
        update_exp_born(p); // update born energy, always do this after updating invblav		
        update_electric_field(p);
    }

    switch (p->hamiltonian)
        {
        case SCMF0:
            update_omega_fields_scmf0(p);
            break;
        case SCMF1:
            update_omega_fields_scmf1(p);
            break;
        default:
            fprintf(stderr, "ERROR: %s:%d Unkown hamiltonian specified %d.\n", __FILE__, __LINE__, p->hamiltonian);
        }
}

//! Set the self interaction terms to the omega fields. This
//! includes the intercation with the external, umbrella and electric fields
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
                    p->omega_field_unified[cell + T_types * p->n_cells_local] =
                        inverse_refbeads * (p->xn[T_types * p->n_types + T_types] * (p->tempfield[cell] - 1.0));
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
                                -inverse_refbeads * p->k_umbrella[T_types] *
                                (p->umbrella_field[cell + T_types * p->n_cells_local] -
                                 p->field_scaling_type[T_types] * p->fields_unified[cell + T_types * p->n_cells_local]); 

                        }
		    // electric field
		    
                          if (p->efieldsolver != -1) {
                            p->omega_field_unified[cell + T_types * p->n_cells_local] +=
                                p->electric_field[cell] * p->charges[T_types];
                          }
                } // cells
        }  // types

// Dielectric contribution
  if (p->efieldsolver != -1) {

     unsigned int  ix, iy, iz;
     unsigned int  ixp, ixm, iyp, iym, izp, izm, cell;

     soma_scalar_t psi[p->nx][p->ny][p->nz]; 

     soma_scalar_t gradpsi2[p->n_cells_local]; 

     soma_scalar_t  deltax = p->Lx/((soma_scalar_t) p->nx);
     soma_scalar_t  deltay = p->Ly/((soma_scalar_t) p->ny);
     soma_scalar_t  deltaz = p->Lz/((soma_scalar_t) p->nz);
 
     soma_scalar_t temp ; // auxiliary variables
     soma_scalar_t constq = 4.0*M_PI/(deltax*deltay*deltaz); // multiplicative constant for Poisson equation


#pragma omp parallel for  
for (ix = 0 ; ix < p->nx ; ix++) {
      for (iy = 0 ; iy < p->ny ; iy++) {
           for (iz = 0 ; iz < p->nz ; iz++) {
              cell = cell_coordinate_to_index(p, ix, iy, iz);
	      psi[ix][iy][iz] = p->electric_field[cell];
	   }
       }
}  

#pragma omp parallel for  
for (ix = 0 ; ix < p->nx ; ix++) {
      ixp = mod((ix+1),p->nx);
      ixm = mod((ix-1),p->nx);
      for (iy = 0 ; iy < p->ny ; iy++) {
         iyp = mod((iy+1),p->ny);
         iym = mod((iy-1),p->ny);
            for (iz = 0 ; iz < p->nz ; iz++) {
               izp = mod((iz+1),p->nz);
               izm = mod((iz-1),p->nz);

               cell = cell_coordinate_to_index(p, ix, iy, iz);

	       temp = (psi[ixp][iy][iz]-psi[ixm][iy][iz])/2./deltax; 
               gradpsi2[cell] = temp*temp;

	       temp = (psi[ix][iyp][iz]-psi[ix][iym][iz])/2./deltay; 
               gradpsi2[cell] += temp*temp;

	       temp = (psi[ix][iy][izp]-psi[ix][iy][izm])/2./deltaz; 
               gradpsi2[cell] += temp*temp;

              } // iz 
       } // iy
} // ix

for (unsigned int type = 0; type < p->n_types; type++) {    /*Loop over all fields according to monotype */
    for (uint64_t cell = 0; cell < p->n_cells_local; cell++) {   /*Loop over all cells, max number of cells is product of nx, ny,nz */
 
//    printf("cell, type, gradpsi2, dl/dN todo  %d %d %f %f %f \n", cell,type,gradpsi2[cell],p->d_invblav[cell + type*p->n_cells_local], gradpsi2[cell]*p->d_invblav[cell + type*p->n_cells_local]);
	    
    p->omega_field_unified[cell + type*p->n_cells_local] += -0.5/constq*gradpsi2[cell]*p->d_invblav[cell + type*p->n_cells_local];

    } // cell 	    
} // type	

} // if efieldsolver

} // end routine


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



void calc_ions(struct Phase *const p)
{

  int ix, iy, iz,cell;
  unsigned int type ;
  soma_scalar_t sumrhoQ = 0.0;                   // total number of A segments 

  for (ix = 0 ; ix < (int) p->nx ; ix++) {
          for (iy = 0 ; iy < (int) p->ny ; iy++) {
                          for (iz = 0 ; iz < (int)  p->nz ; iz++) {
                          cell = cell_coordinate_to_index(p, ix, iy, iz);
			  for (type = 0; type < p->n_types ; type++) {
                             sumrhoQ += p->fields_unified[cell + p->n_cells_local*type]*p->charges[type]; 
			  }
                }
          }
  }	  

  p->Nposions = p->Nnegions = p->Nions;

  if (sumrhoQ > 0.0) { 
      p->Nnegions += sumrhoQ; 
      }
  else {
      p->Nposions += -sumrhoQ;
  }

  soma_scalar_t Nionsdiff = p->Nposions-p->Nnegions;
  soma_scalar_t  netcharge = Nionsdiff + sumrhoQ;

  if (p->info_MPI.sim_rank == 0) {
    fprintf(stdout, "calc_ions: Total charge of fixed beads: %f \n ", sumrhoQ);
    fprintf(stdout, "calc_ions: Total number of added salt ions: %f \n ", p->Nions);
    fprintf(stdout, "calc_ions: Total number of +1 salt ions: %f \n ", p->Nposions);
    fprintf(stdout, "calc_ions: Total number of -1 salt ions: %f \n ", p->Nnegions);
    fprintf(stdout, "check_electro: Net charge:  %f \n ", netcharge);
     fflush(stdout);
  }
 
  assert(fabs(netcharge) < 1.0e-6);
}




void update_electric_field(const struct Phase *const p)
{
    // Update electric potential
    //
    if ((p->efieldsolver == 0)||(p->efieldsolver == 2)) 
    	call_PB(p);
    else if (p->efieldsolver == 1) 
    	call_EN(p);
}  

int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}

void calc_invbls(struct Phase *const p) 
{
unsigned int type;

  p->invbls = (soma_scalar_t *) malloc((p->n_types) * sizeof(soma_scalar_t));
  p->invblav_zero = 0.0; 

  for (type = 0; type < p->n_types ; type++) {
	p->invbls[type] = 1./p->bls[type];
	p->invblav_zero += 1./p->bls[type];
  }
        p->invblav_zero = p->invblav_zero/((soma_scalar_t) p->n_types);
}


void update_invblav(const struct Phase *const p) 

{
int tmpsegsum;
unsigned int cell, type;


    for (cell = 0 ; cell < p->n_cells_local ; cell++) {

                          p->invblav[cell] = 0.0 ; 
			  tmpsegsum = 0;

		          for (type = 0 ; type < p->n_types; type++) {

                             p->invblav[cell] += ((soma_scalar_t) p->fields_unified[cell+p->n_cells_local*type])*p->invbls[type];
                             tmpsegsum +=  p->fields_unified[cell+p->n_cells_local*type];
			  } 
                          if(tmpsegsum != 0)    
                            p->invblav[cell] = p->invblav[cell] / ((soma_scalar_t) tmpsegsum);
			  else if (tmpsegsum == 0) 
                            p->invblav[cell] = p->invblav_zero; // prevents divergence if tmpsegsum = 0
	  }
}


void update_d_invblav(const struct Phase *const p) 

{
unsigned int cell, type;
unsigned int Ntot[p->n_cells_local]; 

#pragma omp parallel for  
 for (cell = 0 ; cell < p->n_cells_local ; cell++) {
   Ntot[cell] = 0;
   for (type = 0; type < p->n_types; type++) {    /*Loop over all fields according to monotype */
         Ntot[cell] += p->fields_unified[cell+p->n_cells_local*type]; // total number of segments in cell  
     } // types
   if(Ntot[cell] == 0) 
	   Ntot[cell] = 1;  // prevents divergence of de/dN when N = 0
 } // cell


for (unsigned int type = 0; type < p->n_types; type++) {    /*Loop over all fields according to monotype */

#pragma omp parallel for
 for (cell = 0 ; cell < p->n_cells_local ; cell++) {

  p->d_invblav[cell+p->n_cells_local*type] = (p->invbls[type] - p->invblav[cell])/((soma_scalar_t) Ntot[cell]);
  }
}

}

void update_exp_born(const struct Phase *const p) 

{
unsigned int cell;

#pragma omp parallel for    
for (i = 0 ; i < p->n_cells_local ; i++) {
   p->exp_born[i] = exp(-1.0/(p->invblav[i]*2.0*p->Born_a)); 
}
}


void update_rhoF(const struct Phase *const p) 

{
unsigned int cell, type;
const soma_scalar_t  deltax = p->Lx/((soma_scalar_t) p->nx);
const soma_scalar_t  deltay = p->Ly/((soma_scalar_t) p->ny);
const soma_scalar_t  deltaz = p->Lz/((soma_scalar_t) p->nz);
const soma_scalar_t  vcell = deltax*deltay*deltaz;

#pragma omp parallel for    
for (i = 0 ; i < p->n_cells_local ; i++) {
   p->rhoF[i] = 0.0; 
   for (type = 0 ; type < p->n_types; type++) {
             p->rhoF[i] += p->fields_unified[i+p->n_cells_local*type]*p->charges[type];
   } 
   p->rhoF[i] = p->rhoF[i] / vcell ; // units of charge per Re^3
   exp_born[i] = exp(-1.0/(p->invblav[i]*2.0*p->Born_a)); 
}
}



