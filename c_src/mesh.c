/* Copyright (C) 2016-2017 Ludwig Schneider
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
#include "mpiroutines.h"

uint64_t cell_to_index_unified(const struct Phase*const p,const uint64_t cell,const unsigned int rtype)
    {
//Unified data layout [type][x][y][z]
    return cell + rtype * p->n_cells;
    }


uint64_t coord_to_index(const struct Phase * p, const soma_scalar_t rx,
			    const soma_scalar_t ry, const soma_scalar_t rz)
    {

    int x, y, z;

    coord_to_cell_coordinate(p, rx, ry, rz, &x, &y, &z);

    return cell_coordinate_to_index(p, x, y, z);
    }

uint64_t coord_to_index_unified(const struct Phase * p, const soma_scalar_t rx,
				    const soma_scalar_t ry, const soma_scalar_t rz, const unsigned int rtype)
    {
    int x, y, z;

    coord_to_cell_coordinate(p, rx, ry, rz, &x, &y, &z);
    const uint64_t cell =cell_coordinate_to_index(p, x, y, z);
    return cell_to_index_unified(p,cell,rtype);
    }

/* #pragma acc routine seq */
/* inline void increment_16_bit_uint(uint16_t*const ptr) */
/* { */
/*   int x = (((size_t)ptr) & 2)>>1 ; */
/*   uint32_t * ptr32 = (uint32_t *) (ptr-x) ; */
/* #pragma acc atomic update */
/*   *ptr32 += (1 << (x*16)) ; */
/* } */

void update_density_fields(const struct Phase *const p)
    {
    static unsigned int last_time_call = 0;
    if (last_time_call == 0 || p->time > last_time_call)
	last_time_call = p->time;
    else			//Quick exit, because the property has already been calculated for the time step.
	return;


    const uint64_t n_indices = p->n_types*p->n_cells;
#pragma acc parallel
    // num_gangs(200) vector_length(128)
	{
#pragma acc loop independent
#pragma omp parallel for
	for (uint64_t index = 0; index < n_indices; index++)	/*Loop over all fields according to monotype */
	    p->fields_32[index] = 0;
	}

	const uint64_t n_polymers = p->n_polymers;

#pragma acc parallel num_gangs(n_polymers) vector_length(128) present(p->poly_arch[0:1])

	    {
#pragma acc loop gang
#pragma omp parallel for
	    for (uint64_t i = 0; i < n_polymers; i++){	/*Loop over polymers */
		const unsigned int N = p->poly_arch[ p->poly_type_offset[p->polymers[i].type]];
#pragma acc loop vector
		for (unsigned int j = 0; j < N; j++) {	/*Loop over monomers */
		    const unsigned int monotype = get_particle_type(
			p->poly_arch[ p->poly_type_offset[p->polymers[i].type]+1+j]);

		    const unsigned int index = coord_to_index_unified(p, p->polymers[i].beads[j].x,
								      p->polymers[i].beads[j].y,
								      p->polymers[i].beads[j].z, monotype);
#pragma acc atomic update
#pragma omp atomic
		    p->fields_32[index] += 1;
		    }
		}
	    }

#pragma acc parallel
	    // num_gangs(200) vector_length(128)
		{
#pragma acc loop independent
		for(uint64_t index = 0; index < n_indices; index++)
		    p->fields_unified[index] = p->fields_32[index];
		}

		/*Share the densityfields -> needed because Hamiltonian may not only quadratic order */
		//Const cast!
		mpi_divergence((struct Phase*const) p);

		if (p->info_MPI.Ncores > 1){
#ifndef ENABLE_MPI_CUDA
#pragma acc update self(p->fields_unified[0:n_indices])
		    MPI_Allreduce(MPI_IN_PLACE, p->fields_unified, n_indices, MPI_UINT16_T, MPI_SUM, p->info_MPI.SOMA_MPI_Comm);
#pragma acc update device(p->fields_unified[0:n_indices])
#else//ENABLE_MPI_CUDA
		    uint16_t * fields_unified = p->fields_unified ;
#pragma acc host_data use_device(fields_unified)
			{
			MPI_Allreduce(MPI_IN_PLACE, fields_unified, n_indices, MPI_UINT16_T, MPI_SUM, p->info_MPI.SOMA_MPI_Comm);
			}
#pragma acc update self(p->fields_unified[0:n_indices])

#endif//ENABLE_MPI_CUDA
		    }
		/* Calculate the added up densities */

		/*Use first type to initialize the fields-> saves set zero routine*/
		soma_scalar_t rescale_density= p->field_scaling_type[0];
#pragma acc parallel loop
#pragma omp parallel for
		for (uint64_t cell = 0; cell < p->n_cells; cell++)
		    p->tempfield[cell] = rescale_density*p->fields_unified[cell];


		for (unsigned int T_types = 1; T_types < p->n_types; T_types++){
		    rescale_density = p->field_scaling_type[T_types];
#pragma acc parallel loop
#pragma omp parallel for
		    for (uint64_t cell = 0; cell < p->n_cells; cell++)
			p->tempfield[cell] += rescale_density*p->fields_unified[T_types * p->n_cells + cell];
		    /*!\todo p->ncells as a temporary variable */
		    }
//#pragma acc update self(p->tempfield[0:p->n_cells])
    }

void update_omega_fields(const struct Phase *const p)
    {
    static unsigned int last_time_call = 0;
    if (last_time_call == 0 || p->time > last_time_call)
	last_time_call = p->time;
    else			//Quick exit, because the property has already been calculated for the time step.
	return;
    switch(p->hamiltonian)
	{
	case SCMF0:
	    update_omega_fields_scmf0(p); break;
	case SCMF1:
	    update_omega_fields_scmf1(p); break;
	default:
	    fprintf(stderr,"ERROR: %s:%d Unkown hamiltonian specified %d.\n",__FILE__,__LINE__,p->hamiltonian);
	}

    }

void update_omega_fields_scmf0(const struct Phase *const p)
    {
    /*Densityfields are shorts and unscaled according to bead_type*/
    /*Tempfields save the complete densities and remain as type soma_scalar_t -> used for insothermal Compressibility*/

    /*Calculate new omegafields on each processor itself */
    const soma_scalar_t inverse_refbeads = 1.0 / p->reference_Nbeads;

    // Compressibility part + external fields
    for (unsigned int T_types = 0; T_types < p->n_types; T_types++)	/*Loop over all fields according to monotype */
#pragma acc parallel loop present(p[:1])
#pragma omp parallel for
	for (uint64_t cell = 0; cell < p->n_cells; cell++){	/*Loop over all cells, max number of cells is product of nx, ny,nz */
	    p->omega_field_unified[cell + T_types*p->n_cells] =
		inverse_refbeads * (p->xn[T_types][T_types] * (p->tempfield[cell] - 1.0));
	    /* the external field is defined such that the energy of a chain of refbeads in this field is x k_B T,
	       thus the normalization per bead */
            if(p->external_field_unified != NULL){
		p->omega_field_unified[cell + T_types*p->n_cells] +=
		    inverse_refbeads * p->external_field_unified[ cell+T_types*p->n_cells];
		}
	    }
    // XN part

    //soma_scalar_t weight;/*can be added if nessassary */

    for (unsigned int T_types = 0; T_types < p->n_types; T_types++) {	/*Loop over all fields according to monotype */
	for (unsigned int S_types = T_types + 1; S_types < p->n_types;
	     S_types++){
            // precalculate the normalization for this type combination
            soma_scalar_t dnorm = -0.5 * inverse_refbeads * p->xn[T_types][S_types];
#pragma acc parallel loop present(p[:1])
#pragma omp parallel for
	    for (uint64_t cell = 0; cell < p->n_cells; cell++) {	/*Loop over all cells, max number of cells is product of nx, ny,nz */
		soma_scalar_t interaction =
		    dnorm *
		    (p->field_scaling_type[T_types]*p->fields_unified[cell+T_types*p->n_cells] - p->field_scaling_type[S_types]*p->fields_unified[cell+S_types*p->n_cells]); /*Added the rescaling cause p->fields are short now*/
		p->omega_field_unified[cell+T_types*p->n_cells] += interaction;
		p->omega_field_unified[cell+S_types*p->n_cells] -= interaction;
		}
	    }
	}
    }

void update_omega_fields_scmf1(const struct Phase *const p)
    {
    const soma_scalar_t inverse_refbeads = 1.0 / p->reference_Nbeads;
    const uint64_t n_cells= p->n_cells;
    const unsigned int n_types = p->n_types;

#pragma acc parallel loop collapse(2) independent present(p[0:1])
#pragma omp parallel for
    for (uint64_t cell = 0; cell < n_cells; cell++)
	{
	for (unsigned int T_types = 0; T_types < n_types; T_types++)
	    {
	    //Overwrite old fields with incompressibility part
	    p->omega_field_unified[cell + T_types*p->n_cells] = inverse_refbeads *(p->xn[T_types][T_types] * (p->tempfield[cell] -1 ));

	    //Add external fields.
	    if( p->external_field_unified != NULL)
		p->omega_field_unified[cell +T_types * p->n_cells] += inverse_refbeads * p->external_field_unified[ cell+T_types*p->n_cells];

	    //Add pairwise interaction for all other types
#pragma acc loop seq
	    for (unsigned int S_types = 0; S_types < n_types;S_types++)
		if( T_types != S_types)
		    {
		    const soma_scalar_t norm =  inverse_refbeads * p->xn[T_types][S_types];
		    const soma_scalar_t interaction =  norm * p->fields_unified[cell+ S_types*p->n_cells] * p->field_scaling_type[S_types];
		    p->omega_field_unified[cell+T_types*p->n_cells] += interaction;
		    }
	    }
	}

    }
