/* Copyright (C) 2016-2017 Ludwig Schneider

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

//! \file generate_positions.c
//! \brief Implementation of generate_positions.h

#include "generate_positions.h"
#include <stdbool.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "polymer.h"
#include "monomer.h"
#include "bond.h"
#include "phase.h"
#include "mc.h"
#include "mesh.h"

//! Helper to get the next not set index in a molecule
//!
//! \private
//! \param already_set Array indicating unset particels
//! \param N Number of particles in molecule
//! \return next index.
int get_next_index(const bool*const already_set,const unsigned int N)
    {
    unsigned int i;
    for(i=0; i < N; i++)
	if(! already_set[i] )
	    break;
    if( i==N )
	return -1;
    return i;
    }

//! Helper to set the neighbour of a particle.
//! \private
//!
//! \param jbead Particle index to set.
//! \param neigh Bonded Neighbour
//! \param bond_type Type of bond to neighbour
//! \param already_set Array of already set particles
//! \param poly Polymer of the particle
//! \param p System
//! \return Errorcode
int set_neighbour(const unsigned int jbead,const Monomer*const neigh,
		  const unsigned int bond_type,bool*const already_set,
		  Polymer*const poly,const struct Phase*const p)
    {
    if(already_set[jbead])
	return 0;
    Monomer dx;
    Monomer new;
    int move_allowed;

    do{
	soma_scalar_t scale = 1.;
	dx.x=dx.y=dx.z=0;
    switch(bond_type)
	{
	case HARMONICVARIABLESCALE:;
	    scale = p->harmonic_normb_variable_scale;
	case HARMONIC: ;
	    soma_normal_vector(&(poly->poly_state),p->args.pseudo_random_number_generator_arg
			       , &(dx.x), &(dx.y), &(dx.z));
	    dx.x /= sqrt(2*p->harmonic_normb*scale);
	    dx.y /= sqrt(2*p->harmonic_normb*scale);
	    dx.z /= sqrt(2*p->harmonic_normb*scale);
	    new.x = neigh->x+dx.x; new.y = neigh->y+dx.y; new.z = neigh->z+dx.z;
	    break;
case STIFF:
	default:
	fprintf(stderr, "ERROR: %s:%d unknow bond type appeared %d\n",
		__FILE__, __LINE__,bond_type);
	new.x=new.y=new.z=0; //Shut up compiler warning
	}
    move_allowed = ! possible_move_area51(p,neigh->x,neigh->y,neigh->z,dx.x,dx.y,dx.z,true);
	}while( move_allowed );

    poly->beads[jbead].x = new.x; poly->beads[jbead].y = new.y;     poly->beads[jbead].z = new.z;
    already_set[jbead] = true;

    //recursively add all connected neighbors
    const int start = get_bondlist_offset(
	p->poly_arch[p->poly_type_offset[poly->type] + jbead + 1]);

    if(start > 0){
	int i = start;
	unsigned int end;
	do{
	    const uint32_t info = p->poly_arch[i++];
	    end = get_end(info);
	    const unsigned int bond_type = get_bond_type(info);
	    const int offset = get_offset(info);

	    const int neighbour_id = jbead + offset;
	    set_neighbour(neighbour_id,&(poly->beads[jbead]),bond_type,already_set,poly,p);

	    }while( end == 0);
	}
    return 0;
    }

int generate_new_beads(struct Phase*const p)
    {
    if( fabs(p->harmonic_normb_variable_scale) < 1e-5 )
	{
	fprintf(stderr,"WARNING: p->harmonic_normb_variable_scale < 1e-5, this may result in unreasonable generated position or even causes divisions by 0.\n");
	}
    for( uint64_t i= 0; i < p->n_polymers; i++)
	{
	Polymer *const poly = &(p->polymers[i]);
	const unsigned int N = p->poly_arch[ p->poly_type_offset[poly->type] ];
	bool *const already_set = (bool*)malloc( N *sizeof(bool));
	if(already_set == NULL)
	    {
	    fprintf(stderr,"ERROR: %s:%d Malloc problem.\n",__FILE__,__LINE__);
	    return -1;
	    }
	memset(already_set,false,N*sizeof(bool));

	int free_index;
	while( (free_index=get_next_index(already_set, N)) >= 0 )
	    {
	    //Set a free monomer
	    soma_scalar_t x,y,z;
	    do{
		x = soma_rng_soma_scalar(&(poly->poly_state),p->args.pseudo_random_number_generator_arg)*p->Lx;
		y = soma_rng_soma_scalar(&(poly->poly_state),p->args.pseudo_random_number_generator_arg)*p->Ly;
		z = soma_rng_soma_scalar(&(poly->poly_state),p->args.pseudo_random_number_generator_arg)*p->Lz;
		}while( p->area51 != NULL &&
			p->area51[coord_to_index(p, x, y, z)] == 1);

	    poly->beads[free_index].x = x; poly->beads[free_index].y = y;
	    poly->beads[free_index].z = z; already_set[free_index] = true;
	    //Set recursively all connected neighbors.
	    const int start = get_bondlist_offset(
		p->poly_arch[p->poly_type_offset[poly->type] + free_index + 1]);
	    if(start > 0){
		int i = start;
		unsigned int end;
		do{
		    const int info = p->poly_arch[i++];
		    end = get_end(info);
		    const unsigned int bond_type = get_bond_type(info);
		    const int offset = get_offset(info);

		    const int neighbour_id = free_index + offset;
		    const unsigned int jbead = neighbour_id;

		    set_neighbour(jbead,&(poly->beads[free_index]),bond_type,already_set,poly,p);

		    }while( end == 0);
		}
	    }

	free(already_set);
//transfer the particle positions after generation to GPU
#pragma acc update device(poly->beads[0:N])
	//Init MSD positions
	memcpy( poly->msd_beads, poly->beads, N*sizeof(Monomer));
	}

    update_density_fields(p);
    memcpy(p->old_fields_unified, p->fields_unified, p->n_cells*p->n_types*sizeof(uint16_t));
    return 0;
    }
