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

//! \file test.c
//! \brief Implementation of test.h

#include "test.h"
#include <stdlib.h>
#include <stdbool.h>
#include "rng.h"
#include "io.h"
#include "init.h"
#include "mesh.h"
#include "mc.h"
#include "independent_sets.h"


int test_read_write_hdf5(const struct Phase * const p)
    {
    int status;
    unsigned int syserror;
    if ((status = write_config_hdf5(p, "/tmp/p1.h5")) != 0) {
	fprintf(stderr, "ERROR: test_read_write_hdf5 %s:%d\n", __FILE__,
		__LINE__);
	return status;
	}

    MPI_Barrier(p->info_MPI.SOMA_MPI_Comm);

    struct Phase phase2;
    struct Phase *const p2 = &phase2;
    p2->info_MPI = p->info_MPI;
    if ((status = read_config_hdf5(p2, "/tmp/p1.h5")) != 0) {
	fprintf(stderr, "ERROR: test_read_write_hdf5 %s:%d\n", __FILE__,
		__LINE__);
	return status;
	}
    init_phase(p2);

    if ((status = write_config_hdf5(p2, "/tmp/p2.h5")) != 0) {
	fprintf(stderr, "ERROR: test_read_write_hdf5 %s:%d\n", __FILE__,
		__LINE__);
	return status;
	}
    free_phase(p2);

    //Maybe the system() approach is not the best. But I think for testing purposes this is fine.
    if ((status = system("h5diff /tmp/p1.h5 /tmp/p2.h5")) != 0) {
	fprintf(stderr, "ERROR: test_read_write_hdf5 %s:%d\n", __FILE__,
		__LINE__);
	return status;
	}
    //Clean up files
    syserror = system("rm -f /tmp/p1.h5");
    if (syserror != 0)
	fprintf(stderr, "ERROR: removing p1.h5 failed\n");
    syserror = system("rm -f /tmp/p2.h5");
    if (syserror != 0)
	fprintf(stderr, "ERROR: removing p2.h5 failed\n");

    if (p->info_MPI.current_core == 0)
	printf("INFO: At t= %d read_write_hdf5 test passed\n", p->time);
    return 0;
    }

int test_particle_types(const struct Phase * const p)
    {
    for (uint64_t i = 0; i < p->n_poly_type; i++)
	{
	const unsigned int type_offset = p->poly_type_offset[i];
	const unsigned int N = p->poly_arch[ type_offset ];
	for(unsigned int j=0; j < N; j++)
	    {
	    const uint32_t info_bl = p->poly_arch[type_offset + 1 +j];
	    const unsigned int type = get_particle_type(info_bl);
	    if (type >= p->n_types)
		{
		fprintf(stderr,
			"ERROR: min. 1 particle has an undefined type. "
			"polytype= %u bead= %u type= %u n_types= %u info_bl=%d \n"
			, (unsigned int)i, j,type, p->n_types,info_bl);
		return type;
		}
	    }
	}
    if (p->info_MPI.current_core == 0)
	printf("INFO: At t= %d particle_type test test passed\n", p->time);
    return 0;
    }

int test_area51_violation(const struct Phase * const p)
    {
    if(p->area51 != NULL){

        for (uint64_t i = 0; i < p->n_polymers; i++)
	    {
            const unsigned int N = p->poly_arch[ p->poly_type_offset[ p->polymers[i].type ] ];
            for (unsigned int j = 0; j < N; j++)
		{
                const uint64_t index = coord_to_index(p,p->polymers[i].beads[j].x,
                                                          p->polymers[i].beads[j].y,
                                                          p->polymers[i].beads[j].z);
                if(p->area51[index] == 1){
                    fprintf(stderr,"ERROR: particle %u %u is in a forbidden area.\n",(unsigned int)i,j);
                    return index;
		    }

		}
	    }
	}
    if (p->info_MPI.current_core == 0)
	printf("INFO: At t= %d area51 violation test passed\n", p->time);
    return 0;
    }

int test_independet_sets(const struct Phase*const p)
    {
    if(p->sets == NULL)
	return 0;

    int ret=0;
    unsigned int divergence = 0;
    for(unsigned int poly_type=0; poly_type < p->n_poly_type; poly_type++)
	{
	struct IndependetSets*const set = &(p->sets[poly_type]);
	unsigned int largest_set=0;
	unsigned int smallest_set = UINT_MAX;
	//printf("PolyType %d:\n",poly_type);
	for(unsigned int set_id =0; set_id < set->n_sets; set_id++)
	    {
	    if(set->set_length[set_id] > largest_set)
		largest_set = set->set_length[set_id];
	    if(set->set_length[set_id] < smallest_set)
		smallest_set = set->set_length[set_id];
	    for(unsigned int i=0; i < set->set_length[set_id]; i++)
		{
		const unsigned int pi = set->sets[ set_id*set->max_member + i];	
		unsigned int n_neigh = 0;
		for(unsigned int j=0; j < set->set_length[set_id]; j++)
		    {
		     const int start = get_bondlist_offset(
		     const unsigned int pj =set->sets[ set_id*set->max_member +j]; 
		     p->poly_arch[p->poly_type_offset[poly_type] + pj + 1]);
		     if(start > 0)
			{
			unsigned int i=start;
			unsigned int end;
			do{
			    const uint32_t info = p->poly_arch[i++];
			    end = get_end(info);
			    const int offset = get_offset(info);
			    const unsigned n_id = pj+offset;
			    if(n_id == pi)
				n_neigh++;
			    }while(end==0);
			}
		    }
		if(p->info_MPI.current_core == 0 && n_neigh >0)
		    fprintf(stderr,"ERROR: mono %d of poly_type %d has neighbors in its set %d.\n",
			    pi,poly_type,set_id);
		ret += n_neigh;
		}
	    //printf("\n");
	    }
	if( largest_set - smallest_set > divergence)
	    divergence = largest_set - smallest_set;
	}

    if( p->info_MPI.current_core == 0 && ret == 0)
	printf("INFO: Test independet sets passed with max divergence %d.\n",divergence);
    return ret;
    }

int test_area51_exact(const struct Phase * const p)
    {
    unsigned int violations = 0;
    if(p->area51 != NULL){
        for (uint64_t i = 0; i < p->n_polymers; i++)
	    {
            const unsigned int N = p->poly_arch[ p->poly_type_offset[ p->polymers[i].type ] ];
            for (unsigned int j = 0; j < N -1 ; j++)
		{
		const Monomer a = p->polymers[i].beads[j];
		const Monomer b = p->polymers[i].beads[j+1];
		Monomer dx;
		dx.x = b.x-a.x;
		dx.y = b.y-a.y;
		dx.z = b.z-a.z;
		const bool ok = possible_move_area51(p,a.x,a.y,a.z,dx.x,dx.y,dx.z,true);
		if( ! ok)
		    violations += 1;
		}
	    }
	}
    MPI_Allreduce(MPI_IN_PLACE,&violations,1,MPI_UNSIGNED,MPI_SUM,p->info_MPI.SOMA_MPI_Comm);
    if (p->info_MPI.current_core == 0)
	{
	if( violations == 0)
	    printf("INFO: At t= %d area51 exact test passed\n", p->time);
	else
	    printf("WARNING: At t= %d area51 exact test **FAILED** with %d violations.\n", p->time,violations);
	}
    return violations;
    }
