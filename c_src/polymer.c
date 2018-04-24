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

//! \file polymer.c
//! \brief Implementation of polymer.h


#include "polymer.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#ifdef _OPENACC
#include <openacc.h>
#endif//_OPENACC
#include "phase.h"

int free_polymer(const struct Phase*const p, Polymer*const poly)
    {
    free(poly->beads);
    free(poly->msd_beads);
    deallocate_rng_state(&(poly->poly_state), p->args.pseudo_random_number_generator_arg);

    if( poly->set_states != NULL)
        {
        for(unsigned int j=0; j < p->max_set_members; j++)
            {
            deallocate_rng_state(poly->set_states+j, p->args.pseudo_random_number_generator_arg);
            free(poly->set_states[j].mt_state);
            free(poly->set_states[j].tt800_state);
            }
        free(poly->set_states);
        free(poly->set_permutation);
        }
    return 0;
    }

int copyin_polymer(struct Phase*const p, Polymer*const poly)
    {
    const unsigned int N = p->poly_arch[ p->poly_type_offset[ poly->type ] ];
#pragma acc enter data copyin(poly->beads[0:N])
    copyin_rng_state( &(poly->poly_state), p->args.pseudo_random_number_generator_arg);

    if(poly->set_permutation !=NULL)
        {
#pragma acc enter data copyin(poly->set_permutation[0:p->max_n_sets])
        }

    if(poly->set_states !=NULL)
        {
#pragma acc enter data copyin(poly->set_states[0:p->max_set_members])
        for(unsigned int j=0; j < p->max_set_members; j++)
            {
            copyin_rng_state(poly->set_states +j, p->args.pseudo_random_number_generator_arg);
            }
        }
    return 0 + 1*N*0;
    }

int copyout_polymer(struct Phase*const p, Polymer*const poly)
    {
    const unsigned int N = p->poly_arch[ p->poly_type_offset[poly->type] ];
#pragma acc exit data copyout(poly->beads[0:N])
    copyout_rng_state(&(poly->poly_state), p->args.pseudo_random_number_generator_arg);

      if(poly->set_permutation !=NULL)
          {
#pragma acc exit data copyout(poly->set_permutation[0:p->max_n_sets])
          }

      if(poly->set_states !=NULL)
          {
          for(unsigned int j=0; j < p->max_set_members; j++)
              {
              copyout_rng_state(poly->set_states+j, p->args.pseudo_random_number_generator_arg);
              }
#pragma acc exit data copyout(poly->set_states[0:p->max_set_members])
          }
      return 0*N;
    }

int reallocate_polymer_mem(struct Phase*const p,uint64_t new_storage)
    {
    const uint64_t heuristics = p->n_polymers_storage * 1.05 + 1;
    if( heuristics > new_storage)
        new_storage = heuristics;

    printf("INFO: @t=%d rank %d is reallocating space for polymers %ld %ld.\n",
           p->time,p->info_MPI.world_rank,new_storage,p->n_polymers_storage);

    struct Polymer*const tmp_poly = (struct Polymer*const)malloc(new_storage*sizeof(struct Polymer));
    if( tmp_poly == NULL)
        {
        fprintf(stderr,"ERROR: %s:%d reallocate malloc %ld\n",__FILE__,__LINE__,new_storage);
        return -1;
        }
    memcpy( tmp_poly, p->polymers, p->n_polymers_storage*sizeof(Polymer));

#ifdef _OPENACC
#pragma acc enter data copyin (tmp_poly[0:new_storage])
    Polymer* tmp_dev = (Polymer*)acc_deviceptr((Polymer*)tmp_poly);
    Polymer*const poly_dev = (Polymer*)acc_deviceptr((Polymer*const)p->polymers);
    acc_memcpy_device( tmp_dev , poly_dev, p->n_polymers_storage*sizeof(Polymer) );
#endif//_OPENACC

#pragma acc exit data delete(p->polymers[0:p->n_polymers_storage])
    free( p->polymers );

    p->n_polymers_storage = new_storage;
    p->polymers = tmp_poly;
#ifdef _OPENACC
    struct Phase*const p_dev = acc_deviceptr(p);

    acc_memcpy_to_device( &(p_dev->polymers), &tmp_dev, sizeof(Polymer*));

#endif//_OPENACC
    return 0;
    }

int push_polymer(struct Phase*const p,const Polymer*const poly)
    {
    assert(poly);
    if(p->n_polymers >= p->n_polymers_storage)
        reallocate_polymer_mem(p,0);
    assert( p->n_polymers < p->n_polymers_storage);

    p->polymers[p->n_polymers] = *poly;

#ifdef _OPENACC
    const unsigned int pos = p->n_polymers;
    Polymer*const polymers_dev = acc_deviceptr( p->polymers);
    acc_memcpy_to_device( polymers_dev + pos, p->polymers + pos, sizeof(Polymer));

#endif//_OPENACC

    copyin_polymer(p,&(p->polymers[p->n_polymers]));

    //Update struct
    p->n_polymers += 1;
#pragma acc update device(p->n_polymers)

    const unsigned int N = p->poly_arch[ p->poly_type_offset[poly->type] ];
    for (unsigned int k = 0; k < N; k++)
        {
        const unsigned int type = get_particle_type(
            p->poly_arch[ p->poly_type_offset[poly->type]+1+k]);
        p->num_bead_type_local[type] += 1;
        p->num_all_beads_local += 1;
        }
#pragma acc update device(p->num_bead_type_local[0:p->n_types])
#pragma acc update device(p->num_all_beads_local[0:1])

    return 0;
    }

int pop_polymer(struct Phase*const p,const uint64_t poly_id,Polymer*const poly)
    {
    if( poly_id >= p->n_polymers)
        {
        fprintf(stderr,"WARNING: Invalid pop attempt of polymer. rank: %d poly_id %ld n_polymers %ld.\n"
                ,p->info_MPI.world_rank,poly_id,p->n_polymers);
        return -1;
        }

    update_self_polymer(p,p->polymers+poly_id,1);
    // Copy out the polymer host
    memcpy( poly, p->polymers+poly_id, sizeof(Polymer) );

    p->n_polymers -= 1;
#pragma acc update device(p->n_polymers)

    //Fill the gap in vector
    memcpy( p->polymers + poly_id, p->polymers + p->n_polymers, sizeof(Polymer));
#ifdef _OPENACC
    Polymer*const polymers_dev = acc_deviceptr(p->polymers);
    acc_memcpy_device( polymers_dev + poly_id, polymers_dev + p->n_polymers, sizeof(Polymer) );
#endif//_OPENACC

    const unsigned int N = p->poly_arch[ p->poly_type_offset[poly->type] ];
    for (unsigned int k = 0; k < N; k++)
        {
        const unsigned int type = get_particle_type(
            p->poly_arch[ p->poly_type_offset[poly->type]+1+k]);
        p->num_bead_type_local[type] -= 1;
        p->num_all_beads_local -= 1;
        }
#pragma acc update device(p->num_bead_type_local[0:p->n_types])
#pragma acc update device(p->num_all_beads_local[0:1])

    copyout_polymer(p, poly);
    return 0;
    }

unsigned int poly_serial_length(const struct Phase*const p,const Polymer*const poly)
    {
    const unsigned int N =p->poly_arch[p->poly_type_offset[poly->type]];

    unsigned int length = 0;
    //Buffer length
    length += sizeof(unsigned int);

    //Type data
    length += sizeof(unsigned int);

    //Rcm data
    length += sizeof( Monomer);

    //Beads data
    length += N*sizeof( Monomer );

    //msd data
    length += N*sizeof( Monomer );

    //poly RNG state
    length += rng_state_serial_length(p);

    if( poly->set_permutation != NULL)
        length += p->max_n_sets * sizeof(unsigned int);

    if( poly->set_states != NULL)
        length += p->max_set_members*rng_state_serial_length(p);

    return length;
    }

int serialize_polymer(const struct Phase*const p,const Polymer*const poly,unsigned char*const buffer)
    {
    const unsigned int N =p->poly_arch[p->poly_type_offset[poly->type]];
    unsigned int position = 0;

    //Buffer length
    const unsigned int length = poly_serial_length(p, poly);
    assert(length != 0);
    memcpy( buffer + position, &length, sizeof(unsigned int));
    position += sizeof(unsigned int);

    //Type data
    memcpy(buffer + position, &(poly->type), sizeof(unsigned int));
    position += sizeof(unsigned int);

    memcpy(buffer + position, &(poly->rcm), sizeof(Monomer) );
    position += sizeof( Monomer );

    //Beads data
    memcpy(buffer + position, poly->beads, N*sizeof(Monomer) );
    position += N*sizeof(Monomer);

    //MSD data
    memcpy(buffer + position , poly->msd_beads, N*sizeof(Monomer) );
    position += N*sizeof(Monomer);

    // Poly state
    position += serialize_rng_state(p, &(poly->poly_state), buffer +position);

    // Set permutation
    if( poly->set_permutation != NULL)
        {
        memcpy(buffer + position, poly->set_permutation, p->max_n_sets * sizeof(unsigned int));
        position += p->max_n_sets * sizeof(unsigned int);
        }

    if( poly->set_states != NULL)
        for(unsigned int i=0; i < p->max_set_members; i++)
            position += serialize_rng_state(p, poly->set_states+i, buffer+position);

    assert( position == length);
    return position;
    }

int deserialize_polymer(const struct Phase*const p, Polymer*const poly,const unsigned char*const buffer)
    {
    unsigned int position = 0;

    //Buffer length
    unsigned int length_tmp;
    memcpy( &length_tmp,buffer + position, sizeof(unsigned int));
    position += sizeof(unsigned int);
    const unsigned int length = length_tmp;
    if(length < sizeof(unsigned int) )
        {
        fprintf(stderr,"ERROR: %s:%d:%d invalid buffer received %d \n",
                __FILE__,__LINE__,p->info_MPI.world_rank,length);
        return length;
        }

    //Type data
    memcpy(&(poly->type),buffer + position, sizeof(unsigned int));
    position += sizeof(unsigned int);
    const unsigned int N =p->poly_arch[p->poly_type_offset[poly->type]];

    //RCM
    memcpy(&(poly->rcm), buffer + position, sizeof(Monomer) );
    position += sizeof(Monomer);

    //Beads data
    poly->beads = (Monomer*)malloc( N*sizeof(Monomer) );
    MALLOC_ERROR_CHECK(poly->beads, N*sizeof(Monomer) );

    memcpy(poly->beads,buffer + position, N*sizeof(Monomer) );
    position += N*sizeof(Monomer);

    //MSD data
    poly->msd_beads = (Monomer*)malloc( N*sizeof(Monomer));
    MALLOC_ERROR_CHECK(poly->msd_beads, N*sizeof(Monomer));

    memcpy(poly->msd_beads, buffer + position, N*sizeof(Monomer) );
    position += N*sizeof(Monomer);

    // Poly state
    position += deserialize_rng_state(p, &(poly->poly_state), buffer +position);

    poly->set_permutation = NULL;
    poly->set_states = NULL;
    // If there is more data in the buffer, this polymer carries set information.
    if( length > position)
        {
        poly->set_permutation = (unsigned int*)malloc( p->max_n_sets * sizeof(unsigned int));
        MALLOC_ERROR_CHECK(poly->set_permutation, p->max_n_sets*sizeof(unsigned int));

        memcpy(poly->set_permutation,buffer + position, p->max_n_sets * sizeof(unsigned int));
        position += p->max_n_sets * sizeof(unsigned int);
        }

    if( length > position)
        {
        poly->set_states = (RNG_STATE*) malloc( p->max_set_members * sizeof(RNG_STATE));
        MALLOC_ERROR_CHECK(poly->set_states, p->max_set_members*sizeof(RNG_STATE));

        for(unsigned int i=0; i < p->max_set_members; i++)
            position += deserialize_rng_state(p, poly->set_states+i, buffer+position);
        }
    else
        {
        assert(poly->set_permutation == NULL);
        }
    if( position != length )
        {
        fprintf(stderr,"ERROR: %s:%d rank %d Deserialization of polymer. "
                " The read buffer size %d, does not coincide with length %d "
                " claimed by the buffer content.\n",
                __FILE__, __LINE__,p->info_MPI.world_rank,position,length);
        return -2;
        }

    return position;
    }

int update_self_polymer(const struct Phase*const p,Polymer*const poly,int rng_update_flag)
    {
    const unsigned int N= p->poly_arch[ p->poly_type_offset[poly->type] ];
    #pragma acc update self(poly->beads[0:N])

    if(poly->set_permutation !=NULL)
        {
#pragma acc update self(poly->set_permutation[0:p->max_n_sets])
        }
    if(rng_update_flag)
        {
          update_self_rng_state(&(poly->poly_state), p->args.pseudo_random_number_generator_arg);
          if(poly->set_states !=NULL)
            {
#pragma acc update self(poly->set_states[0:p->max_set_members])
              for(unsigned int j=0; j < p->max_set_members; j++)
                {
                  update_self_rng_state(poly->set_states +j, p->args.pseudo_random_number_generator_arg);
                }
            }
        }
#pragma acc update self(poly->type)
#pragma acc update self(poly->rcm)
    return 0 + 1*N*0;
    }

int update_polymer_rcm(struct Phase*const p)
    {
    static unsigned int last_time_call = 0;
    if( last_time_call == 0 || p->time > last_time_call)
        last_time_call = p->time;
    else
        return 0;

    const unsigned int n_polymers = p->n_polymers;

    #pragma acc parallel loop
    for(uint64_t npoly=0; npoly < n_polymers; npoly++)
        {
        Polymer *const mypoly = &( p->polymers[npoly] );
        const unsigned int N = p->poly_arch[ p->poly_type_offset[ mypoly->type ]];
        Monomer rcm = make_monomer(0, 0, 0);
        for(unsigned int i=0; i < N; i++)
            {
            rcm.x += mypoly->beads[i].x;
            rcm.y += mypoly->beads[i].y;
            rcm.z += mypoly->beads[i].z;
            }
        rcm.x /= N;
        rcm.y /= N;
        rcm.z /= N;
        mypoly->rcm = rcm;
        }
    return 0;
    }
