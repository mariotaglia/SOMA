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

//! \file polymer.c
//! \brief Implementation of polymer.h

#include "polymer.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#ifdef _OPENACC
#include <openacc.h>
#endif                          //_OPENACC
#include "phase.h"

int reallocate_polymer_mem(struct Phase *const p, uint64_t new_storage)
{
    if (p->present_on_device)
        {
            fprintf(stderr,
                    "ERROR: %s:%d %d reallocate of poly mem, but system is present on device. Call copyout_phase first.\n",
                    __FILE__, __LINE__, p->info_MPI.world_rank);
            return -1;
        }

    const uint64_t heuristics = p->n_polymers_storage * 1.05 + 1;
    if (heuristics > new_storage)
        new_storage = heuristics;

    printf("INFO: @t=%d rank %d is reallocating space for polymers %ld %ld.\n",
           p->time, p->info_MPI.world_rank, new_storage, p->n_polymers_storage);

    struct Polymer *const tmp_poly = (struct Polymer * const)malloc(new_storage * sizeof(struct Polymer));
    if (tmp_poly == NULL)
        {
            fprintf(stderr, "ERROR: %s:%d reallocate malloc %ld\n", __FILE__, __LINE__, new_storage);
            return -1;
        }

    memcpy(tmp_poly, p->polymers, p->n_polymers_storage * sizeof(Polymer));

    free(p->polymers);

    p->n_polymers_storage = new_storage;
    p->polymers = tmp_poly;
    return 0;
}

int push_polymer(struct Phase *const p, const Polymer * const poly)
{
    if (p->present_on_device)
        {
            fprintf(stderr, "ERROR: %s:%d %d PUSH, but system is present on device. Call copyout_phase first.\n",
                    __FILE__, __LINE__, p->info_MPI.world_rank);
            return -1;
        }

    assert(poly);
    if (p->n_polymers >= p->n_polymers_storage)
        reallocate_polymer_mem(p, 0);
    assert(p->n_polymers < p->n_polymers_storage);

    p->polymers[p->n_polymers] = *poly;

    //Update struct
    p->n_polymers += 1;

    const unsigned int N = p->poly_arch[p->poly_type_offset[poly->type]];
    p->num_all_beads_local += N;
    return 0;
}

int pop_polymer(struct Phase *const p, const uint64_t poly_id, Polymer * const poly)
{
    if (p->present_on_device)
        {
            fprintf(stderr, "ERROR: %s:%d POP %d, but system is present on device. Call copyout_phase first.\n",
                    __FILE__, __LINE__, p->info_MPI.world_rank);
            return -1;
        }

    if (poly_id >= p->n_polymers)
        {
            fprintf(stderr, "WARNING: Invalid pop attempt of polymer. rank: %d poly_id %ld n_polymers %ld.\n",
                    p->info_MPI.world_rank, poly_id, p->n_polymers);
            return -1;
        }

    // Copy out the polymer host
    memcpy(poly, p->polymers + poly_id, sizeof(Polymer));

    p->n_polymers -= 1;

    //Fill the gap in vector
    if (poly_id != p->n_polymers)
        {
            memcpy(p->polymers + poly_id, p->polymers + p->n_polymers, sizeof(Polymer));
        }
    const unsigned int N = p->poly_arch[p->poly_type_offset[poly->type]];
    p->num_all_beads_local -= N;

    return 0;
}

int exchange_polymer(struct Phase *const p, const uint64_t poly_i, const uint64_t poly_j)
{
    if (p->present_on_device)
        {
            fprintf(stderr, "ERROR: %s:%d Exchange %d, but system is present on device. Call copyout_phase first.\n",
                    __FILE__, __LINE__, p->info_MPI.world_rank);
            return -1;
        }

    if (poly_i != poly_j)
        {
            if (poly_i >= p->n_polymers)
                {
                    fprintf(stderr, "WARNING: Invalid pop attempt of polymer. rank: %d poly_id %ld n_polymers %ld.\n",
                            p->info_MPI.sim_rank, poly_i, p->n_polymers);
                    return -1;
                }
            if (poly_j >= p->n_polymers)
                {
                    fprintf(stderr, "WARNING: Invalid pop attempt of polymer. rank: %d poly_id %ld n_polymers %ld.\n",
                            p->info_MPI.sim_rank, poly_j, p->n_polymers);
                    return -1;
                }

            Polymer tmp_poly = p->polymers[poly_j];
            p->polymers[poly_j] = p->polymers[poly_i];
            p->polymers[poly_i] = tmp_poly;
        }
    return 0;
}

unsigned int poly_serial_length(const struct Phase *const p, const Polymer * const poly)
{
    const unsigned int N = p->poly_arch[p->poly_type_offset[poly->type]];

    unsigned int length = 0;
    //Buffer length
    length += sizeof(unsigned int);

    //Type data
    length += sizeof(unsigned int);

    //Rcm data
    length += sizeof(Monomer);

    //Tag information
    length += sizeof(uint64_t);

    //Beads data
    length += N * sizeof(Monomer);

    //msd data
    length += N * sizeof(Monomer);

    //poly RNG state
    length += rng_state_serial_length(p);

    //Monomer types (if necessary)
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    length += N * sizeof(uint8_t);
#endif                          //ENABLE_MONOTYPE_CONVERSIONS

    if (poly->set_permutation_offset != UINT64_MAX)
        length += p->max_n_sets * sizeof(unsigned int);

    if (poly->set_states_offset != UINT64_MAX)
        length += p->max_set_members * rng_state_serial_length(p);

    return length;
}

int serialize_polymer(struct Phase *const p, const Polymer * const poly, unsigned char *const buffer)
{
    const unsigned int N = p->poly_arch[p->poly_type_offset[poly->type]];
    unsigned int position = 0;

    //Buffer length
    const unsigned int length = poly_serial_length(p, poly);
    assert(length != 0);
    memcpy(buffer + position, &length, sizeof(unsigned int));
    position += sizeof(unsigned int);

    //Type data
    memcpy(buffer + position, &(poly->type), sizeof(unsigned int));
    position += sizeof(unsigned int);

    //Rcm data
    memcpy(buffer + position, &(poly->rcm), sizeof(Monomer));
    position += sizeof(Monomer);

    //tag data
    memcpy(buffer + position, &(poly->tag), sizeof(uint64_t));
    position += sizeof(uint64_t);

    //Beads data
    memcpy(buffer + position, ((Monomer *) p->ph.beads.ptr) + poly->bead_offset, N * sizeof(Monomer));
    position += N * sizeof(Monomer);

    //MSD data
    memcpy(buffer + position, ((Monomer *) p->ph.msd_beads.ptr) + poly->msd_bead_offset, N * sizeof(Monomer));
    position += N * sizeof(Monomer);

    // Poly state
    position += serialize_rng_state(p, &(poly->poly_state), buffer + position);

    //Monomer types data
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    memcpy(buffer + position, ((uint8_t *) p->ph.monomer_types.ptr) + poly->monomer_type_offset, N * sizeof(uint8_t));
    position += N * sizeof(uint8_t);
#endif                          //ENABLE_MONOTYPE_CONVERSIONS

    // Set permutation
    if (poly->set_permutation_offset != UINT64_MAX)
        {
            memcpy(buffer + position, ((unsigned int *)p->ph.set_permutation.ptr) + poly->set_permutation_offset,
                   p->max_n_sets * sizeof(unsigned int));
            position += p->max_n_sets * sizeof(unsigned int);
        }

    if (poly->set_states_offset != UINT64_MAX)
        for (unsigned int i = 0; i < p->max_set_members; i++)
            position +=
                serialize_rng_state(p, ((RNG_STATE *) p->ph.set_states.ptr) + poly->set_states_offset + i,
                                    buffer + position);

    assert(position == length);
    return position;
}

int deserialize_polymer(struct Phase *const p, Polymer * const poly, const unsigned char *const buffer)
{
    unsigned int position = 0;

    //Buffer length
    unsigned int length_tmp;
    memcpy(&length_tmp, buffer + position, sizeof(unsigned int));
    position += sizeof(unsigned int);
    const unsigned int length = length_tmp;
    if (length < sizeof(unsigned int))
        {
            fprintf(stderr, "ERROR: %s:%d:%d invalid buffer received %d \n",
                    __FILE__, __LINE__, p->info_MPI.world_rank, length);
            return length;
        }

    //Type data
    memcpy(&(poly->type), buffer + position, sizeof(unsigned int));
    position += sizeof(unsigned int);
    const unsigned int N = p->poly_arch[p->poly_type_offset[poly->type]];

    //RCM
    memcpy(&(poly->rcm), buffer + position, sizeof(Monomer));
    position += sizeof(Monomer);

    //tag data
    memcpy(&(poly->tag), buffer + position, sizeof(uint64_t));
    position += sizeof(uint64_t);

    //Beads data
    poly->bead_offset = get_new_soma_memory_offset(&(p->ph.beads), N);
    if (poly->bead_offset == UINT64_MAX)
        {
            fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__, __LINE__,
                    p->info_MPI.world_rank, p->n_polymers);
            return -1;
        }
    memcpy(((Monomer *) p->ph.beads.ptr) + poly->bead_offset, buffer + position, N * sizeof(Monomer));
    position += N * sizeof(Monomer);

    //MSD data
    poly->msd_bead_offset = get_new_soma_memory_offset(&(p->ph.msd_beads), N);
    if (poly->msd_bead_offset == UINT64_MAX)
        {
            fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__, __LINE__,
                    p->info_MPI.world_rank, p->n_polymers);
            return -1;
        }
    memcpy(((Monomer *) p->ph.msd_beads.ptr) + poly->msd_bead_offset, buffer + position, N * sizeof(Monomer));
    position += N * sizeof(Monomer);

    // Poly state
    position += deserialize_rng_state(p, &(poly->poly_state), buffer + position);

    poly->set_permutation_offset = UINT64_MAX;
    poly->set_states_offset = UINT64_MAX;

#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    //Monomer types data
    poly->monomer_type_offset = get_new_soma_memory_offset(&(p->ph.monomer_types), N);
    if (poly->monomer_type_offset == UINT64_MAX)
        {
            fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__, __LINE__,
                    p->info_MPI.world_rank, p->n_polymers);
            return -1;
        }
    memcpy(((uint8_t *) p->ph.monomer_types.ptr) + poly->monomer_type_offset, buffer + position, N * sizeof(uint8_t));
    position += N * sizeof(uint8_t);
#endif                          //ENABLE_MONOTYPE_CONVERSIONS

    // If there is more data in the buffer, this polymer carries set information.
    if (length > position)
        {
            poly->set_permutation_offset = get_new_soma_memory_offset(&(p->ph.set_permutation), p->max_n_sets);
            if (poly->set_permutation_offset == UINT64_MAX)
                {
                    fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__, __LINE__,
                            p->info_MPI.world_rank, p->n_polymers);
                    return -1;
                }
            memcpy(((unsigned int *)p->ph.set_permutation.ptr) + poly->set_permutation_offset, buffer + position,
                   p->max_n_sets * sizeof(unsigned int));
            position += p->max_n_sets * sizeof(unsigned int);
        }

    if (length > position)
        {
            poly->set_states_offset = get_new_soma_memory_offset(&(p->ph.set_states), p->max_set_members);
            if (poly->set_states_offset == UINT64_MAX)
                {
                    fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__, __LINE__,
                            p->info_MPI.world_rank, p->n_polymers);
                    return -1;
                }
            for (unsigned int i = 0; i < p->max_set_members; i++)
                position +=
                    deserialize_rng_state(p, ((RNG_STATE *) p->ph.set_states.ptr) + poly->set_states_offset + i,
                                          buffer + position);
        }
    else
        {
            assert(poly->set_permutation_offset == UINT64_MAX);
        }
    if (position != length)
        {
            fprintf(stderr, "ERROR: %s:%d rank %d Deserialization of polymer. "
                    " The read buffer size %d, does not coincide with length %d "
                    " claimed by the buffer content.\n", __FILE__, __LINE__, p->info_MPI.world_rank, position, length);
            return -2;
        }

    return position;
}

int update_polymer_rcm(struct Phase *const p)
{
    static unsigned int last_time_call = 0;
    if (last_time_call == 0 || p->time > last_time_call)
        last_time_call = p->time;
    else
        return 0;

    const unsigned int n_polymers = p->n_polymers;

#pragma acc parallel loop present(p[0:1])
    for (uint64_t npoly = 0; npoly < n_polymers; npoly++)
        {
            Polymer *const mypoly = &(p->polymers[npoly]);
            const unsigned int N = p->poly_arch[p->poly_type_offset[mypoly->type]];
            Monomer rcm = make_monomer(0, 0, 0);
            Monomer *beads = p->ph.beads.ptr;
            beads += mypoly->bead_offset;

            for (unsigned int i = 0; i < N; i++)
                {
                    rcm.x += beads[i].x;
                    rcm.y += beads[i].y;
                    rcm.z += beads[i].z;
                }
            rcm.x /= N;
            rcm.y /= N;
            rcm.z /= N;
            mypoly->rcm = rcm;
        }
    return 0;
}
