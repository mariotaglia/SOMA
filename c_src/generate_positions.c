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
#include "poly_heavy.h"

//! Helper to get the next not set index in a molecule
//!
//! \private
//! \param already_set Array indicating unset particels
//! \param N Number of particles in molecule
//! \return next index.
int get_next_index(const bool *const already_set, const unsigned int N)
{
    unsigned int i;
    for (i = 0; i < N; i++)
        if (!already_set[i])
            break;
    if (i == N)
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
int set_neighbour(const unsigned int jbead, const Monomer * const neigh,
                  const unsigned int bond_type, bool *const already_set,
                  Polymer * const poly, const struct Phase *const p)
{
    Monomer dx;
    Monomer new;
    int move_allowed;

    Monomer *beads = p->ph.beads.ptr;
    beads += poly->bead_offset;
    do
        {
            soma_scalar_t scale = 1.;
            dx.x = dx.y = dx.z = 0;
            switch (bond_type)
                {
                case HARMONICVARIABLESCALE:;
                    scale = p->harmonic_normb_variable_scale;
                    /* intentionally falls through */
                case HARMONIC:;
                    soma_normal_vector(&(poly->poly_state), p, &(dx.x), &(dx.y), &(dx.z));
                    dx.x /= sqrt(2 * p->harmonic_normb * scale);
                    dx.y /= sqrt(2 * p->harmonic_normb * scale);
                    dx.z /= sqrt(2 * p->harmonic_normb * scale);
                    new.x = neigh->x + dx.x;
                    new.y = neigh->y + dx.y;
                    new.z = neigh->z + dx.z;
                    break;
                case STIFF:
                default:
                    fprintf(stderr, "ERROR: %s:%d unknow bond type appeared %d\n", __FILE__, __LINE__, bond_type);
                    new.x = new.y = new.z = 0;  //Shut up compiler warning
                }
            move_allowed = !possible_move_area51(p, neigh->x, neigh->y, neigh->z, dx.x, dx.y, dx.z, false);
    } while (move_allowed);

    beads[jbead].x = new.x;
    beads[jbead].y = new.y;
    beads[jbead].z = new.z;
    already_set[jbead] = true;
    return 0;
}

int generate_new_beads(struct Phase *const p)
{
    if (fabs(p->harmonic_normb_variable_scale) < 1e-5)
        {
            fprintf(stderr,
                    "WARNING: p->harmonic_normb_variable_scale < 1e-5, this may result in unreasonable generated position or even causes divisions by 0.\n");
        }
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    const soma_scalar_t domain_offset = my_domain * (p->Lx / p->args.N_domains_arg);
    for (uint64_t i = 0; i < p->n_polymers; i++)
        {
            Polymer *const poly = &(p->polymers[i]);
            const unsigned int N = p->poly_arch[p->poly_type_offset[poly->type]];
            Monomer *beads = p->ph.beads.ptr;
            beads += poly->bead_offset;

            unsigned int *chain = (unsigned int *)malloc(N * sizeof(unsigned int));
            if (chain == NULL)
                {
                    fprintf(stderr, "ERROR: %s:%d Malloc problem.\n", __FILE__, __LINE__);
                    return -1;
                }
            memset(chain, 0, N * sizeof(unsigned int));

            int chain_index = 0;
            bool *const already_set = (bool *)malloc(N * sizeof(bool));
            if (already_set == NULL)
                {
                    fprintf(stderr, "ERROR: %s:%d Malloc problem.\n", __FILE__, __LINE__);
                    return -1;
                }
            memset(already_set, false, N * sizeof(bool));

            int free_index;
            free_index = get_next_index(already_set, N);
            while (free_index >= 0)
                {
                    //Set a free monomer
                    soma_scalar_t x, y, z;
                    unsigned int domain_counter = 0;
                    uint64_t idx;
                    do
                        {
                            soma_scalar_t r = soma_rng_soma_scalar(&(poly->poly_state), p);
                            // 128 Tries to place the free bead in the local domain.
                            // If that fails (for example, because the local domain is only area51),
                            // place the bead anywhere. It will be sent later to the correct domain.
                            // This is slower, but still correct + not infinite.
                            if (domain_counter < 128)
                                x = domain_offset + r * (p->Lx / p->args.N_domains_arg);
                            else
                                x = r * p->Lx;
                            r = soma_rng_soma_scalar(&(poly->poly_state), p);
                            y = r * p->Ly;
                            r = soma_rng_soma_scalar(&(poly->poly_state), p);
                            z = r * p->Lz;

                            domain_counter += 1;
                            idx = coord_to_index(p, x, y, z);   // > p->n_cells_local if position is out of the domain
                    } while (p->area51 != NULL && (idx >= p->n_cells_local || p->area51[idx] == 1));

                    beads[free_index].x = x;
                    beads[free_index].y = y;
                    beads[free_index].z = z;
                    already_set[free_index] = true;
                    chain[chain_index] = free_index;
                    unsigned int total_set = 1;

                    unsigned int old_bead = free_index;
                    unsigned int new_bead;
                    unsigned int bond_type;
                    const int old_bead_bondlist_offset =
                        get_bondlist_offset(p->poly_arch[p->poly_type_offset[poly->type] + free_index + 1]);
                    if (old_bead_bondlist_offset != -1)
                        {
                            int bondlist = old_bead_bondlist_offset;
                            unsigned int end = 0;
                            do
                                {
                                    do
                                        {
                                            if (end != 0)
                                                {       //if all the bonds are already set, go back
                                                    chain_index--;
                                                    if (chain_index < 0)
                                                        break;
                                                    bondlist =
                                                        get_bondlist_offset(p->poly_arch
                                                                            [p->poly_type_offset[poly->type] +
                                                                             chain[chain_index] + 1]);
                                                    old_bead = chain[chain_index];
                                                }
                                            const int info = p->poly_arch[bondlist];
                                            end = get_end(info);
                                            bond_type = get_bond_type(info);
                                            const int offset = get_offset(info);
                                            new_bead = old_bead + offset;
                                            bondlist++;
                                    } while (already_set[new_bead]);    //find the first unset neigbour
                                    chain_index++;
                                    set_neighbour(new_bead, &(beads[old_bead]), bond_type, already_set, poly, p);
                                    already_set[new_bead] = true;
                                    total_set++;
                                    bondlist =
                                        get_bondlist_offset(p->poly_arch
                                                            [p->poly_type_offset[poly->type] + new_bead + 1]);
                                    chain[chain_index] = new_bead;
                                    old_bead = new_bead;
                                    end = 0;
                            } while (chain_index > 0);
                        }
                    free_index = get_next_index(already_set, N);
                }
            free(already_set);
            free(chain);
        }                       //loop over polymers
    //We assume that the msd_beads and the beads use the same memory offsets
    memcpy(p->ph.msd_beads.ptr, p->ph.beads.ptr, p->ph.beads.used * p->ph.beads.typelength);
    //transfer the particle positions after generation to GPU
    update_device_polymer_heavy(p, false);
    update_density_fields(p);
    memcpy(p->old_fields_unified, p->fields_unified, p->n_cells_local * p->n_types * sizeof(uint16_t));
    return 0;
}

int generate_monomer_type_array(struct Phase *const p)
{
    (void)p;
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    if (p->info_MPI.sim_rank == 0)
        printf("INFO: Generating monomer type array from poly_arch info, since no data was found in h5-file.\n");
    for (uint64_t i = 0; i < p->n_polymers; i++)
        {
            Polymer *const poly = p->polymers + i;
            const unsigned int N = p->poly_arch[p->poly_type_offset[poly->type]];
            uint8_t *monomer_type = p->ph.monomer_types.ptr;
            monomer_type += poly->monomer_type_offset;

            for (unsigned int j = 0; j < N; j++)
                {
                    monomer_type[j] =
                        get_particle_type_of_poly_arch(p->poly_arch[p->poly_type_offset[poly->type] + 1 + j]);
                }
        }
    if (p->bead_data_read == true)
        {                       //update device and density fields, if bead data is already present. 
            update_device_polymer_heavy(p, false);
            update_density_fields(p);
            memcpy(p->old_fields_unified, p->fields_unified, p->n_cells_local * p->n_types * sizeof(uint16_t));
        }
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    return 0;
}
