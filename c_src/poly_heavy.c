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

/*! \file poly_heavy.c
\brief implementation of poly_heavy.h
*/

#include <string.h>
#include "poly_heavy.h"
#include "soma_memory.h"
#include "phase.h"

int free_polymer_heavy(struct Phase *const p)
{
    int status = 0;
    status |= free_soma_memory(&(p->ph.beads));
    status |= free_soma_memory(&(p->ph.msd_beads));
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    status |= free_soma_memory(&(p->ph.monomer_types));
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    if (p->args.iteration_alg_arg == iteration_alg_arg_SET)
        {
            status |= free_soma_memory(&(p->ph.set_states));
            status |= free_soma_memory(&(p->ph.set_permutation));
        }
    return status;
}

int copyin_polymer_heavy(struct Phase *const p)
{
    int status = 0;
    status |= copyin_soma_memory(&(p->ph.beads));
    status |= copyin_soma_memory(&(p->ph.msd_beads));
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    status |= copyin_soma_memory(&(p->ph.monomer_types));
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    if (p->args.iteration_alg_arg == iteration_alg_arg_SET)
        {
            status |= copyin_soma_memory(&(p->ph.set_states));
            status |= copyin_soma_memory(&(p->ph.set_permutation));
        }
    return status;
}

int copyout_polymer_heavy(struct Phase *const p)
{
    int status = 0;
    status |= copyout_soma_memory(&(p->ph.beads));
    status |= copyout_soma_memory(&(p->ph.msd_beads));
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    status |= copyout_soma_memory(&(p->ph.monomer_types));
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    if (p->args.iteration_alg_arg == iteration_alg_arg_SET)
        {
            status |= copyout_soma_memory(&(p->ph.set_states));
            status |= copyout_soma_memory(&(p->ph.set_permutation));
        }
    return status;
}

int update_device_polymer_heavy(struct Phase *const p, const bool rng_flag)
{
    int status = 0;
    status |= update_device_soma_memory(&(p->ph.beads));
    status |= update_device_soma_memory(&(p->ph.msd_beads));
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    status |= update_device_soma_memory(&(p->ph.monomer_types));
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    if (p->args.iteration_alg_arg == iteration_alg_arg_SET && rng_flag)
        {
            status |= update_device_soma_memory(&(p->ph.set_states));
            status |= update_device_soma_memory(&(p->ph.set_permutation));
        }
    return status;
}

int update_self_polymer_heavy(struct Phase *const p, const bool rng_flag)
{
    int status = 0;
    status |= update_self_soma_memory(&(p->ph.beads));
    status |= update_self_soma_memory(&(p->ph.msd_beads));
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    status |= update_self_soma_memory(&(p->ph.monomer_types));
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    if (p->args.iteration_alg_arg == iteration_alg_arg_SET && rng_flag)
        {
            status |= update_self_soma_memory(&(p->ph.set_states));
            status |= update_self_soma_memory(&(p->ph.set_permutation));
        }
    return status;
}

int consider_compact_polymer_heavy(struct Phase *p, const bool collective)
{
    int compact = p->num_all_beads_local * 1.05 < p->ph.beads.length;
    if (collective)
        MPI_Allreduce(MPI_IN_PLACE, &compact, 1, MPI_INT, MPI_LOR, p->info_MPI.SOMA_comm_sim);

    if (!compact)
        return 0;

    const bool beads_device = p->ph.beads.device_present;
    if (beads_device)
        copyout_soma_memory(&(p->ph.beads));
    const bool msd_beads_device = p->ph.msd_beads.device_present;
    if (msd_beads_device)
        copyout_soma_memory(&(p->ph.msd_beads));
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    const bool monomer_types_device = p->ph.monomer_types.device_present;
    if (monomer_types_device)
        copyout_soma_memory(&(p->ph.monomer_types));
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    const bool states_device = p->ph.set_states.device_present;
    if (states_device && (p->args.iteration_alg_arg == iteration_alg_arg_SET))
        copyout_soma_memory(&(p->ph.set_states));
    const bool permutation_device = p->ph.set_permutation.device_present;
    if (permutation_device && (p->args.iteration_alg_arg == iteration_alg_arg_SET))
        copyout_soma_memory(&(p->ph.set_permutation));

    PolymerHeavy new_ph;
    if (init_soma_memory(&(new_ph.beads), p->num_all_beads_local, sizeof(Monomer)) != 0)
        {
            fprintf(stderr, "ERROR compacting beads\n");
            return -1;
        }
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    if (init_soma_memory(&(new_ph.monomer_types), p->num_all_beads_local, sizeof(uint8_t)) != 0)
        {
            fprintf(stderr, "ERROR compacting beads\n");
            return -1;
        }
#endif                          //ENABLE_MONOTYPE_CONVERSIONS

    if (init_soma_memory(&(new_ph.msd_beads), p->num_all_beads_local, sizeof(Monomer)) != 0)
        {
            fprintf(stderr, "ERROR compacting msd_beads\n");
            return -2;
        }
    if (p->args.iteration_alg_arg == iteration_alg_arg_SET)
        {
            if (init_soma_memory(&(new_ph.set_states), p->n_polymers * p->max_set_members, sizeof(RNG_STATE)) != 0)
                {
                    fprintf(stderr, "ERROR compacting set_states\n");
                    return -3;
                }
            if (init_soma_memory(&(new_ph.set_permutation), p->n_polymers * p->max_n_sets, sizeof(unsigned int)) != 0)
                {
                    fprintf(stderr, "ERROR compacting set_permutation\n");
                    return -4;
                }
        }

    for (uint64_t i = 0; i < p->n_polymers; i++)
        {
            const unsigned int N = p->poly_arch[p->poly_type_offset[p->polymers[i].type]];

            uint64_t new_bead_offset = get_new_soma_memory_offset(&(new_ph.beads), N);
            Monomer *old_beads = p->ph.beads.ptr;
            old_beads += p->polymers[i].bead_offset;
            Monomer *new_beads = new_ph.beads.ptr;
            new_beads += new_bead_offset;

            memcpy(new_beads, old_beads, N * sizeof(Monomer));
            p->polymers[i].bead_offset = new_bead_offset;

#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
            uint64_t new_monomer_type_offset = get_new_soma_memory_offset(&(new_ph.monomer_types), N);
            uint8_t *old_monomer_types = p->ph.monomer_types.ptr;
            old_monomer_types += p->polymers[i].monomer_type_offset;
            uint8_t *new_monomer_types = new_ph.monomer_types.ptr;
            new_monomer_types += new_monomer_type_offset;

            memcpy(new_monomer_types, old_monomer_types, N * sizeof(uint8_t));
            p->polymers[i].monomer_type_offset = new_monomer_type_offset;
#endif                          //ENABLE_MONOTYPE_CONVERSIONS

            uint64_t new_msd_bead_offset = get_new_soma_memory_offset(&(new_ph.msd_beads), N);
            Monomer *old_msd_beads = p->ph.msd_beads.ptr;
            old_msd_beads += p->polymers[i].msd_bead_offset;
            Monomer *new_msd_beads = new_ph.msd_beads.ptr;
            new_msd_beads += new_msd_bead_offset;

            memcpy(new_msd_beads, old_msd_beads, N * sizeof(Monomer));
            p->polymers[i].msd_bead_offset = new_msd_bead_offset;

            if (p->args.iteration_alg_arg == iteration_alg_arg_SET)
                {
                    uint64_t new_state_offset = get_new_soma_memory_offset(&(new_ph.set_states), p->max_set_members);
                    RNG_STATE *old_state = p->ph.set_states.ptr;
                    old_state += p->polymers[i].set_states_offset;
                    RNG_STATE *new_state = new_ph.set_states.ptr;
                    new_state += new_state_offset;

                    memcpy(new_state, old_state, p->max_set_members * sizeof(RNG_STATE));
                    p->polymers[i].set_states_offset = new_state_offset;

                    //No copy for permutation necessary, because it is just temporary memory
                    p->polymers[i].set_permutation_offset =
                        get_new_soma_memory_offset(&(new_ph.set_permutation), p->max_n_sets);
                }

        }

    fprintf(stdout, "World rank %d compacts memory at time step %d. Balance %lu old %lu new %lu\n",
            p->info_MPI.world_rank, p->time, p->num_all_beads_local, p->ph.beads.length, new_ph.beads.length);
    free_polymer_heavy(p);
    p->ph = new_ph;
    if (beads_device)
        copyin_soma_memory(&(p->ph.beads));
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    if (monomer_types_device)
        copyin_soma_memory(&(p->ph.monomer_types));
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    if (msd_beads_device)
        copyin_soma_memory(&(p->ph.msd_beads));
    if (states_device && (p->args.iteration_alg_arg == iteration_alg_arg_SET))
        copyin_soma_memory(&(p->ph.set_states));
    if (permutation_device && (p->args.iteration_alg_arg == iteration_alg_arg_SET))
        copyin_soma_memory(&(p->ph.set_permutation));

    return 0;
}
