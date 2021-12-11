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

//! \file independent_sets.c
//! \brief Implementation of independent_sets.h

#include "independent_sets.h"
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "soma_util.h"
#include "phase.h"
#include "soma_memory.h"
#include "rng.h"

//! Sort the indecies of the set_length in sort_array with the highest at pos=0.
//!
//! \private Used only in the contex of independet set creation
void sort_set_length(const unsigned int n_sets, const unsigned int *const set_length, unsigned int *const sort_array)
{
    //Bubble sort is ok, because it is not time relevant and the data set is small and the set is highly sorted.
    bool elements_swapped;
    do
        {
            elements_swapped = false;
            for (unsigned int i = 1; i < n_sets; i++)
                if (set_length[sort_array[i - 1]] < set_length[sort_array[i]])
                    {
                        elements_swapped = true;
                        unsigned int tmp = sort_array[i];
                        sort_array[i] = sort_array[i - 1];
                        sort_array[i - 1] = tmp;
                    }
    } while (elements_swapped);
}

//! Check, whether a particle can be inserted in a set of independet sets.
//!
//! \private  Used only in the contex of independet set creation
bool try_particle_in_set(const unsigned int set_id, const unsigned int p_id, const struct Phase *const p,
                         const unsigned int N, const unsigned int poly_type, const unsigned int *const set_length,
                         const unsigned int *set_matrix)
{
    bool neighbor_found = false;
    for (unsigned int iTmp = 0; iTmp < set_length[set_id] && neighbor_found == false; iTmp++)
        {
            const unsigned int ibead = set_matrix[set_id * N + iTmp];
            const int start = get_bondlist_offset(p->poly_arch[p->poly_type_offset[poly_type] + ibead + 1]);
            if (start > 0)
                {
                    //BondInfo bn;
                    unsigned int end = 0;
                    for (int i = start; end == 0; i++)
                        {
                            const uint32_t info = p->poly_arch[i];
                            end = get_end(info);
                            const int offset = get_offset(info);
                            const unsigned int neighbour_id = ibead + offset;
                            if (neighbour_id == p_id)
                                neighbor_found = true;
                        }
                }
        }

    return neighbor_found == false;
}

//! Balance the number of set members across the sets, if possible.
//!
//! \private  Used only in the contex of independet set creation
int balance_sets(const struct Phase *const p, const unsigned int N, const unsigned int poly_type,
                 unsigned int *const set_length, unsigned int *const set_matrix, unsigned int recursion_level)
{
    const unsigned int max_recursion = 2;
    if (recursion_level > max_recursion)
        return 0;
    unsigned int n_sets = 0;
    for (unsigned int i = 0; i < N; i++)
        if (set_length[i] > 0)
            n_sets++;
    //Sorted array of set lengths
    unsigned int *sort_array = (unsigned int *)malloc(n_sets * sizeof(unsigned int));
    if (sort_array == NULL)
        {
            fprintf(stderr, "ERROR: malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    for (unsigned int i = 0; i < n_sets; i++)
        sort_array[i] = i;

    bool solvable_unbalanced;
    do
        {
            sort_set_length(n_sets, set_length, sort_array);

            solvable_unbalanced = false;        //Set to true, if a element can be moved in the loop.
            //Try to an element of the largest set to any other set.
            for (unsigned int mono = 0; mono < set_length[sort_array[0]]; mono++)
                {
                    unsigned move_id = set_matrix[sort_array[0] * N + mono];
                    //Try to fit in any of the other sets, start with the smallest.
                    unsigned int new_set = 0;
                    for (new_set = n_sets - 1; new_set > 0; new_set--)
                        //New set has to be smaller than the largest set, even if we remove one element of the largest set.
                        if (set_length[sort_array[new_set]] < set_length[sort_array[0]] - 1)
                            {
                                const bool suitable_set =
                                    try_particle_in_set(sort_array[new_set], move_id, p, N, poly_type, set_length,
                                                        set_matrix);
                                if (suitable_set)
                                    break;
                            }

                    if (new_set > 0)    //Found an element to move and ****BREAK**** for loop, but do not stop while loop.
                        {
                            //reorganise the set_matrix of the largest set.
                            for (unsigned int i = mono; i < set_length[sort_array[0]] - 1; i++)
                                set_matrix[sort_array[0] * N + i] = set_matrix[sort_array[0] * N + i + 1];
                            set_length[sort_array[0]] -= 1;

                            //Insert the element in its new set.
                            set_matrix[sort_array[new_set] * N + set_length[sort_array[new_set]]] = move_id;
                            set_length[sort_array[new_set]] += 1;

                            solvable_unbalanced = true; //
                            break;      //Exit for loop, because we want to move one element at a time.
                        }
                }
    } while (solvable_unbalanced);

    //If we have even after rebalancing with a constant number of sets
    //a bad divergence, we can add a set and try the rebalancing
    //again.
    const unsigned int max_tolerated_divergence = 16;
    sort_set_length(n_sets, set_length, sort_array);
    unsigned int divergence = 0;
    for (unsigned int set = 0; set < n_sets; set++)
        if (set_length[sort_array[0]] - set_length[set] > divergence)
            divergence = set_length[sort_array[0]] - set_length[set];
    if (divergence > max_tolerated_divergence && n_sets < N)
        {
            set_length[sort_array[0]] -= 1;
            set_matrix[n_sets * N + 0] = set_matrix[sort_array[0] * N + set_length[sort_array[0]]];
            set_length[n_sets] += 1;

            balance_sets(p, N, poly_type, set_length, set_matrix, recursion_level + 1);
        }

    free(sort_array);
    return 0;
}

int generate_independet_sets(struct Phase *const p)
{
    int ret = 0;
    switch (p->args.set_generation_algorithm_arg)
        {
        case set_generation_algorithm_arg_FIXEDMINUS_NMINUS_SETS:
            ret = independent_set_fixed(p);
            break;
        case set_generation_algorithm_arg_SIMPLE:
            ret = independent_sets_simple(p);
            break;
        case set_generation_algorithm__NULL:
        default:
            fprintf(stderr, "ERROR: Unknown independent set generation algorithm selected.\n");
            ret = -2;
        }
    return ret;
}

int independent_sets_simple(struct Phase *const p)
{
    int ret = 0;
    struct IndependetSets *const sets = (struct IndependetSets *)malloc(p->n_poly_type * sizeof(IndependetSets));
    if (sets == NULL)
        {
            fprintf(stderr, "ERROR: malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    p->max_set_members = 0;
    unsigned int max_nsets = 0;
    for (unsigned int poly_type = 0; poly_type < p->n_poly_type; poly_type++)
        {
            const unsigned int N = p->poly_arch[p->poly_type_offset[poly_type]];
            unsigned int *const set_length = (unsigned int *)malloc(N * sizeof(unsigned int));
            if (set_length == NULL)
                {
                    fprintf(stderr, "ERROR: malloc %s:%d\n", __FILE__, __LINE__);
                    return -2;
                }
            memset(set_length, 0, N * sizeof(unsigned int));
            unsigned int *const set_matrix = (unsigned int *)malloc(N * N * sizeof(unsigned int));
            if (set_matrix == NULL)
                {
                    fprintf(stderr, "ERROR: malloc %s:%d\n", __FILE__, __LINE__);
                    return -3;
                }

            //Start setting the monomers in the sets.
            for (unsigned int mono = 0; mono < N; mono++)
                {
                    // Find the next possible set to sort in
                    unsigned int set_to_sort_in = 0;
                    bool set_found = false;
                    //Try all sets
                    while (!set_found)
                        {
                            //Check for neighbors in the set
                            const bool suitable_set =
                                try_particle_in_set(set_to_sort_in, mono, p, N, poly_type, set_length, set_matrix);

                            if (suitable_set)
                                set_found = true;
                            set_to_sort_in += 1;
                            assert(set_to_sort_in < N + 1);
                        }
                    set_to_sort_in -= 1;

                    //Sort the monomer in the set
                    set_matrix[set_to_sort_in * N + set_length[set_to_sort_in]] = mono;
                    set_length[set_to_sort_in] += 1;
                }

            balance_sets(p, N, poly_type, set_length, set_matrix, 0);

            //Store the found sets to the final structure and free intermediate arrays.
            unsigned int n_sets = 0;
            unsigned int max_set = 0;
            unsigned int member_sum = 0;
            for (unsigned int i = 0; i < N; i++)
                {
                    if (set_length[i] > 0)
                        n_sets++;
                    if (set_length[i] > max_set)
                        max_set = set_length[i];
                    member_sum += set_length[i];

                }
            assert(member_sum == N);

            sets[poly_type].n_sets = n_sets;
            sets[poly_type].max_member = max_set;
            sets[poly_type].set_length = (unsigned int *)malloc(n_sets * sizeof(unsigned int));
            if (sets[poly_type].set_length == NULL)
                {
                    fprintf(stderr, "ERROR: malloc %s:%d\n", __FILE__, __LINE__);
                    return -4;
                }
            sets[poly_type].sets = (unsigned int *)malloc(n_sets * max_set * sizeof(unsigned int));
            if (sets[poly_type].sets == NULL)
                {
                    fprintf(stderr, "ERROR: malloc %s:%d\n", __FILE__, __LINE__);
                    return -5;
                }
            //Copy data to the new structure
            for (unsigned int iSet = 0; iSet < n_sets; iSet++)
                {
                    sets[poly_type].set_length[iSet] = set_length[iSet];
                    if (set_length[iSet] > p->max_set_members)
                        p->max_set_members = set_length[iSet];
                    for (unsigned int iMember = 0; iMember < set_length[iSet]; iMember++)
                        sets[poly_type].sets[iSet * max_set + iMember] = set_matrix[iSet * N + iMember];
                }
            //Free intermediate arrays
            free(set_length);
            free(set_matrix);
            if (sets[poly_type].n_sets > max_nsets)
                max_nsets = sets[poly_type].n_sets;
        }
    p->sets = sets;
    p->max_n_sets = max_nsets;
    ret = allo_init_memory_for_Polystates(p);
    return ret;
}

int independent_set_fixed(struct Phase *const p)
{
    int ret = 0;
    struct IndependetSets *set_tmp = (struct IndependetSets *)malloc(p->n_poly_type * sizeof(IndependetSets));
    if (set_tmp == NULL)
        {
            fprintf(stderr, "ERROR: malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    p->max_set_members = 0;
    p->max_n_sets = 0;

    for (unsigned int n_poly = 0; n_poly < p->n_poly_type; n_poly++)
        {
            ret = independent_sets_one_polymer(&set_tmp, n_poly, p);
            if (ret != 0)
                {
                    fprintf(stderr, "ERROR: Function independent_sets_one_polymer failed %s:%d\n", __FILE__, __LINE__);
                    return ret;
                }
        }

    p->sets = set_tmp;
    ret = allo_init_memory_for_Polystates(p);
    return ret;
}

int allo_init_memory_for_Polystates(struct Phase *const p)
{
    init_soma_memory(&(p->ph.set_states), p->n_polymers * p->max_set_members, sizeof(RNG_STATE));
    init_soma_memory(&(p->ph.set_permutation), p->n_polymers * p->max_n_sets, sizeof(unsigned int));
    switch (p->args.pseudo_random_number_generator_arg)
        {
        case pseudo_random_number_generator__NULL:
            break;
        case pseudo_random_number_generator_arg_PCG32:
            break;
        case pseudo_random_number_generator_arg_MT:
            reallocate_soma_memory(&(p->rh.mt), p->n_polymers * p->max_n_sets);
            break;
        case pseudo_random_number_generator_arg_TT800:
            reallocate_soma_memory(&(p->rh.tt800), p->n_polymers * p->max_n_sets);
            break;
        }

    for (unsigned int i = 0; i < p->n_polymers; i++)
        {
            Polymer *const poly_tmp = p->polymers + i;

            poly_tmp->set_states_offset = get_new_soma_memory_offset(&(p->ph.set_states), p->max_set_members);
            if (poly_tmp->set_states_offset == UINT64_MAX)
                {
                    fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__, __LINE__,
                            p->info_MPI.world_rank, p->n_polymers);
                    return -1;
                }
            poly_tmp->set_permutation_offset = get_new_soma_memory_offset(&(p->ph.set_permutation), p->max_n_sets);
            if (poly_tmp->set_permutation_offset == UINT64_MAX)
                {
                    fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__, __LINE__,
                            p->info_MPI.world_rank, p->n_polymers);
                    return -1;
                }

            //Init every state in the polymer
            const unsigned int seed = pcg32_random(&(poly_tmp->poly_state.default_state));
            for (unsigned int j = 0; j < p->max_set_members; j++)
                {
                    RNG_STATE *state = p->ph.set_states.ptr;
                    state += poly_tmp->set_states_offset + j;
                    switch (p->args.pseudo_random_number_generator_arg)
                        {
                        case pseudo_random_number_generator__NULL:
                            /* intentionally falls through */
                        case pseudo_random_number_generator_arg_PCG32:
                            state->alternative_rng_offset = UINT64_MAX;
                            break;
                        case pseudo_random_number_generator_arg_MT:
                            state->alternative_rng_offset = get_new_soma_memory_offset(&(p->rh.mt), 1);
                            if (state->alternative_rng_offset == UINT64_MAX)
                                {
                                    fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__,
                                            __LINE__, p->info_MPI.world_rank, p->n_polymers);
                                    return -1;
                                }
                            break;
                        case pseudo_random_number_generator_arg_TT800:
                            state->alternative_rng_offset = get_new_soma_memory_offset(&(p->rh.tt800), 1);
                            if (state->alternative_rng_offset == UINT64_MAX)
                                {
                                    fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__,
                                            __LINE__, p->info_MPI.world_rank, p->n_polymers);
                                    return -1;
                                }
                            break;
                        }
                    seed_rng_state(state, seed, j, p);
                }
        }
    return 0;
}

int independent_sets_one_polymer(struct IndependetSets **const set_tmp_pointer, unsigned int n_poly,
                                 struct Phase *const p)
{
    unsigned int sequence = p->poly_arch[p->poly_type_offset[n_poly]];
    unsigned int max_bond_number = 0;   //the maximal number of bonds of a monomer in the chain
    unsigned int max_bond = 0;  //the monomer which has the maximal number of bonds
    int *bond_number_total;     //total number of bonds of all monomer
    bond_number_total = (int *)malloc(sequence * sizeof(int));
    if (bond_number_total == NULL)
        {
            fprintf(stderr, "Malloc problem %s:%d\n", __FILE__, __LINE__);
            return -4;
        }
    memset(bond_number_total, 0, sequence * sizeof(int));
    unsigned int **bonds_total = (unsigned int **)malloc(sequence * sizeof(unsigned int *));
    if (bonds_total == NULL)
        {
            fprintf(stderr, "Malloc problem %s:%d\n", __FILE__, __LINE__);
            return -4;
        }
    // loop over all monomer to record the bond information
    for (unsigned int monomer_i = p->poly_type_offset[n_poly] + 1;
         monomer_i < sequence + p->poly_type_offset[n_poly] + 1; monomer_i++)
        {
            uint32_t bonds_of_monomer = 0, bond_number = 0;     //the current bond and the number of bonds of the current monomer
            const uint32_t current_poly_arch = p->poly_arch[monomer_i];
            int start_offset_bond = get_bondlist_offset(current_poly_arch);
            //store all bonds information
            if (start_offset_bond < 0)
                {
                    bond_number = 0;
                    bond_number_total[monomer_i - p->poly_type_offset[n_poly] - 1] = bond_number;
                    bonds_total[monomer_i - p->poly_type_offset[n_poly] - 1] = NULL;
                }
            else
                {
                    int end = 0;
                    do
                        {
                            bonds_of_monomer = p->poly_arch[start_offset_bond];
                            end = get_end(bonds_of_monomer);
                            bond_number++;
                            start_offset_bond++;
                    } while (end != 1);
                    bond_number_total[monomer_i - p->poly_type_offset[n_poly] - 1] = bond_number;
                    bonds_total[monomer_i - p->poly_type_offset[n_poly] - 1] =
                        (unsigned int *)malloc(bond_number * sizeof(unsigned int));
                    if (bonds_total[monomer_i - p->poly_type_offset[n_poly] - 1] == NULL)
                        {
                            fprintf(stderr, "ERROR: malloc %s:%d\n", __FILE__, __LINE__);
                            return -1;
                        }
                    memset(&bonds_total[monomer_i - p->poly_type_offset[n_poly] - 1][0], 0,
                           bond_number * sizeof(unsigned int));
                    start_offset_bond = get_bondlist_offset(p->poly_arch[monomer_i]);
                }
            //store the bonds of all monomer to the array bonds_total
            for (unsigned int bond_i = 0; bond_i < bond_number; bond_i++)
                {
                    bonds_of_monomer = p->poly_arch[start_offset_bond];
                    bonds_total[monomer_i - p->poly_type_offset[n_poly] - 1][bond_i] =
                        get_offset(bonds_of_monomer) + monomer_i - p->poly_type_offset[n_poly] - 1;
                    start_offset_bond++;
                }

            if (bond_number > max_bond_number)
                {               //count the number of sets needed
                    max_bond_number = bond_number;
                    max_bond = monomer_i - p->poly_type_offset[n_poly] - 1;
                }
        }
    int *monomer_checked;
    monomer_checked = (int *)malloc(sequence * sizeof(int));
    if (monomer_checked == NULL)
        {
            fprintf(stderr, "Malloc problem %s:%d\n", __FILE__, __LINE__);
            return -4;
        }
    memset(monomer_checked, 0, sequence * sizeof(int));
    unsigned int current_monomer = max_bond;
    unsigned int **independent_sets;
    independent_sets = (unsigned int **)malloc((max_bond_number + 1) * sizeof(unsigned int *));
    if (independent_sets == NULL)
        {
            fprintf(stderr, "Malloc problem %s:%d\n", __FILE__, __LINE__);
            return -4;
        }
    //the monomers stored between offset_set and end_set are the ones whose heighbour are not stored in the sets yet
    unsigned int *end_set = (unsigned int *)malloc((max_bond_number + 1) * sizeof(unsigned int));
    if (end_set == NULL)
        {
            fprintf(stderr, "Malloc problem %s:%d\n", __FILE__, __LINE__);
            return -4;
        }
    unsigned int *offset_set = (unsigned int *)malloc((max_bond_number + 1) * sizeof(unsigned int));
    if (offset_set == NULL)
        {
            fprintf(stderr, "Malloc problem %s:%d\n", __FILE__, __LINE__);
            return -4;
        }
    for (unsigned int set_i = 0; set_i < max_bond_number + 1; set_i++)
        {
            end_set[set_i] = 0;
            offset_set[set_i] = 0;
            independent_sets[set_i] = (unsigned int *)malloc(sequence * sizeof(unsigned int));
            if (independent_sets[set_i] == NULL)
                {
                    fprintf(stderr, "Malloc problem %s:%d\n", __FILE__, __LINE__);
                    return -4;
                }
            memset(&independent_sets[set_i][0], 0, sequence * sizeof(unsigned int));
        }
    unsigned int total_assigned_monomer = 0;
    //loop over all monomer of a poly_type (if the chain consists several molecules)
    for (unsigned int mono_i = 0; mono_i < sequence; mono_i++)
        {
            if (total_assigned_monomer == sequence)
                break;
            if (monomer_checked[mono_i] == -1)
                continue;
            if (mono_i != 0)
                current_monomer = mono_i;
            else
                current_monomer = max_bond;
            monomer_checked[current_monomer] = -1;
            //find the set with least member //only help for if the chain is not separated into single chains
            int start_set = 0;
            for (unsigned int set_i = 0; set_i < max_bond_number + 1; set_i++)
                {
                    if (end_set[set_i] < end_set[start_set])
                        start_set = set_i;
                }

            independent_sets[start_set][end_set[start_set]] = current_monomer;  //write the central monomer into start_set
            end_set[start_set]++;
            total_assigned_monomer++;
            //write the neighbour of the central monomer into remaining sets
            for (int set_i = 1; set_i <= bond_number_total[current_monomer]; set_i++)
                {
                    unsigned int loop_set = set_i + start_set;
                    if (loop_set > max_bond_number)
                        {
                            loop_set = loop_set - max_bond_number - 1;
                            //That should be zero.
                            assert(loop_set == 0);
                        }
                    independent_sets[loop_set][end_set[loop_set]] = bonds_total[current_monomer][set_i - 1];
                    total_assigned_monomer++;
                    monomer_checked[bonds_total[current_monomer][set_i - 1]] = -1;
                    end_set[loop_set]++;
                }

            offset_set[start_set]++;    //neighbour of central monomer are all found now
            unsigned int current_set = start_set + 1;   //the current set we are studying
            if (current_set > max_bond_number)
                current_set = current_set - max_bond_number - 1;
            int chain_finished = 0;
            for (unsigned int set_i = 0; set_i < max_bond_number + 1; set_i++)
                {
                    if (offset_set[set_i] != end_set[set_i])
                        {
                            chain_finished = 1;
                        }
                }
            //start: while not all monomers are checked
            while (chain_finished == 1)
                {
                    chain_finished = 0;
                    if (end_set[current_set] == offset_set[current_set])
                        {
                            current_set++;
                            if (current_set > max_bond_number)
                                current_set = current_set - max_bond_number - 1;
                        }
                    unsigned int writein_set = current_set - 1; //which set to put the new element into (the left one)
                    if (current_set == 0)
                        writein_set = max_bond_number;
                    //loop over all monomers in the current_set, whose neighbours are not studied yet
                    for (unsigned int member_set = offset_set[current_set]; member_set < end_set[current_set];
                         member_set++)
                        {
                            current_monomer = independent_sets[current_set][member_set];
                            for (int bond_i = 0; bond_i < bond_number_total[current_monomer]; bond_i++)
                                {
                                    if (monomer_checked[bonds_total[current_monomer][bond_i]] == -1)
                                        continue;
                                    writein_set =
                                        check_bond_members_of_set(bonds_total, bond_number_total, max_bond_number,
                                                                  writein_set, current_set, offset_set, end_set,
                                                                  independent_sets, bond_i, current_monomer);
                                    independent_sets[writein_set][end_set[writein_set]] =
                                        bonds_total[current_monomer][bond_i];
                                    monomer_checked[bonds_total[current_monomer][bond_i]] = -1;
                                    total_assigned_monomer++;
                                    end_set[writein_set]++;
                                }
                            offset_set[current_set]++;
                        }
                    current_set++;
                    if (current_set > max_bond_number)
                        current_set = current_set - max_bond_number - 1;
                    for (unsigned int set_i = 0; set_i < max_bond_number + 1; set_i++)
                        {
                            if (offset_set[set_i] != end_set[set_i])
                                {
                                    chain_finished = 1;
                                }
                        }
                }
        }                       //end:loop over all monomer of a poly_type
    (*set_tmp_pointer)[n_poly].n_sets = max_bond_number + 1;
    (*set_tmp_pointer)[n_poly].max_member = offset_set[0];
    int start_i = 0;
    for (unsigned int set_i = 1; set_i < max_bond_number + 1; set_i++)
        {
            if (offset_set[set_i] > (*set_tmp_pointer)[n_poly].max_member)
                (*set_tmp_pointer)[n_poly].max_member = offset_set[set_i];
        }
    //restore the sets into inde_set_tmp
    unsigned int *inde_set_tmp =
        (unsigned int *)malloc((*set_tmp_pointer)[n_poly].max_member * (max_bond_number + 1) * sizeof(unsigned int));
    if (inde_set_tmp == NULL)
        {
            fprintf(stderr, "Malloc problem %s:%d\n", __FILE__, __LINE__);
            return -4;
        }
    for (unsigned int set_i = 0; set_i < max_bond_number + 1; set_i++)
        {
            start_i = set_i * (*set_tmp_pointer)[n_poly].max_member;
            for (unsigned int member_i = 0; member_i < end_set[set_i]; member_i++)
                {
                    inde_set_tmp[start_i] = independent_sets[set_i][member_i];
                    start_i++;
                }
        }
    (*set_tmp_pointer)[n_poly].set_length = (unsigned int *)malloc((max_bond_number + 1) * sizeof(unsigned int));
    if ((*set_tmp_pointer)[n_poly].set_length == NULL)
        {
            fprintf(stderr, "Malloc problem %s:%d\n", __FILE__, __LINE__);
            return -4;
        }

    for (unsigned int set_i = 0; set_i < max_bond_number + 1; set_i++)
        {
            (*set_tmp_pointer)[n_poly].set_length[set_i] = end_set[set_i];
        }
    (*set_tmp_pointer)[n_poly].sets = inde_set_tmp;
    if (p->max_n_sets < max_bond_number + 1)
        p->max_n_sets = max_bond_number + 1;

    if (p->max_set_members < (*set_tmp_pointer)[n_poly].max_member)
        p->max_set_members = (*set_tmp_pointer)[n_poly].max_member;

    //free memory
    free(end_set);
    free(offset_set);
    free(bond_number_total);
    for (unsigned int tmp = 0; tmp < sequence; tmp++)
        {
            free(bonds_total[tmp]);
        }
    free(bonds_total);
    free(monomer_checked);

    for (unsigned int tmp = 0; tmp < max_bond_number + 1; tmp++)
        {
            free(independent_sets[tmp]);
        }
    free(independent_sets);
    return 0;
}

unsigned int check_bond_members_of_set(unsigned int **bonds_total, int *bond_number_total, unsigned int max_bond_number,
                                       unsigned int writein_set, unsigned int current_set, unsigned int *offset_set,
                                       unsigned int *end_set, unsigned int **independent_sets, int bond_i,
                                       unsigned int current_monomer)
{
    int found = 0;
    int meet;
    while (found == 0)
        {
            meet = 0;
            for (int neighbour_index = 0; neighbour_index < bond_number_total[bonds_total[current_monomer][bond_i]];
                 neighbour_index++)
                {               //loop over all bonds
                    unsigned int neighbour = bonds_total[bonds_total[current_monomer][bond_i]][neighbour_index];
                    if (neighbour == current_monomer)
                        continue;
                    for (unsigned int member_i = offset_set[writein_set]; member_i < end_set[writein_set]; member_i++)
                        {
                            if (neighbour == independent_sets[writein_set][member_i])
                                {
                                    meet = 1;
                                }
                        }
                }
            if (meet == 1)
                {
                    writein_set++;
                    if (writein_set > max_bond_number)
                        writein_set = writein_set - max_bond_number - 1;
                    if (writein_set == current_set)
                        {
                            writein_set++;
                            if (writein_set > max_bond_number)
                                writein_set = writein_set - max_bond_number - 1;
                        }
                    continue;
                }
            else
                {
                    found = 1;
                    break;
                }
        }
    return writein_set;
}
