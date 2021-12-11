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
#ifndef INDEPENDENT_SETS_H
#define INDEPENDENT_SETS_H
#include "soma_config.h"
struct Phase;

//! \file independent_sets.h
//! Definition of code related to the preparation of independent sets

//! Struct to store the independet set information for a polymer type.
typedef struct IndependetSets {
    unsigned int n_sets;        //!< Number of sets
    unsigned int max_member;    //!< Max number of members per set.
    unsigned int *set_length;   //!< Length of each set
    //! array storing the sets. Flatten array: Access \code
    //! sets[iSet*max_member + iElement]; \endcode Length: n_set*max_member.
    unsigned int *sets;
} IndependetSets;

//! \file independent_sets.h
//! /Brief Functions needed for independent set preparations

//! According to the cmdline argument, this function calls either the "simple algorithm" or the "fixed n sets algorithm".
//! "simple algorithm" is well suited for cases like linear chain while
//! "fixed n set algothrim" is suited for complex crosslinked polymer system.
//! \param p Phase for which the polymers are assigned to sets
//! \return Errorcode
int generate_independet_sets(struct Phase *const p);

//! Generate the independet set information for each poly_type.
//! \param p Phase for which the polymers are assigned to sets
//! \return Errorcode
int independent_sets_simple(struct Phase *const p);

//! Generate independent sets for each poly_type.
//! This algorithm uses n+1 sets to store the particles, with n the number of bonds of the particle with the most bonds.
//! It is much faster then the "simple algorithm" for very long chains.
//! \param p Phase for which the polymers are assigned to sets
//! \return Errorcode
int independent_set_fixed(struct Phase *const p);

//! Private funtion, to be used in combination with independent_set_fixed() to allocate and initiate memory.
//! \param p Phase for which the polymers are assigned to sets
//! \return Errorcode
int allo_init_memory_for_Polystates(struct Phase *const p);

//! Private funtion, to be used in combination with independent_set_fixed().
//! It finds the writein_set that can store the new monomer.
//! \param bonds_total Array to store the bonds of all monomers
//! \param bond_number_total Array to store the number of bonds for all monomers
//! \param max_bond_number Maximal number of bonds of one monomer
//! \param writein_set The set where new monomers are stored
//! \param current_set The set, the member of which is being studied
//! \param offset_set Array which tells how many members of a set is checked
//! \param end_set Array with the number of member in a set
//! \param independent_sets Two-D-Array that stores the independent sets
//! \param bond_i The current bond of the current monomer
//! \param current_monomer The monomer that is being studied
//! \return Set to put new monomer
unsigned int check_bond_members_of_set(unsigned int **bonds_total, int *bond_number_total, unsigned int max_bond_number,
                                       unsigned int writein_set, unsigned int current_set, unsigned int *offset_set,
                                       unsigned int *end_set, unsigned int **independent_sets, int bond_i,
                                       unsigned int current_monomer);

//! Private funtion, to be used in combination with independent_set_fixed().
//! This function calculates independent sets for a single chain.
//! \param set_tmp_pointer The pointer to the array that stores all set information
//! \param n_poly The current poly_type
//! \param p Phase for which the polymers are assigned to sets
//! \return Errorcode
int independent_sets_one_polymer(struct IndependetSets **const set_tmp_pointer, unsigned int n_poly,
                                 struct Phase *const p);
#endif                          //INDEPENDENT_SETS_H
