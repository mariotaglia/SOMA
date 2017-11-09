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

#ifndef INDEPENDENT_SETS_H
#define INDEPENDENT_SETS_H

struct Phase;

//! \file independent_sets.h
//! Definition of code related to the preparation of independent sets

//! Struct to store the independet set information for a polymer type.
typedef struct IndependetSets{
    unsigned int n_sets; //!< Number of sets
    unsigned int max_member; //!< Max number of members per set.
    unsigned int * set_length;//!< Length of each set
    //! array storing the sets. Flatten array: Access \code
    //! sets[iSet*max_member + iElement]; \endcode Length: n_set*max_member.
    unsigned int * sets;
    }IndependetSets;


//! \file independent_sets.h
//! \brief Functions needed for independent set preparations

//! Generate the independet set information for each poly_type.
//!
//! \param p Phase where sets are calculated
//! \return Errorcode
int generate_independet_sets(struct Phase*const p);

#endif//INDEPENDENT_SETS_H
