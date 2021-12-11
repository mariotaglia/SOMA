/* Copyright (C) 2020-2021 Ludwig Schneider

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

//! \file self_documentation.h
//! \brief Define structs and function to enable user friendly documentation of the history and self of simulations.

#ifndef SELF_DOCUMENTATION_H
#define SELF_DOCUMENTATION_H
#include "soma_config.h"
#include <stdio.h>
#include <hdf5.h>

struct Phase;

//! Struct to store the history and self of the running simulation.
//!
//! This info is read of the old configuration file and updated with the next simulation info.
//! Every configuration file written during the simulation will contain this info.
typedef struct SelfDocumentation {
    char **simdoc;              //!< Array of string to self of the simulation. The first dimension counts every simulation start, the second is the string
    unsigned int Ndoc;          //!< length of the self_doc. How many simulation self-documenting strings are stored.
    char *simdoc_internal;      //!< internal memory to hold all the string in a single array
} SelfDocumentation;

//! Function to initialize the SelfDocumentation string.
//!
//! \param p Initialized Phase struct
//! \param filename Filename of HDF5, if NULL is passed a fresh SelfDocumentation is initialized
//! \param sd SelfDocumentation Pointer to struct to init
//! \return Errorcode
int init_self_documentation(struct Phase *p, char *filename, struct SelfDocumentation *sd);

//! Write the self documentation (including the history) to an already opened HDF5 file.
//!
//! If the documentation dataset is present in the file, unlink the old and generate a new one
//! \param sd Initialized SelfDocumentation to write to file
//! \param file_id File handle of the already opened HDF5 file
//! \param plist_id Property identifier to access the file
//! \return Errorcode
int add_self_documentation_to_hdf5(const SelfDocumentation * sd, const hid_t file_id, const hid_t plist_id);

//! Function to release resources allocated for the SelfDocumenation
//!
//! \param sd SelfDocumentation to release
//! \return Errorcode
int free_self_documentation(SelfDocumentation * sd);

//! Fuction to print the self documentation string to stdout.
//!
//! \param sd SelfDocumentation to print out
//! \param f File handle to print the information to
//! \return Errorcode
int print_self_documentation(SelfDocumentation * sd, FILE * f);

#endif                          //SELF_DOCUMENTATION_H
