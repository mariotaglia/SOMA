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

#ifndef SOMA_UTIL_H
#define SOMA_UTIL_H

struct Phase;
struct Info_MPI;
#include "soma_config.h"
#include <stdint.h>
#include "cmdline.h"
#include <math.h>
#include "bond.h"
#include "monomer.h"

//! \file soma_util.h
//! \brief File collecting several helper functions.

//! Enum to select the hamiltonian for the non-bonded interaction
enum Hamiltonian {
    SCMF0 = 0,                  //!< Original SCMF hamiltonian. For details refer doc of update_omega_fields_scmf0().
    SCMF1 = 1,                  //!< Alternative SCMF hamiltonian, especially for more than 2 types. For details refer doc of update_omega_fields_scmf1().
};

//! Enum to classify different cartesian directions
enum CartesionDirection {
    X = 0,                      //!< X direction
    Y = 1,                      //!< Y direction
    Z = 2,                      //!< Z direction
};

//! Enum to specify the different option of mobility modifier functions
enum MobilityEnum {
    DEFAULT_MOBILITY = 0,       //!< Default unchanged mobility. pacc multiplied with 1. Always.
    MULLER_SMITH_MOBILITY = 1,  //!< Mobility reduction factor of the shape \f$ m(r) = \min(1,\sum_\alpha a_\alpha \phi_\alpha(r) + \sum_alpha b_\alpha \phi_alpha^2(r)) \f$. With parameters a and b to be set via xml/hdf5.
    //! Mobility reduction with the shape \f$ M_i(\{\phi_j\}) = 1/2 \left( 1 + \tanh\left(\frac{\phi_{0,i} - \sum_j a_{ij} \phi_j}{\Delta \phi_{i}}\right)\right) \f$
    TANH_MOBILITY = 2,
};
//!  Function to extract the bond_type for poly_arch elements.
//!
//! \param info poly_arch element.
//! \returns Bond type in a bond list part.
//! \warning Call only for bond list part of poly_arch.
#pragma acc routine seq
unsigned int get_bond_type(const uint32_t info);

//!  Function to extract the offset you apply to your Monomer to get the bonded neigbour.
//!
//! \param info poly_arch element.
//! \returns Offset to apply to your Monomer to get the bonded neighbour.
//! \warning Call only for bond list part of poly_arch.
//! \note We use here in32_t instead uint32_t, because for signed
//! values, we need to set shift correct.
#pragma acc routine seq
int get_offset(const int32_t info);
//! Function to extract the end flag of a poly_arch element.
//!
//! \param info poly_arch element.
//! \returns End flag of the bond list. Stop iteration over the bond list, if it is non zero.
//! \warning Call only for bond list part of poly_arch.
#pragma acc routine seq
unsigned int get_end(const uint32_t info);
//! Function to compose a poly_arch element in the bond list region
//!
//! \param offset Offset to the next neighbor.
//! \param bond_type Type of the bond.
//! \param end Singal end of bond list.
//! \returns Element for a poly_arch bond list region.
#pragma acc routine seq
uint32_t get_info(const int offset, const unsigned int bond_type, const unsigned int end);
//! Get the offset for the poly_arch array to start the bond list iteration.
//!
//! \param info_bl poly_arch element in the Monomer region.
//! \return Offset to start the bond list iteration.
//! \warning Call only for Monomer region of poly_arch.
#pragma acc routine seq
int get_bondlist_offset(const int32_t info_bl);
//! Get particle type from a poly_arch element or the  Monomer type polymer heavy struct, depending on where it is stored.
//!
//! \param p Phase
//! \param i local polymer index
//! \param j monomer index
//! \return Particle type of the monomer.
#pragma acc routine seq
unsigned int get_particle_type(const struct Phase *const p, const uint64_t i, const unsigned int j);
//! Get particle type from a poly_arch element or the  Monomer region.
//!
//! \param info_bl poly_arch element in the Monomer region.
//! \return Particle type of the monomer.
//! \warning Call only for Monomer region of poly_arch.
#pragma acc routine seq
unsigned int get_particle_type_of_poly_arch(const uint32_t info_bl);
//! Compose poly_arch element for the
//!
//! \param offset_bl Offset to bondlist to set.
//! \param type Partile type to set.
//! \return Element for poly_arch of the Monomer region.
#pragma acc routine seq
uint32_t get_info_bl(const unsigned int offset_bl, const unsigned int type);

//! \brief after argument parsing of SOMA, this function interpret contradictions and warnings for the user.
//!
//! \param args Argument struct to interpret.
//! \param world_rank World rank of the calling rank.
//! \return Error code
int post_process_args(struct som_args *args, const unsigned int world_rank);

//! \brief Returns the number of bonds of specific type in the poly_arch structure.
//! \param p Phase of the system.
//! \param btype Bondtype to count
//! \warning This does not count the number of bonds in the system, only which bond types are used in the polyarch structure.
//! \return Number of bonds in the poly_arch structure.
//! You may you that function to check if a specific bond type is present in the system.
unsigned int get_number_bond_type(const struct Phase *const p, const enum Bondtype btype);

//! Reseed the random number generators.
//!
//! \param p System to reseed.
//! \param seed New seed for the system.
//! \return Errorcode
int reseed(struct Phase *const p, const unsigned int seed);

//! Macro that helps to check and ana for hdf5 errors.
#define HDF5_ERROR_CHECK(status) if(status < 0){fprintf(stderr, "ERROR: HDF5 %s:%s:%d: %d\n",__func__, __FILE__, __LINE__,(int)status);return status;}

//! Macro to check for errors. An error message is printed, with
//! linenumber and file. In addition, you can specify a specific
//! message that is printed.
#define HDF5_ERROR_CHECK2(status,name) if(status < 0){fprintf(stderr, "ERROR: HDF5 Name:%s %s:%d: %d\n",name, __FILE__, __LINE__,(int)status);return status;}

//! Macro to abort MPI and print a message with you
//! specified if status != 0. \warning This macro includes finalize_MPI() and exit();
#define MPI_ERROR_CHECK(status,msg) if(status != 0){fprintf(stderr, "ERROR: MPI abort Name: %s %s:%d: %d\n",msg, __FILE__, __LINE__,(int)status); ;exit(status);}

//! Macro to check and return error code if malloc failed.
#define MALLOC_ERROR_CHECK( ptr, size ) if( (ptr) == NULL){fprintf(stderr,"MALLOC-ERROR: %s:%d size = %lu\n", __FILE__, __LINE__, (uint64_t) (size)); return -1;}
#pragma acc routine(calc_bond_length) seq
static inline soma_scalar_t calc_bond_length(const soma_scalar_t x_i, const soma_scalar_t x_j, const soma_scalar_t box,
                                             const int mic);
inline soma_scalar_t calc_bond_length(const soma_scalar_t x_i, const soma_scalar_t x_j, const soma_scalar_t box,
                                      const int mic)
{
    soma_scalar_t r = x_i - x_j + 0 * mic * box;        // after "+" to shut up compiler warnings
#if ( ENABLE_MIC == 1)
    if (mic)
        {
#if ( SINGLE_PRECISION == 1)
            const soma_scalar_t img = rintf(r / box);
#else
            const soma_scalar_t img = rint(r / box);
#endif                          //SINGLE_PRECISION
            r -= img * box;
        }
#endif                          //ENABLE_MIC
    return r;
}

#endif                          //SOMA_UTIL_H
