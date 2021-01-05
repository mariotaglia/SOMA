/* Copyright (C) 2016-2021 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016 Marcel Langenberg
   Copyright (C) 2016 Fabien Leonforte
   Copyright (C) 2016 Juan Orozco
   Copyright (C) 2016 Yongzhi Ren
   Copyright (C) 2016 N. Harshavardhan Reddy

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

/*! \file ana.h
  \brief Functions to calculate observables of the configurations.
*/

#ifndef SOMA_ANA_H
#define SOMA_ANA_H

#include "soma_config.h"
#include <hdf5.h>
#include "stdint.h"

struct Phase;

/*! \brief Structure factor type enumerator to indicate the type of
 *  structure factor calculated. It can be time dependent (dynamical) and time independent(static)*/
enum structure_factor_type { DYNAMICAL_STRUCTURE_FACTOR, STATIC_STRUCTURE_FACTOR };

/*!>\brief calculate the end-to-end distance for the polymers of the phase
  \param p Phase configuration to analyze
  \returns global \f$ Re \f$
  \param result pointer to return the result. The length of this array
  should be p->n_poly_types * 4. Order is Re2 Rex Rey Rez for each
  poly type.
*/
void calc_Re(const struct Phase *p, soma_scalar_t * const result);

/*!>\brief calculate the variance between the current density field and the density field at a past time
  \param p Phase configuration to analyze
  \returns global \f$ dvar \f$
  \param dvar pointer to return density variance
*/
void calc_dvar(const struct Phase *p, soma_scalar_t * dvar);

/*!>\brief calculate the gyration radius for the polymers of the phase
  \param p Phase configuration to analyze
  \param result pointer to return the result. The length of this array
  should be p->n_poly_types * 4. Order is Rg2 Rgx Rgy Rgz for each
  poly type.
  \returns global \f$ Rg^2 \f$ and its squared components independently of the polymer type
*/
void calc_Rg(const struct Phase *p, soma_scalar_t * const result);

/*!>\brief calculate the monomer and chain mean square displacement
  \param p Phase configuration to analyze
  \param result array to store the result. Size should be 8*n_poly_types.
  The order is: x-component of MSD, y-component of MSD,z-component of MSD, MSD in 3D and
  x-componenet of mass center MSD, y-componenet of mass center MSD, z-componenet of mass center MSD,
  mass center MSD.
*/
void calc_MSD(const struct Phase *p, soma_scalar_t * const result);

/*! \brief calculate bond anisotropy

  \param p Phase configuration to analyze.
  \param result Array to store the tensor. xx yy zz xy xz yz. Length of array must be 6*n_poly_types.
  \post the memory of where the pointers point is set to the result of calculation.
*/
void calc_anisotropy(const struct Phase *p, soma_scalar_t * const result);

//! \brief calculate the current acceptance ratio.
//! \param p System to analyze
//! \param acc_ratio Pointer to initialized soma_scalar_t where the result is going to be stored.
//! \note After every call the result is reset.
void calc_acc_ratio(struct Phase *const p, soma_scalar_t * const acc_ratio);

//! \brief calculate the non-bonded energy for each particle type
//! \param p System to analyze
//! \param non_bonded_energy Pointer to array of size p->n_types to store the result
void calc_non_bonded_energy(const struct Phase *const p, soma_scalar_t * const non_bonded_energy);

//! \brief calculate the bonded energy for each bond type
//! \param p System to analyze
//! \param bonded_energy Pointer to array of size NUMBER_SOMA_BOND_TYPES to store the result
void calc_bonded_energy(const struct Phase *const p, soma_scalar_t * const bonded_energy);

/*! \brief run the desired analytics

  \param p Configuration phase to analyze.
  \return Errorcode
*/
int analytics(struct Phase *const p);

//! \brief Helper to ouput soma_scalar_t data to a hdf5 file, by extending an \a existing dataset.
//!
//! \param data Pointer to the data to ouput.
//! \param n_data number of data elements.
//! \param name Dataset name
//! \param file_id File specifier for HDF5 output ana file
//! \return Errorcode.
int extent_ana_by_field(const soma_scalar_t * const data, const uint64_t n_data, const char *const name,
                        const hid_t file_id);

#if ( ENABLE_MPI != 1 )
//! Dummy MPI type definition
typedef int MPI_Datatype;
#define MPI_UINT16_T -1
#endif                          //ENABLE_MPI

//! Helper to ouput data to a hdf5 , may be used for the density fields.
//!
//! \note This function does not feature a parallel output of the
//! (density) field. If anyone knows how to fix it, you are welcome to
//! help.
//!
//! \param p Pointer to the state of the system.
//! \param field_pointer Pointer to the field that should be outputted
//! \param field_name Name of the field, eg. "\density_field"
//! \param hdf5_type The H5-Datatype for the output of the field
//! \param mpi_type The MPI_Datatype of the field
//! \param data_size sizeof( datatype ) to get the number of bytes
//! \return Errorcode.
//! \todo The ouput is not yet parallel.
int extent_density_field(const struct Phase *const p, void *const field_pointer,
                         const char *const field_name, hid_t hdf5_type,
                         const MPI_Datatype mpi_type, const size_t data_size);

//! \brief calculate the structure for each poly type
//! \param p System to analyze
//! \param result Pointer to array to store the result
//! \param sf_type Type of structure factor. 0 stands for dynamical, 1 stand for static
//! \return Errorcode.
//#pragma acc routine(calc_structure) seq
int calc_structure(const struct Phase *p, soma_scalar_t * const result, const enum structure_factor_type sf_type);

//! \brief Helper to ouput soma_scalar_t data to a hdf5 file, used for the strcture.
//!
//! \param p Pointer to the state of the system.
//! \param data Pointer to the data to ouput.
//! \param name Dataset name
//! \param file_id File specifier for HDF5 output ana file
//! \param sf_type Type of structure factor. 0 stands for dynamical, 1 stand for static
//! \return Errorcode.
int extent_structure(const struct Phase *p, const soma_scalar_t * const data, const char *const name,
                     const hid_t file_id, const enum structure_factor_type sf_type);

#endif                          //SOMA_ANA_H
