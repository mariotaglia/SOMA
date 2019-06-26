/* Copyright (C) 2016-2019 Ludwig Schneider

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

#ifndef SOMA_POLYTYPE_CONVERSION_H
#define SOMA_POLYTYPE_CONVERSION_H

#include "soma_config.h"
#include "soma_util.h"

//! Top level struct for polymer type conversions.
//! controls the execution frequency. All other fields are only valid, if deltaMC != 0
typedef struct PolyConversion{
    unsigned int deltaMC; //!< control execution frequency of the conversion
    uint8_t * array; //!< Array that contains the reaction start index of the conversion list.

    unsigned int * input_type; //!< Array that contains the input poly type for each reaction (educt)
    unsigned int * output_type;//!< Array that contains the output poly type for each reaction (product)
    unsigned int * reaction_end; //!< Array indicating if this is the last reaction in the list. (boolean)
    unsigned int len_reactions; //!< length of the reaction related arrays input_type, output_type and reaction_end
    }PolyConversion;


//! Helper function to copy the pc data to the device
//! \private
//! \param p Fully CPU initialized Phase struct
//! \return Errorcode
int copyin_poly_conversion(struct Phase* p);

//! Helper function delete the pc data from the device and copy it to the CPU memory
//! \private
//! \param p Fully CPU initialized Phase struct
//! \return Errorcode
int copyout_poly_conversion(struct Phase* p);

//! Helper function to update the host with the pc data
//! \private
//! \param p Fully initialized Phase struct
//! \return Errorcode
int update_self_poly_conversion(const struct Phase* const p);

/*! Helper function to read the conversion array from the config HDF5 file.
    \private
    \param p Phase describing the system
    \param file_id File identifier of open HDF5 file.
    \param plist_id Access properties to use.
    \return Errorcode
*/
int read_poly_conversion_hdf5(struct Phase *const p, const hid_t file_id, const hid_t plist_id);

/*! Helper function to write the polyconversion array to the config HDF5 file.
    \private
    \param p Phase describing the system
    \param file_id File identifier of open HDF5 file.
    \param plist_id Access properties to use.
    \return Errorcode
*/
int write_poly_conversion_hdf5(const struct Phase *const p, const hid_t file_id,const hid_t plist_id);

/*! Helper function to free the CPU memory resources of the pc struct. The function gets automatically called by free_phase().
  \private
  \param p Initialized Phase that is in the process of deallocating its resources.
  \return Errorcode
*/
int free_poly_conversion(struct Phase* p);

#endif//SOMA_POLYTYPE_CONVERSION_H
