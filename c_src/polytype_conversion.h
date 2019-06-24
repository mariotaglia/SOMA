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
    uint16_t * array; //!< Array that contains the reaction start index of the conversion list.

    unsigned int * input_type; //!< Array that contains the input poly type for each reaction (educt)
    unsigned int * output_type;//!< Array that contains the output poly type for each reaction (product)
    unsigned int * reaction_end; //!< Array indicating if this is the last reaction in the list. (boolean)
    BoxConversion * reactions; //!< Array of reactions. The index is stored in the conversion array.
    }PolyConversion;


//! Helper function to copy the pc data to the device
int copyin_poly_conversion(struct PolyConversion* pc);

//! Helper function delete the pc data from the device and copy it to the CPU memory
int copyout_poly_conversion(struct PolyConversion* pc);

//! Helper function to update the host with the pc data
int update_self_poly_conversion(struct PolyConversion* pc);
#endif//SOMA_POLYTYPE_CONVERSION_H
