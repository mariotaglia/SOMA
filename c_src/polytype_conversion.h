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

//! \file polytype_conversion.h
//! \brief File for all declarations of functions and structs for the on-the-fly conversion of polymer types.

//! Lower level for the conversion. Contains the box region, where conversions happen.
//! And describes which polymer types are converted to which others.
typedef struct BoxConversion{
    unsigned int  input_type; //!< polymer type which is converted
    unsigned int output_type; //!< polymer type to which the chains types are converted
    bool stop_iteration; //!< is this the last reaction of this cell?
    }BoxConversion;


//! Top level struct for polymer type conversions.
//! controls the execution frequency. All other fields are only valid, if deltaMC != 0
typedef struct PolyConversion{
    unsigned int deltaMC; //!< control execution frequency of the conversion
    uint16_t * array; //!< Array that contains the reaction start index of the conversion list.
    BoxConversion * reactions; //!< Array of reactions. The index is stored in the conversion array.
    }PolyConversion;


#endif//SOMA_POLYTYPE_CONVERSION_H
