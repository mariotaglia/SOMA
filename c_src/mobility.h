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

#ifndef SOMA_MOBILITY_H
#define SOMA_MOBILITY_H

#include "soma_config.h"
#include "soma_util.h"

//! \file mobility.h
//! \brief Function declaration and struct to modify the particle mobility based on the density composition.

//! Struct to bundle the mobility modifier parameter
typedef struct Mobility {
    enum MobilityEnum type;     //!< Type of the mobility modifier that need to be applied.
    unsigned int param_len;     //!< Length of the parameter array
    soma_scalar_t *param;       //!< Parameter array for the mobility calculation
    unsigned int *poly_type_mc_freq;    //!< Array that contains the execution frequency for different polymer types.
} Mobility;

//! Helper function to copy the mobility data to the device
//! \private
//! \param p Fully CPU initialized Phase struct
//! \return Errorcode
int copyin_mobility(struct Phase *p);

//! Helper function delete the mobility data from the device and copy it to the CPU memory
//! \private
//! \param p Fully CPU initialized Phase struct
//! \return Errorcode
int copyout_mobility(struct Phase *p);

//! Helper function to update the host with the mobility data
//! \private
//! \param p Fully initialized Phase struct
//! \return Errorcode
int update_self_mobility(const struct Phase *const p);

/*! Helper function to read the mobility from the config HDF5 file.
    \private
    \param p Phase describing the system
    \param file_id File identifier of open HDF5 file.
    \param plist_id Access properties to use.
    \return Errorcode
*/
int read_mobility_hdf5(struct Phase *const p, const hid_t file_id, const hid_t plist_id);

/*! Helper function to write the mobility to the config HDF5 file.
    \private
    \param p Phase describing the system
    \param file_id File identifier of open HDF5 file.
    \param plist_id Access properties to use.
    \return Errorcode
*/
int write_mobility_hdf5(const struct Phase *const p, const hid_t file_id, const hid_t plist_id);

/*! Helper function to free the CPU memory resources of the mobility struct. The function gets automatically called by free_phase().
  \private
  \param p Initialized Phase that is in the process of deallocating its resources.
  \return Errorcode
*/
int free_mobility(struct Phase *p);

/*! Actual function to calculate the mobility modifier.

   This function return a probability that needs be multiplied to acceptance function.

   \note this function calculates only the modifier for a single position.
   In order to get the full modification of pacc with respect to detailed balance, both old and new position have to be respected.
   The total modifier is than: \f$ m^* = \sqrt{ m(r) \cdot m(r+\Delta r)} \f$.

   \param p Fully initalized Phase struct
   \param particle_type Type of the moving particle
   \param x x coordinate of the spatial position
   \param y y coordinate of the spatial position
   \param z z coordinate of the spatial position

   \return modifier \f$ m(x,y,z) \in [0,1] \f$
*/
#pragma acc routine(get_mobility_modifier) seq
soma_scalar_t get_mobility_modifier(const struct Phase *const p, const unsigned int particle_type,
                                    const soma_scalar_t x, const soma_scalar_t y, const soma_scalar_t z);

#endif                          //SOMA_MOBILITY_H
