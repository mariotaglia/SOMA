/* This file is part of SOMA.

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

#ifndef SOMA_ELECTRIC_FIELD_H
#define SOMA_ELECTRIC_FIELD_H

#include "soma_config.h"
#include "soma_util.h"
#include "phase.h"


//! Top level struct for electric field implementation.
typedef struct ElectricField{
	soma_scalar_t *eps;				//!< Array that contains the dielectric constants for all particle types.
    soma_scalar_t *eps_arr;         //!< Array that saves information of calculated dielectric constant field.
    uint8_t *electrodes;            //!< Array that contains the electrode positions.
    uint8_t *iter_per_MC;           //!< Value that determines maximum amount of iterations to solve the electric field.
    uint8_t *iter_limit;            //!< Value that determines the upper limit of iterations to solve the electric field.
    soma_scalar_t *thresh_iter;     //!< Value that determines the threshold to stop iterative solution of the electric field.
    soma_scalar_t *Epot;            //!< Array that contains the electric potential field.
    soma_scalar_t *Epot_tmp;        //!< Temporary array that contains the electric potential field after MC step.
    soma_scalar_t *pre_deriv;       //!< Array that contains precomputed derivatives of dielectric constant field.
    soma_scalar_t *H_el_field;      //!< Array that contains electrostatic energy contribution per cell.
    soma_scalar_t H_el;             //!< Sum of electrocstatic contribution.
    
} ElectricField;

/*! Helper function to read the electric field array, electrode array and dielectric constants from the config HDF5 file.
    \private
    \param p Phase describing the system
    \param file_id File identifier of open HDF5 file.
    \param plist_id Access properties to use.
    \return Errorcode
*/
int read_electric_field_hdf5(struct Phase *const p, const hid_t file_id, const hid_t plist_id);

/*! Helper function to write the electric field array, electrode array and dielectric constants to the config HDF5 file.
    \private
    \param p Phase describing the system
    \param file_id File identifier of open HDF5 file.
    \param plist_id Access properties to use.
    \return Errorcode
*/
int write_electric_field_hdf5(const struct Phase *const p, const hid_t file_id, const hid_t plist_id);

/*! Helper function to resolve periodic boundaries for coordinates x,y,z.
    \private
    \param p Phase describing the system
    \param x x-coordinate in 3D representation of field.
    \param y y-coordinate in 3D representation of field.
    \param z z-coordinate in 3D representation of field.
*/
uint64_t cell_to_index(struct Phase *p, const uint64_t x, const uint64_t y, const uint64_t z);

/*! Helper function to compute the dielectric field from the densities of individual particle types (welling2014, eq. 83)
    \private
    \param p Phase describing the system
*/
void calc_dielectric_field(struct Phase *const p);

/*! Helper funtion to precompute derivatives for dielectric field (welling2014, eq. 85)
    \private
    \param p Phase describing the system
*/
void pre_derivatives(struct Phase *const p);

/*! Helper function to iterate solution for derivatives of the electric potential field (using precomputed derivatives for dielectric field)
    \private
    \param p Phase describing the system
    \returns maximum of dE_pot/dr
*/
soma_scalar_t iterate_field(struct Phase *const p);

/*! Helper function to compute the square of the derivative of the electric potential field
    \private
    \param p Phase describing the system
    \returns |Epot/dr|^2 

soma_scalar_t Epot_deriv_sq(struct Phase *const p, uint64_t x, uint64_t y, uint64_t z) */

/*! Helper function to compute the partial derivative of the electric potential field with respect to x
    \private
    \param p Phase describing the system
    \param x x-coordinate in 3D representation of field.
    \param y y-coordinate in 3D representation of field.
    \param z z-coordinate in 3D representation of field.
    \returns dEpot/dx
*/
soma_scalar_t dEpotx(struct Phase *const p, const uint64_t x, const uint64_t y, const uint64_t z);

/*! Helper function to compute the partial derivative of the electric potential field with respect to y
    \private
    \param p Phase describing the system
    \param x x-coordinate in 3D representation of field.
    \param y y-coordinate in 3D representation of field.
    \param z z-coordinate in 3D representation of field.
    \returns dEpot/dy
*/
soma_scalar_t dEpoty(struct Phase *const p, const uint64_t x, const uint64_t y, const uint64_t z);

/*! Helper function to compute the partial derivative of the electric potential field with respect to z
    \private
    \param p Phase describing the system
    \param x x-coordinate in 3D representation of field.
    \param y y-coordinate in 3D representation of field.
    \param z z-coordinate in 3D representation of field.
    \returns dEpot/dz
*/
soma_scalar_t dEpotz(struct Phase *const p, const uint64_t x, const uint64_t y, const uint64_t z);

/*! Main routine, calculates electrostatic energy contribution per cell and total (welling2017, eq. 4)
    \private
    \param p Phase describing the system
    \returns Errorcode
*/
int calc_electric_field_contr(struct Phase *const p);

/*! Helper function to free the CPU memory resources of the pc struct. The function gets automatically called by free_phase().
  \private
  \param p Initialized Phase that is in the process of deallocating its resources.
  \return Errorcode
  */
int free_electric_field(struct Phase *const p);

#endif                          //SOMA_ELECTRIC_FIELD_H