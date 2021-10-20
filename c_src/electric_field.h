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

#include "soma_config.h"
#include "phase.h"


//! Top level struct for electric field implementation.
typedef struct ElectricField{
	soma_scalar_t *eps;				//!< Array that contains the dielectric constants for all particle types.
    soma_scalar_t *eps_arr;         //!< Array that saves information of calculated dielectric constant field.
    soma_scalar_t *eps_arr_tmp;     //!< Temporay array that saves information of calculated dielectric constant field after MC step.
    uint8_t *electrodes;            //!< Array that contains the electrode positions.
    soma_scalar_t *Epot;            //!< Array that contains the electric potential field.
    soma_scalar_t *Epot_tmp;        //!< Temporary array that contains the electric potential field after MC step.
    soma_scalar_t *pre_deriv;       //!< Array that contains precomputed derivatives of dielectric constant field.
    soma_scalar_t *H_el_field;       //!< Array that contains electrostatic energy contribution per cell.
    soma_scalar_t *H_el;            //!< Sum of electrocstatic contribution.
    
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


/*! Helper function to compute the dielectric field from the densities of individual particle types (welling2014, eq. 83)
    \private
    \param p Phase describing the system
*/
void calc_dielectric_field(struct Phase *const p)

/*! Helper funtion to precompute derivatives for dielectric field (welling2014, eq. 85)
    \private
    \param p Phase describing the system
*/
void pre_derivatives(struct PHASE *const p)

/*! Helper function to iterate solution for derivatives of the electric potential field (using precomputed derivatives for dielectric field)
    \private
    \param p Phase describing the system
*/
void iterate_field(struct PHASE *const p)

/*! Helper function to compute the square of the derivative of the electric potential field
    \private
    \param p Phase describing the system
    \returns |Epot/dr|^2 
*/
soma_scalar_t Epot_deriv_sq(struct PHASE *const p, uint64_t x, uint64_t y, uint64_t z)

/*! Main routine, calculates electrostatic energy contribution per cell and total (welling2017, eq. 4)
    \private
    \param p Phase describing the system
    \returns Errorcode
*/
int calc_electric_field_contr(struct PHASE *const p);