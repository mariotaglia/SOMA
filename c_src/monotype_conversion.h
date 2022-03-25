/* Copyright (C) 2016-2022 Gregor Ibbeken, Ludwig Schneider

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

#ifndef SOMA_MONOTYPE_CONVERSION_H
#define SOMA_MONOTYPE_CONVERSION_H

#include "soma_config.h"
#include "soma_util.h"

//! Top level struct for monomer type conversions.
//! controls the execution frequency. All other fields are only valid, if deltaMC != 0
typedef struct MonoConversion {
    unsigned int deltaMC;       //!< control execution frequency of the conversion
    uint8_t *array;             //!< Array that contains the reaction start index of the conversion list.

    unsigned int *input_type;   //!< Array that contains the input mono type for each reaction (educt)
    unsigned int *output_type;  //!< Array that contains the output mono type for each reaction (product)
    unsigned int *reaction_end; //!< Array indicating if this is the last reaction in the list. (boolean)
    soma_scalar_t *rate;        //!< control execution probability of the conversion
    unsigned int block_size;    //!<Contains the block size of all monomer conversions. If it is greater than one, conversion happens in blocks. WARNING: No checks whether blocks fit into the relevant polymers.
    unsigned int *dependency_ntype;     //!<Array that contains the number of  dependency indices
    unsigned int *dependency_type_offset;       //!<Array that contains the start/offset of dependency indices
    unsigned int *dependency_type;      //!<Array that contains the dependency types
    unsigned int len_reactions; //!< length of the reaction related arrays input_type, output_type and reaction_end
    unsigned int len_dependencies;      //!< length of the density dependency array dependency_type (=sum over dependency_ntype) 

} MonoConversion;

//! Helper function to copy the mtc data to the device
//! \private
//! \param p Fully CPU initialized Phase struct
//! \return Errorcode
int copyin_mono_conversion(struct Phase *p);

//! Helper function delete the mtc data from the device and copy it to the CPU memory
//! \private
//! \param p Fully CPU initialized Phase struct
//! \return Errorcode
int copyout_mono_conversion(struct Phase *p);

//! Helper function to update the host with the mtc data
//! \private
//! \param p Fully initialized Phase struct
//! \return Errorcode
int update_self_mono_conversion(const struct Phase *const p);

/*! Helper function to read the conversion array from the config HDF5 file.
    \private
    \param p Phase describing the system
    \param file_id File identifier of open HDF5 file.
    \param plist_id Access properties to use.
    \return Errorcode
*/
int read_mono_conversion_hdf5(struct Phase *const p, const hid_t file_id, const hid_t plist_id);

/*! Helper function to write the monoconversion array to the config HDF5 file.
    \private
    \param p Phase describing the system
    \param file_id File identifier of open HDF5 file.
    \param plist_id Access properties to use.
    \return Errorcode
*/
int write_mono_conversion_hdf5(const struct Phase *const p, const hid_t file_id, const hid_t plist_id);

/*! Helper function to free the CPU memory resources of the mtc struct. The function gets automatically called by free_phase().
  \private
  \param p Initialized Phase that is in the process of deallocating its resources.
  \return Errorcode
*/
int free_mono_conversion(struct Phase *p);

/*! Convert monomer types according to the reaction description of the MonoConversion struct.
  This updates the center of mass of the monomers and chooses between full or partial (with rates) conversions.
  \param p Phase struct describing the simulation
  \return Errorcode
*/
int convert_monotypes(struct Phase *p);

/*! Fully convert monomer types according to the reaction description of the MonoConversion struct.
  \param p Phase struct describing the simulation
  \return Errorcode
*/
int fully_convert_monotypes(struct Phase *p);

/*! Partially Convert monomer types according to the reaction description of the MonoConversion struct.
 This converts the monomer only with a probability given by the rate which may depend (linearly) on the normalized density of some type (for reactions involving multiple types).
  \param p Phase struct describing the simulation
  \return Errorcode
*/
int partially_convert_monotypes(struct Phase *p);

/*! Completely convert a polymer according to the conversion rules if the metropolis criterion is met.
  Rates and dependencies are ignored, conversions take place depending on the energy  difference applied during the conversion.
  \param p Phase struct describing the simulation
  \return Errorcode
*/
int perform_semi_gc_conversions(struct Phase *p);
#endif                          //SOMA_MONOTYPE_CONVERSION_H
