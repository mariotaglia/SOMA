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
#ifndef ANA_INFO_H
#define ANA_INFO_H
#include "soma_config.h"
#include <hdf5.h>
#include "soma_config.h"
struct Phase;
struct global_consts;

//! \file ana_info.h
//! \brief Info needed for output routines.

//! Basic information for analysis.
typedef struct Ana_Info {
    unsigned int delta_mc_Re;   //!< \#mc_sweeps between the ana of Re
    unsigned int delta_mc_Rg;   //!< \#mc_sweeps between the ana of Rg
    unsigned int delta_mc_b_anisotropy; //!< \#mc_sweeps between the ana of b_aniostropy
    unsigned int delta_mc_density_field;        //!< \#mc_sweeps between the ana of density fields
    unsigned int delta_mc_acc_ratio;    //!< \#mc_sweeps between the ana of acc ratios
    unsigned int delta_mc_MSD;  //!< \#mc_sweeps between the ana of MSD
    unsigned int delta_mc_dump; //!< \#mc_sweeps between full dumps of the coords.
    unsigned int delta_mc_density_var;  //!< \#mc_sweeps between the ana of density variance
    unsigned int delta_mc_non_bonded_energy;    //!< \# mc_sweeps between the ana of the non-bonded energy
    unsigned int delta_mc_bonded_energy;        //!< \# mc_sweeps between the ana of the bonded energy
    unsigned int delta_mc_umbrella_field;       //!<\#mc_sweeps between the ana of the string fields
    unsigned int delta_mc_dynamical_structure;  //!<\#mc_sweeps between the ana of the dynamical structure factor
    unsigned int delta_mc_static_structure;     //!<\#mc_sweeps between the ana of the static structure factor
    soma_scalar_t *q_static;    //!< absolute value of wave vector to be calculated for static structure factor
    soma_scalar_t *q_dynamical; //!< absolute value of wave vector to be calculated for dynamical structure factor
    unsigned int q_size_static; //!< number of wave vector to calculate for static structure factor
    unsigned int q_size_dynamical;      //!< number of wave vector to calculate for dynamical structure factor
    char *filename;             //!< filename of the analysis file.
    char *coord_filename;       //!< filename of the configuration files.
    hid_t file_id;              //!< HDF5 file specifier for the ana file. Only valid for current_core == 0. No MPI/IO
} Ana_Info;

//! \brief Initialization of the information needed for analysis routines.
//! \param ana_info (out) returns a fully initialized ana_info
//! \param end_mono (out) if Re is needed, *end_mono stores an array of length (2 * n_poly_types)
//! if Re analysis is turned off (or found to be impossible) *end_mono is NULL.
//! This allocates the array on the heap, and you become the owner of *end_mono.
//! \param global_consts (in) global constants of the simulation as initialized by read_consts_from_config
//! \param field_scaling_type scaling factor according to density, must have gc->n_types elements (simulation ranks have this as p->field_scaling_type)
//! \param filename (in) Filename used for output of the analysis. If
//! pointer to "" or NULL is passed analysis is turned off.
//! \param coord_filename (in) Filename of the readin coord.h5 file.
//! This  used to setup proper dumping filenames. (Passing NULL disables dumping.)
//! \param writer_comm (in) The ranks in this communicator will have write-access to ana_info->file_id.
//! Ranks that pass MPI_COMM_NULL will receive an invalid file_id.
//! All ranks that want write-access must pass the same communicator.
//! \note this call is collective on world,
//! \post initialized p->ana_info
//! \return Errorcode.
int init_ana(Ana_Info * ana_info, unsigned int ** end_mono, const struct global_consts *gc, soma_scalar_t * field_scaling_type,  const char *const filename, const char *const coord_filename, MPI_Comm writer_comm);

//! Release resources connected to the Ana_Info struct after init_ana.
//!
//! \param a Pointer to Ana_Info to release
//! \return errorcode
int close_ana(struct Ana_Info *const a);

#endif                          //ANA_INFO_H
