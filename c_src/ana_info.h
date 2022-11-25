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
#ifndef ANA_INFO_H
#define ANA_INFO_H
#include "soma_config.h"
#include <hdf5.h>
#include "soma_config.h"
struct Phase;

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
    unsigned int delta_mc_Epot_field;           //!<\#mc_sweeps between the ana of the electrostatic potential field
    unsigned int delta_mc_E_field;              //!<\#mc_sweeps between the ana of the electric field
    unsigned int delta_mc_dielectric_field;     //!<\#mc_sweeps between the ana of the dielectric field
    unsigned int delta_mc_amount_iter;          //!<\#mc_sweeps between the ana of the amount of iterations to solve the electric field
    soma_scalar_t *q_static;    //!< absolute value of wave vector to be calculated for static structure factor
    soma_scalar_t *q_dynamical; //!< absolute value of wave vector to be calculated for dynamical structure factor
    unsigned int q_size_static; //!< number of wave vector to calculate for static structure factor
    unsigned int q_size_dynamical;      //!< number of wave vector to calculate for dynamical structure factor
    char *filename;             //!< filename of the analysis file.
    char *coord_filename;       //!< filename of the configuration files.
    hid_t file_id;              //!< HDF5 file specifier for the ana file. Only valid for current_core == 0. No MPI/IO
#if (ENABLE_MPI == 1)
    MPI_Comm inter_domain_communicator; //!< communicator that enables communication between different sim ranks to store density fields
#endif                          //ENABLE_MPI
} Ana_Info;

//! \brief Initialization of the information needed for analysis routines.
//! \param p almost fully initialized Phase pointer.
//! \param filename Filename used for output of the analysis. If
//! pointer to "" or NULL is passed analysis is turned off.
//! \param coord_filename Filename of the readin coord.h5 file.
//! This  used to setup proper dumping filenames. (Passing NULL disables dumping.)
//! \post initialized p->ana_info
//! \return Errorcode.
int init_ana(struct Phase *const p, const char *const filename, const char *const coord_filename);

//! Release resources connected to the Ana_Info struct after init_ana.
//!
//! \param a Pointer to Ana_Info to release
//! \return errorcode
int close_ana(struct Ana_Info *const a);

#endif                          //ANA_INFO_H
