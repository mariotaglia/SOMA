/* Copyright (C) 2016-2019 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016 Marcel Langenberg
   Copyright (C) 2016 Fabien Leonforte
   Copyright (C) 2016 Juan Orozco
   Copyright (C) 2016 Yongzhi Ren

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
#ifndef SOMA_INIT_H
#    define SOMA_INIT_H
#    include "soma_config.h"
/*! \file init.h
  \brief Header file functions required for initialization processes.
*/

//! Forward declaration of the Phase struct. To avoid inclusion of struct.h.
struct Phase;

//!  Print version of SOMA and libraries linked to SOMA
//! \param rank MPI rank of the calling process.
//! \return 0 in successfull non-zero otherwise.
int print_version(const int rank);

//! Initialize the p->msd polymer data with current positions.
//!
//! This will be done after the particles have been initialized,
//! either by reading an input file, or by randomly genrating new
//! particles.
//! \param p System description.
//! \return Errorcode
int init_msd(struct Phase *const p);

//! \brief setup the OpenACC devices according to the commandline arguments.
//! \param p System for which the devices are set.
//! \return Errorcode
int set_openacc_devices(const struct Phase *const p);

#endif                          //SOMA_INIT_H
