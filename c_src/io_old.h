/* Copyright (C) 2016-2021 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016-2017 Marcel Langenberg

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

#ifndef SOMA_IO_OLD_H
#define SOMA_IO_OLD_H

#include "soma_config.h"
#include <hdf5.h>

struct Phase;

//! \file io_old.h
//! \brief File compiling functions for handling of older SOMA file formats

/*! Parser for old style coord.dat files
  \deprecated
  \param p System where read data is stored.
  \param filename Filename to read in.
  \return Errorcode
*/
int read_old_config(struct Phase *p, char *const filename);

/*! Read in old geometry files.
 * \deprecated
 * \param p Phase which is going to get the field added.
 * \param filename Filename containing the info.
 * \return Errorcode
 */
int read_old_geometry(struct Phase *p, const char *filename);

/*! Read beads from configuration file with deprectated version 0
 * \deprecated
 * \param p Phase into which the beads are read
 * \param file_id open HDF5 file
 * \param plist_id parameter for the access of the file.
 * \return Errorcode
 */
int read_beads0(struct Phase *const p, const hid_t file_id, const hid_t plist_id);

#endif                          //SOMA_IO_OLD_H
