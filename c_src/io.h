/* Copyright (C) 2016-2021 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016 Marcel Langenberg

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
#ifndef SOMA_IO_H
#define SOMA_IO_H
/*! \file io.h
  \brief Header file for all functions, that handle with input and ouput
   operations of SOMA.
*/

/*!
  \brief Read a configuration from disk.

  And initializes all memory fields.
  \param p Phase struct which will define the state of the system
  after the call.
  \param filename Relative or absolute path the configuration file to
  read.
  \pre Uninitialied Phase* \a p struct (except the int init_MPI(Phase*) call).
  \post Initialized \a p points to an fully initialized configuration.
  \return Error code. Returns not equal to zero if an error occured.
*/
#include "soma_config.h"
#include <hdf5.h>
struct Phase;

/*!\brief Writes the current configuration to the disk.
  \param p Pointer to a fully initialized configuration.
  \param filename Relative or absolute path the configuration file to
  read.
  \warning The call will overwrite all previous information contained
  in the file "filename" is existing.
  \post Current state is written to the given file.
  \return Error code. Return not equal to zero if an error occured.
 */
int write_config(const struct Phase *const p, const char *const filename);

//! \brief Writes a configuration to disk in hdf5 format using parallel I/O
//!
//! \param p state of the system
//! \param filename Filename of the the ouput file.
//! \warning This function overwrites the file if it exists.
//! \return Errorcode.
int write_config_hdf5(struct Phase *const p, const char *filename);
//! \brief Reads a configuration to disk in hdf5 format using parallel I/O
//!
//! \param p state of the system
//! \param filename Filename of the the ouput file.
//! \return Errorcode.
int read_config_hdf5(struct Phase *const p, const char *filename);

//! \brief Ouput to stdout about the estimated time the program will finish.
//!
//! \param p System state.
//! \param Nsteps Total number of steps the simulation is running.
//! \return Errorcode
int screen_output(struct Phase *const p, const unsigned int Nsteps);

//! \brief Helper function to write HDF5, not parallel splitted data.
//!
//! \private Function is only for internal use.
//! \param ndims Dimensionality of data.
//! \param dims Lenght of the dimensions.
//! \param file_id hdf5-file identifier.
//! \param name Name of the targeted dataset.
//! \param file_type HDF5-type of the data in the file on the disc.
//! \param mem_type HDF5-type of the data in the memory.
//! \param plist_id Identifier for a property list, for accessing the file.
//! \param data pointer to the data to write.
//! \return Errorcode
//! \note Function is MPI-collective.
int write_hdf5(const hsize_t ndims, const hsize_t * const dims, const hid_t file_id,
               const char *const name, const hid_t file_type, const hid_t mem_type,
               const hid_t plist_id, const void *const data);

//! \brief Helper function to read HDF5, not parallel splitted data.
//!
//! \private Function is only for internal use.
//! \param file_id hdf5-file identifier.
//! \param name Name of the targeted dataset.
//! \param mem_type HDF5-type of the data in the memory.
//! \param plist_id Identifier for a property list, for accessing the file.
//! \param data pointer to the data to write.
//! \return Errorcode
//! \note Function is MPI-collective.
int read_hdf5(const hid_t file_id, const char *const name, const hid_t mem_type, const hid_t plist_id,
              void *const data);

#endif                          //SOMA_IO_H
