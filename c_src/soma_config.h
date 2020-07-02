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

//! \file soma_config.h
//! \brief configuration variables for SOMA File is configured by CMAKE.

//! Macro to destinguish between SINGLE_PRECISION and DOUBLE.
//! Automatically set by CMake
#ifndef SOMA_CONFIG_H
#define SOMA_CONFIG_H

#define SINGLE_PRECISION 0

//! Macro to activate MIC support
//! Automatically set by CMake
#define ENABLE_MIC 0

//! Macro to activate domain decomposition support
//! Automatically set by CMake
#define ENABLE_DOMAIN_DECOMPOSITION 1

//! Macro indicating an MPI build
#define ENABLE_MPI 1

#if ( SINGLE_PRECISION == 1)

//! Alias for float or double variables
typedef float soma_scalar_t;
//! Alias for float or double variables in HDF5 memory
#define H5T_SOMA_NATIVE_SCALAR H5T_NATIVE_FLOAT
//! Alias for float or double variables in HDF5 files
#define H5T_SOMA_FILE_SCALAR H5T_IEEE_F32LE
#if ( ENABLE_MPI == 1 )
//! Alias for flow or double variable in MPI
#define MPI_SOMA_SCALAR MPI_FLOAT
#else                           //ENABLE_MPI
//! MPI_DUMMY
#define MPI_SOMA_SCALAR -1
#endif                          //ENABLE_MPI

#else                           //SINGLE_PRECISION

//! Alias for float or double variables
typedef double soma_scalar_t;
//! Alias for float of double variables in HDF5 memory
#define H5T_SOMA_NATIVE_SCALAR H5T_NATIVE_DOUBLE
//! Alias for float of double variables in HDF5 files
#define H5T_SOMA_FILE_SCALAR H5T_IEEE_F64LE
#if ( ENABLE_MPI == 1 )
//! Alias for flow or double variable in MPI
#define MPI_SOMA_SCALAR MPI_DOUBLE
#else                           //ENABLE_MPI
//! MPI_DUMMY
#define MPI_SOMA_SCALAR -1
#endif                          //ENABLE_MPI

#endif                          //SINGLE_PRECISION

//! String containing the git version of SOMA
static const char soma_version[] = "0.6.1-32-ge9ff09e MPI DOUBLE";
//! Returns the version string.
//! \return Pointer to version string
const char *get_soma_version(void);

//! String describing the system info for which SOMA has been compiled.
static const char soma_system_info[] = "Linux-5.6.18-200.fc31.x86_64 x86_64 GNU";
//! Returns the string describing the system SOMA has been compiled.
//! \return Pointer to string.
const char *get_soma_system_info(void);

#endif                          //SOMA_CONFIG_H

//Check for OMP or OPENACC usage
#ifdef _OPENACC
#ifdef _OPENMP
#error "You could either compile with OPENACC or OMP. Not with both."
#endif                          //_OPENMP
#endif                          //_OPENACC
