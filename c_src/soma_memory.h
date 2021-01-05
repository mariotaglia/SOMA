/* Copyright (C) 2016-2021 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016-2017 Marcel Langenberg
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

/*! \file rng_alternative.h
  \brief Definition of alternative pseudo random numbers generators for SOMA.
  PCG32 is the default RNG generation engine, if a different one is used, the definitions for the alternative are in this file.
*/

#ifndef SOMA_MEMORY_H
#define SOMA_MEMORY_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

//! Helper structure to manage larger amounts of memory for multiple arrays or states
//!
//! This struct is intended to store the heavy data of beads, RNG_STATES, permutation arrays etc.
//! The memory can be dynammically appended. It is internal an householder with precached memory.
typedef struct SomaMemory {
    void *ptr;                  //!< Pointer to the memory
    uint64_t length;            //!< Maximum allocated length of the memory
    uint64_t used;              //!< Used memory (< length)
    bool device_present;        //!< Has the memory been allocated on the device ?
    size_t typelength;          //!< Length of each element i.e. sizeof(Monomer) for beads
} SomaMemory;

//! Initialize the memory arrays
//!
//! \param state SomaMemory to inialize
//! \param length Minimum reserved memory after init
//! \param typelength Length of a single element that is stored in the array
//! \return Errorcode
int init_soma_memory(struct SomaMemory *state, const uint64_t length, const size_t typelength);

//! Deallocate all memory of the SomaMemory struct
//!
//! \param state State where the memory shall be deallocated
//! \return Errorcode
int free_soma_memory(struct SomaMemory *state);

//! Register a new element in SomaMemory.
//!
//! \param state initialized state
//! \param n number of allocated slots on the array (used is increased by n)
//! The function reallocates memory automatically if necessary
//! \return offset for the global ptr
uint64_t get_new_soma_memory_offset(struct SomaMemory *state, const uint64_t n);

//! Copyin memory for SomaMemory
//!
//! \param state Initialized state
//! \return Errorcode
int copyin_soma_memory(struct SomaMemory *state);

//! Copyout memory for SomaMemory
//!
//! \param state Initalized state
//! \return Errorcode
int copyout_soma_memory(struct SomaMemory *state);

//! Update device memory for SomaMemory
//!
//! \param state Initalized state
//! \return Errorcode
int update_device_soma_memory(struct SomaMemory *state);

//! Update self memory for SomaMemory
//!
//! \param state Intialized state
//! \return Errorcode
int update_self_soma_memory(struct SomaMemory *state);

//! Reallocate the size of the soma memory
//!
//! \private
//! \param state State to reallocate
//! \param min_increase minimum number by which the array is increased. must be > 1
//! \return Errorcode
int reallocate_soma_memory(struct SomaMemory *state, const uint64_t min_increase);

#endif                          //SOMA_MEMORY_H
