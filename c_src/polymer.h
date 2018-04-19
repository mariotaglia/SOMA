/* Copyright (C) 2016-2017 Ludwig Schneider

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
#ifndef POLYMER_H
#define POLYMER_H

#include "rng.h"
#include "monomer.h"

//! \file polymer.h
//! \brief Code related to the Polymer structures

/*! \brief Polymer information */
typedef struct Polymer{
    //unsigned int N; /*!<\brief number of beads per chain */
    Monomer * beads; /*!<\brief pointer to beads */
    Monomer * msd_beads; //!< \brief bead positions used for MSD calculation. (Not on device.)
    unsigned int type;//!< \brief Type of the polymer architecture.
    RNG_STATE poly_state;   //!< \brief Struct which contains all RNGs
    struct RNG_STATE * set_states; //!< RNG states of independet sets. NULL if not used.
    //! Array to store thr permutation of sets for set iteration. NULL if not used.
    unsigned int*set_permutation;
    }Polymer;

//! \brief Deallocate memory of deep copy elements in Polymer.
//!
//! \param p initialized system.
//! \param poly Polymer to deallocate.
//! \return Errorcode.
int free_polymer(const struct Phase*const p, Polymer*const poly);

//! Copyin deep copy to device memory of a polymer.
//!
//! \param p System
//! \param poly Polymer to copyin.
//! \return Errorcode.
int copyin_polymer(struct Phase*const p, Polymer*const poly);

//! Copyout deep copy from device memory of a polymer.
//!
//! \param p System
//! \param poly Polymer to copyout.
//! \return Errorcode.
int copyout_polymer(struct Phase*const p, Polymer*const poly);

//! If more memory space for polymers is requested than available,
//! this functionallocates more space.
//!
//! \param p System to reallocate memory.
//! \return Errorcode
//! \note This function is expensive to call.
int reallocate_polymer_mem(struct Phase*const p);

//! Push a polymer to the end of the p->polymers array.
//!
//! \param p Fully initialized system.
//! \param poly Pointer to Polymer to insert.
//! If no storage for that polymer is available, a reallocate_polymer_mem()
//! is triggered.
//! \pre The polymer is has to be fully allocated on host memory. But not on GPU memory.
//! \post The polymer is part of the system. The system \a p is now owner of the polymer.
//! \post Device memory for the deep copy of polymer is initialized.
//! \warning If you change properties of the global system, you
//! need to call collective_global_update().
//! \return Errorcode.
int push_polymer(struct Phase*const p,const Polymer*const poly);

//! Extract a polymer from a position.
//!
//! \param p Fully initialized system.
//! \param poly_id index of the polymer to extract from the system.
//! \param poly Pointer to empty space to store the popep polymer in.
//! \pre memory space for the polymer is required.
//! \post \a poly is now the owner of the polymer.
//! \post \a poly is fully deacllocated from device memory.
//! \warning If you change properties of the global system, you
//! need to call collective_global_update().
//! \return Errorcode.
int pop_polymer(struct Phase*const p,const uint64_t poly_id,Polymer*const poly);

//! Obtain the number of bytes, which are necessary to serialize a polymer.
//!
//! Use this function to allocate memory for polymer serialization.
//! \param p System configuration.
//! \param poly Polymer to serialize.
//! \return Number of bytes.
unsigned int poly_serial_length(const struct Phase*const p,const Polymer*const poly);

//! Serialize an Polymer to a raw memory buffer.
//!
//! \param p System
//! \param poly Polymer to serialize.
//! \param buffer Preallocated buffer to store the outcome.
//! \pre Allocation of buffer with return value of poly_state_serial_length() minimum.
//! \note Ownership and allocation status is unchanged.
//! \return Number of written bytes. If < 0 Errorcode.
int serialize_polymer(const struct Phase*const p,const Polymer*const poly,unsigned char*const buffer);

//! Deserialize an Polymer from a raw memory buffer.
//!
//! \param p System.
//! \param poly Polymer to initialize by memory buffer.
//! \param buffer Initialized memory buffer to read.
//! \pre You are owner of \a polymer. And there is no deep
//! copy data allocated. Otherwise, you create memory leaks
//! , because deep copy data is allocated.
//! \post You are owner of the Polymer including deep copy data.
//! \return Number of written bytes. If < 0 Errorcode.
int deserialize_polymer(const struct Phase*const p, Polymer*const poly,unsigned char*const buffer);

//! Update the Self Memory of a given polymer
//!
//! \param p System
//! \param poly Polymer to update
//! \return Errorcode
int update_self_polymer(const struct Phase*const p,Polymer*const poly,const int rng_update_flag);


#endif//POLYMER_H
