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
#ifndef POLYMER_H
#define POLYMER_H
#include "soma_config.h"
#include "rng.h"
#include "monomer.h"

//! \file polymer.h
//! \brief Code related to the Polymer structures

/*! \brief Polymer information */
//! Instead of pointers this struct contains offset counters, that have to be applied to global arrays instead.
//! To invalidate pointers (such as using NULL) the offsets are sets are set to UINT64_MAX
//! To access for example the beads of this polymer use: p->ph.beads[ this->bead_offset + i]
typedef struct Polymer {
    unsigned int type;          //!< \brief Type of the polymer architecture.
    RNG_STATE poly_state;       //!< \brief Struct which contains all RNGs
    Monomer rcm;                //!< center of mass of the polymer
    uint64_t tag;               //!< unique tag of the polymer that remains constant even after sending around and storing in restart files.
    uint64_t bead_offset;       //!< \brief offset to the bead pointer for this polymer.
    uint64_t msd_bead_offset;   //!< \brief offset to the msd bead pointer for this polymer. (Typically not on device).
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    uint64_t monomer_type_offset;       //!< \brief offset to the monomer_types pointer for this polymer.
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    uint64_t set_states_offset; //!< offset to the set_states pointer for this polymer.
    uint64_t set_permutation_offset;    //!< offset to the set_permutation pointer for this polymer.
} Polymer;

//! If more memory space for polymers is requested than available,
//! this functionallocates more space.
//!
//! \param p System to reallocate memory.
//! \param new_storage Suggestion for new storage allocation. If smaller than heuristics, the heuristics is chosen.
//! \return Errorcode
//! \note This function is expensive to call.
int reallocate_polymer_mem(struct Phase *const p, uint64_t new_storage);

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
int push_polymer(struct Phase *const p, const Polymer * const poly);

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
int pop_polymer(struct Phase *const p, const uint64_t poly_id, Polymer * const poly);

//! Exchange the index of two polymers in p.
//!
//! \param p Fully initialized system.
//! \param poly_i the polymer to be exchanged with poly_j
//! \param poly_j the polymer to be exchanged with poly_i
//! \warning If you change properties of the global system, you
//! need to call collective_global_update().
//! \return Errorcode.
int exchange_polymer(struct Phase *const p, const uint64_t poly_i, const uint64_t poly_j);

//! Obtain the number of bytes, which are necessary to serialize a polymer.
//!
//! Use this function to allocate memory for polymer serialization.
//! \param p System configuration.
//! \param poly Polymer to serialize.
//! \return Number of bytes.
unsigned int poly_serial_length(const struct Phase *const p, const Polymer * const poly);

//! Serialize an Polymer to a raw memory buffer.
//!
//! \param p System
//! \param poly Polymer to serialize.
//! \param buffer Preallocated buffer to store the outcome.
//! \pre Allocation of buffer with return value of poly_state_serial_length() minimum.
//! \note Ownership and allocation status is unchanged.
//! \return Number of written bytes. If < 0 Errorcode.
int serialize_polymer(struct Phase *const p, const Polymer * const poly, unsigned char *const buffer);

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
int deserialize_polymer(struct Phase *const p, Polymer * const poly, const unsigned char *const buffer);

//! Update the center of mass of the polymer from its monomer positions.
//!
//! \param p System
//! \return Errorcode
int update_polymer_rcm(struct Phase *const p);

//! Get the domain id to which of the given position
//!
//! \param p reference system
//! \param rcm pointer to rcm
//! \return domain id
unsigned int get_domain_id(const struct Phase *const p, const Monomer * const rcm);
#endif                          //POLYMER_H
