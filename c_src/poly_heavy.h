/* Copyright (C) 2016-2021 Ludwig Schneider
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

#ifndef POLYMER_HEAVY_H
#define POLYMER_HEAVY_H

struct Phase;
#include <stdbool.h>
#include "soma_config.h"
#include "soma_util.h"
#include "soma_memory.h"

/*!
  Struct to combine all heavy data pointer for the polymers of this rank.
  Each individual polymer contains offset values that clarify which area of the memory belongs to the respective polymer.
  It is important to ensure that these areas are initialised to not overlap.

  Advantage of these large memory blocks, instead of small pointer for each polymers is, that they can be more easily allocated and transferred to GPUs and/or other ranks.
  The use of offset integer, instead of pointers, referring to the large memory block is, that we can avoid any kind of deep-copy issues with OpenACC.

  Pointers that are not used, because the current execution mode doesn't need them are set to NULL.
 */
typedef struct PolymerHeavy {
    SomaMemory beads;           //!< Memory handle struct for the beads data
    SomaMemory msd_beads;       //!< Memory handle struct for the msd_beads data (typically not present on device)
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    SomaMemory monomer_types;   //!< Memory handle struct for the monomer type if monomer conversion is switched on.
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    SomaMemory set_states;      //!< Memory handle struct for set states memory
    SomaMemory set_permutation; //!< Memory handle struct for the set permutation info
} PolymerHeavy;

/*! \brief Frees all resources of the PolymerHeavy struct

  \param p Phase that contains a PolymerHeavy to free
  \return error code.
*/
int free_polymer_heavy(struct Phase *const p);

/*! \brief Initalizes and copies all data of the PolymerHeavy struct to device memory

  \param p Phase that contains a PolymerHeavy to initialize
  \return error code.
*/
int copyin_polymer_heavy(struct Phase *const p);

/*! \brief Copies all device data to self and frees devices resources

  \param p Phase that contains a PolymerHeavy to free
  \return error code.
*/
int copyout_polymer_heavy(struct Phase *const p);

/*! \brief Update the device memory with the current self memory of the polymer heavy struct

  \param p Phase that contains a PolymerHeavy to update
  \param rng_flag Update of RNG info required?
  \return error code.
*/
int update_device_polymer_heavy(struct Phase *const p, const bool rng_flag);

/*! \brief Updates the self memory with the current device memory of the polymer heavcy strcut.

  \param p Phase that contains a PolymerHeavy to update
  \param rng_flag Update of RNG info required?
  \return error code.
*/
int update_self_polymer_heavy(struct Phase *const p, const bool rng_flag);

/*! \brief If the allocated memory twice as large as necessary, create a new fresh polymer heavy struct

  \param p Phase that contains a PolymerHeavy to update
  \param collective Is this function collectively called for the simulation SOMA_comm_sim.
  \return error code.
*/
int consider_compact_polymer_heavy(struct Phase *p, const bool collective);
#endif                          // POLYMER_HEAVY_H
