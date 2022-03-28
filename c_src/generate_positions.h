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

#ifndef GENERATE_POSITIONS_H
#define GENERATE_POSITIONS_H
#include "soma_config.h"
struct Phase;

//! \file generate_positions.h
//! \brief Functions needed for the generation of new inital conditions.

//! Generate and overwrite all position of the beads.
//!
//! \param p System with beads to insert.
//! \return Errorcode
//! \note Configuration might not be equilibrated.
//! \warning If your configuration contains rings, the connecting
//! bonds of the rings are not repected. This may result in very
//! unfavorable configurations. For most systems this can be
//! equilibrated, but if your area51 has topological constraints,
//! equilibration might be impossible.
int generate_new_beads(struct Phase *const p);

//! Generate and overwrite all monomer types according to the architectures of the polymers.
//!
//! \param p System with monomer types to insert.
//! \return Errorcode
int generate_monomer_type_array(struct Phase *const p);

#endif                          //GENERATE_POSITIONS_H
