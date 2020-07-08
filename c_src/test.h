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

//! \file test.h
//! \brief Source for various function to check the consistency of the simulation.
#ifndef SOMA_TEST_H
#define SOMA_TEST_H

#include "soma_config.h"
#include <inttypes.h>
#include "server.h"

struct Phase;

//! Test the read and write functionality for a given phase.
//!
//! \param p Phase which is used for the testing.
//! \warning This can only succeed with suffiecient diskspace and of
//! the cmd tool h5diff is installed. /tmp/p1.h5 and /tmp/p2.h5 are
//! going to be accessed.
//! \return Errorcode.
int test_read_write_hdf5(const struct Phase *const p);

//! Test the particle types to be in bounds
//!
//! \param p Phase to test.
//! \return error code (wrong type).
int test_particle_types(const struct Phase *const p);

//! Test, whether the forbidden area51, is violated, by any particle position.
//!
//! \param p System to test.
//! \return Errorcode
int test_area51_violation(const struct Phase *const p);

//! Test, whether the forbidden area51, is violated.
//!
//! In contrast to test_area51_violation() this function checks,
//! whether a molecule penetrates through a forbidden area.
//! \param p System to test.
//! \return Errorcode
int test_area51_exact(const struct Phase *const p);

//! If independet sets are used, test if they are really independet.
//!
//! \param p Phase to check
//! \return Errorcode
int test_independet_sets(const struct Phase *const p);

//! Test if all the polymer chains are in their correct domain.
//!
//! \param p Phase to check
//! \return Errorcode
int test_chains_in_domain(struct Phase *const p);

//! do various tests about the consistency of field-sending to the server
//! \param p Phase of the system
//! \param sim_inf info about the simrank performing the test
//! \param min_cell first cell to be sent (index to fields_unified or omega_fields_unified array)
//! \param max_cell last cell to be sent (inclusive)
//! \param field_size number of cells to be sent
//! \note this call is collective on the Simulation-Communicator. Will abort the run with an error message if inconsistencies are detected.
void test_field_sending_consistency(const struct Phase *p, const struct sim_rank_info *sim_inf, uint64_t min_cell, uint64_t max_cell, uint64_t field_size);
#endif                          //SOMA_TEST_H
