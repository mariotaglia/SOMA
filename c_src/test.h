/* Copyright (C) 2016-2018 Ludwig Schneider

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

struct Phase;

//! Test the read and write functionality for a given phase.
//!
//! \param p Phase which is used for the testing.
//! \warning This can only succeed with suffiecient diskspace and of
//! the cmd tool h5diff is installed. /tmp/p1.h5 and /tmp/p2.h5 are
//! going to be accessed.
//! \return Errorcode.
int test_read_write_hdf5(const struct Phase*const p);

//! Test the particle types to be in bounds
//!
//! \param p Phase to test.
//! \return error code (wrong type).
int test_particle_types(const struct Phase*const p);

//! Test, whether the forbidden area51, is violated, by any particle position.
//!
//! \param p System to test.
//! \return Errorcode
int test_area51_violation(const struct Phase * const p);

//! Test, whether the forbidden area51, is violated.
//!
//! In contrast to test_area51_violation() this function checks,
//! whether a molecule penetrates through a forbidden area.
//! \param p System to test.
//! \return Errorcode
int test_area51_exact(const struct Phase * const p);

//! If independet sets are used, test if they are really independet.
//!
//! \param p Phase to check
//! \return Errorcode
int test_independet_sets(const struct Phase*const p);
#endif//SOMA_TEST_H
