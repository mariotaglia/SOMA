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

#ifndef SOMA_WALLTIME_H
#define SOMA_WALLTIME_H

#include <stdbool.h>

/* \brief check environment if simulation should be stopped.
   Reads environment variable SOMA_WALLTIME_STOP and compares to time(NULL)
   \return true if simulation should be stopped false otherwise
 */
int check_walltime_stop(void);

//! returns true if the iteration of SOMAs main loop needs to be stopped early because
//! of the walltime. If yes, this information is also printed out.
//! collective on world
//! \param no_sync_signal_flag if all ranks should be synchronized for this operation.
//! Note: All ranks must pass the same value for no_sync_signal_flag.
//! \param timestep current timestep (only for printing)
//! \return true if iteration needs to be stop
bool need_walltime_stop(bool no_sync_signal_flag, unsigned int timestep);

#endif
