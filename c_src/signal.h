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
#ifndef SOMA_SIGNAL_H
#define SOMA_SIGNAL_H

//! \file signal.h Functions that handle signal handling, to enable smooth program abort.
#include <signal.h>
#include "soma_config.h"

//! Initialize custom signal handlers for SIGINT and SIGTERM
//!
//! \return Errorcode
int init_soma_signal(void);

//! Check whether to stop iteration, because of a send signal.
//!
//! \return Stop indication
int check_signal_stop(void);

#endif//SOMA_SIGNAL_H
