/* Copyright (C) 2016-2019 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016-2017 Marcel Langenberg
   Copyright (C) 2016 Fabien Leonforte
   Copyright (C) 2016 Juan Orozco
   Copyright (C) 2016 Yongzhi Ren
   Copyright (C) 2016 N. Harshavardhan Reddy

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



#ifndef _ANA_SERVER_H_
#define _ANA_SERVER_H_

#include <stdint.h>
#include <H5Ipublic.h>
#include "soma_config.h"
#include "ana_info.h"
#include "server.h"

int calc_MSD_server(const struct global_consts * gc, const Polymer * polymers,
                    size_t n_polymers, soma_scalar_t * result, const struct comm_with_info * comm, int target);

int analytics_server(const struct global_consts *gc,
                     const Ana_Info * ai, const struct server_info * si, const struct receiver *rcv,
                     unsigned int time);

int extent_ana_by_field(const soma_scalar_t * const data, const uint64_t n_data, const char *const name,
                        const hid_t file_id);
#endif //_ANA_SERVER_H_
