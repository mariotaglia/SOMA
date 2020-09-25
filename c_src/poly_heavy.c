/* Copyright (C) 2016-2019 Ludwig Schneider
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

/*! \file poly_heavy.c
\brief implementation of poly_heavy.h
*/

#include "poly_heavy.h"
#include "soma_memory.h"
#include "phase.h"

int free_polymer_heavy(struct Phase *const p)
{
    int status = 0;
    status |= free_soma_memory(&(p->ph.beads));
    status |= free_soma_memory(&(p->ph.msd_beads));
    if (p->args.iteration_alg_arg == iteration_alg_arg_SET)
        {
            status |= free_soma_memory(&(p->ph.set_states));
            status |= free_soma_memory(&(p->ph.set_permutation));
        }
    return status;
}

int copyin_polymer_heavy(struct Phase *const p)
{
    int status = 0;
    status |= copyin_soma_memory(&(p->ph.beads));
    status |= copyin_soma_memory(&(p->ph.msd_beads));
    if (p->args.iteration_alg_arg == iteration_alg_arg_SET)
        {
            status |= copyin_soma_memory(&(p->ph.set_states));
            status |= copyin_soma_memory(&(p->ph.set_permutation));
        }
    return status;
}

int copyout_polymer_heavy(struct Phase *const p)
{
    int status = 0;
    status |= copyout_soma_memory(&(p->ph.beads));
    status |= copyout_soma_memory(&(p->ph.msd_beads));
    if (p->args.iteration_alg_arg == iteration_alg_arg_SET)
        {
            status |= copyout_soma_memory(&(p->ph.set_states));
            status |= copyout_soma_memory(&(p->ph.set_permutation));
        }
    return status;
}

int update_device_polymer_heavy(struct Phase *const p, const bool rng_flag)
{
    int status = 0;
    status |= update_device_soma_memory(&(p->ph.beads));
    status |= update_device_soma_memory(&(p->ph.msd_beads));
    if (p->args.iteration_alg_arg == iteration_alg_arg_SET && rng_flag)
        {
            status |= update_device_soma_memory(&(p->ph.set_states));
            status |= update_device_soma_memory(&(p->ph.set_permutation));
        }
    return status;
}

int update_self_polymer_heavy(struct Phase *const p, const bool rng_flag)
{
    int status = 0;
    status |= update_self_soma_memory(&(p->ph.beads));
    status |= update_self_soma_memory(&(p->ph.msd_beads));
    if (p->args.iteration_alg_arg == iteration_alg_arg_SET && rng_flag)
        {
            status |= update_self_soma_memory(&(p->ph.set_states));
            status |= update_self_soma_memory(&(p->ph.set_permutation));
        }
    return status;
}
