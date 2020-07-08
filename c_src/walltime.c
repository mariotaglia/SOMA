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

#include "mpi.h"
#include "walltime.h"
#include "err_handling.h"
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <stdio.h>

int check_walltime_stop(void)
{
    char *env = getenv("SOMA_WALLTIME_STOP");
    if (env != NULL)
        {
            const long int walltime = atol(env);
            const time_t ac_time = time(NULL);
            if (walltime != 0)
                if (ac_time > walltime)
                    return 1;
        }
    return 0;
}

bool need_walltime_stop(bool no_sync_signal_flag, unsigned int timestep)
{
    bool stop_iteration = check_walltime_stop();
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    //todo: with this, walltime-stop is always no-sync. We need to define what synched walltime-stop precisely means if server and simranks aren't necessarily on the same timestep.
    (void)no_sync_signal_flag;
    /*
    if (!no_sync_signal_flag)
        {
            MPI_Allreduce(MPI_IN_PLACE, &stop_iteration, 1, MPI_C_BOOL,
                MPI_LOR, MPI_COMM_WORLD);
        }

    if (stop_iteration && world_rank == 0 && !no_sync_signal_flag)
        {
            fprintf(stdout, "One or more ranks have detected the need to do a walltime-stop at timestep %ud, due to synching, all will stop\n",
                timestep);
            return;
        }*/
    if (stop_iteration)
        {
            fprintf(stdout, "Environment to stop iteration at time %ud. "
                            "This was detected by world-rank %d\n",
                            timestep, world_rank);
        }
    return stop_iteration;
}