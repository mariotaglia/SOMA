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

#include "walltime.h"
#include <stdlib.h>
#include <time.h>

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
