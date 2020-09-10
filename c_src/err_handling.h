/* Copyright (C) 2016-2019 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016-2017 Marcel Langenberg

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

//! \file err_handling.h
//! \brief convenience macros/functions for common error handling

#ifndef _ERR_HANDLING_H_
#define _ERR_HANDLING_H_

#define RET_ERR_ON_NULL(ptr, msg) {\
if (ptr == NULL)\
        {\
            fprintf(stderr, "ERROR: %s %s:%d (function %s, %s is NULL)\n", msg, __FILE__, __LINE__, __func__, #ptr);\
            return -1;\
        }\
}

//! prints a mesage from rank 0 only.
//! if MPI is deactivated, it still prints the message, comm is ignored
//! without MPI, comm doesnt even have to be a valid symbol
#if (ENABLE_MPI == 1)
#define WARN_ONCE(comm, msg) {\
    int rank;\
    MPI_Comm_rank(comm, &rank);\
    if (rank == 0){\
        fprintf(stderr, "WARNING: (file %s line %d function %s): %s ", __FILE__, __LINE__, __func__, msg);\
    }\
}
#else
#define WARN_ONCE(comm, msg) {\
    fprintf(stderr, "WARNING: (file %s line %d function %s): %s ", __FILE__, __LINE__, __func__, msg);\
}
#endif

// works like printf (format string + arguments of any type to print),
// but flushes the output and adds additional information to the output:
// rank, line, file, function name
#define DPRINT(...) {int dprint_world_rank;\
	MPI_Comm_rank(MPI_COMM_WORLD, & dprint_world_rank);\
	fprintf(stdout, "worldrank %d line %d of file %s (function %s): ", dprint_world_rank, __LINE__, __FILE__, __func__);\
	fprintf(stdout, __VA_ARGS__);\
	fprintf(stdout,"\n");\
	fflush(stdout);}




#endif //_ERR_HANDLING_H_
