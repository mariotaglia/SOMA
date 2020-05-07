#include "serv_util.h"

#include<mpi.h>
#include<stdio.h>

void nullAbort(const char* file, int line_no, const char* ptr_name, const char *message){
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	fprintf(stderr, "[%s]: %s (%d in file %s) is NULL in world-rank %d -> aborting program\n",
		       	message, ptr_name, line_no, file, world_rank);
	MPI_Abort(MPI_COMM_WORLD, -1);
	exit(-1);
}
