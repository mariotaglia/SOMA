#include "server.h"
#include "serv_util.h"

#include <mpi.h>
#include <stdlib.h>
#include <stdbool.h>

bool need_to_do(int delta, unsigned int time){

	return (delta != 0 && time%delta == 0);

}

int max(int a, int b){

	if (a>b){
		return a;
	}
	return b;
}

int abs(int a){

	if (a > 0){
		return a;
	}

	return -a;
}




int create_servers(MPI_Comm comm_in, int n_servers, struct rank_info *info){

	int rank, size;
	MPI_Comm_rank(comm_in, &rank);
	MPI_Comm_size(comm_in, &size);
	if ( 2*n_servers > size ){
		return -1;
	}

	int is_server = (rank < n_servers);

	MPI_Comm_split(comm_in, is_server, 0, &(info->my_type_comm));
	MPI_Comm_rank(info->my_type_comm, &(info->my_type_rank));
	MPI_Comm_size(info->my_type_comm, &(info->my_type_size));
	
	MPI_Comm_split(comm_in, (info->my_type_rank)%n_servers, !is_server,
		&(info->sim_server_comm));
	MPI_Comm_rank(info->sim_server_comm, &(info->sim_server_rank));

	int sim_server_size;
	MPI_Comm_size(info->sim_server_comm, &sim_server_size);
	info->n_sim_ranks = sim_server_size - 1;

	return 0;
}

bool nothing_to_do(const struct tasks *task, unsigned int time){

	if (need_to_do(task->delta_avg, time)){
		return false;
	}

	if (need_to_do(task->delta_var, time)){
		return false;
	}

	if (need_to_do(task->delta_absmax, time)){
		return false;
	}

	return true;
}

struct receiver{

	//housekeeping to receive the data
	//only to be used by receive_next
	MPI_Request size_req;
	int *size_on_rank;
	int *dspls;
	MPI_Comm comm;
	int commsize;


	//the actual data received
	int size;
	int *data;

};

void init_receiver(struct receiver * r, const struct rank_info *info){

	r->data = NULL;

	r->size_on_rank = (int *)malloc(sizeof(int) * (info->n_sim_ranks + 1));
	ABORT_ON_NULL(r->size_on_rank, "MALLOC ERROR");
	r->dspls = (int *)malloc(sizeof(int) * (info->n_sim_ranks + 1));
	ABORT_ON_NULL(r->dspls, "MALLOC ERROR");

	r->comm = info->sim_server_comm;
	MPI_Comm_size(r->comm, &r->commsize);

}

void free_receiver(struct receiver * r){
	
	free(r->size_on_rank);
	free(r->dspls);
	free(r->data);
}

// blocks until the data for the next timestep is in r.data and ready to be processed
// returns the size of the data in the size-field of the receiver
// requires that the last call on the receiver was either init_receiver or receive_next
void receive_next(struct  receiver* r){

	// Gather the sizes that each rank holds (to set up the first gatherv for the data)
	int my_size = 0; 
	MPI_Igather(&my_size, 1, MPI_INT,
			r->size_on_rank, 1, MPI_INT,
			0, r->comm,
			&(r->size_req));
	MPI_Wait(&(r->size_req), MPI_STATUS_IGNORE);

	r->size = 0;
	for (int i=0; i < r->commsize; i++){
		r->dspls[i] = r->size;
		r->size += r->size_on_rank[i];
	}

	r->data = (int *)realloc(r->data, r->size * sizeof(int));	

	ABORT_ON_NULL(r->data, "REALLOC ERROR");

	MPI_Request data_req;

	// Gatherv the actual data
	MPI_Igatherv(NULL, 0, MPI_INT,
			r->data, r->size_on_rank, r->dspls, MPI_INT,
			0, r->comm,
			&data_req);

	MPI_Wait(&data_req, MPI_STATUS_IGNORE);
}

int run_server(const struct tasks *task, unsigned int steps, const struct rank_info *info,
		int (*write_out)(const char *,...)){

	// must be started on a server
	if (info->sim_server_rank != 0){
		return -1;
	}

	struct receiver recv;
	init_receiver(&recv, info);

	for (unsigned int t=0; t<steps; t++){

		if (nothing_to_do(task, t)){
			continue;
		}

		receive_next(&recv);

		int size = recv.size;
		int *data = recv.data;

		// perform analysis & write to output

		double avg = 0;
		double var = 0;
		double absmax = 0;
		for (long i=0; i<size; i++){
			avg += (double)data[i];
			var +=data[i]*data[i];
			absmax = max(absmax, abs(data[i]));
		}
		// get global values of the observables all servers
		double gl_avg;
		double gl_var;
		double gl_absmax;
		int gl_total;

		MPI_Reduce(&size, &gl_total, 1, MPI_INT,
				MPI_SUM, 0, info->my_type_comm);
		MPI_Reduce(&avg, &gl_avg, 1, MPI_DOUBLE,
				MPI_SUM, 0, info->my_type_comm);
		MPI_Reduce(&var, &gl_var, 1, MPI_DOUBLE,
				MPI_SUM, 0, info->my_type_comm);
		MPI_Reduce(&absmax, &gl_absmax, 1, MPI_DOUBLE,
				MPI_MAX, 0, info->my_type_comm);

		if ( info->my_type_rank == 0){

			gl_var /= gl_total;
			gl_avg /= gl_total;
			gl_var -= gl_avg*gl_avg;

			write_out("%d: ", t, avg);
			if (need_to_do(task->delta_avg, t)){
				write_out("mean=%lf", gl_avg);
			}
			if (need_to_do(task->delta_var, t)){
				write_out("var=%lf", gl_var);
			}
			if (need_to_do(task->delta_absmax, t)){
				write_out("absmax=%lf", gl_absmax);

			}	

			write_out("\n");
		}
		
	}

	free_receiver(&recv);
	return 0;

}
