#ifndef SERVER_H
#define SERVER_H
#include<mpi.h>

// stores information that every rank has after servers are created
struct rank_info{

	// A communicator of exactly one server
	// and 1 or (typically) more simranks
	// associated with this server.
	// The server has rank 0.
	MPI_Comm sim_server_comm;
	// number of simranks (== sim_server_comm_size - 1)
	int n_sim_ranks;
	// my ranknumber in this communicator
	int sim_server_rank;

	// A communicator of every rank
	// of the same type (server or simrank)
	// as this rank
	MPI_Comm my_type_comm;
	// number of servers/simranks
	int my_type_size;
	// my ranknumber in this communicator
	int my_type_rank;

};

// for every observable, hold the delta that indicates how often it should be computed
// 0 means it will never be computed delta>0 means it will be computed every nth timestep,
// starting with timestep 0
struct tasks{
	unsigned int delta_avg;
	unsigned int delta_var;
	unsigned int delta_absmax;
};


// takes a communicator comm_in and creates n_servers server on it
// returns 0 for success, 
int create_servers(MPI_Comm comm_in, int n_servers, struct rank_info *info);


// starts this rank as a server. Will do steps timesteps,
// completing the given tasks at each timestep where it is specified by task
// uses info for communication with simranks
// will fail if the current rank is not a server according to info
int run_server(const struct tasks *task, unsigned int steps, const struct rank_info *info,
		int (*write_out)(const char *,...));


// called by the sim_rank on timesteps where observables are to be calculated
// sends the data to that server that belongs to this rank according to info
// data can be used after return
int send_to_server(int *data, int n_data, struct rank_info *info);


#endif // SERVER_H
