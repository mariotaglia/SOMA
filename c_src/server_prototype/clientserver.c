#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<assert.h>

#include<pthread.h>
#include<stdatomic.h>

#include<limits.h> // for message tags
#include<unistd.h> // for usleep

// messages sent to the server with this tag will be answered by the time-service thread
// with a message that includes the current time step of the server.
#define TIME_TAG 200000000

// messages sent to the server with this tag will be answered with the server-timestep as well,
// but only after the server-timestep has been incremented
#define WAKE_UP_TAG 2000000001


// number of particles per simulation rank
#define N_PARTICLES 5

// send client and server to sleep in each timestep for the specified number of microseconds
#define CLIENT_TIMEOUT 0
#define SERVER_TIMEOUT 100

//information that every rank has (server or client)
struct process_info{

	int world_rank;

	MPI_Comm client_server_comm;
	int client_server_rank;
	int client_num;

	MPI_Comm my_type_comm;
	int my_type_rank;
	int my_type_size;

	bool is_server;

};

// argument to be passed to time_service, which is the start routing of the time-service-thread
struct time_service_args{
	_Atomic(unsigned int) * time_ptr;
	MPI_Comm client_server_comm;
};

// data needed to manage the wait_server function of the client
struct client_blocker{
	unsigned int tserver;
	MPI_Request req;
	MPI_Comm comm;
	int calc_ahead;
	bool has_requested;
};




// splits a communicator old_comm into clients (calculation ranks) and num_serv server ranks.
// returns: 
// in the my_type-communicator: all the processes of the same type (client or server) as the calling process.
// in the client_server-communicator: a communicator with exactly one server that has rank 0, the other ranks are calculation ranks
// assigned to this server.
int create_servers(MPI_Comm old_comm, int num_serv, MPI_Comm *my_type_comm, MPI_Comm *client_server_comm);
void print_info(const struct process_info * info);
void init_info(struct process_info * info, MPI_Comm my_type_comm, MPI_Comm client_server_comm);
void run_server(struct process_info * info, int steps, int calc_ahead);
void run_client(struct process_info * info, int steps, int calc_ahead);

int main(int argc, char **argv){

	// ensure thread-safety of MPI
	int t_safety;
	int init_err = MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &t_safety);
	if (init_err != MPI_SUCCESS || t_safety != MPI_THREAD_MULTIPLE){
		printf("Error initializing MPI, quitting\n");
		MPI_Finalize();
		return -1;
	}

	if ( argc != 2 ){
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		if (rank==0){
			printf("need one argument, number of servers, got %d\n", argc);
		}
		MPI_Finalize();
		return -1;
	}

	int num_serv = atoi(argv[1]);

	
	MPI_Comm my_type_comm, client_server_comm;
	int err = create_servers(MPI_COMM_WORLD, num_serv, &my_type_comm, &client_server_comm);
	if (err){
		printf("error creating servers\n");
		MPI_Finalize();
		return -1;
	}

	struct process_info info;
	init_info(&info, my_type_comm, client_server_comm);

	print_info(&info);

	
	int n = 1000; // number of timesteps
	int calc_ahead = 5; // how many timesteps the clients are allowed to be ahead of the server

	if (info.is_server){
		run_server(&info, n, calc_ahead);
	}
	else{
		run_client(&info, n, calc_ahead);
	}

	printf("simulation ran to the end for rank %d\n", info.client_server_rank);
	
	MPI_Finalize();

	return 0;

}

void *time_service(void * time_service_args){

	struct time_service_args * tsa = time_service_args;
	_Atomic(unsigned int) * time_ptr = tsa->time_ptr;
	MPI_Comm comm = tsa->client_server_comm;

	while (true) {

		MPI_Status status;
		printf("timeservice attempting to receive\n");fflush(stdout);
		MPI_Recv(MPI_BOTTOM, 0, MPI_INT, MPI_ANY_SOURCE, TIME_TAG, comm, &status); 
		printf("timeservice successfully received a message from rank %d\n", status.MPI_SOURCE);fflush(stdout);

		// the server sends this message to quit the thread
		if (status.MPI_SOURCE == 0){
			return 0;
		}
		// other ranks send this message to receive the time
		else {
			unsigned int t = atomic_load(time_ptr);
			printf("timeservice will now message rank %d\n", status.MPI_SOURCE);
			MPI_Send(&t, 1, MPI_UNSIGNED, status.MPI_SOURCE, TIME_TAG, comm);
			printf("timeservice-message to rank %d was sent\n", status.MPI_SOURCE);
		}
	
	}		
}

void run_server(struct process_info * proc_info, int steps, int calc_ahead){


	_Atomic(unsigned int) time = ATOMIC_VAR_INIT(0);
	pthread_t time_service_thread;

	FILE *out_file = fopen("out.dat", "w");
	
	//start the time service
	struct time_service_args tsa;
	tsa.time_ptr = &(time);
	tsa.client_server_comm = proc_info->client_server_comm;
	if (pthread_create(&(time_service_thread), NULL, time_service, &tsa)){
		printf("Error creating thread\n");
	}

	// allocate memory for 1 data-message per timestep and client,
	// for all the clients and for calc_ahead timesteps
	int (*pos)[proc_info->client_num][N_PARTICLES] = malloc(calc_ahead*proc_info->client_num*N_PARTICLES*sizeof(int));
	if (pos==NULL){
		printf("malloc error\n");
	}
	// for every data-message that can be received, have one open receive:
	MPI_Request (*reqs)[proc_info->client_num] = malloc(calc_ahead*proc_info->client_num*sizeof(MPI_Request));
	if (reqs==NULL){
		printf("malloc error\n");
	}
	for (int t=0; t<calc_ahead; t++){
		for (int client=0; client<proc_info->client_num; client++){
			printf("data_server: server creating Irecv with (source,tag)=(%d,%d)\n",client+1,t%calc_ahead);
			MPI_Irecv(pos[t][client],N_PARTICLES, MPI_INT,
					client+1, t%calc_ahead, proc_info->client_server_comm,
					&reqs[t][client]);
		}
	}

	// for every client, have one receive in case they get blocked by being too far ahead of the server
	MPI_Request *block_req= malloc(proc_info->client_num*sizeof(MPI_Request));
	if (block_req==NULL){
		printf("malloc error\n");
	}
	unsigned int *tneeded = malloc(proc_info->client_num*sizeof(unsigned int));
	if (tneeded==NULL){
		printf("malloc error\n");
	}
	for (int client=0; client < proc_info->client_num; client++){
		MPI_Irecv(&tneeded[client], 1, MPI_UNSIGNED,
				client+1, WAKE_UP_TAG, proc_info->client_server_comm,
				&block_req[client]);
	}

	// run server-loop: receive data, calculate and do file-input
	while(time < steps){

		printf("data_server: server is on a new timestep, proceeding to wait for all data for timestep %u (from %d clients)\n", time, proc_info->client_num);
		MPI_Waitall(proc_info->client_num, reqs[time%calc_ahead],
				MPI_STATUSES_IGNORE);
		printf("data_server: all data arrived on server for timestep %d\n", time);fflush(stdout);
		// pos now contains all the data for this timestep

		// do calculation, write to file
		double var = 0;
		double mean = 0;
		for (int c=0; c<proc_info->client_num; c++){
			for (int p=0; p<N_PARTICLES; p++){
				int x = pos[time%calc_ahead][c][p];
				mean += x;
				var += x*x;
			}
		}
		int particle_num = proc_info->client_num *N_PARTICLES;
		var /= particle_num;
		mean /= particle_num;
		var -= mean*mean;
		fprintf(out_file, "t=%u\t|mean=%lf\t|var=%lf\n", time, mean, var);
		fflush(out_file);

		// artificial delay, simulating a long analysis
		usleep(SERVER_TIMEOUT);

		//recreate receives
		for (int client=0; client < proc_info->client_num; client++){

			MPI_Irecv(pos[time%calc_ahead][client],N_PARTICLES, MPI_INT, client+1,
				       	time%calc_ahead, proc_info->client_server_comm,
					&reqs[time%calc_ahead][client]);
		}

		//increment timestep atomically
		atomic_fetch_add(&time, 1);


		//collect distressed clients, if there are any, and create new Irecvs for them
		int blocked_num;
		int index_array[proc_info->client_num];
		MPI_Status statuses[proc_info->client_num];
		MPI_Testsome(proc_info->client_num, block_req,
				&blocked_num, index_array,
				statuses);
		if (blocked_num==0){
			printf("data_server: no clients found struggling at time %d\n", time);
		}
		else{
			printf("data_server: %d clients are struggling at time %d\n", blocked_num, time);
		}
		for (int i=0; i < blocked_num; i++){
			int ind = index_array[i];
			printf("data_server: client (ind+1:%d,source:%d) found to be struggling, asking for time %u, and\
					server time is %u\n", ind+1, statuses[ind].MPI_SOURCE, tneeded[ind],
					time); fflush(stdout);
			assert(tneeded[ind] <= time);
			MPI_Irecv(&tneeded[ind],1,MPI_UNSIGNED,
					ind+1, WAKE_UP_TAG,
					proc_info->client_server_comm,
					&block_req[ind]);

			printf("data_server: server is waking up client %d\n", ind+1);
			MPI_Send(&time, 1, MPI_UNSIGNED,
					ind+1, WAKE_UP_TAG,
					proc_info->client_server_comm);
		}

	}


	// cleanup server
	MPI_Send(MPI_BOTTOM, 0, MPI_INT, 0, TIME_TAG, proc_info->client_server_comm);
	if (pthread_join(time_service_thread, NULL)){
		printf("Error joining time-service-thread");
	}
	free(pos);
	free(reqs);
	free(tneeded);
	free(block_req);
	fclose(out_file);
	
}

void init_client_blocker(struct client_blocker * blocker, int calc_ahead, MPI_Comm client_server_comm){

	// leaves request unitialized, as it will be initialized by Send-functions
	blocker->tserver = 0;
	blocker->has_requested = false;
	blocker->calc_ahead = calc_ahead;
	blocker->comm = client_server_comm;

}

// return when the server is ready to receive the data of timestep tclient
void wait_server(struct client_blocker * blocker, unsigned int tclient){

	if (blocker->has_requested){
		//complete the request to get the timestep of the server
		blocker->has_requested = false;
		MPI_Wait(&(blocker->req), MPI_STATUS_IGNORE);
		//blocker->tserver is now usable again
	}

	if (blocker->tserver + blocker->calc_ahead - 1 > tclient){
		// it is safe to just proceed
		return;
	}
	else if (blocker->tserver + blocker->calc_ahead -1 == tclient){
		// can proceed, but ask the server for a new timestep
		MPI_Ssend(MPI_BOTTOM, 0, MPI_INT, 0,
				TIME_TAG, blocker->comm);
		MPI_Irecv(&(blocker->tserver), 1, MPI_UNSIGNED, 0,
				TIME_TAG, blocker->comm,
				&(blocker->req));
		blocker->has_requested = true;
		return;
	}
	else {
		// block client, because server isn't ready: this is the case that should be avoided.
		//calculate rank just for printing
		int rank;
		MPI_Comm_rank(blocker->comm, &rank);
		printf("client %d is blocking to wait for server at clienttime %u\n", rank, tclient);fflush(stdout);

		// send special message to the server, and block until it is answered.
		unsigned int tneeded = tclient+1-blocker->calc_ahead;
		MPI_Ssend(&tneeded, 1, MPI_UNSIGNED, 0, WAKE_UP_TAG, blocker->comm);
		printf("client %d;tclient %u hangup: send completed\n", rank, tclient);fflush(stdout);
		MPI_Recv(&(blocker->tserver),1, MPI_UNSIGNED, 0,
				WAKE_UP_TAG, blocker->comm, MPI_STATUS_IGNORE);
		printf("client %d;tclient %u hangup: recv completed\n", rank, tclient);fflush(stdout);
		return;
	}
	


}

void run_client(struct process_info * info, int steps, int calc_ahead){

	printf("Hello from client\n");

	struct client_blocker blocker;
	init_client_blocker(&blocker, calc_ahead, info->client_server_comm);

	// initial configuration: N_PARTICLESrandom walkers at position 0
	int pos[N_PARTICLES] = {0};
	srand((unsigned)info->my_type_rank);

	for (unsigned int t=0; t<steps; t++){

		// wait until server is ready to receive data (usually: don't wait at all)
		printf("client %d initiating wait for timestep %d\n", info->client_server_rank, t); fflush(stdout);
		wait_server(&blocker, t);
		printf("client %d finished wait\n", info->client_server_rank); fflush(stdout);

		// start data transfer to server
		MPI_Request req;
		MPI_Issend(pos, N_PARTICLES, MPI_INT, 0, t%calc_ahead,
			       	info->client_server_comm, &req);
		printf("client %d initiated data transfer, expects tag: %d\n", info->client_server_rank, t%calc_ahead); fflush(stdout);

		// simulate: change position of every walker randomly by +-1
		for (int i=0; i<N_PARTICLES; i++){
			int d = 2*(rand()%2) - 1;
			pos[i] += d;
		}

		//artificially increase duration of simulation
		usleep(CLIENT_TIMEOUT);


		// finalize data transfer
		MPI_Wait(&req, MPI_STATUS_IGNORE);
		printf("client %d has finalized data transfer for step %d\n", info->client_server_rank, t);
	
	}

	printf("client %d is exiting his run\n", info->client_server_rank);

}


void simple_gather(const struct process_info * info){

	if (info->is_server){
		printf("server has client_server_rank %d\n", info->client_server_rank);
		int msg;
		for (int i=0; i<info->client_num; i++){
			MPI_Status status;
			MPI_Recv(&msg, 1, MPI_INT, MPI_ANY_SOURCE, 0, info->client_server_comm, &status);
			printf("server %d has received message %d from rank %d\n", info->my_type_rank, msg, status.MPI_SOURCE);
		}
	}
	else{
		int msg = info->client_server_rank;
		MPI_Send(&msg, 1, MPI_INT, 0, 0, info->client_server_comm);

	}
	printf("....\n");

}

void init_info(struct process_info * info, MPI_Comm my_type_comm, MPI_Comm client_server_comm){

	MPI_Comm_rank(MPI_COMM_WORLD, &(info->world_rank));
	MPI_Comm_rank(my_type_comm, &(info->my_type_rank));
	MPI_Comm_rank(client_server_comm, &(info->client_server_rank));
	info->is_server = (info->client_server_rank == 0);
	MPI_Comm_size(my_type_comm, &(info->my_type_size));
	MPI_Comm_size(client_server_comm, &(info->client_num));
	info->client_num -= 1; // do not count server
	info->client_server_comm = client_server_comm;
	info->my_type_comm = my_type_comm;

}

void print_info(const struct process_info * info){

	char * role = info->is_server? "server" : "client";
	printf("I am world rank %d, I am %s %d out of %d, my rank in the cs is %d with %d clients\n ===\n",
			info->world_rank, role, info->my_type_rank, info->my_type_size,
			info->client_server_rank, info->client_num);

}
	


int create_servers(MPI_Comm old_comm, int num_serv, MPI_Comm *my_type_comm, MPI_Comm *client_server_comm){

	int rank, size;
	MPI_Comm_rank(old_comm, &rank);
	MPI_Comm_size(old_comm, &size);
	
	if ( 2*num_serv > size ){
		fprintf(stderr, "too many servers for so few computation ranks");
		return -1;
	}

	int is_server = (rank < num_serv);

	MPI_Comm_split(old_comm, is_server, 0, my_type_comm);
	
	int my_type_rank;
	MPI_Comm_rank(*my_type_comm, &my_type_rank);

	MPI_Comm_split(old_comm, my_type_rank%num_serv, !is_server, client_server_comm);

	//MPI_Comm_set_errhandler(*my_type_comm, MPI_ERRORS_RETURN);
	//MPI_Comm_set_errhandler(*client_server_comm, MPI_ERRORS_RETURN);

	return 0;
}


