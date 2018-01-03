/* Copyright (C) 2016-2018 Ludwig Schneider
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

//! \file mpiroutines.c
//! \brief Implementation of mpiroutines.h


#include"mpiroutines.h"
#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include <assert.h>
#include "phase.h"
#include "mesh.h"

int init_MPI(struct Phase * p)
{

    int test;
    if (MPI_Comm_size(p->info_MPI.SOMA_MPI_Comm, &(test)) != MPI_SUCCESS) {
	fprintf(stderr, "MPI_ERROR: %d (2)\n", test);
	return -2;
    }
    p->info_MPI.Ncores = test;


    if (MPI_Comm_rank(p->info_MPI.SOMA_MPI_Comm, &(test)) != MPI_SUCCESS) {
	fprintf(stderr, "MPI_ERROR (3) %d \n", test);
	return -3;
    }
    p->info_MPI.current_core = test;

    if (MPI_Barrier(p->info_MPI.SOMA_MPI_Comm) != MPI_SUCCESS) {
	fprintf(stderr, "MPI_ERROR (4)\n");
	return -4;
    }

    p->info_MPI.divergence_sec = 0.;
    p->info_MPI.divergence_counter = 0;

    return 0;
}

int finalize_MPI(void)
{
    return MPI_Finalize();
}

int check_status_on_mpi(const struct Phase*const p,int my_status)
    {
    int*const status_array =(int*const)malloc(p->info_MPI.Ncores*sizeof(int));
    if(status_array == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -1;}
    MPI_Allgather(&my_status,1,MPI_INT,
		  status_array,1,
		  MPI_INT,p->info_MPI.SOMA_MPI_Comm);
    for(int i=0; i < p->info_MPI.Ncores; i++)
	if( status_array[i] != 0)
	    return status_array[i];
    free(status_array);
    return 0;
    }

double mpi_divergence(struct Phase*const p)
    {
    if(p->args.load_balance_arg == 0)
	return 0;
    const double start = MPI_Wtime();
    MPI_Barrier(p->info_MPI.SOMA_MPI_Comm);
    const double end = MPI_Wtime();
    p->info_MPI.divergence_sec += (end-start);
    p->info_MPI.divergence_counter += 1;
    return end-start;
    }

int collective_global_update(struct Phase*const p)
    {
    // Global number of polymers
    MPI_Allreduce(&(p->n_polymers), &(p->n_polymers_global), 1,
		  MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_MPI_Comm);
#pragma acc update device(p->n_polymers_global)
    // Total number of beads
    MPI_Allreduce(&(p->num_all_beads_local), &(p->num_all_beads), 1,
		  MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_MPI_Comm);
#pragma acc update device(p->num_all_beads)

    //Beads per type
    MPI_Allreduce(p->num_bead_type_local , p->num_bead_type ,
		  p->n_types, MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_MPI_Comm);
#pragma acc update device(p->num_bead_type[0:p->n_types])

    update_density_fields(p);
    return 0;
    }

int send_polymer_chain(struct Phase*const p, const uint64_t poly_id, const int destination)
    {
    assert(p);
    if( destination >= p->info_MPI.Ncores)
	{
	fprintf(stderr,"ERROR: %s:%d Rank %d requested to send a polymer"
		" to rank %d, but only %d ranks are known.\n",
		__FILE__,__LINE__, p->info_MPI.current_core, destination, p->info_MPI.Ncores);
	return -2;
	}
    if( p->n_polymers == 0)
	{//Send an empty buffer, to signal that no polymer is going to be send.
	unsigned int zero = 0;
	MPI_Send(&zero, 1 , MPI_UNSIGNED, destination, 0 , p->info_MPI.SOMA_MPI_Comm);
	return 1;
	}
    if( poly_id >= p->n_polymers )
	{
	fprintf(stderr,"ERROR: %s:%d Rank %d: Try to send poly_id %ld, but only %ld polymer local.\n",
		__FILE__,__LINE__, p->info_MPI.current_core,poly_id,p->n_polymers);
	return -3;
	}

    Polymer poly; //Memory to store a polymer.

    pop_polymer(p, poly_id, &poly);

    unsigned int buffer_length = poly_serial_length(p, &poly);
    //Communicate this buffer length with the recieving rank.
    MPI_Send(&buffer_length, 1 , MPI_UNSIGNED, destination, 0 , p->info_MPI.SOMA_MPI_Comm);

    unsigned char*const buffer = (unsigned char*)malloc( buffer_length*sizeof(unsigned char));
    MALLOC_ERROR_CHECK(buffer, buffer_length*sizeof(unsigned char));

    const unsigned int bytes_written = serialize_polymer(p, &poly, buffer);
    if(bytes_written != buffer_length)
	{
	fprintf(stderr,"ERROR: %s:%d Rank %d: failed to serialize polymer %d %d\n",
		__FILE__,__LINE__,p->info_MPI.current_core, bytes_written, buffer_length);
	return -3;
	}

    //After serialization the polymer deep memory can be freed.
    free_polymer(p, &poly);

    //Send the buffer to destination.
    MPI_Send( buffer, buffer_length, MPI_UNSIGNED_CHAR, destination, 0,
	      p->info_MPI.SOMA_MPI_Comm);

    free(buffer);

    return 0;
    }

int recv_polymer_chain(struct Phase*const p, const int source)
    {
    assert(p);
    if( source >= p->info_MPI.Ncores)
	{
	fprintf(stderr,"ERROR: %s:%d Rank %d requested to recv a polymer"
		" from rank %d, but only %d ranks are known.\n",
		__FILE__,__LINE__, p->info_MPI.current_core, source, p->info_MPI.Ncores);
	return -2;
	}

    unsigned int buffer_length;
    //Get buffer length from sending rank.
    MPI_Recv(&buffer_length, 1 , MPI_UNSIGNED , source , 0,
	     p->info_MPI.SOMA_MPI_Comm,MPI_STATUS_IGNORE);

    if( buffer_length > 0) //The sending process can signal "no polymer to send."
	{
	//Allocate the memory to recv
	unsigned char*const buffer = (unsigned char*)malloc( buffer_length * sizeof(unsigned char));
	MALLOC_ERROR_CHECK(buffer, buffer_length*sizeof(unsigned char));

	//Recv the polymer buffer.
	MPI_Recv( buffer, buffer_length, MPI_UNSIGNED_CHAR, source, 0,
		  p->info_MPI.SOMA_MPI_Comm,MPI_STATUS_IGNORE);

	//Deserialize and allocate Memory for recved polymer.
	Polymer poly;
	deserialize_polymer(p, &poly, buffer);

	free(buffer);

	//Push the polymer to the system.
	push_polymer(p, &poly);

	//poly is no longer owner of polymer, so no deacllocation needed.

	return 0;
	}
    return 1;
    }

int send_mult_polymers(struct Phase*const p,const int destination,unsigned int Nsends)
    {
    assert(p);
    if( destination >= p->info_MPI.Ncores)
	{
	fprintf(stderr,"ERROR: %s:%d Rank %d requested to send a polymer"
		" to rank %d, but only %d ranks are known.\n",
		__FILE__,__LINE__, p->info_MPI.current_core, destination, p->info_MPI.Ncores);
	return -2;
	}
    Nsends = Nsends > p->n_polymers ? p->n_polymers : Nsends;
    MPI_Send(&Nsends,1,MPI_UNSIGNED, destination, 0 , p->info_MPI.SOMA_MPI_Comm);

    unsigned int buffer_length = 0;
    for(unsigned int i=0; i < Nsends; i++)
	buffer_length += poly_serial_length(p, p->polymers + p->n_polymers - 1 - i);

    //Communicate this buffer length with the recieving rank.
    MPI_Send(&buffer_length, 1 , MPI_UNSIGNED, destination, 0 , p->info_MPI.SOMA_MPI_Comm);

    if(buffer_length > 0)
	{
	unsigned char*const buffer = (unsigned char*)malloc( buffer_length*sizeof(unsigned char));
	MALLOC_ERROR_CHECK(buffer, buffer_length*sizeof(unsigned char));

	unsigned int bytes_written = 0;

	Polymer poly; //Memory to store a polymer.
	for(unsigned int i=0; i < Nsends; i++)
	    {
	    pop_polymer(p, p->n_polymers - 1, &poly);
	    const unsigned int poly_bytes = serialize_polymer(p, &poly, buffer+bytes_written);
	    //After serialization the polymer deep memory can be freed.
	    free_polymer(p, &poly);
	    bytes_written += poly_bytes;
	    }

	if(bytes_written != buffer_length)
	    {
	    fprintf(stderr,"ERROR: %s:%d Rank %d: failed to serialize polymer %d %d\n",
		    __FILE__,__LINE__,p->info_MPI.current_core, bytes_written, buffer_length);
	    return -3;
	    }

	//Send the buffer to destination.
	MPI_Send( buffer, buffer_length, MPI_UNSIGNED_CHAR, destination, 0,
		  p->info_MPI.SOMA_MPI_Comm);

	free(buffer);
	return Nsends;
	}

    return 0;
    }

int recv_mult_polymers(struct Phase*const p, const int source)
    {
    assert(p);
    if( source >= p->info_MPI.Ncores)
	{
	fprintf(stderr,"ERROR: %s:%d Rank %d requested to recv a polymer"
		" from rank %d, but only %d ranks are known.\n",
		__FILE__,__LINE__, p->info_MPI.current_core, source, p->info_MPI.Ncores);
	return -2;
	}

    unsigned int Nsends;
    MPI_Recv(&Nsends, 1 , MPI_UNSIGNED , source , 0,
	     p->info_MPI.SOMA_MPI_Comm,MPI_STATUS_IGNORE);

    unsigned int buffer_length;
    //Get buffer length from sending rank.
    MPI_Recv(&buffer_length, 1 , MPI_UNSIGNED , source , 0,
	     p->info_MPI.SOMA_MPI_Comm,MPI_STATUS_IGNORE);

    if( buffer_length > 0) //The sending process can signal "no polymer to send."
	{
	//Allocate the memory to recv
	unsigned char*const buffer = (unsigned char*)malloc( buffer_length * sizeof(unsigned char));
	MALLOC_ERROR_CHECK(buffer, buffer_length*sizeof(unsigned char));

	//Recv the polymer buffer.
	MPI_Recv( buffer, buffer_length, MPI_UNSIGNED_CHAR, source, 0,
		  p->info_MPI.SOMA_MPI_Comm,MPI_STATUS_IGNORE);

	//Deserialize and allocate Memory for recved polymer.
	Polymer poly;
	unsigned int bytes_read = 0;
	for(unsigned int i=0 ; i < Nsends; i++)
	    {
	    assert(bytes_read < buffer_length);
	    const unsigned int poly_bytes = deserialize_polymer(p, &poly, buffer+bytes_read);
	    bytes_read += poly_bytes;
	    //Push the polymer to the system.
	    push_polymer(p, &poly);
	    }
	free(buffer);

	//poly is no longer owner of polymer, so no deacllocation needed.

	return Nsends;
	}

    return 0;
    }

int load_balance_mpi_ranks(struct Phase*const p)
    {

    if( p->info_MPI.Ncores <= 1)
	return 0;
    double*const waiting_time = (double*const)malloc(p->info_MPI.Ncores*sizeof(double));
    if(waiting_time == NULL){fprintf(stderr,"ERROR: Malloc %s:%d\n",__FILE__,__LINE__); return -1;}
    double div = p->info_MPI.divergence_sec/p->info_MPI.divergence_counter;

    MPI_Allgather( &(div), 1 , MPI_DOUBLE, waiting_time, 1 , MPI_DOUBLE,
		   p->info_MPI.SOMA_MPI_Comm);
    //Reset counters.
    p->info_MPI.divergence_sec = 0;
    p->info_MPI.divergence_counter = 0;

    int arg_max=0,arg_min=0;
    for(int i=0; i < p->info_MPI.Ncores; i++)
	{
	if( waiting_time[i] > waiting_time[arg_max] )
	    arg_max = i;
	if( waiting_time[i] < waiting_time[arg_min] )
	    arg_min = i;
	}

    //Quick exit, if no balance required.
    const double divergences = waiting_time[arg_max]-waiting_time[arg_min];
    free(waiting_time);
    const double seconds_per_step = p->tps_elapsed_time / p->tps_elapsed_steps;
    double p_waiting =  divergences/seconds_per_step;

    //Last tps could be different on different rank, so avoid hangs.
    MPI_Bcast(&p_waiting,1,MPI_DOUBLE,0,p->info_MPI.SOMA_MPI_Comm);

    const double passed_acc = p->args.accepted_load_inbalance_arg /100.;
    const double max_p_waiting = passed_acc >=0 && passed_acc<= 1 ? passed_acc : 0.08;
    char quick_exit = p_waiting < max_p_waiting;
    if( quick_exit )
	return 0;

    unsigned int Nchains = p_waiting*p->n_polymers*0.05 + 1;
    //Limit transferred chains to 5%
    if( Nchains > 0.05 * p->n_polymers )
	Nchains = 0.05 * p->n_polymers +1;

    unsigned int Nsend = UINT_MAX;
    //! \todo Optimize chain to send. (Poly-Type based?)
    if( p->info_MPI.current_core == arg_min )
	Nsend = send_mult_polymers(p, arg_max, Nchains);

    if( p->info_MPI.current_core == arg_max )
	Nsend = recv_mult_polymers(p, arg_min);

    if( arg_min != 0)
	{
	if(p->info_MPI.current_core == arg_min)
	    {
	    MPI_Send(&Nsend, 1 , MPI_UNSIGNED, 0, 0 , p->info_MPI.SOMA_MPI_Comm);
	    }
	if(p->info_MPI.current_core == 0)
	    {
	    MPI_Recv(&Nsend, 1 , MPI_UNSIGNED , arg_min , 0,
		     p->info_MPI.SOMA_MPI_Comm,MPI_STATUS_IGNORE);
	    }
	}

    if( p->info_MPI.current_core == 0)
	printf("INFO: Load balance @t=%d, sending %d chains from rank"
	       " %d to rank %d, because of %f percent waiting.\n",
	       p->time,Nsend,arg_min,arg_max,p_waiting*100.);

    //Synchronize rank, because otherwise a false unbalance would be
    //detected.
    MPI_Barrier(p->info_MPI.SOMA_MPI_Comm);
    return 1;
    }
