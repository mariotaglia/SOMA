/* Copyright (C) 2016-2021 Ludwig Schneider
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
#if ( ENABLE_MPI == 1 )
#include<mpi.h>
#endif                          //ENABLE_MPI
#include<stdio.h>
#include<stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdbool.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include "phase.h"
#include "mesh.h"
#include "polymer.h"
#include "monomer.h"

int init_MPI(struct Phase *p)
{
    //Initialize all variable even without MPI
    p->info_MPI.domain_divergence_sec = 0.;
    p->info_MPI.domain_divergence_counter = 0;

    p->info_MPI.world_rank = 0;
    p->info_MPI.world_size = 1;
    p->info_MPI.sim_size = 1;
    p->info_MPI.sim_rank = 0;
    p->info_MPI.domain_size = 1;
    p->info_MPI.domain_rank = 0;

#if ( ENABLE_MPI == 1 )
    int test;
    if (MPI_Comm_rank(p->info_MPI.SOMA_comm_world, &(test)) != MPI_SUCCESS)
        {
            fprintf(stderr, "MPI_ERROR (5) %d \n", p->info_MPI.world_rank);
            return -5;
        }
    p->info_MPI.world_rank = test;
    if (MPI_Comm_size(p->info_MPI.SOMA_comm_world, &(test)) != MPI_SUCCESS)
        {
            fprintf(stderr, "MPI_ERROR (6) %d \n", p->info_MPI.world_size);
            return -6;
        }
    p->info_MPI.world_size = test;

    if (MPI_Comm_size(p->info_MPI.SOMA_comm_sim, &(test)) != MPI_SUCCESS)
        {
            fprintf(stderr, "MPI_ERROR: %d (2)\n", test);
            return -2;
        }
    p->info_MPI.sim_size = test;

    if (MPI_Comm_rank(p->info_MPI.SOMA_comm_sim, &(test)) != MPI_SUCCESS)
        {
            fprintf(stderr, "MPI_ERROR (3) %d \n", test);
            return -3;
        }
    p->info_MPI.sim_rank = test;

    if (p->info_MPI.sim_size < p->args.N_domains_arg)
        {
            fprintf(stderr, "ERROR: Requested more domains %d than ranks are available %d. %s:%d\n",
                    p->args.N_domains_arg, p->info_MPI.sim_size, __FILE__, __LINE__);
            return -9;
        }
    if (p->info_MPI.sim_size % p->args.N_domains_arg != 0)
        {
            fprintf(stderr, "ERROR: %s:%d invalid number of domains. #ranks mod Ndomains != 0.\n", __FILE__, __LINE__);
            return -7;
        }
    const unsigned int domain_size = p->info_MPI.sim_size / p->args.N_domains_arg;
    const int domain_color = p->info_MPI.sim_rank / domain_size;
    if (MPI_Comm_split(p->info_MPI.SOMA_comm_sim, domain_color, 0, &(p->info_MPI.SOMA_comm_domain)) != MPI_SUCCESS)
        {
            fprintf(stderr, " MPI_ERROR: %s:%d ", __FILE__, __LINE__);
            return -8;
        }
    if (MPI_Comm_size(p->info_MPI.SOMA_comm_domain, &(test)) != MPI_SUCCESS)
        {
            fprintf(stderr, "MPI_ERROR: %d (2)\n", test);
            return -2;
        }
    p->info_MPI.domain_size = test;
    assert((int)domain_size == p->info_MPI.domain_size);

    if (MPI_Comm_rank(p->info_MPI.SOMA_comm_domain, &(test)) != MPI_SUCCESS)
        {
            fprintf(stderr, "MPI_ERROR (3) %d \n", test);
            return -3;
        }
    p->info_MPI.domain_rank = test;

    uint32_t fixed_seed;
    if (!p->args.rng_seed_given || p->args.rng_seed_arg < 0)
        {
            uint32_t time_bytes = time(NULL);
            fixed_seed = time_bytes;
#if __GLIBC__ > 2 || __GLIBC_MINOR__ > 24
#include <sys/random.h>
            uint32_t entropy_bytes;
            if (getentropy(&entropy_bytes, sizeof(uint32_t)) == 0)
                fixed_seed ^= entropy_bytes;
#endif                          //__GLIBC

        }
    else
        fixed_seed = p->args.rng_seed_arg;
    MPI_Bcast(&fixed_seed, 1, MPI_UINT32_T, 0, p->info_MPI.SOMA_comm_sim);
    p->args.rng_seed_arg = fixed_seed;
#endif                          //ENABLE_MPI
    return 0;
}

int finalize_MPI(struct Info_MPI *mpi)
{
#if ( ENABLE_MPI == 1 )
    if (mpi->SOMA_comm_domain != MPI_COMM_NULL)
        MPI_Comm_free(&(mpi->SOMA_comm_domain));
    if (mpi->SOMA_comm_sim != MPI_COMM_NULL)
        MPI_Comm_free(&(mpi->SOMA_comm_sim));

    if (mpi->SOMA_comm_world != MPI_COMM_NULL)
        MPI_Comm_free(&(mpi->SOMA_comm_world));

#ifdef ENABLE_NCCL
    if (mpi->gpu_id >= 0)
        {
            ncclCommDestroy(mpi->SOMA_nccl_world);
            ncclCommDestroy(mpi->SOMA_nccl_sim);
            ncclCommDestroy(mpi->SOMA_nccl_domain);
        }
#endif                          //ENABLE_MPI_CUDA

    return MPI_Finalize();
#else
    return 0;
#endif                          //ENABLE_MPI
}

int check_status_on_mpi(const struct Phase *const p, int my_status)
{
#if ( ENABLE_MPI == 1 )
    int *const status_array = (int *const)malloc(p->info_MPI.sim_size * sizeof(int));
    if (status_array == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d \n", __FILE__, __LINE__);
            return -1;
        }
    MPI_Allgather(&my_status, 1, MPI_INT, status_array, 1, MPI_INT, p->info_MPI.SOMA_comm_sim);
    for (int i = 0; i < p->info_MPI.sim_size; i++)
        if (status_array[i] != 0)
            return status_array[i];
    free(status_array);
    return 0;
#else
    return my_status;
#endif                          //ENABLE_MPI
}

#if ( ENABLE_MPI == 1 )
double mpi_divergence(struct Phase *const p)
{
    if (p->args.load_balance_arg == 0)
        return 0;
    struct timeval start_tv, end_tv;
    gettimeofday(&start_tv, NULL);
    const double start = start_tv.tv_sec + start_tv.tv_usec * 1e-6;
    MPI_Barrier(p->info_MPI.SOMA_comm_domain);
    gettimeofday(&end_tv, NULL);
    const double end = end_tv.tv_sec + end_tv.tv_usec * 1e-6;
    p->info_MPI.domain_divergence_sec += (end - start);
    p->info_MPI.domain_divergence_counter += 1;
    return end - start;
}

int collective_global_update(struct Phase *const p)
{

    // Global number of polymers
    MPI_Allreduce(&(p->n_polymers), &(p->n_polymers_global), 1, MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);
#pragma acc update device(p->n_polymers_global)
    // Total number of beads
    MPI_Allreduce(&(p->num_all_beads_local), &(p->num_all_beads), 1, MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);
#pragma acc update device(p->num_all_beads)

    update_density_fields(p);
    return 0;
}

int send_polymer_chain(struct Phase *const p, const uint64_t poly_id, const int destination, const MPI_Comm comm)
{
    assert(p);
    if (destination >= p->info_MPI.world_size)
        {
            fprintf(stderr, "ERROR: %s:%d World rank %d requested to send a polymer"
                    " to world rank %d, but only world %d ranks are known.\n",
                    __FILE__, __LINE__, p->info_MPI.world_rank, destination, p->info_MPI.world_rank);
            return -2;
        }
    if (p->n_polymers == 0)
        {                       //Send an empty buffer, to signal that no polymer is going to be send.
            unsigned int zero = 0;
            MPI_Send(&zero, 1, MPI_UNSIGNED, destination, 0, comm);
            return 1;
        }
    if (poly_id >= p->n_polymers)
        {
            fprintf(stderr, "ERROR: %s:%d World rank %d: Try to send poly_id %ld, but only %ld polymer local.\n",
                    __FILE__, __LINE__, p->info_MPI.world_rank, poly_id, p->n_polymers);
            return -3;
        }

    Polymer poly;               //Memory to store a polymer.

    pop_polymer(p, poly_id, &poly);

    unsigned int buffer_length = poly_serial_length(p, &poly);
    //Communicate this buffer length with the recieving rank.
    MPI_Send(&buffer_length, 1, MPI_UNSIGNED, destination, 0, comm);

    unsigned char *const buffer = (unsigned char *)malloc(buffer_length * sizeof(unsigned char));
    MALLOC_ERROR_CHECK(buffer, buffer_length * sizeof(unsigned char));

    const unsigned int bytes_written = serialize_polymer(p, &poly, buffer);
    if (bytes_written != buffer_length)
        {
            fprintf(stderr, "ERROR: %s:%d world Rank %d: failed to serialize polymer %d %d\n",
                    __FILE__, __LINE__, p->info_MPI.world_rank, bytes_written, buffer_length);
            return -3;
        }

    //Send the buffer to destination.
    MPI_Send(buffer, buffer_length, MPI_UNSIGNED_CHAR, destination, 0, comm);

    free(buffer);

    return 0;
}

int recv_polymer_chain(struct Phase *const p, const int source, const MPI_Comm comm)
{
    assert(p);
    if (source >= p->info_MPI.world_rank)
        {
            fprintf(stderr, "ERROR: %s:%d World rank %d requested to recv a polymer"
                    " from w. rank %d, but only %d w. ranks are known.\n",
                    __FILE__, __LINE__, p->info_MPI.world_rank, source, p->info_MPI.world_size);
            return -2;
        }

    unsigned int buffer_length;
    //Get buffer length from sending rank.
    MPI_Recv(&buffer_length, 1, MPI_UNSIGNED, source, 0, comm, MPI_STATUS_IGNORE);

    if (buffer_length > 0)      //The sending process can signal "no polymer to send."
        {
            //Allocate the memory to recv
            unsigned char *const buffer = (unsigned char *)malloc(buffer_length * sizeof(unsigned char));
            MALLOC_ERROR_CHECK(buffer, buffer_length * sizeof(unsigned char));

            //Recv the polymer buffer.
            MPI_Recv(buffer, buffer_length, MPI_UNSIGNED_CHAR, source, 0, comm, MPI_STATUS_IGNORE);

            //Deserialize and allocate Memory for recved polymer.
            Polymer poly;
            int err = deserialize_polymer(p, &poly, buffer);
            if (err != (int)buffer_length)
                {
                    fprintf(stderr, "ERROR: %s:%d World rank %d: invalid deserialization %d %d.\n",
                            __FILE__, __LINE__, p->info_MPI.world_rank, err, buffer_length);
                    return -2;
                }

            free(buffer);

            //Push the polymer to the system.
            push_polymer(p, &poly);

            //poly is no longer owner of polymer, so no deacllocation needed.

            return 0;
        }
    return 1;
}

int send_mult_polymers(struct Phase *const p, const int destination, unsigned int Nsends, const MPI_Comm comm)
{
    assert(p);
    if (destination >= p->info_MPI.world_size)
        {
            fprintf(stderr, "ERROR: %s:%d World rank %d requested to send a polymer"
                    " to w. rank %d, but only %d w. ranks are known.\n",
                    __FILE__, __LINE__, p->info_MPI.world_rank, destination, p->info_MPI.world_size);
            return -2;
        }
    Nsends = Nsends > p->n_polymers ? p->n_polymers : Nsends;
    MPI_Send(&Nsends, 1, MPI_UNSIGNED, destination, 0, comm);

    unsigned int buffer_length = 0;
    for (unsigned int i = 0; i < Nsends; i++)
        buffer_length += poly_serial_length(p, p->polymers + p->n_polymers - 1 - i);

    //Communicate this buffer length with the recieving rank.
    MPI_Send(&buffer_length, 1, MPI_UNSIGNED, destination, 1, comm);

    if (buffer_length > 0)
        {
            unsigned char *const buffer = (unsigned char *)malloc(buffer_length * sizeof(unsigned char));
            MALLOC_ERROR_CHECK(buffer, buffer_length * sizeof(unsigned char));

            unsigned int bytes_written = 0;

            Polymer poly;       //Memory to store a polymer.
            for (unsigned int i = 0; i < Nsends; i++)
                {
                    pop_polymer(p, p->n_polymers - 1, &poly);
                    const unsigned int poly_bytes = serialize_polymer(p, &poly, buffer + bytes_written);
                    if (poly_bytes != poly_serial_length(p, &poly))
                        {
                            fprintf(stderr, "ERROR: %s:%d:%d constructing invalid buffer construction %d,%d.\n",
                                    __FILE__, __LINE__, p->info_MPI.world_size,
                                    poly_bytes, poly_serial_length(p, &poly));
                        }
                    unsigned int poly_len;
                    memcpy(&poly_len, buffer + bytes_written, sizeof(unsigned int));
                    unsigned int poly_theo_len = poly_serial_length(p, &poly);
                    if (poly_len != poly_theo_len)
                        {
                            fprintf(stderr, "ERROR: %s:%d:%d invalid buffer length in constructed buffer %d %d\n",
                                    __FILE__, __LINE__, p->info_MPI.world_rank, poly_len, poly_theo_len);
                            return -5;
                        }
                    bytes_written += poly_bytes;
                }

            if (bytes_written != buffer_length)
                {
                    fprintf(stderr, "ERROR: %s:%d World rank %d: failed to serialize polymer %d %d\n",
                            __FILE__, __LINE__, p->info_MPI.world_rank, bytes_written, buffer_length);
                    return -3;
                }

            //Send the buffer to destination.
            int err = MPI_Send(buffer, buffer_length, MPI_UNSIGNED_CHAR, destination, 2, comm);
            if (err != MPI_SUCCESS)
                {
                    fprintf(stderr, "ERROR: %s:%d:%d MPI send not successfull %d %d\n", __FILE__, __LINE__,
                            p->info_MPI.world_rank, err, MPI_SUCCESS);
                    free(buffer);
                    return -4;
                }

            free(buffer);
            return Nsends;
        }

    return 0;
}

//!\private
unsigned char *recv_mult_polymers_core(const int source, const MPI_Comm comm,
                                       unsigned int *const Nsends, unsigned int *const buffer_length)
{
    MPI_Recv(Nsends, 1, MPI_UNSIGNED, source, 0, comm, MPI_STATUS_IGNORE);

    //Get buffer length from sending rank.
    MPI_Recv(buffer_length, 1, MPI_UNSIGNED, source, 1, comm, MPI_STATUS_IGNORE);

    if (*buffer_length > 0)     //The sending process can signal "no polymer to send."
        {
            //Allocate the memory to recv
            unsigned char *const buffer = (unsigned char *)malloc(*buffer_length * sizeof(unsigned char));
            if (buffer == NULL)
                {
                    fprintf(stderr, "MALLOC-ERROR: %s:%d size = %ld\n", __FILE__, __LINE__,
                            *buffer_length * sizeof(unsigned char));
                    return NULL;
                }

            MPI_Status status;
            //Recv the polymer buffer.
            int err = MPI_Recv(buffer, *buffer_length, MPI_UNSIGNED_CHAR, source, 2, comm, &status);

            int recv_count;
            MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &recv_count);
            unsigned char *ret = buffer;
            if ((unsigned int)recv_count != *buffer_length)
                {
                    fprintf(stderr, "ERROR: %s:%d:%d message size differs %d %d\n", __FILE__, __LINE__, source,
                            recv_count, *buffer_length);
                    ret = NULL;
                }
            if (err != MPI_SUCCESS)
                {
                    fprintf(stderr, "ERROR: %s:%d:%d MPI recv not successfull %d %d\n", __FILE__, __LINE__, source, err,
                            MPI_SUCCESS);
                    ret = NULL;
                }
            assert(status.MPI_SOURCE == source);
            assert(status.MPI_TAG == 2);

            return ret;
        }
    return NULL;
}

//!\private
int deserialize_mult_polymers(struct Phase *const p, const unsigned int Nsends,
                              const unsigned int buffer_length, const unsigned char *const buffer)
{
    //Deserialize and allocate Memory for recved polymer.
    Polymer poly;
    unsigned int bytes_read = 0;
    if (p->n_polymers + Nsends > p->n_polymers_storage)
        {
            reallocate_polymer_mem(p, p->n_polymers + Nsends);
        }

    for (unsigned int i = 0; i < Nsends; i++)
        {
            assert(bytes_read < buffer_length);
            int poly_bytes = deserialize_polymer(p, &poly, buffer + bytes_read);
            if (poly_bytes < 0)
                {
                    consider_compact_polymer_heavy(p, false);
                    poly_bytes = deserialize_polymer(p, &poly, buffer + bytes_read);
                }

            bytes_read += poly_bytes;
            //Push the polymer to the system if something has been read
            if (poly_bytes > 0)
                push_polymer(p, &poly);
            else
                {
                    fprintf(stderr, "ERROR: %s:%d rank %d invalid buffer length i=%d pb=%d \t %d %d %d\n",
                            __FILE__, __LINE__, p->info_MPI.world_rank, i, poly_bytes, bytes_read, buffer_length,
                            Nsends);
                    return -1;
                }
        }

    if (bytes_read != buffer_length)
        {
            fprintf(stderr, "ERROR: %s:%d rank %d invalid buffer length %d %d\n",
                    __FILE__, __LINE__, p->info_MPI.world_rank, bytes_read, buffer_length);
            return -1;
        }
    assert(bytes_read == buffer_length);
    return 0;
}

int recv_mult_polymers(struct Phase *const p, const int source, const MPI_Comm comm)
{
    assert(p);
    if (source >= p->info_MPI.world_size)
        {
            fprintf(stderr, "ERROR: %s:%d World rank %d requested to recv a polymer"
                    " from world rank %d, but only %d w. ranks are known.\n",
                    __FILE__, __LINE__, p->info_MPI.world_rank, source, p->info_MPI.world_size);
            return -2;
        }
    unsigned int Nsends;
    unsigned int buffer_length;
    unsigned char *const buffer = recv_mult_polymers_core(source, comm, &Nsends, &buffer_length);

    if (buffer != NULL)
        {
            const int err = deserialize_mult_polymers(p, Nsends, buffer_length, buffer);
            if (err != 0)
                {
                    fprintf(stderr, "ERROR: %s:%d World rank %d: invalid deserialization %d.\n",
                            __FILE__, __LINE__, p->info_MPI.world_rank, err);
                    return err;
                }
            free(buffer);
            return Nsends;
        }
    return 0;
}

int load_balance_mpi_ranks(struct Phase *const p)
{
    if (p->info_MPI.domain_size <= 1)
        return 0;

    double *const waiting_time = (double *const)malloc(p->info_MPI.domain_size * sizeof(double));
    if (waiting_time == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    double div = (p->info_MPI.domain_divergence_sec + 1) / (p->info_MPI.domain_divergence_counter + 1);

    MPI_Allgather(&(div), 1, MPI_DOUBLE, waiting_time, 1, MPI_DOUBLE, p->info_MPI.SOMA_comm_domain);
    //Reset counters.
    p->info_MPI.domain_divergence_sec = 0;
    p->info_MPI.domain_divergence_counter = 0;

    int arg_max = 0, arg_min = 0;
    for (int i = 0; i < p->info_MPI.domain_size; i++)
        {
            if (waiting_time[i] > waiting_time[arg_max])
                arg_max = i;
            if (waiting_time[i] < waiting_time[arg_min])
                arg_min = i;
        }

    //Quick exit, if no balance required.
    const double divergences = waiting_time[arg_max] - waiting_time[arg_min];
    free(waiting_time);
    const double seconds_per_step = p->tps_elapsed_time / p->tps_elapsed_steps;
    double p_waiting = divergences / seconds_per_step;

    //Last tps could be different on different rank, so avoid hangs.
    MPI_Bcast(&p_waiting, 1, MPI_DOUBLE, 0, p->info_MPI.SOMA_comm_domain);

    const double passed_acc = p->args.accepted_load_inbalance_arg / 100.;
    const double max_p_waiting = passed_acc >= 0 && passed_acc <= 1 ? passed_acc : 0.08;
    char quick_exit = p_waiting < max_p_waiting;
    if (quick_exit)
        return 0;

    unsigned int Nchains = p_waiting * p->n_polymers * 0.05 + 1;
    //Limit transferred chains to 5%
    if (Nchains > 0.05 * p->n_polymers)
        Nchains = 0.05 * p->n_polymers + 1;

    unsigned int Nsend = UINT_MAX;
    //! \todo Optimize chain to send. (Poly-Type based?)
    if (p->info_MPI.domain_rank == arg_min)
        {
            //CopyIN/OUT not the optimal solution, but the only option I can see so far
            copyout_phase(p);
            Nsend = send_mult_polymers(p, arg_max, Nchains, p->info_MPI.SOMA_comm_domain);
            consider_compact_polymer_heavy(p, false);
            copyin_phase(p);
        }

    if (p->info_MPI.domain_rank == arg_max)
        {
            //CopyIN/OUT not the optimal solution, but the only option I can see so far
            copyout_phase(p);
            Nsend = recv_mult_polymers(p, arg_min, p->info_MPI.SOMA_comm_domain);
            consider_compact_polymer_heavy(p, false);
            copyin_phase(p);
        }

    if (arg_min != 0)
        {
            if (p->info_MPI.domain_rank == arg_min)
                {
                    MPI_Send(&Nsend, 1, MPI_UNSIGNED, 0, 0, p->info_MPI.SOMA_comm_domain);
                }
            if (p->info_MPI.domain_rank == 0)
                {
                    MPI_Recv(&Nsend, 1, MPI_UNSIGNED, arg_min, 0, p->info_MPI.SOMA_comm_domain, MPI_STATUS_IGNORE);
                }
        }

    if (p->info_MPI.domain_rank == 0)
        printf("INFO world rank %d: Load balance @t=%d, sending %d chains from domain rank"
               " %d to domain rank %d, because of %f percent waiting.\n",
               p->info_MPI.world_rank, p->time, Nsend, arg_min, arg_max, p_waiting * 100.);

    //Synchronize rank, because otherwise a false unbalance would be
    //detected.
    MPI_Barrier(p->info_MPI.SOMA_comm_domain);
    return 1;
}

/*! Extract all polymer chains out of the current system, if their center of mass is in another domain.
    \private
    The serialized buffers of the chains are returned for later communication.
    Depending on the input only neighboring domains or all domains are checked.
    \param p Phase describing the system
    \param domain_lookup_list List of domains, where chains can be send to.
    \param len_domain_list Lenght of the domain list.
    \param n_sends Array of len_domain_list, where the number of chains to be send to the corresponding domain is stored.
    \param buffer_len Array to store the length of the buffers for each domain
    \param buffer 2D Array where all buffers of the chains to be send are stored.
    \param my_domain Domain of the calling rank
    \return Number of chains, which should be transferred, but no valid domain was found. This is an Error if > 0.
*/
int extract_chains_per_domain(struct Phase *const p, const int *domain_lookup_list, const unsigned int len_domain_list,
                              unsigned int *const n_sends, unsigned int *const buffer_len, unsigned char **const buffer,
                              const unsigned int my_domain)
{

    memset(n_sends, 0, len_domain_list * sizeof(unsigned int));
    memset(buffer_len, 0, len_domain_list * sizeof(unsigned int));

    unsigned int chains_missed = 0;
    for (unsigned int i = 0; i < p->n_polymers; i++)
        {
            const unsigned int target_domain = get_domain_id(p, &(p->polymers[i].rcm));
            if (target_domain != my_domain)
                {
                    const int domain_index = domain_lookup_list[target_domain];
                    if (domain_index >= 0)
                        {
                            assert((unsigned int)domain_index < len_domain_list);
                            n_sends[domain_index] += 1;
                            buffer_len[domain_index] += poly_serial_length(p, p->polymers + i);
                        }
                    else
                        chains_missed += 1;
                }
        }

    unsigned int *const buffer_offsets = (unsigned int *const)malloc(len_domain_list * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(buffer_offsets, len_domain_list * sizeof(unsigned int));

    //Allocate buffers
    for (unsigned int i = 0; i < len_domain_list; i++)
        {
            buffer_offsets[i] = 0;
            buffer[i] = (unsigned char *const)malloc(buffer_len[i] * sizeof(unsigned char));
            MALLOC_ERROR_CHECK(buffer[i], buffer_len[i] * sizeof(unsigned char));
        }

    Polymer poly;
    for (unsigned int i = 0; i < p->n_polymers; i++)
        {
            assert(i < p->n_polymers);
            const unsigned int target_domain = get_domain_id(p, &(p->polymers[i].rcm));
            if (my_domain != target_domain)
                {
                    const int domain_index = domain_lookup_list[target_domain];
                    if (domain_index >= 0)
                        {
                            if ((unsigned int)domain_index >= len_domain_list)
                                {
                                    fprintf(stderr, "ERROR %s:%d:%d invalid domain index %d %d\n",
                                            __FILE__, __LINE__, p->info_MPI.world_rank, domain_index, len_domain_list);
                                    return -2;
                                }
                            pop_polymer(p, i, &poly);
                            i -= 1;
                            const unsigned int poly_bytes =
                                serialize_polymer(p, &poly, buffer[domain_index] + buffer_offsets[domain_index]);
                            unsigned int poly_len;
                            memcpy(&poly_len, buffer[domain_index] + buffer_offsets[domain_index],
                                   sizeof(unsigned int));
                            unsigned int poly_theo_len = poly_serial_length(p, &poly);
                            if (poly_len != poly_theo_len)
                                {
                                    fprintf(stderr,
                                            "ERROR: %s:%d:%d invalid buffer length in constructed buffer %d %d\n",
                                            __FILE__, __LINE__, p->info_MPI.world_rank, poly_len, poly_theo_len);
                                    return -5;
                                }
                            buffer_offsets[domain_index] += poly_bytes;
                        }
                }
        }

    //Check buffer for consistency
    for (unsigned int i = 0; i < len_domain_list; i++)
        {
            if (buffer_offsets[i] != buffer_len[i])
                {
                    fprintf(stderr, "ERROR: %s:%d:%d Buffer construction inconsitent %d %d %d\n",
                            __FILE__, __LINE__, p->info_MPI.world_rank, buffer_offsets[i], buffer_len[i], i);
                    return -1;
                }
        }

    if (chains_missed != 0)
        printf
            ("ERROR: %s:%d: World rank %d has to send %d chains to a non neighbor rank. Restart simulation with higher rcm update rate.\n",
             __FILE__, __LINE__, p->info_MPI.world_rank, chains_missed);
    free(buffer_offsets);
    return chains_missed;
}

int send_domain_chains(struct Phase *const p, const bool init)
{
    if (p->args.N_domains_arg < 2)
        return 0;
    update_polymer_rcm(p);

    //Copy IN/OUT is not optimal, but the only option I can see so far.
    copyout_phase(p);

    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;

    const unsigned int len_domain_list = init ? p->args.N_domains_arg - 1 : 2;
    unsigned int *const domain_list = (unsigned int *const)malloc(len_domain_list * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(domain_list, len_domain_list * sizeof(unsigned int));
    int *const domain_lookup_list = (int *const)malloc(p->args.N_domains_arg * sizeof(int));
    MALLOC_ERROR_CHECK(domain_lookup_list, p->args.N_domains_arg * sizeof(int));
    unsigned int *const n_sends = (unsigned int *const)malloc(len_domain_list * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(n_sends, len_domain_list * sizeof(unsigned int));
    unsigned int *const send_len = (unsigned int *const)malloc(len_domain_list * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(send_len, len_domain_list * sizeof(unsigned int));
    unsigned char **const send_buffer = (unsigned char **const)malloc(len_domain_list * sizeof(unsigned char *));
    MALLOC_ERROR_CHECK(send_buffer, len_domain_list * sizeof(unsigned char *));

    unsigned int *const n_recv = (unsigned int *const)malloc(len_domain_list * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(n_recv, len_domain_list * sizeof(unsigned int));
    unsigned int *const recv_len = (unsigned int *const)malloc(len_domain_list * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(recv_len, len_domain_list * sizeof(unsigned int));
    unsigned char **const recv_buffer = (unsigned char **const)malloc(len_domain_list * sizeof(unsigned char *));
    MALLOC_ERROR_CHECK(recv_buffer, len_domain_list * sizeof(unsigned char *));
    if (init)
        {
            for (unsigned int i = 0; i < (unsigned int)p->args.N_domains_arg; i++)
                {
                    if (i < my_domain)
                        domain_list[i] = i;
                    else if (i > my_domain)
                        domain_list[i - 1] = i;
                }
        }
    else
        {
            domain_list[0] = my_domain == 0 ? (unsigned int)p->args.N_domains_arg - 1 : my_domain - 1;
            domain_list[1] = (int)my_domain == p->args.N_domains_arg - 1 ? 0 : my_domain + 1;
        }

    for (unsigned int i = 0; i < (unsigned int)p->args.N_domains_arg; i++)
        domain_lookup_list[i] = -1;     //Default
    for (unsigned int i = 0; i < len_domain_list; i++)
        domain_lookup_list[domain_list[i]] = i;

    const int missed_chains =
        extract_chains_per_domain(p, domain_lookup_list, len_domain_list, n_sends, send_len, send_buffer, my_domain);

    MPI_Request *const req = (MPI_Request * const)malloc(4 * len_domain_list * sizeof(MPI_Request));
    MALLOC_ERROR_CHECK(req, 4 * len_domain_list * sizeof(MPI_Request));
    MPI_Status *const status = (MPI_Status * const)malloc(4 * len_domain_list * sizeof(MPI_Status));
    MALLOC_ERROR_CHECK(status, 4 * len_domain_list * sizeof(MPI_Status));
    //Send out all the Nsends and buffer lengths
    for (unsigned int i = 0; i < len_domain_list; i++)
        {
            const unsigned int comm_domain = domain_list[i];
            const unsigned int comm_to = comm_domain * p->info_MPI.domain_size + p->info_MPI.domain_rank;
            MPI_Isend(n_sends + i, 1, MPI_UNSIGNED, comm_to, 0, p->info_MPI.SOMA_comm_sim,
                      req + (0 * len_domain_list + i));
            MPI_Isend(send_len + i, 1, MPI_UNSIGNED, comm_to, 1, p->info_MPI.SOMA_comm_sim,
                      req + (1 * len_domain_list + i));

            MPI_Irecv(n_recv + i, 1, MPI_UNSIGNED, comm_to, 0, p->info_MPI.SOMA_comm_sim,
                      req + (2 * len_domain_list + i));
            MPI_Irecv(recv_len + i, 1, MPI_UNSIGNED, comm_to, 1, p->info_MPI.SOMA_comm_sim,
                      req + (3 * len_domain_list + i));
        }

    //Finish the communication
    MPI_Waitall(4 * len_domain_list, req, status);
    for (unsigned int i = 0; i < len_domain_list; i++)
        {
            for (unsigned int j = 0; j < 4; j++)
                {
//Reset the requests
                    req[j * len_domain_list + i] = MPI_REQUEST_NULL;
                }
        }

    //Allocate space for recieving the polymers
    for (unsigned int i = 0; i < len_domain_list; i++)
        {
            recv_buffer[i] = (unsigned char *)malloc(recv_len[i] * sizeof(unsigned char));
            MALLOC_ERROR_CHECK(recv_buffer[i], recv_len[i] * sizeof(unsigned char));
        }
    //Send and recv all the buffers
    unsigned int total_chains_recv = 0;
    for (unsigned int i = 0; i < len_domain_list; i++)
        {
            const unsigned int comm_domain = domain_list[i];
            const unsigned int comm_to = comm_domain * p->info_MPI.domain_size + p->info_MPI.domain_rank;
            MPI_Isend(send_buffer[i], send_len[i], MPI_UNSIGNED_CHAR, comm_to, 2, p->info_MPI.SOMA_comm_sim,
                      req + (0 * len_domain_list + i));

            MPI_Irecv(recv_buffer[i], recv_len[i], MPI_UNSIGNED_CHAR, comm_to, 2, p->info_MPI.SOMA_comm_sim,
                      req + (1 + len_domain_list + i));
            total_chains_recv += n_recv[i];
        }
    MPI_Waitall(4 * len_domain_list, req, status);
    for (unsigned int i = 0; i < len_domain_list; i++)
        {
            for (unsigned int j = 0; j < 4; j++)
                {
                    //Deallocate the requests.
                    req[j * len_domain_list + i] = MPI_REQUEST_NULL;
                }
        }

    //Deserialize
    if (p->n_polymers + total_chains_recv > p->n_polymers_storage)
        reallocate_polymer_mem(p, p->n_polymers + total_chains_recv);
    int err = 0;
    for (unsigned int i = 0; i < len_domain_list; i++)
        {
            if (recv_buffer[i] != NULL)
                {
                    err = deserialize_mult_polymers(p, n_recv[i], recv_len[i], recv_buffer[i]);
                    if (err != 0)
                        {
                            fprintf(stderr, "ERROR: %s:%d World rank %d: invalid deserialization %d. "
                                    " i=%d\t %d %d %p\n",
                                    __FILE__, __LINE__, p->info_MPI.world_rank, err,
                                    i, n_recv[i], recv_len[i], recv_buffer[i]);
                            break;
                        }
                }
            else
                {
                    assert(n_recv[i] == 0);
                }
        }
    for (unsigned int i = 0; i < len_domain_list; i++)
        {
            if (recv_buffer[i] != NULL)
                {
                    free(recv_buffer[i]);
                    free(send_buffer[i]);
                }
        }
    free(domain_list);
    free(domain_lookup_list);
    free(req);
    free(status);
    free(n_sends);
    free(send_len);
    free(send_buffer);

    free(n_recv);
    free(recv_len);
    free(recv_buffer);

    consider_compact_polymer_heavy(p, true);
    //Copy IN/OUT is not optimal, but the only option I can see so far.
    copyin_phase(p);

    if (err == 0)
        return missed_chains;
    else
        return -2;
}

unsigned int get_domain_id(const struct Phase *const p, const Monomer * const rcm)
{
    const int fold = rint(rcm->x / p->Lx);
    soma_scalar_t x = rcm->x - fold * p->Lx;
    if (x < 0)
        x += p->Lx;
    const unsigned int target_domain = x / p->Lx * p->args.N_domains_arg;

    assert((int)target_domain < p->args.N_domains_arg);
    return target_domain;
}

#endif                          //ENABLE_MPI
