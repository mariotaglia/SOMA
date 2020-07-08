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

#include "server.h"
#include <stdlib.h>
#include <assert.h>
#include "cmdline.h"
#include "err_handling.h"
#include "serv_util.h"
#include "cmdline.h"
#include "polymer.h"

int complete_comm_with_info(struct comm_with_info * cwi)
{
    if (cwi->comm == MPI_COMM_NULL)
        {
            cwi->size = -1;
            cwi->rank = -1;
            return MPI_SUCCESS;
        }

    int err;
    err = MPI_Comm_rank(cwi->comm, &(cwi->rank));
    if (err != MPI_SUCCESS)
        {
            return err;
        }
    err = MPI_Comm_size(cwi->comm, &(cwi->size));
    if (err != MPI_SUCCESS)
        {
            return err;
        }

    assert(cwi->rank >= 0);
    assert(cwi->rank < cwi->size);

    return MPI_SUCCESS;

}
void m_assert(int expr, const char * str_expr, int line, const char * file, const char * fun,
              const char * mess)
{
    if (!expr)
        {
            int world_rank;
            MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
            fprintf(stderr, "ASSERTION FAILURE on rank %d. Message:\n"
                            "%s\n"
                            "expected (%s) to be true, but was false (function %s, file %s, line %d).\n"
                            "aborting\n",
                    world_rank, mess, str_expr, fun, file, line);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
}

void free_role_info(struct role_info * ri)
{
    free(ri->detailed_info);
    free(ri);
}

bool is_server(unsigned int n_servers, int world_rank, const int * server_ranks)
{
    MESSAGE_ASSERT(n_servers >= 1, "Must have at least one server");
    // get world size for a quick consistency check on the input
    int world_size;
    int status = MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_ERROR_CHECK(status, "failed to get world size");

    for (unsigned int i=0; i < n_servers; i++)
        {
            int serv_rank = server_ranks[i];
            MESSAGE_ASSERT(serv_rank >= 0, "invalid server-rank: server rank is negative.");
            MESSAGE_ASSERT(serv_rank < world_size, "invalid server-rank: server rank number is higher than exists in mpi_comm_world.");

            if (serv_rank == world_rank)
                {
                    return true;
                }
        }

    return false;
}

struct role_info * assign_role(unsigned int n_servers, int n_domains, const int *server_ranks)
{
    struct role_info * ri = malloc(sizeof(struct role_info));
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);


    // consistency checks
    MESSAGE_ASSERT(n_servers*2 <= (unsigned int)world_size, "Can not have more servers than simulation ranks.");
    MESSAGE_ASSERT((unsigned int)n_domains <= (world_size - n_servers), "There can not be more domains than simulation ranks.");
    MESSAGE_ASSERT( n_servers >= 1, "there must be at least one server");
    MESSAGE_ASSERT(n_domains >= 1, "there must be at least one domain");

    ri->is_server = is_server(n_servers, world_rank, server_ranks);

    if (ri->is_server)
        {

            struct server_info * serv_inf = malloc(sizeof(struct server_info));

            // split into server and simranks
            MPI_Comm_split(MPI_COMM_WORLD, 0, 0, &(serv_inf->server_comm.comm));
            complete_comm_with_info(&serv_inf->server_comm);
            assert((unsigned int)serv_inf->server_comm.size == n_servers);

            // assign simranks to servers for poly-sending, the servers get rank 0 in the new communicator
            MPI_Comm_split(MPI_COMM_WORLD, serv_inf->server_comm.rank, 0, &(serv_inf->poly_comm.comm));
            complete_comm_with_info(&serv_inf->poly_comm);
            assert(serv_inf->poly_comm.rank == 0);

            // assign simranks to servers for field-sending, the servers get rank 0
            MPI_Comm_split(MPI_COMM_WORLD, serv_inf->server_comm.rank, 0, &(serv_inf->field_comm.comm));
            complete_comm_with_info(&serv_inf->field_comm);
            assert(serv_inf->field_comm.rank == 0);

            // assign simranks to servers for omega-field-sending, the servers get rank0
            MPI_Comm_split(MPI_COMM_WORLD, serv_inf->server_comm.rank, 0, &(serv_inf->omega_field_comm.comm));
            complete_comm_with_info(&serv_inf->omega_field_comm);
            assert(serv_inf->omega_field_comm.rank == 0);

            ri->detailed_info = serv_inf;

        }
    else // simranks
        {
            struct sim_rank_info * sim_inf = malloc(sizeof(struct sim_rank_info));

            sim_inf->n_domains = n_domains;

            // split into server and simranks
            MPI_Comm_split(MPI_COMM_WORLD, 1, 0, &(sim_inf->sim_comm.comm));
            assert(sim_inf->sim_comm.comm != MPI_COMM_NULL);
            complete_comm_with_info(&sim_inf->sim_comm);
            assert((unsigned int)sim_inf->sim_comm.size == world_size - n_servers);

            // assign simranks to servers for polymer-sending
            MPI_Comm_split(MPI_COMM_WORLD, sim_inf->sim_comm.rank%n_servers, 1, &(sim_inf->poly_comm.comm));
            assert(sim_inf->poly_comm.comm != MPI_COMM_NULL);
            complete_comm_with_info(&sim_inf->poly_comm);
            assert(sim_inf->poly_comm.rank > 0);

            // domain decomposition (servers take no part in this)
            MESSAGE_ASSERT(sim_inf->sim_comm.size % n_domains == 0,
                "all domains must have same # of simranks,"
                "and this is not possible with these parameters " );
            int domain_size = sim_inf->sim_comm.size / n_domains;
            sim_inf->my_domain = sim_inf->sim_comm.rank / domain_size;
            MPI_Comm_split(sim_inf->sim_comm.comm, sim_inf->my_domain, 0, &(sim_inf->domain_comm.comm));
            assert(sim_inf->domain_comm.comm != MPI_COMM_NULL);
            complete_comm_with_info(&sim_inf->domain_comm);
            assert(sim_inf->domain_comm.size = domain_size);

            // assign simranks to servers for field-sending
            // the simrank with domain_rank 0 does all the sending for the domain, the rest do nothing

            int color;
            if (sim_inf->domain_comm.rank == 0)
                    color = sim_inf->my_domain % n_servers;
            else
                    color = MPI_UNDEFINED;
            // splitting with color MPI_UNDEFINED means you get MPI_COMM_NULL,
            // which is what we want for ranks that dont send the fields

            MPI_Comm_split(MPI_COMM_WORLD, color, 1, &(sim_inf->field_comm.comm));
            complete_comm_with_info(&sim_inf->field_comm);
            if (color == MPI_UNDEFINED)
                {
                    assert(sim_inf->field_comm.comm == MPI_COMM_NULL);
                    assert(sim_inf->field_comm.rank == -1);
                    assert(sim_inf->field_comm.size == -1);
                }
            else
                {
                    assert(color >= 0);
                    assert(sim_inf->field_comm.comm != MPI_COMM_NULL);
                    assert(sim_inf->field_comm.rank > 0);
                    assert(sim_inf->field_comm.size >= 2);
                }

            // assign simrank to domain for omega-field-sending.
            // If rank 1 exists in the domain_comm, it does the sending,
            // because rank 0 already sends the density fields.
            // rank 0 will only send the fields if it is the only rank in the domain

            int omega_sender = 0;
            if (sim_inf->domain_comm.size > 1)
                omega_sender = 1;

            if (sim_inf->domain_comm.rank == omega_sender)
                color = sim_inf->my_domain%n_servers;
            else
                color = MPI_UNDEFINED;

            MPI_Comm_split(MPI_COMM_WORLD, color, 1, &(sim_inf->omega_field_comm.comm));
            complete_comm_with_info(&sim_inf->omega_field_comm);

            if (color == MPI_UNDEFINED)
                {
                    assert(sim_inf->omega_field_comm.comm == MPI_COMM_NULL);
                    assert(sim_inf->omega_field_comm.rank == -1);
                    assert(sim_inf->omega_field_comm.size == -1);
                }
            else
                {
                    assert(color >= 0);
                    assert(sim_inf->omega_field_comm.comm != MPI_COMM_NULL);
                    assert(sim_inf->omega_field_comm.rank > 0);
                    assert(sim_inf->omega_field_comm.size >= 2);
                }

            ri->detailed_info = sim_inf;
        }

    return ri;
}

int init_receiver(struct receiver* rcv, const Ana_Info * ai, const struct server_info * si, const struct global_consts * gc, const struct som_args *args)
{

    int n_field_senders = si->poly_comm.size -1;

    if (ai->delta_mc_acc_ratio != 0)
        {
#ifdef _OPENACC
            MESSAGE_ASSERT(-1, "acc-ratio observable not available with openacc-builds");
#endif
            rcv->mv_acc = malloc(2 * sizeof(uint64_t) * si->poly_comm.size);
            RET_ERR_ON_NULL(rcv->mv_acc, "Malloc");
        }
    else
        {
            rcv->mv_acc = NULL;
        }

    // the buffer for serialized polymers and its size will be initialized during receiving, when the size is actually known (it might change over time due to load balancing)
    rcv->poly_size = -1;
    rcv->poly_buffer = NULL;

    rcv->n_polymers = 0;
    rcv->polymers = NULL;

    rcv->n_types = gc->n_types;

    uint64_t n_cells_total = gc->nx * gc->ny * gc->nz;
    uint64_t n_cells_per_domain = n_cells_total / args->N_domains_arg;
    rcv->field_size = n_cells_per_domain * n_field_senders;
    // field size is the number of cells that this server deals with.
    // it may not be the same for all servers.

    rcv->field = malloc(sizeof(uint16_t *) * gc->n_types);
    RET_ERR_ON_NULL(rcv->field, "Malloc");
    rcv->omega_field = malloc(sizeof(uint16_t *) * gc->n_types);
    RET_ERR_ON_NULL(rcv->omega_field, "Malloc");

    for (unsigned int type=0; type < gc->n_types; type++)
        {
            rcv->field[type] = malloc(sizeof(uint16_t) * rcv->field_size);
            RET_ERR_ON_NULL(rcv->field[type], "Malloc");
            rcv->omega_field[type] = malloc(sizeof(uint16_t) * rcv->field_size);
            RET_ERR_ON_NULL(rcv->omega_field, "Malloc");
        }

    //dspls and offset arrays for Igatherv are the same for every receive,
    //create them once
    rcv->fields_recvcount = malloc(sizeof(int) * (n_field_senders + 1));
    RET_ERR_ON_NULL(rcv->fields_recvcount, "Malloc");
    rcv->fields_dspls = malloc(sizeof(int) * (n_field_senders + 1));
    RET_ERR_ON_NULL(rcv->fields_dspls, "Malloc");
    rcv->fields_recvcount[0] = 0;
    rcv->fields_dspls[0] = 0;
    for (int i=1; i <= n_field_senders; i++)
        {
            rcv->fields_recvcount[i] = n_cells_per_domain;
            rcv->fields_dspls[i] = rcv->fields_recvcount[i-1] + rcv->fields_dspls[i-1];
        }
    return 0;
}

void free_receiver(struct receiver* rcv)
{
    free(rcv->poly_buffer);
    for (unsigned int type=0; type < rcv->n_types; type++)
        {
            free(rcv->field[type]);
            free(rcv->omega_field[type]);
        }
    free(rcv->field);
    free(rcv->omega_field);
    free(rcv->fields_recvcount);
    free(rcv->fields_dspls);
}

int receive_from_sim_ranks(const struct server_info * si, const Ana_Info * ai, struct receiver * rcv, unsigned int time, const struct som_args * args, const struct global_consts * gc)
{
    /* note for future maintainers: The order in which non-blocking collectives
     * are posted must be the same in all ranks for MPI to match the data to the right request.
     * The simrank-server communication relies on that fact!
     * If you modify anything here, change send_to_server (send.c) accordingly.
     * see also: https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node126.htm
     */

    // use requests from this simple stack, make it big enough for every possible request
    const unsigned int req_num_max = 4 + 2*rcv->n_types;
    MPI_Request * reqs = malloc(sizeof(MPI_Request) * req_num_max );
    RET_ERR_ON_NULL(reqs, "Malloc");
    for (unsigned int i=0; i<req_num_max; i++)
        reqs[i] = MPI_REQUEST_NULL;
    int req_num=0;


    MPI_Request poly_meta_req[2]; //polymer serialized buffer size request and request for number of polymers is special, we want to wait for it first

    if (no_observables(ai, time))
        {
            return 0;
        }

    int * pol_sz = NULL;
    if (has_poly_obs(ai, time))
        {
            // receive size of serialized polymer-array
            pol_sz = malloc(sizeof(int) * si->poly_comm.size);
            RET_ERR_ON_NULL(pol_sz, "Malloc");
            int sz = 0; //server sends dummy-size
            MPI_Igather(&sz, 1, MPI_INT,
                pol_sz, 1, MPI_INT,
                0, si->poly_comm.comm, &poly_meta_req[0]);
            rcv->n_polymers = 0;
            MPI_Ireduce(MPI_IN_PLACE, &rcv->n_polymers,
                1, MPI_UINT64_T, MPI_SUM, 0, si->poly_comm.comm,
                &poly_meta_req[1]);
        }

    // receive fields

    if (has_field_obs(ai, time))
        {
            // receive fields
            for (unsigned int type=0; type < rcv->n_types; type++)
                {
                    MPI_Igatherv(NULL, 0, MPI_UINT16_T,
                                 rcv->field[type], rcv->fields_recvcount, rcv->fields_dspls, MPI_UINT16_T,
                                 0, si->field_comm.comm, &reqs[req_num++]);
                }
        }

    if (has_omega_field_obs(ai, time))
        {
            // receive omega fields
            for (unsigned int type=0; type < rcv->n_types; type++)
                {
                    MPI_Igatherv(NULL, 0, MPI_UINT16_T,
                                 rcv->omega_field[type], rcv->fields_recvcount, rcv->fields_dspls, MPI_UINT16_T,
                                 0, si->omega_field_comm.comm, &reqs[req_num++]);
                }
        }

    if (need_to_do(ai->delta_mc_acc_ratio, time))
        {
            // server also needs to send something for the gather
            uint64_t dummy_send[2] = {0,0};
            MPI_Igather(dummy_send, 2, MPI_UINT64_T,
                rcv->mv_acc, 2, MPI_UINT64_T,
                0, si->poly_comm.comm,
                &reqs[req_num++]);
        }

    if (has_poly_obs(ai, time))
        {
            int status;
            status = MPI_Waitall(2, poly_meta_req, MPI_STATUSES_IGNORE);
            MPI_ERROR_CHECK(status, "Waiting for the polymer array sizes failed");
            // pol_sz is now the recv_count array for the Igatherv
            // we build the dspls-array from it and find out the total size.
            int *dspls = malloc(sizeof(int) * si->poly_comm.size);
            RET_ERR_ON_NULL(dspls, "Malloc");
            int total_size = 0;
            dspls[0] = 0;
            for (int i=1; i < si->poly_comm.size; i++)
                {
                    dspls[i] = dspls[i-1] + pol_sz[i-1];
                    total_size += pol_sz[i];
                }
            rcv->poly_size = total_size;
            free(rcv->poly_buffer);
            rcv->poly_buffer = malloc(rcv->poly_size * sizeof(unsigned char));


            MPI_Request poly_req;
            MPI_Igatherv(NULL, 0, MPI_UNSIGNED_CHAR, rcv->poly_buffer, pol_sz, dspls, MPI_UNSIGNED_CHAR, 0, si->poly_comm.comm, &poly_req);

            MPI_Wait(&poly_req, MPI_STATUS_IGNORE);


            // si->poly_buffer contains all serialized polymers, we just need to unpack it
            free(rcv->polymers);
            rcv->polymers = malloc(sizeof(Polymer) * rcv->n_polymers);
            RET_ERR_ON_NULL(rcv->polymers, "cannot malloc new memory for polymers on server");

            // this loop is not parallel and that is bad.
            // todo: optimize with new polymer layout.
            size_t position = 0;
            for (size_t i=0; i < rcv->n_polymers; i++)
                {
                    int read = deserialize_poly_server(&(rcv->polymers[i]), rcv->poly_buffer + position, gc, args->pseudo_random_number_generator_arg);
                    if (read < 0)
                        {
                            fprintf(stderr, "failed to deserialize polymer"
                                            " %ld, read was %d "
                                            "(line %d, file %s, function %s)\n",
                                i, read, __LINE__, __FILE__, __func__);
                            MESSAGE_ASSERT(false, "deserializing polymer failed");
                        }
                    position += read;
                }
            free(dspls);
        }

    MPI_Waitall(req_num, reqs, MPI_STATUSES_IGNORE);
    /*
    for (int i=0; i< req_num; i++)
        {
            int fin;
            MPI_Waitany(req_num, reqs, &fin, MPI_STATUS_IGNORE);
        }
    */
    free(reqs);
    free(pol_sz); // freeing null-pointers is a no-op.

    return 0;
}



int receive_field_scaling_type(soma_scalar_t * field_scaling_type, unsigned int n_types, const struct server_info * si)
{
    // receive from rank 1, as rank 1 always exists and is a simrank
    assert(n_types < INT_MAX);
    return MPI_Recv(field_scaling_type, n_types, MPI_SOMA_SCALAR, 1, 0, si->poly_comm.comm, MPI_STATUS_IGNORE);
}

