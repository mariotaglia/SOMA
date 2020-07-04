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
            FILE * fp;
            fp = fopen("/home/julian/output_test/soma.log", "a");
            fprintf(fp, "ASSERTION FAILURE on rank %d. Message:\n"
                            "%s\n"
                            "expected (%s) to be true, but was false (function %s, file %s, line %d).\n"
                            "aborting\n",
                    world_rank, mess, str_expr, fun, file, line);
            fclose(fp);
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
            assert(serv_inf->server_comm.size == n_servers);

            // assign simranks to servers for poly-sending, the servers get rank 0 in the new communicator
            MPI_Comm_split(MPI_COMM_WORLD, serv_inf->server_comm.rank, 0, &(serv_inf->poly_comm.comm));
            complete_comm_with_info(&serv_inf->poly_comm);
            assert(serv_info->poly_comm.rank == 0);

            // assign simranks to servers for field-sending, the servers get rank 0
            MPI_Comm_split(MPI_COMM_WORLD, serv_inf->server_comm.rank, 0, &(serv_inf->field_comm.comm));
            complete_comm_with_info(&serv_inf->field_comm);
            assert(serv_info->field_comm.rank == 0);

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
            assert(sim_inf->sim_comm.size == world_size - n_servers);

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
            assert(sim_inf->domain_comm != MPI_COMM_NULL);
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

            ri->detailed_info = sim_inf;
        }

    return ri;
}




int receive_field_scaling_type(soma_scalar_t * field_scaling_type, unsigned int n_types, const struct server_info * si)
{
    // receive from rank 1, as rank 1 always exists and is a simrank
    assert(n_types < INT_MAX);
    DPRINT("server now receiving the field scaling");
    return MPI_Recv(field_scaling_type, n_types, MPI_SOMA_SCALAR, 1, 0, si->poly_comm.comm, MPI_STATUS_IGNORE);
}

