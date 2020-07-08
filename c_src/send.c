#include "send.h"
#include "mesh.h"
#include "serv_util.h"
#include <stdlib.h>
#include <assert.h>
#include "polymer.h"
#include "server.h"
#include "test.h"
#include "memory.h"
#include "err_handling.h"

void free_sender(struct sender * snd)
{
    MPI_Waitall(snd->req_num, snd->reqs, MPI_STATUSES_IGNORE);
    free(snd->reqs);
    free(snd->poly_buffer);
    for (int type = 0; type < snd->n_types; type++)
        {
            free(snd->fields[type]);
            free(snd->omega_fields[type]);
        }
    free(snd->fields);
    free(snd->omega_fields);
    free(snd->mv_acc);

}
int init_sender (struct sender *snd, const struct sim_rank_info *sim_inf, const struct Phase * p)
{

    // requests are: size of polymer array, number of polymers, polymer array (serialized),
    // acceptance rate, and
    // 2 fields for each monomer type ( density field and omega field )
    // see enum send_requests
    snd->req_num = 4 + 2*p->n_types;
    snd->reqs = malloc(snd->req_num * sizeof(MPI_Request));
    for (int i=0; i<snd->req_num; i++)
        {
            snd->reqs[i] = MPI_REQUEST_NULL;
        }
    snd->n_types = p->n_types;
    snd->fields = malloc( snd->n_types * sizeof(uint16_t *));
    snd->omega_fields = malloc(snd->n_types * sizeof(uint16_t *));
    snd->mv_acc = malloc(2 * sizeof(uint64_t));
    if (snd->reqs == NULL || snd->fields == NULL ||
        snd->omega_fields == NULL || snd->mv_acc == NULL)
        {
            fprintf(stderr, "Malloc error in function %s line %d of file %s\n",
                    __func__ , __LINE__, __FILE__);
            return -1;
        }
    for (int r=0; r < snd->req_num; r++)
        {
            snd->reqs[r] = MPI_REQUEST_NULL;
        }

    // poly_size and buffer will be initialized by polymer serialization (during sending)
    snd->poly_size = -1;
    snd->poly_buffer = NULL;

    // calculate field sizes and starting location for sending
    uint64_t x_min = p->local_nx_low + p->args.domain_buffer_arg;
    uint64_t x_max = p->local_nx_high - p->args.domain_buffer_arg - 1;
    uint64_t min_cell = cell_coordinate_to_index(p, x_min, 0, 0);
    uint64_t max_cell = cell_coordinate_to_index(p, x_max, p->ny - 1, p->nz - 1);
    uint64_t n_cells_to_send = max_cell - min_cell + 1;

    // consistency checks on field sending index calculations
    if (!p->args.skip_tests_flag)
        {
            // aborts on failure
            test_field_sending_consistency(p, sim_inf, min_cell, max_cell, n_cells_to_send);
        }



    snd->n_cells_to_send = n_cells_to_send;
    snd->min_cell = min_cell;

    for (unsigned int r=0; r < p->n_types; r++)
        {
            snd->fields[r] = malloc(n_cells_to_send * sizeof(uint16_t));
            snd->omega_fields[r] = malloc(n_cells_to_send * sizeof(uint16_t));
            if (snd->fields[r] == NULL || snd->omega_fields == NULL)
                {
                    fprintf(stderr, "Malloc error in function %s line %d of file %s\n",
                            __func__ , __LINE__, __FILE__);
                    return -1;
                }
        }

    return 0;
}


void get_open_reqs(struct sender *snd, MPI_Request *open_reqs, int *num_open){
    *num_open = 0;
    for (int r=0; r < snd->req_num; r++)
        {
            if (snd->reqs[r] != MPI_REQUEST_NULL)
                {
                    open_reqs[*num_open] = snd->reqs[r];
                    (*num_open)++;
                }
        }
}

int send_to_server (struct Phase *p, const struct sim_rank_info *sim_inf, struct sender *snd)
{

    /* note for future maintainers: The order in which non-blocking collectives
     * are posted must be the same in all ranks for MPI to match the data to the right request.
     * The simrank-server communication relies on that fact!
     * If you modify anything here, change receive_from_sim_ranks (server.c) accordingly.
     * see also: https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node126.htm
     */

    // on every timestep, test if requests finished
    int err, all_requests_finished;
    err = MPI_Testall(snd->req_num, snd->reqs, &all_requests_finished,
                      MPI_STATUSES_IGNORE);
    if (err != MPI_SUCCESS)
        {
            fprintf(stderr, "error in sending data to server");
            return -1;
        }

    // if no observables are to be calculated on this timestep, exit
    if (no_observables (&p->ana_info, p->time))
        {
            return 0;
        }

    if (!all_requests_finished)
        {
            // this blocks the client until al its requests are completed.
            // this is necessary for correctness but should not be called when SOMA is used with typical data
            // todo: find a good way to log this bad event

            MPI_Request * open_reqs = malloc(snd->req_num * sizeof(MPI_Request));
            int num_open;
            get_open_reqs(snd, open_reqs, &num_open);
            assert(num_open > 0);// if this was 0, all_requests_finished would be lying

            MPI_Waitall(num_open, open_reqs, MPI_STATUSES_IGNORE);
            for (int r=0; r < snd->req_num; r++)
                {
                    snd->reqs[r] = MPI_REQUEST_NULL;
                }
            free(open_reqs);
        }


    if (has_poly_obs(&p->ana_info, p->time))
        {
            // serialize polymers
            err = ser_all_poly (p, &(snd->poly_buffer), &(snd->poly_size));
            if (err)
                {
                    fprintf (stderr, "serializing polymers failed with"
                                     " error code %d "
                                     "(function %s, line %d, file %s)\n",
                             err, __func__, __LINE__, __FILE__);
                    return -1;
                }


            // tell server about size of polymer-array
            // todo: this step can be avoided many times when using the load-balancer-interval
            MPI_Igather (&snd->poly_size, 1, MPI_INT,
                         NULL, 0, MPI_INT,
                         0, sim_inf->poly_comm.comm, &snd->reqs[size_req]);

            // tell server about the number of polymers to be sent
            snd->n_polymers = p->n_polymers;
            MPI_Ireduce (&snd->n_polymers, &snd->n_polymers_recv_buf,
                         1, MPI_UINT64_T, MPI_SUM, 0, sim_inf->poly_comm.comm,
                         &snd->reqs[number_req]);


        }

    // send fields

    if (has_field_obs(&p->ana_info, p->time))
        {
            // send the density fields (not omega)
            if (sim_inf->field_comm.comm != MPI_COMM_NULL)
                {
                    for (unsigned int type=0; type < p->n_types; type++)
                        {
                            const uint64_t index = snd->min_cell + type * p->n_cells_local;
                            memcpy(snd->fields[type], &(p->fields_unified[index]), snd->n_cells_to_send);
                            MPI_Igatherv(snd->fields[type], snd->n_cells_to_send, MPI_UINT16_T,
                                         NULL, NULL, NULL, MPI_INT,
                                         0, sim_inf->field_comm.comm, &(snd->reqs[field_req + 2*type]));
                        }
                }
        }

    if (has_omega_field_obs(&p->ana_info, p->time))
        {
            // send the fields
            if (sim_inf->omega_field_comm.comm != MPI_COMM_NULL)
                {
                    for (unsigned int type=0; type < p->n_types; type++)
                        {
                            const uint64_t index = snd->min_cell + type * p->n_cells_local;
                            memcpy(snd->omega_fields[type], &(p->omega_field_unified[index]), snd->n_cells_to_send);
                            MPI_Igatherv(snd->fields[type], snd->n_cells_to_send, MPI_UINT16_T,
                                         NULL, NULL, NULL, MPI_INT,
                                         0, sim_inf->omega_field_comm.comm,  &(snd->reqs[omega_req + 2*type]));
                        }
                }
        }

    if (need_to_do(p->ana_info.delta_mc_acc_ratio, p->time))
        {
            // send moves and accepts
            snd->mv_acc[0] = p->n_moves;
            snd->mv_acc[1] = p->n_accepts;
            MPI_Igather(snd->mv_acc, 2, MPI_UINT64_T,
                        NULL, 0, MPI_UINT64_T,
                        0, sim_inf->poly_comm.comm, // poly_comm is also used for accepts
                        &snd->reqs[acc_rate_req]);
            p->n_moves = 0;
            p->n_accepts = 0;
        }

    if (has_poly_obs(&p->ana_info, p->time))
        {
            // give server the polymers
            MPI_Igatherv (snd->poly_buffer, snd->poly_size, MPI_UNSIGNED_CHAR, NULL, NULL, NULL, MPI_UNSIGNED_CHAR, 0, sim_inf->poly_comm.comm, &snd->reqs[poly_req]);

        }

    return 0;
}

int send_field_scaling_type(const soma_scalar_t *field_scaling_type, unsigned int n_types, const struct sim_rank_info * sri)
{

    // to make sure every server receives exactly one copy, we use the polymer communicator
    // only rank 1 sends, so every server gets exactly one message

    if (sri->poly_comm.rank != 1)
        return MPI_SUCCESS;

    assert(n_types <= INT_MAX);
    return MPI_Send(field_scaling_type, n_types, MPI_SOMA_SCALAR,
        0, 0, sri->poly_comm.comm);


}