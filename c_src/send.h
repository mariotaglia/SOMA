#ifndef _SEND_H_
#define _SEND_H_

#include <inttypes.h>
#include <mpi.h>
#include <stdbool.h>
#include "phase.h"
#include "server.h"


/*! \brief object used by a simulation rank to manage data sending.
 * \note this is to be initialized with init_sender before polymers and fields can be send
 * using send_to_server. It must be freed with destroy_sender.
 */
struct sender{

    // total number of requests that can be made
    // this is also the size of the reqs, req_is_open and sendbufs array
    int req_num;

    // all created mpi-requests. reqs[i] == MPI_REQUEST_NULL iff request is completed or never made.
    // the poly-size-request, the polymer-request and the request for the acceptance rate are found at the
    // index corresponding to the send_request enum constant. The requests for fields and omega-fields
    // are found at reqs[field_req + 2*types] (fields) and reqs[field_req + 2*types + 1] (omega-field)
    // where type fulfills 0 <= type < p->n_types
    MPI_Request *reqs;

    //true for open, false for completed or never started
    bool *req_is_open;

    // the buffer for the serialized polymers and its size
    uint64_t poly_size;
    unsigned char * poly_buffer;

    // size of the fields to be sent (number of cells, not bytes)
    uint64_t field_size;

    // every element is a pointer to the buffer used for sending one field.
    // The size of the array at omega_fields[i] and fields[i] is field_size
    uint16_t ** fields;
    uint16_t ** omega_fields;

    //indices in the p->fields_unified array that determine the first cell
    // to be send for each monomer type
    uint64_t min_cell;

    // mv_acc[0] == number of moves
    // mv_acc[1] == number of accepts
    uint64_t * mv_acc;

    // number of monomer types
    int n_types;
};


//! meant to index the request-array used by the sender
//! reqs[size_req], reqs[poly_req], reqs[acc_rate_req] refer to size, polymer and acceptance rate requests
//! reqs[field_req + 2*i] refers to the field-request of monomer type i
//! reqs[omega_req + 2*i] refers to the omega_field_request of monomer type i
enum send_requests{
    size_req = 0,
    poly_req = 1,
    acc_rate_req = 2,
    field_req = 3,
    omega_req = 4
};


//! initializes the sender so that send_to_server can be called on it
//! \param snd sender object to be initialized that you are the owner of
//! \param p current configuration (to determine field sizes, number of monomer types etc)
//! \return error code, 0 is for success
int init_sender (struct sender *snd, struct sim_rank_info * sim_inf, const struct Phase * p);

//! frees all buffers associated with the sender, it cant be used afterwards
//! this also waits for all pending requests to complete.
void destroy_sender(struct sender * snd);

//! called by a simrank
//! sends the data that the server needs (according to ana_info and time) to the server that belongs to this simrank.
//! is smart enough to return immediately if server does not need new data, so can and should be called on every timestep
//! by the client
//! \param p phase to be used for sending
//! \return 0 for success, otherwise error code
int send_to_server (struct Phase *p, const struct sim_rank_info *sim_inf, struct sender *snd);

//! called by all simranks, while the servers call receive_field_scaling_type
//! \param field_scaling_type (in) fully initialized field scaling type array (has  (p->field_scaling_type)
//! \param n_types number of monomer types (p->n_types)
//! \param sri mpi information about the simrank
//! \return MPI_SUCCESS or an error code (if mpi is set to return errors)
int send_field_scaling_type(const soma_scalar_t *field_scaling_type, unsigned int n_types, const struct sim_rank_info * sri);

#endif //_SEND_H_
