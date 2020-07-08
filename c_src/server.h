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


#ifndef _SERVER_H_
#define _SERVER_H_

#include <stdbool.h>
#include <mpi.h>
#include <stdio.h>
#include <inttypes.h>
#include "soma_config.h"
#include "phase.h"

#define MESSAGE_ASSERT(expr, mess) m_assert(expr, #expr, __LINE__, __FILE__, __func__, mess)

void m_assert(int expr, const char * str_expr, int line, const char * file, const char * fun,
              const char * mess);

/*!
 * convenience wrapper around an MPI_Comm that also includes the size of the communicator
 * and the rank that this rank has in this communicator.
 * Initialize the comm manually, then call complete_comm_with_info
 */
struct comm_with_info
{
    MPI_Comm comm;
    int size;
    int rank;
};

//! takes a comm_with_info where the communicator is initialized and fills in the other information.
//! MPI_COMM_NULL is valid, in that case size and rank are set to -1
//! \param cwi the comm_with_info to be filled in
//! \return MPI_SUCCESS or an error code (if MPI is set to return on error)
int complete_comm_with_info(struct comm_with_info * cwi);

/*!
 * \brief object to store the role of a rank in the simulation
 * if is_server is true, the rank is a server, and detailed_info should be cast to server_info
 * otherwise, the rank is a simrank, and detailed_info is of type sim_rank_info
 */
struct role_info
{
    bool is_server;
    void * detailed_info;
};

/*!
 * contains information specific to a server.
 * This has the communicators (and satellite information: size and rank)
 * on which the server expects fields and polymers, also a communicator with all the other servers
 */
struct server_info
{
    struct comm_with_info server_comm;
    struct comm_with_info field_comm;
    struct comm_with_info omega_field_comm;
    struct comm_with_info poly_comm;
};

struct sim_rank_info
{
    struct comm_with_info sim_comm;
    struct comm_with_info domain_comm;
    struct comm_with_info poly_comm;
    struct comm_with_info field_comm;
    struct comm_with_info omega_field_comm;

    int n_domains;
    int my_domain;
};

struct receiver
{
    /*! contains the moves and accepts of the simulation ranks since the last time acc_ratio was calculated
     * access: even indices for moves, uneven for accepts, size is the number of polymer-senders of this server times 2.
     */
    uint64_t *mv_acc;

    uint64_t poly_size;
    unsigned char * poly_buffer;

    Polymer * polymers; //! (deserialized) polymer array
    uint64_t n_polymers; //! size of polymer array / number of polymers on this server

    unsigned int n_types; //! number of monomer types in the system (copied from global constants for convenience)
    size_t field_size; //! size of one field for one monomer-type
    uint16_t ** field; //! field[i] is the density field for monomer type i with size field_size
    uint16_t ** omega_field;//! omega_field[i] is the omega-field for monomer type i with size field_size.
    int * fields_recvcount; //! fields_recvcount[i] says how many elements are expected from field rank i (input for Igatherv of fields)
    int * fields_dspls; //! fields_dspls[i] says where to put data from rank i (input for Igatherv of fields)


};

void free_role_info(struct role_info * ri);

//! finds the role of my rank. Every rank must call this at the beginning of the program.
//!
//! \param n_servers
//! \param n_domains
//! \param server_ranks
//! \return
struct role_info * assign_role(unsigned int n_servers, int n_domains, const int *server_ranks);

//! determines if the given world rank is a server or a simulation-rank, given a list of all the servers
//! \param n_servers number of servers (==length of the server_ranks array)
//! \param world_rank world rank of the rank we want to know about.
//! \param server_ranks list of servers
//! \note will perform some consistency checks
//! \return true if world_rank is a server
bool is_server(unsigned int n_servers, int world_rank, const int * server_ranks);

//! receive the field_scaling_type. All servers call this to receive, all simranks call send_field_scaling_type to send
//! \param field_scaling_type (out) Array of n_types that you are the owner of.
//! \param n_types length of the array
//! \param si mpi information about the server
//! \return MPI_SUCCESS or error (if mpi is set to return errors)
int receive_field_scaling_type(soma_scalar_t * field_scaling_type, unsigned int n_types, const struct server_info * si);

//! receive the data from the simulation ranks
//! The new polymers and fields will be found in the receiver
//! \param si MPI information about this server
//! \param ai information about what to analyze in what timestep
//! \param rcv (out) the receiver-object of this server, which will carry the field-information
//! \param time the current timestep for which data is being requested
//! \param args command-line-arguments of the simulation
//! \param gc global constants of the simulation
//! \return 0 for success, otherwise error.
int receive_from_sim_ranks(const struct server_info * si, const Ana_Info * ai, struct receiver * rcv, unsigned int time, const struct som_args * args, const struct global_consts * gc);


//! initialize the receiver object for use in receive_from_simranks
//! \param rcv receiver to be initialized
//! \param ai information about the analysis of the simulation
//! \param si MPI-information about this server-rank
//! \param gc global constants of the simulation
//! \param args command-line-arguments of the simulation
//! \return 0 for success, otherwise failure
int init_receiver(struct receiver* rcv, const Ana_Info * ai, const struct server_info * si, const struct global_consts * gc, const struct som_args *args);
#endif //_SERVER_H_
