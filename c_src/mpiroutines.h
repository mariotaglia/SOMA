/* Copyright (C) 2016-2021 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016 Marcel Langenberg
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
#ifndef SOMA_MPIROUTINES_H
#define SOMA_MPIROUTINES_H

#include "soma_config.h"
#if ( ENABLE_MPI == 1 )
#include <mpi.h>
#endif                          //ENABLE_MPI
#ifdef ENABLE_NCCL
#include <nccl.h>
#endif                          //ENABLE_NCCL
#include <stdint.h>
#include <stdbool.h>
struct Phase;

/*! \file mpiroutines.h
  \brief Header file for functions that require MPI calls for SOMA.
*/

//! \brief All information about the MPI setup.
//! Information about the MPI setup of the simulation.
//! There are three different levels in the MPI hierachy.
//!
//! a) WORLD: This is the coarsed level and contains all MPI ranks
//! assigned to the simulation via mpiexec
//!
//! b) SIMULATION: The world communicator can be split in multiple
//! simulations. Each simulation is basically an indepent instance of
//! SOMA. Each simulation can simulate a different system and is
//! usually supplied by an individual input file and analysis file.
//! SOMA itself does not do any communication between different
//! simulations. For meta-dynamics like the string algorithms this
//! instance is useful. Each simulation can be split into multiple domains.
//!
//! c) DOMAIN: A simulation can be split into different spatial domains.
//! This finest level is used to communicate between the different domains.
//! Each domain can have multiple mpi ranks.
typedef struct Info_MPI {
    int world_rank;             //!< Rank of the world communicator
    int world_size;             //!< Size of the world communicator
    int sim_size;               //!< Size of the simulation communicator
    int sim_rank;               //!< Rank of the simulation communicator
    int domain_size;            //!< Size of a single domain communicator
    int domain_rank;            //!< Rank of a single domain communicator
    int gpu_id;                 //!< ID of the GPU used. if < 0 no GPU used on this rank.
#if ( ENABLE_MPI == 1 )

    MPI_Comm SOMA_comm_world;   //!< Global communicator for 1 simulation
    MPI_Comm SOMA_comm_sim;     /*!< \brief communicator within one conf, SCMF parallelization */
    MPI_Comm SOMA_comm_domain;  /*!< \brief communicator within one domain of a SCMF simulation parallelization */
#ifdef ENABLE_NCCL
    ncclComm_t SOMA_nccl_world; //!< \brief NCCL communicator for world communication
    ncclComm_t SOMA_nccl_sim;   //!< NCCL communicator for sim MPI communicator
    ncclComm_t SOMA_nccl_domain;        //!< NCCL communicator for domain communicator
#endif                          //ENABLE_NCCL
    MPI_Status mpi_status;      //!< Status of the mpi init.
#endif                          //ENABLE_MPI
    //! Store MPI divergence in between domain ranks.
    double domain_divergence_sec;
    //! Counter for the MPI divergence in between domain ranks.
    unsigned int domain_divergence_counter;
} Info_MPI;

/*! \brief Initialize MPI.

  Initialize the MPI-enviroment for any further calls to the MPI routines of SOMA.

  \param p Phase struct, which defines the current state of the
  simulation.
  \note This functions is the only function, which does
  not require a completely initialized Phase struct.
  \return Error code. Error occured if not equal to 0.
  \post All MPI related parts of the Phase struct are initialized.
  \pre Initialized p->info_MPI.SOMA_MPI_Comm.
*/
int init_MPI(struct Phase *p);

//! \brief wrapper for MPI_Finalize
//! \param mpi Input of initialized MPI data.
//! \post the communicators are freed
//! \return Errorcode
int finalize_MPI(struct Info_MPI *mpi);

//! \brief Function to check wheter one MPI-rank passed a non-zero value.
//!
//! Useful for error checking and grace fully exiting MPI.
//! \param p Phase containing the communicator and more.
//! \param my_status Status of the local rank.
//! \return Nonzero value if one MPI-rank passed a non-zero value.
int check_status_on_mpi(const struct Phase *const p, int my_status);

#if ( ENABLE_MPI == 1 )
//! Measure divergence of MPI ranks with an MPI_Barrier call.
//!
//! \param p System which running the simulation. (Reqired for MPI context.)
//! \note this function also updates summed counters in info_MPI
//! \return seconds waiting in Barrier.
double mpi_divergence(struct Phase *const p);

//! Update global properties, which can be combined from local statistics.
//!
//! If you change some local properties, which need to be covered globally,
//! call this function.
//! \param p System to update.
//! (Insert new or completly deleting polymer from the global system
//! is by a good example, since it changes the global number of polymers and
//! the global number of beads and beads per type.)
//! \note This function is MPI collective.
//! \return Errorcode.
int collective_global_update(struct Phase *const p);

//! Send a polymer from one MPI rank to another.
//!
//! \warning No assumptions about the global system are made.
//! If something changes you have to update it afterwards.
//! \param p System.
//! \param poly_id Id of the polymer to send.
//! \param destination Id of the recving MPI rank.
//! \param comm MPI communicator in which the polymer should be sent
//! \return Errorcode.
int send_polymer_chain(struct Phase *const p, const uint64_t poly_id, const int destination, const MPI_Comm comm);

//! Recv a polymer from one MPI rank to another.
//!
//! \warning No assumptions about the global system are made.
//! If something changes you have to update it afterwards.
//! \param p System.
//! \param source Id of the recving MPI rank.
//! \param comm MPI communicator in which the polymer should be recieved
//! \return Errorcode.
int recv_polymer_chain(struct Phase *const p, const int source, const MPI_Comm comm);

//! Send a multiple polymers from one MPI rank to another.
//!
//! \warning No assumptions about the global system are made.
//! If something changes you have to update it afterwards.
//! \param p System.
//! \param Nsends Number of polymers to send
//! \param destination Id of the recving MPI rank.
//! \param comm MPI communicator in which the polymer should be sent
//! \return >= 0 Number of polymers send. else Errorcode
int send_mult_polymers(struct Phase *const p, const int destination, unsigned int Nsends, const MPI_Comm comm);

//! Helper function to obtain the malloc allocated buffer containing polymers from another rank.
//!
//! \param source Source rank
//! \param comm MPI_Communicator for the communication
//! \param Nsends Output parameter for the number of chains
//! \param buffer_length Ouput parameter for the length of allocated buffer
//! \return Pointer to buffer if successful, NULL otherwise
unsigned char *recv_mult_polymers_core(const int source, const MPI_Comm comm,
                                       unsigned int *const Nsends, unsigned int *const buffer_length);

//! Helper function to deserialize and pop in multiple polymer in a buffer.
//!
//! \param p Phase of the system
//! \param Nsends number of polymers in buffer
//! \param buffer_length Length of the buffer
//! \param buffer Pointer to the buffer
//! \return Errorcode
int deserialize_mult_polymers(struct Phase *const p, const unsigned int Nsends,
                              const unsigned int buffer_length, const unsigned char *const buffer);

//! Recv multiple polymers from one MPI rank to another.
//!
//! \warning No assumptions about the global system are made.
//! If something changes you have to update it afterwards.
//! \param p System.
//! \param source Id of the recving MPI rank.
//! \param comm MPI communicator in which the polymer should be sent
//! \return Errorcode.
int recv_mult_polymers(struct Phase *const p, const int source, const MPI_Comm comm);

//! Load balance the MPI ranks.
//!
//! The info_MPI divergence_sec are evaluated and if there is a
//! difference of more than 1/50. s waiting time in the barrier
//! between the ranks, the slowest rank sends a single chain to the
//! fastest rank.
//! The load balancing is done inside an MPI domain
//! \note This function is not designed for frequent calls.
//! \note This function is MPI collective.
//! \param p System.
//! \return Errorcode
int load_balance_mpi_ranks(struct Phase *const p);

//! Distribute all chains to their corresponding rank for a domain decomposition
//!
//! \param p initialied configuration

//! \param init set to true for initial redistribution, chains can be
//! send to all ranks. In the false mode, only neighbor ranks can be
//! addressed. If chains cannot be sent to their corresponding rank an
//! error is issued.
//! \return Number of unsent chains should be 0, otherwise error.
int send_domain_chains(struct Phase *const p, const bool init);

#endif                          //ENABLE_MPI

#endif                          /*SOMA_MPIROUTINES_H */
