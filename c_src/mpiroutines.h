/* Copyright (C) 2016 Ludwig Schneider
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

#include <mpi.h>
#include <stdint.h>
struct Phase;

/*! \file mpiroutines.h
  \brief Header file for functions that require MPI calls for SOMA.
*/

//! \brief All information about the MPI setup.
typedef struct Info_MPI{
    int colour;
    int Ncores;	/*!< \brief total number of cores    */
    int current_core;	/*!< \brief current core index   */
    MPI_Comm  	SOMA_MPI_Comm;	/*!< \brief communicator within one conf, SCMF parallelization */
    MPI_Status     	mpi_status; //!< Status of the mpi init.
    //! Store MPI divergence of all ranks.
    double divergence_sec;
    //! Counter of all mpi divergence calls.
    unsigned int divergence_counter;
    /* ... */
    }Info_MPI;

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

//! \brief Function to check wheter one MPI-rank passed a non-zero value.
//!
//! Useful for error checking and grace fully exiting MPI.
//! \param p Phase containing the communicator and more.
//! \param my_status Status of the local rank.
//! \return Nonzero value if one MPI-rank passed a non-zero value.
int check_status_on_mpi(const struct Phase*const p,int my_status);

//! Measure divergence of MPI ranks with an MPI_Barrier call.
//!
//! \param p System which running the simulation. (Reqired for MPI context.)
//! \note this function also updates summed counters in info_MPI
//! \return seconds waiting in Barrier.
double mpi_divergence(struct Phase*const p);

//! \brief wrapper for MPI_Finalize
//! \return Errorcode
int finalize_MPI(void);

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
int collective_global_update(struct Phase*const p);

//! Send a polymer from one MPI rank to another.
//!
//! \warning No assumptions about the global system are made.
//! If something changes you have to update it afterwards.
//! \param p System.
//! \param poly_id Id of the polymer to send.
//! \param destination Id of the recving MPI rank.
//! \return Errorcode.
int send_polymer_chain(struct Phase*const p, const uint64_t poly_id, const int destination);

//! Recv a polymer from one MPI rank to another.
//!
//! \warning No assumptions about the global system are made.
//! If something changes you have to update it afterwards.
//! \param p System.
//! \param source Id of the recving MPI rank.
//! \return Errorcode.
int recv_polymer_chain(struct Phase*const p, const int source);

//! Send a multiple polymers from one MPI rank to another.
//!
//! \warning No assumptions about the global system are made.
//! If something changes you have to update it afterwards.
//! \param p System.
//! \param Nsends Number of polymers to send
//! \param destination Id of the recving MPI rank.
//! \return >= 0 Number of polymers send. else Errorcode
int send_mult_polymers(struct Phase*const p,const int destination,unsigned int Nsends);

//! Recv multiple polymers from one MPI rank to another.
//!
//! \warning No assumptions about the global system are made.
//! If something changes you have to update it afterwards.
//! \param p System.
//! \param source Id of the recving MPI rank.
//! \return Errorcode.
int recv_mult_polymers(struct Phase*const p, const int source);

//! Load balance the MPI ranks.
//!
//! The info_MPI divergence_sec are evaluated and if there is a
//! difference of more than 1/50. s waiting time in the barrier
//! between the ranks, the slowest rank sends a single chain to the
//! fastest rank.
//! \note This function is not designed for frequent calls.
//! \note This function is MPI collective.
//! \param p System.
//! \return Errorcode
int load_balance_mpi_ranks(struct Phase*const p);
#endif/*SOMA_MPIROUTINES_H*/
