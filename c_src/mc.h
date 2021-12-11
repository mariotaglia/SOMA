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
/*! \file mc.h
  \brief Functions required for any Monte-Carlo move of the system.
*/

#ifndef SOMA_MC_H
#define SOMA_MC_H

#include <stdbool.h>
struct Phase;
#include "soma_config.h"
#include "rng.h"
#include "monomer.h"
#include "polymer.h"
struct Polymer;
//! Main Monte-Carlo function.
//!
//! Automatic selection of the Monte-Carlo algorithm and their execution.
//! \param p System to propagate.
//! \param nsteps Number of timesteps to execute.
//! \return Errorcode
int monte_carlo_propagation(struct Phase *const p, const unsigned int nsteps);

//! Monte-Carlo move of the configuration using the parallel iteration of polymers.
//!
//! \param p Initialized configuration.
//! \param nsteps \#steps to perform with the system.
//! \param tuning_parameter Parameter for ACC kernels. (vector_length)
//! \return Error code. Returns not equal to zero of an error occured.
int mc_polymer_iteration(struct Phase *const p, const unsigned int nsteps, const unsigned int tuning_parameter);

//! Trial Move generation for simple center of mass moves.
//!
//! \param p System description
//! \param poly_type Type of the polymer to move.
//! \param dx Output pointer to dx.
//! \param dy Output pointer to dy.
//! \param dz Output pointer to dz.
//! \param rng_state State of the random number generator.
#pragma acc routine(trial_move_cm) seq
void trial_move_cm(const struct Phase *p, const uint64_t poly_type, soma_scalar_t * const dx,
                   soma_scalar_t * const dy, soma_scalar_t * const dz, RNG_STATE * const rng_state);

//! Calculate the nonbonded energy of a given particle compared to a proposed move.
//!
//! \param p System description
//! \param monomer Original position of the proposed Monomer
//! \param dx Proposed move in X direction.
//! \param dy Proposed move in Y direction.
//! \param dz Proposed move in Z direction.
//! \param iwtype Type of the particle to move.
//! \return Calculated nonbonded energy difference. In case the of an error NAN is returned.
#pragma acc routine(calc_delta_nonbonded_energy) seq
soma_scalar_t calc_delta_nonbonded_energy(const struct Phase *p, const Monomer * const monomer,
                                          const soma_scalar_t dx, const soma_scalar_t dy,
                                          const soma_scalar_t dz, const unsigned int iwtype);

//! Monte-Carlo move: diffusion of the center of mass for entire molecules.
//!
//! \param p Initialized configuration.
//! \param nsteps \#steps to perform with the system.
//! \param tuning_parameter Parameter for ACC kernels. (vector_length)
//! \return Error code. Returns not equal to zero of an error occured.
int mc_center_mass(struct Phase *const p, const unsigned int nsteps, const unsigned int tuning_parameter);

//! Monte-Carlo move of the configuration using the parallel iteration
//! of polymers and independet sets inside each molecule.
//!
//! Via this two level parallelism a better saturation of GPUs can be achieved.
//! (More parallel set the same system.) In addition, the coalesence of memory access is optimized.
//! \param p Initialized configuration.
//! \param nsteps \#steps to perform with the system.
//! \param tuning_parameter Parameter for ACC kernels. (vector_length)
//! \return Error code. Returns not equal to zero of an error occured.
int mc_set_iteration(struct Phase *const p, const unsigned int nsteps, const unsigned int tuning_parameter);

/*! \brief calculate a trial move for the specified bead
  \param p Phase configuration, in which the trial move is proposed
  \param ipoly Local polymer id of the move
  \param ibead Monomer id of the move
  \param dx memory for the resulting x proposal.
  \param dy memory for the resulting x proposal.
  \param dz memory for the resulting x proposal.
  \param iwtype Type of the particle to move.
  \param rng_state State of the RNG.
*/
#pragma acc routine(trial_move) seq
void trial_move(const struct Phase *p, const uint64_t ipoly, const int ibead, soma_scalar_t * dx, soma_scalar_t * dy,
                soma_scalar_t * dz, const unsigned int iwtype, RNG_STATE * const rng_state);

//! Calculate the energy difference for a trial move.
//!
//! \param p configuration
//! \param ipoly polymer of moving bead
//! \param monomer Pointer to old monomer position.
//! \param ibead monomer of moving bead
//! \param dx proposed x move
//! \param dy proposed y move
//! \param dz proposed z move
//! \param iwtype Type of the monomer.
//! \return energy difference of proposed move. In case the of an error NAN is returned.
#pragma acc routine(calc_delta_energy) seq
soma_scalar_t calc_delta_energy(const struct Phase *p, const uint64_t ipoly, const Monomer * monomer,
                                const unsigned int ibead, const soma_scalar_t dx, const soma_scalar_t dy,
                                const soma_scalar_t dz, const unsigned int iwtype);

//! Calculate the bonded energy difference using the NEW2 bond
//! structure for a moved bead.
//!
//! \param p configuration
//! \param monomer Pointer to old Position of the Monomer.
//! \param ipoly polymer of moving bead
//! \param ibead monomer of moving bead
//! \param dx proposed x move
//! \param dy proposed y move
//! \param dz proposed z move
//! \return bonded energy difference of proposed move
#pragma acc routine(calc_delta_bonded_energy) seq
soma_scalar_t calc_delta_bonded_energy(const struct Phase *p, const Monomer * monomer,
                                       const uint64_t ipoly, const unsigned int ibead,
                                       const soma_scalar_t dx, const soma_scalar_t dy, const soma_scalar_t dz);

/*! \brief Calculation of the Metropolis acceptance criteria.
  Random number from rng.h -
  bool as in the macro in <stdbool.h> (valid for C99).

  \param delta_energy: double ( E_new - E_old ), energy change after the intended movement on the random bead
  \param p Phase construct of the simulated system
  \param rng State of the rng to use
  \param modifier Additional factor multiplied to the acceptance probability. Must be in \f$ [0,1] \f$.
  \return true or false according to the Metropolis criteria
*/
#pragma acc routine(som_accept) seq
int som_accept(RNG_STATE * const rng, const struct Phase *const p, const soma_scalar_t delta_energy,
               const soma_scalar_t modifier);

/*! \brief Smart Monte-Carlo (SMC) move.  Calculate the displacement and the energy change from the forces.

  \param p Initialized configuration.
  \param ipoly Polymer of moving bead
  \param ibead Monomer of moving bead
  \param dx proposed displacement calculated from the forces in x
  \param dy proposed displacement calculated from the forces in y
  \param dz proposed displacement calculated from the forces in z
  \param smc_deltaE energy change calculated from the forces acting on the bead before and after the proposed move.
  \param mybead Selected bead for the MC move
  \param myrngstate State of the random number generator, used to generate the random vector.
  \param iwtype Type of the monomer.
  \return displacement of the proposed move (in dx,dy,dz components) and energy change of the proposed SMC move.

  \note only bonded forces are considered.  Calculation of other forces
  (e.g. external fields) can be added in this function.
*/
#pragma acc routine(trial_move_smc) seq
void trial_move_smc(const struct Phase *p, const uint64_t ipoly, const int ibead,
                    soma_scalar_t * dx, soma_scalar_t * dy, soma_scalar_t * dz,
                    soma_scalar_t * smc_deltaE, const Monomer * mybead,
                    RNG_STATE * const myrngstate, const unsigned int iwtype);

#pragma acc routine(propose_smc_move) seq
void propose_smc_move(const struct Phase *p, const uint64_t ipoly, unsigned const int ibead, const unsigned int iwtype,
                      const soma_scalar_t x, const soma_scalar_t y, const soma_scalar_t z,
                      soma_scalar_t rx, soma_scalar_t ry, soma_scalar_t rz, soma_scalar_t * delta_E_bond,
                      soma_scalar_t * dx, soma_scalar_t * dy, soma_scalar_t * dz);

/*! \brief Calculate forces acting on a monomer resulting from all of its bonds.
  \param p Initialized configuration.
  \param ipoly Polymer of moving bead
  \param ibead Monomer of moving bead
  \param x coordinate
  \param y coordinate
  \param z coordinate
  \param fx x component of forces
  \param fy y component of forces
  \param fz z component of forces

  \return forces acting on the provided position (x,y,z) returned as: fx,fy,fz.

  \note coordinates (x,y,z) are passed as a parameter and not taken from
  the configuration to enable the computation of forces after and before the move with this same function.
*/
#pragma acc routine(add_bond_forces) seq
void add_bond_forces(const struct Phase *p, const uint64_t ipoly, unsigned const int ibead,
                     const soma_scalar_t x, const soma_scalar_t y, const soma_scalar_t z,
                     soma_scalar_t * fx, soma_scalar_t * fy, soma_scalar_t * fz);

//! Validates, whether a move for a particle from the old position, by
//! a displacement of dx violates the forbidden area51.

//! Such a violation can be caused by either the final position or by
//! passing through a forbidden area.  This is by sampling the move
//! path with a length of p->max_safe_jump.
//! \param p System description.
//! \param oldx Old position of particle X.
//! \param oldy Old position of particle Y.
//! \param oldz Old position of particle Z.
//! \param dx Displacement of the particle X.
//! \param dy Displacement of the particle Y.
//! \param dz Displacement of the particle Z.
//! \param nonexact Check is done nonexact. Set it to true, if you just want
//!  the end check, but better performance.
//! \warning Normal walls along the coordinate axis are correctly
//! filtered. But if the forbidden area is made up by boxes meeting
//! edge to edge particle can and will pass through.
//! \return True if move is allowed. False otherwise.
#pragma acc routine(possible_move_area51) seq
int possible_move_area51(const struct Phase *p, const soma_scalar_t oldx, const soma_scalar_t oldy,
                         const soma_scalar_t oldz, const soma_scalar_t dx, const soma_scalar_t dy,
                         const soma_scalar_t dz, const int nonexact);

/*! \brief Set iteration function used for one long chain, private function
  \param p Initialized configuration.
  \param nsteps \#steps to perform with the system.
  \param tuning_parameter Parameter for ACC kernels. (vector_length)
  \param nonexact_area51 The exact check of area51
  \param chain_i The index of the chain to be handled
  \return Error code. Returns either pgi error or domain error.
*/
int set_iteration_single_chain(struct Phase *const p, const unsigned int nsteps, const unsigned int tuning_parameter,
                               const int nonexact_area51, uint64_t chain_i);

/*! \brief Set iteration function for all the chains starting from start_chain, private function
  \param p Initialized configuration.
  \param nsteps \#steps to perform with the system.
  \param tuning_parameter Parameter for ACC kernels. (vector_length)
  \param nonexact_area51 The exact check of area51
  \param start_chain the starting index of the chains to be handled with this function
  \return Error code. Returns either pgi error or domain error.
*/
int set_iteration_multi_chain(struct Phase *const p, const unsigned int nsteps, const unsigned int tuning_parameter,
                              const int nonexact_area51, const int start_chain);

/*! \brief Private function used together with set_iteration_multi_chain and set_iteration_single_chain
  \param p Initialized configuration
  \param set_states The set states of the selected polymer
  \param beads ptr to the beads of the current polymer
  \param chain_index Index of the selected polymer
  \param iP Index of the current selected bead
  \param nonexact_area51 The exact check of area51
  \param ibead The selected bead
  \param iwtype The particle type of the selected particle
  \param *accepted_moves_set_ptr The pointer to the number of accepted_moves_set
  \return error_flags[0] indicating domain error
*/
#pragma acc routine(set_iteration_possible_move) seq
int set_iteration_possible_move(const struct Phase *p, RNG_STATE * const set_states, Monomer * const beads,
                                uint64_t chain_index, unsigned int iP, const int nonexact_area51,
                                const unsigned int ibead, const unsigned int iwtype,
                                unsigned int *accepted_moves_set_ptr);
#endif                          //SOMA_MC_H
