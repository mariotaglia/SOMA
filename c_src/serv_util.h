#ifndef _SERV_UTIL_H_
#define _SERV_UTIL_H_

#include <stdbool.h>
#include "phase.h"
#include <mpi.h>


bool need_to_do(unsigned int obs_delta, unsigned int time);

//! determines if no observables need to be calculated on this timestep.
//! dump does not count as an observable
//! \param ai information about the analysis
//! \param t timestep
//! \return true if no obs are to be calculated/written
bool no_observables(const Ana_Info * ai, unsigned int t);

//! determines if any of the observables of the timestep require Polymer-data to be calculated
//! \param ai information about the analysis
//! \param t timestep
//! \return true if polymers are required
bool has_poly_obs(const Ana_Info * ai, unsigned int t);

//! determines if the calculation of in this timestep
//! \param ai information about the analysis
//! \param t timestep
//! \return true if omega fields are required
bool has_omega_field_obs(const Ana_Info * ai, unsigned int t);

//! determines if any observables require the omega field to be calculated in this timestep
//! \param ai information about the analysis
//! \param t timestep
//! \return true if omega fields are required
bool has_field_obs(const Ana_Info * ai, unsigned int t);

#endif //_SERV_UTIL_H_
