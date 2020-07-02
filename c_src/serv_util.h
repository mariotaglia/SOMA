#ifndef _SERV_UTIL_H_
#define _SERV_UTIL_H_

#include <stdbool.h>
#include "phase.h"

bool need_to_do(unsigned int obs_delta, unsigned int time);

bool no_observables(const Ana_Info * ai, unsigned int t);
bool has_poly_obs(const Ana_Info * ai, unsigned int t);
bool has_omega_field_obs(const Ana_Info * ai, unsigned int t);
bool has_field_obs(const Ana_Info * ai, unsigned int t);

#endif //_SERV_UTIL_H_
