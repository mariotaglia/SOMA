//
// Created by julian on 04/04/2020.
//

#ifndef SOMA_ANA_MALLOC_H
#define SOMA_ANA_MALLOC_H

#include <stdint.h> // for uint64_t
#include <stdbool.h>
#include "phase.h"

// for the reduction we need of course two buffers: one for sending, and one for receiving.
// every rank writes into the sendbuffer, and after the in-place-reduction, the root will have the result
// in the sendbuffer, so only the sendbuffer (from now on just the buffer) is the only one that can be
// accessed through this interface

//getters for all the allocated memory (inside the buffer) of the observables
double * get_static_sf_val();
double * get_dynamical_sf_val();
double * get_bonded_energy_val();
double * get_non_bonded_energy_val();
uint64_t * get_moves();
uint64_t * get_accepts();
uint64_t * get_msd_counter();
double * get_msd_vals();
double * get_anisotropy_vals();
uint64_t * get_anisotropy_counter();
uint64_t * get_rg_counter();
double * get_rg_vals();
double * get_dvar_val();
uint64_t * get_re_counter();
double * get_re_vals();

//create and free the buffer,
//note that no pointer to the buffer is returned, instead, it's content is to be accessed with the getters above
int ana_malloc(const struct Phase *p, const bool *needToDo);
void ana_free();

// execute the reduction (in place,
// so that the getters in the root now point to the results of the reduction)
int ana_reduce(int root, MPI_Comm comm);

#endif //SOMA_ANA_MALLOC_H
