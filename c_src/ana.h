//
// Created by julian on 30/03/2020.
//

#include <soma_config.h>
#include <stdint.h>
#include "phase.h"

#ifndef SOMA_ANA_REWORK_H
#define SOMA_ANA_REWORK_H


enum structure_factor_type { DYNAMICAL_STRUCTURE_FACTOR, STATIC_STRUCTURE_FACTOR };

enum obs {Re=0, dvar, Rg, b_anisotropy, acc_ratio, MSD, non_bonded_energy, bonded_energy,
    field, umbrella_field, dynamical_structure, static_structure, dump, /*to make enum iterable */ obs_it_end};

int analytics(struct Phase *const p);
void calc_Rg_local(const struct Phase *p, soma_scalar_t * const local_result, uint64_t * const local_counter);
void calc_Re_local(const struct Phase *p, soma_scalar_t * local_result, uint64_t * local_counter);
void calc_anisotropy_local(struct Phase *const p, soma_scalar_t *const local_result, uint64_t *const local_counter);
void calc_dvar_local(const struct Phase *p, soma_scalar_t * const local_result);
void calc_non_bonded_energy_local(const struct Phase * const p, soma_scalar_t * const non_bonded_energy);
void calc_msd_local(struct Phase *const p, soma_scalar_t *const msd, uint64_t *const counter);
void calc_bonded_energy_local(const struct Phase *p, soma_scalar_t *const bonded_energy);
int calc_structure_local(const struct Phase *p, soma_scalar_t * const result, const enum structure_factor_type sf_type);
void calc_acc_ratio_local(struct Phase *const p, uint64_t *accepts, uint64_t *moves);
int extent_ana_by_field(const soma_scalar_t * const data, const uint64_t n_data, const char *const name,
                        const hid_t file_id);
int extent_density_field(const struct Phase *const p, const void *const field_pointer, const char *const field_name,
                         hid_t hdf5_type, const MPI_Datatype mpi_type, const size_t data_size);
int extent_structure(const struct Phase *p, const soma_scalar_t * const data, const char *const name,
                     const hid_t file_id, const enum structure_factor_type sf_type);

#endif //SOMA_ANA_REWORK_H
