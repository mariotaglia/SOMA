/* Copyright (C) 2016-2021 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016-2017 Marcel Langenberg
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
#include "mc.h"
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "phase.h"
#include "mesh.h"
#include "independent_sets.h"
#include "mobility.h"

void trial_move(const Phase * p, const uint64_t ipoly, const int ibead,
                soma_scalar_t * dx, soma_scalar_t * dy, soma_scalar_t * dz, const unsigned int iwtype,
                RNG_STATE * const rng_state)
{
    //Just to shut up the compiler warning:
    //Any decent compiler optimize it out
    soma_scalar_t scale = ibead + 0 * ipoly;
    scale = p->A[iwtype];

    *dx = scale * (soma_rng_soma_scalar(rng_state, p) - 0.5);
    *dy = scale * (soma_rng_soma_scalar(rng_state, p) - 0.5);
    *dz = scale * (soma_rng_soma_scalar(rng_state, p) - 0.5);
}

void trial_move_cm(const Phase * p, const uint64_t poly_type, soma_scalar_t * const dx, soma_scalar_t * const dy,
                   soma_scalar_t * const dz, RNG_STATE * const rng_state)
{
#ifndef _OPENACC
    assert(p->cm_a);
#endif                          //_OPENACC
    const soma_scalar_t scale = p->cm_a[poly_type];

    *dx = scale * (soma_rng_soma_scalar(rng_state, p) - 0.5);
    *dy = scale * (soma_rng_soma_scalar(rng_state, p) - 0.5);
    *dz = scale * (soma_rng_soma_scalar(rng_state, p) - 0.5);
}

int som_accept(RNG_STATE * const rng, const Phase * const p, const soma_scalar_t delta_energy,
               const soma_scalar_t modifier)
{
    // \todo kBT reqired
    const soma_scalar_t p_acc = exp(-1.0 * delta_energy) * modifier;

    //Use lazy eval.
    if ((p_acc > 1) || (p_acc > soma_rng_soma_scalar(rng, p)))
        {
            return true;
        }
    else
        {
            return false;
        }

}

soma_scalar_t calc_delta_nonbonded_energy(const Phase * p, const Monomer * monomer,
                                          const soma_scalar_t dx, const soma_scalar_t dy, const soma_scalar_t dz,
                                          const unsigned int iwtype)
{
    // Old non-bonded interaction
    const soma_scalar_t xold = monomer->x;
    const soma_scalar_t yold = monomer->y;
    const soma_scalar_t zold = monomer->z;
    const uint64_t cellindex_old = coord_to_index_unified(p, xold, yold, zold, iwtype);
    const uint64_t cellindex_new = coord_to_index_unified(p, xold + dx, yold + dy, zold + dz, iwtype);
    if (cellindex_old > p->n_cells_local * p->n_types || cellindex_new > p->n_cells_local * p->n_types)
        {
#ifdef NAN
#ifdef __PGI
#if __PGIC__ < 20
            return 0 / 0;
#endif                          //__PGIC
#endif                          //_PGI
            return NAN;
#else
            return nan("");
#endif                          //NAN
        }
    const soma_scalar_t energy_old = p->omega_field_unified[cellindex_old];
    // New non-bonded interaction
    const soma_scalar_t energy_new = p->omega_field_unified[cellindex_new];
    const soma_scalar_t energy = energy_new - energy_old;
    return energy;
}

soma_scalar_t calc_delta_energy(const Phase * p, const uint64_t ipoly, const Monomer * const monomer,
                                const unsigned int ibead, const soma_scalar_t dx, const soma_scalar_t dy,
                                const soma_scalar_t dz, const unsigned int iwtype)
{
    const soma_scalar_t delta_bonded_energy = calc_delta_bonded_energy(p, monomer, ipoly, ibead, dx, dy, dz);
    const soma_scalar_t delta_nonbonded_energy = calc_delta_nonbonded_energy(p, monomer, dx, dy, dz, iwtype);

    // non-bonded energy + bonded energy
    soma_scalar_t energy = delta_nonbonded_energy;
    energy += delta_bonded_energy;
    return energy;
}

soma_scalar_t calc_delta_bonded_energy(const Phase * p, const Monomer * monomer,
                                       const uint64_t ipoly, const unsigned int ibead,
                                       const soma_scalar_t dx, const soma_scalar_t dy, const soma_scalar_t dz)
{
    soma_scalar_t delta_energy = 0;
    // loop over bonds of this bead
    const int start = get_bondlist_offset(p->poly_arch[p->poly_type_offset[p->polymers[ipoly].type] + ibead + 1]);
    Monomer *beads = p->ph.beads.ptr;
    beads += p->polymers[ipoly].bead_offset;

    if (start > 0)
        {
            unsigned int end = 0;
            for (int i = start; end == 0; i++)
                {
                    const uint32_t info = p->poly_arch[i];
                    end = get_end(info);
                    const unsigned int bond_type = get_bond_type(info);
                    const int offset = get_offset(info);

                    const int neighbour_id = ibead + offset;
                    const unsigned int jbead = neighbour_id;
                    //printf("    offset=%d jbead=%u  end=%u type=%u\n",neigh->offset,jbead,end,neigh->bond_type);

                    soma_scalar_t scale = 1.;
                    switch (bond_type)
                        {
                        case HARMONICVARIABLESCALE:
                            scale = p->harmonic_normb_variable_scale;
                            /* intentionally falls through */
                        case HARMONIC:
                            //Empty statement, because a statement after a label
                            //has to come before any declaration
                            ;
                            const soma_scalar_t old_rx = calc_bond_length(monomer->x, beads[jbead].x, p->Lx,
                                                                          p->args.bond_minimum_image_convention_flag);
                            const soma_scalar_t new_rx = old_rx + dx;
                            const soma_scalar_t old_ry = calc_bond_length(monomer->y, beads[jbead].y, p->Ly,
                                                                          p->args.bond_minimum_image_convention_flag);
                            const soma_scalar_t new_ry = old_ry + dy;
                            const soma_scalar_t old_rz = calc_bond_length(monomer->z, beads[jbead].z, p->Lz,
                                                                          p->args.bond_minimum_image_convention_flag);
                            const soma_scalar_t new_rz = old_rz + dz;

                            const soma_scalar_t old_r2 = old_rx * old_rx + old_ry * old_ry + old_rz * old_rz;
                            const soma_scalar_t new_r2 = new_rx * new_rx + new_ry * new_ry + new_rz * new_rz;
                            delta_energy += p->harmonic_normb * (new_r2 - old_r2) * scale;
                            break;

                        case STIFF:
#ifndef _OPENACC
                            fprintf(stderr, "ERROR: %s:%d stiff bond not yet implemented.\n", __FILE__, __LINE__);
#endif                          //_OPENACC
                            break;

                        default:
#ifndef _OPENACC
                            fprintf(stderr, "ERROR: %s:%d unknow bond type appeared %d\n",
                                    __FILE__, __LINE__, bond_type);
#endif                          //OPENACC
                            break;
                        }

                }
        }
    return delta_energy;
}

int monte_carlo_propagation(Phase * const p, unsigned int nsteps)
{
    //Update the omega fields for the calculations.
    update_omega_fields(p);
    int ret;
    start_autotuner(&(p->mc_autotuner));
    switch (p->args.iteration_alg_arg)
        {
        case iteration_alg_arg_POLYMER:

            ret = mc_polymer_iteration(p, nsteps, p->mc_autotuner.value);
            break;
        case iteration_alg_arg_SET:
            ret = mc_set_iteration(p, nsteps, p->mc_autotuner.value);
            break;
        case iteration_alg__NULL:
        default:
            fprintf(stderr, "ERROR: Unknown iteration algorithm selected.\n");
            ret = 1;
        }
    end_autotuner(&(p->mc_autotuner));
    if (ret != 0)
        return ret;

    if (p->cm_a)
        {
            start_autotuner(&(p->cm_mc_autotuner));
            ret = mc_center_mass(p, 1, p->cm_mc_autotuner.value);
            end_autotuner(&(p->cm_mc_autotuner));
        }
    if (p->args.autotuner_restart_period_arg > 0 && p->time % p->args.autotuner_restart_period_arg == 0)
        {
            restart_autotuner(&(p->mc_autotuner));
            restart_autotuner(&(p->cm_mc_autotuner));
        }
    update_density_fields(p);
    return ret;
}

int mc_center_mass(Phase * const p, const unsigned int nsteps, const unsigned int tuning_parameter)
{
    assert(p->cm_a);
    //Shutup compiler warning
    unsigned int step = tuning_parameter;
    step = 0;

    int error_flags[1] = { 0 }; //error_flag[0] indicates domain errors
#pragma acc enter data copyin(error_flags[0:1])

    // Loop over the MC scweeps
    for (step = 0; step < nsteps; step++)
        {
            uint64_t n_polymers = p->n_polymers;
            unsigned int n_accepts = 0;
//#pragma acc parallel loop vector_length(tuning_parameter) reduction(+:n_accepts)
#pragma acc parallel loop vector_length(tuning_parameter) present(p[0:1])
#pragma omp parallel for reduction(+:n_accepts)
            for (uint64_t npoly = 0; npoly < n_polymers; npoly++)
                {
                    Polymer *mypoly = &p->polymers[npoly];
                    unsigned int myN = p->poly_arch[p->poly_type_offset[mypoly->type]];
                    RNG_STATE rng_state_local = mypoly->poly_state;
                    RNG_STATE *myrngstate = &rng_state_local;
                    const unsigned int poly_type = mypoly->type;
                    Monomer *beads = p->ph.beads.ptr;
                    beads += mypoly->bead_offset;

                    if (p->cm_a[poly_type] > 0)
                        {
                            //Generate a random displacement for the center of mass.
                            soma_scalar_t dx, dy, dz;
                            trial_move_cm(p, poly_type, &dx, &dy, &dz, myrngstate);

                            soma_scalar_t delta_energy = 0;
                            int move_allowed = 1;
                            soma_scalar_t pacc_modifier = 1;

                            //#pragma acc loop vector reduction(+:delta_energy) reduction(|:move_allowed)
                            //Unfortunately this has to be a seq loop, because the reduction crashes.
#pragma acc loop vector seq
                            for (unsigned int ibead = 0; ibead < myN; ibead++)
                                {
                                    const Monomer mybead = beads[ibead];
                                    const unsigned int iwtype = get_particle_type(p, npoly, ibead);

                                    const int tmp = possible_move_area51(p, mybead.x, mybead.y, mybead.z, dx, dy, dz,
                                                                         p->args.nonexact_area51_flag);
                                    move_allowed &= tmp;

                                    if (tmp)
                                        {
                                            // calculate energy change
                                            delta_energy +=
                                                calc_delta_energy(p, npoly, &mybead, ibead, dx, dy, dz, iwtype);
                                            pacc_modifier *=
                                                get_mobility_modifier(p, iwtype, mybead.x, mybead.y, mybead.z);
                                            pacc_modifier *=
                                                get_mobility_modifier(p, iwtype, mybead.x + dx, mybead.y + dy,
                                                                      mybead.z + dz);
                                        }
                                }
                            if (delta_energy != delta_energy)   // isnan(delta_energy) ) not working with PGI OpanACC
                                {
                                    error_flags[0] = npoly + 1;
                                    move_allowed = 0;
                                }

                            //Take the sqaure root of mobility modifier m (collectively for all beads)
                            pacc_modifier = sqrt(pacc_modifier);

                            //Accept Monte-Carlo call
                            if (move_allowed && som_accept(myrngstate, p, delta_energy, pacc_modifier) == 1)
                                {
#ifndef _OPENACC
                                    n_accepts += 1;
#endif                          //_OPENACC

//#pragma acc loop vector
                                    //See above
#pragma acc loop seq
                                    for (unsigned int ibead = 0; ibead < myN; ibead++)
                                        {
                                            Monomer mybead = beads[ibead];
                                            Monomer *const mybead_ptr = &(beads[ibead]);

                                            mybead.x += dx;
                                            mybead.y += dy;
                                            mybead.z += dz;
                                            *mybead_ptr = mybead;
                                        }
                                }

                            //Copy back the modified RNG state.
                            mypoly->poly_state = rng_state_local;
                        }
                }
        }

#pragma acc exit data copyout(error_flags[0:1])
    if (error_flags[0] != 0)
        {
            fprintf(stderr, "ERROR: Domain error. %d"
                    " A particle has left the buffer domain."
                    " Restart your simulation with larger buffers. %s:%d\n", error_flags[0], __FILE__, __LINE__);
            return error_flags[0];
        }
    return 0;
}

int mc_polymer_iteration(Phase * const p, const unsigned int nsteps, const unsigned int tuning_parameter)
{
    //Shutup compiler warning
    unsigned int step = tuning_parameter;
    step = 0;

    int error_flags[1] = { 0 }; //error_flag[0] indicates domain errors
#pragma acc enter data copyin(error_flags[0:1])
    // Loop over the MC scweeps
    for (step = 0; step < nsteps; step++)
        {
            uint64_t n_polymers = p->n_polymers;
            unsigned int n_accepts = 0;
            const unsigned int gpu_time = p->time;

            //#pragma acc parallel loop vector_length(tuning_parameter) reduction(+:n_accepts)
#pragma acc parallel loop vector_length(tuning_parameter) present(p[0:1])
#pragma omp parallel for reduction(+:n_accepts)
            for (uint64_t npoly = 0; npoly < n_polymers; npoly++)
                {
                    unsigned int accepted_moves_loc = 0;

                    // Rebuild bond information for this chain from bonds, or stay with linear right now?
                    Polymer *mypoly = &p->polymers[npoly];
                    if (gpu_time % p->mobility.poly_type_mc_freq[mypoly->type] != 0)
                        continue;       //EARLY LOOP EXIT FOR MOBILITY CONTRAST

                    Monomer *beads = p->ph.beads.ptr;
                    beads += mypoly->bead_offset;

                    const int mypoly_poly_type_offset = p->poly_type_offset[mypoly->type];
                    unsigned int myN = p->poly_arch[mypoly_poly_type_offset];
                    RNG_STATE *myrngstate = &mypoly->poly_state;        // maybe local copy of rngstate

                    // MC sweep for this chain
#pragma acc loop seq
                    for (unsigned int nmc = 0; nmc < myN; nmc++)
                        {

                            soma_scalar_t dx = 0, dy = 0, dz = 0, delta_energy = 0;
                            unsigned int ibead;

                            // pick a random bead.
                            ibead = soma_rng_uint(myrngstate, p) % myN;
                            const unsigned int iwtype = get_particle_type(p, npoly, ibead);

                            Monomer mybead = beads[ibead];
                            Monomer *mybead_ptr = &(beads[ibead]);

                            // roll normal MC trial move or force biased MC move
                            soma_scalar_t delta_E_bond = 0.0;
                            switch (p->args.move_type_arg)
                                {
                                case move_type_arg_TRIAL:
                                    trial_move(p, npoly, ibead, &dx, &dy, &dz, iwtype, myrngstate);     // normal MC move
                                    delta_E_bond = calc_delta_bonded_energy(p, &mybead, npoly, ibead, dx, dy, dz);
                                    break;
                                case move_type_arg_SMART:
                                    trial_move_smc(p, npoly, ibead, &dx, &dy, &dz, &delta_E_bond, &mybead, myrngstate, iwtype); // force biased move
                                    break;
                                case move_type__NULL:
                                default:
                                    delta_E_bond = 0.0;
                                    break;
                                }

                            const int move_allowed = possible_move_area51(p, mybead.x, mybead.y, mybead.z, dx, dy, dz,
                                                                          p->args.nonexact_area51_flag);
                            if (move_allowed)
                                {
                                    delta_energy = delta_E_bond;
                                    delta_energy += calc_delta_nonbonded_energy(p, &mybead, dx, dy, dz, iwtype);
                                    if (delta_energy != delta_energy)   // isnan(delta_energy) ) not working with PGI OpenACC
                                        {
                                            error_flags[0] = npoly + 1;
                                        }
                                    const soma_scalar_t newx = mybead.x + dx;
                                    const soma_scalar_t newy = mybead.y + dy;
                                    const soma_scalar_t newz = mybead.z + dz;

                                    const soma_scalar_t old_mobility_modifier =
                                        get_mobility_modifier(p, iwtype, mybead.x, mybead.y, mybead.z);
                                    const soma_scalar_t new_mobility_modifier =
                                        get_mobility_modifier(p, iwtype, newx, newy, newz);
                                    const soma_scalar_t mobility_modifier =
                                        sqrt(old_mobility_modifier * new_mobility_modifier);

                                    // MC roll to accept / reject
                                    if (som_accept(myrngstate, p, delta_energy, mobility_modifier) == 1)
                                        {
                                            mybead_ptr->x = newx;
                                            mybead_ptr->y = newy;
                                            mybead_ptr->z = newz;
#ifndef _OPENACC
                                            accepted_moves_loc += 1;
#endif                          //_OPENACC
                                        }
                                }
                        }
#ifndef _OPENACC
                    n_accepts += accepted_moves_loc;
#endif                          //_OPENACC
                }

            p->time += 1;
            p->n_moves += p->num_all_beads_local;
            p->n_accepts += n_accepts;
        }

#pragma acc exit data copyout(error_flags[0:1])
    if (error_flags[0] != 0)
        {
            fprintf(stderr, "ERROR: Domain error. %d"
                    " A particle has left the buffer domain."
                    " Restart your simulation with larger buffers. %s:%d\n", error_flags[0], __FILE__, __LINE__);
            return error_flags[0];
        }
    return 0;
}

int set_iteration_multi_chain(Phase * const p, const unsigned int nsteps, const unsigned int tuning_parameter,
                              const int nonexact_area51, const int start_chain)
{
    int error_flags[2] = { 0 }; // [0] domain error, [1] pgi_bug
#pragma acc enter data copyin(error_flags[0:2])
    for (unsigned int step = 0; step < nsteps; step++)
        {
            const uint64_t n_polymers = p->n_polymers;

            const unsigned int gpu_time = p->time;
            //Shutup compiler warning
            unsigned int n_accepts = tuning_parameter;
            n_accepts = 0;
#pragma acc parallel loop vector_length(tuning_parameter) present(p[0:1]) async
#pragma omp parallel for reduction(+:n_accepts)
            for (uint64_t npoly = start_chain; npoly < n_polymers; npoly++)
                {
                    unsigned int accepted_moves_poly = 0;
                    Polymer *const mypoly = &p->polymers[npoly];

                    const unsigned int poly_type = mypoly->type;
                    if (gpu_time % p->mobility.poly_type_mc_freq[poly_type] != 0)
                        continue;       //EARLY LOOP EXIT FOR MOBILITY CONTRAST

                    //const int mypoly_poly_type_offset = p->poly_type_offset[poly_type];
                    const IndependetSets mySets = p->sets[poly_type];

                    const unsigned int n_sets = mySets.n_sets;
                    const unsigned int *const set_length = mySets.set_length;
                    const unsigned int *const sets = mySets.sets;
                    const unsigned int max_member = mySets.max_member;
                    RNG_STATE *set_states = p->ph.set_states.ptr;
                    set_states += mypoly->set_states_offset;

                    unsigned int *set_permutation = p->ph.set_permutation.ptr;
                    set_permutation += mypoly->set_permutation_offset;

                    Monomer *beads = p->ph.beads.ptr;
                    beads += mypoly->bead_offset;

                    //Generate random permutation of the sets
                    //http://www.wikipedia.or.ke/index.php/Permutation
#pragma acc loop seq
                    for (unsigned int i = 0; i < n_sets; i++)
                        {
                            const unsigned int d = soma_rng_uint(&(mypoly->poly_state), p) % (i + 1);
                            set_permutation[i] = set_permutation[d];
                            set_permutation[d] = i;
                        }

#pragma acc loop seq
                    for (unsigned int iSet = 0; iSet < n_sets; iSet++)
                        {
                            unsigned int accepted_moves_set = 0;
                            const unsigned int set_id = set_permutation[iSet];
                            const unsigned int len = set_length[set_id];
#pragma acc loop vector
                            for (unsigned int iP = 0; iP < len; iP++)
                                {
                                    const unsigned int ibead = sets[set_id * max_member + iP];
                                    const unsigned int iwtype = get_particle_type(p, npoly, ibead);
                                    int error_0 = set_iteration_possible_move(p, set_states, beads, npoly, iP,
                                                                              nonexact_area51, ibead, iwtype,
                                                                              &accepted_moves_set);
                                    error_flags[0] = error_0;
                                }

#ifndef _OPENACC
                            accepted_moves_poly += accepted_moves_set;
#endif                          //_OPENACC
                        }
#ifndef _OPENACC
                    n_accepts += accepted_moves_poly;
#endif                          //_OPENACC
                }
            //p->time += 1;
            p->n_moves += p->num_all_beads_local;
#ifndef _OPENACC
            p->n_accepts += n_accepts;
#endif                          //_OPENACC
        }
    int ret = 0;
#pragma acc exit data copyout(error_flags[0:2])
    if (error_flags[0] != 0)
        {
            fprintf(stderr, "ERROR: Domain error. %d"
                    " A particle has left the buffer domain."
                    " Restart your simulation with larger buffers. %s:%d\n", error_flags[0], __FILE__, __LINE__);
            return error_flags[0];
        }
    ret = error_flags[1];
    return ret;

}

int set_iteration_single_chain(Phase * const p, const unsigned int nsteps, const unsigned int tuning_parameter,
                               const int nonexact_area51, uint64_t chain_i)
{
    int error_flags[2] = { 0 }; // [0] domain error, [1] pgi_bug

#ifndef _OPENACC
    unsigned int accepted_moves_poly = 0;
#endif                          //_OPENACC

    const unsigned int gpu_time = p->time;

#pragma acc enter data copyin(error_flags[0:2])
    for (unsigned int step = 0; step < nsteps; step++)
        {
            unsigned int n_accepts = 0;

            Polymer *const mypoly = &p->polymers[chain_i];

            const unsigned int poly_type = mypoly->type;
            if (gpu_time % p->mobility.poly_type_mc_freq[poly_type] != 0)
                continue;       //EARLY LOOP EXIT FOR MOBILITY CONTRAST

            //const int mypoly_poly_type_offset = p->poly_type_offset[poly_type];
            const IndependetSets mySets = p->sets[poly_type];

            const unsigned int n_sets = mySets.n_sets;
            const unsigned int *const set_length = mySets.set_length;
            const unsigned int *const sets = mySets.sets;
            const unsigned int max_member = mySets.max_member;
            RNG_STATE *set_states = p->ph.set_states.ptr;
            set_states += mypoly->set_states_offset;

            unsigned int *set_permutation = p->ph.set_permutation.ptr;
            set_permutation += mypoly->set_permutation_offset;

            Monomer *beads = p->ph.beads.ptr;
            beads += mypoly->bead_offset;
            //Generate random permutation of the sets
            //http://www.wikipedia.or.ke/index.php/Permutation

            for (unsigned int i = 0; i < n_sets; i++)
                {
                    const unsigned int d = soma_rng_uint(&((&p->polymers[chain_i])->poly_state), p) % (i + 1);
                    set_permutation[i] = set_permutation[d];
                    set_permutation[d] = i;
                }

            for (unsigned int iSet = 0; iSet < n_sets; iSet++)
                {
                    unsigned int accepted_moves_set = 0;
                    const unsigned int set_id = set_permutation[iSet];
                    const unsigned int len = set_length[set_id];

#pragma acc parallel loop vector_length(tuning_parameter) present(p[0:1]) async
#pragma omp parallel for reduction(+:accepted_moves_set)
                    for (unsigned int iP = 0; iP < len; iP++)
                        {
                            const unsigned int ibead = sets[set_id * max_member + iP];
                            const unsigned int iwtype = get_particle_type(p, chain_i, ibead);
                            int error_0 = set_iteration_possible_move(p, set_states, beads, chain_i, iP,
                                                                      nonexact_area51, ibead, iwtype,
                                                                      &accepted_moves_set);
                            error_flags[0] = error_0;
                        }

#ifndef _OPENACC
                    accepted_moves_poly += accepted_moves_set;
#endif                          //_OPENACC
                }
#ifndef _OPENACC
            n_accepts += accepted_moves_poly;
#endif                          //_OPENACC
            p->n_moves += p->num_all_beads_local;
#ifndef _OPENACC
            p->n_accepts += n_accepts;
#endif                          //_OPENACC
        }
    int ret = 0 * tuning_parameter;     //Shutup compiler warning
#pragma acc exit data copyout(error_flags[0:2])
    if (error_flags[0] != 0)
        {
            fprintf(stderr, "ERROR: Domain error. %d"
                    " A particle has left the buffer domain."
                    " Restart your simulation with larger buffers. %s:%d\n", error_flags[0], __FILE__, __LINE__);
            return error_flags[0];
        }

    ret = error_flags[1];
    return ret;
}

int mc_set_iteration(Phase * const p, const unsigned int nsteps, const unsigned int tuning_parameter)
{
    const int nonexact_area51 = p->args.nonexact_area51_flag + 0 * tuning_parameter;    //&Shutup compiler warning.

    //test the order of the polymers and reorder the polymers according to their length if needed
    if (p->time % p->args.set_order_frequency_arg == 0)
        {
            int order = 0;
            uint32_t *const length_poly = (uint32_t * const)malloc(p->n_polymers * sizeof(uint32_t));
            if (length_poly == NULL)
                {
                    fprintf(stderr, "ERROR: %d: %s:%d Malloc error\n", p->info_MPI.world_rank, __FILE__, __LINE__);
                    return -1;
                }
            for (uint64_t poly_i = 0; poly_i < p->n_polymers; poly_i++)
                {
                    Polymer *const this_poly = &p->polymers[poly_i];
                    const unsigned int poly_type = this_poly->type;
                    length_poly[poly_i] = p->poly_arch[p->poly_type_offset[poly_type]];
                }

            for (unsigned int poly_n = p->n_polymers; poly_n > 1; poly_n--)
                {
                    for (unsigned int poly_i = 0; poly_i < poly_n - 1; poly_i++)
                        {
                            if (length_poly[poly_i] > length_poly[poly_i + 1])
                                {
                                    int tmp = length_poly[poly_i];
                                    length_poly[poly_i] = length_poly[poly_i + 1];
                                    length_poly[poly_i + 1] = tmp;
                                    order = 1;
                                    break;
                                }
                        }
                    if (order == 1)
                        break;
                }

            if (order == 1)
                {
                    p->num_long_chain = mc_set_init(p);
                }
            free(length_poly);
        }

    int ret;
    unsigned int num_long_chain = (unsigned int)p->num_long_chain;
    for (unsigned int index = 0; index < num_long_chain; index++)
        {
            ret = set_iteration_single_chain(p, nsteps, tuning_parameter, nonexact_area51, index);
        }
    if (num_long_chain != p->n_polymers)
        {
            ret = set_iteration_multi_chain(p, nsteps, tuning_parameter, nonexact_area51, num_long_chain);
        }
    p->time += 1;
#pragma acc wait
    ret = 0;
    return ret;
}

void trial_move_smc(const Phase * p, const uint64_t ipoly, const int ibead, soma_scalar_t * dx, soma_scalar_t * dy,
                    soma_scalar_t * dz, soma_scalar_t * delta_E_bond, const Monomer * mybead,
                    RNG_STATE * const myrngstate, const unsigned int iwtype)
{
    soma_scalar_t rx, ry, rz;
    soma_normal_vector(myrngstate, p, &rx, &ry, &rz);
    soma_scalar_t x = mybead->x;
    soma_scalar_t y = mybead->y;
    soma_scalar_t z = mybead->z;
    soma_scalar_t E_bond = 0;
    propose_smc_move(p, ipoly, ibead, iwtype, x, y, z, rx, ry, rz, &E_bond, dx, dy, dz);
    *delta_E_bond = E_bond;

}

void propose_smc_move(const Phase * p, const uint64_t ipoly, unsigned const int ibead, const unsigned int iwtype,
                      const soma_scalar_t x, const soma_scalar_t y, const soma_scalar_t z,
                      soma_scalar_t rx, soma_scalar_t ry, soma_scalar_t rz, soma_scalar_t * delta_E_bond,
                      soma_scalar_t * dx, soma_scalar_t * dy, soma_scalar_t * dz)
{
    Monomer *beads = p->ph.beads.ptr;
    beads += p->polymers[ipoly].bead_offset;
    const int start = get_bondlist_offset(p->poly_arch[p->poly_type_offset[p->polymers[ipoly].type] + ibead + 1]);
    soma_scalar_t old_E_B = 0;
    soma_scalar_t new_E_B = 0;
    soma_scalar_t fx = 0.0, fy = 0.0, fz = 0.0;
    soma_scalar_t nfx = 0.0, nfy = 0.0, nfz = 0.0;

    /* Calculate bonded forces & energy before move */
    if (start > 0)
        {
            unsigned int end = 0;
            for (int i = start; end == 0; i++)
                {
                    const uint32_t info = p->poly_arch[i];
                    end = get_end(info);
                    const unsigned int bond_type = get_bond_type(info);
                    const int offset = get_offset(info);

                    const int neighbour_id = ibead + offset;
                    const unsigned int jbead = neighbour_id;
                    soma_scalar_t scale = 1;
                    switch (bond_type)
                        {
                        case HARMONICVARIABLESCALE:
                            scale = p->harmonic_normb_variable_scale;
                            /* intentionally falls through */
                        case HARMONIC:
                            //Empty statement, because a statement after a label
                            //has to come before any declaration
                            ;

                            soma_scalar_t bx_tmp = calc_bond_length(beads[jbead].x, x, p->Lx,
                                                                    p->args.bond_minimum_image_convention_flag);
                            soma_scalar_t by_tmp = calc_bond_length(beads[jbead].y, y, p->Ly,
                                                                    p->args.bond_minimum_image_convention_flag);
                            soma_scalar_t bz_tmp = calc_bond_length(beads[jbead].z, z, p->Lz,
                                                                    p->args.bond_minimum_image_convention_flag);
                            fx += bx_tmp * 2.0 * p->harmonic_normb * scale;
                            fy += by_tmp * 2.0 * p->harmonic_normb * scale;
                            fz += bz_tmp * 2.0 * p->harmonic_normb * scale;
                            soma_scalar_t old_r2_tmp = bx_tmp * bx_tmp + by_tmp * by_tmp + bz_tmp * bz_tmp;
                            old_E_B += p->harmonic_normb * (old_r2_tmp) * scale;

                            break;

                        case STIFF:
#ifndef _OPENACC
                            fprintf(stderr, "ERROR: %s:%d stiff bond not yet implemented.\n", __FILE__, __LINE__);
#endif                          //OPENACC
                            break;

                        default:
#ifndef _OPENACC
                            fprintf(stderr, "ERROR: %s:%d unknow bond type appeared %d\n",
                                    __FILE__, __LINE__, bond_type);
#endif                          //OPENACC
                            break;
                        }
                }
        }

    /* Propose new move */
    const soma_scalar_t R = p->R[iwtype];
    const soma_scalar_t A = p->A[iwtype];
    *dx = A * fx + R * rx;
    *dy = A * fy + R * ry;
    *dz = A * fz + R * rz;

    /* Calculate bonded forces & energy after move */
    if (start > 0)
        {
            unsigned int end = 0;
            for (int i = start; end == 0; i++)
                {
                    const uint32_t info = p->poly_arch[i];
                    end = get_end(info);
                    const unsigned int bond_type = get_bond_type(info);
                    const int offset = get_offset(info);
                    const int neighbour_id = ibead + offset;
                    const unsigned int jbead = neighbour_id;
                    soma_scalar_t scale = 1;
                    switch (bond_type)
                        {
                        case HARMONICVARIABLESCALE:
                            scale = p->harmonic_normb_variable_scale;
                            /* intentionally falls through */
                        case HARMONIC:
                            //Empty statement, because a statement after a label
                            //has to come before any declaration
                            ;

                            soma_scalar_t bx_tmp = calc_bond_length(beads[jbead].x, x + *dx, p->Lx,
                                                                    p->args.bond_minimum_image_convention_flag);
                            soma_scalar_t by_tmp = calc_bond_length(beads[jbead].y, y + *dy, p->Ly,
                                                                    p->args.bond_minimum_image_convention_flag);
                            soma_scalar_t bz_tmp = calc_bond_length(beads[jbead].z, z + *dz, p->Lz,
                                                                    p->args.bond_minimum_image_convention_flag);
                            nfx += bx_tmp * 2.0 * p->harmonic_normb * scale;
                            nfy += by_tmp * 2.0 * p->harmonic_normb * scale;
                            nfz += bz_tmp * 2.0 * p->harmonic_normb * scale;
                            soma_scalar_t new_r2_tmp = bx_tmp * bx_tmp + by_tmp * by_tmp + bz_tmp * bz_tmp;
                            new_E_B += p->harmonic_normb * (new_r2_tmp) * scale;

                            break;

                        case STIFF:
#ifndef _OPENACC
                            fprintf(stderr, "ERROR: %s:%d stiff bond not yet implemented.\n", __FILE__, __LINE__);
#endif                          //OPENACC
                            break;

                        default:
#ifndef _OPENACC
                            fprintf(stderr, "ERROR: %s:%d unknow bond type appeared %d\n",
                                    __FILE__, __LINE__, bond_type);
#endif                          //OPENACC
                            break;
                        }
                }
        }

    *delta_E_bond = new_E_B - old_E_B;
    soma_scalar_t delta_E_SMC = 0.5 * ((nfx + fx) * (*dx) + (nfy + fy) * (*dy) + (nfz + fz) * (*dz)) +
        0.25 * A * ((nfx * nfx) + (nfy * nfy) + (nfz * nfz) - (fx * fx) - (fy * fy) - (fz * fz));
    *delta_E_bond += delta_E_SMC;

}

void add_bond_forces(const Phase * p, const uint64_t ipoly, unsigned const int ibead,
                     const soma_scalar_t x, const soma_scalar_t y, const soma_scalar_t z,
                     soma_scalar_t * fx, soma_scalar_t * fy, soma_scalar_t * fz)
{
    soma_scalar_t v1x = 0.0, v1y = 0.0, v1z = 0.0;
    Monomer *beads = p->ph.beads.ptr;
    beads += p->polymers[ipoly].bead_offset;

    const int start = get_bondlist_offset(p->poly_arch[p->poly_type_offset[p->polymers[ipoly].type] + ibead + 1]);

    if (start > 0)
        {
            unsigned int end = 0;
            for (int i = start; end == 0; i++)
                {
                    const uint32_t info = p->poly_arch[i];
                    end = get_end(info);
                    const unsigned int bond_type = get_bond_type(info);
                    const int offset = get_offset(info);

                    const int neighbour_id = ibead + offset;
                    const unsigned int jbead = neighbour_id;

                    soma_scalar_t scale = 1;
                    switch (bond_type)
                        {
                        case HARMONICVARIABLESCALE:
                            scale = p->harmonic_normb_variable_scale;
                            /* intentionally falls through */
                        case HARMONIC:
                            //Empty statement, because a statement after a label
                            //has to come before any declaration
                            ;
                            soma_scalar_t v1x_tmp = calc_bond_length(beads[jbead].x, x, p->Lx,
                                                                     p->args.bond_minimum_image_convention_flag);
                            soma_scalar_t v1y_tmp = calc_bond_length(beads[jbead].y, y, p->Ly,
                                                                     p->args.bond_minimum_image_convention_flag);
                            soma_scalar_t v1z_tmp = calc_bond_length(beads[jbead].z, z, p->Lz,
                                                                     p->args.bond_minimum_image_convention_flag);
                            v1x += v1x_tmp * 2.0 * p->harmonic_normb * scale;
                            v1y += v1y_tmp * 2.0 * p->harmonic_normb * scale;
                            v1z += v1z_tmp * 2.0 * p->harmonic_normb * scale;
                            break;

                        case STIFF:
#ifndef _OPENACC
                            fprintf(stderr, "ERROR: %s:%d stiff bond not yet implemented.\n", __FILE__, __LINE__);
#endif                          //OPENACC
                            break;

                        default:
#ifndef _OPENACC
                            fprintf(stderr, "ERROR: %s:%d unknow bond type appeared %d\n",
                                    __FILE__, __LINE__, bond_type);
#endif                          //OPENACC
                            break;
                        }
                }
        }
    *fx += v1x;
    *fy += v1y;
    *fz += v1z;
}

inline int possible_move_area51(const Phase * p, const soma_scalar_t oldx, const soma_scalar_t oldy,
                                const soma_scalar_t oldz, soma_scalar_t dx, soma_scalar_t dy, soma_scalar_t dz,
                                const int nonexact)
{
    if (p->area51 == NULL)
        return 1;

    const uint64_t idx = coord_to_index(p, oldx + dx, oldy + dy, oldz + dz);
    if (idx >= p->n_cells_local)
        return 0;               //Domain error. Silent, but delta_E should give the appropriate signaling.

    if (p->area51[idx] != 0)
        return 0;

    if (!nonexact)
        {
            const soma_scalar_t r = sqrt(dx * dx + dy * dy + dz * dz);
            const int num_samples = r / p->max_safe_jump + 1;
            if (num_samples > 0)
                {
                    dx /= num_samples;
                    dy /= num_samples;
                    dz /= num_samples;
                    for (int i = 1; i < num_samples + 1; i++)
                        {
                            const soma_scalar_t jx = oldx + i * dx;
                            const soma_scalar_t jy = oldy + i * dy;
                            const soma_scalar_t jz = oldz + i * dz;
                            const uint64_t idx = coord_to_index(p, jx, jy, jz);
                            if (idx >= p->n_cells_local)
                                return 0;       //Domain error.
                            const unsigned int area = p->area51[idx];
                            if (area != 0)
                                return 0;
                        }
                }
        }

    return 2;
}

int set_iteration_possible_move(const Phase * p, RNG_STATE * const set_states, Monomer * const beads,
                                uint64_t chain_index, unsigned int iP, const int nonexact_area51,
                                const unsigned int ibead, const unsigned int iwtype,
                                unsigned int *accepted_moves_set_ptr)
{

    unsigned int accepted_moves_set = *accepted_moves_set_ptr;
    //local copy of rngstate. For fast updates of state in register.
    RNG_STATE my_state = set_states[iP];
    Monomer mybead = beads[ibead];
    Monomer dx;
    dx.x = dx.y = dx.z = 0;
    soma_scalar_t delta_E_bond = 0.0;
    switch (p->args.move_type_arg)
        {
        case move_type_arg_TRIAL:
            trial_move(p, chain_index, ibead, &dx.x, &dx.y, &dx.z, iwtype, &my_state);  // normal MC move
            delta_E_bond = calc_delta_bonded_energy(p, &mybead, chain_index, ibead, dx.x, dx.y, dx.z);
            break;
        case move_type_arg_SMART:
            trial_move_smc(p, chain_index, ibead, &dx.x, &dx.y, &dx.z, &delta_E_bond, &mybead, &my_state, iwtype);      // force biased move
            break;
        case move_type__NULL:
        default:
            delta_E_bond = 0.0;
            break;
        }
    int error = 0;
    const int move_allowed = possible_move_area51(p, mybead.x, mybead.y, mybead.z, dx.x, dx.y, dx.z, nonexact_area51);
    if (move_allowed)
        {
            soma_scalar_t delta_energy = delta_E_bond;
            delta_energy += calc_delta_nonbonded_energy(p, &mybead, dx.x, dx.y, dx.z, iwtype);
            if (delta_energy != delta_energy)   //isnan(delta_energy) ) not working with PGI OpenACC
                error = chain_index + 1;

            const soma_scalar_t old_mobility_modifier = get_mobility_modifier(p, iwtype, mybead.x, mybead.y, mybead.z);
            const soma_scalar_t new_mobility_modifier =
                get_mobility_modifier(p, iwtype, mybead.x + dx.x, mybead.y + dx.y, mybead.z + dx.z);
            const soma_scalar_t mobility_modifier = sqrt(old_mobility_modifier * new_mobility_modifier);

            // MC roll to accept / reject
            if (som_accept(&my_state, p, delta_energy, mobility_modifier) == 1)
                {
                    Monomer newx;
                    newx.x = mybead.x + dx.x;
                    newx.y = mybead.y + dx.y;
                    newx.z = mybead.z + dx.z;
                    beads[ibead] = newx;
#ifndef _OPENACC
                    accepted_moves_set += 1;
#endif                          //_OPENACC
                }
        }
    //Copy the RNGstate back to global memory
    set_states[iP] = my_state;
    *accepted_moves_set_ptr = accepted_moves_set;
    return error;
}
