/* Copyright (C) 2016-2019 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016-2017 Marcel Langenberg
   Copyright (C) 2016 Fabien Leonforte
   Copyright (C) 2016 Juan Orozco
   Copyright (C) 2016 Yongzhi Ren
   Copyright (C) 2016 N. Harshavardhan Reddy

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

//! \file ana.c
//! \brief Implementation of ana.h


#include "ana.h"
#include <stdbool.h>
#include "phase.h"
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "ana_info.h"
#include "mesh.h"
#include "ana_malloc.h"
#include "io.h"

#define IS_TIME_FOR_OBS(delta) ((delta) != 0 && time%(delta)==0)


// given a timestep and ana_info, it creates a list of tasks to be done in this timestep
void find_tasks(const Ana_Info * const ana_info, unsigned int time, /*out*/ bool * needToDo){

    needToDo[Re] = IS_TIME_FOR_OBS(ana_info->delta_mc_Re);
    needToDo[dvar] = IS_TIME_FOR_OBS(ana_info->delta_mc_density_var);
    needToDo[Rg] = IS_TIME_FOR_OBS(ana_info->delta_mc_Rg);
    needToDo[b_anisotropy] = IS_TIME_FOR_OBS(ana_info->delta_mc_b_anisotropy);
    needToDo[acc_ratio] = IS_TIME_FOR_OBS(ana_info->delta_mc_acc_ratio);
    needToDo[MSD] = IS_TIME_FOR_OBS(ana_info->delta_mc_MSD);
    needToDo[non_bonded_energy] = IS_TIME_FOR_OBS(ana_info->delta_mc_non_bonded_energy);
    needToDo[bonded_energy] = IS_TIME_FOR_OBS(ana_info->delta_mc_bonded_energy);
    needToDo[field] = IS_TIME_FOR_OBS(ana_info->delta_mc_density_field);
    needToDo[umbrella_field] = IS_TIME_FOR_OBS(ana_info->delta_mc_umbrella_field);
    needToDo[dynamical_structure] = IS_TIME_FOR_OBS(ana_info->delta_mc_dynamical_structure);
    needToDo[static_structure] = IS_TIME_FOR_OBS(ana_info->delta_mc_static_structure);
    needToDo[dump] = IS_TIME_FOR_OBS(ana_info->delta_mc_dump);

}

// given a list of tasks ("needToDo" from find_tasks) returns true if it will be necessary to call
// update_self_phase
// This List was created by looking at what calculations in analytics call "update self phase"
bool deviceToHostCopyIsNecessary(const bool *needToDo){
    return needToDo[Re] || needToDo[dvar] || needToDo[Rg] || needToDo[b_anisotropy] ||
    needToDo[MSD] || needToDo[bonded_energy] || needToDo[dump];
}

// decides if *any* task has to be done
bool any(const bool *needToDo){
    for (int i=0; i < obs_it_end; i++){
        if (needToDo[i]){
            return true;
        }
    }
    return false;
}



int analytics(struct Phase *const p)
{
    static unsigned int last_time_call = 0;
    if (last_time_call == 0 || p->time > last_time_call)
        last_time_call = p->time;
    else                        //Quick exit, because the property has already been calculated for the time step.
        return 0;

    //determine what needs to be calculated
    bool needToDo[13]; //there are 13 == SOMA_NUM_OBS observables. todo: this but nicer
    find_tasks(& p->ana_info, p->time, needToDo);

    if (!any(needToDo)){
        // nothing to do on this timestep
        return 0;
    }


    //Copy data to from device to host if necessary
    if (deviceToHostCopyIsNecessary(needToDo)){
        update_self_phase(p, 0);
    }

    //create necessary communication data structures
    ana_malloc(p->n_poly_type,p->n_types,p->ana_info.q_size_dynamical,p->ana_info.q_size_static, needToDo);


    //precomm calc

    if (needToDo[Re]){
        uint64_t *const re_counter = get_re_counter();
        soma_scalar_t *const re_result = get_re_vals();
        if (re_counter == NULL || re_result == NULL)
        {
            fprintf(stderr, "MALLOC ERROR: %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
        calc_Re_local(p, re_result, re_counter);
    }

    if (needToDo[dvar]){
        soma_scalar_t *dvar = get_dvar_val();
        calc_dvar_local(p, dvar);

    }

    if (needToDo[Rg]){
        soma_scalar_t *const Rg = get_rg_vals();
        uint64_t *const counter = get_rg_counter();
        if (Rg == NULL || counter == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d \n", __FILE__, __LINE__);
            return -2;
        }
        calc_Rg_local(p,Rg,counter);
    }

    if (needToDo[b_anisotropy]){
        uint64_t *const counter = get_anisotropy_counter();
        soma_scalar_t *const a = get_anisotropy_vals();
        if (a==NULL || counter==NULL){
            fprintf(stderr, "ERROR: Malloc %s:%d \n", __FILE__, __LINE__);
            return -2;
        }
        calc_anisotropy_local(p, a, counter);
    }

    if (needToDo[acc_ratio]){
        //todo: Openacc-Prüfung einbauen.
        uint64_t * moves = get_moves();
        uint64_t * accepts = get_accepts();
        calc_acc_ratio_local(p, accepts, moves);
    }

    if (needToDo[MSD]) {
        soma_scalar_t *const msd = get_msd_vals();
        uint64_t *const counter =  get_msd_counter();
        if (msd == NULL || counter == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d \n", __FILE__, __LINE__);
            return -2;
        }
        calc_msd_local(p, msd, counter);
    }

    if (needToDo[non_bonded_energy]){
        soma_scalar_t *const nb_energy = get_non_bonded_energy_val();
        if (nb_energy == NULL)
            {
                fprintf(stderr, "ERROR: Malloc %s:%d \n", __FILE__, __LINE__);
                return -2;
            }
#pragma acc update self(p->omega_field_unified[0:p->n_cells_local*p->n_types])
#pragma acc update self(p->fields_unified[0:p->n_cells_local*p->n_types])
            calc_non_bonded_energy_local(p, nb_energy);
    }

    if (needToDo[bonded_energy]) {
        soma_scalar_t *const b_energy =
                get_bonded_energy_val();
        if (b_energy == NULL) {
            fprintf(stderr, "ERROR: Malloc %s:%d \n", __FILE__, __LINE__);
            return -2;

        }
        calc_bonded_energy_local(p, b_energy);
    }

    if (needToDo[dynamical_structure]){
        //update_self_phase(p,0); //update not needed, because the calculation is on the device
        soma_scalar_t *const dynamical_structure_factor = get_dynamical_sf_val();
        if (dynamical_structure_factor == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d \n", __FILE__, __LINE__);
            return -2;
        }
        enum structure_factor_type sf_type = DYNAMICAL_STRUCTURE_FACTOR;
        calc_structure_local(p, dynamical_structure_factor, sf_type);
    }

    if (needToDo[static_structure]){
        //update_self_phase(p,0); //update not needed, because the calculation is on the device
        soma_scalar_t *const static_structure_factor = get_static_sf_val();
        if (static_structure_factor == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d \n", __FILE__, __LINE__);
            return -2;
        }
        enum structure_factor_type sf_type = STATIC_STRUCTURE_FACTOR;
        calc_structure_local(p, static_structure_factor, sf_type);
    }

    //communication
    ana_reduce(0,p->info_MPI.SOMA_comm_sim);

    //postcomm calc & file-IO
    if (p->info_MPI.sim_rank == 0){

        if (needToDo[Re]){
            uint64_t *counter = get_re_counter();
            double *result = get_re_vals();
            for (unsigned int type = 0; type < p->n_poly_type; type++){
                for (unsigned int i = 0; i < 4; i++) {
                    if (counter[type] > 0) {
                        result[type * 4 + i] /= (soma_scalar_t) counter[type];
                    }
                }
            }
            extent_ana_by_field(result, 4*p->n_poly_type, "/Re", p->ana_info.file_id);
        }

        if (needToDo[dvar]){
            double *dvar = get_dvar_val();
            *dvar /= p->args.N_domains_arg;
            extent_ana_by_field(dvar, 1, "/density_var", p->ana_info.file_id);
        }

        if (needToDo[Rg]){
            uint64_t *counter = get_rg_counter();
            double *rg = get_rg_vals();
            for (unsigned int type = 0; type < p->n_poly_type; type++)
                for (unsigned int i = 0; i < 4; i++)
                    if (counter[type] > 0)
                        rg[type * 4 + i] /= (soma_scalar_t) counter[type];

            extent_ana_by_field(rg, 4*p->n_poly_type, "/Rg", p->ana_info.file_id);
        }

        if (needToDo[b_anisotropy]){
            uint64_t *counter = get_anisotropy_counter();
            double *ba = get_anisotropy_vals();

            for (unsigned int type = 0; type < p->n_poly_type; type++) {
                for (unsigned int i = 0; i < 6; i++) {
                    if (counter[type] > 0) {
                        ba[type * 6 + i] /= (soma_scalar_t) counter[type];
                    }
                }
            }
            extent_ana_by_field(ba, 6*p->n_poly_type, "/bond_anisotropy", p->ana_info.file_id);
        }

        if (needToDo[acc_ratio]){
            // ist es Absicht dass die acc_ratio immer zu 0 berechnet wird, wenn MPI disabled ist?
            // im originalen ana.c ist das so.
            uint64_t accepts = *get_accepts();
            uint64_t moves = *get_moves();
            double acc_ratio = (double) accepts / moves;
            extent_ana_by_field(&acc_ratio, 1, "/acc_ratio", p->ana_info.file_id);
        }

        if (needToDo[MSD]){
            uint64_t *counter = get_msd_counter();
            double *result = get_msd_vals();
            //Looping over twice the number of poly types. But loop over half the elements
            // 8/2 = 4, because the norm for first and second half differ.
            for (unsigned int type = 0; type < 2 * p->n_poly_type; type++)
                for (unsigned int i = 0; i < 4; i++)
                    if (counter[type] > 0)
                        result[type * 4 + i] /= (soma_scalar_t) counter[type];
            extent_ana_by_field(result, 8*p->n_poly_type, "/MSD", p->ana_info.file_id);
        }

        if (needToDo[non_bonded_energy]){
            double *nb_energy = get_non_bonded_energy_val();
            extent_ana_by_field(nb_energy, p->n_types, "/non_bonded_energy", p->ana_info.file_id);
        }

        if (needToDo[bonded_energy]){
            double *b_energy = get_bonded_energy_val();
            extent_ana_by_field(b_energy, NUMBER_SOMA_BOND_TYPES, "/bonded_energy", p->ana_info.file_id);
        }

        if (needToDo[dynamical_structure]){
            double * dynamical_sf = get_dynamical_sf_val();
            enum structure_factor_type sf_type = DYNAMICAL_STRUCTURE_FACTOR;
            extent_structure(p, dynamical_sf, "/dynamical_structure_factor", p->ana_info.file_id, sf_type);
        }

        if (needToDo[static_structure]){
            double *static_sf = get_static_sf_val();
            enum structure_factor_type sf_type = STATIC_STRUCTURE_FACTOR;
            extent_structure(p, static_sf, "/static_structure_factor", p->ana_info.file_id, sf_type);
        }
    }


    // other observables that don't fit the pattern

    if (needToDo[dump]){
        const unsigned int len = 1024;
        char filename[1024];        //len
        memset(filename, '\0', len * sizeof(char));
        const unsigned int check = sprintf(filename, "%s_t%u.h5", p->ana_info.coord_filename, p->time);
        if (check >= len)
            fprintf(stderr, "ERROR: %s %d cannot write filename.\n", __FILE__, __LINE__);
        else
            write_config_hdf5(p, filename);
    }

    if (needToDo[umbrella_field]){
        if (p->info_MPI.sim_size == 1)
        {
#pragma acc update self(p->fields_unified[0:p->n_cells*p->n_types])
        }
        extent_density_field(p, p->umbrella_field, "/umbrella_field", H5T_SOMA_NATIVE_SCALAR, MPI_SOMA_SCALAR,
                             sizeof(soma_scalar_t));
    }

    if (needToDo[field]){
        // if we run on one core+gpu only we avoid the stepwise field transfer, so for writing fields
        // we need to update the local field
        if (p->info_MPI.sim_size == 1)
        {
#pragma acc update self(p->fields_unified[0:p->n_cells*p->n_types])
        }

        //Collective IO, not yet.
        extent_density_field(p, p->fields_unified, "/density_field", H5T_NATIVE_UINT16, MPI_UINT16_T,
                             sizeof(uint16_t));
    }

    //flushing (always necessary if we get to this point, as we do a quick exit for no obs being calculated)

    herr_t status = H5Fflush(p->ana_info.file_id, H5F_SCOPE_LOCAL);
    if (status < 0)
    {
        fprintf(stderr, "ERROR: HDF5 %s:%s:%d\n", __func__, __FILE__, __LINE__);
        return -1;
    }


    //freeing resources
    ana_free();

    return 0;

}

void calc_acc_ratio_local(struct Phase *const p, uint64_t *accepts, uint64_t *moves) {

    // copy local results from phase
    *accepts = p->n_accepts;
    *moves = p->n_moves;

    // zero out for next calculation
    p->n_accepts = 0;
    p->n_moves = 0;

}

void calc_msd_local(struct Phase *const p, soma_scalar_t *const msd, uint64_t *const counter) {
    memset(counter, 0, p->n_poly_type * sizeof(uint64_t));
    memset(msd, 0, 6 * p->n_poly_type * sizeof(soma_scalar_t));

    // loop over local chains
    for (uint64_t ipoly = 0; ipoly < p->n_polymers; ipoly++)
    {
        const Polymer *const poly = &(p->polymers[ipoly]);
        const unsigned int type = poly->type;
        const unsigned int N = p->poly_arch[p->poly_type_offset[type]];

        // loop over beads in this chain
        for (unsigned int ibead = 0; ibead < N; ibead++)
        {

            const soma_scalar_t x1 = poly->beads[ibead].x;
            const soma_scalar_t y1 = poly->beads[ibead].y;
            const soma_scalar_t z1 = poly->beads[ibead].z;

            // loop over bonds of this bead
            const int start = get_bondlist_offset(p->poly_arch[p->poly_type_offset[type] + ibead + 1]);
            if (start > 0)
            {
                int i = start;
                int end;
                do
                {
                    const uint32_t bn = p->poly_arch[i++];
                    const int info = bn;
                    end = get_end(info);
                    const int offset = get_offset(info);
                    const int neighbour_id = ibead + offset;
                    const unsigned int jbead = neighbour_id;

                    const soma_scalar_t x2 = p->polymers[ipoly].beads[jbead].x;
                    const soma_scalar_t y2 = p->polymers[ipoly].beads[jbead].y;
                    const soma_scalar_t z2 = p->polymers[ipoly].beads[jbead].z;


                    const int mic_flag = p->args.bond_minimum_image_convention_flag;

                    const soma_scalar_t bx = calc_bond_length(x2,x1,p->Lx,mic_flag);
                    const soma_scalar_t by = calc_bond_length(y2,y1,p->Ly,mic_flag);
                    const soma_scalar_t bz = calc_bond_length(z2,z1,p->Lz,mic_flag);

                    msd[type * 6 + 0] += bx * bx;
                    msd[type * 6 + 1] += by * by;
                    msd[type * 6 + 2] += bz * bz;
                    msd[type * 6 + 3] += bx * by;
                    msd[type * 6 + 4] += by * bz;
                    msd[type * 6 + 5] += bz * bx;
                    counter[type] += 1;
                } while (end == 0);
            }
        }
    }
}

int calc_structure_local(const struct Phase *p, soma_scalar_t *const result, const enum structure_factor_type sf_type) {
    unsigned int error = 0;
    unsigned int q_size, result_tmp_size;
    soma_scalar_t *q_array;

    switch (sf_type)
    {
        case DYNAMICAL_STRUCTURE_FACTOR:
            q_size = p->ana_info.q_size_dynamical;
            q_array = p->ana_info.q_dynamical;
            result_tmp_size = q_size * p->n_types * 4;
            break;
        case STATIC_STRUCTURE_FACTOR:
            q_size = p->ana_info.q_size_static;
            q_array = p->ana_info.q_static;
            result_tmp_size = q_size * p->n_types * 2;
            break;
        default:
            fprintf(stderr, "ERROR: Not correct structure_factor_type. %s:%s:%d\n", __func__, __FILE__, __LINE__);
            return -1;
    }
    uint64_t *const counter = (uint64_t * const)calloc(p->n_poly_type, sizeof(uint64_t));
    if (counter == NULL)
    {
        fprintf(stderr, "MALLOC ERROR: %s:%d\n", __FILE__, __LINE__);
        return -1;
    }
    memset(counter, 0, p->n_poly_type * sizeof(uint64_t));
    memset(result, 0, q_size * p->n_poly_type * p->n_types * p->n_types * sizeof(soma_scalar_t));
    unsigned int n_random_q = p->args.n_random_q_arg;
    soma_scalar_t *const result_tmp =
            (soma_scalar_t * const)malloc(n_random_q * p->n_polymers * result_tmp_size * sizeof(soma_scalar_t));
    if (result_tmp == NULL)
    {
        fprintf(stderr, "MALLOC ERROR: %s:%d\n", __FILE__, __LINE__);
        return -1;
    }
    for (unsigned int index = 0; index < n_random_q * p->n_polymers * result_tmp_size; index++)
    {
        result_tmp[index] = 0;
    }
    soma_scalar_t *const tmp =
            (soma_scalar_t * const)malloc(n_random_q * p->n_polymers * q_size * p->n_types * p->n_types *
                                          sizeof(soma_scalar_t));
    if (tmp == NULL)
    {
        fprintf(stderr, "MALLOC ERROR: %s:%d\n", __FILE__, __LINE__);
        return -1;
    }
    memset(tmp, 0, n_random_q * p->n_polymers * q_size * p->n_types * p->n_types * sizeof(soma_scalar_t));
    enum enum_pseudo_random_number_generator arg_rng_type = p->args.pseudo_random_number_generator_arg;

#pragma acc enter data copyin(result_tmp[0:n_random_q*p->n_polymers*result_tmp_size],q_array[0:q_size])
#pragma acc enter data copyin(tmp[0:n_random_q*p->n_polymers*q_size*p->n_types*p->n_types]) async
#pragma acc parallel loop vector_length(32) present(p[0:1]) async
#pragma omp parallel for
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
    {
        const unsigned int poly_type = p->polymers[poly].type;
        unsigned int poly_length = p->poly_arch[p->poly_type_offset[poly_type]];
        RNG_STATE *const s = &(p->polymers[poly].poly_state);
        //random q generation

#pragma acc loop                //be careful, seq?
        for (unsigned int index_random_q = 0; index_random_q < n_random_q; index_random_q++)
        {
            soma_scalar_t rng1, rng2;
#pragma acc loop seq
            for (unsigned int random_i = 0; random_i < index_random_q; random_i++)
            {
                rng1 = soma_rng_uint(s, arg_rng_type);
                rng2 = soma_rng_uint(s, arg_rng_type);
            }
            rng1 = (uint32_t) soma_rng_uint(s, arg_rng_type) / (soma_scalar_t) soma_rng_uint_max();
            rng2 = (uint32_t) soma_rng_uint(s, arg_rng_type) / (soma_scalar_t) soma_rng_uint_max();
            soma_scalar_t theta = 2 * M_PI * rng1;
            soma_scalar_t phi = acos(1 - 2 * rng2);
            soma_scalar_t unit_q_x = sin(phi) * cos(theta);
            soma_scalar_t unit_q_y = sin(phi) * sin(theta);
            soma_scalar_t unit_q_z = cos(phi);
#pragma acc loop seq
            for (unsigned int mono = 0; mono < poly_length; mono++)
            {
                const unsigned int particle_type =
                        get_particle_type(p->poly_arch[p->poly_type_offset[poly_type] + 1 + mono]);
                // Monomer position at t.
                soma_scalar_t x = p->polymers[poly].beads[mono].x;
                soma_scalar_t y = p->polymers[poly].beads[mono].y;
                soma_scalar_t z = p->polymers[poly].beads[mono].z;
#pragma acc loop seq
                for (unsigned int index_q = 0; index_q < q_size; index_q++)
                {
                    soma_scalar_t qr = q_array[index_q] * (unit_q_x * x + unit_q_y * y + unit_q_z * z);

                    switch (sf_type)
                    {
                        case STATIC_STRUCTURE_FACTOR:
                            result_tmp[index_random_q * p->n_polymers * q_size * p->n_types * 2 +
                                       poly * q_size * p->n_types * 2 + index_q * p->n_types * 2 +
                                       particle_type * 2 + 0] += cos(qr);
                            result_tmp[index_random_q * p->n_polymers * q_size * p->n_types * 2 +
                                       poly * q_size * p->n_types * 2 + index_q * p->n_types * 2 +
                                       particle_type * 2 + 1] += sin(qr);
                            break;
                        case DYNAMICAL_STRUCTURE_FACTOR:
                            ;
                            soma_scalar_t x_0 = p->polymers[poly].msd_beads[mono].x;
                            soma_scalar_t y_0 = p->polymers[poly].msd_beads[mono].y;
                            soma_scalar_t z_0 = p->polymers[poly].msd_beads[mono].z;
                            soma_scalar_t qr_msd =
                                    q_array[index_q] * (unit_q_x * x_0 + unit_q_y * y_0 + unit_q_z * z_0);
                            result_tmp[index_random_q * p->n_polymers * q_size * p->n_types * 4 +
                                       poly * q_size * p->n_types * 4 + index_q * p->n_types * 4 +
                                       particle_type * 4 + 0] += cos(qr);
                            result_tmp[index_random_q * p->n_polymers * q_size * p->n_types * 4 +
                                       poly * q_size * p->n_types * 4 + index_q * p->n_types * 4 +
                                       particle_type * 4 + 1] += sin(qr);
                            result_tmp[index_random_q * p->n_polymers * q_size * p->n_types * 4 +
                                       poly * q_size * p->n_types * 4 + index_q * p->n_types * 4 +
                                       particle_type * 4 + 2] += cos(qr_msd);
                            result_tmp[index_random_q * p->n_polymers * q_size * p->n_types * 4 +
                                       poly * q_size * p->n_types * 4 + index_q * p->n_types * 4 +
                                       particle_type * 4 + 3] += sin(qr_msd);

                            break;
                        default:
                            error = -1;
                    }
                }
            }
        }
    }
    if (error != 0)
    {
        fprintf(stderr, "ERROR: %d unknown structure factor type %s:%d\n", error, __FILE__, __LINE__);
        return error;
    }
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
    {
        const unsigned int poly_type = p->polymers[poly].type;
        counter[poly_type]++;
    }
#pragma acc wait
#pragma acc parallel loop gang vector present(p[0:1])
#pragma omp parallel for
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
    {
        const unsigned int poly_type = p->polymers[poly].type;
        unsigned int poly_length = p->poly_arch[p->poly_type_offset[poly_type]];
#pragma acc loop seq
        for (unsigned int index_random_q = 0; index_random_q < n_random_q; index_random_q++)
        {
#pragma acc loop seq
            for (unsigned int particle_type_i = 0; particle_type_i < p->n_types; particle_type_i++)
            {
#pragma acc loop seq
                for (unsigned int particle_type_j = 0; particle_type_j < p->n_types; particle_type_j++)
                {
#pragma acc loop seq
                    for (unsigned int index_q = 0; index_q < q_size; index_q++)
                    {
                        switch (sf_type)
                        {
                            case STATIC_STRUCTURE_FACTOR:
                                tmp[index_random_q * p->n_polymers * q_size * p->n_types *
                                    p->n_types + poly * q_size * p->n_types * p->n_types +
                                    index_q * p->n_types * p->n_types +
                                    particle_type_i * p->n_types + particle_type_j] +=
                                        (result_tmp
                                         [index_random_q * p->n_polymers * q_size * p->n_types * 2 +
                                          poly * q_size * p->n_types * 2 + index_q * p->n_types * 2 +
                                          particle_type_i * 2 +
                                          0] * result_tmp[index_random_q * p->n_polymers * q_size *
                                                          p->n_types * 2 +
                                                          poly * q_size * p->n_types * 2 +
                                                          index_q * p->n_types * 2 +
                                                          particle_type_j * 2 + 0] +
                                         result_tmp[index_random_q * p->n_polymers * q_size *
                                                    p->n_types * 2 + poly * q_size * p->n_types * 2 +
                                                    index_q * p->n_types * 2 + particle_type_i * 2 +
                                                    1] * result_tmp[index_random_q * p->n_polymers *
                                                                    q_size * p->n_types * 2 +
                                                                    poly * q_size * p->n_types * 2 +
                                                                    index_q * p->n_types * 2 +
                                                                    particle_type_j * 2 +
                                                                    1]) / poly_length;
                                break;
                            case DYNAMICAL_STRUCTURE_FACTOR:
                                tmp[index_random_q * p->n_polymers * q_size * p->n_types *
                                    p->n_types + poly * q_size * p->n_types * p->n_types +
                                    index_q * p->n_types * p->n_types +
                                    particle_type_i * p->n_types + particle_type_j] +=
                                        (result_tmp
                                         [index_random_q * p->n_polymers * q_size * p->n_types * 4 +
                                          poly * q_size * p->n_types * 4 + index_q * p->n_types * 4 +
                                          particle_type_i * 4 +
                                          0] * result_tmp[index_random_q * p->n_polymers * q_size *
                                                          p->n_types * 4 +
                                                          poly * q_size * p->n_types * 4 +
                                                          index_q * p->n_types * 4 +
                                                          particle_type_j * 4 + 2] +
                                         result_tmp[index_random_q * p->n_polymers * q_size *
                                                    p->n_types * 4 + poly * q_size * p->n_types * 4 +
                                                    index_q * p->n_types * 4 + particle_type_i * 4 +
                                                    1] * result_tmp[index_random_q * p->n_polymers *
                                                                    q_size * p->n_types * 4 +
                                                                    poly * q_size * p->n_types * 4 +
                                                                    index_q * p->n_types * 4 +
                                                                    particle_type_j * 4 +
                                                                    3]) / poly_length;
                                break;
                            default:
                                error = -1;
                        }
                    }
                }
            }
        }               //index_random_q
    }                       //poly
    if (error != 0)
    {
        fprintf(stderr, "ERROR: %d unknown structure factor type %s:%d\n", error, __FILE__, __LINE__);
        return error;
    }
#pragma acc exit data copyout(result_tmp[0:n_random_q*p->n_polymers*result_tmp_size],q_array[0:q_size])
#pragma acc exit data copyout(tmp[0:n_random_q*p->n_polymers*q_size*p->n_types*p->n_types])

    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
    {
        unsigned int poly_type = p->polymers[poly].type;
        for (unsigned int index_random_q = 0; index_random_q < n_random_q; index_random_q++)
        {
            for (unsigned int particle_type_i = 0; particle_type_i < p->n_types; particle_type_i++)
            {
                for (unsigned int particle_type_j = 0; particle_type_j < p->n_types; particle_type_j++)
                {
                    for (unsigned int index_q = 0; index_q < q_size; index_q++)
                    {
                        result[index_q * p->n_poly_type * p->n_types * p->n_types +
                               poly_type * p->n_types * p->n_types + particle_type_i * p->n_types +
                               particle_type_j] +=
                                (tmp
                                 [index_random_q * p->n_polymers * q_size * p->n_types * p->n_types +
                                  poly * q_size * p->n_types * p->n_types +
                                  index_q * p->n_types * p->n_types + particle_type_i * p->n_types +
                                  particle_type_j] / (soma_scalar_t) counter[poly_type]) / n_random_q;
                    }
                }
            }
        }
    }
    return 0;
}

void calc_bonded_energy_local(const struct Phase *p, soma_scalar_t *const bonded_energy) {
    memset(bonded_energy, 0, NUMBER_SOMA_BOND_TYPES * sizeof(soma_scalar_t));

    for (unsigned int poly = 0; poly < p->n_polymers; poly++)
    {
        const unsigned int type = p->polymers[poly].type;
        const unsigned int N = p->poly_arch[p->poly_type_offset[type]];

        for (unsigned int mono = 0; mono < N; mono++)
        {

            const int start = get_bondlist_offset(p->poly_arch[p->poly_type_offset[type] + mono + 1]);

            if (start > 0)
            {
                int i = start;
                unsigned int end;
                do
                {
                    const uint32_t info = p->poly_arch[i++];
                    end = get_end(info);
                    const unsigned int bond_type = get_bond_type(info);
                    const int offset = get_offset(info);
                    if (offset > 0)     //Select each bond only once, i<j
                    {
                        const int mono_j = mono + offset;
                        const int mic_flag = p->args.bond_minimum_image_convention_flag;
                        const soma_scalar_t dx = calc_bond_length(p->polymers[poly].beads[mono].x,
                                                                  p->polymers[poly].beads[mono_j].x,
                                                                  p->Lx,
                                                                  mic_flag);
                        const soma_scalar_t dy = calc_bond_length(p->polymers[poly].beads[mono].y,
                                                                  p->polymers[poly].beads[mono_j].y,
                                                                  p->Ly,
                                                                  mic_flag);
                        const soma_scalar_t dz = calc_bond_length(p->polymers[poly].beads[mono].z,
                                                                  p->polymers[poly].beads[mono_j].z,
                                                                  p->Lz,
                                                                  mic_flag);

                        const soma_scalar_t r2 = dx * dx + dy * dy + dz * dz;

                        soma_scalar_t energy = 0;
                        soma_scalar_t scale = 1.;
                        switch (bond_type)
                        {
                            case HARMONICVARIABLESCALE:
                                scale = p->harmonic_normb_variable_scale;
                                /* intentionally falls through */
                            case HARMONIC:
                                energy = p->harmonic_normb * r2 * scale;
                                break;

                            case STIFF:
                                fprintf(stderr,
                                        "ERROR: %s:%d stiff bond not yet implemented.\n",
                                        __FILE__, __LINE__);
                                break;
                            default:
                                fprintf(stderr, "ERROR: %s:%d unknow bond type appeared %d\n",
                                        __FILE__, __LINE__, bond_type);
                                break;
                        }
                        assert(bond_type < NUMBER_SOMA_BOND_TYPES);
                        bonded_energy[bond_type] += energy;
                    }
                } while (end == 0);
            }
        }
    }
}

void calc_non_bonded_energy_local(const struct Phase *const p, soma_scalar_t *const non_bonded_energy) {
    memset(non_bonded_energy, 0, p->n_types * sizeof(soma_scalar_t));

    if (p->info_MPI.domain_rank == 0)   //Only the root domain rank needs to calculate the values.
        {
            const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
            for (unsigned int type = 0; type < p->n_types; type++)
                {
                    for (unsigned int x = my_domain * (p->nx / p->args.N_domains_arg);
                         x < (my_domain + 1) * (p->nx / p->args.N_domains_arg); x++)
                        for (unsigned int y = 0; y < p->ny; y++)
                            for (unsigned int z = 0; z < p->nz; z++)
                                {
                                    const uint64_t cell = cell_coordinate_to_index(p, x, y, z);
                                    non_bonded_energy[type] +=
                                        p->omega_field_unified[cell + type * p->n_cells_local]
                                        * p->fields_unified[cell + type * p->n_cells_local];
                                }
                }
        }
}

void calc_dvar_local(const struct Phase *p, soma_scalar_t *const local_result) {
    uint64_t var = 0.0;
    if (p->info_MPI.domain_rank == 0)   // Only the root domain rank needs to calculate the value
    {
        const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
        for (unsigned int type = 0; type < p->n_types; type++)
            for (unsigned int x = my_domain * (p->nx / p->args.N_domains_arg);
                 x < (my_domain + 1) * (p->nx / p->args.N_domains_arg); x++)
                for (unsigned int y = 0; y < p->ny; y++)
                    for (unsigned int z = 0; z < p->nz; z++)
                    {
                        const uint64_t index = cell_coordinate_to_index(p, x, y, z) + type * p->n_cells_local;
                        const uint64_t value = p->fields_unified[index] - p->old_fields_unified[index];
                        var += value * value;
                    }
    }
    memcpy(p->old_fields_unified, p->fields_unified, p->n_cells_local * p->n_types * sizeof(uint16_t));
    *local_result = var;                //Cast to float
    *local_result /= p->n_cells_local * p->n_types;
}

void calc_anisotropy_local(struct Phase *const p, soma_scalar_t *const local_result, uint64_t *const local_counter) {
    memset( local_counter, 0, p->n_poly_type * sizeof(uint64_t));
    memset(local_result, 0, 6 * p->n_poly_type * sizeof(soma_scalar_t));

    // loop over local chains
    for (uint64_t ipoly = 0; ipoly < p->n_polymers; ipoly++)
    {
        const Polymer *const poly = &(p->polymers[ipoly]);
        const unsigned int type = poly->type;
        const unsigned int N = p->poly_arch[p->poly_type_offset[type]];

        // loop over beads in this chain
        for (unsigned int ibead = 0; ibead < N; ibead++)
        {

            const soma_scalar_t x1 = poly->beads[ibead].x;
            const soma_scalar_t y1 = poly->beads[ibead].y;
            const soma_scalar_t z1 = poly->beads[ibead].z;

            // loop over bonds of this bead
            const int start = get_bondlist_offset(p->poly_arch[p->poly_type_offset[type] + ibead + 1]);
            if (start > 0)
            {
                int i = start;
                int end;
                do
                {
                    const uint32_t bn = p->poly_arch[i++];
                    const int info = bn;
                    end = get_end(info);
                    const int offset = get_offset(info);
                    const int neighbour_id = ibead + offset;
                    const unsigned int jbead = neighbour_id;

                    const soma_scalar_t x2 = p->polymers[ipoly].beads[jbead].x;
                    const soma_scalar_t y2 = p->polymers[ipoly].beads[jbead].y;
                    const soma_scalar_t z2 = p->polymers[ipoly].beads[jbead].z;


                    const int mic_flag = p->args.bond_minimum_image_convention_flag;

                    const soma_scalar_t bx = calc_bond_length(x2,x1,p->Lx,mic_flag);
                    const soma_scalar_t by = calc_bond_length(y2,y1,p->Ly,mic_flag);
                    const soma_scalar_t bz = calc_bond_length(z2,z1,p->Lz,mic_flag);

                    local_result[type * 6 + 0] += bx * bx;
                    local_result[type * 6 + 1] += by * by;
                    local_result[type * 6 + 2] += bz * bz;
                    local_result[type * 6 + 3] += bx * by;
                    local_result[type * 6 + 4] += by * bz;
                    local_result[type * 6 + 5] += bz * bx;
                    local_counter[type] += 1;
                } while (end == 0);
            }
        }
    }

}

void calc_Re_local(const struct Phase *p, soma_scalar_t * local_result, uint64_t * local_counter) {

    memset(local_counter , 0, p->n_poly_type * sizeof(uint64_t));
    memset(local_result, 0, 4 * p->n_poly_type * sizeof(soma_scalar_t));

    for (uint64_t npoly = 0; npoly < p->n_polymers; npoly++)
    {
        const unsigned int type = p->polymers[npoly].type;

        const unsigned int start = p->end_mono[type * 2 + 0];
        const unsigned int end = p->end_mono[type * 2 + 1];
        const soma_scalar_t dx = p->polymers[npoly].beads[start].x - p->polymers[npoly].beads[end].x;
        const soma_scalar_t dy = p->polymers[npoly].beads[start].y - p->polymers[npoly].beads[end].y;
        const soma_scalar_t dz = p->polymers[npoly].beads[start].z - p->polymers[npoly].beads[end].z;

        local_result[type * 4 + 0] += dx * dx + dy * dy + dz * dz;
        local_result[type * 4 + 1] += dx * dx;
        local_result[type * 4 + 2] += dy * dy;
        local_result[type * 4 + 3] += dz * dz;
        local_counter [type] += 1;
    }
}



void calc_Rg_local(const struct Phase *p, soma_scalar_t *const local_result, uint64_t *const local_counter) {
    memset(local_counter, 0, p->n_poly_type * sizeof(uint64_t));
    memset(local_result, 0, 4 * p->n_poly_type * sizeof(soma_scalar_t));

    for (uint64_t ipoly = 0; ipoly < p->n_polymers; ipoly++)
    {
        const unsigned int type = p->polymers[ipoly].type;
        const unsigned int N = p->poly_arch[p->poly_type_offset[type]];
        soma_scalar_t xcm = 0.;
        soma_scalar_t ycm = 0.;
        soma_scalar_t zcm = 0.;
        soma_scalar_t x2 = 0.;
        soma_scalar_t y2 = 0.;
        soma_scalar_t z2 = 0.;
        for (unsigned int ibead = 0; ibead < N; ibead++)
        {
            const soma_scalar_t x1 = p->polymers[ipoly].beads[ibead].x;
            const soma_scalar_t y1 = p->polymers[ipoly].beads[ibead].y;
            const soma_scalar_t z1 = p->polymers[ipoly].beads[ibead].z;
            xcm += x1;
            ycm += y1;
            zcm += z1;
            x2 += x1 * x1;
            y2 += y1 * y1;
            z2 += z1 * z1;
        }
        xcm /= (soma_scalar_t) (N);
        ycm /= (soma_scalar_t) (N);
        zcm /= (soma_scalar_t) (N);

        local_result[type * 4 + 0] += (x2 / (soma_scalar_t) (N) - xcm * xcm) +
                                      (y2 / (soma_scalar_t) (N) - ycm * ycm) + (z2 / (soma_scalar_t) (N) - zcm * zcm);
        local_result[type * 4 + 1] += (x2 / (soma_scalar_t) (N) - xcm * xcm);
        local_result[type * 4 + 2] += (y2 / (soma_scalar_t) (N) - ycm * ycm);
        local_result[type * 4 + 3] += (z2 / (soma_scalar_t) (N) - zcm * zcm);
        local_counter[type] += 1;
    }
}


int extent_density_field(const struct Phase *const p, const void *const field_pointer, const char *const field_name,
                         hid_t hdf5_type, const MPI_Datatype mpi_type, const size_t data_size)
{
    const char *const name = field_name;
    update_density_fields(p);

    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    const unsigned int buffer_size = (p->nx / p->args.N_domains_arg) * p->ny * p->nz;
    const unsigned int ghost_buffer_size = p->args.domain_buffer_arg * p->ny * p->nz;

    if (p->info_MPI.sim_rank == 0)
    {
        herr_t status;

        hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
        HDF5_ERROR_CHECK(plist_id);

        //Open the dataset and space
        hid_t dset = H5Dopen(p->ana_info.file_id, name, H5P_DEFAULT);
        HDF5_ERROR_CHECK(dset);
        hid_t d_space = H5Dget_space(dset);
        HDF5_ERROR_CHECK(d_space);

        // Extent the dimesion of the density field.
        const unsigned int ndims = H5Sget_simple_extent_ndims(d_space);
        if (ndims != 5)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent density field %s.\n",
                    __FILE__, __LINE__, name);
            return -1;
        }
        hsize_t dims[5];    //ndims
        status = H5Sget_simple_extent_dims(d_space, dims, NULL);
        HDF5_ERROR_CHECK(status);

        hsize_t dims_new[5];        //ndims
        dims_new[0] = dims[0] + 1;
        dims_new[1] = p->n_types;
        assert(dims[1] == p->n_types);
        dims_new[2] = p->nx;
        assert(dims[2] == p->nx);
        dims_new[3] = p->ny;
        assert(dims[3] == p->ny);
        dims_new[4] = dims[4];
        assert(dims[4] == p->nz);

        status = H5Dset_extent(dset, dims_new);
        HDF5_ERROR_CHECK(status);

        status = H5Sclose(d_space);
        HDF5_ERROR_CHECK(status);
        status = H5Dclose(dset);
        HDF5_ERROR_CHECK(status);

        //Open the new filespace
        dset = H5Dopen(p->ana_info.file_id, name, H5P_DEFAULT);
        HDF5_ERROR_CHECK(dset);
        hid_t filespace = H5Dget_space(dset);
        HDF5_ERROR_CHECK(filespace);

        hsize_t dims_memspace[5];   //ndims
        dims_memspace[0] = dims_new[0] - dims[0];
        dims_memspace[1] = 1;
        dims_memspace[2] = p->nx / p->args.N_domains_arg;
        dims_memspace[3] = p->ny;
        dims_memspace[4] = p->nz;
        hid_t memspace = H5Screate_simple(ndims, dims_memspace, NULL);

        hsize_t dims_offset[5];     //ndims
        dims_offset[0] = dims[0];
        for (unsigned int i = 1; i < 5; i++)
            dims_offset[i] = 0;

        void *ptr = (void *)malloc(buffer_size * data_size);
        if (ptr == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d %d\n", __FILE__, __LINE__, (int)(buffer_size * data_size));
            return -2;
        }

        for (unsigned int type = 0; type < p->n_types; type++)
        {
            dims_offset[1] = type;
            //memcpy for the root data
            memcpy(ptr, field_pointer + ghost_buffer_size * data_size + p->n_cells_local * type * data_size,
                   buffer_size * data_size);
            for (int i = 0; i < p->args.N_domains_arg; i++)
            {
                filespace = H5Dget_space(dset);
                HDF5_ERROR_CHECK(filespace);

                status =
                        H5Sselect_hyperslab(filespace, H5S_SELECT_SET, dims_offset, NULL, dims_memspace, NULL);
                HDF5_ERROR_CHECK(status);
                status = H5Dwrite(dset, hdf5_type, memspace, filespace, plist_id, ptr);
                HDF5_ERROR_CHECK(status);

                //Recv data from the next rank
                if (i + 1 < p->args.N_domains_arg)
                {
                    const unsigned int rank = (i + 1) * p->info_MPI.domain_size;
#if ( ENABLE_MPI == 1 )         //Safe, because i+1 < p->args.N_domains_arg = 1 for non-mpi simulation
                    MPI_Recv(ptr, buffer_size, mpi_type, rank, (i + 1) + type * p->args.N_domains_arg,
                             p->info_MPI.SOMA_comm_sim, MPI_STATUS_IGNORE);
#endif                          //ENABLEXS
                    dims_offset[2] += p->nx / p->args.N_domains_arg;
                }
            }
            dims_offset[2] = 0;
        }
        free(ptr);
        status = H5Sclose(memspace);
        HDF5_ERROR_CHECK(status);

        status = H5Sclose(filespace);
        HDF5_ERROR_CHECK(status);
        status = H5Dclose(dset);
        HDF5_ERROR_CHECK(status);
    }
    else if (p->info_MPI.domain_rank == 0)
    {                       //Send data to sim rank to write the data
        for (unsigned int type = 0; type < p->n_types; type++)
        {
#if ( ENABLE_MPI == 1 )
            MPI_Send(field_pointer + ghost_buffer_size * data_size + type * p->n_cells_local * data_size,
                     buffer_size, mpi_type, 0, my_domain + type * p->args.N_domains_arg,
                     p->info_MPI.SOMA_comm_sim);
#endif                          //ENABLE_MPI
        }
    }

    return 0;
}
int extent_ana_by_field(const soma_scalar_t * const data, const uint64_t n_data, const char *const name,
                        const hid_t file_id)
{
    herr_t status;
    hid_t dset = H5Dopen(file_id, name, H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset);
    hid_t d_space = H5Dget_space(dset);
    HDF5_ERROR_CHECK(d_space);
    const unsigned int ndims = H5Sget_simple_extent_ndims(d_space);
    if (ndims != 2)
    {
        fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                __FILE__, __LINE__, name);
        return -1;
    }
    hsize_t dims[2];            //ndims

    status = H5Sget_simple_extent_dims(d_space, dims, NULL);
    HDF5_ERROR_CHECK(status);

    hsize_t dims_new[2];        //ndims
    dims_new[0] = dims[0] + 1;
    dims_new[1] = n_data;

    status = H5Dset_extent(dset, dims_new);
    HDF5_ERROR_CHECK(status);

    status = H5Sclose(d_space);
    HDF5_ERROR_CHECK(status);
    hid_t filespace = H5Dget_space(dset);
    HDF5_ERROR_CHECK(filespace);

    hsize_t dims_memspace[2];   //ndims
    dims_memspace[0] = dims_new[0] - dims[0];
    for (unsigned int i = 1; i < ndims; i++)
    {
        dims_memspace[i] = dims[i];
    }

    hid_t memspace = H5Screate_simple(ndims, dims_memspace, NULL);
    hsize_t dims_offset[2];     //ndims
    dims_offset[0] = dims[0];
    for (unsigned int i = 1; i < ndims; i++)
        dims_offset[i] = 0;

    if (dims[0] > 0)
    {
        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, dims_offset, NULL, dims_memspace, NULL);
        HDF5_ERROR_CHECK(status);
    }
    else
        memspace = H5S_ALL;

    status = H5Dwrite(dset, H5T_SOMA_NATIVE_SCALAR, memspace, filespace, H5P_DEFAULT, data);
    HDF5_ERROR_CHECK(status);

    status = H5Sclose(filespace);
    HDF5_ERROR_CHECK(status);
    if (memspace != H5S_ALL)
    {
        status = H5Sclose(memspace);
        HDF5_ERROR_CHECK(status);
    }
    status = H5Dclose(dset);
    HDF5_ERROR_CHECK(status);
    return 0;
}


int extent_structure(const struct Phase *p, const soma_scalar_t * const data, const char *const name,
                     const hid_t file_id, const enum structure_factor_type sf_type)
{
    unsigned int size;

    switch (sf_type)
    {
        case DYNAMICAL_STRUCTURE_FACTOR:
            size = p->ana_info.q_size_dynamical;
            break;
        case STATIC_STRUCTURE_FACTOR:
            size = p->ana_info.q_size_static;
            break;
        default:
            fprintf(stderr, "ERROR: Not correct structure_factor_type. %s:%s:%d\n", __func__, __FILE__, __LINE__);
            return -1;
    }

    herr_t status;
    hid_t dset = H5Dopen(file_id, name, H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset);
    hid_t d_space = H5Dget_space(dset);
    HDF5_ERROR_CHECK(d_space);
    const unsigned int ndims = H5Sget_simple_extent_ndims(d_space);
    if (ndims != 4)
    {
        fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                __FILE__, __LINE__, name);
        return -1;
    }
    hsize_t dims[4];            //ndims, no malloc because!

    status = H5Sget_simple_extent_dims(d_space, dims, NULL);
    HDF5_ERROR_CHECK(status);

    hsize_t dims_new[4];        //ndims
    dims_new[0] = dims[0] + 1;
    dims_new[1] = size;
    assert(dims[1] == size);
    dims_new[2] = p->n_poly_type;
    assert(dims[2] == p->n_poly_type);
    dims_new[3] = p->n_types * p->n_types;
    assert(dims[3] == p->n_types * p->n_types);
    status = H5Dset_extent(dset, dims_new);
    HDF5_ERROR_CHECK(status);

    status = H5Sclose(d_space);
    HDF5_ERROR_CHECK(status);
    hid_t filespace = H5Dget_space(dset);
    HDF5_ERROR_CHECK(filespace);

    hsize_t dims_memspace[4];   //ndims
    dims_memspace[0] = dims_new[0] - dims[0];
    for (unsigned int i = 1; i < ndims; i++)
    {
        dims_memspace[i] = dims[i];
    }

    hid_t memspace = H5Screate_simple(ndims, dims_memspace, NULL);
    hsize_t dims_offset[4];     //ndims
    dims_offset[0] = dims[0];
    for (unsigned int i = 1; i < ndims; i++)
        dims_offset[i] = 0;

    if (dims[0] > 0)
    {
        status = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, dims_offset, NULL, dims_memspace, NULL);
        HDF5_ERROR_CHECK(status);
    }
    else
        memspace = H5P_DEFAULT;

    status = H5Dwrite(dset, H5T_SOMA_NATIVE_SCALAR, memspace, filespace, H5P_DEFAULT, data);
    HDF5_ERROR_CHECK(status);

    status = H5Sclose(filespace);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset);
    HDF5_ERROR_CHECK(status);
    return 0;
}
