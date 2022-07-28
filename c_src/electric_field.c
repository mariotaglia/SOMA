/* Copyright (C) 2016-2021 Ludwig Schneider

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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include "soma_config.h"
#include "phase.h"
#include "mesh.h"
#include "io.h"
#include "electric_field.h"

//! \file electric_field.c
//! \brief Implementation of electric_field.h
int read_electric_field_hdf5(struct Phase *const p, const hid_t file_id, const hid_t plist_id)
{
    p->ef.eps = NULL;
    p->ef.eps_arr = NULL;
    p->ef.electrodes = NULL;
    p->ef.iter_per_MC = 0;
    p->ef.iter_limit = 0;
    p->ef.thresh_iter = 0.0;
    p->ef.Epot = NULL;
    p->ef.Epot_tmp = NULL;
    p->ef.pre_deriv = NULL;
    p->ef.H_el_field = NULL;
    p->ef.H_el = 0.0;
    p->ef.E_field = NULL;
    p->ef.omega_field_el = NULL;
    //p->ef.sqrt_Nbar = 0.0;
    hid_t status;
    
    //Quick exit if no electric field is present in the file
    if (!(H5Lexists(file_id, "/epot_field", H5P_DEFAULT) > 0))
        return 0;

    //Quick exit if no electrodes are present in the file
    if (!(H5Lexists(file_id, "/electrodes", H5P_DEFAULT) > 0))
        return 0;

    //Quick exit if no dielectric constants are specified
    if (!(H5Lexists(file_id, "/parameter/dielectric_constants", H5P_DEFAULT) > 0))
        return 0;

    //Allocate ef.eps_arr
    p->ef.eps_arr =
        (soma_scalar_t *) malloc((p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz *
                           sizeof(soma_scalar_t));
    if (p->ef.eps_arr == NULL)
    {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
    }

    //Allocate ef.Epot_tmp
    p->ef.Epot_tmp =
        (soma_scalar_t *) malloc((p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz *
                           sizeof(soma_scalar_t));
    if (p->ef.Epot_tmp == NULL)
    {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
    }

    //Allocate partial derivatives array
    p->ef.pre_deriv =
        (soma_scalar_t *) malloc((p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz * 6 *
                          sizeof(soma_scalar_t));
    if (p->ef.pre_deriv == NULL)
    {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
    }

    //Allocate ef.H_el_field
    p->ef.H_el_field =
        (soma_scalar_t *) malloc((p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz *
                           sizeof(soma_scalar_t));
    if (p->ef.H_el_field == NULL)
    {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
    }

    //Allocate ef.E_field
    p->ef.E_field =
        (soma_scalar_t *) malloc((p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz *
                           sizeof(soma_scalar_t));
    if (p->ef.E_field == NULL)
    {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
    }

    //Allocate ef.omega_field_el
    p->ef.omega_field_el =
        (soma_scalar_t *) malloc((p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz * p->n_types *
                           sizeof(soma_scalar_t));
    if (p->ef.omega_field_el == NULL)
    {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
    }

    //Read p->ef.eps
    p->ef.eps = (soma_scalar_t * const)malloc(p->n_types * sizeof(soma_scalar_t));
    if (p->ef.eps == NULL)
    {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
    }
    status = read_hdf5(file_id, "/parameter/dielectric_constants", H5T_SOMA_NATIVE_SCALAR, plist_id, p->ef.eps);
    HDF5_ERROR_CHECK2(status, "/parameter/dielectric_constants");

    //Read p->ef.iter_per_MC
    status = read_hdf5(file_id, "/parameter/ef_iter_per_MC", H5T_NATIVE_UINT, plist_id, &(p->ef.iter_per_MC));
    HDF5_ERROR_CHECK2(status, "/parameter/ef_iter_per_MC");

    //Read p->ef.iter_limit
    status = read_hdf5(file_id, "/parameter/ef_iter_limit", H5T_NATIVE_UINT, plist_id, &(p->ef.iter_limit));
    HDF5_ERROR_CHECK2(status, "/parameter/ef_iter_limit");

    //Read p->ef.thresh_iter
    status = read_hdf5(file_id, "/parameter/ef_thresh_iter", H5T_SOMA_NATIVE_SCALAR, plist_id, &(p->ef.thresh_iter));
    HDF5_ERROR_CHECK2(status, "/parameter/ef_thresh_iter");

    /* //Read p->ef.sqrt_Nbar */
    /* status = read_hdf5(file_id, "/parameter/sqrt_Nbar", H5T_SOMA_NATIVE_SCALAR, plist_id, &(p->ef.sqrt_Nbar)); */
    /* HDF5_ERROR_CHECK2(status, "/parameter/sqrt_Nbar"); */

    //Read electrode array
    status = read_field_custom_hdf5(p, (void **) &(p->ef.electrodes), "/electrodes", sizeof(uint8_t), H5T_STD_U8LE, MPI_UINT8_T, file_id, plist_id);
    if (status != 0)
    {
        fprintf(stderr, "ERROR: failed to read electrode field %s:%d.\n", __FILE__, __LINE__);
        return status;
    }                                                                             // Substitute with function in io.c                                                                                                             

    //Read electric potential field
    status = read_field_custom_hdf5(p, (void **) &(p->ef.Epot), "/epot_field", sizeof(soma_scalar_t), H5T_SOMA_NATIVE_SCALAR,
                                            MPI_SOMA_SCALAR, file_id, plist_id);
    if (status != 0)
    {
        fprintf(stderr, "ERROR: failed to read electric potential field %s:%d.\n", __FILE__, __LINE__);
        return status;
    }

    return 0;   
}



int write_electric_field_hdf5(const struct Phase *const p, const hid_t file_id, const hid_t plist_id)
{
    //Quick exit for no electricostatic hamiltonian specified
    if (p->ef.H_el_field == NULL)
        return 0;
    //Write dielectric constants data
    const hsize_t one = 1;
    hsize_t n_types_size = p->n_types;
    hid_t status;
    status =
        write_hdf5(1, &n_types_size, file_id, "/parameter/dielectric_constants", H5T_SOMA_FILE_SCALAR, H5T_SOMA_NATIVE_SCALAR, plist_id,
                   p->ef.eps);
    HDF5_ERROR_CHECK2(status, "/parameter/dielectric_constants");

    //Write amount of iterations per MC step
    status =
        write_hdf5(1, &one, file_id, "/parameter/ef_iter_per_MC", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id, &(p->ef.iter_per_MC)); // H5T_STD_U16LE
    HDF5_ERROR_CHECK2(status, "/parameter/ef_iter_per_MC");

    //Write upper limit of iterations
    status =
        write_hdf5(1, &one, file_id, "/parameter/ef_iter_limit", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id, &(p->ef.iter_limit)); // H5T_STD_U16LE
    HDF5_ERROR_CHECK2(status, "/parameter/ef_iter_limit");

    //Write iteration threshold
    status =
        write_hdf5(1, &one, file_id, "/parameter/ef_thresh_iter", H5T_SOMA_FILE_SCALAR, H5T_SOMA_NATIVE_SCALAR, plist_id, &(p->ef.thresh_iter));
    HDF5_ERROR_CHECK2(status, "/parameter/ef_thresh_iter");

    /* //Write sqrt_Nbar */
    /* status = */
    /*     write_hdf5(1, &one, file_id, "/parameter/sqrt_Nbar", H5T_SOMA_FILE_SCALAR, H5T_SOMA_NATIVE_SCALAR, plist_id, &(p->ef.sqrt_Nbar)); */
    /* HDF5_ERROR_CHECK2(status, "/parameter/sqrt_Nbar"); */

    //Write electrodes field
    status = write_field_custom_hdf5(p, (void **) &(p->ef.electrodes), "/electrodes", H5T_STD_U8LE, H5T_NATIVE_UINT8, file_id, plist_id);
    if (status != 0)
    {
        fprintf(stderr, "ERROR: failed to write electrode field %s:%d.\n", __FILE__, __LINE__);
        return status;
    }

    //Write electric potential field array
    status = write_field_custom_hdf5(p, (void **) &(p->ef.Epot), "/epot_field", H5T_SOMA_FILE_SCALAR, H5T_SOMA_NATIVE_SCALAR, file_id, plist_id);
    if (status != 0)
    {
        fprintf(stderr, "ERROR: failed to write electric potential field %s:%d.\n", __FILE__, __LINE__);
        return status;
    }

    //Write dielectric field array
    status = write_field_custom_hdf5(p, (void **) &(p->ef.eps_arr), "/dielectric_field", H5T_SOMA_FILE_SCALAR, H5T_SOMA_NATIVE_SCALAR, file_id, plist_id);
    if (status != 0)
    {
        fprintf(stderr, "ERROR: failed to write dielectric field %s:%d.\n", __FILE__, __LINE__);
        return status;
    }

    //Write electrostatic hamiltonian field array
    status = write_field_custom_hdf5(p, (void **) &(p->ef.H_el_field), "/H_el_field", H5T_SOMA_FILE_SCALAR, H5T_SOMA_NATIVE_SCALAR, file_id, plist_id);
    if (status != 0)
    {
        fprintf(stderr, "ERROR: failed to write dielectric field %s:%d.\n", __FILE__, __LINE__);
        return status;
    }

    return 0;
}

int copyin_electric_field(struct Phase *p)
{
    if (p->ef.iter_per_MC != 0)
        {
#ifdef _OPENACC
            //The ef struct itself is part of the phase struct and is already present of the device
#pragma acc enter data copyin(p->ef.eps[0:p->n_types])                                                   
#pragma acc enter data copyin(p->ef.eps_arr[0:p->n_cells_local])
#pragma acc enter data copyin(p->ef.electrodes[0:p->n_cells_local])
#pragma acc enter data copyin(p->ef.Epot[0:p->n_cells_local])
#pragma acc enter data copyin(p->ef.Epot_tmp[0:p->n_cells_local])
#pragma acc enter data copyin(p->ef.pre_deriv[0:p->n_cells_local*6])
#pragma acc enter data copyin(p->ef.H_el_field[0:p->n_cells_local])
#pragma acc enter data copyin(p->ef.E_field[0:p->n_cells_local])
#pragma acc enter data copyin(p->ef.omega_field_el[0:(p->n_cells_local*p->n_types)])
#endif                          //_OPENACC
        }
    return 0;
}

int copyout_electric_field(struct Phase *p)
{
    if (p->ef.iter_per_MC != 0)
        {
#ifdef _OPENACC
#pragma acc exit data copyout(p->ef.eps[0:p->n_types])
#pragma acc exit data copyout(p->ef.eps_arr[0:p->n_cells_local])
#pragma acc exit data copyout(p->ef.electrodes[0:p->n_cells_local])
#pragma acc exit data copyout(p->ef.Epot[0:p->n_cells_local])
#pragma acc exit data copyout(p->ef.Epot_tmp[0:p->n_cells_local])
#pragma acc exit data copyout(p->ef.pre_deriv[0:p->n_cells_local*6])
#pragma acc exit data copyout(p->ef.H_el_field[0:p->n_cells_local])
#pragma acc exit data copyout(p->ef.E_field[0:p->n_cells_local])
#pragma acc exit data copyout(p->ef.omega_field_el[0:(p->n_cells_local*p->n_types)])
#endif                          //_OPENACC
        }
    return 0;
}

int update_self_electric_field(const struct Phase *const p)
{
    if (p->ef.iter_per_MC != 0)
        {
#ifdef _OPENACC
#pragma acc update self(p->ef.eps[0:p->n_types])
#pragma acc update self(p->ef.eps_arr[0:p->n_cells_local])
#pragma acc update self(p->ef.electrodes[0:p->n_cells_local])
#pragma acc update self(p->ef.Epot[0:p->n_cells_local])
#pragma acc update self(p->ef.Epot_tmp[0:p->n_cells_local])
#pragma acc update self(p->ef.pre_deriv[0:p->n_cells_local*6])
#pragma acc update self(p->ef.H_el_field[0:p->n_cells_local])
#pragma acc update self(p->ef.E_field[0:p->n_cells_local])
#pragma acc update self(p->ef.omega_field_el[0:(p->n_cells_local*p->n_types)])
#endif                          //_OPENACC
        }
    return 0;
}

uint64_t cell_to_index(struct Phase *const p, const uint64_t x, const uint64_t y, const uint64_t z)         //TODO: Error if domain decomposition active
{
    int64_t xt=x; //changed from uint64_t to int64_t because otherwise xt < 0 always false
    int64_t yt=y;
    int64_t zt=z;
    if (xt >= p->nx) //Wrap back if necessary // formerly: if (xt >= (uint64) p->nx)
      xt -= p->nx;
    if (xt < 0)
      xt += p->nx;
    if (yt >= p->ny) //Wrap back if necessary
      yt -= p->ny;
    if (yt < 0)
      yt += p->ny;
    if (zt >= p->nz) //Wrap back if necessary
      zt -= p->nz;
    if (zt < 0)
      zt += p->nz;
    //Unified data layout [type][x][y][z]
    return xt * p->ny * p->nz + yt * p->nz + zt;
}

soma_scalar_t dEpotx(struct Phase *const p, soma_scalar_t *e_field, const uint64_t x, const uint64_t y, const uint64_t z)
{
    //Resolve periodic boundaries
    int64_t xp = x + 1;
    int64_t xm = x - 1;
    if (xp >= p->nx) xp -= p->nx;
    if (xm < 0) xm += p->nx;
    return (e_field[(xp) * p->ny * p->nz + y * p->nz + z] - e_field[(xm) * p->ny * p->nz + y * p->nz + z]) * 0.5 * p->nx / p->Lx;
}

soma_scalar_t dEpoty(struct Phase *const p, soma_scalar_t *e_field, const uint64_t x, const uint64_t y, const uint64_t z)
{
    int64_t yp = y + 1;
    int64_t ym = y - 1;
    if (yp >= p->ny) yp -= p->ny;
    if (ym < 0) ym += p->ny;
    return (e_field[x * p->ny * p->nz + (yp) * p->nz + z] - e_field[x * p->ny * p->nz + (ym) * p->nz + z]) * 0.5 * p->ny / p->Ly;
}

soma_scalar_t dEpotz(struct Phase *const p, soma_scalar_t *e_field, const uint64_t x, const uint64_t y, const uint64_t z)
{
    int64_t zp = z + 1;
    int64_t zm = z - 1;
    if (zp >= p->nz) zp -= p->nz;
    if (zm < 0) zm += p->nz;
    return (e_field[x * p->ny * p->nz + y * p->nz + (zp)] - e_field[x * p->ny * p->nz + y * p->nz + (zm)]) * 0.5 * p->nz / p->Lz;
}

soma_scalar_t d2Epotx(struct Phase *const p, soma_scalar_t *e_field, const uint64_t x, const uint64_t y, const uint64_t z)
{
    int64_t xp = x + 1;
    int64_t xm = x - 1;
    if (xp >= p->nx) xp -= p->nx;
    if (xm < 0) xm += p->nx;
    return (e_field[(xp) * p->ny * p->nz + y * p->nz + z] - (2 * e_field[(x) * p->ny * p->nz + y * p->nz + z]) + 
            e_field[(xm) * p->ny * p->nz + y * p->nz + z]) * (p->nx / p->Lx) * (p->nx / p->Lx);
}

soma_scalar_t d2Epoty(struct Phase *const p, soma_scalar_t *e_field, const uint64_t x, const uint64_t y, const uint64_t z)
{
    int64_t yp = y + 1;
    int64_t ym = y - 1;
    if (yp >= p->ny) yp -= p->ny;
    if (ym < 0) ym += p->ny;
    return (e_field[x * p->ny * p->nz + (yp) * p->nz + z] - (2 * e_field[x * p->ny * p->nz + (y) * p->nz + z]) +
            e_field[x * p->ny * p->nz + (ym) * p->nz + z]) * (p->ny / p->Ly) * (p->ny / p->Ly);
}

soma_scalar_t d2Epotz(struct Phase *const p, soma_scalar_t *e_field, const uint64_t x, const uint64_t y, const uint64_t z)
{
    int64_t zp = z + 1;
    int64_t zm = z - 1;
    if (zp >= p->nz) zp -= p->nz;
    if (zm < 0) zm += p->nz;
    return (e_field[x * p->ny * p->nz + y * p->nz + (zp)] - (2 * e_field[x * p->ny * p->nz + y * p->nz + (z)]) +
            e_field[x * p->ny * p->nz + y * p->nz + (zm)]) * (p->nz / p->Lz) * (p->nz / p->Lz);
}

/* void calc_sqrt_Nbar(struct Phase *const p) */
/* { */
/*     uint16_t n_mono_t; */
/*     uint64_t sum_chains_t_scaled = 0; */
/*     uint64_t * type_num = (uint64_t*) calloc(0, p->n_poly_type * sizeof(uint64_t)); */

/*     //loop over all polymer chains, save type in type_num */
/*     for (uint64_t c=0; c < p->n_polymers; c++) */
/*         { */
/*             // Polymer * poly = &(p->polymers[c]); */
/*             uint8_t poly_type_c = p->polymers[c].type; */
/*             type_num[poly_type_c] += 1; */
/*         } */

/*     if (p->info_MPI.domain_rank == 0) */
/*     { */
/*         MPI_Reduce(MPI_IN_PLACE, type_num, p->n_poly_type, MPI_UINT8_T, MPI_SUM, 0, p->info_MPI.SOMA_comm_domain); */
/*     } */
/*     else */
/*     { */
/*         MPI_Reduce(type_num, NULL, p->n_poly_type, MPI_UINT8_T, MPI_SUM, 0, p->info_MPI.SOMA_comm_domain); */
/*     } */

/*     for (uint8_t t=0; t < p->n_poly_type; t++) */
/*     { */
/*         n_mono_t = p->poly_arch[p->poly_type_offset[t]]; */
/*         sum_chains_t_scaled += type_num[t] * n_mono_t / p->reference_Nbeads;         */
/*     } */

/*     p->ef.sqrt_Nbar = sum_chains_t_scaled * (1 / (p->Lx * p->Ly * p->Lz)); */

/* #pragma acc update device(p->ef.sqrt_Nbar) */

/*     free(type_num); */
/* } */

void calc_dielectric_field(struct Phase *const p)
{
#pragma acc parallel loop present(p[0:1])
#pragma omp parallel for
    for (uint64_t i = 0; i < p->n_cells; i++)
    {
        if (p->ef.electrodes[i] == 1 || p->area51[i] == 1)
        {
            p->ef.eps_arr[i] = 1.0;
        }
        else
        {
	    // /* Calculation according to welling2014 eq. 83, density dependent expression */
            // soma_scalar_t tmp_frac_n = 1.0;
            // soma_scalar_t tmp_eps_n = 0.0;

            // for (uint8_t n = 0; n < p->n_types; n++)
            // {
            //     tmp_frac_n -= p->fields_unified[i+n*p->n_cells_local] * p->field_scaling_type[n];
            //     tmp_eps_n += p->fields_unified[i+n*p->n_cells_local] * p->field_scaling_type[n] * p->ef.eps[n];
            // }
	    // //eps0 = 1.0
            // p->ef.eps_arr[i] = 1.0 * tmp_frac_n + tmp_eps_n;

            /* Calculation according to welling2014 eq. 82, density independent expression */
            soma_scalar_t tmp_phi_n = 0.0;
            soma_scalar_t tmp_eps_phi_n = 0.0;
            
            for (uint32_t n = 0; n < p->n_types; n++)
            {
                tmp_phi_n += p->fields_unified[i+n*p->n_cells_local] * p->field_scaling_type[n];
                tmp_eps_phi_n += p->fields_unified[i+n*p->n_cells_local] * p->field_scaling_type[n] * p->ef.eps[n];
            }

            p->ef.eps_arr[i] = 1.0 * tmp_eps_phi_n / tmp_phi_n;
        }
    }
}

void pre_derivatives(struct Phase *const p)
{
#pragma acc parallel loop present(p[0:1]) collapse(3) deviceptr(p->ef.eps_arr)
#pragma omp parallel for collapse(3)
    for (uint64_t x=0; x < p->nx; x++)
        for (uint64_t y=0; y < p->ny; y++)
            for (uint64_t z=0; z < p->nz; z++)
            {
                uint64_t i = cell_to_index(p,x,y,z);
                uint64_t pos6 = (x * p->ny * p->nz + y * p->nz + z) * 6;
                if (p->ef.electrodes[i] == 1 || p->area51[i] == 1)
                {
                    p->ef.pre_deriv[pos6 + 0] = 1.0/6.0;
                    p->ef.pre_deriv[pos6 + 1] = 1.0/6.0;
                    p->ef.pre_deriv[pos6 + 2] = 1.0/6.0;
                    p->ef.pre_deriv[pos6 + 3] = 1.0/6.0;
                    p->ef.pre_deriv[pos6 + 4] = 1.0/6.0;
                    p->ef.pre_deriv[pos6 + 5] = 1.0/6.0;
                }
                else
                {
                    soma_scalar_t cr = 1.0/(24.0 * p->ef.eps_arr[i]);

                    soma_scalar_t epsx = cr * ( p->ef.eps_arr[cell_to_index(p,x+1,y,z)] - p->ef.eps_arr[cell_to_index(p,x-1,y,z)] );
                    soma_scalar_t epsy = cr * ( p->ef.eps_arr[cell_to_index(p,x,y+1,z)] - p->ef.eps_arr[cell_to_index(p,x,y-1,z)] );
                    soma_scalar_t epsz = cr * ( p->ef.eps_arr[cell_to_index(p,x,y,z+1)] - p->ef.eps_arr[cell_to_index(p,x,y,z-1)] );
                    
                    p->ef.pre_deriv[pos6 + 0] = epsx + 1.0/6.0;
                    p->ef.pre_deriv[pos6 + 1] =-epsx + 1.0/6.0;
                    p->ef.pre_deriv[pos6 + 2] = epsy + 1.0/6.0;
                    p->ef.pre_deriv[pos6 + 3] =-epsy + 1.0/6.0;
                    p->ef.pre_deriv[pos6 + 4] = epsz + 1.0/6.0;
                    p->ef.pre_deriv[pos6 + 5] =-epsz + 1.0/6.0;
                }
	    }
}

soma_scalar_t iterate_field(struct Phase *const p)
{
    soma_scalar_t max_i = 10.0;
    uint64_t k = 0;

    /* printf("\n area51[2247]: %d\n",p->area51[2247]); */
/*     printf("\n p->ef.electrodes[2247]: %d\n",p->ef.electrodes[2247]); */
/*     printf("\n p->ef.eps_arr[2247]: %.3f \n", p->ef.eps_arr[2247]); */
/* #pragma acc update host (p->ef.eps_arr[0:p->n_cells]) */
/*     printf("\n p->ef.eps_arr[2247] (updated): %.3f \n", p->ef.eps_arr[2247]); */
/*     printf("\n p->ef.pre_deriv[2247*6]: %.3f \n", p->ef.pre_deriv[2247*6]); */
/* #pragma acc update host (p->ef.pre_deriv[0:p->n_cells*6]) */
/*     printf("\n p->ef.pre_deriv[2247*6] (updated): %.3f \n", p->ef.pre_deriv[2247*6]); */
    
/*     printf("\n thresh iter: %.3f\n", p->ef.thresh_iter); */

    /* printf("\n=======\n"); */
/*     printf("\n p->ef.Epot[2247] before iter.: %.3f \n", p->ef.Epot[2247]); */
/* #pragma acc update host (p->ef.Epot[0:p->n_cells]) */
/*     printf("\n p->ef.Epot[2247] (updated) before iter.: %.3f \n", p->ef.Epot[2247]); */

    while(max_i > p->ef.thresh_iter && k < p->ef.iter_limit)
    {
        max_i = p->ef.thresh_iter;
#pragma acc parallel loop present(p[0:1]) collapse(3) deviceptr(p->ef.electrodes,p->ef.Epot,p->ef.Epot_tmp,p->ef.pre_deriv)
#pragma omp parallel for collapse(3) reduction(max: max_i)
        for (uint64_t x=0; x < p->nx; x++)
            for (uint64_t y=0; y < p->ny; y++)
                for (uint64_t z=0; z < p->nz; z++)
                {
                    uint64_t i = cell_to_index(p,x,y,z);
                    /* Keep electric potential field constant at electrode position */
                    if (p->ef.electrodes[i] == 1)
                    {
                        p->ef.Epot_tmp[i] = p->ef.Epot[i];
                    }
                    else
                    {
                        uint64_t pos6 = (x * p->ny * p->nz + y * p->nz + z) * 6;
                        p->ef.Epot_tmp[i] =
                        p->ef.pre_deriv[pos6+0] * p->ef.Epot[cell_to_index(p,x+1,y,z)] + p->ef.pre_deriv[pos6+1] * p->ef.Epot[cell_to_index(p,x-1,y,z)] +
                        p->ef.pre_deriv[pos6+2] * p->ef.Epot[cell_to_index(p,x,y+1,z)] + p->ef.pre_deriv[pos6+3] * p->ef.Epot[cell_to_index(p,x,y-1,z)] +
                        p->ef.pre_deriv[pos6+4] * p->ef.Epot[cell_to_index(p,x,y,z+1)] + p->ef.pre_deriv[pos6+5] * p->ef.Epot[cell_to_index(p,x,y,z-1)];

			/* Option without using pre-calculated derivatives (p->ef.pre_deriv) */
                        /* p->ef.Epot_tmp[i] = (1.0 / 6.0) * ( */
                        /*                      p->ef.Epot[cell_to_index(p,x+1,y,z)] + p->ef.Epot[cell_to_index(p,x-1,y,z)] + */
                        /*                      p->ef.Epot[cell_to_index(p,x,y+1,z)] + p->ef.Epot[cell_to_index(p,x,y-1,z)] + */
                        /*                      p->ef.Epot[cell_to_index(p,x,y,z+1)] + p->ef.Epot[cell_to_index(p,x,y,z-1)]) + */
                        /*                     (1.0 / (24.0 * p->ef.eps_arr[i])) * ( */
                        /*                      (p->ef.eps_arr[cell_to_index(p,x+1,y,z)] - p->ef.eps_arr[cell_to_index(p,x-1,y,z)]) * */
                        /*                      (p->ef.Epot[cell_to_index(p,x+1,y,z)] - p->ef.Epot[cell_to_index(p,x-1,y,z)]) + */
                        /*                      (p->ef.eps_arr[cell_to_index(p,x,y+1,z)] - p->ef.eps_arr[cell_to_index(p,x,y-1,z)]) * */
                        /*                      (p->ef.Epot[cell_to_index(p,x,y+1,z)] - p->ef.Epot[cell_to_index(p,x,y-1,z)]) + */
                        /*                      (p->ef.eps_arr[cell_to_index(p,x,y,z+1)] - p->ef.eps_arr[cell_to_index(p,x,y,z-1)]) * */
                        /*                      (p->ef.Epot[cell_to_index(p,x,y,z+1)] - p->ef.Epot[cell_to_index(p,x,y,z-1)]) ); */

                        
			/* max_i = fmax(max_i, fabs(p->ef.Epot_tmp[i] - p->ef.Epot[i])); */
		    }
                }
	/* additional loop necessary for conv_crit since loop iterations could lead to data dependency using the jacobi stencil */
	/* could be included in former loop, however, might lead to incorrect/inaccurate results - better trade of than looping again? */
	/* OR: just use the difference between new and old solution (fabs(p->ef.Epot_tmp[i] - p->ef.Epot[i]) to trace conversion & omit this loop */
#pragma acc parallel loop present(p[0:1]) collapse(3) reduction(max: max_i) deviceptr(p->ef.Epot_tmp,p->ef.eps_arr)
#pragma omp parallel for collapse(3)
        for (uint64_t x=0; x < p->nx; x++)
            for (uint64_t y=0; y < p->ny; y++)
                for (uint64_t z=0; z < p->nz; z++)
		    {
		      uint64_t j = cell_to_index(p,x,y,z);

			if(p->ef.electrodes[j] != 1)
			{
			    // Laplace equation, welling2017 eq. 5
			    soma_scalar_t conv_crit = (d2Epotx(p,p->ef.Epot_tmp,x,y,z) +
						       d2Epoty(p,p->ef.Epot_tmp,x,y,z) +
						       d2Epotz(p,p->ef.Epot_tmp,x,y,z)) * p->ef.eps_arr[j] +
			    			      (dEpotx(p,p->ef.Epot_tmp,x,y,z) * dEpotx(p,p->ef.eps_arr,x,y,z) +
						       dEpoty(p,p->ef.Epot_tmp,x,y,z) * dEpoty(p,p->ef.eps_arr,x,y,z) +
						       dEpotz(p,p->ef.Epot_tmp,x,y,z) * dEpotz(p,p->ef.eps_arr,x,y,z));

			    if (max_i < conv_crit) max_i = conv_crit;
			}
		    }
	
#pragma acc parallel loop present(p[0:1]) collapse(3) deviceptr(p->ef.Epot,p->ef.Epot_tmp)
#pragma omp parallel for
        for (uint64_t x=0; x < p->nx; x++)
            for (uint64_t y=0; y < p->ny; y++)
                for (uint64_t z=0; z < p->nz; z++)
		    {
		        uint64_t j = cell_to_index(p,x,y,z);

			if(p->ef.electrodes[j] != 1)
			{
			    p->ef.Epot[j] = p->ef.Epot_tmp[j];
			}
		    }

/*	if( (k+1) % 500==0 ) printf("\n max: %.3e\n", max_i); */
/* #pragma acc data deviceptr(p->ef.Epot_tmp, p->ef.Epot) */
/* 	{ */
/* 	soma_scalar_t *tmp = p->ef.Epot_tmp; */
/* 	p->ef.Epot_tmp = p->ef.Epot; */
/* 	p->ef.Epot = tmp; */
/* 	} */
	k += 1;
    }

    /* printf("\n p->ef.Epot[2247] after iter.: %.3f \n", p->ef.Epot[2247]); */
/* #pragma acc update host (p->ef.Epot[0:p->n_cells]) */
/*     printf("\n p->ef.Epot[2247] (updated) after iter.: %.3f \n", p->ef.Epot[2247]); */

/*     soma_scalar_t max_host = 0.0; */
/*     uint64_t err_cells = 0; */
/*     for (uint64_t x=0; x < p->nx; x++) */
/*             for (uint64_t y=0; y < p->ny; y++) */
/*                 for (uint64_t z=0; z < p->nz; z++) */
/* 		{ */
/* 		    uint64_t l = cell_to_index(p,x,y,z); */
/* 		    if(p->ef.electrodes[l] != 1) */
/* 		    { */
/* 		        soma_scalar_t conv_crit = fabs( (d2Epotx(p,p->ef.Epot,x,y,z) + */
/* 						         d2Epoty(p,p->ef.Epot,x,y,z) + */
/* 						         d2Epotz(p,p->ef.Epot,x,y,z)) * p->ef.eps_arr[l] + */
/* 						        (dEpotx(p,p->ef.Epot,x,y,z) * dEpotx(p,p->ef.eps_arr,x,y,z) + */
/* 						         dEpoty(p,p->ef.Epot,x,y,z) * dEpoty(p,p->ef.eps_arr,x,y,z) + */
/* 						         dEpotz(p,p->ef.Epot,x,y,z) * dEpotz(p,p->ef.eps_arr,x,y,z)) ); */
/* 			if(max_host < conv_crit) max_host = conv_crit; */
/* 			if(conv_crit > 1e2) */
/* 		        { */
/* 		            err_cells += 1; */
/* 		            printf("\n %d,%d,%d\n",x,y,z); */
/* 		        } */
/* 		    } */
/* 		} */
/*     printf("\n number of cells in which conv_crit is too high: %d\n", err_cells); */
/*     printf("\n max of Epot field on host: %.3e\n", max_host); */
    
    if ((p->time + 1)  % 400 == 0) printf("MC step: %d, iterations to solve ef: %d \n",p->time+1, k);
    return max_i;
}

int calc_electric_field_contr(struct Phase *const p)
{
    soma_scalar_t max = p->ef.thresh_iter;
    uint64_t k = 0;
    
/*     if (p->time==0) */
/*     { */
/* #pragma acc update device(p->sqrt_Nbar) */
/*         // calc_sqrt_Nbar(p); */
/*     } */

    calc_dielectric_field(p);
#pragma acc wait
    pre_derivatives(p);
#pragma acc wait
    max = iterate_field(p);
    if (max > p->ef.thresh_iter) return -1; 
#pragma acc wait
    // soma_scalar_t sum_H_el = 0.0;
#pragma acc parallel loop present(p[0:1]) collapse(3) deviceptr(p->ef.electrodes,p->ef.Epot,p->ef.E_field,p->ef.eps,p->ef.omega_field_el) //reduction(+:sum_H_el)
#pragma omp parallel for collapse(3) //reduction(+:sum_H_el)
    for (uint64_t x=0; x < p->nx; x++)
        for (uint64_t y=0; y < p->ny; y++)
            for (uint64_t z=0; z < p->nz; z++)
            {
                uint64_t i = cell_to_index(p,x,y,z);
                if (p->ef.electrodes[i] == 1)
                {
                    continue; 
                }
                else
                {
                    soma_scalar_t dEpot_sq = dEpotx(p,p->ef.Epot,x,y,z) * dEpotx(p,p->ef.Epot,x,y,z) +
                                             dEpoty(p,p->ef.Epot,x,y,z) * dEpoty(p,p->ef.Epot,x,y,z) +
                                             dEpotz(p,p->ef.Epot,x,y,z) * dEpotz(p,p->ef.Epot,x,y,z);
                    //p->ef.H_el_field[i] = 0.5 * dEpot_sq * p->ef.eps_arr[i];                  // redundant?
                    //sum_H_el += p->ef.H_el_field[i];                                          // redundant?
                    
                    p->ef.E_field[i] = -(dEpotx(p,p->ef.Epot,x,y,z) + dEpoty(p,p->ef.Epot,x,y,z) + dEpotz(p,p->ef.Epot,x,y,z));

                    for (uint32_t m = 0; m < p->n_types; m++)
                    {   
                        soma_scalar_t phi_sum = 0.0;
                        soma_scalar_t diff_part = 0.0;
                        phi_sum += p->fields_unified[i+m*p->n_cells_local] * p->field_scaling_type[m];

                        for (uint32_t n = 0; n < p->n_types; n++)
                        {
                            if (n != m)
                            {
                                phi_sum += p->fields_unified[i+n*p->n_cells_local] * p->field_scaling_type[n];
                                diff_part += p->fields_unified[i+n*p->n_cells_local] * p->field_scaling_type[n] *
                                            (p->ef.eps[m] - p->ef.eps[n]);
                            }
                        }
                        
                        p->ef.omega_field_el[i+m*p->n_cells_local] = 1.0 * 0.5 * diff_part / (phi_sum * phi_sum * p->sqrt_Nbar) * dEpot_sq;

                        // eps_res = 0.0;
                        // phi_other = 0.0;
                        // for (uint8_t n = 0; n < p->n_types; n++)
                        // {
                        //     if (n != m)
                        //     {
                        //         eps_res += p->ef.eps[n];
                        //         phi_other += p->fields_unified[i+n*p->n_cells_local] * p->field_scaling_type[n];
                        //     }
                        // }

                        // phi_self = p->fields_unified[i+m*p->n_cells_local] * p->field_scaling_type[m];
                        // eps_diff = (p->ef.eps[m] - eps_res) / (p->n_types * p->num_all_beads) *
                        //            (phi_other) / ((phi_self + phi_other) * (phi_self + phi_other)); 
                        
                        // p->ef.omega_field_el[i+m*p->n_cells_local] = eps_diff * dEpot_sq;
                    }
                }

            }

    //tests(p,k);
    //p->ef.H_el = sum_H_el;
#pragma acc update self(p->ef.omega_field_el)
    if (p->ef.omega_field_el != p->ef.omega_field_el)
    {
        return -1;
    }
    //printf("\n p->ef.omega_field_el[2247*2]: %.3e", p->ef.omega_field_el[2247*2]);
    //printf("\n sqrt Nbar = %.3e \n", p->sqrt_Nbar);
    // if domain_decomposition -> allreduce Epot, E_field, eps_arr here.
    return 0;
}

int free_electric_field(struct Phase *const p)
{
    free(p->ef.eps_arr);
    free(p->ef.Epot_tmp);
    free(p->ef.pre_deriv);
    free(p->ef.H_el_field);
    free(p->ef.eps);
    free(p->ef.electrodes);
    free(p->ef.Epot);
    free(p->ef.E_field);
    free(p->ef.omega_field_el);

    return 0;
}

void tests(struct Phase *const p,uint64_t k)
{
    uint32_t counter = 1;
    update_self_electric_field(p);
    if (p->time % counter == 0)
    {
        uint64_t sum_electrodes = 0;
        soma_scalar_t sum_Epot = 0.0;
        soma_scalar_t sum_eps_arr = 0.0;
        soma_scalar_t sum_pre_deriv = 0.0;
        soma_scalar_t sum_omega_field_el = 0.0;
        soma_scalar_t sum_omega_field = 0.0;

#pragma acc update host(p->omega_field_unified)

        for (uint64_t o = 0; o < p->n_cells; o++)
        {
            sum_electrodes += p->ef.electrodes[o];
            sum_Epot += p->ef.Epot[o];
            sum_eps_arr += p->ef.eps_arr[o];
            for (uint8_t m = 0; m < p->n_types; m++)
            {
                sum_omega_field_el += (p->ef.omega_field_el[o+m*p->n_types] * p->ef.omega_field_el[o+m*p->n_types]);
                sum_omega_field += (p->omega_field_unified[o+m*p->n_types] * p->omega_field_unified[o+m*p->n_types]);
            }
            for (uint8_t l = 0; l < 6; l++)
                sum_pre_deriv += p->ef.pre_deriv[o+l];           
        }

        // printf("=====\n");
        // printf("MC STEP: %d\n", p->time);
        // printf("Iterations per MC: %ld\n",p->ef.iter_per_MC);
        // printf("Iteration limit: %ld\n",p->ef.iter_limit);
        // printf("Iteration threshold: %.5f\n",p->ef.thresh_iter);
        // printf("Sum electrode array: %ld\n",sum_electrodes);
        // printf("Sum Epot array: %.3f\n",sum_Epot);
        // printf("Sum epsilon array: %.3f\n",sum_eps_arr);
        // printf("Sum pre deriv array: %.3f\n",sum_pre_deriv);
        // printf("Sum el. hamiltonian: %.3f\n",p->ef.H_el);
        printf("Sum omega_field sq.: \t%.3e\n", sum_omega_field);
        printf("Sum omega_field_el sq.: %.3e\n",sum_omega_field_el);
        printf("=====\n");
    }
}
