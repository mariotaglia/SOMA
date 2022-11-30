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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include "soma_config.h"
#include "phase.h"
#include "mesh.h"
#include "io.h"
#include "electric_field.h"
#include "electric_field_helper.h"

//! \file electric_field.c
//! \brief Implementation of electric_field.h
int read_electric_field_hdf5(struct Phase *const p, const hid_t file_id, const hid_t plist_id)
{
    p->ef.eps = NULL;
    p->ef.eps_arr = NULL;
    p->ef.electrodes = NULL;
    p->ef.iter_limit = 0;
    p->ef.thresh_iter = 0.0;
    p->ef.amt_iter = 0;
    p->ef.Epot = NULL;
    p->ef.Epot_tmp = NULL;
    p->ef.pre_deriv = NULL;
    p->ef.H_el_field = NULL;
    p->ef.H_el = 0.0;
    p->ef.E_field = NULL;
    p->ef.omega_field_el = NULL;
    //p->ef.sqrt_Nbar = 0.0;
    p->ef.stride = 0;
    p->ef.kernel = NULL;
    p->ef.kernel_norm_field = NULL;
    p->ef.kernel_blur = NULL;
    p->ef.kernel_dim = 0;
    p->ef.kernel_rad = 0;
    p->ef.kernel_sigma = 0.0;
    p->ef.conv_nx = 0;
    p->ef.conv_ny = 0;
    p->ef.conv_nz = 0;
    p->ef.n_cells_conv = 0;
    p->ef.el_pos_xy = false;
    p->ef.el_pos_xz = false;
    p->ef.el_pos_yz = false;
    p->ef.x_offset = 0;
    p->ef.y_offset = 0;
    p->ef.z_offset = 0;
    p->ef.electrodes_conv = NULL;
    p->ef.eps_arr_conv = NULL;
    p->ef.pre_deriv_conv = NULL;
    p->ef.Epot_conv = NULL;
    p->ef.Epot_tmp_conv = NULL;

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

    //Read p->ef.iter_limit
    status = read_hdf5(file_id, "/parameter/ef_iter_limit", H5T_NATIVE_UINT32, plist_id, &(p->ef.iter_limit));
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

    //Read kernel dim for convolution
    // p->ef.kernel_dim = 3; // TO DO: READ FROM XML -- MUST BE UNEVEN, default=1!
    status = read_hdf5(file_id, "/parameter/ef_kernel_dim", H5T_NATIVE_UINT8, plist_id, &(p->ef.kernel_dim));
    HDF5_ERROR_CHECK2(status, "/parameter/ef_kernel_dim");

    //Calculate kernel radius
    p->ef.kernel_rad = p->ef.kernel_dim / 2;

    //Read std dev for kernel
    // p->ef.kernel_sigma = 1.0; // TO DO: READ FROM XML, default=1.0!
    status = read_hdf5(file_id, "/parameter/ef_kernel_sigma", H5T_SOMA_NATIVE_SCALAR, plist_id, &(p->ef.kernel_sigma));
    HDF5_ERROR_CHECK2(status, "/parameter/ef_kernel_sigma");

    //Read kernel stride for convolution
    // p->ef.stride = 2; 
    status = read_hdf5(file_id, "/parameter/ef_stride", H5T_NATIVE_UINT8, plist_id, &(p->ef.stride));
    HDF5_ERROR_CHECK2(status, "/parameter/ef_stride");

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

    //Write upper limit of iterations
    status =
        write_hdf5(1, &one, file_id, "/parameter/ef_iter_limit", H5T_STD_U32LE, H5T_NATIVE_UINT32, plist_id, &(p->ef.iter_limit)); // H5T_STD_U16LE
    HDF5_ERROR_CHECK2(status, "/parameter/ef_iter_limit");

    //Write iteration threshold
    status =
        write_hdf5(1, &one, file_id, "/parameter/ef_thresh_iter", H5T_SOMA_FILE_SCALAR, H5T_SOMA_NATIVE_SCALAR, plist_id, &(p->ef.thresh_iter));
    HDF5_ERROR_CHECK2(status, "/parameter/ef_thresh_iter");

    //Write kernel dim
    status =
        write_hdf5(1, &one, file_id, "/parameter/ef_kernel_dim", H5T_STD_U8LE, H5T_NATIVE_UINT8, plist_id, &(p->ef.kernel_dim));
    HDF5_ERROR_CHECK2(status, "/parameter/ef_kernel_dim");

    //Write kernel sigma
    status =
        write_hdf5(1, &one, file_id, "/parameter/ef_kernel_sigma", H5T_SOMA_FILE_SCALAR, H5T_SOMA_NATIVE_SCALAR, plist_id, &(p->ef.kernel_sigma));
    HDF5_ERROR_CHECK2(status, "/parameter/ef_kernel_sigma");

    //Write kernel stride
    status =
        write_hdf5(1, &one, file_id, "/parameter/ef_stride", H5T_STD_U8LE, H5T_NATIVE_UINT8, plist_id, &(p->ef.stride));
    HDF5_ERROR_CHECK2(status, "/parameter/ef_stride");

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

    //TO DO: WRITE CONV VARIABLES/VALUES!!!!!

    return 0;
}

int init_efield(struct Phase *const p)
{
    // determine on which outer surfaces (xy,xz,yz) electrodes are present
    // necessary to resolve non-periodic boundaries
    uint64_t electr_sum_xy_red = 0;
    uint64_t electr_sum_xz_red = 0;
    uint64_t electr_sum_yz_red = 0;

    for (uint64_t x=0; x < p->nx; x++)
        for (uint64_t y=0; y < p->ny; y++)
        {
            if ((x > 0 && x < (p->nx-1)) && (y > 0 && y < (p->ny-1)))
            {
                electr_sum_xy_red += p->ef.electrodes[x * p->ny * p->nz + y * p->nz + 0];
            }
        }

    for (uint64_t x=0; x < p->nx; x++)
        for (uint64_t z=0; z < p->nz; z++)
        {
            if ((x > 0 && x < (p->nx-1)) && (z > 0 && z < (p->nz-1)))
            {
                electr_sum_xz_red += p->ef.electrodes[x * p->ny * p->nz + 0 * p->nz + z];
            }                
        }

    for (uint64_t y=0; y < p->ny; y++)
        for (uint64_t z=0; z < p->nz; z++)
        {
            if ((y > 0 && y < (p->ny-1)) && (z > 0 && z < (p->nz-1)))
            {
                electr_sum_yz_red += p->ef.electrodes[0 * p->ny * p->nz + y * p->nz + z];
            }               
        }
    // verify electrode on opposite surfaces
    if (electr_sum_xy_red > 1) p->ef.el_pos_xy = true;
    if (electr_sum_xz_red > 1) p->ef.el_pos_xz = true;
    if (electr_sum_yz_red > 1) p->ef.el_pos_yz = true;

    if ((p->ef.el_pos_xy && p->ef.el_pos_xz) ||
        (p->ef.el_pos_xy && p->ef.el_pos_yz) || 
        (p->ef.el_pos_xz && p->ef.el_pos_yz) ||
        (!p->ef.el_pos_xy && !p->ef.el_pos_xz && !p->ef.el_pos_yz))
    {
        fprintf(stderr, "ERROR: Electrode position could not be clearly assigned to one plane %s:%d.\n", __FILE__, __LINE__);
        return -1;
    }

    return 0;
}

int init_kernel(struct Phase *const p)
{
    if (p->ef.stride > p->nx || p->ef.stride > p->ny || p->ef.stride > p->nz)
    {
        fprintf(stderr, "ERROR: Kernel stride exceeds box dimension %s:%d.\n", __FILE__, __LINE__);
        return -1;
    }

    // allocate kernel array
    p->ef.kernel = malloc((p->ef.kernel_dim * p->ef.kernel_dim * p->ef.kernel_dim * sizeof(soma_scalar_t)));                          
    if (p->ef.kernel == NULL)                                                                                                         
    {                                                                                                                                 
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);                                                                 
        return -1;                                                                                                                    
    }

    // allocate kernel_blur array always 3x3x3)
    p->ef.kernel_blur = malloc((3 * 3 * 3 * sizeof(soma_scalar_t)));                          
    if (p->ef.kernel_blur == NULL)                                                                                                         
    {                                                                                                                                 
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);                                                                 
        return -1;                                                                                                                    
    }

    // 3D gaussian filter: (1/(sqrt(2*pi)*sigma))^3 * exp(-(x^2+y^2+z^2)/(2*sigma^2))
    soma_scalar_t pi = 3.14159265359;
    soma_scalar_t prefactor = sqrt(2.0 * pi) * sqrt(2.0 * pi) * sqrt(2.0 * pi) *
                              p->ef.kernel_sigma * p->ef.kernel_sigma * p->ef.kernel_sigma;
    soma_scalar_t exp_term = 0.0;
    soma_scalar_t norm_sum = 0.0;
    uint16_t kernel_ind = 0;

    for (int8_t x = -p->ef.kernel_rad; x <= p->ef.kernel_rad; x++)
        for (int8_t y = -p->ef.kernel_rad; y <= p->ef.kernel_rad; y++)
            for (int8_t z = -p->ef.kernel_rad; z <= p->ef.kernel_rad; z++)
            {
                exp_term = x * x + y * y + z * z;

                p->ef.kernel[kernel_ind] = exp(-exp_term / (2 * p->ef.kernel_sigma * p->ef.kernel_sigma)) / prefactor;
                norm_sum += p->ef.kernel[kernel_ind];
                kernel_ind += 1;
            }

    uint16_t kernel_err = 0;
    // normalize kernel
    for (uint16_t i = 0; i < (p->ef.kernel_dim * p->ef.kernel_dim * p->ef.kernel_dim); i++)
    {
        p->ef.kernel[i] /= norm_sum;
        if (p->ef.kernel[i] != p->ef.kernel[i]) kernel_err += 1;
    }

    if (kernel_err != 0)
    {
        fprintf(stderr, "ERROR: Kernel initialization: %d. %s:%d.\n", kernel_err, __FILE__, __LINE__);
    }

    // allocate normalization kernel array
    p->ef.kernel_norm_field = (soma_scalar_t *) malloc(p->nx * p->ny * p->nz * sizeof(soma_scalar_t));         
    if (p->ef.kernel_norm_field == NULL)                                                                                                         
    {                                                                                                                                 
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);                                                                 
        return -1;                                                                                                                    
    }

    // fill normalization field
    for (uint64_t i = 0; i < (p->nx * p->ny * p->nz); i++)
    {
        p->ef.kernel_norm_field[i] = 0.0;
    }

    // fill blur kernel with ones and normalize
    for (uint16_t i = 0; i < (3 * 3 * 3); i++)
    {
        p->ef.kernel_blur[i] = 1.0 / (3 * 3 * 3);
    }

    return 0;
}

int init_convolution(struct Phase *const p)
{
    // DISCLAIMER: Convolution only works for two electrodes opposing electrodes located on outer planes of the box
    // which is resolved in 'ini_efield()'
    
    // get axes in between electrodes to allocate arrays
    if (p->ef.el_pos_yz)
        p->ef.conv_nx = (p->nx-2) / p->ef.stride;
    else
        p->ef.conv_nx = p->nx / p->ef.stride;

    if (p->ef.el_pos_xz)
        p->ef.conv_ny = (p->ny-2) / p->ef.stride;
    else
        p->ef.conv_ny = p->ny / p->ef.stride;

    if (p->ef.el_pos_xy)
        p->ef.conv_nz = (p->nz-2) / p->ef.stride;
    else
        p->ef.conv_nz = p->nz / p->ef.stride;   

    // calculate max. index that can be calculated from convoluted array
    // (max conv cells - offset for index) * stride + (stride - 1)
    uint64_t max_conv_nx = (p->ef.conv_nx - 1) * p->ef.stride + p->ef.stride - 1;
    uint64_t max_conv_ny = (p->ef.conv_ny - 1) * p->ef.stride + p->ef.stride - 1;
    uint64_t max_conv_nz = (p->ef.conv_nz - 1) * p->ef.stride + p->ef.stride - 1;

    // compare to max index of original array (minus electrodes), offset start of indexing
    uint64_t max_nx = p->nx - 1;
    uint64_t max_ny = p->ny - 1;
    uint64_t max_nz = p->nz - 1;

    // in case of electrode, offset index accordingly
    if (p->ef.el_pos_yz)
    {       
        max_conv_nx += 1; 
        max_nx -= 1;
    }
    if (p->ef.el_pos_xz)
    {       
        max_conv_ny += 1; 
        max_ny -= 1;
    }
    if (p->ef.el_pos_xy)
    {       
        max_conv_nz += 1; 
        max_nz -= 1;
    }

    // add additional cell to convoluted array in case max index of original cannot be reached
    if (max_conv_nx < max_nx) p->ef.conv_nx += 1;
    if (max_conv_ny < max_ny) p->ef.conv_ny += 1;
    if (max_conv_nz < max_nz) p->ef.conv_nz += 1;

    // alloc p->ef.eps_arr_conv, p->ef.pre_deriv_conv, p->ef.Epot_conv, p->ef.Epot_tmp_conv, p->ef.electrodes_conv
    // depending on electrode position; include two extra rows of cells to reattach electrodes, else iteration is impossible
    if (p->ef.el_pos_yz)
    {
        p->ef.electrodes_conv =
        (soma_scalar_t *) malloc(((p->ef.conv_nx + 2) / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                                   p->ef.conv_ny * p->ef.conv_nz * sizeof(soma_scalar_t));
        p->ef.eps_arr_conv =
        (soma_scalar_t *) malloc(((p->ef.conv_nx + 2) / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                                   p->ef.conv_ny * p->ef.conv_nz * sizeof(soma_scalar_t));
        p->ef.pre_deriv_conv =
        (soma_scalar_t *) malloc(((p->ef.conv_nx + 2) / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                                   p->ef.conv_ny * p->ef.conv_nz * 6 * sizeof(soma_scalar_t));
        p->ef.Epot_conv =
        (soma_scalar_t *) malloc(((p->ef.conv_nx + 2) / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                                   p->ef.conv_ny * p->ef.conv_nz * sizeof(soma_scalar_t));
        p->ef.Epot_tmp_conv =
        (soma_scalar_t *) malloc(((p->ef.conv_nx + 2) / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                                   p->ef.conv_ny * p->ef.conv_nz * sizeof(soma_scalar_t));
        p->ef.n_cells_conv = ((p->ef.conv_nx + 2) / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                                   p->ef.conv_ny * p->ef.conv_nz;
    }

    if (p->ef.el_pos_xz)
    {
        p->ef.electrodes_conv =
        (soma_scalar_t *) malloc((p->ef.conv_nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                                 (p->ef.conv_ny + 2) * p->ef.conv_nz * sizeof(soma_scalar_t));
        p->ef.eps_arr_conv =
        (soma_scalar_t *) malloc((p->ef.conv_nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                                 (p->ef.conv_ny + 2) * p->ef.conv_nz * sizeof(soma_scalar_t));
        p->ef.pre_deriv_conv =
        (soma_scalar_t *) malloc((p->ef.conv_nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                                 (p->ef.conv_ny + 2) * p->ef.conv_nz * sizeof(soma_scalar_t));
        p->ef.Epot_conv =
        (soma_scalar_t *) malloc((p->ef.conv_nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                                 (p->ef.conv_ny + 2) * p->ef.conv_nz * sizeof(soma_scalar_t));
        p->ef.Epot_tmp_conv =
        (soma_scalar_t *) malloc((p->ef.conv_nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                                 (p->ef.conv_ny + 2) * p->ef.conv_nz * sizeof(soma_scalar_t));
        p->ef.n_cells_conv = (p->ef.conv_nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                             (p->ef.conv_ny + 2) * p->ef.conv_nz;
    }

    if (p->ef.el_pos_xy)
    {
        p->ef.electrodes_conv =
        (soma_scalar_t *) malloc((p->ef.conv_nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                                  p->ef.conv_ny * (p->ef.conv_nz + 2) * sizeof(soma_scalar_t));
        p->ef.eps_arr_conv =
        (soma_scalar_t *) malloc((p->ef.conv_nx  / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                                  p->ef.conv_ny * (p->ef.conv_nz + 2) * sizeof(soma_scalar_t));
        p->ef.pre_deriv_conv =
        (soma_scalar_t *) malloc((p->ef.conv_nx  / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                                  p->ef.conv_ny * (p->ef.conv_nz + 2) * 6 * sizeof(soma_scalar_t));
        p->ef.Epot_conv =
        (soma_scalar_t *) malloc((p->ef.conv_nx  / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                                  p->ef.conv_ny * (p->ef.conv_nz + 2) * sizeof(soma_scalar_t));
        p->ef.Epot_tmp_conv =
        (soma_scalar_t *) malloc((p->ef.conv_nx  / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                                  p->ef.conv_ny * (p->ef.conv_nz + 2) * sizeof(soma_scalar_t));
        p->ef.n_cells_conv = (p->ef.conv_nx  / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) *
                              p->ef.conv_ny * (p->ef.conv_nz + 2);
    }
    
    if (p->ef.electrodes_conv == NULL)
    {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
    }
    if (p->ef.eps_arr_conv == NULL)
    {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
    }
    if (p->ef.pre_deriv_conv == NULL)
    {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
    }
    if (p->ef.Epot_conv == NULL)
    {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
    }
    if (p->ef.Epot_tmp_conv == NULL)
    {
        fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
        return -1;
    }

    // fill Epot_conv & Epot_tmp_conv arrays with 0 - similar to their nonconvoluted counterparts
    for (uint64_t i=0; i < p->ef.n_cells_conv; i++)
    {
        p->ef.Epot_conv[i] = 0.0;
        p->ef.electrodes_conv[i] = 0;
    }

    // fill values for convoluted arrays electrode positions, determine offset to skip these planes during iteration
    // during convolution & iteration the remaining area in between electrodes will be filled;
    // planes consisting of electrodes will be skipped
    
    // soma_scalar_t eps_ave = 0;
    // for (uint8_t t=0; t < p->n_types; t++)
    // {
    //     eps_ave += p->ef.eps[t];
    // }
    // eps_ave /= p->n_types;

    if (p->ef.el_pos_yz)
    {
        p->ef.x_offset = 1;
        p->ef.y_offset = 0;
        p->ef.z_offset = 0;

        for (uint64_t y=0; y < p->ef.conv_ny; y++)
            for (uint64_t z=0; z < p->ef.conv_nz; z++)
            {
                

                p->ef.eps_arr_conv[0 * p->ef.conv_ny * p->ef.conv_nz + y * p->ef.conv_nz + z] = 1.0;                            // mean of all eps?!?!??!
                p->ef.eps_arr_conv[(p->ef.conv_nx + 1) * p->ef.conv_ny * p->ef.conv_nz + y * p->ef.conv_nz + z] = 1.0;

                for (uint8_t i=0; i<6; i++)
                {
                    p->ef.pre_deriv_conv[((0 * p->ef.conv_ny * p->ef.conv_nz + y * p->ef.conv_nz + z) * 6) + i] = 1.0/6.0;
                    p->ef.pre_deriv_conv[(((p->ef.conv_nx + 1) * p->ef.conv_ny * p->ef.conv_nz + y * p->ef.conv_nz + z) * 6) + i] = 1.0/6.0;
                }

                p->ef.Epot_conv[0 * p->ef.conv_ny * p->ef.conv_nz + y * p->ef.conv_nz + z] =
                    p->ef.Epot[0 * p->ny * p->nz + (y * p->ef.stride) * p->nz + (z * p->ef.stride)];
                p->ef.Epot_conv[(p->ef.conv_nx + 1) * p->ef.conv_ny * p->ef.conv_nz + y * p->ef.conv_nz + z] = 
                    p->ef.Epot[(p->nx - 1) * p->ny * p->nz + (y * p->ef.stride) * p->nz + (z * p->ef.stride)];

                p->ef.electrodes_conv[0 * p->ef.conv_ny * p->ef.conv_nz + y * p->ef.conv_nz + z] =
                    p->ef.electrodes[0 * p->ny * p->nz + (y * p->ef.stride) * p->nz + (z * p->ef.stride)];
                p->ef.electrodes_conv[(p->ef.conv_nx + 1) * p->ef.conv_ny * p->ef.conv_nz + y * p->ef.conv_nz + z] = 
                    p->ef.electrodes[(p->nx - 1) * p->ny * p->nz + (y * p->ef.stride) * p->nz + (z * p->ef.stride)];
            }
    }

    if (p->ef.el_pos_xz)
    {
        for (uint64_t x=0; x < p->ef.conv_nx; x++)
            for (uint64_t z=0; z < p->ef.conv_nz; z++)
            {
                p->ef.x_offset = 0;
                p->ef.y_offset = 1;
                p->ef.z_offset = 0;

                p->ef.eps_arr_conv[x * p->ef.conv_ny * p->ef.conv_nz + 0 * p->ef.conv_nz + z] = 1.0;
                p->ef.eps_arr_conv[x * p->ef.conv_ny * p->ef.conv_nz + (p->ef.conv_ny + 1) * p->ef.conv_nz + z] = 1.0;

                for (uint8_t i=0; i<6; i++)
                {
                    p->ef.pre_deriv_conv[((x * p->ef.conv_ny * p->ef.conv_nz + 0 * p->ef.conv_nz + z) * 6) + i] = 1.0/6.0;
                    p->ef.pre_deriv_conv[((x * p->ef.conv_ny * p->ef.conv_nz + (p->ef.conv_ny + 1) * p->ef.conv_nz + z) * 6) + i] = 1.0/6.0;
                }

                p->ef.Epot_conv[x * p->ef.conv_ny * p->ef.conv_nz + 0 * p->ef.conv_nz + z] = 
                    p->ef.Epot[(x * p->ef.stride) * p->ny * p->nz + 0 * p->nz + (z * p->ef.stride)];
                p->ef.Epot_conv[x * p->ef.conv_ny * p->ef.conv_nz + (p->ef.conv_ny + 1) * p->ef.conv_nz + z] = 
                    p->ef.Epot[(x * p->ef.stride) * p->ny * p->nz + (p->ny - 1) * p->nz + (z * p->ef.stride)];

                p->ef.electrodes_conv[x * p->ef.conv_ny * p->ef.conv_nz + 0 * p->ef.conv_nz + z] = 
                    p->ef.electrodes[(x * p->ef.stride) * p->ny * p->nz + 0 * p->nz + (z * p->ef.stride)];
                p->ef.electrodes_conv[x * p->ef.conv_ny * p->ef.conv_nz + (p->ef.conv_ny + 1) * p->ef.conv_nz + z] = 
                    p->ef.electrodes[(x * p->ef.stride) * p->ny * p->nz + (p->ny - 1) * p->nz + (z * p->ef.stride)];
            }
    }

    if (p->ef.el_pos_xy)
    {
        for (uint64_t x=0; x < p->ef.conv_nx; x++)
            for (uint64_t y=0; y < p->ef.conv_ny; y++)
            {
                p->ef.x_offset = 0;
                p->ef.y_offset = 0;
                p->ef.z_offset = 1;

                p->ef.eps_arr_conv[x * p->ef.conv_ny * p->ef.conv_nz + y * p->ef.conv_nz + 0] = 1.0;
                p->ef.eps_arr_conv[x * p->ef.conv_ny * p->ef.conv_nz + y * p->ef.conv_nz + (p->ef.conv_nz + 1)] = 1.0;

                for (uint8_t i=0; i<6; i++)
                {
                    p->ef.pre_deriv_conv[((x * p->ef.conv_ny * p->ef.conv_nz + y * p->ef.conv_nz + 0) * 6) + i] = 1.0/6.0;
                    p->ef.pre_deriv_conv[((x * p->ef.conv_ny * p->ef.conv_nz + y * p->ef.conv_nz + (p->ef.conv_nz + 1)) * 6) + i] = 1.0/6.0;
                }

                p->ef.Epot_conv[x * p->ef.conv_ny * p->ef.conv_nz + y * p->ef.conv_nz + 0] = 
                    p->ef.Epot[(x * p->ef.stride) * p->ny * p->nz + (y * p->ef.stride) * p->nz + 0];
                p->ef.Epot_conv[x * p->ef.conv_ny * p->ef.conv_nz + y * p->ef.conv_nz + (p->ef.conv_nz + 1)] = 
                    p->ef.Epot[(x * p->ef.stride) * p->ny * p->nz + (y * p->ef.stride) * p->nz + (p->nz - 1)];

                p->ef.electrodes_conv[x * p->ef.conv_ny * p->ef.conv_nz + y * p->ef.conv_nz + 0] = 
                    p->ef.electrodes[(x * p->ef.stride) * p->ny * p->nz + (y * p->ef.stride) * p->nz + 0];
                p->ef.electrodes_conv[x * p->ef.conv_ny * p->ef.conv_nz + y * p->ef.conv_nz + (p->ef.conv_nz + 1)] = 
                    p->ef.electrodes[(x * p->ef.stride) * p->ny * p->nz + (y * p->ef.stride) * p->nz + (p->nz - 1)];
            }
    }

    if (p->ef.x_offset == 0 && p->ef.y_offset == 0 && p->ef.z_offset == 0)
    {
        fprintf(stderr, "ERROR: Offset for electrodes could not be determined %s:%d\n", __FILE__, __LINE__);
        return -1;
    }

    // fill kernel normalization field
#pragma acc parallel loop present(p[0:1]) collapse(3)
#pragma omp parallel for collapse(3)
    for (uint64_t xc=0; xc < p->ef.conv_nx; xc++)
        for (uint64_t yc=0; yc < p->ef.conv_ny; yc++)
            for (uint64_t zc=0; zc < p->ef.conv_nz; zc++)
            {
                for (int8_t h = -p->ef.kernel_rad; h <= p->ef.kernel_rad; h++)
                    for (int8_t i = -p->ef.kernel_rad; i <= p->ef.kernel_rad; i++)
                        for (int8_t j = -p->ef.kernel_rad; j <= p->ef.kernel_rad; j++)
                        {
                            // Compute original indices; offset to skip electrode planes
                            int64_t x_i = xc * p->ef.stride + p->ef.x_offset + h;
                            int64_t y_i = yc * p->ef.stride + p->ef.y_offset + i;
                            int64_t z_i = zc * p->ef.stride + p->ef.z_offset + j;
                            // Resolve periodic/non-periodic boundaries, omit saving values in cells of electrodes planes, 
                            // these will be saved subsequently
                            if (p->ef.el_pos_yz)
                            {
                                if (x_i >= p->nx-1) continue;
                                if (x_i <= 0) continue;
                            }
                            else
                            {
                                if (x_i >= p->nx) x_i -= p->nx;
                                if (x_i < 0) x_i += p->nx;    
                            }
                            if (p->ef.el_pos_xz)
                            {
                                if (y_i >= p->ny-1) continue;
                                if (y_i <= 0) continue;
                            }
                            else
                            {
                                if (y_i >= p->ny) y_i -= p->ny;
                                if (y_i < 0) y_i += p->ny;    
                            }
                            if (p->ef.el_pos_xy)
                            {
                                if (z_i >= p->nz-1) continue;
                                if (z_i <= 0) continue;
                            }
                            else
                            {
                                if (z_i >= p->nz) z_i -= p->nz;
                                if (z_i < 0) z_i += p->nz;
                            }
                            
                            p->ef.kernel_norm_field[cell_to_index(p,x_i,y_i,z_i)] += p->ef.kernel[cell_to_index_kernel(p,h,i,j)];
                        }
            }

    // fill electrode planes (2D) of kernel field depending on electrode positions
    if (p->ef.el_pos_yz)
    {
#pragma acc parallel loop present(p[0:1]) collapse(2)
#pragma omp parallel for collapse(2)
    for (uint64_t yc=0; yc < p->ef.conv_ny; yc++)
        for (uint64_t zc=0; zc < p->ef.conv_nz; zc++)
        {
            for (int8_t i = -p->ef.kernel_rad; i <= p->ef.kernel_rad; i++)
                for (int8_t j = -p->ef.kernel_rad; j <= p->ef.kernel_rad; j++)
                {
                    // loop through kernel cell & resolve periodic boundaries
                    int64_t y_i = yc * p->ef.stride + i;
                    int64_t z_i = zc * p->ef.stride + j;
                    if (y_i >= p->ny) y_i -= p->ny;
                    if (y_i < 0) y_i += p->ny;
                    if (z_i >= p->nz) z_i -= p->nz;
                    if (z_i < 0) z_i += p->nz;
                    
                    // first and second electrode cell 
                    uint64_t c_el_1 = cell_to_index(p,0,y_i,z_i);
                    uint64_t c_el_2 = cell_to_index(p,p->nx-1,y_i,z_i);
                    
                    if (p->ef.electrodes[c_el_1] != 1) p->ef.kernel_norm_field[c_el_1] += p->ef.kernel[cell_to_index_kernel(p,0,i,j)];
                    if (p->ef.electrodes[c_el_2] != 1) p->ef.kernel_norm_field[c_el_2] += p->ef.kernel[cell_to_index_kernel(p,0,i,j)];
                }
        }
    }

    if (p->ef.el_pos_xz)
    {
#pragma acc parallel loop present(p[0:1]) collapse(2)
#pragma omp parallel for collapse(2)
    for (uint64_t xc=0; xc < p->ef.conv_nx; xc++)
        for (uint64_t zc=0; zc < p->ef.conv_nz; zc++)
        {
            for (int8_t h = -p->ef.kernel_rad; h <= p->ef.kernel_rad; h++)
                for (int8_t j = -p->ef.kernel_rad; j <= p->ef.kernel_rad; j++)
                {
                    // loop through kernel cell & resolve periodic boundaries
                    int64_t x_i = xc * p->ef.stride + h;
                    int64_t z_i = zc * p->ef.stride + j;
                    if (x_i >= p->nx) x_i -= p->nx;
                    if (x_i < 0) x_i += p->nx;
                    if (z_i >= p->nz) z_i -= p->nz;
                    if (z_i < 0) z_i += p->nz;
                    
                    // first and second electrode cell 
                    uint64_t c_el_1 = cell_to_index(p,x_i,0,z_i);
                    uint64_t c_el_2 = cell_to_index(p,x_i,p->ny-1,z_i);
                    
                    if (p->ef.electrodes[c_el_1] != 1) p->ef.kernel_norm_field[c_el_1] += p->ef.kernel[cell_to_index_kernel(p,h,0,j)];
                    if (p->ef.electrodes[c_el_2] != 1) p->ef.kernel_norm_field[c_el_2] += p->ef.kernel[cell_to_index_kernel(p,h,0,j)];
                }
        }
    }

    if (p->ef.el_pos_xy)
    {
#pragma acc parallel loop present(p[0:1]) collapse(2)
#pragma omp parallel for collapse(2)
    for (uint64_t xc=0; xc < p->ef.conv_nx; xc++)
        for (uint64_t yc=0; yc < p->ef.conv_ny; yc++)
        {
            for (int8_t h = -p->ef.kernel_rad; h <= p->ef.kernel_rad; h++)
                for (int8_t i = -p->ef.kernel_rad; i <= p->ef.kernel_rad; i++)
                {
                    // loop through kernel cell & resolve periodic boundaries
                    int64_t x_i = xc * p->ef.stride + h;
                    int64_t y_i = yc * p->ef.stride + i;
                    if (x_i >= p->nx) x_i -= p->nx;
                    if (x_i < 0) x_i += p->nx;
                    if (y_i >= p->ny) y_i -= p->ny;
                    if (y_i < 0) y_i += p->ny;
                    
                    // first and second electrode cell 
                    uint64_t c_el_1 = cell_to_index(p,x_i,y_i,0);
                    uint64_t c_el_2 = cell_to_index(p,x_i,y_i,p->nz-1);
                    
                    if (p->ef.electrodes[c_el_1] != 1) p->ef.kernel_norm_field[c_el_1] += p->ef.kernel[cell_to_index_kernel(p,h,i,0)];
                    if (p->ef.electrodes[c_el_2] != 1) p->ef.kernel_norm_field[c_el_2] += p->ef.kernel[cell_to_index_kernel(p,h,i,0)];
                }
        }
    }

    return 0;
}


int copyin_electric_field(struct Phase *p)
{
    if (p->ef.iter_limit != 0)
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
            if (p->ef.kernel_dim > 1)
            {
#pragma acc enter data copyin(p->ef.kernel[0:p->ef.kernel_dim*p->ef.kernel_dim*p->ef.kernel_dim])
#pragma acc enter data copyin(p->ef.kernel_norm_field[0:p->n_cells_local])
#pragma acc enter data copyin(p->ef.kernel_blur[0:p->ef.kernel_dim*p->ef.kernel_dim*p->ef.kernel_dim])
#pragma acc enter data copyin(p->ef.electrodes_conv[0:p->ef.n_cells_conv])
#pragma acc enter data copyin(p->ef.eps_arr_conv[0:p->ef.n_cells_conv])
#pragma acc enter data copyin(p->ef.pre_deriv_conv[0:p->ef.n_cells_conv*6])
#pragma acc enter data copyin(p->ef.Epot_conv[0:p->ef.n_cells_conv])
#pragma acc enter data copyin(p->ef.Epot_tmp_conv[0:p->ef.n_cells_conv])
            }
#endif                          //_OPENACC
        }
    return 0;
}

int copyout_electric_field(struct Phase *p)
{
    if (p->ef.iter_limit != 0)
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
            if (p->ef.kernel_dim > 1)
            {
#pragma acc exit data copyout(p->ef.kernel[0:p->ef.kernel_dim*p->ef.kernel_dim*p->ef.kernel_dim])
#pragma acc exit data copyout(p->ef.kernel_norm_field[0:p->n_cells_local])
#pragma acc exit data copyout(p->ef.kernel_blur[0:p->ef.kernel_dim*p->ef.kernel_dim*p->ef.kernel_dim])
#pragma acc exit data copyout(p->ef.eps_arr_conv[0:p->ef.n_cells_conv])
#pragma acc exit data copyout(p->ef.electrodes_conv[0:p->ef.n_cells_conv])
#pragma acc exit data copyout(p->ef.pre_deriv_conv[0:p->ef.n_cells_conv*6])
#pragma acc exit data copyout(p->ef.Epot_conv[0:p->ef.n_cells_conv])
#pragma acc exit data copyout(p->ef.Epot_tmp_conv[0:p->ef.n_cells_conv])
            }
#endif                          //_OPENACC
        }
    return 0;
}

int update_self_electric_field(const struct Phase *const p)
{
    if (p->ef.iter_limit != 0)
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
            if (p->ef.kernel_dim > 1)
            {
#pragma acc update self(p->ef.kernel[0:p->ef.kernel_dim*p->ef.kernel_dim*p->ef.kernel_dim])
#pragma acc update self(p->ef.kernel_norm_field[0:p->n_cells_local])
#pragma acc update self(p->ef.kernel_blur[0:p->ef.kernel_dim*p->ef.kernel_dim*p->ef.kernel_dim])
#pragma acc update self(p->ef.eps_arr_conv[0:p->ef.n_cells_conv])
#pragma acc update self(p->ef.electrodes_conv[0:p->ef.n_cells_conv])
#pragma acc update self(p->ef.pre_deriv_conv[0:p->ef.n_cells_conv*6])
#pragma acc update self(p->ef.Epot_conv[0:p->ef.n_cells_conv])
#pragma acc update self(p->ef.Epot_tmp_conv[0:p->ef.n_cells_conv])
            }
#endif                          //_OPENACC
        }
    return 0;
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
#pragma acc parallel loop present(p[0:1]) collapse(3) //deviceptr(p->ef.eps_arr)
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
    p->ef.amt_iter = 0;

    while(max_i > p->ef.thresh_iter && p->ef.amt_iter < p->ef.iter_limit)
    {
        max_i = p->ef.thresh_iter;

#pragma acc parallel loop present(p[0:1]) collapse(3) //reduction(max: max_i) deviceptr(p->ef.electrodes,p->ef.Epot,p->ef.Epot_tmp,p->ef.pre_deriv)
#pragma omp parallel for collapse(3) //reduction(max: max_i)
        for (uint64_t x=0; x < p->nx; x++)
            for (uint64_t y=0; y < p->ny; y++)
    	        for (uint64_t z=0; z < p->nz; z++)
                    {
                        uint64_t i = cell_to_index(p,x,y,z);
                        //Keep electric potential field constant at electrode position
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

                            //Option without using pre-calculated derivatives (p->ef.pre_deriv)
                            // p->ef.Epot_tmp[i] = (1.0 / 6.0) * ( 
                            //                      p->ef.Epot[cell_to_index(p,x+1,y,z)] + p->ef.Epot[cell_to_index(p,x-1,y,z)] + 
                            //                      p->ef.Epot[cell_to_index(p,x,y+1,z)] + p->ef.Epot[cell_to_index(p,x,y-1,z)] + 
                            //                      p->ef.Epot[cell_to_index(p,x,y,z+1)] + p->ef.Epot[cell_to_index(p,x,y,z-1)]) +
                            //                     (1.0 / (24.0 * p->ef.eps_arr[i])) * (
                            //                     (p->ef.eps_arr[cell_to_index(p,x+1,y,z)] - p->ef.eps_arr[cell_to_index(p,x-1,y,z)]) *
                            //                     (p->ef.Epot[cell_to_index(p,x+1,y,z)] - p->ef.Epot[cell_to_index(p,x-1,y,z)]) +
                            //                     (p->ef.eps_arr[cell_to_index(p,x,y+1,z)] - p->ef.eps_arr[cell_to_index(p,x,y-1,z)]) *
                            //                     (p->ef.Epot[cell_to_index(p,x,y+1,z)] - p->ef.Epot[cell_to_index(p,x,y-1,z)]) +
                            //                     (p->ef.eps_arr[cell_to_index(p,x,y,z+1)] - p->ef.eps_arr[cell_to_index(p,x,y,z-1)]) *
                            //                     (p->ef.Epot[cell_to_index(p,x,y,z+1)] - p->ef.Epot[cell_to_index(p,x,y,z-1)]) );
                        }
                    }
	
#pragma acc parallel loop present(p[0:1]) collapse(3) reduction(max: max_i) //deviceptr(p->ef.Epot_tmp,p->ef.eps_arr)
#pragma omp parallel for collapse(3) reduction(max: max_i)
        for (uint64_t x=0; x < p->nx; x++)
            for (uint64_t y=0; y < p->ny; y++)
                for (uint64_t z=0; z < p->nz; z++)
     		    {
                    uint64_t j = cell_to_index(p,x,y,z);

         			if(p->ef.electrodes[j] != 1)
         			{
         			    // Laplace equation, welling2017 eq. 5
        			    soma_scalar_t conv_crit = 
                            (d2Epotx(p,p->ef.Epot_tmp,x,y,z,p->nx,p->ny,p->nz) +
         					 d2Epoty(p,p->ef.Epot_tmp,x,y,z,p->ny,p->nz) +
         					 d2Epotz(p,p->ef.Epot_tmp,x,y,z,p->ny,p->nz)) * p->ef.eps_arr[j] +
         			    	(dEpotx(p,p->ef.Epot_tmp,x,y,z,p->nx,p->ny,p->nz) * dEpotx(p,p->ef.eps_arr,x,y,z,p->nx,p->ny,p->nz) +
         					 dEpoty(p,p->ef.Epot_tmp,x,y,z,p->ny,p->nz) * dEpoty(p,p->ef.eps_arr,x,y,z,p->ny,p->nz) +
         					 dEpotz(p,p->ef.Epot_tmp,x,y,z,p->ny,p->nz) * dEpotz(p,p->ef.eps_arr,x,y,z,p->ny,p->nz));

        			    if (max_i < conv_crit) max_i = conv_crit;
         			}
                } 
	
#pragma acc parallel loop present(p[0:1]) collapse(3) //deviceptr(p->ef.Epot,p->ef.Epot_tmp)
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

        // if( (p->ef.amt_iter+1) % 500==0 ) printf("\n max: %.3e\n", max_i);

// #pragma acc update host(p->ef)
//         soma_scalar_t *tmp = p->ef.Epot_tmp;
//         p->ef.Epot_tmp = p->ef.Epot;
//         p->ef.Epot = tmp;
// #pragma acc update device(p->ef)
        p->ef.amt_iter += 1;
    }
    
    // if (p->time == 0 || (p->time + 1) % 100 == 0) printf("MC step: %d, iterations to solve ef: %ld \n",p->time+1, p->ef.amt_iter);

    return max_i;
}

void convolution_eps_arr(struct Phase *const p)
{
    uint64_t conv_err = 0;
    uint64_t err_x, err_y, err_z;

    // convoluted indices "xc,yc,zc", "p->ef.conv_nx, etc." only describe the area between electrodes
#pragma acc parallel loop present(p[0:1]) collapse(3)
#pragma omp parallel for collapse(3)
    for (uint64_t xc=0; xc < p->ef.conv_nx; xc++)
        for (uint64_t yc=0; yc < p->ef.conv_ny; yc++)
            for (uint64_t zc=0; zc < p->ef.conv_nz; zc++)
            {
                soma_scalar_t kernel_sum = 0.0;

                for (int8_t h = -p->ef.kernel_rad; h <= p->ef.kernel_rad; h++)
                    for (int8_t i = -p->ef.kernel_rad; i <= p->ef.kernel_rad; i++)
                        for (int8_t j = -p->ef.kernel_rad; j <= p->ef.kernel_rad; j++)
                        {
                            // compute original indices; offset to skip electrode planes
                            int64_t x_i = xc * p->ef.stride + p->ef.x_offset + h;
                            int64_t y_i = yc * p->ef.stride + p->ef.y_offset + i;
                            int64_t z_i = zc * p->ef.stride + p->ef.z_offset + j;
                            // Resolve periodic boundaries, if kernel_dim > 3 make sure not to wrap back on opposing electrode
                            if (p->ef.el_pos_yz)
                            {
                                if (x_i >= p->nx) x_i = p->nx - 1;
                                if (x_i < 0) x_i = 0;    
                            }
                            else
                            {
                                if (x_i >= p->nx) x_i -= p->nx;
                                if (x_i < 0) x_i += p->nx;    
                            }
                            if (p->ef.el_pos_xz)
                            {
                                if (y_i >= p->ny) y_i = p->ny - 1;
                                if (y_i < 0) y_i = 0;
                            }
                            else
                            {
                                if (y_i >= p->ny) x_i -= p->ny;
                                if (y_i < 0) y_i += p->ny;    
                            }
                            if (p->ef.el_pos_xy)
                            {
                                if (z_i >= p->nz) z_i = p->nz - 1;
                                if (z_i < 0) z_i = 0;
                            }
                            else
                            {
                                if (z_i >= p->nz) z_i -= p->nz;
                                if (z_i < 0) z_i += p->nz;
                            }
                            
                            kernel_sum += p->ef.kernel[cell_to_index_kernel(p,h,i,j)] *
                                          p->ef.eps_arr[cell_to_index(p,x_i,y_i,z_i)];
                        }
                
                if (kernel_sum != kernel_sum) 
                {
                    if (conv_err == 0)
                    {
                        err_x = xc;
                        err_y = yc;
                        err_z = zc;
                    }
                    conv_err += 1;
                }
                
                // skip electrode planes by adding offsets to convoluted indices
                uint64_t xc_i = xc + p->ef.x_offset;
                uint64_t yc_i = yc + p->ef.y_offset;
                uint64_t zc_i = zc + p->ef.z_offset;
                p->ef.eps_arr_conv[cell_to_index_conv(p,xc_i,yc_i,zc_i)] = kernel_sum;
            }

    if (conv_err != 0)
    {
        fprintf(stderr, "Convolution errors of eps_arr: %ld. %s:%d. xc=%ld, yc=%ld, zc=%ld.\n", conv_err, __FILE__, __LINE__, err_x,err_y,err_z);
    }
    
    uint64_t eps_conv_err = 0;
    for (uint64_t i=0; i<(p->ef.conv_nx+2)*p->ef.conv_ny*p->ef.conv_nz;i++)
    {
        if(p->ef.eps_arr_conv[i] != p->ef.eps_arr_conv[i]) eps_conv_err += 1;
    }
    if (eps_conv_err != 0)
    {
        fprintf(stderr, "ERROR: Calculating convoluted eps_arr: %ld. %s:%d.", eps_conv_err, __FILE__, __LINE__);
    } 
}

void pre_derivatives_conv(struct Phase *const p)
{
    uint64_t conv_nx = p->ef.conv_nx;
    uint64_t conv_ny = p->ef.conv_ny;
    uint64_t conv_nz = p->ef.conv_nz;
    if (p->ef.el_pos_yz) conv_nx += 2;
    if (p->ef.el_pos_xz) conv_ny += 2;
    if (p->ef.el_pos_xy) conv_nz += 2;

#pragma acc parallel loop present(p[0:1]) collapse(3)
#pragma omp parallel for collapse(3)
    for (uint64_t xc=0; xc < conv_nx; xc++)
        for (uint64_t yc=0; yc < conv_ny; yc++)
            for (uint64_t zc=0; zc < conv_nz; zc++)
            {
                uint64_t i = cell_to_index_conv(p,xc,yc,zc);
                uint64_t pos6 = (xc * (p->ef.conv_ny + p->ef.y_offset * 2) * (p->ef.conv_nz + p->ef.z_offset * 2) +
                                 yc * (p->ef.conv_nz + p->ef.z_offset * 2) + zc) * 6;
                
                if (p->ef.electrodes_conv[i] == 1)
                {
                    for (uint8_t j=0; j < 6; j++)
                        p->ef.pre_deriv_conv[pos6+j] = 1.0/6.0;
                }
                else
                {                    
                    soma_scalar_t cr = 1.0/(24.0 * p->ef.eps_arr_conv[i]);

                    soma_scalar_t epsx = cr * ( p->ef.eps_arr_conv[cell_to_index_conv(p,xc+1,yc,zc)] -
                                                p->ef.eps_arr_conv[cell_to_index_conv(p,xc-1,yc,zc)] );
                    soma_scalar_t epsy = cr * ( p->ef.eps_arr_conv[cell_to_index_conv(p,xc,yc+1,zc)] -
                                                p->ef.eps_arr_conv[cell_to_index_conv(p,xc,yc-1,zc)] );
                    soma_scalar_t epsz = cr * ( p->ef.eps_arr_conv[cell_to_index_conv(p,xc,yc,zc+1)] -
                                                p->ef.eps_arr_conv[cell_to_index_conv(p,xc,yc,zc-1)] );
                    
                    p->ef.pre_deriv_conv[pos6 + 0] = epsx + 1.0/6.0;
                    p->ef.pre_deriv_conv[pos6 + 1] =-epsx + 1.0/6.0;
                    p->ef.pre_deriv_conv[pos6 + 2] = epsy + 1.0/6.0;
                    p->ef.pre_deriv_conv[pos6 + 3] =-epsy + 1.0/6.0;
                    p->ef.pre_deriv_conv[pos6 + 4] = epsz + 1.0/6.0;
                    p->ef.pre_deriv_conv[pos6 + 5] =-epsz + 1.0/6.0;
                }
            }
    
    uint64_t pre_deriv_err = 0;
    for (uint64_t i=0; i<p->ef.n_cells_conv;i++)
    {
        if(p->ef.pre_deriv_conv[i] != p->ef.pre_deriv_conv[i]) pre_deriv_err += 1;
    }
    if (pre_deriv_err != 0)
    {
        fprintf(stderr, "ERROR: Calculating convoluted pre_deriv: %ld. %s:%d.", pre_deriv_err, __FILE__, __LINE__);
    }
}

soma_scalar_t iterate_field_conv(struct Phase *const p)
{
    soma_scalar_t max_i = 10.0;
    p->ef.amt_iter = 0;
    uint64_t conv_nx = p->ef.conv_nx;
    uint64_t conv_ny = p->ef.conv_ny;
    uint64_t conv_nz = p->ef.conv_nz;
    if (p->ef.el_pos_yz) conv_nx += 2;
    if (p->ef.el_pos_xz) conv_ny += 2;
    if (p->ef.el_pos_xy) conv_nz += 2;

    while(max_i > p->ef.thresh_iter && p->ef.amt_iter < p->ef.iter_limit)
    {
        max_i = p->ef.thresh_iter;
        
        // compute new solution, save in Epot_tmp_conv
#pragma acc parallel loop present(p[0:1]) collapse(3)
#pragma omp parallel for collapse(3)
        for (uint64_t xc=0; xc < conv_nx; xc++)
            for (uint64_t yc=0; yc < conv_ny; yc++)
                for (uint64_t zc=0; zc < conv_nz; zc++)
                    {
                        uint64_t i = cell_to_index_conv(p,xc,yc,zc);
                        if (p->ef.electrodes_conv[i] == 1)
                        {
                            p->ef.Epot_tmp_conv[i] = p->ef.Epot_conv[i];
                        }
                        else
                        {

                        uint64_t pos6 = (xc * (p->ef.conv_ny + p->ef.y_offset * 2) * (p->ef.conv_nz + p->ef.z_offset * 2) +
                                         yc * (p->ef.conv_nz + p->ef.z_offset * 2) + zc) * 6;
                        
                        //cell_to_index_conv() resolves non-periodic boundaries
                        p->ef.Epot_tmp_conv[i] =
                            p->ef.pre_deriv_conv[pos6+0] * p->ef.Epot_conv[cell_to_index_conv(p,xc+1,yc,zc)] +
                            p->ef.pre_deriv_conv[pos6+1] * p->ef.Epot_conv[cell_to_index_conv(p,xc-1,yc,zc)] +
                            p->ef.pre_deriv_conv[pos6+2] * p->ef.Epot_conv[cell_to_index_conv(p,xc,yc+1,zc)] +
                            p->ef.pre_deriv_conv[pos6+3] * p->ef.Epot_conv[cell_to_index_conv(p,xc,yc-1,zc)] +
                            p->ef.pre_deriv_conv[pos6+4] * p->ef.Epot_conv[cell_to_index_conv(p,xc,yc,zc+1)] +
                            p->ef.pre_deriv_conv[pos6+5] * p->ef.Epot_conv[cell_to_index_conv(p,xc,yc,zc-1)];
                        }
                    }

        // compute conversion criterion
#pragma acc parallel loop present(p[0:1]) collapse(3) reduction(max: max_i)
#pragma omp parallel for collapse(3) reduction(max: max_i)
        for (uint64_t xc=0; xc < p->ef.conv_nx; xc++)
            for (uint64_t yc=0; yc < p->ef.conv_ny; yc++)
                for (uint64_t zc=0; zc < p->ef.conv_nz; zc++)
                {
                    uint64_t xc_i = xc + p->ef.x_offset;
                    uint64_t yc_i = yc + p->ef.y_offset;
                    uint64_t zc_i = zc + p->ef.z_offset;
                    uint64_t j = cell_to_index_conv(p,xc_i,yc_i,zc_i);

                    uint64_t nx = p->ef.conv_nx + (p->ef.x_offset * 2);
                    uint64_t ny = p->ef.conv_ny + (p->ef.y_offset * 2);
                    uint64_t nz = p->ef.conv_nz + (p->ef.z_offset * 2);
                    // Laplace equation, welling2017 eq. 5

                    soma_scalar_t conv_crit = 
                        (d2Epotx(p,p->ef.Epot_tmp_conv,xc_i,yc_i,zc_i,nx,ny,nz) +
                         d2Epoty(p,p->ef.Epot_tmp_conv,xc_i,yc_i,zc_i,ny,nz) +
                         d2Epotz(p,p->ef.Epot_tmp_conv,xc_i,yc_i,zc_i,ny,nz) ) * p->ef.eps_arr_conv[j] +
                        (dEpotx(p,p->ef.Epot_tmp_conv,xc_i,yc_i,zc_i,nx,ny,nz) * 
                         dEpotx(p,p->ef.eps_arr_conv,xc_i,yc_i,zc_i,nx,ny,nz) +
                         dEpoty(p,p->ef.Epot_tmp_conv,xc_i,yc_i,zc_i,ny,nz) *
                         dEpoty(p,p->ef.eps_arr_conv,xc_i,yc_i,zc_i,ny,nz) +
                         dEpotz(p,p->ef.Epot_tmp_conv,xc_i,yc_i,zc_i,ny,nz) *
                         dEpotz(p,p->ef.eps_arr_conv,xc_i,yc_i,zc_i,ny,nz) );

                    if (max_i < conv_crit) max_i = conv_crit;
                }

        // save new solution to Epot_conv
#pragma acc parallel loop present(p[0:1])
#pragma omp parallel for
        for (uint64_t c=0; c < p->ef.n_cells_conv; c++)
        { 
            p->ef.Epot_conv[c] = p->ef.Epot_tmp_conv[c];
        }

// #pragma acc update host(p->ef)
//         soma_scalar_t *tmp_conv = p->ef.Epot_tmp_conv;
//         p->ef.Epot_tmp_conv = p->ef.Epot_conv;
//         p->ef.Epot_conv = tmp_conv;
// #pragma acc update device(p->ef)

        // if( (p->ef.amt_iter+1) % 500==0 ) printf("max: %.3e\n", max_i);
        p->ef.amt_iter += 1;
    }
    // if (p->time == 0 || (p->time + 1) % 100 == 0) printf("MC step: %d, iterations to solve ef: %ld \n",p->time+1, p->ef.amt_iter);
        
    return max_i;
}

void deconvolution_Epot(struct Phase *const p)
{
    if (p->ef.stride == 1)
    {
#pragma acc parallel loop present(p[0:1]) collapse(3)
#pragma omp parallel for
    for (uint64_t xc=0; xc < p->ef.conv_nx; xc++)
        for (uint64_t yc=0; yc < p->ef.conv_ny; yc++)
            for (uint64_t zc=0; zc < p->ef.conv_nz; zc++)
            {
                int64_t x_i = xc + p->ef.x_offset;
                int64_t y_i = yc + p->ef.y_offset;
                int64_t z_i = zc + p->ef.z_offset;
                uint64_t cell_i = cell_to_index(p,x_i,y_i,z_i);
                p->ef.Epot[cell_i] = p->ef.Epot_conv[cell_i]; 
            }       
    }
    // devonvolution for stride > 1 where empty cells have to be filled
    else
    {
        // empty Epot except electrode planes to fill with updated values 
#pragma acc parallel loop present(p[0:1])
#pragma omp parallel for
        for (uint64_t c=0; c < p->n_cells; c++)
        {
            if (p->ef.electrodes[c] != 1) p->ef.Epot[c] = 0.0;
        }

        // fill empty cells by weighing convoluted values with gauss kernel; omit electrode planes - will be resolved subsequently in 2D
#pragma acc parallel loop present(p[0:1]) collapse(3)
#pragma omp parallel for
        for (uint64_t xc=0; xc < p->ef.conv_nx; xc++)
            for (uint64_t yc=0; yc < p->ef.conv_ny; yc++)
                for (uint64_t zc=0; zc < p->ef.conv_nz; zc++)
                {
                    // kernel loop
                    for (int8_t h=-p->ef.kernel_rad; h <= p->ef.kernel_rad; h++)
                        for (int8_t i=-p->ef.kernel_rad; i <= p->ef.kernel_rad; i++)
                            for (int8_t j=-p->ef.kernel_rad; j <= p->ef.kernel_rad; j++)
                            {
                                // Compute original indices; offset to skip electrode planes
                                int64_t x_i = xc * p->ef.stride + p->ef.x_offset + h;
                                int64_t y_i = yc * p->ef.stride + p->ef.y_offset + i;
                                int64_t z_i = zc * p->ef.stride + p->ef.z_offset + j;
                                // Resolve periodic/non-periodic boundaries, omit indeces within electrode planes
                                if (p->ef.el_pos_yz)
                                {
                                    if (x_i >= p->nx-1) continue;
                                    if (x_i <= 0) continue;
                                }
                                else
                                {
                                    if (x_i >= p->nx) x_i -= p->nx;
                                    if (x_i < 0) x_i += p->nx;
                                }
                                if (p->ef.el_pos_xz)
                                {
                                    if (y_i >= p->ny-1) continue;
                                    if (y_i <= 0) continue;
                                }
                                else
                                {
                                    if (y_i >= p->ny) y_i -= p->ny;
                                    if (y_i < 0) y_i += p->ny;
                                }
                                if (p->ef.el_pos_xy)
                                {
                                    if (z_i >= p->nz-1) continue;
                                    if (z_i <= 0) continue;
                                }
                                else
                                {
                                    if (z_i >= p->nz) z_i -= p->nz;
                                    if (z_i < 0) z_i += p->nz;
                                }
                                
                                uint64_t cell_i = cell_to_index(p,x_i,y_i,z_i);
                                // skip electrode planes (expected to be in area51)
                                // if (p->area51[cell_i] != 1)
                                // {
                                uint64_t xc_i = xc + p->ef.x_offset;
                                uint64_t yc_i = yc + p->ef.y_offset;
                                uint64_t zc_i = zc + p->ef.z_offset;
                                p->ef.Epot[cell_i] += p->ef.Epot_conv[cell_to_index_conv(p,xc_i,yc_i,zc_i)] *
                                                          p->ef.kernel[cell_to_index_kernel(p,h,i,j)];
                                // }
                                
                            }
                }

        // fill empty cells of electrode planes by weighing convoluted values with gauss kernel (only in 2D within electrode plane)
        if (p->ef.el_pos_yz)
        {
#pragma acc parallel loop present(p[0:1]) collapse(2)
#pragma omp parallel for collapse(2)
            for (uint64_t yc=0; yc < p->ef.conv_ny; yc++)
                for (uint64_t zc=0; zc < p->ef.conv_nz; zc++)
                {
                    for (int8_t i = -p->ef.kernel_rad; i <= p->ef.kernel_rad; i++)
                        for (int8_t j = -p->ef.kernel_rad; j <= p->ef.kernel_rad; j++)
                        {
                            // loop through kernel cell & resolve periodic boundaries
                            int64_t y_i = yc * p->ef.stride + i;
                            int64_t z_i = zc * p->ef.stride + j;
                            if (y_i >= p->ny) y_i -= p->ny;
                            if (y_i < 0) y_i += p->ny;
                            if (z_i >= p->nz) z_i -= p->nz;
                            if (z_i < 0) z_i += p->nz;
                            
                            // first and second electrode cell 
                            uint64_t c_el_1 = cell_to_index(p,0,y_i,z_i);
                            uint64_t c_el_2 = cell_to_index(p,p->nx-1,y_i,z_i);
                            
                            if (p->ef.electrodes[c_el_1] != 1) 
                            {
                                p->ef.Epot[c_el_1] += p->ef.Epot_conv[cell_to_index_conv(p,0,yc,zc)] *
                                                                   p->ef.kernel[cell_to_index_kernel(p,0,i,j)];
                            }
                            if (p->ef.electrodes[c_el_2] != 1)
                            {
                                p->ef.Epot[c_el_2] += p->ef.Epot_conv[cell_to_index_conv(p,p->ef.conv_nx+1,yc,zc)] *
                                                                   p->ef.kernel[cell_to_index_kernel(p,0,i,j)];
                            }
                        }
                }
        }

        if (p->ef.el_pos_xz)
        {
#pragma acc parallel loop present(p[0:1]) collapse(2)
#pragma omp parallel for collapse(2)
            for (uint64_t xc=0; xc < p->ef.conv_nx; xc++)
                for (uint64_t zc=0; zc < p->ef.conv_nz; zc++)
                {
                    for (int8_t h = -p->ef.kernel_rad; h <= p->ef.kernel_rad; h++)
                        for (int8_t j = -p->ef.kernel_rad; j <= p->ef.kernel_rad; j++)
                        {
                            // loop through kernel cell & resolve periodic boundaries
                            int64_t x_i = xc * p->ef.stride + h;
                            int64_t z_i = zc * p->ef.stride + j;
                            if (x_i >= p->nx) x_i -= p->nx;
                            if (x_i < 0) x_i += p->nx;
                            if (z_i >= p->nz) z_i -= p->nz;
                            if (z_i < 0) z_i += p->nz;
                            
                            // first and second electrode cell 
                            uint64_t c_el_1 = cell_to_index(p,x_i,0,z_i);
                            uint64_t c_el_2 = cell_to_index(p,x_i,p->ny-1,z_i);
                            
                            if (p->ef.electrodes[c_el_1] != 1) 
                            {
                                p->ef.Epot[c_el_1] += p->ef.Epot_conv[cell_to_index_conv(p,xc,0,zc)] *
                                                                   p->ef.kernel[cell_to_index_kernel(p,h,0,j)];
                            }
                            if (p->ef.electrodes[c_el_2] != 1)
                            {
                                p->ef.Epot[c_el_2] += p->ef.Epot_conv[cell_to_index_conv(p,xc,p->ef.conv_ny+1,zc)] *
                                                                   p->ef.kernel[cell_to_index_kernel(p,h,0,j)];
                            }
                        }
                }
        }

        if (p->ef.el_pos_xy)
        {
#pragma acc parallel loop present(p[0:1]) collapse(2)
#pragma omp parallel for collapse(2)
            for (uint64_t xc=0; xc < p->ef.conv_nx; xc++)
                for (uint64_t yc=0; yc < p->ef.conv_ny; yc++)
                {
                    for (int8_t h = -p->ef.kernel_rad; h <= p->ef.kernel_rad; h++)
                        for (int8_t i = -p->ef.kernel_rad; i <= p->ef.kernel_rad; i++)
                        {
                            // loop through kernel cell & resolve periodic boundaries
                            int64_t x_i = xc * p->ef.stride + h;
                            int64_t y_i = yc * p->ef.stride + i;
                            if (x_i >= p->nx) x_i -= p->nx;
                            if (x_i < 0) x_i += p->nx;
                            if (y_i >= p->ny) y_i -= p->ny;
                            if (y_i < 0) y_i += p->ny;
                            
                            // first and second electrode cell 
                            uint64_t c_el_1 = cell_to_index(p,x_i,y_i,0);
                            uint64_t c_el_2 = cell_to_index(p,x_i,y_i,p->nz-1);
                            
                            if (p->ef.electrodes[c_el_1] != 1) 
                            {
                                p->ef.Epot[c_el_1] += p->ef.Epot_conv[cell_to_index_conv(p,xc,yc,0)] *
                                                                   p->ef.kernel[cell_to_index_kernel(p,h,i,0)];
                            }
                            if (p->ef.electrodes[c_el_2] != 1)
                            {
                                p->ef.Epot[c_el_2] += p->ef.Epot_conv[cell_to_index_conv(p,xc,yc,p->ef.conv_nz+1)] *
                                                                   p->ef.kernel[cell_to_index_kernel(p,h,i,0)];
                            }
                        }
                }
        }

        // normalize values 
#pragma acc parallel loop present(p[0:1])
#pragma omp parallel for
    for (uint64_t c=0; c < p->n_cells; c++)
        {
            if (p->ef.electrodes[c] != 1) p->ef.Epot[c] /= p->ef.kernel_norm_field[c];
        }
    } //else

}

void smooth_Epot(struct Phase *const p)
{
#pragma acc parallel loop present(p[0:1])
#pragma omp parallel for
    for (uint64_t c=0; c < p->n_cells; c++)
    {
        p->ef.Epot_tmp[c] = p->ef.Epot[c];
    }

    // resolve electrode position; skip these planes
    uint64_t nx_cells = p->nx;
    uint64_t ny_cells = p->ny;
    uint64_t nz_cells = p->nz;
    if (p->ef.el_pos_yz) nx_cells -= 2;
    if (p->ef.el_pos_xz) ny_cells -= 2;
    if (p->ef.el_pos_xy) nz_cells -= 2;

#pragma acc parallel loop present(p[0:1]) collapse(3)
#pragma omp parallel for
    for (uint64_t x=p->ef.x_offset; x < nx_cells; x++)
        for (uint64_t y=p->ef.y_offset; y < ny_cells; y++)
            for (uint64_t z=p->ef.z_offset; z < nz_cells; z++)
            {
                soma_scalar_t kernel_blur_sum = 0.0;

                for (int8_t h=-1; h <= 1; h++)
                    for (int8_t i=0-1; i <= 1; i++)
                        for (int8_t j=-1; j <= 1; j++)
                        {
                            int64_t x_i = x + h;
                            int64_t y_i = y + i;
                            int64_t z_i = z + j;
                            // resolve periodic boundaries
                            if(x_i >= p->nx) x_i -= p->nx;
                            if(x_i < p->nx) x_i += p->nx;
                            if(y_i >= p->ny) y_i -= p->ny;
                            if(y_i < p->ny) y_i += p->ny;
                            if(z_i >= p->nz) z_i -= p->nz;
                            if(z_i < p->nz) z_i += p->nz;

                            kernel_blur_sum += p->ef.kernel_blur[cell_to_index_kernel(p,h,i,j)] *
                                               p->ef.Epot_tmp[cell_to_index(p,x_i,y_i,z_i)];
                        }

                p->ef.Epot[cell_to_index(p,x,y,z)] = kernel_blur_sum;
            }
}

int calc_electric_field_contr(struct Phase *const p)
{
    soma_scalar_t max = p->ef.thresh_iter;
    
//     if (p->time==0)
//     {
// #pragma acc update device(p->sqrt_Nbar)
//     // calc_sqrt_Nbar(p);
//     }

    calc_dielectric_field(p);

    if (p->ef.kernel_dim > 1)
    {
        convolution_eps_arr(p);
        pre_derivatives_conv(p);
        max = iterate_field_conv(p);
        deconvolution_Epot(p);
        // if (p->ef.stride > 1) smooth_Epot(p);
    }
    else
    {
        pre_derivatives(p);
        max = iterate_field(p);
    }

    if (max > p->ef.thresh_iter) return -1; 

    // soma_scalar_t sum_H_el = 0.0;
#pragma acc parallel loop present(p[0:1]) collapse(3) //reduction(+:sum_H_el)
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
                    soma_scalar_t dEpot_sq = 
                        dEpotx(p,p->ef.Epot,x,y,z,p->nx,p->ny,p->nz) * dEpotx(p,p->ef.Epot,x,y,z,p->nx,p->ny,p->nz) +
                        dEpoty(p,p->ef.Epot,x,y,z,p->ny,p->nz) * dEpoty(p,p->ef.Epot,x,y,z,p->ny,p->nz) +
                        dEpotz(p,p->ef.Epot,x,y,z,p->ny,p->nz) * dEpotz(p,p->ef.Epot,x,y,z,p->ny,p->nz);
                    //p->ef.H_el_field[i] = 0.5 * dEpot_sq * p->ef.eps_arr[i];                  // redundant?
                    //sum_H_el += p->ef.H_el_field[i];                                          // redundant?
                    
                    p->ef.E_field[i] = -(dEpotx(p,p->ef.Epot,x,y,z,p->nx,p->ny,p->nz) +
                                         dEpoty(p,p->ef.Epot,x,y,z,p->ny,p->nz) + 
                                         dEpotz(p,p->ef.Epot,x,y,z,p->ny,p->nz));
#pragma acc loop seq
                    for (uint32_t m = 0; m < p->n_types; m++)
                    {   
                        soma_scalar_t phi_sum = 0.0;
                        soma_scalar_t diff_part = 0.0;
                        phi_sum += p->fields_unified[i+m*p->n_cells_local] * p->field_scaling_type[m];
#pragma acc loop seq
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

    //tests(p,p->ef.amt_iter);
    //p->ef.H_el = sum_H_el;
// #pragma acc update self(p->ef.omega_field_el)
//     if (p->ef.omega_field_el != p->ef.omega_field_el)
//     {
//         return -1;
//     }
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
    free(p->ef.kernel);
    free(p->ef.kernel_norm_field);
    free(p->ef.kernel_blur);
    free(p->ef.electrodes_conv);
    free(p->ef.eps_arr_conv);
    free(p->ef.pre_deriv_conv);
    free(p->ef.Epot_conv);
    free(p->ef.Epot_tmp_conv);

    return 0;
}

void tests(struct Phase *const p,uint64_t amt_iter)
{
    uint32_t counter = 1;
    update_self_electric_field(p);
    amt_iter += 1;
    if (p->time % counter == 0)
    {
        uint64_t sum_electrodes = 0;
        soma_scalar_t sum_Epot = 0.0;
        soma_scalar_t sum_eps_arr = 0.0;
        soma_scalar_t sum_pre_deriv = 0.0;
        soma_scalar_t sum_omega_field_el = 0.0;
        soma_scalar_t sum_omega_field = 0.0;

#pragma acc update host(p->omega_field_unified[0:(p->n_types*p->n_cells)])

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
