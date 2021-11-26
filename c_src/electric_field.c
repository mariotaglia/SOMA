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
    hid_t status;
    
    //Quick exit if no electric field is present in the file
    if (!(H5Lexists(file_id, "/electric_field", H5P_DEFAULT) > 0))
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
        (soma_scalar_t *) malloc((p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz *
                           6 * sizeof(soma_scalar_t));
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
    // p->ef.iter_per_MC = (uint8_t * const)malloc(sizeof(uint8_t));
    // if (p->ef.iter_per_MC == NULL)
    //     {
    //         fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
    //         return -1;
    //     }
    status = read_hdf5(file_id, "/parameter/ef_iter_per_MC", H5T_NATIVE_UINT, plist_id, &(p->ef.iter_per_MC));
    HDF5_ERROR_CHECK2(status, "/parameter/ef_iter_per_MC");

    //Read p->ef.iter_limit
    // p->ef.iter_limit = (uint8_t * const)malloc(sizeof(uint8_t));
    // if (p->ef.iter_limit == NULL)
    //     {
    //         fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
    //         return -1;
    //     }
    status = read_hdf5(file_id, "/parameter/ef_iter_limit", H5T_NATIVE_UINT, plist_id, &(p->ef.iter_limit));
    HDF5_ERROR_CHECK2(status, "/parameter/ef_iter_limit");

    //Read p->ef.thresh_iter
    // p->ef.thresh_iter = (soma_scalar_t * const)malloc(sizeof(soma_scalar_t));
    // if (p->ef.thresh_iter == NULL)
    //     {
    //         fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
    //         return -1;
    //     }
    status = read_hdf5(file_id, "/parameter/ef_thresh_iter", H5T_SOMA_NATIVE_SCALAR, plist_id, &(p->ef.thresh_iter));
    HDF5_ERROR_CHECK2(status, "/parameter/ef_thresh_iter");

    //Read electrode array
    status = read_field_custom_hdf5(p, (void **) &(p->ef.electrodes), "/electrodes", sizeof(uint8_t), H5T_STD_U8LE, MPI_UINT8_T, file_id, plist_id);
    if (status != 0)
        {
            fprintf(stderr, "ERROR: failed to read electrode field %s:%d.\n", __FILE__, __LINE__);
            return status;
        }                                                                             // Substitute with function in io.c
    

    /*
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    const unsigned int ghost_buffer_size = p->args.domain_buffer_arg * p->ny * p->nz;

    const hsize_t hsize_memspace[3] = { p->nx / p->args.N_domains_arg, p->ny, p->nz };
    hid_t memspace = H5Screate_simple(3, hsize_memspace, NULL);
    hid_t dataset = H5Dopen2(file_id, "/electrodes", H5P_DEFAULT);
    p->ef.electrodes =
        (uint8_t *) malloc((p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz *
                           sizeof(uint8_t));
    if (p->ef.electrodes == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    hid_t dataspace = H5Dget_space(dataset);
    const hsize_t offset[3] = { my_domain * (p->nx / p->args.N_domains_arg), 0, 0 };
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, hsize_memspace, NULL);

    int status;
    if ((status = H5Dread(dataset, H5T_STD_U8LE, memspace, dataspace, plist_id, p->ef.electrodes + ghost_buffer_size)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

#if ( ENABLE_MPI == 1 )
    const int left_neigh_rank =
        (((my_domain - 1) + p->args.N_domains_arg) % p->args.N_domains_arg) * p->info_MPI.domain_size +
        p->info_MPI.domain_rank;
    const int right_neigh_rank =
        (((my_domain + 1) + p->args.N_domains_arg) % p->args.N_domains_arg) * p->info_MPI.domain_size +
        p->info_MPI.domain_rank;

    MPI_Request req[4];
    MPI_Status stat[4];

    uint8_t *ptr = p->ef.electrodes + ghost_buffer_size;                                                                  
    MPI_Isend(ptr, ghost_buffer_size, MPI_UINT8_T, left_neigh_rank, 0, p->info_MPI.SOMA_comm_sim, req + 0);
    ptr = p->ef.electrodes + ((p->nx / p->args.N_domains_arg) * p->ny * p->nz);
    MPI_Isend(ptr, ghost_buffer_size, MPI_UINT8_T, right_neigh_rank, 1, p->info_MPI.SOMA_comm_sim, req + 1);

    ptr = p->ef.electrodes;
    MPI_Irecv(ptr, ghost_buffer_size, MPI_UINT8_T, left_neigh_rank, 1, p->info_MPI.SOMA_comm_sim, req + 2);
    ptr = p->ef.electrodes + ((p->nx / p->args.N_domains_arg) * p->ny * p->nz) + ghost_buffer_size;
    MPI_Irecv(ptr, ghost_buffer_size, MPI_UINT8_T, right_neigh_rank, 0, p->info_MPI.SOMA_comm_sim, req + 3);

    MPI_Waitall(4, req, stat);
    MPI_Barrier(p->info_MPI.SOMA_comm_sim);
#endif                          //ENABLE_MPI

    if ((status = H5Sclose(dataspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    if ((status = H5Sclose(memspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    if ((status = H5Dclose(dataset)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }*/                                                                                                              

    //Read electric potential field
    status = read_field_custom_hdf5(p, (void **) &(p->ef.Epot), "/electric_field", sizeof(soma_scalar_t), H5T_SOMA_NATIVE_SCALAR,
                                            MPI_SOMA_SCALAR, file_id, plist_id);
    if (status != 0)
        {
            fprintf(stderr, "ERROR: failed to read electric potential field %s:%d.\n", __FILE__, __LINE__);
            return status;
        }

    /*
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    const unsigned int ghost_buffer_size = p->args.domain_buffer_arg * p->ny * p->nz;

    const hsize_t hsize_memspace[3] = { p->nx / p->args.N_domains_arg, p->ny, p->nz };
    hid_t memspace = H5Screate_simple(3, hsize_memspace, NULL);
    hid_t dataset = H5Dopen2(file_id, "/electric_field", H5P_DEFAULT);
    p->ef.Epot =
        (soma_scalar_t *) malloc((p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz *
                           sizeof(soma_scalar_t));
    if (p->ef.Epot == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    hid_t dataspace = H5Dget_space(dataset);
    const hsize_t offset[3] = { my_domain * (p->nx / p->args.N_domains_arg), 0, 0 };
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, hsize_memspace, NULL);

    int status;
    if ((status = H5Dread(dataset, H5T_SOMA_NATIVE_SCALAR, memspace, dataspace, plist_id, p->ef.Epot + ghost_buffer_size)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

#if ( ENABLE_MPI == 1 )
    const int left_neigh_rank =
        (((my_domain - 1) + p->args.N_domains_arg) % p->args.N_domains_arg) * p->info_MPI.domain_size +
        p->info_MPI.domain_rank;
    const int right_neigh_rank =
        (((my_domain + 1) + p->args.N_domains_arg) % p->args.N_domains_arg) * p->info_MPI.domain_size +
        p->info_MPI.domain_rank;

    MPI_Request req[4];
    MPI_Status stat[4];

    uint8_t *ptr = p->ef.Epot + ghost_buffer_size;                                                                   
    MPI_Isend(ptr, ghost_buffer_size, MPI_SOMA_SCALAR, left_neigh_rank, 0, p->info_MPI.SOMA_comm_sim, req + 0);
    ptr = p->ef.Epot + ((p->nx / p->args.N_domains_arg) * p->ny * p->nz);
    MPI_Isend(ptr, ghost_buffer_size, MPI_SOMA_SCALAR, right_neigh_rank, 1, p->info_MPI.SOMA_comm_sim, req + 1);

    ptr = p->ef.Epot;
    MPI_Irecv(ptr, ghost_buffer_size, MPI_SOMA_SCALAR, left_neigh_rank, 1, p->info_MPI.SOMA_comm_sim, req + 2);
    ptr = p->ef.Epot + ((p->nx / p->args.N_domains_arg) * p->ny * p->nz) + ghost_buffer_size;
    MPI_Irecv(ptr, ghost_buffer_size, MPI_SOMA_SCALAR, right_neigh_rank, 0, p->info_MPI.SOMA_comm_sim, req + 3);

    MPI_Waitall(4, req, stat);
    MPI_Barrier(p->info_MPI.SOMA_comm_sim);
#endif                          //ENABLE_MPI

    if ((status = H5Sclose(dataspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    if ((status = H5Sclose(memspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    if ((status = H5Dclose(dataset)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;  
        } */
    return 0;   
}



int write_electric_field_hdf5(const struct Phase *const p, const hid_t file_id, const hid_t plist_id)
{
    //Quick exit for no poly conversions
    if (p->ef.iter_per_MC == 0)
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

    //Write electrodes field
    status = write_field_custom_hdf5(p, (void **) &(p->ef.electrodes), "/electrodes", H5T_STD_U8LE, H5T_NATIVE_UINT8, file_id, plist_id);
    if (status != 0)
        {
            fprintf(stderr, "ERROR: failed to write electrode field %s:%d.\n", __FILE__, __LINE__);
            return status;
        }
    /*
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    const unsigned int ghost_buffer_size = p->args.domain_buffer_arg * p->ny * p->nz;

    const hsize_t hsize_dataspace[3] = { p->nx, p->ny, p->nz };
    hid_t dataspace = H5Screate_simple(3, hsize_dataspace, NULL);
    const hsize_t hsize_memspace[3] = { p->nx / p->args.N_domains_arg, p->ny, p->nz };
    hid_t memspace = H5Screate_simple(3, hsize_memspace, NULL);

    hid_t dataset =
        H5Dcreate2(file_id, "/electrodes", H5T_STD_U8LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    dataspace = H5Dget_space(dataset);
    const hsize_t offset[3] = { my_domain * (p->nx / p->args.N_domains_arg), 0, 0 };
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, hsize_memspace, NULL);

#if ( ENABLE_MPI == 1 )
    MPI_Barrier(p->info_MPI.SOMA_comm_world);
#endif                          //ENABLE_MPI

    if ((status =
         H5Dwrite(dataset, H5T_NATIVE_UINT8, memspace, dataspace, plist_id, p->ef.electrodes + ghost_buffer_size)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

    if ((status = H5Sclose(dataspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

    if ((status = H5Sclose(memspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

    if ((status = H5Dclose(dataset)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

    if ((status = H5Gclose(group)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }*/

    //Write electric potential field array
    status = write_field_custom_hdf5(p, (void **) &(p->ef.Epot), "/electric_field", H5T_SOMA_FILE_SCALAR, H5T_SOMA_NATIVE_SCALAR, file_id, plist_id);
    if (status != 0)
        {
            fprintf(stderr, "ERROR: failed to write electric potential field %s:%d.\n", __FILE__, __LINE__);
            return status;
        }
    /*
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    const unsigned int ghost_buffer_size = p->args.domain_buffer_arg * p->ny * p->nz;

    const hsize_t hsize_dataspace[3] = { p->nx, p->ny, p->nz };
    hid_t dataspace = H5Screate_simple(3, hsize_dataspace, NULL);
    const hsize_t hsize_memspace[3] = { p->nx / p->args.N_domains_arg, p->ny, p->nz };
    hid_t memspace = H5Screate_simple(3, hsize_memspace, NULL);

    hid_t dataset =
        H5Dcreate2(file_id, "/electric_field", H5T_SOMA_FILE_SCALAR, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    dataspace = H5Dget_space(dataset);
    const hsize_t offset[3] = { my_domain * (p->nx / p->args.N_domains_arg), 0, 0 };
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, hsize_memspace, NULL);

#if ( ENABLE_MPI == 1 )
    MPI_Barrier(p->info_MPI.SOMA_comm_world);
#endif                          //ENABLE_MPI

    if ((status =
         H5Dwrite(dataset, H5T_SOMA_NATIVE_SCALAR, memspace, dataspace, plist_id, p->ef.Epot + ghost_buffer_size)) < 0)   
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

    if ((status = H5Sclose(dataspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

    if ((status = H5Sclose(memspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

    if ((status = H5Dclose(dataset)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

    if ((status = H5Gclose(group)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }*/
    return 0;
}

uint64_t cell_to_index(struct Phase *const p, const uint64_t x, const uint64_t y, const uint64_t z)         // Missing: Error if domain decomposition active
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

void calc_dielectric_field(struct Phase *const p)
{
    soma_scalar_t tmp_frac_n;
    soma_scalar_t tmp_eps_n;
    soma_scalar_t eps_0 = 1.0;
    uint64_t i;
    
    for (i = 0; i < p-> n_cells; i++)
    {
        tmp_frac_n = 1.0;
        tmp_eps_n = 0.0;

        for (uint8_t k = 0; k < p->n_types; k++)
        {
            tmp_frac_n -= p->fields_unified[i+k*p->n_cells_local] * p->field_scaling_type[k];
            tmp_eps_n += p->fields_unified[i+k*p->n_cells_local] * p->field_scaling_type[k] * p->ef.eps[k];
        }

        p->ef.eps_arr[i] = eps_0 * tmp_frac_n + tmp_eps_n;
    }
}

void pre_derivatives(struct Phase *const p)
{
    uint64_t x,y,z,i,pos6;
    soma_scalar_t cr, epsx, epsy, epsz;
    soma_scalar_t f = 1.0/6.0;
    for (x=0; x < p->nx; x++)
        for (y=0; y < p->ny; y++)
            for (z=0; z < p->nz; z++)
            {
                i = cell_to_index(p,x,y,z);
                pos6 = (x * p->ny * p->nz + y * p->nz + z) * 6;
                cr = 1.0/(24.0 * p->ef.eps_arr[i]);

                epsx = cr * ( p->ef.eps_arr[cell_to_index(p,x+1,y,z)] - p->ef.eps_arr[cell_to_index(p,x-1,y,z)] );
                epsy = cr * ( p->ef.eps_arr[cell_to_index(p,x,y+1,z)] - p->ef.eps_arr[cell_to_index(p,x,y-1,z)] );
                epsz = cr * ( p->ef.eps_arr[cell_to_index(p,x,y,z+1)] - p->ef.eps_arr[cell_to_index(p,x,y,z-1)] );
                
                p->ef.pre_deriv[pos6 + 0] = epsx + f;
                p->ef.pre_deriv[pos6 + 1] =-epsx + f;
                p->ef.pre_deriv[pos6 + 2] = epsy + f;
                p->ef.pre_deriv[pos6 + 3] =-epsy + f;
                p->ef.pre_deriv[pos6 + 4] = epsz + f;
                p->ef.pre_deriv[pos6 + 5] =-epsz + f;
            }
}

soma_scalar_t iterate_field(struct Phase *const p)
{
    uint64_t x,y,z,i,pos6;
    soma_scalar_t dEpot_abs = 0.0;
    soma_scalar_t max = 0.0;

    for (x=0; x < p->nx; x++)
        for (y=0; y < p->ny; y++)
            for (z=0; z < p->nz; z++)
            {
                i = cell_to_index(p,x,y,z);
                //Keep electric potential field constant at electrode position
                if (p->ef.electrodes[i] == 1)
                {
                    p->ef.Epot_tmp[i] = p->ef.Epot[i]; 
                }
                else
                {   
                    pos6 = (x * p->ny * p->nz + y * p->nz + z) * 6;
                    p->ef.Epot_tmp[i] = 
                    p->ef.pre_deriv[pos6+0] * p->ef.Epot[cell_to_index(p,x+1,y,z)] + p->ef.pre_deriv[pos6+1] * p->ef.Epot[cell_to_index(p,x-1,y,z)] +
                    p->ef.pre_deriv[pos6+2] * p->ef.Epot[cell_to_index(p,x,y+1,z)] + p->ef.pre_deriv[pos6+3] * p->ef.Epot[cell_to_index(p,x,y-1,z)] +
                    p->ef.pre_deriv[pos6+4] * p->ef.Epot[cell_to_index(p,x,y,z+1)] + p->ef.pre_deriv[pos6+5] * p->ef.Epot[cell_to_index(p,x,y,z-1)];

                    // p->ef.Epot_tmp[i] = (1.0 / 6.0) * (
                    //                      p->ef.Epot[cell_to_index(p,x+1,y,z)] + p->ef.Epot[cell_to_index(p,x-1,y,z)] + 
                    //                      p->ef.Epot[cell_to_index(p,x,y+1,z)] + p->ef.Epot[cell_to_index(p,x,y-1,z)] +
                    //                      p->ef.Epot[cell_to_index(p,x,y,z+1)] + p->ef.Epot[cell_to_index(p,x,y,z-1)]) +
                    //                     (1.0 / (24.0 * p->ef.eps_arr[i])) * (
                    //                      (p->ef.eps_arr[cell_to_index(p,x+1,y,z)] - p->ef.eps_arr[cell_to_index(p,x-1,y,z)]) *                        
                    //                      (p->ef.Epot[cell_to_index(p,x+1,y,z)] - p->ef.Epot[cell_to_index(p,x-1,y,z)]) +
                    //                      (p->ef.eps_arr[cell_to_index(p,x,y+1,z)] - p->ef.eps_arr[cell_to_index(p,x,y-1,z)]) *                        
                    //                      (p->ef.Epot[cell_to_index(p,x,y+1,z)] - p->ef.Epot[cell_to_index(p,x,y-1,z)]) +
                    //                      (p->ef.eps_arr[cell_to_index(p,x,y,z+1)] - p->ef.eps_arr[cell_to_index(p,x,y,z-1)]) *                        
                    //                      (p->ef.Epot[cell_to_index(p,x,y,z+1)] - p->ef.Epot[cell_to_index(p,x,y,z-1)]) );
                    dEpot_abs = 0.0;

                    dEpot_abs = (d2Epotx(p,p->ef.Epot_tmp,x,y,z) + d2Epoty(p,p->ef.Epot_tmp,x,y,z) +
                                 d2Epotz(p,p->ef.Epot_tmp,x,y,z)) * p->ef.eps_arr[i];
                    //scalar product
                    dEpot_abs += (dEpotx(p,p->ef.Epot_tmp,x,y,z) * dEpotx(p,p->ef.eps_arr,x,y,z) +
                                  dEpoty(p,p->ef.Epot_tmp,x,y,z) * dEpoty(p,p->ef.eps_arr,x,y,z) +
                                  dEpotz(p,p->ef.Epot_tmp,x,y,z) * dEpotz(p,p->ef.eps_arr,x,y,z)) ;
                    if (max <  dEpot_abs) max = dEpot_abs;
                }
            }

    soma_scalar_t *tmp = p->ef.Epot_tmp;
    p->ef.Epot_tmp = p->ef.Epot;
    p->ef.Epot = tmp;

    return max;
}
/*
soma_scalar_t Epot_deriv_sq(struct PHASE *const p, const uint64_t x, const uint64_t y, const uint64_t z)
{
    dEpotx = (p->ef.Epot[(x+1) * p->ny * p->nz + y * p->nz + z] - p->ef.Epot[(x-1) * p->ny * p->nz + y * p->nz + z]) * 0.5 * p->nx / p->Lx;
    dEpoty = (p->ef.Epot[x * p->ny * p->nz + (y+1) * p->nz + z] - p->ef.Epot[x * p->ny * p->nz + (y-1) * p->nz + z]) * 0.5 * p->ny / p->Ly;
    dEpotz = (p->ef.Epot[x * p->ny * p->nz + y * p->nz + (z+1)] - p->ef.Epot[x * p->ny * p->nz + y * p->nz + (z-1)]) * 0.5 * p->nz / p->Lz;
    return(dEpotx * dEpotx + dEpoty * dEpoty + dEpotz * dEpotz);
}*/

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
    //Resolve periodic boundaries
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

int calc_electric_field_contr(struct Phase *const p)
{
    uint64_t x,y,z,i;
    soma_scalar_t t = p->ef.thresh_iter;
    soma_scalar_t dEpot_sq;
    uint64_t k = 0;
    calc_dielectric_field(p);
    pre_derivatives(p);
    
    //Convergence criterion to ensure acceptable values for p->ef.Epot
    while (t >= p->ef.thresh_iter && k < p->ef.iter_limit)
    {                                                   
        t = iterate_field(p);
        k +=1;
        // printf("%.5f\n",t);
    }

    if (k > p->ef.iter_per_MC) return -1; //Check if input for max iterations is exceeded

    p->ef.H_el = 0;

    for (x=0; x < p->nx; x++)
        for (y=0; y < p->ny; y++)
            for (z=0; z < p->nz; z++)
            {
                i = cell_to_index(p,x,y,z);
                // p->ef.H_el_field[i] = 0.5 * Epot_deriv_sq(p,x,y,z) * p->ef.eps_arr[i];
                dEpot_sq = dEpotx(p,p->ef.Epot,x,y,z) * dEpotx(p,p->ef.Epot,x,y,z) + dEpoty(p,p->ef.Epot,x,y,z) * dEpoty(p,p->ef.Epot,x,y,z) +
                            dEpotz(p,p->ef.Epot,x,y,z) * dEpotz(p,p->ef.Epot,x,y,z);
                p->ef.H_el_field[i] = 0.5 * dEpot_sq * p->ef.eps_arr[i];
                p->ef.H_el += p->ef.H_el_field[i]; //log in ana
                // printf("%.3f\n",p->ef.H_el_field[i]);
            }

    // tests(p,k);
    if (p->ef.H_el != p->ef.H_el)
    {
        return -1;
    }
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

    return 0;
}

//copyin?           // used to copy fields to GPU
//copyout? 
//self_update?


void tests(struct Phase *const p,uint64_t k)
{
    uint64_t sum_electrodes = 0;
    soma_scalar_t sum_Epot = 0.0;
    soma_scalar_t sum_eps_arr = 0.0;
    soma_scalar_t sum_pre_deriv = 0.0;
    for (uint64_t o = 0; o < p->n_cells; o++)
    {
        sum_electrodes += p->ef.electrodes[o];
        sum_Epot += p->ef.Epot[o];
        sum_eps_arr += p->ef.eps_arr[o];
        for (uint8_t l = 0; l < 6; l++)
            sum_pre_deriv += p->ef.pre_deriv[o+l];
    }

    printf("=====\n");
    printf("Iterations per MC: %ld\n",p->ef.iter_per_MC);
    printf("Iteration limit: %ld\n",p->ef.iter_limit);
    printf("Iteration threshold: %.5f\n",p->ef.thresh_iter);
    printf("Iterations for current step: %ld\n", k);
    printf("Sum electrode array: %ld\n",sum_electrodes);
    printf("Sum Epot array: %.3f\n",sum_Epot);
    printf("Sum epsilon array: %.3f\n",sum_eps_arr);
    printf("Sum pre deriv array: %.3f\n",sum_pre_deriv);
    printf("Iterations to solve EF: %ld\n", k);
    printf("Sum el contr: %.3f\n",p->ef.H_el);
    printf("=====\n");
}