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
    p->ef.Epot = NULL;
    p->ef.Epot_tmp = NULL;
    p->ef.pre_deriv = NULL;
    p->ef.H_el_field = NULL;
    p->ef.H_el = NULL;
    

    //Quick exit if no electric field is present in the file
    if (!(H5Lexists(file_id, "/electric_field", H5P_DEFAULT) > 0))
        return 0;

    //Quick exit if no electrodes are present in the file
    if (!(H5Lexists(file_id, "/electrodes", H5P_DEFAULT) > 0))
        return 0;

    //Quick exit if no dielectric constants are specified
    if (!(H5Lexists(file_id, "/parameters/dielectric_constants", H5P_DEFAULT) > 0))
        return 0;

    //Allocate ef.eps_arr
    p->ef.eps_arr =
        (soma_scalar_t *) malloc((p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz *
                           sizeof(soma_scalar_t));

    //Allocate ef.Epot_tmp
    p->ef.Epot_tmp =
        (soma_scalar_t *) malloc((p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz *
                           sizeof(soma_scalar_t));

    //Allocate partial derivatives array
    p->ef.pre_deriv =
        (soma_scalar_t *) malloc((p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz * 6
                           sizeof(soma_scalar_t));

    //Allocate ef.H_el_field
    p->ef.H_el_field =
        (soma_scalar_t *) malloc((p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz *
                           sizeof(soma_scalar_t));

    //Allocate ef.H_el
    p->ef.H_el = (soma_scalar_t *) malloc(1.0 * sizeof(soma_scalar_t));

    //Read p->ef.eps
    p->ef.eps = (soma_scalar_t * const)malloc(p->n_types * sizeof(soma_scalar_t));
    if (p->ef.eps == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    status = read_hdf5(file_id, "/parameter/dielectric_constants", H5T_SOMA_NATIVE_SCALAR, plist_id, &(p->ef.eps));
    HDF5_ERROR_CHECK2(status, "/parameter/dielectric_constants");

    //Read electrode array (taken from area51)
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

    uint8_t *ptr = p->area51 + ghost_buffer_size;                                                                   // -> Ludwig? io.c (read_fields_hdf5 is *field opposed to implementation in poly_conversion.c)
    MPI_Isend(ptr, ghost_buffer_size, MPI_UINT8_T, left_neigh_rank, 0, p->info_MPI.SOMA_comm_sim, req + 0);
    ptr = p->area51 + ((p->nx / p->args.N_domains_arg) * p->ny * p->nz);
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
        }


    //Read electric potential field (adjusted from area51)
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
    if ((status = H5Dread(dataset, H5T_STD_U8LE, memspace, dataspace, plist_id, p->ef.Epot + ghost_buffer_size)) < 0)
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

    uint8_t *ptr = p->area51 + ghost_buffer_size;                                                                   //?
    MPI_Isend(ptr, ghost_buffer_size, MPI_UINT8_T, left_neigh_rank, 0, p->info_MPI.SOMA_comm_sim, req + 0);
    ptr = p->area51 + ((p->nx / p->args.N_domains_arg) * p->ny * p->nz);
    MPI_Isend(ptr, ghost_buffer_size, MPI_UINT8_T, right_neigh_rank, 1, p->info_MPI.SOMA_comm_sim, req + 1);

    ptr = p->ef.Epot;
    MPI_Irecv(ptr, ghost_buffer_size, MPI_UINT8_T, left_neigh_rank, 1, p->info_MPI.SOMA_comm_sim, req + 2);
    ptr = p->ef.Epot + ((p->nx / p->args.N_domains_arg) * p->ny * p->nz) + ghost_buffer_size;
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
        }
    return 0;
}



int write_electric_field_hdf5(const struct Phase *const p, const hid_t file_id, const hid_t plist_id)
{
    //Write dielectric constants data
    hsize_t n_types_size = p->n_types;
    status =
        write_hdf5(1, &n_types_size, file_id, "/parameter/dielectric_constants", H5T_SOMA_FILE_SCALAR, H5T_SOMA_NATIVE_SCALAR, plist_id,
                   p->ef.eps);
    HDF5_ERROR_CHECK2(status, "/parameter/dielectric_constants");

    //Write electrodes array
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
         H5Dwrite(dataset, H5T_NATIVE_UINT8, memspace, dataspace, plist_id, p->pc.array + ghost_buffer_size)) < 0)
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
        }

    //Write electric potential field array
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    const unsigned int ghost_buffer_size = p->args.domain_buffer_arg * p->ny * p->nz;

    const hsize_t hsize_dataspace[3] = { p->nx, p->ny, p->nz };
    hid_t dataspace = H5Screate_simple(3, hsize_dataspace, NULL);
    const hsize_t hsize_memspace[3] = { p->nx / p->args.N_domains_arg, p->ny, p->nz };
    hid_t memspace = H5Screate_simple(3, hsize_memspace, NULL);

    hid_t dataset =
        H5Dcreate2(file_id, "/electric_field", H5T_STD_U8LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    dataspace = H5Dget_space(dataset);
    const hsize_t offset[3] = { my_domain * (p->nx / p->args.N_domains_arg), 0, 0 };
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, hsize_memspace, NULL);

#if ( ENABLE_MPI == 1 )
    MPI_Barrier(p->info_MPI.SOMA_comm_world);
#endif                          //ENABLE_MPI

    if ((status =
         H5Dwrite(dataset, H5T_NATIVE_UINT8, memspace, dataspace, plist_id, p->ef.Epot + ghost_buffer_size)) < 0)
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
        }
    return 0;
}


void calc_dielectric_field(struct Phase *const p)
{
    soma_scalar_t tmp_frac_n;
    soma_scalar_t tmp_eps_n;
    soma_scalar_t eps_0 = 1.0;
    uint64_t x,y,z,i;
    for (x = 0; x < p->nx; x++)                
        for (y = 0; y < p->ny; y++)
            for (z = 0; z < p->nz; z++)
            {
                tmp_frac_n = 1.0;
                tmp_eps_n = 0.0;
                i = cell_coordinate_to_index(p,x,y,z);

                for (uint8_t k = 0; k < p->n_types; k++)
                {
                    tmp_frac_n -= p->fields_unified[i+k*p->n_cells_local] * p->field_scaling_type[k];
                    tmp_eps_n += p->fields_unified[i+k*p->n_cells_local] * p->field_scaling_type[k] * p->ef.eps[k];
                }

                p->ef.eps_arr[i] = eps_0 * tmp_frac_n + tmp_eps_n;
            }
}

void pre_derivatives(struct PHASE *const p)
{        
   uint64_t x,y,z,i,pos6;
   soma_scalar_t cr, epsx, epsy, epsz;
   soma_scalar_t f = 1.0/6.0;
   for (x=0; x < p->nx; x++)
      for (y=0; y < p->ny; y++)
         for (z=0; z < p->nz; z++)
            i = cell_coordinate_to_index(p,x,y,z);               //better use direct cell addressing w/o helper function and just one for-loop?
            pos6  =  x * p->ny * p->nz + y * p->nz + z * 6;
            cr = 1.0/(24.0 * p->ef.eps_arr[i]);

            epsx = cr * ( p->ef.eps_arr[cell_coordinate_to_index(p,x+1,y,z)] - p->ef.eps_arr[cell_coordinate_to_index(p,x-1,y,z)] );
            epsy = cr * ( p->ef.eps_arr[cell_coordinate_to_index(p,x,y+1,z)] - p->ef.eps_arr[cell_coordinate_to_index(p,x,y-1,z)] );
            epsz = cr * ( p->ef.eps_arr[cell_coordinate_to_index(p,x,y,z+1)] - p->ef.eps_arr[cell_coordinate_to_index(p,x,y,z-1)] );
            
            p->ef.pre_deriv[pos6 + 0] = epsi + f;
            p->ef.pre_deriv[pos6 + 1] =-epsi + f;
            p->ef.pre_deriv[pos6 + 2] = epsk + f;
            p->ef.pre_deriv[pos6 + 3] =-epsk + f;
            p->ef.pre_deriv[pos6 + 4] = epsj + f;
            p->ef.pre_deriv[pos6 + 5] =-epsj + f;
}

void iterate_field(struct PHASE *const p)
{
    uint64_t x,y,z;
    for (x=0; x < p->nx; x++)
      for (y=0; y < p->ny; y++)
         for (z=0; z < p->nz; z++)
            i = cell_coordinate_to_index(p,x,y,z);
            //Keep electric potential field constant at electrode position
            if (p->ef.electrodes[i] == 1)
            {
                p->ef.Epot_tmp[i] = p->ef.Epot[i]; 
            }
            else
            {
                p->ef.Epot_tmp[i] = p->ef.pre_deriv[pos6+0] * p->ef.Epot[cell_coordinate_to_index(p,x+1,y,z)] + p->ef.pre_deriv[pos6+1] * p->ef.Epot[cell_coordinate_to_index(p,x-1,y,z)] +
                 p->ef.pre_deriv[pos6+2] * p->ef.Epot[cell_coordinate_to_index(p,x,y+1,z)]  + p->ef.pre_deriv[pos6+3] * p->ef.Epot[cell_coordinate_to_index(p,x,y-1,z)] +
                 p->ef.pre_deriv[pos6+4] * p->ef.Epot[cell_coordinate_to_index(p,x,y,z+1)]  + p->ef.pre_deriv[pos6+5] * p->ef.Epot[cell_coordinate_to_index(p,x,y,z-1)];
            }
    memcpy(p->ef.Epot, p->ef.Epot_tmp, p->n_cells * sizeof(soma_scalar_t));         //correct or save directly in p->ef.Epot and omit ef.Epot_tmp?
}

soma_scalar_t Epot_deriv_sq(struct PHASE *const p, uint64_t x, uint64_t y, uint64_t z)
{
    dEpotx = (p->ef.Epot[(x+1) * p->ny * p->nz + y * p->nz + z] - p->ef.Epot[(x-1) * p->ny * p->nz + y * p->nz + z]) * p->nx / p->Lx; // * 0.5 ?!?!
    dEpoty = (p->ef.Epot[x * p->ny * p->nz + (y+1) * p->nz + z] - p->ef.Epot[x * p->ny * p->nz + (y-1) * p->nz + z]) * p->ny / p->Ly;
    dEpotz = (p->ef.Epot[x * p->ny * p->nz + y * p->nz + (z+1)] - p->ef.Epot[x * p->ny * p->nz + y * p->nz + (z-1)]) * p->nz / p->Lz;
    return(dEpotx * dEpotx + dEpoty * dEpoty + dEpotz * dEpotz);
}

int calc_electric_field_contr(struct PHASE *const p)
{
    uint64_t x,y,z,i;
    soma_scalar_t t = 1.0;
    calc_dielectric_field(p);
    pre_derivatives(p);

    //convergence criterion to ensure acceptable values for p->ef.Epot
    while (t <= 0.001)
    {
        for (uint64_t k = 0; k < p->n_cells; k++)
        {
            if (p->ef.electrodes[k] != 1)
                if (t < p->ef.Epot[k]) t = p->ef.Epot[k];               // failsafe?
        }

        iterate_field(p)
    }

    p->ef.H_el = 0;

    for (x=0; x < p->nx; x++)
        for (y=0; y < p->ny; y++)
            for (z=0; z < p->nz; z++)
            {
                i = cell_coordinate_to_index(p,x,y,z);
                p->ef.H_el_field[i] = 0.5 * Epot_deriv_sq(p,x,y,z) * p->ef.eps_arr[i];
                p->ef.H_el += p->ef.H_el_field[i]; //necessary?
            } 

    //free mem?
    return 0;
}

//copyin?
//copyout?
//self_update?
//free?