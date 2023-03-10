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

#include "polytype_conversion.h"

#include <stdlib.h>
#include <hdf5.h>
#include "phase.h"
#include "io.h"
#include "mesh.h"

//! \file polytype_conversion.c
//! \brief Implementation of polytype_conversion.h

int read_poly_conversion_hdf5(struct Phase *const p, const hid_t file_id, const hid_t plist_id)
{

    p->pc.deltaMC = 0;
    p->pc.array = NULL;
    p->pc.len_reactions = 0;
    p->pc.len_dependencies = 0;
    p->pc.rate = NULL;
    p->pc.input_type = NULL;
    p->pc.output_type = NULL;
    p->pc.reaction_end = NULL;
    p->pc.dependency_ntype = NULL;
    p->pc.dependency_type = NULL;
    p->pc.dependency_type_offset = NULL;
    //Quick exit if no poly conversion is present in the file
    if (!(H5Lexists(file_id, "/polyconversion", H5P_DEFAULT) > 0))
        return 0;

    herr_t status;
    unsigned int tmp_deltaMC;
    status = read_hdf5(file_id, "polyconversion/deltaMC", H5T_STD_U32LE, plist_id, &tmp_deltaMC);
    HDF5_ERROR_CHECK(status);

    //Quick exit if no conversion updates are scheduled
    if (tmp_deltaMC == 0)
        return 0;

#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    fprintf(stderr,
            "ERROR: %s: %d, Monotype Conversions are activated so polytype conversions do not work in the current implementation. Switch off the this feature and rerun.",
            __FILE__, __LINE__);
    return -1;
#endif                          //ENABLE_MONOTYPE_CONVERSIONS

    //Read the conversion list before reading the spatial array.
    const hid_t dset_input = H5Dopen(file_id, "/polyconversion/input_type", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_input);
    const hid_t dspace_input = H5Dget_space(dset_input);
    HDF5_ERROR_CHECK(dspace_input);
    const unsigned int ndims_input = H5Sget_simple_extent_ndims(dspace_input);
    if (ndims_input != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "polyconversion/input_type");
            return -1;
        }
    hsize_t dim_input;
    status = H5Sget_simple_extent_dims(dspace_input, &dim_input, NULL);
    HDF5_ERROR_CHECK(status);

    const hid_t dset_output = H5Dopen(file_id, "/polyconversion/output_type", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_output);
    const hid_t dspace_output = H5Dget_space(dset_output);
    HDF5_ERROR_CHECK(dspace_output);
    const unsigned int ndims_output = H5Sget_simple_extent_ndims(dspace_output);
    if (ndims_output != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "polyconversion/output_type");
            return -1;
        }
    hsize_t dim_output;
    status = H5Sget_simple_extent_dims(dspace_output, &dim_output, NULL);
    HDF5_ERROR_CHECK(status);

    const hid_t dset_end = H5Dopen(file_id, "/polyconversion/end", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_end);
    const hid_t dspace_end = H5Dget_space(dset_end);
    HDF5_ERROR_CHECK(dspace_end);
    const unsigned int ndims_end = H5Sget_simple_extent_ndims(dspace_end);
    if (ndims_end != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "polyconversion/end_type");
            return -1;
        }
    hsize_t dim_end;
    status = H5Sget_simple_extent_dims(dspace_end, &dim_end, NULL);
    HDF5_ERROR_CHECK(status);

    if (dim_input != dim_output)
        {
            fprintf(stderr,
                    "ERROR: %s:%d the length of the input and output type for poly conversions is not equal %d %d\n",
                    __FILE__, __LINE__, (int)dim_input, (int)dim_output);
            return -3;
        }
    if (dim_input != dim_end)
        {
            fprintf(stderr,
                    "ERROR: %s:%d the length of the input type and end for poly conversions is not equal %d %d\n",
                    __FILE__, __LINE__, (int)dim_input, (int)dim_end);
            return -3;
        }

    p->pc.input_type = (unsigned int *)malloc(dim_input * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->pc.input_type, dim_input * sizeof(unsigned int));
    p->pc.output_type = (unsigned int *)malloc(dim_output * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->pc.output_type, dim_output * sizeof(unsigned int));
    p->pc.reaction_end = (unsigned int *)malloc(dim_end * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->pc.reaction_end, dim_end * sizeof(unsigned int));

    p->pc.len_reactions = dim_input;

    //Actuablly read the polyconversion list data
    status = H5Dread(dset_input, H5T_STD_U32LE, H5S_ALL, H5S_ALL, plist_id, p->pc.input_type);
    HDF5_ERROR_CHECK(status);
    status = H5Sclose(dspace_input);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_input);
    HDF5_ERROR_CHECK(status);
    status = H5Dread(dset_output, H5T_STD_U32LE, H5S_ALL, H5S_ALL, plist_id, p->pc.output_type);
    HDF5_ERROR_CHECK(status);
    status = H5Sclose(dspace_output);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_output);
    HDF5_ERROR_CHECK(status);
    status = H5Dread(dset_end, H5T_STD_U32LE, H5S_ALL, H5S_ALL, plist_id, p->pc.reaction_end);
    HDF5_ERROR_CHECK(status);
    status = H5Sclose(dspace_end);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_end);
    HDF5_ERROR_CHECK(status);

    //Read the array information
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    const unsigned int ghost_buffer_size = p->args.domain_buffer_arg * p->ny * p->nz;

    const hsize_t hsize_memspace[3] = { p->nx / p->args.N_domains_arg, p->ny, p->nz };
    hid_t memspace = H5Screate_simple(3, hsize_memspace, NULL);
    hid_t dataset = H5Dopen2(file_id, "/polyconversion/array", H5P_DEFAULT);
    p->pc.array =
        (uint8_t *) malloc((p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz *
                           sizeof(uint8_t));
    if (p->pc.array == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    hid_t dataspace = H5Dget_space(dataset);
    const hsize_t offset[3] = { my_domain * (p->nx / p->args.N_domains_arg), 0, 0 };
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, hsize_memspace, NULL);

    if ((status = H5Dread(dataset, H5T_STD_U8LE, memspace, dataspace, plist_id, p->pc.array + ghost_buffer_size)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
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

    uint8_t *ptr = p->pc.array + ghost_buffer_size;
    MPI_Isend(ptr, ghost_buffer_size, MPI_UINT8_T, left_neigh_rank, 0, p->info_MPI.SOMA_comm_sim, req + 0);
    ptr = p->pc.array + ((p->nx / p->args.N_domains_arg) * p->ny * p->nz);
    MPI_Isend(ptr, ghost_buffer_size, MPI_UINT8_T, right_neigh_rank, 1, p->info_MPI.SOMA_comm_sim, req + 1);

    ptr = p->pc.array;
    MPI_Irecv(ptr, ghost_buffer_size, MPI_UINT8_T, left_neigh_rank, 1, p->info_MPI.SOMA_comm_sim, req + 2);
    ptr = p->pc.array + ((p->nx / p->args.N_domains_arg) * p->ny * p->nz) + ghost_buffer_size;
    MPI_Irecv(ptr, ghost_buffer_size, MPI_UINT8_T, right_neigh_rank, 0, p->info_MPI.SOMA_comm_sim, req + 3);

    MPI_Waitall(4, req, stat);
    MPI_Barrier(p->info_MPI.SOMA_comm_sim);
#endif                          //ENABLE_MPI

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

    //If rate is defined in the hdf5, partial conversions are activated and "rate", "n_density_dependencies", "density_dependencies" are read.
    if (!(H5Lexists(file_id, "/polyconversion/rate", H5P_DEFAULT) > 0))
        {
            //Enable the updat only if everything worked fine so far
            p->pc.deltaMC = tmp_deltaMC;
            return 0;
        }

    const hid_t dset_rate = H5Dopen(file_id, "/polyconversion/rate", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_rate);
    const hid_t dspace_rate = H5Dget_space(dset_rate);
    HDF5_ERROR_CHECK(dspace_rate);
    const unsigned int ndims_rate = H5Sget_simple_extent_ndims(dspace_rate);
    if (ndims_rate != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "polyconversion/rate");
            return -1;
        }
    hsize_t dim_rate;
    status = H5Sget_simple_extent_dims(dspace_rate, &dim_rate, NULL);
    HDF5_ERROR_CHECK(status);

    const hid_t dset_ndependency = H5Dopen(file_id, "/polyconversion/n_density_dependencies", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_ndependency);
    const hid_t dspace_ndependency = H5Dget_space(dset_ndependency);
    HDF5_ERROR_CHECK(dspace_ndependency);
    const unsigned int ndims_ndependency = H5Sget_simple_extent_ndims(dspace_ndependency);
    if (ndims_ndependency != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "polyconversion/n_density_dependencies");
            return -1;
        }
    hsize_t dim_ndependency;
    status = H5Sget_simple_extent_dims(dspace_ndependency, &dim_ndependency, NULL);
    HDF5_ERROR_CHECK(status);

    const hid_t dset_dependency = H5Dopen(file_id, "/polyconversion/density_dependencies", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_dependency);
    const hid_t dspace_dependency = H5Dget_space(dset_dependency);
    HDF5_ERROR_CHECK(dspace_dependency);
    const unsigned int ndims_dependency = H5Sget_simple_extent_ndims(dspace_dependency);
    if (ndims_dependency != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "polyconversion/density_dependencies");
            return -1;
        }
    hsize_t dim_dependency;
    status = H5Sget_simple_extent_dims(dspace_dependency, &dim_dependency, NULL);
    HDF5_ERROR_CHECK(status);

    if (dim_input != dim_rate)
        {
            fprintf(stderr,
                    "ERROR: %s:%d the length of the input type and rate for poly conversions is not equal %d %d\n",
                    __FILE__, __LINE__, (int)dim_input, (int)dim_rate);
            return -3;
        }
    if (dim_input != dim_ndependency)
        {
            fprintf(stderr,
                    "ERROR: %s:%d the length of the input type and rate for poly conversions is not equal %d %d\n",
                    __FILE__, __LINE__, (int)dim_input, (int)dim_rate);
            return -3;
        }

    p->pc.rate = (soma_scalar_t *) malloc(dim_rate * sizeof(soma_scalar_t));
    MALLOC_ERROR_CHECK(p->pc.rate, dim_rate * sizeof(soma_scalar_t));
    p->pc.dependency_ntype = (unsigned int *)malloc(dim_ndependency * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->pc.dependency_ntype, dim_ndependency * sizeof(unsigned int));
    p->pc.dependency_type_offset = (unsigned int *)malloc(dim_ndependency * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->pc.dependency_type_offset, dim_ndependency * sizeof(unsigned int));
    p->pc.dependency_type = (unsigned int *)malloc(dim_dependency * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->pc.dependency_type, dim_dependency * sizeof(unsigned int));

    p->pc.len_dependencies = dim_dependency;

    //Read rate and dependencies:
    status = H5Dread(dset_rate, H5T_SOMA_NATIVE_SCALAR, H5S_ALL, H5S_ALL, plist_id, p->pc.rate);
    HDF5_ERROR_CHECK(status);
    status = H5Sclose(dspace_rate);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_rate);
    HDF5_ERROR_CHECK(status);
    status = H5Dread(dset_ndependency, H5T_STD_U32LE, H5S_ALL, H5S_ALL, plist_id, p->pc.dependency_ntype);
    HDF5_ERROR_CHECK(status);
    status = H5Sclose(dspace_ndependency);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_ndependency);
    HDF5_ERROR_CHECK(status);
    p->pc.dependency_type_offset[0] = 0;
    for (unsigned int i = 1; i < dim_ndependency; i++)
        {
            p->pc.dependency_type_offset[i] = p->pc.dependency_type_offset[i - 1] + p->pc.dependency_ntype[i - 1];
        }
    if (p->pc.len_dependencies > 0)
        {
            status = H5Dread(dset_dependency, H5T_STD_U32LE, H5S_ALL, H5S_ALL, plist_id, p->pc.dependency_type);
            HDF5_ERROR_CHECK(status);
        }
    status = H5Sclose(dspace_dependency);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_dependency);
    HDF5_ERROR_CHECK(status);

    //Enable the updat only if everything worked fine so far
    p->pc.deltaMC = tmp_deltaMC;

    return 0;
}

int write_poly_conversion_hdf5(const struct Phase *const p, const hid_t file_id, const hid_t plist_id)
{
    //Quick exit for no poly conversions
    if (p->pc.deltaMC == 0)
        return 0;
    hid_t group = H5Gcreate2(file_id, "/polyconversion", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    herr_t status;
    const hsize_t one = 1;
    status =
        write_hdf5(1, &one, file_id, "/polyconversion/deltaMC", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id,
                   &(p->pc.deltaMC));
    HDF5_ERROR_CHECK(status);

    const hsize_t list_len = p->pc.len_reactions;
    status =
        write_hdf5(1, &list_len, file_id, "/polyconversion/input_type", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id,
                   p->pc.input_type);
    HDF5_ERROR_CHECK(status);
    status =
        write_hdf5(1, &list_len, file_id, "/polyconversion/output_type", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id,
                   p->pc.output_type);
    HDF5_ERROR_CHECK(status);
    status =
        write_hdf5(1, &list_len, file_id, "/polyconversion/end", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id,
                   p->pc.reaction_end);
    HDF5_ERROR_CHECK(status);

    if (p->pc.rate != NULL)
        {
            status =
                write_hdf5(1, &list_len, file_id, "/polyconversion/rate", H5T_SOMA_NATIVE_SCALAR,
                           H5T_SOMA_NATIVE_SCALAR, plist_id, p->pc.rate);
            HDF5_ERROR_CHECK(status);
            status =
                write_hdf5(1, &list_len, file_id, "/polyconversion/n_density_dependencies", H5T_STD_U32LE,
                           H5T_NATIVE_UINT, plist_id, p->pc.dependency_ntype);
            HDF5_ERROR_CHECK(status);
            const hsize_t dep_list_len = p->pc.len_dependencies;
            if (dep_list_len > 0)
                {
                    status =
                        write_hdf5(1, &dep_list_len, file_id, "/polyconversion/density_dependencies", H5T_STD_U32LE,
                                   H5T_NATIVE_UINT, plist_id, p->pc.dependency_type);
                    HDF5_ERROR_CHECK(status);
                }
            else
                {               //no dependencies --> empty dataset, so only create, don't write.
                    herr_t status;
                    const hid_t dataspace = H5Screate_simple(1, &dep_list_len, NULL);
                    HDF5_ERROR_CHECK(dataspace);
                    const hid_t dataset =
                        H5Dcreate2(file_id, "/polyconversion/density_dependencies", H5T_STD_U32LE, dataspace,
                                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
                    HDF5_ERROR_CHECK(dataset);
                    status = H5Sclose(dataspace);
                    HDF5_ERROR_CHECK(status);
                    status = H5Dclose(dataset);
                    HDF5_ERROR_CHECK(status);
                }
        }

    //Write the array to disk
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    const unsigned int ghost_buffer_size = p->args.domain_buffer_arg * p->ny * p->nz;

    const hsize_t hsize_dataspace[3] = { p->nx, p->ny, p->nz };
    hid_t dataspace = H5Screate_simple(3, hsize_dataspace, NULL);
    const hsize_t hsize_memspace[3] = { p->nx / p->args.N_domains_arg, p->ny, p->nz };
    hid_t memspace = H5Screate_simple(3, hsize_memspace, NULL);

    hid_t dataset =
        H5Dcreate2(file_id, "/polyconversion/array", H5T_STD_U8LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

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

    return 0;
}

int copyin_poly_conversion(struct Phase *p)
{
    if (p->pc.deltaMC != 0)
        {
#ifdef _OPENACC
            //The pc struct itself is part of the phase struct and is already present of the device
#pragma acc enter data copyin(p->pc.array[0:p->n_cells_local])
#pragma acc enter data copyin(p->pc.input_type[0:p->pc.len_reactions])
#pragma acc enter data copyin(p->pc.output_type[0:p->pc.len_reactions])
#pragma acc enter data copyin(p->pc.reaction_end[0:p->pc.len_reactions])
            if (p->pc.rate != NULL)
                {
#pragma acc enter data copyin(p->pc.rate[0:p->pc.len_reactions])
#pragma acc enter data copyin(p->pc.dependency_ntype[0:p->pc.len_reactions])
#pragma acc enter data copyin(p->pc.dependency_type_offset[0:p->pc.len_reactions])
#pragma acc enter data copyin(p->pc.dependency_type[0:p->pc.len_dependencies])
                }
#endif                          //_OPENACC
        }
    return 0;
}

int copyout_poly_conversion(struct Phase *p)
{
    if (p->pc.deltaMC != 0)
        {
#ifdef _OPENACC
#pragma acc exit data copyout(p->pc.array[0:p->n_cells_local])
#pragma acc exit data copyout(p->pc.input_type[0:p->pc.len_reactions])
#pragma acc exit data copyout(p->pc.output_type[0:p->pc.len_reactions])
#pragma acc exit data copyout(p->pc.reaction_end[0:p->pc.len_reactions])
            if (p->pc.rate != NULL)
                {
#pragma acc exit data copyout(p->pc.rate[0:p->pc.len_reactions])
#pragma acc exit data copyout(p->pc.dependency_ntype[0:p->pc.len_reactions])
#pragma acc exit data copyout(p->pc.dependency_type_offset[0:p->pc.len_reactions])
#pragma acc exit data copyout(p->pc.dependency_type[0:p->pc.len_dependencies])
                }
#endif                          //_OPENACC
        }
    return 0;
}

int update_self_poly_conversion(const struct Phase *const p)
{
    if (p->pc.deltaMC != 0)
        {
#ifdef _OPENACC
#pragma acc update self(p->pc.array[0:p->n_cells_local])
#pragma acc update self(p->pc.input_type[0:p->pc.len_reactions])
#pragma acc update self(p->pc.output_type[0:p->pc.len_reactions])
#pragma acc update self(p->pc.reaction_end[0:p->pc.len_reactions])
            if (p->pc.rate != NULL)
                {
#pragma acc update self(p->pc.rate[0:p->pc.len_reactions])
#pragma acc update self(p->pc.dependency_ntype[0:p->pc.len_reactions])
#pragma acc update self(p->pc.dependency_type_offset[0:p->pc.len_reactions])
#pragma acc update self(p->pc.dependency_type[0:p->pc.len_dependencies])
                }
#endif                          //_OPENACC
        }
    return 0;
}

int free_poly_conversion(struct Phase *p)
{
    free(p->pc.array);
    free(p->pc.input_type);
    free(p->pc.output_type);
    free(p->pc.reaction_end);
    free(p->pc.rate);
    free(p->pc.dependency_ntype);
    free(p->pc.dependency_type_offset);
    free(p->pc.dependency_type);

    return 0;
}

int convert_polytypes(struct Phase *p)
{
    //Quick exit for
    static unsigned int last_time = 0;
    if (last_time >= p->time)
        return 0;
    last_time = p->time;

    update_polymer_rcm(p);

    if (p->pc.rate == NULL)
        {
            return fully_convert_polytypes(p);
        }
    else
        {
            return partially_convert_polytypes(p);
        }
}

int fully_convert_polytypes(struct Phase *p)
{
    //Iterate all polymers and apply the reaction rules
#pragma acc parallel loop present(p[0:1])
#pragma omp parallel for
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            const Monomer rcm = p->polymers[poly].rcm;
            const uint64_t cell = coord_to_index(p, rcm.x, rcm.y, rcm.z);

            if (p->pc.array[cell] != 0)
                {
                    //Minus 1 because the index in array are shifted by 1
                    int i = p->pc.array[cell] - 1;
                    if (p->polymers[poly].type == p->pc.input_type[i])
                        p->polymers[poly].type = p->pc.output_type[i];
                    for (i++; !p->pc.reaction_end[i - 1]; i++)
                        {
                            if (p->polymers[poly].type == p->pc.input_type[i])
                                p->polymers[poly].type = p->pc.output_type[i];
                        }
                }
        }
    return 0;
}

int partially_convert_polytypes(struct Phase *p)
{
    //Iterate all polymers and apply the reaction rules
#pragma acc parallel loop present(p[0:1])
#pragma omp parallel for
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            Polymer *mypoly = p->polymers + poly;
            const Monomer rcm = mypoly->rcm;
            const uint64_t cell = coord_to_index(p, rcm.x, rcm.y, rcm.z);

            if (p->pc.array[cell] != 0)
                {
                    soma_scalar_t probability = 0.;
                    int i = p->pc.array[cell] - 2;
                    do
                        {
                            i++;
                            if (mypoly->type == p->pc.input_type[i])
                                {
                                    soma_scalar_t norm = 1 - probability;
                                    probability = p->pc.rate[i];
                                    for (unsigned int j = 0; j < p->pc.dependency_ntype[i]; j++)
                                        {
                                            unsigned int type_offset = p->pc.dependency_type_offset[i];
                                            unsigned int dependency_type = p->pc.dependency_type[type_offset + j];
                                            probability *=
                                                p->fields_unified[dependency_type * p->n_cells_local +
                                                                  cell] * p->field_scaling_type[dependency_type];
                                        }
                                    probability /= norm;
                                    soma_scalar_t random_number = soma_rng_soma_scalar(&(mypoly->poly_state), p);
                                    if (random_number < probability)
                                        {
                                            p->polymers[poly].type = p->pc.output_type[i];      //got modifiable lvalue compile error when using mypoly->type = ... and was not able to fix this otherwise.
                                            break;      //to continue with next polymer if conversion has taken place.
                                        }
                                    else
                                        {
                                            probability += (1 - norm);
                                        }
                                }
                    } while (!p->pc.reaction_end[i]);
                }
        }

    return 0;
}
