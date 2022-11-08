/* Copyright (C) 2016-2022 Gregor Ibbeken, Ludwig Schneider

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

#include "monotype_conversion.h"

#include <stdlib.h>
#include <hdf5.h>
#include "phase.h"
#include "io.h"
#include "mesh.h"
#include "mc.h"
#include "monomer.h"

//! \file monotype_conversion.c
//! \brief Implementation of monotype_conversion.h

int read_mono_conversion_hdf5(struct Phase *const p, const hid_t file_id, const hid_t plist_id)
{
#if ( ENABLE_MONOTYPE_CONVERSIONS == 0)
    //Check whether  mono conversion is present in the file:
    (void)p;                    //shut up compiler warnings
    (void)plist_id;
    if ((H5Lexists(file_id, "/monoconversion", H5P_DEFAULT) > 0))
        {
            fprintf(stderr,
                    "ERROR: %s:%d monotype conversions are defined in h5-file but compile option was not chosen for it.\n",
                    __FILE__, __LINE__);
            return -1;
        }

#else                           //ifdef ENABLE_MONOTYPE_CONVERSIONS
    p->mtc.deltaMC = 0;
    p->mtc.array = NULL;
    p->mtc.len_reactions = 0;
    p->mtc.len_dependencies = 0;
    p->mtc.rate = NULL;
    p->mtc.input_type = NULL;
    p->mtc.output_type = NULL;
    p->mtc.reaction_end = NULL;
    p->mtc.block_size = 1;
    p->mtc.dependency_ntype = NULL;
    p->mtc.dependency_type = NULL;
    p->mtc.dependency_type_offset = NULL;
    //Quick exit if no mono conversion is present in the file
    if (!(H5Lexists(file_id, "/monoconversion", H5P_DEFAULT) > 0))
        return 0;

    herr_t status;
    unsigned int tmp_deltaMC;
    status = read_hdf5(file_id, "monoconversion/deltaMC", H5T_STD_U32LE, plist_id, &tmp_deltaMC);
    HDF5_ERROR_CHECK(status);

    //Quick exit if no conversion updates are scheduled
    if (tmp_deltaMC == 0)
        return 0;

    //Read the conversion list before reading the spatial array.
    const hid_t dset_input = H5Dopen(file_id, "/monoconversion/input_type", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_input);
    const hid_t dspace_input = H5Dget_space(dset_input);
    HDF5_ERROR_CHECK(dspace_input);
    const unsigned int ndims_input = H5Sget_simple_extent_ndims(dspace_input);
    if (ndims_input != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "monoconversion/input_type");
            return -1;
        }
    hsize_t dim_input;
    status = H5Sget_simple_extent_dims(dspace_input, &dim_input, NULL);
    HDF5_ERROR_CHECK(status);

    const hid_t dset_output = H5Dopen(file_id, "/monoconversion/output_type", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_output);
    const hid_t dspace_output = H5Dget_space(dset_output);
    HDF5_ERROR_CHECK(dspace_output);
    const unsigned int ndims_output = H5Sget_simple_extent_ndims(dspace_output);
    if (ndims_output != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "monoconversion/output_type");
            return -1;
        }
    hsize_t dim_output;
    status = H5Sget_simple_extent_dims(dspace_output, &dim_output, NULL);
    HDF5_ERROR_CHECK(status);

    const hid_t dset_end = H5Dopen(file_id, "/monoconversion/end", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_end);
    const hid_t dspace_end = H5Dget_space(dset_end);
    HDF5_ERROR_CHECK(dspace_end);
    const unsigned int ndims_end = H5Sget_simple_extent_ndims(dspace_end);
    if (ndims_end != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "monoconversion/end_type");
            return -1;
        }
    hsize_t dim_end;
    status = H5Sget_simple_extent_dims(dspace_end, &dim_end, NULL);
    HDF5_ERROR_CHECK(status);

    unsigned int tmp_block_size;
    status = read_hdf5(file_id, "monoconversion/block_size", H5T_STD_U32LE, plist_id, &tmp_block_size);
    HDF5_ERROR_CHECK(status);

    if (dim_input != dim_output)
        {
            fprintf(stderr,
                    "ERROR: %s:%d the length of the input and output type for mono conversions is not equal %d %d\n",
                    __FILE__, __LINE__, (int)dim_input, (int)dim_output);
            return -3;
        }
    if (dim_input != dim_end)
        {
            fprintf(stderr,
                    "ERROR: %s:%d the length of the input type and end for mono conversions is not equal %d %d\n",
                    __FILE__, __LINE__, (int)dim_input, (int)dim_end);
            return -3;
        }

    p->mtc.input_type = (unsigned int *)malloc(dim_input * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->mtc.input_type, dim_input * sizeof(unsigned int));
    p->mtc.output_type = (unsigned int *)malloc(dim_output * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->mtc.output_type, dim_output * sizeof(unsigned int));
    p->mtc.reaction_end = (unsigned int *)malloc(dim_end * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->mtc.reaction_end, dim_end * sizeof(unsigned int));

    p->mtc.len_reactions = dim_input;

    //Actuablly read the monoconversion list data
    status = H5Dread(dset_input, H5T_STD_U32LE, H5S_ALL, H5S_ALL, plist_id, p->mtc.input_type);
    HDF5_ERROR_CHECK(status);
    status = H5Sclose(dspace_input);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_input);
    HDF5_ERROR_CHECK(status);
    status = H5Dread(dset_output, H5T_STD_U32LE, H5S_ALL, H5S_ALL, plist_id, p->mtc.output_type);
    HDF5_ERROR_CHECK(status);
    status = H5Sclose(dspace_output);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_output);
    HDF5_ERROR_CHECK(status);
    status = H5Dread(dset_end, H5T_STD_U32LE, H5S_ALL, H5S_ALL, plist_id, p->mtc.reaction_end);
    HDF5_ERROR_CHECK(status);
    status = H5Sclose(dspace_end);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_end);
    HDF5_ERROR_CHECK(status);

    p->mtc.block_size = tmp_block_size;

    //Read the array information
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    const unsigned int ghost_buffer_size = p->args.domain_buffer_arg * p->ny * p->nz;

    const hsize_t hsize_memspace[3] = { p->nx / p->args.N_domains_arg, p->ny, p->nz };
    hid_t memspace = H5Screate_simple(3, hsize_memspace, NULL);
    hid_t dataset = H5Dopen2(file_id, "/monoconversion/array", H5P_DEFAULT);
    p->mtc.array =
        (uint8_t *) malloc((p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz *
                           sizeof(uint8_t));
    if (p->mtc.array == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    hid_t dataspace = H5Dget_space(dataset);
    const hsize_t offset[3] = { my_domain * (p->nx / p->args.N_domains_arg), 0, 0 };
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, hsize_memspace, NULL);

    if ((status = H5Dread(dataset, H5T_STD_U8LE, memspace, dataspace, plist_id, p->mtc.array + ghost_buffer_size)) < 0)
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

    uint8_t *ptr = p->mtc.array + ghost_buffer_size;
    MPI_Isend(ptr, ghost_buffer_size, MPI_UINT8_T, left_neigh_rank, 0, p->info_MPI.SOMA_comm_sim, req + 0);
    ptr = p->mtc.array + ((p->nx / p->args.N_domains_arg) * p->ny * p->nz);
    MPI_Isend(ptr, ghost_buffer_size, MPI_UINT8_T, right_neigh_rank, 1, p->info_MPI.SOMA_comm_sim, req + 1);

    ptr = p->mtc.array;
    MPI_Irecv(ptr, ghost_buffer_size, MPI_UINT8_T, left_neigh_rank, 1, p->info_MPI.SOMA_comm_sim, req + 2);
    ptr = p->mtc.array + ((p->nx / p->args.N_domains_arg) * p->ny * p->nz) + ghost_buffer_size;
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
    if (!(H5Lexists(file_id, "/monoconversion/rate", H5P_DEFAULT) > 0))
        {
            //Enable the updat only if everything worked fine so far
            p->mtc.deltaMC = tmp_deltaMC;
            return 0;
        }

    const hid_t dset_rate = H5Dopen(file_id, "/monoconversion/rate", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_rate);
    const hid_t dspace_rate = H5Dget_space(dset_rate);
    HDF5_ERROR_CHECK(dspace_rate);
    const unsigned int ndims_rate = H5Sget_simple_extent_ndims(dspace_rate);
    if (ndims_rate != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "monoconversion/rate");
            return -1;
        }
    hsize_t dim_rate;
    status = H5Sget_simple_extent_dims(dspace_rate, &dim_rate, NULL);
    HDF5_ERROR_CHECK(status);

    const hid_t dset_ndependency = H5Dopen(file_id, "/monoconversion/n_density_dependencies", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_ndependency);
    const hid_t dspace_ndependency = H5Dget_space(dset_ndependency);
    HDF5_ERROR_CHECK(dspace_ndependency);
    const unsigned int ndims_ndependency = H5Sget_simple_extent_ndims(dspace_ndependency);
    if (ndims_ndependency != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "monoconversion/n_density_dependencies");
            return -1;
        }
    hsize_t dim_ndependency;
    status = H5Sget_simple_extent_dims(dspace_ndependency, &dim_ndependency, NULL);
    HDF5_ERROR_CHECK(status);

    const hid_t dset_dependency = H5Dopen(file_id, "/monoconversion/density_dependencies", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_dependency);
    const hid_t dspace_dependency = H5Dget_space(dset_dependency);
    HDF5_ERROR_CHECK(dspace_dependency);
    const unsigned int ndims_dependency = H5Sget_simple_extent_ndims(dspace_dependency);
    if (ndims_dependency != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",
                    __FILE__, __LINE__, "monoconversion/density_dependencies");
            return -1;
        }
    hsize_t dim_dependency;
    status = H5Sget_simple_extent_dims(dspace_dependency, &dim_dependency, NULL);
    HDF5_ERROR_CHECK(status);

    if (dim_input != dim_rate)
        {
            fprintf(stderr,
                    "ERROR: %s:%d the length of the input type and rate for mono conversions is not equal %d %d\n",
                    __FILE__, __LINE__, (int)dim_input, (int)dim_rate);
            return -3;
        }
    if (dim_input != dim_ndependency)
        {
            fprintf(stderr,
                    "ERROR: %s:%d the length of the input type and rate for mono conversions is not equal %d %d\n",
                    __FILE__, __LINE__, (int)dim_input, (int)dim_rate);
            return -3;
        }

    p->mtc.rate = (soma_scalar_t *) malloc(dim_rate * sizeof(soma_scalar_t));
    MALLOC_ERROR_CHECK(p->mtc.rate, dim_rate * sizeof(soma_scalar_t));
    p->mtc.dependency_ntype = (unsigned int *)malloc(dim_ndependency * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->mtc.dependency_ntype, dim_ndependency * sizeof(unsigned int));
    p->mtc.dependency_type_offset = (unsigned int *)malloc(dim_ndependency * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->mtc.dependency_type_offset, dim_ndependency * sizeof(unsigned int));
    p->mtc.dependency_type = (unsigned int *)malloc(dim_dependency * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->mtc.dependency_type, dim_dependency * sizeof(unsigned int));

    p->mtc.len_dependencies = dim_dependency;

    //Read rate and dependencies:
    status = H5Dread(dset_rate, H5T_SOMA_NATIVE_SCALAR, H5S_ALL, H5S_ALL, plist_id, p->mtc.rate);
    HDF5_ERROR_CHECK(status);
    status = H5Sclose(dspace_rate);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_rate);
    HDF5_ERROR_CHECK(status);
    status = H5Dread(dset_ndependency, H5T_STD_U32LE, H5S_ALL, H5S_ALL, plist_id, p->mtc.dependency_ntype);
    HDF5_ERROR_CHECK(status);
    status = H5Sclose(dspace_ndependency);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_ndependency);
    HDF5_ERROR_CHECK(status);
    p->mtc.dependency_type_offset[0] = 0;
    for (unsigned int i = 1; i < dim_ndependency; i++)
        {
            p->mtc.dependency_type_offset[i] = p->mtc.dependency_type_offset[i - 1] + p->mtc.dependency_ntype[i - 1];
        }
    if (p->mtc.len_dependencies > 0)
        {
            status = H5Dread(dset_dependency, H5T_STD_U32LE, H5S_ALL, H5S_ALL, plist_id, p->mtc.dependency_type);
            HDF5_ERROR_CHECK(status);
        }
    status = H5Sclose(dspace_dependency);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_dependency);
    HDF5_ERROR_CHECK(status);

    //Enable the updat only if everything worked fine
    p->mtc.deltaMC = tmp_deltaMC;

#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    return 0;
}

int write_mono_conversion_hdf5(const struct Phase *const p, const hid_t file_id, const hid_t plist_id)
{
    (void)p;                    //shut up compiler warnings
    (void)file_id;
    (void)plist_id;
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    //Quick exit for no mono conversions
    if (p->mtc.deltaMC == 0)
        return 0;
    hid_t group = H5Gcreate2(file_id, "/monoconversion", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    herr_t status;
    const hsize_t one = 1;
    status =
        write_hdf5(1, &one, file_id, "/monoconversion/deltaMC", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id,
                   &(p->mtc.deltaMC));
    HDF5_ERROR_CHECK(status);

    const hsize_t list_len = p->mtc.len_reactions;
    status =
        write_hdf5(1, &list_len, file_id, "/monoconversion/input_type", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id,
                   p->mtc.input_type);
    HDF5_ERROR_CHECK(status);
    status =
        write_hdf5(1, &list_len, file_id, "/monoconversion/output_type", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id,
                   p->mtc.output_type);
    HDF5_ERROR_CHECK(status);
    status =
        write_hdf5(1, &list_len, file_id, "/monoconversion/end", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id,
                   p->mtc.reaction_end);
    HDF5_ERROR_CHECK(status);
    status =
        write_hdf5(1, &one, file_id, "/monoconversion/block_size", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id,
                   &(p->mtc.block_size));
    HDF5_ERROR_CHECK(status);

    if (p->mtc.rate != NULL)
        {
            status =
                write_hdf5(1, &list_len, file_id, "/monoconversion/rate", H5T_SOMA_NATIVE_SCALAR,
                           H5T_SOMA_NATIVE_SCALAR, plist_id, p->mtc.rate);
            HDF5_ERROR_CHECK(status);
            status =
                write_hdf5(1, &list_len, file_id, "/monoconversion/n_density_dependencies", H5T_STD_U32LE,
                           H5T_NATIVE_UINT, plist_id, p->mtc.dependency_ntype);
            HDF5_ERROR_CHECK(status);
            const hsize_t dep_list_len = p->mtc.len_dependencies;
            if (dep_list_len > 0)
                {
                    status =
                        write_hdf5(1, &dep_list_len, file_id, "/monoconversion/density_dependencies", H5T_STD_U32LE,
                                   H5T_NATIVE_UINT, plist_id, p->mtc.dependency_type);
                    HDF5_ERROR_CHECK(status);
                }
            else
                {               //no dependencies --> empty dataset, so only create, don't write.
                    herr_t status;
                    const hid_t dataspace = H5Screate_simple(1, &dep_list_len, NULL);
                    HDF5_ERROR_CHECK(dataspace);
                    const hid_t dataset =
                        H5Dcreate2(file_id, "/monoconversion/density_dependencies", H5T_STD_U32LE, dataspace,
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
        H5Dcreate2(file_id, "/monoconversion/array", H5T_STD_U8LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    dataspace = H5Dget_space(dataset);
    const hsize_t offset[3] = { my_domain * (p->nx / p->args.N_domains_arg), 0, 0 };
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, hsize_memspace, NULL);

#if ( ENABLE_MPI == 1 )
    MPI_Barrier(p->info_MPI.SOMA_comm_world);
#endif                          //ENABLE_MPI

    if ((status =
         H5Dwrite(dataset, H5T_NATIVE_UINT8, memspace, dataspace, plist_id, p->mtc.array + ghost_buffer_size)) < 0)
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

#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    return 0;
}

int copyin_mono_conversion(struct Phase *p)
{
    (void)p;                    //shut up compiler warnings
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    if (p->mtc.deltaMC != 0)
        {
#ifdef _OPENACC
            //The mtc struct itself is part of the phase struct and is already present of the device
#pragma acc enter data copyin(p->mtc.array[0:p->n_cells_local])
#pragma acc enter data copyin(p->mtc.input_type[0:p->mtc.len_reactions])
#pragma acc enter data copyin(p->mtc.output_type[0:p->mtc.len_reactions])
#pragma acc enter data copyin(p->mtc.reaction_end[0:p->mtc.len_reactions])
            if (p->mtc.rate != NULL)
                {
#pragma acc enter data copyin(p->mtc.rate[0:p->mtc.len_reactions])
#pragma acc enter data copyin(p->mtc.dependency_ntype[0:p->mtc.len_reactions])
#pragma acc enter data copyin(p->mtc.dependency_type_offset[0:p->mtc.len_reactions])
#pragma acc enter data copyin(p->mtc.dependency_type[0:p->mtc.len_dependencies])
                }
#endif                          //_OPENACC
        }
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    return 0;
}

int copyout_mono_conversion(struct Phase *p)
{
    (void)p;                    //shut up compiler warnings
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    if (p->mtc.deltaMC != 0)
        {
#ifdef _OPENACC
#pragma acc exit data copyout(p->mtc.array[0:p->n_cells_local])
#pragma acc exit data copyout(p->mtc.input_type[0:p->mtc.len_reactions])
#pragma acc exit data copyout(p->mtc.output_type[0:p->mtc.len_reactions])
#pragma acc exit data copyout(p->mtc.reaction_end[0:p->mtc.len_reactions])
            if (p->mtc.rate != NULL)
                {
#pragma acc exit data copyout(p->mtc.rate[0:p->mtc.len_reactions])
#pragma acc exit data copyout(p->mtc.dependency_ntype[0:p->mtc.len_reactions])
#pragma acc exit data copyout(p->mtc.dependency_type_offset[0:p->mtc.len_reactions])
#pragma acc exit data copyout(p->mtc.dependency_type[0:p->mtc.len_dependencies])
                }
#endif                          //_OPENACC
        }
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    return 0;
}

int update_self_mono_conversion(const struct Phase *const p)
{
    (void)p;                    //shut up compiler warnings
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    if (p->mtc.deltaMC != 0)
        {
#ifdef _OPENACC
#pragma acc update self(p->mtc.array[0:p->n_cells_local])
#pragma acc update self(p->mtc.input_type[0:p->mtc.len_reactions])
#pragma acc update self(p->mtc.output_type[0:p->mtc.len_reactions])
#pragma acc update self(p->mtc.reaction_end[0:p->mtc.len_reactions])
            if (p->mtc.rate != NULL)
                {
#pragma acc update self(p->mtc.rate[0:p->mtc.len_reactions])
#pragma acc update self(p->mtc.dependency_ntype[0:p->mtc.len_reactions])
#pragma acc update self(p->mtc.dependency_type_offset[0:p->mtc.len_reactions])
#pragma acc update self(p->mtc.dependency_type[0:p->mtc.len_dependencies])
                }
#endif                          //_OPENACC
        }
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    return 0;
}

int free_mono_conversion(struct Phase *p)
{
    (void)p;                    //shut up compiler warnings
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    free(p->mtc.array);
    free(p->mtc.input_type);
    free(p->mtc.output_type);
    free(p->mtc.reaction_end);
    free(p->mtc.rate);
    free(p->mtc.dependency_ntype);
    free(p->mtc.dependency_type_offset);
    free(p->mtc.dependency_type);

#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    return 0;
}

int convert_monotypes(struct Phase *p)
{
    (void)p;                    //shut up compiler warnings
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    //Quick exit for
    static unsigned int last_time = 0;
    if (last_time >= p->time)
        return 0;
    last_time = p->time;

    if (p->mtc.rate == NULL)
        {
            return fully_convert_monotypes(p);
        }
    else
        {
            return partially_convert_monotypes(p);
        }
#else                           //ENABLE_MONOTYPE_CONVERSIONS
    return 0;
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
}

int fully_convert_monotypes(struct Phase *p)
{
    (void)p;                    //shut up compiler warnings
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    //Iterate all monomers and apply the reaction rules
#pragma acc parallel loop present(p[0:1])
#pragma omp parallel for
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            const Polymer *polymer = p->polymers + poly;
            const unsigned int N = p->poly_arch[p->poly_type_offset[polymer->type]];
            const unsigned int block_size = p->mtc.block_size;
            for (unsigned int mono = 0; mono < N; mono += block_size)
                {
                    Monomer block_pos = make_monomer(0., 0., 0.);
                    for (unsigned int block_index = 0; block_index < block_size; block_index++)
                        {
                            const Monomer pos = ((Monomer *) p->ph.beads.ptr)[polymer->bead_offset + mono + block_index];       //Read Monomer position
                            block_pos.x += pos.x;
                            block_pos.y += pos.y;
                            block_pos.z += pos.z;
                        }
                    block_pos.x /= block_size;
                    block_pos.y /= block_size;
                    block_pos.z /= block_size;
                    const uint64_t cell = coord_to_index(p, block_pos.x, block_pos.y, block_pos.z);

                    if (p->mtc.array[cell] != 0)
                        {
                            //Minus 1 because the index in array are shifted by 1
                            int i = p->mtc.array[cell] - 1;
                            unsigned int monotype = get_particle_type(p, poly, mono);
                            if (monotype == p->mtc.input_type[i])
                                for (unsigned int block_index = 0; block_index < block_size; block_index++)     //All monomers in the block are switched.
                                    ((uint8_t *) p->ph.monomer_types.ptr)[polymer->monomer_type_offset + mono +
                                                                          block_index] = p->mtc.output_type[i];
                            for (i++; !p->mtc.reaction_end[i - 1]; i++)
                                {
                                    if (monotype == p->mtc.input_type[i])
                                        for (unsigned int block_index = 0; block_index < block_size; block_index++)     //All monomers in the block are switched.
                                            ((uint8_t *) p->ph.monomer_types.ptr)[polymer->monomer_type_offset + mono +
                                                                                  block_index] = p->mtc.output_type[i];
                                }
                        }
                }
        }
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    return 0;
}

int partially_convert_monotypes(struct Phase *p)
{
    (void)p;                    //shut up compiler warnings
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    //Iterate all monomers and apply the reaction rules
#pragma acc parallel loop present(p[0:1])
#pragma omp parallel for
    for (uint64_t poly = 0; poly < p->n_polymers; poly++)
        {
            Polymer *mypoly = p->polymers + poly;
            const unsigned int N = p->poly_arch[p->poly_type_offset[mypoly->type]];

            for (unsigned int block = 0; block < N; block += p->mtc.block_size)
                {
                    Monomer block_rcm = make_monomer(0., 0., 0.);
                    for (unsigned int mono = 0; mono < p->mtc.block_size; mono++)
                        {
                            const Monomer pos = ((Monomer *) p->ph.beads.ptr)[mypoly->bead_offset + block + mono];      //Read Monomer position
                            block_rcm.x += pos.x;
                            block_rcm.y += pos.y;
                            block_rcm.z += pos.z;
                        }
                    block_rcm.x /= p->mtc.block_size;
                    block_rcm.y /= p->mtc.block_size;
                    block_rcm.y /= p->mtc.block_size;

                    const uint64_t cell = coord_to_index(p, block_rcm.x, block_rcm.y, block_rcm.z);

                    if (p->mtc.array[cell] != 0)
                        {
                            soma_scalar_t probability = 0.;
                            int i = p->mtc.array[cell] - 2;
                            //uint8_t monotype = ((uint8_t *) p->ph.monomer_types.ptr)[mypoly->monomer_type_offset + block];
                            unsigned int monotype = get_particle_type(p, poly, block);
                            do
                                {
                                    i++;
                                    if (monotype == p->mtc.input_type[i])
                                        {
                                            soma_scalar_t norm = 1 - probability;
                                            probability = p->mtc.rate[i];
                                            for (unsigned int j = 0; j < p->mtc.dependency_ntype[i]; j++)
                                                {
                                                    unsigned int type_offset = p->mtc.dependency_type_offset[i];
                                                    unsigned int dependency_type =
                                                        p->mtc.dependency_type[type_offset + j];
                                                    probability *=
                                                        p->fields_unified[dependency_type * p->n_cells_local +
                                                                          cell] *
                                                        p->field_scaling_type[dependency_type];
                                                }
                                            probability /= norm;
                                            soma_scalar_t random_number =
                                                soma_rng_soma_scalar(&(mypoly->poly_state), p);
                                            if (random_number < probability)
                                                {
                                                    for (unsigned int mono = 0; mono < p->mtc.block_size; mono++)
                                                        {
                                                            ((uint8_t *) p->ph.monomer_types.ptr)[mypoly->
                                                                                                  monomer_type_offset +
                                                                                                  block + mono] =
                                                                (uint8_t) p->mtc.output_type[i];
                                                        }
                                                    break;      //to continue with next block if conversion has taken place.
                                                }
                                            else
                                                {
                                                    probability += (1 - norm);
                                                }
                                        }
                            } while (!p->mtc.reaction_end[i]);
                        }
                }
        }

#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    return 0;
}
