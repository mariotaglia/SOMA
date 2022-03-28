/* Copyright (C) 2016-2021 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016-2017 Marcel Langenberg

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

/* */

//! \file io.c
//! \brief Implementation of io.h

#include "io.h"
#include <stdio.h>
#include <assert.h>
#if ( ENABLE_MPI == 1 )
#include <mpi.h>
#endif                          //ENABLE_MPI
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#ifdef _OPENMP
#include <omp.h>
#endif                          //_OPENMP
#include "mesh.h"
#include "cmdline.h"
#include "soma_config.h"
#include "polytype_conversion.h"
#include "monotype_conversion.h"
#include "mobility.h"
#include "io_old.h"

#if ( ENABLE_MPI == 1 )
#include <mpi.h>
#endif                          // ( ENABLE_MPI == 1 )
#include "hdf5.h"

int write_hdf5(const hsize_t ndims, const hsize_t * const dims, const hid_t file_id,
               const char *const name, const hid_t file_type, const hid_t mem_type,
               const hid_t plist_id, const void *const data)
{
    herr_t status;
    const hid_t dataspace = H5Screate_simple(ndims, dims, NULL);
    HDF5_ERROR_CHECK2(dataspace, name);
    const hid_t dataset = H5Dcreate2(file_id, name, file_type, dataspace,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_ERROR_CHECK2(dataset, name);

    status = H5Dwrite(dataset, mem_type, H5S_ALL, H5S_ALL, plist_id, data);
    HDF5_ERROR_CHECK2(status, name);

    status = H5Sclose(dataspace);
    HDF5_ERROR_CHECK2(status, name);
    status = H5Dclose(dataset);
    HDF5_ERROR_CHECK2(status, name);
    return status;
}

/*! Helper function to write area51 to the config HDF5 file.
    \private
    \param p Phase describing the system
    \param file_id File identifier of open HDF5 file.
    \param plist_id Access properties to use.
    \return Errorcode
*/
int write_area51_hdf5(const struct Phase *const p, const hid_t file_id, const hid_t plist_id)
{
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    const unsigned int ghost_buffer_size = p->args.domain_buffer_arg * p->ny * p->nz;

    const hsize_t hsize_dataspace[3] = { p->nx, p->ny, p->nz };
    hid_t dataspace = H5Screate_simple(3, hsize_dataspace, NULL);
    const hsize_t hsize_memspace[3] = { p->nx / p->args.N_domains_arg, p->ny, p->nz };
    hid_t memspace = H5Screate_simple(3, hsize_memspace, NULL);

    hid_t dataset = H5Dcreate2(file_id, "/area51", H5T_STD_U8LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    dataspace = H5Dget_space(dataset);
    const hsize_t offset[3] = { my_domain * (p->nx / p->args.N_domains_arg), 0, 0 };
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, hsize_memspace, NULL);

#if ( ENABLE_MPI == 1 )
    MPI_Barrier(p->info_MPI.SOMA_comm_world);
#endif                          //ENABLE_MPI
    int status;
    if ((status =
         H5Dwrite(dataset, H5T_NATIVE_UINT8, memspace, dataspace, plist_id, p->area51 + ghost_buffer_size)) < 0)
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

    return 0;
}

/*! Helper function to write scalar fields (external_field, umbrella_field) to the config HDF5 file.
    \private
    \param p Phase describing the system
    \param file_id File identifier of open HDF5 file.
    \param plist_id Access properties to use.
    \param field Pointer to the scalar field
    \param name Name of the field.
    \return Errorcode
*/
int write_field_hdf5(struct Phase *const p, const hid_t file_id,
                     const hid_t plist_id, const soma_scalar_t * const field, const char *name)
{
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    const unsigned int ghost_buffer_size = p->args.domain_buffer_arg * p->ny * p->nz;

    const hsize_t hsize_dataspace[4] = { p->n_types, p->nx, p->ny, p->nz };
    hid_t dataspace = H5Screate_simple(4, hsize_dataspace, NULL);
    const hsize_t hsize_memspace[4] = { 1, p->nx / p->args.N_domains_arg, p->ny, p->nz };
    hid_t memspace = H5Screate_simple(4, hsize_memspace, NULL);

    hid_t dataset = H5Dcreate2(file_id, name, H5T_SOMA_FILE_SCALAR, dataspace,
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    int status;
    for (unsigned int type = 0; type < p->n_types; type++)
        {
            if ((status = H5Sclose(dataspace)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, status);
                    return status;
                }

            dataspace = H5Dget_space(dataset);
            hsize_t tmp[4];
            H5Sget_simple_extent_dims(dataspace, tmp, NULL);
            const hsize_t offset[4] = { type, my_domain * (p->nx / p->args.N_domains_arg), 0, 0 };

            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, hsize_memspace, NULL);
            const soma_scalar_t *const ptr = field + ghost_buffer_size + type * p->n_cells_local;
            if ((status = H5Dwrite(dataset, H5T_SOMA_NATIVE_SCALAR, memspace, dataspace, plist_id, ptr)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, status);
                    return status;
                }
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

    return 0;
}

/*! Helper function to write the polytypes to disk.

\private
\param p Phase construct to write the polytypes from
\param file_id HDF5 file handle, opened and ready to write.
\param plist_id poperty list to be used for writing the file.
\return Errorcode
*/
int write_poly_type(const struct Phase *const p, hid_t file_id, const hid_t plist_id)
{
    hid_t status;
    //Determine the offset, for each process.
    uint64_t n_polymer_offset = 0;
#if ( ENABLE_MPI == 1 )
    //Cast for MPI_Scan, since some openmpi impl. need a non-const. version.
    MPI_Scan((uint64_t *) & (p->n_polymers), &n_polymer_offset, 1, MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);
    n_polymer_offset -= p->n_polymers;
#endif                          //ENABLE_MPI

    unsigned int *const poly_type = (unsigned int *const)malloc(p->n_polymers_storage * sizeof(unsigned int));
    if (poly_type == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    for (uint64_t i = 0; i < p->n_polymers; i++)
        poly_type[i] = p->polymers[i].type;

    hsize_t hsize_polymers_global = p->n_polymers_global;
    hid_t poly_type_dataspace = H5Screate_simple(1, &(hsize_polymers_global), NULL);
    hsize_t hsize_polymers = p->n_polymers;
    hid_t poly_type_memspace = H5Screate_simple(1, &(hsize_polymers), NULL);

    hid_t poly_type_dataset = H5Dcreate2(file_id, "/poly_type", H5T_STD_U32LE, poly_type_dataspace,
                                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    poly_type_dataspace = H5Dget_space(poly_type_dataset);
    hsize_t hsize_polymers_offset = n_polymer_offset;
    H5Sselect_hyperslab(poly_type_dataspace, H5S_SELECT_SET, &(hsize_polymers_offset), NULL, &(hsize_polymers), NULL);

    if ((status =
         H5Dwrite(poly_type_dataset, H5T_NATIVE_UINT, poly_type_memspace,
                  poly_type_dataspace, plist_id, poly_type)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    free(poly_type);
    if ((status = H5Sclose(poly_type_dataspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    if ((status = H5Sclose(poly_type_memspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    if ((status = H5Dclose(poly_type_dataset)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }
    return status;
}

/*! Helper function to write the poly_tags to disk.

\private
\param p Phase construct to write the polytypes from
\param file_id HDF5 file handle, opened and ready to write.
\param plist_id poperty list to be used for writing the file.
\return Errorcode
*/
int write_poly_tag(const struct Phase *const p, hid_t file_id, const hid_t plist_id)
{
    hid_t status;
    //Determine the offset, for each process.
    uint64_t n_polymer_offset = 0;
#if ( ENABLE_MPI == 1 )
    //Cast for MPI_Scan, since some openmpi impl. need a non-const. version.
    MPI_Scan((uint64_t *) & (p->n_polymers), &n_polymer_offset, 1, MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);
    n_polymer_offset -= p->n_polymers;
#endif                          //ENABLE_MPI

    uint64_t *const poly_tag = (uint64_t * const)malloc(p->n_polymers_storage * sizeof(uint64_t));
    if (poly_tag == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    for (uint64_t i = 0; i < p->n_polymers; i++)
        poly_tag[i] = p->polymers[i].tag;

    hsize_t hsize_polymers_global = p->n_polymers_global;
    hid_t poly_tag_dataspace = H5Screate_simple(1, &(hsize_polymers_global), NULL);
    hsize_t hsize_polymers = p->n_polymers;
    hid_t poly_tag_memspace = H5Screate_simple(1, &(hsize_polymers), NULL);

    hid_t poly_tag_dataset = H5Dcreate2(file_id, "/poly_tag", H5T_STD_U64LE, poly_tag_dataspace,
                                        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    poly_tag_dataspace = H5Dget_space(poly_tag_dataset);
    hsize_t hsize_polymers_offset = n_polymer_offset;
    H5Sselect_hyperslab(poly_tag_dataspace, H5S_SELECT_SET, &(hsize_polymers_offset), NULL, &(hsize_polymers), NULL);

    if ((status =
         H5Dwrite(poly_tag_dataset, H5T_NATIVE_UINT64, poly_tag_memspace, poly_tag_dataspace, plist_id, poly_tag)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    free(poly_tag);
    if ((status = H5Sclose(poly_tag_dataspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    if ((status = H5Sclose(poly_tag_memspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    if ((status = H5Dclose(poly_tag_dataset)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }
    return status;
}

/*! Helper function to write the heavy beads dataset to file. With multiple ranks, MPI-IO is used.

\private
\param p Phase struct intiliazied, ready to write.
\param file_id opened HDF5 file handle to write to.
\param plist_id property list for HDF5 writing.
\return Errorcode
*/
int write_beads(const struct Phase *const p, hid_t file_id, const hid_t plist_id)
{
    hid_t status;
    hsize_t hsize_beads_dataspace[2] = { p->num_all_beads, 3 };
    hid_t beads_dataspace = H5Screate_simple(2, hsize_beads_dataspace, NULL);
    hsize_t hsize_beads_memspace[2] = { p->num_all_beads_local, 3 };
    hid_t beads_memspace = H5Screate_simple(2, hsize_beads_memspace, NULL);

    hid_t beads_dataset = H5Dcreate2(file_id, "/beads", H5T_SOMA_FILE_SCALAR, beads_dataspace,
                                     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    uint64_t bead_offset = 0;
#if ( ENABLE_MPI == 1 )
    //Cast for MPI_Scan, since some openmpi impl. need a non-const. version.
    MPI_Scan((uint64_t *) & (p->num_all_beads_local), &bead_offset, 1, MPI_UINT64_T, MPI_SUM,
             p->info_MPI.SOMA_comm_sim);
    bead_offset -= p->num_all_beads_local;
#endif                          //ENABLE_MPI

    beads_dataspace = H5Dget_space(beads_dataset);
    hsize_t hsize_beads_offset[2] = { bead_offset, 0 };
    H5Sselect_hyperslab(beads_dataspace, H5S_SELECT_SET, hsize_beads_offset, NULL, hsize_beads_memspace, NULL);

    soma_scalar_t *const monomer_data =
        (soma_scalar_t * const)malloc(p->num_all_beads_local * 3 * sizeof(soma_scalar_t));
    if (monomer_data == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    //Transfer intenal bead data to SOMA data format
    uint64_t bead_mem_counter = 0;

    for (uint64_t i = 0; i < p->n_polymers; i++)
        {
            Monomer *beads = p->ph.beads.ptr;
            beads += p->polymers[i].bead_offset;
            const unsigned int N = p->poly_arch[p->poly_type_offset[p->polymers[i].type]];

            for (unsigned int j = 0; j < N; j++)
                {
                    assert(bead_mem_counter < p->num_all_beads_local);
                    monomer_data[bead_mem_counter * 3 + 0] = beads[j].x;
                    monomer_data[bead_mem_counter * 3 + 1] = beads[j].y;
                    monomer_data[bead_mem_counter * 3 + 2] = beads[j].z;

                    bead_mem_counter += 1;
                }
        }
    assert(bead_mem_counter == p->num_all_beads_local);

    if ((status =
         H5Dwrite(beads_dataset, H5T_SOMA_NATIVE_SCALAR, beads_memspace, beads_dataspace, plist_id, monomer_data)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    free(monomer_data);

    if ((status = H5Sclose(beads_dataspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    if ((status = H5Sclose(beads_memspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    if ((status = H5Dclose(beads_dataset)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    hsize_t hsize_mt_dataspace[1] = { p->num_all_beads };
    hid_t mt_dataspace = H5Screate_simple(1, hsize_mt_dataspace, NULL);
    hsize_t hsize_mt_memspace[1] = { p->num_all_beads_local };
    hid_t mt_memspace = H5Screate_simple(1, hsize_mt_memspace, NULL);

    hid_t mt_dataset = H5Dcreate2(file_id, "/monomer_types", H5T_STD_U8LE, mt_dataspace,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    mt_dataspace = H5Dget_space(mt_dataset);
    hsize_t hsize_mt_offset[1] = { bead_offset };       //bead_offset same as above
    H5Sselect_hyperslab(mt_dataspace, H5S_SELECT_SET, hsize_mt_offset, NULL, hsize_mt_memspace, NULL);

    uint8_t *const mt_data = (uint8_t * const)malloc(p->num_all_beads_local * sizeof(uint8_t));
    if (mt_data == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    //Transfer intenal bead data to SOMA data format
    bead_mem_counter = 0;

    for (uint64_t i = 0; i < p->n_polymers; i++)
        {
            uint8_t *mts = p->ph.monomer_types.ptr;
            mts += p->polymers[i].monomer_type_offset;
            const unsigned int N = p->poly_arch[p->poly_type_offset[p->polymers[i].type]];

            for (unsigned int j = 0; j < N; j++)
                {
                    assert(bead_mem_counter < p->num_all_beads_local);
                    mt_data[bead_mem_counter] = mts[j];

                    bead_mem_counter += 1;
                }
        }
    assert(bead_mem_counter == p->num_all_beads_local);

    if ((status = H5Dwrite(mt_dataset, H5T_NATIVE_UINT8, mt_memspace, mt_dataspace, plist_id, mt_data)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    free(mt_data);

    if ((status = H5Sclose(mt_dataspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    if ((status = H5Sclose(mt_memspace)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    if ((status = H5Dclose(mt_dataset)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    return status;
}

int write_config_hdf5(struct Phase *const p, const char *filename)
{
    if (strcmp(filename, "/dev/null") == 0)     //exit if no file should be written
        return 0;
    // Copy polymer data from device to host
    update_self_phase(p, 0);

    herr_t status;

    //Set up file access property list with parallel I/O access
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
#if ( ENABLE_MPI == 1 )
    if (p->info_MPI.sim_size > 1)
        H5Pset_fapl_mpio(plist_id, p->info_MPI.SOMA_comm_sim, MPI_INFO_NULL);
#endif                          //ENABLE_MPI

    //Create a new h5 file and overwrite the content.
    hid_t file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
#if ( ENABLE_MPI == 1 )
    if (p->info_MPI.sim_size > 1)
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif                          //ENABLE_MPI

    //Create a group for all parameter of the simulation
    hid_t parameter_group = H5Gcreate2(file_id, "/parameter", H5P_DEFAULT, H5P_DEFAULT,
                                       H5P_DEFAULT);

    const hsize_t one = 1;
    // Write the version number of the fileformat (hard coded) change if the version changes.
    // A version change make the file incompatible to read/write
    const unsigned int version = 1;
    status = write_hdf5(1, &one, file_id, "/version", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id, &version);
    HDF5_ERROR_CHECK2(status, "/version");

    //Write number of polymers
    status =
        write_hdf5(1, &one, file_id, "/parameter/n_polymers", H5T_STD_U64LE, H5T_NATIVE_UINT64, plist_id,
                   &(p->n_polymers_global));
    HDF5_ERROR_CHECK2(status, "/parameter/n_polymers");

    status =
        write_hdf5(1, &one, file_id, "/parameter/reference_Nbeads", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id,
                   &(p->reference_Nbeads));
    HDF5_ERROR_CHECK2(status, "/parameter/reference_Nbeads");

    status =
        write_hdf5(1, &one, file_id, "/parameter/n_types", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id, &(p->n_types));
    HDF5_ERROR_CHECK2(status, "/parameter/n_types");

    //Number of types
    hsize_t xn_dim[2] = { p->n_types, p->n_types };
    status =
        write_hdf5(2, xn_dim, file_id, "/parameter/xn", H5T_SOMA_FILE_SCALAR, H5T_SOMA_NATIVE_SCALAR, plist_id, p->xn);
    HDF5_ERROR_CHECK2(status, "/parameter/xn");

    //A data
    hsize_t n_types_size = p->n_types;
    status =
        write_hdf5(1, &n_types_size, file_id, "/parameter/A", H5T_SOMA_FILE_SCALAR, H5T_SOMA_NATIVE_SCALAR, plist_id,
                   p->A);
    HDF5_ERROR_CHECK2(status, "/parameter/A");

    //Convert field_scaling_type to density weight for writing
    for (unsigned int i = 0; i < p->n_types; i++)
        p->field_scaling_type[i] /= (p->n_accessible_cells / (soma_scalar_t) p->num_all_beads);
    status =
        write_hdf5(1, &n_types_size, file_id, "/parameter/density_weights", H5T_SOMA_FILE_SCALAR,
                   H5T_SOMA_NATIVE_SCALAR, plist_id, p->field_scaling_type);
    HDF5_ERROR_CHECK2(status, "/parameter/density_weights");
    //Convert the density weights back to field_scaling_type
    for (unsigned int i = 0; i < p->n_types; i++)
        p->field_scaling_type[i] *= (p->n_accessible_cells / (soma_scalar_t) p->num_all_beads);

    //time
    status = write_hdf5(1, &one, file_id, "/parameter/time", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id, &(p->time));
    HDF5_ERROR_CHECK2(status, "/parameter/time");

    status = write_hdf5(1, &one, file_id, "/parameter/hamiltonian",
                        H5T_STD_I32LE, H5T_NATIVE_INT, plist_id, &(p->hamiltonian));
    HDF5_ERROR_CHECK2(status, "parameter/hamiltonian");

    n_types_size = p->n_types;
    status =
        write_hdf5(1, &n_types_size, file_id, "/parameter/k_umbrella", H5T_SOMA_FILE_SCALAR, H5T_SOMA_NATIVE_SCALAR,
                   plist_id, p->k_umbrella);
    HDF5_ERROR_CHECK2(status, "parameter/k_umbrella");

    //Nx Ny Nz
    hsize_t three = 3;
    unsigned int nxyz[3] = { p->nx, p->ny, p->nz };
    status = write_hdf5(1, &three, file_id, "/parameter/nxyz", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id, nxyz);
    HDF5_ERROR_CHECK2(status, "/parameter/nxyz");

    //Lx Ly Lz
    soma_scalar_t lxyz[3] = { p->Lx, p->Ly, p->Lz };
    status =
        write_hdf5(1, &three, file_id, "/parameter/lxyz", H5T_SOMA_FILE_SCALAR, H5T_SOMA_NATIVE_SCALAR, plist_id, lxyz);
    HDF5_ERROR_CHECK2(status, "/parameter/lxyz");

    //p->harmonic_normb_variable_scale
    status =
        write_hdf5(1, &one, file_id, "/parameter/harmonic_normb_variable_scale", H5T_IEEE_F64LE, H5T_SOMA_NATIVE_SCALAR,
                   plist_id, &(p->harmonic_normb_variable_scale));
    HDF5_ERROR_CHECK2(status, "/parameter/harmonic_normb_variable_scale");

    //Polymer architecture
    //Number of polymer type
    status = write_hdf5(1, &one, file_id, "/parameter/n_poly_type", H5T_STD_U32LE,
                        H5T_NATIVE_UINT, plist_id, &(p->n_poly_type));
    HDF5_ERROR_CHECK2(status, "n_poly_types");
    //poly_type_offset
    hsize_t n_poly_type = p->n_poly_type;
    status = write_hdf5(1, &n_poly_type, file_id, "/parameter/poly_type_offset", H5T_STD_I32LE,
                        H5T_NATIVE_INT, plist_id, p->poly_type_offset);
    HDF5_ERROR_CHECK2(status, "poly_type_offset");
    //poly_arch_length
    status = write_hdf5(1, &one, file_id, "/parameter/poly_arch_length", H5T_STD_U32LE,
                        H5T_NATIVE_UINT, plist_id, &(p->poly_arch_length));
    HDF5_ERROR_CHECK2(status, "poly_arch_length");

    //poly_arch
    const hsize_t arch_length = p->poly_arch_length;
    //Warning this is a hidden cast from UINT32 to INT32
    status = write_hdf5(1, &arch_length, file_id, "/parameter/poly_arch", H5T_STD_I32LE,
                        H5T_NATIVE_INT, plist_id, p->poly_arch);
    HDF5_ERROR_CHECK2(status, "poly_arch");

    if (p->cm_a != NULL)
        {
            status = write_hdf5(1, &n_poly_type, file_id, "/parameter/cm_a", H5T_SOMA_FILE_SCALAR,
                                H5T_SOMA_NATIVE_SCALAR, plist_id, p->cm_a);
            HDF5_ERROR_CHECK2(status, "/parameter/cm_a");
        }

    //Close parameter group
    if ((status = H5Gclose(parameter_group)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }
    status = write_poly_type(p, file_id, plist_id);
    HDF5_ERROR_CHECK2(status, "write poly type");

    status = write_poly_tag(p, file_id, plist_id);
    HDF5_ERROR_CHECK2(status, "write poly tag");

    status = write_beads(p, file_id, plist_id);
    HDF5_ERROR_CHECK2(status, "write beads");

    if (p->area51)
        {
            hid_t status = write_area51_hdf5(p, file_id, plist_id);
            if (status != 0)
                {
                    fprintf(stderr, "ERROR: %s:%d cannot write area51.\n", __FILE__, __LINE__);
                    return status;
                }
        }

    if (p->external_field_unified)
        {
            hid_t status = write_field_hdf5(p, file_id, plist_id, p->external_field_unified, "/external_field");
            if (status != 0)
                {
                    fprintf(stderr, "ERROR: %s:%d cannot write external_field.\n", __FILE__, __LINE__);
                    return status;
                }
        }

    if (p->umbrella_field)
        {
            hid_t status = write_field_hdf5(p, file_id, plist_id, p->umbrella_field, "/umbrella_field");
            if (status != 0)
                {
                    fprintf(stderr, "ERROR: %s:%d cannot write umbrella_field.\n", __FILE__, __LINE__);
                    return status;
                }
        }

    status = write_poly_conversion_hdf5(p, file_id, plist_id);
    if (status != 0)
        {
            fprintf(stderr, "ERROR: %s:%d writing the polytype conversion\n", __FILE__, __LINE__);
            return status;
        }

    status = write_mono_conversion_hdf5(p, file_id, plist_id);
    if (status != 0)
        {
            fprintf(stderr, "ERROR: %s:%d writing the monotype conversion\n", __FILE__, __LINE__);
            return status;
        }

    status = write_mobility_hdf5(p, file_id, plist_id);
    if (status != 0)
        {
            fprintf(stderr, "ERROR: %s:%d writing the mobility\n", __FILE__, __LINE__);
            return status;
        }

    H5Pclose(plist_id);
    if ((status = H5Fclose(file_id)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    //Synchronize the closing of the file
    MPI_Barrier(p->info_MPI.SOMA_comm_sim);

    // Variable length data is not supported with parallel IO.
    // So we just reopen the file on one rank and write it anyways.
    if (p->info_MPI.sim_rank == 0)
        {
            hid_t file_id = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
            status = file_id;
            if (status < 0)
                {
                    fprintf(stderr, "ERROR: HDF5 %s:%s:%d\n", __func__, __FILE__, __LINE__);
                    return -1;
                }
            add_self_documentation_to_hdf5(&(p->sd), file_id, H5P_DEFAULT);
            if ((status = H5Fclose(file_id)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, status);
                    return status;
                }
        }
    return status;
}

int read_hdf5(const hid_t file_id, const char *const name, const hid_t mem_type, const hid_t plist_id, void *const data)
{
    herr_t status;
    const hid_t dataset = H5Dopen2(file_id, name, H5P_DEFAULT);
    HDF5_ERROR_CHECK2(dataset, name);
    status = H5Dread(dataset, mem_type, H5S_ALL, H5S_ALL, plist_id, data);
    HDF5_ERROR_CHECK2(status, name);
    status = H5Dclose(dataset);
    HDF5_ERROR_CHECK2(status, name);
    return status;
}

/*! Helper function to read area51 from the config HDF5 file.
    \private
    \param p Phase describing the system
    \param file_id File identifier of open HDF5 file.
    \param plist_id Access properties to use.
    \return Errorcode
*/
int read_area51_hdf5(struct Phase *const p, const hid_t file_id, const hid_t plist_id)
{
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    const unsigned int ghost_buffer_size = p->args.domain_buffer_arg * p->ny * p->nz;

    const hsize_t hsize_memspace[3] = { p->nx / p->args.N_domains_arg, p->ny, p->nz };
    hid_t memspace = H5Screate_simple(3, hsize_memspace, NULL);
    hid_t dataset = H5Dopen2(file_id, "/area51", H5P_DEFAULT);
    p->area51 =
        (uint8_t *) malloc((p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz *
                           sizeof(uint8_t));
    if (p->area51 == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    hid_t dataspace = H5Dget_space(dataset);
    const hsize_t offset[3] = { my_domain * (p->nx / p->args.N_domains_arg), 0, 0 };
    H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, hsize_memspace, NULL);

    int status;
    if ((status = H5Dread(dataset, H5T_STD_U8LE, memspace, dataspace, plist_id, p->area51 + ghost_buffer_size)) < 0)
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

    uint8_t *ptr = p->area51 + ghost_buffer_size;
    MPI_Isend(ptr, ghost_buffer_size, MPI_UINT8_T, left_neigh_rank, 0, p->info_MPI.SOMA_comm_sim, req + 0);
    ptr = p->area51 + ((p->nx / p->args.N_domains_arg) * p->ny * p->nz);
    MPI_Isend(ptr, ghost_buffer_size, MPI_UINT8_T, right_neigh_rank, 1, p->info_MPI.SOMA_comm_sim, req + 1);

    ptr = p->area51;
    MPI_Irecv(ptr, ghost_buffer_size, MPI_UINT8_T, left_neigh_rank, 1, p->info_MPI.SOMA_comm_sim, req + 2);
    ptr = p->area51 + ((p->nx / p->args.N_domains_arg) * p->ny * p->nz) + ghost_buffer_size;
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

/*! Helper function to read scalar fields (external_field, umbrella_field) from the config HDF5 file.
    \private
    \param p Phase describing the system
    \param file_id File identifier of open HDF5 file.
    \param plist_id Access properties to use.
    \param field Pointer to the scalar field
    \param name Name of the field.
    \return Errorcode
*/
int read_field_hdf5(const struct Phase *const p, const hid_t file_id, const hid_t plist_id, soma_scalar_t ** field,
                    const char *name)
{
    const unsigned int my_domain = p->info_MPI.sim_rank / p->info_MPI.domain_size;
    const unsigned int ghost_buffer_size = p->args.domain_buffer_arg * p->ny * p->nz;

    const hsize_t hsize_memspace[4] = { 1, p->nx / p->args.N_domains_arg, p->ny, p->nz };
    hid_t memspace = H5Screate_simple(4, hsize_memspace, NULL);
    hid_t dataset = H5Dopen2(file_id, name, H5P_DEFAULT);
    const uint64_t n_cells_local = (p->nx / p->args.N_domains_arg + 2 * p->args.domain_buffer_arg) * p->ny * p->nz;
    *field = (soma_scalar_t *) malloc(n_cells_local * p->n_types * sizeof(soma_scalar_t));
    if (*field == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    const int left_neigh_rank =
        (((my_domain - 1) + p->args.N_domains_arg) % p->args.N_domains_arg) * p->info_MPI.domain_size +
        p->info_MPI.domain_rank;
    const int right_neigh_rank =
        (((my_domain + 1) + p->args.N_domains_arg) % p->args.N_domains_arg) * p->info_MPI.domain_size +
        p->info_MPI.domain_rank;

    int status;
    for (unsigned int type = 0; type < p->n_types; type++)
        {
            hid_t dataspace = H5Dget_space(dataset);
            const hsize_t offset[4] = { type, my_domain * (p->nx / p->args.N_domains_arg), 0, 0 };
            H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, hsize_memspace, NULL);

            soma_scalar_t *ptr = *field + ghost_buffer_size + (type * n_cells_local);
            if ((status = H5Dread(dataset, H5T_SOMA_NATIVE_SCALAR, memspace, dataspace, plist_id, ptr)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }

#if ( ENABLE_MPI == 1 )
            MPI_Request req[4];
            MPI_Status stat[4];

            ptr = *field + ghost_buffer_size + type * n_cells_local;
            MPI_Isend(ptr, ghost_buffer_size, MPI_SOMA_SCALAR, left_neigh_rank, 0 + 2 * type, p->info_MPI.SOMA_comm_sim,
                      req + 0);
            ptr = *field + ((p->nx / p->args.N_domains_arg) * p->ny * p->nz) + type * n_cells_local;
            MPI_Isend(ptr, ghost_buffer_size, MPI_SOMA_SCALAR, right_neigh_rank, 1 + 2 * type,
                      p->info_MPI.SOMA_comm_sim, req + 1);

            ptr = *field + type * n_cells_local;
            MPI_Irecv(ptr, ghost_buffer_size, MPI_SOMA_SCALAR, left_neigh_rank, 1 + 2 * type, p->info_MPI.SOMA_comm_sim,
                      req + 2);
            ptr = *field + ((p->nx / p->args.N_domains_arg) * p->ny * p->nz) + ghost_buffer_size + type * n_cells_local;
            MPI_Irecv(ptr, ghost_buffer_size, MPI_SOMA_SCALAR, right_neigh_rank, 0 + 2 * type,
                      p->info_MPI.SOMA_comm_sim, req + 3);

            MPI_Waitall(4, req, stat);
            MPI_Barrier(p->info_MPI.SOMA_comm_sim);
#endif                          // ( ENABLE_MPI == 1 )

            if ((status = H5Sclose(dataspace)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }
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

/*! Helper function to read the beads information from the input files of version 1.

\private
\param p Phase struct which is going to be initialized.
\param file_id opened HDF5 file handle to read from.
\param plist_id HDF5 property list to read with.
\return Errorcode.
*/
int read_beads1(struct Phase *const p, const hid_t file_id, const hid_t plist_id)
{
    hid_t status;
    //temporary load poly type to calculate polymer distribution
    unsigned int *poly_type_tmp = (unsigned int *)malloc(p->n_polymers_global * sizeof(unsigned int));
    if (poly_type_tmp == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    status = read_hdf5(file_id, "/poly_type", H5T_NATIVE_UINT, plist_id, poly_type_tmp);
    HDF5_ERROR_CHECK2(status, "poly_type tmp global read");
    p->num_all_beads = 0;
    for (uint64_t i = 0; i < p->n_polymers_global; i++)
        {
            const unsigned int N = p->poly_arch[p->poly_type_offset[poly_type_tmp[i]]];
            p->num_all_beads += N;
        }

    //Distribute the polymers to the different cores according to bead length.
    uint64_t av_beads_per_rank = p->num_all_beads / p->info_MPI.sim_size;
    uint64_t bead_offset = 0;
    uint64_t poly_offset = 0;
    uint64_t my_num_beads = 0;
    uint64_t my_num_poly = 0;

    // distribute the number of polymer to the ranks
    // store only the offset of all previous beads and polymers
    // Loop only the previous ranks and the rank itself
    for (int sim_rank = 0; sim_rank < p->info_MPI.sim_rank + 1; sim_rank++)
        {
            my_num_poly = 0;
            my_num_beads = 0;
            while (my_num_beads < av_beads_per_rank && poly_offset + my_num_poly < p->n_polymers_global)
                {
                    const unsigned int N = p->poly_arch[p->poly_type_offset[poly_type_tmp[poly_offset + my_num_poly]]];
                    my_num_beads += N;
                    my_num_poly += 1;
                }
            bead_offset += my_num_beads;
            poly_offset += my_num_poly;

            assert(bead_offset <= p->num_all_beads);
        }
    //Correct the offset by the number of the current rank
    bead_offset -= my_num_beads;
    poly_offset -= my_num_poly;
    //Make sure all polymers and beads are read
    if (p->info_MPI.sim_rank == p->info_MPI.sim_size - 1)
        {
            my_num_poly = p->n_polymers_global - poly_offset;
            my_num_beads = p->num_all_beads - bead_offset;
        }

    /* printf("DEBUG %s:%d (remove) sim rank: %d bead_offset: %d poly_offset: %d num_beads: %d num_poly: %d\n", __FILE__, */
    /*        __LINE__, p->info_MPI.sim_rank, (int)bead_offset, (int)poly_offset, (int)my_num_beads, (int)my_num_poly); */

    p->n_polymers = my_num_poly;
    p->n_polymers_storage = p->n_polymers;

    p->polymers = (Polymer *) malloc(p->n_polymers_storage * sizeof(Polymer));
    if (p->polymers == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    p->num_all_beads_local = my_num_beads;

    uint64_t *poly_tag = NULL;
    if (H5Lexists(file_id, "/poly_tag", H5P_DEFAULT) > 0)
        {
            poly_tag = (uint64_t *) malloc(p->n_polymers * sizeof(uint64_t));
            if (poly_tag == NULL)
                {
                    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
                    return -1;
                }
            hsize_t hsize_polymers = p->n_polymers;
            hsize_t hsize_polymers_offset = poly_offset;
            hid_t poly_tag_memspace = H5Screate_simple(1, &(hsize_polymers), NULL);
            hid_t poly_tag_dataset = H5Dopen2(file_id, "/poly_tag", H5P_DEFAULT);
            hid_t poly_tag_dataspace = H5Dget_space(poly_tag_dataset);
            H5Sselect_hyperslab(poly_tag_dataspace, H5S_SELECT_SET, &(hsize_polymers_offset), NULL, &(hsize_polymers),
                                NULL);

            if ((status =
                 H5Dread(poly_tag_dataset, H5T_NATIVE_UINT64, poly_tag_memspace, poly_tag_dataspace, plist_id,
                         poly_tag)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }

            if ((status = H5Sclose(poly_tag_memspace)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }

            if ((status = H5Dclose(poly_tag_dataset)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }

            if ((status = H5Sclose(poly_tag_dataspace)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }
        }                       //poly_tag

    init_soma_memory(&(p->ph.beads), p->num_all_beads_local, sizeof(Monomer));
    init_soma_memory(&(p->ph.msd_beads), p->num_all_beads_local, sizeof(Monomer));

    for (uint64_t i = 0; i < p->n_polymers; i++)
        {
            p->polymers[i].type = poly_type_tmp[poly_offset + i];
            if (poly_tag == NULL)       //No tag infor available because older file format
                p->polymers[i].tag = poly_offset + i;
            else
                p->polymers[i].tag = poly_tag[i];

            const unsigned int N = p->poly_arch[p->poly_type_offset[p->polymers[i].type]];
            p->polymers[i].bead_offset = get_new_soma_memory_offset(&(p->ph.beads), N);
            if (p->polymers[i].bead_offset == UINT64_MAX)
                {
                    fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__, __LINE__,
                            p->info_MPI.world_rank, p->n_polymers);
                    return -1;
                }
            p->polymers[i].msd_bead_offset = get_new_soma_memory_offset(&(p->ph.msd_beads), N);
            if (p->polymers[i].msd_bead_offset == UINT64_MAX)
                {
                    fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__, __LINE__,
                            p->info_MPI.world_rank, p->n_polymers);
                    return -1;
                }
        }
    free(poly_type_tmp);

    if (H5Lexists(file_id, "/beads", H5P_DEFAULT) > 0)
        {

            //  Read the bead data into a tempory array that matches the file layout
            soma_scalar_t *const monomer_data =
                (soma_scalar_t * const)malloc(p->num_all_beads_local * 3 * sizeof(soma_scalar_t));
            if (monomer_data == NULL)
                {
                    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
                    return -1;
                }

            hsize_t hsize_beads_memspace[2] = { p->num_all_beads_local, 3 };
            hid_t beads_memspace = H5Screate_simple(2, hsize_beads_memspace, NULL);
            hid_t beads_dataset = H5Dopen2(file_id, "/beads", H5P_DEFAULT);
            hid_t beads_dataspace = H5Dget_space(beads_dataset);
            hsize_t hsize_beads_offset[2] = { bead_offset, 0 };
            H5Sselect_hyperslab(beads_dataspace, H5S_SELECT_SET, hsize_beads_offset, NULL, hsize_beads_memspace, NULL);

            if ((status =
                 H5Dread(beads_dataset, H5T_SOMA_NATIVE_SCALAR, beads_memspace, beads_dataspace, plist_id,
                         monomer_data)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }

            if ((status = H5Sclose(beads_dataspace)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }

            if ((status = H5Sclose(beads_memspace)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }

            if ((status = H5Dclose(beads_dataset)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }

            //Transfer the memory of the temporary array into the polymer structs.
            uint64_t bead_mem_counter = 0;
            for (uint64_t i = 0; i < p->n_polymers; i++)
                {
                    Monomer *beads = p->ph.beads.ptr;
                    beads += p->polymers[i].bead_offset;
                    Monomer *msd_beads = p->ph.msd_beads.ptr;
                    msd_beads += p->polymers[i].msd_bead_offset;

                    const unsigned int N = p->poly_arch[p->poly_type_offset[p->polymers[i].type]];
                    for (unsigned int j = 0; j < N; j++)
                        {
                            assert(bead_mem_counter < p->num_all_beads_local);
                            beads[j].x = monomer_data[bead_mem_counter * 3 + 0];
                            beads[j].y = monomer_data[bead_mem_counter * 3 + 1];
                            beads[j].z = monomer_data[bead_mem_counter * 3 + 2];

                            msd_beads[j].x = monomer_data[bead_mem_counter * 3 + 0];
                            msd_beads[j].y = monomer_data[bead_mem_counter * 3 + 1];
                            msd_beads[j].z = monomer_data[bead_mem_counter * 3 + 2];

                            bead_mem_counter += 1;
                        }
                }
            assert(bead_mem_counter == p->num_all_beads_local);
            free(monomer_data);

            p->bead_data_read = true;
        }
    else
        p->bead_data_read = false;

    /* Read monomer types from file if present: */
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    init_soma_memory(&(p->ph.monomer_types), p->num_all_beads_local, sizeof(uint8_t));
    for (uint64_t i = 0; i < p->n_polymers; i++)
        {
            const unsigned int N = p->poly_arch[p->poly_type_offset[p->polymers[i].type]];
            p->polymers[i].monomer_type_offset = get_new_soma_memory_offset(&(p->ph.monomer_types), N);
            if (p->polymers[i].monomer_type_offset == UINT64_MAX)
                {
                    fprintf(stderr, "ERROR: invalid memory alloc %s:%d rank %d, n_poly %lu\n", __FILE__, __LINE__,
                            p->info_MPI.world_rank, p->n_polymers);
                    return -1;
                }
        }
    if (H5Lexists(file_id, "/monomer_types", H5P_DEFAULT) > 0)
        {

            //  Read the monomer type data into a tempory array that matches the file layout
            uint8_t *const mt_data = (uint8_t * const)malloc(p->num_all_beads_local * sizeof(uint8_t));
            if (mt_data == NULL)
                {
                    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
                    return -1;
                }

            hsize_t hsize_mt_memspace[1] = { p->num_all_beads_local };
            hid_t mt_memspace = H5Screate_simple(1, hsize_mt_memspace, NULL);
            hid_t mt_dataset = H5Dopen2(file_id, "/monomer_types", H5P_DEFAULT);
            hid_t mt_dataspace = H5Dget_space(mt_dataset);
            hsize_t hsize_mt_offset[1] = { bead_offset };       //bead_offset should be the same as monomer_type_offset
            H5Sselect_hyperslab(mt_dataspace, H5S_SELECT_SET, hsize_mt_offset, NULL, hsize_mt_memspace, NULL);

            if ((status = H5Dread(mt_dataset, H5T_NATIVE_UINT8, mt_memspace, mt_dataspace, plist_id, mt_data)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }

            if ((status = H5Sclose(mt_dataspace)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }

            if ((status = H5Sclose(mt_memspace)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }

            if ((status = H5Dclose(mt_dataset)) < 0)
                {
                    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                            p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
                    return status;
                }

            //Transfer the memory of the temporary array into the polymer structs.
            uint64_t mt_mem_counter = 0;
            for (uint64_t i = 0; i < p->n_polymers; i++)
                {
                    uint8_t *mts = p->ph.monomer_types.ptr;
                    mts += p->polymers[i].monomer_type_offset;

                    const unsigned int N = p->poly_arch[p->poly_type_offset[p->polymers[i].type]];
                    for (unsigned int j = 0; j < N; j++)
                        {
                            assert(mt_mem_counter < p->num_all_beads_local);
                            mts[j] = mt_data[mt_mem_counter];
                            mt_mem_counter += 1;
                        }
                }
            assert(mt_mem_counter == p->num_all_beads_local);
            free(mt_data);
            p->mt_data_read = true;

        }
    else
        {
            p->mt_data_read = false;
        }
#else                           //ENABLE_MONOTYPE_CONVERSIONS
    if (H5Lexists(file_id, "/monomer_types", H5P_DEFAULT) > 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n", p->info_MPI.world_rank, __FILE__, __LINE__,
                    (int)status);
            return -1;
        }
    else
        {
            p->mt_data_read = false;
        }
#endif                          //ENABLE_MONOTYPE_CONVERSIONS
    return 0;
}

/*! General helper function to read the beads information from the input file if available.

This function wraps the different actual reading functions and selects the appropriate one according to the file version.
\private
\param p Phase struct which is going to be initialized.
\param file_id opened HDF5 file handle to read from.
\param plist_id HDF5 property list to read with.
\return Errorcode.
*/
int read_beads(struct Phase *p, const hid_t file_id, const hid_t plist_id)
{
    hid_t status;
    unsigned int version = 0;
    if (H5Lexists(file_id, "/version", H5P_DEFAULT) > 0)
        {
            status = read_hdf5(file_id, "/version", H5T_NATIVE_UINT, plist_id, &version);
            HDF5_ERROR_CHECK2(status, "version");
        }
    switch (version)
        {
        case 0:
            return read_beads0(p, file_id, plist_id);
            break;
        case 1:
            return read_beads1(p, file_id, plist_id);
            break;
        default:
            fprintf(stderr, "ERROR: invalid file version %u  %s:%d\n", version, __FILE__, __LINE__);
            return -1;
        }
    return 0;
}

int read_config_hdf5(struct Phase *const p, const char *filename)
{
    herr_t status;
    //Set up file access property list with parallel I/O access
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
#if ( ENABLE_MPI == 1 )
    if (p->info_MPI.sim_size > 1)
        H5Pset_fapl_mpio(plist_id, p->info_MPI.SOMA_comm_sim, MPI_INFO_NULL);
#endif                          //ENABLE_MPI

    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
#if ( ENABLE_MPI == 1 )
    if (p->info_MPI.sim_size > 1)
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
#endif                          //ENABLE_MPI

    status = read_hdf5(file_id, "/parameter/n_polymers", H5T_NATIVE_UINT64, plist_id, &(p->n_polymers_global));
    HDF5_ERROR_CHECK2(status, "/parameter/n_polymers");

    status = read_hdf5(file_id, "/parameter/reference_Nbeads", H5T_NATIVE_UINT, plist_id, &(p->reference_Nbeads));
    HDF5_ERROR_CHECK2(status, "/parameter/reference_Nbeads");

    status = read_hdf5(file_id, "/parameter/n_types", H5T_NATIVE_UINT, plist_id, &(p->n_types));
    HDF5_ERROR_CHECK2(status, "/parameter/n_types");

    p->xn = (soma_scalar_t * const)malloc(p->n_types * p->n_types * sizeof(soma_scalar_t));
    if (p->xn == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    status = read_hdf5(file_id, "/parameter/xn", H5T_SOMA_NATIVE_SCALAR, plist_id, p->xn);
    HDF5_ERROR_CHECK2(status, "/parameter/xn");

    //A array for the diffusivity of the particles
    p->A = (soma_scalar_t * const)malloc(p->n_types * sizeof(soma_scalar_t));
    if (p->A == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    // read p->A
    status = read_hdf5(file_id, "/parameter/A", H5T_SOMA_NATIVE_SCALAR, plist_id, p->A);
    HDF5_ERROR_CHECK2(status, "/parameter/A");

    // read p->time
    status = read_hdf5(file_id, "/parameter/time", H5T_NATIVE_UINT, plist_id, &(p->time));
    HDF5_ERROR_CHECK2(status, "/parameter/time");

    // read nx ny nz
    p->hamiltonian = SCMF0;
    //Don't break old configurations.
    if (H5Lexists(file_id, "/parameter/hamiltonian", H5P_DEFAULT) > 0)
        {
            status = read_hdf5(file_id, "/parameter/hamiltonian", H5T_NATIVE_INT, plist_id, &(p->hamiltonian));
            HDF5_ERROR_CHECK2(status, "/parameter/hamiltonian");
        }

    p->k_umbrella = (soma_scalar_t * const)malloc(p->n_types * sizeof(soma_scalar_t));
    if (p->k_umbrella == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    memset(p->k_umbrella, 0, p->n_types * sizeof(soma_scalar_t));       //Default 0
    if (H5Lexists(file_id, "/parameter/k_umbrella", H5P_DEFAULT) > 0)
        {
            status = read_hdf5(file_id, "/parameter/k_umbrella", H5T_SOMA_NATIVE_SCALAR, plist_id, p->k_umbrella);
            HDF5_ERROR_CHECK2(status, "parameter/k_umbrella");
        }

    unsigned int nxyz[3];
    status = read_hdf5(file_id, "/parameter/nxyz", H5T_NATIVE_UINT, plist_id, nxyz);
    HDF5_ERROR_CHECK2(status, "/parameter/nxyz");
    p->nx = nxyz[0];
    p->ny = nxyz[1];
    p->nz = nxyz[2];

    if (p->nx % p->args.N_domains_arg != 0)
        {
            fprintf(stderr, "ERROR: %s:%d\n\t"
                    "The nx %d number is not divible by the number of domains %d\n",
                    __FILE__, __LINE__, p->nx, p->args.N_domains_arg);
            return -3;
        }

    // read lx ly lz
    soma_scalar_t lxyz[3];
    status = read_hdf5(file_id, "/parameter/lxyz", H5T_SOMA_NATIVE_SCALAR, plist_id, lxyz);
    HDF5_ERROR_CHECK2(status, "/parameter/lxyz");
    p->Lx = lxyz[0];
    p->Ly = lxyz[1];
    p->Lz = lxyz[2];

    //Read in the polymer architectures.
    //Number of polymer type
    status = read_hdf5(file_id, "/parameter/n_poly_type", H5T_NATIVE_UINT, plist_id, &(p->n_poly_type));
    HDF5_ERROR_CHECK2(status, "n_poly_type");
    //poly_type_offset
    p->poly_type_offset = (int *)malloc(p->n_poly_type * sizeof(int));
    if (p->poly_type_offset == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    status = read_hdf5(file_id, "/parameter/poly_type_offset", H5T_NATIVE_INT, plist_id, p->poly_type_offset);
    HDF5_ERROR_CHECK2(status, "poly_type_offset");
    //poly_arch_length
    status = read_hdf5(file_id, "/parameter/poly_arch_length", H5T_NATIVE_UINT, plist_id, &(p->poly_arch_length));
    HDF5_ERROR_CHECK2(status, "poly_arch_length");
    //poly_arch
    p->poly_arch = (uint32_t *) malloc(p->poly_arch_length * sizeof(uint32_t));
    if (p->poly_arch == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    status = read_hdf5(file_id, "/parameter/poly_arch", H5T_NATIVE_INT32, plist_id, p->poly_arch);
    HDF5_ERROR_CHECK2(status, "poly_arch");

    p->cm_a = NULL;             //Default: deactivated.
    //If mobility is specified write it out.
    if (H5Lexists(file_id, "/parameter/cm_a", H5P_DEFAULT) > 0)
        {
            p->cm_a = (soma_scalar_t *) malloc(p->n_poly_type * sizeof(soma_scalar_t));
            if (p->cm_a == NULL)
                {
                    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
                    return -1;
                }
            status = read_hdf5(file_id, "/parameter/cm_a", H5T_SOMA_NATIVE_SCALAR, plist_id, p->cm_a);
            HDF5_ERROR_CHECK2(status, "/parameter/cm_a");
        }

    // read harmonic_normb_variable_scale
    if (H5Lexists(file_id, "/parameter/harmonic_normb_variable_scale", H5P_DEFAULT) > 0)
        {
            status =
                read_hdf5(file_id, "/parameter/harmonic_normb_variable_scale", H5T_SOMA_NATIVE_SCALAR, plist_id,
                          &(p->harmonic_normb_variable_scale));
            HDF5_ERROR_CHECK2(status, "/parameter/harmonic_normb_variable_scale");
        }
    else
        {
            if (p->info_MPI.sim_rank == 0)
                {
                    if (get_number_bond_type(p, HARMONICVARIABLESCALE) != 0)
                        fprintf(stderr,
                                "WARNING: The poly_arch contains HARMONICVARIABLESCALE Bond, but no corresponding value is set. Using 1 instead.\n");
                }
            p->harmonic_normb_variable_scale = 1.;
        }

    status = read_beads(p, file_id, plist_id);
    HDF5_ERROR_CHECK2(status, "total beads read");

    p->area51 = NULL;
    p->external_field_unified = NULL;
    p->umbrella_field = NULL;

    if (H5Lexists(file_id, "/area51", H5P_DEFAULT) > 0)
        {
            hid_t status = read_area51_hdf5(p, file_id, plist_id);
            if (status != 0)
                {
                    fprintf(stderr, "ERROR: failed to read area51 %s:%d.\n", __FILE__, __LINE__);
                    return status;
                }
        }

    p->serie_length = 0;

    if (H5Lexists(file_id, "/external_field", H5P_DEFAULT) > 0)
        {
            hid_t status = read_field_hdf5(p, file_id, plist_id, &(p->external_field_unified), "/external_field");
            if (status != 0)
                {
                    fprintf(stderr, "ERROR: failed to read external_field %s:%d.\n", __FILE__, __LINE__);
                    return status;
                }
            hid_t dataset = H5Dopen2(file_id, "/external_field", H5P_DEFAULT);
            status = dataset;
            HDF5_ERROR_CHECK(status);

            htri_t time_exists = H5Aexists(dataset, "serie_length");
            if (time_exists <= 0)
                {
                    p->serie_length = 1;
                    p->cos_serie = (soma_scalar_t *) malloc(p->serie_length * sizeof(soma_scalar_t));
                    p->sin_serie = (soma_scalar_t *) malloc(p->serie_length * sizeof(soma_scalar_t));
                    p->sin_serie[0] = 0;
                    p->cos_serie[0] = 1;
                    p->period = 1;
                }
            if (time_exists > 0)
                {               //time attr exists

                    hid_t length_attr =
                        H5Aopen_by_name(dataset, "/external_field", "serie_length", H5P_DEFAULT, H5P_DEFAULT);
                    hid_t d_space = H5Aget_space(length_attr);
                    HDF5_ERROR_CHECK(d_space);
                    hsize_t dims[1];    //ndims
                    status = H5Sget_simple_extent_dims(d_space, dims, NULL);
                    HDF5_ERROR_CHECK(status);
                    status = H5Aread(length_attr, H5T_NATIVE_UINT, &(p->serie_length));
                    HDF5_ERROR_CHECK(status);
                    status = H5Aclose(length_attr);
                    HDF5_ERROR_CHECK(status);

                    hid_t period_attr = H5Aopen_by_name(dataset, "/external_field", "period", H5P_DEFAULT, H5P_DEFAULT);
                    d_space = H5Aget_space(period_attr);
                    HDF5_ERROR_CHECK(d_space);
                    status = H5Sget_simple_extent_dims(d_space, dims, NULL);
                    HDF5_ERROR_CHECK(status);
                    status = H5Aread(period_attr, H5T_SOMA_NATIVE_SCALAR, &(p->period));
                    HDF5_ERROR_CHECK(status);
                    status = H5Aclose(period_attr);
                    HDF5_ERROR_CHECK(status);
                    p->cos_serie = (soma_scalar_t *) malloc(p->serie_length * sizeof(soma_scalar_t));
                    if (p->cos_serie == NULL)
                        {
                            fprintf(stderr, "Malloc error: %s:%d .\n", __FILE__, __LINE__);
                            return -1;
                        }

                    hid_t cos_attr = H5Aopen_by_name(dataset, "/external_field", "cos", H5P_DEFAULT, H5P_DEFAULT);
                    d_space = H5Aget_space(cos_attr);
                    HDF5_ERROR_CHECK(d_space);
                    dims[0] = p->serie_length;  //ndims
                    status = H5Sget_simple_extent_dims(d_space, dims, NULL);
                    HDF5_ERROR_CHECK(status);
                    status = H5Aread(cos_attr, H5T_SOMA_NATIVE_SCALAR, p->cos_serie);
                    HDF5_ERROR_CHECK(status);
                    status = H5Aclose(cos_attr);
                    HDF5_ERROR_CHECK(status);

                    p->sin_serie = (soma_scalar_t *) malloc(p->serie_length * sizeof(soma_scalar_t));
                    if (p->sin_serie == NULL)
                        {
                            fprintf(stderr, "Malloc error: %s:%d .\n", __FILE__, __LINE__);
                            return -1;
                        }

                    hid_t sin_attr = H5Aopen_by_name(dataset, "/external_field", "sin", H5P_DEFAULT, H5P_DEFAULT);
                    d_space = H5Aget_space(sin_attr);
                    HDF5_ERROR_CHECK(d_space);
                    dims[0] = p->serie_length;  //ndims
                    status = H5Sget_simple_extent_dims(d_space, dims, NULL);
                    HDF5_ERROR_CHECK(status);
                    status = H5Aread(sin_attr, H5T_SOMA_NATIVE_SCALAR, p->sin_serie);
                    HDF5_ERROR_CHECK(status);
                    status = H5Aclose(sin_attr);
                    HDF5_ERROR_CHECK(status);
                }
            status = H5Dclose(dataset);
            HDF5_ERROR_CHECK(status);
        }

    if (H5Lexists(file_id, "/umbrella_field", H5P_DEFAULT) > 0)
        {
            hid_t status = read_field_hdf5(p, file_id, plist_id, &(p->umbrella_field), "/umbrella_field");
            if (status != 0)
                {
                    fprintf(stderr, "ERROR: failed to read umbrella_field %s:%d.\n", __FILE__, __LINE__);
                    return status;
                }
        }

    status = read_poly_conversion_hdf5(p, file_id, plist_id);
    if (status != 0)
        {
            fprintf(stderr, "ERROR: %s:%d unable to read polytype conversion information.\n", __FILE__, __LINE__);
            return status;
        }

    status = read_mono_conversion_hdf5(p, file_id, plist_id);
    if (status != 0)
        {
            fprintf(stderr, "ERROR: %s:%d unable to read monotype conversion information.\n", __FILE__, __LINE__);
            return status;
        }

    status = read_mobility_hdf5(p, file_id, plist_id);
    if (status != 0)
        {
            fprintf(stderr, "ERROR: %s:%d unable to read mobility information.\n", __FILE__, __LINE__);
            return status;
        }

    p->field_scaling_type = (soma_scalar_t *) malloc(p->n_types * sizeof(soma_scalar_t));
    if (p->field_scaling_type == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    if (H5Lexists(file_id, "/parameter/density_weights", H5P_DEFAULT) > 0)
        {
            hid_t status = read_hdf5(file_id, "/parameter/density_weights", H5T_SOMA_NATIVE_SCALAR, plist_id,
                                     p->field_scaling_type);
            HDF5_ERROR_CHECK2(status, "/parameter/density_weights");
        }
    else
        for (unsigned int i = 0; i < p->n_types; i++)
            p->field_scaling_type[i] = 1;

    if ((status = H5Fclose(file_id)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, (int)status);
            return status;
        }

    return status;
}

int screen_output(struct Phase *const p, const unsigned int Nsteps)
{
    static time_t last_print = 0;
    static unsigned int last_time = 0;
    static double last_sec = 0;
    if (last_time == 0)
        last_time = p->start_time;
    if (last_sec == 0)
        {
            struct timeval last_tv;
            gettimeofday(&last_tv, NULL);
            last_sec = last_tv.tv_sec + last_tv.tv_usec * 1e-6;
        }

    struct timeval now;
    gettimeofday(&now, NULL);

    if (last_print == 0)
        last_print = now.tv_sec;
    const double second = now.tv_sec - p->start_clock.tv_sec + (now.tv_usec - p->start_clock.tv_usec) * 1e-6;
    const unsigned int steps_done = p->time - p->start_time;
    const time_t end = p->start_clock.tv_sec + second * (Nsteps) / (soma_scalar_t) steps_done;

    const double now_sec = now.tv_sec + now.tv_usec * 1e-6;

    const double sec = now_sec - last_sec;
    p->tps_elapsed_time += sec;
    p->tps_elapsed_steps += p->time - last_time;

    if (p->args.screen_output_interval_arg > 0 && now.tv_sec - last_print >= p->args.screen_output_interval_arg)
        {
            const double tps = p->tps_elapsed_steps / p->tps_elapsed_time;
            p->tps_elapsed_time = 1. / tps;
            p->tps_elapsed_steps = 1;

            if (p->info_MPI.sim_rank == 0)
                {
                    fprintf(stdout, "Rank %i:Running for %g [s] | TPS %g | steps-to-go: %u | ETA: %s",
                            p->info_MPI.world_rank, second, tps, Nsteps - steps_done, ctime(&end));
                    fflush(stdout);
                }
            last_print = now.tv_sec;
            last_time = p->time;
            last_sec = now_sec;
        }
    return 0;
}
