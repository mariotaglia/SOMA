/* Copyright (C) 2020-2021 Ludwig Schneider

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

//! \file self_documentation.c
//! \brief Implementation self_documentation.h

#include "self_documentation.h"

#include <stdlib.h>
#include <stdio.h>
#include <hdf5.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <unistd.h>

#include "soma_util.h"
#include "phase.h"

//! Write the current documentation string to a file
//!
//! \param ftmp File handle of the temporary file to write in
//! \param p Phase struct to describe
//! \return Errorcode
int generate_current_documentation_string(FILE * ftmp, struct Phase *p)
{
    time_t now = time(NULL);
    fprintf(ftmp, "Self documentation of SOMA simulation started at %s", ctime(&now));
#if __GLIBC__ > 2 || __GLIBC_MINOR__ > 24
#include <sys/utsname.h>
    struct utsname buffer;
    if (uname(&buffer) < 0)
        fprintf(stderr, "ERROR: accessing uname\n");
    else
        fprintf(ftmp, "Executing machine: %s with %s on %s.\n", buffer.nodename, buffer.sysname, buffer.machine);
#endif                          //__GLIBC

    fprintf(ftmp, "Version info:\n");
    fprintf(ftmp, "\tCompiled on %s with C std %ld.\n", get_soma_system_info(), __STDC_VERSION__);
    fprintf(ftmp, "\tGIT version of SOMA is %s compiled on %s %s.\n", get_soma_version(), __DATE__, __TIME__);

    //HDF5
    unsigned int majnum = 0, minnum = 0, relnum = 0;
    H5get_libversion(&majnum, &minnum, &relnum);
    fprintf(ftmp, "\tHDF5 version is %u.%u.%u\n", majnum, minnum, relnum);
#if ( ENABLE_MPI == 1 )
#ifdef MPI_MAX_LIBRARY_VERSION_STRING
    //MPI
    char mpi_version[MPI_MAX_LIBRARY_VERSION_STRING];
    int length;
    MPI_Get_library_version(mpi_version, &length);
    int mpi_maj = 0, mpi_min = 0;
    MPI_Get_version(&mpi_maj, &mpi_min);
    fprintf(ftmp, "\tMPI version: %s %d.%d\n", mpi_version, mpi_maj, mpi_min);
#else
    fprintf(ftmp, "\tNo MPI lib version available.\n");
#endif                          //mpi_max_library_version_string
#endif                          //( ENABLE_MPI == 1 )

    fprintf(ftmp, "Arguments used to start SOMA.\n");
    cmdline_parser_dump(ftmp, &(p->args));

    fprintf(ftmp, "Selected parameter of the simulation.\n");
    fprintf(ftmp, "\tMC step at beginning %d\n", p->time);
    fprintf(ftmp, "\tPRNG seed %u\n", p->args.rng_seed_arg);
    fprintf(ftmp, "\tHamiltionian used %d\n", p->hamiltonian);
    fprintf(ftmp, "\tInteraction parameter (kN and xN):");
    for (unsigned int i = 0; i < p->n_types; i++)
        for (unsigned int j = i; j < p->n_types; j++)
            fprintf(ftmp, " %dx%d=%f", i, j, p->xn[i * p->n_types + j]);
    fprintf(ftmp, "\n");
    if (p->harmonic_normb_variable_scale != 1.)
        fprintf(ftmp, "\tHarmonic variable scale = %f\n", p->harmonic_normb_variable_scale);
    fprintf(ftmp, "\tDiffusion constants for types:");
    for (unsigned int i = 0; i < p->n_types; i++)
        fprintf(ftmp, " %f", p->A[i]);
    fprintf(ftmp, "\n");
    if (p->cm_a)
        {
            fprintf(ftmp, "\tCenter of mass diffusion constants for molecules:");
            for (unsigned int i = 0; i < p->n_poly_type; i++)
                fprintf(ftmp, " %f", p->cm_a[i]);
            fprintf(ftmp, "\n");
        }
    fprintf(ftmp, "\tField scaling type:");
    for (unsigned int i = 0; i < p->n_types; i++)
        fprintf(ftmp, " %f", p->field_scaling_type[i]);
    fprintf(ftmp, "\n");
    if (p->umbrella_field)
        {
            fprintf(ftmp, "\tStrength of umbrella:");
            for (unsigned int i = 0; i < p->n_types; i++)
                fprintf(ftmp, " %f", p->k_umbrella[i]);
            fprintf(ftmp, "\n");
        }
    fprintf(ftmp, "\tBox dimensions: %f (%d)\t%f (%d)\t%f (%d)\n", p->Lx, p->nx, p->Ly, p->ny, p->Lz, p->nz);
    if (p->bead_data_read)
        fprintf(ftmp, "\tbeads data read: yes\n");
    if (p->area51)
        fprintf(ftmp, "\tArea51 present: yes\n");
    if (p->external_field_unified)
        fprintf(ftmp, "\tExternal field present: yes\n");
    if (p->umbrella_field)
        fprintf(ftmp, "\tUmbrella field present: yes\n");
    if (p->pc.deltaMC)
        fprintf(ftmp, "\tPolytype conversion active %d,\t", p->pc.deltaMC);
    if (p->pc.rate)
        {
            for (unsigned int conv = 0; conv < p->pc.len_reactions; conv++)
                {
                    fprintf(ftmp, "%d --> %d with rate %f", p->pc.input_type[conv],
                            p->pc.output_type[conv], p->pc.rate[conv]);
                    for (unsigned int dd = 0; dd < p->pc.dependency_ntype[conv]; dd++)
                        {
                            fprintf(ftmp, " * phi_%d", p->pc.dependency_type[p->pc.dependency_type_offset[conv] + dd]);
                        }
                    fprintf(ftmp, "\t");

                }
            fprintf(ftmp, "\n");
        }
#if ( ENABLE_MONOTYPE_CONVERSIONS == 1 )
    if (p->mtc.deltaMC)
        fprintf(ftmp, "\tMonotype conversion active %d,\t", p->mtc.deltaMC);
    if (p->mtc.rate)
        {
            for (unsigned int conv = 0; conv < p->mtc.len_reactions; conv++)
                {
                    fprintf(ftmp, "%d --> %d with rate %f", p->mtc.input_type[conv],
                            p->mtc.output_type[conv], p->mtc.rate[conv]);
                    for (unsigned int dd = 0; dd < p->mtc.dependency_ntype[conv]; dd++)
                        {
                            fprintf(ftmp, " * phi_%d",
                                    p->mtc.dependency_type[p->mtc.dependency_type_offset[conv] + dd]);
                        }
                    fprintf(ftmp, "\t");
                }
        }
    if (p->mtc.deltaMC)
        fprintf(ftmp, " with block size %d\n", p->mtc.block_size);
#endif                          //ENABLE_MONOTYPE_CONVERSION
    if (p->mobility.type != DEFAULT_MOBILITY)
        fprintf(ftmp, "\tMobility modification active %d\n", p->mobility.type);
    if (p->serie_length > 1 || (p->serie_length == 1 && p->sin_serie[0] != 0))
        fprintf(ftmp, "\tTime dependent external field active: yes\n");

    return 0;
}

int init_self_documentation(struct Phase *p, char *filename, struct SelfDocumentation *sd)
{
    sd->Ndoc = 0;
    sd->simdoc = NULL;
    sd->simdoc_internal = NULL;

    unsigned int total_size = 0;

    // The documentation string is processed only on simrank 0, and broadcasted to all other ranks
    if (p->info_MPI.sim_rank == 0)
        {
            // A temporary file is created to generate the documentation string.
            // The cmdline parser does not offer a version to print the arguments into a char array, only files.
            // Since, this makes the file required anyway, we use it to generate the documentation string of unknown length conveniently.
            FILE *ftmp = tmpfile();
            if (ftmp == NULL)
                {
                    fprintf(stderr, "ERROR: %s:%d cannot open tmpfile.\n", __FILE__, __LINE__);
                    unsigned int tmp_buffer = 0;
#if (ENABLE_MPI == 1)
                    //Broadcast failure to all MPI ranks to signify failure to avoid hanging.
                    MPI_Bcast(&tmp_buffer, 1, MPI_UNSIGNED, 0, p->info_MPI.SOMA_comm_sim);
#endif                          //ENABLE_MPI
                    return -1;
                }
            generate_current_documentation_string(ftmp, p);

            // obtain file size:
            rewind(ftmp);
            fseek(ftmp, 0, SEEK_END);
            const unsigned int current_size = ftell(ftmp);
            total_size = current_size + 1;
            rewind(ftmp);
            unsigned int size_counter = 0;

            if (filename != NULL)
                {
                    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
                    if (H5Pset_fclose_degree(plist_id, H5F_CLOSE_STRONG) < 0)
                        {
                            fprintf(stderr, "ERROR: setting file access properties %s:%d\n", __FILE__, __LINE__);
                            return -1;
                        }

                    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, plist_id);
                    herr_t status = 0;
                    if (H5Lexists(file_id, "documentation", H5P_DEFAULT))
                        {
                            hid_t dset = H5Dopen(file_id, "documentation", H5P_DEFAULT);
                            HDF5_ERROR_CHECK2(dset, "documentation");
                            hid_t d_space = H5Dget_space(dset);
                            HDF5_ERROR_CHECK2(d_space, "documentation");
#ifndef NDEBUG
                            const hsize_t ndims = H5Sget_simple_extent_ndims(d_space);
                            assert(ndims == 1);
#endif                          //NDEBUG
                            hsize_t dims[1];    //ndims
                            status = H5Sget_simple_extent_dims(d_space, dims, NULL);
                            HDF5_ERROR_CHECK2(status, "documentation");

                            hid_t tid1 = H5Tcopy(H5T_C_S1);
                            status = H5Tset_size(tid1, H5T_VARIABLE);
                            HDF5_ERROR_CHECK2(status, "documentation");

                            char **rdata = (char **)malloc(dims[0] * sizeof(char *));
                            MALLOC_ERROR_CHECK(rdata, dims[0] * sizeof(char *));

                            status = H5Dread(dset, tid1, H5S_ALL, H5S_ALL, H5P_DEFAULT, rdata);
                            HDF5_ERROR_CHECK2(status, "documentation");

                            for (unsigned int i = 0; i < dims[0]; i++)
                                {
                                    total_size += strlen(rdata[i]) + 1;
                                    sd->Ndoc += 1;
                                }
                            sd->simdoc_internal = (char *)malloc(total_size * sizeof(char));
                            MALLOC_ERROR_CHECK(sd->simdoc_internal, total_size * sizeof(char));
                            memset(sd->simdoc_internal, '\0', total_size * sizeof(char));

                            //Copy all old documentations in
                            for (unsigned int i = 0; i < dims[0]; i++)
                                {
                                    memcpy(sd->simdoc_internal + size_counter, rdata[i], strlen(rdata[i]));
                                    size_counter += strlen(rdata[i]) + 1;
                                }

                            status = H5Dvlen_reclaim(tid1, d_space, H5P_DEFAULT, rdata);
                            HDF5_ERROR_CHECK2(status, "documentation");

                            free(rdata);
                            H5Tclose(tid1);
                            H5Sclose(d_space);
                            H5Dclose(dset);
                        }
                    H5Fclose(file_id);
                }

            if (size_counter == 0)      //No previous history. Need to allocate the char buffer
                {
                    assert(sd->simdoc_internal == NULL);
                    assert(total_size == current_size + 1);
                    sd->simdoc_internal = (char *)malloc(total_size * sizeof(char));
                    MALLOC_ERROR_CHECK(sd->simdoc_internal, total_size * sizeof(char));
                    memset(sd->simdoc_internal, '\0', total_size * sizeof(char));
                }
            //copy current doc from tmp file
            const unsigned int result = fread(sd->simdoc_internal + size_counter, 1, current_size, ftmp);
            if (result != current_size)
                {
                    fprintf(stderr, "ERROR: %s:%d reading from tmpfile", __FILE__, __LINE__);
                    unsigned int tmp_buffer = 0;
#if (ENABLE_MPI == 1)
                    //Broadcast failure to all MPI ranks to signify failure to avoid hanging.
                    MPI_Bcast(&tmp_buffer, 1, MPI_UNSIGNED, 0, p->info_MPI.SOMA_comm_sim);
#endif                          //ENABLE_MPI
                    return -1;
                }
            fclose(ftmp);
            sd->Ndoc += 1;
        }
#if (ENABLE_MPI == 1)
    //Broadcast final info to all ranks
    // start with the length of the buffer
    MPI_Bcast(&total_size, 1, MPI_UNSIGNED, 0, p->info_MPI.SOMA_comm_sim);
    if (total_size > 0)         //Failures would have send 0
        {
            //Cast the number of doc strings
            MPI_Bcast(&(sd->Ndoc), 1, MPI_UNSIGNED, 0, p->info_MPI.SOMA_comm_sim);
            if (p->info_MPI.sim_rank != 0)
                {
                    sd->simdoc_internal = (char *)malloc(total_size * sizeof(char));
                    MALLOC_ERROR_CHECK(sd->simdoc_internal, total_size * sizeof(char));
                }
            MPI_Bcast(sd->simdoc_internal, total_size, MPI_CHAR, 0, p->info_MPI.SOMA_comm_sim);
        }
#endif                          //ENABLE_MPI
    sd->simdoc = (char **)malloc(sd->Ndoc * sizeof(char *));
    MALLOC_ERROR_CHECK(sd->simdoc, sd->Ndoc * sizeof(char *));
    sd->simdoc[0] = sd->simdoc_internal;
    for (unsigned int i = 1; i < sd->Ndoc; i++)
        sd->simdoc[i] = sd->simdoc[i - 1] + strlen(sd->simdoc[i - 1]) + 1;

    return 0;
}

int free_self_documentation(SelfDocumentation * sd)
{
    if (sd->Ndoc > 0)
        {
            free(sd->simdoc_internal);
            free(sd->simdoc);
        }
    return 0;
}

int print_self_documentation(SelfDocumentation * sd, FILE * f)
{
    if (sd->Ndoc > 0)
        {
            //Print the current documentation step
            fprintf(f, "%s", sd->simdoc[sd->Ndoc - 1]);
            return 0;
        }
    return -1;
}

int add_self_documentation_to_hdf5(const SelfDocumentation * sd, const hid_t file_id, const hid_t plist_id)
{
    if (sd->Ndoc > 0)
        {
            herr_t status = 0;
            hsize_t size_sd[1] = { sd->Ndoc };
            hid_t dataspace = H5Screate_simple(1, size_sd, NULL);
            HDF5_ERROR_CHECK2(dataspace, "documentation");

            hid_t tid1 = H5Tcopy(H5T_C_S1);
            H5Tset_size(tid1, H5T_VARIABLE);

            //Remove evtl. old dataset of documentation
            if (H5Lexists(file_id, "documentation", H5P_DEFAULT) > 0)
                {
                    status = H5Ldelete(file_id, "documentation", H5P_DEFAULT);
                    HDF5_ERROR_CHECK2(status, "documentation");
                }

            hid_t dataset =
                H5Dcreate2(file_id, "documentation", tid1, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
            HDF5_ERROR_CHECK2(dataset, "documentation");

            status |= H5Dwrite(dataset, tid1, H5S_ALL, H5S_ALL, plist_id, sd->simdoc);
            HDF5_ERROR_CHECK2(status, "documentation");
            status |= H5Sclose(dataspace);
            HDF5_ERROR_CHECK2(status, "documentation");
            status |= H5Dclose(dataset);
            HDF5_ERROR_CHECK2(status, "documentation");
            return status;
        }
    return 0;
}
