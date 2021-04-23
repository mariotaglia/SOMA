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

#include "mobility.h"

#include <stdlib.h>
#include <hdf5.h>
#include "phase.h"
#include "io.h"
#include "mesh.h"

//! \file mobility.c
//! \brief Implementation of mobility.h

int read_mobility_hdf5(struct Phase *const p, const hid_t file_id, const hid_t plist_id)
{
    p->mobility.type = DEFAULT_MOBILITY;
    p->mobility.param_len = 0;
    p->mobility.param = NULL;

    //Default polytype mobility
    p->mobility.poly_type_mc_freq = (unsigned int *)malloc(p->n_poly_type * sizeof(unsigned int));
    MALLOC_ERROR_CHECK(p->mobility.poly_type_mc_freq, p->n_poly_type * sizeof(unsigned int));
    for (unsigned int i = 0; i < p->n_poly_type; i++)
        p->mobility.poly_type_mc_freq[i] = 1;

    //Quick exit if no mobility is present in the file
    if (!(H5Lexists(file_id, "/mobility", H5P_DEFAULT) > 0))
        return 0;

    herr_t status;
    unsigned int tmp_type;
    status = read_hdf5(file_id, "mobility/type", H5T_STD_U32LE, plist_id, &tmp_type);
    HDF5_ERROR_CHECK(status);

    enum MobilityEnum tmp_enum_type = tmp_type;

    //Read the poly_type_mc_freq in any case
    status = read_hdf5(file_id, "mobility/poly_type_mc_freq", H5T_STD_U32LE, plist_id, p->mobility.poly_type_mc_freq);
    HDF5_ERROR_CHECK(status);
    for (unsigned int i = 0; i < p->n_poly_type; i++)
        if (p->mobility.poly_type_mc_freq[i] == 0)
            {
                fprintf(stderr, "ERROR: %s:%d poly_type_mc_freq contains a 0. This is illdefined. min value is 1.\n",
                        __FILE__, __LINE__);
                return -2;
            }

    //Quick exit with default unmodified mobility
    if (tmp_enum_type == DEFAULT_MOBILITY)
        return 0;

    //Read the conversion list before reading the spatial array.
    const hid_t dset_param = H5Dopen(file_id, "/mobility/param", H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset_param);
    const hid_t dspace_param = H5Dget_space(dset_param);
    HDF5_ERROR_CHECK(dspace_param);
    const unsigned int ndims_param = H5Sget_simple_extent_ndims(dspace_param);
    if (ndims_param != 1)
        {
            fprintf(stderr, "ERROR: %s:%d not the correct number of dimensions for %s.\n",
                    __FILE__, __LINE__, "mobility/param");
            return -1;
        }
    hsize_t dim_param;
    status = H5Sget_simple_extent_dims(dspace_param, &dim_param, NULL);
    HDF5_ERROR_CHECK(status);

    switch (tmp_enum_type)      //Verify valid input types
        {
        case DEFAULT_MOBILITY:
            break;
        case MULLER_SMITH_MOBILITY:
            if (dim_param != 2 * p->n_types)
                {
                    fprintf(stderr, "ERROR: %s:%d invalid number of mobility parameter for MULLER_SMITH\n", __FILE__,
                            __LINE__);
                    fprintf(stderr, "\tExpected number of parameters %d but got %d parameters.\n", 2 * p->n_types,
                            (int)dim_param);
                    return -2;
                }
            break;
        case TANH_MOBILITY:
            if (dim_param != 2 * p->n_types + p->n_types * p->n_types)
                {
                    fprintf(stderr, "ERROR: %s:%d invalid number of mobility parameter for TANH\n", __FILE__, __LINE__);
                    fprintf(stderr, "\tExpected number of parameters %d but got %d parameters.\n", 2 * p->n_types,
                            (int)dim_param);
                    return -3;
                }
            break;
        default:
            fprintf(stderr, "ERROR: reading invalid mobility type. Old SOMA version? %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    p->mobility.param = (soma_scalar_t *) malloc(dim_param * sizeof(soma_scalar_t));
    MALLOC_ERROR_CHECK(p->mobility.param, dim_param * sizeof(soma_scalar_t));

    p->mobility.param_len = dim_param;

    //Actuablly read the polyconversion list data
    status = H5Dread(dset_param, H5T_SOMA_FILE_SCALAR, H5S_ALL, H5S_ALL, plist_id, p->mobility.param);
    HDF5_ERROR_CHECK(status);
    status = H5Sclose(dspace_param);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset_param);

    //Enable the updat only if everything worked fine
    p->mobility.type = tmp_enum_type;
    return 0;
}

int write_mobility_hdf5(const struct Phase *const p, const hid_t file_id, const hid_t plist_id)
{
    //Quick exit for default mobility
    if (p->mobility.type == DEFAULT_MOBILITY)
        return 0;

    hid_t group = H5Gcreate2(file_id, "/mobility", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    herr_t status;
    const hsize_t one = 1;
    unsigned int tmp_type = p->mobility.type;
    status = write_hdf5(1, &one, file_id, "/mobility/type", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id, &(tmp_type));
    HDF5_ERROR_CHECK(status);

    const hsize_t n_poly_type = p->n_poly_type;
    status =
        write_hdf5(1, &n_poly_type, file_id, "mobility/poly_type_mc_freq", H5T_STD_U32LE, H5T_NATIVE_UINT, plist_id,
                   p->mobility.poly_type_mc_freq);

    const hsize_t list_len = p->mobility.param_len;
    status =
        write_hdf5(1, &list_len, file_id, "/mobility/param", H5T_SOMA_FILE_SCALAR, H5T_SOMA_NATIVE_SCALAR, plist_id,
                   p->mobility.param);
    HDF5_ERROR_CHECK(status);

    if ((status = H5Gclose(group)) < 0)
        {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
                    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
        }

    return 0;
}

int copyin_mobility(struct Phase *p)
{
#ifdef _OPENACC
#pragma acc enter data copyin(p->mobility.poly_type_mc_freq[0:p->n_poly_type])
    if (p->mobility.type != DEFAULT_MOBILITY)
        {
#pragma acc enter data copyin(p->mobility.param[0:p->mobility.param_len])
        }
#endif                          //_OPENACC
    return p->n_poly_type * 0;
}

int copyout_mobility(struct Phase *p)
{
#ifdef _OPENACC
#pragma acc exit data copyout(p->mobility.poly_type_mc_freq[0:p->n_poly_type])
    if (p->mobility.type != DEFAULT_MOBILITY)
        {
#pragma acc exit data copyout(p->mobility.param[0:p->mobility.param_len])
        }
#endif                          //_OPENACC
    return p->n_poly_type * 0;
}

int update_self_mobility(const struct Phase *const p)
{
#ifdef _OPENACC
#pragma acc update self(p->mobility.poly_type_mc_freq[0:p->n_poly_type])
    if (p->mobility.type != DEFAULT_MOBILITY)
        {
#pragma acc update self(p->mobility.param[0:p->mobility.param_len])
        }
#endif                          //_OPENACC
    return p->n_poly_type * 0;
}

int free_mobility(struct Phase *p)
{
    free(p->mobility.poly_type_mc_freq);
    free(p->mobility.param);
    return 0;
}

soma_scalar_t get_mobility_modifier(const struct Phase *const p, const unsigned int particle_type,
                                    const soma_scalar_t x, const soma_scalar_t y, const soma_scalar_t z)
{
    switch (p->mobility.type)
        {
        default:
            // intentional fall through
        case DEFAULT_MOBILITY:
            return 1;
            break;
        case MULLER_SMITH_MOBILITY:
            ;
            soma_scalar_t modifier = 0;
            for (unsigned int type = 0; type < p->n_types; type++)
                {
                    const uint64_t cellindex = coord_to_index_unified(p, x, y, z, type);
                    const soma_scalar_t phi = p->field_scaling_type[type] * p->fields_unified[cellindex];
                    modifier += p->mobility.param[type] * phi;
                    modifier += p->mobility.param[p->n_types + type] * phi * phi;
                }

            if (modifier == (soma_scalar_t) 0.)
                return 1;
            modifier = ((soma_scalar_t) 1.) / modifier;

            // Fulfill garantied modifier bounds
            if (modifier > 1)
                return 1;
            if (modifier < 0)
                return 1;

            return modifier;
            break;
        case TANH_MOBILITY:
            ;
            const soma_scalar_t *const phi_0 = p->mobility.param;
            const soma_scalar_t *const delta_phi = phi_0 + p->n_types;
            const soma_scalar_t *const a = delta_phi + p->n_types;

            soma_scalar_t a_sum = 0;
            for (unsigned int j = 0; j < p->n_types; j++)
                {
                    const uint64_t cellindex = coord_to_index_unified(p, x, y, z, j);
                    a_sum +=
                        a[particle_type * p->n_types + j] * p->field_scaling_type[j] * p->fields_unified[cellindex];
                }

            return 0.5 * (1 + tanh((phi_0[particle_type] - a_sum) / delta_phi[particle_type]));
            break;
        }

    return 1;
}
