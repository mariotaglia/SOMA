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
#include "ana_server.h"
#include "phase.h" // for global constants. todo: move global constants somewhere else
#include <stdlib.h>
#include <assert.h>
#include "err_handling.h"
#include "string.h"
#include "server.h"
#include "serv_util.h"

int calc_MSD_server(const struct global_consts * gc, const Polymer * polymers,
    size_t n_polymers, soma_scalar_t * result, const struct comm_with_info * comm, int target)
{
    const int result_length = 8*gc->n_poly_type;
    const int counter_length = 2*gc->n_poly_type;
    uint64_t *const counter = (uint64_t *) calloc(counter_length , sizeof(uint64_t));
    RET_ERR_ON_NULL(counter, "CALLOC ERROR");

    memset(result, 0, 8 * gc->n_poly_type * sizeof(soma_scalar_t));
    // Add up local displacement
    for (uint64_t j = 0; j < n_polymers; j++)
        {                       /*Loop over polymers */

            soma_scalar_t tx_c = 0.;
            soma_scalar_t ty_c = 0.;
            soma_scalar_t tz_c = 0.;
            soma_scalar_t mx_c = 0.;
            soma_scalar_t my_c = 0.;
            soma_scalar_t mz_c = 0.;
            const unsigned int type = polymers[j].type;
            const unsigned int N = gc->poly_arch[gc->poly_type_offset[type]];
            for (unsigned int k = 0; k < N; k++)
                {               /*Loop over monomers */
                    const soma_scalar_t tx = polymers[j].msd_beads[k].x - polymers[j].beads[k].x;
                    const soma_scalar_t ty = polymers[j].msd_beads[k].y - polymers[j].beads[k].y;
                    const soma_scalar_t tz = polymers[j].msd_beads[k].z - polymers[j].beads[k].z;

                    // Add up values for chain diffusion
                    tx_c += polymers[j].beads[k].x;
                    ty_c += polymers[j].beads[k].y;
                    tz_c += polymers[j].beads[k].z;
                    mx_c += polymers[j].msd_beads[k].x;
                    my_c += polymers[j].msd_beads[k].y;
                    mz_c += polymers[j].msd_beads[k].z;

                    counter[type * 2 + 0] += 1;
                    result[type * 8 + 0] += tx * tx;
                    result[type * 8 + 1] += ty * ty;
                    result[type * 8 + 2] += tz * tz;
                    result[type * 8 + 3] += (tx * tx + ty * ty + tz * tz);
                }
            counter[type * 2 + 1] += 1;
            result[type * 8 + 4] += (tx_c - mx_c) * (tx_c - mx_c) / (N * N);
            result[type * 8 + 5] += (ty_c - my_c) * (ty_c - my_c) / (N * N);
            result[type * 8 + 6] += (tz_c - mz_c) * (tz_c - mz_c) / (N * N);
            result[type * 8 + 7] +=
                ((tx_c - mx_c) * (tx_c - mx_c) + (ty_c - my_c) * (ty_c - my_c) +
                 (tz_c - mz_c) * (tz_c - mz_c)) / (N * N);
        }

    DPRINT("NOW REDUCING")
    if (comm->rank == target)
        {
            MPI_Reduce(MPI_IN_PLACE, result, result_length, MPI_SOMA_SCALAR,
                MPI_SUM, target, comm->comm);
            MPI_Reduce(MPI_IN_PLACE, counter, counter_length, MPI_UINT64_T,
                MPI_SUM, target, comm->comm);
            //Looping over twice the number of poly types. But loop over half the elements
            // 8/2 = 4, because the norm for first and second half differ.
            for (unsigned int type = 0; type < 2 * gc->n_poly_type; type++)
                for (unsigned int i = 0; i < 4; i++)
                    if (counter[type] > 0)
                        result[type * 4 + i] /= (soma_scalar_t) counter[type];

        }
    else
        {
            soma_scalar_t * result_tmp = calloc(result_length, sizeof(soma_scalar_t));
            RET_ERR_ON_NULL(result_tmp, "Calloc")
            uint64_t * counter_tmp = calloc(counter_length, sizeof(uint64_t));
            RET_ERR_ON_NULL(result_tmp, "Calloc")

            MPI_Reduce(result, result_tmp, result_length, MPI_SOMA_SCALAR,
                MPI_SUM, target, comm->comm);
            MPI_Reduce(counter, counter_tmp, counter_length, MPI_UINT64_T,
                MPI_SUM, target, comm->comm);

            free(result_tmp);
            free(counter_tmp);
        }
    DPRINT("REDUCTION COMPLETE")

    free(counter);

    return 0;
}


int analytics_server(const struct global_consts *gc,
    const Ana_Info * ai, const struct server_info * si, const struct receiver *rcv,
        unsigned int time)
{

    Polymer * polymers = rcv->polymers;
    size_t n_polymers = rcv->n_polymers;


    if (need_to_do(ai->delta_mc_MSD, time))
        {
            DPRINT("server calculating msd")
            int target = 0;
            int msd_len = 8*gc->n_poly_type;
            soma_scalar_t * msd = malloc(msd_len*sizeof(soma_scalar_t));
            int err = calc_MSD_server(gc, polymers, n_polymers,
                msd, &si->server_comm, target);
            if (err)
                {
                    DPRINT("calculating MSD failed")
                    return err;
                }
            if (si->server_comm.rank == target)
                {
                    assert(ai->file_id != -1);
                    extent_ana_by_field(msd, msd_len, "/MSD",ai->file_id);
                }
            free(msd);
            herr_t status = H5Fflush(ai->file_id, H5F_SCOPE_LOCAL);
            if (status < 0)
                {
                    fprintf(stderr, "ERROR: HDF5 %s:%s:%d\n", __func__, __FILE__, __LINE__);
                    return -1;
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
