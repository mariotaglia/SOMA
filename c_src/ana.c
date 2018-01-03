/* Copyright (C) 2016-2018 Ludwig Schneider
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

//! \file ana.c
//! \brief Implementation of ana.h

#include "ana.h"
#include <time.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <stdlib.h>
#include "phase.h"
#include "polymer.h"
#include "mesh.h"
#include "io.h"

void calc_Re(const struct Phase * p, soma_scalar_t *const result)
    {
    uint64_t *const counter = (uint64_t*)malloc( p->n_poly_type*sizeof(uint64_t));
    if(counter == NULL)
	{
	fprintf(stderr,"MALLOC ERROR: %s:%d\n",__FILE__,__LINE__);
	return;
	}
    memset(counter,0, p->n_poly_type*sizeof(uint64_t));
    memset(result, 0 , 4*p->n_poly_type * sizeof(soma_scalar_t));

    for (uint64_t npoly = 0; npoly < p->n_polymers; npoly++) {
	const unsigned int type = p->polymers[npoly].type;

	const unsigned int start = p->end_mono[ type*2 + 0];
	const unsigned int end = p->end_mono[ type*2 + 1];
	const soma_scalar_t dx = p->polymers[npoly].beads[start].x - p->polymers[npoly].beads[end].x;
	const soma_scalar_t dy = p->polymers[npoly].beads[start].y - p->polymers[npoly].beads[end].y;
	const soma_scalar_t dz = p->polymers[npoly].beads[start].z - p->polymers[npoly].beads[end].z;

	result[type*4 + 0 ] += dx*dx + dy*dy + dz*dz;
	result[type*4 + 1 ] += dx*dx ;
	result[type*4 + 2 ] += dy*dy ;
	result[type*4 + 3 ] += dz*dz ;
	counter[type] += 1;
	}

    MPI_Allreduce(MPI_IN_PLACE, result, 4*p->n_poly_type, MPI_SOMA_SCALAR, MPI_SUM,
		  p->info_MPI.SOMA_MPI_Comm);
    MPI_Allreduce(MPI_IN_PLACE, counter, p->n_poly_type, MPI_UINT64_T, MPI_SUM,
		  p->info_MPI.SOMA_MPI_Comm);

    for(unsigned int type=0 ; type < p->n_poly_type; type++)
	for(unsigned int i=0; i < 4; i++)
	    if( counter[type] > 0)
		result[type*4+i] /= (soma_scalar_t) counter[type];

    free(counter);
    }

void calc_dvar(const struct Phase * p, soma_scalar_t *dvar)
    {
    soma_scalar_t var = 0.0;
    for (uint64_t index=0; index < p->n_cells*p->n_types; index++) {
    var += (p->fields_unified[index] - p->old_fields_unified[index])*(p->fields_unified[index] - p->old_fields_unified[index]);
    }
    memcpy(p->old_fields_unified, p->fields_unified, p->n_cells*p->n_types*sizeof(uint16_t));
    var = var/(p->n_cells*p->n_types);
    *dvar = var;
    }

void calc_Rg(const struct Phase *p, soma_scalar_t *const result)
    {
    uint64_t *const counter = (uint64_t*)malloc( p->n_poly_type*sizeof(uint64_t));
    if(counter == NULL)
	{
	fprintf(stderr,"MALLOC ERROR: %s:%d\n",__FILE__,__LINE__);
	return;
	}
    memset(counter,0, p->n_poly_type*sizeof(uint64_t));
    memset(result, 0 , 4*p->n_poly_type * sizeof(soma_scalar_t));

    for (uint64_t ipoly = 0; ipoly < p->n_polymers; ipoly++) {
	const unsigned int type =p->polymers[ipoly].type;
	const unsigned int N = p->poly_arch[ p->poly_type_offset[ type ] ];
	soma_scalar_t xcm = 0.;
	soma_scalar_t ycm = 0.;
	soma_scalar_t zcm = 0.;
	soma_scalar_t x2  = 0.;
	soma_scalar_t y2  = 0.;
	soma_scalar_t z2  = 0.;
	for (unsigned int ibead = 0; ibead < N; ibead++) {
	    const soma_scalar_t x1   = p->polymers[ipoly].beads[ibead].x;
	    const soma_scalar_t y1   = p->polymers[ipoly].beads[ibead].y;
	    const soma_scalar_t z1   = p->polymers[ipoly].beads[ibead].z;
	    xcm += x1;
	    ycm += y1;
	    zcm += z1;
	    x2  += x1*x1;
	    y2  += y1*y1;
	    z2  += z1*z1;
	    }
	xcm /= (soma_scalar_t)(N);
	ycm /= (soma_scalar_t)(N);
	zcm /= (soma_scalar_t)(N);

	result[type*4 + 0] += (x2/(soma_scalar_t)(N) - xcm*xcm) +
	    (y2/(soma_scalar_t)(N) - ycm*ycm) +
	    (z2/(soma_scalar_t)(N) - zcm*zcm);
	result[type*4 + 1] += (x2/(soma_scalar_t)(N) - xcm*xcm);
	result[type*4 + 2] += (y2/(soma_scalar_t)(N) - ycm*ycm);
	result[type*4 + 3] += (z2/(soma_scalar_t)(N) - zcm*zcm);
	counter[type] += 1;
	}

    MPI_Allreduce(MPI_IN_PLACE, result, 4*p->n_poly_type, MPI_SOMA_SCALAR, MPI_SUM,
		  p->info_MPI.SOMA_MPI_Comm);
    MPI_Allreduce(MPI_IN_PLACE, counter, p->n_poly_type, MPI_UINT64_T, MPI_SUM,
		  p->info_MPI.SOMA_MPI_Comm);

    for(unsigned int type=0 ; type < p->n_poly_type; type++)
	for(unsigned int i=0; i < 4; i++)
	    if( counter[type] > 0)
		result[type*4 +i] /= (soma_scalar_t) counter[type];
    free(counter);
    }

void calc_anisotropy(const struct Phase * p, soma_scalar_t *const result)
    {
    uint64_t *const counter = (uint64_t*)malloc( p->n_poly_type*sizeof(uint64_t));
    if(counter == NULL)
	{
	fprintf(stderr,"MALLOC ERROR: %s:%d\n",__FILE__,__LINE__);
	return;
	}
    memset(counter,0, p->n_poly_type*sizeof(uint64_t));
    memset(result, 0 , 6*p->n_poly_type * sizeof(soma_scalar_t));

    // loop over local chains
    for (uint64_t ipoly = 0; ipoly < p->n_polymers; ipoly++){
	const Polymer*const poly = &(p->polymers[ipoly]);
	const unsigned int type = poly->type;
	const unsigned int N = p->poly_arch[ p->poly_type_offset[ type ] ];

	// loop over beads in this chain
	for (unsigned int ibead = 0; ibead < N; ibead++) {

	    const soma_scalar_t x1 = poly->beads[ibead].x;
	    const soma_scalar_t y1 = poly->beads[ibead].y;
	    const soma_scalar_t z1 = poly->beads[ibead].z;

	    // loop over bonds of this bead
	    const int start = get_bondlist_offset(p->poly_arch[p->poly_type_offset[type] + ibead + 1]);
	    if(start > 0){
		int i = start;
		int end;
		do{
		    const uint32_t bn = p->poly_arch[i++];
		    const int info = bn;
		    end = get_end(info);
		    const int offset = get_offset(info);
		    const int neighbour_id = ibead + offset;
		    const unsigned int jbead = neighbour_id;

		    const soma_scalar_t x2 = p->polymers[ipoly].beads[jbead].x;
		    const soma_scalar_t y2 = p->polymers[ipoly].beads[jbead].y;
		    const soma_scalar_t z2 = p->polymers[ipoly].beads[jbead].z;

		    const soma_scalar_t bx = x2 - x1;
		    const soma_scalar_t by = y2 - y1;
		    const soma_scalar_t bz = z2 - z1;

		    result[type*6 + 0] += bx * bx;
		    result[type*6 + 1] += by * by;
		    result[type*6 + 2] += bz * bz;
		    result[type*6 + 3] += bx * by;
		    result[type*6 + 4] += by * bz;
		    result[type*6 + 5] += bz * bx;
		    counter[type] += 1;
		    }while( end == 0);
		}
	    }
	}

    MPI_Allreduce(MPI_IN_PLACE, result, 6*p->n_poly_type, MPI_SOMA_SCALAR, MPI_SUM,
		  p->info_MPI.SOMA_MPI_Comm);
    MPI_Allreduce(MPI_IN_PLACE, counter, p->n_poly_type, MPI_UINT64_T, MPI_SUM,
		  p->info_MPI.SOMA_MPI_Comm);

    for(unsigned int type=0 ; type < p->n_poly_type; type++)
	for(unsigned int i=0; i < 6; i++)
	    if( counter[type] > 0)
		result[type*6 +i] /= (soma_scalar_t) counter[type];
    free(counter);
    }

void calc_MSD(const struct Phase * p, soma_scalar_t *const result)
    {
    uint64_t *const counter = (uint64_t*)malloc( 2*p->n_poly_type*sizeof(uint64_t));
    if(counter == NULL)
	{
	fprintf(stderr,"MALLOC ERROR: %s:%d\n",__FILE__,__LINE__);
	return;
	}
    memset(counter,0, 2*p->n_poly_type*sizeof(uint64_t));
    memset(result, 0 , 8*p->n_poly_type * sizeof(soma_scalar_t));

    // Add up local displacement
    for (uint64_t j = 0; j < p->n_polymers; j++) {	/*Loop over polymers */

        soma_scalar_t tx_c = 0.;
        soma_scalar_t ty_c = 0.;
        soma_scalar_t tz_c = 0.;
        soma_scalar_t mx_c = 0.;
        soma_scalar_t my_c = 0.;
        soma_scalar_t mz_c = 0.;
	const unsigned int type = p->polymers[j].type;
	const unsigned int N = p->poly_arch[ p->poly_type_offset[ type ] ];
	for (unsigned int k = 0; k < N; k++) {	/*Loop over monomers */
	    const soma_scalar_t tx = p->polymers[j].msd_beads[k].x - p->polymers[j].beads[k].x;
	    const soma_scalar_t ty = p->polymers[j].msd_beads[k].y - p->polymers[j].beads[k].y;
	    const soma_scalar_t tz = p->polymers[j].msd_beads[k].z - p->polymers[j].beads[k].z;

            // Add up values for chain diffusion
            tx_c += p->polymers[j].beads[k].x;
            ty_c += p->polymers[j].beads[k].y;
            tz_c += p->polymers[j].beads[k].z;
            mx_c += p->polymers[j].msd_beads[k].x;
            my_c += p->polymers[j].msd_beads[k].y;
            mz_c += p->polymers[j].msd_beads[k].z;

	    counter[type*2 +0] += 1;
	    result[type*8 + 0] += tx * tx;
	    result[type*8 + 1] += ty * ty;
	    result[type*8 + 2] += tz * tz;
	    result[type*8 + 3] += (tx * tx + ty * ty + tz * tz);
	    }
	counter[type*2 + 1] += 1;
        result[type*8 + 4] += (tx_c-mx_c)*(tx_c-mx_c)/(N*N);
        result[type*8 + 5] += (ty_c-my_c)*(ty_c-my_c)/(N*N);
        result[type*8 + 6] += (tz_c-mz_c)*(tz_c-mz_c)/(N*N);
        result[type*8 + 7] += ((tx_c-mx_c)*(tx_c-mx_c) +(ty_c-my_c)*(ty_c-my_c) + (tz_c-mz_c)*(tz_c-mz_c) )/(N*N);
	}

    MPI_Allreduce(MPI_IN_PLACE, result, 8*p->n_poly_type, MPI_SOMA_SCALAR, MPI_SUM,
		  p->info_MPI.SOMA_MPI_Comm);
    MPI_Allreduce(MPI_IN_PLACE, counter, 2*p->n_poly_type, MPI_UINT64_T, MPI_SUM,
		  p->info_MPI.SOMA_MPI_Comm);

    //Looping over twice the number of poly types. But loop over half the elements
    // 8/2 = 4, because the norm for first and second half differ.
    for(unsigned int type=0 ; type < 2*p->n_poly_type; type++)
	for(unsigned int i=0; i < 4; i++)
	    if( counter[type] > 0)
		result[type*4 +i] /= (soma_scalar_t) counter[type];
    free(counter);
    }

void calc_acc_ratio(struct Phase*const p,soma_scalar_t*const acc_ratio)
    {
    // Acceptance ratio does not work with OpenACC builds.
#ifdef _OPENACC
    *acc_ratio = -1;
    return;
#endif
    static unsigned int last_time_call = 0;
    static soma_scalar_t last_acc = 0;
    if (last_time_call == 0 || p->time > last_time_call)
	last_time_call = p->time;
    else{//Quick exit, because the property has already been calculated for the time step.
	*acc_ratio = last_acc;
	return;
	}
    // exchange information among MPI processes
    uint64_t exchange[4] = {p->n_accepts,p->n_moves,p->n_accepts,p->n_moves};
    //Zero out result for next calc
    p->n_accepts = p->n_moves = 0;
    MPI_Reduce(exchange+2,exchange,2,MPI_UINT64_T,MPI_SUM,0,p->info_MPI.SOMA_MPI_Comm);
    if( exchange[1] > 0)
	last_acc = exchange[0]/(soma_scalar_t)exchange[1];
    else
	last_acc = 0;

    *acc_ratio = last_acc;
    }

void calc_non_bonded_energy(const struct Phase*const p, soma_scalar_t*const non_bonded_energy)
    {
    memset(non_bonded_energy,0,p->n_types*sizeof(soma_scalar_t));

    for(unsigned int type=0; type < p->n_types; type++)
	{
	for(uint64_t cell=0; cell < p->n_cells; cell++)
	    {
	    non_bonded_energy[type] +=
		p->omega_field_unified[cell + type*p->n_cells]
		* p->fields_unified[cell + type*p->n_cells ];
	    }
	}
    }

void calc_bonded_energy(const struct Phase*const p, soma_scalar_t*const bonded_energy)
    {
    memset(bonded_energy,0,NUMBER_SOMA_BOND_TYPES*sizeof(soma_scalar_t));

    for(unsigned int poly=0; poly < p->n_polymers; poly++)
	{
	const unsigned int type =p->polymers[poly].type;
	const unsigned int N = p->poly_arch[ p->poly_type_offset[ type ] ];

	for(unsigned int mono=0; mono < N ; mono++)
	    {

	    const int start = get_bondlist_offset(p->poly_arch[p->poly_type_offset[type] + mono + 1]);

	    if(start > 0)
		{
		int i = start;
		unsigned int end;
		do{
		    const uint32_t info = p->poly_arch[i++];
		    end = get_end(info);
		    const unsigned int bond_type = get_bond_type(info);
		    const int offset = get_offset(info);
		    if( offset > 0) //Select each bond only once, i<j
			{
			const int mono_j = mono + offset;

			const soma_scalar_t dx = p->polymers[poly].beads[mono].x - p->polymers[poly].beads[mono_j].x;
			const soma_scalar_t dy = p->polymers[poly].beads[mono].y - p->polymers[poly].beads[mono_j].y;
			const soma_scalar_t dz = p->polymers[poly].beads[mono].z - p->polymers[poly].beads[mono_j].z;
			const soma_scalar_t r2 = dx*dx + dy*dy + dz*dz;

			soma_scalar_t energy = 0;
			soma_scalar_t scale = 1.;
			switch (bond_type) {
			case HARMONICVARIABLESCALE:
			    scale = p->harmonic_normb_variable_scale;
			    /* intentionally falls through */
			case HARMONIC:
			    energy = p->harmonic_normb * r2 *scale;
			    break;

			case STIFF:
			    fprintf(stderr,
				    "ERROR: %s:%d stiff bond not yet implemented.\n",
				    __FILE__, __LINE__);
			    break;
			default:
			    fprintf(stderr, "ERROR: %s:%d unknow bond type appeared %d\n",
				    __FILE__, __LINE__,bond_type);
			    break;
			    }
			assert( bond_type < NUMBER_SOMA_BOND_TYPES);
			bonded_energy[ bond_type ] += energy;
			}
		    }while( end == 0);
		}
	    }
	}
    if( p->info_MPI.current_core == 0)
	MPI_Reduce(MPI_IN_PLACE,bonded_energy,NUMBER_SOMA_BOND_TYPES,MPI_SOMA_SCALAR,MPI_SUM,0,p->info_MPI.SOMA_MPI_Comm);
    else
	MPI_Reduce(bonded_energy, NULL ,NUMBER_SOMA_BOND_TYPES,MPI_SOMA_SCALAR,MPI_SUM,0,p->info_MPI.SOMA_MPI_Comm);
    }

int extent_ana_by_field(const soma_scalar_t*const data,const uint64_t n_data,const char*const name,const hid_t file_id)
    {
    herr_t status;
    hid_t dset = H5Dopen(file_id, name,H5P_DEFAULT);
    HDF5_ERROR_CHECK(dset);
    hid_t d_space = H5Dget_space(dset);
    HDF5_ERROR_CHECK(d_space);
    const unsigned int ndims = H5Sget_simple_extent_ndims(d_space);
    if( ndims != 2)
	{
	fprintf(stderr,"ERROR: %s:%d not the correct number of dimensions to extent the data set for %s.\n",__FILE__,__LINE__,name);
	return -1;
	}
    hsize_t dims[2];//ndims

    status = H5Sget_simple_extent_dims(d_space,dims,NULL);
    HDF5_ERROR_CHECK(status);

    hsize_t dims_new[2];//ndims
    dims_new[0] = dims[0] + 1;
    dims_new[1] = n_data;

    status = H5Dset_extent(dset,dims_new);
    HDF5_ERROR_CHECK(status);

    status = H5Sclose(d_space);
    HDF5_ERROR_CHECK(status);
    hid_t filespace = H5Dget_space(dset);
    HDF5_ERROR_CHECK(filespace);

    hsize_t dims_memspace[2]; //ndims
    dims_memspace[0] = dims_new[0] - dims[0];
    for(unsigned int i = 1; i< ndims;i++){
	dims_memspace[i] = dims[i];
	}

    hid_t memspace = H5Screate_simple(ndims,dims_memspace,NULL);
    hsize_t dims_offset[2];//ndims
    dims_offset[0] = dims[0];
    for(unsigned int i=1; i < ndims; i++)
	dims_offset[i] = 0;

    if(dims[0] > 0){
	status = H5Sselect_hyperslab(filespace,H5S_SELECT_SET,dims_offset,NULL,dims_memspace,NULL);
	HDF5_ERROR_CHECK(status);
	}
    else
	memspace = H5P_DEFAULT;

    status = H5Dwrite(dset,H5T_SOMA_NATIVE_SCALAR,memspace,filespace,H5P_DEFAULT,data);
    HDF5_ERROR_CHECK(status);


    status = H5Sclose(filespace);
    HDF5_ERROR_CHECK(status);
    status = H5Dclose(dset);
    HDF5_ERROR_CHECK(status);
    return 0;
    }

int extent_density_field(const struct Phase*const p)
    {
    char name[] = "/density_field";
    update_density_fields(p);
    //This is a draft verion that outputs only with a single core.
    //! \todo Parallel output
    if(p->info_MPI.current_core == 0 )
    	{
	herr_t status;

	hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	HDF5_ERROR_CHECK(plist_id);

	//Open the dataset and space
	hid_t dset = H5Dopen(p->ana_info.file_id, name,H5P_DEFAULT);
	HDF5_ERROR_CHECK(dset);
	hid_t d_space = H5Dget_space(dset);
	HDF5_ERROR_CHECK(d_space);

	// Extent the dimesion of the density field.
	const unsigned int ndims = H5Sget_simple_extent_ndims(d_space);
	if( ndims != 5)
	    {
	    fprintf(stderr,"ERROR: %s:%d not the correct number of dimensions to extent density field %s.\n",__FILE__,__LINE__,name);
	    return -1;
	    }
	hsize_t dims[5];//ndims
	status = H5Sget_simple_extent_dims(d_space,dims,NULL);
	HDF5_ERROR_CHECK(status);

	hsize_t dims_new[5];//ndims
	dims_new[0] = dims[0] + 1;
	dims_new[1] = p->n_types;
	assert(dims[1] == p->n_types);
	dims_new[2] = p->nx;
	assert(dims[2] == p->nx);
	dims_new[3] = p->ny;
	assert(dims[3] == p->ny);
	dims_new[4] = dims[4];
	assert(dims[4] == p->nz);

	status = H5Dset_extent(dset,dims_new);
	HDF5_ERROR_CHECK(status);

	status = H5Sclose(d_space);
	HDF5_ERROR_CHECK(status);
	status = H5Dclose(dset);
	HDF5_ERROR_CHECK(status);

	//Open the new filespace
	dset = H5Dopen(p->ana_info.file_id, name,H5P_DEFAULT);
	HDF5_ERROR_CHECK(dset);
	hid_t filespace = H5Dget_space(dset);
	HDF5_ERROR_CHECK(filespace);


	hsize_t dims_memspace[5];//ndims
	dims_memspace[0] = dims_new[0] - dims[0];
	dims_memspace[1] = p->n_types;
	dims_memspace[2] = p->nx;
	dims_memspace[3] = p->ny;
	dims_memspace[4] = p->nz;

	hid_t memspace = H5Screate_simple(ndims,dims_memspace,NULL);
	hsize_t dims_offset[5];//ndims
	dims_offset[0] = dims[0];
	for(unsigned int i=1; i < 5; i++)
	    dims_offset[i] = 0;

	status = H5Sselect_hyperslab(filespace,H5S_SELECT_SET,dims_offset,NULL,dims_memspace,NULL);
	HDF5_ERROR_CHECK(status);

	status = H5Dwrite(dset,H5T_NATIVE_UINT16,memspace,filespace,plist_id,p->fields_unified);
	HDF5_ERROR_CHECK(status);

	status = H5Sclose(memspace);
	HDF5_ERROR_CHECK(status);
	status = H5Sclose(filespace);
	HDF5_ERROR_CHECK(status);
	status = H5Dclose(dset);
	HDF5_ERROR_CHECK(status);
	}
    return 0;
    }

int analytics(struct Phase *const p)
    {
    static unsigned int last_time_call = 0;
    if (last_time_call == 0 || p->time > last_time_call)
	last_time_call = p->time;
    else			//Quick exit, because the property has already been calculated for the time step.
	return 0;
    bool written = false;

    // Re calculation
    if(p->ana_info.delta_mc_Re != 0 && p->time % p->ana_info.delta_mc_Re == 0)
	{
	update_self_phase(p);
	soma_scalar_t*const Re=(soma_scalar_t*const)malloc(4*p->n_poly_type*sizeof(soma_scalar_t));
	if(Re == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}
	calc_Re(p,Re);
	if(p->info_MPI.current_core == 0 )
	    extent_ana_by_field(Re, 4*p->n_poly_type, "/Re",p->ana_info.file_id);
	written = true;
	free(Re);
	}
    // dvar calculation
    if(p->ana_info.delta_mc_density_var != 0 && p->time % p->ana_info.delta_mc_density_var == 0)
	{
	update_self_phase(p);
	soma_scalar_t dvar[1];
	calc_dvar(p,dvar);
	if(p->info_MPI.current_core == 0 )
	    extent_ana_by_field(dvar, 1, "/density_var",p->ana_info.file_id);
	}
    // Gyration radius
    if (p->ana_info.delta_mc_Rg != 0 && p->time % p->ana_info.delta_mc_Rg == 0) {
	update_self_phase(p);
    	soma_scalar_t*const Rg = (soma_scalar_t*const)malloc(4*p->n_poly_type*sizeof(soma_scalar_t));
	if(Rg == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}
    	calc_Rg(p, Rg);
    	if (p->info_MPI.current_core == 0)
	    extent_ana_by_field(Rg, 4*p->n_poly_type, "/Rg",p->ana_info.file_id);
	free(Rg);
	written = true;
	}
    // Bond anisotropy
    if (p->ana_info.delta_mc_b_anisotropy != 0 && p->time % p->ana_info.delta_mc_b_anisotropy == 0) {
	update_self_phase(p);
	soma_scalar_t*const a=(soma_scalar_t*const)malloc(6*p->n_poly_type*sizeof(soma_scalar_t));
	if(a == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}
    	calc_anisotropy(p, a);
    	if (p->info_MPI.current_core == 0)
	    extent_ana_by_field(a, 6*p->n_poly_type, "/bond_anisotropy",p->ana_info.file_id);
	written = true;
	free(a);
	}
    // Acceptance ratio
    if (p->ana_info.delta_mc_acc_ratio != 0 && p->time % p->ana_info.delta_mc_acc_ratio == 0) {
	soma_scalar_t acc_ration;
    	calc_acc_ratio(p, &acc_ration);
    	if (p->info_MPI.current_core == 0)
	    extent_ana_by_field(&acc_ration, 1, "/acc_ratio",p->ana_info.file_id);
	written = true;
	}
    // MSD
    if (p->ana_info.delta_mc_MSD != 0 && p->time % p->ana_info.delta_mc_MSD == 0) {
    	soma_scalar_t*const msd = (soma_scalar_t*const)malloc(8*p->n_poly_type*sizeof(soma_scalar_t));
	if(msd == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}
	update_self_phase(p);
    	calc_MSD(p, msd);
    	if (p->info_MPI.current_core == 0)
	    extent_ana_by_field(msd, 8*p->n_poly_type, "/MSD",p->ana_info.file_id);
	written = true;
	free(msd);
	}
    // non-bonded energy calculation
    if(p->ana_info.delta_mc_non_bonded_energy != 0
       && p->time % p->ana_info.delta_mc_non_bonded_energy == 0)
	{
	soma_scalar_t*const nb_energy = (soma_scalar_t*const)malloc(p->n_types*sizeof(soma_scalar_t));
	if(nb_energy == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}
	update_self_phase(p);
	calc_non_bonded_energy(p, nb_energy);
	if(p->info_MPI.current_core == 0)
	    extent_ana_by_field(nb_energy, p->n_types, "/non_bonded_energy",p->ana_info.file_id);
	written = true;
	free(nb_energy);
	}
    // bonded energy calculation
    if(p->ana_info.delta_mc_bonded_energy != 0
       && p->time % p->ana_info.delta_mc_bonded_energy == 0)
	{
	soma_scalar_t*const b_energy = (soma_scalar_t*const)malloc(NUMBER_SOMA_BOND_TYPES*sizeof(soma_scalar_t));
	if(b_energy == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}
	update_self_phase(p);
	calc_bonded_energy(p, b_energy);
	if(p->info_MPI.current_core == 0)
	    extent_ana_by_field(b_energy, NUMBER_SOMA_BOND_TYPES, "/bonded_energy",p->ana_info.file_id);
	written = true;
	free(b_energy);
	}

    //Field
    if (p->ana_info.delta_mc_density_field != 0 && p->time % p->ana_info.delta_mc_density_field == 0) {
	// if we run on one core+gpu only we avoid the stepwise field transfer, so for writing fields
	// we need to update the local field
	if (p->info_MPI.Ncores == 1){
#pragma acc update self(p->fields_unified[0:p->n_cells*p->n_types])
	    }

    	//Collective IO, not yet.
    	extent_density_field(p);
	written = true;
    	}
    //dump
    if (p->ana_info.delta_mc_dump != 0 && p->time % p->ana_info.delta_mc_dump == 0) {
	update_self_phase(p);
	const unsigned int len = 1024;
	char filename[1024];//len
	memset(filename,'\0',len*sizeof(char));
	const unsigned int check = sprintf(filename,"%s_t%u.h5",p->ana_info.coord_filename,p->time);
	if(check >= len)
	    fprintf(stderr,"ERROR: %s %d cannot write filename.\n",__FILE__,__LINE__);
	else
	    write_config_hdf5(p,filename);
	}

    if(written && p->info_MPI.current_core == 0) //If something has been written to hdf5 flush it.
	{
	herr_t status = H5Fflush(p->ana_info.file_id,H5F_SCOPE_LOCAL);
	if( status < 0)
	    {
	    fprintf(stderr, "ERROR: HDF5 %s:%s:%d\n",__func__, __FILE__, __LINE__);
	    return -1;
	    }
	}

    return 0;
    }
