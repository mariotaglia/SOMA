/* Copyright (C) 2017 Ludwig Schneider
   Copyright (C) 2017 De-Wen Sun

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
#include <hdf5.h>
#include <stdlib.h>

#include "soma_util.h"
#include "shear_lib.h"
#include "io.h"


int read_bead_data(const char*const filename,bead_data*const b)
    {
    herr_t status;
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);

    hid_t file_id = H5Fopen(filename,H5F_ACC_RDONLY,plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);

    uint64_t Npolymers;
    status = read_hdf5(file_id,"/parameter/n_polymers",H5T_NATIVE_UINT64,plist_id,&Npolymers);
    HDF5_ERROR_CHECK2(status,"/parameter/n_polymers");

    unsigned int poly_arch_len;
    status = read_hdf5(file_id,"/parameter/poly_arch_length",H5T_NATIVE_UINT,plist_id,&poly_arch_len);
    HDF5_ERROR_CHECK2(status,"poly_arch_length");

    uint32_t*const poly_arch = (uint32_t*) malloc( poly_arch_len * sizeof(uint32_t));
    if(poly_arch == NULL){
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
	}
    status = read_hdf5(file_id,"/parameter/poly_arch",H5T_NATIVE_INT32,plist_id,poly_arch);
    HDF5_ERROR_CHECK2(status,"poly_arch");

    unsigned int n_poly_type;
    status = read_hdf5(file_id,"/parameter/n_poly_type",H5T_NATIVE_UINT,plist_id,&n_poly_type);
    HDF5_ERROR_CHECK2(status,"n_poly_type");
    //poly_type_offset
    int*const poly_type_offset = (int*)malloc(n_poly_type * sizeof(int));
    if(poly_type_offset == NULL){
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
	}
    status = read_hdf5(file_id,"/parameter/poly_type_offset",H5T_NATIVE_INT,plist_id,poly_type_offset);
    HDF5_ERROR_CHECK2(status,"poly_type_offset");


    unsigned int *const poly_type = (unsigned int *const) malloc( Npolymers * sizeof(unsigned int));
    if (poly_type == NULL) {
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
	}
    status = read_hdf5(file_id,"/poly_type",H5T_NATIVE_UINT,plist_id,poly_type);
    HDF5_ERROR_CHECK2(status,"poly_type");

    unsigned int*const number_of_beads= (unsigned int*const) malloc(Npolymers*sizeof(unsigned int));
    if (number_of_beads == NULL) {
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
	}

    unsigned int max_n_beads=0;
    for(unsigned int i=0; i < Npolymers; i++)
	{
	const unsigned int N= poly_arch[ poly_type_offset[ poly_type[i] ] ];
	number_of_beads[i] = N;
	if( N > max_n_beads)
	    max_n_beads = N;
	}

    Monomer *const beads = (Monomer * const) malloc(Npolymers * max_n_beads * sizeof(Monomer));
    if (beads == NULL) {
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
	}

    if(H5Lexists(file_id,"/beads",H5P_DEFAULT) > 0)
    	{
	hid_t monomer_memtype = get_monomer_memtype();
	read_hdf5( file_id, "/beads",monomer_memtype,plist_id,beads);
	if ((status = H5Tclose(monomer_memtype)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    0, __FILE__, __LINE__, status);
	    return status;
	    }
	}
    else
	{
	fprintf(stderr,"ERROR: Given file does not contain bead data\n.");
	return -2;
	}

    free(poly_arch);
    free(poly_type_offset);
    free(poly_type);
    if ((status = H5Fclose(file_id)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		0, __FILE__, __LINE__, status);
	return status;
	}


    //Swap
    b->N_polymers = Npolymers;
    b->number_of_beads = number_of_beads;
    b->max_n_beads = max_n_beads;
    b->beads = beads;
    return 0;
    }

int free_bead_data(bead_data*b)
    {
    free(b->beads);
    free(b->number_of_beads);
    return 0;
    }

int write_bead_data(const char*const filename,bead_data*const b)
    {
    const hid_t memtype = get_monomer_memtype();

    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);

    hid_t file_id = H5Fopen(filename, H5F_ACC_RDWR, plist_id);
    HDF5_ERROR_CHECK2(file_id,"Open file");
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);

    //Open old dataset
    hid_t dataset = H5Dopen(file_id,"/beads",H5P_DEFAULT);
    HDF5_ERROR_CHECK2(dataset,"open beads");

    herr_t status;
    status = H5Dwrite(dataset,memtype,H5S_ALL,H5S_ALL,plist_id, b->beads);
    HDF5_ERROR_CHECK2(status,"write beads");
    status = H5Dclose(dataset);
    HDF5_ERROR_CHECK2(status,"close beads");

    H5Pclose(plist_id);
    status = H5Fclose(file_id);
    HDF5_ERROR_CHECK2(status,"close file");
    return 0;
    }

int apply_shear(bead_data*const b,const unsigned int gradient_direction,
                const unsigned int flow_direction,const soma_scalar_t amplitude)
    {
    uint64_t ciac,gaug;
    Monomer x;
    if( gradient_direction > 2)
    {
     fprintf(stderr,"Invalid gradient direction\n");
     return -1;
    }
    if( flow_direction > 2)
    {
     fprintf(stderr,"Invalid flow direction.\n");
     return -2;
    }
    if( gradient_direction == flow_direction )
     {
      fprintf(stderr,"Ivalid combination of flow and gradient direction.\n");
      return -3;
     }

    for(ciac=0;ciac<b->N_polymers;ciac++)
    {
     for(gaug=0;gaug<b->number_of_beads[ciac];gaug++)
     {

      x=b->beads[get_2d_index(ciac,gaug,b->max_n_beads)];

      if((gradient_direction==0)&&(flow_direction==1))
      {
       x.y+=x.x*amplitude;
      }
      if((gradient_direction==0)&&(flow_direction==2))
      {
       x.z+=x.x*amplitude;
      }

      if((gradient_direction==1)&&(flow_direction==0))
      {
       x.x+=x.y*amplitude;
      }
      if((gradient_direction==1)&&(flow_direction==2))
      {
       x.z+=x.y*amplitude;
      }

      if((gradient_direction==2)&&(flow_direction==0))
      {
       x.x+=x.z*amplitude;
      }
      if((gradient_direction==2)&&(flow_direction==1))
      {
       x.y+=x.z*amplitude;
      }

      b->beads[ get_2d_index(ciac,gaug,b->max_n_beads) ] = x;

     }
    }
    return 0;
    }

unsigned int get_2d_index(const unsigned int poly,const unsigned int mono,
			  const unsigned int max_n_beads)
    {
    return poly*max_n_beads + mono;
    }
