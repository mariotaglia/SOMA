/* Copyright (C) 2016-2017 Ludwig Schneider

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

//! \file ana_info.c
//! \brief Implementation of ana_info.h

#include "ana_info.h"
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "phase.h"

int init_ana(struct Phase * const p,const char*const filename,const char*const coord_filename)
    {
    //******** START EDIT FOR NEW OBSERVABLES HERE**********
    p->ana_info.delta_mc_Re = 0;
    p->ana_info.delta_mc_Rg = 0;
    p->ana_info.delta_mc_b_anisotropy = 0;
    p->ana_info.delta_mc_density_field = 0;
    p->ana_info.delta_mc_acc_ratio = 0;
    p->ana_info.delta_mc_MSD = 0;
    p->ana_info.delta_mc_dump = 0;
    p->ana_info.delta_mc_density_var = 0;
    p->ana_info.delta_mc_non_bonded_energy = 0;
    p->ana_info.delta_mc_bonded_energy = 0;
    p->ana_info.delta_mc_string_field = 0;
    //******** END EDIT FOR NEW OBSERVABLES HERE************
    p->ana_info.filename = NULL;
    p->ana_info.coord_filename = NULL;
    p->ana_info.file_id = -1;
    //Copy no ana conf as backup to device.
#pragma acc update device(p->ana_info)

    if(filename == NULL)
	return 0;
    const unsigned int str_len = strlen(filename);
    if(str_len < 2) //Exclude empty strings.
	return 0;

    p->ana_info.filename = (char*)malloc( (str_len+1)*sizeof(char) );
    if(p->ana_info.filename == NULL)
	{
	fprintf(stderr, "ERROR: Malloc %s:%s:%d\n",__func__, __FILE__, __LINE__);
	return -1;
	}
    strcpy(p->ana_info.filename,filename);

    herr_t status;
    //Set up file access property list with parallel I/O access
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    if (p->info_MPI.sim_size > 1){
	status = H5Pset_fapl_mpio(plist_id, p->info_MPI.SOMA_comm_sim,
				  MPI_INFO_NULL);
	HDF5_ERROR_CHECK(status);
	}

    //Write to the file only for the attr of density_field.
    hid_t file_id_tmp = H5Fopen(p->ana_info.filename, H5F_ACC_RDWR, plist_id);
    status = file_id_tmp;
    if( status < 0){
	fprintf(stderr, "ERROR: HDF5 %s:%s:%d\n",__func__, __FILE__, __LINE__);
	return -1;
    }
    p->ana_info.file_id = file_id_tmp;

    /* H5Pclose(plist_id); */
    /* plist_id = H5Pcreate(H5P_DATASET_XFER); */
    /* if (p->info_MPI.Ncores > 1) */
    /* 	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE); */
    //******** START EDIT FOR NEW OBSERVABLES HERE**********
#ifdef SOMA_NUM_OBS
#error "Namespace violation SOMA_NUM_OBS already defined."
#endif//SOMA_NUM_OBS
//! Number of known observables to SOMA.
//! \private
#define SOMA_NUM_OBS 11
    const char* names[SOMA_NUM_OBS];
    unsigned int*delta_mc[SOMA_NUM_OBS];
    names[0] = "/Re";
    delta_mc[0] = &(p->ana_info.delta_mc_Re);
    names[1] = "/Rg";
    delta_mc[1] = &(p->ana_info.delta_mc_Rg);
    names[2] = "/bond_anisotropy";
    delta_mc[2] = &(p->ana_info.delta_mc_b_anisotropy);
    names[3]="/density_field";
    delta_mc[3] = &(p->ana_info.delta_mc_density_field);
    names[4]="/acc_ratio";
    delta_mc[4] = &(p->ana_info.delta_mc_acc_ratio);
    names[5]="/MSD";
    delta_mc[5] = &(p->ana_info.delta_mc_MSD);
    names[6]="/dump";
    delta_mc[6] = &(p->ana_info.delta_mc_dump);
    names[7] = "/density_var";
    delta_mc[7] = &(p->ana_info.delta_mc_density_var);
    names[8] = "/non_bonded_energy";
    delta_mc[8] = &(p->ana_info.delta_mc_non_bonded_energy);
    names[9] = "/bonded_energy";
    delta_mc[9] = &(p->ana_info.delta_mc_bonded_energy);
    names[10] = "/string_field";
    delta_mc[10] = &(p->ana_info.delta_mc_string_field);
    //******** END EDIT FOR NEW OBSERVABLES HERE************
    for(unsigned int i=0; i< SOMA_NUM_OBS; i++)
	{
	hid_t dataset = H5Dopen2(p->ana_info.file_id,names[i],H5P_DEFAULT);
	status = dataset;
	if( status <  0)
	    fprintf(stderr,"WARNING %s:%d ana element %s not existing (old file version?).\n",__FILE__,__LINE__,names[i]);
	else
	    {
	    hid_t attr = H5Aopen_by_name(dataset,names[i],"DeltaMC",H5P_DEFAULT,H5P_DEFAULT);
	    status = attr;
	    if( status <  0)
		fprintf(stderr,"WARNING %s:%d deltamc for %s not existing (old file version?).\n",__FILE__,__LINE__,names[i]);
	    else
		{
		status = H5Aread(attr,H5T_NATIVE_UINT,delta_mc[i]);
		if(status < 0)
		    fprintf(stderr,"WARNING: %s:%d deltamc for %s unreadable.\n",__FILE__,__LINE__,names[i]);

		if( H5Aclose(attr) <  0)
		    {
		    fprintf(stderr, "ERROR: HDF5 %s:%s:%d\n",__func__, __FILE__, __LINE__);
		    return -1;
		    }
		}
	    if( H5Dclose(dataset) < 0)
		{
		fprintf(stderr, "ERROR: HDF5 %s:%s:%d\n",__func__, __FILE__, __LINE__);
		return -1;
		}
	    }
	}
#undef SOMA_NUM_OBS
    //**************** TEST YOUR OBS EDIT LIKE THIS************************
    /* printf("%d %d %d %d %d %d %d %d\n", */
    /* 	   p->ana_info.delta_mc_Re,p->ana_info.delta_mc_Rg, */
    /* 	   p->ana_info.delta_mc_acc_ratio,p->ana_info.delta_mc_b_anisotropy, */
    /* 	   p->ana_info.delta_mc_density_field,p->ana_info.delta_mc_dump,p->ana_info.delta_mc_non_bonded_energy, */
    /* 	   p->ana_info.delta_mc_bonded_energy); */


    if(p->ana_info.delta_mc_dump > 0)
      {
	if(coord_filename != NULL)
	  {
	    const unsigned int new_len = strlen(coord_filename)+1;
	    p->ana_info.coord_filename = (char*)malloc(new_len * sizeof(char));
	    if (p->ana_info.coord_filename == NULL) {
	      fprintf(stderr,"ERROR: %s %d malloc. So no automatic dumping possible.\n",__FILE__,__LINE__);
	      p->ana_info.delta_mc_dump = 0;
	    }
	    memset(p->ana_info.coord_filename,'\0',new_len*sizeof(char));
	    strcpy(p->ana_info.coord_filename,coord_filename);
	    unsigned int first_dot;
	    for (first_dot = new_len - 1; first_dot != 0; first_dot--)
	      if (p->ana_info.coord_filename[first_dot] == '.' || first_dot == 0)
		{
		  //remove the ending
		  memset(p->ana_info.coord_filename+first_dot,'\0',(new_len-first_dot)*sizeof(char));
		  first_dot = 0;
		  break;
		}
	  }
	else
	  {
	    p->ana_info.delta_mc_dump = 0;
	    fprintf(stderr,"WARNING: %s %d no coord_file specified. So no automatic dumping possible.\n",__FILE__,__LINE__);
	  }
#pragma acc update device(p->ana_info)
      }

    //Check or add the additional attributes for the density fields
    if(p->ana_info.delta_mc_density_field > 0)
	{
	//Exception safety:
	const unsigned int tmd_delta_mc = p->ana_info.delta_mc_density_field;
	p->ana_info.delta_mc_density_field = 0;
#pragma acc update device(p->ana_info)
	const hsize_t three = 3;
	const hsize_t one = 1;

	hid_t dataset = H5Dopen2(p->ana_info.file_id,"/density_field",H5P_DEFAULT);
	status = dataset;
	HDF5_ERROR_CHECK(status);

	//First try to open the attributes. If they exist check for
	//consistence with the current data. If they do not exist
	//create them and fill with data.

	htri_t nxyz_exists = H5Aexists(dataset,"nxyz");
	if( nxyz_exists > 0){//Nxyz attr exists
	    hid_t nxyz_attr = H5Aopen_by_name(dataset,"/density_field","nxyz",H5P_DEFAULT,H5P_DEFAULT);
	    unsigned int nxyz_tmp[3];
	    status = H5Aread(nxyz_attr,H5T_NATIVE_UINT,nxyz_tmp);
	    HDF5_ERROR_CHECK(status);
	    if( (nxyz_tmp[0] != p->nx) || (nxyz_tmp[1] != p->ny) || (nxyz_tmp[2] != p->nz)){
		    fprintf(stderr,"WARNING: existing nxyz attr of density_field"
			    "does not match with the system data. No density field output possible. (Try with fresh ana.h5 file.).\n");
		    return 0;
		}
	    status = H5Aclose(nxyz_attr);
	    HDF5_ERROR_CHECK(status);
	    }
	else{ //Create the nxyz attr
	    hid_t nxyz_dataspace = H5Screate_simple(1, &three, NULL);
	    unsigned int nxyz_tmp[3] = {p->nx,p->ny,p->nz};
	    hid_t nxyz_attr = H5Acreate2(dataset,"nxyz",H5T_STD_U32LE,nxyz_dataspace,H5P_DEFAULT,H5P_DEFAULT);
	    HDF5_ERROR_CHECK(nxyz_attr);
	    status = H5Awrite(nxyz_attr,H5T_NATIVE_UINT,nxyz_tmp);
	    HDF5_ERROR_CHECK(status);
	    status = H5Aclose(nxyz_attr);
	    HDF5_ERROR_CHECK(status);
	    }

	htri_t ntypes_exists = H5Aexists(dataset,"ntypes");
	if( ntypes_exists > 0){//n_types attr exists
	    hid_t ntypes_attr = H5Aopen_by_name(dataset,"/density_field","ntypes",H5P_DEFAULT,H5P_DEFAULT);
	    unsigned int ntypes_tmp;
	    status = H5Aread(ntypes_attr,H5T_NATIVE_UINT,&ntypes_tmp);
	    HDF5_ERROR_CHECK(status);
	    if( (ntypes_tmp != p->n_types) ) {
		fprintf(stderr,"Error: existing ntypes attr of density_field"
			"does not match with the system data. No density field output possible. (Try with fresh ana.h5 file.).\n");
		    return 0;
		}
	    status = H5Aclose(ntypes_attr);
	    HDF5_ERROR_CHECK(status);
	    }
	else{ //Create the nxyz attr
	    hid_t ntypes_dataspace = H5Screate_simple(1, &one, NULL);
	    hid_t ntypes_attr = H5Acreate2(dataset,"ntypes",H5T_STD_U32LE,
					   ntypes_dataspace,H5P_DEFAULT,H5P_DEFAULT);
	    HDF5_ERROR_CHECK(ntypes_attr);
	    status = H5Awrite(ntypes_attr,H5T_NATIVE_UINT,&(p->n_types));
	    HDF5_ERROR_CHECK(status);
	    status = H5Aclose(ntypes_attr);
	    HDF5_ERROR_CHECK(status);
	    }

	htri_t lxyz_exists = H5Aexists(dataset,"lxyz");
	if( lxyz_exists > 0){//lxyz attr exists
	    hid_t lxyz_attr = H5Aopen_by_name(dataset,"/density_field","lxyz",H5P_DEFAULT,H5P_DEFAULT);
	    soma_scalar_t lxyz_tmp[3];
	    status = H5Aread(lxyz_attr,H5T_SOMA_NATIVE_SCALAR,lxyz_tmp);
	    HDF5_ERROR_CHECK(status);
	    if( (lxyz_tmp[0] != p->Lx) || (lxyz_tmp[1] != p->Ly) || (lxyz_tmp[2] != p->Lz)){
		    fprintf(stderr,"WARNING: existing lxyz attr of density_field"
			    "does not match with the system data. No dumping possible.\n");
		    return 0;
		}
	    status = H5Aclose(lxyz_attr);
	    HDF5_ERROR_CHECK(status);
	    }
	else{ //Create the lxyz attr
	    hid_t lxyz_dataspace = H5Screate_simple(1, &three, NULL);
	    soma_scalar_t lxyz_tmp[3] = {p->Lx,p->Ly,p->Lz};
	    hid_t lxyz_attr = H5Acreate2(dataset,"lxyz",H5T_SOMA_FILE_SCALAR,lxyz_dataspace,H5P_DEFAULT,H5P_DEFAULT);
	    HDF5_ERROR_CHECK(lxyz_attr);
	    status = H5Awrite(lxyz_attr,H5T_SOMA_NATIVE_SCALAR,lxyz_tmp);
	    HDF5_ERROR_CHECK(status);
	    status = H5Aclose(lxyz_attr);
	    HDF5_ERROR_CHECK(status);
	    }

	htri_t typescale_exists = H5Aexists(dataset,"typescale");
	if( typescale_exists > 0){//type_scale attr exists
	    hid_t typescale_attr = H5Aopen_by_name(dataset,"/density_field","typescale",H5P_DEFAULT,H5P_DEFAULT);
	    soma_scalar_t*const typescale_tmp =(soma_scalar_t*const)malloc(p->n_types*sizeof(soma_scalar_t));
	    if(typescale_tmp == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}
	    status = H5Aread(typescale_attr,H5T_SOMA_NATIVE_SCALAR,typescale_tmp);
	    HDF5_ERROR_CHECK(status);
	    bool match = true;
	    for(unsigned int i=0; i < p->n_types; i++)
		if( typescale_tmp[i] != p->field_scaling_type[i])
		    match = false;
	    if( !match ){
		fprintf(stderr,"Error: existing typescale attr of density_field"
			"does not match with the system data. No density field output possible. (Try with fresh ana.h5 file.)\n");
		return 0;
		}
	    status = H5Aclose(typescale_attr);
	    HDF5_ERROR_CHECK(status);
	    free(typescale_tmp);
	    }
	else{ //Create the typescale attr
	    const hsize_t ntypes = p->n_types;
	    hid_t typescale_dataspace = H5Screate_simple(1, &ntypes, NULL);
	    hid_t typescale_attr = H5Acreate2(dataset,"typescale",H5T_SOMA_FILE_SCALAR,
					      typescale_dataspace,H5P_DEFAULT,H5P_DEFAULT);
	    HDF5_ERROR_CHECK(typescale_attr);
	    status = H5Awrite(typescale_attr,H5T_SOMA_NATIVE_SCALAR,p->field_scaling_type);
	    HDF5_ERROR_CHECK(status);
	    status = H5Aclose(typescale_attr);
	    HDF5_ERROR_CHECK(status);
	    }

	// Check and set the dimesion of the density field.
	hid_t d_space = H5Dget_space(dataset);
	HDF5_ERROR_CHECK(d_space);
	const unsigned int ndims = H5Sget_simple_extent_ndims(d_space);
	if( ndims != 5)
	    {fprintf(stderr,"ERROR: %s:%d not the correct number of dimensions for a density field.\n",__FILE__,__LINE__);return -3;}
	hsize_t dims[5];//ndims
	status = H5Sget_simple_extent_dims(d_space,dims,NULL);
	HDF5_ERROR_CHECK(status);
	if( dims[1] != p->n_types || dims[2] != p->nx || dims[3] != p->ny || dims[4] != p->nz )
	    {
	    fprintf(stderr,"Error: The density field dimensions do not match! No density field output possible. (Try with fresh ana.h5 file.)\n");
	    return 0;
	    }
	status = H5Dclose(dataset);
	HDF5_ERROR_CHECK(status);

	//If everything was successful set the ana period again.
	p->ana_info.delta_mc_density_field = tmd_delta_mc;
	}


  if(p->ana_info.delta_mc_string_field> 0)
{
  //Exception safety:
  const unsigned int tmd_delta_mc_string = p->ana_info.delta_mc_string_field;
  p->ana_info.delta_mc_string_field = 0;
  #pragma acc update device(p->ana_info)
  const hsize_t three = 3;
  const hsize_t one = 1;

  hid_t dataset = H5Dopen2(p->ana_info.file_id,"/string_field",H5P_DEFAULT);
  status = dataset;
  HDF5_ERROR_CHECK(status);

  //First try to open the attributes. If they exist check for
  //consistence with the current data. If they do not exist
  //create them and fill with data.

  htri_t nxyz_exists = H5Aexists(dataset,"nxyz");
  if( nxyz_exists > 0){//Nxyz attr exists
      hid_t nxyz_attr = H5Aopen_by_name(dataset,"/string_field","nxyz",H5P_DEFAULT,H5P_DEFAULT);
      unsigned int nxyz_tmp[3];
      status = H5Aread(nxyz_attr,H5T_NATIVE_UINT,nxyz_tmp);
      HDF5_ERROR_CHECK(status);
      if( (nxyz_tmp[0] != p->nx) || (nxyz_tmp[1] != p->ny) || (nxyz_tmp[2] != p->nz)){
        fprintf(stderr,"WARNING: existing nxyz attr of string_field"
          "does not match with the system data. No density field output possible. (Try with fresh ana.h5 file.).\n");
        return 0;
    }
      status = H5Aclose(nxyz_attr);
      HDF5_ERROR_CHECK(status);
      }
  else{ //Create the nxyz attr
      hid_t nxyz_dataspace = H5Screate_simple(1, &three, NULL);
      unsigned int nxyz_tmp[3] = {p->nx,p->ny,p->nz};
      hid_t nxyz_attr = H5Acreate2(dataset,"nxyz",H5T_STD_U32LE,nxyz_dataspace,H5P_DEFAULT,H5P_DEFAULT);
      HDF5_ERROR_CHECK(nxyz_attr);
      status = H5Awrite(nxyz_attr,H5T_NATIVE_UINT,nxyz_tmp);
      HDF5_ERROR_CHECK(status);
      status = H5Aclose(nxyz_attr);
      HDF5_ERROR_CHECK(status);
      }

  htri_t ntypes_exists = H5Aexists(dataset,"ntypes");
  if( ntypes_exists > 0){//n_types attr exists
      hid_t ntypes_attr = H5Aopen_by_name(dataset,"/string_field","ntypes",H5P_DEFAULT,H5P_DEFAULT);
      unsigned int ntypes_tmp;
      status = H5Aread(ntypes_attr,H5T_NATIVE_UINT,&ntypes_tmp);
      HDF5_ERROR_CHECK(status);
      if( (ntypes_tmp != p->n_types) ) {
    fprintf(stderr,"Error: existing ntypes attr of string_field"
      "does not match with the system data. No density field output possible. (Try with fresh ana.h5 file.).\n");
        return 0;
    }
      status = H5Aclose(ntypes_attr);
      HDF5_ERROR_CHECK(status);
      }
  else{ //Create the nxyz attr
      hid_t ntypes_dataspace = H5Screate_simple(1, &one, NULL);
      hid_t ntypes_attr = H5Acreate2(dataset,"ntypes",H5T_STD_U32LE,
             ntypes_dataspace,H5P_DEFAULT,H5P_DEFAULT);
      HDF5_ERROR_CHECK(ntypes_attr);
      status = H5Awrite(ntypes_attr,H5T_NATIVE_UINT,&(p->n_types));
      HDF5_ERROR_CHECK(status);
      status = H5Aclose(ntypes_attr);
      HDF5_ERROR_CHECK(status);
      }

  htri_t lxyz_exists = H5Aexists(dataset,"lxyz");
  if( lxyz_exists > 0){//lxyz attr exists
      hid_t lxyz_attr = H5Aopen_by_name(dataset,"/string_field","lxyz",H5P_DEFAULT,H5P_DEFAULT);
      soma_scalar_t lxyz_tmp[3];
      status = H5Aread(lxyz_attr,H5T_SOMA_NATIVE_SCALAR,lxyz_tmp);
      HDF5_ERROR_CHECK(status);
      if( (lxyz_tmp[0] != p->Lx) || (lxyz_tmp[1] != p->Ly) || (lxyz_tmp[2] != p->Lz)){
        fprintf(stderr,"WARNING: existing lxyz attr of string_field"
          "does not match with the system data. No dumping possible.\n");
        return 0;
    }
      status = H5Aclose(lxyz_attr);
      HDF5_ERROR_CHECK(status);
      }
  else{ //Create the lxyz attr
      hid_t lxyz_dataspace = H5Screate_simple(1, &three, NULL);
      soma_scalar_t lxyz_tmp[3] = {p->Lx,p->Ly,p->Lz};
      hid_t lxyz_attr = H5Acreate2(dataset,"lxyz",H5T_SOMA_FILE_SCALAR,lxyz_dataspace,H5P_DEFAULT,H5P_DEFAULT);
      HDF5_ERROR_CHECK(lxyz_attr);
      status = H5Awrite(lxyz_attr,H5T_SOMA_NATIVE_SCALAR,lxyz_tmp);
      HDF5_ERROR_CHECK(status);
      status = H5Aclose(lxyz_attr);
      HDF5_ERROR_CHECK(status);
      }

  htri_t typescale_exists = H5Aexists(dataset,"typescale");
  if( typescale_exists > 0){//type_scale attr exists
      hid_t typescale_attr = H5Aopen_by_name(dataset,"/string_field","typescale",H5P_DEFAULT,H5P_DEFAULT);
      soma_scalar_t*const typescale_tmp =(soma_scalar_t*const)malloc(p->n_types*sizeof(soma_scalar_t));
      if(typescale_tmp == NULL){fprintf(stderr,"ERROR: Malloc %s:%d \n",__FILE__,__LINE__);return -2;}
      status = H5Aread(typescale_attr,H5T_SOMA_NATIVE_SCALAR,typescale_tmp);
      HDF5_ERROR_CHECK(status);
      bool match = true;
      for(unsigned int i=0; i < p->n_types; i++)
    if( typescale_tmp[i] != p->field_scaling_type[i])
        match = false;
      if( !match ){
    fprintf(stderr,"Error: existing typescale attr of string_field"
      "does not match with the system data. No string_field output possible. (Try with fresh ana.h5 file.)\n");
    return 0;
    }
      status = H5Aclose(typescale_attr);
      HDF5_ERROR_CHECK(status);
      free(typescale_tmp);
      }
  else{ //Create the typescale attr
      const hsize_t ntypes = p->n_types;
      hid_t typescale_dataspace = H5Screate_simple(1, &ntypes, NULL);
      hid_t typescale_attr = H5Acreate2(dataset,"typescale",H5T_SOMA_FILE_SCALAR,
                typescale_dataspace,H5P_DEFAULT,H5P_DEFAULT);
      HDF5_ERROR_CHECK(typescale_attr);
      status = H5Awrite(typescale_attr,H5T_SOMA_NATIVE_SCALAR,p->field_scaling_type);
      HDF5_ERROR_CHECK(status);
      status = H5Aclose(typescale_attr);
      HDF5_ERROR_CHECK(status);
      }

  // Check and set the dimesion of the density field.
  hid_t d_space = H5Dget_space(dataset);
  HDF5_ERROR_CHECK(d_space);
  const unsigned int ndims = H5Sget_simple_extent_ndims(d_space);
  if( ndims != 5)
      {fprintf(stderr,"ERROR: %s:%d not the correct number of dimensions for a string_field.\n",__FILE__,__LINE__);return -3;}
  hsize_t dims[5];//ndims
  status = H5Sget_simple_extent_dims(d_space,dims,NULL);
  HDF5_ERROR_CHECK(status);
  if( dims[1] != p->n_types || dims[2] != p->nx || dims[3] != p->ny || dims[4] != p->nz )
      {
      fprintf(stderr,"Error: The string_field dimensions do not match! No string_field output possible. (Try with fresh ana.h5 file.)\n");
      return 0;
      }
  status = H5Dclose(dataset);
  HDF5_ERROR_CHECK(status);

  //If everything was successful set the ana period again.
  p->ana_info.delta_mc_string_field = tmd_delta_mc_string;
}




    p->end_mono = NULL;
    if( p->ana_info.delta_mc_Re > 0)
	{
	const unsigned int tmp = p->ana_info.delta_mc_Re;
	p->ana_info.delta_mc_Re = 0;
	p->end_mono = (unsigned int*)malloc(2*p->n_poly_type*sizeof(unsigned int));
	if(p->end_mono == NULL)
	    {
	    fprintf(stderr,"Malloc error: %s:%d .\n",__FILE__,__LINE__);
	    return -1;
	    }

	hid_t dataset = H5Dopen2(p->ana_info.file_id,"/Re",H5P_DEFAULT);
	hid_t status = dataset;
	HDF5_ERROR_CHECK(status);

	htri_t end_mono_exists = H5Aexists(dataset,"end_mono");
	if( end_mono_exists > 0)
	    {//end_mono attr exists
	    hid_t end_mono_attr = H5Aopen_by_name(dataset,"/Re","end_mono",H5P_DEFAULT,H5P_DEFAULT);

	    hid_t d_space = H5Aget_space(end_mono_attr);
	    HDF5_ERROR_CHECK(d_space);
	    const unsigned int ndims = H5Sget_simple_extent_ndims(d_space);
	    if(ndims != 2)
		{
		fprintf(stderr,"ERROR: %s:%d wrong dimensions for end_mono.\n",__FILE__,__LINE__);
		return -1;
		}
	    hsize_t dims[2];//ndims
	    status = H5Sget_simple_extent_dims(d_space,dims,NULL);
	    HDF5_ERROR_CHECK(status);
	    if( dims[0] != p->n_poly_type )
		{
		fprintf(stderr,"ERROR: %s:%d end_mono for %d polymers, but expected is %d.\n",
			__FILE__,__LINE__,(int)dims[0],p->n_poly_type);
		return -1;
		}
	    if( dims[1] != 2 )
		{
		fprintf(stderr,"ERROR: %s:%d end_mono for %d elements for start/end, but expected is 2.\n",
			__FILE__,__LINE__,(int)dims[1]);
		return -1;
		}

	    status = H5Aread(end_mono_attr,H5T_NATIVE_UINT,p->end_mono);
	    HDF5_ERROR_CHECK(status);

	    status = H5Aclose(end_mono_attr);
	    HDF5_ERROR_CHECK(status);
	    //Sanatize the input values for end_mono
	    for(unsigned int poly_type = 0; poly_type < p->n_poly_type; poly_type++)
		{
		const unsigned int N = p->poly_arch[ p->poly_type_offset[ poly_type ] ];
		if( p->end_mono[poly_type*2 + 0] >= N || p->end_mono[poly_type*2 + 1] >= N)
		    {
		    fprintf(stderr,"ERROR: %s:%d EndMono contains invalid configuration. %d %d %d",__FILE__,__LINE__,p->end_mono[poly_type*2 + 0],p->end_mono[poly_type*2 + 1], N);
		    return -3;
		    }
		}
	    }
	else
	    {
	    fprintf(stderr,"ERROR: analysis of Re requires an end_mono attribute.\n");
	    return -2;
	    }
	status = H5Dclose(dataset);
	HDF5_ERROR_CHECK(status);

	p->ana_info.delta_mc_Re = tmp;
	}

    if( p->ana_info.delta_mc_non_bonded_energy > 0)
	{//Sanitize input
	const unsigned int tmp = p->ana_info.delta_mc_non_bonded_energy;
	p->ana_info.delta_mc_non_bonded_energy = 0;
	hid_t dataset = H5Dopen2(p->ana_info.file_id,"/non_bonded_energy",H5P_DEFAULT);
	hid_t status = dataset;
	HDF5_ERROR_CHECK(status);
	hid_t d_space = H5Dget_space(dataset);
	const unsigned int ndims = H5Sget_simple_extent_ndims(d_space);
	if(ndims != 2)
	    {
	    fprintf(stderr,"ERROR: %s:%d wrong dimensions for non_bonded_energy.\n",__FILE__,__LINE__);
	    return -1;
	    }
	hsize_t dims[2];//ndims
	status = H5Sget_simple_extent_dims(d_space,dims,NULL);
	HDF5_ERROR_CHECK(status);
	if( dims[1] != p->n_types)
	    {
	    fprintf(stderr,"Error: nonbonded energy dataspace is incorrect.\n");
	    return -2;
	    }
	status = H5Dclose(dataset);
	HDF5_ERROR_CHECK(status);

	p->ana_info.delta_mc_non_bonded_energy = tmp;
	}

    if( p->ana_info.delta_mc_bonded_energy > 0)
	{//Sanitize input
	const unsigned int tmp = p->ana_info.delta_mc_bonded_energy;
	p->ana_info.delta_mc_bonded_energy = 0;
	hid_t dataset = H5Dopen2(p->ana_info.file_id,"/bonded_energy",H5P_DEFAULT);
	hid_t status = dataset;
	HDF5_ERROR_CHECK(status);
	hid_t d_space = H5Dget_space(dataset);
	const unsigned int ndims = H5Sget_simple_extent_ndims(d_space);
	if(ndims != 2)
	    {
	    fprintf(stderr,"ERROR: %s:%d wrong dimensions for bonded_energy.\n",__FILE__,__LINE__);
	    return -1;
	    }
	hsize_t dims[2];//ndims
	status = H5Sget_simple_extent_dims(d_space,dims,NULL);
	HDF5_ERROR_CHECK(status);
	if( dims[1] != NUMBER_SOMA_BOND_TYPES)
	    {
	    fprintf(stderr,"Error: bonded energy dataspace is incorrect. (Might be an incompatible version problem.)\n");
	    return -2;
	    }
	status = H5Dclose(dataset);
	HDF5_ERROR_CHECK(status);

	p->ana_info.delta_mc_bonded_energy = tmp;
	}


    if(H5Pclose(plist_id) < 0){
	fprintf(stderr, "ERROR: HDF5 %s:%s:%d\n",__func__, __FILE__, __LINE__);
	return -1;
    }
    //Reopen file with out MPIIO
    if(H5Fclose(p->ana_info.file_id) < 0){
	fprintf(stderr, "ERROR: HDF5 %s:%s:%d\n",__func__, __FILE__, __LINE__);
	return -1;
	}

    if( p->info_MPI.sim_rank == 0)
	{
	file_id_tmp = H5Fopen(p->ana_info.filename, H5F_ACC_RDWR, H5P_DEFAULT);
	status = file_id_tmp;
	if( status < 0){
	    fprintf(stderr, "ERROR: HDF5 %s:%s:%d\n",__func__, __FILE__, __LINE__);
	    return -1;
	    }
	p->ana_info.file_id = file_id_tmp;
	}
    else
	p->ana_info.file_id = -1; //Invalid value, do not use.

#pragma acc update device(p->ana_info)
    return 0;
    }

int close_ana(struct Ana_Info*const a)
    {
    free(a->coord_filename);
    free(a->filename);

    if( a->file_id >= 0)
	{
	//Close the ana file
	if(H5Fclose(a->file_id) < 0){
	    fprintf(stderr, "ERROR: HDF5 %s:%s:%d\n",__func__, __FILE__, __LINE__);
	    return -1;
	    }
	}

    return 0;
    }
