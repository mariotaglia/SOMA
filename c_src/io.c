/* Copyright (C) 2016-2017 Ludwig Schneider
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
#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif//_OPENMP
#include "mesh.h"

int read_old_config(struct Phase * p, char *const filename)
    {
    FILE *inputfile;

    /* open input file */
    if ((inputfile = fopen(filename, "r")) == NULL) {

	fprintf(stderr, "Cannot open coord file %s\n", filename);
	return -1;

	}
    if( p->info_MPI.world_size != 1){//Make sure this is only called for single core processes.
	fprintf(stderr,"ERROR: convert can only be called with a single MPI-core.\n");
	return -11;
	}

    int iread, NX, NY, NZ, S, HYDRO, LEBD, TCHECK, TWRITE, VCHECK, ZA, ZB,
	ZC;
    unsigned int NpolyA, NmonoA, NpolyB, NmonoB, NpolyC, NmonoC;
    soma_scalar_t chiN, deltaN, kappaN, stiff, hN1, gN1, fN1, phia, phib, phic,
	D0, LY, LZ, time, F, f, xiN, gdot, forceN, kNcphi, dtcphi, Neq,
	Nav, Ntdgl, kT, dummy, dt;
    char *Cdummy;
    char cdummy[1024];
    unsigned int i, j;

    Cdummy = &cdummy[0];


    /* read first line of coord.dat */
    iread =
	fscanf(inputfile,
#if (SINGLE_PRECISION == 1)
	       "chiN = %g; deltaN = %g; kappaN = %g; stiff = %g; hN1 = %g; gN1 = %g; fN1 = %g; phia = %g; phib = %g; phic = %g; D0 = %g; LY = %g; LZ = %g; N = %d %d %d; S = %d;\n",
#else//SINGLE_PRECISION
	       "chiN = %lg; deltaN = %lg; kappaN = %lg; stiff = %lg; hN1 = %lg; gN1 = %lg; fN1 = %lg; phia = %lg; phib = %lg; phic = %lg; D0 = %lg; LY = %lg; LZ = %lg; N = %d %d %d; S = %d;\n",
#endif//SINGLE_PRECISION
	       &chiN, &deltaN, &kappaN, &stiff, &hN1, &gN1, &fN1, &phia,
	       &phib, &phic, &D0, &LY, &LZ, &NX, &NY, &NZ, &S);

    if (iread != 17) {
	fprintf(stderr,"Old io ERROR: %s:%d. (%d)\n",__FILE__,__LINE__,iread);
	return -1;
	}


    /* read second line of coord.dat */
    iread =
	fscanf(inputfile,
#if (SINGLE_PRECISION == 1)
	       "time =  %g; F = %g ( %g ); f = %g; xiN = %g; dt*N/xi = %g; gdot = %g; forceN = %g; kNcphi = %g; dtcphi = %g; Neq = %g; Nav = %g; Ntdgl = %g; kT = %g; HYDRO = %d; LEBD = %d; TCHECK = %d; TWRITE = %d; VCHECK = %d;\n",
#else//SINGLE_PRECISION
	       "time =  %lg; F = %lg ( %lg ); f = %lg; xiN = %lg; dt*N/xi = %lg; gdot = %lg; forceN = %lg; kNcphi = %lg; dtcphi = %lg; Neq = %lg; Nav = %lg; Ntdgl = %lg; kT = %lg; HYDRO = %d; LEBD = %d; TCHECK = %d; TWRITE = %d; VCHECK = %d;\n",
#endif//SINGLE_PRECISION
	       &time, &F, &dummy, &f, &xiN, &dt, &gdot, &forceN, &kNcphi,
	       &dtcphi, &Neq, &Nav, &Ntdgl, &kT, &HYDRO, &LEBD, &TCHECK,
	       &TWRITE, &VCHECK);
    if( iread != 19)
	{
	fprintf(stderr,"Old io ERROR: %s:%d. (%d)\n",__FILE__,__LINE__,iread);
	return -2;
	}

    /* read third line of coord.dat - free comment */
    Cdummy = fgets(Cdummy, 1024, inputfile);

    /* read fourth line of coord.dat */
    iread =
	fscanf(inputfile, "%u %u %d  %u %u %d  %u %u %d\n", &NpolyA,
	       &NmonoA, &ZA, &NpolyB, &NmonoB, &ZB, &NpolyC, &NmonoC, &ZC);
    if( iread != 9)
	{
	fprintf(stderr,"Old io ERROR: %s:%d. (%d)\n",__FILE__,__LINE__,iread);
	return -3;
	}


    p->reference_Nbeads = NmonoA;

    /* use the values we need */

    p->nx = NX;
    p->ny = NY;
    p->nz = NZ;
    p->n_cells = NX * NY * NZ;
    p->Lx = D0;
    p->Ly = LY;
    p->Lz = LZ;
    p->time = time;

    p->n_polymers = NpolyA;
    p->n_polymers_storage = p->n_polymers;
    assert(p->info_MPI.world_size == 1);	//Make sure this is only called for single core processes.
    p->n_polymers_global = p->n_polymers;
    p->n_types = 2;

    /* Assuming same diffusion constant for all monomers for this input file, but value can vary! */
    p->A = (soma_scalar_t*) malloc(p->n_types * sizeof(soma_scalar_t));
    if(p->A == NULL){
	fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
	return -4;
	}
    p->A[0] = dt / p->reference_Nbeads;
    p->A[1] = dt / p->reference_Nbeads;

    /* allocate XN interaction matrix */
    p->xn = (soma_scalar_t **) malloc(2 * sizeof(soma_scalar_t *));
    if(p->xn == NULL){
	fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
	return -4;
	}
    p->xn[0] = (soma_scalar_t *) malloc(2 * sizeof(soma_scalar_t));
    if(p->xn[0] == NULL){
	fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
	return -4;
	}
    p->xn[1] = (soma_scalar_t *) malloc(2 * sizeof(soma_scalar_t));
    if(p->xn[1] == NULL){
	fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
	return -4;
	}


    /* set parameters */
    p->xn[0][0] = kappaN;
    p->xn[1][1] = kappaN;
    p->xn[0][1] = deltaN;
    p->xn[1][0] = deltaN;


    //Create the polymer architecture list.
    //The old formate supports only 1 architecture -- a linear chain.
    p->n_poly_type= 1;
    p->poly_type_offset = (int*) malloc( p->n_poly_type*sizeof(int));
    if(p->poly_type_offset == NULL){
	fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
	return -4;
	}
    p->poly_type_offset[0] = 0;
    //Length = 1(N) + N(mono) + 4 (+1,-1,+1-1)
    p->poly_arch_length = 1 + NmonoA + 4;

    p->poly_arch = (uint32_t*) malloc( p->poly_arch_length * sizeof(uint32_t));
    if(p->poly_arch == NULL){
	fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
	return -4;
	}

    //Set up the architecture
    p->poly_arch[0] = NmonoA;

    p->poly_arch[ 0 + 1 ] = get_info_bl(NmonoA +4,0 );
    for(unsigned int i=1; i < NmonoA-1; i++)
	{
	p->poly_arch[ i + 1] = get_info_bl( NmonoA+3, 0);
	}
    p->poly_arch[NmonoA] = get_info_bl( NmonoA +1, 0);

    //Last monomer
    int end,offset,bond_type,info;
    bond_type = HARMONIC;
    end = 1;
    offset = -1;
    info = get_info(offset, bond_type, end);
    p->poly_arch[ NmonoA + 1 ] = info;
    //First monomer
    end = 1;
    offset = 1;
    info = get_info(offset, bond_type, end);
    p->poly_arch[ NmonoA + 2 ] = info;
    //Middle monomers
    end = 0;
    offset = -1;
    info = get_info(offset, bond_type, end);
    p->poly_arch[ NmonoA + 3 ] = info;

    end = 1;
    offset = 1;
    info = get_info(offset, bond_type, end);
    assert( NmonoA + 4 == p->poly_arch_length-1);
    p->poly_arch[ NmonoA + 4 ] = info;

    /* allocate space for polymers */
    p->polymers = (Polymer *) malloc(p->n_polymers_storage * sizeof(Polymer));
    if(p->polymers == NULL){
	fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
	return -4;
	}


    /* read in chains and set up particle data */
    for (i = 0; i < p->n_polymers; i++) {
	p->polymers[i].beads =
	    (Monomer *) malloc(NmonoA * sizeof(Monomer));
	if(p->polymers[i].beads == NULL){
	    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
	    return -4;
	    }
	p->polymers[i].msd_beads = (Monomer*) malloc( NmonoA*sizeof(Monomer) );
	if( p->polymers[i].msd_beads == NULL){
	    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	    return -1;
	    }


	p->polymers[i].type = 0;
	/* read the monomers of this chain */
	for (j = 0; j < NmonoA; j++) {
	    unsigned int tmp;
	    iread =fscanf(inputfile,
#if (SINGLE_PRECISION == 1)
			  "%g %g %g %u\n",
#else//SINGLE_PRECISION
			  "%lg %lg %lg %u\n",
#endif//SINGLE_PRECISION
			  &p->polymers[i].beads[j].x,
			  &p->polymers[i].beads[j].y,
			  &p->polymers[i].beads[j].z, &tmp);
	    if(iread != 4)
		{
		fprintf(stderr,"Old io ERROR: %s:%d. (%d)\n",__FILE__,__LINE__,iread);
		return -6;
		}

	    tmp -= 1;		//First type in conf is 1. First type in program is 0

	    //! \warning Assumed is that all polymers have the same
	    //! type. This includes, that the particle type
	    //! distribution is for all polymers identical. (The last
	    //! specified in the input format.)
	    uint32_t *const type_info = &(p->poly_arch[ p->poly_type_offset[ p->polymers[i].type] + j + 1  ]);
	    const int offset_bl = get_bondlist_offset(*type_info);
	    // Overwrite the type info every single step.
	    *type_info = get_info_bl(offset_bl, tmp);
	    }
	}
    p->bead_data_read = true;
    fclose(inputfile);

    p->area51 = NULL;
    p->external_field_unified = NULL;
    p->umbrella_field = NULL;
    p->hamiltonian = SCMF0;
    p->k_umbrella = (soma_scalar_t*const)malloc(p->n_types*sizeof(soma_scalar_t));
    if( p->k_umbrella == NULL)
	{
	fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
	return -4;
	}
    memset(p->k_umbrella,0,p->n_types*sizeof(soma_scalar_t));

    p->harmonic_normb_variable_scale = 1;
    p->cm_a = NULL; // Deactivate CM movement with the old file format.

    return 0;
    }

int read_old_geometry(struct Phase*p,const char*filename)
    {
    /* check if also an old style geometry needs to be included */
    FILE * geofile = fopen(filename, "rb");
    if (geofile != NULL){
	p->external_field_unified = (soma_scalar_t *) malloc (p->n_types*p->n_cells * sizeof(soma_scalar_t*));
	if (p->external_field_unified == NULL) {
	    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__,
		    __LINE__);
	    return -1;
	    }

	p->area51 = (uint8_t *) malloc( p->n_cells * sizeof(uint8_t ));
	if (p->area51 == NULL) {
	    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__,
		    __LINE__);
	    return -1;
	    }
	int x,y,z;
	soma_scalar_t  wall, field, dielectric, compute_efield, fixed_potential;;

	// in the original definition external fields are symmetric to the two species
	for (unsigned int k = 0; k < p->n_cells; k++){
#if (SINGLE_PRECISION == 1)
	    const int iread = fscanf(geofile,"%d %d %d %g %g %g %g %g", &x, &y, &z, &wall, &field, &dielectric, &compute_efield, &fixed_potential);
#else//SINGLE_PRECISION
	    const int iread = fscanf(geofile,"%d %d %d %lg %lg %lg %lg %lg", &x, &y, &z, &wall, &field, &dielectric, &compute_efield, &fixed_potential);
#endif//SINGLE_PRECISION
	    if (iread != 8) {
		fprintf(stderr,"Geometry io ERROR: %s:%d. (%d)\n",__FILE__,__LINE__,iread);
		return -1;
		}
	    if (wall > 0 ) { p->area51[k] = 0; } else {p->area51[k]=1;}
	    p->external_field_unified[cell_to_index_unified(p,k,0)] = field;
	    p->external_field_unified[cell_to_index_unified(p,k,1)] = -field;
	    }
	fclose(geofile);
	}
    else {
	fprintf(stderr,"WARNING: Passed geometry file %s unreadable. Geometry will be ignored.\n",filename);
	p->area51 = NULL;
	p->external_field_unified = NULL;
	}
    return 0;
    }


#include <mpi.h>
#include "hdf5.h"

int write_hdf5(const hsize_t ndims,const hsize_t*const dims,const hid_t file_id,
	       const char*const name,const hid_t file_type,const hid_t mem_type,
	       const hid_t plist_id,const void*const data)
    {
    herr_t status;
    const hid_t dataspace = H5Screate_simple(ndims, dims, NULL);
    HDF5_ERROR_CHECK2(dataspace,name);
    const hid_t dataset = H5Dcreate2(file_id, name , file_type,dataspace,
				     H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    HDF5_ERROR_CHECK2(dataset,name);

    status = H5Dwrite(dataset, mem_type, H5S_ALL, H5S_ALL,
		      plist_id, data);
    HDF5_ERROR_CHECK2(status,name);

    status = H5Sclose(dataspace);
    HDF5_ERROR_CHECK2(status,name);
    status = H5Dclose(dataset);
    HDF5_ERROR_CHECK2(status,name);
    return status;
    }


int write_config_hdf5(const struct Phase * const p, const char *filename)
    {
    // Copy polymer data from device to host
    update_self_phase(p);

    herr_t status;

    //Set up file access property list with parallel I/O access
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    if (p->info_MPI.sim_size > 1)
	H5Pset_fapl_mpio(plist_id, p->info_MPI.SOMA_comm_sim,
			 MPI_INFO_NULL);

    //Create a new h5 file and overwrite the content.
    hid_t file_id =
	H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    if (p->info_MPI.sim_size > 1)
	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    //Create a group for all parameter of the simulation
    hid_t parameter_group =
	H5Gcreate2(file_id, "/parameter", H5P_DEFAULT, H5P_DEFAULT,
		   H5P_DEFAULT);

    //Write number of polymers
    const hsize_t one = 1;
    status = write_hdf5(1,&one,file_id,"/parameter/n_polymers",H5T_STD_U64LE,H5T_NATIVE_UINT64,plist_id,&(p->n_polymers_global));
    HDF5_ERROR_CHECK2(status,"/parameter/n_polymers");

    status = write_hdf5(1,&one,file_id,"/parameter/reference_Nbeads",H5T_STD_U32LE,H5T_NATIVE_UINT,plist_id,&(p->reference_Nbeads));
    HDF5_ERROR_CHECK2(status,"/parameter/reference_Nbeads");

    status = write_hdf5(1,&one,file_id,"/parameter/n_types",H5T_STD_U32LE,H5T_NATIVE_UINT,plist_id,&(p->n_types));
    HDF5_ERROR_CHECK2(status,"/parameter/reference_Nbeads");

    //Number of types
    hsize_t xn_dim[2] = { p->n_types, p->n_types };
    soma_scalar_t *const tmp_xn =
	(soma_scalar_t *const) malloc(p->n_types * p->n_types * sizeof(soma_scalar_t));
    if (tmp_xn == NULL) {
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
	}
    for (unsigned int i = 0; i < p->n_types; i++)
	for (unsigned int j = 0; j < p->n_types; j++) {
	    tmp_xn[i * p->n_types + j] = p->xn[i][j];
	    }
    status = write_hdf5(2,xn_dim,file_id,"/parameter/xn",H5T_SOMA_FILE_SCALAR,H5T_SOMA_NATIVE_SCALAR,plist_id,tmp_xn);
    HDF5_ERROR_CHECK2(status,"/parameter/xn");
    free(tmp_xn);

    //A data
    hsize_t a_size = p->n_types;
    status = write_hdf5(1,&a_size,file_id,"/parameter/A",H5T_SOMA_FILE_SCALAR,H5T_SOMA_NATIVE_SCALAR,plist_id,p->A);
    HDF5_ERROR_CHECK2(status,"/parameter/A");

    //time
    status = write_hdf5(1,&one,file_id,"/parameter/time",H5T_STD_U32LE,H5T_NATIVE_UINT,plist_id,&(p->time));
    HDF5_ERROR_CHECK2(status,"/parameter/time");

    status = write_hdf5(1,&one,file_id,"/parameter/hamiltonian",
			H5T_STD_I32LE,H5T_NATIVE_INT,plist_id,&(p->hamiltonian));
    HDF5_ERROR_CHECK2(status,"parameter/hamiltonian");

    a_size = p->n_types;
    status = write_hdf5(1, &a_size, file_id, "/parameter/k_umbrella", H5T_SOMA_FILE_SCALAR, H5T_SOMA_NATIVE_SCALAR, plist_id, p->k_umbrella);
    HDF5_ERROR_CHECK2(status, "parameter/k_umbrella");

    //Nx Ny Nz
    hsize_t three = 3;
    unsigned int nxyz[3] = { p->nx, p->ny, p->nz };
    status = write_hdf5(1,&three,file_id,"/parameter/nxyz",H5T_STD_U32LE,H5T_NATIVE_UINT,plist_id,nxyz);
    HDF5_ERROR_CHECK2(status,"/parameter/nxyz");

    //Lx Ly Lz
    soma_scalar_t lxyz[3] = { p->Lx, p->Ly, p->Lz };
    status = write_hdf5(1,&three,file_id,"/parameter/lxyz",H5T_SOMA_FILE_SCALAR,H5T_SOMA_NATIVE_SCALAR,plist_id,lxyz);
    HDF5_ERROR_CHECK2(status,"/parameter/lxyz");

    //p->harmonic_normb_variable_scale
    status = write_hdf5(1,&one,file_id,"/parameter/harmonic_normb_variable_scale",H5T_IEEE_F64LE,H5T_SOMA_NATIVE_SCALAR,plist_id,&(p->harmonic_normb_variable_scale));
    HDF5_ERROR_CHECK2(status,"/parameter/harmonic_normb_variable_scale");

    //Polymer architecture
    //Number of polymer type
    status = write_hdf5(1,&one,file_id,"/parameter/n_poly_type",H5T_STD_U32LE,
			H5T_NATIVE_UINT,plist_id,&(p->n_poly_type));
    HDF5_ERROR_CHECK2(status,"n_poly_types");
    //poly_type_offset
    hsize_t n_poly_type = p->n_poly_type;
    status = write_hdf5(1,&n_poly_type,file_id,"/parameter/poly_type_offset",H5T_STD_I32LE,
			H5T_NATIVE_INT,plist_id,p->poly_type_offset);
    HDF5_ERROR_CHECK2(status,"poly_type_offset");
    //poly_arch_length
    status = write_hdf5(1,&one,file_id,"/parameter/poly_arch_length",H5T_STD_U32LE,
			H5T_NATIVE_UINT,plist_id,&(p->poly_arch_length));
    HDF5_ERROR_CHECK2(status,"poly_arch_length");

    //poly_arch
    const hsize_t arch_length = p->poly_arch_length;
    //Warning this is a hidden cast from UINT32 to INT32
    status = write_hdf5(1,&arch_length,file_id,"/parameter/poly_arch",H5T_STD_I32LE,
			H5T_NATIVE_INT,plist_id,p->poly_arch);
    HDF5_ERROR_CHECK2(status,"poly_arch");

    if(p->cm_a != NULL)
    	{
	status = write_hdf5(1,&n_poly_type,file_id,"/parameter/cm_a",H5T_SOMA_FILE_SCALAR,
			    H5T_SOMA_NATIVE_SCALAR,plist_id,p->cm_a);
	HDF5_ERROR_CHECK2(status,"/parameter/cm_a");
	}


    //Close parameter group
    if ((status = H5Gclose(parameter_group)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		p->info_MPI.world_rank, __FILE__, __LINE__, status);
	return status;
	}

    //Write out the polymers
    //Determine the offset, for each process.
    uint64_t n_polymer_offset;
    //Cast for MPI_Scan, since some openmpi impl. need a non-const. version.
    MPI_Scan((uint64_t *) &(p->n_polymers), &n_polymer_offset, 1,
	     MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);
    n_polymer_offset -= p->n_polymers;

    unsigned int *const poly_type =
	(unsigned int *const) malloc(p->n_polymers_storage * sizeof(unsigned int));
    if (poly_type == NULL) {
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
	}

    for (uint64_t i = 0; i < p->n_polymers; i++)
	poly_type[i] = p->polymers[i].type;

    hsize_t hsize_polymers_global = p->n_polymers_global;
    if( hsize_polymers_global != p->n_polymers_global)
	{
	fprintf(stderr,"\n\n ERROR: %s:%d unhandled uint64 overflow.\n\n\n",__FILE__,__LINE__);
	return -1;
	}
    hid_t poly_type_dataspace =
	H5Screate_simple(1, &(hsize_polymers_global), NULL);
    hsize_t hsize_polymers = p->n_polymers;
    hid_t poly_type_memspace = H5Screate_simple(1, &(hsize_polymers), NULL);

    hid_t poly_type_dataset =
	H5Dcreate2(file_id, "/poly_type", H5T_STD_U32LE, poly_type_dataspace,
		   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    poly_type_dataspace = H5Dget_space(poly_type_dataset);
    hsize_t hsize_polymers_offset = n_polymer_offset;
    H5Sselect_hyperslab(poly_type_dataspace, H5S_SELECT_SET,
			&(hsize_polymers_offset), NULL, &(hsize_polymers),
			NULL);

    if ((status =
	 H5Dwrite(poly_type_dataset, H5T_NATIVE_UINT, poly_type_memspace,
		  poly_type_dataspace, plist_id, poly_type)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		p->info_MPI.world_rank, __FILE__, __LINE__, status);
	return status;
	}

    free(poly_type);
    if ((status = H5Sclose(poly_type_dataspace)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		p->info_MPI.world_rank, __FILE__, __LINE__, status);
	return status;
	}

    if ((status = H5Sclose(poly_type_memspace)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		p->info_MPI.world_rank, __FILE__, __LINE__, status);
	return status;
	}

    if ((status = H5Dclose(poly_type_dataset)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		p->info_MPI.world_rank, __FILE__, __LINE__, status);
	return status;
	}

    unsigned int max_n_beads = 0;
    for (unsigned int i = 0; i < p->n_poly_type; i++)
	{
	const unsigned int N = p->poly_arch[ p->poly_type_offset[i] ];
	if (max_n_beads < N)
	    max_n_beads = N;
	}

    hsize_t hsize_beads_dataspace[2] =
	{ p->n_polymers_global, max_n_beads };
    hid_t beads_dataspace =
	H5Screate_simple(2, hsize_beads_dataspace, NULL);
    hsize_t hsize_beads_memspace[2] = { p->n_polymers, max_n_beads };
    hid_t beads_memspace = H5Screate_simple(2, hsize_beads_memspace, NULL);

    hid_t monomer_filetype = get_monomer_filetype();

    hid_t beads_dataset =
	H5Dcreate2(file_id, "/beads", monomer_filetype, beads_dataspace,
		   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    beads_dataspace = H5Dget_space(beads_dataset);
    hsize_t hsize_beads_offset[2] = { n_polymer_offset, 0 };
    H5Sselect_hyperslab(beads_dataspace, H5S_SELECT_SET,
			hsize_beads_offset, NULL, hsize_beads_memspace,
			NULL);

    hid_t monomer_memtype = get_monomer_memtype();

    //Monomer monomer_data[p->n_polymers][max_n_beads];
    Monomer *const monomer_data =
	(Monomer * const) malloc(p->n_polymers * max_n_beads *
				 sizeof(Monomer));
    if (monomer_data == NULL) {
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
	}
    memset(monomer_data,0,p->n_polymers * max_n_beads * sizeof(Monomer));

    for (uint64_t i = 0; i < p->n_polymers; i++)
	{
	const unsigned int N = p->poly_arch[ p->poly_type_offset[p->polymers[i].type] ];
	memcpy(monomer_data + i * max_n_beads, p->polymers[i].beads,
	       N * sizeof(Monomer));
	}
    /* for(unsigned int j=0; j < p->polymers[i].N; j++) */
    /*    { */
    /*      assert( j < max_n_beads ); */
    /*      monomer_data[i*max_n_beads+j] = p->polymers[i].beads[j]; */
    /*    } */

    if ((status =
	 H5Dwrite(beads_dataset, monomer_memtype, beads_memspace,
		  beads_dataspace, plist_id, monomer_data)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		p->info_MPI.world_rank, __FILE__, __LINE__, status);
	return status;
	}

    free(monomer_data);
    if ((status = H5Tclose(monomer_memtype)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		p->info_MPI.world_rank, __FILE__, __LINE__, status);
	return status;
	}

    if ((status = H5Tclose(monomer_filetype)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		p->info_MPI.world_rank, __FILE__, __LINE__, status);
	return status;
	}

    if ((status = H5Sclose(beads_dataspace)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		p->info_MPI.world_rank, __FILE__, __LINE__, status);
	return status;
	}

    if ((status = H5Sclose(beads_memspace)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		p->info_MPI.world_rank, __FILE__, __LINE__, status);
	return status;
	}

    if ((status = H5Dclose(beads_dataset)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		p->info_MPI.world_rank, __FILE__, __LINE__, status);
	return status;
	}



    if(p->area51)
	{
	const hsize_t hsize_ncells = p->n_cells;
	hid_t n_cells_dataspace =
	    H5Screate_simple(1, &(hsize_ncells), NULL);

	hid_t n_cells_dataset =
	    H5Dcreate2(file_id, "/area51", H5T_STD_U8LE, n_cells_dataspace,
		       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	n_cells_dataspace = H5Dget_space(n_cells_dataset);
	hsize_t area51_offset = p->info_MPI.sim_rank * ( hsize_ncells/p->info_MPI.sim_size);
	if( (unsigned int)p->info_MPI.sim_rank < ( hsize_ncells) % p->info_MPI.sim_size )
	    area51_offset += p->info_MPI.sim_rank;
	else
	    area51_offset += hsize_ncells % p->info_MPI.sim_size;

	hsize_t area51_local = hsize_ncells/p->info_MPI.sim_size;
	if( (unsigned int)p->info_MPI.sim_rank < hsize_ncells % p->info_MPI.sim_size)
	    area51_local += 1;
	hid_t n_cells_memspace = H5Screate_simple(1, &(area51_local), NULL);
	H5Sselect_hyperslab(n_cells_dataspace, H5S_SELECT_SET,
			    &(area51_offset), NULL, &(area51_local),NULL);

	if ((status =
	     H5Dwrite(n_cells_dataset, H5T_NATIVE_UINT8, n_cells_memspace,
		      n_cells_dataspace, plist_id, p->area51 + area51_offset)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
	    return status;
	    }

	if ((status = H5Sclose(n_cells_dataspace)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
	    return status;
	    }

	if ((status = H5Sclose(n_cells_memspace)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
	    return status;
	    }

	if ((status = H5Dclose(n_cells_dataset)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
	    return status;
	    }

	}

    if(p->external_field_unified)
	{
	const hsize_t hsize_ncells = p->n_cells*p->n_types;
	hid_t n_cells_dataspace =
	    H5Screate_simple(1, &(hsize_ncells), NULL);

	hid_t n_cells_dataset =
	    H5Dcreate2(file_id, "/external_field", H5T_SOMA_FILE_SCALAR, n_cells_dataspace,
		       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	n_cells_dataspace = H5Dget_space(n_cells_dataset);
	hsize_t ex_field_offset = p->info_MPI.sim_rank * ( hsize_ncells/p->info_MPI.sim_size);
	if( (unsigned int)p->info_MPI.sim_rank < ( hsize_ncells) % p->info_MPI.sim_size )
	    ex_field_offset += p->info_MPI.sim_rank;
	else
	    ex_field_offset += hsize_ncells % p->info_MPI.sim_size;

	hsize_t ex_field_local = hsize_ncells/p->info_MPI.sim_size;
	if( (unsigned int)p->info_MPI.sim_rank < hsize_ncells % p->info_MPI.sim_size)
	    ex_field_local += 1;
	hid_t n_cells_memspace = H5Screate_simple(1, &(ex_field_local), NULL);
	H5Sselect_hyperslab(n_cells_dataspace, H5S_SELECT_SET,
			    &(ex_field_offset), NULL, &(ex_field_local),NULL);

	if ((status =
	     H5Dwrite(n_cells_dataset, H5T_SOMA_NATIVE_SCALAR, n_cells_memspace,
		      n_cells_dataspace, plist_id, p->external_field_unified+ex_field_offset)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
	    return status;
	    }

	if ((status = H5Sclose(n_cells_dataspace)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
	    return status;
	    }

	if ((status = H5Sclose(n_cells_memspace)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
	    return status;
	    }

	if ((status = H5Dclose(n_cells_dataset)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
	    return status;
	    }

	}
    if(p->umbrella_field)
	{
	const hsize_t hsize_ncells = p->n_cells*p->n_types;
	hid_t n_cells_dataspace =
	    H5Screate_simple(1, &(hsize_ncells), NULL);

	hid_t n_cells_dataset =
	    H5Dcreate2(file_id, "/umbrella_field", H5T_SOMA_FILE_SCALAR, n_cells_dataspace,
		       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	n_cells_dataspace = H5Dget_space(n_cells_dataset);
	hsize_t umbrella_field_offset = p->info_MPI.sim_rank * ( hsize_ncells/p->info_MPI.sim_size);
	if( (unsigned int)p->info_MPI.sim_rank < ( hsize_ncells) % p->info_MPI.sim_size )
	    umbrella_field_offset += p->info_MPI.sim_rank;
	else
	    umbrella_field_offset += hsize_ncells % p->info_MPI.sim_size;

	hsize_t umbrella_field_local = hsize_ncells/p->info_MPI.sim_size;
	if( (unsigned int)p->info_MPI.sim_rank < hsize_ncells % p->info_MPI.sim_size)
	    umbrella_field_local += 1;
	hid_t n_cells_memspace = H5Screate_simple(1, &(umbrella_field_local), NULL);
	H5Sselect_hyperslab(n_cells_dataspace, H5S_SELECT_SET,
			    &(umbrella_field_offset), NULL, &(umbrella_field_local),NULL);

	if ((status =
	     H5Dwrite(n_cells_dataset, H5T_SOMA_NATIVE_SCALAR, n_cells_memspace,
		      n_cells_dataspace, plist_id, p->umbrella_field+umbrella_field_offset)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
	    return status;
	    }

	if ((status = H5Sclose(n_cells_dataspace)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
	    return status;
	    }

	if ((status = H5Sclose(n_cells_memspace)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
	    return status;
	    }

	if ((status = H5Dclose(n_cells_dataset)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
	    return status;
	    }

	}


    H5Pclose(plist_id);
    if ((status = H5Fclose(file_id)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		p->info_MPI.world_rank, __FILE__, __LINE__, status);
	return status;
	}


    return status;
    }

int read_hdf5(const hid_t file_id,const char*const name,const hid_t mem_type,const hid_t plist_id,void*const data)
    {
    herr_t status;
    const hid_t dataset = H5Dopen2(file_id, name, H5P_DEFAULT);
    HDF5_ERROR_CHECK2(dataset,name);
    status = H5Dread(dataset, mem_type, H5S_ALL, H5S_ALL, plist_id, data);
    HDF5_ERROR_CHECK2(status,name);
    status = H5Dclose(dataset);
    HDF5_ERROR_CHECK2(status,name);
    return status;
    }

int read_config_hdf5(struct Phase * const p, const char *filename)
    {
    herr_t status;
    //Set up file access property list with parallel I/O access
    hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
    if (p->info_MPI.sim_size > 1)
	H5Pset_fapl_mpio(plist_id, p->info_MPI.SOMA_comm_sim,
			 MPI_INFO_NULL);

    //Create a new h5 file and overwrite the content.
    hid_t file_id = H5Fopen(filename, H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    if (p->info_MPI.sim_size > 1)
	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    status = read_hdf5(file_id,"/parameter/n_polymers",H5T_NATIVE_UINT64,plist_id,&(p->n_polymers_global));
    HDF5_ERROR_CHECK2(status,"/parameter/n_polymers");


    //Distribute the polymers to the different cores.
    uint64_t n_polymers = p->n_polymers_global / p->info_MPI.sim_size;
    if ((unsigned int) p->info_MPI.sim_rank <
	p->n_polymers_global % p->info_MPI.sim_size)
	n_polymers += 1;
    p->n_polymers = n_polymers;
    p->n_polymers_storage = p->n_polymers;
    uint64_t n_polymer_offset;
    //Cast for MPI_Scan, since some openmpi impl. need a non-const. version.
    MPI_Scan((uint64_t *) &(p->n_polymers), &n_polymer_offset, 1,
	     MPI_UINT64_T, MPI_SUM, p->info_MPI.SOMA_comm_sim);
    n_polymer_offset -= p->n_polymers;

    if (p->info_MPI.sim_rank == p->info_MPI.sim_size - 1)
	assert(p->n_polymers + n_polymer_offset == p->n_polymers_global);

    status = read_hdf5(file_id,"/parameter/reference_Nbeads",H5T_NATIVE_UINT,plist_id,&(p->reference_Nbeads));
    HDF5_ERROR_CHECK2(status,"/parameter/reference_Nbeads");

    status = read_hdf5(file_id,"/parameter/n_types",H5T_NATIVE_UINT,plist_id,&(p->n_types));
    HDF5_ERROR_CHECK2(status,"/parameter/n_types");


    p->xn = (soma_scalar_t **) malloc(p->n_types * sizeof(soma_scalar_t *));
    if (p->xn == NULL) {
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
	}
    for (unsigned int i = 0; i < p->n_types; i++) {
	p->xn[i] = (soma_scalar_t *) malloc(p->n_types * sizeof(soma_scalar_t));
	if (p->xn[i] == NULL) {
	    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	    return -1;
	    }
	}

    //soma_scalar_t tmp_xn[p->n_types][p->n_types];
    soma_scalar_t *const tmp_xn =
	(soma_scalar_t *const) malloc(p->n_types * p->n_types * sizeof(soma_scalar_t));
    if (tmp_xn == NULL) {
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
	}

    status = read_hdf5(file_id,"/parameter/xn",H5T_SOMA_NATIVE_SCALAR,plist_id,tmp_xn);
    HDF5_ERROR_CHECK2(status,"/parameter/xn");

    for (unsigned int i = 0; i < p->n_types; i++)
	for (unsigned int j = 0; j < p->n_types; j++)
	    p->xn[i][j] = tmp_xn[i * p->n_types + j];
    free(tmp_xn);

    //A array for the diffusivity of the particles
    p->A = (soma_scalar_t *const) malloc(p->n_types * sizeof(soma_scalar_t));
    if (p->A == NULL) {
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
	}

    // read p->A
    status = read_hdf5(file_id,"/parameter/A",H5T_SOMA_NATIVE_SCALAR,plist_id,p->A);
    HDF5_ERROR_CHECK2(status,"/parameter/A");

    // read p->time
    status = read_hdf5(file_id,"/parameter/time",H5T_NATIVE_UINT,plist_id,&(p->time));
    HDF5_ERROR_CHECK2(status,"/parameter/time");

    // read nx ny nz
    p->hamiltonian = SCMF0;
    //Don't break old configurations.
    if(H5Lexists(file_id,"/parameter/hamiltonian",H5P_DEFAULT) > 0)
    	{
	status = read_hdf5(file_id,"/parameter/hamiltonian",H5T_NATIVE_INT,
			   plist_id,&(p->hamiltonian));
	HDF5_ERROR_CHECK2(status,"/parameter/hamiltonian");
	}

    p->k_umbrella = (soma_scalar_t*const)malloc(p->n_types * sizeof(soma_scalar_t));
    if (p->k_umbrella == NULL) {
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
	}
    memset(p->k_umbrella,0,p->n_types*sizeof(soma_scalar_t)); //Default 0
    if(H5Lexists(file_id,"/parameter/k_umbrella",H5P_DEFAULT) > 0)
	{
	status = read_hdf5(file_id, "/parameter/k_umbrella", H5T_SOMA_NATIVE_SCALAR, plist_id, p->k_umbrella);
	HDF5_ERROR_CHECK2(status, "parameter/k_umbrella");
	}

    unsigned int nxyz[3];
    status = read_hdf5(file_id,"/parameter/nxyz",H5T_NATIVE_UINT,plist_id,nxyz);
    HDF5_ERROR_CHECK2(status,"/parameter/nxyz");
    p->nx = nxyz[0];
    p->ny = nxyz[1];
    p->nz = nxyz[2];

    // read lx ly lz
    soma_scalar_t lxyz[3];
    status = read_hdf5(file_id,"/parameter/lxyz",H5T_SOMA_NATIVE_SCALAR,plist_id,lxyz);
    HDF5_ERROR_CHECK2(status, "/parameter/lxyz");
    p->Lx = lxyz[0];
    p->Ly = lxyz[1];
    p->Lz = lxyz[2];

    //Read in the polymer architectures.
    //Number of polymer type
    status = read_hdf5(file_id,"/parameter/n_poly_type",H5T_NATIVE_UINT,plist_id,&(p->n_poly_type));
    HDF5_ERROR_CHECK2(status,"n_poly_type");
    //poly_type_offset
    p->poly_type_offset = (int*)malloc(p->n_poly_type * sizeof(int));
    if(p->poly_type_offset == NULL){
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
	}
    status = read_hdf5(file_id,"/parameter/poly_type_offset",H5T_NATIVE_INT,plist_id,p->poly_type_offset);
    HDF5_ERROR_CHECK2(status,"poly_type_offset");
    //poly_arch_length
    status = read_hdf5(file_id,"/parameter/poly_arch_length",H5T_NATIVE_UINT,plist_id,&(p->poly_arch_length));
    HDF5_ERROR_CHECK2(status,"poly_arch_length");
    //poly_arch
    p->poly_arch = (uint32_t*) malloc( p->poly_arch_length * sizeof(uint32_t));
    if(p->poly_arch == NULL){
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
	}
    status = read_hdf5(file_id,"/parameter/poly_arch",H5T_NATIVE_INT32,plist_id,p->poly_arch);
    HDF5_ERROR_CHECK2(status,"poly_arch");

    p->cm_a = NULL; //Default: deactivated.
    //If mobility is specified write it out.
    if(H5Lexists(file_id,"/parameter/cm_a",H5P_DEFAULT) > 0)
    	{
	p->cm_a = (soma_scalar_t*) malloc( p->n_poly_type*sizeof( soma_scalar_t) );
	if(p->cm_a == NULL)
	    {
	    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	    return -1;
	    }
	status = read_hdf5(file_id,"/parameter/cm_a",H5T_SOMA_NATIVE_SCALAR,
			   plist_id,p->cm_a);
	HDF5_ERROR_CHECK2(status,"/parameter/cm_a");
	}

    // read harmonic_normb_variable_scale
    if(H5Lexists(file_id,"/parameter/harmonic_normb_variable_scale",H5P_DEFAULT) > 0)
	{
	status = read_hdf5(file_id,"/parameter/harmonic_normb_variable_scale",H5T_SOMA_NATIVE_SCALAR,plist_id, &(p->harmonic_normb_variable_scale));
	HDF5_ERROR_CHECK2(status,"/parameter/harmonic_normb_variable_scale");
	}
    else
	{
	if( p->info_MPI.sim_rank == 0)
	    {
	    if( get_number_bond_type(p,HARMONICVARIABLESCALE) != 0 )
		fprintf(stderr,"WARNING: The poly_arch contains HARMONICVARIABLESCALE Bond, but no corresponding value is set. Using 1 instead.\n");
	    }
	p->harmonic_normb_variable_scale = 1.;
	}

    //Allocate memory for almost everything of PHASE
    //Polymers (only the local copies).
    p->polymers = (Polymer *) malloc(p->n_polymers_storage * sizeof(Polymer));
    if (p->polymers == NULL) {
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
	}
    //Read in polymer related data, only the locally needed info
    unsigned int *const poly_type =
	(unsigned int *const) malloc(p->n_polymers_storage * sizeof(unsigned int));
    if (poly_type == NULL) {
	fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	return -1;
	}

    hsize_t hsize_polymers = p->n_polymers;
    hsize_t hsize_polymers_offset = n_polymer_offset;
    hid_t poly_type_memspace = H5Screate_simple(1, &(hsize_polymers), NULL);
    hid_t poly_type_dataset = H5Dopen2(file_id, "/poly_type", H5P_DEFAULT);
    hid_t poly_type_dataspace = H5Dget_space(poly_type_dataset);
    H5Sselect_hyperslab(poly_type_dataspace, H5S_SELECT_SET,
			&(hsize_polymers_offset), NULL, &(hsize_polymers),
			NULL);

    if ((status =
	 H5Dread(poly_type_dataset, H5T_NATIVE_UINT, poly_type_memspace,
		 poly_type_dataspace, plist_id, poly_type)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		p->info_MPI.world_rank, __FILE__, __LINE__, status);
	return status;
	}

    if ((status = H5Sclose(poly_type_memspace)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		p->info_MPI.world_rank, __FILE__, __LINE__, status);
	return status;
	}

    if ((status = H5Dclose(poly_type_dataset)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		p->info_MPI.world_rank, __FILE__, __LINE__, status);
	return status;
	}

    if ((status = H5Sclose(poly_type_dataspace)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		p->info_MPI.world_rank, __FILE__, __LINE__, status);
	return status;
	}

    for (uint64_t i = 0; i < p->n_polymers; i++) {
	p->polymers[i].type = poly_type[i];
	}
    free(poly_type);


    //Allocate space for monomers
    for (uint64_t i = 0; i < p->n_polymers; i++)
	{
	const unsigned int N = p->poly_arch[ p->poly_type_offset[ p->polymers[i].type ] ];
	p->polymers[i].beads = (Monomer*) malloc( N*sizeof(Monomer) );
	if( p->polymers[i].beads == NULL){
	    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	    return -1;
	    }
	p->polymers[i].msd_beads = (Monomer*) malloc( N*sizeof(Monomer) );
	if( p->polymers[i].msd_beads == NULL){
	    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	    return -1;
	    }
	}


    if(H5Lexists(file_id,"/beads",H5P_DEFAULT) > 0)
    	{


	//Get the data for the polymers
	//get the dimension
	unsigned int max_n_beads = 0;
	for (unsigned int i = 0; i < p->n_poly_type; i++)
	    {
	    const unsigned int N = p->poly_arch[ p->poly_type_offset[i] ];
	    if (max_n_beads < N)
		max_n_beads = N;
	    }

	Monomer *const monomer_data =
	    (Monomer * const) malloc(p->n_polymers * max_n_beads *
				     sizeof(Monomer));
	if (monomer_data == NULL) {
	    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
	    return -1;
	    }

	hsize_t hsize_beads_memspace[2] = { p->n_polymers, max_n_beads };
	hid_t beads_memspace = H5Screate_simple(2, hsize_beads_memspace, NULL);
	hid_t beads_dataset = H5Dopen2(file_id, "/beads", H5P_DEFAULT);
	hid_t beads_dataspace = H5Dget_space(beads_dataset);
	hsize_t hsize_beads_offset[2] = { n_polymer_offset, 0 };
	H5Sselect_hyperslab(beads_dataspace, H5S_SELECT_SET,
			    hsize_beads_offset, NULL, hsize_beads_memspace,
			    NULL);

	hid_t monomer_memtype = get_monomer_memtype();

	if ((status =
	     H5Dread(beads_dataset, monomer_memtype, beads_memspace,
		     beads_dataspace, plist_id, monomer_data)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
	    return status;
	    }

	for (uint64_t i = 0; i < p->n_polymers; i++)
	    {
	    const unsigned int N = p->poly_arch[ p->poly_type_offset[ p->polymers[i].type ] ];
	    memcpy(p->polymers[i].beads, monomer_data + i * max_n_beads,
		   N * sizeof(Monomer));
	    memcpy(p->polymers[i].msd_beads, monomer_data + i * max_n_beads,
		   N * sizeof(Monomer));
	    }
	free(monomer_data);

	if ((status = H5Tclose(monomer_memtype)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
	    return status;
	    }

	if ((status = H5Sclose(beads_dataspace)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
	    return status;
	    }

	if ((status = H5Sclose(beads_memspace)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
	    return status;
	    }

	if ((status = H5Dclose(beads_dataset)) < 0) {
	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
	    return status;
	    }
	p->bead_data_read = true;
	}
    else
	p->bead_data_read = false;

    p->area51 = NULL;
    p->external_field_unified = NULL;
    p->umbrella_field = NULL;
    p->n_cells = p->nx*p->ny*p->nz;

    if(H5Lexists(file_id,"/area51",H5P_DEFAULT) > 0)
    	{
    	const hsize_t hsize_ncells = p->n_cells;
	hid_t n_cells_memspace = H5Screate_simple(1, &(hsize_ncells), NULL);

    	hid_t n_cells_dataset = H5Dopen2(file_id, "/area51", H5P_DEFAULT);
    	hid_t n_cells_dataspace = H5Dget_space(n_cells_dataset);
	hsize_t dims;
	status = H5Sget_simple_extent_dims(n_cells_dataspace,&dims,NULL);
	HDF5_ERROR_CHECK(status);
	assert(dims == hsize_ncells);

	p->area51 = (uint8_t *) malloc( p->n_cells * sizeof(uint8_t ));
	if (p->area51 == NULL) {
	    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__,
		    __LINE__);
	    return -1;
	    }

    	if ((status =
    	     H5Dread(n_cells_dataset, H5T_STD_U8LE, n_cells_memspace,
    		     n_cells_dataspace, plist_id, p->area51)) < 0) {
    	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
    		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
    	    return status;
    	    }


    	if ((status = H5Sclose(n_cells_dataspace)) < 0) {
    	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
    		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
    	    return status;
    	    }

    	if ((status = H5Sclose(n_cells_memspace)) < 0) {
    	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
    		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
    	    return status;
    	    }

    	if ((status = H5Dclose(n_cells_dataset)) < 0) {
    	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
    		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
    	    return status;
    	    }

    	}

    if(H5Lexists(file_id,"/external_field",H5P_DEFAULT) > 0)
    	{
    	const hsize_t hsize_ncells = p->n_cells*p->n_types;
	hid_t n_cells_memspace = H5Screate_simple(1, &(hsize_ncells), NULL);

    	hid_t n_cells_dataset = H5Dopen2(file_id, "/external_field", H5P_DEFAULT);
    	hid_t n_cells_dataspace = H5Dget_space(n_cells_dataset);
	hsize_t dims;
	status = H5Sget_simple_extent_dims(n_cells_dataspace,&dims,NULL);
	HDF5_ERROR_CHECK(status);
	assert(dims == hsize_ncells);

	assert(p->n_types > 0);
	p->external_field_unified = (soma_scalar_t*)malloc( hsize_ncells*p->n_types*sizeof(soma_scalar_t));
	if (p->external_field_unified == NULL) {
	    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__,
		    __LINE__);
	    return -1;
	    }

    	if ((status =
    	     H5Dread(n_cells_dataset, H5T_SOMA_NATIVE_SCALAR, n_cells_memspace,
    		     n_cells_dataspace, plist_id, p->external_field_unified)) < 0) {
    	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
    		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
    	    return status;
    	    }


    	if ((status = H5Sclose(n_cells_dataspace)) < 0) {
    	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
    		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
    	    return status;
    	    }

    	if ((status = H5Sclose(n_cells_memspace)) < 0) {
    	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
    		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
    	    return status;
    	    }

    	if ((status = H5Dclose(n_cells_dataset)) < 0) {
    	    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
    		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
    	    return status;
    	    }

    	}

    if(H5Lexists(file_id,"/umbrella_field",H5P_DEFAULT) > 0)
        {
        const hsize_t hsize_ncells = p->n_cells*p->n_types;
	hid_t n_cells_memspace = H5Screate_simple(1, &(hsize_ncells), NULL);

        hid_t n_cells_dataset = H5Dopen2(file_id, "/umbrella_field", H5P_DEFAULT);
        hid_t n_cells_dataspace = H5Dget_space(n_cells_dataset);
	hsize_t dims;
	status = H5Sget_simple_extent_dims(n_cells_dataspace,&dims,NULL);
	HDF5_ERROR_CHECK(status);
	assert(dims == hsize_ncells);

	assert(p->n_types > 0);
	p->umbrella_field = (soma_scalar_t*)malloc( hsize_ncells*p->n_types*sizeof(soma_scalar_t));
	if (p->umbrella_field == NULL) {
	    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__,
		    __LINE__);
	    return -1;
	    }

        if ((status =
             H5Dread(n_cells_dataset, H5T_SOMA_NATIVE_SCALAR, n_cells_memspace,
		     n_cells_dataspace, plist_id, p->umbrella_field)) < 0) {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
            }


        if ((status = H5Sclose(n_cells_dataspace)) < 0) {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
            }

        if ((status = H5Sclose(n_cells_memspace)) < 0) {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
            }

        if ((status = H5Dclose(n_cells_dataset)) < 0) {
            fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		    p->info_MPI.world_rank, __FILE__, __LINE__, status);
            return status;
            }

        }


    if ((status = H5Fclose(file_id)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		p->info_MPI.world_rank, __FILE__, __LINE__, status);
	return status;
	}


    return status;
    }

int screen_output(struct Phase*const p,const unsigned int Nsteps)
    {
    static time_t last_print = 0;
    static unsigned int last_time = 0;
    static double last_sec = 0;
    if(last_time == 0) last_time = p->start_time;
    if(last_sec == 0) last_sec = MPI_Wtime();

    const time_t now = time(NULL); if(last_print == 0) last_print = now;
    const unsigned int second = now - p->start_clock;
    const unsigned int steps_done = p->time - p->start_time;
    const time_t end = p->start_clock + second  * (Nsteps)/(soma_scalar_t)steps_done ;
    const double now_sec = MPI_Wtime();
    const double sec =  now_sec - last_sec;
    p->tps_elapsed_time += sec;
    p->tps_elapsed_steps += p->time-last_time;

    if( p->args.screen_output_interval_arg > 0 &&  now - last_print >= p->args.screen_output_interval_arg)
	{
	const double tps = p->tps_elapsed_steps / p->tps_elapsed_time;
	p->tps_elapsed_time = 1./tps;
	p->tps_elapsed_steps = 1;

	if(p->info_MPI.sim_rank == 0)
	    {
	    fprintf(stdout,"Rank %i:Running for %u [s] | TPS %g | steps-to-go: %u | ETA: %s",
		    p->info_MPI.world_rank,second,tps,Nsteps-steps_done,ctime(&end));
	    fflush(stdout);
	    }
	last_print = now;
	last_time = p->time;
	last_sec = now_sec;
	}
    return 0;
    }
