#include <hdf5.h>
#include "readIn_lib.h"
#include "io.h"
#include <math.h>
#include <stdlib.h>
#include "rng.h"

void read_in(const char*const filename, bead_data * const b, polymer_data * const poly)
{
  herr_t status;
  hid_t  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  hid_t  file_id=H5Fopen(filename, H5F_ACC_RDONLY,plist_id);
  H5Pclose(plist_id);
  plist_id=H5Pcreate(H5P_DATASET_XFER);

  float *boxsize = (float*) malloc( 3 * sizeof(float));
  status = read_hdf5(file_id,"/parameter/lxyz",H5T_NATIVE_FLOAT,plist_id,boxsize);
  //HDF5_ERROR_CHECK2(status,"/parameter/lxyz");
  
  uint64_t Npolymers;
  status=read_hdf5(file_id,"/parameter/n_polymers",H5T_NATIVE_UINT64,plist_id,&Npolymers);
  //HDF5_ERROR_CHECK2(status,"/parameter/n_polymers");

  unsigned int poly_arch_len;
  status = read_hdf5(file_id,"/parameter/poly_arch_length",H5T_NATIVE_UINT,plist_id,&poly_arch_len);
  //HDF5_ERROR_CHECK2(status,"poly_arch_length");

  uint32_t*const poly_arch = (uint32_t*) malloc( poly_arch_len * sizeof(uint32_t));
  if(poly_arch == NULL){
    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
    //return -1;
  }
  status = read_hdf5(file_id,"/parameter/poly_arch",H5T_NATIVE_INT32,plist_id,poly_arch);
  //HDF5_ERROR_CHECK2(status,"poly_arch");  
  
  unsigned int n_poly_type;
  status = read_hdf5(file_id,"/parameter/n_poly_type",H5T_NATIVE_UINT,plist_id,&n_poly_type);
  //HDF5_ERROR_CHECK2(status,"n_poly_type");
  
  unsigned int *const poly_type = (unsigned int *const) malloc( Npolymers * sizeof(unsigned int));
  if (poly_type == NULL) {
    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
    //return -1;
  }
  status = read_hdf5(file_id,"/poly_type",H5T_NATIVE_UINT,plist_id,poly_type);
  //HDF5_ERROR_CHECK2(status,"poly_type");
    
  int*const poly_type_offset = (int*)malloc(n_poly_type * sizeof(int));
  if(poly_type_offset == NULL){
    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
    //return -1;
  }
  status = read_hdf5(file_id,"/parameter/poly_type_offset",H5T_NATIVE_INT,plist_id,poly_type_offset);
  //HDF5_ERROR_CHECK2(status,"poly_type_offset");
  
  unsigned int*const number_of_beads= (unsigned int*const) malloc(Npolymers*sizeof(unsigned int));
  if (number_of_beads == NULL) {
    fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
    //return -1;
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
    //return -1;
  }
  
  if(H5Lexists(file_id,"/beads",H5P_DEFAULT) > 0)
    {
      hid_t monomer_memtype = get_monomer_memtype();
      read_hdf5( file_id, "/beads",monomer_memtype,plist_id,beads);
      if ((status = H5Tclose(monomer_memtype)) < 0) {
	fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
		0, __FILE__, __LINE__, status);
	//return status;
      }
    }
  else
    {
      fprintf(stderr,"ERROR: Given file does not contain bead data\n.");
      //return -2;
    }

  H5Pclose(plist_id);
 
  
  poly->n_poly_type=n_poly_type;
  poly->poly_arch=poly_arch;
  poly->poly_type=poly_type;
  poly->poly_type_offset=poly_type_offset;
  poly->boxsize=boxsize;
  poly->poly_arch_length=poly_arch_len;
 
  if ((status = H5Fclose(file_id)) < 0) {
    fprintf(stderr, "ERROR: core: %d HDF5-error %s:%d code %d\n",
	    0, __FILE__, __LINE__, status);
    //return status;
  }
  b->old_N_polymers=Npolymers;
  b->N_polymers=Npolymers;
  b->number_of_beads = number_of_beads;
  b->max_n_beads = max_n_beads;
  b->beads = beads;
  
}


void free_data(bead_data * b,polymer_data*poly)
{
  free(b->beads);
  free(b->number_of_beads); 

  free(poly->poly_arch);
  free(poly->poly_type);
  free(poly->poly_type_offset);
  free(poly->boxsize);
}

int write_out(const char*const filename, bead_data * const b){
  herr_t status;
  
  hid_t plist_id=H5Pcreate(H5P_FILE_ACCESS);
  
  hid_t file_id=H5Fopen(filename,H5F_ACC_RDWR,plist_id);
  HDF5_ERROR_CHECK2(file_id,"open file");
  H5Pclose(plist_id);
  plist_id=H5Pcreate(H5P_DATASET_XFER);

  hid_t dataset=H5Dopen(file_id, "/beads", H5P_DEFAULT);
  HDF5_ERROR_CHECK2(dataset, "open beads");

  const hid_t memtype = get_monomer_memtype();
  status=H5Dwrite(dataset,memtype,H5S_ALL,H5S_ALL,plist_id,b->beads);
  
  HDF5_ERROR_CHECK2(status,"write beads");
  status = H5Dclose(dataset);
  HDF5_ERROR_CHECK2(status,"close beads");
  
  H5Pclose(plist_id);
  status = H5Fclose(file_id);
  HDF5_ERROR_CHECK2(status,"close file");
  return 0;
}



void convert(bead_data * const b,polymer_data*const poly, const char*const outputname)
{ 
  FILE * fp;
  fp = fopen (outputname, "w");

  int sequence=0;
  for (unsigned int i=0;i<b->N_polymers;i++){
    for(unsigned int j=0; j<b->number_of_beads[i]; j++){      
      if (get_particle_type(poly->poly_arch[poly->poly_type_offset[poly->poly_type[i]] + j + 1]) == 0)
	fprintf(fp,"atom %d radius 0.25 name APOLY\n", sequence);               
      else                                                                      
	fprintf(fp,"atom %d radius 0.25 name BPOLY\n", sequence);         
      sequence++;
    }
  }
  int bond_start=poly->poly_type_offset[1];

  for (int beadi=0;beadi<sequence;beadi++){
    int off_bond=0;
    if(poly->poly_arch[bond_start+beadi]>4294967200){
      fprintf(stderr, "ERROR: Function get_bondlist_offset exited the range %s:%d\n", __FILE__, __LINE__);
    }
    off_bond=get_bondlist_offset(poly->poly_arch[bond_start+beadi]);
    while(0==0){   
      int offset=get_offset(poly->poly_arch[off_bond]);
     	
      fprintf(fp,"b %i:%i\n",beadi,offset+beadi);
    
      if (get_end(poly->poly_arch[off_bond])==1) break;
      off_bond++;
    }  
  }
  fprintf(fp,"pbc %f %f %f \n", poly->boxsize[0], poly->boxsize[1], poly->boxsize[2] );
  fprintf(fp,"\ntimestep\n");
  
  for(int i=0;i<sequence;i++){                                            
    fprintf(fp,"%f %f %f\n", b->beads[i].x, b->beads[i].y, b->beads[i].z);
  }
  fclose(fp); 
}


Monomer * feed_back_box(/*const char*const filename,*/bead_data*b,polymer_data*poly)
{
  //  herr_t status;
  
  //  hid_t plist_id=H5Pcreate(H5P_FILE_ACCESS);
  
  //  hid_t file_id=H5Fopen(filename,H5F_ACC_RDWR,plist_id);
  //HDF5_ERROR_CHECK2(file_id,"open file");
  //  H5Pclose(plist_id);
  //  plist_id=H5Pcreate(H5P_DATASET_XFER);

  float boxsize[3];
  boxsize[0]=poly->boxsize[0];
  boxsize[1]=poly->boxsize[1];
  boxsize[2]=poly->boxsize[2];
  int sequence=0;
  int sequence_tmp=0;
  for (unsigned int i=0;i<b->old_N_polymers;i++){
    for(unsigned int j=0; j< b->number_of_beads[i]; j++){
      sequence_tmp++;   
    }
  }
  static Monomer * new_position;
  new_position=malloc(sequence_tmp*sizeof(Monomer));
  
  int tmp=0;
  int tmp1=0;

  for (unsigned int i=0;i<b->old_N_polymers;i++){
    for(unsigned int j=0; j< b->number_of_beads[i]; j++){
      if(b->beads[sequence].x/boxsize[0]<0)
	tmp1=1;
      else
	tmp1=0;
            
      tmp=b->beads[sequence].x/boxsize[0]-tmp1;
  
      new_position[sequence].x=b->beads[sequence].x-boxsize[0]*tmp;

      if(b->beads[sequence].y/boxsize[0]<0)
	tmp1=1;
      else
	tmp1=0;
       
      tmp=b->beads[sequence].y/boxsize[1]-tmp1;
      new_position[sequence].y=b->beads[sequence].y-boxsize[1]*tmp;

      if(b->beads[sequence].z/boxsize[0]<0)
	tmp1=1;
      else
	tmp1=0;
       
      tmp=b->beads[sequence].z/boxsize[2]-tmp1;
      new_position[sequence].z=b->beads[sequence].z-boxsize[2]*tmp;
      sequence++;
    }
  }
  
  return new_position;
}



int change_chi(const char* const filename){
  //herr_t status;

  hid_t dataset_id,dataspace_id;  

  hid_t plist_id=H5Pcreate(H5P_FILE_ACCESS); 
  
  hid_t file_id=H5Fopen(filename,H5F_ACC_RDWR,plist_id);
  HDF5_ERROR_CHECK2(file_id,"open file");
  H5Pclose(plist_id);
  plist_id=H5Pcreate(H5P_DATASET_XFER);

  float tmp[2];
  float tmp_xn[2][2];
  /*status = */read_hdf5(file_id,"/parameter/xn",H5T_NATIVE_FLOAT,plist_id,tmp);

  tmp_xn[1][1]=tmp[0];
  tmp_xn[0][0]=tmp[0];
  tmp_xn[1][0]=20.0;
  tmp_xn[0][1]=20.0;
  printf("%f", tmp_xn[1][0]);
  
  hsize_t    dims_xn[2];
   
  dims_xn[0]=2;
  dims_xn[1]=2;
  dataspace_id = H5Screate_simple(2, dims_xn, NULL);
  /*status=*/H5Ldelete(file_id,"/parameter/xn",H5P_DEFAULT);

  dataset_id = H5Dcreate2(file_id, "/parameter/xn", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /*status=*/H5Dwrite(dataset_id,H5T_NATIVE_FLOAT,H5S_ALL,H5S_ALL,plist_id,&tmp_xn);
  /*status = */H5Dclose(dataset_id);
  /*status = */H5Fclose(file_id);
  
  return 0;
}



int cross_link_new(const char*const filename, bead_data * const b, polymer_data * const poly, unsigned int start_poly, unsigned int end_poly){
  
  //herr_t status;
  hid_t dataset_id, dataspace_id;  
  hid_t plist_id=H5Pcreate(H5P_FILE_ACCESS); 
  
  hid_t file_id=H5Fopen(filename,H5F_ACC_RDWR,plist_id);
  HDF5_ERROR_CHECK2(file_id,"open file");
  H5Pclose(plist_id);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  unsigned int n_crosslink=end_poly-start_poly;
  
  /*change n_polymers*/
  hsize_t     dims_n_polymers[1];  
  unsigned int tmp_n_polymers=b->N_polymers-n_crosslink;
  
  dims_n_polymers[0] = 1;
  
  dataspace_id = H5Screate_simple(1, dims_n_polymers, NULL); 
  /*status=*/H5Ldelete(file_id,"/parameter/n_polymers",H5P_DEFAULT);  
  dataset_id = H5Dcreate2(file_id, "/parameter/n_polymers", H5T_NATIVE_UINT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
  /*status=*/H5Dwrite(dataset_id,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,plist_id,&tmp_n_polymers);
  /*status = */H5Dclose(dataset_id);

  /*change n_poly_type */
 
  hsize_t     dims_n_poly_type[1];
  unsigned int tmp_n_poly_type=/*poly->n_poly_type+*/1;//changed 07.12
  
  dims_n_poly_type[0] = 1;
  
  dataspace_id = H5Screate_simple(1, dims_n_poly_type, NULL);  
  /*status=*/H5Ldelete(file_id,"/parameter/n_poly_type",H5P_DEFAULT);  
  dataset_id = H5Dcreate2(file_id, "/parameter/n_poly_type", H5T_NATIVE_UINT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
  /*status=*/H5Dwrite(dataset_id,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,plist_id,&tmp_n_poly_type);  
  /*status = */H5Dclose(dataset_id);

  /*change array poly_type*/
  hsize_t     dims_poly_type[1];
  dims_poly_type[0] =  tmp_n_polymers;
  unsigned int *tmp_poly_type = (unsigned int * const)malloc(  tmp_n_polymers * sizeof(unsigned int) );
 
  tmp_poly_type[0]=0;  
  /*for(unsigned int i=1;i< tmp_n_polymers;i++){
    tmp_poly_type[i]=poly->poly_type[i-1+n_crosslink];
    }*/
 
  dataspace_id = H5Screate_simple(1, dims_poly_type, NULL);
  /*status=*/H5Ldelete(file_id,"/poly_type",H5P_DEFAULT);
  dataset_id = H5Dcreate2(file_id, "/poly_type", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /*status=*/H5Dwrite(dataset_id,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,plist_id,tmp_poly_type);
  /*status = */H5Dclose(dataset_id);

  /*change number_of_beads*/

  int old_max_n_beads=b-> max_n_beads;

  unsigned int *tmp_number_of_beads[tmp_n_polymers];
  unsigned int new_n_beads=0;
  for(unsigned int i =0;i<=n_crosslink;i++){   
    new_n_beads=b->number_of_beads[i]+new_n_beads;
  }
  if( b-> max_n_beads< new_n_beads)
    b->max_n_beads = new_n_beads;

  tmp_number_of_beads[0]=&new_n_beads; ///really??

  for(unsigned int i=1;i<tmp_n_polymers;i++){
    tmp_number_of_beads[i]=&b->number_of_beads[i+n_crosslink];
  }
 
  /*change size of the beads box*/
  hsize_t     dims_beads[2];

  dims_beads[0] = tmp_n_polymers;  
  dims_beads[1] = b->max_n_beads;
  
  Monomer *tmp_beads = (Monomer * const) malloc(b->max_n_beads*tmp_n_polymers * sizeof(Monomer));

  const hid_t memtype = get_monomer_memtype();

  dataspace_id = H5Screate_simple(2, dims_beads, NULL);
  
  for(unsigned int i=0;i<1;i++){
    for(unsigned int j=0;j<*tmp_number_of_beads[i];j++){
      tmp_beads[i*b->max_n_beads+j].x=b->beads[i*old_max_n_beads+j].x;
      tmp_beads[i*b->max_n_beads+j].y=b->beads[i*old_max_n_beads+j].y;
      tmp_beads[i*b->max_n_beads+j].z=b->beads[i*old_max_n_beads+j].z;
    }
  }

  for(unsigned int i=1;i<tmp_n_polymers;i++){
    for(unsigned int j=0;j<*tmp_number_of_beads[i];j++){
      tmp_beads[i*b->max_n_beads+j].x=b->beads[(i+n_crosslink)*old_max_n_beads+j].x;
      tmp_beads[i*b->max_n_beads+j].y=b->beads[(i+n_crosslink)*old_max_n_beads+j].y;
      tmp_beads[i*b->max_n_beads+j].z=b->beads[(i+n_crosslink)*old_max_n_beads+j].z;
    }    
  }

  /*status=*/H5Ldelete(file_id,"/beads",H5P_DEFAULT);
  dataset_id = H5Dcreate2(file_id, "/beads", memtype,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  /*status=*/H5Dwrite(dataset_id,memtype,H5S_ALL,H5S_ALL,plist_id,tmp_beads);
  /*status =*/ H5Dclose(dataset_id);

  /*change poly_arch*/
  hsize_t     dims_poly_arch[1];
   
  unsigned int Nmono1=poly->poly_arch_length-4-1;
  unsigned int Nmono2=(n_crosslink+1)*Nmono1;
  unsigned int tmp_poly_arch_length = /*1 +Nmono1+*/1+ Nmono2 + 4;//change 07.12
  dims_poly_arch[0]=tmp_poly_arch_length;

  dataspace_id = H5Screate_simple(1, dims_poly_arch, NULL);
   
  uint32_t * tmp_poly_arch = (uint32_t*) malloc( tmp_poly_arch_length * sizeof(uint32_t));
  if(tmp_poly_arch == NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
    return -4;
  }
 
  // tmp_poly_arch[0]=Nmono1;

  // tmp_poly_arch[1] = get_info_bl(1+ Nmono1 +Nmono2+2, get_particle_type(poly->poly_arch[1]));
    
  /* for(unsigned int i=2; i <= Nmono1-2+1; i++)
    {	 
      tmp_poly_arch[ i ] = get_info_bl(1+ Nmono1+Nmono2+3,get_particle_type(poly->poly_arch[i]));
      }*/
  //tmp_poly_arch[Nmono1] = get_info_bl(1+ Nmono1 +Nmono2+1, get_particle_type(poly->poly_arch[Nmono1])); 
   
  tmp_poly_arch[/*Nmono1+1*/0]=Nmono2;
   
  unsigned int tmp_start=1/*Nmono1+2*/;
  
  for(unsigned int j=0;j<=n_crosslink;j++){
    unsigned int tmp_start_short=1;
    tmp_poly_arch[tmp_start] = get_info_bl(/*1+ Nmono1 +*/Nmono2+2, get_particle_type(poly->poly_arch[tmp_start_short]));
     
    for(unsigned int i=tmp_start+1; i <= tmp_start+Nmono1-2; i++)
      {
	tmp_start_short++;
	tmp_poly_arch[ i ] = get_info_bl(/*1+ Nmono1+*/Nmono2+3, get_particle_type(poly->poly_arch[tmp_start_short]));
      }
    tmp_poly_arch[tmp_start+ Nmono1-2+1] = get_info_bl(/*1+ Nmono1 +*/Nmono2+1,get_particle_type(poly->poly_arch[tmp_start_short+1]));    
    tmp_start=tmp_start+Nmono1-2+2;     
  }
   
  //Last monomer
  int end,offset,bond_type,info;
  bond_type = HARMONIC;
  end = 1;
  offset = -1;
  info=get_info(offset, bond_type, end);
  tmp_poly_arch[ /*1+Nmono1+*/Nmono2+ 1 ] = info;

  //First monomer
  end = 1;
  offset = 1;
  info = get_info(offset, bond_type, end);
  tmp_poly_arch[ /*1+Nmono1+*/Nmono2+ 2 ] = info;
  //Middle monomers A
  end = 0;
  offset = -1;
  info = get_info(offset, bond_type, end);
  tmp_poly_arch[ /*1+Nmono1+*/Nmono2+ 3 ] = info;
  
  end = 1;
  offset = 1;
  info = get_info(offset, bond_type, end);
  tmp_poly_arch[/*1+Nmono1+*/Nmono2 + 4] = info;

  /*status=*/H5Ldelete(file_id,"/parameter/poly_arch",H5P_DEFAULT);

  dataset_id = H5Dcreate2(file_id, "/parameter/poly_arch", H5T_NATIVE_INT32,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);       
  /*status=*/H5Dwrite(dataset_id,H5T_NATIVE_INT32,H5S_ALL,H5S_ALL,plist_id,tmp_poly_arch);
  /* status = */H5Dclose(dataset_id);

  /*change poly_type_offset*/
  hsize_t     dims_poly_type_offset[1];
  dims_poly_type_offset[0] = 1/*2*/; //change 07.12
   
  int *tmp_offset =  realloc(poly->poly_type_offset, poly->n_poly_type * sizeof(unsigned int) );
  tmp_offset[0]=0/*Nmono1+1*/;

  dataspace_id = H5Screate_simple(1, dims_poly_type_offset, NULL);
  /*status=*/H5Ldelete(file_id,"/parameter/poly_type_offset",H5P_DEFAULT);
  dataset_id = H5Dcreate2(file_id, "/parameter/poly_type_offset", H5T_NATIVE_INT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  /*status=*/H5Dwrite(dataset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,plist_id,tmp_offset);
  /*status = */H5Dclose(dataset_id);

  /*change poly_arch_len*/
  dataset_id=H5Dopen(file_id, "/parameter/poly_arch_length", H5P_DEFAULT);
  /*status=*/H5Dwrite(dataset_id,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,plist_id,&tmp_poly_arch_length);
  /*status = */H5Dclose(dataset_id);  

  poly->n_poly_type=tmp_n_poly_type;
  poly->poly_arch=tmp_poly_arch;
  poly->poly_type=tmp_poly_type;
  poly->poly_type_offset=tmp_offset;
  poly->poly_arch_length=tmp_poly_arch_length;
  b->N_polymers=tmp_n_polymers;	
  /*status = */H5Fclose(file_id);

  return 0;    
}

int cross_link_monomer2(bead_data * const b, polymer_data * const poly,unsigned int polyN,unsigned int monomerN, unsigned int monoM, unsigned int *poly_arch_len_pointer,int32_t ** poly_arch_tmp_pointer)
{

  int info;    
  int32_t* poly_arch_tmp= *poly_arch_tmp_pointer;
 if(poly_arch_tmp == NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
    return -4;
  }
  unsigned int poly_arch_len=*poly_arch_len_pointer;
  unsigned int poly_arch_lenold=poly_arch_len;

  /*load the poly_arch information of monomerN*/
  
  int offset_poly=poly->poly_type_offset[0];
  for(unsigned int i=0;i<polyN-1;i++){
    offset_poly=offset_poly+b->number_of_beads[i];
  }

  int32_t polyarch=poly_arch_tmp[offset_poly+monomerN];
  if(polyarch>4294967200){
    fprintf(stderr, "ERROR: Function get_bondlist_offset exited the range %s:%d\n", __FILE__, __LINE__);
    return -2;
  }
  int off_bond=get_bondlist_offset(polyarch);
  unsigned int typeN=get_particle_type(poly_arch_tmp[offset_poly+monomerN]);
  int bond_type = HARMONIC;
  do{//count the number of bonds
    unsigned int end=get_end(poly_arch_tmp[off_bond]);
    poly_arch_len++;  
    if(end==1) break;  
    off_bond++;    
  }while(0==0);

  /*load the poly_arch information of monomerM*/
  int offset_polyM=poly->poly_type_offset[0];

  int32_t polyarchM=poly_arch_tmp[offset_polyM+monoM];

  if(polyarchM>4294967200){
    fprintf(stderr, "ERROR: Function get_bondlist_offset exited the range %s:%d\n", __FILE__, __LINE__);
    return -2;
  }
  int off_bondM=get_bondlist_offset(polyarchM);
  unsigned int typeM=get_particle_type(poly_arch_tmp[offset_polyM+monoM]);  
  unsigned int poly_arch_lenM=poly_arch_len+1;
  unsigned int poly_arch_lenoldM=poly_arch_lenM;
  do{//count the number of bonds
    unsigned int end=get_end(poly_arch_tmp[off_bondM]); 
    poly_arch_lenM++;
    if(end==1) break;  
    off_bondM++;    
  }while(0==0);
  
  int bond_exist[poly_arch_len-poly_arch_lenold];//check if the bond already exsists;
  if(polyarch>4294967200||polyarchM>4294967200){
    fprintf(stderr, "ERROR: Function get_bondlist_offset exited the range %s:%d\n", __FILE__, __LINE__);
    return -2;
  }
  off_bond=get_bondlist_offset(polyarch);
  off_bondM=get_bondlist_offset(polyarchM);
  int32_t * tmp_poly_archM = (int32_t*) malloc((poly_arch_lenM+1) * sizeof(int32_t));
  if(tmp_poly_archM == NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
    return -4;
  }
   
  for(unsigned int i=0;i<poly_arch_lenold;i++){    //copy the old poly_arch into the new one
    tmp_poly_archM[i]=poly_arch_tmp[i];   
  }
 
  for(unsigned int i=poly_arch_lenold;i<poly_arch_len;i++){    //copy the old bonds of the new crosslinked monomer into the array, but with end=0
    int offset=get_offset(poly_arch_tmp[off_bond]);
    bond_exist[i-poly_arch_lenold]=offset;
    info=get_info(offset, bond_type, 0);
    tmp_poly_archM[i] = info;
    off_bond++;
  }

  for(unsigned int i=poly_arch_len+1;i<poly_arch_lenM;i++){
    int offset=get_offset(poly_arch_tmp[off_bondM]);
    info=get_info(offset, bond_type, 0);
    tmp_poly_archM[i] = info;
    off_bondM++;
  }

  int bond_match=0;
  int offset_to_monomerM=offset_polyM+monoM-offset_poly-monomerN;
  
  for(unsigned int i=0;i<poly_arch_len-poly_arch_lenold;i++){ //check if the bond already exists
    if(offset_to_monomerM==bond_exist[i])
      bond_match++;
  }  

  if(bond_match==0 && offset_to_monomerM!=0){ //if the bond does not exist and it is not a bond to itself, write this bond into poly_arch
    info=get_info(offset_to_monomerM, bond_type,1);
    tmp_poly_archM[poly_arch_len]=info; 
    tmp_poly_archM[offset_poly+monomerN] = get_info_bl(poly_arch_lenold, typeN);
  
    int offset_to_monomerN=offset_poly+monomerN-offset_polyM-monoM;
    info=get_info(offset_to_monomerN, bond_type,1);
    tmp_poly_archM[poly_arch_lenM]=info;  
    tmp_poly_archM[offset_polyM+monoM] = get_info_bl(poly_arch_lenoldM, typeM);  

    *poly_arch_tmp_pointer = tmp_poly_archM;    //the poly_arch array is updated, but hdf5 file not changed
    *poly_arch_len_pointer=poly_arch_lenM+1;    
    poly->poly_arch_length=poly_arch_lenM+1;
    free(poly_arch_tmp);
  }  
  else{
    free(tmp_poly_archM);
  }  
  return 0;
}


int cross_link_monomer(const char* const filename, bead_data * const b, polymer_data * const poly,unsigned int polyN,unsigned int monomerN, unsigned int polyM,unsigned int monomerM)
{

  herr_t status;
  hid_t dataset_id, dataspace_id;  
  hid_t plist_id=H5Pcreate(H5P_FILE_ACCESS);  
  hid_t  file_id=H5Fopen(filename, H5F_ACC_RDWR,plist_id);
  HDF5_ERROR_CHECK2(file_id,"open file");
  H5Pclose(plist_id);
  plist_id=H5Pcreate(H5P_DATASET_XFER);

  int info;
  unsigned int poly_arch_len;
  status = read_hdf5(file_id,"/parameter/poly_arch_length",H5T_NATIVE_UINT,plist_id,&poly_arch_len);
  HDF5_ERROR_CHECK2(status,"poly_arch_length");
    
  unsigned int poly_arch_lenold=poly_arch_len;
  int32_t*const poly_arch_tmp = (int32_t*) malloc( poly_arch_len * sizeof(int32_t));
  if(poly_arch_tmp == NULL){
    fprintf(stderr, "ERROR: Malloc poly_arch_tmp %s:%d\n", __FILE__, __LINE__);
    return -1;
  }
  status = read_hdf5(file_id,"/parameter/poly_arch",H5T_NATIVE_INT32,plist_id,poly_arch_tmp);
  HDF5_ERROR_CHECK2(status,"poly_arch");
  
  /*load the poly_arch information of monomerN*/

  int offset_poly=poly->poly_type_offset[0];
  for(unsigned int i=0;i<polyN-1;i++){
    offset_poly=offset_poly+b->number_of_beads[i];
  }
  // printf("%i",poly_arch_len);
  int32_t polyarch=poly_arch_tmp[offset_poly+monomerN];
  if(polyarch>4294967200){
    fprintf(stderr, "ERROR: Function get_bondlist_offset exited the range %s:%d\n", __FILE__, __LINE__);
    return -2;
  }
  int off_bond=get_bondlist_offset(polyarch);
  unsigned int typeN=get_particle_type(poly_arch_tmp[offset_poly+monomerN]);
  int bond_type = HARMONIC;
  do{
    unsigned int end=get_end(poly_arch_tmp[off_bond]);
    poly_arch_len++;   
    if(end==1) break;  
    off_bond++;    
  }while(0==0);

  /*load the poly_arch information of monomerM*/
  int offset_polyM=poly->poly_type_offset[0];
  for(unsigned int i=0;i<polyM-1;i++){
    offset_polyM=offset_polyM+b->number_of_beads[i];
  }
  int32_t polyarchM=poly_arch_tmp[offset_polyM+monomerM];

  if(polyarchM>4294967200){
    fprintf(stderr, "ERROR: Function get_bondlist_offset exited the range %s:%d\n", __FILE__, __LINE__);
    return -2;
  }
  int off_bondM=get_bondlist_offset(polyarchM);
  unsigned int typeM=get_particle_type(poly_arch_tmp[offset_polyM+monomerM]);
  
  unsigned int poly_arch_lenM=poly_arch_len+1;
  unsigned int poly_arch_lenoldM=poly_arch_lenM;
  do{
    unsigned int end=get_end(poly_arch_tmp[off_bondM]); 
    poly_arch_lenM++;
    if(end==1) break;  
    off_bondM++;    
  }while(0==0);
  
  int bond_exist[poly_arch_len-poly_arch_lenold];//check if the bond already exsists;
  if(polyarch>4294967200||polyarchM>4294967200){
    fprintf(stderr, "ERROR: Function get_bondlist_offset exited the range %s:%d\n", __FILE__, __LINE__);
    return -2;
  }
  off_bond=get_bondlist_offset(polyarch);
  off_bondM=get_bondlist_offset(polyarchM);
  int32_t * tmp_poly_archM = (int32_t*) malloc((poly_arch_lenM+1) * sizeof(int32_t));
  if(tmp_poly_archM == NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
    return -4;
  }
   
  for(unsigned int i=0;i<poly_arch_lenold;i++){    
    tmp_poly_archM[i]=poly_arch_tmp[i];
  }
   
  for(unsigned int i=poly_arch_lenold;i<poly_arch_len;i++){    
    int offset=get_offset(poly_arch_tmp[off_bond]);
    bond_exist[i-poly_arch_lenold]=offset;
    info=get_info(offset, bond_type, 0);
    tmp_poly_archM[i] = info;
    off_bond++;
  }

  for(unsigned int i=poly_arch_len+1;i<poly_arch_lenM;i++){
    int offset=get_offset(poly_arch_tmp[off_bondM]);
    info=get_info(offset, bond_type, 0);
    tmp_poly_archM[i] = info;
    off_bondM++;
  }

  int bond_match=0;
  int offset_to_monomerM=offset_polyM+monomerM-offset_poly-monomerN;
  
  for(unsigned int i=0;i<poly_arch_len-poly_arch_lenold;i++){
    if(offset_to_monomerM==bond_exist[i])
      bond_match++;
  }  
  if(bond_match==0 && offset_to_monomerM!=0){
    info=get_info(offset_to_monomerM, bond_type,1);
    tmp_poly_archM[poly_arch_len]=info; 
    tmp_poly_archM[offset_poly+monomerN] = get_info_bl(poly_arch_lenold, typeN);
  
    int offset_to_monomerN=offset_poly+monomerN-offset_polyM-monomerM;
    info=get_info(offset_to_monomerN, bond_type,1);
    tmp_poly_archM[poly_arch_lenM]=info;  
    tmp_poly_archM[offset_polyM+monomerM] = get_info_bl(poly_arch_lenoldM, typeM);

    hsize_t     dims_poly_arch[1];
   
    dims_poly_arch[0]=poly_arch_lenM+1;// the poly_arch array is not updated, but hdf5 file is
  
    dataspace_id = H5Screate_simple(1, dims_poly_arch, NULL);

    status=H5Ldelete(file_id,"/parameter/poly_arch",H5P_DEFAULT);

    dataset_id = H5Dcreate2(file_id, "/parameter/poly_arch", H5T_NATIVE_INT32,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
       
    status=H5Dwrite(dataset_id,H5T_NATIVE_INT32,H5S_ALL,H5S_ALL,plist_id,tmp_poly_archM);
    status = H5Dclose(dataset_id);
    
    /*change poly_arch_len*/
     dataset_id=H5Dopen(file_id, "/parameter/poly_arch_length", H5P_DEFAULT);
    unsigned int tmp_poly_arch_len=poly_arch_lenM+1;
    status=H5Dwrite(dataset_id,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,plist_id,&tmp_poly_arch_len);
    status = H5Dclose(dataset_id);
    poly->poly_arch_length=tmp_poly_arch_len;
  }
  /* Close the file. */
  status = H5Fclose(file_id);
  free(tmp_poly_archM);
  free(poly_arch_tmp);
  return 0;
}

int Prob_bonds(const char* const filename,bead_data*b,polymer_data*poly,const double const pzero, Monomer*new_position,const char* const crosslinkfile){

  cross_link_new(filename,b,poly,0,b->N_polymers-1);
  printf("finish crosslinking! \n");
  hid_t dataset_id, dataspace_id;  
  herr_t status;
  hid_t plist_id=H5Pcreate(H5P_FILE_ACCESS);  
  hid_t file_id=H5Fopen(filename,H5F_ACC_RDWR,plist_id);
  HDF5_ERROR_CHECK2(file_id,"open file");
  H5Pclose(plist_id);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  double sum=0;
  double average;
  unsigned int Nref=0;
  read_hdf5(file_id,"/parameter/reference_Nbeads",H5T_NATIVE_UINT32,plist_id,&Nref);
    
  double bfactor=1/sqrt(Nref);
  unsigned int poly_arch_len;
  read_hdf5(file_id,"/parameter/poly_arch_length",H5T_NATIVE_UINT,plist_id,&poly_arch_len);
  poly->poly_arch_length=poly_arch_len;
  
  int32_t* poly_arch_tmp = (int32_t*) malloc( poly_arch_len * sizeof(int32_t));
  if(poly_arch_tmp == NULL){
    fprintf(stderr, "ERROR: Malloc poly_arch_tmp %s:%d\n", __FILE__, __LINE__);
  }
  read_hdf5(file_id,"/parameter/poly_arch",H5T_NATIVE_INT32,plist_id,poly_arch_tmp);
  HDF5_ERROR_CHECK2(status,"poly_arch");

  unsigned int sequence=0;
  
  for(unsigned int i=0;i<b->old_N_polymers;i++){
    sequence=sequence+b->number_of_beads[i];
  }
  
  double gitter_size=0.25; //the size of the box, crosslink is only allowed in neighbour boxes
  int roundr=0,roundw=0,roundh=0;
  double rows=poly->boxsize[0]/gitter_size, width=poly->boxsize[1]/gitter_size,height=poly->boxsize[2]/gitter_size;
  if(rows ==(int) rows)
    roundr=0;
  else
    roundr=1;
   if(width ==(int) width)
    roundw=0;
  else
    roundw=1;
   if(height ==(int) height)
    roundh=0;
  else
    roundh=1;
   
   //start:book memory for particle box
  unsigned int ***index=(unsigned int ***)malloc((int)(rows+roundr)*sizeof(unsigned int **));
  unsigned int ****particle_box=(unsigned int ****)malloc((int)(rows+roundr)*sizeof(unsigned int ***));
  if(particle_box == NULL){
    fprintf(stderr, "ERROR: Malloc poly_arch_tmp %s:%d\n", __FILE__, __LINE__);
  }
  for(unsigned int row_i=0;row_i<rows;row_i++){
    particle_box[row_i]=(unsigned int ***) malloc ((int)(width+roundw)*sizeof(unsigned int **));
    index[row_i]=(unsigned int **)malloc((int)(width+roundw)*sizeof(unsigned int *));
    for(unsigned int width_i=0;width_i<width;width_i++){
      particle_box[row_i][width_i]=(unsigned int **) malloc ((int)(height+roundh)*sizeof(unsigned int *));
      index[row_i][width_i]=(unsigned int *)malloc((int)(height+roundh)*sizeof(unsigned int));   
 
      if(index[row_i][width_i] == NULL){
	  fprintf(stderr, "ERROR: Malloc poly_arch_tmp %s:%d\n", __FILE__, __LINE__);
      }
      for(unsigned int height_i=0;height_i<height;height_i++){
	index[row_i][width_i][height_i]=0;
	particle_box[row_i][width_i][height_i]=(unsigned int *) malloc ((int)(sequence/50)*sizeof(unsigned int));
	if(particle_box[row_i][width_i][height_i] == NULL){
	  fprintf(stderr, "ERROR: Malloc poly_arch_tmp %s:%d\n", __FILE__, __LINE__);
	 }
      }
    }
  }//end:book memory for particle box


  
  int index_offset=0;//start:assign particles to boxes
  for(unsigned int Poly_i=0;Poly_i<b->old_N_polymers;Poly_i++){   
    for(unsigned int Mono_i=0;Mono_i<b->number_of_beads[Poly_i];Mono_i++){   
      particle_box[(int)(new_position[index_offset+Mono_i].x/gitter_size)][(int)(new_position[index_offset+Mono_i].y/gitter_size)][(int)(new_position[index_offset+Mono_i].z/gitter_size)][index[(int)(new_position[index_offset+Mono_i].x/gitter_size)][(int)(new_position[index_offset+Mono_i].y/gitter_size)][(int)(new_position[index_offset+Mono_i].z/gitter_size)]]=index_offset+Mono_i;     
      index[(int)(new_position[index_offset+Mono_i].x/gitter_size)][(int)(new_position[index_offset+Mono_i].y/gitter_size)][(int)(new_position[index_offset+Mono_i].z/gitter_size)]++;
      if(index[(int)(new_position[index_offset+Mono_i].x/gitter_size)][(int)(new_position[index_offset+Mono_i].y/gitter_size)][(int)(new_position[index_offset+Mono_i].z/gitter_size)]>(int)(sequence/50)){
	fprintf(stderr, "ERROR: Malloced memory for particle box not enough %s:%d\n", __FILE__, __LINE__);
	return -4;
      }
    }
    index_offset=index_offset+b->number_of_beads[Poly_i];
  }//end:assign particles to boxes
 
  PCG_STATE state;
  PCG_STATE*const s = &state;
  soma_seed_rng(s, time(NULL), 0);
  
  int ** number_crosslink=(int **)malloc(b->old_N_polymers*sizeof(int*));
  //printf("start initialisation %i\n",b->old_N_polymers);
  for(unsigned int a=0;a<b->old_N_polymers;a++){
    number_crosslink[a]=(int *)malloc(b->old_N_polymers*sizeof(int));

    for(unsigned int bb=0;bb<b->old_N_polymers;bb++){
      number_crosslink[a][bb]=0;
    }
  }
  printf("finish initialisation\n");
  int *bond_polymer=(int*)malloc(b->old_N_polymers*sizeof(int));
  if(bond_polymer == NULL){
    fprintf(stderr, "ERROR: Malloc bond_polymer %s:%d\n", __FILE__, __LINE__);
  }

  double self_crosslink=0;
  double double_crosslink=0;
  double P0=pzero;
  index_offset=0;
  int range=(int)(5*bfactor+1);

  for(unsigned int Poly_i=0;Poly_i<b->old_N_polymers;Poly_i++){
    bond_polymer[Poly_i]=0;
  
    for(unsigned int Mono_i=0;Mono_i<b->number_of_beads[Poly_i];Mono_i++){
      
      for(unsigned int row_i=fmax(0,(int)(new_position[index_offset+Mono_i].x/gitter_size)-range);row_i<fmin((int)poly->boxsize[0]/gitter_size,(int)(new_position[index_offset+Mono_i].x/gitter_size)+range+1);row_i++){

	for(unsigned int width_i=fmax(0,(int)(new_position[index_offset+Mono_i].y/gitter_size)-range);width_i<fmin((int)poly->boxsize[1]/gitter_size,(int)(new_position[index_offset+Mono_i].y/gitter_size)+range+1);width_i++){

	  for(unsigned int height_i=fmax(0,(int)(new_position[index_offset+Mono_i].z/gitter_size)-range);height_i<fmin((int)poly->boxsize[2]/gitter_size,(int)(new_position[index_offset+Mono_i].z/gitter_size)+range+1);height_i++){
	    
	     for(unsigned int index_i=0;index_i<index[(int)(new_position[index_offset+Mono_i].x/gitter_size)][(int)(new_position[index_offset+Mono_i].y/gitter_size)][(int)(new_position[index_offset+Mono_i].z/gitter_size)];index_i++){
	       unsigned int second_tmp=particle_box[row_i][width_i][height_i][index_i];	  
	      double distance_sq=(new_position[index_offset+Mono_i].x-new_position[second_tmp].x)*(new_position[index_offset+Mono_i].x-new_position[second_tmp].x)+(new_position[index_offset+Mono_i].y-new_position[second_tmp].y)*(new_position[index_offset+Mono_i].y-new_position[second_tmp].y)+(new_position[index_offset+Mono_i].z-new_position[second_tmp].z)*(new_position[index_offset+Mono_i].z-new_position[second_tmp].z);
	      float Prob=P0*exp(-1.5/bfactor/bfactor*distance_sq);	
	      unsigned int rng1 = pcg32_random(s);	      
	      float random=rng1/(float)soma_rng_uint_max();
	      if (random<Prob&&(Poly_i!=(int)second_tmp/32||Mono_i!=second_tmp-(int)(second_tmp/32)*32)){
		int old=poly_arch_len;
		cross_link_monomer2(b,poly,Poly_i+1,Mono_i+1,second_tmp+1,&poly_arch_len,&poly_arch_tmp);
		if(old==poly_arch_len){
		  bond_polymer[Poly_i]--;
		  bond_polymer[(int)second_tmp/32]--;
		}
		if(Poly_i==(int)(second_tmp/32))
		  self_crosslink++;
		if((number_crosslink[Poly_i][(int)second_tmp/32]!=0)||(number_crosslink[(int)second_tmp/32][Poly_i]!=0))
		  double_crosslink++;
		number_crosslink[Poly_i][(int)second_tmp/32]=1;	
		number_crosslink[(int)second_tmp/32][Poly_i]=1;	
		bond_polymer[Poly_i]++;
		bond_polymer[(int)second_tmp/32]++;
	      }	     
	     }
	  } 	
	}

      }    
    }	   
    //printf("bond %i, %i\n",Poly_i,bond_polymer[Poly_i]);
    sum=sum+bond_polymer[Poly_i];
    index_offset=index_offset+b->number_of_beads[Poly_i];
  }//end of the loop over all particles  
  average=sum/b->old_N_polymers;
  double self_ave=self_crosslink/b->old_N_polymers;
  double double_ave=double_crosslink/b->old_N_polymers;
  printf("number of polymers %"PRIu64"\n",b->old_N_polymers);
  printf("average %f \n",average);

  FILE *f = fopen(crosslinkfile, "a+");
  if (f == NULL)
    {
      printf("Error opening file!\n");
    }
  fprintf(f, "%f %f %f\n",average,self_crosslink,double_crosslink);
   
  fclose(f);
  hsize_t     dims_poly_arch[1];
   
  dims_poly_arch[0]=poly_arch_len;
  
  dataspace_id = H5Screate_simple(1, dims_poly_arch, NULL);

   H5Ldelete(file_id,"/parameter/poly_arch",H5P_DEFAULT);

  dataset_id = H5Dcreate2(file_id, "/parameter/poly_arch", H5T_NATIVE_INT32,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
       
  H5Dwrite(dataset_id,H5T_NATIVE_INT32,H5S_ALL,H5S_ALL,plist_id,poly_arch_tmp);
  H5Dclose(dataset_id);
   
  /*change poly_arch_len*/
    
  dataset_id=H5Dopen(file_id, "/parameter/poly_arch_length", H5P_DEFAULT);
  H5Dwrite(dataset_id,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,plist_id,&poly_arch_len);
  H5Dclose(dataset_id);
  H5Fclose(file_id);

  free(bond_polymer);
}       


int crosslink1(const char* const filename,const double const pzero,const char* const crosslinkfile)
{
  polymer_data polym;
  polymer_data*const poly=&polym;
  bead_data data;
 
  bead_data*const b=&data;
  data.number_of_beads = NULL;
  data.beads = NULL;
  polym.poly_arch=NULL;
  polym.poly_type=NULL;
  polym.poly_type_offset=NULL;
  polym.boxsize=NULL;
   
  read_in(filename,b,poly);
  printf("finished reading! \n");
   
  Monomer *feed=feed_back_box(/*filename,*/b,poly);
  printf("finish feed back! \n");
   
  Prob_bonds(filename,b,poly,pzero,feed,crosslinkfile);
  printf("finish pbonds! \n");
   
  static int number_of_perkolation=0;
   
  number_of_perkolation=check_perkolation_monomer1(filename,b,poly,0);
  free_data(b,poly);
  return number_of_perkolation;
}

void crosslink(const char* const filename,const double pzero,const char* const crosslinkfile)
{

  polymer_data polym;
  polymer_data*const poly=&polym;
  bead_data data;
 
  bead_data*const b=&data;
  data.number_of_beads = NULL;
  data.beads = NULL;
  polym.poly_arch=NULL;
  polym.poly_type=NULL;
  polym.poly_type_offset=NULL;
  polym.boxsize=NULL;
   
  read_in(filename,b,poly);
  printf("finished reading! \n");
   
  Monomer *feed=feed_back_box(b,poly);
  printf("finish feed back! \n");
   
  Prob_bonds(filename,b,poly,pzero,feed,crosslinkfile);
  printf("finish pbonds! \n");   
   
  free_data(b,poly);
}


