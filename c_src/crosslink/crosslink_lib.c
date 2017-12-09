#include <hdf5.h>
#include "readIn_lib.h"
#include "crosslink_lib.h"
#include "../io.h"
#include <math.h>
#include <stdlib.h>

int cross_link_new(const char*const filename, bead_data * const b, polymer_data * const poly, unsigned int start_poly, unsigned int end_poly){
  
  herr_t status;
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
  status=H5Ldelete(file_id,"/parameter/n_polymers",H5P_DEFAULT);  
  dataset_id = H5Dcreate2(file_id, "/parameter/n_polymers", H5T_NATIVE_UINT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
  status=H5Dwrite(dataset_id,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,plist_id,&tmp_n_polymers);
  status = H5Dclose(dataset_id);

   /*change n_poly_type */
 
  hsize_t     dims_n_poly_type[1];
  unsigned int tmp_n_poly_type=poly->n_poly_type+1;
  
  dims_n_poly_type[0] = 1;
  
  dataspace_id = H5Screate_simple(1, dims_n_poly_type, NULL);  
  status=H5Ldelete(file_id,"/parameter/n_poly_type",H5P_DEFAULT);  
  dataset_id = H5Dcreate2(file_id, "/parameter/n_poly_type", H5T_NATIVE_UINT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
  status=H5Dwrite(dataset_id,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,plist_id,&tmp_n_poly_type);  
  status = H5Dclose(dataset_id);

   /*change array poly_type*/
  hsize_t     dims_poly_type[1];
  dims_poly_type[0] =  tmp_n_polymers;
  unsigned int *tmp_poly_type = (unsigned int * const)malloc(  tmp_n_polymers * sizeof(unsigned int) );
 
  tmp_poly_type[0]=1;  
  for(unsigned int i=1;i< tmp_n_polymers;i++){
    tmp_poly_type[i]=poly->poly_type[i+n_crosslink];
  }
 
  dataspace_id = H5Screate_simple(1, dims_poly_type, NULL);
  status=H5Ldelete(file_id,"/poly_type",H5P_DEFAULT);
  dataset_id = H5Dcreate2(file_id, "/poly_type", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status=H5Dwrite(dataset_id,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,plist_id,tmp_poly_type);
  status = H5Dclose(dataset_id);

  /*change number_of_beads*/

  int old_max_n_beads=b-> max_n_beads;

  unsigned int *tmp_number_of_beads[tmp_n_polymers];
  unsigned int new_n_beads=0;
  for(unsigned int i =0;i<=n_crosslink;i++){   
    new_n_beads=b->number_of_beads[i]+new_n_beads;
  }
  if( b-> max_n_beads< new_n_beads)
    b->max_n_beads = new_n_beads;

  tmp_number_of_beads[0]=new_n_beads;

  for(unsigned int i=1;i<tmp_n_polymers;i++){
    tmp_number_of_beads[i]=b->number_of_beads[i+n_crosslink];
    }
   /*int *tmp_number_of_beads =  realloc(b->number_of_beads, b->N_polymers * sizeof(unsigned int) );*/

  /*change size of the beads box*/
   hsize_t     dims_beads[2];

  dims_beads[0] = tmp_n_polymers;  
  dims_beads[1] = b->max_n_beads;
  
  Monomer *tmp_beads = (Monomer * const) malloc(b->max_n_beads*tmp_n_polymers * sizeof(Monomer));

  const hid_t memtype = get_monomer_memtype();
  
  dataspace_id = H5Screate_simple(2, dims_beads, NULL);
  
  for(unsigned int i=0;i<1;i++){
    for(unsigned int j=0;j<tmp_number_of_beads[i];j++){
      tmp_beads[i*b->max_n_beads+j].x=b->beads[i*old_max_n_beads+j].x;
      tmp_beads[i*b->max_n_beads+j].y=b->beads[i*old_max_n_beads+j].y;
      tmp_beads[i*b->max_n_beads+j].z=b->beads[i*old_max_n_beads+j].z;
    }
  }

  for(unsigned int i=1;i<tmp_n_polymers;i++){
    for(unsigned int j=0;j<tmp_number_of_beads[i];j++){
      tmp_beads[i*b->max_n_beads+j].x=b->beads[(i+n_crosslink)*old_max_n_beads+j].x;
      tmp_beads[i*b->max_n_beads+j].y=b->beads[(i+n_crosslink)*old_max_n_beads+j].y;
      tmp_beads[i*b->max_n_beads+j].z=b->beads[(i+n_crosslink)*old_max_n_beads+j].z;
    }    
  }

  status=H5Ldelete(file_id,"/beads",H5P_DEFAULT);
  dataset_id = H5Dcreate2(file_id, "/beads", memtype,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status=H5Dwrite(dataset_id,memtype,H5S_ALL,H5S_ALL,plist_id,tmp_beads);
  status = H5Dclose(dataset_id);

    /*change poly_arch*/
   hsize_t     dims_poly_arch[1];
   
   unsigned int Nmono1=poly->poly_arch_length-4-1;
   unsigned int Nmono2=(n_crosslink+1)*Nmono1;
   
   //unsigned int NmonoA=63, NmonoB=63;//be careful!
   unsigned int tmp_poly_arch_length = 1 +Nmono1+1+ Nmono2 + 4;
   dims_poly_arch[0]=tmp_poly_arch_length;

   dataspace_id = H5Screate_simple(1, dims_poly_arch, NULL);
   
   uint32_t * tmp_poly_arch = (uint32_t*) malloc( tmp_poly_arch_length * sizeof(uint32_t));
   if(tmp_poly_arch == NULL){
     fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
     return -4;
   }
 
   tmp_poly_arch[0]=Nmono1;

   tmp_poly_arch[1] = get_info_bl(1+ Nmono1 +Nmono2+2, get_particle_type(poly->poly_arch[1]));
    
    for(unsigned int i=2; i <= Nmono1-2+1; i++)
     {	 
       tmp_poly_arch[ i ] = get_info_bl(1+ Nmono1+Nmono2+3,get_particle_type(poly->poly_arch[i]));
     }
   tmp_poly_arch[Nmono1] = get_info_bl(1+ Nmono1 +Nmono2+1, get_particle_type(poly->poly_arch[Nmono1])); 
   
   tmp_poly_arch[Nmono1+1]=Nmono2;
   
   unsigned int tmp_start=Nmono1+2;
  
   for(unsigned int j=0;j<=n_crosslink;j++){
      unsigned int tmp_start_short=1;
     tmp_poly_arch[tmp_start] = get_info_bl(1+ Nmono1 +Nmono2+2, get_particle_type(poly->poly_arch[tmp_start_short]));
     
     for(unsigned int i=tmp_start+1; i <= tmp_start+Nmono1-2; i++)
       {
	 tmp_start_short++;
	 tmp_poly_arch[ i ] = get_info_bl(1+ Nmono1+Nmono2+3, get_particle_type(poly->poly_arch[tmp_start_short]));
       }
       tmp_poly_arch[tmp_start+ Nmono1-2+1] = get_info_bl(1+ Nmono1 +Nmono2+1,get_particle_type(poly->poly_arch[tmp_start_short+1]));    
       tmp_start=tmp_start+Nmono1-2+2;     
       }
   
   //Last monomer
   int end,offset,bond_type,info;
   bond_type = HARMONIC;
   end = 1;
   offset = -1;
   info=get_info(offset, bond_type, end);
   tmp_poly_arch[ 1+Nmono1+Nmono2+ 1 ] = info;
   printf("tmp_poly_arch[1+Nmono1+Nmono2+1]");
    //First monomer
    end = 1;
    offset = 1;
    info = get_info(offset, bond_type, end);
    tmp_poly_arch[ 1+Nmono1+Nmono2+ 2 ] = info;
    //Middle monomers A
    end = 0;
    offset = -1;
    info = get_info(offset, bond_type, end);
    tmp_poly_arch[ 1+Nmono1+Nmono2+ 3 ] = info;
    
    end = 1;
    offset = 1;
    info = get_info(offset, bond_type, end);
    tmp_poly_arch[1+Nmono1+Nmono2 + 4] = info;

    status=H5Ldelete(file_id,"/parameter/poly_arch",H5P_DEFAULT);

    dataset_id = H5Dcreate2(file_id, "/parameter/poly_arch", H5T_NATIVE_INT32,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);       
    status=H5Dwrite(dataset_id,H5T_NATIVE_INT32,H5S_ALL,H5S_ALL,plist_id,tmp_poly_arch);
    status = H5Dclose(dataset_id);

     /*change poly_type_offset*/
    hsize_t     dims_poly_type_offset[1];
    dims_poly_type_offset[0] = 2;
   
    int *tmp_offset =  realloc(poly->poly_type_offset, poly->n_poly_type * sizeof(unsigned int) );
    tmp_offset[1]=Nmono1+1;

    dataspace_id = H5Screate_simple(1, dims_poly_type_offset, NULL);
    status=H5Ldelete(file_id,"/parameter/poly_type_offset",H5P_DEFAULT);
    dataset_id = H5Dcreate2(file_id, "/parameter/poly_type_offset", H5T_NATIVE_INT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status=H5Dwrite(dataset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,plist_id,tmp_offset);
    status = H5Dclose(dataset_id);

    /*change poly_arch_len*/
    dataset_id=H5Dopen(file_id, "/parameter/poly_arch_length", H5P_DEFAULT);
    status=H5Dwrite(dataset_id,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,plist_id,&tmp_poly_arch_length);
    status = H5Dclose(dataset_id);  

    poly->n_poly_type=tmp_n_poly_type;
    poly->poly_arch=tmp_poly_arch;
    poly->poly_type=tmp_poly_type;
    poly->poly_type_offset=tmp_offset;
    poly->poly_arch_length=tmp_poly_arch_length;
    b->N_polymers=tmp_n_polymers;	

    status = H5Fclose(file_id);

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

  uint32_t*const poly_arch_tmp = (uint32_t*) malloc( poly_arch_len * sizeof(uint32_t));
  if(poly_arch_tmp == NULL){

    fprintf(stderr, "ERROR: Malloc poly_arch_tmp %s:%d\n", __FILE__, __LINE__);
    return -1;
  }
  status = read_hdf5(file_id,"/parameter/poly_arch",H5T_NATIVE_INT32,plist_id,poly_arch_tmp);
  HDF5_ERROR_CHECK2(status,"poly_arch");
  
  /*load the poly_arch information of monomerN*/

  int offset_poly=poly->poly_type_offset[1];
  for(int i=0;i<polyN-1;i++){
    offset_poly=offset_poly+b->number_of_beads[i];
  }

  uint32_t polyarch=poly_arch_tmp[offset_poly+monomerN]; 
  int off_bond=get_bondlist_offset(polyarch);
  unsigned int typeN=get_particle_type(poly_arch_tmp[offset_poly+monomerN]);
  int bond_type = HARMONIC;

  do  {
    unsigned int end=get_end(poly_arch_tmp[off_bond]);
    poly_arch_len++;
    if(end==1) break;  
    off_bond++;
    
  }  while(0==0);

  /*load the poly_arch information of monomerM*/
  int offset_polyM=poly->poly_type_offset[1];
  for(int i=0;i<polyM-1;i++){
    offset_polyM=offset_polyM+b->number_of_beads[i];
  }
  uint32_t polyarchM=poly->poly_arch[offset_polyM+monomerM];

  int off_bondM=get_bondlist_offset(polyarchM);
  unsigned int typeM=get_particle_type(poly_arch_tmp[offset_polyM+monomerM]);
  
  unsigned int poly_arch_lenM=poly_arch_len+1;
  unsigned int poly_arch_lenoldM=poly_arch_lenM;
  
  do  {
    unsigned int end=get_end(poly_arch_tmp[off_bondM]); 
    poly_arch_lenM++;
    if(end==1) break;  
    off_bondM++;    
  }  while(0==0);

  
  int bond_exist[poly_arch_len-poly_arch_lenold];//check if the bond already exsists;
  
  off_bond=get_bondlist_offset(polyarch);
  off_bondM=get_bondlist_offset(polyarchM);
  uint32_t * tmp_poly_archM = (uint32_t*) malloc((poly_arch_lenM+1) * sizeof(uint32_t));
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
  
  for(int i=0;i<poly_arch_len-poly_arch_lenold;i++){
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
   
    dims_poly_arch[0]=poly_arch_lenM+1;
  
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
      }
   /* Close the file. */
   status = H5Fclose(file_id);
   free(tmp_poly_archM);
   free(poly_arch_tmp);
   return 0;
}

int Prob_bonds(const char* const filename,bead_data*b,polymer_data*poly){

  herr_t status;

  hid_t plist_id=H5Pcreate(H5P_FILE_ACCESS);
  
  hid_t file_id=H5Fopen(filename,H5F_ACC_RDONLY,plist_id);
  HDF5_ERROR_CHECK2(file_id,"open file");
  H5Pclose(plist_id);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  double sum=0;
  double average;
  unsigned int Nref=0;
  status = read_hdf5(file_id,"/parameter/reference_Nbeads",H5T_NATIVE_UINT32,plist_id,&Nref);
  srand(time(NULL));
  int bond_polymer[b->old_N_polymers];
    
  int first=0;
  status = H5Fclose(file_id);

  double P0=0.01;
  double bfactor=1/sqrt(Nref);

  /*int max_number_bond=10;
  static int **r;
  r = (int **) malloc(b->old_N_polymers * sizeof(int *));
  for(int i=0;i<b->old_N_polymers;i++){
    r[i]=(int*) malloc(max_number_bond * sizeof(int));
  }
  for(int i=0;i<b->old_N_polymers;i++){
    for(int j=0;j<max_number_bond;j++){      
    r[i][j]=-1;
    }
    }*/
  
  for(int Poly_i=0;Poly_i<b->old_N_polymers;Poly_i++){
    int bond_i=0;
    bond_polymer[Poly_i]=0;
    for(int Mono_i=0;Mono_i<b->number_of_beads[Poly_i];Mono_i++){
      int second=0;
      for(int i=0; i<Poly_i;i++){
	second=second+b->number_of_beads[Poly_i];
      }
      for(int Poly_j=Poly_i;Poly_j<b->old_N_polymers;Poly_j++){
	for(int Mono_j=0;Mono_j<b->number_of_beads[Poly_j];Mono_j++){

	  double distance_sq=(b->beads[first+Mono_i].x-b->beads[second+Mono_j].x)*(b->beads[first+Mono_i].x-b->beads[second+Mono_j].x)+(b->beads[first+Mono_i].y-b->beads[second+Mono_j].y)*(b->beads[first+Mono_i].y-b->beads[second+Mono_j].y)+(b->beads[first+Mono_i].z-b->beads[second+Mono_j].z)*(b->beads[first+Mono_i].z-b->beads[second+Mono_j].z);
	  double Prob=P0*exp(-1.5/bfactor/bfactor*distance_sq);        
	  double random= (double)rand()/(double)RAND_MAX;
	 	  
	  if (random<Prob){	  	    	    
	    int crosslink=cross_link_monomer(filename,b,poly,Poly_i+1,Mono_i+1,Poly_j+1,Mono_j+1);	  
	    bond_polymer[Poly_i]++;
	  }
	}
	second=second+b->number_of_beads[Poly_j];		
      }
    }    
    printf("%i ,%i \n",Poly_i,bond_polymer[Poly_i]);
    sum=sum+bond_polymer[Poly_i];
    first=first+b->number_of_beads[Poly_i];    	   
  }
  average=sum/b->old_N_polymers;
  printf("number of polymers %i\n",b->old_N_polymers);
  printf("average %f \n",average);
  return 0;
}

int separate_molecule(bead_data*b,polymer_data*poly)
{
  
  int current_offset=poly->poly_type_offset[1];
  int number_molecule=0;
  int *checked_polymer;
  int number_of_beads=(int16_t) b->number_of_beads[0];
  int current_monomer=poly->poly_type_offset[1]+1+poly_start*number_of_beads;
  checked_polymer= (int*) malloc(b->old_N_polymers*sizeof(int));
  
  for(int a=0;a<b->old_N_polymers;a++){
    checked_polymer[a]=-1;    
  }
  
  checked_polymer[0]=0;

  unsigned int tmp_n_poly_type=poly->n_poly_type;
  
  for(int Poly_i=0;Poly_i<b->old_N_polymers;Poly_i++){   
    if(checked_polymer[Poly_i]==0){
      continue;      
    }
   
    int finish=0;
  
    int exist=0;
    int *loop_member;
    loop_member= (int*) malloc(b->old_N_polymers*number_of_beads*sizeof(int));
    int loop_member_index=0;
    int nodes_bond[100][10];
    int bond_nodes_checked[100];
    int number_of_nodes=0;
    int bond_number_total[100];
    int nodes_offset[100];
    int nodes_list[100];
    int last_node;
      
    current_monomer=current_monomer=poly->poly_type_offset[1]+1+Poly_i*number_of_beads;
    int current_poly_arch=poly->poly_arch[current_monomer];

    while(finish==0){
      for(int loop_index_tmp=0;loop_index_tmp<loop_member_index;loop_index_tmp++){
	if(loop_member[loop_index_tmp]==current_monomer){
	  /*deal with loop*/
	}
      }
      if(exist==1)
	break;
	
      loop_member[loop_member_index]=current_monomer;
      loop_member_index++;
      //printf("current %i,%i,",current_monomer,loop_member_index);
      current_poly_arch=poly->poly_arch[current_monomer];
      int start_offset_bond=get_bondlist_offset(current_poly_arch);
      int start_offset_bondtmp=start_offset_bond;
      int bonds_of_monomer;	    
      int bond_number=0;
      unsigned int crosslink=0;
      int if_end_chain=0;
      unsigned int number_near=0;
      do  {
	bonds_of_monomer=poly->poly_arch[start_offset_bond];
	int end=get_end(bonds_of_monomer);
	nodes_bond[number_of_nodes][bond_number]=get_offset(bonds_of_monomer);
	if(get_offset(bonds_of_monomer)!=1&&get_offset(bonds_of_monomer)!=-1){
	  crosslink=1;		  
	}		
	if(get_offset(bonds_of_monomer)==1||get_offset(bonds_of_monomer)==-1){
	  number_near++;
	}
	bond_number++;
	if(end==1) break;
	start_offset_bond++;

      }while(0==0);
  
      if(number_near==1&&crosslink==0&&get_offset(bonds_of_monomer)==1)
	if_end_chain=1;
	      
      if(number_near==1&&crosslink==0&&get_offset(bonds_of_monomer)==-1)
	if_end_chain=-1;
	      
      if(crosslink==0&&number_near==2){
	//printf("found a middle \n");
	current_monomer++;
	continue;
      }	      
      if(if_end_chain==-1){
	int ii=1;
	printf("found an end ,");
	//printf("number of nodes %i, %i, ", number_of_nodes,ii);
	//printf("%i, %i \n", bond_number_total[number_of_nodes-ii],bond_nodes_checked[number_of_nodes-ii]);
	while(number_of_nodes-ii>=0){
	  if(bond_number_total[number_of_nodes-ii]!=bond_nodes_checked[number_of_nodes-ii])
	    break;
	  else
	    ii++;		  
	}                 
	if(number_of_nodes-ii<0){
	  finish=1;
	  //break;
	}

	if(finish==1){
	  number_of_molecule++;
	  /*change hdf5*/

	  /*change n_poly_type */
 
	  hsize_t     dims_n_poly_type[1];
	  
	  tmp_n_poly_type++;
	  
	  dims_n_poly_type[0] = 1;
  
	  dataspace_id = H5Screate_simple(1, dims_n_poly_type, NULL);  
	  status=H5Ldelete(file_id,"/parameter/n_poly_type",H5P_DEFAULT);  
	  dataset_id = H5Dcreate2(file_id, "/parameter/n_poly_type", H5T_NATIVE_UINT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
	  status=H5Dwrite(dataset_id,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,plist_id,&tmp_n_poly_type);  
	  status = H5Dclose(dataset_id);
	  
	  /*change poly_type_offset*/
	  
	  hsize_t     dims_poly_type_offset[1];
	  dims_poly_type_offset[0] = 1+number_of_molecule;
   
	  int *tmp_offset =  realloc(poly->poly_type_offset, (1+number_of_molecule) * sizeof(unsigned int) );
    tmp_offset[number_of_molecule]=current_offset;
    current_offset=current_offset+???;
    dataspace_id = H5Screate_simple(1, dims_poly_type_offset, NULL);
    status=H5Ldelete(file_id,"/parameter/poly_type_offset",H5P_DEFAULT);
    dataset_id = H5Dcreate2(file_id, "/parameter/poly_type_offset", H5T_NATIVE_INT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status=H5Dwrite(dataset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,plist_id,tmp_offset);
    status = H5Dclose(dataset_id);

    /*change number_of_beads*/

    int old_max_n_beads=b->max_n_beads;

    unsigned int *tmp_number_of_beads[number_of_molecule];
    unsigned int new_n_beads=???;
  
    if( b-> max_n_beads< new_n_beads)
      b->max_n_beads = new_n_beads;
    
    for(unsigned int i=0;i<number_of_molecule;i++){
      tmp_number_of_beads[number_of_molecule-1]=b->number_of_beads[i];
    }
     tmp_number_of_beads[number_of_molecule]=new_n_beads;

     /*change size of the beads box*/
   hsize_t     dims_beads[2];

  dims_beads[0] = number_of_molecule;  
  dims_beads[1] = b->max_n_beads;
  
  Monomer *tmp_beads = (Monomer * const) malloc(b->max_n_beads*number_of_molecule * sizeof(Monomer));

  const hid_t memtype = get_monomer_memtype();
  
  dataspace_id = H5Screate_simple(2, dims_beads, NULL);
  
  /*for(unsigned int i=0;i<1;i++){
    for(unsigned int j=0;j<tmp_number_of_beads[i];j++){
      tmp_beads[i*b->max_n_beads+j].x=b->beads[i*old_max_n_beads+j].x;
      tmp_beads[i*b->max_n_beads+j].y=b->beads[i*old_max_n_beads+j].y;
      tmp_beads[i*b->max_n_beads+j].z=b->beads[i*old_max_n_beads+j].z;
    }
  }

  for(unsigned int i=1;i<tmp_n_polymers;i++){
    for(unsigned int j=0;j<tmp_number_of_beads[i];j++){
      tmp_beads[i*b->max_n_beads+j].x=b->beads[(i+n_crosslink)*old_max_n_beads+j].x;
      tmp_beads[i*b->max_n_beads+j].y=b->beads[(i+n_crosslink)*old_max_n_beads+j].y;
      tmp_beads[i*b->max_n_beads+j].z=b->beads[(i+n_crosslink)*old_max_n_beads+j].z;
    }    
    }*/

  status=H5Ldelete(file_id,"/beads",H5P_DEFAULT);
  dataset_id = H5Dcreate2(file_id, "/beads", memtype,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status=H5Dwrite(dataset_id,memtype,H5S_ALL,H5S_ALL,plist_id,tmp_beads);
  status = H5Dclose(dataset_id);


	}
      
	loop_member_index=nodes_offset[number_of_nodes-ii];
	current_monomer=nodes_list[number_of_nodes-ii]+nodes_bond[number_of_nodes-ii][bond_nodes_checked[number_of_nodes-ii]];
	bond_nodes_checked[number_of_nodes-ii]++;
	number_of_nodes=number_of_nodes-ii+1;
	continue;
      }	      
      if(if_end_chain==1){
	//printf("found a start \n");
	current_monomer++;
	continue;
      }
	      
      if(crosslink==1){
	//printf("found a node,");
	int current_polymer=current_monomer/number_of_beads;
	//printf("current polymer",current_polymer);
	checked_polymer[current_polymer]=0;
	int meetmeet=0;
	int bond_i;
	int bond_ii=0;
	int bond_iii=0;
	for(bond_i=0;bond_i<bond_number;bond_i++){	
	  int meet=0;
	  if(get_offset(poly->poly_arch[start_offset_bondtmp+bond_i])!=-1){
	    for(bond_ii=0;bond_ii<bond_number_total[number_of_nodes-1];bond_ii++){
	      if(-nodes_bond[number_of_nodes-1][bond_ii]==get_offset(poly->poly_arch[start_offset_bondtmp+bond_i])&&last_node+1==loop_member_index){
		meet++;
		meetmeet++;
	      }
	    }
	    if(meet==0||get_offset(poly->poly_arch[start_offset_bondtmp+bond_i])==1){
	      nodes_bond[number_of_nodes][bond_iii]=get_offset(poly->poly_arch[start_offset_bondtmp+bond_i]);
	      bond_iii++;
	    }
	  }
	}
	if(meetmeet==1){
	  for(bond_i=0;bond_i<bond_number;bond_i++){
	    if(get_offset(poly->poly_arch[start_offset_bondtmp+bond_i])==-1){
	      nodes_bond[number_of_nodes][bond_iii]=get_offset(poly->poly_arch[start_offset_bondtmp+bond_i]);
	    }
	  }

	}
      
	nodes_list[number_of_nodes]=current_monomer;
	nodes_offset[number_of_nodes]=loop_member_index;;
	bond_number_total[number_of_nodes]=bond_iii;
	current_monomer=current_monomer+nodes_bond[number_of_nodes][0];
	bond_nodes_checked[number_of_nodes]=1;
	number_of_nodes++;
	last_node=loop_member_index;
      }	     	             
    }
  }
  free(poly_arch);
  return number_loop;

}
