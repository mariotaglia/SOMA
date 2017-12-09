#include <hdf5.h>
#include "readIn_lib.h"
#include "io.h"
#include <math.h>
#include <stdlib.h>
#include "separateMolecule_lib.h"
#include <string.h>



void cross_link_monomer1(int32_t** poly_arch_tmp_pointer,polymer_data*const poly, unsigned int poly_arch_len,unsigned int monomerN, unsigned int monomerM)
{
  unsigned int poly_arch_lenold=poly_arch_len;
  /*load the poly_arch information of monomerN*/
  int32_t* poly_arch_tmp= *poly_arch_tmp_pointer;
  int offset_poly=0;
  int32_t polyarch=poly_arch_tmp[offset_poly+monomerN];

  if(polyarch>4294967200){
    fprintf(stderr, "ERROR: Function get_bondlist_offset exited the range %s:%d\n", __FILE__, __LINE__);
    //return -2;
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
  int offset_polyM=0; 
  int bond_exist[poly_arch_len-poly_arch_lenold];//check if the bond already exsists;
  if(polyarch>4294967200){
    fprintf(stderr, "ERROR: Function get_bondlist_offset exited the range %s:%d\n", __FILE__, __LINE__);
    //return -2;
  }
  
  off_bond=get_bondlist_offset(polyarch);
  int32_t *tmp_poly_archM;
  tmp_poly_archM= malloc((poly_arch_len+1) * sizeof(int32_t));
  if(tmp_poly_archM == NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
  }
 
  memcpy(tmp_poly_archM,poly_arch_tmp,poly_arch_lenold*sizeof(int32_t));
  for(unsigned int i=poly_arch_lenold;i<poly_arch_len;i++){
    int info;
    int offset=get_offset(poly_arch_tmp[off_bond]);
   
    bond_exist[i-poly_arch_lenold]=offset;
    info=get_info(offset, bond_type, 0);
    tmp_poly_archM[i] = info;
    off_bond++;
  }

  int bond_match=0;
  int offset_to_monomerM=offset_polyM+monomerM-offset_poly-monomerN;

  for(unsigned int i=0;i<poly_arch_len-poly_arch_lenold;i++){
    if(offset_to_monomerM==bond_exist[i])
      bond_match++;
  }

  if(bond_match==0 && offset_to_monomerM!=0){
    int info;
    info=get_info(offset_to_monomerM, bond_type,1);
    tmp_poly_archM[poly_arch_len]=info;
    tmp_poly_archM[offset_poly+monomerN] = get_info_bl(poly_arch_lenold, typeN);
     
    poly_arch_len++;   
    *poly_arch_tmp_pointer = tmp_poly_archM;
    poly->poly_arch_length=poly_arch_len;
  }
  
  free(poly_arch_tmp);
}

/////////////
//next one///
////////////

int separate_molecule(const char* const filename,bead_data*b,polymer_data*poly){  
  //int32_t* poly_arch_new=poly->poly_arch;
  int total_length=0;
  int total_beads=0;
  int *checked_monomer;
  unsigned int sequence=0;
  
  for(unsigned int i=0;i<b->old_N_polymers;i++){
    sequence=sequence+b->number_of_beads[i];
  }
  
  int current_monomer=poly->poly_type_offset[1]+1;
  
  checked_monomer= (int*) malloc(sequence*sizeof(int));
  
  for(unsigned int a=0;a<sequence;a++){
    checked_monomer[a]=-1;    
  }
  int molecule_index=-1;
 
  int store_all_crosslink[sequence][8];
  for(unsigned int a=0;a<sequence;a++){   
    for(int b =0;b<8;b++){
      store_all_crosslink[a][b]=0;
    }
  }
  int *loop_member_length;
  loop_member_length= (int*) malloc(sequence*sizeof(int));
  if(loop_member_length == NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
    return -4;
  }  
  for(unsigned int mono_i=0;mono_i<sequence;mono_i++){
  
    if(checked_monomer[mono_i]!=-1)
      continue;
    molecule_index++;
    int molecule_length=0;
    int normal_node=0;
    int finish=0;
    int *loop_member;
    loop_member= (int*) malloc(sequence*sizeof(int));
    if(loop_member == NULL){
      fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
      return -4;
    }
    int loop_member_index=0;
    unsigned int esti=sequence;
    int bond_total[8];
    for(int bbb=0;bbb<8;bbb++){
      bond_total[bbb]=0;
    }
    int checked_bonds[2];
    checked_bonds[0]=0;
    checked_bonds[1]=0;
    int **checked_bonds_node=(int **)malloc(esti*sizeof(int*));
    for(unsigned int aa=0;aa<esti;aa++){
      checked_bonds_node[aa]=(int *)malloc(8*sizeof(int));
    }
    int number_of_nodes=0;
       
    int* bond_nodes_checked= (int*) malloc(esti*sizeof(int));
    
    for(unsigned int a=0;a<esti;a++){
      bond_nodes_checked[a]=0;
      for(int b =0;b<8;b++){
	checked_bonds_node[a][b]=0;
      }
    }

    int* bond_number_total;
    bond_number_total= (int*) malloc(esti*sizeof(int));
    if( bond_number_total== NULL){
      fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
      return -4;
    }
    int* nodes_offset;
    nodes_offset= (int*) malloc(esti*sizeof(int));
    if(nodes_offset == NULL){
      fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
      return -4;
    }
    int* nodes_list;
    nodes_list= (int*) malloc(esti*sizeof(int));
    if(nodes_list == NULL){
      fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
      return -4;
    }
    for(unsigned int bbb=0;bbb<esti;bbb++){
      bond_number_total[bbb]=0;
      nodes_offset[bbb]=0;
      nodes_list[bbb]=0;
    }
    current_monomer=poly->poly_type_offset[1]+1+mono_i;
    int current_poly_arch=poly->poly_arch[current_monomer];
    while(finish==0){
      molecule_length++;
      checked_monomer[current_monomer-poly->poly_type_offset[1]-1]=molecule_index;
      for(int i=0;i<8;i++){
	bond_total[i]=0;
      }
      int side_node=0;
      int exist=0;
      for(int loop_index_tmp=0;loop_index_tmp<loop_member_index&&exist==0;loop_index_tmp++){
	if(loop_member[loop_index_tmp]==current_monomer){
	  exist=1;
	  side_node=1;	  
	}
      }   
      loop_member_index++;
      loop_member[loop_member_index]=current_monomer;      
      current_poly_arch=poly->poly_arch[current_monomer];

      if(current_poly_arch>4294967200){
	fprintf(stderr, "ERROR: Function get_bondlist_offset exited the range %s:%d\n", __FILE__, __LINE__);
	return -2;
      }
      int start_offset_bond=get_bondlist_offset(current_poly_arch);
      int bonds_of_monomer;	    
      int bond_number=0;
      int if_end_chain=0;
      int offset_tmp=0;
      do{
	bonds_of_monomer=poly->poly_arch[start_offset_bond];
	int end=get_end(bonds_of_monomer);
	offset_tmp=get_offset(bonds_of_monomer);
	bond_total[bond_number]=get_offset(bonds_of_monomer);		
	bond_number++;
	if(end==1) break;
	start_offset_bond++;
      }while(0==0);      
      if(bond_number==1){
	  store_all_crosslink[current_monomer-poly->poly_type_offset[1]-1][0]=0;
	  store_all_crosslink[current_monomer-poly->poly_type_offset[1]-1][1]=offset_tmp;
      }    
      if(bond_number==1&&checked_bonds[0]!=0){
	int ii=1;       
	while(number_of_nodes-ii>=0){
	  if(bond_number_total[number_of_nodes-ii]!=bond_nodes_checked[number_of_nodes-ii])
	    break;
	  else{	    
	    ii++;
	  }
	}	
	if(number_of_nodes-ii<0){
	  finish=1;
	  continue;
	}
	current_monomer=nodes_list[number_of_nodes-ii];
	checked_bonds[0]=0;
	number_of_nodes=number_of_nodes-ii;

	for(unsigned int node_tmp=number_of_nodes+1;node_tmp<esti;node_tmp++){
	  for(int tmp=0;tmp<8;tmp++){
	    checked_bonds_node[node_tmp][tmp]=0;	
	  }
	 
	  bond_nodes_checked[node_tmp]=0;
	  bond_number_total[node_tmp]=0;
	}       
	normal_node=1;
	continue;
      }
      if(bond_number==2){
	if(bond_total[0]+bond_total[1]==0){
	  store_all_crosslink[current_monomer-poly->poly_type_offset[1]-1][0]=0;
	  store_all_crosslink[current_monomer-poly->poly_type_offset[1]-1][1]=-2;
	}
	else{
	  store_all_crosslink[current_monomer-poly->poly_type_offset[1]-1][0]=bond_total[0];
	  store_all_crosslink[current_monomer-poly->poly_type_offset[1]-1][1]=bond_total[1];
	}
      }      
      if((bond_number==2&&checked_bonds[0]!=0&&side_node!=1)||(bond_number==1&&checked_bonds[0]==0)){
	//store_all_crosslink[current_monomer-poly->poly_type_offset[1]-1][0]=0;
	int bond_j;
	normal_node=0;
	for(bond_j=0;bond_j<2;bond_j++){
	  int meetj=0; 
	  for(int bond_jj=0;bond_jj<2;bond_jj++){
	    if(bond_total[bond_j]==checked_bonds[bond_jj]){
	      meetj=1;
	      break;
	    }
	  }
	  if(meetj==0){
	    break;
	  }
	}
	current_monomer=bond_total[bond_j]+current_monomer;
	checked_bonds[0]=-bond_total[bond_j];
	continue;	  
      }
        
      if(bond_number==2&&checked_bonds[0]==0){
	side_node=1;
	normal_node=1;
      }
      
      if(bond_number>2||side_node==1){
	for(int bbb=0;bbb<bond_number;bbb++){	  
	  store_all_crosslink[current_monomer-poly->poly_type_offset[1]-1][bbb]=bond_total[bbb];
	  
	}
	int meet_node=0;
	int aa=1;
	int found_the_node=0;
	if(normal_node==0)
	  aa=0;
	for(int i=0;i<number_of_nodes-aa;i++){
	  if(nodes_list[i]==current_monomer){
	    meet_node=i;
	    found_the_node=1;
	    for(int j=0;j<bond_number;j++){
	      checked_bonds_node[number_of_nodes][j]=checked_bonds_node[i][j];
	    }
	    break;
	  }
	}

	int j=0;	  
	while(checked_bonds_node[number_of_nodes][j]!=0){
	  j++;
	}
	bond_nodes_checked[number_of_nodes]=j;
	if(normal_node==0){
	  checked_bonds_node[number_of_nodes][j]=checked_bonds[0];
	  bond_nodes_checked[number_of_nodes]++;
	  j++;
	}

	if(found_the_node==1){
	  bond_nodes_checked[meet_node]=bond_nodes_checked[number_of_nodes];
	  for(int jjj=0;jjj<bond_number;jjj++){
	    checked_bonds_node[meet_node][jjj]=checked_bonds_node[number_of_nodes][jjj];
	  }
	}
	
	if(bond_number==bond_nodes_checked[number_of_nodes]){
	  if_end_chain=1;    
	}
	else{	  	    
	  int bond_i;
	  int bond_ii=0;	 
	  for(bond_i=0;bond_i<bond_number;bond_i++){
	    int meet=0;
	    for(bond_ii=0;bond_ii<bond_number;bond_ii++){
	      if(bond_total[bond_i]==checked_bonds_node[number_of_nodes][bond_ii]){
		meet=1;
		break;
	      }
	    }
	    if(meet==0)
	      break;	  
	  }
	      
	      
	  checked_bonds_node[number_of_nodes][j]=bond_total[bond_i];
	  nodes_list[number_of_nodes]=current_monomer;
	  nodes_offset[number_of_nodes]=loop_member_index;
	
	  normal_node=0;	  
	  bond_nodes_checked[number_of_nodes]++;
	  
	  current_monomer=current_monomer+bond_total[bond_i];
	  
	  if(found_the_node==1){
	    bond_nodes_checked[meet_node]=bond_nodes_checked[number_of_nodes];
	    for(int jjj=0;jjj<bond_number;jjj++){
	      checked_bonds_node[meet_node][jjj]=checked_bonds_node[number_of_nodes][jjj];
	    }  
	  }
	  checked_bonds[0]=-bond_total[bond_i];
	}
      }      	  
      number_of_nodes++;
      bond_number_total[number_of_nodes-1]=bond_number;
      
      if(if_end_chain==1){
	normal_node=1;
	int ii=1;        
	while(number_of_nodes-ii>=0){
	  if(bond_number_total[number_of_nodes-ii]!=bond_nodes_checked[number_of_nodes-ii])
	    break;
	  else{	    
	    ii++;
	  }
	}
	if(number_of_nodes-ii<0){
	  finish=1;
	  break;
	}
	current_monomer=nodes_list[number_of_nodes-ii];
	checked_bonds[0]=0;
	number_of_nodes=number_of_nodes-ii;

	for(unsigned int node_tmp=number_of_nodes+1;node_tmp<esti;node_tmp++){
	  for(int tmp=0;tmp<8;tmp++){
	    checked_bonds_node[node_tmp][tmp]=0;
	  }
	  bond_nodes_checked[node_tmp]=0;
	  bond_number_total[node_tmp]=0;
	}       
	continue;
      }      
    }
    free(checked_bonds_node);
    loop_member_length[molecule_index]=molecule_length;
    free(loop_member);
    free(nodes_list);
    free(nodes_offset);
    free(bond_nodes_checked);
    free(bond_number_total);
 }

  int lastone=molecule_index;
  int start;
  int32_t **tmp_poly_arch=(int32_t **)malloc(molecule_index*sizeof(int32_t*));
  
  int *tmp_poly_offset=(int *)malloc(molecule_index*sizeof(int));
    
  Monomer **tmp_beads=(Monomer ** const)malloc(molecule_index*sizeof(Monomer*));
  int max_length=0;
  for(int tmp_molecule_index=0; tmp_molecule_index<lastone/*molecule_index*/;tmp_molecule_index++){
    printf("now we are at molecule %i\n",tmp_molecule_index);
    int* mapping_list;
    mapping_list= (int*) malloc(sequence*sizeof(int));
    if(mapping_list== NULL){
      fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
      return -4;
    }
     for(unsigned int c=0;c<sequence;c++){
      mapping_list[c]=0;
      }
   
    
    tmp_poly_arch[tmp_molecule_index]= (int32_t*) malloc((total_length+loop_member_length[tmp_molecule_index]*4)*sizeof(int32_t));
    if(tmp_poly_arch[tmp_molecule_index] == NULL){
      fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
      return -4;
    }
    if(tmp_molecule_index!=0){
      memcpy(tmp_poly_arch[tmp_molecule_index],tmp_poly_arch[tmp_molecule_index-1],total_length*sizeof(int32_t));
      free(tmp_poly_arch[tmp_molecule_index-1]);
    }
    start=total_length;

    for(int iii=total_length;iii<total_length+loop_member_length[tmp_molecule_index]*4;iii++){
      tmp_poly_arch[tmp_molecule_index][iii]=0;
    }
  
    tmp_beads[tmp_molecule_index]= (Monomer* const) malloc(loop_member_length[tmp_molecule_index]*sizeof(Monomer));
    if(tmp_beads[tmp_molecule_index] == NULL){
      fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
      return -4;
    }

   
    int tmp_chain_index=start;
    for(unsigned int tmp_monomer_index=0;tmp_monomer_index<sequence;tmp_monomer_index++){
      if(checked_monomer[tmp_monomer_index]==tmp_molecule_index){
	mapping_list[tmp_monomer_index]=tmp_chain_index;      
	tmp_beads[tmp_molecule_index][tmp_chain_index-start].x=b->beads[tmp_monomer_index].x;
	tmp_beads[tmp_molecule_index][tmp_chain_index-start].y=b->beads[tmp_monomer_index].y;
	tmp_beads[tmp_molecule_index][tmp_chain_index-start].z=b->beads[tmp_monomer_index].z;
	tmp_chain_index++;
      }
    }
    tmp_poly_arch[tmp_molecule_index][start]=tmp_chain_index-start;
    //Last monomer
    int end,offset,bond_type,info;
    bond_type = HARMONIC;
    end = 1;
    offset = -1;
    info=get_info(offset, bond_type, end);
    tmp_poly_arch[tmp_molecule_index][tmp_chain_index+1] = info;
 
    //First monomer
    end = 1;
    offset = 1;
    info = get_info(offset, bond_type, end);
    tmp_poly_arch[tmp_molecule_index][tmp_chain_index+2] = info;
    //Middle monomers A
    end = 0;
    offset = -1;
    info = get_info(offset, bond_type, end);
    tmp_poly_arch[tmp_molecule_index][tmp_chain_index+3] = info;
  
    end = 1;
    offset = 1;
    info = get_info(offset, bond_type, end);
    tmp_poly_arch[tmp_molecule_index][tmp_chain_index+4] = info;

    int tmp_iii=1+start;

    for(unsigned int tmp_monomer_index=0;tmp_monomer_index<sequence;tmp_monomer_index++){      
      if(checked_monomer[tmp_monomer_index]==tmp_molecule_index){
	
	if(store_all_crosslink[tmp_monomer_index][0]==0){
	  if(store_all_crosslink[tmp_monomer_index][1]==1){
	    tmp_poly_arch[tmp_molecule_index][tmp_iii]=get_info_bl(tmp_chain_index+2, get_particle_type(poly->poly_arch[tmp_monomer_index+poly->poly_type_offset[1]+1]));
	  }
	  if(store_all_crosslink[tmp_monomer_index][1]==-1){
	    tmp_poly_arch[tmp_molecule_index][tmp_iii]=get_info_bl(tmp_chain_index+1, get_particle_type(poly->poly_arch[tmp_monomer_index+poly->poly_type_offset[1]+1]));
	  }
	  if(store_all_crosslink[tmp_monomer_index][1]==-2){
	    tmp_poly_arch[tmp_molecule_index][tmp_iii]=get_info_bl(tmp_chain_index+3, get_particle_type(poly->poly_arch[tmp_monomer_index+poly->poly_type_offset[1]+1]));
	  }	 
	}
	else{	  
	  int nachbar=0;
	  int which_side=0;
	  int jjj=0;
	  while(store_all_crosslink[tmp_monomer_index][jjj]!=0){	    
	    if(store_all_crosslink[tmp_monomer_index][jjj]==1||store_all_crosslink[tmp_monomer_index][jjj]==-1){
	      which_side=store_all_crosslink[tmp_monomer_index][jjj];
	      nachbar++;
	    }
	    jjj++;
	  }	
	  if(nachbar==1){
	    if(which_side==-1)
	      tmp_poly_arch[tmp_molecule_index][tmp_iii]=get_info_bl(tmp_chain_index+1, get_particle_type(poly->poly_arch[tmp_monomer_index+poly->poly_type_offset[1]+1]));
	    if(which_side==1)
	      tmp_poly_arch[tmp_molecule_index][tmp_iii]=get_info_bl(tmp_chain_index+2, get_particle_type(poly->poly_arch[tmp_monomer_index+poly->poly_type_offset[1]+1]));
	  }
	  
	  if(nachbar==2)	   
	    tmp_poly_arch[tmp_molecule_index][tmp_iii]=get_info_bl(tmp_chain_index+3, get_particle_type(poly->poly_arch[tmp_monomer_index+poly->poly_type_offset[1]+1]));	  
	}
	tmp_iii++;
      }
    }

    poly->poly_arch_length=tmp_chain_index+5;

    int tmp_ii=start; 
    for(unsigned int tmp_monomer_index=0;tmp_monomer_index<sequence;tmp_monomer_index++){  
      if(checked_monomer[tmp_monomer_index]==tmp_molecule_index){	
	if(store_all_crosslink[tmp_monomer_index][0]!=0){
	  //int nachbar=0;
	  int jj=0;
	  while(store_all_crosslink[tmp_monomer_index][jj]!=0){	    
	    if(store_all_crosslink[tmp_monomer_index][jj]!=-1&&store_all_crosslink[tmp_monomer_index][jj]!=1){
	      cross_link_monomer1(&tmp_poly_arch[tmp_molecule_index],poly,poly->poly_arch_length,tmp_ii+1,mapping_list[store_all_crosslink[tmp_monomer_index][jj]+tmp_monomer_index]+1);        
	    }
	    jj++;
	  }
	}
	tmp_ii++;
      }
    }
 
    
    if(max_length<tmp_poly_arch[tmp_molecule_index][start])
      max_length=tmp_poly_arch[tmp_molecule_index][start];
    
    total_beads=total_beads+tmp_poly_arch[tmp_molecule_index][start];
    tmp_poly_offset[tmp_molecule_index]=total_length;
    total_length=poly->poly_arch_length;  
    free(mapping_list);
  }
  free(loop_member_length);
  free(checked_monomer);
  
 Monomer *beads_total= (Monomer* const) malloc(max_length*molecule_index*sizeof(Monomer));
  if(beads_total == NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
    return -4;
  }

  for (int i = 0; i < lastone/*molecule_index*/; i++)
    {
      const unsigned int N = tmp_poly_arch[lastone-1][tmp_poly_offset[i]];
      memcpy(beads_total+max_length*i, tmp_beads[i],N * sizeof(Monomer));
      free(tmp_beads[i]);
    }
  
  herr_t status;
  hid_t dataset_id, dataspace_id;
  hid_t plist_id=H5Pcreate(H5P_FILE_ACCESS);   
  hid_t file_id=H5Fopen(filename,H5F_ACC_RDWR,plist_id);
  HDF5_ERROR_CHECK2(file_id,"open file");
  H5Pclose(plist_id);
  plist_id=H5Pcreate(H5P_DATASET_XFER);

  dataset_id=H5Dopen(file_id, "/parameter/n_polymers", H5P_DEFAULT);   
  status=H5Dwrite(dataset_id,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,plist_id,&molecule_index);
  status = H5Dclose(dataset_id);

  
  hsize_t     dims_poly_archnew[1];
  dims_poly_archnew[0]=total_length;

  dataspace_id = H5Screate_simple(1, dims_poly_archnew, NULL);
   status=H5Ldelete(file_id,"/parameter/poly_arch",H5P_DEFAULT);
  dataset_id = H5Dcreate2(file_id, "/parameter/poly_arch", H5T_NATIVE_INT32,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);       
  status=H5Dwrite(dataset_id,H5T_NATIVE_INT32,H5S_ALL,H5S_ALL,plist_id,tmp_poly_arch[lastone-1]);
  status = H5Dclose(dataset_id);
  
  /*change poly_arch_len*/
  dataset_id=H5Dopen(file_id, "/parameter/poly_arch_length", H5P_DEFAULT);
  status=H5Dwrite(dataset_id,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,plist_id,&total_length);
  status = H5Dclose(dataset_id);

  /*change poly_type*/
  hsize_t     dims_poly_type[1];
  dims_poly_type[0]=molecule_index;

  int * tmp_poly_type=( int *)malloc(molecule_index*sizeof( int));
  for(int c=0;c<molecule_index;c++){
    tmp_poly_type[c]=c;
  }
  
  dataspace_id = H5Screate_simple(1, dims_poly_type, NULL);
  status=H5Ldelete(file_id,"/poly_type",H5P_DEFAULT);
  dataset_id = H5Dcreate2(file_id, "/poly_type", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status=H5Dwrite(dataset_id,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,plist_id,tmp_poly_type);
  status = H5Dclose(dataset_id);
  
  free(tmp_poly_type);

    /*change n_poly_type */
 
  unsigned int tmp_n_poly_type=molecule_index;
  dataset_id=H5Dopen(file_id, "/parameter/n_poly_type", H5P_DEFAULT);
  status=H5Dwrite(dataset_id,H5T_NATIVE_UINT,H5S_ALL,H5S_ALL,plist_id,&tmp_n_poly_type);  
  status = H5Dclose(dataset_id);

  
  /*change poly_type_offset */
  hsize_t     dims_poly_type_offset[1];
  dims_poly_type_offset[0] = molecule_index;
  
  dataspace_id = H5Screate_simple(1, dims_poly_type_offset, NULL);
  status=H5Ldelete(file_id,"/parameter/poly_type_offset",H5P_DEFAULT);
  dataset_id = H5Dcreate2(file_id, "/parameter/poly_type_offset", H5T_NATIVE_INT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status=H5Dwrite(dataset_id,H5T_NATIVE_INT,H5S_ALL,H5S_ALL,plist_id,tmp_poly_offset);
  status = H5Dclose(dataset_id);

  free(tmp_poly_offset);

  /*change beads*/
  hsize_t     dims_beads[2];
  dims_beads[0] =molecule_index;  
  dims_beads[1] =max_length;//change it to max
  const hid_t memtype = get_monomer_memtype();    
  dataspace_id = H5Screate_simple(2, dims_beads, NULL);
  status=H5Ldelete(file_id,"/beads",H5P_DEFAULT);
  dataset_id = H5Dcreate2(file_id, "/beads", memtype,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status=H5Dwrite(dataset_id,memtype,H5S_ALL,H5S_ALL,plist_id,beads_total);
  status = H5Dclose(dataset_id);

  status=H5Fclose(file_id);

  free(beads_total);

  return 0;
}

int independent_set(bead_data*b,polymer_data*poly){
  
  
  int poly_type=0;
  unsigned int sequence=poly->poly_arch[poly->poly_type_offset[poly_type]];
  /*for(unsigned int i=0;i<b->old_N_polymers;i++){
    sequence=sequence+b->number_of_beads[i];
    }*/
  int max_bond_number=0,max_bond=0;

  int* bond_number_total;
  bond_number_total= (int*) malloc(sequence*sizeof(int));
  if( bond_number_total== NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
    return -4;
  }
  memset(bond_number_total,0,sequence*sizeof(int));
  int **bonds_total=(int **)malloc(sequence*sizeof(int*));
  if(bonds_total == NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
    return -4;
  }
  
  for(int current_tmp=poly->poly_type_offset[poly_type]+1;current_tmp<sequence+poly->poly_type_offset[poly_type]+1;current_tmp++){
  
    int bonds_of_monomer=0, bond_number=0;
    int  current_poly_arch=poly->poly_arch[current_tmp];
    int start_offset_bond=get_bondlist_offset(current_poly_arch);
    // printf("current %i %i %i %i %i\n",current_poly_arch,start_offset_bond,current_tmp,sequence,b->old_N_polymers);
    do{
      bonds_of_monomer=poly->poly_arch[start_offset_bond];
      int end=get_end(bonds_of_monomer);		
      bond_number++;
      if(end==1) break;
      start_offset_bond++;
    }while(0==0); 

    bond_number_total[current_tmp-poly->poly_type_offset[poly_type]-1]=bond_number;
    bonds_total[current_tmp-poly->poly_type_offset[poly_type]-1]=(int *)malloc(bond_number*sizeof(int));
    memset(&bonds_total[current_tmp-poly->poly_type_offset[poly_type]-1][0],0,bond_number*sizeof(int));     
	  
    start_offset_bond=get_bondlist_offset(poly->poly_arch[current_tmp]);
    for(int tmp=0;tmp<bond_number;tmp++){      
      bonds_of_monomer=poly->poly_arch[start_offset_bond];
     
      bonds_total[current_tmp-poly->poly_type_offset[poly_type]-1][tmp]=get_offset(bonds_of_monomer)+current_tmp-poly->poly_type_offset[poly_type]-1;
      start_offset_bond++;
    }
    
    if(bond_number>max_bond_number){    
      max_bond_number=bond_number; 
      max_bond=current_tmp-poly->poly_type_offset[poly_type]-1;
    }  
 
  }  
  
  int* monomer_checked;
  monomer_checked= (int*) malloc(sequence*sizeof(int));
  if(monomer_checked== NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
    return -4;
  }
  memset(monomer_checked,0,sequence*sizeof(int));
  //printf("checked %i\n",monomer_checked[51400]);
  int current_monomer=max_bond;
  int** independent_sets;
  independent_sets= (int**) malloc((max_bond_number+1)*sizeof(int*));
  if( independent_sets== NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
    return -4;
  }
  int offset_set[max_bond_number+1]; //only members between offset_set and end_set are unchecked for neighbours  

  int end_set[max_bond_number+1];
  
  for(int aaa=0;aaa<max_bond_number+1;aaa++){
    end_set[aaa]=0;
    offset_set[aaa]=0;
  }
  for(int bb=0;bb<max_bond_number+1;bb++){
    independent_sets[bb]=(int*) malloc(sequence*sizeof(int));
    memset(&independent_sets[bb][0],0,sequence*sizeof(int)); 
  }
  int total_assigned_monomer=0;

  ///loop over all monomer
  for(int mono_i=0;mono_i<sequence;mono_i++){
    if(total_assigned_monomer==sequence)
      break;
    if(monomer_checked[mono_i]==-1)
      continue;
    monomer_checked[mono_i]=-1;
    current_monomer=mono_i;
    independent_sets[0][end_set[0]]=current_monomer;//write the first element in each set
    end_set[0]++;
    //monomer_checked[current_monomer]=-1;
    total_assigned_monomer++;
    for(int aa=1; aa<=bond_number_total[current_monomer];aa++){   
      independent_sets[aa][end_set[aa]]=bonds_total[current_monomer][aa-1];
      total_assigned_monomer++;
      monomer_checked[bonds_total[current_monomer][aa-1]]=-1;
      // if(mono_i==51400)
	// printf("mono %i\n",bonds_total[current_monomer][aa-1]);
      end_set[aa]++;
    }
    offset_set[0]=offset_set[0]+1;//neighbour of inde*[0][offset_set] are all found now 
    int current_set=1; //the one needs to be checked
    int wang=0;
    

    while(total_assigned_monomer<sequence){
      int chain_finished=0;
      for(int i_tmp=0;i_tmp<max_bond_number+1;i_tmp++){
	if(offset_set[i_tmp]!=end_set[i_tmp]){	  
	  chain_finished=1;
	}
      }    
       if(chain_finished==0){
	 break;
       }         
      
      int writein_set=current_set-1; //which set to put the new element into, namely the left one !! no problem, because this is the last current one
      if(current_set==0)
	writein_set=max_bond_number;
     
      if(end_set[current_set]==offset_set[current_set]){
	current_set++;
	if(current_set>max_bond_number)
	  current_set=current_set-max_bond_number-1;
	continue;
      }    
     
      for(int member_set=offset_set[current_set];member_set<end_set[current_set];member_set++){     
	current_monomer=independent_sets[current_set][member_set]; 
	//	if(current_monomer>51394&&current_monomer<51406)
	//  printf("current %i %i %i %i %i \n",current_monomer,bond_number_total[current_monomer],bonds_total[current_monomer][0],bonds_total[current_monomer][1],current_set);
	int number_tmp=0;
	if(monomer_checked[bonds_total[current_monomer][0]]!=-1){
	  
	  int found=0;
	  while(found==0){ 
	    int meet=0;
	    for(int cc=offset_set[writein_set];cc<end_set[writein_set];cc++){
	      for(int neighbour_index=0;neighbour_index<bond_number_total[bonds_total[current_monomer][0]];neighbour_index++){
		int neighbour=bonds_total[bonds_total[current_monomer][0]][neighbour_index];
		if(neighbour==independent_sets[writein_set][cc]){
		  meet=1;
		}
	      }
	    }
	    if(meet==1){
	      writein_set++;
	      if(writein_set>max_bond_number)
		writein_set=writein_set-max_bond_number-1;
	      if(writein_set==current_set){ 
		writein_set++;
		if(writein_set>max_bond_number)
		  writein_set=writein_set-max_bond_number-1;
	      }
	      break;
	    }
	    else 
	      found =1;
	  }
	  independent_sets[writein_set][end_set[writein_set]]=bonds_total[current_monomer][0];
	  monomer_checked[bonds_total[current_monomer][0]]=-1;
	  total_assigned_monomer++;
	  end_set[writein_set]++;
	 
	}
        
	for(int bb=1;bb<bond_number_total[current_monomer];bb++){
	  if(monomer_checked[bonds_total[current_monomer][bb]]==-1)
	    continue;
	  int found=0;
	  while(found==0){ 
	    int meet=0;
	    for(int neighbour_index_of_bb=0;neighbour_index_of_bb<bond_number_total[bonds_total[current_monomer][bb]];neighbour_index_of_bb++){ //loop over all bonds of bb
	      int neighbour_of_bb=bonds_total[bonds_total[current_monomer][bb]][neighbour_index_of_bb];	  
	      for(int cc=offset_set[writein_set];cc<end_set[writein_set];cc++){
		if(neighbour_of_bb==independent_sets[writein_set][cc]){
		  meet=1;
		}
	      }
	      if(meet==1){
		writein_set++;
		if(writein_set>max_bond_number)
		  writein_set=writein_set-max_bond_number-1;
		if(writein_set==current_set){ 
		  writein_set++;
		  if(writein_set>max_bond_number)
		    writein_set=writein_set-max_bond_number-1;
		}
		break;
	      }
	      else
		found=1;
	    }
	  }
	  independent_sets[writein_set][end_set[writein_set]]=bonds_total[current_monomer][bb];
	  monomer_checked[bonds_total[current_monomer][bb]]=-1;
	  total_assigned_monomer++;
	  end_set[writein_set]++;
	}
	offset_set[current_set]++;      
      }
      current_set++;
      if(current_set>max_bond_number)
	current_set=current_set-max_bond_number-1;
    }   //end while
  }//end loop over all monomer
 printf("set %i %i %i\n",offset_set[1]);
   for(int set_tmp=0;set_tmp<1;set_tmp++){
    for(int pi=0;pi<offset_set[set_tmp];pi++){
      for(int pj=0;pj<offset_set[set_tmp];pj++){
	for(int nei_pi_index=0;nei_pi_index<bond_number_total[independent_sets[set_tmp][pi]];nei_pi_index++){
	  // printf("set %i %i %i\n",bonds_total[independent_sets[set_tmp][pi]][nei_pi_index],independent_sets[set_tmp][pj]);
	  if(independent_sets[set_tmp][pj]==bonds_total[independent_sets[set_tmp][pi]][nei_pi_index]){
	    int this=independent_sets[set_tmp][pi];
	    int that=independent_sets[set_tmp][pj];
	    printf("error %i %i !\n",pi,pj);
	    printf("set %i %i %i\n",bond_number_total[this]);
	  }
	}
      }
    }
    }
  return 0;
}

