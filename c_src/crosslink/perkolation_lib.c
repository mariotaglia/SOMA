#include <hdf5.h>
#include "readIn_lib.h"
#include "perkolation_lib.h"
#include "io.h"
#include <math.h>
#include <stdlib.h>

void perkolation(const char* const filename,const char* const perkofile)
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
   
  int perko=check_perkolation_monomer1(filename,b,poly,0);
  printf("number of perkolation: %i.\n",perko);

  FILE *fp;

  fp = fopen(perkofile, "a+");
  if (fp == NULL)
    {
      printf("Error opening file!\n");
    }
  fprintf(fp, "%i\n",perko);

  fclose(fp);
  
  // free_data(b,poly);
}


int check_perkolation_monomer1(const char* const filename,bead_data*b,polymer_data*poly, const int poly_start){
  int wang=0;
  herr_t status;
  hid_t  plist_id=H5Pcreate(H5P_FILE_ACCESS);
  hid_t  file_id=H5Fopen(filename, H5F_ACC_RDONLY,plist_id);
  H5Pclose(plist_id);
  plist_id=H5Pcreate(H5P_DATASET_XFER);
  
  uint32_t*const poly_arch_new = (uint32_t*) malloc( poly->poly_arch_length * sizeof(uint32_t));
  if(poly_arch_new == NULL){fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
    return -1;
  }
  status = read_hdf5(file_id,"/parameter/poly_arch",H5T_NATIVE_INT32,plist_id,poly_arch_new);
  HDF5_ERROR_CHECK2(status,"poly_arch");
  
  static int number_loop=0;
  int bond_max=0;
  int number_of_beads=(int16_t) b->number_of_beads[0];
  unsigned int sequence=0;
  
  for(unsigned int i=0;i<b->old_N_polymers;i++){
    sequence=sequence+b->number_of_beads[i];
  }

  ////new 29.11.2017
  for(int current_tmp=poly->poly_type_offset[1]+1;current_tmp<sequence+poly->poly_type_offset[1]+1;current_tmp++){
    int bonds_of_monomer=0, bond_number=0;
    int32_t  current_poly_arch=poly->poly_arch[current_tmp];
    int start_offset_bond=get_bondlist_offset(current_poly_arch);
    do{
      bonds_of_monomer=poly->poly_arch[start_offset_bond];
      int end=get_end(bonds_of_monomer);		
      bond_number++;
      if(end==1) break;
      start_offset_bond++;
    }while(0==0); 
    if(bond_number>bond_max)
      bond_max=bond_number;
  }
  ////new 29.11.2017


  unsigned int esti=sequence;
  //start: book memory
  int *checked_monomer;
  checked_monomer= (int*) malloc(sequence*sizeof(int));
  memset(checked_monomer,-1,sequence*sizeof(int));

  int *bond_nodes_checked;
  bond_nodes_checked=(int*) malloc(esti*sizeof(int));
  if(bond_nodes_checked == NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
    return -4;
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

  int *bond_total;
  bond_total= (int*) malloc(bond_max*sizeof(int));

  int *loop_member;
  loop_member= (int*) malloc(sequence*sizeof(int));
    
  if(loop_member == NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
    return -4;
  }

  int **checked_bonds_node=(int **)malloc(esti*sizeof(int*));
  for(unsigned int aa=0;aa<esti;aa++){
    checked_bonds_node[aa]=(int *)malloc(bond_max*sizeof(int));  
  }   
  
  if(checked_bonds_node == NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
    return -4;
  }
  //end: book memory
  int current_monomer=poly->poly_type_offset[1]+1+poly_start*number_of_beads;
   
  
  for(unsigned int mono_i=poly_start;mono_i<sequence;mono_i++){ //start: loop over all monomers
    if(checked_monomer[mono_i]==0)
      continue;
  
    current_monomer=poly->poly_type_offset[1]+1+mono_i;	
    //start: if this chain is not crosslinked, then skip it
    int number_of_ends=0;
    int current_monomer_tmp=current_monomer;
    int bond_number_tmp=0,number_of_uncrosslinked_monomer=0;   
    int flag_crosslink=0;
  
    while(number_of_ends<2){
      current_monomer_tmp++;
      
      int32_t current_poly_arch_tmp=poly_arch_new[current_monomer_tmp];
      int start_offset_bond_tmp=get_bondlist_offset(current_poly_arch_tmp);
      number_of_uncrosslinked_monomer++;
      int number_of_neigh=0;
      while(1==1)
	{
	  int end=get_end(poly_arch_new[start_offset_bond_tmp]);	
	  if(get_offset(poly_arch_new[start_offset_bond_tmp])!=1&&get_offset(poly_arch_new[start_offset_bond_tmp])!=-1){
	    flag_crosslink=1;
	  }
	  if(get_offset(poly_arch_new[start_offset_bond_tmp])==1||get_offset(poly_arch_new[start_offset_bond_tmp])==-1)
	    number_of_neigh++;
	  bond_number_tmp++;
	  if(end==1) break;
	  start_offset_bond_tmp++;
	}
      if(number_of_neigh==1)
	number_of_ends++;
    }
    if(flag_crosslink==0){
      mono_i=mono_i+number_of_uncrosslinked_monomer-1;
      continue;
    }
    //end: if this chain is not crosslinked, then skip it
    double Lx=poly->boxsize[0],Ly=poly->boxsize[1];
    int normal_node=0,finish=0;
    int loop_member_index=0;
    int checked_bonds[2]; //saves the path from last monomer to current monomer
    checked_bonds[0]=0;
    checked_bonds[1]=0;
        
    int number_of_nodes=0;
    //start:set all elements of the memory to zero for a new chain
    memset(bond_total,0,bond_max*sizeof(int));
    memset(loop_member,0,sequence*sizeof(int));
    memset(nodes_list,0,esti*sizeof(int));
    memset(nodes_offset,0,esti*sizeof(int));
    memset(bond_number_total,0,esti*sizeof(int));
    memset(bond_nodes_checked,0,esti*sizeof(int));
    for(unsigned int aa=0;aa<esti;aa++){     
      memset(&checked_bonds_node[aa][0],0,bond_max*sizeof(int));     
    }
    //end: set all elements of the memory to zero for a new chain
    int current_poly_arch=poly_arch_new[current_monomer];

    while(finish==0){    //start: if the chain does not end, than stay in this loop   
      checked_monomer[current_monomer-poly->poly_type_offset[1]-1]=0;
      for(int i=0;i<bond_max;i++){
	bond_total[i]=0;
      }
      int side_node=0;
      int exist=0;
      int loop=0;
       
      for(int loop_index_tmp=0;loop_index_tmp<loop_member_index&&exist==0;loop_index_tmp++){//start: loop over all the former elements of the current chain to look for loops
	if(loop_member[loop_index_tmp]==current_monomer){//start: a loop is present, look for percolation
	 
	  loop=1;
	  double sumx=0,sumy=0,sumz=0;
	  exist=1;
	  side_node=1;
	  double *loop_memberx;
	  loop_memberx= (double*) malloc((loop_member_index-loop_index_tmp+1)*sizeof(double));
	  double *loop_membery;
	  loop_membery= (double*) malloc((loop_member_index-loop_index_tmp+1)*sizeof(double));

	  if(loop_memberx == NULL||loop_membery == NULL){
	    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
	    return -4;
	  }
	  memset(loop_memberx,0,(loop_member_index-loop_index_tmp+1)*sizeof(int));
	  memset(loop_membery,0,(loop_member_index-loop_index_tmp+1)*sizeof(int));
	  int jjj=0;
	  for(int index_in_the_loop=loop_index_tmp;index_in_the_loop<loop_member_index-1;index_in_the_loop++){
	    double tmpx=1,tmpy=1;
	    loop_memberx[jjj]=-b->beads[loop_member[index_in_the_loop]-poly->poly_type_offset[1]-1].x+b->beads[loop_member[index_in_the_loop+1]-poly->poly_type_offset[1]-1].x;
	    loop_membery[jjj]=-b->beads[loop_member[index_in_the_loop]-poly->poly_type_offset[1]-1].y+b->beads[loop_member[index_in_the_loop+1]-poly->poly_type_offset[1]-1].y;
	    if(loop_memberx[jjj]<0)
	      tmpx=-1;
	    if(loop_membery[jjj]<0)
	      tmpy=-1;
	    sumx=sumx+loop_memberx[jjj]-(int)(loop_memberx[jjj]/Lx+tmpx/2)*Lx;
	    sumy=sumy+loop_membery[jjj]-(int)(loop_membery[jjj]/Ly+tmpy/2)*Ly;	   
	    jjj++;	   
	  }
	  double tmpx=1,tmpy=1;
	  double distancex=b->beads[current_monomer-poly->poly_type_offset[1]-1].x-b->beads[loop_member[loop_member_index-1]-poly->poly_type_offset[1]-1].x;
	  double distancey=b->beads[current_monomer-poly->poly_type_offset[1]-1].y-b->beads[loop_member[loop_member_index-1]-poly->poly_type_offset[1]-1].y;  
	  if(distancex<0)
	    tmpx=-1;
	  if(distancey<0)
	    tmpy=-1;
	  sumx=sumx+distancex-(int)(distancex/Lx+tmpx/2)*Lx;
	  sumy=sumy+distancey-(int)(distancey/Ly+tmpy/2)*Ly;
	  
	 
	  if(sumx>=poly->boxsize[0]||sumy>=poly->boxsize[1]||sumx<=-poly->boxsize[0]||sumy<=-poly->boxsize[1]){
	    number_loop++;
	    /* for(int index_in_the_loop=loop_index_tmp;index_in_the_loop<loop_member_index;index_in_the_loop++){
	        printf("member %i,%f, %f\n",loop_member[index_in_the_loop]-poly->poly_type_offset[1]-1,b->beads[loop_member[index_in_the_loop]-poly->poly_type_offset[1]-1].x,b->beads[loop_member[index_in_the_loop]-poly->poly_type_offset[1]-1].y);
	    }
	    printf("member %i,%f, %f, %f, %f\n",current_monomer-poly->poly_type_offset[1]-1,b->beads[current_monomer-poly->poly_type_offset[1]-1].x,b->beads[current_monomer-poly->poly_type_offset[1]-1].y,sumx);*/
	     if(sumx>=poly->boxsize[0]||sumx<=-poly->boxsize[0]){
	      printf("found a perkolation in x direction \n");		     
	    }   
	    if(sumy>=poly->boxsize[1]||sumy<=-poly->boxsize[1]){
	      printf("found a perkolation in y direction \n");
	      }	    	    
	    break;
	  }
	  free(loop_memberx);
	  free(loop_membery);
	  if(current_monomer==nodes_list[number_of_nodes-1]){	 
	    number_of_nodes--;	    
	  }	 
	}//end: a loop is present, look for percolation
      }//end: loop over all the former elements of the current chain to look for loops
      loop_member[loop_member_index]=current_monomer;
      loop_member_index++;      
      current_poly_arch=poly_arch_new[current_monomer];
      if(current_poly_arch>4294967200){
	fprintf(stderr, "ERROR: Function get_bondlist_offset exited the range %s:%d\n", __FILE__, __LINE__);
	return -2;
      }     
      //printf("current %i,%i,%i\n",current_monomer,wang,sequence);
      int start_offset_bond=get_bondlist_offset(current_poly_arch);
      int bonds_of_monomer;	    
      int bond_number=0;
      int if_end_chain=0;     
      do{
	bonds_of_monomer=poly_arch_new[start_offset_bond];
	int end=get_end(bonds_of_monomer);
	bond_total[bond_number]=get_offset(bonds_of_monomer);		
	bond_number++;
	if(end==1) break;
	start_offset_bond++;
      }while(0==0);
       
      if(bond_number==1&&checked_bonds[0]!=0){	//start: this is the end of a chain, jump back to the last node or continue with the next chain (if all nodes are checked)
	int ii=1;       
	while(number_of_nodes-ii>=0){
	  if(bond_number_total[number_of_nodes-ii]!=bond_nodes_checked[number_of_nodes-ii])
	    break;
	  else{	    
	    ii++;
	  }
	}	
	if(number_of_nodes-ii<0){//find the last node with unchecked bond
	  finish=1;
	  continue;
	}
	loop_member_index=nodes_offset[number_of_nodes-ii]-1;
	current_monomer=nodes_list[number_of_nodes-ii];
	checked_bonds[0]=0;

	number_of_nodes=number_of_nodes-ii;

	int tmp=fmin(ii+5,(esti-number_of_nodes-1));
	memset(&bond_nodes_checked[number_of_nodes+1],0,tmp*sizeof(int));
	memset(&bond_number_total[number_of_nodes+1],0,tmp*sizeof(int));
	for(unsigned int node_tmp=number_of_nodes+1;node_tmp<number_of_nodes+1+tmp;node_tmp++){
	  memset(&checked_bonds_node[node_tmp][0],0,bond_max*sizeof(int));	 
	  } 
	normal_node=1;
	continue;
      }	//end: this is the end of a chain, jump back to the last node or continue with the next chain (if all nodes are checked)
     
      if((bond_number==2&&checked_bonds[0]!=0&&side_node!=1)||(bond_number==1&&checked_bonds[0]==0)){//start: this is either the first element (with only one neighbour) or a middle one (without other crosslinks)
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
      }   //end: this is either the first element (with only one neighbour) or a middle one (without other crosslinks)
      if(bond_number>=2&&checked_bonds[0]==0){ //decide if the first element of the chain but crosslinked to other chains
	side_node=1;
	normal_node=1;
      }
           
      if(bond_number>2||side_node==1){ //start: this is a node, need to check for more than one bond
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

	if(found_the_node==1){ //if this node is saved before, the previous one needs also to be updated
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
      }      //end: this is a node, need to check for more than one bond     	  
      number_of_nodes++;
      bond_number_total[number_of_nodes-1]=bond_number;

      if(if_end_chain==1){ //jump to the last node with unchecked bonds
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
	loop_member_index=nodes_offset[number_of_nodes-ii]-1;
	current_monomer=nodes_list[number_of_nodes-ii];
	checked_bonds[0]=0;
	number_of_nodes=number_of_nodes-ii;
	int tmpp=fmin(ii+5,(esti-number_of_nodes-1));
	memset(&bond_number_total[number_of_nodes+1],0,tmpp*sizeof(int));
	memset(&bond_nodes_checked[number_of_nodes+1],0,tmpp*sizeof(int));
	for(unsigned int node_tmp=number_of_nodes+1;node_tmp<number_of_nodes+1+tmpp;node_tmp++){
	  memset(&checked_bonds_node[node_tmp][0],0,bond_max*sizeof(int));
	}           
	continue;
      }	         
    }// endt: if the chain does not end, than stay in this loop   
  }//end: loop over all monomers
  for(unsigned int aa=0;aa<esti;aa++){
    free(checked_bonds_node[aa]);     
  }
  free(checked_bonds_node);
  free(bond_nodes_checked);
  free(loop_member);
  free(bond_number_total);
  free(nodes_offset);
  free(nodes_list);
  free(bond_total);

  free(poly_arch_new);
  free(checked_monomer);
  return number_loop;
}

/*int check_configfile(const char* const filename)
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
   
  Monomer *new_position=feed_back_box(b,poly);
  printf("finish feed back! \n");
  int tmp=0;
  //printf(" %f\n",new_position[1].x);
  // for(int Poly_i=0;Poly_i<b->N_polymers;Poly_i++){
  int Poly_i=19185;
  double sumx=0,sumy=0,sumz=0;
  for(int mono_i=0;mono_i<31;mono_i++){
    sumx=sumx+b->beads[Poly_i*32+mono_i+1].x-new_position[Poly_i*32+mono_i].x;
    sumy=sumy+b->beads[Poly_i*32+mono_i+1].y-new_position[Poly_i*32+mono_i].y;
    sumz=sumz+b->beads[Poly_i*32+mono_i+1].z-new_position[Poly_i*32+mono_i].z;
    printf(" %f\n",b->beads[Poly_i*32+mono_i].y);
  }
  printf("tmp %i\n",tmp);
  return 0;

  }*/
