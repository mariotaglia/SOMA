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

//! \file independent_sets.c
//! \brief Implementation of independent_sets.h

#include "independent_sets.h"
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "soma_util.h"
#include "phase.h"
//#include "crosslink/readIn_lib.h"
//! Sort the indecies of the set_length in sort_array with the highest at pos=0.
//!
//! \private Used only in the contex of independet set creation
void sort_set_length(const unsigned int n_sets,const unsigned int *const set_length, unsigned int*const sort_array)
    {
    //Bubble sort is ok, because it is not time relevant and the data set is small and the set is highly sorted.
    bool elements_swapped;
    do{
	elements_swapped = false;
	for(unsigned int i=1; i < n_sets;i++)
	    if( set_length[ sort_array[i-1] ] < set_length[ sort_array[i] ] )
		{
		elements_swapped = true;
		unsigned int tmp= sort_array[i];
		sort_array[i] = sort_array[i-1];
		sort_array[i-1] = tmp;
		}
	}while(elements_swapped);
    }

//! Check, whether a particle can be inserted in a set of independet sets.
//!
//! \private  Used only in the contex of independet set creation
bool try_particle_in_set(const unsigned int set_id,const unsigned int p_id,const struct Phase*const p,const unsigned int N,
			 const unsigned int poly_type,const unsigned int*const set_length,const unsigned int * set_matrix)
    {
    bool neighbor_found = false;
    for(unsigned int iTmp=0; iTmp<set_length[set_id]  && neighbor_found == false;iTmp++)
	{
	const unsigned int ibead = set_matrix[set_id*N + iTmp];
	const int start = get_bondlist_offset(
	    p->poly_arch[p->poly_type_offset[poly_type] + ibead + 1]);
	if(start > 0)
	    {
	    int i = start;
	    //BondInfo bn;
	    unsigned int end;
	    do{
		const uint32_t info = p->poly_arch[i++];
		end = get_end(info);
		const int offset = get_offset(info);
		const unsigned int neighbour_id = ibead + offset;
		if(neighbour_id == p_id)
		    neighbor_found = true;
		}while(end==0);
	    }
	}

    return neighbor_found == false;
    }

//! Balance the number of set members across the sets, if possible.
//!
//! \private  Used only in the contex of independet set creation
int balance_sets(const struct Phase *const p,const unsigned int N,const unsigned int poly_type,unsigned int *const set_length,unsigned int *const set_matrix,unsigned int recursion_level)
    {
    const unsigned int max_recursion=2;
    if(recursion_level > max_recursion)
	return 0;
    unsigned int n_sets=0;
    for(unsigned int i=0; i<N; i++)
	if( set_length[i] > 0)
	    n_sets++;
    //Sorted array of set lengths
    unsigned int *sort_array = (unsigned int*)malloc( n_sets*sizeof(unsigned int) );
    if(sort_array == NULL)
	{
	fprintf(stderr,"ERROR: malloc %s:%d\n",__FILE__,__LINE__);
	return -1;
	}
    for(unsigned int i=0; i < n_sets;i++)
	sort_array[i] = i;

    bool solvable_unbalanced;
    do{
	sort_set_length(n_sets, set_length, sort_array);

	solvable_unbalanced = false; //Set to true, if a element can be moved in the loop.
	//Try to an element of the largest set to any other set.
	for(unsigned int mono=0; mono < set_length[sort_array[0]]; mono++)
	    {
	    unsigned move_id = set_matrix[ sort_array[0] * N + mono];
	    //Try to fit in any of the other sets, start with the smallest.
	    unsigned int new_set = 0;
	    for( new_set = n_sets-1; new_set > 0;new_set--)
		//New set has to be smaller than the largest set, even if we remove one element of the largest set.
		if( set_length[ sort_array[new_set] ] < set_length[sort_array[0]]-1)
		    {
		    const bool suitable_set = try_particle_in_set( sort_array[new_set] , move_id, p, N, poly_type, set_length, set_matrix);
		    if( suitable_set )
			break;
		    }

	    if(new_set > 0)//Found an element to move and ****BREAK**** for loop, but do not stop while loop.
		{
		//reorganise the set_matrix of the largest set.
		for(unsigned int i=mono; i < set_length[ sort_array[0] ]-1; i++)
		    set_matrix[ sort_array[0]*N + i] = set_matrix[ sort_array[0]*N + i+1];
		set_length[ sort_array[0] ] -= 1;

		//Insert the element in its new set.
		set_matrix[ sort_array[new_set]*N + set_length[ sort_array[new_set] ] ] = move_id;
		set_length[ sort_array[new_set] ] += 1;

		solvable_unbalanced = true; //
		break; //Exit for loop, because we want to move one element at a time.
		}
	    }
	}while(solvable_unbalanced);

    //If we have even after rebalancing with a constant number of sets
    //a bad divergence, we can add a set and try the rebalancing
    //again.
    const unsigned int max_tolerated_divergence = 16;
    sort_set_length(n_sets, set_length, sort_array);
    unsigned int divergence = 0;
    for(unsigned int set=0; set < n_sets; set++)
	if( set_length[ sort_array[0] ] - set_length[set] > divergence)
	    divergence = set_length[ sort_array[0] ] - set_length[set];
    if( divergence > max_tolerated_divergence && n_sets < N)
	{
	set_length[ sort_array[0] ] -= 1;
	set_matrix[ n_sets*N + 0] = set_matrix[ sort_array[0]*N + set_length[ sort_array[0]] ];
	set_length[ n_sets ] += 1;

	balance_sets(p, N, poly_type, set_length, set_matrix, recursion_level+1);
	}

    free(sort_array);
    return 0;
    }

int generate_independet_sets(struct Phase*const p)
{
  int res=0; 
  
  switch( p->args.set_generation_algorithm_arg)
    {
    case set_generation_algorithm_arg_FIXEDMINUS_NMINUS_SETS:
      res=independent_set_fixed(p);
      break;
    case set_generation_algorithm_arg_SIMPLE:  
      res=independent_sets_simple(p);
      break;
    case set_generation_algorithm__NULL:
    default:
      fprintf(stderr,"ERROR: Unknown independent set generation algorithm selected.\n");
    }
  return 0;
}


int independent_sets_simple(struct Phase* const p)
    {
    struct IndependetSets*const sets = (struct IndependetSets*)malloc( p->n_poly_type * sizeof(IndependetSets) );
    p->max_set_members = 0;
    if(sets == NULL)
	{
	fprintf(stderr,"ERROR: malloc %s:%d\n",__FILE__,__LINE__);
	return -1;
	}
    unsigned int max_nsets = 0;
    for(unsigned int poly_type = 0; poly_type < p->n_poly_type; poly_type++)
	{
	const unsigned int N = p->poly_arch[ p->poly_type_offset[poly_type] ];
	unsigned int*const set_length = (unsigned int*)malloc( N*sizeof(unsigned int));
	if(set_length == NULL)
	    {
	    fprintf(stderr,"ERROR: malloc %s:%d\n",__FILE__,__LINE__);
	    return -2;
	    }
	memset(set_length,0,N*sizeof(unsigned int));
	unsigned int*const set_matrix = (unsigned int*)malloc( N*N*sizeof(unsigned int));
	if(set_matrix == NULL)
	    {
	    fprintf(stderr,"ERROR: malloc %s:%d\n",__FILE__,__LINE__);
	    return -3;
	    }

	//Start setting the monomers in the sets.
	for(unsigned int mono = 0; mono < N ;mono++)
	    {
	    // Find the next possible set to sort in
	    unsigned int set_to_sort_in = 0;
	    bool set_found = false;
	    //Try all sets
	    while(!set_found)
		{
		//Check for neighbors in the set
		const bool suitable_set = try_particle_in_set(set_to_sort_in, mono, p, N, poly_type, set_length, set_matrix);

		if(suitable_set)
		    set_found = true;
		set_to_sort_in += 1;
		assert(set_to_sort_in < N+1);
		}
	    set_to_sort_in -= 1;

	    //Sort the monomer in the set
	    set_matrix[ set_to_sort_in*N + set_length[set_to_sort_in]] = mono;
	    set_length[set_to_sort_in] += 1;
	    }

	balance_sets(p,N,poly_type,set_length,set_matrix,0);

	//Store the found sets to the final structure and free intermediate arrays.
	unsigned int n_sets=0;
	unsigned int max_set=0;
	unsigned int member_sum=0;
	for(unsigned int i=0; i<N;i++)
	    {
	    if(set_length[i] > 0)
		n_sets++;
	    if(set_length[i] > max_set)
		max_set = set_length[i];
	    member_sum += set_length[i];

	    }
	assert(member_sum == N);

	sets[poly_type].n_sets= n_sets;
	sets[poly_type].max_member = max_set;
	sets[poly_type].set_length = (unsigned int*)malloc( n_sets * sizeof(unsigned int));
	if(sets[poly_type].set_length == NULL)
	    {
	    fprintf(stderr,"ERROR: malloc %s:%d\n",__FILE__,__LINE__);
	    return -4;
	    }
	sets[poly_type].sets = (unsigned int*)malloc( n_sets*max_set * sizeof(unsigned int));
	if(sets[poly_type].sets == NULL)
	    {
	    fprintf(stderr,"ERROR: malloc %s:%d\n",__FILE__,__LINE__);
	    return -5;
	    }
	//Copy data to the new structure
	for(unsigned int iSet=0; iSet < n_sets; iSet++)
	    {
	    sets[poly_type].set_length[iSet] = set_length[iSet];
	    if(set_length[iSet] > p->max_set_members)
		p->max_set_members =set_length[iSet];
	    for(unsigned int iMember=0;iMember<set_length[iSet];iMember++)
		sets[poly_type].sets[ iSet*max_set + iMember] = set_matrix[iSet*N + iMember];
	    }
	//Free intermediate arrays
	free(set_length);
	free(set_matrix);
	if(sets[poly_type].n_sets > max_nsets)
	  max_nsets = sets[poly_type].n_sets;
	}
    p->sets = sets;
    p->max_n_sets = max_nsets;

    //Allocate and init memory for polymer states
    for(unsigned int i=0; i < p->n_polymers;i++)
	{
	Polymer*const poly = p->polymers+i;

	poly->set_permutation = (unsigned int*)malloc( p->max_n_sets * sizeof(unsigned int));
	if( poly->set_permutation == NULL )
	    {
	    fprintf(stderr,"ERROR: Malloc %s:%d\n",__FILE__,__LINE__);
	    return -1;
	    }

	poly->set_states = (struct RNG_STATE*)malloc(p->max_set_members * sizeof(struct RNG_STATE));
	if( poly->set_states == NULL )
	    {
	    fprintf(stderr,"ERROR: Malloc %s:%d\n",__FILE__,__LINE__);
	    return -1;
	    }

	//Init every state in the polymer
	const unsigned int seed = soma_rng_uint(&(poly->poly_state),pseudo_random_number_generator_arg_PCG32);
	for(unsigned int j=0; j < p->max_set_members; j++)
	    {
	    struct RNG_STATE*const state = &(poly->set_states[j]);
	    allocate_rng_state(state, p->args.pseudo_random_number_generator_arg);
	    seed_rng_state(state, seed, j, p->args.pseudo_random_number_generator_arg);
	    }
	}

    return 0;
    }



int independent_set_fixed(struct Phase* const poly){
  
  struct IndependetSets*const set_tmp = (struct IndependetSets*)malloc(1* sizeof(IndependetSets) ); 
  unsigned int sequence=poly->num_all_beads;
  unsigned int max_bond_number=0,max_bond=0;

  int* bond_number_total;
  bond_number_total= (int*) malloc(sequence*sizeof(int));
  if( bond_number_total== NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
    return -4;
  }
  memset(bond_number_total,0,sequence*sizeof(int));
  unsigned int **bonds_total=(unsigned int **)malloc(sequence*sizeof(unsigned int*));
  if(bonds_total == NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
    return -4;
  }
  
  for(int current_tmp=poly->poly_type_offset[1]+1;current_tmp<sequence+poly->poly_type_offset[1]+1;current_tmp++){
  
    uint32_t bonds_of_monomer=0, bond_number=0;
    uint32_t  current_poly_arch=poly->poly_arch[current_tmp];
    int start_offset_bond=get_bondlist_offset(current_poly_arch);
    do{
      bonds_of_monomer=poly->poly_arch[start_offset_bond];
      int end=get_end(bonds_of_monomer);		
      bond_number++;
      if(end==1) break;
      start_offset_bond++;
    }while(0==0); 

    bond_number_total[current_tmp-poly->poly_type_offset[1]-1]=bond_number;
    bonds_total[current_tmp-poly->poly_type_offset[1]-1]=(unsigned int *)malloc(bond_number*sizeof(unsigned int));
    memset(&bonds_total[current_tmp-poly->poly_type_offset[1]-1][0],0,bond_number*sizeof(unsigned int));     
	  
    start_offset_bond=get_bondlist_offset(poly->poly_arch[current_tmp]);
    for(int tmp=0;tmp<bond_number;tmp++){      
      bonds_of_monomer=poly->poly_arch[start_offset_bond];     
      bonds_total[current_tmp-poly->poly_type_offset[1]-1][tmp]=get_offset(bonds_of_monomer)+current_tmp-poly->poly_type_offset[1]-1;
      start_offset_bond++;
    }
    
    if(bond_number>max_bond_number){    
      max_bond_number=bond_number; 
      max_bond=current_tmp-poly->poly_type_offset[1]-1;
    }  
 
  }  
  
  int* monomer_checked;
  monomer_checked= (int*) malloc(sequence*sizeof(int));
  if(monomer_checked== NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
    return -4;
  }
  memset(monomer_checked,0,sequence*sizeof(int));
  unsigned int current_monomer=max_bond;
  unsigned int** independent_sets;
  independent_sets= (unsigned int**) malloc((max_bond_number+1)*sizeof(unsigned int*));
  if( independent_sets== NULL){
    fprintf(stderr,"Malloc problem %s:%d\n",__FILE__,__LINE__);
    return -4;
  }
  //unsigned int offset_set[max_bond_number+1]; //only members between offset_set and end_set are unchecked for neighbours  
  //unsigned int end_set[max_bond_number+1];

  unsigned int* end_set= (unsigned int*) malloc((max_bond_number+1)*sizeof(unsigned int));
  unsigned int* offset_set= (unsigned int*) malloc((max_bond_number+1)*sizeof(unsigned int));
  for(int aaa=0;aaa<max_bond_number+1;aaa++){
    end_set[aaa]=0;
    offset_set[aaa]=0;
  }
  for(int bb=0;bb<max_bond_number+1;bb++){
    independent_sets[bb]=(unsigned int*) malloc(sequence*sizeof(unsigned int));
    memset(&independent_sets[bb][0],0,sequence*sizeof(unsigned int)); 
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
      end_set[aa]++;
    }
    offset_set[0]=offset_set[0]+1;//neighbour of inde*[0][offset_set] are all found now 
    int current_set=1; //the one needs to be checked
 
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
	      unsigned int neighbour_of_bb=bonds_total[bonds_total[current_monomer][bb]][neighbour_index_of_bb];	  
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
  set_tmp[0].n_sets=max_bond_number+1;  
  int ccc=0;
  for(int aaa=0;aaa<max_bond_number+1;aaa++){
    if(offset_set[aaa]>set_tmp[0].max_member)
      set_tmp[0].max_member=offset_set[aaa];
  }

  unsigned int * inde_set_tmp= (unsigned int*) malloc(set_tmp[0].max_member*(max_bond_number+1)*sizeof(unsigned int));
  for(int aaa=0;aaa<max_bond_number+1;aaa++){
    ccc=aaa*set_tmp[0].max_member;
    for(int bbb=0;bbb<offset_set[aaa];bbb++){
      inde_set_tmp[ccc]=independent_sets[aaa][bbb];
      ccc++;
    }
  }
  set_tmp[0].set_length = (unsigned int*)malloc((max_bond_number+1) * sizeof(unsigned int));
  set_tmp[0].set_length=offset_set;
  /* for(int aaa=0;aaa<max_bond_number+1;aaa++){
    set_tmp[0].set_length[aaa]=offset_set[aaa]; 
    }*/
  set_tmp[0].sets=inde_set_tmp;
  poly->sets=set_tmp;

  //Allocate and init memory for polymer states
  for(unsigned int i=0; i < poly->n_polymers;i++)
    {
      Polymer*const poly_tmp = poly->polymers+i;

      poly_tmp->set_permutation = (unsigned int*)malloc( poly->max_n_sets * sizeof(unsigned int));
      if( poly_tmp->set_permutation == NULL )
	{
	  fprintf(stderr,"ERROR: Malloc %s:%d\n",__FILE__,__LINE__);
	  return -1;
	}

      poly_tmp->set_states = (struct RNG_STATE*)malloc(poly->max_set_members * sizeof(struct RNG_STATE));
      if( poly_tmp->set_states == NULL )
	{
	  fprintf(stderr,"ERROR: Malloc %s:%d\n",__FILE__,__LINE__);
	  return -1;
	}

      //Init every state in the polymer
      const unsigned int seed = soma_rng_uint(&(poly_tmp->poly_state),pseudo_random_number_generator_arg_PCG32);
      for(unsigned int j=0; j < poly->max_set_members; j++)
	{
	  struct RNG_STATE*const state = &(poly_tmp->set_states[j]);
	  allocate_rng_state(state, poly->args.pseudo_random_number_generator_arg);
	  seed_rng_state(state, seed, j, poly->args.pseudo_random_number_generator_arg);
	}
    }
  return 0;
}

