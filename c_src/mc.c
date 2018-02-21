/* Copyright (C) 2016-2017 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016-2017 Marcel Langenberg
   Copyright (C) 2016 Fabien Leonforte
   Copyright (C) 2016 Juan Orozco
   Copyright (C) 2016 Yongzhi Ren

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
#include "mc.h"
#include <math.h>
#include <stdbool.h>
#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include "phase.h"
#include "mesh.h"
#include "independent_sets.h"

void trial_move(const Phase * p, const uint64_t ipoly, const int ibead,
		soma_scalar_t *dx, soma_scalar_t *dy, soma_scalar_t *dz,const unsigned int iwtype,const enum enum_pseudo_random_number_generator arg_rng_type,RNG_STATE*const rng_state)
{
  //Just to shut up the compiler warning:
  //Any decent compiler optimize it out
  soma_scalar_t scale=ibead + 0*ipoly;
  scale = p->A[iwtype];

  *dx = scale * (soma_rng_soma_scalar(rng_state,arg_rng_type) - 0.5);
  *dy = scale * (soma_rng_soma_scalar(rng_state,arg_rng_type) - 0.5);
  *dz = scale * (soma_rng_soma_scalar(rng_state,arg_rng_type) - 0.5);
}

void trial_move_cm(const Phase * p, const uint64_t poly_type,soma_scalar_t *const dx, soma_scalar_t *const dy, soma_scalar_t *const dz,
		   const enum enum_pseudo_random_number_generator arg_rng_type,RNG_STATE*const rng_state)
{
#ifndef _OPENACC
  assert( p->cm_a );
#endif//_OPENACC
  const soma_scalar_t scale= p->cm_a[poly_type];

  *dx = scale * (soma_rng_soma_scalar(rng_state,arg_rng_type) - 0.5);
  *dy = scale * (soma_rng_soma_scalar(rng_state,arg_rng_type) - 0.5);
  *dz = scale * (soma_rng_soma_scalar(rng_state,arg_rng_type) - 0.5);
}


bool som_accept(RNG_STATE *const rng,  enum enum_pseudo_random_number_generator rng_type , soma_scalar_t delta_energy)
{
  // \todo kBT reqired
  const soma_scalar_t p_acc = exp(-1.0 * delta_energy );

  //Use lazy eval.
  if ((p_acc > 1) || (p_acc > soma_rng_soma_scalar(rng, rng_type)))
    {
      return true;
    }
  else
    {
      return false;
    }
}

soma_scalar_t calc_delta_nonbonded_energy(const Phase * p,const Monomer*const monomer,
					  const soma_scalar_t dx, const soma_scalar_t dy,const soma_scalar_t dz,
					  const unsigned int iwtype)
{
  // Old non-bonded interaction
  const soma_scalar_t xold = monomer->x;
  const soma_scalar_t yold = monomer->y;
  const soma_scalar_t zold = monomer->z;
  const uint64_t cellindex_old = coord_to_index_unified(p, xold, yold, zold, iwtype);
  const soma_scalar_t energy_old = p->omega_field_unified[cellindex_old];

  const uint64_t cellindex_new = coord_to_index_unified(p, xold+dx, yold+dy, zold+dz, iwtype);

  // New non-bonded interaction
  const soma_scalar_t energy_new = p->omega_field_unified[cellindex_new];
  const soma_scalar_t energy = energy_new - energy_old;
  return energy;
}

soma_scalar_t calc_delta_energy(const Phase * p, const uint64_t ipoly,const Monomer*const monomer,
				const unsigned int ibead,const soma_scalar_t dx,const soma_scalar_t dy,
				const soma_scalar_t dz,const unsigned int iwtype)
{
  const  soma_scalar_t delta_nonbonded_energy = calc_delta_nonbonded_energy(p,monomer,dx,dy,dz,iwtype);
  const soma_scalar_t delta_bonded_energy = calc_delta_bonded_energy(p, monomer,ipoly, ibead, dx, dy, dz);

  // non-bonded energy + bonded energy
  soma_scalar_t energy = delta_nonbonded_energy;
  energy += delta_bonded_energy;
  return energy;
}

soma_scalar_t calc_delta_bonded_energy(const Phase * const p,const Monomer*const monomer,
				       const uint64_t ipoly,const unsigned int ibead,
				       const soma_scalar_t dx,const soma_scalar_t dy, const soma_scalar_t dz)
{
  soma_scalar_t delta_energy = 0;
  // loop over bonds of this bead
  //printf("== %d \n",ibead) ;
  //printf("   poly type = %d \n",p->polymers[ipoly].type) ;
  const int start = get_bondlist_offset(p->poly_arch[p->poly_type_offset[p->polymers[ipoly].type] + ibead + 1]);

  //printf("   start = %d \n",start) ;

  if(start > 0){
    int i = start;
    //BondInfo bn;
    unsigned int end;
    do{
      const uint32_t info = p->poly_arch[i++];
      end = get_end(info);
      const unsigned int bond_type = get_bond_type(info);
      const int offset = get_offset(info);

      const int neighbour_id = ibead + offset;
      const unsigned int jbead = neighbour_id;
      //printf("    offset=%d jbead=%u  end=%u type=%u\n",neigh->offset,jbead,end,neigh->bond_type);

      soma_scalar_t scale = 1.;
      switch (bond_type) {
      case HARMONICVARIABLESCALE:
	scale = p->harmonic_normb_variable_scale;
	/* intentionally falls through */
      case HARMONIC:
	//Empty statement, because a statement after a label
	//has to come before any declaration
	;
	const soma_scalar_t old_rx = calc_bond_length(monomer->x,p->polymers[ipoly].beads[jbead].x,p->Lx,p->args.bond_minimum_image_convention_flag);
	const soma_scalar_t new_rx = old_rx + dx;
	const soma_scalar_t old_ry = calc_bond_length(monomer->y,p->polymers[ipoly].beads[jbead].y,p->Ly,p->args.bond_minimum_image_convention_flag);
	const soma_scalar_t new_ry = old_ry + dy;
	const soma_scalar_t old_rz = calc_bond_length(monomer->z,p->polymers[ipoly].beads[jbead].z,p->Lz,p->args.bond_minimum_image_convention_flag);
	const soma_scalar_t new_rz = old_rz + dz;

	const soma_scalar_t old_r2 =
	  old_rx * old_rx + old_ry * old_ry + old_rz * old_rz;
	const soma_scalar_t new_r2 =
	  new_rx * new_rx + new_ry * new_ry + new_rz * new_rz;
	delta_energy += p->harmonic_normb * (new_r2 - old_r2) *scale;
	break;

      case STIFF:
#ifndef _OPENACC
	fprintf(stderr,
		"ERROR: %s:%d stiff bond not yet implemented.\n",
		__FILE__, __LINE__);
#endif//_OPENACC
	break;



      default:
#ifndef _OPENACC
	fprintf(stderr, "ERROR: %s:%d unknow bond type appeared %d\n",
		__FILE__, __LINE__,bond_type);
#endif//OPENACC
	break;
      }

    }while( end == 0);
  }
  return delta_energy;
}

int monte_carlo_propagation(Phase*const p,unsigned int nsteps)
{
  //Update the omega fields for the calculations.

  update_omega_fields(p);

  int ret;
  start_autotuner(&(p->mc_autotuner));
  switch( p->args.iteration_alg_arg)
    {
    case iteration_alg_arg_POLYMER:
      ret = mc_polymer_iteration(p,nsteps,p->mc_autotuner.value);
      break;
    case iteration_alg_arg_SET:
      ret = mc_set_iteration(p,nsteps,p->mc_autotuner.value);
      break;
    case iteration_alg__NULL:
    default:
      fprintf(stderr,"ERROR: Unknown iteration algorithm selected.\n");
      ret = 1;
    }
  end_autotuner(&(p->mc_autotuner));

  if( p-> cm_a )
    {
      start_autotuner(&(p->cm_mc_autotuner));
      mc_center_mass(p, 1, p->cm_mc_autotuner.value);
      end_autotuner(&(p->cm_mc_autotuner));
    }
  if( p->args.autotuner_restart_period_arg > 0 && p->time % p->args.autotuner_restart_period_arg == 0)
    {
      restart_autotuner( &(p->mc_autotuner) );
      restart_autotuner( &(p->cm_mc_autotuner) );
    }

  update_density_fields(p);
  return ret;
}

int mc_center_mass(Phase*const p, const unsigned int nsteps,const unsigned int tuning_parameter)
{
  assert(p->cm_a);
  //Shutup compiler warning
  unsigned int step=tuning_parameter;step=0;
  // Loop over the MC scweeps
  for (step = 0; step < nsteps; step++)
    {
      uint64_t n_polymers = p->n_polymers ;
      unsigned int n_accepts = 0;
      //#pragma acc parallel loop vector_length(tuning_parameter) reduction(+:n_accepts)
#pragma acc parallel loop vector_length(tuning_parameter)
#pragma omp parallel for reduction(+:n_accepts)
      for (uint64_t npoly = 0; npoly < n_polymers; npoly++)
	{
	  Polymer * mypoly = &p->polymers[npoly];
	  unsigned int myN = p->poly_arch[p->poly_type_offset[mypoly->type]];
	  RNG_STATE rng_state_local = mypoly->poly_state;
	  RNG_STATE * myrngstate = &rng_state_local;
	  const unsigned int poly_type = mypoly->type;
	  if( p->cm_a[poly_type] > 0)
	    {
	      enum enum_pseudo_random_number_generator arg_rng_type;
	      arg_rng_type = p->args.pseudo_random_number_generator_arg;


	      //Generate a random displacement for the center of mass.
	      soma_scalar_t dx,dy,dz;
	      trial_move_cm(p,poly_type,&dx,&dy,&dz,arg_rng_type,myrngstate);

	      soma_scalar_t delta_energy = 0;
	      int move_allowed = 1;

	      //#pragma acc loop vector reduction(+:delta_energy) reduction(|:move_allowed)
	      //Unfortunately this has to be a seq loop, because the reduction crashes.
#pragma acc loop vector seq
	      for (unsigned int ibead = 0; ibead < myN; ibead++)
		{
		  const Monomer  mybead = mypoly->beads[ibead];
		  const unsigned int iwtype =get_particle_type(p->poly_arch[ p->poly_type_offset[poly_type]+1+ibead]);

		  const int tmp = possible_move_area51(p, mybead.x,mybead.y,mybead.z, dx,dy,dz,p->args.nonexact_area51_flag);
		  move_allowed &= tmp;

		  if ( tmp  ){
		    // calculate energy change
		    delta_energy += calc_delta_energy(p, npoly,&mybead, ibead, dx, dy, dz,iwtype);
		  }
		}

	      //Accept Monte-Carlo call
	      if (move_allowed && som_accept(myrngstate, arg_rng_type ,delta_energy) == 1)
		{
#ifndef _OPENACC
		  n_accepts += 1;
#endif//_OPENACC

		  //#pragma acc loop vector
		  //See above
#pragma acc loop seq
		  for (unsigned int ibead = 0; ibead < myN; ibead++)
		    {
		      Monomer  mybead = mypoly->beads[ibead];
		      Monomer*const mybead_ptr = &mypoly->beads[ibead];

		      mybead.x += dx;
		      mybead.y += dy;
		      mybead.z += dz;
		      *mybead_ptr = mybead;
		    }
		}

	      //Copy back the modified RNG state.
	      mypoly->poly_state = rng_state_local;
	    }
	}
      /* p->time += 1; */
      /* p->n_moves += p->num_all_beads_local; */
      /* p->n_accepts += n_accepts; */
    }

  return 0;
}

int mc_polymer_iteration(Phase * const p, const unsigned int nsteps,const unsigned int tuning_parameter)
{
  //Shutup compiler warning
  unsigned int step=tuning_parameter;step=0;


  // Loop over the MC scweeps
  for (step = 0; step < nsteps; step++)
    {
      uint64_t n_polymers = p->n_polymers ;
      unsigned int n_accepts = 0;

      //#pragma acc parallel loop vector_length(tuning_parameter) reduction(+:n_accepts)
#pragma acc parallel loop vector_length(tuning_parameter)
#pragma omp parallel for reduction(+:n_accepts)
      for (uint64_t npoly = 0; npoly < n_polymers; npoly++)
	{
	  unsigned int accepted_moves_loc = 0;

	  // Rebuild bond information for this chain from bonds, or stay with linear right now?
	  Polymer * mypoly = &p->polymers[npoly];

	  unsigned int myN = p->poly_arch[p->poly_type_offset[mypoly->type]];
	  RNG_STATE * myrngstate = &mypoly->poly_state; // maybe local copy of rngstate
	  enum enum_pseudo_random_number_generator arg_rng_type;
	  arg_rng_type = p->args.pseudo_random_number_generator_arg;

	  // MC sweep for this chain
#pragma acc loop seq
	  for (unsigned int nmc = 0; nmc < myN; nmc++) {

	    soma_scalar_t dx=0, dy=0, dz=0, delta_energy=0;
	    unsigned int ibead;

	    // pick a random bead.
	    ibead = soma_rng_uint( myrngstate, arg_rng_type ) % myN;
	    const unsigned int iwtype =get_particle_type(
							 p->poly_arch[ p->poly_type_offset[mypoly->type]+1+ibead]);

	    Monomer  mybead = mypoly->beads[ibead];
	    Monomer* mybead_ptr = &mypoly->beads[ibead];

	    // roll normal MC trial move or force biased MC move
	    soma_scalar_t smc_deltaE;
	    switch(p->args.move_type_arg){
	    case move_type_arg_TRIAL:
	      trial_move(p, npoly, ibead, &dx, &dy, &dz,iwtype,arg_rng_type,myrngstate);	// normal MC move
	      smc_deltaE=0.0;
	      break;
	    case move_type_arg_SMART:
	      trial_move_smc(p, npoly, ibead, &dx, &dy, &dz, &smc_deltaE, &mybead, myrngstate,arg_rng_type,iwtype);	// force biased move
	      break;
	    case move_type__NULL:
	    default:
	      smc_deltaE=0.0;
	      break;
	    }
	    soma_scalar_t newx = mybead.x+dx;
	    soma_scalar_t newy = mybead.y+dy;
	    soma_scalar_t newz = mybead.z+dz;
	    const int move_allowed = possible_move_area51(p, mybead.x,mybead.y,mybead.z, dx,dy,dz,p->args.nonexact_area51_flag);
	    if ( move_allowed  )
	      {
		delta_energy = calc_delta_energy(p, npoly,&mybead, ibead, dx, dy, dz,iwtype);
		delta_energy+=smc_deltaE;

		// MC roll to accept / reject
		if (som_accept(myrngstate, arg_rng_type ,delta_energy) == 1)
		  {
		    mybead_ptr->x = newx;
		    mybead_ptr->y = newy;
		    mybead_ptr->z = newz;
#ifndef _OPENACC
		    accepted_moves_loc += 1;
#endif//_OPENACC
		  }
	      }
	  }
#ifndef _OPENACC
	  n_accepts += accepted_moves_loc;
#endif//_OPENACC
	}

      p->time += 1;
      p->n_moves += p->num_all_beads_local;
      p->n_accepts += n_accepts;
    }

  return 0;
}
void set_iteration_multi_chain(Phase * const p, const unsigned int nsteps,const unsigned int tuning_parameter,const enum enum_pseudo_random_number_generator my_rng_type,const int nonexact_area51, const int start_chain){
  for(unsigned int step = 0; step < nsteps; step++)
    {
      const uint64_t n_polymers = p->n_polymers ;
      unsigned int n_accepts=0;
#pragma acc parallel loop vector_length(tuning_parameter) async
#pragma omp parallel for reduction(+:n_accepts)
      for (uint64_t npoly = start_chain; npoly < n_polymers; npoly++)
	{
	  unsigned int accepted_moves_poly = 0;
	  const uint32_t*const poly_arch = p->poly_arch;
	  const int*const poly_type_offset = p->poly_type_offset;

	  Polymer *const mypoly = &p->polymers[npoly];
	  const unsigned int poly_type = mypoly->type;

	  const unsigned int myN = p->poly_arch[p->poly_type_offset[mypoly->type]];
	  accepted_moves_poly += 0*myN; // Shutup compiler warning
	  const IndependetSets mySets= p->sets[mypoly->type];

	  //Thanks to the PGI compiler, I have to local copy everything I need.
	  const unsigned int n_sets = mySets.n_sets;
	  const unsigned int*const set_length = mySets.set_length;
	  const unsigned int* const sets = mySets.sets;
	  const unsigned int max_member = mySets.max_member;
	  RNG_STATE * const set_states = mypoly->set_states;
	  unsigned int*const set_permutation = mypoly->set_permutation;

	  //Generate random permutation of the sets
	  //http://www.wikipedia.or.ke/index.php/Permutation
#pragma acc loop seq
	  for(unsigned int i=0; i < n_sets; i++)
	    {
	      const unsigned int d = soma_rng_uint(&(mypoly->poly_state),my_rng_type) % (i+1) ;
	      set_permutation[i] = set_permutation[d];
	      set_permutation[d] = i;
	    }

#pragma acc loop seq
	  for(unsigned int iSet=0; iSet < n_sets; iSet++)
	    {
	      unsigned int accepted_moves_set = 0;
	      const unsigned int set_id = set_permutation[iSet];
	      const unsigned int len = set_length[set_id];
#pragma acc loop vector
	      for(unsigned int iP=0; iP < len; iP++)
		{
		  const unsigned int ibead = sets[ set_id*max_member + iP];		  
		  const unsigned int iwtype =get_particle_type(poly_arch[ poly_type_offset[poly_type]+1+ibead]);
		  //local copy of rngstate. For fast updates of state in register.
		  //assert( iP < p->max_set_members );
		  RNG_STATE my_state = set_states[iP];
		  Monomer  mybead = mypoly->beads[ibead];
		  Monomer dx;
		  dx.x=dx.y=dx.z=0;
		  soma_scalar_t smc_deltaE=0;
		  switch(p->args.move_type_arg){
		  case move_type_arg_TRIAL:
		    trial_move(p, npoly, ibead, &dx.x, &dx.y, &dx.z,iwtype,my_rng_type,&my_state);	// normal MC move
		    smc_deltaE=0.0;
		    break;
		  case move_type_arg_SMART:
		    trial_move_smc(p, npoly, ibead, &dx.x, &dx.y, &dx.z, &smc_deltaE, &mybead, &my_state,my_rng_type,iwtype);	// force biased move
		    break;
		  case move_type__NULL:
		  default:
		    smc_deltaE=0.0;
		    break;
		  }
		  const int move_allowed = possible_move_area51(p, mybead.x,mybead.y,mybead.z, dx.x,dx.y,dx.z,nonexact_area51);
		  if ( move_allowed  )
		    {
		      // calculate energy change
		      const soma_scalar_t delta_energy =
			calc_delta_energy(p, npoly,&mybead, ibead, dx.x, dx.y, dx.z,iwtype) + smc_deltaE;

		      // MC roll to accept / reject
		      if (som_accept(&my_state, my_rng_type ,delta_energy) == 1)
			{
			  Monomer newx;
			  newx.x = mybead.x + dx.x;
			  newx.y = mybead.y + dx.y;
			  newx.z = mybead.z + dx.z;
			  mypoly->beads[ibead] = newx;
#ifndef _OPENACC
			  accepted_moves_set += 1;
#endif//_OPENACC
			}
		    }
		  //Copy the RNGstate back to global memory
		  mypoly->set_states[iP] = my_state;    
		}
#ifndef _OPENACC
	      accepted_moves_poly += accepted_moves_set;
#endif//_OPENACC
	    }
#ifndef _OPENACC
	  n_accepts += accepted_moves_poly;
#endif//_OPENACC
	}
      //p->time += 1;
      p->n_moves += p->num_all_beads_local;
#ifndef _OPENACC
      p->n_accepts += n_accepts;
#endif//_OPENACC
    } 
}

void set_iteration_single_chain(Phase * const p, const unsigned int nsteps,const unsigned int tuning_parameter,const enum enum_pseudo_random_number_generator my_rng_type,const int nonexact_area51,uint64_t chain_i){
  for(unsigned int step = 0; step < nsteps; step++)
    {
      unsigned int n_accepts=0;

      unsigned int accepted_moves_poly = 0;
      const IndependetSets mySets= p->sets[(&p->polymers[chain_i])->type];

      //Thanks to the PGI compiler, I have to local copy everything I need.
      const unsigned int n_sets = mySets.n_sets;
      const unsigned int*const set_length_tmp = mySets.set_length;
      const unsigned int* const sets_tmp = mySets.sets;
      const unsigned int max_member = mySets.max_member;
      RNG_STATE * const set_states_tmp = (&p->polymers[chain_i])->set_states;
      unsigned int*const set_permutation = (&p->polymers[chain_i])->set_permutation;

      //Generate random permutation of the sets
      //http://www.wikipedia.or.ke/index.php/Permutation
     
      for(unsigned int i=0; i < n_sets; i++)
	{
	  const unsigned int d = soma_rng_uint(&((&p->polymers[chain_i])->poly_state),my_rng_type) % (i+1) ;
	  set_permutation[i] = set_permutation[d];
	  set_permutation[d] = i;
	}

      for(unsigned int iSet=0; iSet < n_sets; iSet++)
	{
	  unsigned int accepted_moves_set = 0;
	  const unsigned int set_id = set_permutation[iSet];		    
	  const unsigned int len = set_length_tmp[set_id];

#pragma acc parallel loop vector_length(tuning_parameter) async				  
#pragma omp parallel for reduction(+:n_accepts)
	  for(unsigned int iP=0; iP < len; iP++)
	    {	
	   const uint32_t*const poly_arch = p->poly_arch; 
	   const int*const poly_type_offset = p->poly_type_offset; 
	   Polymer *const mypoly = &p->polymers[chain_i];
	   const unsigned int poly_type = mypoly->type; 
	   const IndependetSets mySets= p->sets[poly_type];
	   const unsigned int*const set_length = mySets.set_length;
	   const unsigned int* const sets = mySets.sets;
	   RNG_STATE * const set_states = mypoly->set_states;

	     
	     
	      const unsigned int myN = p->poly_arch[p->poly_type_offset[mypoly->type]]; 
	      accepted_moves_poly += 0*myN; // Shutup compiler warning 
	      
	      const unsigned int ibead = sets[ set_id*max_member + iP];
	      const unsigned int iwtype =get_particle_type(poly_arch[ poly_type_offset[poly_type]+1+ibead]);

	      RNG_STATE my_state = set_states[iP];
	      Monomer  mybead = mypoly->beads[ibead];
	      
	      Monomer dx;
	      dx.x=dx.y=dx.z=0;
	      soma_scalar_t smc_deltaE=0;
	      switch(p->args.move_type_arg){
	      case move_type_arg_TRIAL:
		trial_move(p, 0, ibead, &dx.x, &dx.y, &dx.z,iwtype,my_rng_type,&my_state);	// normal MC move
		smc_deltaE=0.0;
		break;
	      case move_type_arg_SMART:
		trial_move_smc(p, 0, ibead, &dx.x, &dx.y, &dx.z, &smc_deltaE, &mybead, &my_state,my_rng_type,iwtype);	// force biased move
		break;
	      case move_type__NULL:
	      default:
		smc_deltaE=0.0;
		break;
	      }
	      const int move_allowed = possible_move_area51(p, mybead.x,mybead.y,mybead.z, dx.x,dx.y,dx.z,nonexact_area51);
	      if ( move_allowed  )
		{
		  // calculate energy change
		  const soma_scalar_t delta_energy =
		    calc_delta_energy(p, 0,&mybead, ibead, dx.x, dx.y, dx.z,iwtype) + smc_deltaE;
		  // MC roll to accept / reject
		  if (som_accept(&my_state, my_rng_type ,delta_energy) == 1)
		    {
		      Monomer newx;
		      newx.x = mybead.x + dx.x;
		      newx.y = mybead.y + dx.y;
		      newx.z = mybead.z + dx.z;
		      mypoly->beads[ibead] = newx;
#ifndef _OPENACC
		      accepted_moves_set += 1;
#endif//_OPENACC
		    }
		}	      
	      //Copy the RNGstate back to global memory
	      mypoly->set_states[iP] = my_state;
	    }
		    
#ifndef _OPENACC
	  accepted_moves_poly += accepted_moves_set;
#endif//_OPENACC
	}
#ifndef _OPENACC
      n_accepts += accepted_moves_poly;
#endif//_OPENACC	  
      //p->time += 1;
      p->n_moves += p->num_all_beads_local;
#ifndef _OPENACC
      p->n_accepts += n_accepts;
#endif//_OPENACC
    } 
}
	      

int mc_set_iteration(Phase * const p, const unsigned int nsteps,const unsigned int tuning_parameter){
  // We need to check that, otherwise there is a problem in some iterations.
  // But all occuring errors, if any, are reported to the users.
  const enum enum_pseudo_random_number_generator my_rng_type = p->args.pseudo_random_number_generator_arg;
  const int nonexact_area51=p->args.nonexact_area51_flag  + 0*tuning_parameter; //&Shutup compiler warning.
 
  //reorder the polymers according to their length
  int num_long_chain=p->num_long_chain;
  for(int index=0;index<num_long_chain;index++){
  // printf("number %i\n",num_long_chain);
    set_iteration_single_chain(p,nsteps,tuning_parameter,my_rng_type,nonexact_area51,index);
    }
  if(num_long_chain!=p->n_polymers){
    set_iteration_multi_chain(p,nsteps,tuning_parameter,my_rng_type,nonexact_area51,num_long_chain);
  }
  p->time += 1;
#pragma acc wait

  int ret = 0;
  return ret;
}

void trial_move_smc(const Phase * p, const uint64_t ipoly, const int ibead, soma_scalar_t *const dx, soma_scalar_t *const dy, soma_scalar_t *const dz,
		    soma_scalar_t * smc_deltaE,const Monomer *const mybead, RNG_STATE *const myrngstate, const enum enum_pseudo_random_number_generator arg_rng_type,const unsigned int iwtype)
{
  soma_scalar_t x=mybead->x;
  soma_scalar_t y=mybead->y;
  soma_scalar_t z=mybead->z;

  /* R calculated from A according to: Rossky, Doll and Friedman, J.Chem.Phys 69(10)1978 */
  const soma_scalar_t A=p->A[iwtype];
  const soma_scalar_t R=p->R[iwtype];

  /* calculate forces in current position */
  soma_scalar_t fx=0.0; soma_scalar_t fy=0.0; soma_scalar_t fz=0.0;
  add_bond_forces(p,ipoly,ibead,x,y,z,&fx,&fy,&fz);

  /* generate a normal distributed random vector */
  soma_scalar_t rx, ry, rz;
  soma_normal_vector(myrngstate, arg_rng_type, &rx, &ry, &rz);

  /* combine the random offset with the forces, to obtain Brownian motion */
  *dx = A*fx + rx*R;
  *dy = A*fy + ry*R;
  *dz = A*fz + rz*R;

  /* calculate proposed position */
  x+=*dx;
  y+=*dy;
  z+=*dz;

  /* calculate forces in the proposed position */
  soma_scalar_t nfx=0.0; soma_scalar_t nfy=0.0; soma_scalar_t nfz=0.0;
  add_bond_forces(p,ipoly,ibead,x,y,z,&nfx,&nfy,&nfz);

  /* calculate additional terms for scm energy change */
  *smc_deltaE = 0.0;
  *smc_deltaE += 0.5*((nfx+fx)*(*dx) +
		      (nfy+fy)*(*dy) +
		      (nfz+fz)*(*dz));

  *smc_deltaE += 0.25*A*((nfx*nfx)+(nfy*nfy)+(nfz*nfz) -
			 (fx*fx)-(fy*fy)-(fz*fz));

}

void add_bond_forces(const Phase * p, const uint64_t ipoly, unsigned const int ibead,
                     const soma_scalar_t x, const soma_scalar_t y, const soma_scalar_t z,
                     soma_scalar_t *fx, soma_scalar_t *fy, soma_scalar_t *fz)
{
  soma_scalar_t v1x=0.0,v1y=0.0,v1z=0.0;

  const int start = get_bondlist_offset(p->poly_arch[p->poly_type_offset[p->polymers[ipoly].type] + ibead + 1]);

  if( start > 0)
    {
      int i = start;
      unsigned int end;
      do{
	const uint32_t info = p->poly_arch[i++];
	end = get_end(info);
	const unsigned int bond_type = get_bond_type(info);
	const int offset = get_offset(info);

	const int neighbour_id = ibead + offset;
	const unsigned int jbead = neighbour_id;

	soma_scalar_t scale = 1;
	switch (bond_type)
	  {
	  case HARMONICVARIABLESCALE:
	    scale = p->harmonic_normb_variable_scale;
	    /* intentionally falls through */
	  case HARMONIC:
	    //Empty statement, because a statement after a label
	    //has to come before any declaration
	    ;
	    soma_scalar_t v1x_tmp = calc_bond_length(p->polymers[ipoly].beads[jbead].x,x,p->Lx,p->args.bond_minimum_image_convention_flag);
	    soma_scalar_t v1y_tmp = calc_bond_length(p->polymers[ipoly].beads[jbead].y,y,p->Ly,p->args.bond_minimum_image_convention_flag);
	    soma_scalar_t v1z_tmp = calc_bond_length(p->polymers[ipoly].beads[jbead].z,z,p->Lz,p->args.bond_minimum_image_convention_flag);
	    v1x += v1x_tmp*2.0*p->harmonic_normb *scale;
	    v1y += v1y_tmp*2.0*p->harmonic_normb *scale;
	    v1z += v1z_tmp*2.0*p->harmonic_normb *scale;
	    break;
		    
	  case STIFF:
#ifndef _OPENACC
	    fprintf(stderr,
		    "ERROR: %s:%d stiff bond not yet implemented.\n",
		    __FILE__, __LINE__);
#endif//OPENACC
	    break;


	  default:
#ifndef _OPENACC
	    fprintf(stderr, "ERROR: %s:%d unknow bond type appeared %d\n",
		    __FILE__, __LINE__,bond_type);
#endif//OPENACC
	    break;
	  }
      }while( end == 0);
    }
  *fx += v1x;
  *fy += v1y;
  *fz += v1z;
}

inline int possible_move_area51(const Phase*p,const soma_scalar_t oldx,const soma_scalar_t oldy,const soma_scalar_t oldz, soma_scalar_t dx,soma_scalar_t dy,soma_scalar_t dz,const int nonexact)
{
  if( p->area51 == NULL)
    return 1;

  if( p->area51[ coord_to_index(p, oldx+dx, oldy+dy, oldz+dz) ] != 0)
    return 0;

  if(!nonexact)
    {
      const soma_scalar_t r = sqrt(dx*dx + dy*dy + dz*dz);
      const int num_samples = r/p->max_safe_jump ;
      if(num_samples > 0)
	{
	  dx /= num_samples;     dy /= num_samples;     dz /= num_samples;
	  for(int i = 1 ; i < num_samples+1; i++)
	    {
	      const soma_scalar_t jx = oldx+i*dx;
	      const soma_scalar_t jy = oldy+i*dy;
	      const soma_scalar_t jz = oldz+i*dz;
	      const unsigned int index = p->area51[coord_to_index(p, jx, jy, jz)];
	      if(  index != 0 )
		return 0;
	    }
	}
    }

  return 2;
}
