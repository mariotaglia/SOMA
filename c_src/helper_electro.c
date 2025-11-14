
#include "helper_electro.h"
#include "mesh.h"
#include <stdbool.h>
#include <assert.h>
#include "mpiroutines.h"
#include "soma_config.h"
#include "cmdline.h"

/*  Helper routines for efield calculation */

void calc_ions(struct Phase *const p)
  {

  soma_scalar_t sumrhoQ;                   // total number of charges
  soma_scalar_t sumrhoQtmp = 0.0;                   // total number of charges
  unsigned int polytype ;

  for (uint64_t i = 0; i < p->n_polymers; i++)   { ; // loope over pol chains
    polytype = p->polymers[i].type; // polymer type
    const unsigned N = p->poly_arch[p->poly_type_offset[polytype]];
    unsigned int iN;   
    for (iN = 0; iN < N ; iN++) {  // loop over monomers
         const unsigned int type = get_particle_type(p, i, iN);
//	fprintf(stdout, "Hola %d %d type:%d \n ", i, iN, type); 
         sumrhoQtmp += p->charges[type]; 
    }
  }

#if ( ENABLE_MPI == 1 )
    MPI_Allreduce(&sumrhoQtmp, &sumrhoQ, 1,  MPI_SOMA_SCALAR, MPI_SUM, p->info_MPI.SOMA_comm_sim);
#endif                          //ENABLE_MPI

  p->Nposions = p->Nnegions = p->Nions;

  if (sumrhoQ > 0.0) { 
      p->Nnegions += sumrhoQ; 
      }
  else {
      p->Nposions += -sumrhoQ;
  }

  soma_scalar_t Nionsdiff = p->Nposions-p->Nnegions;
  soma_scalar_t  netcharge = Nionsdiff + sumrhoQ;

  if (p->info_MPI.sim_rank == 0) {
    fprintf(stdout, "calc_ions: Total charge of fixed beads: %f \n ", sumrhoQ);
    fprintf(stdout, "calc_ions: Total number of added salt ions: %f \n ", p->Nions);
    fprintf(stdout, "calc_ions: Total number of +1 salt ions: %f \n ", p->Nposions);
    fprintf(stdout, "calc_ions: Total number of -1 salt ions: %f \n ", p->Nnegions);
    fprintf(stdout, "check_electro: Net charge:  %f \n ", netcharge);
     fflush(stdout);
  }
 
  assert(fabs(netcharge) < 1.0e-6);
}

void update_electric_field(const struct Phase *const p)
{
    // Update electric potential
    //
//    if ((p->args.efieldsolver_arg == efieldsolver_arg_PB)||(p->args.efieldsolver_arg == efieldsolver_arg_NP)) 
//    	call_PB(p);

    if (p->args.efieldsolver_arg == efieldsolver_arg_PB) 
    	call_PB(p);

    if (p->args.efieldsolver_arg == efieldsolver_arg_NP)
    	call_J(p);

    if (p->args.efieldsolver_arg == efieldsolver_arg_EN) 
    	call_EN(p);
}  

void calc_invbls(struct Phase *const p) 
{
unsigned int type;

  p->invbls = (soma_scalar_t *) malloc((p->n_types) * sizeof(soma_scalar_t));
  p->invblav_zero = 0.0; 

  for (type = 0; type < p->n_types ; type++) {
	p->invbls[type] = 1./p->bls[type];
	p->invblav_zero += 1./p->bls[type];
  }
        p->invblav_zero = p->invblav_zero/((soma_scalar_t) p->n_types);
}


void update_invblav(const struct Phase *const p) // Updates invblav = average of inverse Bjerrum length

{
unsigned int tmpsegsum[p->n_cells];
unsigned int cell, type;


#pragma data create(tmpsegsumacc[0:p->n_cells])  
{

#pragma acc parallel loop present(p[:1]) 
#pragma omp parallel for
for (cell = 0 ; cell < p->n_cells ; cell++) {
     p->invblav[cell] = 0.0 ; 
     tmpsegsum[cell] = 0;
}

for (type = 0 ; type < p->n_types; type++) {

#pragma acc parallel loop present(p[:1])
#pragma omp parallel for
    for (cell = 0 ; cell < p->n_cells ; cell++) {

       p->invblav[cell] += ((soma_scalar_t) p->fields_unified[cell+p->n_cells*type])*p->invbls[type];
       tmpsegsum[cell] +=  p->fields_unified[cell+p->n_cells*type];
 } 
}

#pragma acc parallel loop present(p[:1])
#pragma omp parallel for
    for (cell = 0 ; cell < p->n_cells ; cell++) {
           if(tmpsegsum[cell] != 0)    
                p->invblav[cell] = p->invblav[cell] / ((soma_scalar_t) tmpsegsum [cell]);
	   else if (tmpsegsum[cell] == 0) 
                p->invblav[cell] = p->invblav_zero; // prevents divergence if tmpsegsum = 0
	  }
    }
}

void update_d_invblav(const struct Phase *const p) // Updates d_invblav = derivative of average inverse Bjerrum length respect to N_i

{
unsigned int cell, type;
unsigned int tmpsegsum[p->n_cells];

#pragma data create(tmpsegsumacc[0:p->n_cells])  
{

#pragma acc parallel loop present(p[:1]) 
#pragma omp parallel for
for (cell = 0 ; cell < p->n_cells ; cell++) {
     tmpsegsum[cell] = 0;
}

for (type = 0 ; type < p->n_types; type++) {
#pragma acc parallel loop present(p[:1])
#pragma omp parallel for
    for (cell = 0 ; cell < p->n_cells ; cell++) {
       tmpsegsum[cell] +=  p->fields_unified[cell+p->n_cells*type];
 } 
}


#pragma acc parallel loop present(p[:1]) 
#pragma omp parallel for
for (cell = 0 ; cell < p->n_cells ; cell++) {
	if(tmpsegsum[cell] == 0) 
	   tmpsegsum[cell] = 1;  // prevents divergence of de/dN when N = 0
}


for (type = 0 ; type < p->n_types; type++) {
#pragma acc parallel loop present(p[:1])
#pragma omp parallel for
    for (cell = 0 ; cell < p->n_cells ; cell++) {
       p->d_invblav[cell+p->n_cells_local*type] = (p->invbls[type] - p->invblav[cell])/((soma_scalar_t) tmpsegsum[cell]);
  } 
 }
} // pragma create

}

void update_exp_born(const struct Phase *const p) // Updates exp_born = exp(-u_B)

{
unsigned int cell;
soma_scalar_t borntmp;

#pragma acc parallel loop present(p[:1])
#pragma omp parallel for    
for (cell = 0 ; cell < p->n_cells_local ; cell++) {

   borntmp = 1.0/(p->invblav[cell]*2.0*p->Born_pol);
   borntmp += -1.0/(p->invblav_zero*2.0*p->Born_pol); // shift to prevent very large numbers...
   p->exp_born_pol[cell] = exp(-borntmp); 
}

#pragma acc parallel loop present(p[:1])
#pragma omp parallel for    
for (cell = 0 ; cell < p->n_cells_local ; cell++) {

   borntmp = 1.0/(p->invblav[cell]*2.0*p->Born_pos);
   borntmp += -1.0/(p->invblav_zero*2.0*p->Born_pos); // shift to prevent very large numbers...
   p->exp_born_pos[cell] = exp(-borntmp); 
}

#pragma acc parallel loop present(p[:1])
#pragma omp parallel for    
for (cell = 0 ; cell < p->n_cells_local ; cell++) {

   borntmp = 1.0/(p->invblav[cell]*2.0*p->Born_neg);
   borntmp += -1.0/(p->invblav_zero*2.0*p->Born_neg); // shift to prevent very large numbers...
   p->exp_born_neg[cell] = exp(-borntmp); 
}
}

void update_rhoF(const struct Phase *const p) // Updates rhoF (charge per Re^3) 

{
unsigned int cell, type;

#pragma acc parallel loop present(p[:1]) 
#pragma omp parallel for
for (cell = 0 ; cell < p->n_cells ; cell++) {
     p->rhoF[cell] = 0;
}

for (type = 0 ; type < p->n_types; type++) {
#pragma acc parallel loop present(p[:1])
#pragma omp parallel for    
    for (cell = 0 ; cell < p->n_cells_local ; cell++) {
        p->rhoF[cell] += ((soma_scalar_t) p->fields_unified[cell+p->n_cells_local*type])*p->charges[type]/p->vcell;
   } 
}
}

void update_NB(const struct Phase *const p) // Updates NB (|charge| per Re^3, including counter ions) 

{
unsigned int cell, type;

#pragma acc parallel loop present(p[:1]) 
#pragma omp parallel for
for (cell = 0 ; cell < p->n_cells ; cell++) {
     p->NB[cell] = 0;
}

for (type = 0 ; type < p->n_types; type++) {
#pragma acc parallel loop present(p[:1])
#pragma omp parallel for    
    for (cell = 0 ; cell < p->n_cells_local ; cell++) {
        p->NB[cell] += ((soma_scalar_t) p->fields_unified[cell+p->n_cells_local*type])*fabs(p->charges[type]);
   } 
}

#pragma acc parallel loop present(p[:1])
#pragma omp parallel for    
    for (cell = 0 ; cell < p->n_cells_local ; cell++) {
        p->NB[cell] += p->Born_pol*(p->nneg_field[cell]/p->Born_neg + p->npos_field[cell]/p->Born_pos)*p->vcell; 
   } 
}



void calc_born_S(struct Phase *const p) // calculates uB+ + uB- from exp_born

{
unsigned int cell;

#pragma acc parallel loop present(p[:1])
#pragma omp parallel for    
  for (cell = 0 ; cell < p->n_cells ; cell++) {
          p->born_Sc[cell] = -log(p->exp_born_pos[cell]) -log(p->exp_born_neg[cell]);
  }
}

void calc_J_umbrella(const struct Phase *const p) // calculates J fluxes and put it into umbrella field for export

{
  soma_scalar_t  eps[p->n_cells]; // c/ceq
  soma_scalar_t  cions[p->n_cells]; // c/ceq
  soma_scalar_t  Jx[p->n_cells]; // c/ceq
  soma_scalar_t  Jy[p->n_cells]; // c/ceq
  soma_scalar_t  Jz[p->n_cells]; // c/ceq
  unsigned int ix,iy,iz,i,cell; 
  unsigned int ixm,iym,izm,ixp,iyp,izp, cellm, cellp;
 
  int mod(int a, int b); // modulus

    for (cell = 0 ; cell < p->n_cells ; cell++) {
	      cions[cell] = p->npos_field[cell]; // save solution in non-eq
    }


   // Calc ions in equilibrium
   call_EN(p);

    for (cell = 0 ; cell < p->n_cells ; cell++) {
	      eps[cell] = cions[cell]/p->npos_field[cell]; // ratio non-eq/eq
    }


// Calculation of ion fluxes

  
    for (ix = 0 ; ix < p->nx ; ix++) {

     ixp = mod((ix+1),p->nx);
     ixm = mod((ix-1),p->nx);
 
     for (iy = 0 ; iy < p->ny ; iy++) {

        iyp = mod((iy+1),p->ny);
        iym = mod((iy-1),p->ny);

	for (iz = 1 ; iz < p->nz-1 ; iz++) {
      
	izp = iz+1;       	
	izm = iz-1;       	
 
// Jz
	 cellm = cell_coordinate_to_index(p, ix, iy, izm);
	 cellp = cell_coordinate_to_index(p, ix, iy, izp);
	 cell = cell_coordinate_to_index(p, ix, iy, iz);
         Jz[cell] =  -(p->npos_field[cell])*(eps[cellp]-eps[cellm])/p->deltaz;
// Jx
	 cellm = cell_coordinate_to_index(p, ixm, iy, iz);
	 cellp = cell_coordinate_to_index(p, ixp, iy, iz);
	 cell = cell_coordinate_to_index(p, ix, iy, iz);
         Jx[cell] =  -(p->npos_field[cell])*(eps[cellp]-eps[cellm])/p->deltax;
// Jy
	 cellm = cell_coordinate_to_index(p, ix, iym, iz);
	 cellp = cell_coordinate_to_index(p, ix, iyp, iz);
	 cell = cell_coordinate_to_index(p, ix, iy, iz);
         Jy[cell] =  -(p->npos_field[cell])*(eps[cellp]-eps[cellm])/p->deltay;
			    
          }
     }
  }


    for (ix = 0 ; ix < p->nx ; ix++) {
     for (iy = 0 ; iy < p->ny ; iy++) {
	 iz = 0;       	
	 cell = cell_coordinate_to_index(p, ix, iy, iz);
         Jz[cell] = 0.0;
         Jx[cell] = 0.0;
         Jy[cell] = 0.0;

	 iz = p->nz-1;       	
	 cell = cell_coordinate_to_index(p, ix, iy, iz);
         Jz[cell] = 0.0;
         Jx[cell] = 0.0;
         Jy[cell] = 0.0;
     }
  }


// save to umbrella field and restore non-eq solution

  for (cell = 0 ; cell < p->n_cells ; cell++) {

      p->umbrella_field[cell] = Jz[cell];
      p->umbrella_field[cell+p->n_cells] = sqrt(Jz[cell]*Jz[cell]+Jx[cell]*Jx[cell]+Jy[cell]*Jy[cell]);
      p->npos_field[cell] = cions[cell];
  }

}


