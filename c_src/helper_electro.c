
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
         const unsigned int type = get_particle_type(p, polytype, iN); 
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
//    if ((p->args.efieldsolver_arg == efieldsolver_arg_PB)||(p->args.efieldsolver_arg == efieldsolver_arg_PH)) 
//    	call_PB(p);

    if (p->args.efieldsolver_arg == efieldsolver_arg_PB) 
    	call_PB(p);

    if (p->args.efieldsolver_arg == efieldsolver_arg_PH)
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



void calc_exp_noneq(const struct Phase *const p) // Updates exp_noneq = exp(-mu_plus(z))

{
const soma_scalar_t alfa = p->args.noneq_ratio_arg;
soma_scalar_t temp;
unsigned int cell,ix,iy,iz;

  if (p->info_MPI.sim_rank == 0) {
    fprintf(stdout, "calc_exp_noneq: c(L)/c(z) ratio: %f \n ", alfa);
     fflush(stdout);
  }
 

  for (ix = 0 ; ix < p->nx ; ix++) {
     for (iy = 0 ; iy < p->ny ; iy++) {
	for (iz = 0 ; iz < p->nz ; iz++) {
        cell = cell_coordinate_to_index(p, ix, iy, iz);

        temp = 1.0 + ((soma_scalar_t) iz)/(((soma_scalar_t) p->nz)-1.)*(1.0-alfa)/alfa; 
        p->exp_noneq[cell] = temp*temp; 
        }
     }
   }        
}

void calc_born_S(const struct Phase *const p) // calculates uB+ + uB- from exp_born

{
unsigned int cell;

#pragma acc parallel loop present(p[:1])
#pragma omp parallel for    
  for (cell = 0 ; cell < p->n_cells ; cell++) {
          p->born_Sc[cell] = -log(p->exp_born_pos[cell]) -log(p->exp_born_neg[cell]);
  }
}


