
#include "donnan.h"
#include <assert.h>
#include "mpiroutines.h"
#include "soma_config.h"
#include "float.h"

int call_EN(const struct Phase *const p);


int call_EN(const struct Phase *const p)
{
unsigned int i;
soma_scalar_t  Qpos, Qneg, Qposnew, Qnegnew; 
soma_scalar_t iterror_i = DBL_MAX;
soma_scalar_t iterror_o = DBL_MAX;

const soma_scalar_t  vall = p->Lx*p->Ly*p->Lz ; //total volume of simulation box
 
const soma_scalar_t maxiterror_i = 1.e-5 ; // maximum relative iteration error for Qpos and Qneg
const soma_scalar_t maxiterror_o = 1.e-5 ; // maximum relative iteration error for ion densities

const soma_scalar_t v_pos = 4.0/3.0*M_PI*p->Born_pos*p->Born_pos*p->Born_pos; // volume of cations, units Re^3
const soma_scalar_t v_neg = 4.0/3.0*M_PI*p->Born_neg*p->Born_neg*p->Born_neg; // volume of anions, units Re^3
const soma_scalar_t rho0 = p->num_all_beads/vall; // density, beads/Re^3
const soma_scalar_t kappa0 = p->xn[0]/p->reference_Nbeads; // kappa for ions, use the A-A value

soma_scalar_t exprep_pol_pos[p->n_cells], exprep_pol_neg[p->n_cells]; // contribution from ion-polymer repulsions
soma_scalar_t exprep_ion_pos[p->n_cells], exprep_ion_neg[p->n_cells]; // contribution from ion-ion repulsions
								      //
soma_scalar_t exprep_ion_pos_old[p->n_cells], exprep_ion_neg_old[p->n_cells]; // auxiliary fields for iteration

soma_scalar_t exprep_pos[p->n_cells], exprep_neg[p->n_cells]; // total repulsions
soma_scalar_t tmpsumdens[p->n_cells]; // auxiliary field
	
Qneg = vall;
Qpos = vall;


// auxiliary fields for steric interactions
#pragma acc data copyin(omega_rep_pol)  
#pragma acc parallel loop present(p[:1]) 
#pragma omp parallel for    
    for (i = 0 ; i < p->n_cells ; i++) {
    exprep_pol_pos[i] = exp(-p->omega_rep_pol[i]*rho0*v_pos);  
    exprep_pol_neg[i] = exp(-p->omega_rep_pol[i]*rho0*v_neg);  
    exprep_ion_pos[i] = 1; // initialize to one, then iterate
    exprep_ion_neg[i] = 1; // initialize to one, then iterate
     }

////////////////////////////////////////////////////////////////////////////////////////
// Double iteration loop -- 
// Inner loop converges Q
// Outter loop converges rho+ and rho-, which is necessary because of ion-ion repulsions
///////////////////////////////////////////////////////////////////////////////////////


while (iterror_o > maxiterror_o) { // outer
  while (iterror_i > maxiterror_i) { // inner

    Qnegnew = 0.0;
    Qposnew = 0.0; 

 if (p->Nnegions > p->Nposions) {    	

#pragma acc data copyin(Qneg) copyin(Qpos)  
#pragma acc parallel loop present(p[:1]) 
#pragma omp parallel for    
    for (i = 0 ; i < p->n_cells ; i++) {
       
        exprep_pos[i] = exprep_pol_pos[i]*exprep_ion_pos[i];
        exprep_neg[i] = exprep_pol_neg[i]*exprep_ion_neg[i];

	p->electric_field[i] = p->rhoF[i];
	p->electric_field[i] += 
	  sqrt(p->rhoF[i]*p->rhoF[i] + 4.*p->Nposions*p->Nnegions/Qpos/Qneg*p->exp_born_pos[i]*p->exp_born_neg[i]*exprep_pos[i]*exprep_neg[i]);
	p->electric_field[i] = p->electric_field[i] / (2.0*p->Nnegions/Qneg*p->exp_born_neg[i]*exprep_neg[i]) ;
	p->electric_field[i] = log(p->electric_field[i]);
     }
 }

 else {
#pragma acc data copyin(Qneg) copyin(Qpos)  
#pragma acc parallel loop present(p[:1]) 
#pragma omp parallel for  
    for (i = 0 ; i < p->n_cells ; i++) {

        exprep_pos[i] = exprep_pol_pos[i]*exprep_ion_pos[i];
        exprep_neg[i] = exprep_pol_neg[i]*exprep_ion_neg[i];

	p->electric_field[i] = -p->rhoF[i];
	p->electric_field[i] += 
	  sqrt(p->rhoF[i]*p->rhoF[i] + 4.*p->Nposions*p->Nnegions/Qpos/Qneg*p->exp_born_pos[i]*p->exp_born_neg[i]*exprep_pos[i]*exprep_neg[i]);
	p->electric_field[i] = p->electric_field[i] / (2.0*p->Nposions/Qpos*p->exp_born_pos[i]*exprep_pos[i]) ;
	p->electric_field[i] = -log(p->electric_field[i]);
     }
 }    

#pragma acc data copy(Qposnew) copy(Qnegnew) 
#pragma acc parallel loop reduction (+:Qposnew) reduction (+:Qnegnew)  
#pragma omp parallel for reduction (+:Qposnew) reduction (+:Qnegnew)
    for (i = 0 ; i < p->n_cells ; i++) {
        Qposnew += exp(-p->electric_field[i])*p->exp_born_pos[i]*exprep_pos[i];
        Qnegnew += exp(p->electric_field[i])*p->exp_born_neg[i]*exprep_neg[i];
    }

        Qposnew = Qposnew*p->vcell;
        Qnegnew = Qnegnew*p->vcell;

        iterror_i = fabs(Qpos-Qposnew)/Qpos + fabs(Qneg-Qnegnew)/Qneg;
	Qpos = Qposnew;
	Qneg = Qnegnew;

        printf("Qpos, Qneg, error, %.3e %.3e %.3e \n", Qpos, Qneg, iterror_i);
} // inner iteration loop, Q is converged

iterror_i = DBL_MAX;

// calculate ion densities
#pragma acc parallel loop present(p[:1])   
#pragma omp parallel for  
    for (i = 0 ; i < p->n_cells ; i++) {
        p->npos_field[i] = exp(-p->electric_field[i])*p->exp_born_pos[i]*exprep_pos[i]/Qpos*p->Nposions; 
        p->nneg_field[i] = exp(p->electric_field[i])*p->exp_born_neg[i]*exprep_neg[i]/Qneg*p->Nnegions  ; 
     }

/*    for (i = 0 ; i < p->n_cells ; i++) {
    printf("i pos neg born rhoF %d %.3e %.3e %.3e %.3e \n",i, p->npos_field[i], p->nneg_field[i], p->exp_born[i], p->rhoF[i]);
    }
*/

// recalculate auxiliary fields 
#pragma acc parallel loop present(p[:1]) 
#pragma omp parallel for    
    for (i = 0 ; i < p->n_cells ; i++) {
    tmpsumdens[i] =  p->npos_field[i]*v_pos+p->nneg_field[i]*v_neg;
    exprep_ion_pos_old[i] = exprep_ion_pos[i]; 
    exprep_ion_neg_old[i] = exprep_ion_neg[i]; 
    exprep_ion_pos[i] = exp(-kappa0*rho0*v_pos*tmpsumdens[i]); 
    exprep_ion_neg[i] = exp(-kappa0*rho0*v_neg*tmpsumdens[i]); 
    }

        iterror_o = 0;
// Outer loop iteration
    for (i = 0 ; i < p->n_cells ; i++) {
        iterror_o += fabs((exprep_ion_pos[i]-exprep_ion_pos_old[i])/exprep_ion_pos[i]);
        iterror_o += fabs((exprep_ion_neg[i]-exprep_ion_neg_old[i])/exprep_ion_neg[i]);
    }

        printf("error outer, %.3e \n", iterror_o);
} // outer loop iteration, exprep_ion is converged

  exit(0); // still need to implement ion-pol repulsions in omega_field
           // copyout in last loop?
  return(0);
}

