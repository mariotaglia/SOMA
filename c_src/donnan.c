
#include "donnan.h"
#include <assert.h>
#include "mpiroutines.h"
#include "soma_config.h"
#include "float.h"

int call_EN(const struct Phase *const p);


int call_EN(const struct Phase *const p)
{
unsigned int i;
soma_scalar_t  Qpos, Qneg, Qpos_new, Qneg_new; 
soma_scalar_t iterror_1 = DBL_MAX;
soma_scalar_t iterror_2 = DBL_MAX;
const soma_scalar_t alfa = 0.25; // mixing factor in iteration

const soma_scalar_t  vall = p->Lx*p->Ly*p->Lz ; //total volume of simulation box
 
const soma_scalar_t maxiterror_1 = 1.e-5 ; // maximum relative iteration error for Qpos and Qneg
const soma_scalar_t maxiterror_2 = 1.e-5 ; // maximum relative iteration error for ion densities

const soma_scalar_t v_pos = 4.0/3.0*M_PI*p->Born_pos*p->Born_pos*p->Born_pos; // volume of cations, units Re^3
const soma_scalar_t v_neg = 4.0/3.0*M_PI*p->Born_neg*p->Born_neg*p->Born_neg; // volume of anions, units Re^3

const soma_scalar_t v_salt = (v_pos+v_neg)/2.0;


const soma_scalar_t rho0 = p->num_all_beads/vall; // density, beads/Re^3
const soma_scalar_t kappa0 = p->xn[0]/p->reference_Nbeads; // kappa for ions, use the A-A value

soma_scalar_t *exprep_pol = (soma_scalar_t *) malloc(p->n_cells * sizeof(soma_scalar_t));
    if (exprep_pol == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

soma_scalar_t *exprep_ion = (soma_scalar_t *) malloc(p->n_cells * sizeof(soma_scalar_t));
    if (exprep_ion == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

soma_scalar_t *exprep_ion_new = (soma_scalar_t *) malloc(p->n_cells * sizeof(soma_scalar_t));
    if (exprep_ion_new == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }


soma_scalar_t *tmpsumdens = (soma_scalar_t *) malloc(p->n_cells * sizeof(soma_scalar_t));
    if (tmpsumdens == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }


// initial guess for Q     
Qneg = vall;
Qpos = vall;


// exprep_pol does not change during iteration
#pragma acc parallel loop present(p[:1])   
#pragma omp parallel for  
    for (i = 0 ; i < p->n_cells ; i++) {
        exprep_pol[i] = exp(-p->omega_rep_pol[i]*rho0*v_salt);  
    }

 
// initial guess for exprep_ion
#pragma acc parallel loop present(p[:1])   
#pragma omp parallel for  
    for (i = 0 ; i < p->n_cells ; i++) {
        exprep_ion[i] = 1; // assume no ion-ion repulsions initially
    }


while ((iterror_1 > maxiterror_1) || (iterror_2 > maxiterror_2)) {

    Qneg_new = 0.0;
    Qpos_new = 0.0; 

 if (p->Nnegions > p->Nposions) {    	

#pragma acc data copyin(Qneg) copyin(Qpos)  
#pragma acc parallel loop present(p[:1]) 
#pragma omp parallel for    
    for (i = 0 ; i < p->n_cells ; i++) {
	p->electric_field[i] = p->rhoF[i];
	p->electric_field[i] += 
        sqrt(p->rhoF[i]*p->rhoF[i] + 4.*p->Nposions*p->Nnegions/Qpos/Qneg*p->exp_born_pos[i]*p->exp_born_neg[i]);
	p->electric_field[i] = p->electric_field[i] / (2.0*p->Nnegions/Qneg*p->exp_born_neg[i]) ;
	p->electric_field[i] = log(p->electric_field[i]);
     }
 }

 else {
#pragma acc data copyin(Qneg) copyin(Qpos)  
#pragma acc parallel loop present(p[:1]) 
#pragma omp parallel for  
    for (i = 0 ; i < p->n_cells ; i++) {
	p->electric_field[i] = -p->rhoF[i];
	p->electric_field[i] += 
	sqrt(p->rhoF[i]*p->rhoF[i] + 4.*p->Nposions*p->Nnegions/Qpos/Qneg*p->exp_born_pos[i]*p->exp_born_neg[i]);
	p->electric_field[i] = p->electric_field[i] / (2.0*p->Nposions/Qpos*p->exp_born_pos[i]) ;
	p->electric_field[i] = -log(p->electric_field[i]);
     }
 }    

#pragma acc data copy(Qpos_new) copy(Qneg_new) 
#pragma acc parallel loop reduction (+:Qpos_new) reduction (+:Qneg_new)  
#pragma omp parallel for reduction (+:Qpos_new) reduction (+:Qneg_new)
    for (i = 0 ; i < p->n_cells ; i++) {
        exprep_pol[i] = exp(-p->omega_rep_pol[i]*rho0*v_salt);  

        Qpos_new += exp(-p->electric_field[i])*p->exp_born_pos[i]*exprep_pol[i]*exprep_ion[i];
        Qneg_new += exp( p->electric_field[i])*p->exp_born_neg[i]*exprep_pol[i]*exprep_ion[i];
    }

        Qpos_new = Qpos_new*p->vcell;
        Qneg_new = Qneg_new*p->vcell;

// ion densities
#pragma acc parallel loop present(p[:1])   
#pragma omp parallel for  
    for (i = 0 ; i < p->n_cells ; i++) {
        p->npos_field[i] = exp(-p->electric_field[i])*p->exp_born_pos[i]*exprep_pol[i]*exprep_ion[i]/Qpos_new*p->Nposions; 
        p->nneg_field[i] = exp( p->electric_field[i])*p->exp_born_neg[i]*exprep_pol[i]*exprep_ion[i]/Qneg_new*p->Nnegions; 
        tmpsumdens[i] =  (p->npos_field[i]+p->nneg_field[i])*v_salt;
        exprep_ion_new[i] = exp(-kappa0*rho0*v_salt*tmpsumdens[i]);
    }

// Calculate error norms

// Q    
iterror_1 = fabs(Qpos-Qpos_new)/Qpos + fabs(Qneg-Qneg_new)/Qneg;
Qpos = Qpos_new*alfa + Qpos*(1.-alfa);
Qneg = Qneg_new*alfa + Qneg*(1.-alfa);

// ion-ion repulsion Boltzmann factor
iterror_2 = 0.0;
#pragma acc data copy(iterror_2)
#pragma acc parallel loop reduction (+:iterror_2)
    for (i = 0 ; i < p->n_cells ; i++) {
        iterror_2 += fabs((exprep_ion[i]-exprep_ion_new[i])/exprep_ion[i]);
    }

#pragma acc parallel loop present(p[:1])   
#pragma omp parallel for  
    for (i = 0 ; i < p->n_cells ; i++) {
        exprep_ion[i] = exprep_ion_new[i]*alfa + exprep_ion[i]*(1.-alfa);
    }


//        printf("Error 1, Error 2 : %.3e  %.3e \n", iterror_1, iterror_2);

} // iteration loop


//        printf("Converged \n");
  return(0);
}

