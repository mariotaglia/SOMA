
#include "constantJ.h"
#include <assert.h>
#include "mpiroutines.h"
#include "soma_config.h"
#include "float.h"

int call_CJ(const struct Phase *const p);


int call_CJ(const struct Phase *const p)
{
unsigned int i;

/* Constant J solver for 1D calculations
 * Using a starting concentration, integrates de NP equation without iteration
 * Iteration is used at the end to enforce the fixed number of ions
 * See notes for the derivation 
*/

  soma_scalar_t *cion = (soma_scalar_t *) malloc(NEQ * sizeof(soma_scalar_t)); //solution for c+
    if (cion == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
  soma_scalar_t *bornS = (soma_scalar_t *) malloc(NEQ * sizeof(soma_scalar_t)); // bornS
    if (cion == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }


    const soma_scalar_t invdz = 1.0/p->deltaz ;



soma_scalar_t  Qpos, Qneg, Qposnew, Qnegnew; 
soma_scalar_t iterror = DBL_MAX;

const soma_scalar_t  vall = p->Lx*p->Ly*p->Lz ;
 
const soma_scalar_t maxiterror = 1e-5 ; // maximum relative iteration error for Qpos and Qneg

Qneg = vall;
Qpos = vall;

while (iterror > maxiterror) {

    Qnegnew = 0.0;
    Qposnew = 0.0; 

 if (p->Nnegions > p->Nposions) {    	

#pragma acc data copyin(Qneg) copyin(Qpos)  
#pragma acc parallel loop present(p[:1]) 
#pragma omp parallel for    
    for (i = 0 ; i < p->n_cells ; i++) {
	p->electric_field[i] = p->rhoF[i];
	p->electric_field[i] += sqrt(p->rhoF[i]*p->rhoF[i] + 4.*p->Nposions*p->Nnegions/Qpos/Qneg*p->exp_born_pos[i]*p->exp_born_neg[i]);
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
	p->electric_field[i] += sqrt(p->rhoF[i]*p->rhoF[i] + 4.*p->Nposions*p->Nnegions/Qpos/Qneg*p->exp_born_pos[i]*p->exp_born_neg[i]);
	p->electric_field[i] = p->electric_field[i] / (2.0*p->Nposions/Qpos*p->exp_born_pos[i]) ;
	p->electric_field[i] = -log(p->electric_field[i]);
     }
 }    

#pragma acc data copy(Qposnew) copy(Qnegnew) 
#pragma acc parallel loop reduction (+:Qposnew) reduction (+:Qnegnew)  
#pragma omp parallel for reduction (+:Qposnew) reduction (+:Qnegnew)
    for (i = 0 ; i < p->n_cells ; i++) {
        Qposnew += exp(-p->electric_field[i])*p->exp_born_pos[i];
        Qnegnew += exp(p->electric_field[i])*p->exp_born_neg[i];
    }

        Qposnew = Qposnew*p->vcell;
        Qnegnew = Qnegnew*p->vcell;

        iterror = fabs(Qpos-Qposnew)/Qpos + fabs(Qneg-Qnegnew)/Qneg;
	Qpos = Qposnew;
	Qneg = Qnegnew;

//        printf("Qpos, Qneg, error, %.3e %.3e %.3e \n", Qpos, Qneg, iterror);
}

#pragma acc parallel loop present(p[:1])   
#pragma omp parallel for  
    for (i = 0 ; i < p->n_cells ; i++) {
        p->npos_field[i] = exp(-p->electric_field[i])*p->exp_born_pos[i]/Qpos*p->Nposions; 
        p->nneg_field[i] = exp(p->electric_field[i])*p->exp_born_neg[i]/Qneg*p->Nnegions  ; 
     }

/*    for (i = 0 ; i < p->n_cells ; i++) {
    printf("i pos neg born rhoF %d %.3e %.3e %.3e %.3e \n",i, p->npos_field[i], p->nneg_field[i], p->exp_born[i], p->rhoF[i]);
    }
*/



  return(0);
}

