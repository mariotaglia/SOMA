
#include "constantJ.h"
#include <assert.h>
#include "mpiroutines.h"
#include "soma_config.h"
#include "float.h"
#include "helper_electro.h"

int call_CJ(struct Phase *const p);


int call_CJ(struct Phase *const p)
{
unsigned int i;

/* Constant J solver for 1D calculations
 * Using a starting concentration, integrates de NP equation without iteration
 * Iteration is used at the end to enforce the fixed number of ions
 * See notes for the derivation 
*/

  soma_scalar_t *cion = (soma_scalar_t *) malloc(p->nz * sizeof(soma_scalar_t)); //solution for c+
    if (cion == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
  soma_scalar_t *diffBs = (soma_scalar_t *) malloc(p->nz * sizeof(soma_scalar_t)); // gradient bornS
    if (diffBs == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }


    const soma_scalar_t invdz = 1.0/p->deltaz ;
    const soma_scalar_t  vall = p->Lx*p->Ly*p->Lz ;
    const soma_scalar_t Jpos = p->args.noneq_curr_arg/p->deltax/p->deltay; // J pos is the flux of cations, units Re^-4. 
						       // The input is the total current of cations, units Re^-2 

// Calc Born_S -> puts it into p->born_Sc
  calc_born_S(p);

// Calc gradient of Born_S 
    for (unsigned int i = 0 ; i < p->n_cells-1 ; i++) {
      diffBs[i] = p->born_Sc[i+1] - p->born_Sc[i];
     }
 
// initial guess
cion[0] = p->Nposions/vall ; // initial guess for ion concentration in first cell

// iterate until the total number of cations in the system is ok

soma_scalar_t iterror = DBL_MAX;
const soma_scalar_t maxiterror = 1e-5 ; // maximum relative iteration error for Nposions

soma_scalar_t Aeq, Beq, Ceq;

while (iterror > maxiterror) {

// not paralelizable c(i+1)<-c(i)
    for (unsigned int i = 0 ; i < p->n_cells-1 ; i++) {

	    Aeq = 1.0;
	    Beq = p->rhoF[i+1];
	    Ceq = -exp(invdz*Jpos/cion[i] - diffBs[i])*(cion[i]*(cion[i]+p->rhoF[i]));

	    /*
	    printf("invdz gradBs %.3e %.3e \n", invdz, gradBs[i]);
	    printf("Ceq 1 2 %.3e %.3e \n", -Jpos*0.5*(p->rhoF[i+1]+cion[i]+p->rhoF[i]), (-invdz+0.5*gradBs[i])*(cion[i]*cion[i]+cion[i]*p->rhoF[i]));
            printf("A B C %d %.3e %.3e %.3e \n",i, Aeq, Beq, Ceq);
            */
	    
            cion[i+1] = (-Beq + sqrt(Beq*Beq - 4.0*Aeq*Ceq))/2.0/Aeq;

//	    printf("%f %f \n", cion[i+1], Jpos);
	    }


    soma_scalar_t sumions = 0.0;
    for (unsigned int i = 0 ; i < p->n_cells ; i++) {
	    sumions += cion[i]; }
    
    sumions = sumions*p->vcell;

    iterror = fabs(p->Nposions-sumions)/p->Nposions;
    cion[0] = cion[0]*p->Nposions/sumions; // decrease or increase cion to achieve sumions = p->Nposions

//    printf("cion, sumions, iterror  %.3e %.3e %.3e \n", cion[0], sumions, iterror);
   
}	
//    printf("CONVERGED \n");

// need to check if solution is equal to equilibrium for J = 0

//exit(0);

#pragma acc data copyin(cions)
#pragma acc parallel loop present(p[:1])
#pragma omp parallel for
    for (unsigned int i = 0 ; i < p->n_cells ; i++) {
   
    p->npos_field[i]=cion[i];
    p->nneg_field[i]=cion[i]+p->rhoF[i];
    p->electric_field[i] = log(p->nneg_field[i])-log(p->exp_born_neg[i]);

    }

/*
    for (unsigned int i = 0 ; i < p->n_cells ; i++) {
    printf("i + - e %d %.3e %.3e %.3e  \n ", i, p->npos_field[i], p->nneg_field[i],  p->electric_field[i]); 
    }
    exit(0);
*/

free(cion);
free(diffBs);
return 0;

}
