
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
soma_scalar_t iterror = DBL_MAX;

const soma_scalar_t  deltax = p->Lx/((soma_scalar_t) p->nx);
const soma_scalar_t  deltay = p->Ly/((soma_scalar_t) p->ny);
const soma_scalar_t  deltaz = p->Lz/((soma_scalar_t) p->nz);
const soma_scalar_t  vcell = deltax*deltay*deltaz;
const soma_scalar_t  vall = p->Lx*p->Ly*p->Lz ;
 
const soma_scalar_t maxiterror = 1e-5 ; // maximum relative iteration error for Qpos and Qneg

Qneg = vall;
Qpos = vall;

    printf("Hola \n");

while (iterror > maxiterror) {
    Qnegnew = 0.0;
    Qposnew = 0.0; 

 if (p->Nnegions > p->Nposions) {    	
//#pragma omp parallel for    
    for (i = 0 ; i < p->n_cells_local ; i++) {
	p->electric_field[i] = p->rhoF[i];
	p->electric_field[i] += sqrt(p->rhoF[i]*p->rhoF[i] + 4.*p->Nposions*p->Nnegions/Qpos/Qneg*p->exp_born[i]*p->exp_born[i]);
	p->electric_field[i] = p->electric_field[i] / (2.0*p->Nnegions/Qneg*p->exp_born[i]) ;
	p->electric_field[i] = log(p->electric_field[i]);
     }
 }

 else {
//#pragma omp parallel for  
    for (i = 0 ; i < p->n_cells_local ; i++) {
	p->electric_field[i] = -p->rhoF[i];
	p->electric_field[i] += sqrt(p->rhoF[i]*p->rhoF[i] + 4.*p->Nposions*p->Nnegions/Qpos/Qneg*p->exp_born[i]*p->exp_born[i]);
	p->electric_field[i] = p->electric_field[i] / (2.0*p->Nposions/Qpos*p->exp_born[i]) ;
	p->electric_field[i] = -log(p->electric_field[i]);
     }
 }    

//#pragma omp parallel for reduction (+:Qposnew) 
    for (i = 0 ; i < p->n_cells_local ; i++) {
        Qposnew += exp(-p->electric_field[i])*p->exp_born[i];
    }

//#pragma omp parallel for reduction (+:Qnegnew)
    for (i = 0 ; i < p->n_cells_local ; i++) {
        Qnegnew += exp(p->electric_field[i])*p->exp_born[i];
    }

        Qposnew = Qposnew*vcell;
        Qnegnew = Qnegnew*vcell;

        iterror = fabs(Qpos-Qposnew)/Qpos + fabs(Qneg-Qnegnew)/Qneg;
	Qpos = Qposnew;
	Qneg = Qnegnew;

//        printf("Qpos, Qneg, error, %.3e %.3e %.3e \n", Qpos, Qneg, iterror);
}

//#pragma omp parallel for  
    for (i = 0 ; i < p->n_cells_local ; i++) {
        p->npos_field[i] = exp(-p->electric_field[i])*p->exp_born[i]/Qpos*p->Nposions ; 
        p->nneg_field[i] = exp(p->electric_field[i])*p->exp_born[i]/Qneg*p->Nnegions  ; 
     }

/*    for (i = 0 ; i < p->n_cells_local ; i++) {
    printf("i pos neg born rhoF %d %.3e %.3e %.3e %.3e \n",i, p->npos_field[i], p->nneg_field[i], p->exp_born[i], p->rhoF[i]);
    }
*/



  return(0);
}

