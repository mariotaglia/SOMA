
#include "donnan.h"
#include <assert.h>
#include "mpiroutines.h"
#include "soma_config.h"
#include "float.h"

int call_EN(const struct Phase *const p);
int call_NO(const struct Phase *const p);


int call_EN(const struct Phase *const p)
{
unsigned int i, type;
soma_scalar_t  rhoQ[p->n_cells_local]; // total number of charges
soma_scalar_t  Qpos, Qneg, Qposnew, Qnegnew; 
soma_scalar_t iterror = DBL_MAX;

const soma_scalar_t  deltax = p->Lx/((soma_scalar_t) p->nx);
const soma_scalar_t  deltay = p->Ly/((soma_scalar_t) p->ny);
const soma_scalar_t  deltaz = p->Lz/((soma_scalar_t) p->nz);
const soma_scalar_t  vcell = deltax*deltay*deltaz;
const soma_scalar_t  vall = p->Lx*p->Ly*p->Lz ;
 
const soma_scalar_t maxiterror = 1e-6*vall ; // maximum iteration error for Qpos and Qneg

Qneg = vall;
Qpos = vall;

while (iterror > maxiterror) {
    Qnegnew = 0.0;
    Qposnew = 0.0; 

 if (p->Nnegions > p->Nposions) {    	
#pragma omp parallel for    
    for (i = 0 ; i < p->n_cells_local ; i++) {
        rhoQ[i] = 0.0; 
        for (type = 0 ; type < p->n_types; type++) {
                   rhoQ[i] += p->fields_unified[i+p->n_cells_local*type]*p->charges[type];
        } 
        rhoQ[i] = rhoQ[i] / vcell ; // units of charge per Re^3


	p->electric_field[i] = rhoQ[i];
	p->electric_field[i] += sqrt(rhoQ[i]*rhoQ[i] + 4.*p->Nposions*p->Nnegions/Qpos/Qneg) ;
	p->electric_field[i] = p->electric_field[i] / (2.0*p->Nnegions/Qneg) ;
	p->electric_field[i] = log(p->electric_field[i]);
     }
 }

 else {
#pragma omp parallel for  
    for (i = 0 ; i < p->n_cells_local ; i++) {
        rhoQ[i] = 0.0; 
        for (type = 0 ; type < p->n_types; type++) {
                   rhoQ[i] += p->fields_unified[i+p->n_cells_local*type]*p->charges[type];
        } 
        rhoQ[i] = rhoQ[i] / vcell ; // units of charge per Re^3

	p->electric_field[i] = -rhoQ[i];
	p->electric_field[i] += sqrt(rhoQ[i]*rhoQ[i] + 4.*p->Nposions*p->Nnegions/Qpos/Qneg) ;
	p->electric_field[i] = p->electric_field[i] / (2.0*p->Nposions/Qpos) ;
	p->electric_field[i] = -log(p->electric_field[i]);
     }
 }    

#pragma omp parallel for reduction (+:Qposnew) 
    for (i = 0 ; i < p->n_cells_local ; i++) {
        Qposnew += exp(-p->electric_field[i]);
    }

#pragma omp parallel for reduction (+:Qnegnew)
    for (i = 0 ; i < p->n_cells_local ; i++) {
        Qnegnew += exp(p->electric_field[i]);
    }


        Qposnew = Qposnew*vcell;
        Qnegnew = Qnegnew*vcell;

        iterror = fabs(Qpos-Qposnew) + fabs(Qneg-Qnegnew);
	Qpos = Qposnew;
	Qneg = Qnegnew;

//        printf("Qpos, Qneg, error, %.3e %.3e %.3e \n", Qpos, Qneg, iterror);
}
	
  return(0);
}



int call_NO(const struct Phase *const p)
{
unsigned int i;

for (i = 0 ; i < p->n_cells_local ; i++) 
	p->electric_field[i] = 0.0;
  return(0);
}
