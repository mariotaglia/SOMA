
#include "donnan.h"
#include <assert.h>
#include "mpiroutines.h"
#include "soma_config.h"
#include "float.h"
#include "mesh.h"

int call_EN(const struct Phase *const p);


int call_EN(const struct Phase *const p)
{
unsigned int i, ix, iy, iz, cell;
soma_scalar_t  Qpos, Qneg, Qposnew, Qnegnew; 
soma_scalar_t iterror = DBL_MAX;

const soma_scalar_t falfa = p->args.noneq_ratio_arg*p->deltaz*((soma_scalar_t) p->nx)*((soma_scalar_t) p->ny) ;

const soma_scalar_t  vall = p->Lx*p->Ly*p->Lz ;
const soma_scalar_t maxiterror = 1e-5 ; // maximum relative iteration error for Qpos and Qneg

soma_scalar_t cpos_av[p->nz], mupos[p->nz], expmupos_sc[p->n_cells];

Qneg = vall;
Qpos = vall;

/* Calculate cpos_av from previous iteration */

#pragma acc data copyout(cpos_av[:])
#pragma acc parallel loop present(p[:1]) 
#pragma omp parallel for    
for (iz = 0 ; iz < p->nz ; iz++) {

cpos_av[iz] = 0.0;

   for (ix = 0 ; ix < p->nx ; ix++) {
     for (iy = 0 ; iy < p->ny ; iy++) {
      
     cell = cell_coordinate_to_index(p, ix, iy, iz);
     cpos_av[iz] += p->npos_field[cell];
    
     }
  }
}


/* Calculate mu pos, do not paralelize */

mupos[0] = 0.0;
mupos[1] = falfa/cpos_av[0];

/* Integration based on trapezoid rule */

for (iz = 2 ; iz < p->nz ; iz++) {
   mupos[iz] = mupos[iz-1] + 0.5*falfa*(1.0/cpos_av[iz]+1.0/cpos_av[iz-1]) ; 
}

//#pragma acc data copyin(mupos)
//#pragma acc parallel loop present(p[:1]) 
//#pragma omp parallel for    
for (ix = 0 ; ix < p->nx ; ix++) {
 for (iy = 0 ; iy < p->ny ; iy++) {
   for (iz = 0 ; iz < p->nz ; iz++) {
       cell = cell_coordinate_to_index(p, ix, iy, iz);
       expmupos_sc[cell] = exp(mupos[iz]);
     }
   }
}

//printf("cpos, expmupos, %.3e %.3e \n", cpos_av[p->nz-1], exp(mupos[p->nz-1]));

while (iterror > maxiterror) {

    Qnegnew = 0.0;
    Qposnew = 0.0; 

 if (p->Nnegions > p->Nposions) {    	

#pragma acc data copyin(Qneg) copyin(Qpos) 
#pragma acc parallel loop present(p[:1]) 
#pragma omp parallel for    
    for (i = 0 ; i < p->n_cells ; i++) {
	p->electric_field[i] = p->rhoF[i];
	p->electric_field[i] += sqrt(p->rhoF[i]*p->rhoF[i] + 4.*p->Nposions*p->Nnegions/Qpos/Qneg*p->exp_born_pos[i]*p->exp_born_neg[i]*expmupos_sc[i]);
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
	p->electric_field[i] += sqrt(p->rhoF[i]*p->rhoF[i] + 4.*p->Nposions*p->Nnegions/Qpos/Qneg*p->exp_born_pos[i]*p->exp_born_neg[i]*expmupos_sc[i]);
	p->electric_field[i] = p->electric_field[i] / (2.0*p->Nposions/Qpos*p->exp_born_pos[i]*expmupos_sc[i]) ;
	p->electric_field[i] = -log(p->electric_field[i]);
     }
 }    

#pragma acc data copy(Qposnew) copy(Qnegnew) 
#pragma acc parallel loop reduction (+:Qposnew) reduction (+:Qnegnew)  
#pragma omp parallel for reduction (+:Qposnew) reduction (+:Qnegnew)
    for (i = 0 ; i < p->n_cells ; i++) {
        Qposnew += exp(-p->electric_field[i])*p->exp_born_pos[i]*expmupos_sc[i];
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
        p->npos_field[i] = exp(-p->electric_field[i])*p->exp_born_pos[i]/Qpos*p->Nposions*expmupos_sc[i] ; 
        p->nneg_field[i] = exp(p->electric_field[i])*p->exp_born_neg[i]/Qneg*p->Nnegions  ; 
     }

/*    for (i = 0 ; i < p->n_cells ; i++) {
    printf("i pos neg born rhoF %d %.3e %.3e %.3e %.3e \n",i, p->npos_field[i], p->nneg_field[i], p->exp_born[i], p->rhoF[i]);
    }
*/




  return(0);
}

