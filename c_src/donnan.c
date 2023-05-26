
#include "donnan.h"
#include <assert.h>
#include "mpiroutines.h"
#include "soma_config.h"
#include "float.h"

int call_EN(const struct Phase *const p);


int call_EN(const struct Phase *const p)
{
unsigned int i;
soma_scalar_t  Qpos, Qneg; 
soma_scalar_t iterror = DBL_MAX;
soma_scalar_t tmp_npos_field[p->n_cells];
soma_scalar_t tmp_nneg_field[p->n_cells];

const soma_scalar_t  vall = (p->Lx*p->Ly*p->Lz) ;
 
const soma_scalar_t maxiterror = 1e-5 ; // maximum relative iteration error for Qpos and Qneg

Qneg = vall;
Qpos = vall;


#pragma acc data copyin(tmp_npos_field) copyin(tmp_nneg_field)
#pragma acc parallel loop present(p[:1])
#pragma omp parallel for
    for (i = 0 ; i < p->n_cells ; i++) {
        tmp_npos_field[i] = p->Nposions/(p->Lx*p->Ly*p->Lz);
        tmp_nneg_field[i] = p->Nnegions/(p->Lx*p->Ly*p->Lz);
     }


while (iterror > maxiterror) {

#pragma acc data copyin(Qneg) copyin(Qpos)  
#pragma acc parallel loop present(p[:1]) 
#pragma omp parallel for    
    for (i = 0 ; i < p->n_cells ; i++) {
        p->electric_field[i] = 0.5*log(Qneg/Qpos*p->exp_noneq[i]*p->exp_born_pos[i]/p->exp_born_neg[i]);
        p->tempfield_ions[i] = p->alfaions*p->field_scaling_type[0]*(tmp_npos_field[i]+tmp_nneg_field[i])*p->vcell;
    }

    Qnew = 0;
    Qpos = 0;

#pragma acc data copy(Qpos) copy(Qneg) 
#pragma acc parallel loop reduction (+:Qpos) reduction (+:Qneg)  
#pragma omp parallel for reduction (+:Qpos) reduction (+:Qneg)
    for (i = 0 ; i < p->n_cells ; i++) {
        Qpos += exp(-p->electric_field[i])*p->exp_born_pos[i]*p->exp_noneq[i]
&		*exp(-p->alfaions*(p->tempfield_ions[i]+p->tempfield[i]-1.))*vcell;
        Qneg += exp(p->electric_field[i])*p->exp_born_neg[i]
&                *exp(-p->alfaions*(p->tempfield_ions[i]+p->tempfield[i]-1.))*vcell;		;
    }

#pragma acc parallel loop present(p[:1])   
#pragma omp parallel for  
    for (i = 0 ; i < p->n_cells ; i++) {
        p->npos_field[i] = p->Nposions/Qpos*exp(-p->electric_field[i])*p->exp_born_pos[i]*p->exp_noneq[i]
	p->npos_field[i] *= exp(-p->alfaions*(p->tempfield_ions[i]+p->tempfield[i]-1.)) ; 
        p->nneg_field[i] = p->Nnegions/Qneg*exp( p->electric_field[i])*p->exp_born_neg[i]  ; 
	p->nneg_field[i] *= exp(-p->alfaions*(p->tempfield_ions[i]+p->tempfield[i]-1.)) ; 
     }


    iterror = 0;
    for (i = 0 ; i < p->n_cells ; i++) {
        iterror += fabs(p->npos_field[i]-tmp_npos_field[i]);
        iterror += fabs(p->nneg_field[i]-tmp_nneg_field[i]);
    }

        printf("Qpos, Qneg, error, %.3e %.3e %.3e \n", Qpos, Qneg, iterror);
}

  return(0);
}

