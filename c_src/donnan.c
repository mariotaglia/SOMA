
#include "donnan.h"
#include <assert.h>
#include "mpiroutines.h"
#include "soma_config.h"

int call_donnan(const struct Phase *const p);


int call_donnan(const struct Phase *const p)
{
unsigned int i, type;
soma_scalar_t  rhoQ[p->n_cells_local]; // total charge density

for (i = 0 ; i < p->n_cells_local ; i++) {
        rhoQ[i] = 0.0 ; 
        for (type = 0 ; type < p->n_types; type++) {
                   rhoQ[i] += p->fields_unified[i+p->n_cells_local*type]*p->charges[type];
        } 

   	p->electric_field[i] = rhoQ[i] + sqrt(rhoQ[i]*rhoQ[i] + 4.*p->Nions*p->Nions) ;
	p->electric_field[i] = p->electric_field[i] / (2.0*p->Nions) ;
	p->electric_field[i] = log(p->electric_field[i]);
}

  return(0);
}
