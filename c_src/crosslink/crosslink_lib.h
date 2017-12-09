#ifndef crosslink_lib
#define crosslink_lib

//#include "../struct.h"
#include "../monomer.h"
#include "../soma_util.h"

//Connect all polymers to one molecule
int cross_link_new(const char*const filename, bead_data * const b, polymer_data * const poly, unsigned int start, unsigned int end);

//Connect monomers
//poly->poly_arch not updated!!! Need to be reloaded in the next function.
int cross_link_monomer(const char*const FILE, bead_data * const b, polymer_data * const poly,unsigned int polyN,unsigned int monomerN,unsigned int polyM, unsigned int monomerM);
//Connect monomers according to their distance
int  Prob_bonds(const char* const filename,bead_data*b,polymer_data*poly);

int separate_molecule(bead_data*b,polymer_data*poly);

#endif
