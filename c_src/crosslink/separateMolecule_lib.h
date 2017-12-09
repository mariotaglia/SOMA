#ifndef separateMolecule_lib
#define separateMolecule_lib

//#include "../struct.h"
#include "monomer.h"
#include "soma_util.h"

int separate_molecule(const char* const filename,bead_data*b,polymer_data*poly);


void cross_link_monomer1(int32_t** poly_arch_tmp_pointer,polymer_data*const poly, unsigned int poly_arch_len,unsigned int monomerN, unsigned int monomerM);

int independent_set(bead_data*b,polymer_data*poly);
#endif
