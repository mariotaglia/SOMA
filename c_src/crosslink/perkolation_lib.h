#ifndef perkolation_lib
#define perkolation_lib

//#include "../struct.h"
#include "monomer.h"
#include "soma_util.h"

int check_perkolation_monomer(const char* const filename,bead_data*b,polymer_data*poly, const int* const feed,  const int poly_start);

//int check_perkulation(const char* const filename,bead_data*b,polymer_data*poly,const int** const pbond, const int* const feed);


//Count the number of perkolation
int check_perkolation_monomer1(const char* const filename,bead_data*b,polymer_data*poly, const int poly_start);

void perkolation(const char* const filename,const char* const perkofile);
int check_configfile(const char* const filename);
#endif
