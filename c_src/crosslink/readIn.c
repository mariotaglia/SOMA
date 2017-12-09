#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>
#include "soma_config.h"
#include "readIn_lib.h"
#include "separateMolecule_lib.h"
#include "rng.h"


int main(int argc, char *argv[])
{
  char *filename = argv[1];
  int choice=strtod (argv[2], NULL);

  if(choice==0){
    char *perkofile=argv[3];
    perkolation(filename,perkofile);
    printf("perkolation!\n");
  }

  if(choice==1){    
    double pzero=strtod (argv[3], NULL);
    char *crosslinkfile = argv[4];
    crosslink(filename,pzero,crosslinkfile);
    printf("crosslinked!\n");
  }

  if(choice==2){
    change_chi(filename);
    printf("changed chi!\n");
  }

  if(choice==4){
    polymer_data polym;
    polymer_data*const poly=&polym;
    bead_data data;
 
    bead_data*const b=&data;
    data.number_of_beads = NULL;
    data.beads = NULL;
    polym.poly_arch=NULL;
    polym.poly_type=NULL;
    polym.poly_type_offset=NULL;
    polym.boxsize=NULL;
   
    read_in(filename,b,poly);
    printf("finished reading! \n");
    independent_set(b,poly);
  }


   return 0;
}
