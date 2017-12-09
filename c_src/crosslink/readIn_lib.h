#ifndef readIn_lib
#define readIn_lib

//#include "../struct.h"
#include "monomer.h"
#include "soma_util.h"


/*! file shear_lib.h
  \brief File declaring functions to apply instantaneous shear to coordinates.
*/

//! Struct combining relevant parameter to store bead data.
typedef struct
{
  //! Max number of beads per polymer
  unsigned int max_n_beads;
  //! Number of polymers
  uint64_t N_polymers;
    
  uint64_t old_N_polymers;
  //! Array, which holds the number of beads for every polymer.
  unsigned int*number_of_beads;
  //Two D array containing the bead information.
  Monomer * beads;
}bead_data;

typedef struct
{
  float*boxsize;
  
  uint32_t*poly_arch;

  unsigned int*poly_type;
  
  int*poly_type_offset;

  unsigned int n_poly_type;

  unsigned int poly_arch_length;
}polymer_data;


//! Read bead data from a hdf5 file.
//! \param filename Filename to read from
//! \param struct to fill with bead data
void read_in(const char*const filename,bead_data*const b,polymer_data*const poly);
int write_out(const char*const filename, bead_data*const b);

//! Free dynamic arrays of bead data
//! \param b Bead data to free.
void free_data(bead_data*b, polymer_data*poly);

int change_chi(const char* const filename);

// Convert filename(h5) into outputname(vtf)
void convert(bead_data*const b, polymer_data*const poly, const char*const outputname);
// Apply periodic bondary condition
Monomer * feed_back_box(/*const char*const filename,*/bead_data*b,polymer_data*poly);

//Connect all polymers to one molecule
int cross_link_new(const char*const filename, bead_data * const b, polymer_data * const poly, unsigned int start, unsigned int end);

//Connect monomers
//Update the hdf5 file, can be used independently
//poly->poly_arch not updated!!! Need to be reloaded in the next function.
int cross_link_monomer(const char*const FILE, bead_data * const b, polymer_data * const poly,unsigned int polyN,unsigned int monomerN,unsigned int polyM, unsigned int monomerM);

//Connect monomers
//Does not update hdf5 file, to be used in other functions(Prob_bonds)
int cross_link_monomer2(bead_data * const b, polymer_data * const poly,unsigned int polyN,unsigned int monomerN,unsigned int monoM, unsigned int *poly_arch_len_pointer,int32_t ** poly_arch_tmp_pointer);

//Crosslink the network according to their distance
int Prob_bonds(const char* const filename,bead_data*b,polymer_data*poly,const double pzero,Monomer*new_position,const char* const crosslinkfile);


//Crosslink the network and return the number of perkulation
int crosslink1(const char* const filename,const double pzero,const char* const crosslinkfile);

//Crosslink the network and return zero
void crosslink(const char* const filename,const double pzero,const char* const crosslinkfile);

#endif
