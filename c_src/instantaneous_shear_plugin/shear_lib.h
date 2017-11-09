/* Copyright (C) 2017 Ludwig Schneider
   Copyright (C) 2017 De-Wen Sun

   This file is part of SOMA.

   SOMA is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   SOMA is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with SOMA.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SHEAR_LIB_H
#define SHEAR_LIB_H

#include "monomer.h"

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
        //! Array, which holds the number of beads for every polymer.
        unsigned int*number_of_beads;
        //Two D array containing the bead information.
        Monomer*beads;
    }bead_data;

//! Get index to acces 2d index to acces matrix of beads
//! \param poly polymer ID
//! \param mono monomer ID
//! \param max_n_beads
unsigned int get_2d_index(const unsigned int poly,const unsigned int mono,
                          const unsigned int max_n_beads);

//! Applies a instantaneous shear to the passed coordinates.
//! \param b Bead data to manipulate
//! \param gradient_direction Direction of the shear gradient (0,1,2)
//! \param flow_direction Direction of the shear flow (0,1,2)
//! \param amplitude shear amplitude
int apply_shear(bead_data*const b,const unsigned int gradient_direction,
                const unsigned int flow_direction,const soma_scalar_t amplitude);

//! Read bead data from a hdf5 file.
//! \param filename Filename to read from
//! \param struct to fill with bead data
int read_bead_data(const char*const filename,bead_data*const b);

//! Write bad data to hdf5 file.
//! \param filename Filename to read from.
//! \param struct with data to write to hdf5.
int write_bead_data(const char*const filename,bead_data*const b);

//! Free dynamic arrays of bead data
//! \param b Bead data to free.
int free_bead_data(bead_data*b);

#endif//SHEAR_LIB_H
