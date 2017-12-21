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

#include <stdio.h>
#include <stdlib.h>
#include <hdf5.h>

#include "soma_config.h"
#include "shear_lib.h"

/*! file shear.c
  \brief Helper program to apply instantaneous shear to a configuration.

  The program reads the "/beads" dataset of a hdf5 file modifies the
  position according to a specified instantaneous shear and writes
  them back in the configuration.
*/

int main(int argc, char *argv[])
    {
    printf("Component of SOMA Copyright (C) 2016-2017 Ludwig Schneider, Ulrich\n"
           "Welling, Marcel Langenberg, Fabien Leonforte, Juan Orozco, Yongzhi\n"
	   "Ren, De-Wen Sun. This program comes with ABSOLUTELY NO WARRANTY; see\n"
	   "GNU Lesser General Public License v3 for details. This is free\n"
	   "software, and you are welcome to redistribute it under certain\n"
	   "conditions; see GNU Lesser General Public License v3 for details. \n\n"
	   "You are using SOMA please cite: \n * Daoulas, Kostas Ch. and Müller, "
	   "Marcus , J. Chem.Phys.2006, 125,18\n");

    if(argc != 5)
	{
	fprintf(stderr,"Usage: filename.h5 shear_gradient_direction(x=0,y=1,z=2) shear_flow_direction(x=0,y=1,z=2) shear_amplitude\n");
	return -1;
	}
    const char*const filename = argv[1];
    const unsigned int gradient_dir = atoi(argv[2]);
    const unsigned int flow_dir = atoi(argv[3]);
    const soma_scalar_t amplitude = atof(argv[4]);

    bead_data data;
    bead_data*const b=&data;
    data.number_of_beads = NULL;
    data.beads = NULL;

    int return_value=0;

    if( return_value == 0)
	{
	const int read = read_bead_data(filename, b);
	if( read != 0)
	    {
	    fprintf(stderr,"ERROR: Unsuccessful read of data.\n");
	    return_value= -2;
	    }
	}

    if( return_value == 0)
	{
	const int shear = apply_shear(b, gradient_dir, flow_dir, amplitude);
	if( shear != 0)
	    {
	    fprintf(stderr,"ERROR: Unsuccessful attempt to shear beads.\n");
	    return_value= -3;
	    }
	}

    if( return_value == 0)
	{
	const int write = write_bead_data(filename, b);
	if( write != 0)
	    {
	    fprintf(stderr,"ERROR: Unsuccessful attempt to write data to file.\n");
	    return_value = -4;
	    }
	}

    free_bead_data(b);
    return return_value;
    }
