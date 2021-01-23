/* Copyright (C) 2016-2021 Ludwig Schneider

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

//! \file monomer.c
//! \brief Implementation of monomer.h

#include "monomer.h"

hid_t get_monomer_memtype(void)
{
    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(Monomer));
    H5Tinsert(memtype, "x", HOFFSET(Monomer, x), H5T_SOMA_NATIVE_SCALAR);
    H5Tinsert(memtype, "y", HOFFSET(Monomer, y), H5T_SOMA_NATIVE_SCALAR);
    H5Tinsert(memtype, "z", HOFFSET(Monomer, z), H5T_SOMA_NATIVE_SCALAR);

    return memtype;
}

hid_t get_monomer_filetype(void)
{
    hid_t memtype = H5Tcreate(H5T_COMPOUND, 4 * 8);
    H5Tinsert(memtype, "x", 0, H5T_SOMA_FILE_SCALAR);
    H5Tinsert(memtype, "y", 8, H5T_SOMA_FILE_SCALAR);
    H5Tinsert(memtype, "z", 16, H5T_SOMA_FILE_SCALAR);
    return memtype;
}
