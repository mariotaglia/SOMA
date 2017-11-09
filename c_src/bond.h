/* Copyright (C) 2016-2017 Ludwig Schneider

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
#ifndef BOND_H
#define BOND_H

/*!\file bond.h
  \brief Definition of Bond related code pieces.
*/

/*! \brief Bond type enumerator to indicate the different bond
 *  types. Matches the bondDict in the ConfGen.py.*/
enum Bondtype{
    HARMONIC=0, /*!<\brief Harmonic bond with a single spring const. */
/*!\deprecated stiff is not implemented, at least for now */
    STIFF=1,    /*!<\brief Stiff bonds.*/
    HARMONICVARIABLESCALE=2 /*!<\brief Harmonic bond with a single spring const. But with an additional scaling factor for influence of Photoswitches.*/
    };

//! Number of implemented bonds, definition in init.c
extern const unsigned int NUMBER_SOMA_BOND_TYPES;


#endif//BOND_H
