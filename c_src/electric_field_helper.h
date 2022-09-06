/* This file is part of SOMA.                                                                                                                                    
                                                                                                                                                                 
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

// File necessary because helper functions access struct Phase p which is not defined yet when 'electric_field.h' is read.
#include "soma_config.h"
#include "phase.h"

#ifndef SOMA_ELECTRIC_FIELD_HELPER_H
#define SOMA_ELECTRIC_FIELD_HELPER_H

#pragma acc routine(cell_to_index) seq
inline uint64_t cell_to_index(struct Phase *const p, const uint64_t x, const uint64_t y, const uint64_t z)
{
    int64_t xt=x; //changed from uint64_t to int64_t because otherwise xt < 0 always false                                                                                                      
    int64_t yt=y;
    int64_t zt=z;
    if (xt >= p->nx) //Wrap back if necessary // formerly: if (xt >= (uint64) p->nx)                                                                                                            
      xt -= p->nx;
    if (xt < 0)
      xt += p->nx;
    if (yt >= p->ny) //Wrap back if necessary                                                                                                                                                   
      yt -= p->ny;
    if (yt < 0)
      yt += p->ny;
    if (zt >= p->nz) //Wrap back if necessary                                                                                                                                                   
      zt -= p->nz;
    if (zt < 0)
      zt += p->nz;
    //Unified data layout [type][x][y][z]                                                                                                                                                       
    return xt * p->ny * p->nz + yt * p->nz + zt;
}

/*! Helper function to compute the partial derivative of the electric potential field with respect to x                                                                                        
    \private                                                                                                                                                                                   
    \param p Phase describing the system                                                                                                                                                       
    \param x x-coordinate in 3D representation of field.                                                                                                                                       
    \param y y-coordinate in 3D representation of field.                                                                                                                                       
    \param z z-coordinate in 3D representation of field.                                                                                                                                       
    \returns dEpot/dx */
#pragma acc routine(dEpotx) seq
inline soma_scalar_t dEpotx(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y, const uint64_t z)
{
    //Resolve periodic boundaries                                                                                                                                                               
    int64_t xp = x + 1;
    int64_t xm = x - 1;
    if (xp >= p->nx) xp -= p->nx;
    if (xm < 0) xm += p->nx;
    return (e_field[(xp) * p->ny * p->nz + y * p->nz + z] - e_field[(xm) * p->ny * p->nz + y * p->nz + z]) * 0.5 * p->nx / p->Lx;
}

/*! Helper function to compute the partial derivative of the electric potential field with respect to y                                                                                        
    \private                                                                                                                                                                                   
    \param p Phase describing the system                                                                                                                                                       
    \param x x-coordinate in 3D representation of field.                                                                                                                                       
    \param y y-coordinate in 3D representation of field.                                                                                                                                       
    \param z z-coordinate in 3D representation of field.                                                                                                                                       
    \returns dEpot/dy */
#pragma acc routine(dEpoty) seq
inline soma_scalar_t dEpoty(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y, const uint64_t z)
{
    int64_t yp = y + 1;
    int64_t ym = y - 1;
    if (yp >= p->ny) yp -= p->ny;
    if (ym < 0) ym += p->ny;
    return (e_field[x * p->ny * p->nz + (yp) * p->nz + z] - e_field[x * p->ny * p->nz + (ym) * p->nz + z]) * 0.5 * p->ny / p->Ly;
}

/*! Helper function to compute the partial derivative of the electric potential field with respect to z                                                                                        
    \private                                                                                                                                                                                   
    \param p Phase describing the system                                                                                                                                                       
    \param x x-coordinate in 3D representation of field.                                                                                                                                       
    \param y y-coordinate in 3D representation of field.                                                                                                                                       
    \param z z-coordinate in 3D representation of field.                                                                                                                                       
    \returns dEpot/dz */
#pragma acc routine(dEpotz) seq
inline soma_scalar_t dEpotz(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y, const uint64_t z)
{
    int64_t zp = z + 1;
    int64_t zm = z - 1;
    if (zp >= p->nz) zp -= p->nz;
    if (zm < 0) zm += p->nz;
    return (e_field[x * p->ny * p->nz + y * p->nz + (zp)] - e_field[x * p->ny * p->nz + y * p->nz + (zm)]) * 0.5 * p->nz / p->Lz;
}

/*! Helper function to compute the second partial derivative of the electric potential field with respect to x                                                                                 
    \private                                                                                                                                                                                   
    \param p Phase describing the system                                                                                                                                                       
    \param x x-coordinate in 3D representation of field.                                                                                                                                       
    \param y y-coordinate in 3D representation of field.                                                                                                                                       
    \param z z-coordinate in 3D representation of field.                                                                                                                                       
    \returns d2Epot/dx^2 */
#pragma acc routine(d2Epotx) seq
inline soma_scalar_t d2Epotx(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y, const uint64_t z)
{
    int64_t xp = x + 1;
    int64_t xm = x - 1;
    if (xp >= p->nx) xp -= p->nx;
    if (xm < 0) xm += p->nx;
    return (e_field[(xp) * p->ny * p->nz + y * p->nz + z] - (2 * e_field[(x) * p->ny * p->nz + y * p->nz + z]) +
            e_field[(xm) * p->ny * p->nz + y * p->nz + z]) * (p->nx / p->Lx) * (p->nx / p->Lx);
}

/*! Helper function to compute the second partial derivative of the electric potential field with respect to y                                                                                 
    \private                                                                                                                                                                                   
    \param p Phase describing the system                                                                                                                                                       
    \param x x-coordinate in 3D representation of field.                                                                                                                                       
    \param y y-coordinate in 3D representation of field.                                                                                                                                       
    \param z z-coordinate in 3D representation of field.                                                                                                                                       
    \returns d2Epot/dy^2 */
#pragma acc routine(d2Epoty) seq
inline soma_scalar_t d2Epoty(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y, const uint64_t z)
{
    int64_t yp = y + 1;
    int64_t ym = y - 1;
    if (yp >= p->ny) yp -= p->ny;
    if (ym < 0) ym += p->ny;
    return (e_field[x * p->ny * p->nz + (yp) * p->nz + z] - (2 * e_field[x * p->ny * p->nz + (y) * p->nz + z]) +
            e_field[x * p->ny * p->nz + (ym) * p->nz + z]) * (p->ny / p->Ly) * (p->ny / p->Ly);
}

/*! Helper function to compute the second partial derivative of the electric potential field with respect to z                                                                                 
    \private                                                                                                                                                                                   
    \param p Phase describing the system                                                                                                                                                       
    \param x x-coordinate in 3D representation of field.                                                                                                                                       
    \param y y-coordinate in 3D representation of field.                                                                                                                                       
    \param z z-coordinate in 3D representation of field.                                                                                                                                       
    \returns d2Epot/dz^2 */
#pragma acc routine(d2Epotz) seq
inline soma_scalar_t d2Epotz(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y, const uint64_t z)
{
    int64_t zp = z + 1;
    int64_t zm = z - 1;
    if (zp >= p->nz) zp -= p->nz;
    if (zm < 0) zm += p->nz;
    return (e_field[x * p->ny * p->nz + y * p->nz + (zp)] - (2 * e_field[x * p->ny * p->nz + y * p->nz + (z)]) +
            e_field[x * p->ny * p->nz + y * p->nz + (zm)]) * (p->nz / p->Lz) * (p->nz / p->Lz);
}


#endif                          //SOMA_ELECTRIC_FIELD_HELPER_H     
