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

/*! Helper function to compute cell description (x,y,z) to index of the 1D arrays
    \private
    \param p Phase describing the system
    \param x x-coordinate in 3D representation of field.
    \param y y-coordinate in 3D representation of field.
    \param z z-coordinate in 3D representation of field.
    \returns index */
#pragma acc routine(cell_to_index) seq
static inline uint64_t cell_to_index(struct Phase *const p, const int64_t x, const int64_t y, const int64_t z);
inline uint64_t cell_to_index(struct Phase *const p, const int64_t x, const int64_t y, const int64_t z)
{
    int64_t xt = x;
    int64_t yt = y;
    int64_t zt = z;
    // Wrap back if necessary; add padding in case non-periodic boundaries
    if (p->ef.el_pos_yz)
    {
        if (xt >= p->nx) xt = p->nx-1;
        if (xt < 0) xt = 0;
    }
    else
    {
        if (xt >= p->nx) xt -= p->nx;
        if (xt < 0) xt += p->nx;
    }
    if (p->ef.el_pos_xz)
    {
        if (yt >= p->ny) yt = p->ny-1;
        if (yt < 0) yt = 0;
    }
    else
    {
        if (yt >= p->ny) yt -= p->ny;
        if (yt < 0) yt += p->ny;   
    }
    if (p->ef.el_pos_xy)
    {
        if (zt >= p->nz) zt = p->nz-1;
        if (zt < 0) zt = 0;   
    }   
    else
    {
        if (zt >= p->nz) zt -= p->nz;
        if (zt < 0) zt += p->nz;
    }
    return xt * p->ny * p->nz + yt * p->nz + zt;
}

/*! Helper function to compute cell description (x,y,z) to index of the convoluted 1D arrays
 *  Will automatically add electrode planes by considering offsets
    \private
    \param p Phase describing the system
    \param x x-coordinate in 3D representation of field.
    \param y y-coordinate in 3D representation of field.
    \param z z-coordinate in 3D representation of field.
    \returns index */
#pragma acc routine(cell_to_index_conv) seq
static inline uint64_t cell_to_index_conv(struct Phase *const p, const int64_t x, const int64_t y, const int64_t z);
inline uint64_t cell_to_index_conv(struct Phase *const p, const int64_t x, const int64_t y, const int64_t z)
{
    int64_t xt = x;
    int64_t yt = y;
    int64_t zt = z;
    // Wrap back if necessary; consider offset for electrode planes & non-periodic boundaries
    // in order to not wrap back on opposing electrode
    if (p->ef.el_pos_yz)
    {
        if (xt >= p->ef.conv_nx + p->ef.x_offset * 2) xt = (p->ef.conv_nx + p->ef.x_offset * 2) - 1;
        if (xt < 0) xt = 0;
    }
    else
    {
        if (xt >= p->ef.conv_nx + p->ef.x_offset * 2) xt -= p->ef.conv_nx + p->ef.x_offset * 2;
        if (xt < 0) xt += p->ef.conv_nx;
    }
    if (p->ef.el_pos_xz)
    {
        if (yt >= p->ef.conv_ny + p->ef.y_offset * 2) yt = (p->ef.conv_ny + p->ef.y_offset * 2) - 1;
        if (yt < 0) yt = 0;
    }
    else
    {
        if (yt >= p->ef.conv_ny + p->ef.y_offset * 2) yt -= p->ef.conv_ny + p->ef.y_offset * 2;
        if (yt < 0) yt += p->ef.conv_ny;
    }
    if (p->ef.el_pos_xy)
    {
        if (zt >= p->ef.conv_nz + p->ef.z_offset * 2) zt = (p->ef.conv_nz + p->ef.z_offset * 2) - 1;
        if (zt < 0) zt = 0;
    }
    else
    {
        if (zt >= p->ef.conv_nz + p->ef.z_offset * 2) zt -= p->ef.conv_nz + p->ef.z_offset * 2;
        if (zt < 0) zt += p->ef.conv_nz;
    }
    // addition of "y_offset * 2, etc." necessary to add potential planes of electrodes
    return xt * (p->ef.conv_ny + p->ef.y_offset * 2) * (p->ef.conv_nz + p->ef.z_offset * 2) + yt * (p->ef.conv_nz + p->ef.z_offset * 2) + zt;
}

/*! Helper function to compute cell description (x,y,z) to index of the kernel
    \private
    \param p Phase describing the system
    \param x x-coordinate in 3D representation of field.
    \param y y-coordinate in 3D representation of field.
    \param z z-coordinate in 3D representation of field.
    \returns index */
#pragma acc routine(cell_to_index_kernel) seq
static inline uint16_t cell_to_index_kernel(struct Phase *const p, const int16_t x, const int16_t y, const int16_t z);
inline uint16_t cell_to_index_kernel(struct Phase *const p, const int16_t x, const int16_t y, const int16_t z)
{
    int16_t ki_x = x;
    int16_t ki_y = y;
    int16_t ki_z = z;
    ki_x += p->ef.kernel_rad;
    ki_y += p->ef.kernel_rad;
    ki_z += p->ef.kernel_rad;
    //Unified data layout [type][x][y][z]
    return ki_x * p->ef.kernel_dim * p->ef.kernel_dim + ki_y * p->ef.kernel_dim + ki_z;
}

/*! Helper function to compute the partial derivative of a field with respect to x
    \private
    \param p Phase describing the system
    \param x x-coordinate in 3D representation of field.
    \param y y-coordinate in 3D representation of field.
    \param z z-coordinate in 3D representation of field.
    \param nx_cells maximum amount of cells along x axis.
    \param ny_cells maximum amount of cells along y axis.
    \param nz_cells maximum amount of cells along z axis.
    \returns dEpot/dx */
#pragma acc routine(dEpotx) seq
static inline soma_scalar_t dEpotx(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y, 
                                   const uint64_t z, const int64_t nx_cells, const int64_t ny_cells, const int64_t nz_cells);
inline soma_scalar_t dEpotx(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y, 
                            const uint64_t z, const int64_t nx_cells, const int64_t ny_cells, const int64_t nz_cells)
{
    //Resolve periodic boundaries
    int64_t xp = x + 1;
    int64_t xm = x - 1;
    if (xp >= nx_cells) xp -= nx_cells;
    if (xm < 0) xm += nx_cells;
    return (e_field[(xp) * ny_cells * nz_cells + y * nz_cells + z] - e_field[(xm) * ny_cells * nz_cells + y * nz_cells + z]) * 0.5 * p->nx / p->Lx;
}

/*! Helper function to compute the partial derivative of a field with respect to y
    \private
    \param p Phase describing the system
    \param x x-coordinate in 3D representation of field.
    \param y y-coordinate in 3D representation of field.
    \param z z-coordinate in 3D representation of field.
    \param ny_cells maximum amount of cells along y axis.
    \param nz_cells maximum amount of cells along z axis.
    \returns dEpot/dy */
#pragma acc routine(dEpoty) seq
static inline soma_scalar_t dEpoty(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y,
                                   const uint64_t z, const int64_t ny_cells, const int64_t nz_cells);
inline soma_scalar_t dEpoty(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y,
                            const uint64_t z, const int64_t ny_cells, const int64_t nz_cells)
{
    int64_t yp = y + 1;
    int64_t ym = y - 1;
    if (yp >= ny_cells) yp -= ny_cells;
    if (ym < 0) ym += ny_cells;
    return (e_field[x * ny_cells * nz_cells + (yp) * nz_cells + z] - e_field[x * ny_cells * nz_cells + (ym) * nz_cells + z]) * 0.5 * p->ny  / p->Ly;
}

/*! Helper function to compute the partial derivative of a field with respect to z
    \private
    \param p Phase describing the system
    \param x x-coordinate in 3D representation of field.
    \param y y-coordinate in 3D representation of field.
    \param z z-coordinate in 3D representation of field.
    \param ny_cells maximum amount of cells along y axis.
    \param nz_cells maximum amount of cells along z axis.
    \returns dEpot/dz */
#pragma acc routine(dEpotz) seq
static inline soma_scalar_t dEpotz(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y,
                                   const uint64_t z, const int64_t ny_cells, const int64_t nz_cells);
inline soma_scalar_t dEpotz(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y, 
                            const uint64_t z, const int64_t ny_cells, const int64_t nz_cells)
{
    int64_t zp = z + 1;
    int64_t zm = z - 1;
    if (zp >= nz_cells) zp -= nz_cells;
    if (zm < 0) zm += nz_cells;
    return (e_field[x * ny_cells * nz_cells + y * nz_cells + (zp)] - e_field[x * ny_cells * nz_cells + y * nz_cells + (zm)]) * 0.5 * p->nz / p->Lz;
}

/*! Helper function to compute the second partial derivative of a field with respect to x
    \private
    \param p Phase describing the system
    \param x x-coordinate in 3D representation of field.
    \param y y-coordinate in 3D representation of field.
    \param z z-coordinate in 3D representation of field.
    \param nx_cells maximum amount of cells along x axis.
    \param ny_cells maximum amount of cells along y axis.
    \param nz_cells maximum amount of cells along z axis.
    \returns d2Epot/dx^2 */
#pragma acc routine(d2Epotx) seq
static inline soma_scalar_t d2Epotx(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y, 
                                    const uint64_t z, const int64_t nx_cells, const int64_t ny_cells, const int64_t nz_cells);
inline soma_scalar_t d2Epotx(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y,
                             const uint64_t z, const int64_t nx_cells, const int64_t ny_cells, const int64_t nz_cells)
{
    int64_t xp = x + 1;
    int64_t xm = x - 1;
    if (xp >= nx_cells) xp -= nx_cells;
    if (xm < 0) xm += nx_cells;
    return (e_field[(xp) * ny_cells * nz_cells + y * nz_cells + z] - (2 * e_field[(x) * ny_cells * nz_cells + y * nz_cells + z]) +
            e_field[(xm) * ny_cells * nz_cells + y * nz_cells + z]) * (p->nx / p->Lx) * (p->nx / p->Lx);
}

/*! Helper function to compute the second partial derivative of a field with respect to y
    \private
    \param p Phase describing the system
    \param x x-coordinate in 3D representation of field.
    \param y y-coordinate in 3D representation of field.
    \param z z-coordinate in 3D representation of field.
    \param ny_cells maximum amount of cells along y axis.
    \param nz_cells maximum amount of cells along z axis.
    \returns d2Epot/dy^2 */
#pragma acc routine(d2Epoty) seq
static inline soma_scalar_t d2Epoty(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y,
                                    const uint64_t z, const int64_t ny_cells, const int64_t nz_cells);
inline soma_scalar_t d2Epoty(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y,
                             const uint64_t z, const int64_t ny_cells, const int64_t nz_cells)
{
    int64_t yp = y + 1;
    int64_t ym = y - 1;
    if (yp >= ny_cells) yp -= ny_cells;
    if (ym < 0) ym += ny_cells;
    return (e_field[x * ny_cells * nz_cells + (yp) * nz_cells + z] - (2 * e_field[x * ny_cells * nz_cells + (y) * nz_cells + z]) +
            e_field[x * ny_cells * nz_cells + (ym) * nz_cells + z]) * (p->ny / p->Ly) * (p->ny / p->Ly);
}

/*! Helper function to compute the second partial derivative of a field with respect to z
    \private
    \param p Phase describing the system
    \param x x-coordinate in 3D representation of field.
    \param y y-coordinate in 3D representation of field.
    \param z z-coordinate in 3D representation of field.
    \param ny_cells maximum amount of cells along y axis.
    \param nz_cells maximum amount of cells along z axis.
    \returns d2Epot/dz^2 */
#pragma acc routine(d2Epotz) seq
static inline soma_scalar_t d2Epotz(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y,
                                    const uint64_t z, const int64_t ny_cells, const int64_t nz_cells);
inline soma_scalar_t d2Epotz(struct Phase *const p, soma_scalar_t * e_field, const uint64_t x, const uint64_t y,
                             const uint64_t z, const int64_t ny_cells, const int64_t nz_cells)
{
    int64_t zp = z + 1;
    int64_t zm = z - 1;
    if (zp >= nz_cells) zp -= nz_cells;
    if (zm < 0) zm += nz_cells;
    return (e_field[x * ny_cells * nz_cells + y * nz_cells + (zp)] - (2 * e_field[x * ny_cells * nz_cells + y * nz_cells + (z)]) +
            e_field[x * ny_cells * nz_cells + y * nz_cells + (zm)]) * (p->nz / p->Lz) * (p->nz / p->Lz);
}


#endif                          //SOMA_ELECTRIC_FIELD_HELPER_H     
