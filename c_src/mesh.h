/* Copyright (C) 2016-2021 Ludwig Schneider
   Copyright (C) 2016 Ulrich Welling
   Copyright (C) 2016 Marcel Langenberg
   Copyright (C) 2016 Fabien Leonforte
   Copyright (C) 2016 Juan Orozco
   Copyright (C) 2016 Yongzhi Ren

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
#ifndef SOMA_MESH_H
#define SOMA_MESH_H
/*! \file mesh.h
  \brief Functions related to the mesh of the density fields.
*/

/*! calculate the 3D cell index from 3 spatial coordinates
  \param p Phase configuration.
  \param rx spatial position x.
  \param ry spatial position y.
  \param rz spatial position z.
  \param x pointer to resulting x index
  \param y pointer to resulting y index
  \param z pointer to resulting z index
*/

#include "soma_config.h"
#include "phase.h"

#pragma acc routine(coord_to_cell_coordinate) seq
static inline void coord_to_cell_coordinate(const struct Phase *p, const soma_scalar_t rx, const soma_scalar_t ry,
                                            const soma_scalar_t rz, int *x, int *y, int *z);
inline void coord_to_cell_coordinate(const struct Phase *p, const soma_scalar_t rx, const soma_scalar_t ry,
                                     const soma_scalar_t rz, int *x, int *y, int *z)
{

    soma_scalar_t px, py, pz;

    // Fold coordinate back into the box
    px = rx - p->Lx * (int)(rx * p->iLx);
    py = ry - p->Ly * (int)(ry * p->iLy);
    pz = rz - p->Lz * (int)(rz * p->iLz);

    // Assure correct symmetry at coordinate = 0
    if (px < 0)
        px = p->Lx + px;
    if (py < 0)
        py = p->Ly + py;
    if (pz < 0)
        pz = p->Lz + pz;

    // Calculate index
    *x = (int)(px * p->iLx * p->nx);
    *y = (int)(py * p->iLy * p->ny);
    *z = (int)(pz * p->iLz * p->nz);
}

/*! calculate field array index from 3D cell coordinates

  Unified data layout [type][x][y][z]
  \param p Phase configuration
  \param x cell index x
  \param y cell index y
  \param z cell index z
  \return memory index to access a field p.fields[type][index]
*/
static inline uint64_t cell_coordinate_to_index(const struct Phase *p, const int x, const int y, const int z);
#pragma acc routine(cell_coordinate_to_index) seq
inline uint64_t cell_coordinate_to_index(const struct Phase *p, const int x, const int y, const int z)
{
    int xt = x;
#if ( ENABLE_DOMAIN_DECOMPOSITION == 1 )
    //For performance reasons compile the domain decomposition only if necessary
    if (p->args.N_domains_arg > 1)
        {
            if (xt >= p->local_nx_high) //Wrap back if necessary
                xt -= p->nx;
            if (xt < p->local_nx_low)
                xt += p->nx;
            if (xt < p->local_nx_low || xt >= p->local_nx_high)
                {
                    return UINT64_MAX;  //Error, requested index out of local bounds
                }
            xt -= p->local_nx_low;
        }
#endif                          //ENABLE_MIC
    //Unified data layout [type][x][y][z]
    return xt * p->ny * p->nz + y * p->nz + z;
}

/*! \brief Calculate an index for any given position and any density field of p.
  This function should ALWAYS be called to get an index of a field.
  Issues:

 - What happens, if the position is out of the box?
 - How can an error be signaled?
 - For performance issues, we could think of a macro instead.

 \param p Initialized Phase struct.
 \param rx X-coordinate of the position in question.
 \param ry Y-coordinate of the position in question.
 \param rz Z-coordinate of the position in question.
 \return Index for the field index referencing.
*/
#pragma acc routine(coord_to_index) seq
static inline uint64_t coord_to_index(const struct Phase *const p, const soma_scalar_t rx, const soma_scalar_t ry,
                                      const soma_scalar_t rz);
inline uint64_t coord_to_index(const struct Phase *p, const soma_scalar_t rx, const soma_scalar_t ry,
                               const soma_scalar_t rz)
{
    int x, y, z;
    coord_to_cell_coordinate(p, rx, ry, rz, &x, &y, &z);
    return cell_coordinate_to_index(p, x, y, z);
}

//! Get the unified index from the cell_index and the type
//!
//! Unified data layout [type][x][y][z]
//! \param p Phase of system (p->n_cells needed)
//! \param cell cell index
//! \param rtype type you request.
//! \return Calculated unified index.
#pragma acc routine(cell_to_index_unified) seq
static inline uint64_t cell_to_index_unified(const struct Phase *const p, const uint64_t cell,
                                             const unsigned int rtype);
inline uint64_t cell_to_index_unified(const struct Phase *const p, const uint64_t cell, const unsigned int rtype)
{
    if (cell == UINT64_MAX)     //Error cell out of bounds
        return UINT64_MAX;
    //Unified data layout [type][x][y][z]
    return cell + rtype * p->n_cells_local;
}

/*! \brief Calculate an index for any given position and type and the density field / external field of p.
  This function should ALWAYS be called to get an index of a field.
  Issues:

 - What happens, if the position is out of the box?
 - How can an error be signaled?
 - For performance issues, we could think of a macro instead.

 \param p Initialized Phase struct.
 \param rx X-coordinate of the position in question.
 \param ry Y-coordinate of the position in question.
 \param rz Z-coordinate of the position in question.
 \param rtype Type of the accessed particle.
 \return Index for the field index referencing.
*/
#pragma acc routine(coord_to_index_unified) seq
static inline uint64_t coord_to_index_unified(const struct Phase *const p, const soma_scalar_t rx,
                                              const soma_scalar_t ry, const soma_scalar_t rz, unsigned int rtype);
inline uint64_t coord_to_index_unified(const struct Phase *p, const soma_scalar_t rx, const soma_scalar_t ry,
                                       const soma_scalar_t rz, const unsigned int rtype)
{
    int x, y, z;

    coord_to_cell_coordinate(p, rx, ry, rz, &x, &y, &z);
    const uint64_t cell = cell_coordinate_to_index(p, x, y, z);
    return cell_to_index_unified(p, cell, rtype);
}

/*! Update the density fields \f$ \rho \f$ of the system.
 The field update fist sets all entries to zero. Then performs a loop over all polymers eg. monomers.
 The field itself is chosen according to the monomer type. Be shure that there is no additional type set otherwise
 the adressing fails.
 \param p Phase configuration to update the  densityfields.
 \return Errorcode
 \note this call requires MPI collectives so call it on every MPI-core.
 */
int update_density_fields(const struct Phase *const p);

/*! Update the omega \f$ \omega\f$-fields of the system.

  The function requires updated density fields, as a consequence
  update_density_fields is called. So the MPI collective
  characteristics are inherited.
  \note In accordance to the hamiltonian specified via the cmd line interface,
  the correct hamiltonian is selected.

 \param p Phase configuration to update the omega-fields.
 \note this call requires MPI collectives so call it on every MPI-core.
*/
void update_omega_fields(const struct Phase *const p);

//! Function to calculate the omega fields for hamiltonian scmf0.
//!
//! \param p Phase configuration to update.
//! \f[\frac{H}{k_BT } = \frac{\rho_0}{N_{ref}}\int_V D(\{r\}) \sum_i^{n_{types}}\frac{k_i}{2} \left(\phi_i(r)-s_i(r)\right)^2 + \sum_i^{n_{types}} \phi_i(r) f_i(r) + \frac{\kappa N_{ref}}{2} (\sum_{i} ^{n_{types}} \phi_i(r) - 1)^2 - \sum_{i<j}^{n_{types}} \frac{\chi_{i,j}N_{ref}}{4} (\phi_i(r) - \phi_j(r))^2 \f]
//! \f[  \frac{\omega_k(c)}{k_BT} = \frac{N}{k_BT \rho_0 \Delta L^3}\frac{\partial H(c)}{k_BT \partial \phi_k} = \sum_i^{n_{types}} k_i\left(\phi_i(c)-s_i(c)\right) + f_k(c) + \kappa (\sum_i^{n_{types}} \phi_i(c)-1) - \sum_{i< k}^{n_{types}} \frac{\chi_{i,k}}{2} (\phi_k(c)-\phi_i(c))   \f]
//! \note In case of the use of the density_weights every \f$\phi_i\f$ is multiplied with its weight \f$w_i\f$. For this Hamiltonian weights differing from 1 are not recommended -- use SCMF1 instead.
void update_omega_fields_scmf0(const struct Phase *const p);

//! Function to calculate the omega fields for hamiltonian scmf1.
//!
//! \param p Phase configuration to update.
//! \f[\frac{H}{k_BT } = \frac{\rho_0}{N_{ref}}\int_V D(\{r\}) \sum_i^{n_{types}} \frac{k_i}{2} \left(\phi_i(r)-s_i(r)\right)^2+ \sum_i^{n_{types}} \phi_i(r) f_i(r) + \frac{\kappa N_{ref}}{2} \left(\sum_{i} ^{n_{types}} \phi_i(r) - 1\right)^2 + \sum_{i< j}^{n_{types}} \chi_{i,j}N_{ref} \phi_i(r) \phi_j(r) \f]
//! \f[ \omega_k(c) = \frac{N_{ref}}{k_BT\rho_0 \Delta L^3}\frac{\partial H(c)}{\partial \phi_k} = \sum_i^{n_{types}} k_i\left(\phi_i(c)-s_i(c)\right)+ f_k(c) + \kappa \left(\sum_i^{n_{types}} \phi_i(c)-1\right) + \sum_{i< k}^{n_{types}} \chi_{i,k} \phi_i(c)   \f]
//! \note In case of the use of the density_weights every \f$\phi_i\f$ is multiplied with its weight \f$w_i\f$.
void update_omega_fields_scmf1(const struct Phase *const p);

#endif                          //SOMA_MESH_H
