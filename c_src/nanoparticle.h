/// p->nanoparticle pointer to nanoparticle position (similar to polymer[bead].x.y.z)pointer to nanoparticle position (similar to polymer[bead].x.y.z)o/* Copyright (C) 2016-2019 Ludwig Schneider

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
//#ifndef POLYMER_H
//#define POLYMER_H
#include "soma_config.h"
#include "rng.h"
#include "phase.h"
//! \file nanpoarticle.h
//! \brief Code related to the nanoparticles

/*! \brief nanpoarticle information */
typedef struct Nanoparticle {
    soma_scalar_t x;            /*!<\brief pointer to x position */
    soma_scalar_t y;            /*!<\brief pointer to y position */
    soma_scalar_t z;            /*!<\brief pointer to z position */
    soma_scalar_t radius;       /*!<\brief Radius of nanoparticle */
    soma_scalar_t *field;       /*!<\brief Density field created by THIS nanoparticle */
    //  RNG_STATE nanoparticle_state;       //!< \brief Struct which contains all RNGs
    //  struct RNG_STATE *set_states;       //!< RNG states of independet sets. NULL if not used.
} Nanoparticle;

/*! \brief Init nanoparticle data
  \param p Initialized configuration.
  \return error value
*/
int init_nanoparticle(struct Phase *p);
/*! \brief Calculate the density field contribution of ALL nanoparticles into p->nanoparticle_field
  \param p Initialized configuration.
  \return error value
*/
int calc_np_field(struct Phase *p);
int calc_np_field_total(struct Phase *p);

/*! \brief Calculate the density field of nanoparticle np, saved into its struct
  \param p Initialized configuration.
  \param np nanoparticle.
  \return error value
*/

int calc_my_np_field(struct Phase *p, Nanoparticle * np);
/*! \brief Map the density of a box nanoparticle to the discrete density field
  \param p Initialized configuration.
  \param np nanoparticle.
  \return error value
*/
int add_my_np_field(struct Phase *p, Nanoparticle * np);

int box_to_grid(struct Phase *p, Nanoparticle * np);
