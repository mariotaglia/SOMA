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
#include "soma_config.h"
//#include "rng.h"
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
/*! \brief Calc p->nanoparticle_field containing all nanoparticle densities
  \param p Initialized configuration.
  \return error value
*/
int calc_np_field_total(struct Phase *p);
/*! \brief Calculate the density field of nanoparticle np, saved into its struct
  \param p Initialized configuration.
  \param np nanoparticle.
  \return error value
*/
int calc_my_np_field(struct Phase *p, Nanoparticle * np);
/*! \brief Add density field of nanoparticle np to p->nanoparticle_field
  \param p Initialized configuration.
  \param np Nanoparticle.
  \return error value
*/
int add_my_np_field(struct Phase *p, Nanoparticle * np);
/*! \brief Map density field of box nanoparticle to discrete grid
  \param p Initialized configuration.
  \param np Nanoparticle.
  \return error value
*/
int box_to_grid(struct Phase *p, Nanoparticle * np);
/*! \brief switch area51 at nps location to given value
  \param p Initialized configuration.
  \param switch_value value
  \return error value
*/
int nanoparticle_area51_switch(struct Phase *p, int switch_value);
