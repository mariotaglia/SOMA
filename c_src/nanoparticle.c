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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "soma_config.h"
#include "rng.h"
#include "phase.h"

//! \file nanpoarticle.h
//! \brief Code related to the nanoparticles

/*! \brief nanpoarticle information */
/* typedef struct Nanoparticle { */
/*   soma_scalar_t x;             /\*!<\brief pointer to x position *\/ */
/*   soma_scalar_t y;             /\*!<\brief pointer to y position *\/ */
/*   soma_scalar_t z;             /\*!<\brief pointer to z position *\/ */
/*   soma_scalar_t radius;  /\*!<\brief Radius of nanoparticle *\/ */
/*   //  RNG_STATE nanoparticle_state;       //!< \brief Struct which contains all RNGs */
/*   //  struct RNG_STATE *set_states;       //!< RNG states of independet sets. NULL if not used. */
/* } Nanoparticle; */

int init_nanoparticle(struct Phase *p)
{
    p->nanoparticles = (Nanoparticle *) malloc(p->n_nanoparticles * sizeof(Nanoparticle));
    if (p->nanoparticles == NULL)
        {
            fprintf(stderr, "Malloc problem %s:%d\n", __FILE__, __LINE__);
            return -4;
        }
    p->nanoparticles[0].x = 0.5;        //test wall
    p->nanoparticles[0].y = 0.5;
    p->nanoparticles[0].z = 0.5;
    p->nanoparticles[0].radius = 0.25;
    p->nanoparticles[0].field = (soma_scalar_t *) malloc(p->n_cells * sizeof(soma_scalar_t));
    if (p->nanoparticles[0].field == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }

    return 0;
}

int calc_np_field_total(struct Phase *p)
{
    memset(p->nanoparticle_field, 0, p->n_cells * sizeof(soma_scalar_t));
    for (uint64_t inp = 0; inp < p->n_nanoparticles; inp++)
        {                       //later something like mpi_allreduce
            calc_my_np_field(p, &p->nanoparticles[0]);
            add_my_np_field(p, &p->nanoparticles[0]);
        }
    return 0;
}

int calc_my_np_field(struct Phase *p, Nanoparticle * np)
{
    box_to_grid(p, np);
    return 0;
}

int box_to_grid(struct Phase *p, Nanoparticle * np)
{                               //map box nanoparticle to discrete mesh
    soma_scalar_t *xfield = calloc(p->nx, sizeof(soma_scalar_t));
    soma_scalar_t dl = p->Lx / p->nx;
    soma_scalar_t xlo = (np->x - np->radius);
    soma_scalar_t xhi = (np->x + np->radius);
    if (2 * np->radius < dl)
        {
            xfield[(int)(np->x / dl)] = 2 * np->radius / dl;
        }
    else
        {
            for (uint64_t x = 0; x < p->nx; x++)
                {
                    if ((soma_scalar_t) x * dl >= xlo && (soma_scalar_t) (x + 1) * dl <= xhi)
                        xfield[x] = 1;
                    if ((soma_scalar_t) x * dl < xlo && (soma_scalar_t) (x + 1) * dl > xlo)
                        xfield[x] = 1 - fmod(xlo, dl) / dl;
                    if ((soma_scalar_t) (x + 1) * dl > xhi && (soma_scalar_t) (x) * dl < xhi)
                        xfield[x] = fmod(xhi, dl) / dl;
                }
        }
    for (uint64_t x = 0; x < p->nx; x++)        //lazy 1d copy
        for (uint64_t y = 0; y < p->ny; y++)
            for (uint64_t z = 0; z < p->nz; z++)
                np->field[x * p->ny * p->nz + y * p->nz + z] = xfield[x];
    return 0;
}

int add_my_np_field(struct Phase *p, Nanoparticle * np)
{
    for (uint64_t x = 0; x < p->nx; x++)
        for (uint64_t y = 0; y < p->ny; y++)
            for (uint64_t z = 0; z < p->nz; z++)
                p->nanoparticle_field[x * p->ny * p->nz + y * p->nz + z] +=
                    np->field[x * p->ny * p->nz + y * p->nz + z];
    return 0;
}
