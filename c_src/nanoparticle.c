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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "soma_config.h"
#include "phase.h"

int init_nanoparticle(struct Phase *p)
{
    p->nanoparticles = (Nanoparticle *) malloc(p->n_nanoparticles * sizeof(Nanoparticle));
    if (p->nanoparticles == NULL)
        {
            fprintf(stderr, "Malloc problem %s:%d\n", __FILE__, __LINE__);
            return -4;
        }
    p->nanoparticles[0].x = 2;  //test wall
    p->nanoparticles[0].y = 0.5;
    p->nanoparticles[0].z = 0.5;
    p->nanoparticles[0].radius = 0.2;
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
{
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
    /* WARNING this is only for the case of 1 nano particle. otherwise addition to existing field */
    memcpy(p->nanoparticle_field, np->field, p->n_cells * sizeof(soma_scalar_t));
    return 0;
}

int nanoparticle_area51_switch(struct Phase *p, int switch_value)
{

    if (p->area51 == NULL)
        {
            /* WARNING: no domain decomposition in this case */
            printf
                ("No area51 found.. creating for initial nanoparticles\n WARNING: This wont work with domain decomposition\n");
            p->area51 = (uint8_t *) malloc(p->n_cells * sizeof(uint8_t));
        }
    if (p->area51 == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
    for (uint64_t cell = 0; cell < p->n_cells; cell++)
        if (p->nanoparticle_field[cell] > 1e-8)
            p->area51[cell] = switch_value;
}
