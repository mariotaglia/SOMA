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
#include "mesh.h"
int calc_np_field_total(struct Phase *p)
{
    soma_scalar_t *tempfield = calloc(p->nx, sizeof(soma_scalar_t));
    memset(p->nanoparticle_field, 0, p->n_cells * sizeof(soma_scalar_t));
    for (uint64_t inp = 0; inp < p->n_nanoparticle; inp++)
        {
            calc_my_np_field(p, &p->nanoparticles[inp], tempfield);
            add_my_np_field(p, &p->nanoparticles[inp], tempfield);
        }

    //       if (p->umbrella_field)
    //           memcpy(p->umbrella_field, p->nanoparticle_field, p->n_cells * sizeof(soma_scalar_t));   //for debug purposes
    free(tempfield);
#pragma acc update device(p->nanoparticle_field[p->n_cells])
    return 0;
}

int calc_my_np_field(struct Phase *p, Nanoparticle * np, soma_scalar_t * tempfield)
{
    box_to_grid(p, np, tempfield);
    return 0;
}

int box_to_grid(struct Phase *p, Nanoparticle * np, soma_scalar_t * tempfield)
{

    soma_scalar_t dl = p->Lx / p->nx;
    soma_scalar_t xlo = (np->x - np->radius);
    soma_scalar_t xhi = (np->x + np->radius);
    while (xlo < 0)
        xlo += p->Lx;
    while (xhi < 0)
        xhi += p->Lx;
    xlo = fmod(xlo, p->Lx);
    xhi = fmod(xhi, p->Lx);
    int clo = floor(xlo / dl);
    int chi = floor(xhi / dl);
    if (clo == chi)
        tempfield[clo] = 2 * np->radius / dl;
    else
        {
            if (clo < chi)
                for (int i = clo + 1; i < chi; i++)
                    tempfield[i] = 1;
            else
                for (int i = clo + 1; i <= (int)(chi + p->nx); i = (i + 1))
                    tempfield[i % p->nx] = 1;
            tempfield[clo] = ((clo + 1.0) * dl - xlo) / dl;
            tempfield[chi] = (xhi - chi * dl) / dl;
        }
    /* for (uint64_t x = 0; x < p->nx; x++) */
    /*   printf("%lf\n",tempfield[x]); */
    return 0;
}

int add_my_np_field(struct Phase *p, Nanoparticle * np, soma_scalar_t * tempfield)
{
    for (uint64_t x = 0; x < p->nx; x++)
        for (uint64_t y = 0; y < p->ny; y++)
            for (uint64_t z = 0; z < p->nz; z++)
                p->nanoparticle_field[x * p->ny * p->nz + y * p->nz + z] += tempfield[x] * np->interaction;
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
    return 0;
}

int move_nanoparticle(Nanoparticle * np, soma_scalar_t displacement)
{
    np->x += displacement;
    return 0;
}

int resize_nanoparticle(Nanoparticle * np, soma_scalar_t factor)
{
    np->radius *= factor;
    return 0;
}

int test_nanoparticle(struct Phase *p, Nanoparticle * np)
{
    int test = 0;
    soma_scalar_t *tempfield = calloc(p->nx, sizeof(soma_scalar_t));
    for (uint64_t i = 0; i < p->nx * 10; i++)
        {
            memset(tempfield, 0.0, p->nx * sizeof(soma_scalar_t));

            move_nanoparticle(np, -p->Lx / p->nx / 5);

            resize_nanoparticle(np, 1.001);
            calc_my_np_field(p, np, tempfield);

            soma_scalar_t ref_vol = np->radius * 2 * np->interaction;
            soma_scalar_t vol = 0;

            for (uint64_t x = 0; x < p->nx; x++)
                vol += tempfield[x] * p->Lx / p->nx;

            if (vol - ref_vol > 1e-6)
                {
                    test += 1;
                }

        }
    for (uint64_t i = 0; i < p->nx * 10; i++)
        {
            move_nanoparticle(np, p->Lx / p->nx / 5);
            resize_nanoparticle(np, 1.0 / 1.001);
        }
    if (test == 0)
        printf("Nanoparticle test passed.\n");
    else
        printf("Nanoparticle test failed.\n");
    free(tempfield);
    return test;
}
