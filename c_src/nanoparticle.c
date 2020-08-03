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
#include <math.h>
#include "soma_config.h"
#include "phase.h"
#include "mesh.h"
#include "mc.h"

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
    soma_scalar_t d=fmod(np->x,dl)/dl;
    soma_scalar_t dd2=(d-0.5)*(d-0.5);
    soma_scalar_t fnp=1.2;
    
    if (clo < chi)
      for (int i = clo + 1; i < chi; i++)
	tempfield[i] = 1* fnp;
    else
      for (int i = clo + 1; i <= (int)(chi + p->nx); i = (i + 1))
	tempfield[i % p->nx] = 1* fnp;
    soma_scalar_t dlo= ((clo + 1.0) * dl - xlo) / dl;
    soma_scalar_t dhi=(xhi - chi * dl) / dl;



    /* ///const edge */
    /* soma_scalar_t xa= 2.61393661e-01; */
    /* soma_scalar_t xb=-5.07971418e-01; */
    /* soma_scalar_t xc= 4.90827324e-01; */
    /* soma_scalar_t xd=-1.15725804e-01; */
    /* soma_scalar_t xe= 1.55394607e-01; */
    /* soma_scalar_t xf= 7.16087153e-01; */
    /* soma_scalar_t xg= 1.27522649e-05; */


    /* ///const  edge small */
    /* soma_scalar_t xa= 2.79761581e-01; */
    /* soma_scalar_t xb=-6.04518082e-01; */
    /* soma_scalar_t xc= 5.73541725e-01; */
    /* soma_scalar_t xd=-2.18959190e-01; */
    /* soma_scalar_t xe= 1.87293972e-01; */
    /* soma_scalar_t xf= 7.83705283e-01; */
    /* soma_scalar_t xg= 1.69980393e-04; */

    
    /* ///const  bulk small */
    /* soma_scalar_t xa= 2.38969144e-01; */
    /* soma_scalar_t xb=-4.98099364e-01; */
    /* soma_scalar_t xc= 4.35923284e-01; */
    /* soma_scalar_t xd=-1.69777167e-01; */
    /* soma_scalar_t xe= 2.40942369e-01; */
    /* soma_scalar_t xf= 7.54415562e-01; */
    /* soma_scalar_t xg=-6.61919857e-04; */

    /* soma_scalar_t xa=5.58533470e-02; */
    /* soma_scalar_t xb=-7.28120871e-02; */
    /* soma_scalar_t xc=8.38215045e-02; */
    /* soma_scalar_t xd=5.97644130e-02; */
    /* soma_scalar_t xe=6.43003438e-02; */
    /* soma_scalar_t xf=8.10989400e-01; */
    /* soma_scalar_t xg=-3.46970772e-05; */
    
    /* soma_scalar_t xa=-2.33275278e-01; */
    /* soma_scalar_t xb=6.51175036e-01; */
    /* soma_scalar_t xc=-5.26773121e-01; */
    /* soma_scalar_t xd=2.37455615e-01; */
    /* soma_scalar_t xe=6.79208152e-02; */
    /* soma_scalar_t xf=8.04816319e-01; */
    /* soma_scalar_t xg=1.79962597e-04; */

    /* const edge small newest */
    /* soma_scalar_t xa=0.2682037467621129; */
    /* soma_scalar_t xb=-0.5208632422052761; */
    /* soma_scalar_t xc=0.500753988738143; */
    /* soma_scalar_t xd=-0.11691499827239908; */
    /* soma_scalar_t xe=0.1467827089162869; */
    /* soma_scalar_t xf=0.7218849608800436; */
    /* soma_scalar_t xg=3.428351380292647e-05; */

    /* const bulk edge small newest */


    /* soma_scalar_t xa=-0.29619148614586593; */
    /* soma_scalar_t xb=0.8117401630556317; */
    /* soma_scalar_t xc=-0.657459689453385; */
    /* soma_scalar_t xd=0.2617042904010601; */
    /* soma_scalar_t xe=0.08054724396460784; */
    /* soma_scalar_t xf=0.800372671630164; */
    /* soma_scalar_t xg=0.00024082204615293888; */


    //const bulk field scal var
    soma_scalar_t xa=(-0.29619148614586593*3  -0.3226686627648815)/4.0;
    soma_scalar_t xb=(0.8117401630556317*3    +0.8749964661358357)/4.0;
    soma_scalar_t xc=(-0.657459689453385*3    -0.7221984713707235)/4.0;
    soma_scalar_t xd=(0.2617042904010601 *3   +0.33881015200606035)/4.0;
    soma_scalar_t xe=(0.08054724396460784 *3  -0.009165689033214073)/4.0;
    soma_scalar_t xf=(0.800372671630164    *3 +0.8408752811677045)/4.0;
    soma_scalar_t xg=(0.00024082204615293888*3-3.0333153610555586e-05)/4.0;
      




    
    soma_scalar_t flo=xg+xf*pow(dlo,1)+xe*pow(dlo,2)+xd*pow(dlo,3)+xc*pow(dlo,4)+xb*pow(dlo,5)+xa*pow(dlo,6);
    soma_scalar_t fhi=xg+xf*pow(dhi,1)+xe*pow(dhi,2)+xd*pow(dhi,3)+xc*pow(dhi,4)+xb*pow(dhi,5)+xa*pow(dhi,6);
    tempfield[clo] = flo*fnp;
    tempfield[chi] = fhi*fnp;
    printf("%lf\t\%lf\n%lf\t%lf\n",dlo,flo,dhi,fhi);
    /* if(tempfield[chi]<1e-6) */
    /*   tempfield[chi-1]*=np->interaction; */
    
    FILE *f = fopen("wall", "w");
    if (f == NULL)
      {
        printf("Error opening file!\n");
        exit(1);
      }
    for(int i =0;i<p->nx;i++)
      fprintf(f,"%lf\n",tempfield[i]);
    fclose(f);
    
    return 0;
}

int add_my_np_field(struct Phase *p, Nanoparticle * np, soma_scalar_t * tempfield)
{
  for (uint64_t x = 0; x < p->nx; x++)
      for (uint64_t y = 0; y < p->ny; y++)
        for (uint64_t z = 0; z < p->nz; z++){
	      p->nanoparticle_field[x * p->ny * p->nz + y * p->nz + z] += tempfield[x] ;        
        }
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
int init_nanoparticle_rng(struct Phase *p, Nanoparticle *np)
{
  allocate_rng_state(&(np->nanoparticle_rng_state), p->args.pseudo_random_number_generator_arg);
  seed_rng_state(&(np->nanoparticle_rng_state), p->args.rng_seed_arg,
                 0, p->args.pseudo_random_number_generator_arg);
  return 0;
}
int nanoparticle_mc_move(struct Phase *p, Nanoparticle *np){
  int accept=0;
  soma_scalar_t scale=0.01;
  soma_scalar_t delta_energy=0.0;
  soma_scalar_t e0=0.0;
  soma_scalar_t e1=0.0;

  calc_np_field_total(p);
  for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
    p->tempfield[cell] = p->field_scaling_type[0] * p->fields_unified[cell] + p->nanoparticle_field[cell];
  self_omega_field(p);
  calc_non_bonded_energy( p, &e0);

  soma_scalar_t dx = scale * (soma_rng_soma_scalar(&np->nanoparticle_rng_state, p->args.pseudo_random_number_generator_arg) - 0.5);
  np->x+=dx;
  
  calc_np_field_total(p);
  for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
    p->tempfield[cell] = p->field_scaling_type[0] * p->fields_unified[cell] + p->nanoparticle_field[cell];
  self_omega_field(p);
  //  printf("%lf\n",p->omega_field_unified[38 * p->ny * p->nz ]);
  calc_non_bonded_energy( p, &e1);
  delta_energy=e1-e0;
  if(som_accept(&np->nanoparticle_rng_state, p->args.pseudo_random_number_generator_arg, delta_energy)==1)
    {
      accept=1;
    }
  else
    {
      np->x-=dx;
      calc_np_field_total(p);
      for (uint64_t cell = 0; cell < p->n_cells_local; cell++)
        p->tempfield[cell] = p->field_scaling_type[0] * p->fields_unified[cell] + p->nanoparticle_field[cell];
      self_omega_field(p);
      

    }
  printf("%i\t%lf\n",accept,delta_energy);
  return accept;
}
