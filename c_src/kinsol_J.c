/* Routines for solving a set of non-linear coupled equations
 * Finds x that satisfy F(x) = 0
 * Based on kinsol  6.5 example "inKrylovDemo_ls.c"  */

#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <kinsol/kinsol.h>               /* access to KINSOL func., consts.      */
#include <sunlinsol/sunlinsol_spgmr.h>   /* access to SPGMR SUNLinearSolver      */
#include <sunlinsol/sunlinsol_spbcgs.h>  /* access to SPBCGS SUNLinearSolver     */
#include <sunlinsol/sunlinsol_sptfqmr.h> /* access to SPTFQMR SUNLinearSolver    */
#include <sunlinsol/sunlinsol_spfgmr.h>  /* access to SPFGMR SUNLinearSolver     */
#include <sundials/sundials_dense.h>     /* use generic dense solver in precond. */
#include <sundials/sundials_types.h>     /* defs. of realtype, sunindextype      */
#include  "helper_electro.h"
#include "mesh.h"
#include "helper_electro.h"

#ifdef _OPENMP

#include <omp.h>
#include <nvector/nvector_openmp.h>    /* access to OpenMP N_Vector            */
#define  NVITH NV_Ith_OMP

#else

#include <nvector/nvector_serial.h>    /* access to Serial N_Vector            */
#define  NVITH NV_Ith_S

#endif



int mod(int a, int b); // modulus

#include "kinsol_soma.h" 

/* Problem Constants */

/* Linear Solver Loop Constants */

#define USE_SPGMR   0
#define USE_SPBCGS  1
#define USE_SPTFQMR 2
#define USE_SPFGMR  3


/* Functions Called by the KINSOL Solver */

static int funcJ(N_Vector cc, N_Vector fval, void *user_data);


/* Template for preconditioner, currently not in use */

static int PrecSetupJ(N_Vector cc, N_Vector cscale,
                       N_Vector fval, N_Vector fscale,
                       void *user_data);

static int PrecSolveJ(N_Vector cc, N_Vector cscale,
                       N_Vector fval, N_Vector fscale,
                       N_Vector vv, void *user_data);



/* Private Helper Functions */

static Phase *AllocUserData(void);
static void SetInitialProfilesJ(N_Vector cc);
static realtype SetScaleJ(const struct Phase *const p);
static int check_flag(void *flagvalue, const char *funcname, int opt);

int iters;

 /*
 *--------------------------------------------------------------------
 * MAIN ROUTINE
 *--------------------------------------------------------------------
 */

int call_J(struct Phase *const p)
{
  static realtype *ccx; // last solution

  unsigned int ix,iy,iz,i,cell, cellp, cellm; 
  int globalstrategy, linsolver;
  realtype fnormtol, scsteptol; // tolerances
  N_Vector cc, sc, constraints;
  static int flagsolved = 1; // turn to 0 after first solution 
  static realtype scale; 
  int num_threads;
  int flag, maxl, maxlrst, mset;
  void *kmem;
  SUNLinearSolver LS;
  Phase *data;
  soma_scalar_t cions[p->n_cells]; // concentration
  soma_scalar_t sumions;
  const soma_scalar_t alfa = p->args.noneq_ratio_arg;
  realtype fnorm;
  soma_scalar_t current0, currentL ;

  soma_scalar_t  eps[p->n_cells]; // c/ceq

/* Kinsol runs on CPU only, update fields */
#pragma acc update self(p->exp_born_pos[0:p->n_cells])
#pragma acc update self(p->exp_born_neg[0:p->n_cells])

// Calc Born_S

  calc_born_S(p);   

  int NEQ; //<- Number of equations 
  NEQ = (int) p->nx*p->ny*(p->nz-2); /* the concentration is fixed near electrodes */

  iters = 0; // number of iterations

  /* Create the SUNDIALS context object for this simulation. */
  SUNContext sunctx = NULL;
  SUNContext_Create(NULL, &sunctx);

  cc = sc = constraints = NULL;
  kmem = NULL;
  LS = NULL;
  data = NULL; 

  /* Allocate memory, and set problem data, initial values, tolerances */
  globalstrategy = KIN_NONE ; /* KIN_NONE = basic Newton iteration
				KIN_LINESEARCH = Newton with globalization
				KIN_FP = fixed point interaction
				KIN_PICARD = Picard interaction */

  /* Set the number of threads to use */
  num_threads = 1;     /* default value*/
#ifdef _OPENMP
  num_threads = p->args.omp_threads_arg;
#endif

  data = AllocUserData(); 
  if (check_flag((void *)data, "AllocUserData", 2)) return(1);

  *data = *p; // Pointer to phase information

//  InitUserData(data);

  /* Create serial vectors of length NEQ */

#ifdef _OPENMP
cc = N_VNew_OpenMP(NEQ, num_threads, sunctx);
if (check_flag((void *)cc, "N_VNew_OpenMP", 0)) return(1);
sc = N_VNew_OpenMP(NEQ, num_threads, sunctx);
if (check_flag((void *)sc, "N_VNew_OpenMP", 0)) return(1);
constraints = N_VNew_OpenMP(NEQ, num_threads, sunctx);
if (check_flag((void *)constraints, "N_VNew_OpenMP", 0)) return(1);
#else
cc = N_VNew_Serial(NEQ, sunctx);
if (check_flag((void *)cc, "N_VNew_Serial", 0)) return(1);
sc = N_VNew_Serial(NEQ, sunctx);
if (check_flag((void *)sc, "N_VNew_Serial", 0)) return(1);
constraints = N_VNew_Serial(NEQ, sunctx);
if (check_flag((void *)constraints, "N_VNew_Serial", 0)) return(1);
#endif

N_VConst(1.0, constraints);  // constrains c >= 0

  linsolver = 1  ; // linear solver, use 0 = SPGMR, 1 = SPBCGS, 2 = SPTFQMR, 3 = SPFGMR

    /* Allocate ccx */
   if (flagsolved) {
	   ccx = (realtype*)malloc(NEQ*sizeof(realtype));
	   if (ccx == NULL) return(1);
    } 

    /* Initial guess */

    /* (Re-)Initialize user data */

   fnormtol = 1e-5;   
   scsteptol = 1e-10; 


   // Calc ions in equilibrium
  
   call_EN(p);

   if (flagsolved)  {   
   // initial guess
        for (ix = 0 ; ix < p->nx ; ix++) {
         for (iy = 0 ; iy < p->ny ; iy++) {
     	  for (iz = 1 ; iz < p->nz-1 ; iz++) {
              i = iz + (p->nz-2)*iy + (p->nz-2)*p->ny*ix - 1 ;
              cell = cell_coordinate_to_index(p, ix, iy, iz);
              NVITH(cc,i) = 1.0 ;
	  }
          }
        } 
   }
   else {
       for (ix = 0 ; ix < p->nx ; ix++) {
         for (iy = 0 ; iy < p->ny ; iy++) {
     	  for (iz = 1 ; iz < p->nz-1 ; iz++) {
              i = iz + (p->nz-2)*iy + (p->nz-2)*p->ny*ix - 1 ;
              cell = cell_coordinate_to_index(p, ix, iy, iz);
              NVITH(cc,i) = ccx[i] ;
	  }
          }
        } 




   }



    /* Set scale vector */
    if (flagsolved) scale = SetScaleJ(p);
    N_VConst(scale, sc);

    /* Call KINCreate/KINInit to initialize KINSOL:
       A pointer to KINSOL problem memory is returned and stored in kmem. */
    kmem = KINCreate(sunctx);
    if (check_flag((void *)kmem, "KINCreate", 0)) return(1);

    /* Vector cc passed as template vector. */
    flag = KINInit(kmem, funcJ, cc);
    if (check_flag(&flag, "KINInit", 1)) return(1);

    flag = KINSetUserData(kmem, data);
    if (check_flag(&flag, "KINSetUserData", 1)) return(1);

    flag = KINSetConstraints(kmem, constraints);  // CONSTRAINTS NO NEEDED
    if (check_flag(&flag, "KINSetConstraints", 1)) return(1);

    flag = KINSetFuncNormTol(kmem, fnormtol);
    if (check_flag(&flag, "KINSetFuncNormTol", 1)) return(1);
    flag = KINSetScaledStepTol(kmem, scsteptol);
    if (check_flag(&flag, "KINSetScaledStepTol", 1)) return(1);

    /* Attach a linear solver module */
    switch(linsolver) {

    /* (a) SPGMR */
    case(USE_SPGMR):

      /* Print header */
//      printf(" -------");
//      printf(" \n| SPGMR |\n");
//      printf(" -------\n");

      /* Create SUNLinSol_SPGMR object with right preconditioning and the
         maximum Krylov dimension maxl */
      maxl = 1000;

//      LS = SUNLinSol_SPGMR(cc, SUN_PREC_NONE, maxl, sunctx);
//      if(check_flag((void *)LS, "SUNLinSol_SPGMR", 0)) return(1); 

      LS = SUNLinSol_SPGMR(cc, SUN_PREC_RIGHT, maxl, sunctx);
      if(check_flag((void *)LS, "SUNLinSol_SPGMR", 0)) return(1);

      /* Attach the linear solver to KINSOL */
      flag = KINSetLinearSolver(kmem, LS, NULL);
      if (check_flag(&flag, "KINSetLinearSolver", 1)) return 1;

      /* Set the maximum number of restarts */
      maxlrst = 1000;
      flag = SUNLinSol_SPGMRSetMaxRestarts(LS, maxlrst);
      if (check_flag(&flag, "SUNLinSol_SPGMRSetMaxRestarts", 1)) return(1);

      break;

    /* (b) SPBCGS */
    case(USE_SPBCGS):

      /* Print header */
//      printf(" --------");
//      printf(" \n| SPBCGS |\n");
//      printf(" --------\n");

      /* Create SUNLinSol_SPBCGS object and the
         maximum Krylov dimension maxl */
      maxl = 1000;

//      LS = SUNLinSol_SPBCGS(cc, SUN_PREC_NONE, maxl, sunctx);
//      if(check_flag((void *)LS, "SUNLinSol_SPBCGS", 0)) return(1); 

      LS = SUNLinSol_SPBCGS(cc, SUN_PREC_RIGHT, maxl, sunctx);
      if(check_flag((void *)LS, "SUNLinSol_SPBCGS", 0)) return(1); 

      /* Attach the linear solver to KINSOL */
      flag = KINSetLinearSolver(kmem, LS, NULL);
      if (check_flag(&flag, "KINSetLinearSolver", 1)) return 1;

      /* Set the maximum number of restarts */
      maxlrst = 10;
      flag = SUNLinSol_SPGMRSetMaxRestarts(LS, maxlrst);
      if (check_flag(&flag, "SUNLinSol_SPGMRSetMaxRestarts", 1)) return(1);

      break;

    /* (c) SPTFQMR */
    case(USE_SPTFQMR):

      /* Print header */
//      printf(" ---------");
//      printf(" \n| SPTFQMR |\n");
//      printf(" ---------\n");

      /* Create SUNLinSol_SPTFQMR object with right preconditioning and the
         maximum Krylov dimension maxl */
      maxl = 1000;

/*      LS = SUNLinSol_SPTFQMR(cc, SUN_PREC_NONE, maxl, sunctx);
      if(check_flag((void *)LS, "SUNLinSol_SPTFQMR", 0)) return(1);
*/

      LS = SUNLinSol_SPTFQMR(cc, SUN_PREC_RIGHT, maxl, sunctx);
      if(check_flag((void *)LS, "SUNLinSol_SPTFQMR", 0)) return(1);

      /* Attach the linear solver to KINSOL */
      flag = KINSetLinearSolver(kmem, LS, NULL);
      if (check_flag(&flag, "KINSetLinearSolver", 1)) return 1;

      break;

    /* (d) SPFGMR */
    case(USE_SPFGMR):

      /* Print header */
 //     printf(" -------");
 //     printf(" \n| SPFGMR |\n");
 //     printf(" -------\n");

      /* Create SUNLinSol_SPFGMR object with right preconditioning and the
         maximum Krylov dimension maxl */
      maxl = 1000;
/*      LS = SUNLinSol_SPFGMR(cc, SUN_PREC_NONE, maxl, sunctx);
      if(check_flag((void *)LS, "SUNLinSol_SPFGMR", 0)) return(1); */

      LS = SUNLinSol_SPFGMR(cc, SUN_PREC_RIGHT, maxl, sunctx);
      if(check_flag((void *)LS, "SUNLinSol_SPFGMR", 0)) return(1);


      /* Attach the linear solver to KINSOL */
      flag = KINSetLinearSolver(kmem, LS, NULL);
      if (check_flag(&flag, "KINSetLinearSolver", 1)) return 1;

      /* Set the maximum number of restarts */
      maxlrst = 100;
      flag = SUNLinSol_SPGMRSetMaxRestarts(LS, maxlrst);
      if (check_flag(&flag, "SUNLinSol_SPGMRSetMaxRestarts", 1)) return(1);

      break;

    }

    /* Set preconditioner functions*/
    flag = KINSetPreconditioner(kmem, PrecSetupJ, PrecSolveJ);
    if (check_flag(&flag, "KINSetPreconditioner", 1)) return(1);

    mset = 1; // maximum number of iterations before recalc diagonal preconditioner

    flag = KINSetMaxSetupCalls(kmem, mset);
    if (check_flag(&flag, "KINSetMaxSetupCalls", 1)) return(1);
    
    /* Call KINSol */

    flag = KINSol(kmem,           /* KINSol memory block */
		  cc,             /* initial guess on input; solution vector */
		  globalstrategy, /* global strategy choice */
		  sc,             /* scaling vector, for the variable cc */
		  sc);            /* scaling vector for function values fval */



    if (check_flag(&flag, "KINSol", 1)) return(1);

        KINGetFuncNorm(kmem, &fnorm);
//        printf("flag %d \n", flag);
    if (((flag == 0)||(flag == 1)||(flag == 2))&&(!isnan(fnorm))) {  // converged
							       //

        p->aviter += iters;
        p->countiter++;
 
/* Save solution */
        // Save profile  
        soma_scalar_t avpsi = 0; //average psi
        for (i = 0 ; i < NEQ ; i++) {
        	ccx[i] = NVITH(cc,i);
                flagsolved = 0;
         } // converged
	   
    } else {  // did not converged
        if (p->info_MPI.sim_rank == 0) 
             fprintf(stdout, "Kinsol failed to converge last step \n");
    }

/* Calc electric field */

// recover c from kinsol

// Transform from ix, iy, iz to kinsol's index: (the calculation box is smaller in the z direction than the simulation box)
// index = iz + (nz-2)*iy + (nz-2)*ny*ix - 1



  for (ix = 0 ; ix < p->nx ; ix++) {
	  for (iy = 0 ; iy < p->ny ; iy++) {
			  for (iz = 1 ; iz < p->nz-1 ; iz++) {
                          i = iz + (p->nz-2)*iy + (p->nz-2)*p->ny*ix - 1 ;
			  cell = cell_coordinate_to_index(p, ix, iy, iz);
	                  cions[cell] = NVITH(cc,i)*p->npos_field[cell]/p->npos_field[0]; // cions not normalized 
		   }
           }
   }


// fill borders
  iz = 0;   
#pragma omp parallel for  
  for (ix = 0 ; ix < p->nx ; ix++) {
	  for (iy = 0 ; iy < p->ny ; iy++) {
               cell = cell_coordinate_to_index(p, ix, iy, iz);
	       cions[cell] = alfa*p->npos_field[cell]/p->npos_field[0];
           }
   }

  iz = p->nz-1;   
#pragma omp parallel for  
  for (ix = 0 ; ix < p->nx ; ix++) {
	  for (iy = 0 ; iy < p->ny ; iy++) {
               cell = cell_coordinate_to_index(p, ix, iy, iz);
	       cions[cell] = p->npos_field[cell]/p->npos_field[0]; 
           }
   }


/* Normalize c and save ion densities */

sumions = 0.0;
#pragma omp parallel for reduction(+:sumions) 
    for (cell = 0 ; cell < p->n_cells ; cell++) {
         sumions += cions[cell]; 
    }

// Normalize to match totalnumber and divide by cell volume
#pragma omp parallel for  
    for (cell = 0 ; cell < p->n_cells ; cell++) {
           cions[cell] = cions[cell] * p->Nposions / sumions / p->vcell;
    }

#pragma omp parallel for  
    for (cell = 0 ; cell < p->n_cells ; cell++) {
	      eps[cell] = cions[cell]/p->npos_field[cell];
    }


// Calculation of ion currents at electrode

current0 = 0.0;
currentL = 0.0;

  for (ix = 0 ; ix < p->nx ; ix++) {
	  for (iy = 0 ; iy < p->ny ; iy++) {

            iz = 0; 	  
	    cellm = cell_coordinate_to_index(p, ix, iy, iz);
            iz = 1; 	  
	    cell = cell_coordinate_to_index(p, ix, iy, iz);


	    current0 -= (p->npos_field[cell])*(eps[cell]-eps[cellm]);
	    current0 -= (p->npos_field[cellm])*(eps[cell]-eps[cellm]);
			    

	    iz = p->nz-2; 	  
	    cell = cell_coordinate_to_index(p, ix, iy, iz);
            iz = p->nz-1; 	  
	    cellp = cell_coordinate_to_index(p, ix, iy, iz);

	    currentL -= (p->npos_field[cellp])*(eps[cellp]-eps[cell]);
	    currentL -= (p->npos_field[cell])*(eps[cellp]-eps[cell]);
       
 
          }
   }

  current0 = current0 * p->deltax*p->deltay/p->deltaz;
  currentL = currentL * p->deltax*p->deltay/p->deltaz;


// Save non-eq ion densities
#pragma omp parallel for  
    for (cell = 0 ; cell < p->n_cells ; cell++) {
           p->nneg_field[cell] = cions[cell];
           p->npos_field[cell] = cions[cell];
    }


/* Now calculate electric field */
#pragma omp parallel for  
    for (cell = 0 ; cell < p->n_cells ; cell++) {
           p->electric_field[cell] = log(p->nneg_field[cell])-log(p->exp_born_neg[cell]); 
    }


// print    
        printf("Transport converged, flag %d, iters %d, norm %.3e, normtol %.3e, I(0) %.3e, I(L) %.3e \n", flag, iters, fnorm, fnormtol, current0, currentL);

    
	/* Free memory */

    KINFree(&kmem);
    SUNLinSolFree(LS);

  N_VDestroy(constraints);
  N_VDestroy(cc);
  N_VDestroy(sc);
/*  FreeUserData(data); */

  SUNContext_Free(&sunctx);
  return(0);
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY KINSOL
 *--------------------------------------------------------------------
 */

static int funcJ(N_Vector cc, N_Vector fval, void *user_data)
{

#include <assert.h>

  unsigned int ix, iy, iz, cell, i;
  int ixp ,ixm, iyp, iym, izp, izm;
  struct Phase *const p = user_data;
  const soma_scalar_t alfa = p->args.noneq_ratio_arg;

  int NEQ; //<- Number of equations 
  NEQ = (int) p->nx*p->ny*(p->nz-2); /* the concentration is fixed near electrodes */

  soma_scalar_t  res[p->nx][p->ny][p->nz]; // residual Poisson Eq.
  soma_scalar_t  c[p->nx][p->ny][p->nz]; // ion concetration
  soma_scalar_t  eps[p->nx][p->ny][p->nz]; // auxiliary field
 

  iters++;	   

// born_S
soma_scalar_t  born_S[p->nx][p->ny][p->nz];
#pragma omp parallel for  
  for (ix = 0 ; ix < p->nx ; ix++) {
	  for (iy = 0 ; iy < p->ny ; iy++) {
	        for (iz = 0 ; iz < p->nz ; iz++) {
	        cell = cell_coordinate_to_index(p, ix, iy, iz);
	        born_S[ix][iy][iz] = p->born_Sc[cell]; 
	        }
          }
   }



// c from npos_ions
for (ix = 0 ; ix < p->nx ; ix++) {
	  for (iy = 0 ; iy < p->ny ; iy++) {
			  for (iz = 0 ; iz <  p->nz ; iz++) {
                          cell = cell_coordinate_to_index(p, ix, iy, iz);
                          c[ix][iy][iz] = p->npos_field[cell]/p->npos_field[0];
			  }
		     }
	 	}


// epsilon from kinsol's input

// Transform from ix, iy, iz to kinsol's index: (the calculation box is smaller in the z direction than the simulation box)
// index = iz + (nz-2)*iy + (nz-2)*ny*ix - 1

  for (ix = 0 ; ix < p->nx ; ix++) {
	  for (iy = 0 ; iy < p->ny ; iy++) {
			  for (iz = 1 ; iz < p->nz-1 ; iz++) {
                          i = iz + (p->nz-2)*iy + (p->nz-2)*p->ny*ix - 1 ;
	                  eps[ix][iy][iz] = NVITH(cc,i);
	           }
           }
   }

// fill borders
  iz = 0;   
#pragma omp parallel for  
  for (ix = 0 ; ix < p->nx ; ix++) {
	  for (iy = 0 ; iy < p->ny ; iy++) {
	       eps[ix][iy][iz] = alfa;  
           }
   }

  iz = p->nz-1;
#pragma omp parallel for  
  for (ix = 0 ; ix < p->nx ; ix++) {
	  for (iy = 0 ; iy < p->ny ; iy++) {
	       eps[ix][iy][iz] = 1.0; 
           }
   }

/*
  for (ix = 0 ; ix < p->nx ; ix++) {
  for (iy = 0 ; iy < p->ny ; iy++) {
  for (iz = 0 ; iz < p->nz ; iz++) {
        printf("i %d %d %d %f %f \n", ix, iy, iz, c[ix][iy][iz], born_S[ix][iy][iz]);
  }
  }
  }
*/

//soma_scalar_t slope = (1-alfa)/p->Lz;

// DO NOT PARALELIZE HERE  
  for (ix = 0 ; ix < p->nx ; ix++) {

     ixp = mod((ix+1),p->nx);
     ixm = mod((ix-1),p->nx);
 
     for (iy = 0 ; iy < p->ny ; iy++) {

        iyp = mod((iy+1),p->ny);
        iym = mod((iy-1),p->ny);

#pragma omp parallel for  
	for (iz = 1 ; iz < p->nz-1 ; iz++) {
      
	izp = iz+1;       	
	izm = iz-1;       	
     
        cell = cell_coordinate_to_index(p, ix, iy, iz); // cell in simulation box

        res[ix][iy][iz] = 0.0;
        res[ix][iy][iz] += 0.5*((c[ixp][iy][iz]+c[ix][iy][iz])*(eps[ixp][iy][iz]-eps[ix][iy][iz]))/(p->deltax*p->deltax);
        res[ix][iy][iz] += 0.5*(-(c[ix][iy][iz]+c[ixm][iy][iz])*(eps[ix][iy][iz]-eps[ixm][iy][iz]))/(p->deltax*p->deltax);
        res[ix][iy][iz] += 0.5*((c[ix][iyp][iz]+c[ix][iy][iz])*(eps[ix][iyp][iz]-eps[ix][iy][iz]))/(p->deltay*p->deltay);
        res[ix][iy][iz] += 0.5*(-(c[ix][iy][iz]+c[ix][iym][iz])*(eps[ix][iy][iz]-eps[ix][iym][iz]))/(p->deltay*p->deltay);
        res[ix][iy][iz] += 0.5*((c[ix][iy][izp]+c[ix][iy][iz])*(eps[ix][iy][izp]-eps[ix][iy][iz]))/(p->deltaz*p->deltaz);
        res[ix][iy][iz] += 0.5*(-(c[ix][iy][iz]+c[ix][iy][izm])*(eps[ix][iy][iz]-eps[ix][iy][izm]))/(p->deltaz*p->deltaz);
	
        }
    }
  }

// Transform from ix, iy, iz to kinsol's index: (the calculation box is smaller in the z direction than the simulation box)
// index = iz + (nz-2)*iy + (nz-2)*ny*ix - 1



// DO NOT PARALELIZE #pragma omp parallel for  
  for (ix = 0 ; ix < p->nx ; ix++) {
     for (iy = 0 ; iy < p->ny ; iy++) {
	for (iz = 1 ; iz < p->nz-1 ; iz++) {

         i = iz + (p->nz-2)*iy + (p->nz-2)*p->ny*ix - 1 ;
         NVITH(fval,i) = res[ix][iy][iz];
//         printf("func: cell, res %d %d %f %f \n", p->iter, i, res[ix][iy][iz], c[ix][iy][iz]);
		     }
	 	}
	  }

 

// DEBUG print norm 
soma_scalar_t norma = 0;
        for (ix = 0 ; ix < p->nx ; ix++) {
               for (iy = 0 ; iy < p->ny ; iy++) {
                  for (iz = 1 ; iz <  p->nz-1 ; iz++) {
                  cell = cell_coordinate_to_index(p, ix, iy, iz);
              			  norma += fabs(res[ix][iy][iz]); 
                     }
                }
         }
  printf("func: iter, norma, res(nx,ny,nz): %d %f %f \n ", iters, norma, res[0][0][0]); 

  
//  printf("func: Nposions, Nnegions: %f, %f \n ", p->Nposions, p->Nnegions);
//  printf("func: Number of Equations: %d \n", NEQ);

//  printf("func: iter, norm %d %.3e \n", iter, norma);

//  exit(1);
  return(0);
}
  

/*
 * Set initial conditions in cc
 */


static realtype SetScaleJ(const struct Phase *const p)
{
   realtype scale;

   scale = 1.0 ;
           
   return(scale);
   }

static void SetInitialProfilesJ(N_Vector cc)
{ 
  N_VConst(1.0, cc);  
}
 
static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr,
            "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr,
              "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1);
    }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr,
            "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1);
  }

  return(0);
}

/*
 * Allocate memory for data structure of type Phase
 */

static Phase *AllocUserData(void)
{
  Phase *data;
  
  data = malloc(sizeof *data);
//  printf("data size %lu \n", sizeof *data);
  return(data);
}



/*
 * Preconditioner setup routine. Generate and preprocess P.
 */

static int PrecSetupJ(N_Vector cc, N_Vector cscale,
                       N_Vector fval, N_Vector fscale,
                       void *user_data)  {

#include <assert.h>

  unsigned int ix, iy, iz, i;
  int ixp ,ixm, iyp, iym, izp, izm, cell;
  struct Phase *const p = user_data;
  soma_scalar_t c[p->nx][p->ny][p->nz]; // concentration
  soma_scalar_t lc[p->nx][p->ny][p->nz]; // concentration
  int NEQ;
  NEQ = (int) p->nx*p->ny*(p->nz-2); /* the concentration is fixed near electrodes */
  const soma_scalar_t alfa = p->args.noneq_ratio_arg;
  soma_scalar_t  eps[p->nx][p->ny][p->nz]; // auxiliary field

// born_S
soma_scalar_t  born_S[p->nx][p->ny][p->nz];
#pragma omp parallel for  
  for (ix = 0 ; ix < p->nx ; ix++) {
	  for (iy = 0 ; iy < p->ny ; iy++) {
	        for (iz = 0 ; iz < p->nz ; iz++) {
	        cell = cell_coordinate_to_index(p, ix, iy, iz);
	        born_S[ix][iy][iz] = p->born_Sc[cell]; 
	        }
          }
   }



// c from npos_ions
for (ix = 0 ; ix < p->nx ; ix++) {
	  for (iy = 0 ; iy < p->ny ; iy++) {
			  for (iz = 0 ; iz <  p->nz ; iz++) {
                          cell = cell_coordinate_to_index(p, ix, iy, iz);
                          c[ix][iy][iz] = p->npos_field[cell]/p->npos_field[0];
			  }
		     }
	 	}


// epsilon from kinsol's input

// Transform from ix, iy, iz to kinsol's index: (the calculation box is smaller in the z direction than the simulation box)
// index = iz + (nz-2)*iy + (nz-2)*ny*ix - 1

  for (ix = 0 ; ix < p->nx ; ix++) {
	  for (iy = 0 ; iy < p->ny ; iy++) {
			  for (iz = 1 ; iz < p->nz-1 ; iz++) {
                          i = iz + (p->nz-2)*iy + (p->nz-2)*p->ny*ix - 1 ;
	                  eps[ix][iy][iz] = NVITH(cc,i);
	           }
           }
   }

// fill borders
  iz = 0;   
#pragma omp parallel for  
  for (ix = 0 ; ix < p->nx ; ix++) {
	  for (iy = 0 ; iy < p->ny ; iy++) {
	       eps[ix][iy][iz] = alfa;  
           }
   }

  iz = p->nz-1;
#pragma omp parallel for  
  for (ix = 0 ; ix < p->nx ; ix++) {
	  for (iy = 0 ; iy < p->ny ; iy++) {
	       eps[ix][iy][iz] = 1.0; 
           }
   }

/*
  for (ix = 0 ; ix < p->nx ; ix++) {
  for (iy = 0 ; iy < p->ny ; iy++) {
  for (iz = 0 ; iz < p->nz ; iz++) {
        printf("i %d %d %d %f %f \n", ix, iy, iz, c[ix][iy][iz], born_S[ix][iy][iz]);
  }
  }
  }
*/

//soma_scalar_t slope = (1-alfa)/p->Lz;


/// Calculate diagonal preconditioner, temp_prec_field

// Transform from ix, iy, iz to kinsol's index: (the calculation box is smaller in the z direction than the simulation box)
// index = iz + (nz-2)*iy + (nz-2)*ny*ix - 1
// i = iz + (p->nz-2)*iy + (p->nz-2)*p->ny*ix - 1 ;

  for (ix = 0 ; ix < p->nx ; ix++) {

     ixp = mod((ix+1),p->nx);
     ixm = mod((ix-1),p->nx);
 
     for (iy = 0 ; iy < p->ny ; iy++) {

        iyp = mod((iy+1),p->ny);
        iym = mod((iy-1),p->ny);

	for (iz = 1 ; iz < p->nz-1 ; iz++) {
        
        	izp = iz+1;
 		izm = iz-1;

		 i = iz + (p->nz-2)*iy + (p->nz-2)*p->ny*ix - 1 ;

                 p->temp_prec_field[i] = 0.0;
       		 p->temp_prec_field[i] += -0.5*(c[ixp][iy][iz]+2*c[ix][iy][iz]+c[ixm][iy][iz])/(p->deltax*p->deltax); 
        	 p->temp_prec_field[i] += -0.5*(c[ix][iyp][iz]+2*c[ix][iy][iz]+c[ix][iym][iz])/(p->deltay*p->deltay); 
        	 p->temp_prec_field[i] += -0.5*(c[ix][iy][izp]+2*c[ix][iy][iz]+c[ix][iy][izm])/(p->deltaz*p->deltaz); 
 
	}
      }
    }

return(0);

}	

/*
 * Preconditioner solve routine
 */

static int PrecSolveJ(N_Vector cc, N_Vector cscale,
                       N_Vector fval, N_Vector fscale,
                       N_Vector vv, void *user_data)
{
  unsigned int i;
  struct Phase *const p = user_data;
  int NEQ;
  NEQ = (int) p->nx*p->ny*(p->nz-2); /* the concentration is fixed near electrodes */

  soma_scalar_t vvin[NEQ]; 
  soma_scalar_t vvout[NEQ]; 

  for (i = 0 ; i < NEQ ; i++) {
	  vvin[i] = NVITH(vv,i); 
   }

  for (i = 0 ; i < NEQ ; i++) {
          vvout[i] = vvin[i]/p->temp_prec_field[i]; // Diagonal precond.
	  NVITH(vv,i) = vvout[i]; 
   }

  return(0);
}

