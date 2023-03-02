/* Routines for solving a set of non-linear coupled equations
 * Finds x that satisfy F(x) = 0
 * Based on kinsol  6.5 example "inKrylovDemo_ls.c"  */

#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <kinsol/kinsol.h>               /* access to KINSOL func., consts.      */
#include <nvector/nvector_openmp.h>    /* access to OpenMP N_Vector            */
#include <sunlinsol/sunlinsol_spgmr.h>   /* access to SPGMR SUNLinearSolver      */
#include <sunlinsol/sunlinsol_spbcgs.h>  /* access to SPBCGS SUNLinearSolver     */
#include <sunlinsol/sunlinsol_sptfqmr.h> /* access to SPTFQMR SUNLinearSolver    */
#include <sunlinsol/sunlinsol_spfgmr.h>  /* access to SPFGMR SUNLinearSolver     */
#include <sundials/sundials_dense.h>     /* use generic dense solver in precond. */
#include <sundials/sundials_types.h>     /* defs. of realtype, sunindextype      */
#include "mesh.h"

int mod(int a, int b); // modulus

#include "kinsol_soma.h" 

#ifdef _OPENMP
#include <omp.h>
#endif

/* Problem Constants */

/* Linear Solver Loop Constants */

#define USE_SPGMR   0
#define USE_SPBCGS  1
#define USE_SPTFQMR 2
#define USE_SPFGMR  3


/* Functions Called by the KINSOL Solver */

static int func(N_Vector cc, N_Vector fval, void *user_data);


/* Template for preconditioner, currently not in use */

/*static int PrecSetupBD(N_Vector cc, N_Vector cscale,
                       N_Vector fval, N_Vector fscale,
                       void *user_data);

static int PrecSolveBD(N_Vector cc, N_Vector cscale,
                       N_Vector fval, N_Vector fscale,
                       N_Vector vv, void *user_data);
*/


/* Private Helper Functions */

static Phase *AllocUserData(void);
static void SetInitialProfiles(N_Vector cc);
static realtype SetScale(void);
static int check_flag(void *flagvalue, const char *funcname, int opt);

int iter = 0;
realtype fnorm;
 
 /*
 *--------------------------------------------------------------------
 * MAIN ROUTINE
 *--------------------------------------------------------------------
 */

int call_PB(const struct Phase *const p)
{
  static realtype *ccx; // last solution

  int i;
  int globalstrategy, linsolver;
  realtype fnormtol, scsteptol; // tolerances
  N_Vector cc, sc, constraints;
  static int flagsolved = 1; // turn to 0 after first solution 
  static realtype scale; 
  int num_threads;
  int flag, maxl, maxlrst;
  void *kmem;
  SUNLinearSolver LS;
  Phase *data;

  int NEQ; //<- Number of equations 
  NEQ = (int) p->n_cells_local - 1; /* Due to PBC the set of equations is no longer LI, so psi(nx,ny,nz) can
				       be fixed to zero (see notes) */

  iter = 0; // number of iterations

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
//  printf("Number of threads %d \n", num_threads);
		  

  data = AllocUserData(); 
  if (check_flag((void *)data, "AllocUserData", 2)) return(1);

  *data = *p; // Pointer to phase information

//  InitUserData(data);

  /* Create serial vectors of length NEQ */
  cc = N_VNew_OpenMP(NEQ, num_threads, sunctx);
  if (check_flag((void *)cc, "N_VNew_OpenMP", 0)) return(1);
  sc = N_VNew_OpenMP(NEQ, num_threads, sunctx);
  if (check_flag((void *)sc, "N_VNew_OpenMP", 0)) return(1);

/* Template for using constrainst, not in use 
 constraints = N_VNew_Serial(NEQ, sunctx);
 if (check_flag((void *)constraints, "N_VNew_Serial", 0)) return(1);
 N_VConst(0, constraints); */

  linsolver = 1; // linear solver, use 0 = SPGMR, 1 = SPBCGS, 2 = SPTFQMR, 3 = SPFGMR

    /* Allocate ccx */
   if (flagsolved) {
	   ccx = (realtype*)malloc(NEQ*sizeof(realtype));
	   if (ccx == NULL) return(1);
    } 


    /* Initial guess */


    /* (Re-)Initialize user data */

   fnormtol = 1e-5;   
   scsteptol = 1e-13; 
 
   if (p->efieldsolver == 2) {   // homogeneous as initial guess, then reuse last solution
        if (flagsolved) {
           if (p->info_MPI.sim_rank == 0) 
              fprintf(stdout, "Set initial guess for electrostatics \n");
	   SetInitialProfiles(cc);
        } 
        else {
	    // Recover profile, need to implement in a function   
            for (i = 0 ; i < NEQ ; i++) {
		   NV_Ith_OMP(cc,i) = ccx[i]; 
            }
        }
   }	 
   else if (p->efieldsolver == 0) { // EN as initial
	call_EN(p); 

        for (i = 0 ; i < NEQ ; i++) {
           NV_Ith_OMP(cc,i) = p->electric_field[i] - p->electric_field[NEQ]; // sets efield to zero in the last cell
	}
   }
   else {
       fprintf(stderr, "Invalid value of p->efieldsolver in kinsol_soma.c\n");
       assert(0);
   } 


    /* Set scale vector */
    if (flagsolved) scale = SetScale();
    N_VConst(scale, sc);

    /* Call KINCreate/KINInit to initialize KINSOL:
       A pointer to KINSOL problem memory is returned and stored in kmem. */
    kmem = KINCreate(sunctx);
    if (check_flag((void *)kmem, "KINCreate", 0)) return(1);

    /* Vector cc passed as template vector. */
    flag = KINInit(kmem, func, cc);
    if (check_flag(&flag, "KINInit", 1)) return(1);

    flag = KINSetUserData(kmem, data);
    if (check_flag(&flag, "KINSetUserData", 1)) return(1);

//    flag = KINSetConstraints(kmem, constraints);  // CONSTRAINTS NO NEEDED
//    if (check_flag(&flag, "KINSetConstraints", 1)) return(1);

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
      LS = SUNLinSol_SPGMR(cc, SUN_PREC_NONE, maxl, sunctx);
      if(check_flag((void *)LS, "SUNLinSol_SPGMR", 0)) return(1);

      /* Attach the linear solver to KINSOL */
      flag = KINSetLinearSolver(kmem, LS, NULL);
      if (check_flag(&flag, "KINSetLinearSolver", 1)) return 1;

      /* Set the maximum number of restarts */
      maxlrst = 100;
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
      LS = SUNLinSol_SPBCGS(cc, SUN_PREC_NONE, maxl, sunctx);
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
      LS = SUNLinSol_SPTFQMR(cc, SUN_PREC_NONE, maxl, sunctx);
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
      LS = SUNLinSol_SPFGMR(cc, SUN_PREC_NONE, maxl, sunctx);
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

    /* Set preconditioner functions, not in use*/
/*    flag = KINSetPreconditioner(kmem, PrecSetupBD, PrecSolveBD);
    if (check_flag(&flag, "KINSetPreconditioner", 1)) return(1);
*/
    
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
//        printf("Elec. converged, flag %d, iters %d, norm %.3e, normtol %.3e \n", flag, iter, fnorm, fnormtol);
        /* Save solution */
        // Save profile, need to implement in a function   
        soma_scalar_t avpsi = 0; //average psi
        for (i = 0 ; i < NEQ ; i++) {
        	ccx[i] = NV_Ith_OMP(cc,i);
                avpsi += ccx[i]; 
        	p->electric_field[i] = ccx[i];
        }
        p->electric_field[p->n_cells_local-1] = 0.0;
        avpsi = avpsi / ((soma_scalar_t) p->n_cells_local); 

        for (i = 0 ; i < (int) p->n_cells_local ; i++) {
        	p->electric_field[i] += -avpsi;
        }
        flagsolved = 0;
    } // converge
    else {  // did not converged
        if (p->info_MPI.sim_rank == 0) 
             fprintf(stdout, "Kinsol failed to converge last step, restart initial guess \n");
    }

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

static int func(N_Vector cc, N_Vector fval, void *user_data)
{
//  realtype xx, yy, delx, dely, *cxy, *rxy, *fxy, dcyli, dcyui, dcxli, dcxri;
//  int jx, jy, is, idyu, idyl, idxr, idxl;

#include <assert.h>

  unsigned int ix, iy, iz, cell;
  int ixp ,ixm, iyp, iym, izp, izm;
  const struct Phase *const p = user_data;
  int NEQ = (int) p->n_cells_local - 1; /* Due to PBC the set of equations is no longer LI, so psi(nx,ny,nz) can
				       be fixed to zero (see notes) */


  const soma_scalar_t  deltax = p->Lx/((soma_scalar_t) p->nx);
  const soma_scalar_t  deltay = p->Ly/((soma_scalar_t) p->ny);
  const soma_scalar_t  deltaz = p->Lz/((soma_scalar_t) p->nz);
  const soma_scalar_t  vcell = deltax*deltay*deltaz;
  soma_scalar_t  res[p->nx][p->ny][p->nz]; // residual Poisson Eq.
  soma_scalar_t  rhoposion[p->n_cells_local]; // number density of positive ions in 3D lattice
  soma_scalar_t  rhonegion[p->n_cells_local]; // number density of negative ions in 3D lattice
  soma_scalar_t  rhoQ[p->nx][p->ny][p->nz]; // total charge density
  soma_scalar_t  psi[p->nx][p->ny][p->nz]; // electrostatic potential 3D lattice, units of kBT/|e|
  soma_scalar_t  psic[p->n_cells_local]; // electrostatic potential 3D lattice, units of kBT/|e|
  soma_scalar_t  sumposions = 0 , sumnegions = 0 ;
  soma_scalar_t  sumrhoQ = 0 ; // sum of rhoQ, for debug only
  
  soma_scalar_t constq = 4.0*M_PI; // multiplicative constant for Poisson equation
  soma_scalar_t invbl[p->nx][p->ny][p->nz]; // inverse of local Bjerrum length

  iter++;	   

// psi from kinsol's input

#pragma omp parallel for  
  for (ix = 0 ; ix < p->nx ; ix++) {
	  for (iy = 0 ; iy < p->ny ; iy++) {
			  for (iz = 0 ; iz <  p->nz ; iz++) {
                          cell = cell_coordinate_to_index(p, ix, iy, iz);
                          if ((int) cell < NEQ) {
				  psi[ix][iy][iz] = NV_Ith_OMP(cc,cell); 
				  psic[cell] = NV_Ith_OMP(cc,cell); 
			  }
		     }
	 	}
	  }

psi[p->nx-1][p->ny-1][p->nz-1] = 0.0; // choice of zero of electrostatic potential due to PBC
psic[p->n_cells_local] = 0.0; // choice of zero of electrostatic potential due to PBC

// Pos ion 


if (p->Nposions > 0) {
#pragma omp parallel for reduction(+:sumposions) 
     for (cell = 0 ; cell < p->n_cells_local ; cell++) {
            rhoposion[cell] = exp(-psic[cell])*p->exp_born[cell];
	    sumposions += rhoposion[cell]; 
      }

// Normalize to match totalnumber and divide by cell volume

#pragma omp parallel for  
    for (cell = 0 ; cell < p->n_cells_local ; cell++) {
            rhoposion[cell] = rhoposion[cell] * p->Nposions / sumposions / vcell; 
    }
} // if

else {

#pragma omp parallel for  
    for (cell = 0 ; cell < p->n_cells_local ; cell++) {
             rhoposion[cell] = 0.0;
    }
}

// Neg ion 
if (p->Nnegions > 0) {

#pragma omp parallel for reduction(+:sumnegions) 
    for (cell = 0 ; cell < p->n_cells_local ; cell++) {
             rhonegion[cell] = exp(psic[cell])*p->exp_born[cell];
	     sumnegions += rhonegion[cell]; 
    }

// Normalize to match totalnumber and divide by cell volume
#pragma omp parallel for  
    for (cell = 0 ; cell < p->n_cells_local ; cell++) {
           rhonegion[cell] = rhonegion[cell] * p->Nnegions / sumnegions / vcell; 
    }
} // if 
else {

#pragma omp parallel for  
     for (cell = 0 ; cell < p->n_cells_local ; cell++) {
           rhonegion[cell] = 0.0;
     }
}

// rhoQ from unified fields 

#pragma omp parallel for reduction(+:sumrhoQ) 
    for (ix = 0 ; ix < p->nx ; ix++) {
	  for (iy = 0 ; iy < p->ny ; iy++) {
		  for (iz = 0 ; iz < p->nz ; iz++) {
                          cell = cell_coordinate_to_index(p, ix, iy, iz);

			  rhoQ[ix][iy][iz] = p->rhoF[cell];
                          rhoQ[ix][iy][iz] += rhoposion[cell]-rhonegion[cell]; 
			  sumrhoQ += rhoQ[ix][iy][iz];	  
		     }
	 	}
	  }

// recast inverse of average bjerrum length into x,y,z coordinates

    for (ix = 0 ; ix < p->nx ; ix++) {
	  for (iy = 0 ; iy < p->ny ; iy++) {
		  for (iz = 0 ; iz < p->nz ; iz++) {
                          cell = cell_coordinate_to_index(p, ix, iy, iz);
                          invbl[ix][iy][iz] = p->invblav[cell] ;
//				  printf("func: cell, res %d %f \n", cell, p->invblav[cell]);
		     }
	 	}
	  }

/// Calculate residual from Poisson's equation

#pragma omp parallel for  
  for (ix = 0 ; ix < p->nx ; ix++) {

     ixp = mod((ix+1),p->nx);
     ixm = mod((ix-1),p->nx);
 
     for (iy = 0 ; iy < p->ny ; iy++) {

        iyp = mod((iy+1),p->ny);
        iym = mod((iy-1),p->ny);

	for (iz = 0 ; iz < p->nz ; iz++) {
        
	izp = mod((iz+1),p->nz);
	izm = mod((iz-1),p->nz);
        cell = cell_coordinate_to_index(p, ix, iy, iz);

	res[ix][iy][iz] = rhoQ[ix][iy][iz]*constq;

        res[ix][iy][iz] += 0.5*((invbl[ixp][iy][iz]+invbl[ix][iy][iz])*(psi[ixp][iy][iz]-psi[ix][iy][iz]))/(deltax*deltax); 
        res[ix][iy][iz] += 0.5*(-(invbl[ix][iy][iz]+invbl[ixm][iy][iz])*(psi[ix][iy][iz]-psi[ixm][iy][iz]))/(deltax*deltax); 
        res[ix][iy][iz] += 0.5*((invbl[ix][iyp][iz]+invbl[ix][iy][iz])*(psi[ix][iyp][iz]-psi[ix][iy][iz]))/(deltay*deltay); 
        res[ix][iy][iz] += 0.5*(-(invbl[ix][iy][iz]+invbl[ix][iym][iz])*(psi[ix][iy][iz]-psi[ix][iym][iz]))/(deltay*deltay); 
        res[ix][iy][iz] += 0.5*((invbl[ix][iy][izp]+invbl[ix][iy][iz])*(psi[ix][iy][izp]-psi[ix][iy][iz]))/(deltaz*deltaz); 
        res[ix][iy][iz] += 0.5*(-(invbl[ix][iy][iz]+invbl[ix][iy][izm])*(psi[ix][iy][iz]-psi[ix][iy][izm]))/(deltaz*deltaz); 

	res[ix][iy][iz] = -res[ix][iy][iz];
                          
			  if ((int) cell < NEQ) { 
				  NV_Ith_OMP(fval,cell) = -res[ix][iy][iz];
//				  printf("func: cell, res %d %f \n", cell, invbl[ix][iy][iz]);
			  }     
		     }
	 	}
	  }


/* DEBUG print norm 
soma_scalar_t norma = 0;
        for (ix = 0 ; ix < (int) p->nx ; ix++) {
               for (iy = 0 ; iy < (int) p->ny ; iy++) {
                  for (iz = 0 ; iz < (int)  p->nz ; iz++) {
                  cell = cell_coordinate_to_index(p, ix, iy, iz);
			  if (cell < NEQ)  
              			  norma += res[ix][iy][iz]; 
                          }
                     }
                }

  printf("func: iter, norma, res(nx,ny,nz): %d %f %f \n ", iter, norma, res[p->nx-1][p->ny-1][p->nz-1]); 
*/
  assert(fabs(sumrhoQ) < 1.0e-5);

//  printf("func: Nposions, Nnegions: %f, %f \n ", p->Nposions, p->Nnegions);
//  printf("func: Number of Equations: %d \n", NEQ);

//  printf("func: iter, norm %d %.3e \n", iter, norma);

//  exit(1);
  return(0);
}
  

/*
 * Set initial conditions in cc
 */

static realtype SetScale(void)
{ 
   realtype scale;
   soma_scalar_t constq = 4.0*M_PI; // multiplicative constant for Poisson equation
   
   scale = 1./constq/1. ;  // the factor 1. is an estimation for the average charge per lattice site
	   
   return(scale);
   }   
 
static void SetInitialProfiles(N_Vector cc)
{ 
  N_VConst(0.0, cc);  
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



