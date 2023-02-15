/* Routines for solving a set of non-linear coupled equations
 * Finds x that satisfy F(x) = 0
 * Based on kinsol  6.5 example "inKrylovDemo_ls.c"  */

#include <stdlib.h>
#include <math.h>

#include <kinsol/kinsol.h>               /* access to KINSOL func., consts.      */
#include <nvector/nvector_serial.h>      /* access to serial N_Vector            */
#include <sunlinsol/sunlinsol_spgmr.h>   /* access to SPGMR SUNLinearSolver      */
#include <sunlinsol/sunlinsol_spbcgs.h>  /* access to SPBCGS SUNLinearSolver     */
#include <sunlinsol/sunlinsol_sptfqmr.h> /* access to SPTFQMR SUNLinearSolver    */
#include <sunlinsol/sunlinsol_spfgmr.h>  /* access to SPFGMR SUNLinearSolver     */
#include <sundials/sundials_dense.h>     /* use generic dense solver in precond. */
#include <sundials/sundials_types.h>     /* defs. of realtype, sunindextype      */

#include "mesh.h"

#include "kinsol_soma.h" 


/* helpful macros */

#ifndef MAX
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

/* Problem Constants */

#define NUM_SPECIES     6  /* must equal 2*(number of prey or predators)
                              number of prey = number of predators       */

#define PI       RCONST(3.1415926535898)   /* pi */

#define MX          5              /* MX = number of x mesh points */
#define MY          5              /* MY = number of y mesh points */
#define NSMX        (NUM_SPECIES * MX)
// #define NEQ         (NSMX * MY)    /* number of equations in the system */
#define AA          RCONST(1.0)    /* value of coefficient AA in above eqns */
#define EE          RCONST(10000.) /* value of coefficient EE in above eqns */
#define GG          RCONST(0.5e-6) /* value of coefficient GG in above eqns */
#define BB          RCONST(1.0)    /* value of coefficient BB in above eqns */
#define DPREY       RCONST(1.0)    /* value of coefficient dprey above */
#define DPRED       RCONST(0.5)    /* value of coefficient dpred above */
#define ALPHA       RCONST(1.0)    /* value of coefficient alpha above */
#define AX          RCONST(1.0)    /* total range of x variable */
#define AY          RCONST(1.0)    /* total range of y variable */
#define FTOL        RCONST(1.e-7)  /* ftol tolerance */
#define STOL        RCONST(1.e-13) /* stol tolerance */
#define THOUSAND    RCONST(1000.0) /* one thousand */
#define ZERO        RCONST(0.)     /* 0. */
#define ONE         RCONST(1.0)    /* 1. */
#define TWO         RCONST(2.0)    /* 2. */
#define PREYIN      RCONST(1.0)    /* initial guess for prey concentrations. */
#define PREDIN      RCONST(30000.0)/* initial guess for predator concs.      */

/* Linear Solver Loop Constants */

#define USE_SPGMR   0
#define USE_SPBCGS  1
#define USE_SPTFQMR 2
#define USE_SPFGMR  3

/* User-defined vector access macro: IJ_Vptr */

/* IJ_Vptr is defined in order to translate from the underlying 3D structure
   of the dependent variable vector to the 1D storage scheme for an N-vector.
   IJ_Vptr(vv,i,j) returns a pointer to the location in vv corresponding to
   indices is = 0, jx = i, jy = j.    */

#define IJ_Vptr(vv,i,j)   (&NV_Ith_S(vv, i*NUM_SPECIES + j*NSMX))

/* Type : UserData
   contains preconditioner blocks, pivot arrays, and problem constants 

typedef struct {
  realtype **P[MX][MY];
  sunindextype *pivot[MX][MY];
  realtype **acoef, *bcoef;
  N_Vector rates;
  realtype *cox, *coy;
  realtype ax, ay, dx, dy;
  realtype uround, sqruround;
  int mx, my, ns, np;
} *UserData;  */

/* Functions Called by the KINSOL Solver */

static int func(N_Vector cc, N_Vector fval, void *user_data);

/*static int PrecSetupBD(N_Vector cc, N_Vector cscale,
                       N_Vector fval, N_Vector fscale,
                       void *user_data);

static int PrecSolveBD(N_Vector cc, N_Vector cscale,
                       N_Vector fval, N_Vector fscale,
                       N_Vector vv, void *user_data);
*/
/* Private Helper Functions */

static Phase *AllocUserData(void);
/*static void InitUserData(UserData data);
static void FreeUserData(UserData data); */

 static void SetInitialProfiles(N_Vector cc, N_Vector sc, int NEQ, const struct Phase *const p);
/* static void PrintHeader(int globalstrategy, int maxl, int maxlrst,
                        realtype fnormtol, realtype scsteptol,
			int linsolver);
*/
// static void PrintOutput(N_Vector cc);
// static void PrintFinalStats(void *kmem, int linsolver);

/*
static void WebRate(realtype xx, realtype yy, realtype *cxy, realtype *ratesxy,
                    void *user_data);
*/
// static realtype DotProd(int size, realtype *x1, realtype *x2);
 static int check_flag(void *flagvalue, const char *funcname, int opt);



 soma_scalar_t scale;
 int iter = 0;
 soma_scalar_t norma;  // sum of residuals

 
 /*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

int call_kinsol(const struct Phase *const p)
{
  int globalstrategy, linsolver;
  realtype fnormtol, scsteptol;
  N_Vector cc, sc, constraints;
  int flag, maxl, maxlrst;
  void *kmem;
  SUNLinearSolver LS;
  Phase *data;


  int NEQ; //<- Number of equations 
  NEQ = (int) p->n_cells_local - 1; /* Due to PBC the set of equations is no longer LI, so phi(nx,ny,nz) can
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
  globalstrategy = KIN_NONE; /* KIN_NONE = basic Newton iteration
				KIN_LINESEARCH = Newton with globalization
				KIN_FP = fixed point interaction
				KIN_PICARD = Picard interaction */

  data = AllocUserData(); 
  if (check_flag((void *)data, "AllocUserData", 2)) return(1);

  *data = *p; // Pointer to phase information

//  InitUserData(data); */

  /* Create serial vectors of length NEQ */
  cc = N_VNew_Serial(NEQ, sunctx);
  if (check_flag((void *)cc, "N_VNew_Serial", 0)) return(1);
  sc = N_VNew_Serial(NEQ, sunctx);
  if (check_flag((void *)sc, "N_VNew_Serial", 0)) return(1);


/* CHECK WHAT NEXT LINE DOES OJO 
  data->rates = N_VNew_Serial(NEQ, sunctx);
  if (check_flag((void *)data->rates, "N_VNew_Serial", 0)) return(1);
*/

// constraints = N_VNew_Serial(NEQ, sunctx);
//  if (check_flag((void *)constraints, "N_VNew_Serial", 0)) return(1);
//  N_VConst(ZERO, constraints);

  fnormtol=FTOL; scsteptol=STOL;

  linsolver = 1; // linear solver, use 0 = SPGMR, 1 = SPBCGS, 2 = SPTFQMR, 3 = SPFGMR

    /* (Re-)Initialize user data */
    SetInitialProfiles(cc, sc, NEQ, p);

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

      /* Create SUNLinSol_SPBCGS object with right preconditioning and the
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

    /* Set preconditioner functions */
/*    flag = KINSetPreconditioner(kmem, PrecSetupBD, PrecSolveBD);
    if (check_flag(&flag, "KINSetPreconditioner", 1)) return(1);
*/
    
    /* Print out the problem size, solution parameters, initial guess. */
    /*PrintHeader(globalstrategy, maxl, maxlrst, fnormtol, scsteptol, linsolver);*/
    


    /* Call KINSol and print output concentration profile */

    flag = KINSol(kmem,           /* KINSol memory block */
		  cc,             /* initial guess on input; solution vector */
		  globalstrategy, /* global strategy choice */
		  sc,             /* scaling vector, for the variable cc */
		  sc);            /* scaling vector for function values fval */
    if (check_flag(&flag, "KINSol", 1)) return(1);

    printf("Electrostatic converged in %d iters, with norm %.3e \n", iter, norma*scale);
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



/* Readability definitions used in other routines below */
/*
#define acoef  (data->acoef)
#define bcoef  (data->bcoef)
#define cox    (data->cox)
#define coy    (data->coy)
*/
/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY KINSOL
 *--------------------------------------------------------------------
 */

static int func(N_Vector cc, N_Vector fval, void *user_data)
{
//  realtype xx, yy, delx, dely, *cxy, *rxy, *fxy, dcyli, dcyui, dcxli, dcxri;
//  int jx, jy, is, idyu, idyl, idxr, idxl;

  int ix, iy, iz, cell;	
  int ixp ,ixm, iyp, iym, izp, izm;
  const struct Phase *const p = user_data;

  int NEQ = (int) p->n_cells_local - 1; /* Due to PBC the set of equations is no longer LI, so phi(nx,ny,nz) can
				       be fixed to zero (see notes) */


  soma_scalar_t  rhoA[p->nx][p->ny][p->nz]; // number of A segments in 3D lattice 
  soma_scalar_t  sumrhoA = 0;                   // total number of A segments 
  soma_scalar_t  deltax = p->Lx/((soma_scalar_t) p->nx);
  soma_scalar_t  deltay = p->Ly/((soma_scalar_t) p->ny);
  soma_scalar_t  deltaz = p->Lz/((soma_scalar_t) p->nz);
  soma_scalar_t  res[p->nx][p->ny][p->nz]; // residual Poisson Eq.
  soma_scalar_t  rhoposion[p->nx][p->ny][p->nz]; // number of positive ions in 3D lattice
  soma_scalar_t  rhonegion[p->nx][p->ny][p->nz]; // number of negative ions in 3D lattice
  soma_scalar_t  rhoQ[p->nx][p->ny][p->nz]; // total charge density
  soma_scalar_t  phi[p->nx][p->ny][p->nz]; // electrostatic potential 3D lattice, units of kBT/|e|
  soma_scalar_t  sumposions = 0 , sumnegions = 0 ;
  soma_scalar_t  sumrhoQ = 0 ; // sum of rhoQ, for debug only
  
  soma_scalar_t constq = 4.0*PI*p->Bjerrum/(deltax*deltay*deltaz); // multiplicative constant for Poisson equation

  iter++;	   

// phi from kinsol's input

  for (ix = 0 ; ix < (int) p->nx ; ix++) {
	  for (iy = 0 ; iy < (int) p->ny ; iy++) {
			  for (iz = 0 ; iz < (int)  p->nz ; iz++) {
                          cell = cell_coordinate_to_index(p, ix, iy, iz);
                          if (cell < NEQ) phi[ix][iy][iz] = NV_Ith_S(cc,cell); 
		     }
	 	}
	  }

phi[p->nx-1][p->ny-1][p->nz-1] = 0.0; // choice of zero of electrostatic potential due to PBC

// Pos ion and neg ion

  for (ix = 0 ; ix < (int) p->nx ; ix++) {
	  for (iy = 0 ; iy < (int) p->ny ; iy++) {
			  for (iz = 0 ; iz < (int)  p->nz ; iz++) {
                          rhoposion[ix][iy][iz] = exp(-phi[ix][iy][iz]);
			  sumposions += rhoposion[ix][iy][iz]; 

			  rhonegion[ix][iy][iz] = exp(phi[ix][iy][iz]);
			  sumnegions += rhonegion[ix][iy][iz]; 

		     }
	 	}
	  }

// Normalize to match totalnumber

  for (ix = 0 ; ix < (int) p->nx ; ix++) {
	  for (iy = 0 ; iy < (int) p->ny ; iy++) {
			  for (iz = 0 ; iz < (int)  p->nz ; iz++) {
                          rhoposion[ix][iy][iz] = rhoposion[ix][iy][iz] * p->Nposions / sumposions; 
                          rhonegion[ix][iy][iz] = rhonegion[ix][iy][iz] * p->Nnegions / sumnegions; 

		     }
	 	}
	  }

// rhoA from unified fields 
  for (ix = 0 ; ix < (int) p->nx ; ix++) {
	  for (iy = 0 ; iy < (int) p->ny ; iy++) {
			  for (iz = 0 ; iz < (int)  p->nz ; iz++) {
                          cell = cell_coordinate_to_index(p, ix, iy, iz);
                          rhoA[ix][iy][iz] = p->fields_unified[cell]; // density of A segments because no n_type offset is used                           
	                  sumrhoA += rhoA[ix][iy][iz];
//                          printf("ix, iy, iz, cell, %d, %d, %d, %d, %f \n", ix, iy, iz, cell, rhoA[ix][iy][iz]);
		     }
	 	}
	  }


// rhoQ from all contributions 
  for (ix = 0 ; ix < (int) p->nx ; ix++) {
	  for (iy = 0 ; iy < (int) p->ny ; iy++) {
			  for (iz = 0 ; iz < (int)  p->nz ; iz++) {
                          rhoQ[ix][iy][iz] = rhoA[ix][iy][iz] * p->Acharge; 
                          rhoQ[ix][iy][iz] += rhoposion[ix][iy][iz]-rhonegion[ix][iy][iz]; 
			  sumrhoQ += rhoQ[ix][iy][iz];	  
		     }
	 	}
	  }


/// Calculate residual from Poisson's equation

  norma = 0.0; 

  for (ix = 0 ; ix < (int) p->nx ; ix++) {

	  ixp = (ix+1) % (p->nx);
          ixm = (ix-1) % (p->nx);
 
	  for (iy = 0 ; iy < (int) p->ny ; iy++) {

		  iyp = (iy+1)%(p->ny);
		  iym = (iy-1)%(p->ny);

	          for (iz = 0 ; iz < (int)  p->nz ; iz++) {
        
			  izp = (iz+1)%(p->nz);
			  izm = (iz-1)%(p->nz);
                  	  cell = cell_coordinate_to_index(p, ix, iy, iz);

			  res[ix][iy][iz] = rhoQ[ix][iy][iz]*constq;

			  res[ix][iy][iz] += (phi[ixp][iy][iz]-2.*phi[ix][iy][iz]+phi[ixm][iy][iz])/(deltax*deltax);	  
			  res[ix][iy][iz] += (phi[ix][iyp][iz]-2.*phi[ix][iy][iz]+phi[ix][iym][iz])/(deltay*deltay);	  
			  res[ix][iy][iz] += (phi[ix][iy][izp]-2.*phi[ix][iy][iz]+phi[ix][iy][izm])/(deltaz*deltaz);	  
		  
			  if (cell < NEQ) { 
				  NV_Ith_S(fval,cell) = -res[ix][iy][iz];
				  norma += fabs(-res[ix][iy][iz]);
//				  printf("func: cell, res %d %f \n", cell, rhoA[ix][iy][iz]);
			  }     
		     }
	 	}
	  }

//  printf("func: sumrhoA: %f \n ", sumrhoA);
//  printf("func: sumrhoQ: %f \n ", sumrhoQ);
//  printf("func: Bjerrum lenght is: %f \n ", p->Bjerrum);
//  printf("func: Nposions, Nnegions: %f, %f \n ", p->Nposions, p->Nnegions);
//  printf("func: Number of Equations: %d \n", NEQ);

//  printf("func: iter, norm %d %.3e \n", iter, norma);

//  exit(1);
  return(0);
}
  

/*  data = (UserData)user_data;
  delx = data->dx;
  dely = data->dy; */



  /* Loop over all mesh points, evaluating rate array at each point*/
//  for (jy = 0; jy < MY; jy++) {

//    yy = dely*jy;

    /* Set lower/upper index shifts, special at boundaries. */
/*    idyl = (jy != 0   ) ? NSMX : -NSMX;
    idyu = (jy != MY-1) ? NSMX : -NSMX;

    for (jx = 0; jx < MX; jx++) {

      xx = delx*jx;
*/
      /* Set left/right index shifts, special at boundaries. */
/*      idxl = (jx !=  0  ) ?  NUM_SPECIES : -NUM_SPECIES;
      idxr = (jx != MX-1) ?  NUM_SPECIES : -NUM_SPECIES;

      cxy = IJ_Vptr(cc,jx,jy);
      rxy = IJ_Vptr(data->rates,jx,jy);
      fxy = IJ_Vptr(fval,jx,jy);
*/
      /* Get species interaction rate array at (xx,yy) */
/*      WebRate(xx, yy, cxy, rxy, user_data);

      for(is = 0; is < NUM_SPECIES; is++) {
*/
        /* Differencing in x direction */
/*        dcyli = *(cxy+is) - *(cxy - idyl + is) ;
        dcyui = *(cxy + idyu + is) - *(cxy+is);
*/
        /* Differencing in y direction */
/*        dcxli = *(cxy+is) - *(cxy - idxl + is);
        dcxri = *(cxy + idxr +is) - *(cxy+is);
*/
        /* Compute the total rate value at (xx,yy) */
/*        fxy[is] = (coy)[is] * (dcyui - dcyli) +
          (cox)[is] * (dcxri - dcxli) + rxy[is];
*/
//      } /* end of is loop */

//    } /* end of jx loop */

//  } /* end of jy loop */



//}

/*
 * Set initial conditions in cc
 */

static void SetInitialProfiles(N_Vector cc, N_Vector sc, int NEQ, const struct Phase *const p)
{ 
   int i, ix, iy, iz;
   int cell; 
   soma_scalar_t deltax = p->Lx/((soma_scalar_t) p->nx);
   soma_scalar_t deltay = p->Ly/((soma_scalar_t) p->ny);
   soma_scalar_t deltaz = p->Lz/((soma_scalar_t) p->nz);
   soma_scalar_t constq = 4.0*PI*p->Bjerrum/(deltax*deltay*deltaz); // multiplicative constant for Poisson equation
   soma_scalar_t  sumrhoA = 0;                   // total number of A segments

   for (ix = 0 ; ix < (int) p->nx ; ix++) {
          for (iy = 0 ; iy < (int) p->ny ; iy++) {
                          for (iz = 0 ; iz < (int)  p->nz ; iz++) {
                          cell = cell_coordinate_to_index(p, ix, iy, iz);
                          sumrhoA += p->fields_unified[cell]; // density of A segments because no n_type offset is used                           
     			  }
                }
          }


   scale = 1./constq/sumrhoA*((soma_scalar_t) p->n_cells_local);


// Initial guess for electrostatic potential is phi = 0 everywhere

for (i = 0 ; i < NEQ ; i++) {
   NV_Ith_S(cc,i) = 0.0; // Initial Guess
   NV_Ith_S(sc,i) = scale;   }
}
   
   /*static void SetInitialProfiles(N_Vector cc, N_Vector sc)
{
  int i, jx, jy;
  realtype *cloc, *sloc;
  realtype  ctemp[NUM_SPECIES], stemp[NUM_SPECIES];

  printf("ENTER SET INITIAL PROFILES \n");

  for (i = 0; i < NUM_SPECIES/2; i++) {
    ctemp[i] = PREYIN;
    stemp[i] = ONE;
  }
  for (i = NUM_SPECIES/2; i < NUM_SPECIES; i++) {
    ctemp[i] = PREDIN;
    stemp[i] = RCONST(0.00001);
  }

  for (jy = 0; jy < MY; jy++) {
    for (jx = 0; jx < MX; jx++) {
      cloc = IJ_Vptr(cc,jx,jy);
      sloc = IJ_Vptr(sc,jx,jy);
      for (i = 0; i < NUM_SPECIES; i++) {
        cloc[i] = ctemp[i];
        sloc[i] = stemp[i];
      }
    }
  }
  for (i = 0 ; i <= NSMX * MY ; i++) {
	  printf("Value %d ,  cc = %f \n", i, cc[i]); 
	  printf("Value %d ,  cloc = %f \n", i, cloc[i]); 
  }
}
*/
/*
 * Print first lines of output (problem description)
 */
/*
static void PrintHeader(int globalstrategy, int maxl, int maxlrst,
                        realtype fnormtol, realtype scsteptol,
			int linsolver)
{
  printf("\nPredator-prey test problem --  KINSol (serial version)\n\n");
  printf("Mesh dimensions = %d X %d\n", MX, MY);
  printf("Number of species = %d\n", NUM_SPECIES);
  printf("Total system size = %d\n\n", NEQ);
  printf("Flag globalstrategy = %d (0 = None, 1 = Linesearch)\n",
         globalstrategy);

  switch(linsolver) {

  case(USE_SPGMR):
    printf("Linear solver is SPGMR with maxl = %d, maxlrst = %d\n",
	   maxl, maxlrst);
    break;

  case(USE_SPBCGS):
    printf("Linear solver is SPBCGS with maxl = %d\n", maxl);
    break;

  case(USE_SPTFQMR):
    printf("Linear solver is SPTFQMR with maxl = %d\n", maxl);
    break;

  case(USE_SPFGMR):
    printf("Linear solver is SPFGMR with maxl = %d, maxlrst = %d\n",
	   maxl, maxlrst);
    break;

  }

  printf("Preconditioning uses interaction-only block-diagonal matrix\n");
  printf("Positivity constraints imposed on all components \n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("Tolerance parameters:  fnormtol = %Lg   scsteptol = %Lg\n",
         fnormtol, scsteptol);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("Tolerance parameters:  fnormtol = %g   scsteptol = %g\n",
         fnormtol, scsteptol);
#else
  printf("Tolerance parameters:  fnormtol = %g   scsteptol = %g\n",
         fnormtol, scsteptol);
#endif

  printf("\nInitial profile of concentration\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At all mesh points:  %Lg %Lg %Lg   %Lg %Lg %Lg\n",
         PREYIN, PREYIN, PREYIN,
         PREDIN, PREDIN, PREDIN);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At all mesh points:  %g %g %g   %g %g %g\n",
         PREYIN, PREYIN, PREYIN,
         PREDIN, PREDIN, PREDIN);
#else
  printf("At all mesh points:  %g %g %g   %g %g %g\n",
         PREYIN, PREYIN, PREYIN,
         PREDIN, PREDIN, PREDIN);
#endif
}
*/
/*
 * Print sampled values of current cc
 */
/*
static void PrintOutput(N_Vector cc)
{
  int is, jx, jy;
  realtype *ct;

  jy = 0; jx = 0;
  ct = IJ_Vptr(cc,jx,jy);
  printf("\nAt bottom left:");
*/
  /* Print out lines with up to 6 values per line */
/*
for (is = 0; is < NUM_SPECIES; is++){
    if ((is%6)*6 == is) printf("\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf(" %Lg",ct[is]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf(" %g",ct[is]);
#else
    printf(" %g",ct[is]);
#endif
  }

  jy = MY-1; jx = MX-1;
  ct = IJ_Vptr(cc,jx,jy);
  printf("\n\nAt top right:");
*/
  /* Print out lines with up to 6 values per line */
/*  for (is = 0; is < NUM_SPECIES; is++) {
    if ((is%6)*6 == is) printf("\n");
#if defined(SUNDIALS_EXTENDED_PRECISION)
    printf(" %Lg",ct[is]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
    printf(" %g",ct[is]);
#else
    printf(" %g",ct[is]);
#endif
  }
  printf("\n\n");
}
*/
/*
 * Print final statistics contained in iopt
 */
/*
static void PrintFinalStats(void *kmem, int linsolver)
{
  long int nni, nfe, nli, npe, nps, ncfl, nfeSG;
  int flag;

  flag = KINGetNumNonlinSolvIters(kmem, &nni);
  check_flag(&flag, "KINGetNumNonlinSolvIters", 1);
  flag = KINGetNumFuncEvals(kmem, &nfe);
  check_flag(&flag, "KINGetNumFuncEvals", 1);

  flag = KINGetNumLinIters(kmem, &nli);
  check_flag(&flag, "KINGetNumLinIters", 1);
  flag = KINGetNumPrecEvals(kmem, &npe);
  check_flag(&flag, "KINGetNumPrecEvals", 1);
  flag = KINGetNumPrecSolves(kmem, &nps);
  check_flag(&flag, "KINGetNumPrecSolves", 1);
  flag = KINGetNumLinConvFails(kmem, &ncfl);
  check_flag(&flag, "KINGetNumLinConvFails", 1);
  flag = KINGetNumLinFuncEvals(kmem, &nfeSG);
  check_flag(&flag, "KINGetNumLinFuncEvals", 1);

  printf("Final Statistics.. \n");
  printf("nni    = %5ld    nli   = %5ld\n", nni, nli);
  printf("nfe    = %5ld    nfeSG = %5ld\n", nfe, nfeSG);
  printf("nps    = %5ld    npe   = %5ld     ncfl  = %5ld\n", nps, npe, ncfl);

  if (linsolver < 3) printf("\n=========================================================\n\n");

}
*/
/*
 * Check function return value...
 *    opt == 0 means SUNDIALS function allocates memory so check if
 *             returned NULL pointer
 *    opt == 1 means SUNDIALS function returns a flag so check if
 *             flag >= 0
 *    opt == 2 means function allocates memory so check if returned
 *             NULL pointer
 */

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

