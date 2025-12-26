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

static int funcJD(N_Vector cc, N_Vector fval, void *user_data);



static int PrecSetupJD(N_Vector cc, N_Vector cscale,
                       N_Vector fval, N_Vector fscale,
                       void *user_data);

static int PrecSolveJD(N_Vector cc, N_Vector cscale,
                       N_Vector fval, N_Vector fscale,
                       N_Vector vv, void *user_data);


static int jactimes(N_Vector v, N_Vector Jv, N_Vector cc, booleantype *new_u,
                    void *user_data);


/* Private Helper Functions */

static Phase *AllocUserData(void);
static void SetInitialProfilesJD(N_Vector cc);
static realtype SetScaleJD(const struct Phase *const p);
static int check_flag(void *flagvalue, const char *funcname, int opt);
int itersJD;


 /*
 *--------------------------------------------------------------------
 * MAIN ROUTINE
 *--------------------------------------------------------------------
 */

int call_JD(struct Phase *const p)
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
  const soma_scalar_t alfa = p->args.noneq_ratio_arg; // in this case, electrostatic potential difference in kBT/e
  realtype fnorm;
  soma_scalar_t current0, currentL;


  soma_scalar_t *psiC = (soma_scalar_t *) malloc(p->n_cells * sizeof(soma_scalar_t));
    if (psiC == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }


/* Kinsol runs on CPU only, update fields */
#pragma acc update self(p->exp_born_pos[0:p->n_cells])
#pragma acc update self(p->exp_born_neg[0:p->n_cells])

// Calc Born_S

  calc_born_S(p);   

  int NEQ; //<- Number of equations 
  NEQ = (int) p->nx*p->ny*p->nz-1; /* electrostatic potential in all system, fixed reference */

  itersJD = 0; // number of iterations

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

N_VConst(0.0, constraints);  // no constrains c

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
   // initial guess, electrostatic potential equal to equilibrium, so \delta psi = psi' = 0
   // right solution for alfa -> 0 
   //
        for (i = 0 ; i < NEQ ; i++) {
              NVITH(cc,i) = 0.0 ;
        }


   }  else {
       for (i = 0 ; i < NEQ ; i++)  {
              NVITH(cc,i) = ccx[i] ;
	  }
   }


    /* Set scale vector */
    if (flagsolved) scale = SetScaleJD(p);
    N_VConst(scale, sc);

    /* Call KINCreate/KINInit to initialize KINSOL:
       A pointer to KINSOL problem memory is returned and stored in kmem. */
    kmem = KINCreate(sunctx);
    if (check_flag((void *)kmem, "KINCreate", 0)) return(1);

    /* Vector cc passed as template vector. */
    flag = KINInit(kmem, funcJD, cc);
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

      LS = SUNLinSol_SPBCGS(cc, SUN_PREC_NONE, maxl, sunctx);
      if(check_flag((void *)LS, "SUNLinSol_SPBCGS", 0)) return(1); 

//      LS = SUNLinSol_SPBCGS(cc, SUN_PREC_RIGHT, maxl, sunctx);
//      if(check_flag((void *)LS, "SUNLinSol_SPBCGS", 0)) return(1); 

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


//      LS = SUNLinSol_SPTFQMR(cc, SUN_PREC_RIGHT, maxl, sunctx);
//      if(check_flag((void *)LS, "SUNLinSol_SPTFQMR", 0)) return(1);

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
//      LS = SUNLinSol_SPFGMR(cc, SUN_PREC_NONE, maxl, sunctx);
//      if(check_flag((void *)LS, "SUNLinSol_SPFGMR", 0)) return(1); 

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


   /* ------------------------------------
   * Set Jacobian vector product function
   * ------------------------------------ */

    flag = KINSetJacTimesVecFn(kmem, jactimes);
    if (check_flag(&flag, "KINSetJacTimesVecFn", 1)) return(1);

    /* Set preconditioner functions*/
//    flag = KINSetPreconditioner(kmem, PrecSetupJD, PrecSolveJD);
//    if (check_flag(&flag, "KINSetPreconditioner", 1)) return(1);

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

        p->aviter += itersJD;
        p->countiter++;
 
/* Save solution */
        // Save profile  

	for (i = 0 ; i < NEQ ; i++) {
        	ccx[i] = NVITH(cc,i);
                flagsolved = 0;
         } // converged
	   
    } else {  // did not converged
        if (p->info_MPI.sim_rank == 0) 
             fprintf(stdout, "Kinsol failed to converge last step \n");
    }


/* Calc electric field */

// recover electric field from kinsol

  for (i = 0 ; i < NEQ ; i++) {
        psiC[i] = NVITH(cc,i); // note that psi'[NEQ+1] = 0.0 
        p->electric_field[i] += p->electric_field[i] + psiC[i];  
   }

   psiC[p->n_cells-1] = 0.0;

// Calculation of ion currents


//for (iz = 0 ; iz < p->nz-1 ; iz++) { // DEBUG
// iz = 0; // no debug
current0 = 0.0;
  for (ix = 0 ; ix < p->nx ; ix++) {
     for (iy = 0 ; iy < p->ny ; iy++) {

	    cellm = cell_coordinate_to_index(p, ix, iy, iz);
	    cell = cell_coordinate_to_index(p, ix, iy, iz+1);


	    current0 -= (p->npos_field[cell]+p->npos_field[cellm])*(psiC[cell]-psiC[cellm]);
			    
          } // ix
   } //iy
  current0 = current0 * p->deltax*p->deltay/p->deltaz/2.0;
//} // iz -- DEBUG

currentL = 0.0;
  for (ix = 0 ; ix < p->nx ; ix++) {
     for (iy = 0 ; iy < p->ny ; iy++) {

	    cellm = cell_coordinate_to_index(p, ix, iy, iz);
	    cell = cell_coordinate_to_index(p, ix, iy, iz+1);

	    currentL -= (p->npos_field[cell]+p->npos_field[cellm])*(psiC[cell]-psiC[cellm]);
			    
          } // ix
   } //iy
  currentL = currentL * p->deltax*p->deltay/p->deltaz/2.0;

  printf("check: iz, current: %d  %.3e %.3e \n", iz, current0, currentL); // DEBUG
  p->current=current0; // store to save in ana file


// print    
//        printf("Transport converged, flag %d, iters %d, norm %.3e, normtol %.3e, I(0) %.3e \n", flag, itersJD, fnorm, fnormtol, current);

    
	/* Free memory */

// update electric field



    free(psiC);


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

static int funcJD(N_Vector cc, N_Vector fval, void *user_data)
{

#include <assert.h>

  int ix, iy, iz, cell, i;
  int ixp ,ixm, iyp, iym, izp, izm;
  unsigned int iixp ,iixm, iiyp, iiym, iizp, iizm;
  struct Phase *const p = user_data;
  const soma_scalar_t alfa = p->args.noneq_ratio_arg;

  int NEQ; //<- Number of equations 
  NEQ = (int) p->nx*p->ny*p->nz-1; /* the concentration is fixed near electrodes */


  soma_scalar_t *vvin = (soma_scalar_t *) malloc(NEQ * sizeof(soma_scalar_t));
    if (vvin == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }


  soma_scalar_t  res; // residual Poisson Eq.
  soma_scalar_t  psizm, psizp; // auxiliary for PBC

  soma_scalar_t *psiC = (soma_scalar_t *) malloc(p->n_cells * sizeof(soma_scalar_t));
    if (psiC == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }



  itersJD++;	   

// recover difference electrostatic potential, psi' = psi - psieq
  
for (i = 0 ; i < NEQ ; i++) {
        psiC[i] = NVITH(cc,i); 
}
psiC[p->n_cells-1] = 0.0;

soma_scalar_t norma = 0;
// DO NOT PARALELIZE HERE  
  for (ix = 0 ; ix < p->nx ; ix++) {

     ixp = mod((ix+1),p->nx);
     ixm = mod((ix-1),p->nx);
 
     for (iy = 0 ; iy < p->ny ; iy++) {

        iyp = mod((iy+1),p->ny);
        iym = mod((iy-1),p->ny);

#pragma omp parallel for  
	for (iz = 0 ; iz < p->nz ; iz++) {
 
	izp = mod((iz+1),p->nz);
        izm = mod((iz-1),p->nz);

                 i = iz + p->nz*iy + p->nz*p->ny*ix ;
		 iixp = iz + p->nz*iy + p->nz*p->ny*ixp ;
		 iixm = iz + p->nz*iy + p->nz*p->ny*ixm ;
		 iiyp = iz + p->nz*iyp + p->nz*p->ny*ix ;
		 iiym = iz + p->nz*iym + p->nz*p->ny*ix ;
		 iizp = izp + p->nz*iy + p->nz*p->ny*ix ;
		 iizm = izm + p->nz*iy + p->nz*p->ny*ix ;

	psizp = psiC[iizp] + floor((soma_scalar_t)(iz+1)/(soma_scalar_t)p->nz)*alfa; 
	psizm = psiC[iizm] + floor((soma_scalar_t)(iz-1)/(soma_scalar_t)p->nz)*alfa; 
     
	//printf("iz: %d %f %f \n ", c[ix][iy][iz], psi[ix][iy][iz]); 

        res = 0.0;

        res += ((p->npos_field[iixp]+p->npos_field[i])*(psiC[iixp]-psiC[i]))/(p->deltax*p->deltax);
	res += (-(p->npos_field[i]+p->npos_field[iixm])*(psiC[i]-psiC[iixm]))/(p->deltax*p->deltax);

        res += ((p->npos_field[iiyp]+p->npos_field[i])*(psiC[iiyp]-psiC[i]))/(p->deltay*p->deltay);
	res += (-(p->npos_field[i]+p->npos_field[iiym])*(psiC[i]-psiC[iiym]))/(p->deltay*p->deltay);

        res += ((p->npos_field[iizp]+p->npos_field[i])*(psizp-psiC[i]))/(p->deltaz*p->deltaz);
	res += (-(p->npos_field[i]+p->npos_field[iizm])*(psiC[i]-psizm))/(p->deltaz*p->deltaz);

        if (i < NEQ) {NVITH(fval,i) = res;} 
        norma += fabs(res); 
        }
    }
  }

// DEBUG print norm 
        for (ix = 0 ; ix < p->nx ; ix++) {
               for (iy = 0 ; iy < p->ny ; iy++) {
                  for (iz = 0 ; iz <  p->nz ; iz++) {
                  cell = cell_coordinate_to_index(p, ix, iy, iz);
//  printf("check: iz, %.3e %.3e \n", iz,  c[ix][iy][iz], psi[ix][iy][iz]); // DEBUG
                     }
                }
         }
  printf("func: iter, norma: %d %f %f %f %f \n ", itersJD, norma, psiC[0]); 
  
//  printf("func: Nposions, Nnegions: %f, %f \n ", p->Nposions, p->Nnegions);
//  printf("func: Number of Equations: %d \n", NEQ);

//  printf("func: iter, norm %d %.3e \n", iter, norma);

//  exit(1);

  free(psiC);
  return(0);
}
  

/*
 * Set initial conditions in cc
 */


static realtype SetScaleJD(const struct Phase *const p)
{
   realtype scale;

   scale = 1.0 ;
           
   return(scale);
   }

static void SetInitialProfilesJD(N_Vector cc)
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

static int PrecSetupJD(N_Vector cc, N_Vector cscale,
                       N_Vector fval, N_Vector fscale,
                       void *user_data)  {

  unsigned int ix, iy, iz, cell, i;
  unsigned int ixp ,ixm, iyp, iym, izp, izm;
  unsigned int iixp ,iixm, iiyp, iiym, iizp, iizm;
  struct Phase *const p = user_data;
  const soma_scalar_t alfa = p->args.noneq_ratio_arg;

  int NEQ; //<- Number of equations 
  NEQ = (int) p->nx*p->ny*p->nz-1; /* the concentration is fixed near electrodes */

  soma_scalar_t *c = (soma_scalar_t *) malloc(p->n_cells * sizeof(soma_scalar_t));
    if (c == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }


// c from npos_ions
for (i = 0 ; i < p->n_cells ; i++) {
      c[i] =  p->npos_field[i];
}


/// Calculate diagonal preconditioner, temp_prec_field

  for (ix = 0 ; ix < p->nx ; ix++) {

     ixp = mod((ix+1),p->nx);
     ixm = mod((ix-1),p->nx);
 
     for (iy = 0 ; iy < p->ny ; iy++) {

        iyp = mod((iy+1),p->ny);
        iym = mod((iy-1),p->ny);

	for (iz = 0 ; iz < p->nz ; iz++) {
 
	izp = mod((iz+1),p->nz);
        izm = mod((iz-1),p->nz);

		 i = iz + p->nz*iy + p->nz*p->ny*ix ;
		 iixp = iz + p->nz*iy + p->nz*p->ny*ixp ;
		 iixm = iz + p->nz*iy + p->nz*p->ny*ixm ;
		 iiyp = iz + p->nz*iyp + p->nz*p->ny*ix ;
		 iiym = iz + p->nz*iym + p->nz*p->ny*ix ;
		 iizp = izp + p->nz*iy + p->nz*p->ny*ix ;
		 iizm = izm + p->nz*iy + p->nz*p->ny*ix ;

                 p->temp_prec_field[i] = 0.0;
       		 p->temp_prec_field[i] += -(c[iixp]+2*c[i]+c[iixm])/(p->deltax*p->deltax); 
        	 p->temp_prec_field[i] += -(c[iiyp]+2*c[i]+c[iiym])/(p->deltay*p->deltay); 
        	 p->temp_prec_field[i] += -(c[iizp]+2*c[i]+c[iizm])/(p->deltaz*p->deltaz); 
 
	}
      }
    }

free(c);
return(0);

}	

/*
 * Preconditioner solve routine
 */

static int PrecSolveJD(N_Vector cc, N_Vector cscale,
                       N_Vector fval, N_Vector fscale,
                       N_Vector vv, void *user_data)
{
  unsigned int i;
  struct Phase *const p = user_data;
  int NEQ;
  NEQ = (int) p->nx*p->ny*p->nz-1; 

  soma_scalar_t *vvin = (soma_scalar_t *) malloc(NEQ * sizeof(soma_scalar_t));
    if (vvin == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }
  soma_scalar_t *vvout = (soma_scalar_t *) malloc(NEQ * sizeof(soma_scalar_t));
    if (vvout == NULL)
        {
            fprintf(stderr, "ERROR: Malloc %s:%d\n", __FILE__, __LINE__);
            return -1;
        }


  for (i = 0 ; i < NEQ ; i++) {
	  vvin[i] = NVITH(vv,i); 
   }

  for (i = 0 ; i < NEQ ; i++) {
          vvout[i] = vvin[i]/p->temp_prec_field[i]; // Diagonal precond.
	  NVITH(vv,i) = vvout[i]; 
   }

  free(vvin);
  free(vvout);
  return(0);
}


/*
 * Jacobian vector product function
 */

static int jactimes(N_Vector v, N_Vector Jv, N_Vector cc, booleantype *new_u,
                    void *user_data)
{

  unsigned int ix, iy, iz, cell, i;
  unsigned int ixp ,ixm, iyp, iym, izp, izm;
  int iixp ,iixm, iiyp, iiym, iizp, iizm, j;
  struct Phase *const p = user_data;
  const soma_scalar_t alfa = p->args.noneq_ratio_arg;

  int NEQ; //<- Number of equations 
  NEQ = (int) p->nx*p->ny*p->nz-1; /* the concentration is fixed near electrodes */

  soma_scalar_t  tmp;  

  for (ix = 0 ; ix < p->nx ; ix++) {

     ixp = mod((ix+1),p->nx);
     ixm = mod((ix-1),p->nx);
 
     for (iy = 0 ; iy < p->ny ; iy++) {

        iyp = mod((iy+1),p->ny);
        iym = mod((iy-1),p->ny);

	for (iz = 0 ; iz < p->nz ; iz++) {
 
	izp = mod((iz+1),p->nz);
        izm = mod((iz-1),p->nz);

                 i = iz + p->nz*iy + p->nz*p->ny*ix ;


	if (i != NEQ) {

		 iixp = iz + p->nz*iy + p->nz*p->ny*ixp ;
		 iixm = iz + p->nz*iy + p->nz*p->ny*ixm ;
		 iiyp = iz + p->nz*iyp + p->nz*p->ny*ix ;
		 iiym = iz + p->nz*iym + p->nz*p->ny*ix ;
		 iizp = izp + p->nz*iy + p->nz*p->ny*ix ;
		 iizm = izm + p->nz*iy + p->nz*p->ny*ix ;

        // fij for j = i
         tmp  = 0.0;
         tmp  += -(p->npos_field[iixp]+2*p->npos_field[i]+p->npos_field[iixm])/(p->deltax*p->deltax); 
         tmp  += -(p->npos_field[iiyp]+2*p->npos_field[i]+p->npos_field[iiym])/(p->deltay*p->deltay); 
         tmp  += -(p->npos_field[iizp]+2*p->npos_field[i]+p->npos_field[iizm])/(p->deltaz*p->deltaz); 
         j = i; 
	 tmp = tmp*NVITH(v,j);
	 NVITH(Jv,i) = tmp; 

 	 // fij for j = x+1,y,z 
         tmp  = (p->npos_field[iixp]+p->npos_field[i])/(p->deltax*p->deltax); 
         j = iz + p->nz*iy + p->nz*p->ny*ixp ;
	 tmp = tmp*NVITH(v,j);
	 if (j != NEQ) { NVITH(Jv,i) += tmp; } 

        // fij for j = x-1,y,z 
         tmp  = (p->npos_field[iixm]+p->npos_field[i])/(p->deltax*p->deltax); 
         j = iz + p->nz*iy + p->nz*p->ny*ixm ;
	 tmp = tmp*NVITH(v,j);
	 if (j != NEQ) { NVITH(Jv,i) += tmp; } 

        // fij for j = x,y+1,z 
         tmp  = (p->npos_field[iiyp]+p->npos_field[i])/(p->deltay*p->deltay); 
         j = iz + p->nz*iyp + p->nz*p->ny*ix ;
	 tmp = tmp*NVITH(v,j);
	 if (j != NEQ) { NVITH(Jv,i) += tmp; } 

        // fij for j = x,y-1,z 
         tmp  = (p->npos_field[iiym]+p->npos_field[i])/(p->deltay*p->deltay); 
         j = iz + p->nz*iym + p->nz*p->ny*ix ;
	 tmp = tmp*NVITH(v,j);
	 if (j != NEQ) { NVITH(Jv,i) += tmp; } 

        // fij for j = x,y,z+1 
         tmp  = (p->npos_field[iizp]+p->npos_field[i])/(p->deltaz*p->deltaz); 
         j = izp + p->nz*iy + p->nz*p->ny*ix ;
	 tmp = tmp*NVITH(v,j);
	 if (j != NEQ) { NVITH(Jv,i) += tmp; } 

	// fij for j = x,y,z-1 
         tmp  = (p->npos_field[iizm]+p->npos_field[i])/(p->deltaz*p->deltaz); 
         j = izm + p->nz*iy + p->nz*p->ny*ix ;
	 tmp = tmp*NVITH(v,j);
	 if (j != NEQ) { NVITH(Jv,i) += tmp; } 

        } // i =! NEQ

	}
      }
    }

  return(0);
}


