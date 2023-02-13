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

/* static UserData AllocUserData(void);
static void InitUserData(UserData data);
static void FreeUserData(UserData data); */

 static void SetInitialProfiles(N_Vector cc, N_Vector sc);
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
/*  UserData data; */
  int flag, maxl, maxlrst;
  void *kmem;
  SUNLinearSolver LS;
  void *userdata;

  int NEQ; //<- Number of equations 

  NEQ = (int) p->n_cells_local;

  fprintf(stdout, "call_kinsol: Number of cells: %lu \n ", p->n_cells_local);
  fprintf(stdout, "call_kinsol: Number of equations: %d \n ", NEQ);
  fprintf(stdout, "call_kinsol: Bjerrum lenght is: %f \n ", p->Bjerrum);
  fprintf(stdout, "call_kinsol: Nposions, Nnegions: %d, %d \n ", p->Nposions, p->Nnegions);

  userdata = p; // Pointer to phase information

  /* Create the SUNDIALS context object for this simulation. */
  SUNContext sunctx = NULL;
  SUNContext_Create(NULL, &sunctx);

  cc = sc = constraints = NULL;
  kmem = NULL;
  LS = NULL;
/*  data = NULL; */

  /* Allocate memory, and set problem data, initial values, tolerances */
  globalstrategy = KIN_NONE;

/*  data = AllocUserData(); */

/*  if (check_flag((void *)data, "AllocUserData", 2)) return(1);
  InitUserData(data); */

  /* Create serial vectors of length NEQ */
  cc = N_VNew_Serial(NEQ, sunctx);
  if (check_flag((void *)cc, "N_VNew_Serial", 0)) return(1);
  sc = N_VNew_Serial(NEQ, sunctx);
  if (check_flag((void *)sc, "N_VNew_Serial", 0)) return(1);


/* CHECK WHAT NEXT LINE DOES OJO 
  data->rates = N_VNew_Serial(NEQ, sunctx);
  if (check_flag((void *)data->rates, "N_VNew_Serial", 0)) return(1);
*/

  constraints = N_VNew_Serial(NEQ, sunctx);
  if (check_flag((void *)constraints, "N_VNew_Serial", 0)) return(1);
  N_VConst(TWO, constraints);

  fnormtol=FTOL; scsteptol=STOL;

  linsolver = 0; // linear solver, use 0 = SPGMR, 1 = SPBCGS, 2 = SPTFQMR, 3 = SPFGMR

    /* (Re-)Initialize user data */
    SetInitialProfiles(cc, sc);

    /* Call KINCreate/KINInit to initialize KINSOL:
       A pointer to KINSOL problem memory is returned and stored in kmem. */
    kmem = KINCreate(sunctx);
    if (check_flag((void *)kmem, "KINCreate", 0)) return(1);

    /* Vector cc passed as template vector. */
    flag = KINInit(kmem, func, cc);
    if (check_flag(&flag, "KINInit", 1)) return(1);

    flag = KINSetUserData(kmem, userdata);
    if (check_flag(&flag, "KINSetUserData", 1)) return(1);
    flag = KINSetConstraints(kmem, constraints);
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
      printf(" -------");
      printf(" \n| SPGMR |\n");
      printf(" -------\n");

      /* Create SUNLinSol_SPGMR object with right preconditioning and the
         maximum Krylov dimension maxl */
      maxl = 15;
      LS = SUNLinSol_SPGMR(cc, SUN_PREC_RIGHT, maxl, sunctx);
      if(check_flag((void *)LS, "SUNLinSol_SPGMR", 0)) return(1);

      /* Attach the linear solver to KINSOL */
      flag = KINSetLinearSolver(kmem, LS, NULL);
      if (check_flag(&flag, "KINSetLinearSolver", 1)) return 1;

      /* Set the maximum number of restarts */
      maxlrst = 2;
      flag = SUNLinSol_SPGMRSetMaxRestarts(LS, maxlrst);
      if (check_flag(&flag, "SUNLinSol_SPGMRSetMaxRestarts", 1)) return(1);

      break;

    /* (b) SPBCGS */
    case(USE_SPBCGS):

      /* Print header */
      printf(" --------");
      printf(" \n| SPBCGS |\n");
      printf(" --------\n");

      /* Create SUNLinSol_SPBCGS object with right preconditioning and the
         maximum Krylov dimension maxl */
      maxl = 15;
      LS = SUNLinSol_SPBCGS(cc, SUN_PREC_RIGHT, maxl, sunctx);
      if(check_flag((void *)LS, "SUNLinSol_SPBCGS", 0)) return(1);

      /* Attach the linear solver to KINSOL */
      flag = KINSetLinearSolver(kmem, LS, NULL);
      if (check_flag(&flag, "KINSetLinearSolver", 1)) return 1;

      break;

    /* (c) SPTFQMR */
    case(USE_SPTFQMR):

      /* Print header */
      printf(" ---------");
      printf(" \n| SPTFQMR |\n");
      printf(" ---------\n");

      /* Create SUNLinSol_SPTFQMR object with right preconditioning and the
         maximum Krylov dimension maxl */
      maxl = 25;
      LS = SUNLinSol_SPTFQMR(cc, SUN_PREC_RIGHT, maxl, sunctx);
      if(check_flag((void *)LS, "SUNLinSol_SPTFQMR", 0)) return(1);

      /* Attach the linear solver to KINSOL */
      flag = KINSetLinearSolver(kmem, LS, NULL);
      if (check_flag(&flag, "KINSetLinearSolver", 1)) return 1;

      break;

    /* (d) SPFGMR */
    case(USE_SPFGMR):

      /* Print header */
      printf(" -------");
      printf(" \n| SPFGMR |\n");
      printf(" -------\n");

      /* Create SUNLinSol_SPFGMR object with right preconditioning and the
         maximum Krylov dimension maxl */
      maxl = 15;
      LS = SUNLinSol_SPFGMR(cc, SUN_PREC_RIGHT, maxl, sunctx);
      if(check_flag((void *)LS, "SUNLinSol_SPFGMR", 0)) return(1);

      /* Attach the linear solver to KINSOL */
      flag = KINSetLinearSolver(kmem, LS, NULL);
      if (check_flag(&flag, "KINSetLinearSolver", 1)) return 1;

      /* Set the maximum number of restarts */
      maxlrst = 2;
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

    /*
    printf("\n\nComputed equilibrium species concentrations:\n");
    PrintOutput(cc);*/
    
    /* Print final statistics and free memory */
    /*PrintFinalStats(kmem, linsolver); */

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

/*
 * System function for predator-prey system
 */

static int func(N_Vector cc, N_Vector fval, void *user_data)
{
//  realtype xx, yy, delx, dely, *cxy, *rxy, *fxy, dcyli, dcyui, dcxli, dcxri;
//  int jx, jy, is, idyu, idyl, idxr, idxl;


  const struct Phase *const p = user_data;

  fprintf(stdout, "FUNC2: Bjerrum lenght is: %f \n ", p->Bjerrum);
  fprintf(stdout, "FUNC2: Nposions, Nnegions: %d, %d \n ", p->Nposions, p->Nnegions);

  exit(0);

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



  return(0);
}

/*
 * Set initial conditions in cc
 */

static void SetInitialProfiles(N_Vector cc, N_Vector sc)
{
  int i, jx, jy;
  realtype *cloc, *sloc;
  realtype  ctemp[NUM_SPECIES], stemp[NUM_SPECIES];

  /* Initialize arrays ctemp and stemp used in the loading process */
  for (i = 0; i < NUM_SPECIES/2; i++) {
    ctemp[i] = PREYIN;
    stemp[i] = ONE;
  }
  for (i = NUM_SPECIES/2; i < NUM_SPECIES; i++) {
    ctemp[i] = PREDIN;
    stemp[i] = RCONST(0.00001);
  }

  /* Load initial profiles into cc and sc vector from ctemp and stemp. */
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
}

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
