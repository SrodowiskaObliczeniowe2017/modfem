/************************************************************************
File pdh_conv_diff.h - internal information for problem dependent module
                     for solving convection-diffusion equations

Contains declarations of constants and interface routines:

  pdr_time_integration - time integration driver (in time integration) 
  pdr_adapt - to enforce adaptation strategy for a given problem (in adapt)

in pds_conv_diff_util.c :
  pdr_select_problem - to select the proper problem   
  pdr_bc_val - to get the boundary condition values
  pdr_ic_val - to get the initial condition values
  pdr_time_i_params - to return parameters of time discretization
  pdr_time_d_params - to return parameters of time discretization
  pdr_set_time_i_params - to change parameters of time discretization
  pdr_set_time_d_params - to change parameters of time discretization
  pdr_lins_i_params - to return parameters of linear equations solver
  pdr_lins_d_params - to return parameters of linear equations solver
  pdr_conv_diff_error_test - to compute error norm for test examples
  pdr_conv_diff_post_process - to privide simple interactive post-processing
  pdr_conv_diff_average_sol - to compute the average of solution over element

  pdr_get_mat_data - to get material data from material structures

  pdr_exact_sol - to return values and derivatives at a point
	for functions used as exact solutions for test problems

  pdr_slope_limit - to limit the slope of linear solution

  pdr_zzhu_error - to get Zienkiewicz-Zhu error estimator
  pdr_zzhu_error_small - to get Zienkiewicz-Zhu error estimator based on small patches

// SHOULD BE CHANGED ?
  pdr_fa_coeff - to return coefficients for face integrals
  pdr_get_bc_type - to get BC type given BC flag from mesh data structure
                    !!! according to some adopted convention !!!
  pdr_pde_coeff - to return coefficients of the original pdes
  pdr_bc_diri_coeff - to return DIRICHLET boundary coeficients
                     (routine substitutes only non-zero values)
  pdr_bc_neum_coeff - to return NEUMANN boundary coeficients
                     (routine substitutes only non-zero values)
  pdr_bc_mixed_coeff - to return ROBIN (mixed) boundary coeficients

------------------------------  			
History:
	02.2002 - Krzysztof Banas, initial version		
*************************************************************************/

#ifndef _pdh_conv_diff_
#define _pdh_conv_diff_

#include <stdio.h>

#ifdef __cplusplus
extern "C" 
{
#endif

/* USED AND IMPLEMENTED PROBLEM DEPENDENT DATA STRUCTURES AND INTERFACES */

#include "../include/pdh_conv_diff_problem.h"	
/* types and functions related to boundary conditions handling */
  //#include "../include/pdh_conv_diff_bc.h"	
/* types and functions related to materials handling */
  //#include "../include/pdh_conv_diff_materials.h"	
  //#include "../include/pdh_conv_diff_weakform.h"	

/*** CONSTANTS ***/


/* GLOBAL VARIABLES */
extern int pdv_conv_diff_nr_problems;  /* the number of problems  */
extern int pdv_conv_diff_current_problem_id;  /* ID of the current problem */
extern pdt_conv_diff_problem pdv_conv_diff_problems[PDC_CONV_DIFF_MAX_NUM_PROB];
                                                       /* array of problems */
  /* pdv_conv_diff_current_problem_id < pdv_conv_diff_nr_problems */
  /* pdv_conv_diff_nr_problems <= PDC_CONV_DIFF_MAX_NUM_PROB */

#ifdef PARALLEL
extern int pdv_exchange_table_index;
#endif


/*** INTERNAL PROCEDURES ***/

/**-----------------------------------------------------------
pdr_time_integration - time integration driver. 
------------------------------------------------------------*/  
extern void pdr_time_integration(
  int Problem_id,       /* in: problem data structure to be used */
  char* Work_dir,
  FILE *interactive_input, 
  FILE *interactive_output
);


/**--------------------------------------------------------
pdr_adapt - to enforce adaptation strategy for a given problem 
---------------------------------------------------------*/
extern int pdr_adapt( /* returns: >0 - success, <=0 - failure */
  int Problem_id,       /* in: problem data structure to be used */
  char* Work_dir,
  FILE *Interactive_input, 
  FILE *Interactive_output
  );

/**--------------------------------------------------------
pdr_zzhu_error - to get Zienkiewicz-Zhu error estimator
----------------------------------------------------------*/
double pdr_zzhu_error(
/* returns  - Zienkiewicz-Zhu error for the whole mesh */
  int Field_id,    /* in: approximation field ID  */
  FILE *interactive_output
);

/**--------------------------------------------------------
pdr_zzhu_error_small - to get Zienkiewicz-Zhu error based on small patches
----------------------------------------------------------*/
double pdr_zzhu_error_small( 
/* returns  - Zienkiewicz-Zhu error for the whole mesh */
  int Field_id,    /* in: approximation field ID  */
  FILE *interactive_output
			     );

/**--------------------------------------------------------
  pdr_select_problem - to select the proper problem   
---------------------------------------------------------*/
extern pdt_conv_diff_problem* pdr_select_problem( 
                                 /* returns pointer to the chosen problem */
			         /* (if input is not valid it returns */
		                 /* the pointer to the current problem) */
  int Problem_id  /* in: problem ID or 
		     0 (PDC_CURRENT_PROBLEM_ID) for the current problem */
  );


/**-------------------------------------------------------------------------
  pdr_bc_val - to get the boundary condition values
---------------------------------------------------------------------------*/
extern int pdr_bc_val( /* returns: 1 - success, <=0 - failure */
	int Problem_id,	/* in: data structure to be used */
	int Num,       	/* in: index in the array of boundary values */
	double *Val     /* out: vector of boundary values */
	);

/**-------------------------------------------------------------------------
  pdr_ic_val - to get the initial condition values
---------------------------------------------------------------------------*/
extern int pdr_ic_val( /* returns: 1 - success, <=0 - failure */
	int Problem_id,	/* in: data structure to be used */
	int Num,       	/* in: index in the array of initial values */
	double *Val     /* out: vector of initial values */
	);

/**--------------------------------------------------------
  pdr_time_i_params - to return parameters of time discretization
---------------------------------------------------------*/
extern int pdr_time_i_params( /* returns: integer time integration parameter */
	int Problem_id,	/* in: data structure to be used  */
	int Num         /* in: parameter number in control structure */
	);

/**--------------------------------------------------------
  pdr_time_d_params - to return parameters of time discretization
---------------------------------------------------------*/
extern double pdr_time_d_params( /* returns: real time integration parameter */
	int Problem_id,	/* in: data structure to be used  */
	int Num         /* in: parameter number in control structure */
	);

/**--------------------------------------------------------
  pdr_set_time_i_params - to change parameters of time discretization
---------------------------------------------------------*/
extern void pdr_set_time_i_params( 
        int Problem_id,	     /* in: data structure to be used  */
	int Num,             /* in: parameter number in control structure */
	int Value            /* in: parameter value */
	);

/**--------------------------------------------------------
  pdr_set_time_d_params - to change parameters of time discretization
---------------------------------------------------------*/
extern void pdr_set_time_d_params( 
        int Problem_id,	     /* in: data structure to be used  */
	int Num,             /* in: parameter number in control structure */
	double Value         /* in: parameter value */
	);

/**--------------------------------------------------------
  pdr_lins_i_params - to return parameters of linear equations solver
---------------------------------------------------------*/
extern int pdr_lins_i_params( /* returns: integer linear solver parameter */
	int Problem_id,	/* in: data structure to be used  */
	int Num         /* in: parameter number in control structure */
	);

/**--------------------------------------------------------
  pdr_lins_d_params - to return parameters of linear equations solver
---------------------------------------------------------*/
extern double pdr_lins_d_params( /* returns: real linear solver parameter */
	int Problem_id,	/* in: data structure to be used  */
	int Num         /* in: parameter number in control structure */
	);


/**--------------------------------------------------------
  pdr_conv_diff_error_test - to compute error norm for test examples
----------------------------------------------------------*/
double pdr_conv_diff_error_test( /* returns H1 norm of error for */
	         	/* for several known exact solutions */
  int Field_id,    /* in: approximation field ID  */
  FILE *Interactive_output
  );

/**--------------------------------------------------------
  pdr_conv_diff_post_process - to privide simple interactive post-processing
----------------------------------------------------------*/
double pdr_conv_diff_post_process(
  int Field_id,    /* in: approximation field ID  */
  FILE *Interactive_input,
  FILE *Interactive_output
  );

/**--------------------------------------------------------
pdr_average_sol_el - to compute the average of solution over element
----------------------------------------------------------*/
double pdr_average_sol_el( /* returns: the average of solution over element */
        int Problem_id,      /* in: data structure to be used  */
        int El          /* in: element number */
				  );


/**--------------------------------------------------------
pdr_conv_diff_exact_sol - to return values and derivatives at a point
	for functions used as exact solutions for test problems
----------------------------------------------------------*/
extern int pdr_conv_diff_exact_sol(   /* returns: 1-found exact solution, 0-not found */
	int Flag,    /* in: flag, e.g. material number */
	double X,	/* in: coordinates of a point */
	double Y,
	double Z,
	double Time,	/* in: time instant */
	double *Exact,	/* out: values of solution and derivatives */
	double *Exact_x,
	double *Exact_y,
	double *Exact_z,
	double *Lapl	/* out: Laplacian of function */
	);

/**--------------------------------------------------------
pdr_slope_limit - to limit the slope of linear solution
---------------------------------------------------------*/
int pdr_slope_limit( /* returns: >0 - success, <=0 - failure */
	int Problem_id	/* in: data structure to be used  */
	);

/**--------------------------------------------------------
pdr_fa_coeff - to return coefficients for face integrals
----------------------------------------------------------*/
void pdr_fa_coeff(
	int Face,	/* in: face number */
        int BC_flag,	/* in: boundary condition flag */
	int Elem,	/* in: element number */
	int Mat_num,	/* in: material number */
	double Hsize,	/* in: size of an element */
	double Time,    /* in: time instant */
	int Pdeg,	/* in: local degree of polynomial */
	double* Xcoor,	/* in: global coordinates of a point */
	double* Vec_norm, /* in: unit normal to the boundary */
	double* U_val, 	/* in: computed solution */
	double* U_x, 	/* in: gradient of computed solution */
	double* U_y,   	/* in: gradient of computed solution */
	double* U_z,   	/* in: gradient of computed solution */
	double* Anx,	/* out: diffusion coefficients */
	double* Any, 	
	double* Anz, 	
	double* Bn,	/* out: convection coefficients */
	double* Fval,	/* out: rhs  Dirichlet data f */
	double* Gval,	/* out: rhs  Neumann data g */
	double* Qn,	/* out: rhs normal fluxes */
	double* Vel_norm /* out: velocity normal to the boundary */
	);

/**-----------------------------------------------------------
  pdr_get_bc_type - to get BC type given BC flag from mesh data structure
                    !!! according to some adopted convention !!!
------------------------------------------------------------*/
int pdr_get_bc_type( /* returns: bc type for a face: */
		/*  	1 (PDC_INTERIOR) - interior face */
		/*  	2 (PDC_BC_DIRI) - Dirichlet boundary face */
		/*  	3 (PDC_BC_NEUM) - Neumann boundary face */
		/*  	4 (PDC_BC_MIXED) - Robin boundary face */
   int Fa_bc	/* in: BC flag (returned by mesh data structure) */
   );

/**--------------------------------------------------------
  pdr_pde_coeff - to return coefficients of the original pdes
----------------------------------------------------------*/
void pdr_pde_coeff(
	int Mat_num,	/* in: material number */
	double* Xcoor,	/* in: global coordinates of a point */
	double  Time,	/* in: time instant */
	double* U_val, 	/* in: computed solution */
	double* U_x, 	/* in: gradient of computed solution */
	double* U_y,   	/* in: gradient of computed solution */
	double* U_z,   	/* in: gradient of computed solution */
	double* Axx, 	/* out: diffusion coefficients */		
	double* Axy,		
	double* Axz,		
	double* Ayx,
	double* Ayy,
	double* Ayz,
	double* Azx,
	double* Azy,
	double* Azz,
	double* Bx,	/* out: convection coefficients */	
	double* By,	
	double* Bz,	
	double* Tx,	/* out: convection coefficients */	
	double* Ty,	
	double* Tz,	
	double* Cval,	/* out: reaction coefficients */
	double* Mval,	/* out: mass matrix coefficients */
	double* Lval,	/* out: rhs time coefficients */
	double* Sval,	/* out: rhs coefficients without derivatives */	
	double* Qx,	/* out: rhs coefficients with derivatives */	
	double* Qy,	/* out: rhs coefficients with derivatives */	
	double* Qz	/* out: rhs coefficients with derivatives */	
	);

/**--------------------------------------------------------
pdr_bc_diri_coeff - to return DIRICHLET boundary coeficients
                     (routine substitutes only non-zero values)
----------------------------------------------------------*/
void pdr_bc_diri_coeff(
	int Face,       /* in: index in an array of BC coefficients */
	int Mat_num,	/* in: material number */
	double* Xcoor,	/* in: global coordinates of a point */
	double* Vec_norm,/* in: unit normal to the boundary */
	double  Time,	/* in: time instant */
	double* U_val, 	/* in: computed solution */
	double* U_x, 	/* in: gradient of computed solution */
	double* U_y,   	/* in: gradient of computed solution */
	double* U_z,   	/* in: gradient of computed solution */
	double* Fval	/* out: BC DIRICHLET coefficients */	
	);

/**--------------------------------------------------------
pdr_bc_neum_coeff - to return NEUMANN boundary coeficients
                     (routine substitutes only non-zero values)
----------------------------------------------------------*/
void pdr_bc_neum_coeff(
	int Face,	/* in: face number */
	int Mat_num,	/* in: material number */
	double* Xcoor,	/* in: global coordinates of a point */
	double* Vec_norm,/* in: unit normal to the boundary */
	double  Time,	/* in: time instant */
	double* U_val, 	/* in: computed solution */
	double* U_x, 	/* in: gradient of computed solution */
	double* U_y,   	/* in: gradient of computed solution */
	double* U_z,   	/* in: gradient of computed solution */
	double* Gval	/* out: BC NEUMANN coefficients */	
	);

/**--------------------------------------------------------
pdr_bc_mixed_coeff - to return ROBIN (mixed) boundary coeficients
----------------------------------------------------------*/
void pdr_bc_mixed_coeff(
	int Face,	/* in: face number */
	int Mat_num,	/* in: material number */
	double* Xcoor,	/* in: global coordinates of a point */
	double* Vec_norm,/* in: unit normal to the boundary */
	double  Time,	/* in: time instant */
	double* U_val, 	/* in: computed solution */
	double* U_x, 	/* in: gradient of computed solution */
	double* U_y,   	/* in: gradient of computed solution */
	double* U_z,   	/* in: gradient of computed solution */
	double* Fval,	/* out: BC mixed coefficients */	
	double* Kr	/* out: proportionality coefficient */	
	);

#ifdef __cplusplus
}
#endif

#endif
