#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

#include "lsh_ns_supg_ext_intf.h"

#ifdef __cplusplus
//extern "C" {
#endif

/**-----------------------------------------------------------
   lsr_ns_supg_ext_init - to create a new solver instance, read its control
   parameters and initialize its data structure
   ------------------------------------------------------------*/
int lsr_ns_supg_ext_init( /* returns: >0 - solver ID, <0 - error code */
	int Solver_id,   /* in: solver ID (used to identify the subproblem) */
	int Parallel,      /* parameter specifying sequential (LSC_SEQUENTIAL) */
	/* or parallel (LSC_PARALLEL) execution */
	int* Max_num_levels_p,  /* in: number of levels for multigrid: */
	/*     1 - enforce single level solver */
	/*     >1 - enforce the number of levels */
	/* out: actual number of levels !!! */
	char* Filename,  /* in: name of the file with control parameters */
	int Max_iter, /* maximal number of iterations, -1 for values from Filename */
	int Error_type, /* type of error norm (stopping criterion), -1 for Filename*/
	double Error_tolerance, /* value for stopping criterion, -1.0 for Filename */
	int Monitoring_level /* Level of output, -1 for Filename */
	){
		return lsr_ns_supg_ext_petsc_init(Solver_id, Parallel, Max_num_levels_p, Filename, Max_iter, Error_type, Error_tolerance, Monitoring_level);
	}

/**--------------------------------------------------------
  lsr_ns_supg_ext_create_matrix - to allocate space for a global system matrix
  ---------------------------------------------------------*/
int lsr_ns_supg_ext_create_matrix(
	/* returns: >=0 - success code, <0 - error code */
	int Solver_id,   /* in: solver ID (used to identify the subproblem) */
	int Level_id,    /* in: level ID */
	int Nrblocks,    /* in: number of DOF blocks */
	int Nrdof_glob,  /* in: total number of DOFs */
	int Max_sm_size, /* in: maximal size of the stiffness matrix */
	int* Nrdofbl,	   /* in: list of numbers of dofs in a block */
	int* Posglob,	   /* in: list of global numbers of first dof */
	int* Nroffbl,	   /* in: list of numbers of off diagonal blocks */
	int** L_offbl	   /* in: list of lists of off diagonal blocks */
	)
{
	return lsr_ns_supg_ext_petsc_create_matrix(Solver_id, Level_id, Nrblocks, Nrdof_glob, Max_sm_size, Nrdofbl, Posglob, Nroffbl, L_offbl);
}


/**--------------------------------------------------------
lsr_ns_supg_ext_create_precon - to create preconditioner blocks corresponding
to small subdomains of neighboring elements
---------------------------------------------------------*/
int lsr_ns_supg_ext_create_precon( /* returns:   >0 number of diagonal blocks */
	/*	       <=0 - error */
	int Solver_id,         /* in: solver ID (used to identify the subproblem) */
	int Level_id           /* in: level ID */
	){

	return 0;
}

/**--------------------------------------------------------
  lsr_ns_supg_ext_clear_matrix - to initialize block structure of system matrix
  ---------------------------------------------------------*/
int lsr_ns_supg_ext_clear_matrix(
	/* returns: >=0 - success code, <0 - error code */
	int Solver_id,   /* in: solver ID (used to identify the subproblem) */
	int Level_id,    /* in: level ID */
	int Comp_type    /* in: indicator for the scope of computations: */
	/*   NS_SUPG_EXT_SOLVE - solve the system */
	/*   NS_SUPG_EXT_RESOLVE - resolve for the new rhs vector */
	){
	return lsr_ns_supg_ext_petsc_clear_matrix(Solver_id, Level_id, Comp_type);
}

/**-----------------------------------------------------------
  lsr_ns_supg_ext_assemble_local_sm - to assemble entries to the global stiffness matrix
  and the global load vector using the provided local
  stiffness matrix and load vector
	  ------------------------------------------------------------*/
int lsr_ns_supg_ext_assemble_local_sm(
	/* returns: >=0 - success code, <0 - error code */
	int Solver_id,         /* in: solver ID (used to identify the subproblem) */
	int Level_id,          /* in: level ID */
	int Comp_type,         /* in: indicator for the scope of computations: */
	/*   NS_SUPG_EXT_SOLVE - solve the system */
	/*   NS_SUPG_EXT_RESOLVE - resolve for the new rhs vector */
	int Nr_dof_bl,         /* in: number of global dof blocks */
	/*     associated with the local stiffness matrix */
	int* L_bl_id,          /* in: list of dof blocks' IDs */
	int* L_bl_nrdof,       /* in: list of blocks' numbers of dof */
	double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
	double* Rhs_vect,      /* in: rhs vector */
	char* Rewr_dofs         /* in: flag to rewrite or sum up entries */
	/*   'T' - true, rewrite entries when assembling */
	/*   'F' - false, sum up entries when assembling */
	){
		return lsr_ns_supg_ext_petsc_assemble_local_sm(Solver_id, Level_id, Comp_type, Nr_dof_bl, L_bl_id, L_bl_nrdof, Stiff_mat, Rhs_vect, Rewr_dofs);
}




int lsr_ns_supg_ext_show_matrix(
	/* returns: >=0 - success code, <0 - error code */
	int Solver_id,         /* in: solver ID (used to identify the subproblem) */
	int Level_id          /* in: level ID */

	){
	return 0;
}


/**--------------------------------------------------------
  lsr_ns_supg_ext_fill_precon - to prepare preconditioner by factorizing the stiffness matrix,
  either only diagonal blocks or block ILU(0)
  ---------------------------------------------------------*/
int lsr_ns_supg_ext_fill_precon(
	/* returns: >=0 - success code, <0 - error code */
	int Solver_id,         /* in: solver ID (used to identify the subproblem) */
	int Level_id           /* in: level ID */
	){

	return lsr_ns_supg_ext_petsc_fill_precon(Solver_id, Level_id);
}

/**--------------------------------------------------------
  lsr_ns_supg_ext_solve - to solve a system of equations, given previously constructed
  system matrix, preconditioner
  ---------------------------------------------------------*/
int lsr_ns_supg_ext_solve( /* returns: convergence indicator: */
	/* 1 - convergence */
	/* 0 - noconvergence */
	/* <0 - error code */
	int Solver_id,  /* in: solver ID */
	int Ndof, 	/* in: 	the number of degrees of freedom */
	int Ini_zero,   /* in:  indicator whether initial guess is zero (0/1) */
	double* X, 	/* in: 	the initial guess */
	/* out:	the iterated solution */
	double* B,	/* in:  the rhs vector, if NULL take rhs */
	/*      from block data structure */
	int* Nr_iter, 	/* in:	the maximum iterations to be performed */
	/* out:	actual number of iterations performed */
	double* Toler, 	/* in:	tolerance level for chosen measure */
	/* out:	the final value of convergence measure */
	int Monitor,	/* in:	flag to determine monitoring level */
	/*	0 - silent run, 1 - warning messages */
	/*	2 - 1+restart data, 3 - 2+iteration data */
	double* Conv_rate /* out: convergence rate */
	){

	return lsr_ns_supg_ext_petsc_solve(Solver_id, Ndof, Ini_zero, X, B, Nr_iter, Toler, Monitor, Conv_rate);
}

/**--------------------------------------------------------
  lsr_ns_supg_ext_free_matrix - to free space for solver data structure (system
  matrix, preconditioner and possibly other)
  ---------------------------------------------------------*/
int lsr_ns_supg_ext_free_matrix(
	int Solver_id   /* in: solver ID (used to identify the subproblem) */
	){
	return lsr_ns_supg_ext_petsc_free_matrix(Solver_id);
}

/**--------------------------------------------------------
  lsr_ns_supg_ext_compreres - to compute the residual of the left
  preconditioned system of equations, v = M^-1 * ( b - Ax )
  ( used also to compute the product v = M^-1 * Ax)
  ---------------------------------------------------------*/
extern void lsr_ns_supg_ext_compreres(
   int Solver_id,      /* in: pointer to solver data structure to be passed
   to data structure dependent routines */
   int Subsystem_id,  /* in: subsystem data structure to be used */
   int Level_id,	/* in: index of the mesh (level) */
   int Control,	/* in: indicator whether to compute residual (1)
				  or matrix-vector product (0)  */
   int Ini_zero,	/* in: flag for zero input vector X */
   int Ndof, 	/* in: number of unknowns (components of X and V) */
   double* X, 	/* in: input vector */
   double* B,	/* in: the rhs vector, if NULL take rhs */
   /*     from block data structure */
   double* V 	/* out: output vector, V = M^-1 * ( B - A*X ) */
   ){

	lsr_ns_supg_ext_petsc_compreres(Solver_id, Subsystem_id, Level_id, Control, Ini_zero, Ndof, X, B, V);
}

double lsr_ns_supg_ext_compres(
	/* returns: residual norm */
	int Solver_id,   /* in: solver ID (used to identify the subproblem) */
	double* X, /* in: input vector */
	int Ndof /* in: input vector size*/
	)
{
	return lsr_ns_supg_ext_petsc_compres(Solver_id, X, Ndof);
}

void lsr_ns_supg_ext_compres_vector(
	int Solver_id,   /* in: solver ID (used to identify the subproblem) */
	int Ndof, 	/* in: number of unknowns (components of x) */
	double* X, 	/* in: input vector (may be NULL if Ini_zero==0) */
	double* B,	/* in:  the rhs vector, if NULL take rhs */
					/*      from block data structure ( B is not taken */
					/*      into account if Use_rhs!=1) */
	double* V 	/* out: v = b-Ax */
	)
{
	return lsr_ns_supg_ext_petsc_compres_vector(Solver_id, Ndof, X, B, V);
}

#ifdef __cplusplus
//}
#endif
