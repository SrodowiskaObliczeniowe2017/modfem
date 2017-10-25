
#ifndef LAS_NS_SUPG_PETSC_HPP
#define LAS_NS_SUPG_PETSC_HPP

#include <map>
#include <petscksp.h>
#include <time.h>
#include <math.h>
#include "uth_intf.h"
#include "uth_system.h"

#include "AMGSolverStructure.hpp"
#include "ResidualBasedErrorEvaluator.hpp"
#include "KnownSolutionErrorEvaluator.hpp"
#include "las_petsc_intf.hpp"
extern "C" {
#include "../../lsd_mkb/lah_intf.h"
}
#include "ns_supg/MatrixUtil.hpp"
#include "ns_supg/ApproximateInverseOpt.hpp"
#include "ns_supg/VectorTransformationUtil.hpp"
#include "ns_supg/SchurComplement.hpp"


extern struct AMGSolverData amg_solver_data;

int lsr_ns_supg_ext_petsc_init( /* returns: >0 - solver ID, <0 - error code */
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
		);


int lsr_ns_supg_ext_petsc_create_matrix(
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
		);

int lsr_ns_supg_ext_petsc_assemble_local_sm(
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
		);

int lsr_ns_supg_ext_petsc_assemble_local_mm(
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
		);

int lsr_ns_supg_ext_petsc_fill_precon(
	/* returns: >=0 - success code, <0 - error code */
	int Solver_id,         /* in: solver ID (used to identify the subproblem) */
	int Level_id           /* in: level ID */
	);

int lsr_ns_supg_ext_petsc_solve( /* returns: convergence indicator: */
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
	);

void lsr_ns_supg_ext_petsc_compreres (
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
	);

int lsr_ns_supg_ext_petsc_free_matrix(int Solver_id);

int lsr_ns_supg_ext_petsc_clear_matrix(int Solver_id, int Level_id, int Comp_type);

double lsr_ns_supg_ext_petsc_compres(int Solver_id, double* X, int Ndof);

void lsr_ns_supg_ext_petsc_compres_vector(int Solver_id, int Ndof, double* X, double* B, double* V);

//for testing
void lsr_ns_supg_ext_petsc_get_system(Mat* Avv, Mat* Avp, Mat* Apv, Mat* App, Vec* v, Vec* p);

#endif
