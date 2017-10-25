#ifndef LAS_PETSC_INTF_HPP
#define LAS_PETSC_INTF_HPP

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../../lsd_mkb/amg_ext/amg_ext.h"
#include "AMGSolverStructure.hpp"
#include "ResidualBasedErrorEvaluator.hpp"
#include "KnownSolutionErrorEvaluator.hpp"
extern "C"{
	#include "GhostBlockExchanger.h"
}


struct AMGSolverData{
	Mat A;
	Vec b;
	AMGSolverStructure* solver;
	int* posglobs;
	struct AMGAlgorithmData amgInitData;
};

extern struct AMGAlgorithmData amg_algorithm_data;
extern struct AMGSolverData amg_solver_data;

extern "C" {
void amg_project_solution_to_level(int level_id);

void amg_project_solution_from_level(int level_id);
}

extern "C" int lar_petsc_allocate_SM_and_LV( // returns: matrix index in itv_matrices array
  int Max_SM_size, /* maximal size of element stiffness matrix */
  int Nrblocks,    /* in: number of DOF blocks */
  int Nrdof_glob,  /* in: total number of DOFs */
  int* Nrdofbl,	   /* in: list of numbers of dofs in a block */
  int* Posglob,	   /* in: list of global numbers of first dof */
  int* Nroffbl,	   /* in: list of numbers of off diagonal blocks */
  // if Nroffbl[iblock]<=0 - iblock is a ghost block with no rows associated
  int** L_offbl	   /* in: list of lists of off diagonal blocks */
	);

extern "C" int lar_petsc_initialize_SM_and_LV(int Matrix_id, int Comp_type);

/*------------------------------------------------------------
  lar_petsc_assemble_SM_and_LV - to assemble entries to the global stiffness matrix
                           and the global load vector using the provided local
                           stiffness matrix and load vector
------------------------------------------------------------*/
extern "C" int lar_petsc_assemble_SM_and_LV(
                         /* returns: >=0 - success code, <0 - error code */
  int Matrix_id,   /* in: matrix ID */
  int Comp_type,         /* in: indicator for the scope of computations: */
                         /*   ITC_SOLVE - solve the system */
                         /*   ITC_RESOLVE - resolve for the new rhs vector */
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

extern "C" int lar_petsc_fill_preconditioner(int Matrix_id);

extern "C" int lar_petsc_free_SM_and_LV(int Matrix_id);

extern "C" void lar_petsc_perform_BJ_or_GS_iterations(
  int Matrix_id,   /* in: matrix ID */
  int Use_rhs,	/* in: 0 - no rhs, 1 - with rhs */
  int Ini_zero,	/* in: flag for zero initial guess */
  int Nr_prec,  /* in: number of preconditioner iterations */
  int Ndof,	/* in: number of unknowns (components of v*) */
  double* V,	/* in,out: vector of unknowns updated */
		/* during the loop over subdomains */
  double* B	/* in:  the rhs vector, if NULL take rhs */
		/*      from block data structure */
	);
#endif
