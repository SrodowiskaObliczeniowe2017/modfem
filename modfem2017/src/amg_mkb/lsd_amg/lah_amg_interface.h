#ifndef LAH_AMG_INTERFACE_H
#define LAH_AMG_INTERFACE_H

#include "../../lsd_mkb/amg_ext/lah_petsc_interface.h"
#include "../lad_petsc/las_petsc_intf.hpp"


extern "C" int lar_allocate_SM_and_LV_petsc(
  int Max_SM_size,
  int Nrblocks,
  int Nrdof_glob,
  int* Nrdofbl,
  int* Posglob,
  int* Nroffbl,
  int** L_offbl);

extern "C" int lar_initialize_SM_and_LV_petsc(int Matrix_id, int Comp_type);

extern "C" double lar_get_storage_petsc(int Matrix_id);

extern "C" int lar_assemble_SM_and_LV_petsc(
  int Matrix_id,
  int Comp_type,
  int Nr_dof_bl,
  int* L_bl_id,
  int* L_bl_nrdof,
  double* Stiff_mat,
  double* Rhs_vect,
  char* Rewr_dofs);

extern "C" int lar_allocate_preconditioner_petsc(int Matrix_id);

extern "C" int lar_fill_preconditioner_petsc(int Matrix_id);

extern "C" int lar_free_preconditioner_petsc(int Matrix_id);

extern "C" int lar_free_SM_and_LV_petsc(int Matrix_id);

extern "C" void lar_compute_residual_petsc(int Matrix_id, int Use_rhs, int Ini_zero, int Ndof, double* X, double* B, double* V);

extern "C" void lar_perform_BJ_or_GS_iterations_petsc(
  int Matrix_id,
  int Use_rhs,
  int Ini_zero,
  int Nr_prec,
  int Ndof,
  double* V,
  double* B);

extern "C" void lar_perform_rhsub_petsc(int Matrix_id, int Ndof,	double* V, double* B);

extern "C" int lar_block_print_matrix_petsc(int Matrix_id);

#endif
