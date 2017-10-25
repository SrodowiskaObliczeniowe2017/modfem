#include "lah_amg_interface.h"

/*---------------------------------------------------------
  lar_allocate_SM_and_LV - to allocate space for stiffness matrix and load vector
---------------------------------------------------------*/
int lar_allocate_SM_and_LV_petsc( // returns: matrix index in itv_matrices array
  int Max_SM_size, /* maximal size of element stiffness matrix */
  int Nrblocks,    /* in: number of DOF blocks */
  int Nrdof_glob,  /* in: total number of DOFs */
  int* Nrdofbl,	   /* in: list of numbers of dofs in a block */
  int* Posglob,	   /* in: list of global numbers of first dof */
  int* Nroffbl,	   /* in: list of numbers of off diagonal blocks */
  // if Nroffbl[iblock]<=0 - iblock is a ghost block with no rows associated
  int** L_offbl
	)
{
	return lar_petsc_allocate_SM_and_LV(Max_SM_size, Nrblocks, Nrdof_glob, Nrdofbl, Posglob, Nroffbl,
		 	L_offbl);
}



int lar_initialize_SM_and_LV_petsc(int Matrix_id, int Comp_type)
{
      return lar_petsc_initialize_SM_and_LV(Matrix_id, Comp_type);
}

double lar_get_storage_petsc(int Matrix_id)
{
	return 0;
}

/*------------------------------------------------------------
  lar_assemble_SM_and_LV - to assemble entries to the global stiffness matrix
                           and the global load vector using the provided local
                           stiffness matrix and load vector
------------------------------------------------------------*/
int lar_assemble_SM_and_LV_petsc(
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
  )
{
	return lar_petsc_assemble_SM_and_LV(Matrix_id, Comp_type, Nr_dof_bl, L_bl_id, L_bl_nrdof, Stiff_mat, Rhs_vect, Rewr_dofs);
}

int lar_allocate_preconditioner_petsc(int Matrix_id)
{
	return 0;
}

int lar_fill_preconditioner_petsc(int Matrix_id)
{
	return lar_petsc_fill_preconditioner(Matrix_id);
}

int lar_free_preconditioner_petsc(int Matrix_id)
{
	return 0;
}

int lar_free_SM_and_LV_petsc(int Matrix_id)
{
	return lar_petsc_free_SM_and_LV(Matrix_id);
}

void lar_compute_residual_petsc(int Matrix_id, int Use_rhs, int Ini_zero, int Ndof, double* X, double* B, double* V)
{
  return;
}

void lar_perform_BJ_or_GS_iterations_petsc(
  int Matrix_id,   /* in: matrix ID */
  int Use_rhs,	/* in: 0 - no rhs, 1 - with rhs */
  int Ini_zero,	/* in: flag for zero initial guess */
  int Nr_prec,  /* in: number of preconditioner iterations */
  int Ndof,	/* in: number of unknowns (components of v*) */
  double* V,	/* in,out: vector of unknowns updated */
		/* during the loop over subdomains */
  double* B	/* in:  the rhs vector, if NULL take rhs */
		/*      from block data structure */
	)
{
	return lar_petsc_perform_BJ_or_GS_iterations(Matrix_id, Use_rhs, Ini_zero, Nr_prec, Ndof, V, B);
}


void lar_perform_rhsub_petsc(int Matrix_id, int Ndof,	double* V, double* B)
{
	return;
}


int lar_block_print_matrix_petsc(int Matrix_id){
	return 0;
}

/*
void print_crs_matrix(itt_matrices *matrix)
{
	//itt_matrices *matrix = &itv_matrices[matrix_id];
	//printf("TOTAL COUNT %d\n", count_entries(matrix));
	int block_index;
	for(block_index = 1; block_index <= matrix->Nrblocks; block_index++)
	{
		itt_blocks* current_block = matrix->Block[block_index];

		if(current_block->Lngb != NULL)
		{
			int block_row_number, block_column_number;
			block_row_number = current_block->Posg;
			int i,j;
			//diagonal block
			for(i = 0; i<current_block->Ndof; i++)
			{
				block_column_number = current_block->Posg;
				for(j = 0; j<current_block->Ndof; j++)
				{
					//PetscErrorCode ierr;
					//ierr = MatSetValues(A,1,&block_row_number,1,&block_column_number,
						//	&(current_block->Dia[i + current_block->Ndof*j]),INSERT_VALUES);CHKERRABORT(PETSC_COMM_WORLD, ierr);

					printf("row %d column %d, value %12.3le\n", block_row_number, block_column_number,
							current_block->Dia[i + current_block->Ndof*j]);
					block_column_number++;
				}
				block_row_number++;
			}
			//printf("\n");
			//aux block
			int current_aux_block_number;
			for(current_aux_block_number = 1; current_aux_block_number <= current_block->Lngb[0];
					current_aux_block_number++)
			{
				itt_blocks* current_aux_block = matrix->Block[current_block->Lngb[current_aux_block_number]];
				block_row_number = current_block->Posg;
				for(i = 0; i<current_block->Ndof; i++)
				{
					block_column_number = current_aux_block->Posg;
					for(j = 0; j<current_aux_block->Ndof; j++)
					{
						//PetscErrorCode ierr;
						//ierr = MatSetValues(A,1,&block_row_number,1,&block_column_number,
								//&(current_block->Aux[current_aux_block_number - 1][i + current_block->Ndof*j]),
								//INSERT_VALUES);CHKERRABORT(PETSC_COMM_WORLD, ierr);
						printf("row %d column %d, value %12.3le\n", block_row_number, block_column_number,
								current_block->Aux[current_aux_block_number - 1][i + current_block->Ndof*j]);
						block_column_number++;
					}
					block_row_number++;
				}
			}

			//rhs
			for(i = 0; i<current_block->Ndof; i++)
			{
				//PetscErrorCode ierr;
				int global_colum_number = current_block->Posg + i;
				//ierr = VecSetValues(b,1,&global_colum_number,&(current_block->Rhs[i]),INSERT_VALUES); CHKERRABORT(PETSC_COMM_WORLD, ierr);
				printf("row %d value %12.3le \n", i, current_block->Rhs[i]);
			}
		}
		printf("\n");
	}
	printf("KONIEC SUMOWANIA\n");
	//PetscErrorCode ierr;
	//ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	//ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	//ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	//ierr = VecView(b,PETSC_VIEWER_STDOUT_WORLD); CHKERRABORT(PETSC_COMM_WORLD, ierr);
}
*/
