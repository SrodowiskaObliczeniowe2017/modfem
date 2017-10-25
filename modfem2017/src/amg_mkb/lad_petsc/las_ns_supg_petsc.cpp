#include "las_ns_supg_petsc.hpp"

double get_timeClock(){
	return time_clock();
}

struct AMGSUPGSolverData{
	Mat Avv;
	Mat Avp;
	Mat Apv;
//	Mat App;
	Mat schur_complement;
	Mat velocity_inverse;
	KSP preonly_solver;
	Vec bvp;
	Mat Vmm;
	Mat commutator;
//	Vec bpv;
//	AMGSolverStructure* solver = NULL;
	int* posglobs;
} amg_ns_supg_solver_data {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};

Vec solution_v;
Vec solution_p;
Vec result_v;
Vec result_p;
Vec zero_vec;
double total_time = 0;
double schur_time = 0;
double vv_time = 0;

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
	int Monitoring_level /* Level of output, -1 for Filename */)
	{
		return 0;
	}

void create_ns_supg_matrix(Mat* matrix, PetscInt rows_number, PetscInt columns_number)
{
	  PetscErrorCode ierr;
	  ierr = MatCreate(PETSC_COMM_WORLD, matrix);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = MatSetSizes(*matrix,PETSC_DECIDE,PETSC_DECIDE,rows_number,columns_number);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = MatSetFromOptions(*matrix);CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

void create_ns_supg_vector(Vec* vector, PetscInt rows_number)
{
	  PetscErrorCode ierr;
	  ierr = VecCreate(PETSC_COMM_WORLD,vector);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = VecSetSizes(*vector,PETSC_DECIDE,rows_number);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = VecSetFromOptions(*vector);CHKERRABORT(PETSC_COMM_WORLD, ierr);
}


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
	)
	{
	  amg_ns_supg_solver_data.posglobs = (int*)malloc(sizeof(int) * Nrblocks);
	  memcpy(amg_ns_supg_solver_data.posglobs, Posglob,sizeof(int)* Nrblocks);

	  PetscBool initialized;
	  PetscInitialized(&initialized);
	  if(!initialized)
		  PetscInitializeNoArguments();

	  // ERROR in 3.7.4 - need 2 arguments - KB
	  //PetscOptionsView(PETSC_VIEWER_STDOUT_WORLD);

	  create_ns_supg_matrix(&(amg_ns_supg_solver_data.Avv), Nrdof_glob*3/4, Nrdof_glob*3/4);
	  create_ns_supg_matrix(&(amg_ns_supg_solver_data.Avp), Nrdof_glob*3/4, Nrdof_glob*1/4);
	  create_ns_supg_matrix(&(amg_ns_supg_solver_data.Apv), Nrdof_glob*1/4, Nrdof_glob*3/4);
//	  create_ns_supg_matrix(&(amg_ns_supg_solver_data.App), Nrdof_glob*1/4, Nrdof_glob*1/4);

	  create_ns_supg_vector(&(amg_ns_supg_solver_data.bvp), Nrdof_glob * 3 / 4);
//	  create_ns_supg_vector(&(amg_ns_supg_solver_data.bpv), Nrdof_glob * 1 / 4);

	  int block_index;
	  int max_nonzeros_in_row = 0;
	  for(block_index=0;block_index<Nrblocks;block_index++)
	  {
	      int diagonal_block_size = Nrdofbl[block_index];
	      int off_diagonal_block_sizes = 0;
	      if(Nroffbl[block_index]>0)
	      {
	    	  int i;
	    	  for(i=0;i<Nroffbl[block_index];i++)
	    	  {
	    		  int offdiagonal_block_id = L_offbl[block_index][i];
	    		  off_diagonal_block_sizes += Nrdofbl[offdiagonal_block_id];
	    	  }
	      }
	      max_nonzeros_in_row = std::max(max_nonzeros_in_row, diagonal_block_size + off_diagonal_block_sizes);
	  }

	  PetscErrorCode ierr;
	  //TODO: set the second parameter
	  ierr = MatMPIAIJSetPreallocation(amg_ns_supg_solver_data.Avv,max_nonzeros_in_row*3/4,NULL,0,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = MatSeqAIJSetPreallocation(amg_ns_supg_solver_data.Avv,max_nonzeros_in_row*3/4,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = MatSetUp(amg_ns_supg_solver_data.Avv); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = MatMPIAIJSetPreallocation(amg_ns_supg_solver_data.Avp,max_nonzeros_in_row*1/4,NULL,0,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = MatSeqAIJSetPreallocation(amg_ns_supg_solver_data.Avp,max_nonzeros_in_row*1/4,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = MatSetUp(amg_ns_supg_solver_data.Avp); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = MatMPIAIJSetPreallocation(amg_ns_supg_solver_data.Apv,max_nonzeros_in_row*3/4,NULL,0,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = MatSeqAIJSetPreallocation(amg_ns_supg_solver_data.Apv,max_nonzeros_in_row*3/4,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = MatSetUp(amg_ns_supg_solver_data.Apv); CHKERRABORT(PETSC_COMM_WORLD, ierr);

	  create_ns_supg_matrix(&(amg_ns_supg_solver_data.Vmm), Nrdof_glob*3/4, Nrdof_glob*3/4);
	  ierr = MatMPIAIJSetPreallocation(amg_ns_supg_solver_data.Vmm,max_nonzeros_in_row*3/4,NULL,0,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = MatSeqAIJSetPreallocation(amg_ns_supg_solver_data.Vmm,max_nonzeros_in_row*3/4,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	  ierr = MatSetUp(amg_ns_supg_solver_data.Vmm); CHKERRABORT(PETSC_COMM_WORLD, ierr);

	  int* nrdofbl_app = (int*)malloc(sizeof(int)*Nrblocks);
	  for(int i = 0; i<Nrblocks; i++)
		  nrdofbl_app[i] = 1;
	  int* posglob_app = (int*)malloc(sizeof(int)*Nrblocks);
	  for(int i = 0; i<Nrblocks; i++)
	  {
		  posglob_app[i] = Posglob[i] / 4;
	  }
	  lar_allocate_SM_and_LV(LAC_STORAGE_PETSC,99,Nrblocks,Max_sm_size,Nrblocks,0,nrdofbl_app,posglob_app,Nroffbl,L_offbl);

//	  lar_petsc_allocate_SM_and_LV(Max_sm_size,Nrblocks,Nrblocks,nrdofbl_app,posglob_app,Nroffbl,L_offbl);
	  free(nrdofbl_app);
	  free(posglob_app);


//	  ierr = MatMPIAIJSetPreallocation(amg_ns_supg_solver_data.App,max_nonzeros_in_row*1/4,NULL,0,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSeqAIJSetPreallocation(amg_ns_supg_solver_data.App,max_nonzeros_in_row*1/4,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSetUp(amg_ns_supg_solver_data.App); CHKERRABORT(PETSC_COMM_WORLD, ierr);

	  return 0;
	}

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
	)
	{

	PetscErrorCode ierr;
	int i;
	int j;
	int all_blocks_height = 0;
	for(i = 0; i<Nr_dof_bl; i++)
		all_blocks_height += L_bl_nrdof[i];

	int row_offset = 0;
	for(i = 0; i<Nr_dof_bl; i++)
	{
		int column_offset = 0;
		for(j = 0; j<Nr_dof_bl; j++)
		{
			int k;
			int l;

			//Avv
			int j_index = amg_ns_supg_solver_data.posglobs[L_bl_id[j]] * 3 / 4;
			for(k = 0; k<3; k++)
			{
				int i_index = amg_ns_supg_solver_data.posglobs[L_bl_id[i]] * 3 / 4;
				for(l = 0; l<3; l++)
				{
					double value = Stiff_mat[(column_offset+k)*all_blocks_height + row_offset +l];
					//printf("(%d %d) v %lf ",i_index,j_index,value);
					ierr = MatSetValues(amg_ns_supg_solver_data.Vmm,1,&i_index,1,&j_index,&value,ADD_VALUES); CHKERRABORT(PETSC_COMM_WORLD, ierr);
					i_index++;
				}
				j_index++;
				//printf("\n");
			}
			column_offset += L_bl_nrdof[j];
		}
		row_offset += L_bl_nrdof[i];
	}

  return(1);
}

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
	)
	{
	if(Comp_type != 0 && Comp_type != 1){
		lsr_ns_supg_ext_petsc_assemble_local_mm(Solver_id, Level_id, Comp_type, Nr_dof_bl,
				L_bl_id, L_bl_nrdof, Stiff_mat, Rhs_vect, Rewr_dofs);
		return 1;
	}

	PetscErrorCode ierr;
	int i;
	int j;
	int all_blocks_height = 0;
	for(i = 0; i<Nr_dof_bl; i++)
		all_blocks_height += L_bl_nrdof[i];

	double* app_local_mat = (double*)malloc(sizeof(double)*Nr_dof_bl*Nr_dof_bl);
	double* app_local_rhs = (double*)malloc(sizeof(double)*Nr_dof_bl);
	int* app_local_L_bl_nrdof = (int*)malloc(sizeof(int)*Nr_dof_bl);
	for(int i = 0; i<Nr_dof_bl; i++)
		app_local_L_bl_nrdof[i] = 1;


//
	int posglob_short_column_buffer[20];
	int posglob_column_buffer[20*3];


	if(all_blocks_height > 80 || Nr_dof_bl > 20){
		printf("Nr_dof_bl max is 20 not %d\n", Nr_dof_bl);
		CHKERRABORT(PETSC_COMM_WORLD, 17);
	}

	int column_index = 0;
	for(i = 0; i<Nr_dof_bl; i++){
		int posglob = amg_ns_supg_solver_data.posglobs[L_bl_id[i]]*3/4;
		posglob_short_column_buffer[column_index/3] = posglob / 3;
		for(j = 0; j<3; j++){
			posglob_column_buffer[column_index++] = posglob + j;
		}

	}

	double local_avv[9*20*20];
	double local_avp[3*20*20];
	double local_apv[3*20*20];
	double local_app[1*20*20];

	for(int i = 0; i<Nr_dof_bl; i++){
		for(int j = 0; j<Nr_dof_bl; j++){
			for(int k = 0; k<3; k++){
				for(int l = 0; l<3; l++){
					local_avv[(i*3+k)*all_blocks_height*3/4 + j*3 + l] = Stiff_mat[(i*4+k) + (j*4 + l)*all_blocks_height];
				}
			}
			for(int k = 0; k<3; k++){
				for(int l = 0; l<1; l++){
					local_avp[(i*3+k)*all_blocks_height*1/4 + j + l] = Stiff_mat[(i*4+k) + (j*4 + l + 3)*all_blocks_height];
				}
			}
			for(int k = 0; k<1; k++){
				for(int l = 0; l<3; l++){
					local_apv[(i+k)*all_blocks_height*3/4 + j*3 + l] = Stiff_mat[(i*4+k+3) + (j*4 + l)*all_blocks_height];
				}
			}
			for(int k = 0; k<1; k++){
				for(int l = 0; l<1; l++){
					local_app[(i+k)*all_blocks_height*1/4 + j + l] = Stiff_mat[(i*4+k+3) + (j*4 + l + 3)*all_blocks_height];
				}
			}
		}
	}

	int row_index = 0;
//	printf("%d ", posglob_column_buffer[row_index]);
	for(i = 0; i<Nr_dof_bl; i++){
//		int posglob = amg_solver_data.posglobs[L_bl_id[i]];
//		if(first_ghost_block == -1 || (amg_solver_data.posglobs[L_bl_id[i]]< first_ghost_block + posglob_offset
//				&& amg_solver_data.posglobs[L_bl_id[i]] >= posglob_offset)){

			ierr = MatSetValues(amg_ns_supg_solver_data.Avv,3,posglob_column_buffer + row_index,
					column_index,posglob_column_buffer,local_avv + (all_blocks_height*3/4)*row_index,ADD_VALUES);
			CHKERRABORT(PETSC_COMM_WORLD, ierr);

			ierr = MatSetValues(amg_ns_supg_solver_data.Avp,3,posglob_column_buffer + row_index,
					column_index/3,posglob_short_column_buffer,local_avp + (all_blocks_height*1/4)*row_index,ADD_VALUES);
			CHKERRABORT(PETSC_COMM_WORLD, ierr);

			ierr = MatSetValues(amg_ns_supg_solver_data.Apv,1,posglob_short_column_buffer + row_index/3,
					column_index,posglob_column_buffer,local_apv + (all_blocks_height*3/4)*row_index/3,ADD_VALUES);
			CHKERRABORT(PETSC_COMM_WORLD, ierr);

			ierr = VecSetValues(amg_ns_supg_solver_data.bvp,3,posglob_column_buffer + row_index,
					&(Rhs_vect[row_index + row_index/3]),ADD_VALUES);
			CHKERRABORT(PETSC_COMM_WORLD, ierr);
//		}

		row_index += 3;
	}

	int row_offset = 0;
	for(i = 0; i<Nr_dof_bl; i++)
	{
		int column_offset = 0;
		for(j = 0; j<Nr_dof_bl; j++)
		{
			//App
			int j_index = amg_ns_supg_solver_data.posglobs[L_bl_id[j]] * 1 / 4;
			int i_index = amg_ns_supg_solver_data.posglobs[L_bl_id[i]] * 1 / 4;

			double value = Stiff_mat[(column_offset+3)*all_blocks_height + row_offset + 3];
			//printf("(%d %d) v %lf ",i_index,j_index,value);
			//ierr = MatSetValues(amg_ns_supg_solver_data.App,1,&i_index,1,&j_index,&value,ADD_VALUES); CHKERRABORT(PETSC_COMM_WORLD, ierr);
			app_local_mat[column_offset/4*all_blocks_height/4 + row_offset/4] = value;
			column_offset += L_bl_nrdof[j];
		}
		int r;
		for(r = 0; r<1; r++)
		{
			//printf("rhs %d wartosc %lf\n",r + posglobs[L_bl_id[i]] * 1 / 4,Rhs_vect[r+row_offset + 3]);
			int row = r + amg_ns_supg_solver_data.posglobs[L_bl_id[i]] * 1 / 4;
//			ierr = VecSetValues(amg_ns_supg_solver_data.bpv,1,&row,&(Rhs_vect[r+row_offset+3]),ADD_VALUES); CHKERRABORT(PETSC_COMM_WORLD, ierr);
			app_local_rhs[row_offset/4] = Rhs_vect[r+row_offset+3];
		}
		row_offset += L_bl_nrdof[i];
	}

  lar_petsc_assemble_SM_and_LV(0,0,Nr_dof_bl,L_bl_id,app_local_L_bl_nrdof,app_local_mat,app_local_rhs,NULL);
  free(app_local_mat);
  free(app_local_rhs);
  free(app_local_L_bl_nrdof);
  return(1);
}

int lsr_ns_supg_ext_petsc_solve(
	int Solver_id,
	int Ndof,
	int Ini_zero,
	double* X,
	double* B,
	int* Nr_iter,
	double* Toler,
	int Monitor,
	double* Conv_rate
	){
	return 0;
}

int lsr_ns_supg_ext_petsc_fill_precon(
	/* returns: >=0 - success code, <0 - error code */
	int Solver_id,         /* in: solver ID (used to identify the subproblem) */
	int Level_id           /* in: level ID */
	){

//	  lar_petsc_fill_preconditioner(0);

	PetscErrorCode ierr;

	ierr = MatAssemblyBegin(amg_ns_supg_solver_data.Avv,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatAssemblyEnd(amg_ns_supg_solver_data.Avv,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatAssemblyBegin(amg_ns_supg_solver_data.Avp,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatAssemblyEnd(amg_ns_supg_solver_data.Avp,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatAssemblyBegin(amg_ns_supg_solver_data.Apv,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatAssemblyEnd(amg_ns_supg_solver_data.Apv,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatAssemblyBegin(amg_solver_data.A,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatAssemblyEnd(amg_solver_data.A,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecAssemblyBegin(amg_solver_data.b); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecAssemblyEnd(amg_solver_data.b); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecAssemblyBegin(amg_ns_supg_solver_data.bvp); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecAssemblyEnd(amg_ns_supg_solver_data.bvp); CHKERRABORT(PETSC_COMM_WORLD, ierr);

//	PetscViewer viewer;
//	PetscViewerASCIIOpen(PETSC_COMM_WORLD, "OptAvv.data", &viewer);
//	MatView(amg_ns_supg_solver_data.Avv, viewer);
//	PetscViewerASCIIOpen(PETSC_COMM_WORLD, "OptAvp.data", &viewer);
//	MatView(amg_ns_supg_solver_data.Avp, viewer);
//	PetscViewerASCIIOpen(PETSC_COMM_WORLD, "OptApv.data", &viewer);
//	MatView(amg_ns_supg_solver_data.Apv, viewer);
//	PetscViewerASCIIOpen(PETSC_COMM_WORLD, "OptApp.data", &viewer);
//	MatView(amg_solver_data.A, viewer);
//
//	PetscViewerASCIIOpen(PETSC_COMM_WORLD, "OptV.data", &viewer);
//	VecView(amg_ns_supg_solver_data.bvp, viewer);
//	PetscViewerASCIIOpen(PETSC_COMM_WORLD, "OptP.data", &viewer);
//	VecView(amg_solver_data.b, viewer);
//	PetscViewerDestroy(&viewer);
//
//
//	CHKERRABORT(PETSC_COMM_WORLD, 17);

//	ierr = MatAssemblyBegin(amg_ns_supg_solver_data.Vmm,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	ierr = MatAssemblyEnd(amg_ns_supg_solver_data.Vmm,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);


//	create_exact_schur_complement(&(amg_ns_supg_solver_data.schur_complement), amg_ns_supg_solver_data.Avp, amg_ns_supg_solver_data.Avv,
//			amg_ns_supg_solver_data.Apv,amg_solver_data.A);
//
//	amg_ns_supg_solver_data.velocity_inverse =
//			create_velocity_exact_inverse(&amg_ns_supg_solver_data.velocity_inverse, amg_ns_supg_solver_data.Avv, amg_ns_supg_solver_data.bvp);
//	draw_matrix(amg_ns_supg_solver_data.velocity_inverse, 400);
//	MatView(amg_ns_supg_solver_data.Avv, PETSC_VIEWER_STDOUT_WORLD);
//	MatView(amg_ns_supg_solver_data.Avv, PETSC_VIEWER_DRAW_WORLD);

	clock_t begin = clock();
//	amg_ns_supg_solver_data.velocity_inverse = create_velocity_ilu_inverse(&(amg_ns_supg_solver_data.preonly_solver), amg_ns_supg_solver_data.Avv);
//	amg_ns_supg_solver_data.velocity_inverse = create_velocity_mass_inverse(amg_ns_supg_solver_data.Vmm, amg_ns_supg_solver_data.bvp);
	//155 140
//	ierr = MatCreateSchurComplementPmat(
//			amg_ns_supg_solver_data.Avv,
//			amg_ns_supg_solver_data.Avp,
//			amg_ns_supg_solver_data.Apv,
//			amg_solver_data.A,
//			MAT_SCHUR_COMPLEMENT_AINV_DIAG,MAT_INITIAL_MATRIX,&(amg_ns_supg_solver_data.schur_complement));

//	create_s_complement_ilu(amg_ns_supg_solver_data.Avv, amg_ns_supg_solver_data.Avp, amg_ns_supg_solver_data.Apv, amg_solver_data.A,&(amg_ns_supg_solver_data.schur_complement));
//	create_s_complement_Vmm(amg_ns_supg_solver_data.Avv, amg_ns_supg_solver_data.Avp, amg_ns_supg_solver_data.Apv, amg_solver_data.A,&(amg_ns_supg_solver_data.schur_complement));
	//Mat restricted_matrix = restrict_matrix(amg_ns_supg_solver_data.velocity_inverse, 1000);
//	amg_ns_supg_solver_data.commutator = getApproximateInverseCommutator(amg_ns_supg_solver_data.Avv, amg_ns_supg_solver_data.Avp, 20);
//	print_sym_difference(amg_ns_supg_solver_data.Avv);
//	print_sym_difference(amg_solver_data.A);
	amg_ns_supg_solver_data.velocity_inverse = create_s_complement_spai(amg_ns_supg_solver_data.Avv, amg_ns_supg_solver_data.Avp, amg_ns_supg_solver_data.Apv, amg_solver_data.A,
			NULL, &(amg_ns_supg_solver_data.schur_complement));
//	MatDestroy(&restricted_matrix);
//	VecView(amg_ns_supg_solver_data.bvp, PETSC_VIEWER_STDOUT_WORLD);
	//change_pressure_for_boundary_matrix();
	CHKERRABORT(PETSC_COMM_WORLD, ierr);
	int size;
	MatGetSize(amg_ns_supg_solver_data.schur_complement,&size,NULL);
	MatScale(amg_ns_supg_solver_data.schur_complement, -1.0);
//
//	Vec scaling;
//	VecDuplicate(amg_solver_data.b, &scaling);
//	MatGetDiagonal(amg_ns_supg_solver_data.schur_complement, scaling);
//	int pin_index;
//	double xxxx;
//	for(pin_index = 0; pin_index<size; pin_index++){
//		VecGetValues(scaling,1,&pin_index,&xxxx);
//		if(xxxx > 1.0e+7 || xxxx < -1.0e+7)
//			break;
//	}
//
//	VecSet(scaling, 1.0);
//
//	double scaling_value = 1.0e-7;
//	VecSetValues(scaling,1,&pin_index,&scaling_value,INSERT_VALUES);
//	MatDiagonalScale(amg_ns_supg_solver_data.schur_complement, scaling, NULL);

	amg_solver_data.solver = new AMGSolverStructure(amg_ns_supg_solver_data.schur_complement, NULL, size, NULL, false);
//	MatView(amg_ns_supg_solver_data.schur_complement, PETSC_VIEWER_STDOUT_WORLD);

	amg_solver_data.solver->CreateAMGSolverLevels(amg_solver_data.amgInitData.coarsening_algorithm,
			amg_solver_data.amgInitData.interpolation_algorithm,
			amg_solver_data.amgInitData.strength_threshold,
			amg_solver_data.amgInitData.levels_number);

//	scaling_value = 1.0e+7;
//	VecSetValues(scaling,1,&pin_index,&scaling_value,INSERT_VALUES);
//	MatDiagonalScale(amg_ns_supg_solver_data.schur_complement, scaling, NULL);
//	VecDestroy(&scaling);


	ierr = VecCreate(PETSC_COMM_WORLD,&solution_v);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecSetSizes(solution_v,PETSC_DECIDE,size*3);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecSetFromOptions(solution_v);CHKERRABORT(PETSC_COMM_WORLD, ierr);

	ierr = VecCreate(PETSC_COMM_WORLD,&solution_p);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecSetSizes(solution_p,PETSC_DECIDE,size);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecSetFromOptions(solution_p);CHKERRABORT(PETSC_COMM_WORLD, ierr);

	VecDuplicate(solution_v, &result_v);
	VecDuplicate(solution_p, &result_p);

	ierr = VecCreate(PETSC_COMM_WORLD,&zero_vec);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecSetSizes(zero_vec,PETSC_DECIDE,size);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecSetFromOptions(zero_vec);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	VecSet(zero_vec, 0.0);

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	printf("Subsystems creation time %lf \n",elapsed_secs);

//		MatScale(Avv,-1.0);
//		MatScale(Apv,-1.0);
//		MatScale(Avp,-1.0);
//		MatScale(App,-1.0);
//		VecScale(bvp,-1.0);
//		VecScale(bpv,-1.0);
	//solve(Ndof, V);
	return 0;
}

int lsr_ns_supg_ext_petsc_free_matrix(int Solver_id)
{
	MatDestroy(&(amg_ns_supg_solver_data.Avv));
	MatDestroy(&(amg_ns_supg_solver_data.Avp));
	MatDestroy(&(amg_ns_supg_solver_data.Apv));
	MatDestroy(&(amg_ns_supg_solver_data.Vmm));
	MatDestroy(&(amg_ns_supg_solver_data.commutator));
//	MatDestroy(&(amg_ns_supg_solver_data.App));
	VecDestroy(&(amg_ns_supg_solver_data.bvp));
//	VecDestroy(&(amg_ns_supg_solver_data.bpv));
	lar_petsc_free_SM_and_LV(0);
	VecDestroy(&zero_vec);
	VecDestroy(&solution_v);
	VecDestroy(&solution_p);
	VecDestroy(&result_v);
	VecDestroy(&result_p);
	return 0;
}

int lsr_ns_supg_ext_petsc_clear_matrix(int Solver_id, int Level_id, int Comp_type)
{
	lar_petsc_initialize_SM_and_LV(0,0);
	MatZeroEntries(amg_ns_supg_solver_data.Avv);
	MatZeroEntries(amg_ns_supg_solver_data.Avp);
	MatZeroEntries(amg_ns_supg_solver_data.Apv);
//	MatZeroEntries(amg_ns_supg_solver_data.App);
	MatDestroy(&(amg_ns_supg_solver_data.schur_complement));
	//velocity inveres
	if(amg_ns_supg_solver_data.preonly_solver != NULL)
	{
		KSPDestroy(&(amg_ns_supg_solver_data.preonly_solver));
	}
	if(amg_ns_supg_solver_data.velocity_inverse != NULL)
	{
		MatDestroy(&(amg_ns_supg_solver_data.velocity_inverse));
	}
//	else
//	{
//		MatZeroEntries(velocity_inverse);
//	}

	VecSet(amg_ns_supg_solver_data.bvp, 0.0);
//	VecSet(amg_ns_supg_solver_data.bpv, 0.0);
//	delete amg_ns_supg_solver_data.solver;
//	amg_ns_supg_solver_data.solver = NULL;
	if(total_time != 0){
		printf("TOTAL TIME: %lf\n", total_time);
		printf("SCHUR TIME: %lf\n", schur_time);
		printf("VV TIME: %lf\n", vv_time);
	}
	total_time = 0;
	schur_time = 0;
	vv_time = 0;
	return 0;
}

void simplest_preconditioner_execution(int Ndof, double* V, Vec u_v, Vec v_residual, Vec bvp_residual, Vec bpv_residual)
{
	int size;
	MatGetSize(amg_ns_supg_solver_data.Avv,&size,NULL);
	GaussSidelMethod* gaussSidelMethod = new GaussSidelMethod();
	gaussSidelMethod->PreSmoothing(amg_ns_supg_solver_data.Avv, bvp_residual);
	ErrorEvaluator* velocity_error_evaluator =
			new ResidualBasedErrorEvaluator(amg_ns_supg_solver_data.Avv, bvp_residual, 1.0e-8, 1.0e-8, size, (int)sqrt(size)/40);
	//new ResidualBasedErrorEvaluator(Avv, bvp_residual, 1.0e-8, 1.0e-7, size);
	clock_t begin = clock();
	int iteration_nr = 0;
	PetscReal norm;
	do
	{
		gaussSidelMethod->Smooth(u_v);
	}
	while(!velocity_error_evaluator->Stop(u_v,iteration_nr++));
	delete velocity_error_evaluator;
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//printf("Total time of Dvv solution: %lf\n",elapsed_secs);

	amg_solver_data.solver->ChangeSystem(bpv_residual,(int)sqrt(size)/80);
//	MatGetSize(schur_complement,&size,NULL);
//	ErrorEvaluator* amg_error_evaluator = new ResidualBasedErrorEvaluator(schur_complement, bpv_residual, 1.0e-7, 1.0e-8, size);
//	AMGSolverStructure* solver = new AMGSolverStructure(schur_complement, bpv_residual, size, amg_error_evaluator, false);
//	clock_t begin = clock();
//	solver->CreateAMGSolverLevels();
//	clock_t end = clock();
//	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//time_t start = time(NULL);
	begin = clock();
	Vec u_p = amg_solver_data.solver->RunVCycle();
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//	printf("Schur: %lf \n",elapsed_secs);

	convert_vec_to_x(u_v, u_p, V, Ndof);
	//VecView(u_v, PETSC_VIEWER_STDOUT_WORLD);
	//VecView(u_p, PETSC_VIEWER_STDOUT_WORLD);
	VecDestroy(&u_p);
	delete gaussSidelMethod;
}

void simplest_preconditioner(int Ndof, double* V, Vec bvp_residual, Vec bpv_residual)
{
	Vec u_v;
	Vec v_residual;
	VecDuplicate(amg_ns_supg_solver_data.bvp,&u_v);
	VecDuplicate(amg_ns_supg_solver_data.bvp,&v_residual);
	VecSet(u_v, 0.0);

	simplest_preconditioner_execution(Ndof, V, u_v, v_residual, bvp_residual, bpv_residual);

	VecDestroy(&u_v);
	VecDestroy(&v_residual);
}

void simpler_preconditioner_iteration(double gaussIterationsError, Vec* bvp_residual, Vec u_v, Vec v_residual, Vec intermediate_v,
		Vec v_tmp, Vec* bpv_residual, Vec u_p, Vec intermediate_p, Vec p_tmp, bool simpler)
{
	Mat X = amg_ns_supg_solver_data.velocity_inverse;
////	MatView(amg_ns_supg_solver_data.velocity_inverse, PETSC_VIEWER_STDOUT_WORLD);
//	VecView(*bvp_residual,PETSC_VIEWER_STDOUT_WORLD);
//	VecView(*bpv_residual,PETSC_VIEWER_STDOUT_WORLD);
	if(simpler)
	{
//		MatSolve(X, *bvp_residual, intermediate_v);
		MatMult(X, *bvp_residual, intermediate_v);

		//VecView(intermediate_v,PETSC_VIEWER_STDOUT_WORLD);
		MatMult(amg_ns_supg_solver_data.Apv, intermediate_v, intermediate_p);
		VecAYPX(intermediate_p,-1.0,*bpv_residual);
	}

	int size;
	MatGetSize(amg_solver_data.A,&size,NULL);
//
//	ErrorEvaluator* error_evaluator = new ResidualBasedErrorEvaluator(schur_complement, intermediate_p, 1.0e-8, 1.0e-8, size);
//	AMGSolverStructure* solver = new AMGSolverStructure(schur_complement, intermediate_p, size, error_evaluator, false);
//	//clock_t begin = clock();
//	solver->CreateAMGSolverLevels();
	//clock_t end = clock();
	//double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//printf("Total time of level creation: %lf\n",elapsed_secs);
	//time_t start = time(NULL);
	//begin = clock();



	clock_t begin = clock();
//	for(int i = 0; i<(int)sqrt(size)/80; i++)
//	{
//		Vec u_p_hat = amg_solver_data.solver->RunVCycle();
//		VecDestroy(&u_p_hat);
//	}
//	Vec u_p_hat = amg_solver_data.solver->RunVCycle();
//	VecSet(u_p_hat, 0.0);
	Vec u_p_hat;
	if(simpler)
	{
		amg_solver_data.solver->ChangeSystem(intermediate_p, 1);
		u_p_hat = amg_solver_data.solver->RunVCycle();
	}
	//VecSet(u_p_hat,0.0);
//	VecView(u_p_hat, PETSC_VIEWER_STDOUT_WORLD);
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//printf("Schur1: %lf \n",elapsed_secs);
	if(simpler)
	{
		MatMult(amg_ns_supg_solver_data.Avp, u_p_hat, intermediate_v);
		VecAYPX(intermediate_v, -1.0, *bvp_residual);
	}
	else
	{
		VecCopy(*bvp_residual, intermediate_v);
	}
	//VecView(intermediate_v, PETSC_VIEWER_STDOUT_WORLD);
	GaussSidelMethod* gaussSideMethod = new GaussSidelMethod();
	gaussSideMethod->PreSmoothing(amg_ns_supg_solver_data.Avv, intermediate_v);
	//VecSet(u_v,0.0);
	MatGetSize(amg_ns_supg_solver_data.Avv,&size,NULL);
	ResidualBasedErrorEvaluator* velocity_error_evalutaor =
			new ResidualBasedErrorEvaluator(amg_ns_supg_solver_data.Avv, intermediate_v, 0, 0, size, 4);
	PetscReal norm;
	double vv_start = get_timeClock();
	int iteration_nr = 0;
//	MatView(amg_ns_supg_solver_data.Avv, PETSC_VIEWER_STDOUT_WORLD);
//	VecView(intermediate_v, PETSC_VIEWER_STDOUT_WORLD);
	do
	{
		gaussSideMethod->Smooth(u_v);
//		VecView(u_v, PETSC_VIEWER_STDOUT_WORLD);
	}
	while(!velocity_error_evalutaor->Stop(u_v,iteration_nr++));
//	//		VecSetValues(b,1,&i,&value,INSERT_VALUES);
//	const PetscInt indx[] = {0,1,2,3,4,5};
//	const PetscScalar values[] = {0.0,1.0,2.0,3.0,4.0,5.0};
//	VecSetValues(u_v, 6, indx, values, INSERT_VALUES);
	vv_time += (get_timeClock() - vv_start);
	delete velocity_error_evalutaor;

	MatMult(amg_ns_supg_solver_data.Apv, u_v, intermediate_p);
	VecAYPX(intermediate_p,-1.0,*bpv_residual);
	//VecView(intermediate_p, PETSC_VIEWER_STDOUT_WORLD);

	if(simpler)
	{
		MatMult(amg_solver_data.A, u_p_hat, p_tmp);
		VecAXPY(intermediate_p, -1.0, p_tmp);
	}

	//change_pressure_for_boundary_vec(intermediate_p);
	amg_solver_data.solver->ChangeSystem(intermediate_p, 1);
//	MatView(amg_ns_supg_solver_data.schur_complement, PETSC_VIEWER_STDOUT_WORLD);
//	VecView(intermediate_p, PETSC_VIEWER_STDOUT_WORLD);
//	error_evaluator = new ResidualBasedErrorEvaluator(schur_complement, intermediate_p, 1.0e-8, 1.0e-8, size);
//	solver = new AMGSolverStructure(schur_complement, intermediate_p, size, error_evaluator, false);
//	//clock_t begin = clock();
//	solver->CreateAMGSolverLevels();
	//clock_t end = clock();
	//double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//printf("Total time of level creation: %lf\n",elapsed_secs);
	//time_t start = time(NULL);
	//begin = clock();
	//MatView(schur_complement, PETSC_VIEWER_STDOUT_WORLD);
	VecDestroy(&u_p);
	begin = clock();
	double schur_start = get_timeClock();
	u_p = amg_solver_data.solver->RunVCycle();
	schur_time += (get_timeClock() - schur_start);
	double pressure_dump = 1.0;
//	VecView(u_p, PETSC_VIEWER_STDOUT_WORLD);
	VecScale(u_p, pressure_dump);
	end = clock();
	elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//printf("Schur2: %lf \n",elapsed_secs);
	//VecView(u_p, PETSC_VIEWER_STDOUT_WORLD);
	//end = clock();
	//elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	//printf("Total time of solving: %lf %.2f\n",elapsed_secs, (double)(time(NULL) - start));
	//VecView(u_p,PETSC_VIEWER_STDOUT_WORLD);

	MatMult(amg_ns_supg_solver_data.Avp, u_p, intermediate_v);

	//VecView(intermediate_v,PETSC_VIEWER_STDOUT_WORLD);

	MatMult(X, intermediate_v, v_tmp);
//	MatSolve(X, intermediate_v, v_tmp);

	//VecView(v_tmp,PETSC_VIEWER_STDOUT_WORLD);
	VecCopy(v_tmp, intermediate_v);
	//VecView(intermediate_v, PETSC_VIEWER_STDOUT_WORLD);
	//VecView(u_v, PETSC_VIEWER_STDOUT_WORLD);
	//change_pressure_for_boundary_vec(intermediate_v);
	VecAYPX(intermediate_v, -1.0, u_v);
	//VecView(intermediate_v, PETSC_VIEWER_STDOUT_WORLD);

	VecCopy(intermediate_v,*bvp_residual);
	//commutator
//	MatMult(amg_ns_supg_solver_data.commutator, u_p, *bpv_residual);
		VecCopy(u_p,*bpv_residual);
	//commutator


	if(simpler)
	{
		VecAYPX(*bpv_residual, 1.0, u_p_hat);
		VecDestroy(&u_p_hat);
	}

	//VecView(intermediate_v, PETSC_VIEWER_STDOUT_WORLD);
	//VecView(u_p, PETSC_VIEWER_STDOUT_WORLD);
	delete gaussSideMethod;
}

void simpler_preconditioner(int Ndof, double* X, Vec bvp_residual, Vec bpv_residual)
{
	Vec u_v;
	Vec u_p = NULL;
	Vec intermediate_v;
	Vec intermediate_p;
	Vec v_residual;
	Vec v_tmp;
	Vec p_tmp;
	VecDuplicate(amg_ns_supg_solver_data.bvp,&v_residual);
	VecDuplicate(amg_ns_supg_solver_data.bvp,&intermediate_v);
	VecDuplicate(amg_solver_data.b,&intermediate_p);
	VecDuplicate(amg_ns_supg_solver_data.bvp,&u_v);
	VecDuplicate(amg_ns_supg_solver_data.bvp,&v_tmp);
	VecDuplicate(amg_solver_data.b,&p_tmp);

	//MatCreateSchurComplementPmat(Avv,Avp,Apv,App,MAT_SCHUR_COMPLEMENT_AINV_DIAG,MAT_INITIAL_MATRIX,&schur_complement);
	//MatView(schur_complement, PETSC_VIEWER_STDOUT_WORLD);
	simpler_preconditioner_iteration(0.0000001, &bvp_residual, u_v, v_residual, intermediate_v, v_tmp, &bpv_residual, u_p, intermediate_p, p_tmp, false);

	//simpler_preconditioner_iteration(0.0001, &uv, u_v, v_residual, intermediate_v, vv_diagonal, &up, u_p, intermediate_p, schur_complement);

	VecDestroy(&intermediate_v);
	VecDestroy(&intermediate_p);
	VecDestroy(&v_residual);
	VecDestroy(&u_v);
	VecDestroy(&v_tmp);
	VecDestroy(&p_tmp);

	convert_vec_to_x(bvp_residual, bpv_residual, X, Ndof);
}

void mult_velocity_pressure(Vec current_v, Vec current_p, Vec result_v, Vec result_p)
{
	Vec tmp_v, tmp_p;
	VecDuplicate(current_v,&tmp_v);
	VecDuplicate(current_p,&tmp_p);

	MatMult(amg_ns_supg_solver_data.Avv,current_v,tmp_v);
	MatMultAdd(amg_ns_supg_solver_data.Avp,current_p,tmp_v,result_v);
	MatMult(amg_ns_supg_solver_data.Apv,current_v, tmp_p);
	MatMultAdd(amg_solver_data.A,current_p,tmp_p,result_p);

	VecDestroy(&tmp_v);
	VecDestroy(&tmp_p);
}

void compute_residuals(Vec* residual_bv, Vec* residual_bp, Vec current_v, Vec current_p)
{
	//VecView(bvp,PETSC_VIEWER_STDOUT_WORLD);
	VecDuplicate(amg_ns_supg_solver_data.bvp,residual_bv);
	VecDuplicate(amg_solver_data.b,residual_bp);

	mult_velocity_pressure(current_v, current_p, *residual_bv, *residual_bp);
	VecAYPX(*residual_bv,-1.0,amg_ns_supg_solver_data.bvp);
	VecAYPX(*residual_bp,-1.0,amg_solver_data.b);
}

//void solve(int nrdof, double* v)
//{
//	Mat A;
//	MatCreate(PETSC_COMM_WORLD,&A);
//	MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,nrdof,nrdof);
//	MatSetFromOptions(A);
//	MatMPIAIJSetPreallocation(A,100,NULL,100,NULL);
//	MatSeqAIJSetPreallocation(A,100,NULL);
//	MatSetUp(A);
//
//	PetscInt ncols;
//	const PetscInt *cols;
//	const PetscScalar *vals;
//	int i_index, j_index;
//	double value;
//
//	Vec b;
//	VecCreate(PETSC_COMM_WORLD,&b);
//	VecSetSizes(b,PETSC_DECIDE,nrdof);
//	VecSetFromOptions(b);
//
//	for(int i = 0; i<nrdof * 3 / 4; i++)
//	{
//		MatGetRow(Avv,i,&ncols,&cols,&vals);
//		for(int j = 0; j<ncols; j++)
//		{
//			i_index = i;
//			j_index = cols[j];
//			value = vals[j];
//			MatSetValues(A,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		}
//		MatRestoreRow(Avv,i,&ncols,&cols,&vals);
//		MatGetRow(Avp,i,&ncols,&cols,&vals);
//		for(int j = 0; j<ncols; j++)
//		{
//			i_index = i;
//			j_index = cols[j] + nrdof * 3 / 4;
//			value = vals[j];
//			MatSetValues(A,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		}
//		MatRestoreRow(Avp,i,&ncols,&cols,&vals);
//
//		VecGetValues(bvp,1,&i,&value);
//		VecSetValues(b,1,&i,&value,INSERT_VALUES);
//	}
//	for(int i = 0; i<nrdof * 1 / 4; i++)
//	{
//		MatGetRow(Apv,i,&ncols,&cols,&vals);
//		for(int j = 0; j<ncols; j++)
//		{
//			i_index = i + nrdof * 3 / 4;
//			j_index = cols[j];
//			value = vals[j];
//			MatSetValues(A,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		}
//		MatRestoreRow(Apv,i,&ncols,&cols,&vals);
//		MatGetRow(App,i,&ncols,&cols,&vals);
//		for(int j = 0; j<ncols; j++)
//		{
//			i_index = i + nrdof * 3 / 4;
//			j_index = cols[j] + nrdof * 3 / 4;
//			value = vals[j];
//			MatSetValues(A,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		}
//		MatRestoreRow(App,i,&ncols,&cols,&vals);
//
//		i_index = i + nrdof * 3 / 4;
//		VecGetValues(bpv,1,&i,&value);
//		VecSetValues(b,1,&i_index,&value,INSERT_VALUES);
//	}
//
//	MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
//	MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);
//
//	VecAssemblyBegin(b);
//	VecAssemblyEnd(b);
//
//	Vec solution, residual;
//	VecDuplicate(b, &solution);
//	VecDuplicate(b, &residual);
//	VecSet(solution, 0);
//
//	KSP ksp;
//	KSPCreate(PETSC_COMM_WORLD,&ksp);
//
//
//	KSPSetOperators(ksp,A,A);
//	KSPSetFromOptions(ksp);
//
//
//	KSPSolve(ksp,b,solution);
//	PetscScalar* values;
//	VecGetArray(solution,&values);
//	for(int i = 0; i<nrdof; i+=4)
//	{
//		//v[i] = values[i - i/4];
//		printf("%lf",values[i - i/4]);
//		//v[i+1] = values[i+1 - i/4];
//		printf(" %lf",values[i+1 - i/4]);
//		//v[i+2] = values[i+2 - i/4];
//		printf(" %lf",values[i+2 - i/4]);
//		//v[i+3] = values[nrdof * 3 / 4 + i/4];
//		printf(" %lf\n",values[nrdof * 3 / 4 + i/4]);
//	}
//	MatDestroy(&A);
//	VecDestroy(&b);
//	VecDestroy(&solution);
//}

//TODO: parallel version + vector reusage!
void get_residual_vectors(double* X, int Ndof, Vec* v_residual, Vec* p_residual){
	Vec v;
	Vec p;
	Vec vv_residual;
	Vec vp_residual;
	Vec pv_residual;
	Vec pp_residual;
	PetscErrorCode ierr;
	ierr = VecCreate(PETSC_COMM_WORLD,&v);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecSetSizes(v,PETSC_DECIDE,Ndof * 3 / 4);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecSetFromOptions(v);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecCreate(PETSC_COMM_WORLD,&p);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecSetSizes(p,PETSC_DECIDE,Ndof * 1 / 4);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecSetFromOptions(p);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	VecDuplicate(v, &vp_residual);
	VecDuplicate(v, &vv_residual);
	VecDuplicate(p, &pv_residual);
	VecDuplicate(p, &pp_residual);
	convert_x_to_vec(v,p,X,Ndof);

	ierr = VecAssemblyBegin(v); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecAssemblyEnd(v); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecAssemblyBegin(p); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecAssemblyEnd(p); CHKERRABORT(PETSC_COMM_WORLD, ierr);


	MatResidual(amg_ns_supg_solver_data.Avv,amg_ns_supg_solver_data.bvp,v,vv_residual);
//	VecView(vv_residual, PETSC_VIEWER_STDOUT_WORLD);
	MatResidual(amg_ns_supg_solver_data.Avp,amg_ns_supg_solver_data.bvp,p,vp_residual);
//	VecView(vp_residual, PETSC_VIEWER_STDOUT_WORLD);
	MatResidual(amg_ns_supg_solver_data.Apv,amg_solver_data.b,v,pv_residual);
//	VecView(pv_residual, PETSC_VIEWER_STDOUT_WORLD);
	MatResidual(amg_solver_data.A,amg_solver_data.b,p,pp_residual);
//	VecView(pp_residual, PETSC_VIEWER_STDOUT_WORLD);

	VecAXPY(vv_residual,1.0,vp_residual);
	VecAXPY(pp_residual,1.0,pv_residual);

	VecAXPY(vv_residual,-1.0,amg_ns_supg_solver_data.bvp);
	VecAXPY(pp_residual,-1.0,amg_solver_data.b);
//
//	VecView(vv_residual, PETSC_VIEWER_STDOUT_WORLD);
//	VecView(pp_residual, PETSC_VIEWER_STDOUT_WORLD);

	*v_residual = vv_residual;
	*p_residual = pp_residual;

	VecDestroy(&v);
	VecDestroy(&p);
	VecDestroy(&vp_residual);
	VecDestroy(&pv_residual);
}

double lsr_ns_supg_ext_petsc_compres(int Solver_id, double* X, int Ndof){

	Vec v_residual, p_residual;
	get_residual_vectors(X, Ndof, &v_residual, &p_residual);
	double r1, r2;
	VecNorm(v_residual,NORM_2,&r1);
	VecNorm(p_residual,NORM_2,&r2);
	VecDestroy(&v_residual);
	VecDestroy(&p_residual);
	return sqrt(r1*r1 + r2*r2);
}

void lsr_ns_supg_ext_petsc_compres_vector(
	int Solver_id,   /* in: solver ID (used to identify the subproblem) */
	int Ndof, 	/* in: number of unknowns (components of x) */
	double* X, 	/* in: input vector (may be NULL if Ini_zero==0) */
	double* B,	/* in:  the rhs vector, if NULL take rhs */
					/*      from block data structure ( B is not taken */
					/*      into account if Use_rhs!=1) */
	double* V 	/* out: v = b-Ax */
	)
{
	Vec v_residual, p_residual;
	get_residual_vectors(X, Ndof, &v_residual, &p_residual);
	convert_vec_to_x(v_residual, p_residual, V, Ndof);
	VecDestroy(&v_residual);
	VecDestroy(&p_residual);
}

void lsr_ns_supg_ext_petsc_compreres (int Solver_id, int Subsystem_id, int Level_id, int Control, int Ini_zero,	int Ndof,
		double* X, double* B, double* V)
{	  PetscErrorCode ierr;
//	  Ndof = 8;
//	  ierr = MatCreate(PETSC_COMM_WORLD,&amg_ns_supg_solver_data.Avv);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSetSizes(amg_ns_supg_solver_data.Avv,PETSC_DECIDE,PETSC_DECIDE,Ndof*3/4,Ndof*3/4);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSetFromOptions(amg_ns_supg_solver_data.Avv);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//
//	  ierr = MatCreate(PETSC_COMM_WORLD,&amg_ns_supg_solver_data.Avp);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSetSizes(amg_ns_supg_solver_data.Avp,PETSC_DECIDE,PETSC_DECIDE,Ndof*3/4,Ndof*1/4);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSetFromOptions(amg_ns_supg_solver_data.Avp);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//
//	  ierr = MatCreate(PETSC_COMM_WORLD,&amg_ns_supg_solver_data.Apv);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSetSizes(amg_ns_supg_solver_data.Apv,PETSC_DECIDE,PETSC_DECIDE,Ndof*1/4,Ndof*3/4);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSetFromOptions(amg_ns_supg_solver_data.Apv);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//
//	  ierr = MatCreate(PETSC_COMM_WORLD,&amg_solver_data.A);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSetSizes(amg_solver_data.A,PETSC_DECIDE,PETSC_DECIDE,Ndof*1/4,Ndof*1/4);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSetFromOptions(amg_solver_data.A);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//
//	  ierr = VecCreate(PETSC_COMM_WORLD,&amg_ns_supg_solver_data.bvp);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = VecSetSizes(amg_ns_supg_solver_data.bvp,PETSC_DECIDE,Ndof * 3 / 4);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = VecSetFromOptions(amg_ns_supg_solver_data.bvp);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//
//	  ierr = VecCreate(PETSC_COMM_WORLD,&amg_solver_data.b);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = VecSetSizes(amg_solver_data.b,PETSC_DECIDE,Ndof * 1 / 4);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = VecSetFromOptions(amg_solver_data.b);CHKERRABORT(PETSC_COMM_WORLD, ierr);
////
//	  ierr = MatMPIAIJSetPreallocation(amg_ns_supg_solver_data.Avv,100,NULL,100,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSeqAIJSetPreallocation(amg_ns_supg_solver_data.Avv,100,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSetUp(amg_ns_supg_solver_data.Avv); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatMPIAIJSetPreallocation(amg_ns_supg_solver_data.Avp,100,NULL,100,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSeqAIJSetPreallocation(amg_ns_supg_solver_data.Avp,100,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSetUp(amg_ns_supg_solver_data.Avp); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatMPIAIJSetPreallocation(amg_ns_supg_solver_data.Apv,100,NULL,100,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSeqAIJSetPreallocation(amg_ns_supg_solver_data.Apv,100,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSetUp(amg_ns_supg_solver_data.Apv); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatMPIAIJSetPreallocation(amg_solver_data.A,100,NULL,100,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSeqAIJSetPreallocation(amg_solver_data.A,100,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSetUp(amg_solver_data.A); CHKERRABORT(PETSC_COMM_WORLD, ierr);

//	int i_index, j_index;
//	double value;
//	i_index = 0; j_index = 0; value = 2;
//	MatSetValues(amg_ns_supg_solver_data.Avv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//	i_index = 0; j_index = 1; value = 1;
//	MatSetValues(amg_ns_supg_solver_data.Avv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//	i_index = 0; j_index = 2; value = 1;
//	MatSetValues(amg_ns_supg_solver_data.Avv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//	i_index = 1; j_index = 0; value = 1;
//	MatSetValues(amg_ns_supg_solver_data.Avv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//	i_index = 1; j_index = 1; value = 2;
//	MatSetValues(amg_ns_supg_solver_data.Avv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//	i_index = 1; j_index = 2; value = 1;
//	MatSetValues(amg_ns_supg_solver_data.Avv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//	i_index = 2; j_index = 0; value = 1;
//	MatSetValues(amg_ns_supg_solver_data.Avv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//	i_index = 2; j_index = 1; value = 1;
//	MatSetValues(amg_ns_supg_solver_data.Avv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//	i_index = 2; j_index = 2; value = 2;
//	MatSetValues(amg_ns_supg_solver_data.Avv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//
//	i_index = 0; j_index = 0; value = 1;
//	MatSetValues(amg_ns_supg_solver_data.Avp,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//	i_index = 1; j_index = 0; value = 1;
//	MatSetValues(amg_ns_supg_solver_data.Avp,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//	i_index = 2; j_index = 0; value = 1;
//	MatSetValues(amg_ns_supg_solver_data.Avp,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//
//	i_index = 0; j_index = 0; value = 1;
//	MatSetValues(amg_ns_supg_solver_data.Apv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//	i_index = 0; j_index = 1; value = 1;
//	MatSetValues(amg_ns_supg_solver_data.Apv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//	i_index = 0; j_index = 2; value = 1;
//	MatSetValues(amg_ns_supg_solver_data.Apv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//
//	i_index = 0; j_index = 0; value = 1;
//	MatSetValues(amg_solver_data.A,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//
//	i_index = 0; value = 11;
//	VecSetValues(amg_ns_supg_solver_data.bvp,1,&i_index,&value,INSERT_VALUES);
//	i_index = 1; value = 12;
//	VecSetValues(amg_ns_supg_solver_data.bvp,1,&i_index,&value,INSERT_VALUES);
//	i_index = 2; value = 13;
//	VecSetValues(amg_ns_supg_solver_data.bvp,1,&i_index,&value,INSERT_VALUES);
//	i_index = 0; value = 10;
//	VecSetValues(amg_solver_data.b,1,&i_index,&value,INSERT_VALUES);
//
//	    int i_index, j_index;
//	    double value;
//		i_index = 0; j_index = 0; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Avv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//
//
//		for(int i = 1; i<6; i++)
//		{
//			for(int j = 0; j<6; j++){
//				double value;
//				if(i==j){
//					value = -10;
//					MatSetValues(amg_ns_supg_solver_data.Avv,1,&i,1,&j,&value,INSERT_VALUES);
//				}
//				else{
//					value = 1;
//					MatSetValues(amg_ns_supg_solver_data.Avv,1,&i,1,&j,&value,INSERT_VALUES);
//				}
//			}
//		}
//
//		i_index = 0; j_index = 0; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Apv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 0; j_index = 1; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Apv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 0; j_index = 2; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Apv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 0; j_index = 3; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Apv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 0; j_index = 4; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Apv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 0; j_index = 5; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Apv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 1; j_index = 4; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Apv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 1; j_index = 5; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Apv,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//
//		i_index = 1; j_index = 0; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Avp,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 1; j_index = 1; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Avp,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 2; j_index = 0; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Avp,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 2; j_index = 1; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Avp,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 3; j_index = 0; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Avp,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 3; j_index = 1; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Avp,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 4; j_index = 0; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Avp,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 4; j_index = 1; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Avp,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 5; j_index = 0; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Avp,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 5; j_index = 1; value = 1;
//		MatSetValues(amg_ns_supg_solver_data.Avp,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//
//
//		i_index = 0; j_index = 0; value = -10;
//		MatSetValues(amg_solver_data.A,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 0; j_index = 1; value = 1;
//		MatSetValues(amg_solver_data.A,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 1; j_index = 0; value = 1;
//		MatSetValues(amg_solver_data.A,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//		i_index = 1; j_index = 1; value = 1;
//		MatSetValues(amg_solver_data.A,1,&i_index,1,&j_index,&value,INSERT_VALUES);
//
//
//		i_index = 0; value = 1;
//		VecSetValues(amg_ns_supg_solver_data.bvp,1,&i_index,&value,INSERT_VALUES);
//		i_index = 1; value = 2;
//		VecSetValues(amg_ns_supg_solver_data.bvp,1,&i_index,&value,INSERT_VALUES);
//		i_index = 2; value = 3;
//		VecSetValues(amg_ns_supg_solver_data.bvp,1,&i_index,&value,INSERT_VALUES);
//		i_index = 3; value = 4;
//		VecSetValues(amg_ns_supg_solver_data.bvp,1,&i_index,&value,INSERT_VALUES);
//		i_index = 4; value = 5;
//		VecSetValues(amg_ns_supg_solver_data.bvp,1,&i_index,&value,INSERT_VALUES);
//		i_index = 5; value = 6;
//		VecSetValues(amg_ns_supg_solver_data.bvp,1,&i_index,&value,INSERT_VALUES);
//		i_index = 0; value = 7;
//		VecSetValues(amg_solver_data.b,1,&i_index,&value,INSERT_VALUES);
//		i_index = 1; value = 8;
//		VecSetValues(amg_solver_data.b,1,&i_index,&value,INSERT_VALUES);
//		lsr_ns_supg_ext_petsc_fill_precon(0,0);
//		MatView(amg_ns_supg_solver_data.Avv, PETSC_VIEWER_STDOUT_WORLD);
//		MatView(amg_ns_supg_solver_data.Avp, PETSC_VIEWER_STDOUT_WORLD);
//		MatView(amg_ns_supg_solver_data.Apv, PETSC_VIEWER_STDOUT_WORLD);
//		MatView(amg_solver_data.A, PETSC_VIEWER_STDOUT_WORLD);


	//zzzzzzzzzzzzzzzzzzzzzzzz
	double start_time = get_timeClock();
	//in fact if initial is zero its copied form the X anyway
	if(Ini_zero)
	{
		VecSet(solution_v,0.0);
		VecSet(solution_p,0.0);
	}
	else
	{
		convert_x_to_vec(solution_v, solution_p, X, Ndof);
	}

	ierr = VecAssemblyBegin(solution_v); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecAssemblyEnd(solution_v); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecAssemblyBegin(solution_p); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecAssemblyEnd(solution_p); CHKERRABORT(PETSC_COMM_WORLD, ierr);



//	MatView(amg_ns_supg_solver_data.Avv, PETSC_VIEWER_STDOUT_WORLD);
//	MatView(amg_ns_supg_solver_data.Avp, PETSC_VIEWER_STDOUT_WORLD);
//	MatView(amg_ns_supg_solver_data.Apv, PETSC_VIEWER_STDOUT_WORLD);
//	MatView(amg_solver_data.A, PETSC_VIEWER_STDOUT_WORLD);
//
//	VecView(amg_ns_supg_solver_data.bvp, PETSC_VIEWER_STDOUT_WORLD);
//	VecView(amg_solver_data.b, PETSC_VIEWER_STDOUT_WORLD);

	/*solve(Ndof, V);
	-ksp_rtol 0.0000001
	-ksp_monitor

	MatView(Avv, PETSC_VIEWER_STDOUT_WORLD);
	MatView(Avp, PETSC_VIEWER_STDOUT_WORLD);
	MatView(Apv, PETSC_VIEWER_STDOUT_WORLD);
	MatView(App, PETSC_VIEWER_STDOUT_WORLD);
*/
//	VecView(amg_ns_supg_solver_data.bvp, PETSC_VIEWER_STDOUT_WORLD);
//	VecView(bpv, PETSC_VIEWER_STDOUT_WORLD);

	bool simplest = false;

//	if(Control){
//		Vec conv;
//		VecCreate(PETSC_COMM_WORLD,&conv);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//		VecSetSizes(conv,PETSC_DECIDE,Ndof);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//		VecSetFromOptions(conv);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//		convert_x_to_vec_simple(conv, X, Ndof);
//		PetscViewer viewer;
//		PetscViewerASCIIOpen(PETSC_COMM_WORLD, "/home/damian/comparison/SUPG_SOLUTIONX.data", &viewer);
//		VecView(conv, viewer);
//		VecDestroy(&conv);
//	}

	if(Control)
	{
		Vec residual_bv, residual_bp;
		compute_residuals(&residual_bv, &residual_bp, solution_v, solution_p);
//		VecView(residual_bv,PETSC_VIEWER_STDOUT_WORLD);
		//VecView(residual_bp,PETSC_VIEWER_STDOUT_WORLD);

		printf("ReS %.15lf \n", lsr_ns_supg_ext_petsc_compres(0,X,Ndof));
		if(simplest)
			simplest_preconditioner(Ndof, V, residual_bv, residual_bp);
		else
			simpler_preconditioner(Ndof, V, residual_bv, residual_bp);
		VecDestroy(&residual_bv);
		VecDestroy(&residual_bp);
	}
	else
	{

		clock_t begin = clock();

		mult_velocity_pressure(solution_v, solution_p, result_v, result_p);
		clock_t mid = clock();
		double elapsed = double(mid - begin) / CLOCKS_PER_SEC;
//		printf("Total time of iteration %lf \n",elapsed);
		if(simplest)
			simplest_preconditioner(Ndof, V, result_v, result_p);
		else
			simpler_preconditioner(Ndof, V, result_v, result_p);
		clock_t end = clock();
		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//		printf("Total time of iteration %lf \n",elapsed_secs);
	}
	total_time += (get_timeClock() - start_time);
//	for(int i = 0; i<8; i++)
//		printf("%lf ", V[i]);
//	printf("\n");
}

//for testing
void lsr_ns_supg_ext_petsc_get_system(Mat* Avv, Mat* Avp, Mat* Apv, Mat* App, Vec* v, Vec* p){
	(*Avv) = amg_ns_supg_solver_data.Avv;
	(*Avp) = amg_ns_supg_solver_data.Avp;
	(*Apv) = amg_ns_supg_solver_data.Apv;
	(*App) = amg_solver_data.A;
	(*v) = amg_ns_supg_solver_data.bvp;
	(*p) = amg_solver_data.b;
}
