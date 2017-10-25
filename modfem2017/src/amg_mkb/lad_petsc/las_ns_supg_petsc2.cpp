//#include "las_ns_supg_petsc.hpp"
//
//
//struct AMGSUPGSolverData{
//	Mat A;
//	Vec b;
////	Vec bpv;
////	AMGSolverStructure* solver = NULL;
//	int* posglobs;
//	KSP preonly_solver = NULL;
//	Mat velocity_inverse = NULL;
//} amg_ns_supg_solver_data;
//
//int lsr_ns_supg_ext_petsc_init( /* returns: >0 - solver ID, <0 - error code */
//	int Solver_id,   /* in: solver ID (used to identify the subproblem) */
//	int Parallel,      /* parameter specifying sequential (LSC_SEQUENTIAL) */
//	/* or parallel (LSC_PARALLEL) execution */
//	int* Max_num_levels_p,  /* in: number of levels for multigrid: */
//	/*     1 - enforce single level solver */
//	/*     >1 - enforce the number of levels */
//	/* out: actual number of levels !!! */
//	char* Filename,  /* in: name of the file with control parameters */
//	int Max_iter, /* maximal number of iterations, -1 for values from Filename */
//	int Error_type, /* type of error norm (stopping criterion), -1 for Filename*/
//	double Error_tolerance, /* value for stopping criterion, -1.0 for Filename */
//	int Monitoring_level /* Level of output, -1 for Filename */)
//	{
//		return 0;
//	}
//
//void create_ns_supg_matrix(Mat* matrix, PetscInt rows_number, PetscInt columns_number)
//{
//	  PetscErrorCode ierr;
//	  ierr = MatCreate(PETSC_COMM_WORLD, matrix);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSetSizes(*matrix,PETSC_DECIDE,PETSC_DECIDE,rows_number,columns_number);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSetFromOptions(*matrix);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//}
//
//void create_ns_supg_vector(Vec* vector, PetscInt rows_number)
//{
//	  PetscErrorCode ierr;
//	  ierr = VecCreate(PETSC_COMM_WORLD,vector);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = VecSetSizes(*vector,PETSC_DECIDE,rows_number);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = VecSetFromOptions(*vector);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//}
//
//int velocity_size;
//Vec solution;
//Vec result;
//
//int lsr_ns_supg_ext_petsc_create_matrix(
//	/* returns: >=0 - success code, <0 - error code */
//	int Solver_id,   /* in: solver ID (used to identify the subproblem) */
//	int Level_id,    /* in: level ID */
//	int Nrblocks,    /* in: number of DOF blocks */
//	int Nrdof_glob,  /* in: total number of DOFs */
//	int Max_sm_size, /* in: maximal size of the stiffness matrix */
//	int* Nrdofbl,	   /* in: list of numbers of dofs in a block */
//	int* Posglob,	   /* in: list of global numbers of first dof */
//	int* Nroffbl,	   /* in: list of numbers of off diagonal blocks */
//	int** L_offbl	   /* in: list of lists of off diagonal blocks */
//	)
//	{
//	velocity_size = 3*Nrdof_glob/4;
//	  amg_ns_supg_solver_data.posglobs = (int*)malloc(sizeof(int) * Nrblocks);
//	  memcpy(amg_ns_supg_solver_data.posglobs, Posglob,sizeof(int)* Nrblocks);
//
//	  PetscBool initialized;
//	  PetscInitialized(&initialized);
//	  if(!initialized)
//		  PetscInitialize(NULL,NULL,NULL,NULL);
//
//	  // ERROR in 3.7.4 - need 2 arguments - KB
//	  //PetscOptionsView(PETSC_VIEWER_STDOUT_WORLD);
//
//	  create_ns_supg_matrix(&(amg_ns_supg_solver_data.A), Nrdof_glob, Nrdof_glob);
//	  create_ns_supg_vector(&(amg_ns_supg_solver_data.b), Nrdof_glob);
//
//
//	  int block_index;
//	  int max_nonzeros_in_row = 0;
//	  for(block_index=0;block_index<Nrblocks;block_index++)
//	  {
//	      int diagonal_block_size = Nrdofbl[block_index];
//	      int off_diagonal_block_sizes = 0;
//	      if(Nroffbl[block_index]>0)
//	      {
//	    	  int i;
//	    	  for(i=0;i<Nroffbl[block_index];i++)
//	    	  {
//	    		  int offdiagonal_block_id = L_offbl[block_index][i];
//	    		  off_diagonal_block_sizes += Nrdofbl[offdiagonal_block_id];
//	    	  }
//	      }
//	      max_nonzeros_in_row = std::max(max_nonzeros_in_row, diagonal_block_size + off_diagonal_block_sizes);
//	  }
//
//	  PetscErrorCode ierr;
//	  //TODO: set the second parameter
//	  ierr = MatMPIAIJSetPreallocation(amg_ns_supg_solver_data.A,max_nonzeros_in_row,NULL,0,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSeqAIJSetPreallocation(amg_ns_supg_solver_data.A,max_nonzeros_in_row,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	  ierr = MatSetUp(amg_ns_supg_solver_data.A); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//
//	  return 0;
//	}
//
//int lsr_ns_supg_ext_petsc_assemble_local_sm(
//	/* returns: >=0 - success code, <0 - error code */
//	int Solver_id,         /* in: solver ID (used to identify the subproblem) */
//	int Level_id,          /* in: level ID */
//	int Comp_type,         /* in: indicator for the scope of computations: */
//	/*   NS_SUPG_EXT_SOLVE - solve the system */
//	/*   NS_SUPG_EXT_RESOLVE - resolve for the new rhs vector */
//	int Nr_dof_bl,         /* in: number of global dof blocks */
//	/*     associated with the local stiffness matrix */
//	int* L_bl_id,          /* in: list of dof blocks' IDs */
//	int* L_bl_nrdof,       /* in: list of blocks' numbers of dof */
//	double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
//	double* Rhs_vect,      /* in: rhs vector */
//	char* Rewr_dofs         /* in: flag to rewrite or sum up entries */
//	/*   'T' - true, rewrite entries when assembling */
//	/*   'F' - false, sum up entries when assembling */
//	)
//	{
//
//	//	printf("%d\n",first_ghost_block);
//
//		PetscErrorCode ierr;
//		int i;
//		int j;
//		int all_blocks_height = 0;
//		for(i = 0; i<Nr_dof_bl; i++)
//			all_blocks_height += L_bl_nrdof[i];
//
//		int posglob_column_buffer[1000];
//
//		if(all_blocks_height > 1000){
//			printf("POSGLOB BUFFER MAX IS 1000");
//			CHKERRABORT(PETSC_COMM_WORLD, 17);
//		}
//
//		int column_index = 0;
//		for(i = 0; i<Nr_dof_bl; i++){
//			int posglob = amg_ns_supg_solver_data.posglobs[L_bl_id[i]];
//			for(j = 0; j<L_bl_nrdof[i]; j++){
//				posglob_column_buffer[column_index++] = posglob + j;
//			}
//		}
//
//		//transpose
//		for(i = 0; i<all_blocks_height; i++){
//			for(j = i+1; j<all_blocks_height; j++){
//				double tmp = Stiff_mat[i*all_blocks_height + j];
//				Stiff_mat[i*all_blocks_height + j] = Stiff_mat[j*all_blocks_height + i];
//				Stiff_mat[j*all_blocks_height + i] = tmp;
//			}
//		}
//
//		int row_index = 0;
//		for(i = 0; i<Nr_dof_bl; i++){
//	//		int posglob = amg_solver_data.posglobs[L_bl_id[i]];
//				ierr = MatSetValues(amg_ns_supg_solver_data.A,L_bl_nrdof[i],posglob_column_buffer + row_index,
//						column_index,posglob_column_buffer,Stiff_mat + all_blocks_height*row_index,ADD_VALUES);
//				CHKERRABORT(PETSC_COMM_WORLD, ierr);
//
//				ierr = VecSetValues(amg_ns_supg_solver_data.b,L_bl_nrdof[i],posglob_column_buffer + row_index,
//						&(Rhs_vect[row_index]),ADD_VALUES);
//				CHKERRABORT(PETSC_COMM_WORLD, ierr);
//
//			row_index += L_bl_nrdof[i];
//		}
//
////	if(Comp_type != 0 && Comp_type != 1){
////		return 1;
////	}
////
////	PetscErrorCode ierr;
////	int i;
////	int j;
////	int all_blocks_height = 0;
////	for(i = 0; i<Nr_dof_bl; i++)
////		all_blocks_height += L_bl_nrdof[i];
////
////	int row_offset = 0;
////	for(i = 0; i<Nr_dof_bl; i++)
////	{
////		int column_offset = 0;
////		for(j = 0; j<Nr_dof_bl; j++)
////		{
////			int k;
////			int l;
////
////			//Avv
////			int j_index = amg_ns_supg_solver_data.posglobs[L_bl_id[j]] * 3 / 4;
////			for(k = 0; k<3; k++)
////			{
////				int i_index = amg_ns_supg_solver_data.posglobs[L_bl_id[i]] * 3 / 4;
////				for(l = 0; l<3; l++)
////				{
////					double value = Stiff_mat[(column_offset+k)*all_blocks_height + row_offset +l];
////					//printf("(%d %d) v %lf ",i_index,j_index,value);
////					ierr = MatSetValues(amg_ns_supg_solver_data.A,1,&i_index,1,&j_index,&value,ADD_VALUES); CHKERRABORT(PETSC_COMM_WORLD, ierr);
////					i_index++;
////				}
////				j_index++;
////				//printf("\n");
////			}
////
////			//Avp
////			j_index = velocity_size + amg_ns_supg_solver_data.posglobs[L_bl_id[j]] * 1 / 4;
////			for(k = 0; k<1; k++)
////			{
////				int i_index = amg_ns_supg_solver_data.posglobs[L_bl_id[i]] * 3 / 4;
////				for(l = 0; l<3; l++)
////				{
////					double value = Stiff_mat[(column_offset+k+3)*all_blocks_height + row_offset +l];
////					//printf("(%d %d) v %lf ",i_index,j_index,value);
////					ierr = MatSetValues(amg_ns_supg_solver_data.A,1,&i_index,1,&j_index,&value,ADD_VALUES); CHKERRABORT(PETSC_COMM_WORLD, ierr);
////					i_index++;
////				}
////				j_index++;
////				//printf("\n");
////			}
////
////			//Apv
////			j_index = amg_ns_supg_solver_data.posglobs[L_bl_id[j]] * 3 / 4;
////			for(k = 0; k<3; k++)
////			{
////				int i_index = velocity_size + amg_ns_supg_solver_data.posglobs[L_bl_id[i]] * 1 / 4;
////				for(l = 0; l<1; l++)
////				{
////					double value = Stiff_mat[(column_offset+k)*all_blocks_height + row_offset + l + 3];
////					//printf("(%d %d) v %lf ",i_index,j_index,value);
////					ierr = MatSetValues(amg_ns_supg_solver_data.A,1,&i_index,1,&j_index,&value,ADD_VALUES); CHKERRABORT(PETSC_COMM_WORLD, ierr);
////					i_index++;
////				}
////				j_index++;
////				//printf("\n");
////			}
////
////			//App
////			j_index = velocity_size + amg_ns_supg_solver_data.posglobs[L_bl_id[j]] * 1 / 4;
////			for(k = 0; k<1; k++)
////			{
////				int i_index = velocity_size + amg_ns_supg_solver_data.posglobs[L_bl_id[i]] * 1 / 4;
////				for(l = 0; l<1; l++)
////				{
////					double value = Stiff_mat[(column_offset+k+3)*all_blocks_height + row_offset + l + 3];
////					//printf("(%d %d) v %lf ",i_index,j_index,value);
////					ierr = MatSetValues(amg_ns_supg_solver_data.A,1,&i_index,1,&j_index,&value,ADD_VALUES); CHKERRABORT(PETSC_COMM_WORLD, ierr);
////					i_index++;
////				}
////				j_index++;
////				//printf("\n");
////			}
////			column_offset += L_bl_nrdof[j];
////		}
////		int r;
////		for(r = 0; r<3; r++)
////		{
////			//printf("rhs %d wartosc %lf\n",r + posglobs[L_bl_id[i]] * 3 / 4,Rhs_vect[r+row_offset]);
////			int row = r + amg_ns_supg_solver_data.posglobs[L_bl_id[i]] * 3 / 4;
////			ierr = VecSetValues(amg_ns_supg_solver_data.b,1,&row,&(Rhs_vect[r+row_offset]),ADD_VALUES); CHKERRABORT(PETSC_COMM_WORLD, ierr);
////		}
////		for(r = 0; r<1; r++)
////		{
////			//printf("rhs %d wartosc %lf\n",r + posglobs[L_bl_id[i]] * 1 / 4,Rhs_vect[r+row_offset + 3]);
////			int row = velocity_size + r + amg_ns_supg_solver_data.posglobs[L_bl_id[i]] * 1 / 4;
////			ierr = VecSetValues(amg_ns_supg_solver_data.b,1,&row,&(Rhs_vect[r+row_offset+3]),ADD_VALUES); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//////			app_local_rhs[row_offset/4] = Rhs_vect[r+row_offset+3];
////		}
////		row_offset += L_bl_nrdof[i];
////	}
//
//  return(1);
//}
//
//int lsr_ns_supg_ext_petsc_solve(
//	int Solver_id,
//	int Ndof,
//	int Ini_zero,
//	double* X,
//	double* B,
//	int* Nr_iter,
//	double* Toler,
//	int Monitor,
//	double* Conv_rate
//	){
//	return 0;
//}
//
//void scale_system(){
//
//	int size;
//	MatGetSize(amg_ns_supg_solver_data.A, &size, NULL);
//	Vec scaling_vec;
//	create_ns_supg_vector(&scaling_vec, size);
//
//	double value = -1;
//	VecSet(scaling_vec, 1.0);
//	for(int i = size * 3 / 4; i<size; i++){
//		VecSetValues(scaling_vec, 1, &i, &value, INSERT_VALUES);
//	}
//	VecAssemblyBegin(amg_ns_supg_solver_data.b);
//	VecAssemblyEnd(amg_ns_supg_solver_data.b);
//	VecPointwiseMult(amg_ns_supg_solver_data.b, amg_ns_supg_solver_data.b,scaling_vec);
//	MatDiagonalScale(amg_ns_supg_solver_data.A, scaling_vec,NULL);
//}
//
//int lsr_ns_supg_ext_petsc_fill_precon(
//	/* returns: >=0 - success code, <0 - error code */
//	int Solver_id,         /* in: solver ID (used to identify the subproblem) */
//	int Level_id           /* in: level ID */
//	){
//
////	  lar_petsc_fill_preconditioner(0);
//
//	PetscErrorCode ierr;
//	amg_solver_data.amgInitData = amg_algorithm_data;
//
//	ierr = MatAssemblyBegin(amg_ns_supg_solver_data.A,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	ierr = MatAssemblyEnd(amg_ns_supg_solver_data.A,MAT_FINAL_ASSEMBLY); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	ierr = VecAssemblyBegin(amg_ns_supg_solver_data.b); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	ierr = VecAssemblyEnd(amg_ns_supg_solver_data.b); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//
//
//	clock_t begin = clock();
//
//	VecDuplicate(amg_ns_supg_solver_data.b, &solution);
//	VecDuplicate(amg_ns_supg_solver_data.b, &result);
//
//	amg_ns_supg_solver_data.velocity_inverse = create_velocity_ilu_inverse(&(amg_ns_supg_solver_data.preonly_solver), amg_ns_supg_solver_data.A);
//
////	scale_system();
////	print_positive_to_negative(amg_ns_supg_solver_data.A);
//
////	int size;
////	MatGetSize(amg_ns_supg_solver_data.A, &size, NULL);
////	amg_solver_data.solver = new AMGSolverStructure(amg_ns_supg_solver_data.A, NULL, size, NULL, false);
//////	MatView(amg_ns_supg_solver_data.schur_complement, PETSC_VIEWER_STDOUT_WORLD);
////
////	amg_solver_data.solver->CreateAMGSolverLevels(amg_solver_data.amgInitData.coarsening_algorithm,
////			amg_solver_data.amgInitData.interpolation_algorithm,
////			amg_solver_data.amgInitData.strength_threshold,
////			amg_solver_data.amgInitData.levels_number);
//
//
//
//	clock_t end = clock();
//	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
//	printf("Subsystems creation time %lf \n",elapsed_secs);
//
//	return 0;
//}
//
//int lsr_ns_supg_ext_petsc_free_matrix(int Solver_id)
//{
//	MatDestroy(&(amg_ns_supg_solver_data.A));
//	VecDestroy(&(amg_ns_supg_solver_data.b));
//	VecDestroy(&solution);
//	VecDestroy(&result);
//	return 0;
//}
//
//int lsr_ns_supg_ext_petsc_clear_matrix(int Solver_id, int Level_id, int Comp_type)
//{
//	MatZeroEntries(amg_ns_supg_solver_data.A);
//	VecSet(amg_ns_supg_solver_data.b, 0.0);
//	delete amg_solver_data.solver;
//	if(amg_ns_supg_solver_data.preonly_solver != NULL)
//	{
//		KSPDestroy(&(amg_ns_supg_solver_data.preonly_solver));
//	}
////	if(amg_ns_supg_solver_data.velocity_inverse != NULL)
////	{
////		MatDestroy(&(amg_ns_supg_solver_data.velocity_inverse));
////	}
//
//	return 0;
//}
//
//double lsr_ns_supg_ext_petsc_compres(int Solver_id, double* X, int Ndof){
//	Vec solution;
//
//	Vec residual;
//
//	PetscErrorCode ierr;
//	ierr = VecCreate(PETSC_COMM_WORLD,&solution);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	ierr = VecSetSizes(solution,PETSC_DECIDE,Ndof);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	ierr = VecSetFromOptions(solution);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//
//	VecDuplicate(solution, &residual);
//	convert_x_to_vec_simple(solution,X,Ndof);
//
//	ierr = VecAssemblyBegin(solution); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	ierr = VecAssemblyEnd(solution); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//
//	MatResidual(amg_ns_supg_solver_data.A,amg_ns_supg_solver_data.b,solution,residual);
//
//	double r;
//	VecNorm(residual,NORM_2,&r);
//
//	VecDestroy(&solution);
//	VecDestroy(&residual);
//	return r;
//}
//
//void compute_residuals(Vec* residual, Vec current_solution)
//{
//	//VecView(bvp,PETSC_VIEWER_STDOUT_WORLD);
//	VecDuplicate(amg_ns_supg_solver_data.b,residual);
//	MatMult(amg_ns_supg_solver_data.A, current_solution, *residual);
//	VecAYPX(*residual,-1.0,amg_ns_supg_solver_data.b);
//}
//
//void preconditioner(int Ndof, double* V, Vec vector){
//
//
////	VecScale(vector, -1.0);
////	amg_solver_data.solver->ChangeSystem(vector,1);
////	Vec result = amg_solver_data.solver->RunVCycle();
//	Vec result;
//	VecDuplicate(vector, &result);
//	MatSolve(amg_ns_supg_solver_data.velocity_inverse, vector, result);
//	convert_vec_to_x_simple(result, V, Ndof);
//	VecDestroy(&result);
//
//
////	Vec x;
////	VecDuplicate(vector, &x);
////	GaussSidelMethod* gaussSidelMethod = new GaussSidelMethod();
////	gaussSidelMethod->PreSmoothing(amg_ns_supg_solver_data.A, vector);
////	ErrorEvaluator* velocity_error_evaluator =
////			new ResidualBasedErrorEvaluator(amg_ns_supg_solver_data.A, vector, 0, 0, -1, 10);
////	//new ResidualBasedErrorEvaluator(Avv, bvp_residual, 1.0e-8, 1.0e-7, size);
////	clock_t begin = clock();
////	int iteration_nr = 0;
////	PetscReal norm;
////	do
////	{
////		gaussSidelMethod->Smooth(x);
////	}
////	while(!velocity_error_evaluator->Stop(x,iteration_nr++));
////	delete velocity_error_evaluator;
////
////	convert_vec_to_x(x, V, Ndof);
////	VecDestroy(&x);
//}
//
//void lsr_ns_supg_ext_petsc_compreres (int Solver_id, int Subsystem_id, int Level_id, int Control, int Ini_zero,	int Ndof,
//		double* X, double* B, double* V)
//{	  PetscErrorCode ierr;
//
//
//	double total_time = 0;
//	//zzzzzzzzzzzzzzzzzzzzzzzz
//	double start_time = time_clock();
//	//in fact if initial is zero its copied form the X anyway
//	if(Ini_zero)
//	{
//		VecSet(solution,0.0);
//	}
//	else
//	{
//		convert_x_to_vec_simple(solution, X, Ndof);
//	}
//
//	ierr = VecAssemblyBegin(solution); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	ierr = VecAssemblyEnd(solution); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//
//	if(Control){
//		Vec conv;
//		VecDuplicate(solution, &conv);
//		convert_x_to_vec_simple(conv, X, Ndof);
//		PetscViewer viewer;
//		PetscViewerASCIIOpen(PETSC_COMM_WORLD, "/home/damian/comparison/SUPG_SOLUTIONX.data", &viewer);
//		VecView(conv, viewer);
//		VecDestroy(&conv);
//	}
//
//	if(Control)
//	{
//		Vec residual;
//		compute_residuals(&residual, solution);
////		VecView(residual_bv,PETSC_VIEWER_STDOUT_WORLD);
//		//VecView(residual_bp,PETSC_VIEWER_STDOUT_WORLD);
//
//		printf("ReS %.15lf \n", lsr_ns_supg_ext_petsc_compres(0,X,Ndof));
//		preconditioner(Ndof, V, residual);
//		VecDestroy(&residual);
//	}
//	else
//	{
//
//		clock_t begin = clock();
//		MatMult(amg_ns_supg_solver_data.A, solution, result);
//		clock_t mid = clock();
//		double elapsed = double(mid - begin) / CLOCKS_PER_SEC;
//		preconditioner(Ndof, V, result);
//		clock_t end = clock();
//		double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
////		printf("Total time of iteration %lf \n",elapsed_secs);
//	}
//	total_time += (time_clock() - start_time);
////	for(int i = 0; i<8; i++)
////		printf("%lf ", V[i]);
////	printf("\n");
//}
//
//
