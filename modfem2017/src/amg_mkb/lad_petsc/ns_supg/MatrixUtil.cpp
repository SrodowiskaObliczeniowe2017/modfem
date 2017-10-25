#include "MatrixUtil.hpp"

Mat create_exact_schur_complement(Mat* schur_complement, Mat Avp, Mat Avv, Mat Apv, Mat App)
{
	MatDestroy(schur_complement);
	int rows, cols;
	MatGetSize(Avp,&rows,&cols);
	Mat t;
	MatCreateDense(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,rows,cols,NULL,&t);
	IS isr, isc; MatFactorInfo info;
	MatGetOrdering(Avv,MATORDERINGNATURAL,&isr,&isc);
	Mat Avv_dense;
	MatConvert(Avv, MATSEQDENSE,MAT_INITIAL_MATRIX,&Avv_dense);
	MatLUFactor(Avv_dense,isr,isc,&info);
	Mat Avp_dense;
	MatConvert(Avp, MATSEQDENSE,MAT_INITIAL_MATRIX,&Avp_dense);
	MatMatSolve(Avv_dense,Avp_dense,t);
	MatMatMult(Apv,t,MAT_INITIAL_MATRIX,1.0,schur_complement);

	Mat App_dense;
	MatConvert(App, MATSEQDENSE,MAT_INITIAL_MATRIX,&App_dense);
	MatAXPY(*schur_complement,-1.0,App_dense,DIFFERENT_NONZERO_PATTERN);
	MatScale(*schur_complement,-1.0);

	MatDestroy(&t);
	MatDestroy(&Avv_dense);
	MatDestroy(&Avp_dense);
	MatDestroy(&App_dense);

	return *schur_complement;
}

//no longer valid - factored matrix instead!
Mat create_velocity_exact_inverse(Mat* velocity_inverse, Mat Avv, Vec bvp)
{
	MatDestroy(velocity_inverse);
	Vec tmp_v;
	VecDuplicate(bvp,&tmp_v);
	int rows, cols;
	MatGetSize(Avv,&rows,&cols);
	Mat B;
	MatCreateDense(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,rows,cols,NULL,&B);
	VecSet(tmp_v, 1.0);
	MatDiagonalSet(B,tmp_v,INSERT_VALUES);
	Mat X;
	MatCreateDense(PETSC_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE,rows,cols,NULL,&X);
	IS isr, isc; MatFactorInfo info;
	MatGetOrdering(Avv,MATORDERINGNATURAL,&isr,&isc);
	Mat Avv_dense;
	MatConvert(Avv, MATSEQDENSE,MAT_INITIAL_MATRIX,&Avv_dense);
	MatLUFactor(Avv_dense,isr,isc,&info);
	MatMatSolve(Avv_dense, B, X);

	MatDestroy(&B);
	MatDestroy(&Avv_dense);
	VecDestroy(&tmp_v);
	return X;
}

Mat create_velocity_ilu_inverse(KSP* preonly_solver, Mat Avv)
{
	Mat ilu_inverse;
	PC ilu_preconditioner;
	KSPCreate(PETSC_COMM_WORLD,preonly_solver);
	KSPSetOperators(*preonly_solver,Avv,Avv);

	KSPGetPC(*preonly_solver,&ilu_preconditioner);
//	PCSetType(ilu_preconditioner,PCSPAI);
	PCSetType(ilu_preconditioner,PCILU);
	PCFactorSetLevels(ilu_preconditioner,0);
	//PCFactorSetUseInPlace(ilu_preconditioner);

	/* PCFactorSetShiftType(prec,MAT_SHIFT_POSITIVE_DEFINITE); */

	KSPSetFromOptions(*preonly_solver);
	KSPSetType(*preonly_solver,KSPPREONLY);
	KSPSetUp(*preonly_solver);

	PCFactorGetMatrix(ilu_preconditioner,&ilu_inverse);

	return ilu_inverse;
}

Mat create_velocity_mass_inverse(Mat Vmm, Vec bvp)
{
  Mat velocity_inverse;
  PetscErrorCode ierr;
  int size;
  MatGetSize(Vmm,&size, NULL);
  ierr = MatCreate(PETSC_COMM_WORLD, &velocity_inverse);CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSetSizes(velocity_inverse,PETSC_DECIDE,PETSC_DECIDE,size,size);CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSetFromOptions(velocity_inverse);CHKERRABORT(PETSC_COMM_WORLD, ierr);

  ierr = MatMPIAIJSetPreallocation(velocity_inverse,1,NULL,0,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSeqAIJSetPreallocation(velocity_inverse,1,NULL);CHKERRABORT(PETSC_COMM_WORLD, ierr);
  ierr = MatSetUp(velocity_inverse); CHKERRABORT(PETSC_COMM_WORLD, ierr);

  Vec mass_diagonal;
  VecDuplicate(bvp,&mass_diagonal);CHKERRABORT(PETSC_COMM_WORLD, ierr);

  MatGetDiagonal(Vmm,mass_diagonal);
  ierr = VecReciprocal(mass_diagonal);CHKERRABORT(PETSC_COMM_WORLD, ierr);
  VecView(mass_diagonal, PETSC_VIEWER_STDOUT_WORLD);
  MatDiagonalSet(velocity_inverse,mass_diagonal,INSERT_VALUES);


  	MatAssemblyBegin(velocity_inverse,MAT_FINAL_ASSEMBLY);
  	MatAssemblyEnd(velocity_inverse,MAT_FINAL_ASSEMBLY);
  	VecDestroy(&mass_diagonal);
	return velocity_inverse;
}

void draw_matrix(Mat mat, double threshold){

	Mat restricted_matrix = restrict_matrix(mat, threshold);
	MatView(restricted_matrix, PETSC_VIEWER_DRAW_WORLD);
	MatDestroy(&restricted_matrix);
}

Mat restrict_matrix(Mat mat, double threshold){
	const int* columns;
	const int* rows;
	const double* values;
	int columns_number;
	Mat pattern;
	PetscErrorCode ierr;

	int size;
	MatGetSize(mat, &size, NULL);
	ierr = MatCreate(PETSC_COMM_WORLD, &pattern);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatSetSizes(pattern,PETSC_DECIDE,PETSC_DECIDE,size,size);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatSetFromOptions(pattern);CHKERRABORT(PETSC_COMM_WORLD, ierr);

	int range_begin, range_end;
	MatGetOwnershipRange(mat,&range_begin,&range_end);

	int nnz[size];
	for(int i = range_begin; i<range_end; i++){
		MatGetRow(mat, i, &columns_number, &columns, &values);
		nnz[i] = 0;
		for(int j = 0; j<columns_number; j++){
			if((values[j] > threshold) || (values[j] < -threshold))
				nnz[i]++;
		}
		MatRestoreRow(mat, i, &columns_number, &columns, &values);
	}

	ierr = MatSeqAIJSetPreallocation(pattern,-1,nnz);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatSetUp(pattern); CHKERRABORT(PETSC_COMM_WORLD, ierr);



	double one = 1.0;
	for(int i = range_begin; i<range_end; i++){
		MatGetRow(mat, i, &columns_number, &columns, &values);
		for(int j = 0; j<columns_number; j++)
			if((values[j] > threshold) || (values[j] < -threshold))
				MatSetValues(pattern,1,&i,1,&(columns[j]),&one,INSERT_VALUES);
		MatRestoreRow(mat, i, &columns_number, &columns, &values);
	}

	MatAssemblyBegin(pattern, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(pattern, MAT_FINAL_ASSEMBLY);
	MatView(pattern, PETSC_VIEWER_DRAW_WORLD);
	return pattern;
}

void print_sym_difference(Mat mat){
	double max = 0;
	double value1, max_diff_v1;
	double value2, max_diff_v2;
	int size;
	MatGetSize(mat, &size, NULL);
	for(int i = 0; i<size; i++){
		for(int j = i+1; j<size; j++){
			MatGetValues(mat,1,&i,1,&j,&value1);
			MatGetValues(mat,1,&j,1,&i,&value2);
			if(fabs(value1 - value2) > max){
				max = fabs(value1 - value2);
				max_diff_v1 = value1;
				max_diff_v2 = value2;
			}
		}
	}
	printf("MAX DIFFERENCE %.10lf for %.10lf and %.10lf\n", max, max_diff_v1, max_diff_v2);
}

void print_positive_to_negative(Mat mat){
	const int* rows;
	const double* values;
	int columns_number;
	int size;
	MatGetSize(mat,&size, NULL);
	int positive = 0;
	int negative = 0;
	for(int i = 0; i< size; i++){
		MatGetRow(mat, i, &columns_number, NULL, &values);
		for(int j = 0; j<columns_number; j++){
			if(values[j] > 0)
				positive++;
			else
				negative++;
		}
		MatRestoreRow(mat, i, &columns_number, NULL, &values);
	}
	printf("POSITIVE: %d NEGATIVE %d\n", positive, negative);
}
