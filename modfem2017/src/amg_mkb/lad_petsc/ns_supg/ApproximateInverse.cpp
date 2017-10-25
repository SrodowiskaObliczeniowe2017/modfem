#include "ApproximateInverse.hpp"

double t4, t5, t6;
double x4, x5, x6;

void getRows(int column, int range_begin, int range_end, int* columns, int band, int* result_band){

	*result_band = band;
	if(range_end - range_begin < band){
		*result_band = range_end - range_begin;
		band = range_end - range_begin;
	}

	for(int i = 0; i<band; i++){
		columns[i] = column - band / 2 + i;
	}

	while(columns[0] < range_begin)
		for(int i = 0; i<band; i++)
			columns[i]++;

	while(columns[band - 1] > range_end - 1)
		for(int i = 0; i<band; i++)
			columns[i]--;
}

void getRows(Mat Avv, int column, int range_begin, int range_end, int* columns, int band, int* result_band){
	const int* columns_x;
	int columns_number;
	int diagonal_location = 0;

	MatGetRow(Avv,column,&columns_number,&columns_x,NULL);

	if(columns_number <= band){
		(*result_band) = columns_number;
		for(int i = 0; i<columns_number; i++)
			columns[i] = columns_x[i];
		return;
	}
	(*result_band) = band;

	while(columns_x[diagonal_location] != column)
		diagonal_location++;

	if((band - 1) / 2 > diagonal_location){
		for(int i = 0; i<band; i++)
			columns[i] = columns_x[i];
	}
	else if(diagonal_location + (band - 1) / 2 >= columns_number) {
		for(int i = columns_number - 1 ; i> columns_number - 1 - band; i--)
			columns[band - 1 + i - (columns_number - 1)] = columns_x[i];
	}
	else{
		for(int i = diagonal_location - (band - 1) / 2; i <= diagonal_location + (band - 1) / 2; i++)
			columns[i-(diagonal_location - (band - 1) / 2)] = columns_x[i];

	}
	MatRestoreRow(Avv,column,&columns_number,&columns_x,NULL);
}

void getSystem(double* matrix, int band, int size, double* system){
	for(int i = 0; i<band; i++){
		for(int j = i; j<band; j++){
			system[i*band + j] = 0;
			for(int k = 0; k<size; k++){
				system[i*band + j] += matrix[k*band + i]*matrix[k*band + j];
			}
			system[j*band + i] = system[i*band + j];
		}
	}
}

std::map<std::pair<int,int>, double> pairMap;

int xCount = 0;
int yCount = 0;
int columns_buffer[500];
double values_buffer[500];
double getDifferentColumnsProductSymmetric(Mat Avv, int size, int i, int j){
	double product = 0;
//	t4 = time_clock();
	int columns_number;
	const int* columns;
	const double* values;
	MatGetRow(Avv,i,&columns_number,&columns,&values);

	int columns_buffer_number = columns_number;
	memcpy ( columns_buffer, columns, columns_number *sizeof(PetscInt));
	memcpy ( values_buffer, values, columns_number *sizeof(PetscScalar));
//	for(int k = 0; k<columns_number; k++){
//		columns_buffer[k] = columns[k];
//		values_buffer[k] = values[k];
//	}
	MatRestoreRow(Avv,i,&columns_number,&columns,&values);
//	x4 += time_clock() - t4;
//	t5 = time_clock();
	MatGetRow(Avv,j,&columns_number,&columns,&values);
	int columns_buffer_index = 0;
	int columns_index = 0;
	while(columns_index < columns_number && columns_buffer_index < columns_buffer_number){
		if(columns[columns_index] < columns_buffer[columns_buffer_index]){
			columns_index++;
		}
		else if(columns[columns_index] > columns_buffer[columns_buffer_index]){
			columns_buffer_index++;
		}
		else{
			product += values[columns_index]*values_buffer[columns_buffer_index];
			columns_index++;
			columns_buffer_index++;
		}
	}
	MatRestoreRow(Avv,j,&columns_number,&columns,&values);
//	x5 += time_clock() - t5;

	return product;
}

double getDifferentColumnsProduct(Mat Avv, int size, int i, int j){

	PetscErrorCode ierr;
	int* rows;
	const int* columns;
	double* values;
	int columns_number;
	double product = 0;
//	t4 = time_clock();
	MatGetRow(Avv,i,&columns_number,&columns,NULL);
	int* columns_buffer = (int*) malloc(sizeof(int)* columns_number);
	int columns_buffer_number = columns_number;
	for(int k = 0; k<columns_number; k++)
		columns_buffer[k] = columns[k];
	MatRestoreRow(Avv,i,&columns_number,&columns,NULL);

	MatGetRow(Avv,j,&columns_number,&columns,NULL);
	int columns_buffer_index = 0;
	int columns_index = 0;
	int result_index = 0;
	while(columns_index < columns_number && columns_buffer_index < columns_buffer_number){
		if(columns[columns_index] < columns_buffer[columns_buffer_index]){
			columns_index++;
		}
		else if(columns[columns_index] > columns_buffer[columns_buffer_index]){
			columns_buffer_index++;
		}
		else{
			columns_buffer[result_index] = columns_buffer[columns_buffer_index];
			result_index++;
			columns_index++;
			columns_buffer_index++;
		}
	}
	MatRestoreRow(Avv,j,&columns_number,&columns,NULL);

	double* values_i = (double* )malloc(sizeof(double) * result_index);
	double* values_j = (double* )malloc(sizeof(double) * result_index);
	ierr = MatGetValues(Avv,result_index,columns_buffer,1,&i,values_i); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatGetValues(Avv,result_index,columns_buffer,1,&j,values_j); CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	x4 += time_clock() - t4;
//	t5 = time_clock();
	for(int i = 0; i<result_index; i++)
		product += values_i[i]*values_j[i];
//	x5 += time_clock() - t5;
	free(values_i);
	free(values_j);
	free(columns_buffer);
	return product;
}

double getColumnsProduct(Mat Avv, int size, int i, int j, bool symmetric){

	if(i > j)
	{
		int swap = i;
		i = j;
		j = swap;
	}
	std::pair<int, int> pp(i,j);
//	t6 = time_clock();
	std::map<std::pair<int, int>,double>::iterator it = pairMap.find(pp);
	if(it != pairMap.end()){
		yCount++;
		double v = it->second;
//		x6 += time_clock() - t6;
		return v;
	}
	xCount++;
	PetscErrorCode ierr;
	int* rows; //= (int*)malloc(sizeof(int) * size);
	const int* columns;// = (int*)malloc(sizeof(int) * size);
	const double* values;// = (double* )malloc(sizeof(double) * size);
	double values_buffer[200];
	int columns_number;
	double product = 0;

	if(i == j){
//		t4 = time_clock();

		if(symmetric)
		{
			MatGetRow(Avv,i,&columns_number,&columns,&values);
			for(int k = 0; k<columns_number; k++)
				product += values[k]*values[k];
		}
		else{
			MatGetRow(Avv,i,&columns_number,&columns,NULL);
			ierr = MatGetValues(Avv,columns_number,columns,1,&i,values_buffer); CHKERRABORT(PETSC_COMM_WORLD, ierr);
			for(int k = 0; k<columns_number; k++)
				product += values_buffer[k]*values_buffer[k];
		}
//		x4 += time_clock() - t4;
//		t5 = time_clock();

//		x5 += time_clock() - t5;
		if(symmetric){
			MatRestoreRow(Avv,i,&columns_number,&columns,&values);
		}
		else{
			MatRestoreRow(Avv,i,&columns_number,&columns,NULL);
		}

	}
	else{
		if(symmetric){
			product = getDifferentColumnsProductSymmetric(Avv, size, i, j);
		}
		else{
			product = getDifferentColumnsProduct(Avv, size, i, j);
		}
	}
	pairMap.insert(std::pair<std::pair<int, int>, double>(pp,product));
	return product;
}

void getSystem(Mat Avv, int *columns, int band, int size, double* system){

	for(int i = 0; i<band; i++){
		for(int j = i; j<band; j++){
			system[i*band + j] = getColumnsProduct(Avv, size, columns[i], columns[j], false);
			system[j*band + i] = system[i*band + j];
		}
	}
}

void getCommutatorSystem(Mat Apv_vp, int *columns, int band, int size, double* system, Vec tmp1, Vec tmp2){
	for(int i = 0; i<band; i++){
		for(int j = i; j<band; j++){
			MatGetColumnVector(Apv_vp,tmp1,columns[i]);
			MatGetColumnVector(Apv_vp,tmp2,columns[j]);
			VecDot(tmp1,tmp2,&(system[i*band + j]));
			system[j*band + i] = system[i*band + j];
		}
	}
}

void luDecomposition(double* matrix, int size){
	for(int k = 0; k<size-1; k++){
		for(int i = k+1; i<size; i++){
			matrix[k + i*size] /= matrix[k + k*size];
			for(int j = k+1; j<size; j++){
				matrix[j + i*size] -= matrix[k + i*size]*matrix[j + k*size];
			}
		}
	}
}

void bsfr(double* matrix, int size, double* rhs){
	for(int i = 0; i<size; i++){
		for(int j = 0; j<i; j++)
			rhs[i] -= matrix[j + i*size]*rhs[j];
		//rhs[i] /= matrix[i + i*size];
	}

	for(int i = size -1; i>-1; i--){
		for(int j = size -1; j> i; j--)
			rhs[i] -= matrix[j + i*size]*rhs[j];
		rhs[i] /= matrix[i + i*size];
	}
}

void fillApproximateInversePreallocation(int* nnz, int band, int size){
	for(int i = 0; i<size; i++){
		if(band*2 > size)
			nnz[i] = size;
		else
			nnz[i] = 2*band;
	}

//	for(int i = 0; i<band; i++){
//		nnz[i] *= 2;
//		nnz[size - i - 1] *= 2;
//	}

}

Mat getApproximateInverse2(Mat Avv, int band, Mat restricted_matrix){

	PetscErrorCode ierr;
	Mat sparse_approximate_inverse;
	int size;
	MatGetSize(Avv, &size, NULL);
	if(band > size)
		band = size;
	if(restricted_matrix == NULL){
		ierr = MatCreate(PETSC_COMM_WORLD, &sparse_approximate_inverse);CHKERRABORT(PETSC_COMM_WORLD, ierr);
		ierr = MatSetSizes(sparse_approximate_inverse,PETSC_DECIDE,PETSC_DECIDE,size,size);CHKERRABORT(PETSC_COMM_WORLD, ierr);
		ierr = MatSetFromOptions(sparse_approximate_inverse);CHKERRABORT(PETSC_COMM_WORLD, ierr);
		int* nnz = (int*)malloc(sizeof(int) * size);
		fillApproximateInversePreallocation(nnz, band, size);
		ierr = MatSeqAIJSetPreallocation(sparse_approximate_inverse,-1,nnz);CHKERRABORT(PETSC_COMM_WORLD, ierr);
		ierr = MatSetUp(sparse_approximate_inverse); CHKERRABORT(PETSC_COMM_WORLD, ierr);
		free(nnz);
	}
	else{
		MatDuplicate(restricted_matrix,MAT_SHARE_NONZERO_PATTERN, &sparse_approximate_inverse);
	}


	int range_begin, range_end;
	MatGetOwnershipRange(Avv,&range_begin,&range_end);

	double* least_squares_system = (double*)malloc(sizeof(double)*band*band);
	int* rows = (int*)malloc(sizeof(int) * size);
	int* columns = (int*)malloc(sizeof(int) * band);
	double* rhs;// = (double* )malloc(sizeof(double) * band);
	double* values = (double* )malloc(sizeof(double) * band);
	for(int i = 0; i<size; i++)
		rows[i] = i;

	int result_band;
	pairMap.clear();

	for(int i = range_begin; i<range_end; i++){
		if(restricted_matrix == NULL){
			getRows(i,range_begin, range_end, columns, band, &result_band);
		}
		else{
			getRows(restricted_matrix, i,range_begin, range_end, columns, band, &result_band);
		}


		ierr = MatGetValues(Avv,1,&i,result_band,columns,values); CHKERRABORT(PETSC_COMM_WORLD, ierr);

		getSystem(Avv, columns, result_band, size, least_squares_system);
		rhs = values;
		luDecomposition(least_squares_system, result_band);
		bsfr(least_squares_system, result_band, rhs);
		ierr = MatSetValues(sparse_approximate_inverse,result_band,columns,1,&i,rhs, INSERT_VALUES); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	}
	pairMap.clear();
	free(least_squares_system);
	free(rows);
	free(columns);
	free(values);


	MatAssemblyBegin(sparse_approximate_inverse,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(sparse_approximate_inverse,MAT_FINAL_ASSEMBLY);
//	MatView(Avv, PETSC_VIEWER_DRAW_WORLD);
//	MatView(sparse_approximate_inverse, PETSC_VIEWER_DRAW_WORLD);
//	MatView(sparse_approximate_inverse, PETSC_VIEWER_STDOUT_WORLD);
//	draw_matrix(sparse_approximate_inverse, 10);
	return sparse_approximate_inverse;
}

void fill_vec_to_interpolate_to(Vec vec, Vec tmp_vec, Mat Apv_vv, Mat Avp, int column){
	MatGetColumnVector(Avp,tmp_vec,column);
//	MatView(Apv_vv, PETSC_VIEWER_STDOUT_WORLD);
//	VecView(tmp_vec, PETSC_VIEWER_STDOUT_WORLD);
	MatMult(Apv_vv,tmp_vec,vec);
//	MatMultTranspose(Avp,tmp_vec2, vec);
}

void fill_rhs(Vec vec_to_interpolate_to, Vec tmp_vec, Mat Apv_vp, Mat Avp, int column, double* result){
	MatGetColumnVector(Apv_vp,tmp_vec,column);
//	MatMultTranspose(Avp, tmp_vec, tmp_vec2);
	VecDot(vec_to_interpolate_to,tmp_vec,result);
}

Mat getApproximateInverseCommutator(Mat Avv, Mat Avp, int band){

	Mat Apv_vv, Apv_vp;

	MatTransposeMatMult(Avp,Avv,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&Apv_vv);
	MatTransposeMatMult(Avp,Avp,MAT_INITIAL_MATRIX, PETSC_DEFAULT,&Apv_vp);
//	MatView(Apv_vv, PETSC_VIEWER_STDOUT_WORLD);
//	MatView(Apv_vp, PETSC_VIEWER_STDOUT_WORLD);

	PetscErrorCode ierr;
	Mat sparse_approximate_inverse;
	int size;
	MatGetSize(Avp, NULL, &size);
	if(band > size)
		band = size;
	ierr = MatCreate(PETSC_COMM_WORLD, &sparse_approximate_inverse);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatSetSizes(sparse_approximate_inverse,PETSC_DECIDE,PETSC_DECIDE,size,size);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatSetFromOptions(sparse_approximate_inverse);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	int* nnz = (int*)malloc(sizeof(int) * size);
	fillApproximateInversePreallocation(nnz, band, size);
	ierr = MatSeqAIJSetPreallocation(sparse_approximate_inverse,-1,nnz);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	free(nnz);

	ierr = MatSetUp(sparse_approximate_inverse); CHKERRABORT(PETSC_COMM_WORLD, ierr);


	int range_begin, range_end;
	range_begin = 0;
	range_end = size;

	double* least_squares_system = (double*)malloc(sizeof(double)*band*band);
	int* columns = (int*)malloc(sizeof(int) * band);
	double* rhs;// = (double* )malloc(sizeof(double) * band);
	Vec vec_to_interpolate_to, tmp_vec, tmp_vec2;
	ierr = VecCreate(PETSC_COMM_WORLD,&vec_to_interpolate_to);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecSetSizes(vec_to_interpolate_to,PETSC_DECIDE,size);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecSetFromOptions(vec_to_interpolate_to);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecCreate(PETSC_COMM_WORLD,&tmp_vec);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecSetSizes(tmp_vec,PETSC_DECIDE,size*3);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecSetFromOptions(tmp_vec);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	VecDuplicate(vec_to_interpolate_to, &tmp_vec2);

	rhs = (double* )malloc(sizeof(double) * band);
	double* values = (double* )malloc(sizeof(double) * size);

	int result_band;
	pairMap.clear();

	for(int i = range_begin; i<range_end; i++){
		getRows(i,range_begin, range_end, columns, band, &result_band);
		fill_vec_to_interpolate_to(vec_to_interpolate_to, tmp_vec, Apv_vv, Avp, i);
		for(int j = 0; j<result_band; j++){
			fill_rhs(vec_to_interpolate_to, tmp_vec2, Apv_vp, Avp, columns[j], rhs+j);
		}
		getSystem(Apv_vp, columns, result_band, size, least_squares_system);
//		getCommutatorSystem(Apv_vp, columns, result_band, size, least_squares_system, tmp_vec2, vec_to_interpolate_to);
		luDecomposition(least_squares_system, result_band);
		bsfr(least_squares_system, result_band, rhs);
		ierr = MatSetValues(sparse_approximate_inverse,result_band,columns,1,&i,rhs, INSERT_VALUES); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	}
	pairMap.clear();
	VecDestroy(&vec_to_interpolate_to);
	VecDestroy(&tmp_vec);
	VecDestroy(&tmp_vec2);
	free(least_squares_system);
	free(columns);
	free(rhs);
	free(values);
	MatDestroy(&Apv_vv);
	MatDestroy(&Apv_vp);

	MatAssemblyBegin(sparse_approximate_inverse,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(sparse_approximate_inverse,MAT_FINAL_ASSEMBLY);
//	MatView(Avv, PETSC_VIEWER_DRAW_WORLD);
//	MatView(sparse_approximate_inverse, PETSC_VIEWER_DRAW_WORLD);
//	MatView(sparse_approximate_inverse, PETSC_VIEWER_STDOUT_WORLD);
//	draw_matrix(sparse_approximate_inverse, 10);
	return sparse_approximate_inverse;
}
