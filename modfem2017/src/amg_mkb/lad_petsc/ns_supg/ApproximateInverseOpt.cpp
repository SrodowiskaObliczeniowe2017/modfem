#include "ApproximateInverseOpt.hpp"

double* row_dot_products;
Mat dot_product_matrix;

#define MAX_ROW_SIZE 100

double get_time(){
	return time_clock();
}

int max(int a, int b){
	if(a > b)
		return a;
	return b;
}

int min(int a, int b){
	if(a > b)
		return b;
	return a;
}

void fill_row_dot_products(Mat Avv, double* row_dot_products, int index, int first_column, int last_column){

	const int* columns;
	const double* values;
	int columns_number;
	MatGetRow(Avv,last_column,&columns_number,&columns,&values);

	double main_values[MAX_ROW_SIZE];
	int main_columns[MAX_ROW_SIZE];
	int main_columns_number = columns_number;

	memcpy ( main_columns, columns, main_columns_number *sizeof(PetscInt));
	memcpy ( main_values, values, main_columns_number *sizeof(PetscScalar));

	MatRestoreRow(Avv,last_column,&columns_number,&columns,&values);

	for(int i = first_column; i<= last_column; i++){
		MatGetRow(Avv,i,&columns_number,&columns,&values);

		double row_dot_product = 0;
		int main_columns_index = 0;
		int columns_index = 0;
		while(columns_index < columns_number && main_columns_index < main_columns_number){
			if(columns[columns_index] < main_columns[main_columns_index]){
				columns_index++;
			}
			else if(columns[columns_index] > main_columns[main_columns_index]){
				main_columns_index++;
			}
			else{
				row_dot_product += values[columns_index]*main_values[main_columns_index];
				columns_index++;
				main_columns_index++;
			}
		}
		row_dot_products[index] = row_dot_product;
		index++;
		MatRestoreRow(Avv,i,&columns_number,&columns,&values);
	}


}

int getAuxilaryMin(int i, int band){
	int x = min(i, band);
	return x*(x+1)/2;
}

void create_row_dot_products(Mat Avv, int band, int size){

	int row_dot_products_size = (1+band)*band / 2 + (size - band)*band;
	row_dot_products = (double*)malloc(row_dot_products_size*sizeof(double));
	for(int i = 0; i<size; i++){
		int index = getAuxilaryMin(i, band) + max((i-band)*band, 0);
		fill_row_dot_products(Avv, row_dot_products, index, max(0, i - band + 1), i);
	}
}

double get_row_dot_product(int i, int j, int band){
	if(j > i){
		int swap = j;
		j = i;
		i = swap;
	}
	int index = getAuxilaryMin(i, band) + max((i-band)*band, 0) + j - max(0, i - band + 1);
	return row_dot_products[index];
}

void fillApproximateInversePreallocation2(int* nnz, int band, int size){
	for(int i = 0; i<size; i++){
		if(band*2 > size)
			nnz[i] = size;
		else
			nnz[i] = 2*band;
	}
}

void getRows2(int column, int range_begin, int range_end, int* columns, int band){

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

void getRowsPattern(Mat Avv, int column, int range_begin, int range_end, int* columns, int* band){

	const int* columns_x;
	int columns_number;
	int diagonal_location = 0;

	MatGetRow(Avv,column,&columns_number,&columns_x,NULL);

	if(columns_number <= (*band)){
		(*band) = columns_number;
		for(int i = 0; i<columns_number; i++)
			columns[i] = columns_x[i];
		return;
	}

	while(columns_x[diagonal_location] != column)
		diagonal_location++;

	if(((*band) - 1) / 2 > diagonal_location){
		for(int i = 0; i<(*band); i++)
			columns[i] = columns_x[i];
	}
	else if(diagonal_location + ((*band) - 1) / 2 >= columns_number) {
		for(int i = columns_number - 1 ; i> columns_number - 1 - (*band); i--)
			columns[(*band) - 1 + i - (columns_number - 1)] = columns_x[i];
	}
	else{
		for(int i = diagonal_location - ((*band) - 1) / 2; i <= diagonal_location + ((*band) - 1) / 2; i++)
			columns[i-(diagonal_location - ((*band) - 1) / 2)] = columns_x[i];

	}
	MatRestoreRow(Avv,column,&columns_number,&columns_x,NULL);
}

double tt3 =0;

void solve_lapack(double* matrix, int size, double *rhs){
	double x = get_time();
	int pivots[MAX_ROW_SIZE];
	int info;
	dgetrf_( &size, &size, matrix, &size, pivots, &info);
	char trans = 'C';
	int nrhs = 1;
	dgetrs_(&trans,&size,&nrhs,matrix,&size,pivots,rhs,&size,&info);
	tt3 += get_time() -x;
}

void getSystem2(int *columns, int band, double* system){

	for(int i = 0; i<band; i++){
		for(int j = i; j<band; j++){
//			MatGetValues(dot_product_matrix,1,&(columns[i]),1,&(columns[j]),&(system[i*band + j]));
			system[i*band + j] = get_row_dot_product(columns[i], columns[j], band);
			system[j*band + i] = system[i*band + j];
		}
	}
}

double tt2 = 0;

void getSystem3(int *columns, int band, double* system){

	double x = get_time();
//	for(int i = 0; i<band; i++){
//		for(int j = i; j<band; j++){
			MatGetValues(dot_product_matrix,band,columns,band,columns,system);
//			system[i*band + j] = get_row_dot_product(columns[i], columns[j], band);
//			system[j*band + i] = system[i*band + j];
//		}
//	}

	tt2 += get_time() - x;
}



Mat getApproximateInverseOpt(Mat Avv, int band){

	double tt = get_time();
	tt2 = 0;
	tt3 = 0;
	MatTransposeMatMult(Avv,Avv,MAT_INITIAL_MATRIX,PETSC_DEFAULT,&dot_product_matrix);
//	MatView(dot_product_matrix, PETSC_VIEWER_STDOUT_WORLD);;
//	MatTranspose(Avv,MAT_REUSE_MATRIX,&Avv);
//
	PetscErrorCode ierr;
	Mat sparse_approximate_inverse;
	int size;
	MatGetSize(Avv, &size, NULL);
	if(band > size)
		band = size;
	if(band % 2 == 0)
		band--;

	ierr = MatCreate(PETSC_COMM_WORLD, &sparse_approximate_inverse);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatSetSizes(sparse_approximate_inverse,PETSC_DECIDE,PETSC_DECIDE,size,size);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatSetFromOptions(sparse_approximate_inverse);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	int* nnz = (int*)malloc(sizeof(int) * size);
//	fillApproximateInversePreallocation2(nnz, band, size);
//	ierr = MatSeqAIJSetPreallocation(sparse_approximate_inverse,-1,nnz);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	ierr = MatSetUp(sparse_approximate_inverse); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	MatDuplicate(Avv,MAT_SHARE_NONZERO_PATTERN,&sparse_approximate_inverse);
//	free(nnz);
//
//	create_row_dot_products(Avv, band, size);
	printf("FIRST PART %lf\n", get_time() - tt);
	int range_begin, range_end;
	MatGetOwnershipRange(Avv,&range_begin,&range_end);

	double* least_squares_system = (double*)malloc(sizeof(double)*band*band);
	int* columns = (int*)malloc(sizeof(int) * band);
	double* rhs = (double* )malloc(sizeof(double) * band);
	int original_band = band;
	for(int i = range_begin; i<range_end; i++){
		band = original_band;
//		getRows2(i,range_begin, range_end, columns, band);
		getRowsPattern(Avv, i, range_begin, range_end, columns, &band);
		ierr = MatGetValues(Avv,1,&i,band,columns,rhs); CHKERRABORT(PETSC_COMM_WORLD, ierr);
		getSystem3(columns, band, least_squares_system);
//		luDecomposition2(least_squares_system, band);
//		bsfr2(least_squares_system, band, rhs);
		solve_lapack(least_squares_system, band, rhs);
		ierr = MatSetValues(sparse_approximate_inverse,band,columns,1,&i,rhs, INSERT_VALUES); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	}
	printf("SECOND PART %lf\n", tt2);
	printf("Third PART %lf\n", tt3);
	free(least_squares_system);
	free(columns);
	free(rhs);
//	free(row_dot_products);

	MatDestroy(&dot_product_matrix);
	MatAssemblyBegin(sparse_approximate_inverse,MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(sparse_approximate_inverse,MAT_FINAL_ASSEMBLY);
//	MatView(Avv, PETSC_VIEWER_DRAW_WORLD);
//	MatView(sparse_approximate_inverse, PETSC_VIEWER_DRAW_WORLD);
//	MatView(sparse_approximate_inverse, PETSC_VIEWER_STDOUT_WORLD);
//	draw_matrix(sparse_approximate_inverse, 10);

//	MatTranspose(Avv,MAT_REUSE_MATRIX,&Avv);
	return sparse_approximate_inverse;


}
