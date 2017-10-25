#include "petscmat.h"

int rows[] = {0,0,0,0,0,0,0,1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,5,5,6,6,6,6,6,6,6};
int columns[] = {0,1,2,3,4,5,6,0,1,2,3,4,5,6,0,1,2,3,4,5,6,0,1,2,3,4,5,6,0,1,2,3,4,5,6,0,1,2,3,4,5,6,0,1,2,3,4,5,6};

int size = 7;

void prepare_simple_matrix(Mat* matrix, double* values)
{

	MatCreate(PETSC_COMM_WORLD, matrix);
	MatSetSizes(*matrix, PETSC_DECIDE, PETSC_DECIDE, size, size);
	MatSetFromOptions(*matrix);
	MatMPIAIJSetPreallocation(*matrix,30,NULL,30,NULL);
	MatSeqAIJSetPreallocation(*matrix,30,NULL);
	MatSetUp(*matrix);

	int i;
	for(i = 0; i<49; i++)
		MatSetValues(*matrix, 1, &rows[i], 1, &columns[i], &values[i], INSERT_VALUES);

	MatAssemblyBegin(*matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*matrix, MAT_FINAL_ASSEMBLY);
	MatView(*matrix, PETSC_VIEWER_STDOUT_WORLD);
}

void prepare_simple_rhs(Vec* vec, double* values)
{
	VecCreate(PETSC_COMM_WORLD,vec);
	VecSetSizes(*vec,PETSC_DECIDE,size);
	VecSetFromOptions(*vec);
	double zero = 0.0;
	double one = 1.0;
	int i;
	for (i=0; i<size -1; i++) {
	  VecSetValues(*vec,1,&i,&zero,INSERT_VALUES);
	}
	VecSetValues(*vec,1,&i,&one,INSERT_VALUES);
	VecView(*vec, PETSC_VIEWER_STDOUT_WORLD);
}

void prepare_testcase_1_matrix(Mat* matrix, int size)
{
	MatCreate(PETSC_COMM_WORLD, matrix);
	MatSetSizes(*matrix, PETSC_DECIDE, PETSC_DECIDE, size, size);
	MatSetFromOptions(*matrix);
	MatMPIAIJSetPreallocation(*matrix,30,NULL,30,NULL);
	MatSeqAIJSetPreallocation(*matrix,30,NULL);
	MatSetUp(*matrix);

	int row;
	int column;
	double value;

	for(int i = 1; i<size - 1; i++)
	{
		row = i;
		column = i - 1;
		value = -1;
		MatSetValues(*matrix, 1, &row, 1, &column, &value, INSERT_VALUES);
		column++;
		value = 2;
		MatSetValues(*matrix, 1, &row, 1, &column, &value, INSERT_VALUES);
		column++;
		value = -1;
		MatSetValues(*matrix, 1, &row, 1, &column, &value, INSERT_VALUES);
	}
	row = 0;
	column = 0;
	value = 1.0;
	MatSetValues(*matrix, 1, &row, 1, &column, &value, INSERT_VALUES);
	row = size - 1;
	column = size - 1;
	MatSetValues(*matrix, 1, &row, 1, &column, &value, INSERT_VALUES);

	MatAssemblyBegin(*matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*matrix, MAT_FINAL_ASSEMBLY);
	//MatView(*matrix, PETSC_VIEWER_STDOUT_WORLD);
}

void prepare_testcase_1_rhs(Vec* vec, int size)
{
	VecCreate(PETSC_COMM_WORLD,vec);
	VecSetSizes(*vec,PETSC_DECIDE,size);
	VecSetFromOptions(*vec);
	double zero = 0.0;
	double one = 1.0;
	int i;
	for (i=0; i<size -1; i++) {
	  VecSetValues(*vec,1,&i,&zero,INSERT_VALUES);
	}
	VecSetValues(*vec,1,&i,&one,INSERT_VALUES);
	//VecView(*vec, PETSC_VIEWER_STDOUT_WORLD);
}

void prepare_testcase_1_solution(Vec* vec, int size)
{
	VecCreate(PETSC_COMM_WORLD,vec);
	VecSetSizes(*vec,PETSC_DECIDE,size);
	VecSetFromOptions(*vec);
	double value = 0.0;
	double increase = 1.0/(size-1);
	int i;
	for (i=0; i<size; i++) {
	  VecSetValues(*vec,1,&i,&value,INSERT_VALUES);
	  value += increase;
	}
	//VecView(*vec, PETSC_VIEWER_STDOUT_WORLD);
}

void prepare_matrix_general_crs(Mat* matrix, double* values, int* rows, int* columns, int matrix_height, int matrix_width, int entries_size)
{
	MatCreate(PETSC_COMM_WORLD, matrix);
	MatSetSizes(*matrix, PETSC_DECIDE, PETSC_DECIDE, matrix_height, matrix_width);
	MatSetFromOptions(*matrix);
	MatMPIAIJSetPreallocation(*matrix,30,NULL,30,NULL);
	MatSeqAIJSetPreallocation(*matrix,30,NULL);
	MatSetUp(*matrix);

	int i;
	for(i = 0; i<entries_size; i++)
		MatSetValues(*matrix, 1, &rows[i], 1, &columns[i], &values[i], INSERT_VALUES);

	MatAssemblyBegin(*matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*matrix, MAT_FINAL_ASSEMBLY);
	MatView(*matrix, PETSC_VIEWER_STDOUT_WORLD);
}

void prepare_matrix_general_crs(Mat* matrix, double* values, int* rows, int* columns, int matrix_size, int entries_size)
{
	prepare_matrix_general_crs(matrix, values, rows, columns, matrix_size, matrix_size, entries_size);
}

void prepare_matrix_inverse_example_1(Mat* matrix){
	double values[] = {-4, 1, 1, 1, -4, 1, 1, 1, 1, -4, 1, 1, 1, 1, -4, 1, 1, 1, 1, -4, 1, 1, 1, -4};
	int rows[] = {0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5};
	int cols[] = {0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4, 1, 2, 3, 4, 5, 2, 3, 4, 5, 3, 4, 5};
	prepare_matrix_general_crs(matrix, values, rows, cols, 6, 24);
//	-4 1 1 0 0 0
//	1 -4 1 1 0 0
//	1 1 -4 1 1 0
//	0 1 1 -4 1 1
//	0 0 1 1 -4 1
//	0 0 0 1 1 -4
}
void prepare_matrix_inverse_example_2(Mat* matrix){
//	double values[] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
//			0, 0, 0, 0, 0, 1};
	double values[] = {5, 3, -1, 2, 4, 6, 3, -4, 2, -5, -2, 1, 1, 4, 3, -2, -2, 1, 8, -3, 1, 7, 3, 3, -3, 5, 3, -1, 4, 8,
			5, 7, 3, 3, 8, 1};
	int rows[] = {0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5};
	int cols[] = {0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
	prepare_matrix_general_crs(matrix, values, rows, cols, 6, 36);
	//5 3 -1 2 4 6
	//3 -4 2 -5 -2 1
	//1 4 3 -2 -2 1
	//8 -3 1 7 3 3
	//-3 5 3 -1 4 8
	//5 7 3 3 8 1
}

void prepare_matrix_inverse_example_3(Mat* matrix){
//	double values[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	double values[] = {1, 2, 3, 4, 6, 7, 9, 1, 5, 3, 8, 7};
	int rows[] = {0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5};
	int cols[] = {0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1};
	prepare_matrix_general_crs(matrix, values, rows, cols, 6, 2, 12);
//	1 2
//	3 4
//	6 7
//	9 1
//	5 3
//	8 7

//	1 3 6 9 5 8
//	2 4 7 1 3 7
}

void prepare_matrix_inverse_example_4(Mat* matrix){
//	double values[] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
	double values[] = {5, 4, 2, -2, 5, 9, 1, 3,-5, 3, 2, 5,-1, 2, 5, 8};
	int rows[] = {0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
	int cols[] = {0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
	prepare_matrix_general_crs(matrix, values, rows, cols, 4, 4, 16);
//	5 4 2 -2
//	5 9 1 3
//	-5 3 2 5
//	-1 2 5 8
}

void prepare_matrix_inverse_example_5(Mat* matrix){
	double values[] = {2, 3, 1, 4, 5, 8, 7, 3, 2, 1, 9, 0, 1, 3, 2, 8, 7, 2, 1, 5, 4, 2, 1, 4, 5};
	int rows[] = {0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4};
	int cols[] = {0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4, 0, 1, 2, 3, 4};
	prepare_matrix_general_crs(matrix, values, rows, cols, 5, 5, 25);
	//2 3 1 4 5
	//8 7 3 2 1
	//9 0 1 3 2
	//8 7 2 1 5
	//4 2 1 4 5
}

void prepare_matrix_inverse_full_pattern(Mat* matrix){
	double values[] = {0, 0, 0, 0, 0, 0,
					  1, 1, 1, 1, 1, 1,
					  2, 2, 2, 2, 2, 2,
					  3, 3, 3, 3, 3, 3,
					  4, 4, 4, 4, 4, 4,
					  5, 5, 5, 5, 5, 5};
	int rows[] = {0, 0, 0, 0, 0, 0,
				  1, 1, 1, 1, 1, 1,
				  2, 2, 2, 2, 2, 2,
				  3, 3, 3, 3, 3, 3,
				  4, 4, 4, 4, 4, 4,
				  5, 5, 5, 5, 5, 5};
	int cols[] = {0, 1, 2, 3, 4, 5,
				  0, 1, 2, 3, 4, 5,
				  0, 1, 2, 3, 4, 5,
				  0, 1, 2, 3, 4, 5,
				  0, 1, 2, 3, 4, 5,
				  0, 1, 2, 3, 4, 5};

	prepare_matrix_general_crs(matrix, values, rows, cols, 6, 36);
}
