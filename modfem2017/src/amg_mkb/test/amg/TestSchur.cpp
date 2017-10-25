/*
 * JacomiMethodTest.cpp
 *
 *  Created on: Jun 5, 2015
 *      Author: damian
 */

#define BOOST_TEST_MODULE SchurTest
#include <boost/test/included/unit_test.hpp>
#include <petscksp.h>

FILE* utv_log_out = stdout;

Mat Avv, Avp, Apv, App;
Mat schur_complement;
Mat preconditioner;
Mat identity;

void prepare_matrix(int i_start, int j_start, Mat* matrix)
{
	/*double values[] = { 3.0, 4.0, 1.0, 6.0,
						2.0, 3.0, 2.0, 4.0,
						3.0, 3.0, 2.0, 2.0,
						2.0, 1.0, 3.0, 4.0};*/

	double values[] = { 2.0, 0.0, 1.0, 6.0,
						0.0, 2.0, 2.0, 4.0,
						3.0, 3.0, 2.0, 2.0,
						2.0, 1.0, 3.0, 4.0};
	int size = 2;
	MatCreate(PETSC_COMM_WORLD, matrix);
	MatSetSizes(*matrix, PETSC_DECIDE, PETSC_DECIDE, size, size);
	MatSetFromOptions(*matrix);
	MatMPIAIJSetPreallocation(*matrix,30,NULL,30,NULL);
	MatSeqAIJSetPreallocation(*matrix,30,NULL);
	MatSetUp(*matrix);

	int i, j;
	for(i = 0; i<2; i++)
		for(j = 0; j<2; j++)
			MatSetValues(*matrix, 1, &i, 1, &j, &values[(i+i_start)*4 + j+j_start], INSERT_VALUES);

	MatAssemblyBegin(*matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*matrix, MAT_FINAL_ASSEMBLY);
	MatView(*matrix, PETSC_VIEWER_STDOUT_WORLD);
}

void prepare_identity(int size, Mat* matrix)
{
	MatCreate(PETSC_COMM_WORLD, matrix);
	MatSetSizes(*matrix, PETSC_DECIDE, PETSC_DECIDE, size, size);
	MatSetFromOptions(*matrix);
	MatMPIAIJSetPreallocation(*matrix,30,NULL,30,NULL);
	MatSeqAIJSetPreallocation(*matrix,30,NULL);
	MatSetUp(*matrix);

	int i;
	double value = 1;
	for(i = 0; i<size; i++)
		MatSetValues(*matrix, 1, &i, 1, &i, &value, INSERT_VALUES);

	MatAssemblyBegin(*matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*matrix, MAT_FINAL_ASSEMBLY);
	MatView(*matrix, PETSC_VIEWER_STDOUT_WORLD);
}

BOOST_AUTO_TEST_CASE( SchurTest )
{
	PetscInitialize(NULL,NULL,NULL,NULL);
	prepare_matrix(0,0,&Avv);
	prepare_matrix(0,2,&Avp);
	prepare_matrix(2,0,&Apv);
	prepare_matrix(2,2,&App);
	prepare_identity(2, &preconditioner);
//MatCreateSchurComplementPmat
	//MatCreateSchurComplement(Avv,preconditioner,Avp,Apv,App,&schur_complement);
	MatCreateSchurComplementPmat(Avv,Avp,Apv,App,MAT_SCHUR_COMPLEMENT_AINV_DIAG,MAT_INITIAL_MATRIX,&schur_complement);

	prepare_identity(2, &identity);
	Mat result;
	MatMatMult(identity,schur_complement,
					MAT_INITIAL_MATRIX,PETSC_DEFAULT,&result);


	MatView(schur_complement, PETSC_VIEWER_STDOUT_WORLD);
	printf("OK\n");
}
