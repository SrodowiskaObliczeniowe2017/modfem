/*
 * JacomiMethodTest.cpp
 *
 *  Created on: Jun 5, 2015
 *      Author: damian
 */

#define BOOST_TEST_MODULE JacobiMethodTest
#include <boost/test/included/unit_test.hpp>
#include "../../amg/JacobiMethod.hpp"
#include "../../amg/GaussSidelMethod.hpp"
#include "petscmat.h"
#include "petscvec.h"

FILE* utv_log_out = stdout;


void prepare_matrix(Mat* matrix, double* values)
{
	int size = 3;
	MatCreate(PETSC_COMM_WORLD, matrix);
	MatSetSizes(*matrix, PETSC_DECIDE, PETSC_DECIDE, size, size);
	MatSetFromOptions(*matrix);
	MatMPIAIJSetPreallocation(*matrix,30,NULL,30,NULL);
	MatSeqAIJSetPreallocation(*matrix,30,NULL);
	MatSetUp(*matrix);

	int row_numbers[] = {0,0,0,1,1,1,2,2,2};
	int column_numbers[] = {0,1,2,0,1,2,0,1,2};

	int i;
	for(i = 0; i<9; i++)
		MatSetValues(*matrix, 1, &row_numbers[i], 1, &column_numbers[i], &values[i], INSERT_VALUES);

	MatAssemblyBegin(*matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*matrix, MAT_FINAL_ASSEMBLY);
	MatView(*matrix, PETSC_VIEWER_STDOUT_WORLD);
}

void prepare_rhs(Vec* rhs, double* values)
{
	VecCreate(PETSC_COMM_WORLD, rhs);
	VecSetSizes(*rhs, PETSC_DECIDE, 3);
	VecSetFromOptions(*rhs);
	int row_numbers[] = {0,1,2};
	int i;

	for(i = 0; i<3; i++)
		VecSetValues(*rhs, 1, &row_numbers[i], &values[i], INSERT_VALUES);
	VecView(*rhs, PETSC_VIEWER_STDOUT_WORLD);
}

void prepare_system(Mat* matrix, Vec* rhs, Vec* x, double* matrix_values, double* rhs_values)
{
	prepare_matrix(matrix, matrix_values);
	prepare_rhs(rhs, rhs_values);

	VecCreate(PETSC_COMM_WORLD, x);
	VecSetSizes(*x, PETSC_DECIDE, 3);
	VecSetFromOptions(*x);
	VecSet(*x, 0.0);
}

BOOST_AUTO_TEST_CASE( JacobiMethodTest )
{
	PetscInitialize(NULL,NULL,NULL,NULL);
	Mat matrix;
	Vec rhs,x;
	double matrix_values[] = {5.0,-2.0,3.0,-3.0,9.0,1.0,2.0,-1.0,-7.0};
	double rhs_values[] = {-1.0,2.0,3.0};
	prepare_system(&matrix, &rhs, &x, matrix_values, rhs_values);
	JacobiMethod* jacobiMethod = new JacobiMethod();
	jacobiMethod->PreSmoothing(matrix, rhs);
	int i;
	for(i = 0; i<10; i++)
	{
		jacobiMethod->Smooth(x);
	}
	printf("\n Vector after  10 jacobi method iterations \n");
		VecView(x, PETSC_VIEWER_STDOUT_WORLD);
	printf("\n");
	PetscScalar values[] = {0.186119, 0.331232, -0.422711};
	PetscScalar computed_solution[3];
	PetscInt indices[] = {0,1,2};
	VecGetValues(x,3,indices,computed_solution);

    BOOST_CHECK_CLOSE(computed_solution[0], values[0], 0.001);
    BOOST_CHECK_CLOSE(computed_solution[1], values[1], 0.001);
    BOOST_CHECK_CLOSE(computed_solution[2], values[2], 0.001);
}

BOOST_AUTO_TEST_CASE( GaussSidelMethodTest )
{
	PetscInitialize(NULL,NULL,NULL,NULL);
	Mat matrix;
	Vec rhs,x;
	double matrix_values[] = {4.0,1.0,-1.0,2.0,7.0,1.0,1.0,-3.0,12.0};
	double rhs_values[] = {3.0,19.0,31.0};
	prepare_system(&matrix, &rhs, &x, matrix_values, rhs_values);
	GaussSidelMethod* gsMethod = new GaussSidelMethod();
	gsMethod->PreSmoothing(matrix, rhs);
	int i;
	for(i = 0; i<1; i++)
	{
		gsMethod->Smooth(x);
	}
	printf("\n Vector after  10 gs method iterations \n");
		VecView(x, PETSC_VIEWER_STDOUT_WORLD);
	printf("\n");
	PetscScalar values[] = {1.0, 2.0, 3.0};
	PetscScalar computed_solution[3];
	PetscInt indices[] = {0,1,2};
	VecGetValues(x,3,indices,computed_solution);

    BOOST_CHECK_CLOSE(computed_solution[0], values[0], 0.001);
    BOOST_CHECK_CLOSE(computed_solution[1], values[1], 0.001);
    BOOST_CHECK_CLOSE(computed_solution[2], values[2], 0.001);
}
