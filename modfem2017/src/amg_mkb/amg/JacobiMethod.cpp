/*
 * JacobiMethod.cpp
 *
 *  Created on: Jun 5, 2015
 *      Author: damian
 */

#include "JacobiMethod.hpp"

JacobiMethod::JacobiMethod()
{
	// TODO Auto-generated constructor stub

}

JacobiMethod::~JacobiMethod()
{

	PetscErrorCode ierr;

	ierr = VecDestroy(&diagonal); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecDestroy(&bdiagonal); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecDestroy(&temporary_result); CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

void JacobiMethod::Smooth(Vec x)
{
	PetscErrorCode ierr;

	ierr = MatMult(scaledMatrix, x, temporary_result); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecAXPY(x, -1.0, temporary_result); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecAXPY(x, 1.0, bdiagonal); CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

void JacobiMethod::SetRhs(Vec rhs)
{
	throw new std::runtime_error("Not implemented");
}

Vec JacobiMethod::GetResidual(Vec x)
{
	throw new std::runtime_error("Not implemented");
}

void JacobiMethod::PreSmoothing(Mat matrix, Vec b)
{
	PetscErrorCode ierr;

	//same layout of diagonal as b and x vectors
	ierr = VecDuplicate(b,&diagonal); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatGetDiagonal(matrix, diagonal); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	//inversed diagonal
	ierr = VecReciprocal(diagonal); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatDiagonalScale(matrix,diagonal,NULL); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	scaledMatrix = matrix;


	ierr = VecDuplicate(b,&bdiagonal); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecPointwiseMult(bdiagonal, diagonal, b); CHKERRABORT(PETSC_COMM_WORLD, ierr);

	ierr = VecDuplicate(b,&temporary_result); CHKERRABORT(PETSC_COMM_WORLD, ierr);
}
