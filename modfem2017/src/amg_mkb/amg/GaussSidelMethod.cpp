/*
 * GaussSidelMethod.cpp
 *
 *  Created on: Jul 5, 2015
 *      Author: damian
 */

#include "GaussSidelMethod.hpp"

GaussSidelMethod::GaussSidelMethod() {

}

GaussSidelMethod::~GaussSidelMethod() {

}

void GaussSidelMethod::Smooth(Vec x)
{
	//MatView(matrix,PETSC_VIEWER_STDOUT_WORLD);
	//MatFactorInfoInitialize(MatFactorInfo *info)
	PetscErrorCode ierr = MatSOR(matrix,rhs,1.0,SOR_LOCAL_FORWARD_SWEEP,0.0,1,1,x);
	CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

void GaussSidelMethod::PreSmoothing(Mat matrix, Vec b)
{
	this->matrix = matrix;
	this->rhs = b;
	//MatView(matrix, PETSC_VIEWER_STDOUT_WORLD);
}

void GaussSidelMethod::SetRhs(Vec rhs)
{
	this->rhs = rhs;
}

Vec GaussSidelMethod::GetResidual(Vec x)
{
	Vec residual;
	Vec temporary_result;

	VecDuplicate(rhs,&residual);
	VecDuplicate(rhs,&temporary_result);
	PetscErrorCode ierr;

	ierr = MatMult(matrix, x, temporary_result); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = VecAXPY(rhs, -1.0, temporary_result); CHKERRABORT(PETSC_COMM_WORLD, ierr);

}
