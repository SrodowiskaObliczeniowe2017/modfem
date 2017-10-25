/*
 * DirectSolver.cpp
 *
 *  Created on: Jul 27, 2015
 *      Author: damian
 */

#include "DirectSolver.hpp"

DirectSolver::DirectSolver() {


}

DirectSolver::~DirectSolver() {
	KSPDestroy(&ksp);
}

Vec DirectSolver::GetResidual(Vec x)
{
	throw new std::runtime_error("Not implemented");
}

void DirectSolver::Smooth(Vec x)
{
	//KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
	KSPSolve(ksp,rhs,x);


	//PetscInt iterations_number;
	//KSPGetIterationNumber(ksp,&iterations_number);
}

void DirectSolver::PreSmoothing(Mat matrix, Vec b)
{
	this->matrix = matrix;
	this->rhs = b;
	PetscErrorCode ierr;

	KSPCreate(PETSC_COMM_WORLD,&ksp);
	KSPSetOperators(ksp,matrix,matrix);

	//PC pc;
	//KSPGetPC(ksp,&pc);
	//PCSetType(pc,PCNONE);

	ierr = KSPSetFromOptions(ksp);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	KSPSetUp(ksp);CHKERRABORT(PETSC_COMM_WORLD, ierr);
}

void DirectSolver::SetRhs(Vec rhs)
{
	this->rhs = rhs;
}
