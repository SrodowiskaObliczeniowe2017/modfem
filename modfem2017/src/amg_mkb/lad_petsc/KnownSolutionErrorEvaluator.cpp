/*
 * KnownSolutionErrorEvaluator.cpp
 *
 *  Created on: Aug 10, 2015
 *      Author: damian
 */

#include "KnownSolutionErrorEvaluator.hpp"

KnownSolutionErrorEvaluator::KnownSolutionErrorEvaluator(int size) : ErrorEvaluator() {
	this->size = size;
	VecCreate(PETSC_COMM_WORLD,&exact_solution);
	VecSetSizes(exact_solution,PETSC_DECIDE,size);
	VecSetFromOptions(exact_solution);
	double value = 0.0;
	double increase = 1.0/(size-1);
	int i;
	for (i=0; i<size; i++) {
	  VecSetValues(exact_solution,1,&i,&value,INSERT_VALUES);
	  value += increase;
	}

	VecDuplicate(exact_solution,&solution_difference);
	VecAssemblyBegin(exact_solution);
	VecAssemblyEnd(exact_solution);
	VecAssemblyBegin(solution_difference);
	VecAssemblyEnd(solution_difference);
}

KnownSolutionErrorEvaluator::~KnownSolutionErrorEvaluator() {

}

bool KnownSolutionErrorEvaluator::Stop(Vec x, int iteration_nr)
{
	double l2_norm;
	VecCopy(exact_solution,solution_difference);
	VecAXPY(solution_difference, -1.0, x);
	VecNorm(solution_difference,NORM_2,&l2_norm);

	printf("L2 norm: %lf\n", l2_norm);

	return l2_norm < 0.0001;
}

