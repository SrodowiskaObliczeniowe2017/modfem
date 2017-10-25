/*
 * ResidualBasedErrorEvaluator.cpp
 *
 *  Created on: Aug 25, 2015
 *      Author: damian
 */

#include "ResidualBasedErrorEvaluator.hpp"

ResidualBasedErrorEvaluator::ResidualBasedErrorEvaluator(Mat mat, Vec rhs, double difference_threshold,
		double exact_solution_threshold, int size, int max_iterations) : ErrorEvaluator(){
	VecDuplicate(rhs, &residual);
	this->mat = mat;
	this->rhs = rhs;
	this->residual = residual;
	this->difference_threshold = difference_threshold;
	this->exact_solution_threshold = exact_solution_threshold;
	this->first_iteration = true;
	this->max_iterations = max_iterations;
}

ResidualBasedErrorEvaluator::~ResidualBasedErrorEvaluator() {
	VecDestroy(&residual);
}

bool ResidualBasedErrorEvaluator::Stop(Vec x, int iteration_nr)
{
	if(max_iterations != -1 && iteration_nr >= max_iterations)
		return true;
	if(max_iterations != -1)
		return false;
	MatResidual(mat,rhs,x,residual);
	PetscReal norm;
	VecNorm(residual,NORM_2,&norm);
//	printf("%20.15lf\n", norm);
//	PetscReal rhsNorm;
//	VecNorm(rhs,NORM_2,&rhsNorm);
	if(first_iteration)
	{
		first_iteration = false;
		previous_residual = norm;
		return norm < exact_solution_threshold;
	}
	bool stop = (fabs(previous_residual - norm) < difference_threshold) || (norm < exact_solution_threshold);
	previous_residual = norm;
	return stop;
}
