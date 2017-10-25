/*
 * ResidualBasedErrorEvaluator.hpp
 *
 *  Created on: Aug 25, 2015
 *      Author: damian
 */

#ifndef SRC_AMG_MKB_LAD_PETSC_RESIDUALBASEDERROREVALUATOR_HPP_
#define SRC_AMG_MKB_LAD_PETSC_RESIDUALBASEDERROREVALUATOR_HPP_

#include <petscmat.h>
#include <petscvec.h>
#include "ErrorEvaluator.hpp"

class ResidualBasedErrorEvaluator : public ErrorEvaluator{
public:
	ResidualBasedErrorEvaluator(Mat mat, Vec rhs, double difference_threshold, double exact_solution_threshold, int size, int max_iterations);
	virtual ~ResidualBasedErrorEvaluator();
	virtual bool Stop(Vec x, int iteration_nr);
private:
	Vec residual;
	Mat mat;
	Vec rhs;
	double difference_threshold;
	double previous_residual;
	double exact_solution_threshold;
	bool first_iteration;
	int max_iterations;
};

#endif /* SRC_AMG_MKB_LAD_PETSC_RESIDUALBASEDERROREVALUATOR_HPP_ */
