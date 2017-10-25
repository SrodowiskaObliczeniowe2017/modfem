/*
 * KnownSolutionErrorEvaluator.hpp
 *
 *  Created on: Aug 10, 2015
 *      Author: damian
 */

#ifndef SRC_AMG_MKB_LAD_PETSC_KNOWNSOLUTIONERROREVALUATOR_HPP_
#define SRC_AMG_MKB_LAD_PETSC_KNOWNSOLUTIONERROREVALUATOR_HPP_
#include "ErrorEvaluator.hpp"

class KnownSolutionErrorEvaluator : public ErrorEvaluator{
public:
	KnownSolutionErrorEvaluator(int size);
	virtual ~KnownSolutionErrorEvaluator();
	virtual bool Stop(Vec x, int iteration_nr);
private:
	int size;
	Vec exact_solution;
	Vec solution_difference;
};

#endif /* SRC_AMG_MKB_LAD_PETSC_KNOWNSOLUTIONERROREVALUATOR_HPP_ */
