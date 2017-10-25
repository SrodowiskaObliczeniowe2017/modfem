/*
 * DirectSolver.hpp
 *
 *  Created on: Jul 27, 2015
 *      Author: damian
 */

#ifndef SRC_AMG_MKB_AMG_DIRECTSOLVER_HPP_
#define SRC_AMG_MKB_AMG_DIRECTSOLVER_HPP_

#include "Smoother.hpp"
#include <stdexcept>
#include <petscksp.h>

class DirectSolver : public Smoother {
public:
	DirectSolver();
	virtual void Smooth(Vec x);
	virtual void PreSmoothing(Mat matrix, Vec b);
	virtual void SetRhs(Vec rhs);
	virtual Vec GetResidual(Vec x);
	virtual ~DirectSolver();
private:
	Mat matrix;
	Vec rhs;
	KSP ksp;
};

#endif /* SRC_AMG_MKB_AMG_DIRECTSOLVER_HPP_ */
