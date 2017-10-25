/*
 * GaussSidelMethod.h
 *
 *  Created on: Jul 5, 2015
 *      Author: damian
 */

#ifndef SRC_AMG_MKB_AMG_GAUSSSIDELMETHOD_HPP_
#define SRC_AMG_MKB_AMG_GAUSSSIDELMETHOD_HPP_

#include "Smoother.hpp"
#include <petscksp.h>

class GaussSidelMethod : public Smoother {

public:
	GaussSidelMethod();
	virtual void Smooth(Vec x);
	virtual void PreSmoothing(Mat matrix, Vec b);
	virtual ~GaussSidelMethod();
	virtual void SetRhs(Vec rhs);
	virtual Vec GetResidual(Vec x);
private:
	Mat matrix;
	Vec rhs;
	KSP ksp;
	PC pc;
};

#endif /* SRC_AMG_MKB_AMG_GAUSSSIDELMETHOD_HPP_ */
