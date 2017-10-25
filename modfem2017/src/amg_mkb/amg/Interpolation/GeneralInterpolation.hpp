/*
 * GeneralInterpolation.hpp
 *
 *  Created on: Feb 3, 2017
 *      Author: damian
 */

#ifndef SRC_AMG_MKB_AMG_INTERPOLATION_GENERALINTERPOLATION_HPP_
#define SRC_AMG_MKB_AMG_INTERPOLATION_GENERALINTERPOLATION_HPP_

#include "InterpolationStrategy.hpp"

class GeneralInterpolation : public InterpolationStrategy{
public:
	GeneralInterpolation();
	virtual ~GeneralInterpolation();
	Mat GetMatrixFromCoarseToFine();
};

#endif /* SRC_AMG_MKB_AMG_INTERPOLATION_GENERALINTERPOLATION_HPP_ */
