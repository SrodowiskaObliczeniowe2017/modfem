/*
 * JacobiMethod.hpp
 *
 *  Created on: Jun 5, 2015
 *      Author: damian
 */

#ifndef JACOBIMETHOD_HPP_
#define JACOBIMETHOD_HPP_

#include "Smoother.hpp"
#include <exception>
#include <stdexcept>

class JacobiMethod : public Smoother {
	public:
		JacobiMethod();
		virtual ~JacobiMethod();
		virtual void Smooth(Vec x);
		virtual void PreSmoothing(Mat matrix, Vec b);
		virtual void SetRhs(Vec rhs);
		virtual Vec GetResidual(Vec x);

	private:
		Mat scaledMatrix;
		Vec diagonal;
		Vec bdiagonal;
		Vec temporary_result;
};

#endif /* SRC_AMG_MKB_AMG_JACOBIMETHOD_HPP_ */
