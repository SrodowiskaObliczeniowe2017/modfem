#ifndef SMOOTHER_H
#define SMOOTHER_H

#include <petscvec.h>
#include <petscmat.h>

class Smoother
{
	public:
		virtual ~Smoother() {}
		virtual void Smooth(Vec x) = 0;
		virtual void PreSmoothing(Mat matrix, Vec b) = 0;
		virtual void SetRhs(Vec rhs) = 0;
		virtual Vec GetResidual(Vec x) = 0;

};


#endif
