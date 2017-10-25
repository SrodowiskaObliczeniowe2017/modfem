#ifndef ERROR_EVALUATOR_H
#define ERROR_EVALUATOR_H

#include <petscvec.h>
#include <petscmat.h>

class ErrorEvaluator
{
	public:
		ErrorEvaluator() {}
		virtual ~ErrorEvaluator() {}
		virtual bool Stop(Vec x, int iteration_nr) = 0;
	protected:
		Mat mat;
};


#endif
