#ifndef CFSPLITTER_H
#define CFSPLITTER_H

#include <petscmat.h>
#include <list>

class CFSplitter
{
	public:
		virtual void MakeCFSplitting() = 0;
		virtual Mat GetMatrixFromCoarseToFine() = 0;
		virtual std::list<int>* GetNumbersOfCoarseRows() = 0;
		virtual ~CFSplitter();
	protected:
		Mat mat;
		CFSplitter(Mat mat);
};

#endif
