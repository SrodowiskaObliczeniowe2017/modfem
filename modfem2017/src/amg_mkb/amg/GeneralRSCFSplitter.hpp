/*
 * GeneralRSCFSplitter.hpp
 *
 *  Created on: Feb 3, 2017
 *      Author: damian
 */

#ifndef SRC_AMG_MKB_AMG_GENERALRSCFSPLITTER_HPP_
#define SRC_AMG_MKB_AMG_GENERALRSCFSPLITTER_HPP_

#include "RSCFSplitter.hpp"

class GeneralRSCFSplitter : public RSCFSplitter {
public:
	GeneralRSCFSplitter(Mat mat, double strength_threshold, InterpolationStrategy* strategy);
	virtual ~GeneralRSCFSplitter();

protected:
	virtual PetscScalar GetExtremalValue(const PetscScalar* values, PetscInt size, PetscInt row_index,
			const PetscInt* columns, int first_row_in_range, int range_end);
	virtual bool IsStrongDependenceWithinRange(PetscScalar value, PetscScalar row_min, PetscInt row_index, PetscInt column_index,
			PetscInt first_row_in_range, PetscInt range_end, double strength_threshold);
};

#endif /* SRC_AMG_MKB_AMG_GENERALRSCFSPLITTER_HPP_ */
