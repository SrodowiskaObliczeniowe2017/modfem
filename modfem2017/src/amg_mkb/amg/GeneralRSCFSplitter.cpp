/*
 * GeneralRSCFSplitter.cpp
 *
 *  Created on: Feb 3, 2017
 *      Author: damian
 */

#include "GeneralRSCFSplitter.hpp"

GeneralRSCFSplitter::GeneralRSCFSplitter(Mat mat, double strength_threshold, InterpolationStrategy* strategy) :
	RSCFSplitter(mat,strength_threshold,strategy) {
	// TODO Auto-generated constructor stub

}

GeneralRSCFSplitter::~GeneralRSCFSplitter() {
	// TODO Auto-generated destructor stub
}


PetscScalar GeneralRSCFSplitter::GetExtremalValue(const PetscScalar* values, PetscInt size, PetscInt row_index,
		const PetscInt* columns, int first_row_in_range, int range_end)
{
	return getExtremalValue(values, size, row_index, columns, first_row_in_range, range_end);
}

bool GeneralRSCFSplitter::IsStrongDependenceWithinRange(PetscScalar value, PetscScalar row_min, PetscInt row_index, PetscInt column_index,
		PetscInt first_row_in_range, PetscInt range_end, double strength_threshold)
{
	return isGeneralStrongDependenceWithinRange(value, row_min, row_index, column_index, first_row_in_range, range_end, strength_threshold);
}
