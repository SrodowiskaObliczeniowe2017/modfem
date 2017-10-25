/*
 * DirectInterpolation.hpp
 *
 *  Created on: Jan 22, 2016
 *      Author: damian
 */

#ifndef SRC_AMG_MKB_AMG_INTERPOLATION_DIRECTINTERPOLATION_HPP_
#define SRC_AMG_MKB_AMG_INTERPOLATION_DIRECTINTERPOLATION_HPP_
#include "InterpolationStrategy.hpp"
#include "../AMGMatrixUtilityFunctions.hpp"
class DirectInterpolation : public InterpolationStrategy{
public:
	DirectInterpolation(Mat mat, struct row_info* row_info_array, struct influenced_info* influenced_info_array,
			PetscInt first_row_in_range, PetscInt range_end, double strength_threshold);
	DirectInterpolation();
	//void PrepareFromCoarseToFineMat(Mat* from_coarse_to_fine_mat, int coarse_rows_number, int coarse_columns_number);
	virtual ~DirectInterpolation();
	Mat GetMatrixFromCoarseToFine();
private:
	int SetCoarseRowsNumber();
	void Handle3RS(RowData* rowData);
};

#endif /* SRC_AMG_MKB_AMG_INTERPOLATION_DIRECTINTERPOLATION_HPP_ */
