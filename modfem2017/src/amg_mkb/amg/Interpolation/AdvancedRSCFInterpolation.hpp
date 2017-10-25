/*
 * AdvancedRSCFInterpolation.hpp
 *
 *  Created on: Jan 23, 2016
 *      Author: damian
 */

#ifndef SRC_AMG_MKB_AMG_INTERPOLATION_ADVANCEDRSCFINTERPOLATION_HPP_
#define SRC_AMG_MKB_AMG_INTERPOLATION_ADVANCEDRSCFINTERPOLATION_HPP_

#include "InterpolationStrategy.hpp"
#include <stdexcept>

#define MAX_ROW_SIZE 1200

class AdvancedRSCFInterpolation : public InterpolationStrategy {
public:
	AdvancedRSCFInterpolation(Mat mat, struct row_info* row_info_array, struct influenced_info* influenced_info_array,
			PetscInt first_row_in_range, PetscInt range_end, double strength_threshold);
	Mat GetMatrixFromCoarseToFine();
	virtual ~AdvancedRSCFInterpolation();

};

#endif /* SRC_AMG_MKB_AMG_INTERPOLATION_ADVANCEDRSCFINTERPOLATION_HPP_ */
