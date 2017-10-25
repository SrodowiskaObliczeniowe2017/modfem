/*
 * InterpolationStrategy.hpp
 *
 *  Created on: Jan 22, 2016
 *      Author: damian
 */

#ifndef SRC_AMG_MKB_AMG_INTERPOLATION_INTERPOLATIONSTRATEGY_HPP_
#define SRC_AMG_MKB_AMG_INTERPOLATION_INTERPOLATIONSTRATEGY_HPP_
#include <petscksp.h>
#include "../AMGRows.hpp"
#include "../RowData.hpp"
#include "../RowDataCollector.hpp"
#include "../AMGMatrixUtilityFunctions.hpp"

class InterpolationStrategy {
public:
	InterpolationStrategy(Mat mat, struct row_info* row_info_array, struct influenced_info* influenced_info_array,
			PetscInt first_row_in_range, PetscInt range_end, double strength_threshold);
	InterpolationStrategy();
	void InitStrategy(Mat mat, struct row_info* row_info_array, struct influenced_info* influenced_info_array,
			PetscInt first_row_in_range, PetscInt range_end, double strength_threshold);
	virtual ~InterpolationStrategy();

	virtual Mat GetMatrixFromCoarseToFine() = 0;


protected:
	void PrepareFromCoarseToFineMat(Mat* from_coarse_to_fine_mat, int coarse_rows_number, int coarse_columns_number, RowData* rowData);
	bool IsInCSet(PetscInt row_index, RowData* rowData, PetscInt first_row_in_range, PetscInt range_end,
			struct influenced_info* influenced_info_array);
	PetscInt GetColumnIndex(PetscInt row_index, RowData* rowData, PetscInt first_row_in_range, PetscInt range_end,
			struct influenced_info* influenced_info_array);
	Mat mat;
	struct row_info* row_info_array;
	struct influenced_info* influenced_info_array;
	PetscInt first_row_in_range;
	PetscInt range_end;
	double strength_threshold;
};

#endif /* SRC_AMG_MKB_AMG_INTERPOLATION_INTERPOLATIONSTRATEGY_HPP_ */
