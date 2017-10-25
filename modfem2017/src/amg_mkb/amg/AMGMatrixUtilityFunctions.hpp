
#ifndef SRC_AMG_MKB_AMG_MATRIX_UTILITY_FUNCTIONS_HPP_
#define SRC_AMG_MKB_AMG_MATRIX_UTILITY_FUNCTIONS_HPP_
#include <petscksp.h>
#include "AMGRows.hpp"
#include "RowData.hpp"
#include "RowDataCollector.hpp"
#include "AMGMatrixUtilityFunctions.hpp"
#include <stdexcept>

extern int max_unfiltered_row_size;
extern double value_filter_threshold;


typedef struct {
	int columns_number;
	int columns[6000];
	double values[6000];
	PetscInt unfiltered_columns_number;
	PetscInt* unfiltered_columns;
	PetscScalar* unfiltered_values;
} filtered_row;

bool isInCSet(PetscInt row_index, RowData* rowData, PetscInt first_row_in_range, PetscInt range_end,
		struct influenced_info* influenced_info_array);

PetscInt getColumnIndex(PetscInt row_index, RowData* rowData, PetscInt first_row_in_range, PetscInt range_end,
		struct influenced_info* influenced_info_array);

void increase_nonzero_number(int column_number_in_coarse, int row_number,
		int from_coarse_to_fine_ownership_column_begin, int from_coarse_to_fine_ownership_column_end,
		int* nondiagonal_nonzero_numbers, int* diagonal_nonzero_numbers);

bool isStrongDependenceWithinRange(PetscScalar value, PetscScalar row_min, PetscInt row_index, PetscInt column_index,
		PetscInt first_row_in_range, PetscInt range_end, double strength_threshold);

bool isStrongDependence(PetscScalar value, PetscScalar row_min, PetscInt row_index, PetscInt column_index, double strength_threshold);

bool isInFSetOutOfRange(PetscInt row_index, RowData* rowData, PetscInt first_row_in_range, PetscInt range_end,
		struct influenced_info* influenced_info_array);

PetscScalar getMin(const PetscScalar* values, PetscInt size, PetscInt row_index,const PetscInt* columns, int first_row_in_range, int range_end);

bool isGeneralStrongDependenceWithinRange(PetscScalar value, PetscScalar row_min, PetscInt row_index, PetscInt column_index,
		PetscInt first_row_in_range, PetscInt range_end, double strength_threshold);

bool isGeneralStrongDependence(PetscScalar value, PetscScalar row_min, PetscInt row_index, PetscInt column_index, double strength_threshold);

PetscScalar getExtremalValue(const PetscScalar* values, PetscInt size, PetscInt row_index,const PetscInt* columns, int first_row_in_range, int range_end);

void getFilteredRow(Mat mat, int row_number, filtered_row *row);
void restoreFilteredRow(Mat mat, int row_number, filtered_row *row);

#endif
