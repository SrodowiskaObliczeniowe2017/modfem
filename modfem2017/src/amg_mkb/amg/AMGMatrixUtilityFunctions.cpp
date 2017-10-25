#include "AMGMatrixUtilityFunctions.hpp"

//double value_filter_threshold = 1.0e-12;
double value_filter_threshold = 0;
//double value_filter_threshold = 0;
int max_unfiltered_row_size = 4000;

bool isInCSet(PetscInt row_index, RowData* rowData, PetscInt first_row_in_range, PetscInt range_end,
		struct influenced_info* influenced_info_array)
{
	if(row_index - first_row_in_range >= 0 && row_index < range_end)
		return influenced_info_array[row_index - first_row_in_range].row_info->set == CSET;
	return rowData->GetRowData(row_index) != -1;

}

bool isInFSetOutOfRange(PetscInt row_index, RowData* rowData, PetscInt first_row_in_range, PetscInt range_end,
		struct influenced_info* influenced_info_array)
{
	if(row_index < first_row_in_range || row_index >= range_end)
		return rowData->GetRowData(row_index) == -1;
	return false;

}

PetscInt getColumnIndex(PetscInt row_index, RowData* rowData, PetscInt first_row_in_range, PetscInt range_end,
		struct influenced_info* influenced_info_array)
{
	if(row_index - first_row_in_range >= 0 && row_index < range_end)
		return influenced_info_array[row_index - first_row_in_range].row_number_in_coarse;
	return rowData->GetRowData(row_index);
}

void increase_nonzero_number(int column_number_in_coarse, int row_number,
		int from_coarse_to_fine_ownership_column_begin, int from_coarse_to_fine_ownership_column_end,
		int* nondiagonal_nonzero_numbers, int* diagonal_nonzero_numbers)
{
	if(column_number_in_coarse < from_coarse_to_fine_ownership_column_begin ||
			column_number_in_coarse >= from_coarse_to_fine_ownership_column_end)
	{
		nondiagonal_nonzero_numbers[row_number]++;
	}
	else
	{
		diagonal_nonzero_numbers[row_number]++;
	}
}

bool isStrongDependenceWithinRange(PetscScalar value, PetscScalar row_min, PetscInt row_index, PetscInt column_index,
		PetscInt first_row_in_range, PetscInt range_end, double strength_threshold)
{
	//printf("strong? %lf %d %d\n", value, row_index, column_index);
//	if(value > 0)
//		value = -value;
	return isStrongDependence(value, row_min, row_index, column_index, strength_threshold) &&
			(column_index >= first_row_in_range && column_index < range_end);
}

bool isStrongDependence(PetscScalar value, PetscScalar row_min, PetscInt row_index, PetscInt column_index, double strength_threshold)
{
	return (-value > strength_threshold*-row_min) && (column_index != row_index);
}

//TODO: consider min evaluation boundaries
PetscScalar getMin(const PetscScalar* values, PetscInt size, PetscInt row_index,const PetscInt* columns, int first_row_in_range, int range_end)
{
	PetscScalar min_value = values[0];
//	PetscScalar max_value = values[0];
	if(columns[0] == row_index && size > 1 /*&& columns[1] >= first_row_in_range && columns[1] < range_end*/){
		min_value = values[1];
//		max_value = values[1];
	}
	for(PetscInt i = 0; i<size; i++)
	{
		if(columns[i] != row_index /*&& columns[i] >= first_row_in_range && columns[i] < range_end*/)
		{
			min_value = std::min(min_value, values[i]);
//			max_value = std::max(max_value, values[i]);
		}
	}
//	min_value = std::min(min_value, -max_value);
	return min_value;
}

bool isGeneralStrongDependenceWithinRange(PetscScalar value, PetscScalar row_min, PetscInt row_index, PetscInt column_index,
		PetscInt first_row_in_range, PetscInt range_end, double strength_threshold)
{
	//printf("strong? %lf %d %d\n", value, row_index, column_index);
//	if(value > 0)
//		value = -value;
	return isGeneralStrongDependence(value, row_min, row_index, column_index, strength_threshold) &&
			(column_index >= first_row_in_range && column_index < range_end);
}

bool isGeneralStrongDependence(PetscScalar value, PetscScalar row_min, PetscInt row_index, PetscInt column_index, double strength_threshold)
{
	return (fabs(value) > strength_threshold*fabs(row_min)) && (column_index != row_index);
}

//TODO: consider min evaluation boundaries
PetscScalar getExtremalValue(const PetscScalar* values, PetscInt size, PetscInt row_index,const PetscInt* columns, int first_row_in_range, int range_end)
{
	PetscScalar extremal_value = values[0];
//	PetscScalar max_value = values[0];
	if(columns[0] == row_index && size > 1 /*&& columns[1] >= first_row_in_range && columns[1] < range_end*/){
		extremal_value = values[1];
//		max_value = values[1];
	}
	for(PetscInt i = 0; i<size; i++)
	{
		if(columns[i] != row_index /*&& columns[i] >= first_row_in_range && columns[i] < range_end*/)
		{
			extremal_value = std::max(extremal_value, fabs(values[i]));
//			max_value = std::max(max_value, values[i]);
		}
	}
//	min_value = std::min(min_value, -max_value);
	return extremal_value;
}

void getFilteredRow(Mat mat, int row_number, filtered_row *row){
	MatGetRow(mat,row_number,&(row->unfiltered_columns_number),
			(const PetscInt**)&(row->unfiltered_columns),
			(const PetscScalar**)&(row->unfiltered_values));
	if(row->unfiltered_columns_number > max_unfiltered_row_size)
		throw new std::runtime_error("Not implemented");
	//row->columns_number = row->unfiltered_columns_number;
//	bool pin = false;
	int columns_number = 0;
	for(int i = 0; i<row->unfiltered_columns_number; i++){
		if(row->unfiltered_values[i] > value_filter_threshold || row->unfiltered_values[i] < -value_filter_threshold){
			row->columns[columns_number] = row->unfiltered_columns[i];
			row->values[columns_number++] = row->unfiltered_values[i];
		}

//		if(row->unfiltered_columns[i] == row_number){
//			if(fabs(row->unfiltered_values[i]) > 1e+7){
//				pin = true;
//			}
//		}
	}
//	if(pin){
//		for(int i = 0; i<columns_number; i++)
//			row->values[i] /= 1000000;
//	}

	row->columns_number = columns_number;
}

void restoreFilteredRow(Mat mat, int row_number, filtered_row *row){
	MatRestoreRow(mat, row_number, &(row->unfiltered_columns_number),
			(const PetscInt**)&(row->unfiltered_columns),
			(const PetscScalar**)&(row->unfiltered_values));
}


