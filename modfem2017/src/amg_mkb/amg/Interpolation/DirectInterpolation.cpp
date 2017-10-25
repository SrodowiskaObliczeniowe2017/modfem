/*
 * DirectInterpolation.cpp
 *
 *  Created on: Jan 22, 2016
 *      Author: damian
 */

#include "DirectInterpolation.hpp"

DirectInterpolation::DirectInterpolation(Mat mat, struct row_info* row_info_array, struct influenced_info* influenced_info_array,
		PetscInt first_row_in_range, PetscInt range_end, double strength_threshold) : InterpolationStrategy(mat, row_info_array, influenced_info_array,
				first_row_in_range, range_end, strength_threshold)
{

}

DirectInterpolation::DirectInterpolation() {

}

DirectInterpolation::~DirectInterpolation() {
	// TODO Auto-generated destructor stub
}


PetscScalar getSum(const PetscScalar* values, PetscInt size, PetscInt row_index,const PetscInt* columns)
{
	PetscScalar sum = 0;
	for(PetscInt i = 0; i<size; i++)
	{
		if(columns[i] != row_index)
		{
			sum += values[i];
		}
	}
	return sum;
}

int DirectInterpolation::SetCoarseRowsNumber()
{

	PetscInt ownership_range_size = range_end - first_row_in_range;
	PetscInt row_number_in_coarse = 0;
	for(PetscInt i = 0; i<ownership_range_size; i++)
	{
		if(influenced_info_array[i].row_info->set == CSET)
		{
			influenced_info_array[i].row_number_in_coarse = row_number_in_coarse++;
		}
	}

	int rank = 0;
	int size = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	int coarse_columns_numbers[size];

	MPI_Allgather(&row_number_in_coarse, 1, MPI_INT, coarse_columns_numbers, 1, MPI_INT, PETSC_COMM_WORLD);

	for(int i = 0; i<size-1; i++)
	{
		coarse_columns_numbers[i+1] += coarse_columns_numbers[i];
	}
	int coarse_column_number_offset = 0;
	if(rank > 0)
		coarse_column_number_offset = coarse_columns_numbers[rank - 1];

	for(PetscInt i = 0; i<ownership_range_size; i++)
	{
		if(influenced_info_array[i].row_info->set == CSET)
		{
			influenced_info_array[i].row_number_in_coarse += coarse_column_number_offset;
		}
	}
//	printf("coarse rows %d\n",row_number_in_coarse);

	return coarse_columns_numbers[size-1];
}

void DirectInterpolation::Handle3RS(RowData* rowData){

	int rank = 0;
	int size = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &size);

	int potential_strong_ff[size];
	for(int i = 0; i<size; i++)
		potential_strong_ff[i] = 0;

	PetscInt ownership_range_size = range_end - first_row_in_range;
	PetscInt row_index;
	filtered_row row;

	for(PetscInt i = 0; i<ownership_range_size; i++)
	{
		row_index = i + first_row_in_range;
		if(influenced_info_array[i].row_info->set == FSET)
		{
			getFilteredRow(mat,row_index,&row);
			double min = getMin(row.values, row.columns_number, row_index, row.columns, first_row_in_range, range_end);
			PetscInt strongInterpolating[row.columns_number];
			PetscInt strongNonInterpolatingOutOfRange[row.columns_number];
			PetscInt strongInterpolatingSize = 0;
			PetscInt strongNonInterpolatingOutOfRangeSize = 0;
			for(PetscInt j = 0; j<row.columns_number; j++)
			{
				if(row.columns[j] != row_index && isInFSetOutOfRange(row.columns[j], rowData, first_row_in_range, range_end, influenced_info_array) &&
						isStrongDependence(row.values[j], min, row_index, row.columns[j], strength_threshold))
				{
					strongNonInterpolatingOutOfRange[strongNonInterpolatingOutOfRangeSize++] = row.columns[j];
					potential_strong_ff[rowData->GetRowDataOwner(row.columns[j])]++;
					printf("STRONG FF %d %d owner %d\n", row_index, row.columns[j], rowData->GetRowDataOwner(row.columns[j]));
				}
				if(row.columns[j] != row_index && isInCSet(row.columns[j], rowData, first_row_in_range, range_end, influenced_info_array) &&
						isStrongDependence(row.values[j], min, row_index, row.columns[j], strength_threshold))
				{
					strongInterpolating[strongInterpolatingSize++] = row.columns[j];
				}
			}

			if(strongNonInterpolatingOutOfRangeSize > 0)
			{
//				for(int j = 0; i<strongNonInterpolatingOutOfRangeSize; j++)
//					potential_strong_ff[rowData->GetRowDataOwner(row.columns[strongNonInterpolatingOutOfRange[j])] +=  (2 + strongInterpolatingSize);
			}
			restoreFilteredRow(mat,row_index,&row);
		}
	}
}

Mat DirectInterpolation::GetMatrixFromCoarseToFine()
{
	int coarse_columns_number = SetCoarseRowsNumber();
	RowDataCollector* collector = new RowDataCollector();
	RowData* rowData = collector->Collect(mat, influenced_info_array);
	//Handle3RS(rowData);
	MPI_Barrier(PETSC_COMM_WORLD);

	PetscErrorCode ierr;
	Mat from_coarse_to_fine_mat;
	Vec coarse_rhs;
	PetscInt coarse_rows_number;
	MatGetSize(mat,&coarse_rows_number,NULL);

	PetscInt ownership_range_size = range_end - first_row_in_range;

	PrepareFromCoarseToFineMat(&from_coarse_to_fine_mat, coarse_rows_number, coarse_columns_number, rowData);

	const PetscInt* columns;
	const PetscScalar* values;
	PetscInt columns_number;
	filtered_row row;

	PetscScalar one = 1.0;
	PetscInt row_index;

	PetscInt* interpolatory;
	PetscInt* strong_non_interpolatory;
	PetscInt* weak_non_interpolatory;

	for(PetscInt i = 0; i<ownership_range_size; i++)
	{
		row_index = i + first_row_in_range;
		if(influenced_info_array[i].row_info->set == FSET)
		{
//			MatGetRow(mat,i + first_row_in_range,&columns_number,&columns,&values);
			getFilteredRow(mat,i + first_row_in_range,&row);
			columns_number = row.columns_number;
			columns = row.columns;
			values = row.values;
			//TODO: Change sum ranges
			PetscScalar row_sum = getSum(values, columns_number, i+first_row_in_range, columns);
			double min = getMin(values, columns_number, i+first_row_in_range,columns, first_row_in_range, range_end);
			PetscScalar interpolating_sum = 0;
			PetscScalar value_on_diagonal;
			for(PetscInt j = 0; j<columns_number; j++)
			{
				//TODO: change
				if(columns[j] != i+first_row_in_range && isInCSet(columns[j], rowData, first_row_in_range, range_end, influenced_info_array) &&
						isStrongDependence(values[j], min, i + first_row_in_range, columns[j], strength_threshold))
				{
					//mf_log_info("%s %lf %d","adding to interpolating sum ", values[j], columns[j]);
					interpolating_sum += values[j];
					//mf_log_info("%s %lf","interpolating sum ", interpolating_sum);
				}
				else if(columns[j] == (i+first_row_in_range))
				{
					//mf_log_info("%s %lf","value on diagonal ", values[j]);
					value_on_diagonal = values[j];
				}

			}

			PetscScalar value;
			for(PetscInt j = 0; j<columns_number; j++)
			{

				//TODO: change
				if(columns[j] != i+first_row_in_range && isInCSet(columns[j], rowData, first_row_in_range, range_end, influenced_info_array) &&
						isStrongDependence(values[j], min, i + first_row_in_range, columns[j], strength_threshold))
				{

					//printf("%lf %lf %lf %lf\n", row_sum, interpolating_sum, values[j], value_on_diagonal);
					value = -(row_sum/interpolating_sum)*values[j]/value_on_diagonal;
					int global_row_number_in_coarse = getColumnIndex(columns[j], rowData, first_row_in_range, range_end,
							influenced_info_array);
					//printf("%d global %d\n", global_row_number_in_coarse, rank);
	
					ierr = MatSetValues(from_coarse_to_fine_mat,1,&row_index,1,&global_row_number_in_coarse,&value,INSERT_VALUES);
					CHKERRABORT(PETSC_COMM_WORLD, ierr);
				}
			}
			restoreFilteredRow(mat,i + first_row_in_range,&row);
//			MatRestoreRow(mat,i + first_row_in_range,&columns_number,&columns,&values);
		}
		else
		{
			//printf("%d is C\n", row_index);
			int global_row_number_in_coarse = influenced_info_array[i].row_number_in_coarse;
			//printf("%d global\n", global_row_number_in_coarse);

			ierr = MatSetValues(from_coarse_to_fine_mat,1,&row_index,1,&(global_row_number_in_coarse),&one,INSERT_VALUES);

			CHKERRABORT(PETSC_COMM_WORLD, ierr);
		}
	}

	delete[] row_info_array;
	delete[] influenced_info_array;

	ierr = MatAssemblyBegin(from_coarse_to_fine_mat,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatAssemblyEnd(from_coarse_to_fine_mat,MAT_FINAL_ASSEMBLY);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	//MatView(from_coarse_to_fine_mat, PETSC_VIEWER_STDOUT_WORLD);

	return from_coarse_to_fine_mat;
}
