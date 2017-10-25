/*
 * AdvancedRSCFInterpolation.cpp
 *
 *  Created on: Jan 23, 2016
 *      Author: damian
 */

#include "AdvancedRSCFInterpolation.hpp"

AdvancedRSCFInterpolation::AdvancedRSCFInterpolation(Mat mat, struct row_info* row_info_array, struct influenced_info* influenced_info_array,
		PetscInt first_row_in_range, PetscInt range_end, double strength_threshold) : InterpolationStrategy(mat, row_info_array, influenced_info_array,
				first_row_in_range, range_end, strength_threshold)
{
}

AdvancedRSCFInterpolation::~AdvancedRSCFInterpolation()
{

}

bool isInFSetWeak(PetscInt row_index, PetscInt column_index, PetscInt first_row_in_range, PetscInt range_end,
		struct influenced_info* influenced_info_array, double value, double row_min, double strength_threshold)
{
	if(influenced_info_array[column_index].row_info->set == FSET)
		return !isStrongDependenceWithinRange(value, row_min, row_index, column_index, first_row_in_range, range_end,
				strength_threshold);
	return false;
}


Mat AdvancedRSCFInterpolation::GetMatrixFromCoarseToFine()
{
	PetscInt ownership_range_size = range_end - first_row_in_range;
	PetscInt column_number_in_coarse = 0;
	for(PetscInt i = 0; i<ownership_range_size; i++)
	{
		if(influenced_info_array[i].row_info->set == CSET)
		{
			influenced_info_array[i].row_number_in_coarse = column_number_in_coarse++;
		}
	}

	PetscErrorCode ierr;
	Mat from_coarse_to_fine_mat;
	Vec coarse_rhs;
	PetscInt coarse_rows_number;
	MatGetSize(mat,&coarse_rows_number,NULL);
	PetscInt coarse_columns_number = column_number_in_coarse;
	printf("coarse rows %d\n",column_number_in_coarse);

	PrepareFromCoarseToFineMat(&from_coarse_to_fine_mat, coarse_rows_number, coarse_columns_number, NULL);

	const PetscInt* columns;
	const PetscScalar* values;
	PetscInt columns_number;
	filtered_row row;


	PetscScalar one = 1.0;
	PetscInt row_index;

	PetscInt* interpolatory;
	PetscInt* strong_non_interpolatory;
	PetscInt* weak_non_interpolatory;

	double master_row_value_buffer[MAX_ROW_SIZE];
	int master_row_column_buffer[MAX_ROW_SIZE];
	int master_row_size;

	RowData* rowData = NULL;
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
			double row_min = getMin(values, columns_number, i+first_row_in_range, columns, first_row_in_range, range_end);
			double weigth_scale = 0;
			for(PetscInt j = 0; j<columns_number; j++)
			{
				if(columns[j] == i + first_row_in_range)
				{
					weigth_scale += values[j];
				}
				else if(isInFSetWeak(i+first_row_in_range, columns[j], first_row_in_range, range_end,
						influenced_info_array, values[j], row_min, strength_threshold))
				{
					weigth_scale += values[j];
				}
			}
			weigth_scale = - 1.0 / weigth_scale;

			PetscInt strongInterpolating[columns_number];
			PetscInt strongNonInterpolating[columns_number];
			PetscInt strongInterpolatingSize = 0;
			PetscInt strongNonInterpolatingSize = 0;

			for(PetscInt j = 0; j<columns_number; j++)
			{
				if(isStrongDependenceWithinRange(values[j],row_min,row_index,
						columns[j],first_row_in_range, range_end, strength_threshold))
				{
					if(influenced_info_array[columns[j] - first_row_in_range].row_info->set == CSET)
					{
						strongInterpolating[strongInterpolatingSize++] = columns[j];
					}
					//
					else if(influenced_info_array[columns[j] - first_row_in_range].row_info->set == FSET)
					{
						strongNonInterpolating[strongNonInterpolatingSize++] = columns[j];
					}
				}
			}

			if(columns_number > MAX_ROW_SIZE)
			{
				printf("%d %d columns number, max row size\n",columns_number, MAX_ROW_SIZE );
				throw new std::runtime_error("max row buffer size too small");
			}

			memcpy(master_row_value_buffer, values, columns_number * sizeof(PetscScalar));
			memcpy(master_row_column_buffer, columns, columns_number * sizeof(PetscInt));
			master_row_size = columns_number;
			restoreFilteredRow(mat, i + first_row_in_range, &row);
			//MatRestoreRow(mat,i + first_row_in_range,&columns_number,&columns,&values);

			for(int master_row_column_index = 0; master_row_column_index < master_row_size; master_row_column_index++)
			{
				if(master_row_column_buffer[master_row_column_index] != i+first_row_in_range && isInCSet(master_row_column_buffer[master_row_column_index],
																			rowData, first_row_in_range, range_end, influenced_info_array))
				{
					double value;
					MatGetValues(mat,1,&row_index,1,&master_row_column_buffer[master_row_column_index],&value);
					for(int j = 0; j<strongNonInterpolatingSize; j++)
					{
//						MatGetRow(mat,strongNonInterpolating[j],&columns_number,&columns,&values);
						getFilteredRow(mat,strongNonInterpolating[j],&row);
						columns_number = row.columns_number;
						columns = row.columns;
						values = row.values;
						double denominator = 0;
						bool change = false;
						for(int k = 0; k<columns_number; k++)
						{
							for(int l = 0; l<strongInterpolatingSize; l++)
							{
//								printf("strong non interpolating column nbr %lf\n",strongNonInterpolating[j]);
								if(columns[k] == strongInterpolating[l])
								{
//									if(values[k] < 0){
										denominator += values[k];
										change = true;
//										printf("DENOMINATOR ATER CHANGE  %lf\n",denominator);
//									}
								}
							}
						}
						double master_row_value;
						double strong_non_interpolating_row_value;
						MatGetValues(mat,1,&row_index,1,&strongNonInterpolating[j],&master_row_value);
						MatGetValues(mat,1,&strongNonInterpolating[j],1,&master_row_column_buffer[master_row_column_index],&strong_non_interpolating_row_value);
						//printf("DENOMINATOR %lf %d %d\n", denominator, change, strongNonInterpolating[j]);
//						double kkk = abs(denominator);
//						if(abs(denominator) < 0.001 && change)
//							denominator = 0.001;
//						if(strong_non_interpolating_row_value > 0)
//							strong_non_interpolating_row_value = 0;
						value += master_row_value*strong_non_interpolating_row_value / denominator;
//						printf("denominator %lf\n", denominator);
						restoreFilteredRow(mat, strongNonInterpolating[j], &row);
//						MatRestoreRow(mat,strongNonInterpolating[j],&columns_number,&columns,&values);
					}
					value *= weigth_scale;
//					printf("%lf\n",value);

					int global_row_number_in_coarse = getColumnIndex(master_row_column_buffer[master_row_column_index], rowData, first_row_in_range, range_end,
												influenced_info_array);
					if(value != value){
						printf("STOP\n");
					}
					ierr = MatSetValues(from_coarse_to_fine_mat,1,&row_index,1,&global_row_number_in_coarse,&value,INSERT_VALUES);
										CHKERRABORT(PETSC_COMM_WORLD, ierr);
				}
			}
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
	return from_coarse_to_fine_mat;
}


