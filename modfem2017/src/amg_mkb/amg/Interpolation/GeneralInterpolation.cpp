/*
 * GeneralInterpolation.cpp
 *
 *  Created on: Feb 3, 2017
 *      Author: damian
 */

#include "GeneralInterpolation.hpp"

GeneralInterpolation::GeneralInterpolation() {
	// TODO Auto-generated constructor stub

}

GeneralInterpolation::~GeneralInterpolation() {
	// TODO Auto-generated destructor stub
}

void getStrongCountAndValuesSums(const int* columns1, const double* values1, int array1_size,
		const int* columns2, const double* values2, int array2_size, int* count, double* values_sum, double* absolute_values_sum){
	int index1 = 0;
	int index2 = 0;
	int cnt = 0;
	double v_sum = 0;
	double abs_v_sum = 0;
	while(index1 < array1_size && index2 < array2_size){
		if(columns1[index1] < columns2[index2])
			index1++;
		else if(columns1[index1] > columns2[index2])
			index2++;
		else{
			v_sum += values2[index2];
			abs_v_sum += fabs(values2[index2]);
			index1++; index2++;
			cnt++;
		}
	}
	*count = cnt;
	*values_sum = v_sum;
	*absolute_values_sum = abs_v_sum;
}

Mat GeneralInterpolation::GetMatrixFromCoarseToFine(){
	RowData* rowData = NULL;
	int coarse_columns_number = 0;
	PetscInt ownership_range_size = range_end - first_row_in_range;
	for(PetscInt i = 0; i<ownership_range_size; i++)
	{
		if(influenced_info_array[i].row_info->set == CSET)
		{
			influenced_info_array[i].row_number_in_coarse = coarse_columns_number++;
		}
	}

	PetscErrorCode ierr;
	Mat from_coarse_to_fine_mat;
	Vec coarse_rhs;
	PetscInt coarse_rows_number;
	MatGetSize(mat,&coarse_rows_number,NULL);

	PrepareFromCoarseToFineMat(&from_coarse_to_fine_mat, coarse_rows_number, coarse_columns_number, NULL);

	const PetscInt* columns;
	const PetscScalar* values;
	PetscInt columns_number;
	filtered_row row;

	const PetscInt* helper_columns;
	const PetscScalar* helper_values;
	PetscInt helper_columns_number;
	filtered_row helper_row;

	PetscScalar one = 1.0;
	PetscInt row_index;

	PetscInt* interpolatory;
	PetscInt* strong_non_interpolatory;
	PetscInt* weak_non_interpolatory;

	int nonzeros_number_in_interpolating_columns[max_unfiltered_row_size];
	int strong_interpolating_indices[max_unfiltered_row_size];
	double values_sum_in_interpolating_columns[max_unfiltered_row_size];
	double absolute_values_sum_in_interpolating_columns[max_unfiltered_row_size];

	for(PetscInt i = 0; i<ownership_range_size; i++)
	{
		row_index = i + first_row_in_range;
		if(influenced_info_array[i].row_info->set == FSET)
		{
			int nonzeros_number_in_interpolating_columns_size = 0;
			int strong_interpolating_indices_size = 0;
			int values_sum_in_interpolating_columns_size = 0;

			getFilteredRow(mat,i + first_row_in_range,&row);
			columns_number = row.columns_number;
			columns = row.columns;
			values = row.values;
			double extremal = getExtremalValue(values, columns_number, i+first_row_in_range, columns, first_row_in_range, range_end);
			for(PetscInt j = 0; j<columns_number; j++)
			{
				if(columns[j] != i+first_row_in_range && isInCSet(columns[j], rowData, first_row_in_range, range_end, influenced_info_array) &&
						isGeneralStrongDependence(values[j], extremal, i + first_row_in_range, columns[j], strength_threshold)){
					strong_interpolating_indices[strong_interpolating_indices_size++] = columns[j];
				}
			}

			for(PetscInt j = 0; j<columns_number; j++)
			{
				if(columns[j] != i+first_row_in_range && !isInCSet(columns[j], rowData, first_row_in_range, range_end, influenced_info_array))
				{
					getFilteredRow(mat,columns[j],&helper_row);
					helper_columns_number = helper_row.columns_number;
					helper_columns = helper_row.columns;
					helper_values = helper_row.values;
					int nonzeros_number;
					double values_sum;
					double absolute_values_sum;
					getStrongCountAndValuesSums(strong_interpolating_indices, values, strong_interpolating_indices_size,
							helper_columns, helper_values, helper_columns_number, &nonzeros_number, &values_sum, &absolute_values_sum);

					nonzeros_number_in_interpolating_columns[nonzeros_number_in_interpolating_columns_size] = nonzeros_number;
					values_sum_in_interpolating_columns[nonzeros_number_in_interpolating_columns_size] = values_sum;
					absolute_values_sum_in_interpolating_columns[nonzeros_number_in_interpolating_columns_size] = absolute_values_sum;
					nonzeros_number_in_interpolating_columns_size++;
					restoreFilteredRow(mat,i + first_row_in_range,&helper_row);
				}
			}

			double denominator = 0;
			int non_c_index = 0;
			for(PetscInt k = 0; k<columns_number; k++)
			{
				if(columns[k] == i+first_row_in_range){
					denominator += values[k];
				}
				else if(!isInCSet(columns[k], rowData, first_row_in_range, range_end, influenced_info_array)){
					bool isStrong = isGeneralStrongDependence(values[k], extremal, i + first_row_in_range, columns[k], strength_threshold);
					if(nonzeros_number_in_interpolating_columns_size == 0 || nonzeros_number_in_interpolating_columns[non_c_index] == 0){
						if(isStrong)
						{
							throw new std::runtime_error("Should never happen");
						}
						else{
							denominator -= fabs(values[k]);
						}
					}
					else if(!isStrong){
						double xi = -values_sum_in_interpolating_columns[non_c_index]/absolute_values_sum_in_interpolating_columns[non_c_index];
						if(xi >= 0.5 && values[k] < 0){
							denominator -= values[k];
						}

					}
					else{
						double eta;
						int row = i + first_row_in_range;
						MatGetValues(mat,1,&(columns[k]),1,&row,&eta);
						eta = fabs(eta);
						double xi = -values_sum_in_interpolating_columns[non_c_index]/absolute_values_sum_in_interpolating_columns[non_c_index];
						eta *= nonzeros_number_in_interpolating_columns[non_c_index]/absolute_values_sum_in_interpolating_columns[non_c_index];
						if(eta < 0.75 && xi >= 0.5 && values[k] < 0){
							denominator -= values[k];
						}
						else if(eta > 2 && xi >= 0.5 && values[k] < 0){
							denominator += 0.5*values[k];
						}
					}

					non_c_index++;
				}
			}

			for(PetscInt j = 0; j<columns_number; j++)
			{
				if(columns[j] != i+first_row_in_range && isInCSet(columns[j], rowData, first_row_in_range, range_end, influenced_info_array) &&
					isGeneralStrongDependence(values[j], extremal, i + first_row_in_range, columns[j], strength_threshold))
				{
					double numerator = values[j];
					non_c_index = 0;
					for(PetscInt k = 0; k<columns_number; k++){
						if(!isInCSet(columns[k], rowData, first_row_in_range, range_end, influenced_info_array) && (columns[k] != i + first_row_in_range)){
							if(nonzeros_number_in_interpolating_columns_size == 0 || nonzeros_number_in_interpolating_columns[non_c_index] == 0){
								continue;
							}

							bool isStrong = isGeneralStrongDependence(values[k], extremal, i + first_row_in_range, columns[k], strength_threshold);
							int column = columns[j];
							double v;
							MatGetValues(mat,1,&(columns[k]),1,&column,&v);
							v = fabs(v);
							double interpolation_value = v*values[k]/absolute_values_sum_in_interpolating_columns[non_c_index];
							double xi = -values_sum_in_interpolating_columns[non_c_index]/absolute_values_sum_in_interpolating_columns[non_c_index];
							if(!isStrong){
								if(xi >= 0.5 && values[k] < 0){
									numerator += 2*interpolation_value;
								}
								else{
									numerator += interpolation_value;
								}
							}
							else{
								double eta = v*nonzeros_number_in_interpolating_columns[non_c_index]/absolute_values_sum_in_interpolating_columns[non_c_index];
								if(eta < 0.75 && xi >= 0.5 && values[k] < 0){
									numerator += 2*interpolation_value;
								}
								else if(eta > 2 && xi >= 0.5 && values[k] < 0){
									numerator += 0.5*interpolation_value;
								}
								else{
									numerator+=interpolation_value;
								}
							}
							non_c_index++;
						}
					}
					double value = numerator/denominator;
					int global_row_number_in_coarse = getColumnIndex(columns[j], rowData, first_row_in_range, range_end,
							influenced_info_array);
					MatSetValues(from_coarse_to_fine_mat,1,&row_index,1,&global_row_number_in_coarse, &value, INSERT_VALUES);
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
	//MatView(from_coarse_to_fine_mat, PETSC_VIEWER_STDOUT_WORLD);

	return from_coarse_to_fine_mat;
}
