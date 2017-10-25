/*
 * InterpolationStrategy.cpp
 *
 *  Created on: Jan 22, 2016
 *      Author: damian
 */

#include "InterpolationStrategy.hpp"


InterpolationStrategy::InterpolationStrategy(Mat mat, struct row_info* row_info_array, struct influenced_info* influenced_info_array,
		PetscInt first_row_in_range, PetscInt range_end, double strength_threshold) : mat(mat), row_info_array(row_info_array),
		influenced_info_array(influenced_info_array), first_row_in_range(first_row_in_range), range_end(range_end), strength_threshold(strength_threshold){
	this->influenced_info_array = influenced_info_array;

}

InterpolationStrategy::InterpolationStrategy(){

}

InterpolationStrategy::~InterpolationStrategy() {
	// TODO Auto-generated destructor stub
}

void InterpolationStrategy::InitStrategy(Mat mat, struct row_info* row_info_array, struct influenced_info* influenced_info_array,
		PetscInt first_row_in_range, PetscInt range_end, double strength_threshold){
	this->mat = mat;
	this->row_info_array = row_info_array;
	this->influenced_info_array = influenced_info_array;
	this->first_row_in_range = first_row_in_range;
	this->range_end = range_end;
	this->strength_threshold = strength_threshold;
}


void InterpolationStrategy::PrepareFromCoarseToFineMat(Mat* from_coarse_to_fine_mat, int coarse_rows_number, int coarse_columns_number, RowData* rowData)
{

	PetscErrorCode ierr;
	const PetscInt* columns;
	const PetscScalar* values;
	PetscInt columns_number;
	filtered_row row;

	int mat_local_rows_number, mat_local_cols_number;
	MatGetLocalSize(mat,&mat_local_rows_number, &mat_local_cols_number);

	ierr = MatCreate(PETSC_COMM_WORLD,from_coarse_to_fine_mat);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatSetSizes(*from_coarse_to_fine_mat,mat_local_rows_number,PETSC_DECIDE,coarse_rows_number, coarse_columns_number);CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatSetFromOptions(*from_coarse_to_fine_mat);CHKERRABORT(PETSC_COMM_WORLD, ierr);

	int from_coarse_to_fine_ownership_column_begin;
	int from_coarse_to_fine_ownership_column_end;

	int rank, world_size;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &world_size);

	//info: in accordance with petsc split ownership range formula
	from_coarse_to_fine_ownership_column_begin = coarse_columns_number / world_size * rank + std::min(coarse_columns_number % world_size, rank);
	int local_size = PETSC_DECIDE;
	PetscSplitOwnership(PETSC_COMM_WORLD,&local_size,&coarse_columns_number);
	from_coarse_to_fine_ownership_column_end = from_coarse_to_fine_ownership_column_begin + local_size;
	//TODO: workaround for seq
	//VecGetOwnershipRange
//	PetscInt localColumns;
//	PetscInt localRows;
//	MatGetLocalSize(*from_coarse_to_fine_mat,&localRows,&localColumns);
//	printf("LOcALS %d %d\n", localRows, localColumns);
//	ierr = MatGetOwnershipRangeColumn(*from_coarse_to_fine_mat,
//			&from_coarse_to_fine_ownership_column_begin,
//			&from_coarse_to_fine_ownership_column_end); CHKERRABORT(PETSC_COMM_WORLD, ierr);


	int* diagonal_nonzero_numbers = (int* )calloc(coarse_rows_number, sizeof(int));
	int* nondiagonal_nonzero_numbers = (int *)calloc(coarse_rows_number, sizeof(int));

	for(PetscInt i = 0; i<range_end - first_row_in_range; i++)
	{

		nondiagonal_nonzero_numbers[i] = 0;
		diagonal_nonzero_numbers[i] = 0;
		if(influenced_info_array[i].row_info->set == CSET)
		{
			increase_nonzero_number(influenced_info_array[i].row_number_in_coarse,i,
					from_coarse_to_fine_ownership_column_begin, from_coarse_to_fine_ownership_column_end,
					nondiagonal_nonzero_numbers, diagonal_nonzero_numbers);
		}
		else
		{
			getFilteredRow(mat,i + first_row_in_range,&row);
			columns_number = row.columns_number;
			columns = row.columns;
//			MatGetRow(mat,i + first_row_in_range,&columns_number,&columns,NULL);
			for(int j = 0; j<columns_number; j++)
			{
				if(columns[j] >= first_row_in_range && columns[j] < range_end){
					if(influenced_info_array[columns[j] - first_row_in_range].row_info->set == CSET)
					{
						increase_nonzero_number(influenced_info_array[columns[j] - first_row_in_range].row_number_in_coarse,i,
								from_coarse_to_fine_ownership_column_begin, from_coarse_to_fine_ownership_column_end,
								nondiagonal_nonzero_numbers, diagonal_nonzero_numbers);
					}
				}
				else{
					int columnIndex = rowData->GetRowData(columns[j]);
					if(columnIndex != -1 ){
						if(columnIndex >= from_coarse_to_fine_ownership_column_begin && columnIndex < from_coarse_to_fine_ownership_column_end){
							diagonal_nonzero_numbers[i]++;
						}
						else{
							nondiagonal_nonzero_numbers[i]++;
						}
					}
				}
			}
//			MatRestoreRow(mat,i + first_row_in_range,&columns_number,&columns,NULL);
			restoreFilteredRow(mat, i + first_row_in_range, &row);
		}
//		if(nondiagonal_nonzero_numbers[i] == 0 && diagonal_nonzero_numbers[i] == 0){
//			printf("ZERO ROW %d %d\n", i, columns_number);
//			getFilteredRow(mat,i + first_row_in_range,&row);
//			for(int i=0; i<row.columns_number; i++){
//				printf("%d ", row.values[i]);
//			}
//			printf("\n");
//			restoreFilteredRow(mat, i + first_row_in_range, &row);
//		}
	}

	ierr = MatMPIAIJSetPreallocation(*from_coarse_to_fine_mat,0,diagonal_nonzero_numbers ,0,nondiagonal_nonzero_numbers );
	CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatSeqAIJSetPreallocation(*from_coarse_to_fine_mat,0,diagonal_nonzero_numbers);
	CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatSetUp(*from_coarse_to_fine_mat);CHKERRABORT(PETSC_COMM_WORLD, ierr);
//	getchar();
	free(diagonal_nonzero_numbers);
	free(nondiagonal_nonzero_numbers);
}

