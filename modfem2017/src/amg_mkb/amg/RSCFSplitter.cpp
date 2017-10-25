/*
 * RSCFSplitter.cpp
 *
 *  Created on: Jun 6, 2015
 *      Author: damian
 */

#include "RSCFSplitter.hpp"

RSCFSplitter::RSCFSplitter(Mat mat, double strength_threshold, InterpolationStrategy* strategy) : CFSplitter(mat)
{
	this->strength_threshold = strength_threshold;
	this->strategy = strategy;
}

RSCFSplitter::~RSCFSplitter()
{
}

PetscScalar RSCFSplitter::GetExtremalValue(const PetscScalar* values, PetscInt size, PetscInt row_index,
		const PetscInt* columns, int first_row_in_range, int range_end)
{
	return getMin(values, size, row_index, columns, first_row_in_range, range_end);
}

bool RSCFSplitter::IsStrongDependenceWithinRange(PetscScalar value, PetscScalar row_min, PetscInt row_index, PetscInt column_index,
		PetscInt first_row_in_range, PetscInt range_end, double strength_threshold)
{
	return isStrongDependenceWithinRange(value, row_min, row_index, column_index, first_row_in_range, range_end, strength_threshold);
}

void RSCFSplitter::MakeCFSplitting()
{
	//mf_log_info("%s","CF Splitting - beginning ");
	PetscErrorCode ierr;
	ierr = MatGetOwnershipRange(mat,&first_row_in_range,&range_end); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	//mf_log_info("%s %d %d","range: ", first_row_in_range, range_end);

	PetscInt ownership_range_size = range_end - first_row_in_range;

	const PetscInt* columns;
	const PetscScalar* values;
	PetscInt columns_number;
	filtered_row row;


	row_info_array = new struct row_info[ownership_range_size];
	influenced_info_array = new struct influenced_info[ownership_range_size];


	for(PetscInt i = first_row_in_range; i<range_end; i++)
	{
		row_info_array[i - first_row_in_range].local_row_number = i - first_row_in_range;
		influenced_info_array[i - first_row_in_range].row_info = &row_info_array[i - first_row_in_range];
	}

	//compute initial influence

	for(PetscInt i = first_row_in_range; i<range_end; i++)
	{

		//mf_log_info("%s %d","Numebr of row: ", i);

//		MatGetRow(mat,i,&columns_number,&columns,&values);
		getFilteredRow(mat,i,&row);
		columns_number = row.columns_number;
		columns = row.columns;
		values = row.values;

		PetscScalar row_min = GetExtremalValue(values, columns_number, i, columns, first_row_in_range, range_end);

		//mf_log_info("%s %lf","Min value: ", row_min);
		for(PetscInt j = 0; j<columns_number; j++)
		{
			//mf_log_info("%s %d %d %lf","(row, column, value) ", i, columns[j], values[j]);
			PetscScalar value = values[j];
			if(IsStrongDependenceWithinRange(value, row_min, i, columns[j], first_row_in_range, range_end, strength_threshold))
			{
				//mf_log_info("%s %d %d %lf","(row, column, value) is strong ", i, columns[j], values[j]);
				influenced_info_array[columns[j] - first_row_in_range].influenced_number++;
				ReplaceAfterInfluencedNumberChange(influenced_info_array[columns[j] - first_row_in_range].row_info, true);

			}
		}
		restoreFilteredRow(mat, i, &row);
//		MatRestoreRow(mat,i,&columns_number,&columns,&values);
	}
//	getchar();

	int c_rows_count = 0;
	//setting c from dirichlet bc
//	if(dirichlet_bc != NULL)
//		for(PetscInt i = first_row_in_range; i<range_end; i++)
//		{
//			if(dirichlet_bc[i - first_row_in_range])
//			{
//				row_info_array[i - first_row_in_range].set = CSET;
//			}
//		}

	//
	//c-f splitting P1
	//mf_log_info("%s ","Reducing the heap");
	int rows_on_heap = ownership_range_size;
	while(rows_on_heap)
	{
		//mf_log_info("%s %d","Heap size", rows_on_heap);
		row_info* highest_influence_row_info = &row_info_array[0];
		row_info* last_heap_row_info = &row_info_array[rows_on_heap - 1];

		PetscInt row_number = highest_influence_row_info->local_row_number + first_row_in_range;
		if(highest_influence_row_info->set == DEFAULT)
		{
			//mf_log_info("%s %d","Setting c on ", highest_influence_row_info->local_row_number);
			highest_influence_row_info->set = CSET;
			c_rows_count++;
			PetscInt row_index = highest_influence_row_info->local_row_number + first_row_in_range;
//			MatGetRow(mat,row_index,&columns_number,&columns,&values);
			getFilteredRow(mat,row_index,&row);
			columns_number = row.columns_number;
			columns = row.columns;
			values = row.values;


			PetscInt neighboring[columns_number];
			PetscScalar row_min = GetExtremalValue(values, columns_number, row_index, columns, first_row_in_range, range_end);
			//set to F all strong neighbors that are not C and are in range
			PetscInt column_index = 0;
			for(PetscInt j = 0; j<columns_number; j++)
			{
				if(IsStrongDependenceWithinRange(values[j],row_min,row_index,columns[j], first_row_in_range, range_end, strength_threshold))
				{
					neighboring[column_index++] = columns[j];
					row_info* neighboring_row_info = influenced_info_array[columns[j] - first_row_in_range].row_info;
					if(neighboring_row_info->set == DEFAULT)
					{
						neighboring_row_info->set = FSET;
					}
				}
			}
			PetscInt neighboring_number = column_index;

			restoreFilteredRow(mat, highest_influence_row_info->local_row_number + first_row_in_range, &row);
//			MatRestoreRow(mat,highest_influence_row_info->local_row_number + first_row_in_range,&columns_number,&columns,&values);

			for(PetscInt j = 0; j<neighboring_number; j++)
			{

//				MatGetRow(mat,neighboring[j],&columns_number,&columns,&values);
				getFilteredRow(mat,neighboring[j],&row);
				columns_number = row.columns_number;
				columns = row.columns;
				values = row.values;

				PetscScalar row_min = GetExtremalValue(values, columns_number, neighboring[j], columns, first_row_in_range, range_end);
				//printf("A neighobor %d with min %lf\n",neighboring[j],row_min);
				for(PetscInt k = 0; k<columns_number; k++)
				{
					if(IsStrongDependenceWithinRange(values[k],row_min,neighboring[j],columns[k], first_row_in_range, range_end, strength_threshold))
					{
						row_info* neighbor_neighboring_row_info = influenced_info_array[columns[k] - first_row_in_range].row_info;
						if(neighbor_neighboring_row_info->set == DEFAULT)
						{
							influenced_info_array[columns[k] - first_row_in_range].influenced_number++;
							//mf_log_info("%s %d %d","Incrementing (on, to) ", columns[k] - first_row_in_range, influenced_info_array[columns[k] - first_row_in_range].influenced_number);
							ReplaceAfterInfluencedNumberChange(influenced_info_array[columns[k] - first_row_in_range].row_info, false);
						}
					}
				}
//				MatRestoreRow(mat,neighboring[j],&columns_number,&columns,&values);
				restoreFilteredRow(mat, neighboring[j], &row);
			}

		}
		else
		{
			//just remove highest influence from the heap
		}

		if(rows_on_heap != 1)
		{
			replaceRowInfos(highest_influence_row_info, last_heap_row_info);
		}

		rows_on_heap--;
		makeHeap(rows_on_heap);
	}
//	printf("After the first step: %d of C rows\n", c_rows_count);
	//c-f splitting P2 (strong F-F dependencies)
	//mf_log_info("%s  ", "C F splitting II - removing FF dependencies");

	for(PetscInt i = 0; i<ownership_range_size; i++)
	{
		if(row_info_array[i].set == FSET)
		{
//			MatGetRow(mat,row_info_array[i].local_row_number + first_row_in_range,&columns_number,&columns,&values);
			getFilteredRow(mat,row_info_array[i].local_row_number + first_row_in_range,&row);
			columns_number = row.columns_number;
			columns = row.columns;
			values = row.values;
			PetscScalar row_min = GetExtremalValue(values, columns_number, row_info_array[i].local_row_number + first_row_in_range, columns, first_row_in_range, range_end);

			PetscInt strongInterpolating[columns_number];
			PetscInt strongNonInterpolating[columns_number];
			PetscInt strongInterpolatingSize = 0;
			PetscInt strongNonInterpolatingSize = 0;


			for(PetscInt j = 0; j<columns_number; j++)
			{
				if(IsStrongDependenceWithinRange(values[j],row_min,row_info_array[i].local_row_number + first_row_in_range,columns[j],
						first_row_in_range, range_end, strength_threshold))
				{

					//to moze byc wyzej
					//mf_log_info("%s %d %d","Strong dependence ", row_info_array[i].local_row_number, columns[j]);
					if(influenced_info_array[columns[j] - first_row_in_range].row_info->set == CSET)
					{
						strongInterpolating[strongInterpolatingSize++] = columns[j];
					}
					//
					else if(influenced_info_array[columns[j] - first_row_in_range].row_info->set == FSET)
					{
						strongNonInterpolating[strongNonInterpolatingSize++] = columns[j];
					}
					else
					{
						mf_log_err("%s","Error in C-F splitting");
						throw new std::runtime_error("Some variables neither in c nor in f");
						MPI_Abort(PETSC_COMM_WORLD, -1);
					}
				}
			}
			restoreFilteredRow(mat, row_info_array[i].local_row_number + first_row_in_range, &row);
//			MatRestoreRow(mat,row_info_array[i].local_row_number + first_row_in_range,&columns_number,&columns,&values);
			bool addToCSet = strongNonInterpolatingSize > 0;
			for(PetscInt j = 0; j<strongNonInterpolatingSize; j++)
			{
				addToCSet = true;
				//mf_log_info("%s %d %d","Taking ",row_info_array[i].local_row_number, strongNonInterpolating[j]);
//				MatGetRow(mat,strongNonInterpolating[j],&columns_number,&columns,&values);
				getFilteredRow(mat,strongNonInterpolating[j],&row);
				columns_number = row.columns_number;
				columns = row.columns;
				values = row.values;

				PetscScalar row_min = GetExtremalValue(values, columns_number, strongNonInterpolating[j], columns, first_row_in_range, range_end);
				for(PetscInt k = 0; k<columns_number; k++)
				{
					//for(PetscInt z = 0; z<columns_number; z++)
					//{
					//	printf("z %d %lf\n",columns[z], values[z]);
					//}
					//mf_log_info("%s %d %d %lf","Taking ",strongNonInterpolating[j], columns[k], values[k]);
					if(IsStrongDependenceWithinRange(values[k],row_min,strongNonInterpolating[j],columns[k],
							first_row_in_range, range_end, strength_threshold))
					{
						for(PetscInt l = 0; l<strongInterpolatingSize; l++)
						{
							if(columns[k] == strongInterpolating[l])
							{
								//mf_log_info("%s %d %d","Cancelling adding to C ", columns[k], strongInterpolating[l]);
								addToCSet = false;
								break;
							}
						}
					}
				}

//				MatRestoreRow(mat,strongNonInterpolating[j],&columns_number,&columns,&values);
				restoreFilteredRow(mat, strongNonInterpolating[j], &row);
				if(addToCSet)
				{
					break;
				}
			}
			if(addToCSet)
			{
				row_info_array[i].set = CSET;
			}
		}
	}

	for(PetscInt i = first_row_in_range; i<range_end; i++)
	{
		if(row_info_array[i - first_row_in_range].set == DEFAULT)
		{
			throw new std::runtime_error("Some variables neither in c nor in f");
			MPI_Abort(PETSC_COMM_WORLD, -1);
		}
	}
}

struct row_info* RSCFSplitter::getHeapParent(int index)
{
	return &row_info_array[(index+1)/2 - 1];
}

void RSCFSplitter::makeHeap(int size)
{
	if(size == 0)
		return;

	int current_row_info_index = 0;
	while(true)
	{
		//mf_log_info("%s %d","current row info ", current_row_info_index);
		int left_sibling = current_row_info_index*2 + 1;
		int right_sibling = left_sibling + 1;

		if(left_sibling >= size)
			return;

		int current_influenced = influenced_info_array[row_info_array[current_row_info_index].local_row_number].influenced_number;
		int left_influenced = influenced_info_array[row_info_array[left_sibling].local_row_number].influenced_number;
		int right_influenced = -1;

		if(right_sibling < size)
			right_influenced = influenced_info_array[row_info_array[right_sibling].local_row_number].influenced_number;

		int higher_influence_index = left_sibling;
		if(right_influenced > left_influenced)
			higher_influence_index = right_sibling;

		if(left_influenced > current_influenced || right_influenced > current_influenced)
		{
			replaceRowInfos(&row_info_array[current_row_info_index], &row_info_array[higher_influence_index]);
			current_row_info_index = higher_influence_index;
		}
		else
		{
			return;
		}
	}
}

void RSCFSplitter::replaceRowInfos(struct row_info* row_info1, struct row_info* row_info2)
{
	influenced_info_array[row_info1->local_row_number].row_info = row_info2;
	influenced_info_array[row_info2->local_row_number].row_info = row_info1;

	Set set = row_info1->set;
	PetscInt tmp_local_row_number = row_info1->local_row_number;

	row_info1->set = row_info2->set;
	row_info1->local_row_number = row_info2->local_row_number;

	row_info2->set = set;
	row_info2->local_row_number = tmp_local_row_number;
}

void RSCFSplitter::ReplaceAfterInfluencedNumberChange(struct row_info* row_info, bool switch_root)
{
	while(true)
	{
		int index_in_heap = row_info - row_info_array;
		//mf_log_info("%s %d","index in heap ", index_in_heap);
		if(index_in_heap == 0)
			return;
		if(!switch_root && index_in_heap < 3)
			return;

		struct row_info* parent_row_info = getHeapParent(index_in_heap);
		if(influenced_info_array[row_info->local_row_number].influenced_number >
					influenced_info_array[parent_row_info->local_row_number].influenced_number)
		{
			replaceRowInfos(row_info, parent_row_info);
			row_info = parent_row_info;
		}
		else
		{
			return;
		}
	}
}

void RSCFSplitter::GetSetsWithInfluencedNumebr(std::map<int, int> *c_set,std::map<int, int> *f_set)
{
	for(PetscInt i = 0; i<range_end-first_row_in_range; i++)
	{
		if(row_info_array[i].set == CSET)
		{
			(*c_set)[row_info_array[i].local_row_number] = influenced_info_array[row_info_array[i].local_row_number].influenced_number;
		}
		else
		{
			(*f_set)[row_info_array[i].local_row_number] = influenced_info_array[row_info_array[i].local_row_number].influenced_number;
		}
	}

}

std::list<int>* RSCFSplitter::GetNumbersOfCoarseRows()
{
	if(influenced_info_array != NULL)
	{
		std::list<int>* numbersOfCoarseRows = new std::list<int>();
		for(int i = 0; i<range_end - first_row_in_range; i++)
		{
			if(influenced_info_array[i].row_info->set == CSET)
				numbersOfCoarseRows->push_back(i);
		}
		return numbersOfCoarseRows;
	}
	return new std::list<int>();
}

Mat RSCFSplitter::GetMatrixFromCoarseToFine()
{
	//TODO: interpolation strategy option
//	InterpolationStrategy* strategy = new DirectInterpolation(mat, row_info_array, influenced_info_array,
//		first_row_in_range, range_end, strength_threshold);
//	InterpolationStrategy* strategy = new AdvancedRSCFInterpolation(mat, row_info_array, influenced_info_array,
//			first_row_in_range, range_end, strength_threshold);
	strategy->InitStrategy(mat, row_info_array, influenced_info_array,first_row_in_range, range_end, strength_threshold);
	Mat mat = strategy->GetMatrixFromCoarseToFine();
	delete strategy;
	return mat;
}

