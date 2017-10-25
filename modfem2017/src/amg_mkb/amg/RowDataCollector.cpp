/*
 * RowDataCollector.cpp
 *
 *  Created on: Oct 17, 2015
 *      Author: damian
 */

#include "RowDataCollector.hpp"
#include "AMGRows.hpp"

void print_buffer(PetscInt* exchangeBuffer, int request_size, int rank, bool send){
	char s[500]{0};
	printf("REQUEST SIZE %d\n", request_size);
	for(int i = 0; i<request_size; i++)
		sprintf(s,"%s %d ",s,exchangeBuffer[i]);
	printf("rank: %d send: %d request: %s\n",rank, send, s);
}

RowData* RowDataCollector::Collect(Mat mat, struct influenced_info* influenced_info_array)
{
	int rank = 0;
	int size = 0;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	RowData* required_row_data = GetRequiredRowData(mat, influenced_info_array);
	int row_data_request_sizes[size];
	required_row_data->GetRowDataRequestSizes(row_data_request_sizes);
	//for(int i = 0; i<size; i++)
		//printf("request sizes %d %d %d\n",i,row_data_request_sizes[i],rank);

	// could be optimized
	int localMaxExchangeBufferSize = GetExchangeBufferSize(row_data_request_sizes, size);
	int globalMaxExchangeBufferSize;
	MPI_Allreduce(&localMaxExchangeBufferSize, &globalMaxExchangeBufferSize, 1,
            MPI_INT, MPI_MAX, PETSC_COMM_WORLD);

	//printf("local exchange buffer size %d\n",localMaxExchangeBufferSize);
	//printf("global exchange buffer size %d\n",globalMaxExchangeBufferSize);
	// could also be optimized
	PetscInt* exchangeBuffer = new PetscInt[globalMaxExchangeBufferSize];
	PetscInt range_begin, range_end;
	MatGetOwnershipRange(mat,&range_begin,&range_end);
	MPI_Status status;
	for(int i = 0; i<size; i++){
		if(i < rank){
			int request_size = required_row_data->FillExchangeRequest(exchangeBuffer, i);
//			print_buffer(exchangeBuffer,request_size, i,true);
			MPI_Send(exchangeBuffer, request_size, MPI_INT, i, 0, PETSC_COMM_WORLD);
			MPI_Recv(exchangeBuffer, request_size, MPI_INT, i, 0, PETSC_COMM_WORLD, &status);
//			print_buffer(exchangeBuffer,request_size, i,false);
			required_row_data->FillExchangeResponse(exchangeBuffer,i);
		}
		if(i > rank){
			MPI_Recv(exchangeBuffer, globalMaxExchangeBufferSize, MPI_INT, i, 0, PETSC_COMM_WORLD, &status);
			int received_count;
			MPI_Get_count(&status,MPI_INT,&received_count);
			FillExchangeBuffer(exchangeBuffer, influenced_info_array, received_count, range_begin);
			MPI_Send(exchangeBuffer, received_count, MPI_INT, i, 0, PETSC_COMM_WORLD);
		}
	}


	for(int i = 0; i<size; i++){
		if(i > rank){
			int request_size = required_row_data->FillExchangeRequest(exchangeBuffer, i);
//			print_buffer(exchangeBuffer,request_size, i,true);
			MPI_Send(exchangeBuffer, request_size, MPI_INT, i, 0, PETSC_COMM_WORLD);
			MPI_Recv(exchangeBuffer, request_size, MPI_INT, i, 0, PETSC_COMM_WORLD, &status);
//			print_buffer(exchangeBuffer,request_size, i,false);
			required_row_data->FillExchangeResponse(exchangeBuffer,i);
		}
		if(i < rank){
			MPI_Recv(exchangeBuffer, globalMaxExchangeBufferSize, MPI_INT, i, 0, PETSC_COMM_WORLD, &status);
			int received_count;
			MPI_Get_count(&status,MPI_INT,&received_count);
			FillExchangeBuffer(exchangeBuffer, influenced_info_array, received_count, range_begin);
			MPI_Send(exchangeBuffer, received_count, MPI_INT, i, 0, PETSC_COMM_WORLD);
		}
	}
	required_row_data->Print(rank);
	if(globalMaxExchangeBufferSize > 0)
		delete exchangeBuffer;
	return required_row_data;
}

void RowDataCollector::FillExchangeBuffer(PetscInt* exchangeBuffer, struct influenced_info* influenced_info_array,
		int received_count, PetscInt range_begin)
{
	for(int i = 0; i<received_count; i++)
	{
		struct influenced_info* influenced_info = &(influenced_info_array[exchangeBuffer[i] - range_begin]);
		if(influenced_info->row_info->set == CSET){
			exchangeBuffer[i] = influenced_info->row_number_in_coarse;
		}
		else{
			exchangeBuffer[i] = -1;
		}
	}
}

int RowDataCollector::GetExchangeBufferSize(int* row_data_request_sizes, int world_size)
{
	int max_request_size = 0;
	for(int i = 0; i<world_size; i++)
		max_request_size = std::max(max_request_size, row_data_request_sizes[i]);

	return max_request_size;
}

RowData* RowDataCollector::GetRequiredRowData(Mat mat, struct influenced_info* influenced_info_array)
{
	int size;
	MPI_Comm_size(PETSC_COMM_WORLD, &size);
	const PetscInt* ranges = new PetscInt[size+1];
	MatGetOwnershipRanges(mat,&ranges);
	RowData* row_data = new RowData(size, ranges);

	PetscInt first_row_in_range;
	PetscInt range_end;
	PetscErrorCode ierr;
	ierr = MatGetOwnershipRange(mat,&first_row_in_range,&range_end); CHKERRABORT(PETSC_COMM_WORLD, ierr);
	PetscInt ownership_range_size = range_end - first_row_in_range;

	const PetscInt* columns;
	PetscInt columns_number;
	filtered_row row;

	for(PetscInt i = 0; i<ownership_range_size; i++)
	{
		if(influenced_info_array[i].row_info->set == FSET)
		{
//			MatGetRow(mat,i + first_row_in_range,&columns_number,&columns,NULL);
			getFilteredRow(mat,i + first_row_in_range,&row);
			columns_number = row.columns_number;
			columns = row.columns;
			for(PetscInt j = 0; j<columns_number; j++)
			{
				if(columns[j] - first_row_in_range < 0 || columns[j] - first_row_in_range >= ownership_range_size)
				{
//					printf("%d is required my range: %d %d\n",columns[j], first_row_in_range, range_end);
					row_data->AddRowDataRequest(columns[j]);
				}
			}
			restoreFilteredRow(mat, i + first_row_in_range, &row);
//			MatRestoreRow(mat,i + first_row_in_range,&columns_number,&columns,NULL);
		}
	}

	return row_data;
}

RowDataCollector::RowDataCollector() {


}

RowDataCollector::~RowDataCollector() {

}

