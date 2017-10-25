/*
 * RowData.cpp
 *
 *  Created on: Oct 17, 2015
 *      Author: damian
 */

#include "RowData.hpp"

RowData::RowData(int world_size,const PetscInt *ranges) {
	//for(int i = 0; i<world_size+1; i++)
	//	printf("ranges %d %d\n", ranges[i], world_size);
	this->world_size = world_size;
	this->ranges = ranges;
	this->row_data = new std::map<PetscInt, PetscInt>*[world_size];
	for(int i = 0; i<world_size; i++)
	{
		row_data[i] = new std::map<PetscInt, PetscInt>();
	}
}

void RowData::Print(int rank)
{
//	char buf[1000]{0};
//	//sprintf(buf,"from %d\n",rank);
//	for(int i = 0; i<world_size; i++)
//	{
//		for (std::map<PetscInt,PetscInt>::iterator it=row_data[i]->begin(); it!=row_data[i]->end(); ++it)
//		{
//			sprintf(buf,"%s %d -> %d",buf,it->first,it->second);
//		}
//		sprintf(buf,"%s\n",buf);
//	}
//	printf("%s",buf);
}

void RowData::AddRowDataRequest(PetscInt column_index){
	for(int i = 0; i<world_size; i+=1){
		if(ranges[i] <= column_index && ranges[i+1] > column_index){
			(*row_data[i])[column_index] = -1;
		}
	}
}

PetscInt RowData::GetRowData(PetscInt column_index){
	for(int i = 0; i<world_size; i+=1){
		if(ranges[i] <= column_index && ranges[i+1] > column_index){
			return (*row_data[i])[column_index];
		}
	}
}

int RowData::GetRowDataOwner(PetscInt column_index){
	for(int i = 0; i<world_size; i+=1){
		if(ranges[i] <= column_index && ranges[i+1] > column_index){
			return i;
		}
	}
}

void RowData::GetRowDataRequestSizes(int* row_data_request_sizes){
	for(int i = 0; i<world_size; i++){
		row_data_request_sizes[i] = row_data[i]->size();
	}
}

int RowData::FillExchangeRequest(PetscInt* exchangeBuffer, int rank)
{
	int index = 0;
	for (std::map<PetscInt,PetscInt>::iterator it=row_data[rank]->begin(); it!=row_data[rank]->end(); ++it)
	{
		exchangeBuffer[index++] = it->first;
	}
	return row_data[rank]->size();
}

void RowData::FillExchangeResponse(PetscInt* exchangeBuffer, int rank)
{
	//size not needed
	int index = 0;
	for (std::map<PetscInt,PetscInt>::iterator it=row_data[rank]->begin(); it!=row_data[rank]->end(); ++it)
	{
		it->second = exchangeBuffer[index++];
	}
}

RowData::~RowData() {
	for(int i = 0; i<world_size; i++){
		delete row_data[i];
	}
	delete[] row_data;
	delete[] ranges;
}

