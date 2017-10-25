/*
 * RowData.hpp
 *
 *  Created on: Oct 17, 2015
 *      Author: damian
 */

#ifndef SRC_AMG_MKB_AMG_ROWDATA_HPP_
#define SRC_AMG_MKB_AMG_ROWDATA_HPP_

#include <map>
#include <petscsys.h>


class RowData {
public:
	RowData(int world_size,const PetscInt *ranges);
	virtual ~RowData();
	PetscInt GetRowData(PetscInt column_index);
	int GetRowDataOwner(PetscInt column_index);
	void AddRowDataRequest(PetscInt column_index);
	void GetRowDataRequestSizes(int* row_data_request_sizes);
	int FillExchangeRequest(PetscInt* exchangeBuffer, int rank);
	void FillExchangeResponse(PetscInt* exchangeBuffer, int rank);
	void Print(int rank);
private:
	std::map<PetscInt, PetscInt>** row_data;
	int world_size;
	const PetscInt *ranges;
};

#endif /* SRC_AMG_MKB_AMG_ROWDATA_HPP_ */
