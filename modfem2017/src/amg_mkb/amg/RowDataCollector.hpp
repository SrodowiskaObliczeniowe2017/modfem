/*
 * RowDataCollector.hpp
 *
 *  Created on: Oct 17, 2015
 *      Author: damian
 */

#ifndef SRC_AMG_MKB_AMG_ROWDATACOLLECTOR_HPP_
#define SRC_AMG_MKB_AMG_ROWDATACOLLECTOR_HPP_
#include "AMGRows.hpp"
#include <map>
#include <petscmat.h>
#include <algorithm>
#include <mpi.h>
#include "RowData.hpp"
#include "AMGMatrixUtilityFunctions.hpp"




class RowDataCollector {
public:
	RowData* Collect(Mat mat, struct influenced_info* influenced_info_array);
	RowDataCollector();
	virtual ~RowDataCollector();
private:
	RowData* GetRequiredRowData(Mat mat, struct influenced_info* influenced_info_array);
	int GetExchangeBufferSize(int* row_data_request_sizes, int world_size);
	void FillExchangeBuffer(PetscInt* exchangeBuffer, struct influenced_info* influenced_info_array, int received_count, PetscInt range_begin);
};

#endif /* SRC_AMG_MKB_AMG_ROWDATACOLLECTOR_HPP_ */
