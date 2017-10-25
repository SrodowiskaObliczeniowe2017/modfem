/*
 * AMGRows.hpp
 *
 *  Created on: Oct 18, 2015
 *      Author: damian
 */

#ifndef SRC_AMG_MKB_AMG_AMGROWS_HPP_
#define SRC_AMG_MKB_AMG_AMGROWS_HPP_
#include <petscsys.h>

enum Set
{
	DEFAULT = 0,
	CSET = 1,
	FSET = 2
};

struct row_info
{
	Set set;
	PetscInt local_row_number;

	row_info() : set(DEFAULT) {};
};

struct influenced_info
{
	int influenced_number;
	struct row_info* row_info;
	PetscInt row_number_in_coarse;

	influenced_info() : row_info(NULL), influenced_number(0), row_number_in_coarse(-1) {};
};


#endif /* SRC_AMG_MKB_AMG_AMGROWS_HPP_ */
