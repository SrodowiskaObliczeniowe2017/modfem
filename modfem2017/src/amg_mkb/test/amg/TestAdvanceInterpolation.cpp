#define BOOST_TEST_MODULE AdvanceInterpolationTest
#include <boost/test/included/unit_test.hpp>
#include "../../amg/Interpolation/AdvancedRSCFInterpolation.hpp"
#include "../util/util.cpp"

FILE* utv_log_out = stdout;

double values[] =   {1000.0,100.0,10.0,2000.0,200.0,20.0,
					1001.0,101.0,11.0,2001.0,201.0,21.0,
					1002.0,102.0,12.0,2002.0,202.0,22.0,
					1003.0,103.0,13.0,2003.0,203.0,23.0,
					1004.0,104.0,14.0,2004.0,204.0,24.0,
					1005.0,105.0,15.0,2005.0,205.0,25.0};

int rows_numbers[] = {0,0,0,0,0,0,
				 1,1,1,1,1,1,
				 2,2,2,2,2,2,
				 3,3,3,3,3,3,
				 4,4,4,4,4,4,
				 5,5,5,5,5,5,
};

int columns_numbers[] = {0,1,2,3,4,5,
   			  0,1,2,3,4,5,
			  0,1,2,3,4,5,
			  0,1,2,3,4,5,
			  0,1,2,3,4,5,
			  0,1,2,3,4,5
};

double values2[] =   {2.0,-2.0,4.0,-4.0,
					  5.0,-5.0,3.0,-3.0,
					-7.0,7.0,1.0,-1.0,
					2.0,2.0,-6.0,-6.0};

int rows_numbers2[] = {0,0,0,0,
				 1,1,1,1,
				 2,2,2,2,
				 3,3,3,3};

int columns_numbers2[] = {0,1,2,3,
   			  0,1,2,3,
			  0,1,2,3,
			  0,1,2,3};


/*
BOOST_AUTO_TEST_CASE( AdvancedInterpolationTest1 )
{
	PetscInitialize(NULL,NULL,NULL,NULL);
	Mat matrix;
	prepare_matrix_general_crs(&matrix, values, rows_numbers, columns_numbers, 6, 36);
	struct row_info* row_info_array = new row_info[6];
	row_info_array[0].set = CSET;
	row_info_array[1].set = FSET;
	row_info_array[2].set = FSET;
	row_info_array[3].set = CSET;
	row_info_array[4].set = FSET;
	row_info_array[5].set = FSET;
	struct influenced_info* influened_info_array = new influenced_info[6];
	influened_info_array[0].row_number_in_coarse = 0;
	influened_info_array[3].row_number_in_coarse = 1;
	influened_info_array[0].row_info = row_info_array;
	influened_info_array[1].row_info = row_info_array+1;
	influened_info_array[2].row_info = row_info_array+2;
	influened_info_array[3].row_info = row_info_array+3;
	influened_info_array[4].row_info = row_info_array+4;
	influened_info_array[5].row_info = row_info_array+5;
	AdvancedRSCFInterpolation* interpolation = new AdvancedRSCFInterpolation(matrix, row_info_array,influened_info_array,
			0,6,0.5);

	Mat coarse_to_fine_mat = interpolation->GetMatrixFromCoarseToFine();
	MatView(coarse_to_fine_mat, PETSC_VIEWER_STDOUT_WORLD);

}*/

BOOST_AUTO_TEST_CASE( AdvancedInterpolationTest2 )
{
	PetscInitialize(NULL,NULL,NULL,NULL);
	Mat matrix;
	prepare_matrix_general_crs(&matrix, values2, rows_numbers2, columns_numbers2, 4, 16);
	struct row_info* row_info_array = new row_info[4];
	row_info_array[0].set = FSET;
	row_info_array[1].set = CSET;
	row_info_array[2].set = CSET;
	row_info_array[3].set = FSET;
	struct influenced_info* influened_info_array = new influenced_info[4];
	influened_info_array[1].row_number_in_coarse = 0;
	influened_info_array[2].row_number_in_coarse = 1;
	influened_info_array[0].row_info = row_info_array;
	influened_info_array[1].row_info = row_info_array+1;
	influened_info_array[2].row_info = row_info_array+2;
	influened_info_array[3].row_info = row_info_array+3;
	AdvancedRSCFInterpolation* interpolation = new AdvancedRSCFInterpolation(matrix, row_info_array,influened_info_array,
			0,4,0.5);

	Mat coarse_to_fine_mat = interpolation->GetMatrixFromCoarseToFine();
	MatView(coarse_to_fine_mat, PETSC_VIEWER_STDOUT_WORLD);

}

