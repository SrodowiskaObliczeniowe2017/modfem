/*
 * JacomiMethodTest.cpp
 *
 *  Created on: Jun 5, 2015
 *      Author: damian
 */

#define BOOST_TEST_MODULE RSCFSplitterTest
#include <boost/test/included/unit_test.hpp>
#include "../../amg/RSCFSplitter.hpp"
#include "petscmat.h"

FILE* utv_log_out = stdout;

int row_numbers[] = {0,0,0,0,0,0,0,1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,5,5,6,6,6,6,6,6,6};
int column_numbers[] = {0,1,2,3,4,5,6,0,1,2,3,4,5,6,0,1,2,3,4,5,6,0,1,2,3,4,5,6,0,1,2,3,4,5,6,0,1,2,3,4,5,6,0,1,2,3,4,5,6};
double values[] = {-2.0,6.0,-6.0,8.0,-10.0,12.0,-10.0,
					6.0,-10.0,12.0,-10.0,8.0,-6.0,6.0,
					-6.0,12.0,6.0,-6.0,8.0,-10.0,12.0,
					8.0,-10.0,-6.0,-10.0,8.0,-6.0,6.0,
					-10.0,8.0,8.0,8.0,6.0,-6.0,8.0,
					12.0,-6.0,-10.0,-6.0,-6.0,8.0,-6.0,
					-10.0,6.0,12.0,6.0,8.0,-6.0,-6.0};

double values2[] = {2.0,-2.0,6.0,-8.0,-10.0,-12.0,-14.0,
					-2.0,4.0,-2.0,-4.0,-6.0,-8.0,-10.0,
					-6.0,-2.0,6.0,-2.0,-4.0,-6.0,-8.0,
					-8.0,-4.0,-2.0,8.0,-6.0,-8.0,-10.0,
					-10.0,-6.0,-4.0,-6.0,10.0,-12.0,-14.0,
					-12.0,-8.0,-6.0,-8.0,-12.0,12.0,-6.0,
					-14.0,-10.0,-8.0,-10.0,-14.0,-6.0,14.0};

double interpolation_matrix[] = {9.2307692308, 10.7692307692,
								3.5555555556, 4.4444444444,
								2, 2.6666666667,
								2.1111111111, 2.6388888889,
								2.4, 2.8,
								1.0, 0.0,
								0.0, 1.0};


void prepare_matrix(Mat* matrix, double* values)
{
	int size = 7;
	MatCreate(PETSC_COMM_WORLD, matrix);
	MatSetSizes(*matrix, PETSC_DECIDE, PETSC_DECIDE, size, size);
	MatSetFromOptions(*matrix);
	MatMPIAIJSetPreallocation(*matrix,30,NULL,30,NULL);
	MatSeqAIJSetPreallocation(*matrix,30,NULL);
	MatSetUp(*matrix);

	int i;
	for(i = 0; i<49; i++)
		MatSetValues(*matrix, 1, &row_numbers[i], 1, &column_numbers[i], &values[i], INSERT_VALUES);

	MatAssemblyBegin(*matrix, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*matrix, MAT_FINAL_ASSEMBLY);
	MatView(*matrix, PETSC_VIEWER_STDOUT_WORLD);
}

BOOST_AUTO_TEST_CASE( RSCFSplitterTestSplitting )
{
	PetscInitialize(NULL,NULL,NULL,NULL);
	Mat matrix;
	prepare_matrix(&matrix, values);

	RSCFSplitter* splitter = new RSCFSplitter(matrix,0.5);
	splitter->MakeCFSplitting();

	std::map<int,int> c_map;
	std::map<int,int> f_map;
	splitter->GetSetsWithInfluencedNumebr(&c_map,&f_map);

	for (std::map<int,int>::iterator it=c_map.begin(); it!=c_map.end(); ++it)
	    printf("%d belongs to c with influenced nr %d\n", it->first, it->second);

	BOOST_CHECK(c_map[0] == 6);
	BOOST_CHECK(c_map[5] == 5);

	printf("\n");
	for (std::map<int,int>::iterator it=f_map.begin(); it!=f_map.end(); ++it)
	    printf("%d belongs to f with influenced nr %d\n", it->first, it->second);

	BOOST_CHECK(f_map[1] == 2);
	BOOST_CHECK(f_map[2] == 3);
	BOOST_CHECK(f_map[3] == 3);
	BOOST_CHECK(f_map[4] == 2);
	BOOST_CHECK(f_map[6] == 2);

}

BOOST_AUTO_TEST_CASE( RSCFSplitterTestCoarse )
{
	PetscInitialize(NULL,NULL,NULL,NULL);
	Mat matrix;
	prepare_matrix(&matrix, values2);

	RSCFSplitter* splitter = new RSCFSplitter(matrix,0.5);
	splitter->MakeCFSplitting();
	Mat mat = splitter->GetMatrixFromCoarseToFine();
	MatView(mat,PETSC_VIEWER_STDOUT_WORLD);

	double values[14];
	for(int i = 0; i<7; i++)
		for(int j = 0; j<2; j++)
		{
			MatGetValues(mat,1,&i,1,&j,values+i*2+j);
			BOOST_CHECK_CLOSE(values[i*2+j], interpolation_matrix[i*2+j], 0.001);
		}
}
