#define BOOST_TEST_MODULE JacobiMethodTest
#include <boost/test/included/unit_test.hpp>
#include "../../lad_petsc/AMGSolverStructure.hpp"
#include "../util/util.cpp"
#include "../../lad_petsc/KnownSolutionErrorEvaluator.cpp"

FILE* utv_log_out = stdout;

double values[] =   {1.0,0.0,0.0,0.0,0.0,0.0,0.0,
					-1.0,2.0,-1.0,0.0,0.0,0.0,0.0,
					0.0,-1.0,2.0,-1.0,0.0,0.0,0.0,
					0.0,0.0,-1.0,2.0,-1.0,0.0,0.0,
					0.0,0.0,0.0,-1.0,2.0,-1.0,0.0,
					0.0,0.0,0.0,0.0,-1.0,2.0,-1.0,
					0.0,0.0,0.0,0.0,0.0,0.0,1.0};

double vec[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};

BOOST_AUTO_TEST_CASE( AMGSolverStructureTest )
{
    clock_t begin = clock();

	int size = 100;
	PetscInitialize(NULL,NULL,NULL,NULL);
	Mat matrix;
	Vec rhs;
	prepare_testcase_1_matrix(&matrix, size);
	prepare_testcase_1_rhs(&rhs, size);

	VecAssemblyBegin(rhs);
	VecAssemblyEnd(rhs);

	ErrorEvaluator* error_evaluator = new KnownSolutionErrorEvaluator(size);
	AMGSolverStructure* solver = new AMGSolverStructure(matrix, rhs, size, error_evaluator, false);

	solver->CreateAMGSolverLevels(0,0,0.5,10);
	solver->RunVCycle();

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	printf("TOTAL TIME: %lf\n",elapsed_secs);
	delete solver;
}

BOOST_AUTO_TEST_CASE( AMGSolverStructureChangeSystemTest )
{
    clock_t begin = clock();

	int size = 10;
	PetscInitialize(NULL,NULL,NULL,NULL);
	Mat matrix;
	Vec rhs;
	prepare_testcase_1_matrix(&matrix, size);
	prepare_testcase_1_rhs(&rhs, size);

	VecAssemblyBegin(rhs);
	VecAssemblyEnd(rhs);

	ErrorEvaluator* error_evaluator = new ResidualBasedErrorEvaluator(matrix, rhs, 1.0e-7, 1.0e-8, size, -1);
	AMGSolverStructure* solver = new AMGSolverStructure(matrix, rhs, size, error_evaluator, false);

	solver->CreateAMGSolverLevels(0,0,0.5,3);
	Vec result = solver->RunVCycle();

	for(int i = 0; i<10; i++)
	{
		double value;
		VecGetValues(result,1,&i,&value);
		BOOST_CHECK_CLOSE( value, 1.0/9.0 * i, 0.0001 );
	}

	Vec rhs2;
	VecDuplicate(rhs, &rhs2);

	for(int i = 0; i<5; i++)
	{
		double value = i + 1;
		VecSetValues(rhs2,1,&i,&value,INSERT_VALUES);
		int index = 9 - i;
		VecSetValues(rhs2,1,&index,&value,INSERT_VALUES);
	}

	solver->ChangeSystem(rhs2,-1);
	result = solver->RunVCycle();

	double solution[] = { -1.000, -15.000, -27.000, -36.000, -41.000, -41.000, -36.000, -27.000, -15.000, -1.000};
	for(int i = 0; i<10; i++)
	{
		double value;
		VecGetValues(result,1,&i,&value);
		BOOST_CHECK_CLOSE( value, solution[i], 0.0001 );
	}

	Vec rhs3;
	VecDuplicate(rhs, &rhs3);
	VecCopy(rhs, rhs3);

	solver->ChangeSystem(rhs3,-1);
	result = solver->RunVCycle();

	for(int i = 0; i<10; i++)
	{
		double value;
		VecGetValues(result,1,&i,&value);
		BOOST_CHECK_CLOSE( value, -1.0/9.0 * i, 0.0001 );
	}

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	printf("TOTAL TIME: %lf\n",elapsed_secs);
	delete solver;
}
