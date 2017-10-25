#define BOOST_TEST_MODULE ApproximateIverseTest

#include "../../lad_petsc/ns_supg/ApproximateInverse.hpp"
#include "../../lad_petsc/ns_supg/ApproximateInverseOpt.hpp"
#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include "util.cpp"

FILE* utv_log_out = stdout;

BOOST_AUTO_TEST_CASE( ApproximateInverseTest1 )
{
	PetscInitialize(NULL,NULL,NULL,NULL);
	Mat example_matrix, pattern_matrix;
	prepare_matrix_inverse_example_1(&example_matrix);
	prepare_matrix_inverse_full_pattern(&pattern_matrix);
	Mat result = getApproximateInverse2(example_matrix, 8, pattern_matrix);
	MatView(result, PETSC_VIEWER_STDOUT_WORLD);

//	row 0: (0, -0.331321)  (1, -0.154259)  (2, -0.171026)  (3, -0.114688)  (4, -0.0838364)  (5, -0.0496311)
//	row 1: (0, -0.154259)  (1, -0.395708)  (2, -0.221328)  (3, -0.207243)  (4, -0.128102)  (5, -0.0838364)
//	row 2: (0, -0.171026)  (1, -0.221328)  (2, -0.462777)  (3, -0.251509)  (4, -0.207243)  (5, -0.114688)
//	row 3: (0, -0.114688)  (1, -0.207243)  (2, -0.251509)  (3, -0.462777)  (4, -0.221328)  (5, -0.171026)
//	row 4: (0, -0.0838364)  (1, -0.128102)  (2, -0.207243)  (3, -0.221328)  (4, -0.395708)  (5, -0.154259)
//	row 5: (0, -0.0496311)  (1, -0.0838364)  (2, -0.114688)  (3, -0.171026)  (4, -0.154259)  (5, -0.331321)
}

//BOOST_AUTO_TEST_CASE( ApproximateInverseCommutatorTest1 )
//{
//	PetscInitialize(NULL,NULL,NULL,NULL);
//	Mat Avv, Avp;
//	prepare_matrix_inverse_example_2(&Avv);
//	prepare_matrix_inverse_example_3(&Avp);
//	Mat result = getApproximateInverseCommutator(Avv, Avp, 2);
//	MatView(result, PETSC_VIEWER_STDOUT_WORLD);
//
////	row 0: (0, 14.5455)  (1, 6.76486)
////	row 1: (0, -4.95455)  (1, 3.8514)
//}

BOOST_AUTO_TEST_CASE( ApproximateInverseTest3 )
{
	PetscInitialize(NULL,NULL,NULL,NULL);
	Mat example_matrix;
	prepare_matrix_inverse_example_4(&example_matrix);
	Mat result = getApproximateInverse2(example_matrix, 10, NULL);
	MatView(result, PETSC_VIEWER_STDOUT_WORLD);

//	-0.050  0.069 -0.198  0.085
//	 0.064  0.057  0.139 -0.092
//	 0.293 -0.180  0.096  0.080
//	-0.205  0.107 -0.120  0.109
}


BOOST_AUTO_TEST_CASE( ApproximateInverseTest4 )
{
	PetscInitialize(NULL,NULL,NULL,NULL);
	Mat example_matrix;
	prepare_matrix_inverse_example_5(&example_matrix);
	Mat result = getApproximateInverseOpt(example_matrix, 5);
	int i,j;
	double value;
	i = 0; j = 0;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - 0.364, 0.001 );
	i = 0; j = 1;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - -0.091, 0.001 );
	i = 0; j = 2;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - 0.273, 0.001 );
	i = 0; j = 3;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - 0.091, 0.001 );
	i = 0; j = 4;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - -0.545, 0.001 );
	i = 1; j = 0;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - 1.727, 0.001 );
	i = 1; j = 1;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - -0.182, 0.001 );
	i = 1; j = 2;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - 0.545, 0.001 );
	i = 1; j = 3;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - 0.182, 0.001 );
	i = 1; j = 4;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - -2.091, 0.001 );
	i = 2; j = 0;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - -5.642, 0.001 );
	i = 2; j = 1;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - 1.051, 0.001 );
	i = 2; j = 2;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - -2.216, 0.001 );
	i = 2; j = 3;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - -0.614, 0.001 );
	i = 2; j = 4;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - 6.932, 0.001 );
	i = 3; j = 0;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - 1.483, 0.001 );
	i = 3; j = 1;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - -0.074, 0.001 );
	i = 3; j = 2;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - 0.534, 0.001 );
	i = 3; j = 3;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - -0.114, 0.001 );
	i = 3; j = 4;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - -1.568, 0.001 );
	i = 4; j = 0;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - -1.040, 0.001 );
	i = 4; j = 1;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - -0.006, 0.001 );
	i = 4; j = 2;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - -0.420, 0.001 );
	i = 4; j = 3;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - 0.068, 0.001 );
	i = 4; j = 4;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_CHECK_SMALL( value - 1.341, 0.001 );
//	MatView(result, PETSC_VIEWER_STDOUT_WORLD);

//	 0.364 -0.091  0.273  0.091 -0.545
//	 1.727 -0.182  0.545  0.182 -2.091
//	-5.642  1.051 -2.216 -0.614  6.932
//	 1.483 -0.074  0.534 -0.114 -1.568
//	-1.040 -0.006 -0.420  0.068  1.341
}


BOOST_AUTO_TEST_CASE( ApproximateInverseTest5 )
{
	PetscInitialize(NULL,NULL,NULL,NULL);
	Mat example_matrix;
	prepare_matrix_inverse_example_5(&example_matrix);
	Mat result = getApproximateInverseOpt(example_matrix, 3);
	int i,j;
	double value;
	i = 0; j = 0;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_REQUIRE_CLOSE( value, -0.0327702, 0.01 );
	i = 0; j = 1;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_REQUIRE_CLOSE( value, -0.101567, 0.01 );
	i = 1; j = 0;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_REQUIRE_CLOSE( value, 0.0111948, 0.01 );
	i = 1; j = 1;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_REQUIRE_CLOSE( value, -0.151639, 0.01 );
	i = 1; j = 2;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_REQUIRE_CLOSE( value, -0.217164, 0.01 );
	i = 2; j = 0;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_REQUIRE_CLOSE( value, 0.147161, 0.01 );
	i = 2; j = 1;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_REQUIRE_CLOSE( value, 0.915734, 0.01 );
	i = 2; j = 2;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_REQUIRE_CLOSE( value, 0.585097, 0.01 );
	i = 2; j = 3;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_REQUIRE_CLOSE( value, 0.125, 0.01 );
	i = 2; j = 4;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_REQUIRE_CLOSE( value, -0.09375, 0.01 );
	i = 3; j = 2;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_REQUIRE_CLOSE( value, 0.0171064, 0.01 );
	i = 3; j = 3;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_REQUIRE_CLOSE( value, -0.240385, 0.01 );
	i = 3; j = 4;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_REQUIRE_CLOSE( value, 0.0841346, 0.01 );
	i = 4; j = 3;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_REQUIRE_CLOSE( value, 0.182692, 0.01 );
	i = 4; j = 4;
	MatGetValues(result,1,&i,1,&j,&value); BOOST_REQUIRE_CLOSE( value, 0.0360577, 0.01 );

//	MatView(result, PETSC_VIEWER_STDOUT_WORLD);
//	row 0: (0, -0.0327702)  (1, -0.101567)
//	row 1: (0, 0.0111948)  (1, -0.151639)  (2, -0.217164)
//	row 2: (0, 0.147161)  (1, 0.915734)  (2, 0.585097)  (3, 0.125)  (4, -0.09375)
//	row 3: (2, 0.0171064)  (3, -0.240385)  (4, 0.0841346)
//	row 4: (3, 0.182692)  (4, 0.0360577)
}
