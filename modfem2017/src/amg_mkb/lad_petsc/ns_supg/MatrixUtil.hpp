#ifndef NS_SUPG_MATRIX_UTIL_H
#define NS_SUPG_MATRIX_UTIL_H

#include <petscksp.h>

Mat create_exact_schur_complement(Mat* schur_complement, Mat Avp, Mat Avv, Mat Apv, Mat App);
Mat create_velocity_exact_inverse(Mat* velocity_inverse, Mat Avv, Vec bvp);
Mat create_velocity_ilu_inverse(KSP* preonly_solver, Mat Avv);
Mat create_velocity_mass_inverse(Mat Vmm, Vec bvp);
void draw_matrix(Mat mat, double threshold);
Mat restrict_matrix(Mat mat, double threshold);
void print_sym_difference(Mat mat);
void print_positive_to_negative(Mat mat);


#endif
