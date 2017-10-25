#ifndef NS_SUPG_VECTOR_UTIL_H
#define NS_SUPG_VECTOR_UTIL_H

#include <petscvec.h>

void convert_vec_to_x(Vec solution_v, Vec solution_p, double* x, int Ndof);
void convert_x_to_vec(Vec solution_v, Vec solution_p, double* x, int Ndof);
void convert_vec_to_x(Vec solution, double* x, int Ndof);
void convert_x_to_vec(Vec solution, double* x, int Ndof);
void convert_vec_to_x_simple(Vec solution, double* x, int Ndof);
void convert_x_to_vec_simple(Vec solution, double* x, int Ndof);


#endif
