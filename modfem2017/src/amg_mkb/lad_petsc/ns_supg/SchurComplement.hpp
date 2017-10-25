#ifndef NS_SUPG_SCHUR_H
#define NS_SUPG_SCHUR_H

#include <petscksp.h>
#include "ApproximateInverseOpt.hpp"

void create_s_complement(Mat Avv, Mat Avp, Mat Apv, Mat App, Mat* S, Vec bvp);
void create_s_complement_Vmm(Mat Avv, Mat Avp, Mat Apv, Mat App, Mat* S, Mat Vmm, Vec bvp);
Mat create_s_complement_spai(Mat Avv, Mat Avp, Mat Apv, Mat App, Mat restricted_matrix, Mat* S);
void create_s_complement_block(Mat Avv, Mat Avp, Mat Apv, Mat App, Mat* S);

#endif
