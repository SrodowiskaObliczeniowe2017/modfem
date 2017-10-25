#ifndef APP_INVERSE_H
#define APP_INVERSE_H

#include <petscksp.h>
#include "MatrixUtil.hpp"
#include <time.h>
#include <math.h>
#include "../../../include/uth_system.h"
#include <map>

Mat getApproximateInverse(Mat Avv);
Mat getApproximateInverse2(Mat Avv, int band, Mat restricted_matrix);
Mat getApproximateInverseCommutator(Mat Avv, Mat Avp, int band);

#endif
