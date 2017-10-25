#ifndef APP_INVERSE_OPT_H
#define APP_INVERSE_OPT_H

#include <petscksp.h>
#include "MatrixUtil.hpp"
#include <time.h>
#include <math.h>
#include "../../../include/uth_system.h"
#include <map>
#include "lin_alg_intf.h"

Mat getApproximateInverseOpt(Mat Avv, int band);

#endif
