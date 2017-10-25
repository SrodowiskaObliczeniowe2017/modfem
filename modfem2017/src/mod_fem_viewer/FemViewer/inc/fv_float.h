/*
 * fv_float.h
 *
 *  Created on: 22 maj 2014
 *      Author: dwg
 */

#ifndef FV_FLOAT_H_
#define FV_FLOAT_H_

#include"MathHelper.h"

#define USE_DB_PREC 1





#if USE_DB_PREC == 0
typedef float mfvFloat_t;
typedef float ScalarValueType;
typedef FemViewer::fvmath::CVec3f ScalarValueType3;
    #define FLOAT_CONVERT(p) p##f

#else
typedef double mfvFloat_t;
typedef double ScalarValueType;
typedef FemViewer::fvmath::CVec3d ScalarValueType3;
    #define FLOAT_CONVERT(p) p
#endif

#define COORD_TYPE_FLOAT 1
#if COORD_TYPE_FLOAT == 1
typedef float CoordType;
#else
typedef mfvFloat_t CoordType;
#endif

#define FV_FLOAT_MAX FLOAT_CONVERT(3.2e37)
#define FV_FLOAT_COMPARE FLOAT_CONVERT(3.1e37)
#endif /* FV_FLOAT_H_ */
