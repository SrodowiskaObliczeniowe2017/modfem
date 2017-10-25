/*
 * mic.h
 *
 *  Created on: 16 lis 2015
 *      Author: pmaciol
 */

#ifndef MOD_FEM_VIEWER_FEMVIEWER_MIC_MIC_H_
#define MOD_FEM_VIEWER_FEMVIEWER_MIC_MIC_H_


extern int init_platorm(int argc,char** argv);

enum { mic_CPU, mic_XEONPHI };

typedef enum { mic_NoSegment, mic_DensityZero, mic_ColorHavebreakpoints } vol_integr_t;
typedef enum { mic_Riemann, mic_Trapezoid, mic_Full } vol_rend_integr_t;

#endif /* MOD_FEM_VIEWER_FEMVIEWER_MIC_MIC_H_ */
