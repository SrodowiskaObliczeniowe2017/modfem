/*
 * types.h
 *
 *  Created on: 10 lut 2015
 *      Author: pmaciol
 */

#ifndef TYPES_H_
#define TYPES_H_

#ifndef __OPENCL_VERSION__
#include <CL/cl.h>
#include <CL/cl_platform.h>
typedef cl_float3 float3;
typedef cl_float2 float2;
typedef cl_int3 int3;
typedef cl_int2 int2;

#endif

typedef struct _gridinfo_t {
	double min[3];
	double max[3];
	double cellSize[3];
	int   resolution[3];
 } grid_t;

 struct gridaccel_s
 {
 #ifdef __OPENCL_VERSION__
  __constant gridinfo_t * gridinfo;
  __constant int * grid;
  __constant int * cells;
 #else
  cl_mem gridinfo;
  cl_mem grid;
  cl_mem cells;
 #endif
 };

struct frameconfig_s
{
	float3 eyepoint;
	float3 lookat;
	float3 upvector;
	int2   angles;
};

struct pixel_s
{
	float3 color;
	float fill;
	unsigned int lock;
};

struct reconstruct_s
{
	int pixelNumber;
	float importance;
	unsigned int randomSeed;
};

struct ray_s
{
	float3 origin;
	float3 direction;
	struct reconstruct_s reconstructInfo;
};

struct elem_info {
	float vts[3*6];
	float pls[4*5];
	float min_v;
	float max_v;
};

typedef elem_info elem_info_t;

struct coeffs_info {
	short pdeg;
	short base;
	int   nr_shap;
	float coeffs[128];
};

typedef coeffs_info coeffs_info_t;

/*
 * This structure is passed once between the host and kernel for
 * data initialized before the run.
 */
typedef struct raycaster_context {
  cl_float4 dir_top_start;
  cl_float4 dir_bottom_start;
  cl_float4 dir_bottom_stop;
  cl_float4 eyeP;
  cl_float4 lightP;
  cl_int2   window_size;
  cl_float  epsilon;
  cl_int maxIterations;
  cl_int stride;
  cl_int pad[3];
} rcc_t;

typedef struct {
  /*  camera parameters */
  float fov;
  float aspect;
  float translate[2];
  float zoom;
  float curRotation[16];

  /*  rendering parameters */
  int shadows; /*  whether self-shadowing is being rendered */
} camera_t;








#endif /* TYPES_H_ */
