/*
 * ocl.h
 *
 *  Created on: 15 cze 2014
 *      Author: dwg
 */

#ifndef OCL_H_
#define OCL_H_


#include<CL/cl.h>
#include<stdio.h>
#include<string>
#include<vector>

#ifdef __cplusplus
extern "C" {
#endif

/* A structure for information about device type
 * and it's count
 */
typedef struct {
	cl_device_type type; /* A type of device */
	cl_uint   num;/* The number of given type */
	cl_int	  gl_index; /* The index of device with GL support. otherwise -1 */
} device_info_t;

/* Info about process it's hoat name and devices */
typedef struct __host_info_id {

	//cl_int gl_support; /* Index of the device supported GL, otherwise -1*/
	std::string platform; /* A name of OpenCL platform */
	std::vector<device_info_t> devices;

	void clear() {
		//gl_support = -1;
		platform.clear();
		devices.clear();
	}
} host_info, * host_info_id;



/* stuff for raytracing */



/* Initialize OpenCL. It returns the number of awailabe
 * OpnecL devices for 1st supported platform.
 * -1 means that no OpenCL support or some erros appears
 **/
extern cl_int initOpenCL(host_info_id host,cl_device_type devType,cl_bool supportGL);


extern cl_int loadDataCL();
extern cl_int displayCL();
extern cl_int runCL();
extern cl_int loadField(void *Ptr);
extern cl_int bindOpenCL(cl_uint GLpbo,cl_int w,cl_int h);
extern cl_int setDevice();

/* Get information of all process */

/* Shut down OpenCL */
extern void shutdownOpenCL(void);

#ifdef __cplusplus
}
#endif

namespace FemViewer {


} // end namespace


#endif /* OCL_H_ */
