/*
 * ocl.cpp
 *
 *  Created on: 15 cze 2014
 *      Author: Paweł Macioł
 */
#include "Log.h"
#include "types.h"
#include "fv_threads.h"
#include "fv_helper_routines.h"
#include "fv_dictstr.h"
#include "fv_config.h"
#include "ocl.h"
#include "Mesh.h"
#include "VtxAccumulator.h"
#include "Field.h"
#include "ModelControler.h"
#include "ViewManager.h"
#include <memory>
#include <vector>
#include <limits>
#include <string>
#include <fstream>
#include <iostream>
// OpenCL C++ wraper
//#define CL_USE_DEPRECATED_OPENCL_1_1_APIS
//#undef CL_VERSION_1_2
#define  __CL_ENABLE_EXCEPTIONS
// This header is needed becauseof 1.1 OpenCL by Nvidia
#include "cl.hpp"
#include <CL/cl_gl.h>
#ifndef WIN32
#include <GL/glx.h>
#endif

#define LOCAL_SIZE_X 16
#define LOCAL_SIZE_Y 16
#define GL_SHARING_EXTENSION "cl_khr_gl_sharing"

cl_int error;
cl_int pixcount;
// Grid of threads sizes
size_t gwsize[2];
size_t lwsize[2] = {16,16};
// Pointer to array of pixel structures
pixel_s * image = NULL;

static rcc_t rctx;
static float N[3];                     /* Normal to the image plane */
static float T[3];                     /* Tangent of the image plane */
static float B[3];                     /* Binormal of the image plane */
static camera_t rc;           		   /* Misc screen parameters */

static int first_frame   = 1;
static int current_frame = 0;
static int image_size;
static unsigned char *mapped_addr[2]; // Buffers for frame images

// Timer
static fv_timer ocl_timer;

extern int gridcount, cellcount;
extern int img_width;
extern int img_height;



using namespace FemViewer;
/* Paths to sources of kernels */
const char* kernels[] = {
		"kernel_grid.cl"
		/*"/home/dwg/development/mod_fem_code/src/mod_fem_viewer/FemViewer/ocl/*/"kernel_zero.cl"
		/*"/home/dwg/development/mod_fem_code/src/mod_fem_viewer/FemViewer/ocl/*/"kernels.cl"

};

/* Registered types od compute devices */
static const
cl_device_type _device_types[] = {
		CL_DEVICE_TYPE_CPU,
		CL_DEVICE_TYPE_GPU,
		CL_DEVICE_TYPE_ACCELERATOR
};

static void calculateView(RContext* mrc, camera_t * cam);
static void initWorkload(RContext* mrc);
//static void *next_frame(size_t **size);
//static void *last_frame(int *size);


size_t roundUp(int group_size, int global_size)
{
    int r = global_size % group_size;
    if(r == 0)
    {
        return global_size;
    } else
    {
        return global_size + group_size - r;
    }
}

typedef cl::Platform 		oclPlatform;
typedef cl::Device   		oclDevice;
typedef cl::Context  		oclContext;
typedef cl::CommandQueue	oclCmdQ;
typedef cl::Program		 	oclProgram;
typedef cl::Kernel	 		oclKernel;
typedef cl::Event		 	oclEvent;
typedef cl::Error			oclError;
typedef cl::BufferGL		oclBufferGL;
typedef cl::Buffer			oclBuffer;

// Bufers
oclBuffer outBuffers[2];
oclBuffer elDataBuffer;
oclBuffer coeffBuffer;
oclBuffer gridInfoBuffer;
oclBuffer gridBuffer;
oclBuffer cellsBuffer;
// Events
oclEvent mapEvt, unmapEvt;
oclEvent evt;
std::vector<oclEvent> events;
// Kernels
oclKernel kern_priray;

//onst char * extGLstring = "cl_khr_sharing";



class ocl
{
public:
	enum eErr {
		TEST_OK		= CL_SUCCESS,
		TEST_FAILED = -100,
	};

	static StrDict errMsgs[];
	static void    printPlatformInfo();
	static void    printDeviceInfo();
	// External handle functions
	friend cl_int  initOpenCL(host_info_id host,cl_device_type devType,cl_bool supportGL);
	friend cl_int  loadDataCL();
	friend cl_int  displayCL();
	friend cl_int  bindOpenCL(cl_uint GLpbo, cl_int w, cl_int h);
	friend void    shutdownOpenCL(void);

private:
	// Internal containers for platforms and devices
	typedef std::vector<cl::Platform> vPlatforms;
	typedef std::vector<cl::Device>	  vDevices;
	static vPlatforms _platforms;
	static vDevices   _devices;
	// Self pointer
	static ocl * _selfp;
	static Lock _lock;

	//
	static void   lock() { _lock.lock(); }
	static void unlock() { _lock.unlock(); }

	static cl_int invokeTestKernel(void);
	static cl_int integrateWithOpenGL(void);
  public:

	static ocl*& instance() { return _selfp; }
	// Dtr
	~ocl();
	void clean();

	// Returns the number of supported devices for current platform
	static cl_uint countDevices(const oclPlatform& pl,const cl_device_type devType);
	static cl_uint isDeviceSupported(cl_device_type devType);
	static cl_int  isExtensionSupported(const char * extString,const cl_device_type type);

	cl_int  loadProgram(const std::string& srcPath,bool quiet = true);
	static void* next_frame(size_t **size);
	static void* last_frame(int *size);
 protected:
	cl_int 		_currPlatform; // Index of platform in use
	cl_int 		_currDevice;   // index of device in use
	oclContext	_context;
	oclCmdQ		_queue;
	oclProgram	_program;
	oclKernel	_kernel;
	oclEvent 	_event;
	cl_int	  	_err;
	oclBufferGL _bufferGL;
	cl_uint		_gridsize[2];

 private:
	// Ctr
	explicit ocl(cl_uint devType = CL_DEVICE_TYPE_DEFAULT,
			cl_bool = CL_FALSE, cl_bool = CL_FALSE);
	oclContext createContext(cl_uint devType, cl_bool supportGL = false);

	// Not implemented
	ocl(const ocl&);
	ocl& operator=(const ocl&);
};


StrDict ocl::errMsgs[] = {
		{TEST_FAILED, "Test failed at startup!\n"}
};

std::vector<cl::Platform> ocl::_platforms;
std::vector<cl::Device>   ocl::_devices;
ocl * ocl::_selfp = nullptr;
Lock ocl::_lock;



/*---------------------------------------------------------------------
 * Invoke kernel dummy for start working
 *--------------------------------------------------------------------*/
cl_int ocl::invokeTestKernel(void)
{
	//mfp_debug("invokekernel\n");
	static const char srckernel[] =
	"__kernel void                       \n"
	"VecAdd(const __global float * M,    \n"
	"       const __global float * V,      \n"
	"             __global float * U,      \n"
    "        const unsigned int     N)      \n"
	"{                                   \n"
	"    size_t id =  get_global_id(0);  \n"
	"    if (id < N) {                   \n"
	"      U[id] = M[id] * V[id]; }         \n"
	"}                                   \n"
	;
	const cl_uint width = 10000U;
	const cl_uint mem_size = width * sizeof(cl_float);
	const cl_float3 values = { 0.5f, 2.0f, 1.0f };
	float M[width], V[width], U[width];
	fvFillArray(M, (cl_int)width, &values.s[0]);
	fvFillArray(V, (cl_int)width, &values.s[1]);
	// Get a list of supported platforms
	cl::Platform::get(& ocl::_platforms);

	// Pick first platform
	//cl_context_properties cprops[] = {
	//		CL_CONTEXT_PLATFORM,
	//		(cl_context_properties)(ocl::_platforms[ocl::instance()->_currPlatform])(), 0};
//	GLXContext glCtx = glXGetCurrentContext();
//	cl_context_properties cprops[] = {
//		CL_CONTEXT_PLATFORM,(cl_context_properties)(ocl::_platforms[ocl::instance()->_currPlatform])(),
//	    CL_GL_CONTEXT_KHR, (intptr_t)glCtx,
//		CL_GLX_DISPLAY_KHR, (intptr_t)glXGetCurrentDisplay(),
//		0
//	};
	cl_int status;
	cl::Context context = ocl::instance()->createContext(CL_DEVICE_TYPE_GPU);
	//mfp_log_debug("Status %d\n",status);

	// Query the set of devices attached to the context
	std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();

	// Create and program from source
	cl::Program::Sources sources(1, std::make_pair(srckernel, 0));
	cl::Program program(context, sources);
	//mfp_log_debug("Before building CL program\n");
	try {
		cl_int err = program.build(devices, "-Werror -cl-fast-relaxed-math");
	}
	catch (cl::Error err ) {
		std::string log( program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]));
		throw cl::Error(err.err(), log.c_str());
	}

	// Allocate memory and fill in arrays
	// Create buffer for M and copy host contents
	cl::Buffer MBuffer = cl::Buffer(
	            context,
	            CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
	            mem_size,
	            (void *) &M[0]);
	// Create buffer for V and copy host contents
	cl::Buffer VBuffer = cl::Buffer(
	            context,
	            CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
	            mem_size,
	            (void *) &V[0]);
	// Create buffer for U vector result
	cl::Buffer UBuffer = cl::Buffer(
	            context,
	            CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR,
	            mem_size,
	            (void *) &U[0]);

	// Create kernel object
	cl::Kernel kernel(program, "VecAdd");

	// Set kernel arguments
	kernel.setArg(0, MBuffer);
	kernel.setArg(1, VBuffer);
	kernel.setArg(2, UBuffer);
    kernel.setArg(3, width);

	// Create command queue
	cl::CommandQueue queue(context, devices[0], 0);

	// Do the work
	queue.enqueueNDRangeKernel(
	            kernel,
	            cl::NullRange,
	            cl::NDRange(width),
	            cl::NullRange);

	// Map UBuffer to host pointer. This enforces a sync with
	// the host backing space, remember we choose GPU device.
	float * output = (float *) queue.enqueueMapBuffer(
	            UBuffer,
	            CL_TRUE, // block
	            CL_MAP_READ,
	            0,
	            mem_size);

	// validate results
	for (cl_uint i = 0; i < width; ++i) {
		if (values.s[2] != U[i])
			throw cl::Error(TEST_FAILED, getString(errMsgs,FV_SIZEOF_ARRAY(errMsgs),TEST_FAILED));
	}

	// Finally release our hold on accessing the memory
	cl_int err = queue.enqueueUnmapMemObject(
	            UBuffer,
	            (void *) output);
	return(err);
}

cl_int ocl::integrateWithOpenGL(void)
{
	// Whether current platform support GL-CL operation
	if (_selfp->isExtensionSupported(GL_SHARING_EXTENSION,CL_DEVICE_TYPE_DEFAULT) > 0) {
		// Get platforms
		oclPlatform::get(& ocl::_platforms);
		// Specify properties for context
		cl_context_properties properties[] = {
		   CL_GL_CONTEXT_KHR, (cl_context_properties)glXGetCurrentContext(),
		   CL_GLX_DISPLAY_KHR, (cl_context_properties)glXGetCurrentDisplay(),
		   CL_CONTEXT_PLATFORM, (cl_context_properties)(ocl::_platforms[ocl::_selfp->_currPlatform])(),
		   0
		};
		//_cr cl::Context context(CL_DEVICE_TYPE_GPU, cprops);
		// Query the set of devices attched to the context
		// Query the set of devices attched to the context
			//std::vector<cl::Device> devices = context.getInfo<CL_CONTEXT_DEVICES>();
	}
	return(1);
}

ocl::ocl(cl_uint devType,cl_bool initGL, cl_bool quiet)
: _currPlatform(-1)
, _currDevice(-1)
{
	//mfp_log_debug("Ctr for %u\n",devType);
	// Check for OpenCL support
	if (ocl::_platforms.empty()) {
		//mfp_debug("Getting platforms\n");
		_err = oclPlatform::get(&_platforms);
	}

	// Find first platform which support this type of device
	for (cl_int ip(0); ip < _platforms.size(); ++ip) {
		if (!quiet) printPlatformInfo();
		_err = _platforms.at(ip).getDevices(devType, &ocl::_devices);
		if (initGL) {
			cl_int idev = isExtensionSupported(GL_SHARING_EXTENSION,devType);
			//mfp_debug("Here\n");
			if (idev > 0) {
				_currPlatform = ip; break;
			}
		} else {
			_currPlatform = ip;
		}
	}

	if (_currPlatform < 0) throw oclError(_currPlatform, "Specific OpenCL Platform not found");
	// For the found platform get specific first device
//	GLXContext glCtx = glXGetCurrentContext();
//	assert(glCtx != nullptr);
//	cl_context_properties cpsGL[] = {
//		CL_CONTEXT_PLATFORM,(cl_context_properties)(_platforms[_currPlatform])(),
// 	    CL_GL_CONTEXT_KHR, (cl_context_properties)glCtx,
//		CL_GLX_DISPLAY_KHR, (cl_context_properties)glXGetCurrentDisplay(),
//		0
//	};
	_context = createContext(devType,initGL);
	//_context = oclContext(_devices, NULL, NULL, NULL, &_err);
	// Query the set of devices attached to the context
	//mfp_debug("Searching for type: %u\n",devType);
	_devices = _context.getInfo<CL_CONTEXT_DEVICES>();
	_currDevice = 0;


	//mfp_log_debug("After creating ocl context\n");
	// Setup a context and a command queue
	_queue = oclCmdQ(_context, ocl::_devices.at(_currDevice), 0, &_err);

	_gridsize[0] = LOCAL_SIZE_X;
	_gridsize[1] = LOCAL_SIZE_Y;

}

ocl::~ocl()
{
	ocl::_devices.clear();
	ocl::_platforms.clear();
}

void ocl::clean()
{
	ocl::_devices.clear();
}

oclContext ocl::createContext(cl_uint devType, cl_bool supportGL)
{
	oclContext ctx;
	if (supportGL) {
		#ifdef WIN32
		cl_context_properties props[] = {
			CL_GL_CONTEXT_KHR, (cl_context_properties)wglGetCurrentContext(),
			CL_WGL_HDC_KHR, (cl_context_properties)wglGetCurrentDC(),
			CL_CONTEXT_PLATFORM, (cl_context_properties)cpPlatform,
			0
		};
		ctx = oclContext(devType, props, NULL, NULL, &_err);
		#else // UNIX
		cl_context_properties props[] = {
			CL_CONTEXT_PLATFORM,(cl_context_properties)(_platforms[_currPlatform])(),
			CL_GL_CONTEXT_KHR, (cl_context_properties)glXGetCurrentContext(),
			CL_GLX_DISPLAY_KHR, (cl_context_properties)glXGetCurrentDisplay(),
			0
		};
		ctx = oclContext(devType, props, NULL, NULL, &_err);
		#endif
	}
	else {
		cl_context_properties props[] = {
			CL_CONTEXT_PLATFORM, (cl_context_properties)(_platforms[_currPlatform])(),
		    0
		 };
		 ctx = oclContext(devType, props, NULL, NULL, &_err);
	}
    return ctx;
}

cl_uint ocl::countDevices(const oclPlatform& pl,const cl_device_type devType)
{
	if (!ocl::_devices.empty()) ocl::_devices.clear();
	cl_uint num = 0;
	cl_int err = pl.getDevices(devType, &ocl::_devices);
	if (err == CL_SUCCESS) {
		num = static_cast<cl_uint>(ocl::_devices.size());
		ocl::_devices.clear();
	}
	return num;
}

cl_uint ocl::isDeviceSupported(cl_device_type devType)
{
	cl_uint numDevs = 0;
	for (size_t i = 0; i < _devices.size(); ++i) {
		if (devType == _devices[i].getInfo<CL_DEVICE_TYPE>()) ++numDevs;
	}
	return numDevs;
}

// Find first device supported given extension
cl_int ocl::isExtensionSupported(const char * extString,const cl_device_type type)
{
	// Check if devices of current platform support given extension
	cl_int index(-1), i(0);
	for (vDevices::iterator it = ocl::_devices.begin(); it != ocl::_devices.end(); ++it, ++i) {
		// Check if given extension is on the list
		if ((*it).getInfo<CL_DEVICE_EXTENSIONS>().find(
			extString) == std::string::npos)  continue;
		cl_device_type type_ = it->getInfo<CL_DEVICE_TYPE>();
		if (type_ == type) {
			index = i;
			break;
		}
	}
	return(index);
}

cl_int ocl::loadProgram(const std::string& srcPath,bool quiet)
{
	// Load OpenCL source
	std::ifstream srcfile(srcPath);
	std::string src(std::istreambuf_iterator<char>(srcfile),
			std::istreambuf_iterator<char>(0));
	try {
		// Compile the device program
		_program = oclProgram(_context, oclProgram::Sources(1,
			std::make_pair(src.c_str(), src.size())));
		_err = _program.build(_devices, "-Werror");
	} catch(oclError err) {
		_err = err.err();
		std::cerr  << "ERROR: "
				<< err.what()
				<< "("
				<< err.err()
				<< ")"
				<< std::endl;
	}

	if (! quiet) {
		std::cout << "Build Status: " << _program.getBuildInfo<CL_PROGRAM_BUILD_STATUS>(_devices[_currDevice]) << std::endl;
		std::cout << "Build Options:\t" << _program.getBuildInfo<CL_PROGRAM_BUILD_OPTIONS>(_devices[_currDevice]) << std::endl;
		std::cout << "Build Log:\t " << _program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(_devices[_currDevice]) << std::endl;
	}

	return _err;

}
//} // end namespace

void ocl::printPlatformInfo()
{
	std::string name;
	for (vPlatforms::iterator it = _platforms.begin(); it != _platforms.end(); ++it)
	{
		(*it).getInfo(CL_PLATFORM_NAME, &name);
		std::cout << "Platform Name:\t\t" << name << std::endl;

		(*it).getInfo(CL_PLATFORM_VENDOR, &name);
		std::cout << "Platform Vendor:\t\t" << name << std::endl;

		(*it).getInfo(CL_PLATFORM_ICD_SUFFIX_KHR, &name);
		std::cout << "Platform ICD-suffix:\t\t" << name << std::endl;
	}
}

void ocl::printDeviceInfo()
{
	std::string name;
	for (vDevices::iterator it = _devices.begin(); it != _devices.end(); ++it)
	{
		(*it).getInfo(CL_DEVICE_NAME, &name);
		std::cout << "Name:\t\t" << name << std::endl;

		(*it).getInfo(CL_DEVICE_VENDOR, &name);
		std::cout << "Vendor:\t\t" << name << std::endl;

		(*it).getInfo(CL_DEVICE_VERSION, &name);
		std::cout << "Version:\t\t" << name << std::endl;
	}
}

cl_int initOpenCL(host_info_id host,cl_device_type devType,cl_bool supportGL)
{
	ocl::lock();

	cl_int err;

	if (! ocl::_selfp) {
		//mfp_debug("creating OpenCL engine\n");
		atexit(shutdownOpenCL);
		// Get all platforms
		cl::Platform::get(&ocl::_platforms);
		if (ocl::_platforms.empty()) return (-1);
	}

	try
	{
		// Reserve memory for host_info objects
		host = new host_info[ocl::_platforms.size()];
		device_info_t dev_info;
		int count = 0;
		// Scan each platform
		//mfp_debug("There is %d platforms\n",(int)ocl::_platforms.size());
		for (auto pl : ocl::_platforms) {
			//mfp_debug("CL_platform %d\n",count);
			// Get a name of first supported platform
			host[count].platform = pl.getInfo<CL_PLATFORM_NAME>();

			cl_int total = 0;
			device_info_t dev;
			for (int i = 0; i < FV_SIZEOF_ARRAY(_device_types); ++i) {
				//mfp_debug("here\n");
				dev.type = _device_types[i];
				dev.num = ocl::isDeviceSupported(_device_types[i]);
				//mfp_debug("After\n");
				if (dev.num > 0) {
					// Check if GL/CL is supported
					dev.gl_index = ocl::isExtensionSupported(GL_SHARING_EXTENSION,_device_types[i]);
					host[count].devices.push_back(dev);
					total += dev.num;
				}
			}
			++count;
		}
		//mfp_debug("befiorw new for ocl\n");
		// Clear before creating
		if (ocl::instance() == nullptr) ocl::instance() = new ocl(devType,false,true);


    		// Collect info about CL devices

		// Inovke test
		//mfp_debug("before test\n");
		err = ocl::invokeTestKernel();
		if (err != CL_SUCCESS) throw oclError(ocl::TEST_FAILED, "Test failed!");
		//mfp_log_info("OpenCL test passed\n");
		//err = ocl::integrateWithOpenGL();
		err = static_cast<cl_int>(ocl::_platforms.size());
	} catch (oclError& ex) {
		mfp_log_err("CL: %d %s\n",ex.err(),ex.what());
		err = -1;
		delete [] host;
		shutdownOpenCL();
	}
	ocl::unlock();
	return err;
}

cl_int loadDataCL()
{
	try
	{
		//mfp_log_debug("In loadCL\n");
		// Handle to singleton device
		ocl*& OCL(ocl::instance());
		assert(OCL != NULL);
		ViewManager& vm = ViewManagerInst();
		ModelCtrl* mc = &ModelCtrlInst();
		initWorkload(mc->RenderingContext());
		// Init image size
		tilesize[0] = vm.GetWidth();
		tilesize[1] = vm.GetHeight();
		pixcount = tilesize[0]*tilesize[1];
		image_size = pixcount * 4;
		// Init grid thrs sizes and number of pixels
		gwsize[0] = tilesize[0];// / IMG_WIDTH_CONSTRAINT;
		gwsize[1] = tilesize[1];

//		cl_context_properties cprops[] = {
//									CL_CONTEXT_PLATFORM, (cl_context_properties)(OCL->_platforms[0])(), 0 };
//		cl::Context context(CL_DEVICE_TYPE_GPU, cprops);
		oclContext& context = OCL->_context;
		// Init out buffers
		outBuffers[0] = cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, image_size, NULL);
		outBuffers[1] = cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_ALLOC_HOST_PTR, image_size, NULL);

		cl_image_format imgform = {CL_RGBA,CL_UNORM_INT8};


		// Allocate data for
		size_t ms_els  = mc->elData.size() * sizeof(elem_info);
		size_t ms_cofs = mc->coeffsData.size() * sizeof(coeffs_info_t);
		elDataBuffer = cl::Buffer(context, CL_MEM_READ_WRITE, ms_els, mc->elData.data());
		coeffBuffer  = cl::Buffer(context, CL_MEM_READ_WRITE, ms_cofs, mc->coeffsData.data());
		// Allocate buffers for grid structure
		// Grid info
		gridInfoBuffer = cl::Buffer(context,CL_MEM_READ_ONLY,sizeof(grid_t), &mc->gridData);
		// Grid
		gridBuffer = cl::Buffer(context,CL_MEM_READ_ONLY,sizeof(cl_int)*gridcount,&mc->C_ptr);
		// Cells
		cellsBuffer = cl::Buffer(context,CL_MEM_READ_WRITE,sizeof(cl_int)*cellcount,&mc->L_ptr);
		// Buffer for rctx


		// Pixels
		//cl::Buffer PixelsBuffer = cl::Buffer(context, CL_MEM_READ_WRITE, sizeof(struct pixel_s)*pixcount,NULL);
		//
		image = new pixel_s [pixcount];
		delete [] image;

		// Load kernels
		std::ifstream srcfile(kernels[1]);
		std::string src(std::istreambuf_iterator<char>(srcfile),
								std::istreambuf_iterator<char>(0));

		// Compile the device program
		oclProgram prog(context, cl::Program::Sources(1,
						std::make_pair(src.c_str(), src.size())));
		prog.build(ocl::_devices, "-Werror");

		kern_priray = oclKernel(prog, "MinMaxRenge");
/*
		cl::Buffer CoeffBuffer = cl::Buffer(
				context,
				CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				mem_CoeffSize,
				(void *) coeffs.data());

		cl::Buffer MnMxBuffer = cl::Buffer(
				context,
				CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR,
			    2 * sizeof(cl_float),
			    (void *) &MinMax[0]);

		// Load OpenCL kernel
		std::ifstream srcfile(kernels[1]);
		std::string src(std::istreambuf_iterator<char>(srcfile),
							std::istreambuf_iterator<char>(0));
		*/
/*

	// Fill in with indices
	std::vector<cl_int> indices(totalSize,0);
	FemViewer::Mesh::FillArrayWithIndicesOfElemNodes(mesh_,indices.data());
	assert(indices.size() == totalSize);
	const cl_uint mem_IndSize = sizeof(cl_int) * totalSize;

	// Fill in coordinates
	std::vector<CoordType> vertices(3*mesh_->GetMaxNodeId(),cl_max);
	FemViewer::Mesh::FillArrayWithElemCoords(mesh_,vertices.data());
	const cl_uint mem_VertSize = sizeof(cl_float) * vertices.size();

	// Fill in solution coefficient
	cl_uint nCoeffs= field_->CountCoeffs();
	assert(nCoeffs > 0 && nCoeffs <= mesh_->GetElems().size()*450);
	std::vector<cl_float> coeffs(nCoeffs,-10e200);
	cl_uint coeffSize = FemViewer::Field::FillArrayWithCoeffs(field_,coeffs.data());
	assert(coeffSize == coeffs.size());
	const cl_uint mem_CoeffSize = sizeof(cl_float) * coeffSize;

	cl_context_properties cprops[] = {
				CL_CONTEXT_PLATFORM, (cl_context_properties)(OCL->_platforms[0])(), 0};
	cl::Context context(CL_DEVICE_TYPE_GPU, cprops);

	cl::Buffer IndicesBuffer = cl::Buffer(
				context,
				CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
				mem_IndSize,
				(void *) indices.data());

	cl::Buffer CoordBuffer = cl::Buffer(
			context,
			CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			mem_VertSize,
			(void *) vertices.data());

	cl::Buffer CoeffBuffer = cl::Buffer(
			context,
			CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
			mem_CoeffSize,
			(void *) coeffs.data());

	cl::Buffer MnMxBuffer = cl::Buffer(
			context,
			CL_MEM_WRITE_ONLY | CL_MEM_USE_HOST_PTR,
		    2 * sizeof(cl_float),
		    (void *) &MinMax[0]);


	kernel.setArg(0, CoordBuffer);
	kernel.setArg(1, IndicesBuffer);
	kernel.setArg(2, CoeffBuffer);
	kernel.setArg(3, MnMxBuffer);
	kernel.setArg(4, mesh_->GetElems().size());

	// Create command queue
	cl::CommandQueue queue(context, OCL->_devices[0], 0);
	queue.enqueueWriteBuffer(IndicesBuffer,CL_TRUE,0,mem_IndSize, indices.data());
	queue.enqueueWriteBuffer(CoordBuffer,CL_TRUE,0,mem_VertSize, vertices.data());
	queue.enqueueWriteBuffer(CoeffBuffer,CL_TRUE,0,mem_CoeffSize, coeffs.data());

	// Do the work
	queue.enqueueNDRangeKernel(
		            kernel,
		            cl::NullRange,
		            cl::NDRange(mesh_->GetElems().size()),
		            cl::NullRange);

	queue.finish();
	queue.enqueueReadBuffer(MnMxBuffer,CL_TRUE,0,2*sizeof(cl_float), &MinMax[0]);
*/mfp_log_debug("leaving loadCL");
	} catch(oclError err) {
		error = err.err();
		std::cerr  << "ERROR: "
				<< err.what()
				<< "("
				<< err.err()
				<< ")"
				<< std::endl;
		return -1;
	}
	return 0;
}

cl_int displayCL()
{
	return 0;
}

cl_int runCL()
{
	mfp_log_debug("In runCL\n");
	for (int f = 0; f < nr_frames-1; f++) {
	    (void)ocl::next_frame(NULL);
	}
	(void)ocl::last_frame(NULL);

	double delta = ocl_timer.get_duration_sec();
	printf ("%d Frames took %f seconds. Rate = %f Mpixels/sec, %f frames/second\n",
		          nr_frames, delta,
		          (double) (img_width * img_height) * (double) (nr_frames) /
		          (1000000.0 * (double) (delta)), nr_frames/delta);
	return 0;
}

cl_int bindOpenCL(cl_uint GLpbo, cl_int width, cl_int height)
{
	if (glIsBuffer(GLpbo) == GL_FALSE || ocl::instance() == nullptr) return(-1);

	cl_int err;
	ocl::_selfp->_bufferGL = oclBufferGL(ocl::_selfp->_context,
			CL_MEM_WRITE_ONLY, GLpbo, &err);
	if (err != CL_SUCCESS) {
		std::cerr << "Error in binding GL-CL\n";
		throw oclError(err,"ERROR in binding GL context\n");
	}

	ocl::_selfp->_gridsize[0] = roundUp(LOCAL_SIZE_X, width);
	ocl::_selfp->_gridsize[1] = roundUp(LOCAL_SIZE_Y, height);

	return(1);
}

void shutdownOpenCL(void)
{
	using namespace FemViewer;
	mfp_debug("OpenCL: shutdown\n");
	ocl::lock();
	if (ocl::_selfp != nullptr) { delete ocl::_selfp; ocl::_selfp = nullptr; }
	ocl::unlock();
}


static void calculateView(RContext* mrc, camera_t * cam)
{
  float tmpEye[4];
  float tmpLightP[4];

  /*
   * First apply the view transformations to the initial eye, look at,
   * and up.  These will be used later to determine the basis.
   */
  int i, j;
  float mag;

  /* eye starts on the unit sphere */
  float eyeStart[4] = { mrc->camera.pos.x, mrc->camera.pos.y, mrc->camera.pos.z, 1.0f };

  /* initially look at the origin */
  float lookatStart[4] = { mrc->camera.target.x, mrc->camera.target.y, mrc->camera.target.z, 1.0f };

  /* up is initially along the y-axis */
  float upStart[4] = { mrc->camera.up.x, mrc->camera.up.y, mrc->camera.up.z, 0.0f };

  /* point light location */
  static float lookAt[4], up[4];


  /* translate the eye and look at points */
  eyeStart[0]    += cam->translate[0];
  eyeStart[1]    += cam->translate[1];
  eyeStart[2]    += cam->zoom;
  lookatStart[0] += cam->translate[0];
  lookatStart[1] += cam->translate[1];
  lookatStart[2] += cam->zoom;

  /* rotate eye, lookat, and up by multiplying them with the current rotation matrix */
  for (i = 0; i < 4; i++) {
    tmpEye[i] = 0.0f;
    lookAt[i] = 0.0f;
    up[i] = 0.0f;

    for (j = 0; j < 4; j++) {
      tmpEye[i] += cam->curRotation[i * 4 + j] * eyeStart[j];
      lookAt[i] += cam->curRotation[i * 4 + j] * lookatStart[j];
      up[i] += cam->curRotation[i * 4 + j] * upStart[j];
    }
  }


  /* Now we construct the basis: */
  /*   N = (look at) - (eye)     */
  /*   T = up                    */
  /*   B = N x T                 */

  /* find and normalize N = (lookat - eye) */
  for (i = 0; i < 3; i++) N[i] = lookAt[i] - tmpEye[i];
  mag = 1.0f / sqrt (N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
  for (i = 0; i < 3; i++) N[i] *= mag;

  /* find and normalize T = up */
  for (i = 0; i < 3; i++) T[i] = up[i];
  mag = 1.0f / sqrt (T[0] * T[0] + T[1] * T[1] + T[2] * T[2]);
  for (i = 0; i < 3; i++) T[i] *= mag;

  /* find B = N x T (already unit length) */
  B[0] = N[1] * T[2] - N[2] * T[1];
  B[1] = N[2] * T[0] - N[0] * T[2];
  B[2] = N[0] * T[1] - N[1] * T[0];

  /* move the light a little bit up and to the right of the eye. */
  for (i = 0; i < 3; i++) {
    tmpLightP[i] = tmpEye[i] - B[i] * 0.5f;
    tmpLightP[i] += T[i] * 0.5f;
  }

  rctx.eyeP.s[0] = tmpEye[0];
  rctx.eyeP.s[1] = tmpEye[1];
  rctx.eyeP.s[2] = tmpEye[2];
  rctx.eyeP.s[3] = tmpEye[3];

  rctx.lightP.s[0] = tmpLightP[0];
  rctx.lightP.s[1] = tmpLightP[1];
  rctx.lightP.s[2] = tmpLightP[2];
  rctx.lightP.s[3] = tmpLightP[3];
}

static void initWorkload(RContext* mrc)
{
	mfp_log_debug("initWorkload\n");
  float alpha;                  		/* height for aspect ratio */
  float beta;                   		/* width for aspect ratio */
  Matrix<float> w2c = Matrix<float>::LookAt(mrc->camera.pos,
		  mrc->camera.target, mrc->camera.up);

  /* Set initial view parametrs */
  rc.translate[0]    = w2c[12];
  rc.translate[1]    = w2c[13];
  rc.zoom            = w2c[14];
  for (int i=0;i<16;++i) {
	  rc.curRotation[i] = w2c[i];
}
//  rc.curRotation[0]  = 0.0F;
//  rc.curRotation[1]  = 1.0F;
//  rc.curRotation[2]  = 0.0F;
//  rc.curRotation[3]  = 0.0F;
//  rc.curRotation[4]  = -1.0F;
//  rc.curRotation[5]  = 0.0F;
//  rc.curRotation[6]  = 0.0F;
//  rc.curRotation[7]  = 0.0F;
//  rc.curRotation[8]  = 0.0F;
//  rc.curRotation[9]  = 0.0F;
//  rc.curRotation[10] = 1.0F;
//  rc.curRotation[11] = 0.0F;
//  rc.curRotation[12] = 0.0F;
//  rc.curRotation[13] = 0.0F;
//  rc.curRotation[14] = 0.0F;
//  rc.curRotation[15] = 1.0F;
  rc.shadows = 0;

  rc.fov = mrc->camera.fov;
  rc.aspect = ((float) img_width) / ((float) img_height);
  rctx.window_size.s[0] = (float) img_width;
  rctx.window_size.s[1] = (float) img_height;
  rctx.maxIterations = 100;
  rctx.epsilon = 0.001f;
  rctx.stride = img_width;

  calculateView (mrc,&rc);

  beta = tan ((rc.fov * M_PI / 180.0f) / 2.0f); /*find height */
  alpha = beta * rc.aspect;     /*find width */

  /* rendering region: upper left corner */
  rctx.dir_top_start.s[0] = -alpha * T[0] - beta * B[0] + N[0];
  rctx.dir_top_start.s[1] = -alpha * T[1] - beta * B[1] + N[1];
  rctx.dir_top_start.s[2] = -alpha * T[2] - beta * B[2] + N[2];
  mfp_log_debug("Top start: %f %f %f\n",rctx.dir_top_start.s[0] ,rctx.dir_top_start.s[1],rctx.dir_top_start.s[2]);
  /* rendering region: lower left corner */
  rctx.dir_bottom_start.s[0] = -alpha * T[0] + beta * B[0] + N[0];
  rctx.dir_bottom_start.s[1] = -alpha * T[1] + beta * B[1] + N[1];
  rctx.dir_bottom_start.s[2] = -alpha * T[2] + beta * B[2] + N[2];
  mfp_log_debug("Bottom start: %f %f %f\n",rctx.dir_bottom_start.s[0] ,rctx.dir_bottom_start.s[1],rctx.dir_bottom_start.s[2]);
  /* rendering region: lower right corner */
  rctx.dir_bottom_stop.s[0] = alpha * T[0] + beta * B[0] + N[0];
  rctx.dir_bottom_stop.s[1] = alpha * T[1] + beta * B[1] + N[1];
  rctx.dir_bottom_stop.s[2] = alpha * T[2] + beta * B[2] + N[2];
  mfp_log_debug("Bottom stop: %f %f %f\n",rctx.dir_bottom_stop.s[0] ,rctx.dir_bottom_stop.s[1],rctx.dir_bottom_stop.s[2]);
}

void* ocl::next_frame(size_t **size)
{
	cl_int err;
	int    next;				/* next frame index */
	void * frame_buffer;              	/* image frame buffer to display */
	//cl_event map_event;
	//cl_event unmap_event = NULL;
	oclEvent map_event;
	oclEvent unmap_event; // z def is NULL
	/*-----------------------------------------------------------------------------------
	 * Advance to the parameterization to next frame
	 *-----------------------------------------------------------------------------------
	 */
	// nextMu();

	/* ----------------------------------------------------------------------------------
   * The rendering pipeline when double buffering. This ordering of requests is
   * design to support both in-order and out-of-order queues.
   *
   *          Compute buffer 0
   *                 |
   *      .--------->|
   *      |          v
   *      |     Map buffer 0
   *      |          |
   *      |          v
   *      |    Unmap buffer 1 (if mapped)
   *      |          |
   *      |          v
   *      |   Compute buffer 1
   *      |          |------------------> return buffer 0
   *      |          v
   *      |     Map buffer 1
   *      |          |
   *      |          v
   *      |    Unmap buffer 0
   *      |          |
   *      |          v
   *      |   Compute buffer 0
   *      |          |------------------> return buffer 1
   *      `----------'
   *
   * ----------------------------------------------------------------------------------
   */
  /* If this is the first frame, prime the pipeline by computing the first buffer.
   */
	  if (first_frame) {
	    first_frame = 0;
//	    cluRunKernel (clu, kernel, &compute_event, 3,
//			  sizeof(currMu), currMu,
//			  sizeof(cl_mem), &outBuffer[current_frame],
//			  sizeof(struct julia_context), &jc);
//	    nextMu();
	    _selfp->_queue.enqueueNDRangeKernel(kern_priray,cl::NullRange,cl::NDRange(image_size),cl::NullRange,&events,&evt);
	    events = { evt };
	    // Tu cos update
	  }

	  next = 1 - current_frame;
	  mapped_addr[current_frame] = (cl_uchar*)_selfp->_queue.enqueueMapBuffer(outBuffers[current_frame],false,CL_MAP_READ, 0,
			  image_size, &events, &mapEvt, &_selfp->_err);

	 // release event
	 //evt = NULL



//	  mapped_addr[current_frame] = clEnqueueMapBuffer (_selfp->_queue.enqueueMapBuffer(), outBuffer[current_frame],
//							   CL_FALSE, CL_MAP_READ, 0, img_size,
//							   1, &compute_event, &map_event, &err);
//	  CLU_CHECK_ERROR ("clEnqueueMapBuffer", err);
//	  CLU_CHECK_ERROR ("clReleaseEvent", clReleaseEvent(compute_event));
//
	  //Unmap the next frame if it was previously mapped.
	  if (mapped_addr[next]) {
		  _selfp->_queue.enqueueUnmapMemObject(outBuffers[next],mapped_addr[next],&events,&unmapEvt);
		  //CLU_CHECK_ERROR ("clEnqueueUnmapMemObject",
		  //	     clEnqueueUnmapMemObject (commands, outBuffer[next], mapped_addr[next], 0, NULL, &unmap_event));
		  mapped_addr[next] = NULL;
		  //cluSetKernelDependency(clu, kernel, 1, &unmap_event);
	  }
//
//	  /* Wait for the map to complete */
//	  CLU_CHECK_ERROR ("clWaitForEvents", clWaitForEvents(1, &map_event));
//	  CLU_CHECK_ERROR ("clReleaseEvent", clReleaseEvent(map_event));
//
//
//	  /* Generate the next frame
//	   */
//	  cluRunKernel (clu, kernel, &compute_event, 3,
//			sizeof(currMu), currMu,
//			sizeof(cl_mem), &outBuffer[next],
//			sizeof(struct julia_context), &jc);
//
//	  if (unmap_event) {
//	    CLU_CHECK_ERROR ("clReleaseEvent", clReleaseEvent(unmap_event));
//	  }
//
//	  frame_buffer = (void *)mapped_addr[current_frame];

	  /*------------------------------------------------------------------------------------
	   * Return the resulting frame buffer pointer and size information/
	   *------------------------------------------------------------------------------------
	   */
	  current_frame = next;

	  return (frame_buffer);
}


void *ocl::last_frame(int *size)
{
	cl_int err;
	int    frame_size = 0;		/* default frame size */
	char * frame_buffer;              	/* image frame buffer to display */

  /* Compute the frame if no frame is already in the pipeline
   */
  if (first_frame) {
    first_frame = 0;
//    cluRunKernel (clu, kernel, &compute_event, 3,
//		  sizeof(currMu), currMu,
//		  sizeof(cl_mem), &outBuffer[current_frame],
//		  sizeof(struct julia_context), &jc);
  }
  /* Map the computed frame into the host address space
   */
//  frame_buffer = clEnqueueMapBuffer (commands, outBuffer[current_frame],
//				     CL_TRUE, CL_MAP_READ, 0, img_size,
//				     1, &compute_event,
//				     NULL, &err);
//  CLU_CHECK_ERROR ("clEnqueueMapBuffer", err);
//  CLU_CHECK_ERROR ("clReleaseEvent", clReleaseEvent(compute_event));
//  CLU_CHECK_ERROR ("clFinish", clFinish(commands));
//
//  compute_event = NULL;

  if (size) *size = frame_size;
  return (frame_buffer);
}





