/*
 * AppImpl.cpp
 *
 *  Created on: 4 sty 2014
 *      Author: Paweł Macioł
 */
#include "Log.h"
#include "types.h"
#include "fv_config.h"
#include "fv_timer.h"
#include "fv_txt_utls.h"
#include "ocl.h"
#include "AppImpl.hpp"
#include "ViewManager.h"
#include "ModelControler.h"
#include "View.h"

// GUI interface
#if defined(WIN32) && defined(VIEWER_WIN_GUI)
# include "win/Common.h"
# include "win/WinApp.h"
#elif defined(VIEWER_FGLT_GUI)
# include "glt/GlutWindow.h"
#elif defined(VIEWER_WX_GUI)
# include "wx/wxFemViewerApp.h"
#else
# error "Unknown GUI-type"
#endif

cl_int numelems;
struct triangle_s * triangles=NULL;
char** mtllist=NULL;
char* mtllib=NULL;
cl_int* grid_data=NULL;
cl_int* cell_data=NULL;
grid_t gridinfo;
size_t imgsize[2];

extern struct pixel_s * image;

extern cl_mem pixels;
extern cl_event shade;
extern cl_event extraevent;
extern cl_command_queue cq;
extern cl_int pixcount;
extern const bool initgl;


int texnum=0;


#ifdef DEBUG
extern unsigned long waitfor_trace;
extern unsigned long waitfor_count;
extern unsigned long waitfor_light;
extern unsigned long waitfor_shade;
#endif

namespace FemViewer {



ModuleThread::ModuleThread(int argc,char** argv)
: _shutdown(0)
, _argc(argc)
, _argv(argv)
, _running(false)
{
	//mfp_debug("ModuleThread ctr\n");
}

ModuleThread::~ModuleThread()
{
	//mfp_debug("ModuleThread dtr\n");
	//Thread::Terminate();
}



void ModuleThread::Start()
{
	//mfp_log_debug("Starting new thread procedure\n");
	if (!this->_handle) Thread::Start();
	else { // tell to refresh
		//mfp_debug("refresh screen\n");
		GUIEngine::Refresh();
	}
	// Wait for Notify from slave
	//mfp_log_debug("Waiting for wakeup from slave\n");
//#ifndef WIN32
//	usleep(10);
//#endif
	Thread::WaitForNotify();
	_running = true;
	//mfp_log_debug("OK, I have been wakeuped\n");
}

int ModuleThread::Run()
{
	//mfp_log_debug("ModuleThread Run\n");
	int result = run_application(this->_argc,this->_argv,initgl,this);
	//mfp_log_debug("Ending start thread\n");
	Thread::Notify();
	this->_handle = 0;
	return result;
}

int ModuleThread::Update()
{

}

void ModuleThread::Stop()
{
	mfp_log_debug("Waiting for slave");
	_shutdown = 1;
	GUIEngine::ShutDown();
	Thread::Stop();
	_running = false;
}



int run_application(int argc, char** argv, bool initgl, Thread* parent)
{
	// Init GUI and OpenGL
	//mfp_log_debug("Main core of application");
	//RenderManager<GLCore>::Init(argc,argv);
	ModelCtrl & mc = ModelCtrlInst();
	ViewManager & vm = ViewManagerInst();


	if (initgl) {
		static int nr_run;
		//mfp_log_info("Initializing %d display device...",nr_run++);
		//GLCore::init(argc,argv);
		//mc.InitData();
		vm.GLSupport() = GUIEngine::Init(argc,argv,parent);
		if (!mc.Do()) return -1;
		GUIEngine::GetInstance()->OnRun();
		//mfp_log_info("Initializing display device...done");
	}
	else {
		// Correct width/height
		mfp_log_info("Initializing data...");

		//t.start();
		mc.InitData();

		gridcount = mc.gridData.resolution[0] * mc.gridData.resolution[1] * mc.gridData.resolution[2] + 1;
		cellcount = mc.C_ptr[gridcount-1];
		printf("created %d cells, with %d entries\n",gridcount-1,cellcount);
		printf("grid resolution: %d %d %d\n",mc.gridData.resolution[0], mc.gridData.resolution[1], mc.gridData.resolution[2]);
		img_width  = IMG_WIDTH_CONSTRAINT  * ((img_width  + (IMG_WIDTH_CONSTRAINT-1))  / IMG_WIDTH_CONSTRAINT);
		img_height = IMG_HEIGHT_CONSTRAINT * ((img_height + (IMG_HEIGHT_CONSTRAINT-1)) / IMG_HEIGHT_CONSTRAINT);
		vm.SetWidth(img_width);
		vm.SetHeight(img_height);
		int image_size = img_width * img_height * 4;
		//initWorkload(mc.RenderingContext());
		host_info_id host = NULL;
		int size = initOpenCL(host,CL_DEVICE_TYPE_GPU,false);
		loadDataCL();
		//runCL();
		// move to atexit
		shutdownOpenCL();
		if (host) delete [] host;
	}

	//if (initgl)		GUIEngine::GetInstance()->OnRun();

	// Init renderer
	//mfp_log_info("Initializing data...done");



	return EXIT_SUCCESS;
}




} // end namespace





