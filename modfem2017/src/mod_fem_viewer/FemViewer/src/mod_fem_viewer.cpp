#include<stdlib.h>
#include<assert.h>

/* Internal configuration */
#include "Log.h"
#include "mod_fem_viewer.h"
#include "fv_config.h"
#include "fv_exception.h"
//#include "fv_threads.h"

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

#include "AppImpl.hpp"

using namespace FemViewer;

#ifndef _USE_FV_LIB

FILE * out = stdout;
typedef AppImpl<AppEngine> app_core;

int main(int argc, char **argv)

#else

typedef AppImpl<GUIEngine *> app_core;

enum {
	APP_OFF = 0,
	APP_INIT,

};


femviewer_export
int init_mod_fem_viewer(int argc, char **argv, FILE *out)
#endif
{
	static ModuleThread app_thread(argc,argv);
	int result(-1);

	try {
		// Init logger
		Log::Init(out);
		// Starrting thread of graphics
		mfp_log_info("Launch graphic module....\n");
		app_thread.Start();
		mfp_log_info("...Launch graphic module: Done\n");
		//result = app_core::run(argc,argv);
	}
	catch (const fv_exception &ex) {
		mfp_log_err("%s\n",ex.what());
	}
	catch (...) {
		mfp_log_err("Some errors appear in graphics module!\n");
	}
	return result;
}

