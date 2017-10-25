/*
 * FemViewerApp.cpp
 *
 *  Created on: 26-01-2012
 *      Author: dwg
 */
#include "fv_config.h"
#include "defs.h"
#include "../../utils/fv_exception.h"
#include "FemViewerApp.h"
#include "ModelControler.h"
#ifdef _USE_FV_IO_MGR
#include"../../include/IOManager.h"
#endif
#ifdef _USE_FV_EX_MOD
#include"../../include/ApproxModule.h"
#endif


/* Graphic GUI */
#ifdef _USE_FV_WX
#include "../wx/wxFemViewerApp.h"
#else
#include "GlutWindow.h"
#endif

namespace FemViewer {



FemViewerApp::FemViewerApp(int argc, char** argv, int posX, int posY,
						   int width, int height, const char* title)
:  _pwnd(NULL)
  ,_pmodel(NULL)
  ,_running(0)
{
	try
	{
		#ifdef _USE_FV_WX
		// Start application
		WxGUI::wxFemViewerApp::Init(argc,argv,posX,posY,width,height,title);
		//_pwnd = static_cast<WindowFV*>(&WxGUI::wxGetApp().GetWindow());
		#else
		GlutGUI::GlutWindow::Init(argc,argv);
		_pwnd = new GlutGUI::GlutWindow(posX,posY,width,height,title);
		#endif

		#ifdef _USE_FV_IO_MGR
		if (IOMgr_Initialize(argc,(const char**)argv) != FV_SUCCESS){
			throw( fv_exception("I/O manager initializing error!"));
		}
   	   #endif

   	   #ifdef _USE_FV_EXT_MOD
		if (init_approx_modules() < 0) {
			throw( fv_exception("Approximation module error!"));
		}
   	   #endif

		//_pmodel = ModelControler::Init(this);
		//FV_ASSERT(_pmodel!=NULL);
	}
	catch(fv_exception& excpt)
	{
		std::cout << "Exception: " << excpt << std::endl;
		throw(-1);
	}
}

FemViewerApp::~FemViewerApp(void)
{
	destructor_debug("FemViewerApp")
#ifdef _USE_FV_IO_MGR
	// Delete the filemanager
	IOMgr_Destroy();
#endif

#ifdef _USE_FV_EXT_MOD
	// Delete the approxmodule
	destroy();
#endif
	/// Cleanup
#ifdef _USE_FV_WX
	WxGUI::wxGetApp().OnExit();
	wxEntryCleanup();
#else
	delete static_cast<GlutGUI::GlutWindow*>(_pwnd);
#endif
	_pwnd = NULL;
	_running = 0;
}

int FemViewerApp::Run(void)
{
	return 0;
	_running++;
#ifdef _USE_FV_WX
	return WxGUI::wxGetApp().OnRun();
#else
	return (static_cast<GlutGUI::GlutWindow*>(_pwnd))->MainLoop();
#endif

}

bool FemViewerApp::Close()
{
	_pwnd->ShutDown();
	return true;
}

} // end namespace

