#include "WindowFV.h"
#include "ViewManager.h"
#include "ModelControler.h"
#include "Log.h"

namespace FemViewer {

	//Window* Window::_self = NULL;
    bool mfvWindow::init(int argc,char **argv)
    {
    	//static mfvWindow wnd;
    	ModelCtrlInst().Do();
    	return true;
    }

	mfvWindow::mfvWindow() : m_pview(ViewManagerInst()), m_pmodel(ModelCtrlInst())
	{
		mfp_debug("mfvWindow ctr\n");
	}

	mfvWindow::~mfvWindow(void)
	{
		mfp_debug("mfvWindow dtr\n");
	}



} // end namespace FemViewer
