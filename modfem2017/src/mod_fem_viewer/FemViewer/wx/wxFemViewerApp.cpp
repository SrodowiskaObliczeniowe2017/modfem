#include "AppImpl.hpp"
#include "ViewManager.h"
#include "ModelControler.h"
#include "wxFemViewerApp.h"
#include "wxMainWindow.h"

//#include "../../utils/fv_exception.h"
#ifdef _USE_FV_IO_MGR
#include "../../include/IOManager.h"
#endif

#ifdef _USE_FV_EXT_MOD
#include "../../include/ApproxModule.h"
#endif
#include "fv_compiler.h"
#include "fv_config.h"
//#include "fv_txt_utls.h"
#include "Log.h"
//#include <wx/log.h>
#include<ostream>
//#include<X11/Xlib.h>



namespace FemViewer {
namespace WxGUI		{
	
IMPLEMENT_APP_NO_MAIN(wxFemViewerApp)


Thread* wxFemViewerApp::m_parent = nullptr;

//wxFemViewerApp * wxFemViewerApp::GetInstance()
//{
//	if (IsInit()) return static_cast<wxFemViewerApp*>(wxApp::GetInstance());
//	return nullptr;
//}

int wxFemViewerApp::Init(int argc,char **argv,Thread* parent)
{
	wxDISABLE_DEBUG_SUPPORT();
	//gtk_disable_setlocale ();
	//int result = XInitThread ();
	if (parent) m_parent = parent;
	argc = 1; // ignore options, only program name
	//wxApp::SetInstance( new wxFemViewerApp(parent->_shutdown))
	bool result = wxEntryStart(argc,argv);
	//bool result = wxInitialize(argc, argv);
	assert(result);
	//mfp_log_debug("Init");
	wxTheApp->OnInit();
	//wxTheApp->OnRun();
	return true;
}


int wxFemViewerApp::Init2(int argc,char **argv)
{
	wxDISABLE_DEBUG_SUPPORT();
	//mfp_log_debug("init2");
	//if( !glutGet( GLUT_INIT_STATE ) )
	//	glutInit( &argc,argv );

	//ModelCtrlInst().Do();
	return 0;
}

bool wxFemViewerApp::IsInit()
{
	return !!(wxApp::GetInstance() != NULL);
}

void wxFemViewerApp::Refresh()
{
	if (IsInit()) {
		wxTheApp->GetTopWindow()->Raise();
		//mfp_debug("Sending event\n");
		wxCommandEvent evt(wxEVT_MENU, IDM_FILE_RELOAD);
		wxQueueEvent(wxTheApp->GetTopWindow(), evt.Clone());
		//mfp_debug("evt was send\n");
	}
}



void wxFemViewerApp::ShutDown()
{
	//mfp_log_debug("ShutDown");
	wxGetApp().m_shutdown = true;
}



wxFemViewerApp::wxFemViewerApp()
: m_shutdown(false)
, m_pframe(nullptr)
, m_pcanvas(nullptr)
{
	//mfp_log_debug("ctr\n");
	// Some posts in the internet says that it has to be placed here?
	//XInitThreads();
}

wxFemViewerApp::~wxFemViewerApp(void)
{
}

bool wxFemViewerApp::OnInit()
{
	//mfp_log_debug("OnInit");

	if (!wxApp::OnInit())
		return false;

	wxInitAllImageHandlers();

	// Create main window
	m_pframe = new wxMainWindow(imageWidth,imageHeight,wxT(MOD_FEM_VIEWER_NAME));


//	m_pcanvas = new wxRenderCanvas(m_pWnd, wxID_ANY, gl_attrib, m_pWnd->GetPosition(),
//		m_pWnd->GetClientSize(),wxDEFAULT_FRAME_STYLE );
//	if (!m_pcanvas) throw "Can't initialize OpenGL machine!";


	//m_pWnd->Refresh();
	//mfp_log_debug("Show iwndow now!\n");
	//m_pWnd->Show();
	//wnd_ptr->DoInit();

	//mfp_log_debug("Window is on top\n");
//#if wxUSE_IMAGE
//    wxInitAllImageHandlers();
//#endif
    //m_pWnd->Show(true);
    //SetTopWindow(m_pWnd);
	return true;
}

int wxFemViewerApp::OnRun()
{
	//mfp_log_debug("OnRun\n");
	//m_pWnd->DoInit();
	assert(true == m_pframe->IsShownOnScreen());

	m_pframe->Update();
	//m_pWnd->Show(true);
	m_pframe->Raise();
	assert(m_pframe->IsEnabled());
	SetTopWindow(m_pframe);
	wxGetApp().GetTopWindow()->SetFocus();
	ViewManagerInst().Init(ModelCtrlInst().Boundary());
	//assert(m_pWnd->IsActive());
	//m_pWnd->StartRenderThread();
	//mfp_log_debug("Start rendering loop\n");
	assert(m_parent);
	m_parent->Notify();
	wxApp::OnRun();
	return OnExit();
}

//int wxFemViewerApp::Update()
//{
//	wxCommandEvent event(wxEVT_MENU, IDM_FILE_RELOAD);
//	GetInstance()->QueueEvent(&event);
//	return 1;
//}



int wxFemViewerApp::OnExit()
{
  //  mfp_log_debug("OnExit");
	int result = wxApp::OnExit();
	wxApp:CleanUp();
    m_parent = nullptr;
    //wxUninitialize();
   //	mfp_debug("Po entrycleanup\n");
   	wxEntryCleanup();
	// Returns
	return result;
}

}// end namespace WxGUI
}// end namespace FemViewer


