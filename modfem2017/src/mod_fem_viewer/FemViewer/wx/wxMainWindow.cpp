#include "fv_inc.h"
#include "wxFemViewerApp.h"
#include "wxMainWindow.h"
#include "wxOpenGL.h"
#include "wxApproxModWnd.h"
#include "wxLegendDialog.h"
#include "wxScreenSettings.h"
#include "wxSolCalcDialog.h"
#include "wx/busyinfo.h"
#include "wx/image.h"
#include "wx/toolbar.h"
#include "wx/config.h"
//#include "wx/menu.h"
#include "wx/docview.h"
#include "fv_config.h"
#include "fv_exception.h"
#include "fv_txt_utls.h"
#include "fv_timer.h"
#include "BaseField.h"
#include "ViewManager.h"
#include "ModelControler.h"
#include "AppImpl.hpp"

#include "Log.h"

#include <sstream>
#include <algorithm>


// define this to use XPMs everywhere (by default, BMPs are used under Win)
// BMPs use less space, but aren't compiled into the executable on other platforms
#ifdef __WXMSW__
# define USE_XPM_BITMAPS 0
#else
# define USE_XPM_BITMAPS 1
#endif

#if USE_XPM_BITMAPS && defined(__WXMSW__) && !wxUSE_XPM_IN_MSW
# error You need to enable XPM support to use XPM bitmaps with toolbar!
#endif // USE_XPM_BITMAPS

// If this is 1, the sample will test an extra toolbar identical to the
// main one, but not managed by the frame. This can test subtle differences
// in the way toolbars are handled, especially on Mac where there is one
// native, 'installed' toolbar.
#define USE_UNMANAGED_TOOLBAR 0
#define MAX_HISTORY_FILES 	  8

// Define this as 0 for the platforms not supporting controls in toolbars
#define USE_CONTROLS_IN_TOOLBAR 1

// The toolbar
const int ID_TOOLBAR = 500;
static const long TOOLBAR_STYLE = wxTB_FLAT | wxTB_DOCKABLE | wxTB_TEXT;

// Resources
#if !defined(__WXMSW__) && !defined(__WXPM__)
    #include "../res/variable.xpm"
#endif

#if USE_XPM_BITMAPS
    #include "../res/xaxis.xpm"
    #include "../res/yaxis.xpm"
    #include "../res/zaxis.xpm"
    #include "../res/free.xpm"
#endif // USE_XPM_BITMAPS

enum {
	IDT_XDOF,IDT_YDOF,IDT_ZDOF,IDT_FDOF
};

namespace FemViewer {
namespace WxGUI		{
wxDEFINE_EVENT(wxEVT_COMMAND_MYTHREAD_COMPLETED, wxThreadEvent);
wxDEFINE_EVENT(wxEVT_COMMAND_MYTHREAD_UPDATE, wxThreadEvent);

wxBEGIN_EVENT_TABLE(wxMainWindow, wxFrame)
	EVT_CHAR(wxMainWindow::OnChar)
	#ifndef _USE_FV_LIB
	EVT_MENU(IDM_FILE_OPENMESH,wxMainWindow::OnMenuFileOpen)
	EVT_MENU(IDM_FILE_OPENFIELD,wxMainWindow::OnMenuFileOpen)
	#endif
	EVT_MENU(IDM_FILE_REFRESH,wxMainWindow::OnMenuFileRefresh)
	EVT_MENU(IDM_FILE_RELOAD,wxMainWindow::OnMenuFileReload)
	EVT_MENU(IDM_FILE_CLEAR,wxMainWindow::OnMenuFileReset)
	#ifndef _USE_FV_LIB
	EVT_MENU_RANGE(wxID_FILE,wxID_FILE9, wxMainWindow::OnMenuFileHistory)
	#endif
	EVT_MENU(IDM_FILE_EXIT,wxMainWindow::OnMenuFileQuit)
	EVT_MENU(IDM_VIEW_PERSPECTIVE,wxMainWindow::OnMenuViewProjection)
	EVT_MENU(IDM_VIEW_ORTOGRAPHIC,wxMainWindow::OnMenuViewProjection)
	EVT_MENU(IDM_VIEW_TOP,wxMainWindow::OnMenuChangeView)
	EVT_MENU(IDM_VIEW_BOTTOM,wxMainWindow::OnMenuChangeView)
	EVT_MENU(IDM_VIEW_FRONT,wxMainWindow::OnMenuChangeView)
	EVT_MENU(IDM_VIEW_BACK,wxMainWindow::OnMenuChangeView)
	EVT_MENU(IDM_VIEW_LEFT,wxMainWindow::OnMenuChangeView)
	EVT_MENU(IDM_VIEW_RIGHT,wxMainWindow::OnMenuChangeView)
	EVT_MENU(IDM_VIEW_DEFAULT,wxMainWindow::OnMenuChangeView)
	EVT_MENU(IDM_VIEW_FULL,wxMainWindow::OnMenuViewMode)
	EVT_MENU(IDM_VIEW_BOUNDINGBOX,wxMainWindow::OnMenuViewMode)
	EVT_MENU(IDM_VIEW_FAST,wxMainWindow::OnMenuViewMode)
	EVT_MENU(IDM_VIEW_NEWVIEW,wxMainWindow::OnMenuViewNew)
	EVT_MENU(IDM_VIEW_NEXTVIEW,wxMainWindow::OnMenuViewNext)
	EVT_MENU(IDM_VIEW_PREVIOUSVIEW,wxMainWindow::OnMenuViewPrevious)
	EVT_MENU(IDM_VIEW_DUMPCURRENTVIEW,wxMainWindow::OnMenuViewDumpCurr)
	EVT_MENU(IDM_VIEW_DUMPALLVIEWS,wxMainWindow::OnMenuViewDumpAll)
	EVT_MENU(IDM_CONFIGURE_AXES,wxMainWindow::OnMenuConfigAxes)
	EVT_MENU(IDM_CONFIGURE_GRID,wxMainWindow::OnMenuConfigGrid)
	EVT_MENU(IDM_CONFIGURE_CUT_PLANE,wxMainWindow::OnMenuConfigCutPlane)
	EVT_MENU(IDM_CONFIGURE_BACKGROUNDCOLOR, wxMainWindow::OnMenuConfigBkgColor)
	EVT_MENU(IDM_CONFIGURE_LIGHT,wxMainWindow::OnMenuConfigLight)
	EVT_MENU(IDM_CONFIGURE_RESET,wxMainWindow::OnMenuConfigReset)
	EVT_MENU(IDM_CONFIGURE_SAVECONFIGURATION,wxMainWindow::OnMenuConfigSave)
	EVT_MENU(IDM_CONFIGURE_EDITLEGEND,wxMainWindow::OnMenuConfigLegendEdit)
	EVT_MENU(IDM_CONFIGURE_CHANGEAPROXIMATIONMODULE,wxMainWindow::OnMenuConfigModuleApprox)
	EVT_MENU(IDM_RENDER_SETSOLUTIONFORMULA,wxMainWindow::OnMenuRenderSolutionSettings)
	EVT_MENU(IDM_RENDER_WIREFRAME,wxMainWindow::OnMenuRenderDraw)
	EVT_MENU(IDM_RENDER_SHADED,wxMainWindow::OnMenuRenderDraw)
	EVT_MENU(IDM_RENDER_CONTOURLINES,wxMainWindow::OnMenuRenderDraw)
	EVT_MENU(IDM_RENDER_NUM_VERTEX,wxMainWindow::OnMenuRenderDraw)
	EVT_MENU(IDM_RENDER_NUM_ELEM,wxMainWindow::OnMenuRenderDraw)
	EVT_MENU(IDM_RENDER_GL,wxMainWindow::OnMenuRenderDraw)
	EVT_MENU(IDM_RENDER_GL_CL,wxMainWindow::OnMenuRenderDraw)

	//EVT_MENU(IDM_RDRAW_FLOODED,wxMainWindow::OnMenuRenderDraw)
	EVT_MENU(IDM_RENDER_TRIM,wxMainWindow::OnMenuRenderDraw)
	EVT_MENU(IDM_RENDER_GRID,wxMainWindow::OnMenuRenderDraw)
	EVT_MENU(IDM_RENDER_RASTR,wxMainWindow::OnMenuRenderDraw)
	EVT_MENU(IDM_RENDER_RAYS,wxMainWindow::OnMenuRenderDraw)
	EVT_MENU(IDM_RENDER_SETSECTION,wxMainWindow::OnMenuRenderCutPlaneParams)
	EVT_MENU(IDM_RENDER_DUMPSCREEN,wxMainWindow::OnMenuRenderDumpScreen)
	EVT_MENU(IDM_HELP_ABOUT,wxMainWindow::OnMenuHelpAbout)
	EVT_IDLE(wxMainWindow::OnIdle)
	EVT_CLOSE(wxMainWindow::OnCloseWindow)
	EVT_THREAD(WORKER_EVENT_UPDATE, wxMainWindow::OnThreadUpdate)
	EVT_THREAD(WORKER_EVENT_COMPLETE, wxMainWindow::OnThreadCompletion)

	EVT_UPDATE_UI(IDM_RENDER_TRIM, wxMainWindow::OnUpdateMenuItemCutPlane)
	//EVT_SIZE(wxMainWindow::OnResize)
wxEND_EVENT_TABLE()



RenderThread::~RenderThread()
{
   wxCriticalSectionLocker enter(wxGetApp().m_critsect);
   //_handler->m_pThread = NULL;
}

void* RenderThread::Entry()
{
	//mfp_log_debug("entry\n");
	//m_pGLCanvas->SeCurrent(m_pGLCanvas->)
	//wxPaintDC dc(m_pGLCanvas);
	static fv_timer t;
	static int i, nb_frame;
	//_handler->Show(true);

	//OpenGL params
//	#ifdef __WXMSW__
//	int *gl_attrib = NULL;
//	#else
//	int gl_attrib[20] = { WX_GL_RGBA, WX_GL_MIN_RED, 1, WX_GL_MIN_GREEN, 1,
//			 WX_GL_MIN_BLUE, 1, WX_GL_DEPTH_SIZE, 1,
//			 WX_GL_DOUBLEBUFFER,
//	#  if defined(__WXMAC__) || defined(__WXCOCOA__)
//			 GL_NONE };
//	#  else
//		  None };
//	#  endif
//	#endif
	//wxRenderCanvas * canvas = new wxRenderCanvas(_handler, wxID_ANY, gl_attrib, _handler->GetPosition(),
	//			_handler->GetClientSize(), wxDEFAULT_FRAME_STYLE);
	//_handler->m_canvas = canvas;
	//wxASSERT(_handler->m_canvas != NULL);
	// Init GL objects like axes grid
	// This should be placed in other place not here
	//canvas->m_modelCtr->Do(ModelCtrl::INIT); // INIT - default
	//canvas->m_viewCtr->Init(canvas->GetSize().x,canvas->GetSize().y);
	wxClientDC dc(_handler->m_canvas);
	_handler->m_canvas->SetCurrent(*(_handler->m_canvas->m_context));
	//mfp_log_debug("before loop\n");
	t.start();
	double dLastTime  = 0.0;
	while(!TestDestroy())
	{
		//wxLogError("in loop\n");


		wxMutexGuiEnter();
//		wxClientDC dc(_handler->m_canvas);
		mfp_log_debug("loop %d\n",i++);
		//_handler->m_canvas->SetCurrent(*(_handler->m_canvas->m_context));
		mfp_log_debug("Before dipplaycllback\n");
		//_handler->m_pVMngr->DisplayCallback();
		_handler->m_canvas->PaintNow();
		//_handler->m_canvas->SwapBuffers();
		mfp_log_debug("after swapbuffers\n");
		//_handler->m_canvas->SetCurrent(*(_handler->m_canvas->m_dummyconetxt));
		//m_pGLCanvas->PaintNow();
		wxMutexGuiLeave();
		//Sleep(50);
		//m_pGLCanvas->SetCurrent(*m_pGLCanvas->m_context);
		//ViewManagerInst().DisplayCallback();
		//m_pGLCanvas->SwapBuffers();
		//m_pGLCanvas->SetCurrent(*m_pGLCanvas->m_dummyconetxt);
		//wxCommandEvent event(wxEVT_ERASE_BACKGROUND, EVT_MY_REFRESH_CANVAS);
		//m_pGLCanvas->GetEventHandler()->AddPendingEvent(event);
		//
		Sleep(_iRefreshDelay);
		//wxQueueEvent(_handler,new wxThreadEvent(wxEVT_COMMAND_MYTHREAD_UPDATE));
#if 1
	nb_frame++;
	double currentTime = t.get_time_ms();
	double delta = currentTime - dLastTime;
	if (delta > 1000.0 ) {
		// If last prinf() was more than 1 sec ago
		printf("render %d %.3f ms/frame, %.1f FPS\n",nb_frame, delta / double(nb_frame),
				double(nb_frame) * delta  * 0.001);
		nb_frame = 0;
		dLastTime = currentTime;
	}

#endif
	}
	//mfp_log_debug("Leaving\n");
	wxQueueEvent(_handler,new wxThreadEvent(WORKER_EVENT_COMPLETE));
	return (void *)0;
}


wxMainWindow::wxMainWindow(const int width, const int height,
						   const wxChar* title)
: wxFrame(NULL, wxID_ANY, wxString(title))
{
	//mfp_log_debug("wxMainWindow ctr\n");
	m_pMCtrl = &ModelCtrlInst();
	m_pVMngr = &ViewManagerInst();
	m_gsettings = m_pVMngr->GetSettings();

	// Create statusbar & specify its properties
	int widths[] = { -1, 100,100, 100, 100 };
	m_status = CreateStatusBar();

	wxASSERT( m_status!=NULL );
	m_status->SetFieldsCount( WXSIZEOF(widths), widths );

	// Set icon
	SetIcon( wxICON( variable ) );
/*
	// Set visual configuration
	m_gridOn = m_pVMngr->GetSettings()->bIsGridOn;
	m_axesOn = m_pVMngr->GetSettings()->bIsAxesOn;
	m_colorMapOn = m_pVMngr->GetSettings()->bIsLegendOn;
	m_edgesOn = m_pVMngr->GetSettings()->bEdgeOn;
	m_shadingOn = m_pVMngr->GetSettings()->bShadingOn;
	m_contourLinesOn = m_pVMngr->GetSettings()->bIsovalueLineOn;
	m_numVertsOn = m_pVMngr->GetSettings()->bShowNumVertices;
	m_numElemsOn = m_pVMngr->GetSettings()->bShowNumElems;
	m_orthoProjOn = m_pVMngr->GetSettings()->bIsOrthoOn;
	m_bvhGridOn = m_pVMngr->GetSettings()->bIsBVHGridDraw;
	m_bvhOctreeOn = m_pVMngr->GetSettings()->bIsBVHOctreeDraw;
	m_rayTraceOn = m_pVMngr->GetSettings()->bIsRayTracing;
	m_useCutPlane = m_pVMngr->GetSettings()->bSliceModel;
	m_sliceModelOn = m_pVMngr->GetSettings()->bShowCutPlane[0];
*/
	#ifndef _USE_FV_LIB
	// Read configs
	wxConfigBase * pConfig = wxConfigBase::Get();
	pConfig->SetPath( wxT("FileHistory") );
	m_fhistory = new ::wxFileHistory();
	m_fhistory->Load( *pConfig );
	#endif

   	// Menu File
	m_filemenu = new wxMenu;
	#ifndef _USE_FV_LIB
	m_filemenu->Append(IDM_FILE_OPENMESH,wxT("Open &mesh file...\tCtrl+M"));
	m_filemenu->Append(IDM_FILE_OPENFIELD,wxT("Open &field file...\tCtrl+F"));
	m_filemenu->AppendSeparator();
	#endif
	m_filemenu->Append(IDM_FILE_REFRESH,wxT("Refresh field data"));
	m_filemenu->Append(IDM_FILE_RELOAD, wxT("Reload mesh and field data"));
	m_filemenu->Append(IDM_FILE_CLEAR,wxT("&Reset model"));
	#ifndef _USE_FV_LIB
	m_filemenu->AppendSeparator();
	m_fhistory->AddFilesToMenu(m_filemenu);
	#endif
	m_filemenu->AppendSeparator();
	m_filemenu->Append(IDM_FILE_EXIT,wxT("&Quit"));

	/// Menu View
	m_viewmenu = new wxMenu;
	m_viewmenu->AppendRadioItem(IDM_VIEW_PERSPECTIVE,wxT("Perspective"));
	m_viewmenu->AppendRadioItem(IDM_VIEW_ORTOGRAPHIC,wxT("Parallel"));
	m_viewmenu->AppendSeparator();
	m_viewmenu->AppendRadioItem(IDM_VIEW_TOP,wxT("&Top"));
	m_viewmenu->AppendRadioItem(IDM_VIEW_BOTTOM,wxT("Botto&m"));
	m_viewmenu->AppendRadioItem(IDM_VIEW_FRONT,wxT("&Front"));
	m_viewmenu->AppendRadioItem(IDM_VIEW_BACK,wxT("&Back"));
	m_viewmenu->AppendRadioItem(IDM_VIEW_LEFT,wxT("&Left"));
	m_viewmenu->AppendRadioItem(IDM_VIEW_RIGHT,wxT("Rig&ht"));
	m_viewmenu->AppendRadioItem(IDM_VIEW_DEFAULT,wxT("Defa&ult"));
	m_viewmenu->AppendSeparator();
	m_viewmenu->AppendRadioItem(IDM_VIEW_FULL,wxT("Full mode"));
	m_viewmenu->AppendRadioItem(IDM_VIEW_BOUNDINGBOX,wxT("Bounding box mode"));
	m_viewmenu->AppendRadioItem(IDM_VIEW_FAST,wxT("Fast mode"));
	m_viewmenu->AppendSeparator();
	m_viewmenu->Append(IDM_VIEW_NEWVIEW,wxT("New add"));
	m_viewmenu->Append(IDM_VIEW_NEXTVIEW,wxT("Next"));
	m_viewmenu->Append(IDM_VIEW_PREVIOUSVIEW,wxT("Previous"));
	m_viewmenu->AppendSeparator();
	m_viewmenu->Append(IDM_VIEW_DUMPCURRENTVIEW,wxT("Dump current"));
	m_viewmenu->Append(IDM_VIEW_DUMPALLVIEWS,wxT("Dump all"));

	// Menu Configuration
	m_cfgmenu = new wxMenu;
	m_cfgmenu->Append(IDM_CONFIGURE_AXES,wxT("Axes"),wxT("Displaying axexs"),true);
	m_cfgmenu->Append(IDM_CONFIGURE_GRID,wxT("Grid"),wxT("Displaying grid"), true);
	m_cfgmenu->Append(IDM_CONFIGURE_CUT_PLANE,wxT("Cut-plane"),wxT("On/Off display of cut-planes"),true);
	m_cfgmenu->Append(IDM_CONFIGURE_BACKGROUNDCOLOR,wxT("Background color"));
	m_cfgmenu->Append(IDM_CONFIGURE_LIGHT,wxT("Lighting"));
	m_cfgmenu->AppendSeparator();
	m_cfgmenu->Append(IDM_CONFIGURE_RESET,wxT("Defaults"));
	m_cfgmenu->Append(IDM_CONFIGURE_SAVECONFIGURATION,wxT("Save"));
	m_cfgmenu->AppendSeparator();
	m_cfgmenu->Append(IDM_CONFIGURE_EDITLEGEND,wxT("Edit legend..."));
	m_cfgmenu->AppendSeparator();
	m_cfgmenu->Append(IDM_CONFIGURE_CHANGEAPROXIMATIONMODULE,wxT("Use external approximation..."));

	/// Menu Render
	m_rendermenu = new wxMenu;
	m_rendermenu->Append(IDM_RENDER_SETSOLUTIONFORMULA,wxT("Solution settings..."));
	m_rendermenu->AppendSeparator();
	m_rendermenu->Append(IDM_RENDER_WIREFRAME,wxT("Draw element edges"),wxT("Turn On/Off edges of a mesh on the screen"),true);
	m_rendermenu->Append(IDM_RENDER_SHADED,wxT("Draw field"),wxT("Turn On/Off a solid solution on the screen"), true);
	m_rendermenu->Append(IDM_RENDER_CONTOURLINES,wxT("Draw contours"),wxT("Turn On/Off displaying of contours"), true);
	m_rendermenu->AppendSeparator();
	m_rendermenu->Append(IDM_RENDER_NUM_VERTEX,wxT("Draw vertex id"),wxT("Turn On/Off displaying of id of vertices"), true);
	m_rendermenu->Append(IDM_RENDER_NUM_ELEM,wxT("Draw element id"),wxT("Turn On/Off displaying of id of elements"), true);
	m_rendermenu->AppendSeparator();
	m_rendermenu->AppendCheckItem(IDM_RENDER_TRIM, wxT("Show cutted"),wxT("Turn On/Off slice of a model"));
	m_rendermenu->Append(IDM_RENDER_SETSECTION, wxT("Cut plane settings..."),wxT("Enter to cutting settings"));
	m_rendermenu->AppendSeparator();
	m_rendermenu->Append(IDM_RENDER_GRID, wxT("Draw grid"), wxT("Show grid acceleration structure"), true);
	m_rendermenu->AppendSeparator();
	m_rendermenu->AppendRadioItem(IDM_RENDER_GL, wxT("Use rasterization"), wxT("Switch to OpenGL rasterization"));
	m_rendermenu->AppendRadioItem(IDM_RENDER_GL_CL, wxT("Use tracing rays"), wxT("Switch to Ray Tracing with OpenCL"));
	m_rendermenu->AppendSeparator();
	m_rendermenu->Append(IDM_RENDER_DUMPSCREEN,wxT("Save screen..."));
	m_rendermenu->AppendSeparator();

	// Menu Help
	m_helpmenu = new wxMenu;
	m_helpmenu->Append(IDM_HELP_ABOUT,wxT("&About..."));

	// Set menu toolbar
	m_menubar = new wxMenuBar(wxMB_DOCKABLE);
	m_menubar->Append(m_filemenu,	wxT("&File"));
	m_menubar->Append(m_viewmenu,	wxT("&View"));
	m_menubar->Append(m_cfgmenu, 	wxT("&Configuration"));
	m_menubar->Append(m_rendermenu, wxT("&Render"));
	m_menubar->Append(m_helpmenu, 	wxT("&Help"));

	// Check appropriate menu items
	m_menubar->Check(IDM_CONFIGURE_AXES, m_gsettings->bIsAxesOn);
	m_menubar->Check(IDM_CONFIGURE_GRID, m_gsettings->bIsGridOn);
	m_menubar->Check(IDM_CONFIGURE_CUT_PLANE, m_gsettings->bShowCutPlane[0]);
	m_menubar->Check(IDM_RENDER_WIREFRAME, m_gsettings->bEdgeOn);
	m_menubar->Check(IDM_RENDER_SHADED, m_gsettings->bShadingOn);
	m_menubar->Check(IDM_RENDER_CONTOURLINES, m_gsettings->bIsovalueLineOn);
	m_menubar->Check(IDM_RENDER_NUM_VERTEX, m_gsettings->bShowNumVertices);
	m_menubar->Check(IDM_RENDER_NUM_ELEM, m_gsettings->bShowNumElems);
	m_menubar->Check(IDM_RENDER_TRIM, m_gsettings->bShowCutPlane[0]);
	m_menubar->Check(IDM_RENDER_GRID, m_gsettings->bIsBVHGridDraw);

	SetMenuBar(m_menubar);
	//m_menubar->Check(IDM_RENDER_RAYS, m_pVMngr->GetSettings()->bIsRayTracing);
	// Create toolbar
	//CreateWndToolBar();
	// Resize window to required sizes
	SetClientSize(wxSize(width,height));
	CenterOnScreen();
	Show(true);
	Raise();
	//DoInit();
	//OpenGL params
	int attribs[] = {WX_GL_RGBA, WX_GL_DOUBLEBUFFER, WX_GL_DEPTH_SIZE, 16, 0};

	m_canvas = new wxRenderCanvas(this, wxID_ANY, attribs, GetPosition(),
			GetClientSize(), wxDEFAULT_FRAME_STYLE);

	if (!m_canvas) throw fv_exception("Can't initialize OpenGL machine!");
	//DoInit();
	//Update();
	//FV_CHECK_ERROR_GL();

	//StartRenderThread();
	//SetProperties();
	m_status->SetStatusText( wxT("Graphic module initialized"), 0 );
}

wxMainWindow::~wxMainWindow()
{
	//mfp_log_debug("Dtr +++++++++++++++++++++++++++++++\n");

	// Destroy OpenGL machine
	wxDELETE( m_canvas );

    /*if( m_searchTool && !m_searchTool->GetToolBar() )
    {
        GetToolBar()->AddTool(m_searchTool);
    }*/

    #ifndef _USE_FV_LIB
	wxConfigBase * pConfig = wxConfigBase::Get();
	if( pConfig != NULL )
	{
		// TDO: write configs
		m_fhistory->Save(*pConfig);
		wxDELETE(m_fhistory)
	}
	#endif

	wxDELETE( m_status );
	//mfp_debug("after delet status\n");
	//wxWindow::Destroy();
}


void wxMainWindow::DoInit(int cmd)
{
	//mfp_log_debug("In wxMainWindow::DoIni\n");
	//mfp_debug("wxMainWindow: DoInit\n");
	#ifndef _USE_FV_LIB
	wxCommandEvent event;
	OnMenuConfigModuleApprox(event);
	#else
	try {
		// First, read mesh data
		m_pMCtrl->Do((ModelCtrl::eCommnad)cmd);
		FV_CHECK_ERROR_GL();
		UpdateStatusBar();
	} catch(const fv_exception& ex) {
		std::cout << ex << std::endl;
	} catch(...)
	{
		std::cout<< "Something is wrong!\n";
	}
	#endif
}

void wxMainWindow::Update()
{
	int nno = ModelCtrlInst().GetCurrentMesh()->GetNumNodes();
	int ned = ModelCtrlInst().GetCurrentMesh()->GetNumEdges();
	int nfa = ModelCtrlInst().GetCurrentMesh()->GetNumFaces();
	int nel = ModelCtrlInst().GetCurrentMesh()->GetNumElems();
	m_status->SetStatusText( wxString::Format("no: %d", nno), 1);
	m_status->SetStatusText( wxString::Format("ed: %d", ned), 2);
	m_status->SetStatusText( wxString::Format("fa: %d", nfa), 3);
	m_status->SetStatusText( wxString::Format("el: %d", nel), 4);

}

void wxMainWindow::StartRenderThread()
{
//	m_pThread = new RenderThread(this);
//	wxASSERT(m_pThread != NULL);
//	if (m_pThread->Create() != wxTHREAD_NO_ERROR  ||
//		m_pThread->Run() != wxTHREAD_NO_ERROR)
//	{
//		wxLogError(wxT("Can't create thread!"));
//		delete m_pThread; m_pThread = NULL;
//	}l501

}


void wxMainWindow::OnThreadUpdate(wxThreadEvent& WXUNUSED(event))
{
	//wxCriticalSectionLocker enter(m_pThreadCS);
	//mfp_log_debug("Update from thread\n");
	m_canvas->Refresh();
	//wxMessageOutputDebug().Printf("MYFRAME: MyThread update...\n");
}

void wxMainWindow::OnThreadCompletion(wxThreadEvent& WXUNUSED(event))
{
	//mfp_log_debug("End thread\n");
	wxMessageOutputDebug().Printf("MYFRAME: MyThread exited!\n");
}




void wxMainWindow::DoPauseThread()
{
	// anytime we access the m_pThread pointer we must ensure that it won't
	// be modified in the meanwhile; since only a single thread may be
	// inside a given critical section at a given time, the following code
	// is safe:
//	wxCriticalSectionLocker enter(m_pThreadCS);
//	if (m_pThread) // does the thread still exist?
//	{
//		// without a critical section, once reached this point it may happen
//		// that the OS scheduler gives control to the MyThread::Entry() function,
//		// which in turn may return (because it completes its work) making
//		// invalid the m_pThread pointer
//		if (m_pThread->Pause() != wxTHREAD_NO_ERROR )
//			wxLogError("Can't pause the thread!");
//	}
}

//void wxMainWindow::RefreshScreen(void) 
//{
//	wxASSERT( m_canvas!=NULL );
//	m_canvas->wxWindow::Refresh(false);InitGLVisual
//}

///----Private methods --------------------------------------------
void wxMainWindow::OnChar(wxKeyEvent& event)
{
    #if wxUSE_UNICODE
    int key = event.GetUnicodeKey();
    #else
    int key = event.GetKeyCode();
    #endif

    switch( key )
    {
        case WXK_ESCAPE:
            wxTheApp->ExitMainLoop();
            return;

        default:
			//getViewManager().KeyboardCallback(
            event.Skip();
            return;
    }

    this->Refresh(false);
}


// Menu
#ifndef _USE_FV_LIB
void wxMainWindow::OnMenuFileOpen(wxCommandEvent& event)
{
	wxMessageDialog *dlg(NULL);
	wxString s_extdef, wxErrMsg;
	wxString title(_T("Open "));
	title += (event.GetId() == IDM_OPEN_MESH) ? wxString(_T("mesh file")) : wxString(_T("field file"));
	wxString wxPath = wxFileSelector(title,
								   wxEmptyString,
								   wxEmptyString,
								   wxString::Format
			                       (
			                    		   _T("Input file (*.dat;*dmp)|*.dat;*.dmp|All files (%s)|%s"),
			                               wxFileSelectorDefaultWildcardStr,
			                               wxFileSelectorDefaultWildcardStr
			                       ),
			                       wxFileSelectorDefaultWildcardStr,
			                       wxCHANGE_DIR,
			                       this);
	if (!wxPath)   return;
	else
    {
		try
		{
			#ifndef __WXDEBUG__
			wxBusyInfo infoRead(_T("Reading file, please wait..."), this);
			#endif
			//_proj.Init(path.mb_str(), file_type);
			_fhistory->AddFileToHistory(wxPath);
			_fhistory->AddFilesToMenu();
			m_canvas->Refresh(true);
		}
		catch (fv_exception& ex)
		{
			std::ostringstream oss;
			oss << ex;
			wxString reson;
			reson.Printf(_T("%s"),oss.str() );
			(void)wxMessageBox(reson,_T("Error"),wxICON_ERROR | wxOK);

		}
	}

	try 
	{
		//GetFilePath(wxPath,event.GetId());

	}
	catch(std::runtime_error & Ex)
	{
		wxErrMsg.Printf(_T("%s"), Ex.what());
		dlg = new wxMessageDialog(NULL, wxErrMsg, _T("Error"), wxOK | wxICON_ERROR);
	}
	catch(...)
	{
		wxErrMsg = _T("Error accrued while loading file ");
		wxErrMsg.append(wxPath);
		dlg = new wxMessageDialog(NULL, wxErrMsg, _T("Error"), wxOK | wxICON_ERROR);
    }

	if(dlg) {
		dlg->ShowModal();
		dlg->Destroy();
	}
	
}
#endif

void wxMainWindow::OnMenuFileReset(wxCommandEvent& WXUNUSED(event))
{
	DoInit(ModelCtrl::CLEAR);
	wxFrame::SetTitle(wxT("FemViewer"));
	Update();
	m_canvas->Refresh(true);

}

/*void wxMainWindow::OnMenuFileDump(wxCommandEvent& WXUNUSED(event))
{
	if(m_pVMngr->MenuCallback("general_dump_characteristics")){
		m_canvas->Refresh(false);
	}
}*/



void wxMainWindow::OnMenuFileHistory(wxCommandEvent& event)
{
//	const int idFile = event.GetId() - wxID_FILE;
//	if(idFile < MAX_HISTORY_FILES)
//	{
//		wxString wxPath = _fhistory->GetHistoryFile(idFile);
//		_proj.setFileName(wxPath.mb_str());
//	} else {
//		event.Skip();
//	}

	m_canvas->Refresh(false);
}


// This function refresh field data
void wxMainWindow::OnMenuFileRefresh(wxCommandEvent& WXUNUSED(event))
{
	//m_pVMngr->MenuCallback(IDM_FILE_REFRESH);
	this->m_pMCtrl->Do(ModelCtrl::UPDATE | ModelCtrl::INIT_FIELD);
	this->Update();
	m_canvas->Refresh(false);
}

void wxMainWindow::OnMenuFileReload(wxCommandEvent& WXUNUSED(event))
{
	//m_pVMngr->MenuCallback(IDM_FILE_REFRESH);
	mfp_debug("reseive signal to update");
	m_pMCtrl->Do(ModelCtrl::CLEAR);
	m_pMCtrl->Do(ModelCtrl::INIT);
	this->Update();
	m_canvas->Refresh(false);
	wxGetApp().m_parent->Notify();
}

void wxMainWindow::OnMenuFileQuit(wxCommandEvent& WXUNUSED(event))
{
	//mfp_log_debug("wxMainWindow: Quiting!\n");
	// This will be chaught in CLOSE_EVT
	//m_pMCtrl->Clear();
	//m_canvas->CleanUp();
	Close(true);
}

void wxMainWindow::OnMenuViewProjection(wxCommandEvent& event)
{
	if(event.GetId() == IDM_VIEW_PERSPECTIVE)
		m_pVMngr->MenuCallback(IDM_VIEW_PERSPECTIVE);
	else if(event.GetId() == IDM_VIEW_ORTOGRAPHIC)
		m_pVMngr->MenuCallback(IDM_VIEW_ORTOGRAPHIC);
	m_canvas->Refresh(false);

}

void wxMainWindow::OnMenuChangeView(wxCommandEvent& event)
{
	int evt_id = event.GetId();
	switch(evt_id)
	{
	case IDM_VIEW_TOP:
	case IDM_VIEW_BOTTOM:
	case IDM_VIEW_FRONT:
	case IDM_VIEW_BACK:
	case IDM_VIEW_LEFT:
	case IDM_VIEW_RIGHT:
	case IDM_VIEW_DEFAULT:
	case IDM_VIEW_FULL:
	case IDM_VIEW_FAST:
		if(m_pVMngr->MenuCallback(evt_id))
				m_canvas->Refresh(false);
		break;
	

	default:
		break;
	}

}

void wxMainWindow::OnMenuViewMode(wxCommandEvent& event)
{
	const int evt_id = event.GetId();
	switch(evt_id)
	{
	case IDM_VIEW_FULL:
	case IDM_VIEW_BOUNDINGBOX:
	case IDM_VIEW_FAST:
		if(m_pVMngr->MenuCallback(evt_id))
				m_canvas->Refresh(false);
		return;
	default:
		return;
	}

}

void wxMainWindow::OnMenuViewNew(wxCommandEvent& event)
{
	if(m_pVMngr->MenuCallback(IDM_VIEW_NEWVIEW)) {
		m_canvas->Refresh(false);
	}
}

void wxMainWindow::OnMenuViewNext(wxCommandEvent& event)
{
	if(m_pVMngr->MenuCallback(IDM_VIEW_NEXTVIEW)) {
		m_canvas->Refresh(false);
	}
}

void wxMainWindow::OnMenuViewPrevious(wxCommandEvent& event)
{
	if(m_pVMngr->MenuCallback(IDM_VIEW_PREVIOUSVIEW)) {
		m_canvas->Refresh(false);
	}
}

void wxMainWindow::OnMenuViewDumpCurr(wxCommandEvent& event)
{
	if(m_pVMngr->MenuCallback(IDM_VIEW_DUMPCURRENTVIEW)) {
		m_canvas->Refresh(false);
	}
}

void wxMainWindow::OnMenuViewDumpAll(wxCommandEvent& event)
{
	if(m_pVMngr->MenuCallback(IDM_VIEW_DUMPALLVIEWS)) {
		m_canvas->Refresh(false);
	}
}

void wxMainWindow::OnMenuConfigAxes(wxCommandEvent& event)
{
	if(m_pVMngr->MenuCallback(IDM_CONFIGURE_AXES)) {
		m_canvas->Refresh(false);
	}
}

void wxMainWindow::OnMenuConfigGrid(wxCommandEvent& event)
{
	if(m_pVMngr->MenuCallback(IDM_CONFIGURE_GRID)) {
		m_canvas->Refresh(false);
	}
}

void wxMainWindow::OnMenuConfigCutPlane(wxCommandEvent& event)
{
	m_gsettings->bShowCutPlane[0] = event.IsChecked();
	if(m_pVMngr->MenuCallback(IDM_CONFIGURE_CUT_PLANE)) {
		m_canvas->Refresh(false);
	}
}

void wxMainWindow::OnMenuConfigBkgColor(wxCommandEvent& event)
{
	if(m_pVMngr->MenuCallback(IDM_CONFIGURE_BACKGROUNDCOLOR)) {
		m_canvas->Refresh(true);
	}
}

void wxMainWindow::OnMenuConfigLight(wxCommandEvent& event)
{
	if(m_pVMngr->MenuCallback(IDM_CONFIGURE_LIGHT)) {
		m_canvas->Refresh(false);
	}
}

void wxMainWindow::OnMenuConfigReset(wxCommandEvent& event)
{
	if(m_pVMngr->MenuCallback(IDM_VIEW_DEFAULT)) {
		m_canvas->Refresh(false);
	}
}

void wxMainWindow::OnMenuConfigSave(wxCommandEvent& event)
{
	if(m_pVMngr->MenuCallback(IDM_CONFIGURE_SAVECONFIGURATION)) {
		m_canvas->Refresh(true);
	}
}

void wxMainWindow::OnMenuConfigLegendEdit(wxCommandEvent& event)
{
	// Temporary store current legend configuration
	Legend alegend(m_pVMngr->GetLegend());
	wxLegendDialog *dlg = new wxLegendDialog(wxString(wxT("Edit legend")),alegend);

    if (dlg->ShowModal() == wxID_OK) {
    	// Now model points to changed legend
    	m_pMCtrl->Do(ModelCtrl::UPDATE | ModelCtrl::INIT_LEGEND);
    	m_pMCtrl->Do();
    }
    dlg->Destroy();
	m_canvas->Refresh(false);
}

void wxMainWindow::OnMenuConfigModuleApprox(wxCommandEvent & WXUNUSED(event))
{
//	mod_approx& mod = m_Project.GetModuleType();
//	wxApproxModWnd* dlg = new wxApproxModWnd(this,mod);
//	dlg->CenterOnParent();
//	if(wxID_OK == dlg->ShowModal())
//	{
//		m_Project.set_module();
//	}
//
//	dlg->Destroy();
	m_canvas->Refresh(false);
}

void wxMainWindow::OnMenuRenderSolutionSettings(wxCommandEvent& WXUNUSED(event))
{
	SolutionData lcurrSol( m_pMCtrl->GetCurrentField()->GetSolution() );
	//int nr_sols = lcurSolModelCtrlInst()->GetCurrentField()->GetNrSol();
	//int nr_equs = ModelCtrlInst()->GetCurrentField()->GetNreq();
	//int cur_sol = ModelCtrlInst()->GetCurrentField()->GetCurrSol();
	//const std::string& formula = ModelCtrlInst()->GetCurrentField()->GetFromula();

	//wxSolCalcDialog* dlg = new wxSolCalcDialog(nr_equs,nr_sols,cur_sol,formula);
	wxSolCalcDialog* dlg = new wxSolCalcDialog(lcurrSol.nrEquetions,lcurrSol.nrSolutions,lcurrSol.currSolution,lcurrSol.sFormula);
	if(dlg->ShowModal() == wxID_OK)
	{
		lcurrSol.currSolution = dlg->GetCurrSol();
		lcurrSol.sFormula     = dlg->GetFunction();
		if (m_pMCtrl->GetCurrentSolution() = lcurrSol) {
			//mfp_log_debug("Solution Changed\n");
			//ModelCtrlInst()->GetCurrentField()->GetCurrSol() = dlg->GetCurrSol();
			//ModelCtrlInst()->GetCurrentField()->GetFromula() = dlg->GetFunction();

			#ifndef __WXDEBUG__
			wxBusyInfo infoCalc(_T("Calculation, please wait..."), this);
			#endif
			this->DoInit(ModelCtrl::UPDATE | ModelCtrl::INIT_SOL);
			m_canvas->Refresh(true);
		}
	}

	dlg->Destroy();
}


void wxMainWindow::OnMenuRenderDraw(wxCommandEvent & event)
{
	const bool check = event.IsChecked() ? true : false;
	const unsigned int key(event.GetId());

	int cmd = ModelCtrl::UPDATE;
	assert(m_pVMngr);
	switch(key)
	{
		case IDM_RENDER_WIREFRAME: // Draw mesh edges
			m_gsettings->bEdgeOn = check;
			cmd |= ModelCtrl::INIT_GL;
			break;
		case IDM_RENDER_SHADED: // Draw field
			m_gsettings->bShadingOn = check;
			cmd |= ModelCtrl::INIT_GL;
			break;
		case IDM_RENDER_CONTOURLINES: // Draw contours
			m_gsettings->bIsovalueLineOn = check;
			cmd |= ModelCtrl::INIT_GL;
			break;
		case IDM_RENDER_NUM_VERTEX:
			m_gsettings->bShowNumVertices = check;
			//mfp_log_debug("Setting show num vertices to %d\n",m_gsettings->bShowNumVertices);
			cmd |= ModelCtrl::INIT_DISPLAY;
			break;
		case IDM_RENDER_NUM_ELEM:
			m_gsettings->bShowNumElems = check;
			//mfp_log_debug("Setting show num elements to %d\n",m_gsettings->bShowNumElems);
			cmd |= ModelCtrl::INIT_DISPLAY;
			break;
		case IDM_RENDER_GL:
			m_gsettings->eRenderType = RASTERIZATION_GL;
			break;
		case IDM_RENDER_GL_CL:
			m_gsettings->eRenderType = RAYTRACE_GL_CL;
			break;
		//case IDM_RDRAW_FLOODED: // Draw colored contours
		//	ViewManagerInst()GetSettings()->bLineCoolored = check;
		//	return;
		case IDM_RENDER_TRIM: // Draw cutted
			//mfp_debug("This is happedmnd\n");
			m_gsettings->bSliceModel = check;
			m_gsettings->bShowCutPlane[0] = check;
			m_gsettings->bDisplaySlices = check;
			//m_pMCtrl->SetCutPlaneChange();
			cmd |= ModelCtrl::INIT_PLANE;
			break;
		case IDM_RENDER_GRID:
			m_pVMngr->GetSettings()->bIsBVHGridDraw = check;
			cmd |= ModelCtrl::INIT_ACCEL;
			break;
		//case IDM_RENDER_RASTR:
		//case IDM_RENDER_RAYS:
		//	m_pVMngr->GetSettings()->bIsRayTracing = check;
		//	break;
		default:
			return;
	}
	
	m_pMCtrl->Do(cmd);
	m_canvas->Refresh(true);

}



void wxMainWindow::OnMenuRenderDumpScreen(wxCommandEvent &event)
{
	wxFileDialog dialog(this,
		_T("Specyfy a file name to save"),
		wxEmptyString,
		_T("screen.png"),
        _T("PNG files (*.png)|*.png|All files (*.*)|*.*"),
		wxFD_SAVE|wxFD_OVERWRITE_PROMPT);
	dialog.SetFilterIndex(0);
	if (dialog.ShowModal() == wxID_OK)
    {
        //wxLogMessage(_T("%s, filter %d"),
          //           dialog.GetPath().c_str(), dialog.GetFilterIndex());
		ScreenShot( dialog.GetPath());
    }
	m_pVMngr->AddStatusMsg("Saving screen");
	m_canvas->Refresh(true);

}

void wxMainWindow::OnMenuRenderSettings(wxCommandEvent &event)
{
	wxScreenSetts dlg(ViewManagerInst(), this);
    dlg.ShowModal();
}

void wxMainWindow::OnMenuRenderCutPlaneParams(wxCommandEvent& event)
{
	// Local vector of cut-planes
	std::vector<CutPlane> vCutPlanes(m_pMCtrl->GetCutPlanes());

	wxLogDebug("Before reading: use: %d display: %d",m_gsettings->bSliceModel,m_gsettings->bShowCutPlane[0]);
	wxPlaneCutDialog * dlg = new wxPlaneCutDialog(this,
			m_pMCtrl->GetCurrentMesh()->GetMeshBBox3D(),
			vCutPlanes,
			m_gsettings->bSliceModel,
			m_gsettings->bShowCutPlane[0]);
	int oper = 0;
    if (dlg->ShowModal() == wxID_OK)
	{
    	m_gsettings->bSliceModel = dlg->GetUseAll();
    	m_gsettings->bShowCutPlane[0] = dlg->GetDisplayAll();
    	wxLogDebug("Ok is presed with %d %d",m_gsettings->bSliceModel,m_gsettings->bShowCutPlane[0]);

    	if (! dlg->GetCutPlanes().empty())
    	{
    		wxLogDebug("After reading: use: %d display: %d",m_gsettings->bSliceModel,m_gsettings->bShowCutPlane[0]);

    		std::vector<CutPlane>& dlg_pls = dlg->GetCutPlanes();
    		std::vector<CutPlane>::iterator it;

    		// First check old data if there are changes
    		// If not, remove them
    		for (CutPlane & cp : vCutPlanes)
    		{
    			it = std::find_if(dlg_pls.begin(),
    				          dlg_pls.end(),
    				          [cp](CutPlane & pl) { return pl == cp; });
    			if (it != dlg_pls.end()) {
    				// Check for visibility
    				cp.Show(it->IsVisible());
    			}
    			else {
    				// Mark for delete
    				cp.SetChanged(false);
    			}
    		}
    		// Erase unused items
			it = std::remove_if(vCutPlanes.begin(),vCutPlanes.end(),[](CutPlane & cp) { return (cp.IsPlaneChanged() == false); });
			vCutPlanes.erase(it,vCutPlanes.end());

			// Append new
			for (CutPlane cp : dlg_pls)
			{
				it = std::find_if(dlg_pls.begin(),
								  dlg_pls.end(),
								  [cp](CutPlane & pl) { return pl == cp; });
				if (it != vCutPlanes.end())
					vCutPlanes.push_back(cp);
			}

			oper = ModelCtrl::INIT_PLANE;
			m_pMCtrl->SetCutPlaneChange();
    	}
    	else {
    		vCutPlanes.clear();
    	}
	}
    dlg->Destroy();
    m_pMCtrl->Do(ModelCtrl::UPDATE | oper);
}

void wxMainWindow::OnUpdateMenuItemCutPlane(wxUpdateUIEvent& WXUNUSED(event))
{
	//mfp_log_debug("diplay = %d\n",m_displayCutPlane);

	wxMenuItem *item = FindItemInMenuBar(IDM_RENDER_TRIM);
	if (item && item->IsCheckable())
	{
		wxLogDebug("OnupdateMenuItem state: %d",m_gsettings->bShowCutPlane[0]);
		item->Check(m_gsettings->bShowCutPlane[0]);
	}
}

void wxMainWindow::OnMenuHelpAbout(wxCommandEvent& WXUNUSED(event))
{
	wxMessageBox(wxT("Paweł Macioł - Finite Element Viewer"));
}

void wxMainWindow::OnToolLeftClick(wxCommandEvent & event)
{
	//static const char * cmds[] = {"x_dof","y_dof","z_dof","f_dof"};
	int cmd = 0;
	//const int id = event.GetId();
	if (event.GetId() == IDT_XDOF)
	{
		cmd = IDM_VIEW_RIGHT;
	} else if (event.GetId() == IDT_YDOF)
	{
		cmd = IDM_VIEW_TOP;
	} else if (event.GetId() == IDT_ZDOF)
	{
		cmd = IDM_VIEW_FRONT;
	} else if (event.GetId() == IDT_FDOF)
	{
		cmd = IDM_VIEW_DEFAULT;
	} else {
		event.Skip();
		return;
	}

	m_pVMngr->MenuCallback(cmd);
	
}


void wxMainWindow::ScreenShot(const wxString& path)
{
	int w, h;
	m_canvas->GetClientSize(&w, &h);

	const size_t imgDepth = 3;
	const size_t imgSize = imgDepth * size_t(w) * size_t(h);

	typedef unsigned char byte;

	byte* pixbuf = static_cast<byte*>(malloc(imgSize));

	int tl_x, tl_y;

	//m_canvas->DoGetPosition(&tl_x, &tl_y);
	m_canvas->GetPosition(&tl_x, &tl_y);

	const int bl_y = GetClientSize().GetHeight() - (tl_y + h);

	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadBuffer(GL_BACK_LEFT);
	glReadPixels(tl_x, bl_y, w, h, GL_RGB, GL_UNSIGNED_BYTE, pixbuf);

	for(int y = 0; y < h / 2; y++)
	{
		const int swapY = h - y - 1;
		for(int x = 0; x < w; x++)
		{
			const int offset = int(imgDepth) * (x + y * w);
			const int swapOffset = int(imgDepth) * (x + swapY * w);

			// Swap R, G and B of the 2 pixels
			std::swap(pixbuf[offset + 0], pixbuf[swapOffset + 0]);
			std::swap(pixbuf[offset + 1], pixbuf[swapOffset + 1]);
			std::swap(pixbuf[offset + 2], pixbuf[swapOffset + 2]);
		}
	}

	wxImage img(w, h, pixbuf);
	img.SaveFile(path,wxBITMAP_TYPE_PNG );

	//if(pixbuf) free(pixbuf);
	mfp_log_info("Screen sacved succesfuly\n");

}

void wxMainWindow::OnIdle(wxIdleEvent &event)
{
	//static int i;
	//mfp_log_debug("OnIdlw in wxamin %d\n",i++);
    if (wxGetApp().m_shutdown) Close(true);
    //m_canvas->PaintNow();
    event.RequestMore();
}

void wxMainWindow::OnCloseWindow(wxCloseEvent& event)
{
	/*{
		wxCriticalSectionLocker enter(m_pThreadCS);
		if (m_pThread) // does the thread still exist?
		{
			wxMessageOutputDebug().Printf("wxMainWindow: deleting thread");
			if (m_pThread->Delete() != wxTHREAD_NO_ERROR )
				wxLogError("Can't delete the thread!");
		}
	}

	while (1)
	{
		{ // was the ~MyThread() function executed?
			wxCriticalSectionLocker enter(m_pThreadCS);
			if (!m_pThread) break;
		}
		// wait for thread completion
		wxThread::This()->Sleep(1);
	}*/
	this->m_canvas->CleanUp();
	if (event.CanVeto()) {
	//	mfp_debug("pppppppppppppppppppppppppppppppppppppp\n");
		this->m_canvas->CleanUp();
		//Show(false);
		//event.Veto();
	}
	this->Destroy();
}

//void wxMainWindow::OnResizeWindow(wxSizeEvent& event)
//{
//	ViewManagerInst()ReshapeCallback();
//}
/*
wxImage wxMainWindow::DumpImage()
{
	const RenderingOptions& renderingOptions = m_appData.GetRenderingOptions();

	const int width = renderingOptions.OutputWidth();
	const int height = renderingOptions.OutputHeight();
	
	const size_t bytesPerPixel = 3;	// RGB
	const size_t imageSizeInBytes = bytesPerPixel * size_t(width) * size_t(height);
	
	// Allocate with malloc, because the data will be managed by wxImage
	BYTE* pixels = static_cast<BYTE*>(malloc(imageSizeInBytes));

	// glReadPixels takes the lower-left corner, while GetViewportOffset gets the top left corner
	const wxPoint topLeft = GetViewportOffset();
	const wxPoint lowerLeft(topLeft.x, GetClientSize().GetHeight() - (topLeft.y + height));

	// glReadPixels can align the first pixel in each row at 1-, 2-, 4- and 8-byte boundaries. We
	// have allocated the exact size needed for the image so we have to use 1-byte alignment
	// (otherwise glReadPixels would write out of bounds)
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadPixels(lowerLeft.x, topLeft.y, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);

	// glReadPixels reads the given rectangle from bottom-left to top-right, so we must
	// reverse it
	for(int y = 0; y < height / 2; y++)
	{
		const int swapY = height - y - 1;
		for(int x = 0; x < width; x++)
		{
			const int offset = int(bytesPerPixel) * (x + y * width);
			const int swapOffset = int(bytesPerPixel) * (x + swapY * width);

			// Swap R, G and B of the 2 pixels
			std::swap(pixels[offset + 0], pixels[swapOffset + 0]);
			std::swap(pixels[offset + 1], pixels[swapOffset + 1]);
			std::swap(pixels[offset + 2], pixels[swapOffset + 2]);
		}
	}
	
	return wxImage(width, height, pixels);
}
*/

/*void wxMainWindow::SetProperties()
{
	m_menubar->Check(IDM_CONFIGURE_AXES, m_axesOn);
	m_menubar->Check(IDM_CONFIGURE_GRID, m_gridOn);
	m_menubar->Check(IDM_RENDER_WIREFRAME, m_edgesOn);
	m_menubar->Check(IDM_RENDER_SHADED, m_shadingOn);
	m_menubar->Check(IDM_RENDER_CONTOURLINES, m_contourLinesOn);
	//mfp_debug("Before checking 2\n");
	m_menubar->Check(IDM_RENDER_NUM_VERTEX, m_numVertsOn);
	m_menubar->Check(IDM_RENDER_NUM_ELEM, m_numElemsOn);
	m_menubar->Check(IDM_RENDER_TRIM,	  m_sliceModelOn);
	//mfp_debug("Before checking 2.5\n");
	m_menubar->Check(IDM_RENDER_GRID, m_bvhGridOn);
	//m_menubar->Check(IDM_RENDER_RAYS, m_pVMngr->GetSettings()->bIsRayTracing);
}*/

void wxMainWindow::CreateWndToolBar()
{
//    // Delete and recreate the toolbar
//    wxToolBarBase *toolBar = GetToolBar();
//    long style = toolBar ? toolBar->GetWindowStyle() : TOOLBAR_STYLE;
//
//    if (toolBar && m_searchTool && m_searchTool->GetToolBar() == NULL)
//    {
//        // see ~wxMainWindow()
//        toolBar->AddTool(m_searchTool);
//    }
//
//    m_searchTool = NULL;
//
//    delete toolBar;
//
//    SetToolBar(NULL);
//
//    style &= ~(wxTB_HORIZONTAL | wxTB_VERTICAL | wxTB_BOTTOM | wxTB_RIGHT | wxTB_HORZ_LAYOUT);
//	// On top
//	style |= wxTB_TOP;
//	// With tips
//    style &= ~wxTB_NO_TOOLTIPS;
//
//    if ( (style & wxTB_TEXT) && !(style & wxTB_NOICONS))
//        style |= wxTB_HORZ_LAYOUT;
//
//    toolBar = CreateToolBar(style, ID_TOOLBAR);
//
//    PopulateToolbar(toolBar);

}


enum {
	IDT_XDOF, IDT_YDOF, IDT_ZDOF, IDT_FREEDOF
};

void wxMainWindow::PopulateToolbar(wxToolBarBase* toolBar)
{
    // Set up toolbar
    enum {
        Tool_xaxis,
        Tool_yaxis,
		Tool_zaxis,
		#ifdef __WXMSW__
        Tool_faxis,
		#else
        Tool_free,
		#endif
        Tool_Max
    };

	// Declare table of toolbar bitmaps
    wxBitmap toolBarBitmaps[Tool_Max];

	#if USE_XPM_BITMAPS
    	#define INIT_TOOL_BMP(bmp) \
        toolBarBitmaps[Tool_##bmp] = wxBitmap(bmp##_xpm)
	#else // !USE_XPM_BITMAPS
    	#define INIT_TOOL_BMP(bmp) \
        toolBarBitmaps[Tool_##bmp] = wxBITMAP(bmp)
	#endif // USE_XPM_BITMAPS/!USE_XPM_BITMAPS

    INIT_TOOL_BMP(xaxis);
    INIT_TOOL_BMP(yaxis);
    INIT_TOOL_BMP(zaxis);
	#ifdef __WXMSW__
    INIT_TOOL_BMP(faxis);
	#else
    INIT_TOOL_BMP(free);
	#endif

    int w = toolBarBitmaps[Tool_xaxis].GetWidth(),
        h = toolBarBitmaps[Tool_xaxis].GetHeight();


    // this call is actually unnecessary as the toolbar will adjust its tools
    // size to fit the biggest icon used anyhow but it doesn't hurt neither
    toolBar->SetToolBitmapSize(wxSize(w, h));

	#ifdef __WXMSW__
	int width = 24;
	#else
	int width = 16;
	#endif
	int currentX = 5;

	toolBar->AddTool(IDT_XDOF, toolBarBitmaps[0], wxNullBitmap, false, currentX, -1, (wxObject *) NULL, wxT("X dof"));
	currentX += width + 5;
	toolBar->AddTool(IDT_YDOF, toolBarBitmaps[1], wxNullBitmap, false, currentX, -1, (wxObject *) NULL, wxT("Y dof"));
	currentX += width + 5;
	toolBar->AddTool(IDT_ZDOF, toolBarBitmaps[2], wxNullBitmap, false, currentX, -1, (wxObject *) NULL, wxT("Z dof"));
	currentX += width + 5;
	toolBar->AddTool(IDT_FREEDOF, toolBarBitmaps[3], wxNullBitmap, false, currentX, -1, (wxObject *) NULL, wxT("Free dof"));
	currentX += width + 5;
	toolBar->AddSeparator();

   // This is needed!
    toolBar->Realize();
}

void wxMainWindow::SetTitle(const wxString& title)
{
	wxString tmpTile(wxT("FemViewer - "));
	tmpTile.append(wxString(m_pMCtrl->GetCurrentMesh()->Name().c_str(),wxConvUTF8));
	wxFrame::SetTitle(tmpTile);
}

void wxMainWindow::UpdateStatusBar()
{
	int nno =0,ned =0,nfa =0,nel =0;
	if (m_pMCtrl->GetCurrentMesh()) {
		nno = m_pMCtrl->GetCurrentMesh()->GetNumNodes();
		ned = m_pMCtrl->GetCurrentMesh()->GetNumEdges();
		nfa = m_pMCtrl->GetCurrentMesh()->GetNumFaces();
		nel = m_pMCtrl->GetCurrentMesh()->GetNumElems();
	}
	m_status->SetStatusText( wxString::Format(wxT("no: %d"), nno), 1);
	m_status->SetStatusText( wxString::Format(wxT("ed: %d"), ned), 2);
	m_status->SetStatusText( wxString::Format(wxT("fa: %d"), nfa), 3);
	m_status->SetStatusText( wxString::Format(wxT("el: %d"), nel), 4);
}
}
}
