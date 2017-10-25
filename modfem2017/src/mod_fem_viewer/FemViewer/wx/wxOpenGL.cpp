#include "fv_inc.h"
#include "wxFemViewerApp.h"
#include "wxOpenGL.h"
#include "ModelControler.h"
#include "ViewManager.h"
#include "Log.h"
#include "fv_txt_utls.h"
#include <wx/wx.h>

//#include <GL/glu.h>
//#include <GL/gl.h>

#ifdef __WXDEBUG__
#include<iostream>
#endif

namespace FemViewer {
namespace WxGUI		{


//void* RenderThread::Entry()
//{
//	//mfp_log_debug("entry\n");
//	//m_pGLCanvas->SeCurrent(m_pGLCanvas->)
//	//wxPaintDC dc(m_pGLCanvas);
//	static int i;
//	wxClientDC dc(m_pGLCanvas);
//
//	while(!TestDestroy())
//	{
//		//wxLogError("in loop\n");
//
//		mfp_log_debug("loop %d\n",i++);
//		wxMutexGuiEnter();
//		m_pGLCanvas->SetCurrent(*(m_pGLCanvas->m_context));
//		m_pGLCanvas->m_viewCtr->DisplayCallback();
//		m_pGLCanvas->SwapBuffers();
//		m_pGLCanvas->SetCurrent(*(m_pGLCanvas->m_dummyconetxt));
//		//m_pGLCanvas->PaintNow();
//		wxMutexGuiLeave();
//		Sleep(50);
//		//m_pGLCanvas->SetCurrent(*m_pGLCanvas->m_context);
//		//ViewManagerInst().DisplayCallback();
//		//m_pGLCanvas->SwapBuffers();
//		//m_pGLCanvas->SetCurrent(*m_pGLCanvas->m_dummyconetxt);
//		//wxCommandEvent event(wxEVT_ERASE_BACKGROUND, EVT_MY_REFRESH_CANVAS);
//		//m_pGLCanvas->GetEventHandler()->AddPendingEvent(event);
//		//
//		//wxQueueEvent(m_pGLCanvas,new wxThreadEvent(MY_RENDER_THREAD_UPDATE));
//	}
//	mfp_log_debug("Leaving\n");
//	wxQueueEvent(m_pGLCanvas,new wxThreadEvent(MY_RENDER_THREAD_COMPLETE));
//	return (void *)0;
//}



BEGIN_EVENT_TABLE(wxRenderCanvas, wxWindow)
    EVT_SIZE(wxRenderCanvas::OnSize)
    EVT_PAINT(wxRenderCanvas::OnPaint)
 	EVT_ERASE_BACKGROUND(wxRenderCanvas::OnEraseBackground)
    EVT_IDLE(wxRenderCanvas::OnIdle)
 	EVT_MOTION(wxRenderCanvas::OnMouseMove)
	EVT_MOUSEWHEEL(wxRenderCanvas::OnMouseWheel)
	EVT_LEFT_DOWN(wxRenderCanvas::OnLeftBtnDown)
	EVT_LEFT_UP(wxRenderCanvas::OnLeftBtnUp)
	EVT_MIDDLE_DOWN(wxRenderCanvas::OnMiddleBtnDown)
	EVT_MIDDLE_UP(wxRenderCanvas::OnMiddleBtnUp)
	EVT_RIGHT_DOWN(wxRenderCanvas::OnRightBtnDown)
	EVT_RIGHT_UP(wxRenderCanvas::OnRightBtnUp)
	EVT_WINDOW_CREATE(wxRenderCanvas::OnWindowCreate)
    //EVT_MOUSE_EVENTS(wxRenderCanvas::OnMouse)
	//EVT_CHAR(wxRenderCanvas::OnChar);
	EVT_ENTER_WINDOW(wxRenderCanvas::OnEnterWindow )
	EVT_MENU(EVT_MY_REFRESH_CANVAS, wxRenderCanvas::OnRefresh)
	//EVT_THREAD(MY_RENDER_THREAD_UPDATE, wxRenderCanvas::OnThreadUpdate)
	//EVT_THREAD(MY_RENDER_THREAD_COMPLETE, wxRenderCanvas::OnThreadComplete)
END_EVENT_TABLE()

wxRenderCanvas::wxRenderCanvas(
		wxFrame* parent,
		wxWindowID id,
		int* gl_attrib,
		const wxPoint& pos,
		const wxSize& size,
		long style,
		const wxString& name
		)
: wxGLCanvas(static_cast<wxWindow*>(parent),
             id, gl_attrib, pos, size,
		     style | wxFULL_REPAINT_ON_RESIZE, name)
, m_rendering(false)
, m_binitGL(false)
, m_ndraging(0)
, m_frame_cnt(0)
, m_timer()
{
	//mfp_log_debug("-----------------------Ctr\n");

	//XInitThreads();

	// Setup GL context
	m_context = new wxGLContext(this);
	wxASSERT(m_context);

	// Set dummy context for release state
	//m_dummyconetxt = new wxGLContext(this);
	//wxASSERT(m_dummyconetxt);

	//SetCurrent(*m_context);
	SetCurrent(*m_context);
	int result = glutGet( GLUT_INIT_STATE );

	//mfp_debug("init glut = %d\n",result);
	if (result == 0) {
		//mfp_debug("initializing glut");
		int argc = 1;
		char* argv[1] = { wxTheApp->argv[0].char_str() };
		//mfp_debug("Application name: %s \n",argv[0]);
		GLCore::init(argc,argv);
		//GLCore::instance().initGLLists(NUM_TOTAL_LISTS);
	}

	m_modelCtr = & ModelCtrlInst();
	m_viewCtr  = & ViewManagerInst();

}

wxRenderCanvas::~wxRenderCanvas()
{
	//mfp_log_debug("Destructor: ~wxRendercanvas");
	//wxDELETE(m_context);
	//mfp_log_debug("After deleting context\n");
	//wxDELETE(m_dummyconetxt);
}

void wxRenderCanvas::PaintNow()
{
	static double dLastTime = 0.0;
	//static int i;
	//mfp_log_debug("PainNow %d\n",i++);
	//wxClientDC dc(this);
	//SetCurrent(*m_context);
	m_viewCtr->DisplayCallback();
	FV_CHECK_ERROR_GL();
	SwapBuffers();

#if 0
	m_frame_cnt++;
	double currentTime = m_timer.get_time_ms();
	double delta = currentTime - dLastTime;
	if (delta > 1000.0 ) {
		// If last prinf() was more than 1 sec ago
		printf("%d %.3f ms/frame, %.1f FPS\n",m_frame_cnt, delta / double(m_frame_cnt),
				double(m_frame_cnt) * delta  * 0.001);
		m_frame_cnt = 0;
		dLastTime = currentTime;
	}

#endif
}

void wxRenderCanvas::CleanUp()
{
	//mfp_debug("ClenaUp z CLERA z MODCTRL\n");
	SetCurrent( *m_context );
	m_modelCtr->Do(ModelCtrl::CLEAR);
	m_viewCtr->GLSupport() = false;
	FV_CHECK_ERROR_GL();
}

int wxRenderCanvas::GetWidth()
{
	wxClientDC dc(this);
	wxSize size = dc.GetSize();
	return size.GetWidth();
}

int wxRenderCanvas::GetHeight()
{
	wxClientDC dc(this);
	wxSize size = dc.GetSize();
	return size.GetHeight();
}

wxBitmap wxRenderCanvas::GetScreenShots()
{
	// Read the OpenGL image into a pixel array
	GLint view[4];
	glGetIntegerv(GL_VIEWPORT, view);
	void* pixels = malloc(3 * view[2] * view[3]); // must use malloc
	glPixelStorei(GL_PACK_ALIGNMENT, 1);
	glReadBuffer(GL_BACK_LEFT);
	glReadPixels(0, 0, view[2], view[3], GL_RGB, GL_UNSIGNED_BYTE, pixels);
	// Put the image into a wxImage
	wxImage image((int) view[2], (int) view[3]);
	image.SetData((unsigned char*) pixels);
	image = image.Mirror(false);
	return wxBitmap(image);
}

void wxRenderCanvas::OnRefresh(wxCommandEvent &WXUNUSED(event)) {
	Refresh();
}

void wxRenderCanvas::OnThreadUpdate(wxThreadEvent &WXUNUSED(event)) {
	//mfp_log_debug("Render thread update\n");
	Update();
}

void wxRenderCanvas::OnThreadComplete(wxThreadEvent &WXUNUSED(event)) {
	//mfp_log_debug("Render thread complete\n");
}

void wxRenderCanvas::OnPaint(wxPaintEvent& WXUNUSED(event))
{
	static int i;
	//mfp_log_debug("in paint %d\n",i++);
	//if(!IsShownOnScreen()) return;
	wxPaintDC dc(this);
	//wxMutexGuiEnter();
	PaintNow();
	//wxMutexGuiLeave();

	//{
		//wxCriticalSectionLocker lock(wxGetApp().m_pcritSec);
		//SetCurrent(*m_context);
		//m_viewCtr->DisplayCallback();
		//SwapBuffers();
	//}
    //SetCurrent(*m_dummyconetxt);
}


void wxRenderCanvas::OnIdle(wxIdleEvent& event)
{
	//static int i;
	//mfp_log_debug("OnIdle %d\n",i++);
	//return;
	//if (wxGetApp().m_shutdown) m_parent->Close(true);
	//m_viewCtr->DisplayCallback();
	wxClientDC dc(this);
	//wxMutexGuiEnter();
	PaintNow();
	//wxMutexGuiLeave();
	//event.RequestMore();
}
int wxRenderCanvas::Modifires(wxMouseEvent& event)
{
	int modiffier = 0;
	if(event.ShiftDown()) modiffier |= 0x1;
	if(event.ControlDown()) modiffier |= 0x2;
	if(event.AltDown()) modiffier |= 0x4;
	return modiffier;
}

void wxRenderCanvas::OnSize(wxSizeEvent &event)
{
	static int i;
	if (!IsShownOnScreen())
	 	return;
	int w,h;
	//wxGLCanvas::OnSize(event );
	SetCurrent(*m_context);
	GetClientSize( &w, &h );
	//mfp_log_debug("OnSize %d %d %d\n",i++,w,h);
	glViewport( 0, 0, event.GetSize().x, event.GetSize().y );
	FV_CHECK_ERROR_GL();
	m_viewCtr->ReshapeCallback(w,h);
	
	Update();
}

void wxRenderCanvas::OnEraseBackground(wxEraseEvent& WXUNUSED(event))
{
    // Do nothing, to avoid flashing on MSW
	//m_VMgr.DisplayCallback();
	//mfp_log_debug("Onerasebackground\n");
	Update();
}

void wxRenderCanvas::OnMouseMove(wxMouseEvent& event)
{
	if( m_ndraging ) {
		m_viewCtr->MouseMotionCallback( event.GetX(), event.GetY() );
		Refresh(false);
	}
}

void wxRenderCanvas::OnLeftBtnDown(wxMouseEvent& event)
{
	int modif = Modifires( event );
	m_viewCtr->MouseButtonCallback(0,0,modif,event.GetX(),event.GetY());
	Refresh(false);
	m_ndraging = 1;
}

void wxRenderCanvas::OnLeftBtnUp(wxMouseEvent& event)
{
	int modif = Modifires( event );
	m_viewCtr->MouseButtonCallback(0,1,modif,event.GetX(),event.GetY());
	Refresh(false);
	m_ndraging = 0;
}

void wxRenderCanvas::OnMiddleBtnDown(wxMouseEvent& event)
{
	int modif = Modifires( event );
	m_viewCtr->MouseButtonCallback(1,0,modif,event.GetX(),event.GetY());
	Refresh(false);
	m_ndraging = 1;
}

void wxRenderCanvas::OnMiddleBtnUp(wxMouseEvent& event)
{
	int modif = Modifires( event );
	m_viewCtr->MouseButtonCallback( 1,1,modif,event.GetX(),event.GetY() );
	Refresh(false);
	m_ndraging = 0;
}

void wxRenderCanvas::OnRightBtnDown(wxMouseEvent& event)
{
	int modif = Modifires( event );
	m_viewCtr->MouseButtonCallback(2,0,modif,event.GetX(),event.GetY());
	Refresh(false);
	m_ndraging = 1;
}

void wxRenderCanvas::OnRightBtnUp(wxMouseEvent& event)
{
	int modif = Modifires(event);
	m_viewCtr->MouseButtonCallback( 2,1,modif,event.GetX(),event.GetY() );
	Refresh(false);
	m_ndraging = 0;
}

//void wxRenderCanvas::OnMouse(wxMouseEvent& event)
//{
//
//	if(event.Dragging()) {
//		m_VMgr.MouseMotionCallback(event.GetX(), event.GetY());
//	} else {
//
//		int modif = 0, button;
//		if(event.ShiftDown()) modif |= 0x1;
//		if(event.ControlDown()) modif |= 0x2;
//		if(event.AltDown()) modif |= 0x4;
//
//		if(event.ButtonUp()) {
//			button = event.GetButton();
//			m_VMgr.MouseButtonCallback(button,1,modif,event.GetX(),event.GetY());
//		} else if(event.ButtonDown()) {
//			button = event.GetButton();
//			m_VMgr.MouseButtonCallback(button,0,modif,event.GetX(),event.GetY());
//		}
//	}
//	Refresh(false);
//}

void wxRenderCanvas::OnMouseWheel(wxMouseEvent& event)
{
	int rot = event.GetWheelRotation();
	if(rot < 0) {
		m_viewCtr->MouseButtonCallback(3,0,0,0,0);
		//Refresh(false);
	} else if(rot > 0) {
		m_viewCtr->MouseButtonCallback(4,0,0,0,0);
		//Refresh(false);
	} else {
		m_viewCtr->MouseButtonCallback(4,1,0,0,0);

	}
	Refresh(false);
}


void wxRenderCanvas::OnEnterWindow(wxMouseEvent& WXUNUSED(event))
{
	SetFocus();
}

void wxRenderCanvas::OnWindowCreate(wxWindowCreateEvent& event)
{
	mfp_log_debug("======tworzenie okna===\n");
}

//void wxRenderCanvas::Redraw(const std::string& title)
//{
//	m_VMgr.SetTitle(title);
//	m_VMgr.GetGraphicData().Reset();
//	m_pModel->SetGoemetry(m_VMgr.GetGraphicData().GetRootObject());
//	m_pModel->Redraw();
//}
//void wxRenderCanvas::GetViewportOffset() const
//{
//	SetCurrent(*m_glRC);
//	wxGLCanvas::GetPosition(
//}

void wxRenderCanvas::InitGL()
{
	//mfp_log_debug("INIT OPENGL!!!!\n");
//	SetCurrent(*m_context);
//	int result = glutGet( GLUT_INIT_STATE );
//
//	//mfp_debug("init glut = %d\n",result);
//	if (result == 0) {
//		mfp_debug("initializing glut");
//		int argc = 1;
//		char* argv[1] = { wxTheApp->argv[0].char_str() };
//		GLCore::init(argc,argv);
//	}


	// naruszenie ochrony pamieci?????
	//mfp_log_debug("x: %d y: %d\n",GetSize().x,GetSize().y);
	//m_viewCtr->Init(GetSize().x,GetSize().y);
	//return NULL;

//	RenderThread* pThread = new RenderThread(this);
//	wxASSERT_MSG(pThread!=NULL,wxT("Can't allocate thread object"));
//
//	if (pThread->Create() != wxTHREAD_NO_ERROR )
//	{
//		wxLogError(wxT("Can't create thread!"));
//		delete pThread;
//		pThread = NULL;
//	}

	return;
}


}// end namespace WxGUI
}// end namespace FemViewer
