#ifndef _WX_OPENGL_H_
#define _WX_OPENGL_H_

#include "wx/thread.h"
#include "wx/glcanvas.h"
#include "fv_timer.h"

namespace FemViewer {
class ViewManager;
class ModelCtrl;

namespace WxGUI		{


#define EVT_MY_REFRESH_CANVAS 20000

enum {
	MY_RENDER_THREAD_UPDATE = wxID_HIGHEST + 2000,
	MY_RENDER_THREAD_COMPLETE,
};

class wxMainWindow;
class wxRenderCanvas;

/*class RenderThread
		: public wxThread
{
  public:
	RenderThread(wxRenderCanvas *pcanvas,unsigned delay = 16)
    : wxThread(wxTHREAD_DETACHED)
    , _Handler(pcanvas)
    , m_iRefreshDelay(delay)
  	{ ; }

  protected:
	virtual void* Entry();
	wxMainWindow *_Handler;
	//wxWindowID 	m_iGLCanvasID;
	unsigned 		m_iRefreshDelay;
};
*/

//class wxMainWindow;

class wxRenderCanvas : public wxGLCanvas
{
	friend class RenderThread;
  public:
	wxRenderCanvas(wxFrame* parent,
			 	   wxWindowID id = wxID_ANY,
			 	   int* gl_attrib = NULL,
			 	   const wxPoint& pos = wxDefaultPosition,
			 	   const wxSize& size = wxDefaultSize,
			 	   long style = 0,
			 	   const wxString& name = wxT("FemViewer"));
	~wxRenderCanvas();
	void CleanUp();
	void InitGL();
	int GetWidth();
	int GetHeight();
	wxBitmap GetScreenShots();
	void PaintNow();

  //protected:
	void OnRefresh(wxCommandEvent& event);
	void OnThreadUpdate(wxThreadEvent& event);
	void OnThreadComplete(wxThreadEvent& event);
	void OnPaint(wxPaintEvent& event);
	void OnIdle(wxIdleEvent& event);
	int  Modifires(wxMouseEvent& event);
	void OnSize(wxSizeEvent& event);
	void OnEraseBackground(wxEraseEvent& event);
	void OnMouseMove(wxMouseEvent& event);
	void OnLeftBtnDown(wxMouseEvent& event);
	void OnLeftBtnUp(wxMouseEvent& event);
	void OnMiddleBtnDown(wxMouseEvent& event);
	void OnMiddleBtnUp(wxMouseEvent& event);
	void OnRightBtnDown(wxMouseEvent& event);
	void OnRightBtnUp(wxMouseEvent& event);
	//void OnMouse(wxMouseEvent& event);
	void OnMouseWheel(wxMouseEvent& event);
	

	void OnEnterWindow(wxMouseEvent& event);
	void OnWindowCreate(wxWindowCreateEvent& event);
	/*void GetViewportOffset() const;*/
	


  private:
	bool			m_rendering;
	bool			m_binitGL;
	int 			m_ndraging;
	wxGLContext		*m_context, *m_dummyconetxt;
	ModelCtrl		*m_modelCtr;
	ViewManager		*m_viewCtr;
	int 			m_frame_cnt;
	fv_timer		m_timer;
	
	wxDECLARE_NO_COPY_CLASS(wxRenderCanvas);
	wxDECLARE_EVENT_TABLE();

};

}// end namespace WxGUI
}// end namespace FemViewer


#endif /* _WX_OPENGL_H_
*/

