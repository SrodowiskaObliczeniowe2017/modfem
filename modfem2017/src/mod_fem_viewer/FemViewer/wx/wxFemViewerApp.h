#ifndef _FEM_VIEWER_APP_H_
#define _FEM_VIEWER_APP_H_

#include "wx/wxprec.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#if !wxUSE_THREADS
    #error "This sample requires thread support!"
#endif // wxUSE_THREADS

#include "wx/thread.h"
#include "wx/glcanvas.h"


class Thread;

namespace FemViewer {
namespace WxGUI 	{

class RenderThread;
class wxMainWindow;
class wxRenderCanvas;


class wxFemViewerApp : public wxGLApp //,public wxThreadHelper
{
  public:
	static Thread * m_parent;
	static int  Init  (int argc,char **argv,Thread* parent= nullptr);
	static int  Init2 (int argc,char **argv);
	static bool IsInit();
	static void Refresh();
	//static wxFemViewerApp * GetInstance();
	static void ShutDown();

	wxFemViewerApp(void);
	virtual ~wxFemViewerApp(void);
	virtual bool OnInit();
	virtual int  OnRun();
	//virtual int  Update();
	virtual int  OnExit();

	bool m_shutdown;
	wxCriticalSection m_critsect;

  protected:
	/// Private members
	wxMainWindow 	* m_pframe;
	wxRenderCanvas 	* m_pcanvas;
	RenderThread	* m_pThread;

  private:
	/// Block use of these
	DECLARE_NO_COPY_CLASS(wxFemViewerApp);
};


DECLARE_APP(wxFemViewerApp)

}// end namespace WxGUI

typedef WxGUI::wxFemViewerApp GUIEngine;

}// end namespace FemViewer

#endif /* _FEM_VIEWER_APP_H_
*/
