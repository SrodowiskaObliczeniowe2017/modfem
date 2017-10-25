#include"Common.h"
#include"WinGL.h"
#include"ModelControler.h"
#include"Log.h"
// namespaces
namespace FemViewer {
namespace Win       {

#define GL_WIDTH  800
#define GL_HEIGHT 600

#undef SubclassWindow
WinGL::WinGL(FemViewer::ViewManager& pVmgr,HWND hWnd) 
: m_hDC(0), m_hRC(0), m_Vmgr(pVmgr)
{
	ATLVERIFY(this->SubclassWindow( hWnd ));
	CreateContext(32,24,8);

	RECT rcWnd, rcGL = {0, 0, GL_WIDTH, GL_HEIGHT};
    GetWindowRect(&rcWnd);

    AdjustWindowRectEx(&rcGL, 
                       GetWindowLong(GWL_STYLE), 
                       GetMenu() != NULL, 
                       GetWindowLong(GWL_EXSTYLE)
					   );

	SetWindowPos(NULL, rcWnd.left, rcWnd.top, rcGL.right - rcGL.left, rcGL.bottom - rcGL.top, SWP_NOZORDER | SWP_SHOWWINDOW);
	
	this->InitGL(GL_WIDTH,GL_HEIGHT);
	//this->InitData();
	logd("Constructor WinGL");
}

WinGL::~WinGL(void)
{
	if(m_hWnd) {
		UnsubclassWindow();
		DestroyContext();
	}
	logd("Destructor WinGL");
}

void WinGL::InitGL(int width,int height)
{
	ATLVERIFY(wglMakeCurrent(m_hDC,m_hRC));
	m_Vmgr.InitGL(width,height);
	ATLVERIFY(wglMakeCurrent( NULL, NULL ));
}

void WinGL::InitData()
{
	ATLVERIFY(wglMakeCurrent(m_hDC,m_hRC));
	ModelCtrlInst().Do(ModelCtrl::UPDATE);
	ATLVERIFY(wglMakeCurrent( NULL, NULL ));
}

int WinGL::FindPixelFormat(HDC hdc, int colorBits, int depthBits, int stencilBits)
{
    int currMode;                               // pixel format mode ID
    int bestMode = 0;                           // return value, best pixel format
    int currScore = 0;                          // points of current mode
    int bestScore = 0;                          // points of best candidate
    PIXELFORMATDESCRIPTOR pfd;

    // search the available formats for the best mode
    bestMode = 0;
    bestScore = 0;
    for(currMode = 1; ::DescribePixelFormat(hdc, currMode, sizeof(pfd), &pfd) > 0; ++currMode)
    {
        // ignore if cannot support opengl
        if(!(pfd.dwFlags & PFD_SUPPORT_OPENGL))
            continue;

        // ignore if cannot render into a window
        if(!(pfd.dwFlags & PFD_DRAW_TO_WINDOW))
            continue;

        // ignore if cannot support rgba mode
        if((pfd.iPixelType != PFD_TYPE_RGBA) || (pfd.dwFlags & PFD_NEED_PALETTE))
            continue;

        // ignore if not double buffer
        if(!(pfd.dwFlags & PFD_DOUBLEBUFFER))
            continue;

        // try to find best candidate
        currScore = 0;

        // colour bits
        if(pfd.cColorBits >= colorBits) ++currScore;
        if(pfd.cColorBits == colorBits) ++currScore;

        // depth bits
        if(pfd.cDepthBits >= depthBits) ++currScore;
        if(pfd.cDepthBits == depthBits) ++currScore;

        // stencil bits
        if(pfd.cStencilBits >= stencilBits) ++currScore;
        if(pfd.cStencilBits == stencilBits) ++currScore;

        // alpha bits
        if(pfd.cAlphaBits > 0) ++currScore;

        // check if it is best mode so far
        if(currScore > bestScore)
        {
            bestScore = currScore;
            bestMode = currMode;
        }
    }

    return bestMode;
}

bool WinGL::SetPixelFormat(HDC hdc, int colorBits, int depthBits, int stencilBits)
{
    PIXELFORMATDESCRIPTOR pfd;

    // find out the best matched pixel format
    int pixelFormat = FindPixelFormat(hdc, colorBits, depthBits, stencilBits);
    if(pixelFormat == 0) return false;

    // set members of PIXELFORMATDESCRIPTOR with given mode ID
    ::DescribePixelFormat(hdc, pixelFormat, sizeof(pfd), &pfd);

    // set the fixel format
    if(!::SetPixelFormat(hdc, pixelFormat, &pfd)) return false;

    return true;
}

bool WinGL::CreateContext(int colorBits,int depthBits,int stencilBits)
{
	 // retrieve a handle to a display device context
   m_hDC = ::GetDC(m_hWnd);
   ATLASSERT(m_hDC != NULL);

   // set pixel format
   if(!SetPixelFormat(m_hDC, colorBits, depthBits, stencilBits))
   {
        ::MessageBox(0, _T("Cannot set a suitable pixel format."), _T("Error"), MB_ICONEXCLAMATION | MB_OK);
        ::ReleaseDC(m_hWnd, m_hDC);  // remove device contex
		m_hDC = NULL;
        return FALSE;
   }

    // create a new OpenGL rendering context
    m_hRC = ::wglCreateContext(m_hDC);
    //::wglMakeCurrent(hdc, hglrc);
	 ATLVERIFY(wglMakeCurrent(m_hDC,m_hRC));
    //::ReleaseDC(m_hWnd, m_hDC);
	//m_hDC = NULL;

    return TRUE;
}

void WinGL::DestroyContext()
{
    if(!m_hDC || !m_hRC) return;

    // delete DC and RC
    ::wglMakeCurrent(NULL, NULL);
    ::wglDeleteContext(m_hRC);
	::ReleaseDC(m_hWnd, m_hDC);
    m_hDC = NULL;
    m_hRC = NULL;
}

}// end namespace Win
}// end namespace Femviwer