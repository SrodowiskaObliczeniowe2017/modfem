#ifndef WinGL_H_
#define WinGL_H_

#include"WinApp.h"
#include"ViewManager.h"

namespace FemViewer {
namespace Win		{ 

class WinGL : public CWindowImpl<WinGL, CWindow, CControlWinTraits>
{
public:

	WinGL(FemViewer::ViewManager& pVmgr, HWND hWnd);
	~WinGL(void);

	void Render(void);
	BEGIN_MSG_MAP(WinGL)
		MESSAGE_HANDLER(WM_PAINT, OnPaint)
		MESSAGE_HANDLER(WM_LBUTTONDOWN, OnLBttnDown)
		MESSAGE_HANDLER(WM_LBUTTONUP, OnLBttnUp)
		MESSAGE_HANDLER(WM_RBUTTONDOWN, OnRBttnDown)
		MESSAGE_HANDLER(WM_RBUTTONUP, OnRBttnUp)
		MESSAGE_HANDLER(WM_MOUSEMOVE, OnMouseMove)
		MESSAGE_HANDLER(WM_SIZE, OnSize)
		MESSAGE_HANDLER(WM_CHAR, OnChar )
    END_MSG_MAP()

private:
	enum {
		MBTN_DOWN = 0,
		MBTN_UP   = 1
	};

	LRESULT OnPaint(UINT nMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);
	LRESULT OnLBttnDown(UINT nMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);
	LRESULT OnLBttnUp(UINT nMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);
	LRESULT OnRBttnDown(UINT nMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);
	LRESULT OnRBttnUp(UINT nMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);
	LRESULT OnMouseMove(UINT nMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);
	LRESULT OnSize(UINT nMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);
	LRESULT OnChar(UINT nMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled);
	
	void InitGL(int width,int height);
	void InitData();
	bool SetPixelFormat(HDC hdc, int colorBits, int depthBits, int stencilBits);
	int  FindPixelFormat(HDC hdc, int colorBits, int depthBits, int stencilBits);
	bool CreateContext(int colorBits,int depthBits,int stencilBits1);
	void DestroyContext();
	void ResetView(void) { if(IsWindow()) Invalidate(FALSE); }
	void Refresh(BOOL bflag = FALSE) { Invalidate(bflag); }
	unsigned  GetModifires() const;
	HDC			 m_hDC;
	HGLRC		 m_hRC;
	FemViewer::ViewManager& m_Vmgr;
public:
	LRESULT OnMenu(UINT nMsg, WPARAM wParam, LPARAM lParam);
};

inline void WinGL::Render(void)
{
	//::DrawScene();
	Invalidate(TRUE);
}

inline LRESULT WinGL::OnPaint(UINT nMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled)
{
	nMsg; wParam; lParam; bHandled;

    BOOL bRet = wglMakeCurrent( m_hDC, m_hRC );
    ATLASSERT(bRet); bRet;
	m_Vmgr.DisplayCallback();
    ATLASSERT(SwapBuffers( m_hDC ));
    ATLASSERT(wglMakeCurrent( NULL, NULL ));
    ATLASSERT(ValidateRect( NULL ));
    return 0L;
}

inline LRESULT WinGL::OnLBttnDown(UINT nMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled)
{
	nMsg; bHandled;

	int posx = LOWORD(lParam);
	int posy = HIWORD(lParam);
	m_Vmgr.MouseButtonCallback(0,MBTN_DOWN,GetModifires(),posx,posy);
    SetCapture();
	Refresh();
    return 0L;
}

inline LRESULT WinGL::OnLBttnUp(UINT nMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled)
{
	nMsg; wParam; lParam; bHandled;

	UINT nFlags = (UINT) wParam;
	int posx = LOWORD(lParam);
	int posy = HIWORD(lParam);
	int modif = 0;
	m_Vmgr.MouseButtonCallback(0,MBTN_UP,GetModifires(),posx,posy);
    ReleaseCapture();
	Refresh();
    return 0L;
}

inline LRESULT WinGL::OnRBttnDown(UINT nMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled)
{
	nMsg; bHandled;

	int posx = LOWORD(lParam);
	int posy = HIWORD(lParam);
	BOOL bRet = wglMakeCurrent( m_hDC, m_hRC );
    ATLASSERT(bRet); bRet;
	m_Vmgr.MouseButtonCallback(1,MBTN_DOWN,GetModifires(),posx,posy);
	ATLASSERT(wglMakeCurrent( NULL, NULL ));
    SetCapture();
	Refresh();
    return 0L;
}

inline LRESULT WinGL::OnRBttnUp(UINT nMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled)
{
	nMsg; wParam; lParam; bHandled;

	int posx = LOWORD(lParam);
	int posy = HIWORD(lParam);
	m_Vmgr.MouseButtonCallback(1,MBTN_UP,GetModifires(),posx,posy);
    ReleaseCapture();
	Refresh();
    return 0L;
}

inline LRESULT WinGL::OnMouseMove(UINT nMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled)
{
	nMsg; bHandled;

	int posx = LOWORD(lParam);
	int posy = HIWORD(lParam);
	m_Vmgr.MouseMotionCallback(posx,posy);
	Refresh();
    return 0L;
}

inline LRESULT WinGL::OnSize(UINT nMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled)
{
	nMsg; wParam; bHandled;
	int cx = LOWORD(lParam);
    int cy = HIWORD(lParam);
	::SendMessage(MainApp::GetInstance()->GetStatusBar(),WM_SIZE,0,0);
    if (!m_hWnd || !cy) return 0L;
    HDC hDC = ::GetDC(m_hWnd);
    ATLASSERT(hDC);
	ATLASSERT(m_hRC);
    ATLVERIFY(::wglMakeCurrent( hDC, m_hRC ));
	m_Vmgr.ReshapeCallback(cx,cy);
    ATLVERIFY(ReleaseDC(hDC));
    return 0L;
}

inline LRESULT WinGL::OnChar(UINT nMsg, WPARAM wParam, LPARAM lParam, BOOL& bHandled)
{
    nMsg; lParam; bHandled;

    switch (wParam)
    {
        case ' ': ResetView();
	}
    return 0L;
}

inline LRESULT WinGL::OnMenu(UINT nMsg, WPARAM wParam, LPARAM lParam)
{
	nMsg, wParam, lParam;
	int wmId = LOWORD(wParam);
	m_Vmgr.MenuCallback(wmId);
	Refresh();
	return 0L;
}

inline unsigned WinGL::GetModifires() const
{
	unsigned int modif = 0;
	static bool btn_state;
	if(GetAsyncKeyState(VK_LSHIFT) & 0x8000)
		modif |= 0x1;
	if(GetAsyncKeyState(VK_LCONTROL) & 0x8000)
		modif |= 0x2;
	if(GetAsyncKeyState(VK_RMENU) & 0x8000)
		modif |= 0x4;
/*	modif |= (1 << (btn_state ? 1 : 0));
	btn_state = !!( ::GetKeyState(VK_CONTROL) & 0x8000 );
	modif |= (1 << (btn_state ? 2 : 0));
	btn_state = !!( ::GetKeyState(VK_MENU) & 0x8000 );
	modif |= (1 << (btn_state ? 4 : 0))*/;
	return modif;
}

} // end namsepace Win
} // end namespace FemViewer

#endif /* WinGL_H_
*/