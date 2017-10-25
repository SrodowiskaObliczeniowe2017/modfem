#ifndef WinApp_H_
#define WinApp_H_

//#define WINDOWS_LEAN_AND_MEAN
//#include<Windows.h>
#include"Handle.h"

namespace FemViewer {
namespace Win		{

void printError(LPTSTR lpszFunction);
/* Forwward declaretion */
class WinGL;
class DialogWindow;
// Main class of the application
class MainApp : public Handle<HWND>
{
public:
	static volatile bool m_bRun;
	static MainApp * GetInstance();
	BOOL IsInit() { return ::IsWindow(static_cast<Handle<HWND> >(*GetInstance())); }
	BOOL Init(int argc,char** argv);
	void ShutDown();
	int  Run();
	static void SetStatusBarInfo(HWND hwndSB, int part, unsigned style, char *fmt, ...);

	~MainApp();
	BOOL Create();
	void Update();
	void Destroy();
	WinGL*& GetWindowGL();
	HWND& GetStatusBar();  
	DialogWindow*& GetChildDlg(); 
protected:
	static MainApp m_application;
	HINSTANCE      m_hInstance;
	HWND		   m_hStatusBar;
	WNDCLASSEX     m_wndClassEx;
	WinGL*         m_pMainWnd;
	DialogWindow*  m_pChildDialog;
private:
	MainApp();
	MainApp(const MainApp& rhs);
	MainApp& operator=(const MainApp& rhs);
};

inline WinGL*& MainApp::GetWindowGL() 
{ 
	return m_pMainWnd; 
}

inline HWND& MainApp::GetStatusBar()
{
	return m_hStatusBar;
}

inline DialogWindow*& MainApp::GetChildDlg()
{
	return m_pChildDialog;
}


} // endnamespace Win

typedef Win::MainApp AppEngine;

} // end namespace FemViewer
#endif