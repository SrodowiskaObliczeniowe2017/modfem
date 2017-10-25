#include"defs.h"
#include<Windows.h>
#include"Common.h"
#include"WinApp.h"
#include"WinGL.h"
#include"DialogWindow.h"
#include"ViewManager.h"
#include"ModelControler.h"
#include"win/resource.h"
#include"Log.h"
#include"fv_txt_utls.h"
#include"fv_assert.h"


#include<cstdio>
#include<cstdarg>
#include<iostream>
#include<process.h>
#include<commctrl.h>  
#include <strsafe.h>

#pragma comment( lib, "comctl32.lib" )
// namespaces
using namespace std;
namespace FemViewer {
namespace Win		{
// defines
#define MAX_LOADSTRING 100


// Global Variables:
const TCHAR AppName[] = TEXT("FemViewer");



// Forword procedure declaretion
INT_PTR CALLBACK About(HWND hDlg, UINT msg, WPARAM wParam, LPARAM lParam);
INT_PTR CALLBACK ApproxDlgProc(HWND hDlg, UINT msg, WPARAM wParam, LPARAM lParam);
LRESULT CALLBACK WndProcedure(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam);

volatile bool MainApp::m_bRun = true;
MainApp MainApp::m_application;

// Get handle to main application
MainApp * MainApp::GetInstance() 
{ 
	return &m_application;
}

void printError(LPTSTR lpszFunction) 
{ 
    // Retrieve the system error message for the last-error code

    LPVOID lpMsgBuf;
    LPVOID lpDisplayBuf;
    DWORD dw = GetLastError(); 

    FormatMessage(
        FORMAT_MESSAGE_ALLOCATE_BUFFER | 
        FORMAT_MESSAGE_FROM_SYSTEM |
        FORMAT_MESSAGE_IGNORE_INSERTS,
        NULL,
        dw,
        MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT),
        (LPTSTR) &lpMsgBuf,
        0, NULL );

    // Display the error message and exit the process

    lpDisplayBuf = (LPVOID)LocalAlloc(LMEM_ZEROINIT, 
        (lstrlen((LPCTSTR)lpMsgBuf) + lstrlen((LPCTSTR)lpszFunction) + 40) * sizeof(TCHAR)); 
    StringCchPrintf((LPTSTR)lpDisplayBuf, 
        LocalSize(lpDisplayBuf) / sizeof(TCHAR),
        TEXT("%s failed with error %d: %s"), 
        lpszFunction, dw, lpMsgBuf); 
    MessageBox(NULL, (LPCTSTR)lpDisplayBuf, TEXT("Error"), MB_OK); 

    LocalFree(lpMsgBuf);
    LocalFree(lpDisplayBuf);
    //ExitProcess(dw); 
}

// Init function
BOOL MainApp::Init(int argc,char **argv)
{
	if( !glutGet( GLUT_INIT_STATE ) ) glutInit( &argc,argv );
	int err = TRUE;
	if(IsInit() == FALSE)
	{
		INITCOMMONCONTROLSEX commonCtrls;
        commonCtrls.dwSize = sizeof(commonCtrls);
        commonCtrls.dwICC = ICC_UPDOWN_CLASS;
        ::InitCommonControlsEx( &commonCtrls );
		// Init ModelControler
		ModelCtrlInst().Do(ModelCtrl::INIT);
	    cout << "StarTING NEW THREAD\n";
		if(Create())
		{
			GetInstance()->Update();
			err = MainApp::Run();
		}
	} 
	return err;
}

void MainApp::ShutDown()
{
	GetInstance()->Destroy();
	//if( m_papplication )
	//{
	//	delete m_papplication;
	//	m_papplication = NULL;
	//}
}

// The main constr.
MainApp::MainApp()
: m_hStatusBar(NULL)
, m_pMainWnd(NULL)
{

	// Initializing the application using the application member variable
	m_hInstance = ::GetModuleHandle(NULL);
	m_wndClassEx.cbSize        = sizeof(WNDCLASSEX);
	m_wndClassEx.style         = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
	m_wndClassEx.lpfnWndProc   = WndProcedure;
	m_wndClassEx.cbClsExtra    = 0;
	m_wndClassEx.cbWndExtra    = 0;
	m_wndClassEx.hInstance     = m_hInstance;
	m_wndClassEx.hIcon         = LoadIcon(NULL, IDI_APPLICATION);
	m_wndClassEx.hCursor       = LoadCursor(NULL, IDC_ARROW);
	m_wndClassEx.hbrBackground = NULL; //cause WM_ERASEBKGND //(HBRUSH)(COLOR_WINDOW+1)
	m_wndClassEx.lpszMenuName  = MAKEINTRESOURCE(IDR_MAIN_MENU);
	m_wndClassEx.lpszClassName = AppName;
	m_wndClassEx.hIconSm       = LoadIcon(NULL, IDI_APPLICATION);
	logd("Constructor MainApp");
}

// The main destr.
MainApp::~MainApp()
{
	::UnregisterClass(AppName, m_hInstance);
	//m_pMainWnd->DestroyWindow();
	if(m_pMainWnd) delete m_pMainWnd;
	m_hInstance = NULL;
	m_pMainWnd  = NULL;
	logd("Destructor MainApp");
}

BOOL MainApp::Create()
{
	if(!::RegisterClassEx(&m_wndClassEx)) return FALSE;

	_handle = ::CreateWindowEx(WS_EX_APPWINDOW | WS_EX_WINDOWEDGE,
							  AppName,
							  AppName,
							  WS_OVERLAPPEDWINDOW | WS_CLIPCHILDREN,   
							  CW_USEDEFAULT,
							  CW_USEDEFAULT,
							  CW_USEDEFAULT,
							  CW_USEDEFAULT,
							  NULL,
							  NULL,
							  m_hInstance,
							  NULL);
	if(!_handle) return FALSE;

	
	::ShowWindow(_handle, SW_SHOWDEFAULT);
	::UpdateWindow(_handle);

	m_pMainWnd = new WinGL( ViewManagerInst(), _handle );
	FV_ASSERT( m_pMainWnd != NULL );

	return TRUE;
}

void MainApp::Update()
{
	// Update fields with mesh parameters
	char buff[256];
	int nno = static_cast<int>(ModelCtrlInst().GetCurrentMesh()->GetNumNodes());
	sprintf(buff,"no: %d",nno);
	MainApp::SetStatusBarInfo(m_hStatusBar,0,0,buff);
	int ned = static_cast<int>(ModelCtrlInst().GetCurrentMesh()->GetNumEdges());
	sprintf(buff,"ed: %d",ned);
	MainApp::SetStatusBarInfo(m_hStatusBar,1,0,buff);
	int nfa = static_cast<int>(ModelCtrlInst().GetCurrentMesh()->GetNumFaces());
	sprintf(buff,"fa: %d",nfa);
	MainApp::SetStatusBarInfo(m_hStatusBar,2,0,buff);
	int nel = static_cast<int>(ModelCtrlInst().GetCurrentMesh()->GetNumElems());
	sprintf(buff,"el: %d",nel);
	MainApp::SetStatusBarInfo(m_hStatusBar,3,0,buff);
	
}

void MainApp::Destroy()
{
	if(::IsWindow(_handle)) ::PostMessage(_handle,WM_CLOSE,0,0);
}


int MainApp::Run()
{
	MSG msg;
    do 
	{ 
        if( PeekMessage( &msg, NULL, 0, 0, PM_REMOVE ) != 0 )
		{
			TranslateMessage( &msg );
			DispatchMessage( &msg );
		}
		else
		{
			// No messages 2 process?
			(GetInstance()->GetWindowGL())->Render();
		}
		 
	} while ( WM_QUIT != msg.message );

	if(GetInstance()->GetWindowGL() != NULL)
	{
		delete GetInstance()->GetWindowGL();
		GetInstance()->GetWindowGL() = NULL;
	}

	return (int) msg.wParam;
}

void MainApp::SetStatusBarInfo(HWND hwndSB, int part, unsigned style, char *fmt, ...)
{
    char tmpbuf[128];
    va_list argp;
    va_start(argp, fmt);
 
    _vsnprintf(tmpbuf, sizeof(tmpbuf), fmt, argp);
 
    va_end(argp);
 
    //cannot use PostMessage, as the panel type is not set correctly
	SendMessage(hwndSB, SB_SETTEXT, (WPARAM)part | style, (LPARAM)(tmpbuf));
}


LRESULT CALLBACK WndProcedure(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	//bool lresult = false;        // return value
	//static IController* lpctrl(NULL);
	//
	//if(uMsg == WM_NCCREATE)
	//{
	//	lpctrl = (IController*)(((CREATESTRUCT*)lParam)->lpCreateParams);
	//	lpctrl->Atach(hWnd);

	//	::SetWindowLongPtr(hWnd, GWLP_USERDATA, (LONG_PTR)lpctrl);

	//	return ::DefWindowProc(hWnd, uMsg, wParam, lParam);
	//} else {
	//	lpctrl = (IController*)::GetWindowLongPtr(hWnd, GWLP_USERDATA);
	//	if(lpctrl)
	//		lresult = lpctrl->Command(uMsg,wParam,lParam);
	//}

	//if(!lresult)
	//	return ::DefWindowProc(hWnd, uMsg, wParam, lParam);

	//return 0L;	
	int wmId, wmEvent;
    PAINTSTRUCT ps;
    HDC hdc;
	RECT rectGL;
	static WinGL*& GLWnd = MainApp::GetInstance()->GetWindowGL();
	switch (uMsg)
    {
	case WM_CREATE:
		{
			HWND& hStBr = MainApp::GetInstance()->GetStatusBar();
		    // Create status bar
			hStBr = ::CreateWindowEx(0, STATUSCLASSNAME, 0, WS_CHILD | WS_VISIBLE | SBARS_SIZEGRIP,
                                 0, 0, 0, 0, hWnd, (HMENU)IDC_STATUSBAR, ::GetModuleHandle(0), 0);

			int Statwidths[] = { 100,200,300,400,-1};
			::SendMessage(hStBr, SB_SETPARTS, (WPARAM)(sizeof(Statwidths)/sizeof(int)), (LPARAM)Statwidths);

		    //::SendMessage(hStBr, (UINT) SB_SETTEXT, (WPARAM)(INT) 1 | 0, (LPARAM) (LPSTR) TEXT("Ready to use"));

			MainApp::SetStatusBarInfo(hStBr,4,0,_T("Ready to use"));
		}
		break;
    case WM_COMMAND:
        wmId    = LOWORD(wParam);
        wmEvent = HIWORD(wParam);
        // Parse the menu selections:
        switch(wmId)
		{
		case IDM_HELP_ABOUT:
			::DialogBox(GetModuleHandle(0), MAKEINTRESOURCE(IDD_DIALOGABOUT), hWnd, About);
				//::DialogBox((HINSTANCE)::GetWindowLongPtr(_hwnd, GWLP_HINSTANCE), MAKEINTRESOURCE(IDD_DIALOGABOUT), _hwnd, About);
            break;
		case IDM_CONFIGURE_CHANGEAPROXIMATIONMODULE:
			{
				nodule_params params = { ModelCtrlInst().GetMeshModuleType(), ModelCtrlInst().GetApproximationType() };
				DialogWindow*& dlgPtr= MainApp::GetInstance()->GetChildDlg(); 
				dlgPtr = new DialogApprox((void *)&params, GetModuleHandle(0), IDD_DIALOG_APPROX, hWnd, ApproxDlgProc);
				//dlgPtr->Create();
				FV_FREE_PTR(dlgPtr);
			} break;
		case IDM_FILE_EXIT:
			// cleanup
				::DestroyWindow(hWnd);
				break;
		case IDM_FILE_OPENMESH:
			{
				OPENFILENAME ofn;       // common dialog box structure
				char szFile[260] = "";       // buffer for file name  
				ZeroMemory( & ofn, sizeof( ofn ) );
				ofn.lStructSize = sizeof( ofn );
				ofn.lpstrFilter = "Input mesh data(*.dat; *.dmp)\0*.dat;0*.dmp\0All files\0*.*\0";
				ofn.lpstrTitle = TEXT("Opem mesh input file");
				ofn.nMaxFile = MAX_PATH;
				ofn.lpstrFile = szFile;
				ofn.lpstrDefExt = "dat";
				ofn.Flags = OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
				if( GetOpenFileName( & ofn ) == TRUE )
				{
					//WczytajPlik( sNazwaPliku, hEdit );
				}
			
			} break;
		case IDM_FILE_OPENFIELD:
			{
				OPENFILENAME ofn;       // common dialog box structure
				char szFile[260] = "";       // buffer for file name  
				ZeroMemory( & ofn, sizeof( ofn ) );
				ofn.lStructSize = sizeof( ofn );
				ofn.lpstrFilter = "Input field data(*.dmp; *.dat)\0*.dmp;\0*.dat\0All files\0*.*\0";
				ofn.lpstrTitle = TEXT("Opem field input file");
				ofn.nMaxFile = MAX_PATH;
				ofn.lpstrFile = szFile;
				ofn.lpstrDefExt = "dmp";
				ofn.Flags = OFN_FILEMUSTEXIST | OFN_HIDEREADONLY;
				if( GetOpenFileName( & ofn ) == TRUE )
				{
					//WczytajPlik( sNazwaPliku, hEdit );
				}
			
			} break;
		case IDM_RENDER_DUMPSCREEN:
			{
				OPENFILENAME ofn;       // common dialog box structure
				char szFile[260] = "";       // buffer for file name  
				ZeroMemory( & ofn, sizeof( ofn ) );
				ofn.lStructSize = sizeof( ofn );
				ofn.lpstrFilter = "Ouput image file (*.png)\0*.png;\0All files\0*.*\0";
				ofn.lpstrTitle = TEXT("Save screen to file");
				ofn.nMaxFile = MAX_PATH;
				ofn.lpstrFile = szFile;
				ofn.lpstrDefExt = "png";
				ofn.Flags = OFN_HIDEREADONLY | OFN_LONGNAMES;
				if( GetSaveFileName( & ofn ) == TRUE )
				{
					//WczytajPlik( sNazwaPliku, hEdit );
				}
			} break;
		case IDM_FILE_RELOAD:
		case IDM_FILE_REFRESH:
		case IDM_FILE_CLEAR:
		case IDM_VIEW_ORTOGRAPHIC:
		case IDM_VIEW_PERSPECTIVE:
		case IDM_VIEW_TOP:
		case IDM_VIEW_BOTTOM:
		case IDM_VIEW_FRONT:
		case IDM_VIEW_BACK:
		case IDM_VIEW_LEFT:
		case IDM_VIEW_RIGHT:
		case IDM_VIEW_DEFAULT:
		case IDM_VIEW_FULL:
		case IDM_VIEW_FAST:
		case IDM_VIEW_BOUNDINGBOX:
		case IDM_VIEW_NEWVIEW:
		case IDM_VIEW_NEXTVIEW:
		case IDM_VIEW_PREVIOUSVIEW:
		case IDM_CONFIGURE_AXES:
		case IDM_CONFIGURE_GRID:
		case IDM_CONFIGURE_LEGEND:
		case IDM_CONFIGURE_BACKGROUNDCOLOR:
		case IDM_CONFIGURE_LIGHT:
		case IDM_CONFIGURE_EDITLEGEND:
		case IDM_CONFIGURE_RESET:
		case IDM_CONFIGURE_SAVECONFIGURATION:
		case IDM_RENDER_SETSOLUTIONFORMULA:
		case IDM_RENDER_SETSECTION:
		case IDM_RENDER_WIREFRAME:
		case IDM_RENDER_SHADED:
		case IDM_RENDER_CONTOURLINES:
		case IDM_RENDER_TRIM:
		case ID_VIEW_DUMPCURRENTVIEW:
		case IDM_VIEW_DUMPCURRENTVIEW:
		case IDM_VIEW_DUMPALLVIEWS:
			GLWnd->OnMenu(uMsg,wParam,lParam);
			break;
		default:
			return ::DefWindowProc(hWnd,uMsg,wParam,lParam);
		}
		break;
	case WM_CLOSE:
		IO::log("otrzymalem komunokat");
		// do cleanup
		::DestroyWindow(hWnd); //?
		break;
	//case WM_SIZE:
	//	//SetWindowText(hWnd,_T("Zmiana rozmiaru"));
	//	::SendMessage(GetInstance().GetStatusBar(),WM_SIZE,0,0);
	//	break;
    case WM_ERASEBKGND:
        // do nothing
        break;
    case WM_PAINT:
        hdc = BeginPaint(hWnd, &ps);
        // TODO: Add any drawing code here...
        EndPaint(hWnd, &ps);
        break;
    case WM_DESTROY:
		{
          ::PostQuitMessage(0);
		 // MainApp*& pApp = MainApp::GetInstance();
		 // if( pApp != NULL )
		 // {
			//delete pApp;
			//pApp = NULL;
		 // }
		}
        break;
    default:
        return ::DefWindowProc(hWnd, uMsg, wParam, lParam);
    }
    return 0L;
}


// Message handler for about box.
INT_PTR CALLBACK About(HWND hDlg, UINT msg, WPARAM wParam, LPARAM lParam)
{
    UNREFERENCED_PARAMETER(lParam);
    switch (msg)
    {
    case WM_INITDIALOG:
		
        return (INT_PTR)TRUE;

    case WM_COMMAND:
        if (LOWORD(wParam) == IDOK)
        {
            ::EndDialog(hDlg, IDOK);
            return (INT_PTR)TRUE;
        }
        break;
    }
    return (INT_PTR)FALSE;
}

// message handler for Approximation dialog
INT_PTR CALLBACK ApproxDlgProc(HWND hDlg, UINT msg, WPARAM wParam, LPARAM lParam)
{
	UNREFERENCED_PARAMETER(lParam);
    switch (msg)
    {
    case WM_INITDIALOG:
		MainApp::GetInstance()->GetChildDlg()->InitControls();
        return (INT_PTR)TRUE;

    case WM_COMMAND:
		{
			switch( LOWORD(wParam) )
            {
			case IDOK:
				// Do test for correct values of controls
			case IDCANCEL:
				::EndDialog(hDlg, IDOK);
				return (INT_PTR)TRUE;
			}
		
		} break;
    }
    return (INT_PTR)FALSE;
}

}
}

