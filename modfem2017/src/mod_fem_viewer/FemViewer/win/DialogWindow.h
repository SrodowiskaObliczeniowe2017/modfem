#ifndef DIALOG_WINDOW_H_
#define DIALOG_WINDOW_H_

#include"Common.h"
#include"Controls.h"

namespace FemViewer {
namespace Win		{



class DialogWindow {
public:
	DialogWindow(HINSTANCE hInst, WORD id, HWND hParent, DLGPROC Proc);
	virtual ~DialogWindow();
	virtual BOOL Create();
	virtual void InitControls(){}
	VOID Show(int cmdShow=SW_SHOWDEFAULT);
	HWND GetHandle() { return m_handle; } 
	VOID SetPosition(int x, int y) { this->m_posx = x; m_posy = y; }
	
protected:
	HWND m_handle;
	HWND m_parent;
	WORD m_id;
	int  m_posx, m_posy;
	int m_width, m_height;
	HINSTANCE m_hInstance;
	DLGPROC m_pfnDialogMsg;
};


class DialogApprox : public DialogWindow {
public:
	DialogApprox(void * params, HINSTANCE hInst, WORD id, HWND hParent, DLGPROC Proc)
    : DialogWindow(hInst,id,hParent,Proc)
	{ 
		InitControls(params);
	}
	~DialogApprox() {}

	BOOL Create(void * params);
	void InitControls(void * params);
private:
	

	// Mesh radio-butons
	RadioButton m_rbtnMesh2D;
	RadioButton m_rbtnMesh3D;
	RadioButton m_rbtnMeshTR;
	// Approximation ardio-buttons
	RadioButton m_rbtnApproxDG;
	RadioButton m_rbtnApproxSTD;
};

} // end namespace Win
} // end namespace FemViewer

#endif /* DIALOG_WINDOW_H_
*/