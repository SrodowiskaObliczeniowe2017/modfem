#include "DialogWindow.h"
#include "resource.h"

namespace FemViewer {
namespace Win       {

DialogWindow::DialogWindow(HINSTANCE hInst,WORD id,HWND hParent, DLGPROC Proc)
: m_handle(0)
, m_parent(hParent)
, m_id(id)
, m_posx(0)
, m_posy(0)
, m_width(100)
, m_height(100)
, m_hInstance(hInst)
, m_pfnDialogMsg(Proc)
{
}

DialogWindow::~DialogWindow()
{
}

BOOL DialogWindow::Create()
{
	if(DialogBox(GetModuleHandle(0), MAKEINTRESOURCE(m_id), m_handle, m_pfnDialogMsg)>0) return TRUE;
	return FALSE;
}

BOOL DialogApprox::Create(void * params)
{
	BOOL result = TRUE;
	if(DialogWindow::Create() == IDOK);

	return result;
}

void DialogApprox::InitControls(void * params)
{
	m_rbtnMesh2D.set(this->m_handle, IDC_RADIO_MESH2D);
	m_rbtnMesh3D.set(this->m_handle, IDC_RADIO_MESH3D);
	m_rbtnMeshTR.set(this->m_handle, IDC_RADIO_MESHTR);
	m_rbtnMesh3D.check();
	m_rbtnApproxDG.set(this->m_handle, IDC_RADIO_APPROXDG);
	m_rbtnApproxSTD.set(this->m_handle, IDC_RADIO_APPROXSTD);
	m_rbtnApproxDG.check();
	BOOL result = DialogWindow::Create();

}

}
}