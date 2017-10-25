#ifndef _WX_SCREEN_SETTINGS_H_
#define _WX_SCREEN_SETTINGS_H_

#include <wx/wx.h>
//#include <wx/propdlg.h>
//#include <wx/bookctrl.h>
//#include <wx/defs.h>

//class ViewManager;
#include "../inc/ViewManager.h"
using namespace FemViewer;


class wxScreenSetts : public wxDialog
{
public:
	wxScreenSetts( ViewManager& vmgr, wxWindow* parent);
	~wxScreenSetts();

	//wxPanel* CreateGeneralSettingsPage(wxWindow* parent);
 //   wxPanel* CreateMeshSettingsPage(wxWindow* parent);
	//wxPanel* CreateFieldSettingsPage(wxWindow* parent);



private:
	void OnBtnApply(wxCommandEvent& event);
	void OnBtnClose(wxCommandEvent& event);


private:
	ViewManager& _viewMgr;
	/*wxImageList* _imageList;*/
	enum {
		IDS_AXES,
		IDS_GRID,
		IDS_LEGEND,
		IDS_MESH,
		IDS_FIELD,
	};

	enum {
		btnAxes,
		btnGrid,
		btnLegend,
		btnMesh,
		btnField,
		btnMax,
	};

	struct BtnInfo
    {
        int flag;
        const char *name;
    };

	static const BtnInfo _btnInfo[btnMax];

	wxCheckBox *buttons[btnMax];


	DECLARE_EVENT_TABLE()

    //wxDECLARE_NO_COPY_CLASS(wxScreenSetts)
	wxScreenSetts(const wxScreenSetts &);
	wxScreenSetts& operator=(const wxScreenSetts&);

};

#endif /* _WX_SCREEN_SETTINGS_H_
*/
