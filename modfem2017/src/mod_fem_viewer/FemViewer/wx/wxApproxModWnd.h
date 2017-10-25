#ifndef _APPROX_MOD_WND_H_
#define _APPROX_MOD_WND_H_

#include "ApproxModule.h"

class wxApproxModWnd : public wxDialog {

public:
	// Constructor
	wxApproxModWnd(wxWindow *parent,mod_approx& option);

	// Types
	enum {
		IDR_INTRENAL,
		IDR_EXTERNAL_DG_PRSIM,
		IDR_EXTERNAL_STD_PRISM,
		IDR_EXTERNAL_STD_HYBRID,
		IDR_MAX
	};
	
private:
	// Functions
	void OnEvent(wxCommandEvent & event);
	void OnOk(wxCommandEvent & event);
	void OnCancel(wxCommandEvent & event);

	

	struct BtnInfo {
		int id;
		wxString name;
	};

	static const BtnInfo btnInfo[IDR_MAX];

	// Members
	wxRadioButton * btnOpts[IDR_MAX];

	int &mode;
	mod_approx opt;
	
	DECLARE_EVENT_TABLE()
	//wxDECLARE_NO_COPY_CLASS(wxApproxModWnd);
	wxApproxModWnd(const wxApproxModWnd&);
	wxApproxModWnd& operator=(const wxApproxModWnd&);
};


#endif /* _APPROX_MOD_WND_H_
*/
