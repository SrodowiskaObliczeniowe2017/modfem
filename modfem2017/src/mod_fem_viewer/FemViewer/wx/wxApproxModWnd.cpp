// For compilers that support precompilation, includes "wx/wx.h".
#include "wx/wxprec.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#include <wx/defs.h>
#endif

#include "wxApproxModWnd.h"


/* static */
const wxApproxModWnd::BtnInfo wxApproxModWnd::btnInfo[] =
{
	{ wxApproxModWnd::IDR_INTRENAL, wxT("Use internal approximation (dg/std) for prism") },
	{ wxApproxModWnd::IDR_EXTERNAL_DG_PRSIM, wxT("Use external approximation (dg only) for prism meshes") },
	{ wxApproxModWnd::IDR_EXTERNAL_STD_PRISM, wxT("Use external approximation (std only) for prism meshes") },
	{ wxApproxModWnd::IDR_EXTERNAL_STD_HYBRID, wxT("Use external approximation (std only) for hybrid meshes") },
};

BEGIN_EVENT_TABLE(wxApproxModWnd, wxDialog)
	EVT_RADIOBUTTON(IDR_INTRENAL, wxApproxModWnd::OnEvent)
	EVT_RADIOBUTTON(IDR_EXTERNAL_DG_PRSIM, wxApproxModWnd::OnEvent)
    EVT_RADIOBUTTON(IDR_EXTERNAL_STD_PRISM, wxApproxModWnd::OnEvent)
	EVT_RADIOBUTTON(IDR_EXTERNAL_STD_HYBRID, wxApproxModWnd::OnEvent)
	EVT_BUTTON(wxID_OK, wxApproxModWnd::OnOk)
	EVT_BUTTON(wxID_CANCEL, wxApproxModWnd::OnCancel)
END_EVENT_TABLE()

wxApproxModWnd::wxApproxModWnd(wxWindow * parent, mod_approx & option)
: wxDialog(parent, wxID_ANY, wxString(wxT("Approximation module options")),
		   wxDefaultPosition, wxSize(200,100), wxDEFAULT_DIALOG_STYLE|wxRESIZE_BORDER),
		   mode((int&)option),
		   opt(option)
{
	wxASSERT((option >=0) && (option < IDR_MAX));

	// vertical sizer
	wxBoxSizer *const sizerTop = new wxBoxSizer(wxVERTICAL);
	sizerTop->AddSpacer(10);

	// static box
	wxStaticBox *const staticBox = new wxStaticBox(this, wxID_ANY, wxT("Choose approximation module"));

	// static sizer
	wxStaticBoxSizer *const sizer1 = new wxStaticBoxSizer(staticBox,wxVERTICAL);

	sizerTop->Add(sizer1, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);

	btnOpts[0] = new wxRadioButton(this, btnInfo[0].id, wxString(btnInfo[0].name), wxDefaultPosition, wxDefaultSize, wxRB_GROUP);
	btnOpts[1] = new wxRadioButton(this, btnInfo[1].id, wxString(btnInfo[1].name));
	btnOpts[2] = new wxRadioButton(this, btnInfo[2].id, wxString(btnInfo[2].name));
	btnOpts[3] = new wxRadioButton(this, btnInfo[3].id, wxString(btnInfo[3].name));

	sizer1->Add(btnOpts[0], 0, wxALL, 5);
	sizer1->AddSpacer(3);
	sizer1->Add(btnOpts[1], 0, wxALL, 5);
	sizer1->AddSpacer(3);
	sizer1->Add(btnOpts[2], 0, wxALL, 5);
	sizer1->AddSpacer(3);
	sizer1->Add(btnOpts[3], 0, wxALL, 5);
	sizer1->AddSpacer(3);
	//for(int i(0);i<IDR_MAX;i++)
	//{
	//	if(!i) 
	//		btnOpts[i] = new wxRadioButton(this, btnInfo[i].id, btnInfo[i].descr, wxDefaultPosition, wxDefaultSize, wxRB_GROUP);
	//	else
	//		btnOpts[i] = new wxRadioButton(this, btnInfo[i].id, btnInfo[i].descr);

	//	sizer1->Add(btnOpts[i], 0, wxALL, 5);
	//}

	// Update current options
	btnOpts[mode]->SetValue(true);

	//sizerTop->Add(sizer1, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);
	//sizerTop->AddStretchSpacer(0,0,5);

	// bottom sizer
	wxBoxSizer *const sizerBottom = new wxBoxSizer(wxHORIZONTAL);

	// Add OK & CANCEL btns
	sizerBottom->Add( new wxButton(this,wxID_OK,wxT("OK")), 0, wxALL, 10);
	sizerBottom->Add( new wxButton(this,wxID_CANCEL,wxT("Cancel")), 0, wxALL, 10);

	sizerTop->Add(sizerBottom, 0, wxALIGN_CENTER, 5);

	SetSizer(sizerTop);
	sizerTop->Fit(this);
	sizerTop->SetSizeHints(this);


}

void wxApproxModWnd::OnEvent(wxCommandEvent & event)
{
	opt = static_cast<mod_approx>(event.GetId());
}

void wxApproxModWnd::OnOk(wxCommandEvent& WXUNUSED(event))
{
	mode = opt;
	EndModal(wxID_OK);
}

void wxApproxModWnd::OnCancel(wxCommandEvent& WXUNUSED(event))
{
	EndModal(wxID_CANCEL);
}
