#include "wx/wxprec.h"

#ifdef __BORLANDC__
#pragma hdrstop
#endif

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#include "wx/defs.h"

#define wxAPPLY					0x00000020
#define wxCLOSE					0x00000040
#include "wxScreenSettings.h"
#include "wx/bookctrl.h"

#ifdef __WXDEBUG__
#include<iostream>
using namespace std;
#endif

const wxScreenSetts::BtnInfo wxScreenSetts::_btnInfo[] = 
{
	{ wxScreenSetts::IDS_AXES,		"&Axes on" },
	{ wxScreenSetts::IDS_GRID,		"&Grid on" },
	{ wxScreenSetts::IDS_LEGEND,	"&Legend on" },
	{ wxScreenSetts::IDS_MESH,		"&Mesh on" },
	{ wxScreenSetts::IDS_FIELD,		"&Field on" },
};

BEGIN_EVENT_TABLE(wxScreenSetts, wxDialog)
    EVT_BUTTON(wxID_APPLY, wxScreenSetts::OnBtnApply)
    EVT_BUTTON(wxID_CLOSE, wxScreenSetts::OnBtnClose)
END_EVENT_TABLE()

//BEGIN_EVENT_TABLE(wxScreenSetts, wxDialog)
//  EVT_BUTTON    ( wxID_OK, wxScreenSetts::OnBtnOK )
//  //EVT_BUTTON    ( wxID_HELP, wxMgrViewDialogLegend::OnHelp )
//  //EVT_RADIOBOX  (ID_RADIOBOX, wxMgrViewDialogLegend::OnRadioColGrad)
//
//  //EVT_CHECKBOX ( MenuColorCheckBox, wxMgrViewDialogLegend::OnMenuColorCheckBox)
//
//END_EVENT_TABLE()

//IMPLEMENT_CLASS(wxScreenSetts, wxPropertySheetDialog)

//BEGIN_EVENT_TABLE(wxScreenSetts, wxPropertySheetDialog)
//END_EVENT_TABLE()

wxScreenSetts::wxScreenSetts(ViewManager &vmgr, wxWindow *parent)
: wxDialog(parent, wxID_ANY, wxT("Settings"), wxDefaultPosition, wxDefaultSize,
		   wxDEFAULT_DIALOG_STYLE | wxRESIZE_BORDER),
  _viewMgr(vmgr)//, _imageList(NULL)
{
#ifdef __WXDEBUG__
	cout<< "Constructor: wxScreenSetts\n";
#endif

	wxSizer * const sizerTop = new wxBoxSizer(wxVERTICAL);

	wxSizer * const sizerHor = new wxBoxSizer(wxHORIZONTAL);
    // this sizer allows to configure the messages shown in the message box
    wxSizer * const sizerScreen = new wxStaticBoxSizer(wxVERTICAL, this, wxString("&Screen",wxConvUTF8));

	int flag = true;

	for(int n= 0; n< btnMax; n++)
	{
		buttons[n] = new wxCheckBox(this, wxID_ANY, wxString(_btnInfo[n].name,wxConvUTF8));
		sizerScreen->Add(buttons[n], flag ? wxSizerFlags().Centre().Left().Border(3) : wxSizerFlags().Centre().Expand().Border(6));
		flag = false;
	}

	sizerHor->Add(sizerScreen, wxSizerFlags(3).Center().Border());


    sizerTop->Add(sizerHor, wxSizerFlags(1).Expand().Border());


    // this one is for configuring the buttons
    //wxFlexGridSizer * const sizerBtns = new wxFlexGridSizer(2, 5, 5);
    //sizerBtns->AddGrowableCol(1);

    //sizerBtns->Add(new wxStaticText(this, wxID_ANY, "Button(s)"));
    //sizerBtns->Add(new wxStaticText(this, wxID_ANY, "Custom label"));

    //for ( int n = 0; n < Btn_Max; n++ )
    //{
    //    m_buttons[n] = new wxCheckBox(this, wxID_ANY, ms_btnInfo[n].name);
    //    sizerBtns->Add(m_buttons[n], wxSizerFlags().Centre().Left());

    //    m_labels[n] = new wxTextCtrl(this, wxID_ANY);
    //    sizerBtns->Add(m_labels[n], wxSizerFlags(1).Centre().Expand());

    //    m_labels[n]->Connect(wxEVT_UPDATE_UI,
    //                         wxUpdateUIEventHandler(
    //                             TestMessageBoxDialog::OnUpdateLabelUI),
    //                         NULL,
    //                         this);
    //}

    //wxSizer * const
    //    sizerBtnsBox = new wxStaticBoxSizer(wxVERTICAL, this, "&Buttons");
    //sizerBtnsBox->Add(sizerBtns, wxSizerFlags(1).Expand());
    //sizerTop->Add(sizerBtnsBox, wxSizerFlags().Expand().Border());


    //// icon choice
    //const wxString icons[] = {
    //    "&Information", "&Question", "&Warning", "&Error"
    //};

    //m_icons = new wxRadioBox(this, wxID_ANY, "&Icon:",
    //                         wxDefaultPosition, wxDefaultSize,
    //                         WXSIZEOF(icons), icons);
    //sizerTop->Add(m_icons, wxSizerFlags().Expand().Border());


    //// miscellaneous other stuff
    //wxSizer * const
    //    sizerFlags = new wxStaticBoxSizer(wxHORIZONTAL, this, "&Other flags");

    //m_chkNoDefault = new wxCheckBox(this, wxID_ANY, "Make \"No\" &default");
    //m_chkNoDefault->Connect(wxEVT_UPDATE_UI,
    //                        wxUpdateUIEventHandler(
    //                            TestMessageBoxDialog::OnUpdateNoDefaultUI),
    //                        NULL,
    //                        this);
    //sizerFlags->Add(m_chkNoDefault, wxSizerFlags(1).Border());

    //m_chkCentre = new wxCheckBox(this, wxID_ANY, "Centre on &parent");
    //sizerFlags->Add(m_chkCentre, wxSizerFlags(1).Border());

    //sizerTop->Add(sizerFlags, wxSizerFlags().Expand().Border());

    // finally buttons to show the resulting message box and close this dialog
    sizerTop->Add(CreateStdDialogButtonSizer(wxAPPLY | wxCLOSE),
                  wxSizerFlags().Right().Border());

    SetSizerAndFit(sizerTop);

	buttons[btnAxes]->SetValue(_viewMgr.GetSettings()->bIsAxesOn ? true : false);
	buttons[btnGrid]->SetValue(_viewMgr.GetSettings()->bIsGridOn ? true : false);
	buttons[btnLegend]->SetValue(_viewMgr.GetSettings()->bIsLegendOn ? true : false);
	buttons[btnMesh]->SetValue(_viewMgr.GetSettings()->bEdgeOn ? true : false);
	buttons[btnField]->SetValue(_viewMgr.GetSettings()->bShadingOn ? true : false);
}


wxScreenSetts::~wxScreenSetts()
{
	//delete _imageList;
}

//wxPanel* wxScreenSetts::CreateGeneralSettingsPage(wxWindow *parent)
//{
//	wxPanel* panel = new wxPanel(parent, wxID_ANY);
//
//	wxBoxSizer *topSizer = new wxBoxSizer( wxVERTICAL );
//    wxBoxSizer *item0 = new wxBoxSizer( wxVERTICAL );
//
//    // General setting tab
//
//	// Axes On/Off
//    wxBoxSizer* itemSizer3 = new wxBoxSizer( wxHORIZONTAL );
//    wxCheckBox* checkBox3 = new wxCheckBox(panel, ID_DRAW_AXES, _("&Axes On/Off"), wxDefaultPosition, wxDefaultSize);
//	checkBox3->SetValue( _viewMgr.GetSettings()->bIsAxesOn ? true : false );
//    itemSizer3->Add(checkBox3, 0, wxALL|wxALIGN_CENTER_VERTICAL, 5);
//    item0->Add(itemSizer3, 0, wxGROW|wxALL, 0);
//
//	// Grid On/Off
//	wxBoxSizer* itemSizer2 = new wxBoxSizer( wxHORIZONTAL );
//    wxCheckBox* checkBox2 = new wxCheckBox(panel, ID_DRAW_AXES, _("&Grid On/Off"), wxDefaultPosition, wxDefaultSize);
//    itemSizer2->Add(checkBox2, 0, wxALL|wxALIGN_CENTER_VERTICAL, 5);
//	checkBox2->SetValue( _viewMgr.GetSettings()->bIsGridOn ? true : false );
//    item0->Add(itemSizer2, 0, wxGROW|wxALL, 0);
//
//	// Legend On/Off
//	wxBoxSizer* itemSizer1 = new wxBoxSizer( wxHORIZONTAL );
//    wxCheckBox* checkBox1 = new wxCheckBox(panel, ID_DRAW_AXES, _("&Legend On/Off"), wxDefaultPosition, wxDefaultSize);
//    itemSizer1->Add(checkBox1, 0, wxALL|wxALIGN_CENTER_VERTICAL, 5);
//	//checkBox2->SetValue( _viewMgr.GetSettings()->bIsGridOn ? true : false );
//    item0->Add(itemSizer1, 0, wxGROW|wxALL, 0);
//
//
///*
//    // 
//
//    wxString autoSaveLabel = _("&Auto-save every");
//    wxString minsLabel = _("mins");
//
//    wxBoxSizer* itemSizer12 = new wxBoxSizer( wxHORIZONTAL );
//    wxCheckBox* checkBox12 = new wxCheckBox(panel, ID_AUTO_SAVE, autoSaveLabel, wxDefaultPosition, wxDefaultSize);
//
////#if wxUSE_SPINCTRL
//    wxSpinCtrl* spinCtrl12 = new wxSpinCtrl(panel, ID_AUTO_SAVE_MINS, wxEmptyString,
//        wxDefaultPosition, wxSize(40, wxDefaultCoord), wxSP_ARROW_KEYS, 1, 60, 1);
////#endif
//
//    itemSizer12->Add(checkBox12, 0, wxALL|wxALIGN_CENTER_VERTICAL, 5);
////#if wxUSE_SPINCTRL
//    itemSizer12->Add(spinCtrl12, 0, wxALL|wxALIGN_CENTER_VERTICAL, 5);
////#endif
//    itemSizer12->Add(new wxStaticText(panel, wxID_STATIC, minsLabel), 0, wxALL|wxALIGN_CENTER_VERTICAL, 5);
//    item0->Add(itemSizer12, 0, wxGROW|wxALL, 0);
//
//    //// TOOLTIPS
//
//    wxBoxSizer* itemSizer8 = new wxBoxSizer( wxHORIZONTAL );
//    wxCheckBox* checkBox6 = new wxCheckBox(panel, ID_SHOW_TOOLTIPS, _("Show &tooltips"), wxDefaultPosition, wxDefaultSize);
//    itemSizer8->Add(checkBox6, 0, wxALL|wxALIGN_CENTER_VERTICAL, 5);
//    item0->Add(itemSizer8, 0, wxGROW|wxALL, 0);
//*/
//    topSizer->Add( item0, 1, wxGROW|wxALIGN_CENTRE|wxALL, 5 );
//
//    panel->SetSizerAndFit(topSizer);
//	return panel;
//}
//
//wxPanel* wxScreenSetts::CreateMeshSettingsPage(wxWindow *parent)
//{
//	wxPanel* panel = new wxPanel(parent, wxID_ANY);
//	return panel;
//}
//
//wxPanel* wxScreenSetts::CreateFieldSettingsPage(wxWindow *parent)
//{
//	wxPanel* panel = new wxPanel(parent, wxID_ANY);
//	return panel;
//}

void wxScreenSetts::OnBtnApply(wxCommandEvent& event)
{
	_viewMgr.GetSettings()->bIsAxesOn = buttons[btnAxes]->IsChecked() ? true : false;
	_viewMgr.GetSettings()->bIsGridOn = buttons[btnGrid]->IsChecked() ? true : false;
	_viewMgr.GetSettings()->bIsLegendOn = buttons[btnLegend]->IsChecked() ? true : false;
	_viewMgr.GetSettings()->bEdgeOn = buttons[btnMesh]->IsChecked() ? true : false;
	_viewMgr.GetSettings()->bShadingOn = buttons[btnField]->IsChecked() ? true : false;
	


    //switch ( m_icons->GetSelection() )
    //{
    //    case 0: style |= wxICON_INFORMATION; break;
    //    case 1: style |= wxICON_QUESTION; break;
    //    case 2: style |= wxICON_WARNING; break;
    //    case 3: style |= wxICON_ERROR; break;
    //}

    //if( m_chkCentre->IsChecked() )
    //    style |= wxCENTRE;

    //if ( m_chkNoDefault->IsEnabled() && m_chkNoDefault->IsChecked() )
    //    style |= wxNO_DEFAULT;


    //wxMessageDialog dlg(this, m_textMsg->GetValue(), "Test Message Box",
    //                    style);
    //if ( !m_textExtMsg->IsEmpty() )
    //    dlg.SetExtendedMessage(m_textExtMsg->GetValue());

    //if ( style & wxYES_NO )
    //{
    //    if ( style & wxCANCEL )
    //    {
    //        dlg.SetYesNoCancelLabels(m_labels[Btn_Yes]->GetValue(),
    //                                 m_labels[Btn_No]->GetValue(),
    //                                 m_labels[Btn_Cancel]->GetValue());
    //    }
    //    else
    //    {
    //        dlg.SetYesNoLabels(m_labels[Btn_Yes]->GetValue(),
    //                           m_labels[Btn_No]->GetValue());
    //    }
    //}
    //else
    //{
    //    if ( style & wxCANCEL )
    //    {
    //        dlg.SetOKCancelLabels(m_labels[Btn_Ok]->GetValue(),
    //                              m_labels[Btn_Cancel]->GetValue());
    //    }
    //    else
    //    {
    //        dlg.SetOKLabel(m_labels[Btn_Ok]->GetValue());
    //    }
    //}

    //dlg.ShowModal();
}

void wxScreenSetts::OnBtnClose(wxCommandEvent& event)
{
	EndModal(wxID_CANCEL);
}
