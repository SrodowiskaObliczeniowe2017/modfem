#include "wxLegendDialog.h"

#ifdef __WXDEBUG__
    #include <iostream>
    using namespace std;
#endif

BEGIN_EVENT_TABLE(wxLegendDialog, wxDialog)
    EVT_BUTTON( ResetPanels, wxLegendDialog::OnPanelsResetToDefaults )
    EVT_BUTTON( PanelNewRefresh, wxLegendDialog::OnPanelNewRefresh )
    EVT_BUTTON( PanelNewActive, wxLegendDialog::OnPanelNewActive )
  //EVT_SIZE
  //EVT_SIZING
END_EVENT_TABLE()

wxLegendDialog::wxLegendDialog(const wxString & title, FemViewer::Legend & legend)
: wxDialog
    (
        NULL, wxID_ANY, title,
        wxDefaultPosition,
        //size,
        //wxDefaultSize,
        wxSize(400,500),
        //style
        wxCAPTION | wxRESIZE_BORDER
    ),
    activeLegend(legend)
{

    //stworzenie nowego plutna
    newLegend = activeLegend;

    // (legendy oraz przyciski)
    //wxBoxSizer *vbTopTop = new wxBoxSizer(wxVERTICAL);
    //nadrzedny
    wxBoxSizer* topsizer = new wxBoxSizer(wxVERTICAL);//vbTopTop

    //this->SetSizer(topsizer);
    //topsizer->SetSizeHints(this);


    //glowny
    wxBoxSizer *mainsizer = new wxBoxSizer(wxHORIZONTAL);//hbTop
    //topsizer->Add(mainsizer);
    topsizer->Add(mainsizer, wxSizerFlags().Expand());


    //pobranie wielkosci mainsizer
    //odjecia tyle co na info
    //potem po polowie dla Tworzenia pozostalych

    wxSize size;
    size =  mainsizer->GetSize();

    int in;
    in = 40;

    size.GetWidth();
    size.GetHeight();

    //ramka aktywnej legendy
    wxStaticBox * activeTxtBox = new wxStaticBox (this, wxID_ANY, wxT("Current legend") );
    wxStaticBoxSizer * vbActive = new wxStaticBoxSizer (activeTxtBox, wxVERTICAL);

    //mainsizer->Add(vbActive);
    mainsizer->Add( vbActive, wxSizerFlags(1).Expand() );

    //ranka nowej legendy
    wxStaticBox * newTxtBox = new wxStaticBox (this, wxID_ANY, wxT("New legend") );
    wxStaticBoxSizer * vbNew = new wxStaticBoxSizer (newTxtBox, wxVERTICAL);

    //mainsizer->Add(vbNew);
    mainsizer->Add( vbNew, wxSizerFlags(1).Expand() );



    //m_canvas_active = new wxLegendCanvas( this, activeLegend, wxSize( (size.GetWidth()-in)/2 , size.GetHeight()-10) );
    m_canvas_active = new wxLegendCanvas( this, activeLegend, wxDefaultSize );

    //m_canvas_new = new wxLegendCanvas( this, newLegend, wxSize( (size.GetWidth()-in)/2 , size.GetHeight()-10) );
    m_canvas_new = new wxLegendCanvas( this, newLegend, wxDefaultSize );


    vbActive->Add(m_canvas_active, wxSizerFlags(1).Expand() );

    //vbNew->Add(m_canvas_new, 0, wxALIGN_LEFT, 5);
    vbNew->Add(m_canvas_new, wxSizerFlags(1).Expand());


    //ramka informacyjna
    wxBoxSizer *infosizer = new wxBoxSizer(wxVERTICAL);//vbInfo

    mainsizer->Add(infosizer, wxSizerFlags().Expand());


    //informacje aktywnej
    wxStaticBox * activeInfoTxtBox = new wxStaticBox (this, -1, wxT("Current settings") );
    wxStaticBoxSizer * vbActiveInfo = new wxStaticBoxSizer (activeInfoTxtBox, wxVERTICAL);

    //informcje nowej
    wxStaticBox * newInfoTxtBox = new wxStaticBox (this, -1, wxT("New settings") );
    wxStaticBoxSizer * vbNewInfo = new wxStaticBoxSizer (newInfoTxtBox, wxVERTICAL);


    //infosizer->Add(vbActiveInfo, 0, wxBORDER, 5);
    infosizer->Add(vbActiveInfo, wxSizerFlags().Border(0));
    //infosizer->Add(vbNewInfo, 0, wxBORDER, 5);
    infosizer->Add(vbNewInfo, wxSizerFlags().Border(0));

    //************************************************
    //  Poczaek info
    //************************************************


    wxString choices2[] = {wxT("Color"), wxT("Gradient")};

    wxStaticText *stActiveMin = new wxStaticText(this, wxID_ANY, _T("Min"));
    wxStaticText *stActiveMax = new wxStaticText(this, wxID_ANY, _T("Max"));
    wxStaticText *stActiveStart = new wxStaticText(this, wxID_ANY, _T("Start"));
    wxStaticText *stActiveEnd = new wxStaticText(this, wxID_ANY, _T("End"));
    wxStaticText *stActiveRengNr = new wxStaticText(this, wxID_ANY, wxT("L."));
    wxStaticText *stActiveLegType = new wxStaticText(this, wxID_ANY, wxT("Type"));
    wxStaticText *stActiveOutBlck = new wxStaticText(this, wxID_ANY, wxT("Color black outside"));


    wxStaticText *stNewMin = new wxStaticText(this, wxID_ANY, wxT("Min"));
    wxStaticText *stNewMax = new wxStaticText(this, wxID_ANY, wxT("Max"));
    wxStaticText *stNewStart = new wxStaticText(this, wxID_ANY, wxT("Start"));
    wxStaticText *stNewEnd = new wxStaticText(this, wxID_ANY, wxT("End"));
    wxStaticText *stNewRengNr = new wxStaticText(this, wxID_ANY, wxT("L."));
    wxStaticText *stNewLegType = new wxStaticText(this, wxID_ANY, wxT("Type"));
    wxStaticText *stNewOutBlck = new wxStaticText(this, wxID_ANY, wxT("Color black outside"));


        //min
    tcActiveMin = new wxTextCtrl(this, wxID_ANY);
        //max
    tcActiveMax = new wxTextCtrl(this, wxID_ANY);
        //start
    tcActiveStart = new wxTextCtrl(this, wxID_ANY);
        //end
    tcActiveEnd = new wxTextCtrl(this, wxID_ANY);
        //reng nr
    //tcActiveRengNr = new wxTextCtrl(this, -1);

    scActiveRengNr = new wxSpinCtrl( this, ID_SPINCTRL, wxEmptyString );
    //scActiveRengNr = new wxSpinCtrl( panel, ID_SPINCTRL, wxEmptyString, wxPoint(200, 160), wxSize(80, wxDefaultCoord) );
    //scActiveRengNr->SetRange(-10,30);
    //scActiveRengNr->SetValue(15);

        //leg type
    //wxTextCtrl *tcActiveLegType = new wxTextCtrl(this, -1);
    rbActiveLegType = new wxRadioBox( this, ID_RADIOBOX, _T(""), wxDefaultPosition, wxDefaultSize, WXSIZEOF(choices2), choices2, 1, wxRA_SPECIFY_COLS );
    //wxRA_SPECIFY_ROWS
        //out black
    //wxTextCtrl *tcActiveOutBlck = new wxTextCtrl(this, -1);
    cbActiveOutBlck = new wxCheckBox( this, MenuColorCheckBox, wxT(""), wxDefaultPosition, wxDefaultSize );
    cbActiveOutBlck->SetValue( legend.GetOutsideBlack() );

    tcActiveMin->Disable();
    tcActiveMax->Disable();
    tcActiveStart->Disable();
    tcActiveEnd->Disable();
    scActiveRengNr->Disable();
    rbActiveLegType->Disable();
    cbActiveOutBlck->Disable();


    //do testow formy wylaczyc





       //min
    //tcNewMin = new wxTextCtrl(this, -1);
    tcNewMin = new wxTextCtrl(this, -1, wxT(""), wxDefaultPosition, wxDefaultSize, 0,
        wxTextValidator(
        //wxFILTER_ALPHA
        wxFILTER_NUMERIC
        )
    );



        //mx
    tcNewMax = new wxTextCtrl(this, wxID_ANY);
        //start
    tcNewStart = new wxTextCtrl(this, wxID_ANY);
        //end
    tcNewEnd = new wxTextCtrl(this, wxID_ANY);
        //reng nr
    //tcNewRengNr = new wxTextCtrl(this, -1);
    scNewRengNr = new wxSpinCtrl( this, ID_SPINCTRL, wxEmptyString );

    scNewRengNr->SetRange(1, 50);
        //leg type

    //wxTextCtrl *tcNewLegType = new wxTextCtrl(this, -1);
    rbNewLegType = new wxRadioBox( this, ID_RADIOBOX, _T(""), wxDefaultPosition, wxDefaultSize, WXSIZEOF(choices2), choices2, 1, wxRA_SPECIFY_COLS );

        //out black

    //wxTextCtrl *tcNewOutBlck = new wxTextCtrl(this, -1);
    cbNewOutBlck = new wxCheckBox( this, MenuColorCheckBox, wxT(""), wxDefaultPosition, wxDefaultSize );

    //stworzenie wnetrza
        //aktywnego
    //wxButton *btnActiveInfo = new wxButton(this, wxID_EXIT, wxT("ActiveInfo"), wxDefaultPosition, wxSize(140, 199));
    //wxFlexGridSizer *fgActiveInfo = new wxFlexGridSizer(6, 2, 9, 25);
    wxFlexGridSizer *fgActiveInfo = new wxFlexGridSizer(6, 2, 1, 1);
    //wxFlexGridSizer *fgActiveInfo = new wxFlexGridSizer(6, 2, 5, 20);

    fgActiveInfo->Add(stActiveMin, wxSizerFlags().Align(wxALIGN_CENTER) );
    fgActiveInfo->Add(tcActiveMin, wxSizerFlags().Expand());

    fgActiveInfo->Add(stActiveMax, wxSizerFlags().Align(wxALIGN_CENTER));
    fgActiveInfo->Add(tcActiveMax, wxSizerFlags().Expand());

    fgActiveInfo->Add(stActiveStart, wxSizerFlags().Align(wxALIGN_CENTER));
    fgActiveInfo->Add(tcActiveStart, wxSizerFlags().Expand());

    fgActiveInfo->Add(stActiveEnd, wxSizerFlags().Align(wxALIGN_CENTER));
    fgActiveInfo->Add(tcActiveEnd, wxSizerFlags().Expand());

    fgActiveInfo->Add(stActiveRengNr, wxSizerFlags().Align(wxALIGN_CENTER));
    //fgActiveInfo->Add(tcActiveRengNr, 1, wxEXPAND);
    fgActiveInfo->Add(scActiveRengNr, wxSizerFlags().Expand());


    fgActiveInfo->Add(stActiveLegType, wxSizerFlags().Align(wxALIGN_CENTER));
    //fgActiveInfo->Add(tcActiveLegType, 1, wxEXPAND);
    fgActiveInfo->Add(rbActiveLegType, wxSizerFlags().Expand());


    ////fgActiveInfo->Add(tcActiveOutBlck, 1, wxEXPAND);
    //fgActiveInfo->Add(stActiveOutBlck, 1, wxALIGN_CENTER_VERTICAL );
    //fgActiveInfo->Add(cbActiveOutBlck, 1, wxEXPAND);

    fgActiveInfo->Add(cbActiveOutBlck, 1, wxALIGN_CENTER_VERTICAL);
    fgActiveInfo->Add(stActiveOutBlck, 1, wxEXPAND );

        //nowego
    //wxButton *btnNewInfo = new wxButton(this, wxID_EXIT, wxT("NewInfo"), wxDefaultPosition, wxSize(140, 199));

    //wxFlexGridSizer *fgNewInfo = new wxFlexGridSizer(6, 2, 5, 20);
    wxFlexGridSizer *fgNewInfo = new wxFlexGridSizer(6, 2, 1, 1);

    fgNewInfo->Add(stNewMin, 1, wxALIGN_CENTER);
    fgNewInfo->Add(tcNewMin, 1, wxEXPAND);

    fgNewInfo->Add(stNewMax, 1, wxALIGN_CENTER);
    fgNewInfo->Add(tcNewMax, 1, wxEXPAND);


    fgNewInfo->Add(stNewStart, 1, wxALIGN_CENTER);
    fgNewInfo->Add(tcNewStart, 1, wxEXPAND);

    fgNewInfo->Add(stNewEnd, 1, wxALIGN_CENTER);
    fgNewInfo->Add(tcNewEnd, 1, wxEXPAND);

    fgNewInfo->Add(stNewRengNr, 1, wxALIGN_CENTER);
    //fgNewInfo->Add(tcNewRengNr, 1, wxEXPAND);
    fgNewInfo->Add(scNewRengNr, 1, wxEXPAND);


    fgNewInfo->Add(stNewLegType, 1, wxALIGN_CENTER);
    fgNewInfo->Add(rbNewLegType, 1, wxEXPAND);

    fgNewInfo->Add(cbNewOutBlck, 1, wxALIGN_CENTER);
    fgNewInfo->Add(stNewOutBlck, 1, wxEXPAND);




    tcNewMin->Disable();
    tcNewMax->Disable();


    //wypelnienie wnetrza
        //altywnego
    //vbActiveInfo->Add(btnActiveInfo, 0, wxBORDER, 5);
    vbActiveInfo->Add(fgActiveInfo, 0, wxBORDER, 5);
        //nowego
    //vbNewInfo->Add(btnNewInfo, 0, wxBORDER, 5);
    vbNewInfo->Add(fgNewInfo, 0, wxBORDER, 5);


    //************************************************
    //
    //************************************************
    btnReset   = new wxButton(this, ResetPanels, wxT("Reset"));
    btnRefresh = new wxButton(this, PanelNewRefresh, wxT("Refresh"));
    btnActive  = new wxButton(this, PanelNewActive, wxT("Active"));

    wxBoxSizer *hb3Buttons = new wxBoxSizer(wxHORIZONTAL);

    //hb2Buttons->Add(btnRefresh);
    //hb2Buttons->Add(btnActive);
    hb3Buttons->Add(btnReset,  wxSizerFlags(1).Expand());
    hb3Buttons->Add(btnRefresh, wxSizerFlags(1).Expand());
    hb3Buttons->Add(btnActive, wxSizerFlags(1).Expand());

    //infosizer->Add(hb2Buttons, 0, wxBORDER, 10);
    infosizer->Add(hb3Buttons, wxSizerFlags().Border(0).Expand());


/*
    btnRefresh = new wxButton(this, PanelNewRefresh, wxT("Refresh"));
    btnActive =  new wxButton(this, PanelNewActive, wxT("Active"));

    wxBoxSizer *hb2Buttons = new wxBoxSizer(wxHORIZONTAL);

    hb2Buttons->Add(btnRefresh);
    hb2Buttons->Add(btnActive);

    vbInfo->Add(hb2Buttons, 0, wxBORDER, 10);

*/

    //************************************************
    //
    //************************************************

    wxSizer * sizerOkHelp = CreateButtonSizer(wxOK | wxHELP);
    //topsizer->Add(sizerOkHelp, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);
    topsizer->Add(sizerOkHelp, wxSizerFlags().Border(5).Bottom());

    this->SetSizer(topsizer);
    topsizer->SetSizeHints(this);
    //topsizer->SetSizerAndFit(this);

    //ustawienie wartosci w aktive
    SetActiveValues(activeLegend);
    SetNewValues(newLegend);

}

wxLegendDialog::wxLegendDialog(const wxString& title, long style, FemViewer::Legend& legend, const wxSize& size)
: wxDialog(NULL, -1, title, wxDefaultPosition, size, style), activeLegend(legend)
{

    //hb - horizontal box sizer
    //vb - vertical box sizer

    //fg - wxFlexGridSizer

    //tc - wxTextCtrl

    //st, stx - wxStaticText


    //btn -  button (wxButton)

    //rb - wxRadioBox
    //cb - wxCheckBox

    //sc, sctr - -wxSpinCtrl


    wxBoxSizer *vbTopTop = new wxBoxSizer(wxVERTICAL);

    //horyzontalne - ogolne (legenda, legenda, informacje)
    wxBoxSizer *hbTop = new wxBoxSizer(wxHORIZONTAL);

    wxBoxSizer *vbInfo = new wxBoxSizer(wxVERTICAL);


    wxStaticBox * activeTxtBox = new wxStaticBox (this, -1, wxT("Active legend settings") );
    wxStaticBoxSizer * vbActive = new wxStaticBoxSizer (activeTxtBox, wxVERTICAL);

    wxStaticBox * newTxtBox = new wxStaticBox (this, -1, wxT("New legend settings") );
    wxStaticBoxSizer * vbNew = new wxStaticBoxSizer (newTxtBox, wxVERTICAL);


    wxStaticBox * activeInfoTxtBox = new wxStaticBox (this, -1, wxT("Aktywna") );
    wxStaticBoxSizer * vbActiveInfo = new wxStaticBoxSizer (activeInfoTxtBox, wxVERTICAL);

    wxStaticBox * newInfoTxtBox = new wxStaticBox (this, -1, wxT("Nowa") );
    wxStaticBoxSizer * vbNewInfo = new wxStaticBoxSizer (newInfoTxtBox, wxVERTICAL);


    hbTop->Add(vbActive);
    hbTop->Add(vbNew);

    vbTopTop->Add(hbTop);

    this->SetSizer(vbTopTop);

    //stworzenie aktywnego plutna
    //wxButton *btnActive = new wxButton(this, wxID_EXIT, wxT("Active"), wxDefaultPosition, wxSize(95, 420));
    m_canvas_active = new wxLegendCanvas( this, activeLegend, wxSize(95, 420) );

    //stworzenie nowego plutna
    newLegend = activeLegend;
    //wxButton *btnNew = new wxButton(this, wxID_EXIT, wxT("New"), wxDefaultPosition, wxSize(95, 420));
    m_canvas_new = new wxLegendCanvas( this, newLegend, wxSize(95, 420) );


    //dodanie do aktywnego
    //vbActive->Add(btnActive, 0, wxBORDER, 5);
    vbActive->Add(m_canvas_active, 0, wxALIGN_LEFT, 5);

    //vbActive->Add(btnActive, 0, wxEXPAND, 5);

    //vbActive->AddStretchSpacer();
    //vbActive->AddSpacer(5);

    //dodanie do nwego
    //vbNew->Add(btnNew, 0, wxBORDER, 5);
    vbNew->Add(m_canvas_new, 0, wxALIGN_LEFT, 5);
    //sizerLeftRight->Add(m_canvas_new, 0, wxALIGN_LEFT, 5);

    hbTop->Add(vbInfo);

    vbInfo->Add(vbActiveInfo, 0, wxBORDER, 5);
    vbInfo->Add(vbNewInfo, 0, wxBORDER, 5);

    wxString choices2[] = {wxT("Color"), wxT("Gradient")};

    wxStaticText *stActiveMin = new wxStaticText(this, -1, wxT("Min"));
    wxStaticText *stActiveMax = new wxStaticText(this, -1, wxT("Max"));
    wxStaticText *stActiveStart = new wxStaticText(this, -1, wxT("Start"));
    wxStaticText *stActiveEnd = new wxStaticText(this, -1, wxT("End"));
    wxStaticText *stActiveRengNr = new wxStaticText(this, -1, wxT("L."));
    wxStaticText *stActiveLegType = new wxStaticText(this, -1, wxT("Typ"));
    wxStaticText *stActiveOutBlck = new wxStaticText(this, -1, wxT("Kolor czarny na zewn�trz"));


    wxStaticText *stNewMin = new wxStaticText(this, -1, wxT("Min"));
    wxStaticText *stNewMax = new wxStaticText(this, -1, wxT("Max"));
    wxStaticText *stNewStart = new wxStaticText(this, -1, wxT("Start"));
    wxStaticText *stNewEnd = new wxStaticText(this, -1, wxT("End"));
    wxStaticText *stNewRengNr = new wxStaticText(this, -1, wxT("L."));
    wxStaticText *stNewLegType = new wxStaticText(this, -1, wxT("Typ"));
    wxStaticText *stNewOutBlck = new wxStaticText(this, -1, wxT("Kolor czarny na zewn�trz"));

        //min
    tcActiveMin = new wxTextCtrl(this, -1);
        //max
    tcActiveMax = new wxTextCtrl(this, -1);
        //start
    tcActiveStart = new wxTextCtrl(this, -1);
        //end
    tcActiveEnd = new wxTextCtrl(this, -1);
        //reng nr
    //tcActiveRengNr = new wxTextCtrl(this, -1);

    scActiveRengNr = new wxSpinCtrl( this, ID_SPINCTRL, wxEmptyString );
    //scActiveRengNr = new wxSpinCtrl( panel, ID_SPINCTRL, wxEmptyString, wxPoint(200, 160), wxSize(80, wxDefaultCoord) );
    //scActiveRengNr->SetRange(-10,30);
    //scActiveRengNr->SetValue(15);

        //leg type
    //wxTextCtrl *tcActiveLegType = new wxTextCtrl(this, -1);
    rbActiveLegType = new wxRadioBox( this, ID_RADIOBOX, _T(""), wxDefaultPosition, wxDefaultSize, WXSIZEOF(choices2), choices2, 1, wxRA_SPECIFY_COLS );
    //wxRA_SPECIFY_ROWS
        //out black
    //wxTextCtrl *tcActiveOutBlck = new wxTextCtrl(this, -1);
    cbActiveOutBlck = new wxCheckBox( this, MenuColorCheckBox, wxT(""), wxDefaultPosition, wxDefaultSize );
    cbActiveOutBlck->SetValue( legend.GetOutsideBlack() );

    tcActiveMin->Disable();
    tcActiveMax->Disable();
    tcActiveStart->Disable();
    tcActiveEnd->Disable();
    scActiveRengNr->Disable();
    rbActiveLegType->Disable();
    cbActiveOutBlck->Disable();

    //tcNewMin->Disable();
    //tcNewMax->Disable();

        //min
    //tcNewMin = new wxTextCtrl(this, -1);
    tcNewMin = new wxTextCtrl(this, -1, wxT(""), wxDefaultPosition, wxDefaultSize, 0,
        wxTextValidator(
        //wxFILTER_ALPHA
        wxFILTER_NUMERIC
        )
    );


        //mx
    tcNewMax = new wxTextCtrl(this, -1);
        //start
    tcNewStart = new wxTextCtrl(this, -1);
        //end
    tcNewEnd = new wxTextCtrl(this, -1);
        //reng nr
    //tcNewRengNr = new wxTextCtrl(this, -1);
    scNewRengNr = new wxSpinCtrl( this, ID_SPINCTRL, wxEmptyString );

    scNewRengNr->SetRange(1, 50);
        //leg type

    //wxTextCtrl *tcNewLegType = new wxTextCtrl(this, -1);
    rbNewLegType = new wxRadioBox( this, ID_RADIOBOX, _T(""), wxDefaultPosition, wxDefaultSize, WXSIZEOF(choices2), choices2, 1, wxRA_SPECIFY_COLS );

        //out black

    //wxTextCtrl *tcNewOutBlck = new wxTextCtrl(this, -1);
    cbNewOutBlck = new wxCheckBox( this, MenuColorCheckBox, wxT(""), wxDefaultPosition, wxDefaultSize );

    //stworzenie wnetrza
        //aktywnego
    //wxButton *btnActiveInfo = new wxButton(this, wxID_EXIT, wxT("ActiveInfo"), wxDefaultPosition, wxSize(140, 199));
    //wxFlexGridSizer *fgActiveInfo = new wxFlexGridSizer(6, 2, 9, 25);
    wxFlexGridSizer *fgActiveInfo = new wxFlexGridSizer(6, 2, 1, 1);
    //wxFlexGridSizer *fgActiveInfo = new wxFlexGridSizer(6, 2, 5, 20);

    fgActiveInfo->Add(stActiveMin, 1, wxALIGN_CENTER);
    fgActiveInfo->Add(tcActiveMin, 1, wxEXPAND);

    fgActiveInfo->Add(stActiveMax, 1, wxALIGN_CENTER);
    fgActiveInfo->Add(tcActiveMax, 1, wxEXPAND);

    fgActiveInfo->Add(stActiveStart, 1, wxALIGN_CENTER);
    fgActiveInfo->Add(tcActiveStart, 1, wxEXPAND);

    fgActiveInfo->Add(stActiveEnd, 1, wxALIGN_CENTER);
    fgActiveInfo->Add(tcActiveEnd, 1, wxEXPAND);

    fgActiveInfo->Add(stActiveRengNr, 1, wxALIGN_CENTER);
    //fgActiveInfo->Add(tcActiveRengNr, 1, wxEXPAND);
    fgActiveInfo->Add(scActiveRengNr, 1, wxEXPAND);


    fgActiveInfo->Add(stActiveLegType, 1, wxALIGN_CENTER);
    //fgActiveInfo->Add(tcActiveLegType, 1, wxEXPAND);
    fgActiveInfo->Add(rbActiveLegType, 1, wxEXPAND);


    ////fgActiveInfo->Add(tcActiveOutBlck, 1, wxEXPAND);
    //fgActiveInfo->Add(stActiveOutBlck, 1, wxALIGN_CENTER_VERTICAL );
    //fgActiveInfo->Add(cbActiveOutBlck, 1, wxEXPAND);

    fgActiveInfo->Add(cbActiveOutBlck, 1, wxALIGN_CENTER_VERTICAL);
    fgActiveInfo->Add(stActiveOutBlck, 1, wxEXPAND );

        //nowego
    //wxButton *btnNewInfo = new wxButton(this, wxID_EXIT, wxT("NewInfo"), wxDefaultPosition, wxSize(140, 199));

    //wxFlexGridSizer *fgNewInfo = new wxFlexGridSizer(6, 2, 5, 20);
    wxFlexGridSizer *fgNewInfo = new wxFlexGridSizer(6, 2, 1, 1);

    fgNewInfo->Add(stNewMin, 1, wxALIGN_CENTER);
    fgNewInfo->Add(tcNewMin, 1, wxEXPAND);

    fgNewInfo->Add(stNewMax, 1, wxALIGN_CENTER);
    fgNewInfo->Add(tcNewMax, 1, wxEXPAND);


    fgNewInfo->Add(stNewStart, 1, wxALIGN_CENTER);
    fgNewInfo->Add(tcNewStart, 1, wxEXPAND);

    fgNewInfo->Add(stNewEnd, 1, wxALIGN_CENTER);
    fgNewInfo->Add(tcNewEnd, 1, wxEXPAND);

    fgNewInfo->Add(stNewRengNr, 1, wxALIGN_CENTER);
    //fgNewInfo->Add(tcNewRengNr, 1, wxEXPAND);
    fgNewInfo->Add(scNewRengNr, 1, wxEXPAND);


    fgNewInfo->Add(stNewLegType, 1, wxALIGN_CENTER);
    fgNewInfo->Add(rbNewLegType, 1, wxEXPAND);

    fgNewInfo->Add(cbNewOutBlck, 1, wxALIGN_CENTER);
    fgNewInfo->Add(stNewOutBlck, 1, wxEXPAND);




    //wypelnienie wnetrza
        //altywnego
    //vbActiveInfo->Add(btnActiveInfo, 0, wxBORDER, 5);
    vbActiveInfo->Add(fgActiveInfo, 0, wxBORDER, 5);
        //nowego
    //vbNewInfo->Add(btnNewInfo, 0, wxBORDER, 5);
    vbNewInfo->Add(fgNewInfo, 0, wxBORDER, 5);


    //dodanie horyzontalnego
    //w ktorym sa dwa przyciski

    btnRefresh = new wxButton(this, PanelNewRefresh, wxT("Refresh"));
    btnActive =  new wxButton(this, PanelNewActive, wxT("Active"));

    wxBoxSizer *hb2Buttons = new wxBoxSizer(wxHORIZONTAL);

    hb2Buttons->Add(btnRefresh);
    hb2Buttons->Add(btnActive);

    vbInfo->Add(hb2Buttons, 0, wxBORDER, 10);

    wxSizer * sizerOkHelp = CreateButtonSizer(wxOK | wxHELP);
    vbTopTop->Add(sizerOkHelp, 0, wxALIGN_CENTER_HORIZONTAL|wxALL, 5);



    //ustawienie wartosci w aktive
    SetActiveValues(activeLegend);
    SetNewValues(newLegend);

    //ustawienie wartosci w new


}

wxLegendDialog::~wxLegendDialog()
{

}

void wxLegendDialog::MsgBoxErr(const wxString & m, const wxString & t)
{
    wxMessageBox(m, t, wxOK | wxICON_ERROR);
}

void wxLegendDialog::OnPanelsResetToDefaults( wxCommandEvent & event )
{
	newLegend.Set(newLegend.GetMin(),newLegend.GetMax(),newLegend.GetMin(),newLegend.GetMax());
	SetNewValues(newLegend);
	m_canvas_new->Refresh(true);
	wxCommandEvent evt;
	OnPanelNewActive(evt);
}

void wxLegendDialog::OnPanelNewRefresh( wxCommandEvent &event )
{

    #ifdef __WXDEBUG__
        cout << "wxLegendDialog::OnPanelNewRefresh" << endl;
    #endif

    //sprawdzenie poprawnosci pod kontem struktury danych (czyli czy sa liczbami double)

    double d;
	FemViewer::LegendType lt;
    double min, max, start, end;
    unsigned int rn;

    //wxString s;
    //s.ToDouble()




    //test min
    if ( !tcNewMin->GetValue().ToDouble(&d) )
    {
        //blad
        #ifdef __WXDEBUG__
            cout << "blad podania wartosci min" << endl;
        #endif
        MsgBoxErr( wxT("blad podania wartosci min") );
        return;
    }

    min = d;
    #ifdef __WXDEBUG__
        cout << "min: " << d << endl;
    #endif

    //Test max
    //1
    //wxString s;
    //s = tcNewMax->GetValue();
    //if ( !s.ToDouble(&d) )
    //2
    if ( !tcNewMax->GetValue().ToDouble(&d) )
    {
        //blad
        #ifdef __WXDEBUG__
            cout << "blad podania wartosci max" << endl;
        #endif
        MsgBoxErr( wxT("blad podania wartosci max") );
        return;
    }

    max = d;
    #ifdef __WXDEBUG__
        cout << "max: " << d << endl;
    #endif


    if ( !tcNewStart->GetValue().ToDouble(&d) )
    {
        //blad
        #ifdef __WXDEBUG__
            cout << "blad podania wartosci start" << endl;
        #endif
        MsgBoxErr( wxT("blad podania wartosci start") );
        return;
    }

    start = d;
    #ifdef __WXDEBUG__
        cout << "start: " << d << endl;
    #endif


    if ( !tcNewEnd->GetValue().ToDouble(&d) )
    {
        //blad
        #ifdef __WXDEBUG__
            cout << "blad podania wartosci end" << endl;
        #endif
        MsgBoxErr( wxT("blad podania wartosci end") );
        return;
    }

    end = d;
    #ifdef __WXDEBUG__
        cout << "end: " << d << endl;
    #endif



    #ifdef __WXDEBUG__
        cout << "liczba przedzialow: " << scNewRengNr->GetValue() << endl;
    #endif
    rn = scNewRengNr->GetValue();


    #ifdef __WXDEBUG__
        cout << "typ legendy: " << rbNewLegType->GetSelection()<< endl;
    #endif

    switch( rbNewLegType->GetSelection() )
    {
        case 0:
			lt = FemViewer::SOLID;
            break;
        case 1:
			lt = FemViewer::GRADIENT;
            break;
    }



    bool b;

    b = cbNewOutBlck->IsChecked();

    #ifdef __WXDEBUG__
        cout << "out side black: ";
        if (b)
        {
            cout << "true" << endl;
        }
        else
        {
            cout << "false" << endl;
        }
    #endif


    newLegend.Set(
            min,//-1.0,
            max,//2.0,
            start,//0.0,//od - od do podzial na 4 przedzialy
            end,//5.0,//do
            rn,//5,//20,
            lt,//COLOR,
            b//false
            //true
            );

    m_canvas_new->Refresh(true);

    //ustawienie nowych wartosci legendy

    //odrysownie nowej legendy

    //podanie wartosci z nowej legendy

}

void wxLegendDialog::OnPanelNewActive( wxCommandEvent &event )
{

    activeLegend = newLegend;

    SetActiveValues(activeLegend);

    m_canvas_active->Refresh(true);

}

void wxLegendDialog::SetActiveValues(const FemViewer::Legend & l)
{
    wxString s;

    s.Clear();
    s << l.GetMin();

    tcActiveMin->SetValue(s);

    s.Clear();
    s << l.GetMax();
    tcActiveMax->SetValue(s);

    s.Clear();
    s << l.GetStart();
    tcActiveStart->SetValue(s);

    s.Clear();
    s << l.GetEnd();
    tcActiveEnd->SetValue(s);

    scActiveRengNr->SetValue(l.GetRangesNumber());

    switch(l.GetLegendType())
    {
	case FemViewer::SOLID:
            rbActiveLegType->SetSelection(0);
            break;
	case FemViewer::GRADIENT:
            rbActiveLegType->SetSelection(1);
            break;
    }

    cbActiveOutBlck->SetValue(l.GetOutsideBlack());

}

void wxLegendDialog::SetNewValues(const FemViewer::Legend & l)
{

    wxString s;

    s.Clear();
    s << l.GetMin();
    tcNewMin->SetValue(s);

    s.Clear();
    s << l.GetMax();
    tcNewMax->SetValue(s);

    s.Clear();
    s << l.GetStart();
    tcNewStart->SetValue(s);

    s.Clear();
    s << l.GetEnd();
    tcNewEnd->SetValue(s);
    std::cout<<"start = " << l.GetStart() << " end = " << l.GetEnd()  << std::endl;

    scNewRengNr->SetValue(l.GetRangesNumber());

    switch(l.GetLegendType())
    {
	case FemViewer::SOLID:
            rbNewLegType->SetSelection(0);
            break;
	case FemViewer::GRADIENT:
            rbNewLegType->SetSelection(1);
            break;
    }

    cbNewOutBlck->SetValue(l.GetOutsideBlack());

}
