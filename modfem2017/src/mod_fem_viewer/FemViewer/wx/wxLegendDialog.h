#ifndef _WX_LEGEND_DIALOG_H_
#define _WX_LEGEND_DIALOG_H_

#include "../inc/Legend.h"
#include <wx/spinctrl.h>
#include "wxLegendCanvas.h"

class wxLegendDialog: public wxDialog
{

    public:

		wxLegendDialog(const wxString& title, long style, FemViewer::Legend& legend, const wxSize& size = wxSize(300,400));

		wxLegendDialog(const wxString& title, FemViewer::Legend& legend);

		~wxLegendDialog();

    protected:

		FemViewer::Legend& activeLegend;
		FemViewer::Legend  newLegend;

        wxLegendCanvas *m_canvas_active;
        wxLegendCanvas *m_canvas_new;

        //active
        wxTextCtrl *tcActiveMin;
        wxTextCtrl *tcActiveMax;
        wxTextCtrl *tcActiveStart;
        wxTextCtrl *tcActiveEnd;
        //wxTextCtrl *tcActiveRengNr;
        wxSpinCtrl *scActiveRengNr;
        wxRadioBox *rbActiveLegType;
        wxCheckBox *cbActiveOutBlck;


        //new
        wxTextCtrl *tcNewMin;
        wxTextCtrl *tcNewMax;
        wxTextCtrl *tcNewStart;
        wxTextCtrl *tcNewEnd;
        //wxTextCtrl *tcNewRengNr;
        wxSpinCtrl *scNewRengNr;
        wxRadioBox *rbNewLegType;
        wxCheckBox *cbNewOutBlck;

        // Buttons
        wxButton *btnReset, *btnRefresh, *btnActive;

        enum
        {
            ID_RADIOBOX,
            MenuColorCheckBox,
            ID_SPINCTRL,
            ResetPanels,
            PanelNewRefresh,
            PanelNewActive
        };

        //!Ustawia wartosci w czesci informacyjnej active
		void SetActiveValues(const FemViewer::Legend & l);
        //!Ustawia wartosci w czesci informacyjnej new
		void SetNewValues(const FemViewer::Legend & l);

		//! This resets legend to defaults
		void OnPanelsResetToDefaults( wxCommandEvent &event);
        //! Obs�uga przycisku odswierzajacego nowa liste legendy na podstawie podanych paramterow
        void OnPanelNewRefresh( wxCommandEvent &event );
        //! Obs�uga przycisku aktywujaca nowa liste legendy jako aktualna
        void OnPanelNewActive( wxCommandEvent &event );

        void MsgBoxErr(const wxString & m, const wxString & t = wxT("Error") );

        DECLARE_EVENT_TABLE()

};

#endif /* _WX_LEGEND_DIALOG_H_
*/
