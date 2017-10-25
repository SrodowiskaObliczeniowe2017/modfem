#ifndef _SOL_CALC_DIALOG_H_
#define _SOL_CALC_DIALOG_H_

#include <wx/wx.h>
#include <wx/defs.h>
#include <wx/spinctrl.h>
#include <string>
#include <vector>
#include "CutPlane.h"


namespace FemViewer {
class BBox3D;

namespace WxGUI		{

//! Okno dialogowe dajeca mozliwosc ustalenia funkcji
class wxSolCalcDialog : public wxDialog
{
	public:
		// Constructors
		wxSolCalcDialog(
			int & nreq,	     // the number of components of each solution
			int & nr_sol,		 // the number of solutions
			int & curr_sol,// the indef current solution - default: first solution
			const std::string& formula // the formula expresion
		);

        // Gets value from range [nr_sol_start, nr_sol_start + nr_sol - 1]
        // 0 - nreq-1
		int GetCurrSol(void){return _spinbutton->GetValue(); };

		// Sets values
		void SetCurrSol(int val);

		// Gets function's formula
		std::string GetFunction(void);
		
		// Updates function's description field
		void SetFunction(const std::string &str);

	private:

		// Tne number of components in solution
		int _nreq;
		// The number of solutions
		int _nr_sol;
		// The current solution's index
		int _curr_sol;

        wxTextCtrl   * _spintext;
        wxSpinButton * _spinbutton;
        wxTextCtrl   * _textctrl;

        wxString	   _txtStd;  //(wxT("v0"));
        std::string   _formula;

        void OnSpinButtonUpdate( wxSpinEvent &event );
        void OnButtonTextReset( wxCommandEvent& event );
        void OnButtonTextFuncTest( wxCommandEvent& event );

        void OnOk( wxCommandEvent& event );
        void OnHelp( wxCommandEvent& event );

        enum
        {
            ID_SPIN = 1000,
            ID_TEXT_RESET,
            ID_TEXT_FUNC_TEST
        };

        DECLARE_EVENT_TABLE()

		inline void ClumpToRange(int val)
		{
			if(val < 0) _curr_sol = 0;
			else if(val > (_nr_sol - 1)) _curr_sol = _nr_sol -1;
			else _curr_sol = val;
		}

};

//! Okno dialogowe pozwalaj�ce na mo�liwo�� modyfikacjiparamert�w p�aszczyzny
class wxPlaneCutDialog : public wxDialog
{
	DECLARE_CLASS(wxPlaneCutDialog)
	DECLARE_EVENT_TABLE()

	public:

		wxPlaneCutDialog();
	    wxPlaneCutDialog(wxWindow * parent, const BBox3D& bbox, std::vector<CutPlane>& ctpls,
	    		const bool use, const bool display);
        ~wxPlaneCutDialog();

        void SetUseAll(bool use) { m_useAll = use; }
        bool GetUseAll() const { return m_useAll; }

        void SetDisplayAll(bool display) { m_displayAll = display; }
        bool GetDisplayAll() const { return m_displayAll; }

        void SetCutPlanes(const std::vector<CutPlane> planes) { m_planes = planes; }
        std::vector<CutPlane>& GetCutPlanes() { return m_planes; }


    protected:

        enum {
        	ID_USE_ALL = 100,
        	ID_DISP_ALL,
        	ID_NUM_PLANES,
        	ID_USE_ITEM = 200,
        	ID_DISP_ITEM = 300,
        };

        static const int max_planes = 10;
        const BBox3D& m_bbox;
        bool 	m_useAll;
        bool 	m_displayAll;
        std::vector<CutPlane> m_planes;
        int  	m_numPls;

    private:
        // Private controls
        wxSpinCtrl *m_scNumPlsl;
        wxCheckBox *m_cbUseAll, *m_cbDisplayAll;
        wxCheckBox *m_cbUseItem[max_planes];
        wxCheckBox *m_cbDisplayItem[max_planes];
        wxTextCtrl *m_tcParams[max_planes][4];
        wxFlexGridSizer *m_fgSizer;
        wxBoxSizer *m_topSizer;

    private:

    	//void SetProperties();
        void Init();
        void DoLayout();
        void OnUseAll(wxCommandEvent& evt);
        void OnDisplayAll(wxCommandEvent& evt);
        void OnNumPlanes(wxSpinEvent& evt);
        void OnOk( wxCommandEvent& event );

        void HideRow(size_t index);
        void ShowRow(size_t index);
        void OnSize(wxSizeEvent& evt);


};

}// end namespace WxGUI
}// end namespace FemViewer
#endif /* _SOL_CALC_DIALOG_H_
*/


