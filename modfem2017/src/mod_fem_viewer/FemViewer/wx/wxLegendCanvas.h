#ifndef WXMGRVIEWDIALOGLEGENDCANVAS_H_INCLUDED
#define WXMGRVIEWDIALOGLEGENDCANVAS_H_INCLUDED

#include <wx/wx.h>
#include <wx/dc.h>
//#include "wxMGRViewDialogLegend.h"
//#include "../LegendValue.h"
#include "../inc/Legend.h"

/*! \class wxMgrViewDialogLegendCanvas
    \brief Legend class.

    Klasa dziedziczonapo wxScrolledWindow sluzaca do rysowania pozwalajaca na przewijanie.
    Ma zaimplemenetowana funkcje rysowania na podstawie przekazanej legendy.
    Odczytuje jej parametry i na ich podstawie rysuje.
*/
class wxLegendCanvas: public wxScrolledWindow
{
	public:
		wxLegendCanvas(wxDialog* parent, FemViewer::Legend& legend, const wxSize& size = wxDefaultSize);
		~wxLegendCanvas();

		FemViewer::Legend& legend;

	protected:

	private:

        wxRadioBox  *m_radioColGrad;//wybor miedzy kolorem a gradientem

        void OnPaint(wxPaintEvent &event);
        void OnPaint2(wxPaintEvent &event);
        void OnPaint3(wxPaintEvent &event);

		void GradientFillLinear(wxDC &dc, const wxRect& rect, FemViewer::colormap_value& initial, FemViewer::colormap_value& dest, wxDirection nDirection);

        DECLARE_EVENT_TABLE()

        //std::string doubleToString(double inValue);

};


#endif // WXMGRVIEWDIALOGLEGENDCANVAS_H_INCLUDED
