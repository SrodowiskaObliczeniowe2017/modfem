#include "wxLegendCanvas.h"
#include "wxLegendDialog.h"
#include "../inc/Color.h"

#ifdef __WXDEBUG__
    #include <iostream>
    using namespace std;
#endif

BEGIN_EVENT_TABLE(wxLegendCanvas, wxScrolledWindow)
    //EVT_PAINT  (wxLegendCanvas::OnPaint)
    //EVT_PAINT  (wxLegendCanvas::OnPaint2)
    EVT_PAINT  (wxLegendCanvas::OnPaint3)
END_EVENT_TABLE()

wxLegendCanvas::wxLegendCanvas(wxDialog*parent,  FemViewer::Legend& legend, const wxSize& size)
//: wxScrolledWindow(parent, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxHSCROLL | wxVSCROLL | wxNO_FULL_REPAINT_ON_RESIZE)
: wxScrolledWindow(parent, wxID_ANY, wxDefaultPosition, wxSize(290, 420), wxSUNKEN_BORDER|wxHSCROLL|wxVSCROLL), legend(legend)
//: wxScrolledWindow(parent, wxID_ANY, wxDefaultPosition, size, wxSUNKEN_BORDER|wxHSCROLL|wxVSCROLL), legend(legend)

//: wxScrolledWindow(parent, wxID_ANY, wxDefaultPosition, wxDefaultSize, wxNO_FULL_REPAINT_ON_RESIZE )
{ }

wxLegendCanvas::~wxLegendCanvas()
{ }

void wxLegendCanvas::OnPaint3(wxPaintEvent &event)
{
    #ifdef __WXDEBUG__
        cout << "wxLegendCanvas::OnPaint3" << endl;
    #endif

    wxString s;

    wxPaintDC pdc(this);

    wxDC &dc = pdc ;

    dc.Clear();

    if (!legend.IsValid())
        return;

    wxSize size;
    size =  this->GetSize();

    int height = size.GetHeight() - 15; // rozmiar okna

    //liczba elementow, liczba granic
    int lGr = legend.GetLiczbaGranic();
    int lEl = lGr - 1;

    int deltaHeightMod = height % lEl;
    height = height - deltaHeightMod;
    int deltaHeight = height / lEl;

    const int UpText = 9;

    wxRect rLeg(10, 5, 20, deltaHeight-1);
   	FemViewer::colormap_value a, b;



    //const int sbx = 25;//startBorderX
    const int sbx = 5;//startBorderX
    const int ebx = 50;//endBorderX
    const int stx = 60;//startTextX

	FemViewer::ColorRGB rgbA, rgbB;
    double dA, dB;



    for (int i = 0; i < lEl; ++i)
    {
        if ( i > 0 )
            rLeg.Offset(0, rLeg.height+1);


        legend.GetRange3(i, dA, rgbA, dB, rgbB);

        b.value = dA;
        b.SetColor(rgbA);

        a.value = dB;
        a.SetColor(rgbB);

        GradientFillLinear(dc, rLeg, b, a, wxDOWN);

        dc.DrawLine(sbx, rLeg.GetTop()-1, ebx, rLeg.GetTop()-1);

        s.Clear();
        s << b.value;
        //s << 1;
        dc.DrawText(s, stx, rLeg.GetTop() - UpText);


        if (i == (lEl-1))
        {
            dc.DrawLine(sbx, rLeg.GetBottom()+1, ebx, rLeg.GetBottom()+1);

            s.Clear();
            s << a.value;
            //s << 1;
            dc.DrawText(s, stx, rLeg.GetBottom() - UpText);
        }


    }



}
/*
void wxLegendCanvas::OnPaint2(wxPaintEvent &event)
{

    #ifdef __WXDEBUG__
        cout << "wxLegendCanvas::OnPaint2" << endl;
    #endif


    //wxString v;
    //v.Clear();

    wxString s;

    wxPaintDC pdc(this);

    wxDC &dc = pdc ;

    dc.Clear();

    if (legend.GetColors().size() < 1)
        return;

    wxSize size;
    size =  this->GetSize();

    #ifdef __WXDEBUG__
        cout << "wysokosc:" << size.GetHeight() << endl;
        cout << "szerokosc:" << size.GetWidth() << endl;
    #endif

    unsigned int lEl, lGr;

    //z = legend.GetRangesCount();//liczba elementow

    lEl = legend.GetRangesCount();
    lGr = legend.GetBordersCount();

    #ifdef __WXDEBUG__
        cout << "liczba elementow: " << lEl << endl;
        cout << "liczba granic: " << lGr << endl;
    #endif

    int height = size.GetHeight(); // rozmiar okna
    //size.GetWidth()
    //size.GetHeight()

    height = height - 10;

    #ifdef __WXDEBUG__
        cout << "wysokosc(420-10):" << height << endl;
    #endif

    int deltaHeight;//przedzial wysokosci dla elementu
    int deltaHeightMod;

    deltaHeightMod = height % lEl;
    //deltaHeightMod = height % lGr;

    #ifdef __WXDEBUG__
        cout << height << " % "<< lEl << " = " << deltaHeightMod << endl;
    #endif

    #ifdef __WXDEBUG__
        cout << height << " - "<< deltaHeightMod << " = " << height - deltaHeightMod << endl;
    #endif

    height = height - deltaHeightMod;

    #ifdef __WXDEBUG__
        cout << height << " / "<< lEl << " = " << height / lEl << endl;
    #endif

    deltaHeight = height / lEl;
    //deltaHeight = height / lGr;

    const int UpText = 9;

    //wxRect rLeg(sX, sY, sX+x, sY+y);
    wxRect rLeg(0, 0, 20, deltaHeight-1);

    FemViewer::LegendValue  a, b;

    rLeg.Offset(10, 5);

    //const int sbx = 25;//startBorderX
    const int sbx = 5;//startBorderX
    const int ebx = 50;//endBorderX
    const int stx = 60;//startTextX

	FemViewer::ColorRGB rgbA, rgbB;
    double dA, dB;


    /*
    for (unsigned int i = lEl; i > 0; --i)
    {

        if ( i < lEl)
            rLeg.Offset(0, rLeg.height+1);

        legend.GetRange(i, dA, rgbA, dB, rgbB);


        a.value = dA;
        a.SetColor(rgbA);

        b.value = dB;
        b.SetColor(rgbB);

        GradientFillLinear(dc, rLeg, b, a, wxDOWN);

        dc.DrawLine(sbx, rLeg.GetTop()-1, ebx, rLeg.GetTop()-1);

        v.Clear();
        v << b.value;

        dc.DrawText(v, stx, rLeg.GetTop() - UpText);

        if (i == 1)
        {
            v.Clear();
            v << a.value;

            dc.DrawLine(sbx, rLeg.GetBottom()+1, ebx, rLeg.GetBottom()+1);
            //dc.DrawText(wxT("aaa"), stx, rLeg.GetBottom() - UpText);
            dc.DrawText(v, stx, rLeg.GetBottom() - UpText);
        }



    }//for (unsigned int i = lEl; i > 0; --i)
    * /


    for (unsigned int i = 0; i < lEl; ++i)
    {
        if ( i > 0 )
            rLeg.Offset(0, rLeg.height+1);


        legend.GetRange(i, dA, rgbA, dB, rgbB);
        //legend.GetRange(lEl-i, dA, rgbA, dB, rgbB);

        b.value = dA;
        b.SetColor(rgbA);

        a.value = dB;
        a.SetColor(rgbB);

        GradientFillLinear(dc, rLeg, b, a, wxDOWN);

        dc.DrawLine(sbx, rLeg.GetTop()-1, ebx, rLeg.GetTop()-1);

        s.Clear();
        s << b.value;
        dc.DrawText(s, stx, rLeg.GetTop() - UpText);


        if (i == (lEl-1))
        {
            dc.DrawLine(sbx, rLeg.GetBottom()+1, ebx, rLeg.GetBottom()+1);

            s.Clear();
            s << a.value;
            dc.DrawText(s, stx, rLeg.GetBottom() - UpText);
        }


    }



}
*/
/*
void wxLegendCanvas::OnPaint(wxPaintEvent &event)
{

    #ifdef __WXDEBUG__
        cout << "wxLegendCanvas::OnPaint" << endl;
    #endif

    wxString v;

    v.Clear();

    wxPaintDC pdc(this);

    wxDC &dc = pdc ;

    dc.Clear();

    if (legend.legendValues.size() < 1)
        return;

    wxSize size;
    size =  this->GetSize();

    #ifdef __WXDEBUG__
        cout << "wysokosc:" << size.GetHeight() << endl;
        cout << "szerokosc:" << size.GetWidth() << endl;
    #endif


    //int z, lEl, lGr;//liczba elementow, liczba granic
    unsigned int z, lEl, lGr;

    ////int lEl = legend.legendValues.size();//liczba elementow
    z = legend.legendValues.size();//liczba elementow
    //z = legend.GetRangesCount();
    //z = 2;
    //z = 4;
    //z = 5;
    //z = 7;
    //z = 8;
    //z = 10;
    //z = 12;
    //z = 15;


    lEl = z - 1;
    lGr = z;

    #ifdef __WXDEBUG__
        cout << "liczba elementow: " << lEl << endl;
        cout << "liczba granic: " << lGr << endl;
    #endif

    int height = size.GetHeight(); // rozmiar okna
    //size.GetWidth()
    //size.GetHeight()

    height = height - 10;

    #ifdef __WXDEBUG__
        cout << "wysokosc(420-10):" << height << endl;
    #endif

    //const int sLeft = 10;
    //const int sTop = 5;
    //const int legendHeight = 20;

    int deltaHeight;//przedzial wysokosci dla elementu
    int deltaHeightMod;

    deltaHeightMod = height % lEl;

    #ifdef __WXDEBUG__
        cout << height << " % "<< lEl << " = " << deltaHeightMod << endl;
    #endif

    #ifdef __WXDEBUG__
        cout << height << " - "<< deltaHeightMod << " = " << height - deltaHeightMod << endl;
    #endif

    height = height - deltaHeightMod;

    #ifdef __WXDEBUG__
        cout << height << " / "<< lEl << " = " << height / lEl << endl;
    #endif

    deltaHeight = height / lEl;

    //int sX = sLeft;
    //int sY = sTop;

    //int x = legendHeight;
    //int y = deltaHeight;

    const int UpText = 9;

    //wxRect rLeg(sX, sY, sX+x, sY+y);
    wxRect rLeg(0, 0, 20, deltaHeight-1);

    FemViewer::LegendValue  a, b;

    //a.red = 1.0;
    //a.green = 1.0;
    //a.blue = 1.0;

    //b.red = 1.0;
    //b.green = 1.0;
    //b.blue = 1.0;

    a.red = 0.9f;
    a.green = 0.5f;
    a.blue = 0.0f;

    b.red = 0.0f;
    b.green = 0.0f;
    b.blue = 1.0f;

    //dc.DrawRectangle(10, 5, 100, 410);
    //dc.DrawRectangle(10, 5, 100, 410-deltaHeightMod);

    rLeg.Offset(10, 5);

    //const int sbx = 25;//startBorderX
    const int sbx = 5;//startBorderX
    const int ebx = 50;//endBorderX
    const int stx = 60;//startTextX

    //for (unsigned int i = 0; i<=lEl; ++i)



    for (unsigned int i = lEl; i > 0; --i)
    {
        if ( i < lEl)
            rLeg.Offset(0, rLeg.height+1);

		if (legend.legendType == FemViewer::GRADIENT)
        {
            //jesli gradient
            b = legend.legendValues[i];
            a = legend.legendValues[i-1];
        }
		else if (legend.legendType == FemViewer::SOLID)
        {
            //b = legend.legendValues[i];
            //a = legend.legendValues[i-1];
            b = legend.legendValues[i-1];
            a = legend.legendValues[i-1];
        }
        //jesli jednolity kolor zaczynajac od pierwszego
        //b = legend.legendValues[i];
        //a = legend.legendValues[i];

        GradientFillLinear(dc, rLeg, b, a, wxDOWN);

        dc.DrawLine(sbx, rLeg.GetTop()-1, ebx, rLeg.GetTop()-1);

        v.Clear();
        v << b.value;

        //dc.DrawText(wxT("aaa"), stx, rLeg.GetTop() - UpText);
        dc.DrawText(v, stx, rLeg.GetTop() - UpText);

        //if (i == lEl-1)
        //{
        //    dc.DrawLine(sbx, rLeg.GetBottom(), ebx, rLeg.GetBottom());
        //    dc.DrawText(wxT("aaa"), stx, rLeg.GetBottom() - UpText);
        //}

        if (i == 1)
        {
            v.Clear();
            v << a.value;

            dc.DrawLine(sbx, rLeg.GetBottom()+1, ebx, rLeg.GetBottom()+1);
            //dc.DrawText(wxT("aaa"), stx, rLeg.GetBottom() - UpText);
            dc.DrawText(v, stx, rLeg.GetBottom() - UpText);
        }
    }

    /*
    for (unsigned int i = 0; i<lEl; ++i)
    {

        #ifdef __WXDEBUG__
            cout << i+1 << ", " << endl;
        #endif

        if ( i > 0)
            rLeg.Offset(0, rLeg.height-1);

        #ifdef __WXDEBUG__
            cout << "top: " << rLeg.GetTop();
            cout << " bottom: " << rLeg.GetBottom() << endl;
        #endif

        //dc.DrawRectangle(rLeg);

        GradientFillLinear(dc, rLeg, a, b, wxDOWN);

        dc.DrawLine(5, rLeg.GetTop(), 50, rLeg.GetTop());

        dc.DrawText(wxT("aaa"), 60, rLeg.GetTop() - UpText);

        if (i == lEl-1)
        {
            dc.DrawLine(5, rLeg.GetBottom(), 50, rLeg.GetBottom());
            dc.DrawText(wxT("aaa"), 60, rLeg.GetBottom() - UpText);
        }

        //tam gdzie top wpisujemy liczbe z wartoscia

        //dc.DrawRectangle(10, 5 + i*deltaHeight, legendHeight, deltaHeight);
    }
    * /
    /*
    for (unsigned int i = lEl-1; i > 0; --i)
    {
        if (i < lEl-1)
            rLeg.Offset(0, rLeg.height);
        dc.DrawRectangle(rLeg);
    }
    * /

    //testowe rysowanie prostokatu
    //dc.DrawRectangle(rLeg);

    //patla rysujaca prostokaty



}
*/
//GradientFillLinear(const wxRect& rect, const wxColour& initialColour, const wxColour& destColour, wxDirection nDirection = wxEAST)
//wxDOWN
void wxLegendCanvas::GradientFillLinear(wxDC &dc, const wxRect& rect, FemViewer::colormap_value & initial, FemViewer::colormap_value & dest, wxDirection nDirection)
{

    dc.GradientFillLinear(rect, wxColour( initial.RedAsUChar(), initial.GreenAsUChar(), initial.BlueAsUChar() ),
                                    wxColour( dest.RedAsUChar(), dest.GreenAsUChar(), dest.BlueAsUChar() ), nDirection);

}
