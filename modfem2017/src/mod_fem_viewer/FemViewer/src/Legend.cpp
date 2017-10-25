//#include "defs.h"
//#include "fv_config.h"
//#include "fv_compiler.h"
//#include "fv_inc.h"
//#include "fv_assert.h"
//#include "Common.h"
#include "defs.h"
#include "fv_txt_utls.h"
#include "fv_exception.h"
#include "Legend.h"
#include "Color.h"
#include "MathHelper.h"
#include "Log.h"

#include <cmath>
#include <limits>
#include <fstream>

using std::numeric_limits;

using namespace std;
namespace FemViewer {

	Legend::Legend()
	: _hBar()
	, _min(0.0), _max(1.0), _start(0.0), _end(1.0)
	, _nranges(_ncont)
	, _type(_dftype)
	, _outsideBlack(false)
	, _ok(false)
	, _colors()
	, color()
	{
		//mfp_log_debug("default ctr\n");
	}

	Legend::Legend(const Legend & ref)
	: _hBar(ref._hBar)
	, _min(ref._min)
	, _max(ref._max)
	, _start(ref._start)
	, _end(ref._end)
	, _nranges(ref._nranges)
	, _type(ref._type)
	, _outsideBlack(ref._outsideBlack)
	, _ok(ref._ok)
	, color(ref.color)
	{
		mfp_log_debug("Copy ctr\n");
		_colors.clear();
		_colors.resize(ref._colors.size());
		std::copy (ref._colors.begin(), ref._colors.end(), _colors.begin());
	}

	Legend::Legend(const std::string& path)
	{
		std::ifstream inFile(path.c_str());
		std::string name;
		inFile >> name;
		if (name != "linear" || name != "texture")
			throw fv_exception("Incorrect format in colormap config file!\n");
		inFile >> _min;
		inFile >> _max;
		if (_min > _max) std::swap(_min, _max);
		inFile >> _nranges; // remember true _nranges  = _nranges-1;
		assert(_nranges >= 2);
		colormap_value cv;
		if (name == "linear") {
			_colors.reserve(_nranges);
			double dlt = (_max - _min)/static_cast<double>(_nranges-1);

			for (unsigned i=0;i<_nranges;++i)
			{
				cv.value = _min + i*dlt;
				inFile >> cv.rgb.r;
				inFile >> cv.rgb.g;
				inFile >> cv.rgb.b;
				if (i==0) cv.type == BORDER | INSIDE | START | LVT_COLOR;
				else if (i==(_nranges-1)) cv.type == BORDER | INSIDE | END | LVT_COLOR;
				else cv.type = cv.type == BORDER | INSIDE | LVT_COLOR;
				_colors.push_back(cv);
			}
		}
		else { // texture
			_colors.reserve(2*_nranges);
			double dlt = (_max - _min)/static_cast<double>(_nranges-1);
			for (unsigned i=0;i<_nranges;++i)
			{
				cv.clear();
				cv.value = _min + i*dlt;
				if (i==0) cv.type == BORDER | INSIDE | START | LVT_COLOR;
				else if (i==(_nranges-1)) cv.type == BORDER | INSIDE | END | LVT_COLOR;
				else cv.type = cv.type == BORDER | INSIDE | LVT_COLOR;
				_colors.push_back(cv);
				cv.clear();
				inFile >> cv.rgb.r;
				inFile >> cv.rgb.g;
				inFile >> cv.rgb.b;
				cv.type = LVT_COLOR;
				_colors.push_back(cv);
			}
		}
	}

	Legend::~Legend()
	{
		mfp_log_debug("Dtr\n");
		_colors.clear();
	}

	Legend & Legend::operator=(const Legend & ref)
	{
		mfp_log_debug("Legend::operator=\n");
		_hBar = ref._hBar;
		_min = ref._min;
		_max = ref._max;
		_start = ref._start;
		_end = ref._end;
		_nranges = ref._nranges;
		_type = ref._type;
		_outsideBlack = ref._outsideBlack;
		_ok = ref._ok;
		color = ref.color;
		_colors.clear();
		_colors.resize(ref._colors.size());
		std::copy (ref._colors.begin(), ref._colors.end(), _colors.begin());
		return *this;
	}

	bool Legend::operator==(const Legend& rhs) const
	{
		unsigned int res = 0x0;;
		res |= (1 << !(_max == rhs._max));
		res |= (2 << !(_min == rhs._min));
		res |= (4 << !(_start == rhs._start));
		res |= (8 << !(_end == rhs._end));
		res |= (16 << !(_nranges == rhs._nranges));
		res |= (32 << !(_outsideBlack == rhs._outsideBlack));
		return (res == 0) ? true : false;
	}

	void Legend::Reset()
	{
		_hBar.reset();
		_min = _start = 0.0; //
		_max = _end = 1.0;//LARGE_F;
		_nranges = 10;
		_type = GRADIENT;
		_outsideBlack = false;
		_ok = true;
		_hBar->Create(this);
	}

	void Legend::Create()
	{
//		mfp_log_debug("Legend: create\n");
		if(!_hBar) {
	//		mfp_log_debug("Create new ColorBar");
			_hBar = std::shared_ptr<ColorBar>(new ColorBar);
		}
		//mfp_debug("Before horizontal create\n");
		_hBar->Create(this);
	}
	void Legend::Set(double minValue,double maxValue,
					 double startValue,double endValue,
					 unsigned int rangesNumber,
					 LegendType legendType,
					 bool outsideBlack)
	{

		//FV_ASSERT(maxv >= minv);
		//FV_ASSERT(endv >= startv);
		//mfp_log_debug("Set\n");
		_min = minValue;
		_max = maxValue;
		_start = startValue;
		_end = endValue;
		_nranges = rangesNumber;
		_type = legendType;
		_outsideBlack = outsideBlack;
		_colors.clear();

		if (_type == SOLID) {
			FillStartEndColor();
			AddMinMaxColor();
		}
		else {
			FillStartEndGradient();
			AddMinMaxGradient();
		}
		//mfp_debug("before create\n");
		//Create();
		_ok = true;
	}


	void Legend::FillStartEndColor(void)
	{
		const float s(1.f), v(1.f);
		const double deltaStartEnd((_end - _start) / double(_nranges));
		const double deltaColor(240.0 / double(_nranges));
		float h(deltaColor / 2.0);

		ColorRGB sRgb;
		ColorHSV sHsv;
		colormap_value val;
		val.value = _start;
		for (unsigned int i = 0; i <= _nranges; ++i)
		{
			val.clear();
			val.value = _start + i * deltaStartEnd;

			// Set type
			if (i == 0) {
				val.type = INSIDE | BORDER | START;
			} else if (i == (_nranges) ) {
				val.type = INSIDE | BORDER | END;
			} else {
				val.type = INSIDE | BORDER;
			}
			_colors.push_back(val);
			//jesli mniejszy od liczbaEl-1 to dodac kolor
			if (i <= _nranges )
			{
				//podanie tyko koloru

				//reset poprzedniej wartosci
				val.clear();

				//ustawienie nowych kolorow
				sHsv.H = h;
				sHsv.S = s;
				sHsv.V = v;
				sRgb   = ConvertHSVToRGB(sHsv);
				//std::cout << sRgb;

				val.rgb.r = sRgb.R;
				val.rgb.g = sRgb.G;
				val.rgb.b = sRgb.B;
				//val.value = -1.0; //this is done by val.Reset()
				val.type  = LVT_COLOR;
				_colors.push_back(val);

				h += deltaColor;
			}
		}
	}

	void Legend::FillStartEndGradient(void)
	{
		const double deltaStartEnd((_end - _start) / double(_nranges));
		const double deltaColor(240.0 / double(_nranges));//+1); // from red to blue over hue
		//float h = 0.0f;//deltaColor / 2.0f;

		ColorHSV sHsv = {240.0f, 1.0f, 1.0f};
		ColorRGB sRgb = ConvertHSVToRGB(sHsv);
		
		colormap_value val;
		val.value = _start;
		val.rgb.r = sRgb.R;
		val.rgb.g = sRgb.G;
		val.rgb.b = sRgb.B;
		val.type  = INSIDE | BORDER | LVT_COLOR | START;

		_colors.push_back(val);
		//std::cout << "0" << " = {" << int(sRgb.R * 255) << ", " << int(sRgb.G * 255) << ", " << int(sRgb.B * 255) << "}" << std::endl;
		//std::cout << "H = " << sHsv.H << " S = " << sHsv.S << " V = " << sHsv.V << std::endl;

		for (unsigned i = _nranges - 1; i > 0; --i)
		{
			sHsv.H -= deltaColor;
			sRgb    = ConvertHSVToRGB(sHsv);

			//std::cout << i << " = {" << int(sRgb.R * 255) << ", " << int(sRgb.G * 255) << ", " << int(sRgb.B * 255) << "}" << std::endl;
			//std::cout<< "H = " << sHsv.H << " S = " << sHsv.S << " V = " << sHsv.V << std::endl;
			//zapis wartosci i kolorow
			val.value += deltaStartEnd;
			val.rgb.r = sRgb.R;
			val.rgb.g = sRgb.G;
			val.rgb.b = sRgb.B;
			val.type = INSIDE | BORDER | LVT_COLOR;
			_colors.push_back(val);
			//h += deltaColor;
		}
		sHsv.H -= deltaColor;
		sRgb = ConvertHSVToRGB(sHsv);
		//hsv = ConvertRGBToHSV(sRgb);
		val.value = _end;
		val.rgb.r = sRgb.R;
		val.rgb.g = sRgb.G;
		val.rgb.b = sRgb.B;
		val.type  = INSIDE | BORDER | LVT_COLOR | END;
		_colors.push_back(val);
		//std::cout << "20" << " = {" << int(sRgb.R * 255) << ", " << int(sRgb.G * 255) << ", " << int(sRgb.B * 255) << "}" << std::endl;
		//std::cout << "H = " << sHsv.H << " S = " << sHsv.S << " V = " << sHsv.V << std::endl;
	}

	void Legend::MinMaxPos(int & minPos, int & maxPos)
	{
		// min
		if ((fvmath::Compare(_min, _start)) || (fvmath::Compare(_min, _end))) {
			minPos = 0;
		}
		else {
			if ((_min > _start) && (_min < _end)) {
				minPos = 0;
			}
			else if (_min < _start)	{
				minPos = -1;
			}
			else if (_min > _end) {
				minPos = 1;
			}
		}

		// max
		if ((fvmath::Compare(_max, _start)) || (fvmath::Compare(_max, _end))) {
			maxPos = 0;
		}
		else {
			if ((_max > _start) && (_max < _end) ) {
				minPos = 0;
			}
			else if (_max < _start) {
				minPos = -1;
			}
			else if (_max > _end) {
				minPos = 1;
			}
		}
	}

	void Legend::AddMinMaxColor(void)
	{
		int minPos, maxPos;
		colormap_value val;

		MinMaxPos(minPos, maxPos);

		if ((minPos  < 0) && (maxPos < 0)) {
			// Both are below _start

			// Set color
			if (! _outsideBlack) {
				// First color
				val.rgb.r = _colors[1].rgb.r;
				val.rgb.g = _colors[1].rgb.g;
				val.rgb.b = _colors[1].rgb.b;
			}
			val.type = LVT_COLOR;
			_colors.insert(_colors.begin(), val);

			// max
			val.clear();
			val.value = _max;
			val.type  = OUTSIDE;
			_colors.insert(_colors.begin(), val);

			// Color
			val.clear();
			if (! _outsideBlack) {
				// First color
				val.rgb.r = _colors[1].rgb.r;
				val.rgb.g = _colors[1].rgb.g;
				val.rgb.b = _colors[1].rgb.b;
			}
			val.type = LVT_COLOR;
			_colors.insert(_colors.begin(), val);

			// min
			val.clear();
			val.value = _min;
			val.type = OUTSIDE;
			_colors.insert(_colors.begin(), val);
		}
		else if ((minPos  > 0) && (maxPos > 0)) {
			// Both are over _end:
			// add color then min, next color then max.

			// Add min value
			val.value = _min;
			val.type = OUTSIDE;
			_colors.push_back(val);

			// Add color for min
			val.clear();
			if (!_outsideBlack) {
				// Fisrt color
				val.rgb.r = _colors[_colors.size()-2].rgb.r;
				val.rgb.g = _colors[_colors.size()-2].rgb.g;
				val.rgb.b = _colors[_colors.size()-2].rgb.b;
			}
			val.type = LVT_COLOR;
			_colors.push_back(val);

			// Add max value
			val.clear();
			val.value = _max;
			val.type = OUTSIDE;
			_colors.push_back(val);
			// Add color max
			val.clear();
			if (!_outsideBlack) {
				// First color
				val.rgb.r = _colors[_colors.size()-2].rgb.r;
				val.rgb.g = _colors[_colors.size()-2].rgb.g;
				val.rgb.b = _colors[_colors.size()-2].rgb.b;
			}
			val.type = LVT_COLOR;
			_colors.push_back(val);


		}
		else {
			// Min is under _start
			if (minPos  < 0) {
				// Add at the begining color then min
				if (! _outsideBlack) {
					// First color
					val.rgb.r = _colors[1].rgb.r;
					val.rgb.g = _colors[1].rgb.g;
					val.rgb.b = _colors[1].rgb.b;
				}
				val.type = LVT_COLOR;
				_colors.insert(_colors.begin(), val);

				//min
				val.clear();
				val.value = _min;
				val.type  = OUTSIDE;
				_colors.insert(_colors.begin(), val);
			}
			// Max is over the end
			if (maxPos  > 0) {
				// Add max value then color
				val.clear();
				val.value = _max;
				val.type = OUTSIDE;
				_colors.push_back(val);

				// Add color
				val.clear();
				if (!_outsideBlack) {
					// Fisrt color
					val.rgb.r = _colors[_colors.size()-1].rgb.r;
					val.rgb.g = _colors[_colors.size()-1].rgb.g;
					val.rgb.b = _colors[_colors.size()-1].rgb.b;
				}
				val.type = LVT_COLOR;
				_colors.push_back(val);
			}
		}
	}

	void Legend::AddMinMaxGradient(void)
	{
		int minPos, maxPos;
		colormap_value val;
		MinMaxPos(minPos, maxPos);
		//mfp_debug("AddMinMaxGradient\n");
		if ((minPos < 0) && (maxPos < 0)) {
			// Both are below
			// Put max and min at the beginning of values with the first color enabled
			val.value = _max;
			if (!_outsideBlack) {
				// First color
				val.rgb.r = _colors[0].rgb.r;
				val.rgb.g = _colors[0].rgb.g;
				val.rgb.b = _colors[0].rgb.b;
			}
			val.type  = OUTSIDE;
			_colors.insert(_colors.begin(), val);

			// The same for min
			val.clear();
			val.value = _min;

			// 1st color already in
			val.rgb.r = _colors[0].rgb.r;
			val.rgb.g = _colors[0].rgb.g;
			val.rgb.b = _colors[0].rgb.b;
			val.type = OUTSIDE;
			_colors.insert(_colors.begin(), val);
		}
		else if ((minPos  > 0 ) && ( maxPos > 0)) {
			// Both are over the end
			// so, put them into values at the end and set last color
			// First for min
			val.value = _min;
			if (!_outsideBlack) {
				//last col
				val.rgb.r = _colors[_colors.size()-1].rgb.r;
				val.rgb.g = _colors[_colors.size()-1].rgb.g;
				val.rgb.b = _colors[_colors.size()-1].rgb.b;
			}
			val.type = OUTSIDE;
			_colors.push_back(val);

			// For max
			val.clear();
			val.value = _max;
			val.rgb.r = _colors[_colors.size()-1].rgb.r;
			val.rgb.g = _colors[_colors.size()-1].rgb.g;
			val.rgb.b = _colors[_colors.size()-1].rgb.b;
			val.type = OUTSIDE;
			_colors.push_back(val);
		}
		else {
			// When minimal value is below the start value
			if (minPos  < 0) {
				// Put a new color at the beginning
				// add color
				val.value = _min;
				if (!_outsideBlack) {
					//first color
					val.rgb.r = _colors[0].rgb.r;
					val.rgb.g = _colors[0].rgb.g;
					val.rgb.b = _colors[0].rgb.b;
				}
				val.type = OUTSIDE;
				_colors.insert(_colors.begin(), val);
			}

			if (maxPos  > 0) {
				//dodajemy na koncu
				val.clear();
				val.value = _max;
				if (!_outsideBlack) {
					//pierwszy kolor
					val.rgb.r = _colors[_colors.size()-1].rgb.r;
					val.rgb.g = _colors[_colors.size()-1].rgb.g;
					val.rgb.b = _colors[_colors.size()-1].rgb.b;
				}
				val.type = OUTSIDE;
				_colors.push_back(val);
			}
		}

		//mfp_debug("Ater minmaxgradient\n");
	}

bool Legend::IsBorder2(double value) const
{

    for (int i = 0; i < GetLiczbaGranic(); ++i)
    {

		if ( fvmath::Compare(value, GetWartoscGranicy(i)) )
        {
            return true;
        }

    }
    return false;

}
//
//bool Legend::IsBorder3(double value) const
//{
//
//    for (int i = 0; i < GetLiczbaGranic(); ++i)
//    {
//
//        //if ( value == GetWartoscGranicy(i) )
//        if ( value == GetWartoscGranicy(i) )
//        {
//            return true;
//        }
//
//    }
//    return false;
//
//}

int Legend::IsBorderNr(double value) const
{

    for (int i = 0; i < GetLiczbaGranic(); ++i)
    {

        //if ( value == GetWartoscGranicy(i) )
        if ( fvmath::Compare(value, GetWartoscGranicy(i)) )
        {
            return i;
        }

    }
    return -1;

}


//dla podanej wartosci zwraca mozliwe przedzialy w jakich sie znajduje
//dla granicy moga to byc dwa lewy lub pracy
//dla wartosci miedzy granicami jest to tylko jeden przedzial
//void Legend::GetRanges(const double & d, std::vector<unsigned int> & vUi)
std::vector<unsigned int> Legend::GetRanges(const double & d) const
{

    //bool border;
    int bnr;
    std::vector<unsigned int> vUi;

    vUi.clear();

    //sprawdzenie czy jest granica
    //border = IsBorder2(d);
    bnr = IsBorderNr(d);

    //#ifdef DEBUG4
    //    cout << "jest granica nr:" << bnr << endl;
    //#endif

    if (bnr >= 0)
    {
        //granica


        if (bnr == 0)//jesli granica zerowa
        {
            vUi.push_back(0);
        }
        else if (bnr == (GetLiczbaGranic()-1) )//jesli granica ostatnia
        {
            vUi.push_back(bnr-1);
        }
        else//jesli granica wewnetrzna
        {
            vUi.push_back(bnr-1);
            vUi.push_back(bnr);
        }

    }
    else
    {
        //wartosc wewnetrzna

        //jesli nie jest granica

        //petla po przedzialach legendy i zwrocenie numeru przedzialu w ktorym sie zawiera
        for (int i = 0; i < GetLiczbaGranic()-1; ++i)
        {
            if
            (
                ( d > GetWartoscGranicy(i) )
                &&
                ( d < GetWartoscGranicy(i+1) )
            )
            //if ( MathFunction::compare(value, GetWartoscGranicy(i)) )
            {
                vUi.push_back(i);
				break;
            }
        }//for (unsigned int i = 0; i < GetLiczbaGranic()-1; ++i)
    }//else

    return vUi;

}

size_t Legend::GetRangesSize(void) const
{

    if (_type == SOLID)
    {
        return ((_colors.size()-1)/2);
    }
    else if (_type == GRADIENT)
    {
        return (_colors.size()-1);
    }

    return(0);

}

void Legend::GetRangeValues(unsigned int i, double & v1, double & v2) const
{

    if (_type == SOLID)
    {
        //indekssy:
        //0 - 0, 2
        //1 - 2, 4
        //2 - 4, 6
        v1 = _colors[i*2].value;
        v2 = _colors[(i+1)*2].value;
    }
    else if (_type == GRADIENT)
    {
        //indeksy: i, i+1
        v1 = _colors[i].value;
        v2 = _colors[i+1].value;
    }

}

int Legend::GetRanges(const double & v1, const double & v2) const
{

    double sd, ed;

    for (unsigned int i = 0; i < GetRangesSize(); ++i)
    {

        GetRangeValues(i, sd, ed);

        if ( v1 < v2)
        {
            if ( (v1 >= sd) && ( v2 <= ed ) )
                return i;
        }
        else
        {
            if ( (v2 >= sd) && ( v1 <= ed ) )
                return i;
        }

    }
    //petla po przedzialach

        //jesli obie wartosci zawieraja sie w podanym przedziale to zwroceniejego numeru

    return -1;

}

bool Legend::GetColor(const double & value, const unsigned int & range, ColorRGB & c) const
{
	//sprawdzenie jesli po za zasiegiem to false

    #ifdef DEBUG2
	std::cout << "Legend::GetColor" << endl;
    #endif

    if (_type == SOLID)
    {
        c = _colors[range*2+1].GetColor();

        #ifdef DEBUG4
            cout << "wybrany kolor: " << c.AsString() << endl;
        #endif
    }
    else if (_type == GRADIENT)
    {
        double s, e, interp;
        //wartosc koloru poczatku
        //wartosc koloru konca

        s = _colors[range].value;
        e = _colors[range+1].value;

        interp = this->GetInterp(s, e, value);

        c.SetInterp(
            _colors[range].GetColor(),
            _colors[range+1].GetColor(),
            interp
            );

        //wartosc iterpolacji interpolacaja do granicy range do range+!
        //na podstawie wartoci kolorow oraz interpolacji wartosc koloru

    }
    return true;

}

std::vector<ColorRGB> Legend::GetColors(const std::vector<double> & vValues)
//void Legend::GetColors(const std::vector<double> & vValues, std::vector<rgb> & vRgb)
{

    std::vector<ColorRGB> vRgb;

    vRgb.clear();

    //sprawdzic czy podane wartosci mieszcza sie w jednym przedziale

    //wektor wektorow wartosci mozliwych przedzialow
    //rozmiar = 3

    std::vector< std::vector<unsigned int> > vv;

    std::vector<unsigned int> vIlosc;


    vv.resize( vValues.size() );

    for (unsigned int i = 0; i < vv.size(); ++i)
    {
        vv[i] = this->GetRanges(vValues[i]);
    }

    double max;

    #ifdef DEBUG4
        //cout << "jest granica nr:" << bnr << endl;

        for (unsigned int i = 0; i < vv.size(); ++i)
        {
            for (unsigned int j = 0; j < vv[i].size(); ++j)
            {
                cout << vv[i][j]<< ", ";
            }
            cout << endl;
        }

    #endif

    //znalezienie maksymalnej wartosci przedzialu do ktorego naleza wartosci
    max = vv[0][0];

    for (unsigned int i = 0; i < vv.size(); ++i)
    {
        for (unsigned int j = 0; j < vv[i].size(); ++j)
        {
            if (vv[i][j] > max)
                max = vv[i][j];
        }
    }

    vIlosc.resize(max+1);

    for (unsigned int i = 0; i < vv.size(); ++i)
    {
        for (unsigned int j = 0; j < vv[i].size(); ++j)
        {
            //if (vv[i][j] > max)
            //    max = vv[i][j];
            vIlosc[vv[i][j]] = vIlosc[vv[i][j]] + 1;
        }
    }

    #ifdef DEBUG4
        //cout << "jest granica nr:" << bnr << endl;

        for (unsigned int i = 0; i < vIlosc.size(); ++i)
        {
            cout << i << ": "<< vIlosc[i] << endl;
        }

    #endif

    //sprawdzenie ile jest elementow rownych ilosci elementow
    unsigned int x, xIndex;
    x = 0;
    for (unsigned int i = 0; i < vIlosc.size(); ++i)
    {
        if ( vValues.size() == vIlosc[i])
        {
            x++;
            xIndex = i;
        }
    }


    #ifdef DEBUG4
        cout << "liczba przedzialow wspolnych dla wszytskich wartosci: " << x << endl;
        if (x==1)
        {
            cout << "OK - dla przedzialu: " << xIndex << endl;
        }
        else
        {
            cout << "NIE OK - nie mo�na jednoznacznie okre�li� przedzia�u" << endl;
        }
    #endif


    if (x != 1)
        return vRgb;

    ColorRGB rgbTmp;

    //petla po  wartosciach wejsciowych
    for (unsigned int i = 0; i < vValues.size(); ++i)
    {

        if (this->GetColor(vValues[i], xIndex, rgbTmp))
        {
            //jesli pobrany dobrze kolor
            vRgb.push_back(rgbTmp);
        }
        else
        {
            //jesli blad pobrania koloru to podanie czarnego zeby nie bylo problemu z wyswietlaniem
        }



    }


    #ifdef DEBUG4
        cout << "wybrane kolory: " << endl;
        for (unsigned int i = 0; i < vRgb.size(); ++i)
        {
            cout << vRgb[i].AsString() << endl;
        }
    #endif


    //pobranie wszytskich wmozliwych warto�ci

    //zlikwidowanie powtorek

    //sprawdzenie ktora wartosc znajduje sie we wszytskich i czy jest to tylko jedna
    //dla kazdej przedzialu w w ilu elementach sie znajduje oraz czy jest tylko jeden ktory rowny jest max liczbie

    return vRgb;

}


//bool Legend::GetColors(GraphElement2<double> & ge2) const
//{
//    //dla pojedynczego elementu
//    std::vector<std::vector<unsigned int> > vv;
//    std::vector<unsigned int> vIlosc;
//    unsigned int max;
//
//
//    vv.resize( ge2.values.size() );
//
//    for (size_t i = 0; i < vv.size(); ++i)
//    {
//        vv[i] = this->GetRanges(ge2.values[i]);
//    }
//
//
//    #ifdef DEBUG2
//        cout << endl;
//        for (unsigned int i = 0; i < vv.size(); ++i)
//        {
//            for (unsigned int j = 0; j < vv[i].size(); ++j)
//            {
//                cout << vv[i][j] << ", ";
//            }
//            cout << endl;
//        }
//        cout << endl;
//    #endif
//
//    //znalezienie maksymalnej wartosci przedzialu do ktorego naleza wartosci
//    max = vv[0][0];
//
//    for (size_t i = 0; i < vv.size(); ++i)
//    {
//        for (size_t j = 0; j < vv[i].size(); ++j)
//        {
//            if (vv[i][j] > max)
//                max = vv[i][j];
//        }
//    }
//
//    vIlosc.resize(max+1);
//
//    for (unsigned int i = 0; i < vv.size(); ++i)
//    {
//        for (unsigned int j = 0; j < vv[i].size(); ++j)
//        {
//            //if (vv[i][j] > max)
//            //    max = vv[i][j];
//            vIlosc[vv[i][j]] = vIlosc[vv[i][j]] + 1;
//        }
//    }
//
//    //sprawdzenie ile jest elementow rownych ilosci elementow
//    unsigned int x, xIndex;
//    x = 0;
//    for (unsigned int i = 0; i < vIlosc.size(); ++i)
//    {
//        if ( ge2.values.size() == vIlosc[i])
//        {
//            x++;
//            xIndex = i;
//        }
//    }
//
//    #ifdef DEBUG2
//        cout << "x: " << x << endl;
//    #endif
//
//    //jesli nie udalo sie ustalic wspolnego przedzialu
//    if (x == 0)
//    {
//        return false;
//    }
//
//    //jesl byl tylko jeden przedzial - to wczensiej wybrany
//    //if (x == 1)
//    //    xIndex = ;
//
//
//    //jesli zostaly dwa przedzialy poniewaz wszytskie wartosci sa
//    //granicami lub bardzo blisko granicy to ponowne wybranie ale mniejszego
//    /*
//    else if (x == 2)
//    {
//        for (unsigned int i = 0; i < vIlosc.size(); ++i)
//        {
//            if ( ge2.values.size() == vIlosc[i])
//            {
//                xIndex = i;
//                break;
//            }
//        }
//    }
//    else
//    {
//        return false;
//    }
//    */
//
//    #ifdef DEBUG22
//        cout << "xIndex: " << xIndex << endl;
//    #endif
//
//    //rgb rgbTmp;
//
//
//    ge2.colors.clear();
//    ge2.colors.resize( ge2.values.size()  );
//
//
//    //petla po  wartosciach wejsciowych
//    for (unsigned int i = 0; i < ge2.values.size(); ++i)
//    {
//
//        this->GetColor( ge2.values[i], xIndex, ge2.colors[i] );
//
//        //if (this->GetColor(vValues[i], xIndex, rgbTmp))
//        //{
//            //jesli pobrany dobrze kolor
//        //    vRgb.push_back(rgbTmp);
//        //}
//        //else
//        //{
//            //jesli blad pobrania koloru to podanie czarnego zeby nie bylo problemu z wyswietlaniem
//        //}
//
//    }
//
//    return true;/
//}


//bool Legend::GetColors(std::vector< GraphElement2<double> > & vGe2) const
//{
//
//    for (unsigned int i = 0; i < vGe2.size(); ++i)
//    {
//        GetColors(vGe2[i]);
//            //return false;
//    }
//    return true;
//}


// Dla podanej wartosci zwraca gorna oraz dolna wartosc przedzialu - nowa wersja
bool Legend::GetRangesBorder(const double & value, double & s, double &e) const
{

    for (int i(0); i < GetLiczbaGranic()-1; ++i)
    {

        if
        (
            ( value > GetWartoscGranicy(i) )
            &&
            ( value < GetWartoscGranicy(i+1) )
        )
        //if ( MathFunction::compare(value, GetWartoscGranicy(i)) )
        {
            s = GetWartoscGranicy(i);
            e = GetWartoscGranicy(i+1);
            return true;
        }

    }
    return false;


}

void Legend::GetRange3(const unsigned int & range, double & sValue, ColorRGB & sColor, double & eValue, ColorRGB & eColor)
{

    if (_type == SOLID)
    {
        this->GetRangeColor(range, sValue, sColor, eValue, eColor);
    }
    else if (_type == GRADIENT)
    {
        this->GetRangeGradient(range, sValue, sColor, eValue, eColor);
    }

}

void Legend::GetRangeColor(const unsigned int & range, double & sValue, ColorRGB & sColor, double & eValue, ColorRGB & eColor)
{
    unsigned int r = 2 * range;

    sColor = _colors[r+1].GetColor();
    eColor = _colors[r+1].GetColor();

    sValue = _colors[r].value;
    eValue = _colors[r+2].value;
}

void Legend::GetRangeGradient(const unsigned int & range, double & sValue, ColorRGB & sColor, double & eValue, ColorRGB & eColor)
{

    //dla 0 : 0, 1

    //dla 1 : 1, 2


    //sColor = _colors[range].GetColor();
    //eColor = _colors[range+1].GetColor();

    sValue = _colors[range].value;
    eValue = _colors[range+1].value;

    unsigned int r = range;

    if ( (_colors[r].type & OUTSIDE) || (_colors[r+1].type & OUTSIDE))
    {

        if ( _colors[r].type & OUTSIDE )
        {
            sColor = _colors[r].GetColor();
            eColor = _colors[r].GetColor();
        }
        else if (_colors[r+1].type & OUTSIDE)
        {
            sColor = _colors[r+1].GetColor();
            eColor = _colors[r+1].GetColor();
        }
        //gradient czy kolor
        /*
        if (_type == SOLID)
        {
            sColor = _colors[r+1].GetColor();
            eColor = _colors[r+1].GetColor();
        }
        else if (_type == GRADIENT)
        {
            sColor = _colors[r+1].GetColor();
            eColor = _colors[r+1].GetColor();
        }
        */

    }
    else
    {
        sColor = _colors[r].GetColor();
        eColor = _colors[r+1].GetColor();
    }

}

/*
unsigned int Legend::GetIndex(double value)
{
    unsigned int r = 0;

    if ( _colors.size() < 2)
        return r;

    for (unsigned int i = 0; i < _colors.size()-1; ++i)
    {
        if ( (value >= _colors[i].value) && (value <= _colors[i+1].value) )
        {
            r = i+1;
            break;
        }
    }

    return r;

}*/

void Legend::GetColor(double value, colormap_value &lv)
{
    if ( _colors.size() < 2) return;
    lv.clear();
    std::vector<colormap_value>::const_iterator it = _colors.begin();
    std::vector<colormap_value>::const_iterator it_e = _colors.end() -1;
    //--it_e;
    if (it->value >= value) {
    	lv = *it;
    }
    else if (it_e->value <= value) {
    	lv = *it_e;
    }
    else {
    	bool find(false);
    	//std::cout.precision(numeric_limits<double>::digits10 + 1);
    	//mfp_debug("To sie wykonuje>>>>>\n");
    	std::vector<colormap_value>::const_iterator it_n(it);
    	for (++it_n; it != it_e; it = it_n, ++it_n) {
    	//for (unsigned int i = 0; i < _colors.size()-1; ++i)//  why is -1 ???? -> because of indexMax++
    	//{
    		//bool res1( !fvmath::Less( value, it->value));
    		//bool res2( !fvmath::Less)
			//std::cout<< " i = " << i << " " << res1 << " " << res2 << "val = " << value << " lv " << _colors[i+1].value<<" dif= " <<(_colors[i+1].value -value + fvmath::epsild) << std::endl;
			if ((it->value < value) && (it_n->value > value))
			//if ( (value >= legendValues[i].value) && (value <= legendValues[i+1].value) )
			{
				double delta   	  = it_n->value - it->value;
				double deltaValue = value - it->value;
				float interp 	  = float(deltaValue / delta);

				lv.value = value;
				lv.rgb.r = it->rgb.r + (it_n->rgb.r - it->rgb.r) * interp;
				lv.rgb.g = it->rgb.g + (it_n->rgb.g - it->rgb.g) * interp;
				lv.rgb.b = it->rgb.b + (it_n->rgb.b - it->rgb.b) * interp;
				find = !find;
				break;
			}
    	}

    	if (!find) {
    		//std::cout.precision(numeric_limits<double>::digits10 + 1);
    		std::cout<< "pocz: " << _colors[0].value << " kon: " << _colors[_colors.size()-1].value << std::endl;
    		std::cout<< "nie znalazłem dla val: " << value << std::endl;
    		std::cout<< "liczba granic: " << _colors.size() << std::endl;
    	}
    }// else
}

void Legend::GetColor(double value)
{

    GetColor(value, color);

    return;

}

double Legend::GetValue(const ColorRGB& rgb)
{

	// For given color in RGB mode get it HSV value
	ColorHSV hsv = ConvertRGBToHSV(rgb);
	//std::cout << rgb;
	//std::cout << static_cast<HSV>(hsv);
	// Get the range where H is foll in
	// We must get rab=ge and scale it over this range to get real value
	const float delta = 240.0f / (float)_nranges;
	float begin = 240.0f;
	float end = begin - delta;
	unsigned int index, next;
	for (unsigned int i = 0; i < _nranges; i++)
	{
		if (hsv.H <= begin && hsv.H >= end) {
			index = i;
			next  = i + 1;
			break;
		}

		begin = end;
		end -= delta;
	}

	double interp = (double)((begin - hsv.H) / delta);

	return (_colors[index].value + (_colors[next].value - _colors[index].value) * interp);
}

/*double Legend::GetInterp(double value)
{

    //sprawdzenie czy podana wartosc zawiera sie w zakresie
    if ( (value < legendValues[0].value) || (value > legendValues[legendValues.size()-1].value) )
        return -1.0;

    if ( legendValues.size() < 2)
        return -1.0;

    double delta;//wartoci przedzialu
    double deltaValue;//od poczatu przedzialu do podanej wartosci
    double interp;// wartosc interpolacji z zakresu 0-1

    unsigned int indexMin, indexMax;

    bool find = false;

    for (unsigned int i = 0; i < legendValues.size()-1; ++i)
    {
        if ( (value >= legendValues[i].value) && (value <= legendValues[i+1].value) )
        {
            indexMin = i;
            indexMax = i+1;
            find = true;
            break;
        }
    }

    if (!find)
        return -1.0;

    delta = legendValues[indexMax].value - legendValues[indexMin].value;

    deltaValue = value - legendValues[indexMin].value;

    interp = deltaValue/delta;

    return interp;

}
*/

double Legend::GetInterp(double a, double b, double value) const
{

    //if (value < a )
    //    return -1.0;

    //if (value > b )
    //    return -1.0;

    double delta = 0.0;//wartoci przedzialu
    double deltaValue;//od poczatu przedzialu do podanej wartosci
    double interp = 0.0;// wartosc interpolacji z zakresu 0-1

    delta = b - a;
    deltaValue = value - a;
	
	if(delta != 0.0)
		interp = deltaValue/delta;

    return interp;

}


double Legend::GetValue(double a, double b, double interp) const
{

    //this->x = a.x + (b.x-a.x) * interp;
    return ( a + (b-a) * interp );

}
/*
void Legend::GetInterval(double valueStart, double valueEnd)
{

    colorRange.clear();

    unsigned int rangeStart, rangeEnd;

    bool findStart = false;
    bool findEnd = false;

	//int var=0;
	//int vars=0;

	
    for (unsigned int i = 0; i < legendValues.size()-1; ++i)
    {
		//printf("ws= %.16lf lvi= %.16lf lvi1= %.16lf\n",valueStart,legendValues[i].value,legendValues[i+1].value);
        //if ( (valueStart >= legendValues[i].value) && (fabs(valueStart - legendValues[i+1].value) >= 0.000000000000001 ))
			if ( (valueStart >= legendValues[i].value) && (valueStart <= legendValues[i+1].value) )
        {
            rangeStart = i+1;
			//vars++;
            findStart = true;
        }

		//printf("we= %.16lf lvi= %.16lf lvi1= %.16lf\n",valueEnd,legendValues[i].value,legendValues[i+1].value);
		if ( (valueEnd >= legendValues[i].value) && (valueEnd <= legendValues[i+1].value ))
        {
			
            rangeEnd = i+1;
			//var++;
            findEnd = true;
        }

        if ( findStart && findEnd )
            break;
    }

	//if(!var){
	//	printf("stop: liczba= %16.lf",1e-15);
		//std::cout<< "vs 
	//}

	//if(!vars){
	//	printf("stop: liczba= %16.lf",1e-15);
	//	//std::cout<< "vs 
	//}

    //if ( rangeStart == rangeEnd)
    //{
    //    GetColor(valueStart, color);
    //    colorRange.push_back(color);
    //    GetColor(valueEnd, color);
    //    colorRange.push_back(color);
    //}

    LegendValue colorTmp;

    GetColor(valueStart, colorTmp);
    colorRange.push_back(colorTmp);

    //petla po progach miedzy watosciami

	//std::cout << "start: " << rangeStart << std::endl;
	//std::cout << "end: " << rangeEnd << std::endl;

    if (rangeStart < rangeEnd)
    {
		//cout << "rangeend: " << endl;
        for (unsigned int j = rangeStart; j < rangeEnd; ++j)
        {
            GetColor(legendValues[j].value, colorTmp);
            colorRange.push_back(colorTmp);
        }
    }

    if (rangeStart > rangeEnd)
    {
        //cout << "a" << endl;
        for (unsigned int j = rangeStart; j > rangeEnd; j--)
        {
            //cout << "j: " << j << endl;
            GetColor(legendValues[j-1].value, colorTmp);
            colorRange.push_back(colorTmp);
        }
    }

    GetColor(valueEnd, colorTmp);
    colorRange.push_back(colorTmp);

    return;

}
*/
/*
bool Legend::DCompare(double x, double y)
//bool equal(double x, double y)
{
    const double EPS = 0.000001;
    const double f = fabs(x - y);

    if ( f <= EPS)
        {return true;}
    else
        {return false;}
}
*/

/*bool Legend::IsBorder(double value)
{


    for (unsigned int i = 0; i < _colors.size()-1; ++i)
    {
			if (fvmath::Compare(value, _colors[i].value))
        //if (compare(value, legendValues[i].value))
            return true;
    }
    return false;

    //return (MathFunction::compare(value, legendValues[i].value) ? true : false);

}*/



/*
void Legend::GetColor(double value)
{

    color.value = 0;
    color.red = 0;
    color.green = 0;
    color.blue = 0;

    if ( legendValues.size() < 2)
        return;

    double min, max;
    double delta;//wartoci przedzialu
    double deltaValue;//od poczatu przedzialu do podanej wartosci
    double interp;// wartosc interpolacji z zakresu 0-1

    unsigned int indexMin, indexMax;

    bool find = false;

    for (unsigned int i = 0; i < legendValues.size()-1; ++i)
    {
        if ( (value >= legendValues[i].value) && (value <= legendValues[i+1].value) )
        {
            indexMin = i;
            indexMax = i+1;
            find = true;
            break;
        }
    }

    if (!find)
        return;

    delta = legendValues[indexMax].value - legendValues[indexMin].value;

    deltaValue = value - legendValues[indexMin].value;

    interp = deltaValue/delta;


    #ifdef DEBUG
        cout << value << endl;
        cout << interp << endl;
    #endif

    color.value = value;
    #ifdef DEBUG
        cout << "wyliczenie red: " << endl;
        cout << legendValues[indexMin].red << " + (" << legendValues[indexMax].red << " - " << legendValues[indexMin].red << ") * " << interp << endl;
        cout << "wyliczenie green: " << endl;
        cout << legendValues[indexMin].green << " + (" << legendValues[indexMax].green << " - " << legendValues[indexMin].green << ") * " << interp << endl;
        cout << "wyliczenie blue: " << endl;
        cout << legendValues[indexMin].blue << " + (" << legendValues[indexMax].blue << " - " << legendValues[indexMin].blue << ") * " << interp << endl;
    #endif
    color.red = legendValues[indexMin].red + (legendValues[indexMax].red - legendValues[indexMin].red) * interp;
    color.green = legendValues[indexMin].green + (legendValues[indexMax].green - legendValues[indexMin].green) * interp;
    color.blue = legendValues[indexMin].blue + (legendValues[indexMax].blue - legendValues[indexMin].blue) * interp;

    //sprawdzenie w ktorym przedziale jest wartosc
    //przygotowanie do interpolacji
    //wykonanie interpolacji
    //zapis podanej wartosci oraz kolorow jej odpowiadajacych

    return;
}
*/



/*
void Legend::SetMinMax2(double min, double max)
{

    #ifdef FV_DEBUG
        std::cout << std::endl << "Legend::SetMinMax2" << endl;
    #endif

    unsigned int rangeNumber;
    rangeNumber = 4;

    //***************************************

    double delta, deltaPart;
    //, division;
    unsigned int colorNumber;//liczba kolorow potrzebna dla podanej ilosci przedzialow w zaleznosci od rodzaju legendy
    unsigned int colorDiv;//liczba przez jaka nalezy podzilic 1 alby osiagnac dana liczbe kolorow

    //jajco idzie brzegiem lasu

    delta = max - min;

    deltaPart = delta/rangeNumber;

    float colorDelta;//przedzial miedzy dwoma granicami

    //this->_type = SOLID;
    this->_type = GRADIENT;

    //jesli kazdy przedzial jest
    if (this->_type == SOLID)
    {
        //jednilitym kolorem
        colorNumber = rangeNumber;
        colorDiv = colorNumber-1;

        #ifdef DEBUG
            cout << "Kolor" << endl;
        #endif
    }
    else if (this->_type == GRADIENT)
    {
        //w zakresie barwy od do
        colorNumber = rangeNumber+1;
        colorDiv = colorNumber-1;

        #ifdef DEBUG
            cout << "Gradient" << endl;
        #endif

    }

    #ifdef DEBUG
        cout << "6 / " << colorDiv << " = " << 6.0 / colorDiv << endl;
    #endif

    colorDelta = 6.0f / colorDiv;

    if (this->_type == SOLID)
    {

        for (unsigned int i = 0; i < rangeNumber; ++i)
        {

            if (i < rangeNumber-1)
            {
            #ifdef DEBUG
        cout << i * colorDelta  << " - " << i * colorDelta  << ":" << min + i * deltaPart  << " - " << min + (i+1) * deltaPart << endl;
            #endif
            }
            else
            {
            #ifdef DEBUG
        cout << i * colorDelta  << " - " << i * colorDelta  << ":" << min + i * deltaPart  << " - " << max << endl;
            #endif
            }
        }



    }
    else if (this->_type == GRADIENT)
    {

        for (unsigned int i = 0; i <= rangeNumber; ++i)
        {
            #ifdef DEBUG
                cout << i * colorDelta << endl;
            #endif
        }

        for (unsigned int i = 0; i < rangeNumber; ++i)
        {
            #ifdef DEBUG
cout << i * colorDelta << " - " << (i+1) * colorDelta << " : " << min + i * deltaPart  << " - " << min + (i+1) * deltaPart << endl;
            #endif
        }
    }

    //this->_type = SOLID;

}

*/
} // end namespace FemViewer
