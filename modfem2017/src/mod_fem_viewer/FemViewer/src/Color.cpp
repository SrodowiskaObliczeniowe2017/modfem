#include "Color.h"
#include <cstdio>
#include <cmath>

template<typename Ty>
inline Ty min_of(const Ty& a, const Ty& b, const Ty& c)
{
	Ty _min = (a < b) ? a : b;
	return ( _min < c) ? _min : c;
}

template<typename Ty>
inline Ty max_of(const Ty& a, const Ty& b, const Ty& c)
{
	Ty _min = (a > b) ? a : b;
	return ( _min > c) ? _min : c;
}

namespace FemViewer {

	ColorRGB::ColorRGB(const ColorRGB& col)
	{
		R = col.R; G = col.G; B = col.B;
		this->Clump2float();
	}

	ColorRGB& ColorRGB::operator =(const ColorRGB& col)
	{
		this->R = col.R;
		this->G = col.G;
		this->B = col.B;
		return(*this);	
	}

	void ColorRGB::Clump2float()
	{	
		this->R = Clump2Type<float>(this->R,0.0f,1.0f);
		this->G = Clump2Type<float>(this->G,0.0f,1.0f);
		this->B = Clump2Type<float>(this->B,0.0f,1.0f);
	}

	void ColorRGB::outcolor(const char * name) const
	{
		(void)printf("%s (%f, %f, %f)\n",name,R,G,B);
	}

	ColorRGB ConvertHSVToRGB( const ColorHSV& hsv)
	{
		float hue, sat, val;

		hue = hsv.H;
		sat = hsv.S;
		val = hsv.V;

		ColorRGB sRgb;

		float red, grn, blu;
		float i, f, p, q, t;

		if(val==0) {
				red = 0.0f;
				grn = 0.0f;
				blu = 0.0f;
		} else {
				hue/=60.0f;
				i = floor(hue);
				f = hue-i;
				p = val*(1-sat);
				q = val*(1-sat*f);
				t = val*(1-sat*(1-f));
				if (i==0) {red=val; grn=t; blu=p;}
				else if (i==1) {red=q; grn=val; blu=p;}
				else if (i==2) {red=p; grn=val; blu=t;}
				else if (i==3) {red=p; grn=q; blu=val;}
				else if (i==4) {red=t; grn=p; blu=val;}
				else if (i==5) {red=val; grn=p; blu=q;}
		}
		sRgb.R = red;
		sRgb.G = grn;
		sRgb.B = blu;
		return sRgb;

	}
	
	ColorHSV ConvertRGBToHSV(const ColorRGB& rgb)
	{
		float _min = min_of<float>(rgb.R,rgb.G,rgb.B);
		float _max = max_of<float>(rgb.R,rgb.G,rgb.B);
		float _dlt = _max - _min;

		ColorHSV _hsv;
		_hsv.V = _max;

		if ( _max != 0.0f)
			_hsv.S = _dlt / _max;		// s
		else {
			// r = g = b = 0		// s = 0, v is undefined
			_hsv.S = 0.0f;
			_hsv.H = -1.0f;
			return _hsv;
		}

		if ( rgb.R == _max )
			_hsv.H = ( rgb.G - rgb.B ) / _dlt;		// between yellow & magenta
		else if( rgb.G == _max )
			_hsv.H = 2.0f + ( rgb.B - rgb.R ) / _dlt;	// between cyan & yellow
		else
			_hsv.H = 4.0f + ( rgb.R - rgb.G ) / _dlt;	// between magenta & cyan

		_hsv.H  *= 60.0f;				// degrees

		if( _hsv.H < 0.0f )
			_hsv.H += 360.0f;

		return _hsv;
	}


} // end namespace FemViewer
