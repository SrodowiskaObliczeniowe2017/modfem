#ifndef _COLOR_H_
#define _COLOR_H_

#include<iostream>
template<typename Type>
inline Type Clump2Type(Type el, Type Min, Type Max)
{
	if( el < Min) return Min;
	else if( el > Max) return Max;
	else return el;
}

namespace FemViewer {

struct RGBColor {
	union { struct { float r, g, b; }; float v[3]; };
};

struct RGBAColor : public RGBColor{
	float a;
};
struct ColorHSV {
	float H, S, V;
};

class HSV : public ColorHSV {
public:
	HSV(float h = 0.f,float s = 0.f,float v=0.f) {
		H = h; S = s; V = v;
	}
	HSV(const ColorHSV& hsv) { H = hsv.H; S = hsv.S; V = hsv.V; }
	friend std::ostream& operator<< (std::ostream& os,const HSV hsv);
};

inline std::ostream& operator<< (std::ostream& os,const HSV hsv) {
	return os << "Hue = " << hsv.H << " Sat = " << hsv.S << " Val = " << hsv.V << std::endl;
}

class ColorRGB {
	public:
		float R, G, B;

		explicit ColorRGB(float r = 0.0f, float g = 0.0f, float b = 0.0f)
			: R(r), G(g), B(b) 
		{
			this->Clump2float();
		}

		explicit ColorRGB(const float* col)
		{
			R = *col++;
			G = *col++;
			B = *col;
		}

		operator float*()  {
			return &this->R;
		}

		ColorRGB(const ColorRGB& refCol);
		virtual ~ColorRGB(){}

		void SetInterp(const ColorRGB& a, const ColorRGB& b, const double interp)
        {
            R = static_cast<float>(a.R+ (b.R-a.R) * interp);
            G = static_cast<float>(a.G+ (b.G-a.G) * interp);
            B = static_cast<float>(a.B+ (b.B-a.B) * interp);
        }


		ColorRGB& operator= (const ColorRGB& col);

		void outcolor(const char * name) const;

		friend std::ostream& operator<< (std::ostream& os, const ColorRGB& rhs);

	protected:
		virtual void Clump2float();

	};

	inline std::ostream& operator<< (std::ostream& os, const ColorRGB& rhs)
	{
		return os << "RGB = ("<< int(rhs.R * 255) << ", " << int(rhs.G * 255) << ", " << int(rhs.B * 255) << ")\n";
	}

	class ColorRGBA : public ColorRGB {
	public:
		float A;

		ColorRGBA(float r = 0.0f, float g = 0.0f, float b = 0.0f, float a = 1.0f)
			: ColorRGB(r,g,b), A(a)
		{}

		ColorRGBA(const ColorRGBA& rhs) : ColorRGB(rhs), A(rhs.A)
		{}

		ColorRGBA& operator=(const ColorRGBA& rhs) {
			if (this == &rhs) return *this;
			ColorRGB::operator=(rhs);
			A = rhs.A;
		}
	};

	ColorRGB ConvertHSVToRGB(const ColorHSV& hsv);
	ColorHSV ConvertRGBToHSV(const ColorRGB& rgb);

} // end namespace FemViewer

#endif /* _COLOR_H_
*/
