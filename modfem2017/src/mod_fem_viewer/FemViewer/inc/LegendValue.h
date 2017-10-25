#ifndef _LEGEND_VALUE_H_
#define _LEGEND_VALUE_H_

#include <iostream>
#include <sstream>
#include <string>
#include <string.h>

#include "Color.h"

namespace FemViewer {

/* Types of legend items */
enum LegendValueType
{
    NON = 0, //nieoznaczony
    BORDER = 1,
    INSIDE = 2,//elementy bedace granica wewnatrz start end
    OUTSIDE = 4,//elementy bedace granica po za start end

    START = 8,
    END = 16,

    LVT_COLOR = 32

};

/* Item of a Legend */
typedef struct colormap_value {
	/* Value */
    double value;
    /* Color */
    RGBColor rgb;
    /* Type */
    int type;
    /*!
    *    A constructor.
    */
    explicit colormap_value(void) { this->clear(); };
    colormap_value(const colormap_value& ref)
    : value(ref.value)
    , rgb(ref.rgb)
    , type(ref.type)
    {}
    void operator=(const colormap_value& rhs) {
    	value = rhs.value;
    	rgb.r = rhs.rgb.r;
    	rgb.g = rhs.rgb.g;
    	rgb.b = rhs.rgb.b;
    	type  = rhs.type;
    }
    void clear(void) {
        //value = -1.0;
        //rgb.r = 0.0f;
		//rgb.g = 0.0f;
		//rgb.b = 0.0f;
		//type  = NON;
		memset(this,0,sizeof(*this));
    }

    //zwraca wartosci z zakresu 0-1
    //float RedAsFloat(void){return red;};

    //zwraca warto�ci z zakresu 0-255
    //unsigned int RedAsUInt(void){return (unsigned int)(red*255);};

    unsigned char RedAsUChar(void){return static_cast<unsigned char>(rgb.r*255.99);};
    unsigned char GreenAsUChar(void){return static_cast<unsigned char>(this->rgb.g*255.99);};
    unsigned char BlueAsUChar(void){return static_cast<unsigned char>(this->rgb.b*255.99);};

    //! Pobiera dane koloru wartosci legendy jako kolor
	ColorRGB GetColor(void) const
    {
        return ColorRGB(rgb.r, rgb.g, rgb.b);
    };

    //! Ustawia kolor
    void SetColor( const ColorRGB & c )
    {
        rgb.r = c.R;
        rgb.g = c.G;
        rgb.b = c.B;
    };

    /*! \fn std::string TypeAsString(void)
    *   \brief Zwraca nazwe typu elementu
    *
    *   Zwraca nazwe typu elementu.
    *
    *   \return nazwa typu elementu
    */
    std::string TypeAsString(void)
    {
        std::stringstream oss;

        if ( this->type == NON)
            oss << "NON ";

        if ( this->type & BORDER)
            oss << "BORDER ";

        if ( this->type & INSIDE)
            oss << "INSIDE ";

        if ( this->type & OUTSIDE)
            oss << "OUTSIDE ";

        if ( this->type & LVT_COLOR)
            oss << "COLOR ";

        if ( this->type & START)
            oss << "START ";

        if ( this->type & END)
            oss << "END ";

        return oss.str();
    };

    //! Zwraca obiekt jako string
    std::string AsString(void)
    {
        std::stringstream oss;

        oss << "war: " << value << " r: " << (int)RedAsUChar() << " g: " << (int)GreenAsUChar() << " b: "<< (int)BlueAsUChar() << " t: " << TypeAsString();

        return oss.str();
    };
} colormap_value;






} // end namespace FemViewer
#endif /* _LEGEND_VALUE_H_
*/
