/*
 * Light.h
 *
 *  Created on: 14 sie 2014
 *      Author: dwg
 */

#ifndef LIGHT_H_
#define LIGHT_H_

#include "Log.h"
#include "Color.h"
#include "MathHelper.h"

namespace FemViewer {

using namespace fvmath;

class Light {

public:
		enum {
			Flat = 0,
			Camera,
			Fixed
		};

 public:

	explicit Light();
    explicit Light(const Light& rhs);
	Light& operator=(const Light& rhs);
	~Light(){
		//mfp_debug("Light detr");
	}

		  int& Type() 	   { return m_type; }
	const int& Type() const { return m_type; }
	      ColorRGBA& Color()       { return m_color; }
	const ColorRGBA& Color() const { return m_color; }
		  CVec3f& AmbientIntensity()       { return m_ambientIntensity; }
	const CVec3f& AmbientIntensity() const { return m_ambientIntensity; }
		  CVec3f& DiffuseIntensity()       { return m_diffuseIntensity; }
	const CVec3f& DiffuseIntensity() const { return m_diffuseIntensity; }
		  CVec3f& Position()       { return m_position; }
    const CVec3f& Position() const { return m_position; }

 private:

    int m_type;
	ColorRGBA m_color;
	CVec3f m_ambientIntensity;
	CVec3f m_diffuseIntensity;
	CVec3f m_position;
};


} //

#endif /* LIGHT_H_ */
