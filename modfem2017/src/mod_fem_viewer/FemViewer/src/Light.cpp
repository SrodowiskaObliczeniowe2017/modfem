/*
 * Light.cpp
 *
 *  Created on: 14 sie 2014
 *      Author: dwg
 */

#include "Light.h"
namespace FemViewer {

Light::Light()
: m_type(Fixed)
, m_color(1.0f, 1.0f, 1.0f, 1.0f)
, m_ambientIntensity(0.3f, 0.3f, 0.3f)
, m_diffuseIntensity(0.7f, 0.7f, 0.7f)
, m_position(3.f, 3.f, 3.f)
{
	//mfp_log_debug("Light default ctr");
}

Light::Light(const Light& rhs)
: m_type(rhs.m_type)
, m_color(rhs.m_color)
, m_ambientIntensity(rhs.m_ambientIntensity)
, m_diffuseIntensity(rhs.m_diffuseIntensity)
, m_position(rhs.m_position)
{
	//mfp_log_debug("Light copy ctr");
}

Light& Light::operator=(const Light& rhs)
{
	//mfp_log_debug("Light operator=");
	m_type = rhs.m_type;
	m_color = rhs.m_color;
	m_ambientIntensity = rhs.m_ambientIntensity;
	m_diffuseIntensity = rhs.m_diffuseIntensity;
	m_position = rhs.m_position;
	return *this;
}

}// end namespace



