#ifndef _Homogenous_H_
#define _Homogenous_H_

#include "../../include/fv_limits.h"
#include <float.h>
#include <iostream>
#include <cmath>

class Homogenous
{
public:
	// components
	float x, y, z, w;

	// default constructor
	Homogenous(): x(0.0), y(0.0), z(0.0), w(1.0) {}

	// special constructor
	Homogenous(float x, float y, float z = 0.0, float w = 1.0): x(x), y(y), z(z), w(w) {}
	
	// add
	Homogenous operator + (const Homogenous &hc) const
	{
		float x1 = x, y1 = y, z1 = z;
		if (w)
		{
			x1 /= w;
			y1 /= w;
			z1 /= w;
		}
		
		float x2 = hc.x, y2 = hc.y, z2 = hc.z;
		if (hc.w)
		{
			x2 /= hc.w;
			y2 /= hc.w;
			z2 /= hc.w;
		}

		return Homogenous(x1 + x2, y1 + y2, z1 + z2, (w == 0.0 || hc.w == 0.0) ? 0.0f : 1.0f);
	}

	// sub
	Homogenous operator - (const Homogenous &hc) const
	{
		float x1 = x, y1 = y, z1 = z;
		if (w)
		{
			x1 /= w;
			y1 /= w;
			z1 /= w;
		}
		
		float x2 = hc.x, y2 = hc.y, z2 = hc.z;
		if (hc.w)
		{
			x2 /= hc.w;
			y2 /= hc.w;
			z2 /= hc.w;
		}

		return Homogenous(x1 - x2, y1 - y2, z1 - z2, (w == 0.0 || hc.w == 0.0) ? 0.0f : 1.0f);
	}

	// dot product
	float operator * (const Homogenous &hc) const
	{
		if (w == 0.0 || hc.w == 0.0)
			return FLT_MAX;

		float x1 = x, y1 = y, z1 = z;
		x1 /= w;
		y1 /= w;
		z1 /= w;
	
		float x2 = hc.x, y2 = hc.y, z2 = hc.z;
		x2 /= hc.w;
		y2 /= hc.w;
		z2 /= hc.w;

		return x1 * x2 + y1 * y2 + z1 * z2;
	}

	// cross product
	Homogenous operator % (const Homogenous &hc) const
	{
		float x1 = x, y1 = y, z1 = z;
		if (w)
		{
			x1 /= w;
			y1 /= w;
			z1 /= w;
		}
		
		float x2 = hc.x, y2 = hc.y, z2 = hc.z;
		if (hc.w)
		{
			x2 /= hc.w;
			y2 /= hc.w;
			z2 /= hc.w;
		}

		return Homogenous(y1 * z2 - z1 * y2, z1 * x2 - x1 * z2, x1 * y2 - y1 * x2, (w == 0.0 || hc.w == 0.0) ? 0.0f : 1.0f);
	}

	// add
	Homogenous& operator += (const Homogenous &hc)
	{
		float x1 = x, y1 = y, z1 = z;
		if (w)
		{
			x1 /= w;
			y1 /= w;
			z1 /= w;
		}
		
		float x2 = hc.x, y2 = hc.y, z2 = hc.z;
		if (hc.w)
		{
			x2 /= hc.w;
			y2 /= hc.w;
			z2 /= hc.w;
		}

		x = x1 + x2;
		y = y1 + y2;
		z = z1 + z2;
		w = (w == 0.0 || hc.w == 0.0) ? 0.0f : 1.0f;

		return *this;
	}

	// sub
	Homogenous& operator -= (const Homogenous &hc)
	{
		float x1 = x, y1 = y, z1 = z;
		if (w)
		{
			x1 /= w;
			y1 /= w;
			z1 /= w;
		}
		
		float x2 = hc.x, y2 = hc.y, z2 = hc.z;
		if (hc.w)
		{
			x2 /= hc.w;
			y2 /= hc.w;
			z2 /= hc.w;
		}

		x = x1 - x2;
		y = y1 - y2;
		z = z1 - z2;
		w = (w == 0.0 || hc.w == 0.0) ? 0.0f : 1.0f;

		return *this;
	}

	// scaling from right
	Homogenous& operator *= (float lambda_xyz)
	{
		x *= lambda_xyz;
		y *= lambda_xyz;
		z *= lambda_xyz;
		return *this;
	}

	// scaling from right
	Homogenous& operator /= (float lambda_xyz)
	{
		x /= lambda_xyz;
		y /= lambda_xyz;
		z /= lambda_xyz;
		return *this;
	}

	// cross product
	Homogenous& operator %= (const Homogenous &hc)
	{
		float x1 = x, y1 = y, z1 = z;
		x = y1*hc.z - z1*hc.y;
		y = z1*hc.x - x1*hc.z;
		z = x1*hc.y - y1*hc.x;
		return(*this);
	}

	// length
	float Magnitude() const
	{
		if (w == 0.0)
			return FLT_MAX;

		float w2 = w * w;
		return sqrt((x * x + y * y + z * z) / w2);
	}

	// convert to unit vector
	void Normalize()
	{
		float length = Magnitude();
		x /= length;
		y /= length;
		z /= length;
	}
};

// scaling from left
inline Homogenous operator * (float lhs, const Homogenous &rhs)
{
	return Homogenous(lhs*rhs.x, lhs*rhs.y, lhs*rhs.z, lhs*rhs.w); 
}

// scaling from right
inline Homogenous operator * (const Homogenous &lhs, float rhs)
{
	return Homogenous(lhs.x*rhs, lhs.y*rhs, lhs.z*rhs, lhs.w*rhs); 
}

// output
inline std::ostream& operator << (std::ostream &lhs, const Homogenous &hc)
{
	return lhs << hc.x << " " << hc.y << " " << hc.z << " " << hc.w;
}

// input
//inline std::istream& operator >> (std::istream &lhs, Homogenous &hc);


#endif