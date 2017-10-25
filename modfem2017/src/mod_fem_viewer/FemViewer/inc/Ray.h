/*
 * Ray.h
 *
 *  Created on: 28 sty 2014
 *      Author: Paweł Macioł
 */

#ifndef _Ray_h__
#define _Ray_h__

#include "MathHelper.h"
#include <stdint.h> //uint64_t

namespace FemViewer {

using namespace fvmath;

template<typename T>
class IsectInfo
{
  public:
	IsectInfo() {}
	CVec3<T> P, N;
	CVec3<T> dPdu, dPdv;
	T s, t;
	CVec3<T> worldToLocal(const Vec3<T> &v) {
		CVec3<T> sn(dPdu);
		(void)Normalize(sn);
		CVec3<T> tn = N * sn;
		return CVec3<T>(Dot(v, sn), Dot(v, tn), Dot(v, N));
	}
};

template<typename T>
class Ray {
  public:
	bool debug;
	mutable Vec3f color;
	CVec3<T> orig, dir;		/// ray orig and dir
	mutable T tmin, tmax;	/// ray min and max distances
	unsigned triangleId;	/// used with triangle mesh (id of the intersected triangle)
	CVec3<T> invdir; 		/// precomputed for ray-box intersection
	int sign[3];			/// precomputed for ray-box intersection
	uint64_t id;			/// unique ray id
	Ray(Vec3<T> orig, Vec3<T> dir, T near = T(0), T far = std::numeric_limits<T>::max());
	CVec3<T> operator () (const T &t) const { return orig + dir * t; }
};

}// end namespace


#endif /* _Ray_h__ */
