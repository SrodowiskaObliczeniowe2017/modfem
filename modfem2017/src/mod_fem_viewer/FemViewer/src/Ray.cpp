/*
 * Ray.cpp
 *
 *  Created on: 28 sty 2014
 *      Author: Paweł Macioł
 */
#include "Ray.h"

namespace FemViewer {
using namespace fvmath;

int64_t uniqueRayId = 0;

template<typename T>
Ray<T>::Ray(Vec3<T> orig, Vec3<T> dir, T near, T far) : debug(false),
	orig(orig), dir(dir), tmin(near), tmax(far)
{
	this->id = uniqueRayId;
	__sync_fetch_and_add(&uniqueRayId, 1);
	invdir = T(1) / dir;
	sign[0] = (invdir.v[0] < 0);
	sign[1] = (invdir.v[1] < 0);
	sign[2] = (invdir.v[2] < 0);
}

template<> Ray<float>::Ray(Vec3<float> orig, Vec3<float> dir, float near, float far);

}// end namespace FemViewer



