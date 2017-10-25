#include "MathHelper.h"
#include <iostream>

namespace FemViewer {
namespace fvmath {




//template double Normalize<double>(Vec3d& v);


template<typename T>
bool Compare(const T& arg1,const T& arg2)
{
	#ifdef FV_DEBUG
	std::cout << "Compare: " << arg1 << " " << arg2 << std::endl;
	#endif
	return( arg1 == arg2);
}

template<>
bool Compare<float>(const float& arg1,const float& arg2)
{
	//#ifdef FV_DEBUG
	//std::cout << "Compare: " << arg1 << " " << arg2 << std::endl;
	//#endif

	const float epsilf = 0.000001f;
	const float delta = fabsf(arg1 - arg2);
	return ( epsilf >= delta) ? true : false;
}

template<>
bool Compare<double>(const double& arg1,const double& arg2)
{
    const double epsil = 0.0000000000001;
	const double delta = fabs(arg1 - arg2);
	return ( epsil >= delta) ? true : false;
}

template<>
bool Compare<CVec3f>(const CVec3f& v1,const CVec3f& v2)
{
	bool result = ( (v1.x == v2.x) &&
					(v1.y == v2.y) &&
					(v1.z == v2.z) );
	return result;
}


template<>
bool Compare<Vec3d>(const Vec3d& v1,const Vec3d& v2)
{
	bool result = ( (v1.x == v2.x) &&
					(v1.y == v2.y) &&
					(v1.z == v2.z) );
	return result;
}

template<class T1, class T2>
bool Compare2(const T1& arg1, const T2& arg2)
{
	return (arg1.x == arg2 && arg1.y == arg2 && arg1.z == arg2);
}


template<>
bool Compare2<Vec3d,double>(const Vec3d& v,const double& val)
{
	return (v.x == val && v.y == val && v.z == val);
}

const float epsil = 0.000001f;
const double epsild = 0.00000000000001;
const double ref = 1.4142135623730950;

template<> bool Less(const double& arg1,const double& arg2) {
	double tmp = arg1 - arg2 + epsild;
	return (tmp < 0.0);
}

template<> bool Less(const float& arg1,const float& arg2) {
	double tmp = arg1 - arg2 + epsil;
	return (tmp < 0.f);
}
} // end namespace fvmath
} // end namespace FemViewer
