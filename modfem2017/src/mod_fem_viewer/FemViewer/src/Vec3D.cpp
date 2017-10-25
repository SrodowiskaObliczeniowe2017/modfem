#include "Vec3D.h"
#include <cstdio>
#include <cmath>
#include <stdexcept>

#define min_(a,b) (((a < b)) ? (a) : (b))
#define max_(a,b) (((a > b)) ? (a) : (b))

namespace FemViewer {

Vec3D operator -(const Vec3D& vec_)
{
	return Vec3D(-vec_[0],-vec_[1],-vec_[2]);
}
bool operator ==(const Vec3D& vec_1, const Vec3D& vec_2)
{
	return ( vec_1.x[0] == vec_2.x[0] && 
		     vec_1.x[1] == vec_2.x[1] && 
			 vec_1.x[2] == vec_2.x[2] );
}

bool operator !=(const Vec3D& vec_1, const Vec3D& vec_2)
{
	return ( vec_1.x[0] != vec_2.x[0] || 
		     vec_1.x[1] != vec_2.x[1] || 
			 vec_1.x[2] != vec_2.x[2] );
}

float Vec3D::quadNorm() const
{
	return(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
}

float Vec3D::Norm() const
{
	return sqrt(quadNorm());
}

void Vec3D::normalize()
{
	float length = Norm();
	if( length > 0.0f){
		x[0] /= length;
		x[1] /= length;
		x[2] /= length;
	}else
		throw std::runtime_error("Can't normalize a null 3D vector!!!");
	
}

Vec3D Vec3D::mini(const Vec3D& vec_1,const Vec3D& vec_2)
{
  return Vec3D(min_(vec_1[0],vec_2[0]),min_(vec_1[1],vec_2[1]),min_(vec_1[2],vec_2[2])) ;
}

Vec3D Vec3D::maxi(const Vec3D& vec_1,const Vec3D& vec_2)
{
  return Vec3D(max_(vec_1[0],vec_2[0]),max_(vec_1[1],vec_2[1]),max_(vec_1[2],vec_2[2])) ;
}

void Vec3D::out(const char *vname) const
{
	(void)printf("%s (%f, %f, %f)\n",vname,x[0],x[1],x[2]);
}

} // end namespace FemViewer
