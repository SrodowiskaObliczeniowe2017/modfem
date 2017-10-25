#ifndef _VECTOR_3D_H_
#define _VECTOR_3D_H_

#include <cassert>
#include <cmath>
#include "../utils/fv_assert.h"
#include "Point3D.h"
namespace FemViewer {

	class Vec3D {
private:
	float x[3];
public:
	Vec3D() { x[0] = x[1] = x[2] = 0.0f; }
	~Vec3D(){}
	Vec3D(float v_) { x[0] = x[1] = x[2] = v_; } 
	Vec3D(float x_, float y_, float z_ = 0.0f) {
		x[0] = x_; x[1] = y_, x[2] = z_; }
	Vec3D(const float *x_) {
		x[0] = x_[0]; x[1] = x_[1]; x[2] = x_[2]; } 
	Vec3D(const Vec3D& vec_){ 
		x[0] = vec_.x[0]; x[1] = vec_.x[1]; x[2] = vec_.x[2]; }
	template<class T>
	Vec3D(const Point3D<T>& p)
	{
		x[0] = (float)p.x;
		x[1] = (float)p.y;
		x[2] = (float)p.z;
	}
	inline float operator [](int idx) const {
		assert((idx>=0) && (idx<3));
		return x[idx];
	}

	inline float& operator [](int idx){
		assert((idx>=0) && (idx<3));
		return x[idx];
	}
	
	inline const float* coord()const { return &x[0]; }
	inline void Set(float _x, float _y, float _z)
	{ x[0] = _x; x[1] = _y; x[2] = _z; }
	inline float _x()const { return x[0]; } 
	inline float _y()const { return x[1]; } 
	inline float _z()const { return x[2]; }
	inline bool IsZero() const;

	Vec3D& operator  =(const Vec3D& vec_);
	Vec3D& operator +=(const Vec3D& vec_);
	Vec3D& operator -=(const Vec3D& vec_);
	Vec3D& operator *=(const float v_);
	Vec3D& operator /=(const float v_);
	Vec3D& operator  +(const float v_);
	Vec3D& operator  -(const float v_);
	Vec3D& operator %=(const Vec3D& vec_);

	Vec3D getOXY() const;
	Vec3D getOYZ() const;
	Vec3D getOZX() const;

	Vec3D operator +(const Vec3D& vec_) const;
	Vec3D operator -(const Vec3D& vec_) const;
	Vec3D operator -() const;
	float getDistance(const Vec3D& vec) const;
	float getAngle(const Vec3D& vec) const;
		
	friend Vec3D operator *(const Vec3D& vec_, const float rhv_);
	friend Vec3D operator *(const float lhv_, const Vec3D& vec_);
	friend Vec3D operator /(const Vec3D& vec_, const float rhv_);
	friend Vec3D operator %(const Vec3D& vec_1, const Vec3D& vec_2);
	friend float dotProd(const Vec3D& vec_1, const Vec3D& vec_2);
	friend Vec3D operator -(const Vec3D& vec_);
	friend bool operator ==(const Vec3D& vec_1, const Vec3D& vec_2);
	friend bool operator !=(const Vec3D& vec_1, const Vec3D& vec_2);
	friend float getAngle(const Vec3D& vec_1, const Vec3D& vec_2);

	float quadNorm() const;
	float Norm() const;
	void normalize();

	static Vec3D mini(const Vec3D& vec_1, const Vec3D& vec_2);
	static Vec3D maxi(const Vec3D& vec_1, const Vec3D& vec_2);
	void out(const char* vname = "") const;
	
	
};

inline bool Vec3D::IsZero() const
{
	return( this->Norm() == 0.0f);
		    
}
inline Vec3D& Vec3D::operator =(const Vec3D& vec_)
{
	x[0] = vec_.x[0];
	x[1] = vec_.x[1];
	x[2] = vec_.x[2];
	return(*this);
}

inline Vec3D& Vec3D::operator +=(const Vec3D& vec_)
{
	x[0] += vec_.x[0]; x[1] += vec_.x[1]; x[2] += vec_.x[2];
	return (*this);
}

inline Vec3D& Vec3D::operator -=(const Vec3D& vec_)
{
	x[0] -= vec_.x[0]; x[1] -= vec_.x[1]; x[2] -= vec_.x[2];
	return (*this);
}

inline Vec3D& Vec3D::operator *=(const float v_)
{
	x[0] *= v_;
	x[1] *= v_;
	x[2] *= v_; 
	return(*this);
}

inline Vec3D& Vec3D::operator /=(const float v_)
{
	x[0] /= v_;
	x[1] /= v_;
	x[2] /= v_; 
	return(*this);
}

inline Vec3D& Vec3D::operator +(const float v_)
{
	x[0] += v_;
	x[1] += v_;
	x[2] += v_;
	return (*this);
}

inline Vec3D& Vec3D::operator -(const float v_)
{
	return this->operator+(-v_);
}

inline Vec3D& Vec3D::operator %=(const Vec3D& vec_)
{
	float tmp_xyz[3];
	tmp_xyz[0] = x[0]; tmp_xyz[1] = x[1]; tmp_xyz[2] = x[2];
	x[0] = tmp_xyz[1]*vec_.x[2] - tmp_xyz[2]*vec_.x[1];
	x[1] = tmp_xyz[2]*vec_.x[0] - tmp_xyz[0]*vec_.x[2];
	x[2] = tmp_xyz[0]*vec_.x[1] - tmp_xyz[1]*vec_.x[0];
	return(*this);
}

inline Vec3D Vec3D::getOXY() const
{
	return Vec3D(x[0],x[1],0.0f);
}

inline Vec3D Vec3D::getOYZ() const
{
	return Vec3D(0.0f,x[1],x[2]);
}

inline Vec3D Vec3D::getOZX() const
{
	return Vec3D(x[0],0.0f,x[2]);
}

inline Vec3D Vec3D::operator +(const Vec3D& vec_) const
{
	return Vec3D(x[0]+vec_.x[0],x[1]+vec_.x[1],x[2]+vec_.x[2]);
}

inline Vec3D Vec3D::operator -(const Vec3D& vec_) const
{
	return Vec3D(x[0]-vec_.x[0],x[1]-vec_.x[1],x[2]-vec_.x[2]);
}

inline Vec3D Vec3D::operator -() const
{
	return Vec3D(-x[0], -x[1], -x[2]);
}


inline float Vec3D::getDistance(const Vec3D& vec) const
{
	return (*this - vec).Norm();
}

inline float Vec3D::getAngle(const Vec3D& vec) const
{
	//FV_ASSERT(Norm() > 0.0f);
	//FV_ASSERT(vec.Norm() > 0.0f);
	float dot = dotProd(*this,vec);
	return acos(dot/(Norm()*vec.Norm()));

}


inline Vec3D operator *(const Vec3D& vec_, const float rhv_)
{ 
	return Vec3D(vec_.x[0]*rhv_,vec_.x[1]*rhv_,vec_.x[2]*rhv_);
}

inline Vec3D operator *(const float lhv_, const Vec3D& vec_)
{ 
	return Vec3D(vec_.x[0]*lhv_,vec_.x[1]*lhv_,vec_.x[2]*lhv_);
}

inline Vec3D operator /(const Vec3D& vec_, const float rhv_)
{
	return Vec3D(vec_.x[0]/rhv_,vec_.x[1]/rhv_,vec_.x[2]/rhv_);
}

inline Vec3D operator %(const Vec3D& vec_1, const Vec3D& vec_2)
{
	//printf("rezultat: %f\n",vec_1.x[2]*vec_2.x[0]-vec_2.x[0]*vec_2.x[2]);
	return Vec3D(vec_1.x[1]*vec_2.x[2]-vec_1.x[2]*vec_2.x[1],
		         vec_1.x[2]*vec_2.x[0]-vec_1.x[0]*vec_2.x[2],
				 vec_1.x[0]*vec_2.x[1]-vec_1.x[1]*vec_2.x[0]);
}

inline float dotProd(const Vec3D &vec_1, const Vec3D &vec_2)
{
	return(vec_1.x[0]*vec_2.x[0] + vec_1.x[1]*vec_2.x[1] + vec_1.x[2]*vec_2.x[2]);
}

//inline float getAngle(const Vec3D& vec_1, const Vec3D& vec_2)
//{
//	FV_ASSERT(vec_1.Norm() > 0.0f);
//	FV_ASSERT(vec_2.Norm() > 0.0f);
//	float dot = dotProd(vec_1,vec_2);
//	return acos(dot/(vec_1.Norm()*vec_2.Norm()));
//
//}

} // end namespace FemViewer

#endif /* _VECTOR_3D_H_ */
