/*
 * point_T.hh
 *
 *  Created on: 12-01-2012
 *      Author: Paweł Macioł
 */

#ifndef POINT_T_hh__
#define POINT_T_hh__

#include"vector_T.hh"

template<typename T>
class CPoint3 : public CVecT<T,3>
{
public:
	CPoint3() {}

	CPoint3(T x_, T y_, T z_)
    {
        (*this)(0) = x_;
        (*this)(1) = y_;
        (*this)(2) = z_;
    }

	CPoint3(const CVecT<T,3>& v_) : CVecT<T,3><T,3>(v_)
    {}

    const T& x() const { return this->at(0); }
    const T& y() const { return this->at(1); }
    const T& z() const { return this->at(2); }
    T& x() { return (*this)(0); }
    T& y() { return (*this)(1); }
    T& z() { return (*this)(2); }
};

typedef CPoint3<float > Point3f;
typedef CPoint3<double> Point3d;

#endif /* POINT_T_hh__ */
