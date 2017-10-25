#ifndef _MATRIX_H_
#define _MATRIX_H_

#include<vector>
#include<iostream>
#include"fv_inc.h"
#include"Vec3D.h"
#include"Homogenous.h"
#include"MathHelper.h"

namespace FemViewer {
using namespace fvmath;

template<typename T> class Matrix;

template<typename T>
std::ostream& operator << (std::ostream& strm,const Matrix<T>& mt);

template<typename T = float>
class Matrix {
public:
	static Matrix Identity;
private:
	bool isIdentity;
public:

		//
		// T = (matrix[0] matrix[4] matrix[ 8] matrix[12]
		//      matrix[1] matrix[5] matrix[ 9] matrix[13]
		//		matrix[2] matrix[6] matrix[10] matrix[14]
		//		matrix[3] matrix[7] matrix[11] matrix[15])
		//
		std::vector<T> matrix;
		
		// default constructor
		Matrix()
		{
			SetIdentity();
		}

		Matrix(T v) : matrix(16,v)
		{ }

		// copy constructor
		Matrix(const Matrix<T> &Tr) : matrix(Tr.matrix)
		{ }

		// assingment operator
		Matrix& operator = (const Matrix<T> &Tr)
		{
			if (this != &Tr)
				matrix = Tr.matrix;
			return *this;
		}

		bool IsIdentity() const { return isIdentity; }
		void SetIdentity()
		{
			matrix.resize(16, 0);
			matrix[0] = matrix[5] = matrix[10] = matrix[15] = T(1);
			isIdentity = true;
		}

			  T& operator [](size_t idx)       { return matrix[idx]; }
		const T& operator [](size_t idx) const { return matrix[idx]; }

		// add
		Matrix operator + (const Matrix<T> &Tr) const
		{
			Matrix result(*this);
			for(size_t i = 0; i < 16; i++)
				result.matrix[i] += Tr.matrix[i];
			return result;
		}

		// sub
		Matrix operator - (const Matrix<T> &Tr) const
		{
			Matrix result(*this);
			for(size_t i = 0; i < 16; i++)
				result.matrix[i] -= Tr.matrix[i];
			return result;
		}

		// multiplicate
		Matrix operator * (const Matrix<T>& Tr) const
		{
			Matrix result;
			result.matrix[0] = matrix[0]*Tr.matrix[0] + matrix[4]*Tr.matrix[1] + matrix[8]*Tr.matrix[2] + matrix[12]*Tr.matrix[3];
			result.matrix[1] = matrix[1]*Tr.matrix[0] + matrix[5]*Tr.matrix[1] + matrix[9]*Tr.matrix[2] + matrix[13]*Tr.matrix[3];
			result.matrix[2] = matrix[2]*Tr.matrix[0] + matrix[6]*Tr.matrix[1] + matrix[10]*Tr.matrix[2] + matrix[14]*Tr.matrix[3];
			result.matrix[3] = matrix[3]*Tr.matrix[0] + matrix[7]*Tr.matrix[1] + matrix[11]*Tr.matrix[2] + matrix[15]*Tr.matrix[3];

			result.matrix[4] = matrix[0]*Tr.matrix[4] + matrix[4]*Tr.matrix[5] + matrix[8]*Tr.matrix[6] + matrix[12]*Tr.matrix[7];
			result.matrix[5] = matrix[1]*Tr.matrix[4] + matrix[5]*Tr.matrix[5] + matrix[9]*Tr.matrix[6] + matrix[13]*Tr.matrix[7];
			result.matrix[6] = matrix[2]*Tr.matrix[4] + matrix[6]*Tr.matrix[5] + matrix[10]*Tr.matrix[6] + matrix[14]*Tr.matrix[7];
			result.matrix[7] = matrix[3]*Tr.matrix[4] + matrix[7]*Tr.matrix[5] + matrix[11]*Tr.matrix[6] + matrix[15]*Tr.matrix[7];

			result.matrix[8] = matrix[0]*Tr.matrix[8] + matrix[4]*Tr.matrix[9] + matrix[8]*Tr.matrix[10] + matrix[12]*Tr.matrix[11];
			result.matrix[9] = matrix[1]*Tr.matrix[8] + matrix[5]*Tr.matrix[9] + matrix[9]*Tr.matrix[10] + matrix[13]*Tr.matrix[11];
			result.matrix[10] = matrix[2]*Tr.matrix[8] + matrix[6]*Tr.matrix[9] + matrix[10]*Tr.matrix[10] + matrix[14]*Tr.matrix[11];
			result.matrix[11] = matrix[3]*Tr.matrix[8] + matrix[7]*Tr.matrix[9] + matrix[11]*Tr.matrix[10] + matrix[15]*Tr.matrix[11];

			result.matrix[12] = matrix[0]*Tr.matrix[12] + matrix[4]*Tr.matrix[13] + matrix[8]*Tr.matrix[14] + matrix[12]*Tr.matrix[15];
			result.matrix[13] = matrix[1]*Tr.matrix[12] + matrix[5]*Tr.matrix[13] + matrix[9]*Tr.matrix[14] + matrix[13]*Tr.matrix[15];
			result.matrix[14] = matrix[2]*Tr.matrix[12] + matrix[6]*Tr.matrix[13] + matrix[10]*Tr.matrix[14] + matrix[14]*Tr.matrix[15];
			result.matrix[15] = matrix[3]*Tr.matrix[12] + matrix[7]*Tr.matrix[13] + matrix[11]*Tr.matrix[14] + matrix[15]*Tr.matrix[15];

			return result;
		}

		// add to this
		Matrix& operator += (const Matrix<T> &Tr)
		{
			for(size_t i = 0; i < 16; i++)
				matrix[i] += Tr.matrix[i];
			return *this;
		}

		// sub from this
		Matrix& operator -= (const Matrix<T> &Tr)
		{
			for(size_t i = 0; i < 16; i++)
				matrix[i] -= Tr.matrix[i];
			return *this;
		}

		// multiplicate this by a constant
		Matrix& operator *= (T constant)
		{
			for(size_t i = 0; i < 16; i++)
				matrix[i] *= constant;
			return *this;
		}

		// divide this by a constant
		Matrix& operator /= (T constant)
		{
			for(size_t i = 0; i < 16; i++)
				matrix[i] /= constant;
			return *this;
		}

		// multiplicate by a D�scartes coordinate, i.e., transform a D�scartes coordinate into another one
		Vec3D operator * (const Vec3D& dc) const
		{
			float w = matrix[3] * dc[0] + matrix[7] * dc[1] + matrix[11] * dc[2] + matrix[15];
			// note, if w = 0, then transformation matrix is singular
			// so, be careful :)
			return Vec3D(
				(matrix[0] * dc[0] + matrix[4] * dc[1] + matrix[ 8] * dc[2] + matrix[12]) / w,
				(matrix[1] * dc[0] + matrix[5] * dc[1] + matrix[ 9] * dc[2] + matrix[13]) / w,
				(matrix[2] * dc[0] + matrix[6] * dc[1] + matrix[10] * dc[2] + matrix[14]) / w);
		}

		CVec3<T> operator * (const Vec3<T>& dc) const
		{
			T w = this->matrix[3] * dc.x + this->matrix[7] * dc.y + this->matrix[11] * dc.z + this->matrix[15];
			T inv_w = 1 / w;
			// note, if w = 0, then transformation matrix is singular
			// so, be careful :)
			return CVec3<T>(
				(this->matrix[0] * dc.x + this->matrix[4] * dc.y + this->matrix[ 8] * dc.z + this->matrix[12]) * inv_w,
				(this->matrix[1] * dc.x + this->matrix[5] * dc.y + this->matrix[ 9] * dc.z + this->matrix[13]) * inv_w,
				(this->matrix[2] * dc.x + this->matrix[6] * dc.y + this->matrix[10] * dc.z + this->matrix[14]) * inv_w);
		}

		// multiplicate by a homogeneous coordinate, i.e., transform a homogeneous coordinate into another one
		Homogenous operator * (const Homogenous &hc) const
		{
			float w = matrix[3] * hc.x + matrix[7] * hc.y + matrix[11] * hc.z + matrix[15];
			// note, if w = 0, then transformation matrix is singular
			// so, be careful :)
			return Homogenous(
				(matrix[0] * hc.x + matrix[4] * hc.y + matrix[ 8] * hc.z + matrix[12]) / w,
				(matrix[1] * hc.x + matrix[5] * hc.y + matrix[ 9] * hc.z + matrix[13]) / w,
				(matrix[2] * hc.x + matrix[6] * hc.y + matrix[10] * hc.z + matrix[14]) / w);
		}


		void make3x3(T* out) const
		{
			int i(0);
			do {
				if (i == 3 || i == 7) i++;
				*out++ = matrix[i++];
			} while (i<11);
		}

		void LoadIdentity()
		{
			for (size_t i = 0; i < 16; i++)
				matrix[i] = (i % 5) ? 0.0f : 1.0f;
		}

		void LoadNullMatrix()
		{
			for (size_t i = 0; i < 16; i++)
				matrix[i] = 0.0f;
		}

		void SetScalefactors(const T& sx, const T&sy, const T& sz)
		{
			matrix[0] = sx;
			matrix[5] = sy;
			matrix[10] = sz;
		}

		void SetTransformationFactors(const Vec3<T>& coef)
		{
			matrix[12] = coef.x;
			matrix[13] = coef.y;
			matrix[14] = coef.z;
		}

		// Is euqal to 0 the determinant of the transformation matrix?
		bool IsSingular() const;

		// perform the transformation
		void Apply() const
		{
			glMultMatrixf(&matrix[0]);
		}

		// virtual destructor
		virtual ~Matrix()
		{
			matrix.clear();
		}

		// rotate about 
		void rotateAbout(const float *coord, const float angle);

		// scale
		void scale(const Vec3D& vec);

		// translate
		void translate(const Vec3D& vec);

		// build
		Matrix Invers() const;

		static Matrix LookAt(const fvmath::Vec3<T>& eye,
				      const fvmath::Vec3<T>& center,
				      const fvmath::Vec3<T>& up)
		{
			fvmath::CVec3<T> d = center - eye;
			Normalize(d);
			fvmath::CVec3<T> u(up);
			Normalize(u);
			fvmath::CVec3<T> s = d * u;
			u = s * d;
			// Not needed because of earlier
			//Normalize(u);

			Matrix<T> Result;
			Result[0] = s.x;
			Result[4] = s.y;
		    Result[8] = s.z;
			Result[1] = u.x;
			Result[5] = u.y;
			Result[9] = u.z;
			Result[2] =-d.x;
			Result[6] =-d.y;
			Result[10] =-d.z;
			Result[12] =-Dot(s, eye);
			Result[13] =-Dot(u, eye);
			Result[14] = Dot(d, eye);

			return Result;
		}

		// dump data to screen
		friend std::ostream& operator << <>(std::ostream& strm,const Matrix& mt);
	};

   	template<typename T>
 	Matrix<T> Matrix<T>::Identity;

	// multiplicate the transformation matrix by a constant from left
	template<typename T>
	inline Matrix<T> operator * (T lhs, const Matrix<T> &rhs)
	{

		Matrix<T> result(rhs);
		for(typename std::vector<T>::iterator i = result.matrix.begin(); i < result.matrix.end(); i++)
			*i *= lhs;
		return result;
	}

	// multiplicate the transformation matrix by a constant from right
	template<typename T>
	inline Matrix<T> operator * (const Matrix<T> &lhs, T rhs)
	{
		Matrix<T> result(lhs);
		for(typename std::vector<T>::iterator i = result.matrix.begin(); i < result.matrix.end(); i++)
			*i *= rhs;
		return result;
	}

//	// output to stream
//	template<typename T>
//	inline std::ostream& operator << (std::ostream& lhs, const Matrix<T> &transformation)
//	{
//		lhs << transformation.matrix[0];
//		for(size_t i = 1; i < 16; i++)
//			lhs << " " << transformation.matrix[i];
//		return lhs;
//	}

	// input from stream
	template<typename T>
	inline std::istream& operator >> (std::istream& lhs, Matrix<T>& transformation)
	{
		lhs >> transformation.matrix[0];
		for(size_t i = 1; i < 16; i++)
			lhs >> transformation.matrix[i];
		return lhs;
	}

	template<typename T>
	std::ostream & operator << (std::ostream& strm,const Matrix<T>& mt)
	{
		strm << "matrix: \n"
		<< "m[0] = "<< mt.matrix[0] << " m[4] = " << mt.matrix[4] << " m[ 8] = " << mt.matrix[ 8] << " m[12] = " << mt.matrix[12] << std::endl
		<< "m[1] = "<< mt.matrix[1] << " m[5] = " << mt.matrix[5] << " m[ 9] = " << mt.matrix[ 9] << " m[13] = " << mt.matrix[13] << std::endl
		<< "m[2] = "<< mt.matrix[2] << " m[6] = " << mt.matrix[6] << " m[10] = " << mt.matrix[10] << " m[14] = " << mt.matrix[14] << std::endl
		<< "m[3] = "<< mt.matrix[3] << " m[7] = " << mt.matrix[7] << " m[11] = " << mt.matrix[11] << " m[15] = " << mt.matrix[15] << std::endl;
		return strm;
	}


	typedef Matrix<float> Matrixf;
	typedef Matrix<double> Matrixd;

} // end namespace FemViewer
#endif 
