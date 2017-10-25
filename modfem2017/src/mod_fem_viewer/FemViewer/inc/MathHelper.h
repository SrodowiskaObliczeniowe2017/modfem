#ifndef _MATH_HELPER_H_
#define _MATH_HELPER_H_

#include <limits>
#include <cmath>
#include <ostream>
#include <map>
#include <vector>
#include "fv_compiler.h"
#include <stdint.h>
#if defined(WIN32) && defined(_MSC_VER)
# pragma warning( disable : 4201 )
#endif


namespace FemViewer {
namespace fvmath {

#define AVEC3(T)    union { struct { T x,y,z; }; T v[3]; }
#define AVEC4(T)    union { struct { T x,y,z,w; }; T v[4]; }
#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

static const float kFloatMaximumValue = std::numeric_limits<float>::max();
static const float kInfinity = std::numeric_limits<float>::infinity();

template<typename T>
struct Vec3 {
	AVEC3(T);
};

template<typename T>
struct Vec4 {
	AVEC4(T);
};

typedef Vec4<float> Vec4f;


template<typename T>
class Ray3 {
	Vec3<T> org, dir;
};



template<typename T> class CVec3;

template<typename T>
std::ostream& operator << (std::ostream& os, const CVec3<T>& rhs);

template<typename T>
class CVec3 : public Vec3<T> {
public:
	CVec3(T x_= T(0), T y_=T(0), T z_=T(0)) //: Vec3<T>::x(x_), Vec3<T>::y(y_), Vec3<T>::z(z_)
	{ Vec3<T>::x = x_; Vec3<T>::y = y_; Vec3<T>::z = z_; }
	template<typename U> CVec3(const U v_[]);
	template<typename U> CVec3(U x_, U y_, U z_);
	template<typename U> CVec3(const Vec3<U>& v_);
	template<typename U = T> CVec3<T>& operator=(const Vec3<U>& v_);
	template<typename U>
	void operator+=(const fvmath::Vec3<U>& v_);
	template<typename U>
	operator U*() { return static_cast<U *>(this->v); }
	T operator[](int idx) const { return this->v[idx]; }
	T& operator[](int idx) { return this->v[idx]; }
	template<typename U>
	void operator +=(const U* rhs);
	void operator *=(const T& alpha);
	void operator -=(const Vec3<T>& rhs);
	void operator /=(const Vec3<T>& rhs);

	friend std::ostream& operator << <>(std::ostream& os, const CVec3& rhs);
};
template<typename T> template<typename U>
inline CVec3<T>::CVec3(const U v_[]) {
	for(int i(0);i<3;i++) Vec3<T>::v[i] = static_cast<T>(v_[i]);
}

template<> template<>
inline CVec3<float>::CVec3(const float v_[]) {
	for(int i(0);i<3;i++) Vec3<float>::v[i] = v_[i];
}

template<typename T> template<typename U>
inline CVec3<T>::CVec3(U x_,U y_, U z_)
{
	Vec3<T>::x = static_cast<T>(x_);
	Vec3<T>::y = static_cast<T>(y_);
	Vec3<T>::z = static_cast<T>(z_);
}



template<typename T> template<typename U>
inline CVec3<T>::CVec3(const Vec3<U>& v_)
{
	for(int i(0);i<3;i++) Vec3<T>::v[i] = v_.v[i];
}

template<typename T> template<typename U>
inline CVec3<T>& CVec3<T>::operator=(const Vec3<U>& v_)
{
	Vec3<T>::x = static_cast<T>(v_.x);
	Vec3<T>::y = static_cast<T>(v_.y);
	Vec3<T>::z = static_cast<T>(v_.z);

	return *this;
}


template<> template<>
inline CVec3<float>& CVec3<float>::operator=(const Vec3<float>& v) {
	for (unsigned i(0);i<3;++i) this->v[i] = v.v[i];
	return *this;
}

template<typename T> template<typename U>
inline void CVec3<T>::operator+=(const Vec3<U>& v_)
{
	Vec3<T>::x += static_cast<T>(v_.x);
	Vec3<T>::y += static_cast<T>(v_.y);
	Vec3<T>::z += static_cast<T>(v_.z);
}

template<> template<>
inline void CVec3<float>::operator+=(const Vec3<float>& v_)
{
	Vec3<float>::x += v_.x;
	Vec3<float>::y += v_.y;
	Vec3<float>::z += v_.z;
}


template<typename T> template<typename U>
inline void CVec3<T>::operator+=(const U* rhs)
{
	for (int i(0);i<3;++i) Vec3<T>::v[i] = static_cast<T>(rhs[i]);
}

template<typename T>
inline void CVec3<T>::operator*=(const T& alpha)
{
	for (int i(0);i<3;++i) Vec3<T>::v[i] *= alpha;
}

template<typename T>
inline void CVec3<T>::operator-=(const Vec3<T>& rhs)
{
	for (int i(0);i<3;++i) Vec3<T>::v[i] -= rhs.v[i];
}

template<typename T>
inline void CVec3<T>::operator/=(const Vec3<T>& rhs)
{
	for (int i(0);i<3;++i) Vec3<T>::v[i] /= rhs.v[i];
}

template<typename T>
std::ostream& operator <<(std::ostream& os, const CVec3<T>& rhs)
{
	return os << "x = " << rhs.x << " y = " << rhs.y << " z = " << rhs.z;
}

typedef Vec3<int> Vec3i;
typedef Vec3<uint32_t> Vec3u;
typedef Vec3<float> Vec3f;
typedef Vec3<double> Vec3d;
typedef Vec3d Point3d;
typedef CVec3<int> CVec3i;
typedef CVec3<uint32_t> CVec3u;
typedef CVec3<float> CVec3f;
typedef CVec3<double> CVec3d;

template<typename T> inline T Norm(const Vec3<T>& v)
{ return sqrt(v.x*v.x+v.y*v.y+v.z*v.z); }

template<typename T> T Normalize(Vec3<T>& v) {
	T len = Norm<T>(v);
	if (len == 0) len = 1;
	T inv = 1 / len;
		//throw("Can't normalize null vector!");
	v.x *= inv;
	v.y *= inv;
	v.z *= inv;

	return(len);
}

template<typename T>
inline T Normalize(Vec3<T> *v) {
	return Normalize(*v);
}

template<typename T> inline CVec3<T> operator+(const Vec3<T>& v, const Vec3<T>& u )
{ return CVec3<T>( v.x+u.x, v.y+u.y, v.z+u.z ); }

template<typename T> inline CVec3<T> operator-( const Vec3<T>& v, const Vec3<T>& u )
{
	Vec3<T> res;
	for (int i(0);i<3;++i) res.v[i] = v.v[i] - u.v[i];
	return res;
}

/*
template<typename T> inline CVec3<T> operator-( const Vec3<T> *v, const Vec3<T> *u )
{
	return CVec3<T>(*v - *u);
}
*/

template<typename T> T Dot(const Vec3<T>& v, const Vec3<T>& u)
{ return (v.x*u.x+v.y*u.y+v.z*u.z); }

template<typename T> inline CVec3<T> operator*( const Vec3<T>& v, const Vec3<T>& u)
{ return CVec3<T>( v.y*u.z - v.z*u.y,
                   v.z*u.x - v.x*u.z,
                   v.x*u.y - v.y*u.x );
}

template<typename T> inline void Cross( const Vec3<T>* v, const Vec3<T>* u, Vec3<T>* out)
{
	out->x = (v->y*u->z - v->z*u->y);
	out->y = (v->z*u->x - v->x*u->z);
    out->z = (v->x*u->y - v->y*u->x);
    Normalize(out);
}

template<typename T> inline CVec3<T> operator *(T alfa, const Vec3<T>& v)
{ return CVec3<T>(alfa*v.x,alfa*v.y,alfa*v.z); }

template<typename T> inline CVec3<T> operator *(const Vec3<T>& v,T alfa) {
	return CVec3<T>(alfa*v.x,alfa*v.y,alfa*v.z);
}

template<typename T> inline CVec3<T> operator / (const T& l,const Vec3<T>& rhs){
	CVec3<T> w;
	for (int i = 0; i<3; ++i) w.v[i] = l / rhs.v[i];
	return w;
}

template<typename T> inline Vec3<T> operator / (const Vec3<T>& lhs,const Vec3<T>& rhs) {
	Vec3<T> res;
	for (int i = 0; i<3; ++i) res.v[i] = lhs.v[i] / rhs.v[i];
	return res;
}

template<typename T> inline bool operator ==(const CVec3<T>& lhs,const CVec3<T>& rhs) {
	for (unsigned int i(0); i < 3; ++i) {
		if (lhs[i] != rhs[i]) return false;
	}
	return true;
}

template<typename T>
inline CVec3<T> MulM3x3V(const T mat3[9], const Vec3<T>& v) {
	CVec3<T> row0(mat3[0],mat3[3],mat3[6]);
	CVec3<T> row1(mat3[1],mat3[4],mat3[7]);
	CVec3<T> row2(mat3[2],mat3[5],mat3[8]);

	return CVec3<T>(Dot(row0,v),Dot(row1,v),Dot(row2,v));
}


#undef AVEC3
#undef AVEC4

template<typename T>
inline bool ValueWithin(const T& v, const T& v1, const T& v2)
{ return (v >= v1 && v <= v2) || (v >= v2 && v <= v1); }

template<typename T> bool Compare(const T& arg1,const T& arg2);

template<>
bool Compare<float>(const float& arg1,const float& arg2);
template<>
bool Compare<double>(const double& arg1,const double& arg2);
template<>
bool Compare<CVec3f>(const CVec3f& v1,const CVec3f& v2);

template<typename T> bool Less(const T& arg1, const T& arg2);
template<> bool Less(const double& arg1,const double& arg2);

template<typename T>
inline T max(T arg1, T arg2) {
	return (arg1 > arg2) ? arg1 : arg2;
}

template<typename T>
inline T min(T arg1, T arg2) {
	return (arg1 < arg2) ? arg1 : arg2;
}

template<typename T>
T DotProd(const T v1[], const T v2[]);



extern const float  epsil;
extern const double ref;

template<typename T> inline const T abs(const T& arg)
{ return (arg < T(0)) ? -arg : arg; }

template<typename T> inline T clamp(const T &v, const T &lo, const T &hi)
{ return max(lo, min(v, hi)); }

// 32-bit calculations

/* Find integer log base 2 of a 32-bit IEEE float */
/* Conditions: f > 0.0 && finite(f) && isnormal(f) */
inline int log2(const float f)
{
	int c;
	c = *(unsigned int*) &f;
	return ((c >> 23) - 127);
}



/* Round up to the next highest power of 2 of 32-bit IEE float*/
inline int npow2(const float f)
{
	unsigned int r;
	if (f > 1.0f) {
		const unsigned int t = 1U << ((*(unsigned int *)&f >> 23) - 0x7f);
		r = t << (unsigned int)((t < (unsigned int)f));
	} else {
		r = 1;
	}

	return(r);
}

inline int npow2(unsigned int v)
{
	return npow2((float)v);
}



// CCW-direction p1,p2,p3
template<typename T,typename U>
CVec3<T> GetNormal2Plane(const U* A,const U* B,const U* C)
{
	CVec3<T> vAC( C[0] - A[0], C[1] - A[1], C[2] - A[2]);
	CVec3<T> vAB( B[0] - A[0], B[1] - A[1], B[2] - A[2]);
	CVec3<T> nv((vAC * vAB).v);
	Normalize(nv);
	return nv;
}

/* Routines for ray triangle test */
/* Muller algorithm */
template<typename T>
static int testIntersection(const Ray3<T>& ray,
		                    const CVec3<T>& A, const CVec3<T>& B, const CVec3<T>& C,
					        T& t, T& u, T& v,
					        bool singleSide = true)
{
	// Make triangle
	CVec3<T> AB = B - A;
	CVec3<T> AC = C - A;
	CVec3<T> pvec = ray.dir * AC;
	T det = Dot(AB, pvec);
	// test
	if (singleSide) {
		if (det < epsil)
			return 0;
		T tvec = ray.orig - A;
		u = Dot(tvec, pvec);
		if (u < 0 || u > det)
			return 0;
		CVec3<T> qvec = tvec - AB;
		v = Dot(ray.dir, qvec);
		if (v < 0 || u + v > det)
			return 0;
		t = Dot(AC, qvec);
		T inv_det = 1 / det;
		t *= inv_det;
		v *= inv_det;
		u *= inv_det;
	} else {
		if (det < -epsil && det > epsil)
			return 0;
		T inv_det = 1 / det;
		T tvec = ray.orig - A;
		u = Dot(tvec, pvec) * inv_det;
		if (u < 0 || u > 1)
			return 0;
		CVec3<T> qvec = tvec - AB;
		v = Dot(ray.dir, qvec) * inv_det;
		if (v < 0 || u + v > 1)
			return 0;
		t = t = Dot(AC, qvec) * inv_det;
	}
	return 1;
}

inline int roundup(const int value,const int div)
{
	int ret = value % div ? (value / div + 1)*div : value;
	return ret;
}

template <typename T>
inline bool is_near(T v1, T v2, T margin) {
	return fabs(v1-v2) < margin;
}

// Convert value to rgb color
template<typename TValue>
CVec3<TValue> convertToRGB(TValue f)
{
	const TValue dx = 0.8;
	f = (f < 0) ? 0 : (f > 1) ? 1 : f;
	TValue g = (6 - 2*dx)*f + dx;
	TValue R = max(0, (3 -fabs(g-4)-fabs(g-5))/2);
	TValue G = max(0, (4 -fabs(g-2)-fabs(g-4))/2);
	TValue B = max(0, (3 -fabs(g-1)-fabs(g-2))/2);

	return CVec3<TValue>(R, G, B);

}

typedef unsigned int indextype;

template< typename TVertex >
struct ComparatorSpec
{
	bool operator ()(const TVertex& lh,const TVertex& rh)
	{
		return lh.info < rh.info;
	}
};


template< typename TIndex,
          typename TVertex,
		  typename TComp = ComparatorSpec<TVertex>>
bool getSimilarVertexIndex(
	TVertex & packed,
	std::map<TVertex,TIndex,TComp> & VertexToOutIndex,
	TIndex & result)
{
	typename std::map<TVertex,TIndex,TComp>::iterator it = VertexToOutIndex.find(packed);
	if ( it == VertexToOutIndex.end() ){
		return false;
	}else{
		result = it->second;
		return true;
	}
}



template< typename TIndex,
          typename TVertex,
		  typename TComp = ComparatorSpec<TVertex>>
void removeDuplicates(std::vector< TVertex >& inout_vertices,
			  std::vector< TIndex >& inout_indices)
{
	// Temporary container
	std::vector<TVertex> out_vertices;
	out_vertices.reserve(inout_vertices.size());

	assert(inout_indices.empty()==false);

	TIndex index = 0;
	TIndex num = inout_indices.size();
	typename std::map<TVertex,TIndex,TComp>  VertexToOutIndex;

	// For each input index
	for (TIndex i=0; i < num; i++)
	{
		// Get index value
		TIndex vt_index = inout_indices[i];
		// Get vertex from input vector
		TVertex v = TVertex(inout_vertices[ vt_index ]);

		// Check whether vertex is present already
		TIndex ii;
		bool found = getSimilarVertexIndex(v,VertexToOutIndex,ii);

		if (found) {
			inout_indices[i] = ii;
		}
		else {
			//mfp_debug("Size of vertices: %u",out_vertices.size());
			out_vertices.push_back(v);
			TIndex newindex = out_vertices.size() - 1;
			inout_indices[i] = newindex;
			VertexToOutIndex[v] = newindex;
		}
	}
	// Resize
	//std::vector<Vertex>(out_vertices).swap(inout_vertices);
	inout_vertices = out_vertices;
}

		
} // end namespace fvmath
} // end namespace FemViewer

#endif /* _MATH_HELPER_H_
*/
