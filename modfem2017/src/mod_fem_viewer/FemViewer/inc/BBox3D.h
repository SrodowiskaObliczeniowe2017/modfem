#ifndef _BOUNDING_BOX_H_
#define _BOUNDING_BOX_H_

#define NOMINMAX 1
#include "defs.h"
#include "Matrix.h"
#include "Vec3D.h"
#include "MathHelper.h"
#include "Ray.h"

#include <string>
#include <ostream>
#include <limits>

namespace FemViewer {

// Forward declarations
template<typename T> class Ray;
class BBox3D;
BBox3D operator*(const Matrix<float>& pMatrix,
                 const BBox3D& pBBox3D);

template<typename T>
struct AAbb {
	union {
		struct { CVec3<T> mn, mx; };
		CVec3<T> bounds[2];
		T coords[6];
		struct {
			T x_min, y_min, z_min;
			T x_max, y_max, z_max;
		};
	};
	AAbb()
	: x_min(std::numeric_limits<T>::max())
	, y_min(std::numeric_limits<T>::max())
	, z_min(std::numeric_limits<T>::max())
	, x_max(-std::numeric_limits<T>::max())
	, y_max(-std::numeric_limits<T>::max())
	, z_max(-std::numeric_limits<T>::max())
	{}
	AAbb(T x0,T y0,T z0,T x1,T y1,T z1)
	: x_min(x0),y_min(y0),z_min(z0)
	, x_max(x1),y_max(y1),z_max(z1)
	{}
	AAbb(const Vec3<T>& l,const Vec3<T>& u) : mn(l), mx(u) {}
	AAbb(const AAbb& rh) : mn(rh.mn), mx(rh.mx) {}
	void operator=(const AAbb& rh) { mn = rh.mn; mx = rh.mx; }
};

typedef AAbb<float> AAbbf;
typedef AAbb<double> AAbbd;


// Class to manipulate and display a bounding box
class BBox3D : public AAbbf
{
  bool _aInitialized;
public:

  BBox3D();
  BBox3D(const BBox3D& pBBox3D);
  explicit BBox3D(const Vec3D&   pPoint);
  explicit BBox3D(const fvmath::Vec3f& pCenter, const float size);
  explicit BBox3D(const fvmath::Vec3f& pLower,const fvmath::Vec3f& pUpper);
  explicit BBox3D(const Vec3D& pLower, const Vec3D& pUpper);
  template<class T>
  explicit BBox3D(const fvmath::Vec3<T>& l,const fvmath::Vec3<T>& u);
  ~BBox3D        () {;}
  Vec3f & operator [](size_t idx) { return bounds[idx]; }
  const Vec3f & operator [](size_t idx) const { return bounds[idx]; }
  void Reset();

  void         dumpCharacteristics         (std::ostream&      pOstream,
                                            const std::string& pIndentation,
                                            const Matrix<float>&   pTransformation);

  void         render                      () const;

  //Vec3D		getCenter                   () const;
  CVec3f    getCenter() const;
  template<class T>
  fvmath::CVec3<T> center() const;
  float     getCircumscribedSphereRadius() const;

  float		Xmin() const { return x_min; }
  float		Ymin() const { return y_min; }
  float		Zmin() const { return z_min; }
  float		Xmax() const { return x_max; }
  float		Ymax() const { return y_max; }
  float     Zmax() const { return z_max; }



  void		minX(float x) { x_min = x; }
  void		maxX(float x) { x_max = x; }
  void		minY(float y) { y_min = y; }
  void		maxY(float y) { y_max = y; }
  void		minZ(float z) { z_min = z; }
  void		maxZ(float z) { z_max = z; }

  Vec3D   lower() const {
	  return Vec3D(x_min,y_min,z_min);
  }

  Vec3D   upper() const {
	  return Vec3D(x_max,y_max,z_max);
  }
  template<class T>
  void    move(const fvmath::Vec3<T>& mv);

  bool	    IsInside(const Vec3D& pt) const
  {
	  return( (pt._x() >= (x_min - 0.001)) && (pt._x() <= (x_max + 0.001)) &&
		      (pt._y() >= (y_min - 0.001)) && (pt._y() <= (y_max + 0.001)) &&
			  (pt._z() >= (z_min - 0.001)) && (pt._z() <= (z_max + 0.001)) );
  }
  // m
  template<class T>
  bool IsInside(const fvmath::Vec3<T>& pt) const
  {
	  return( (pt.x >= (x_min - 0.001)) && (pt.x <= (x_max + 0.001)) &&
			  (pt.y >= (y_min - 0.001)) && (pt.y <= (y_max + 0.001)) &&
			  (pt.z >= (z_min - 0.001)) && (pt.z <= (z_max + 0.001)) );
  }

  BBox3D&  operator= (const BBox3D& pBBox3D);
  BBox3D   operator+ (const BBox3D& pBBox3D) const;
  BBox3D   operator+ (const Vec3D&  pPoint) const;
  BBox3D&  operator+=(const BBox3D& pBBox3D);
  BBox3D&  operator+=(const Vec3D&  pPoint      );
  template<class T>
  inline BBox3D&  operator+=(const fvmath::Vec3<T>&  pPoint)
  {
	  return *this;
  }

  template<class T>
  inline BBox3D&  operator+=(const fvmath::Vec3<T>*  pPoint)
  {
	  this->bounds[0] = fv_min(pPoint[0],this->bounds[0]);
	  this->bounds[1] = fv_min(pPoint[1],this->bounds[1]);
	  this->bounds[2] = fv_min(pPoint[2],this->bounds[2]);
	  this->bounds[3] = fv_max(pPoint[0],this->bounds[3]);
	  this->bounds[4] = fv_max(pPoint[1],this->bounds[4]);
	  this->bounds[5] = fv_max(pPoint[2],this->bounds[5]);
	  return *this;
  }

  inline uint8_t majorAxis(float *pSize) const {
	  pSize[0] = x_max - x_min;
	  pSize[1] = y_max - y_min;
	  pSize[2] = z_max - z_min;
	  uint8_t major = 0;
	  if (pSize[0] < pSize[1]) major = 1;
	  if (pSize[1] < pSize[1]) major = 2;
	  return major;
  }

  inline float maxDim() const {

	  if(!_aInitialized) return 0.0f;
	  const float ux = fvmath::abs(x_max-x_min);
	  const float uy = fvmath::abs(y_max-y_min);
	  const float uz = fvmath::abs(z_max-z_min);

	  return fv_max(fv_max(ux,uy),uz);
  }

  inline float minDim() const {
	  if(!_aInitialized) return 0.0f;
  	  const float ux = fvmath::abs(x_max-x_min);
  	  const float uy = fvmath::abs(y_max-y_min);
  	  const float uz = fvmath::abs(z_max-z_min);

  	  return fv_min(fv_min(ux,uy),uz);
    }
  inline float * data(void) { return coords; }
  inline bool isInitialized() const { return _aInitialized; }
  inline void setInitialized(const float flg=true) { _aInitialized = flg; }
  //template<typename T>
  bool intersect(const Ray<float>& r) const;

public:
  class Edge {
  public:
	  Edge() : _begin(0), _end(0) {}
	  Edge(unsigned index0, unsigned index1) : _begin(index0), _end(index1){}
	  unsigned int begin() const { return _begin; }
	  unsigned int end() const { return _end; }

  private:
		  unsigned int _begin, _end;
  };

  template<typename T = float>
  static uint8_t  InitVertices(const BBox3D* pBBox3D,CVec3<T> *pVertices);
private:

  friend BBox3D operator*(const Matrix<float>& pMatrix,
                          const BBox3D& pBBox3D);
  friend BBox3D operator*(const BBox3D& pBBox3D, const float zoom);

  friend std::ostream& operator <<(std::ostream& os, const BBox3D& pBBox);

};

template<class T>
inline BBox3D::BBox3D(const fvmath::Vec3<T>& l,const fvmath::Vec3<T>& u)
: AAbb(l,u)
, _aInitialized(true)
{
	FV_ASSERT(x_max >= x_min);
	FV_ASSERT(y_max >= y_min);
	FV_ASSERT(z_max >= z_min);
}

template<class T>
fvmath::CVec3<T> BBox3D::center() const
{
	fvmath::CVec3<T> lcen;

	if (_aInitialized) {
		lcen.x = 0.5 * (x_min + x_max);
		lcen.y = 0.5 * (y_min + y_max);
		lcen.z = 0.5 * (z_min + z_max);
	}

	return lcen;
}

template<class T>
inline void BBox3D::move(const fvmath::Vec3<T>& mv)
{
	if (_aInitialized) {
		x_min += mv.x;
		x_max += mv.x;
		y_min += mv.y;
		y_max += mv.y;
		z_min += mv.z;
		z_max += mv.z;
	}
}

template<>
inline BBox3D& BBox3D::operator+=(const fvmath::Vec3f& pPoint)
{
	 if (_aInitialized) {
	    x_min = fv_min(x_min, pPoint.x);
	    x_max = fv_max(x_max, pPoint.x);
	    y_min = fv_min(y_min, pPoint.y);
	    y_max = fv_max(y_max, pPoint.y);
	    z_min = fv_min(z_min, pPoint.z);
	    z_max = fv_max(z_max, pPoint.z);
	  }
	  else {
	    x_min = x_max = pPoint.x;
	    y_min = y_max = pPoint.y;
	    z_min = z_max = pPoint.z;
	    _aInitialized = true;
	  }

	  FV_ASSERT(x_max >= x_min);
	  FV_ASSERT(y_max >= y_min);
	  FV_ASSERT(z_max >= z_min);

	  return *this;
}

template<>
inline BBox3D& BBox3D::operator+=(const fvmath::Vec3d& pPoint)
{
	float x = static_cast<float>(pPoint.x);
	float y = static_cast<float>(pPoint.y);
	float z = static_cast<float>(pPoint.z);
	if (_aInitialized) {
		x_min = fv_min(x_min, x);
		x_max = fv_max(x_max, x);
		y_min = fv_min(y_min, y);
		y_max = fv_max(y_max, y);
		z_min = fv_min(z_min, z);
		z_max = fv_max(z_max, z);
	 }
	 else {
		x_min = x_max = x;
		y_min = y_max = y;
		z_min = z_max = z;
		_aInitialized = true;
	  }

	  FV_ASSERT(x_max >= x_min);
	  FV_ASSERT(y_max >= y_min);
	  FV_ASSERT(z_max >= z_min);

	  return *this;

}
//template<typename T>
inline bool BBox3D::intersect(const Ray<float>& r) const
{
	float tmin, tmax, tymin, tymax, tzmin, tzmax;
	tmin = (bounds[r.sign[0]].x - r.orig.x) * r.invdir.x;
	tmax = (bounds[1 - r.sign[0]].x - r.orig[0]) * r.invdir.x;
	tymin = (bounds[r.sign[1]].y - r.orig.y) * r.invdir.y;
	tymax = (bounds[1 - r.sign[1]].y - r.orig.y) * r.invdir.y;
	if ((tmin > tymax) || (tymin > tmax))
		return false;
	if (tymin > tmin)
		tmin = tymin;
	if (tymax < tmax)
		tmax = tymax;
	tzmin = (bounds[r.sign[2]].z - r.orig.z) * r.invdir.z;
	tzmax = (bounds[1-r.sign[2]].z - r.orig.z) * r.invdir.z;
	if ((tmin > tzmax) || (tzmin > tmax))
		return false;
	if (tzmin > tmin)
		tmin = tzmin;
	if (tzmax < tmax)
		tmax = tzmax;
	if (tmin > r.tmin) r.tmin = tmin;
	if (tmax < r.tmax) r.tmax = tmax;

	return true;
}

template<typename T>
uint8_t BBox3D::InitVertices(const BBox3D *pBBox3D,CVec3<T>* pVertices)
{
	pVertices[0] = CVec3<T>(pBBox3D->Xmin(),pBBox3D->Ymin(),pBBox3D->Zmin());
	pVertices[1] = CVec3<T>(pBBox3D->Xmax(),pBBox3D->Ymin(),pBBox3D->Zmin());
	pVertices[2] = CVec3<T>(pBBox3D->Xmax(),pBBox3D->Ymax(),pBBox3D->Zmin());
	pVertices[3] = CVec3<T>(pBBox3D->Xmin(),pBBox3D->Ymax(),pBBox3D->Zmin());
	pVertices[4] = CVec3<T>(pBBox3D->Xmin(),pBBox3D->Ymin(),pBBox3D->Zmax());
	pVertices[5] = CVec3<T>(pBBox3D->Xmax(),pBBox3D->Ymin(),pBBox3D->Zmax());
	pVertices[6] = CVec3<T>(pBBox3D->Xmax(),pBBox3D->Ymax(),pBBox3D->Zmax());
	pVertices[7] = CVec3<T>(pBBox3D->Xmin(),pBBox3D->Ymax(),pBBox3D->Zmax());
	return(8);
}

//template<> bool BBox3D::intersect(const Ray<float>& r) const;

// OpenMP version of getting bounding box of model volume
class Mesh;
extern const BBox3D::Edge edges[];
extern BBox3D getModelBoundingBox(void);
extern const BBox3D& getModelBoundingBoxd_omp(Mesh* pmesh);
extern void min_sub_tab(long start_point, long n_points, double *l_mnmx);

//template<> bool BBox3D::intersect(const Ray<float>& r) const;

} // end namespace FemViewer
#endif // _BOUNDING_BOX_H_
