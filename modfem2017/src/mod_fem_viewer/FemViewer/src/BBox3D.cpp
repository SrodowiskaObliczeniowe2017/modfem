#include "fv_assert.h"
#include "fv_timer.h"
#include "defs.h"
#include "Log.h"
#include "BBox3D.h"
#include "Matrix.h"
#include "Ray.h"
#include "Mesh.h"
#include "mmh_intf.h"
#include <cmath>
#include <ostream>
#include <string>
#include <algorithm>
#include <GL/gl.h>
#include <GL/glu.h>
#include <omp.h>


namespace FemViewer
{

const BBox3D::Edge edges[] = {
		BBox3D::Edge(0, 1),
		//Edge(1, 3),
		BBox3D::Edge(1, 2),
		//Edge(3, 2),
		BBox3D::Edge(2, 3),
		//Edge(0, 2),
		BBox3D::Edge(0, 3),
		BBox3D::Edge(0, 4),
		BBox3D::Edge(1, 5),
		//Edge(3, 7),
		BBox3D::Edge(2, 6),
		//Edge(2, 6),
		BBox3D::Edge(3, 7),
		BBox3D::Edge(4, 5),
		//Edge(5, 7),
		BBox3D::Edge(5, 6),
		//Edge(6, 7),
		BBox3D::Edge(7, 6),
		//Edge(4, 6)
		BBox3D::Edge(4, 7)
};


BBox3D::BBox3D()
: AAbb(), _aInitialized(false)
{}

BBox3D::BBox3D(const BBox3D& pBBox3D)
: AAbb(pBBox3D), _aInitialized(pBBox3D._aInitialized)
{
  FV_ASSERT(x_max >= x_min || !_aInitialized);
  FV_ASSERT(y_max >= y_min || !_aInitialized);
  FV_ASSERT(z_max >= z_min || !_aInitialized);
}

/*BBox3D::BBox3D(const Vec3D& pPoint)
: _aInitialized  (true),
    x_max         (pPoint[0]),
    y_max         (pPoint[1]),
    z_max         (pPoint[2]),
    x_min         (pPoint[0]),
    y_min         (pPoint[1]),
    z_min         (pPoint[2])
{
	FV_ASSERT(x_max >= x_min);
	FV_ASSERT(y_max >= y_min);
	FV_ASSERT(z_max >= z_min);
}*/

BBox3D::BBox3D(const fvmath::Vec3f& Center,const float halfsize)
: _aInitialized(true)
{
	FV_ASSERT(halfsize > 0);
	const CVec3f offset(halfsize,halfsize,halfsize);
	this->mn = Center - offset;
	this->mx = Center + offset;
}

BBox3D::BBox3D(const fvmath::Vec3f& pLower,const fvmath::Vec3f& pUpper)
: AAbb(pLower,pUpper), _aInitialized(true)
{
	FV_ASSERT(x_max >= x_min);
	FV_ASSERT(y_max >= y_min);
	FV_ASSERT(z_max >= z_min);
}

BBox3D::BBox3D(const Vec3D& pLower,const Vec3D& pUpper)
: AAbb(pLower._x(),pLower._y(),pLower._z(),
	   pUpper._x(),pUpper._y(),pUpper._z())
, _aInitialized(true)
{
	FV_ASSERT(x_max >= x_min);
	FV_ASSERT(y_max >= y_min);
	FV_ASSERT(z_max >= z_min);
}


void BBox3D::Reset()
{
	x_min = FLT_MAX;
	y_min = FLT_MAX;
	z_min = FLT_MAX;
	x_max = -FLT_MAX;
    y_max = -FLT_MAX;
	z_max = -FLT_MAX;
	_aInitialized = false;
}


// Dump in ASCII the caracteristics of the BBox3D
void BBox3D::dumpCharacteristics(std::ostream&       pOstream,
                                      const std::string&  pIndentation,
                                      const Matrix<float>&    pTransformation)
{
  pOstream << pIndentation << "BBox3D " << std::endl;
  std::string lIndentation = pIndentation + " ";

#ifdef FV_DUMP_MEMORY_USAGE

  pOstream << lIndentation << "Memory used by the BBox3D = " << sizeof(*this) << std::endl;

#endif // #ifdef FD_DUMP_MEMORY_USAGE

  if (_aInitialized) {

    // Compute the transformed BBox3D
    BBox3D lTransformed = pTransformation * (*this);

    pOstream << lIndentation << "[MinX, MaxX] = [" << lTransformed.x_min << ", " << lTransformed.x_max << "]" << std::endl;
    pOstream << lIndentation << "[MinY, MaxY] = [" << lTransformed.y_min << ", " << lTransformed.y_max << "]" << std::endl;
    pOstream << lIndentation << "[MinZ, MaxZ] = [" << lTransformed.z_min << ", " << lTransformed.z_max << "]" << std::endl;
  }
  else {
    pOstream << lIndentation << "Empty" << std::endl;
  }
}

void BBox3D::render() const
{

  if (_aInitialized) {

    FV_ASSERT(x_max >= x_min);
    FV_ASSERT(y_max >= y_min);
    FV_ASSERT(z_max >= z_min);

    // We push the attributes on the openGL attribute stack; so
    // we dont disturb the current values of color (current) or line width (line)
    glPushAttrib(GL_CURRENT_BIT | GL_LINE_BIT | GL_LIGHTING_BIT);
    // Warning: performance hit; when a lot of push(GL_LIGHTING_BIT) calls are made
    glDisable(GL_LIGHTING);
    glLineWidth(2);

    glBegin(GL_LINES);
    {
      glColor3f(0.7f, 0.5f, 0.1f);

      glVertex3f(x_min,y_min,z_min); glVertex3f(x_max,y_min,z_min);
      glVertex3f(x_min,y_min,z_min); glVertex3f(x_min,y_max,z_min);
      glVertex3f(x_min,y_min,z_min); glVertex3f(x_min,y_min,z_max);
      glVertex3f(x_max,y_min,z_min); glVertex3f(x_max,y_max,z_min);
      glVertex3f(x_max,y_min,z_min); glVertex3f(x_max,y_min,z_max);
      glVertex3f(x_max,y_max,z_min); glVertex3f(x_min,y_max,z_min);
      glVertex3f(x_max,y_max,z_min); glVertex3f(x_max,y_max,z_max);
      glVertex3f(x_min,y_max,z_min); glVertex3f(x_min,y_max,z_max);
      glVertex3f(x_min,y_min,z_max); glVertex3f(x_max,y_min,z_max);
      glVertex3f(x_min,y_min,z_max); glVertex3f(x_min,y_max,z_max);
      glVertex3f(x_max,y_min,z_max); glVertex3f(x_max,y_max,z_max);
      glVertex3f(x_min,y_max,z_max); glVertex3f(x_max,y_max,z_max);
    }
    glEnd();

    // Revert the lighting state and the line state
    glPopAttrib();
  }
}

/*Vec3D BBox3D::getCenter() const
{
	Vec3D lCenter;

	if(_aInitialized) {

		lCenter = Vec3D(0.5f*(x_min + x_max),
                        0.5f*(y_min + y_max),
                        0.5f*(z_min + z_max));
	}
	return lCenter;
}*/

CVec3f BBox3D::getCenter() const
{
	CVec3f lCenter;
	if (_aInitialized) {
		lCenter.x = .5f*(x_min + x_max);
		lCenter.y = .5f*(y_min + y_max);
		lCenter.z = .5f*(z_min + z_max);
	}

	return lCenter;
}

float BBox3D::getCircumscribedSphereRadius() const
{
  float lRadius = 0.0f;

  if (_aInitialized) {
    FV_ASSERT(x_max >= x_min);
    FV_ASSERT(y_max >= y_min);
    FV_ASSERT(z_max >= z_min);

    float lDX = 0.5f*(x_max - x_min);
    float lDY = 0.5f*(y_max - y_min);
    float lDZ = 0.5f*(z_max - z_min);
    lRadius = sqrt(lDX*lDX + lDY*lDY + lDZ*lDZ);
  }
  return lRadius;
}

BBox3D& BBox3D::operator=(const BBox3D& pBBox3D)
{
  if (&pBBox3D != this) {

    _aInitialized = pBBox3D._aInitialized;
    x_max        = pBBox3D.x_max;
    y_max        = pBBox3D.y_max;
    z_max        = pBBox3D.z_max;
    x_min        = pBBox3D.x_min;
    y_min        = pBBox3D.y_min;
    z_min        = pBBox3D.z_min;

    FV_ASSERT(x_max >= x_min || !_aInitialized);
    FV_ASSERT(y_max >= y_min || !_aInitialized);
    FV_ASSERT(z_max >= z_min || !_aInitialized);
  }

  return *this;
}

BBox3D BBox3D::operator+(const BBox3D& pBBox3D) const
{
  BBox3D lBBox3D;

  if (_aInitialized && pBBox3D._aInitialized) {
    lBBox3D.x_max = std::max(x_max, pBBox3D.x_max);
    lBBox3D.y_max = std::max(y_max, pBBox3D.y_max);
    lBBox3D.z_max = std::max(z_max, pBBox3D.z_max);
    lBBox3D.x_min = std::min(x_min, pBBox3D.x_min);
    lBBox3D.y_min = std::min(y_min, pBBox3D.y_min);
    lBBox3D.z_min = std::min(z_min, pBBox3D.z_min);
    lBBox3D._aInitialized = true;

    FV_ASSERT(lBBox3D.x_max >= lBBox3D.x_min);
    FV_ASSERT(lBBox3D.y_max >= lBBox3D.y_min);
    FV_ASSERT(lBBox3D.z_max >= lBBox3D.z_min);
  }
  else if (_aInitialized) {
    lBBox3D.x_max = x_max;
    lBBox3D.y_max = y_max;
    lBBox3D.z_max = z_max;
    lBBox3D.x_min = x_min;
    lBBox3D.y_min = y_min;
    lBBox3D.z_min = z_min;
    lBBox3D._aInitialized = true;

    FV_ASSERT(lBBox3D.x_max >= lBBox3D.x_min);
    FV_ASSERT(lBBox3D.y_max >= lBBox3D.y_min);
    FV_ASSERT(lBBox3D.z_max >= lBBox3D.z_min);
  }
  else if (pBBox3D._aInitialized) {
    lBBox3D.x_max = pBBox3D.x_max;
    lBBox3D.y_max = pBBox3D.y_max;
    lBBox3D.z_max = pBBox3D.z_max;
    lBBox3D.x_min = pBBox3D.x_min;
    lBBox3D.y_min = pBBox3D.y_min;
    lBBox3D.z_min = pBBox3D.z_min;
    lBBox3D._aInitialized = true;

    FV_ASSERT(lBBox3D.x_max >= lBBox3D.x_min);
    FV_ASSERT(lBBox3D.y_max >= lBBox3D.y_min);
    FV_ASSERT(lBBox3D.z_max >= lBBox3D.z_min);
  }

  return lBBox3D;
}

BBox3D BBox3D::operator+(const Vec3D& pPoint) const
{
  BBox3D lBBox3D;

  if (_aInitialized) {
    lBBox3D.x_max = std::max(x_max, pPoint[0]);
    lBBox3D.y_max = std::max(y_max, pPoint[1]);
    lBBox3D.z_max = std::max(z_max, pPoint[2]);
    lBBox3D.x_min = std::min(x_min, pPoint[0]);
    lBBox3D.y_min = std::min(y_min, pPoint[1]);
    lBBox3D.z_min = std::min(z_min, pPoint[2]);
  }
  else {
    lBBox3D.x_max = pPoint[0];
    lBBox3D.y_max = pPoint[1];
    lBBox3D.z_max = pPoint[2];
    lBBox3D.x_min = pPoint[0];
    lBBox3D.y_min = pPoint[1];
    lBBox3D.z_min = pPoint[2];
  }

  lBBox3D._aInitialized = true;

  FV_ASSERT(lBBox3D.x_max >= lBBox3D.x_min);
  FV_ASSERT(lBBox3D.y_max >= lBBox3D.y_min);
  FV_ASSERT(lBBox3D.z_max >= lBBox3D.z_min);

  return lBBox3D;
}

BBox3D& BBox3D::operator+=(const BBox3D& pBBox3D)
{
  if (_aInitialized && pBBox3D._aInitialized) {
    x_max = std::max(x_max, pBBox3D.x_max);
    y_max = std::max(y_max, pBBox3D.y_max);
    z_max = std::max(z_max, pBBox3D.z_max);
    x_min = std::min(x_min, pBBox3D.x_min);
    y_min = std::min(y_min, pBBox3D.y_min);
    z_min = std::min(z_min, pBBox3D.z_min);

    FV_ASSERT(x_max >= x_min);
    FV_ASSERT(y_max >= y_min);
    FV_ASSERT(z_max >= z_min);
  }
  else if (pBBox3D._aInitialized) {
    x_max = pBBox3D.x_max;
    y_max = pBBox3D.y_max;
    z_max = pBBox3D.z_max;
    x_min = pBBox3D.x_min;
    y_min = pBBox3D.y_min;
    z_min = pBBox3D.z_min;
    _aInitialized = true;

    FV_ASSERT(x_max >= x_min);
    FV_ASSERT(y_max >= y_min);
    FV_ASSERT(z_max >= z_min);
  }

  return *this;
}


BBox3D& BBox3D::operator+=(const Vec3D& pPoint)
{
  if (_aInitialized) {
    x_max = std::max(x_max, pPoint[0]);
    y_max = std::max(y_max, pPoint[1]);
    z_max = std::max(z_max, pPoint[2]);
    x_min = std::min(x_min, pPoint[0]);
    y_min = std::min(y_min, pPoint[1]);
    z_min = std::min(z_min, pPoint[2]);
  }
  else {
    x_max = pPoint[0];
    y_max = pPoint[1];
    z_max = pPoint[2];
    x_min = pPoint[0];
    y_min = pPoint[1];
    z_min = pPoint[2];
    _aInitialized = true;
  }

  FV_ASSERT(x_max >= x_min);
  FV_ASSERT(y_max >= y_min);
  FV_ASSERT(z_max >= z_min);

  return *this;
}




//template<typename T>
//bool BBox3D::intersect(const Ray<T>& r) const
//{
//	T tmin, tmax, tymin, tymax, tzmin, tzmax;
//	tmin = (bounds[r.sign[0]].x - r.orig.x) * r.invdir.x;
//	tmax = (bounds[1 - r.sign[0]].x - r.orig[0]) * r.invdir.x;
//	tymin = (bounds[r.sign[1]].y - r.orig.y) * r.invdir.y;
//	tymax = (bounds[1 - r.sign[1]].y - r.orig.y) * r.invdir.y;
//	if ((tmin > tymax) || (tymin > tmax))
//		return false;
//	if (tymin > tmin)
//		tmin = tymin;
//	if (tymax < tmax)
//		tmax = tymax;
//	tzmin = (bounds[r.sign[2]].z - r.orig.z) * r.invdir.z;
//	tzmax = (bounds[1-r.sign[2]].z - r.orig.z) * r.invdir.z;
//	if ((tmin > tzmax) || (tzmin > tmax))
//		return false;
//	if (tzmin > tmin)
//		tmin = tzmin;
//	if (tzmax < tmax)
//		tmax = tzmax;
//	if (tmin > r.tmin) r.tmin = tmin;
//	if (tmax < r.tmax) r.tmax = tmax;
//
//	return true;
//}
//
//template<> bool BBox3D::intersect(const Ray<float>& r) const;

// Transform all the corners of the BBox3D
// and modify the current BBox3D to include all the
// transformed points. This is overkill if we only use translations
// and scalings, but it will work for any linear transformation
// if we add support to let's say glrotate.
BBox3D operator*(const Matrix<float>& pMatrix, const BBox3D& pBBox3D)
{
  BBox3D lBBox3D;

  if (pBBox3D._aInitialized) {

    lBBox3D += pMatrix * Vec3D(pBBox3D.x_min, pBBox3D.y_min, pBBox3D.z_min);
    lBBox3D += pMatrix * Vec3D(pBBox3D.x_max, pBBox3D.y_min, pBBox3D.z_min);
    lBBox3D += pMatrix * Vec3D(pBBox3D.x_max, pBBox3D.y_max, pBBox3D.z_min);
    lBBox3D += pMatrix * Vec3D(pBBox3D.x_min, pBBox3D.y_max, pBBox3D.z_min);
    lBBox3D += pMatrix * Vec3D(pBBox3D.x_min, pBBox3D.y_min, pBBox3D.z_max);
    lBBox3D += pMatrix * Vec3D(pBBox3D.x_max, pBBox3D.y_min, pBBox3D.z_max);
    lBBox3D += pMatrix * Vec3D(pBBox3D.x_max, pBBox3D.y_max, pBBox3D.z_max);
    lBBox3D += pMatrix * Vec3D(pBBox3D.x_min, pBBox3D.y_max, pBBox3D.z_max);
  }

  return lBBox3D;
}

BBox3D operator*(const BBox3D& pBBox3D, const float zoom)
{
	BBox3D lBBox3D;

	float zoom1;
	if(zoom == 0.0f) zoom1 = 1.0f;
	else zoom1= 1/ zoom;

  if (pBBox3D._aInitialized) {

	  lBBox3D += pBBox3D.x_min * zoom1; 
	  lBBox3D += pBBox3D.x_max * zoom1;
	  lBBox3D += pBBox3D.y_min * zoom1;
	  lBBox3D += pBBox3D.y_max * zoom1;
	  lBBox3D += pBBox3D.z_min * zoom1;
	  lBBox3D += pBBox3D.z_max * zoom1;
  }

  return lBBox3D;
}

//operator double(const BBox3D& pBox3D) const
//{
//	double * bounds = new double[6];
//
//	if(pBox3D._aInitialized) {
//		bounds[0] = static_cast<double>(pBox3D.x_min);
//		bounds[1] = static_cast<double>(pBox3D.y_min);
//		bounds[2] = static_cast<double>(pBox3D.z_min);
//		bounds[3] = static_cast<double>(pBox3D.x_max);
//		bounds[4] = static_cast<double>(pBox3D.y_max);
//		bounds[5] = static_cast<double>(pBox3D.z_max);
//	} else {
//		bounds[0] = DBL_MIN;
//		bounds[1] = DBL_MIN;
//		bounds[2] = DBL_MIN;
//		bounds[3] = DBL_MAX;
//		bounds[4] = DBL_MAX;
//		bounds[5] = DBL_MAX;
//	}
//	return(bounds);
//}

std::ostream& operator << (std::ostream& os, const BBox3D& pBBox)
{
	if (!pBBox.isInitialized())
		os << "Not initialized\n";
	else
		os << "Bounding box boundaries:\n"
		   << "min=(" << pBBox.x_min << ", " << pBBox.y_min<< ", " <<pBBox.z_min << ")\n"
		   << "max=(" << pBBox.x_max << ", " << pBBox.y_max<< ", " <<pBBox.z_max << ")\n";
	return os;
}

BBox3D getModelBoundingBox(void)
{
	BBox3D lbbox;
	fvmath::Vec3d lcoord;
	int mesh_id = MMC_CUR_MESH_ID;
	int nno = mmr_get_max_node_id(MMC_CUR_MESH_ID);
	int i;
	fv_timer t;
	t.start();
	#pragma omp parallel default(none) private(lcoord,mesh_id,nno,i) shared(lbbox)
	{
		#pragma omp single
		{
			printf("Program jest wykonywany na %d watkach.\n", omp_get_num_threads());
		}
		#pragma omp parallel for schedule(static)
		for (i = 1; i <= nno; ++i) {
			if (mmr_node_status(mesh_id,i)) {
				mmr_node_coor(MMC_CUR_MESH_ID,i,lcoord.v);
				lbbox += lcoord;
			}
		}
	}
	t.stop();
	std::cout << "End with " << t.get_time_us() << std::endl;

	return lbbox;
}

void min_sub_tab(long start_point, long n_points, double *l_mnmx)
{
  double mn[3], mx[3], node[3];
  int sid = 0;

  for (int i = 0; i < n_points; ++i) {
	  if (mmr_node_coor(MMC_CUR_MESH_ID,start_point++,node) < 0) continue;
	  if (!i) {
		  mn[0] = mx[0] = node[0];
		  mn[1] = mx[1] = node[1];
		  mn[2] = mx[2] = node[2];
	  } else {
		  mn[0] = MIN(mn[0],node[0]);
		  mx[0] = MAX(mx[0],node[0]);
		  mn[1] = MIN(mn[1],node[1]);
		  mx[1] = MAX(mx[1],node[1]);
		  mn[2] = MIN(mn[2],node[2]);
		  mx[2] = MAX(mx[2],node[2]);
	  }
  }
  // Store in array
  l_mnmx[0] = mn[0];
  l_mnmx[1] = mn[1];
  l_mnmx[2] = mn[2];
  l_mnmx[3] = mx[0];
  l_mnmx[4] = mx[1];
  l_mnmx[5] = mx[2];
  //printf("elem %d start idx %d\n",n_points,start_point);
  //printf("BBox model: min = {%lf, %lf, %lf}\n",xmn,ymn,zmn);
  //printf("BBox model: max = {%lf, %lf, %lf}\n",xmx,ymx,zmx);
}



} // end namespace FemViewer
