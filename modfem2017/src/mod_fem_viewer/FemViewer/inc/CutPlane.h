#ifndef _CUT_PLANE_H_
#define _CUT_PLANE_H_

#include "Enums.h"
#include "Vec3D.h"
#include "BBox3D.h"
#include "Plane.h"
#include "fv_txt_utls.h"

#include "Point3D.h"
#include "GraphicElem.hpp"
#include "MathHelper.h"
#include "Log.h"
#include <set>
#include <vector>
#include <ostream>
#include <cmath>


namespace FemViewer {

  // Forward class declaration
  class BBox3D;
   
  class CutPlane : public Plane
  {
  public:

	  static int Defaults(const BBox3D& bbox,std::vector<CutPlane>& out_planes);
	  static CutPlane SetPlane(fvmath::CVec3d normal,fvmath::CVec3d point);

  public:

    typedef std::pair<int,int> elem_idx;

    struct comp {
      bool operator()(const elem_idx& lhs,const elem_idx rhs) const {
	if (lhs.first < rhs.first) return lhs.second < rhs.second;
	return lhs.second < rhs.second;
      }
    };

    typedef std::set<elem_idx,comp> elem_indices;
    typedef elem_indices::iterator elem_indices_itr;

    // Default position is z = 0.5
    CutPlane(double nx=0.0,double ny=0.0, double nz=1.0, double d=-0.05, bool normalize = true)
    {
      Set(nx,ny,nz,d,normalize);
      _active = false;
      _display = false;
      _isPlaneChanged = false;
      _isNormalized = normalize;
    }

    CutPlane(double p[],bool normalize = true)
      {
	Set(p[0],p[1],p[2],p[3],normalize);
      }
		
  CutPlane(const CutPlane& rhs) : Plane(rhs)
      {
	if (this != &rhs) {
	  _active = rhs._active;
	  _display = rhs._display;
	  _isPlaneChanged = rhs._isPlaneChanged;
	  _isNormalized = rhs._isNormalized;
	}
      }

    CutPlane& operator=(const CutPlane& rhs)
      {
	if (this != &rhs) {
	  for (int i(0); i<4; ++i) p.n[i] = rhs.p.n[i];
	  _active = rhs._active;
	  _display = rhs._display;
	  _isPlaneChanged = rhs._isPlaneChanged;
	  _isNormalized = rhs._isNormalized;
	  _vplaneCorrners = rhs._vplaneCorrners;
	  _scale = rhs._scale;
	  _elementIndices = rhs._elementIndices;
	}
	return *this;
      }

    virtual ~CutPlane(){}



    
    void Enable(bool status = true) { _active = status; }
    bool IsActive() const { return _active; }
    void Show(bool status = true) { _display = status; }
    void Hide() { _display = false; }
    bool IsVisible() const { return _display; }
    void SetChanged(bool state = true) { _isPlaneChanged = state; }
    bool IsPlaneChanged() const { return _isPlaneChanged; }
    bool IsNormalized() const { return _isNormalized; }

    elem_indices& GetElementIndices() { return _elementIndices; }
    const elem_indices& GetElementIndices() const { return _elementIndices; }

    void Set(double a, double b, double c, double d, bool normalize = true)
    {
      fvmath::CVec3d v(a,b,c);
      double l = normalize ? Normalize(v) : 1.0;
      p.a = v.x; p.b = v.y; p.c = v.z;
      p.d = d / l;
    }

    /*void Set(double a, double b, double c, Point3D<double> pt)
      {
      Set(a,b,c,0.0);
      p.d = -(p.a*pt.x + p.b*pt.y + p.c*pt.z);
      }*/

    bool operator () (const Plane& rhPl)
    {
      bool res = false;
      if(p.n[0] != rhPl.GetParams()[0] || 
	 p.n[1] != rhPl.GetParams()[1] ||
	 p.n[2] != rhPl.GetParams()[2] ||
	 p.n[3] != rhPl.GetParams()[3] ) res = true;
      
      if(res) { 
	p.a = GetParams()[0];
	p.b = GetParams()[1];
	p.c = GetParams()[2];
	p.d = GetParams()[3];
      }

      return res;
    }

    bool operator==(const CutPlane& rhs);


		// -3 all in opozite side
		// -2 to 2  plane cuts the triangle
		//  3 all in the same saide
		/*int CheckWithTraingle(const double triang[], int type=3) const
		{
			int ret(0);
			for (int i(0),j(0);i<9;i+=3)
			{
				ret += CheckLocation(triang[i],triang[i+1],triang[i+2]);
			}
			return( ret);
		}*/
		
		// Check if quad is intersected by plane
		/*inline int CheckWithQuad(const double quad[]) const
		{
			int ret = CheckWithTraingle(quad);
			ret += CheckLocation(quad[9],quad[10],quad[11]);
			return ret;
		}*/

		/*inline bool CheckWithFigure(const double fig[], int& type) const
		{
			int ret = (type == 3) ? CheckWithTraingle(fig) :  CheckWithQuad(fig) ;
			return (type == 3) ? (type=ret, !(ret == 3 || ret == -3)) : (type=ret, !(ret == 4 || ret == -4));
		}*/


    inline int IntersectWithLine(double* x2,double* x1,double& u) const {
      return Plane::IntersectWithLine(this->p.n,x2,x1,u);
    }


    //int CheckElementLocation(const double x[], std::vector<GraphElement2<double> >& vels, int size);

    
    
    int CheckOrientation(const double v[]);
    
    int InversPlane(const double v[]);

    bool InitGeometryOutlines(const BBox3D& bbox);

    bool IsValid(const AAbbf& bbox) const {
      return IntersectBBox3D(bbox) == 0 ? true : false;
    }

    template<typename T>
      int IntersectBBox3D(const AAbb<T>& bbox) const;

    bool InitOutlines(const BBox3D& bbox);
    inline void Draw() const {
      if (_vplaneCorrners.empty() || !_display) return;
      drawPlane(_vplaneCorrners[0].v,_scale.v,_vplaneCorrners.size());
    }

    void AddIndex(const elem_idx&& idx);
    bool IsElement(const elem_idx&& idx) const;

    unsigned GetMaxNUmVertices(int max_div);
    unsigned GetMaxNumIndices(int max_div);

    static unsigned CutElement(const CutPlane* pl,const double vts[18],const int nodes[7],
			       const int index, int div, std::vector<Node_t>& vertices,std::vector<unsigned>& triagles);
  private:

    inline Vec3D ToVec3D(const double& a_, 
			 const double& b_, 
			 const double& c_) const;
    friend std::ostream& operator << (std::ostream& os, const CutPlane& rhs);

    bool _active;
    bool _display;
		bool _isPlaneChanged;
		bool _isNormalized;

		std::vector<CVec3f> _vplaneCorrners;
		CVec3f _scale;

		elem_indices _elementIndices;

	};



	inline Vec3D CutPlane::ToVec3D(const double& a, 
							       const double& b, 
							       const double& c) const
	{
		float _a = static_cast<float>(a);
		float _b = static_cast<float>(b);
		float _c = static_cast<float>(c);
		return Vec3D(_a, _b, _c);
	}

	// Returns:
	// -1 - box is inside
	//  0 - plane intersect box
	//  1 - box is outside the plane
	template<typename T>
	  inline int CutPlane::IntersectBBox3D(const AAbb<T>& bbox) const
	  {
		using namespace fvmath;
		CVec3<T> c = (bbox.mx + bbox.mn) * T(0.5);
		CVec3<T> h = (bbox.mx - bbox.mn) * T(0.5);
		T nx = std::fabs(p.n[0]);
		T ny = std::fabs(p.n[1]);
		T nz = std::fabs(p.n[2]);
		T e  = h.x*nx+ h.y*ny + h.z*nz;
		T s  = c.x*p.a + c.y*p.b + c.z*p.c + p.d;
		if (s - e > T(0.0)) return 1;
		if (s + e < T(0.0)) return -1;
		return 0;
	  }

	

	
} // end namespace
#endif /* _CUT_PLANE_H_
*/
