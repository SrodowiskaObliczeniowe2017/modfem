#ifndef _PLANE_H_
#define _PLANE_H_

#include "fv_float.h"
#include "MathHelper.h"

#include <cstring>
//#include<iostream>


namespace FemViewer {

class Plane {

public:
  typedef fvmath::Vec3d  vec3;
  typedef fvmath::CVec3d cvec3;
  typedef fvmath::Vec4<CoordType>  vec4;

  static double DistanceToPlane(double x[3], double n[3], double p0[3]);

  static double Evaluate(const double n[4],const double x[3]);

  static int    Check(const double n[4],const double x[3]); 

  static int    IntersectWithLine(const double pn[4], double p2[3], double p1[3], double & u);

  static int    GetCuttedLine(double pn[4], double p2[3], double p1[3], double p0[3]);

  static int    CreatePlanesForPrism(const vec3 el_coords[6],vec4 el_planes[5]);

  static int    CreatePlanesForTetra(const vec3 el_coords[4],vec4 el_planes[4]);

  static void   Extract(const vec3* pt,const vec3* n,vec4* prams);

protected:

  typedef union {
    struct { double a, b, c, d; };
    double n[4];
  } params;
	
  params p;

public:
  Plane()
  {
    Reset();
  }

  Plane(const Plane& rhs)
  {
    if (this != &rhs)
      {
	memcpy(&p,&rhs.p,sizeof(params));
      }
  }

  void operator=(const Plane& rhs)
    {
      memcpy(p.n, rhs.p.n, sizeof(this->p));
    }

  void Reset()
  {
    memset(p.n, 0x0, sizeof(p));
  }

  virtual ~Plane() {}

  virtual void SetPlane(double n[3], double p[3]);

	double * GetParams() { return p.n; }
	const double * GetParams() const { return p.n; }

	double EvaluatePlane(const double x1,const double x2,const double x3) const
	{
		return x1*p.a + x2*p.b + x3*p.c + p.d;
	}

	inline double Calculate(const double *v) {
		return (p.a*v[0] + p.b*v[1] +p.c*v[2] + p.d);
	}

	inline int CheckLocation(const double& x1,const double& x2,const double& x3) const
	{
		double res = EvaluatePlane(x1,x2,x3);
		//std::cout <<"Checkloc::x= ("<< x1 << ", " << x2 << ", " << x3 <<")\n";
		//std::cout <<"CheckLoc::res= " << res << std::endl;
		if(fvmath::Compare(res,0.0)) return( 0);
		if(res < 0.0) return(-1);
		return( 1);
	}

	inline int CheckLocation(const double x[]) const
	{
		double res = EvaluatePlane(x[0],x[1],x[2]);
		if(fvmath::Compare(res,0.0)) return( 0);
		if(res < 0.0) return(-1);
		return( 1);
	}

     

};

inline double Plane::Evaluate(const double n[4],const double x[3])
{
	return (n[0]*x[0]+n[1]*x[1]+n[2]*x[2]+n[3]);
}

inline double Plane::DistanceToPlane(double x[3], double n[3], double p0[3])
{
#define abs_val(x) ((x)>0 ? (x) : -(x))
	return abs_val( n[0]*(x[0]-p0[0]) + n[1]*(x[1]-p0[1]) + n[2]*(x[2]-p0[2]) );
}

} // end naespace 

#endif /* _PLANE_H_
*/
