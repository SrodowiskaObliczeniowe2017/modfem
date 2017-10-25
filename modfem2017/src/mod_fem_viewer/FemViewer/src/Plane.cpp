#include "Plane.h"
#include "BBox3D.h"

#include <iostream>
#include<omp.h>

namespace FemViewer {

#define FV_PLANE_TOL 1.0e-06
#define FV_LARGE_DOUBlE  10.0e+06


int Plane::Check(const double n[4],const double x[3])
{
	double res = Evaluate(n,x);
	if (fvmath::Compare(res,0.0)) return( 0);
	if (res > 0.0)				  return( 1);
								  return(-1);
}

// Given a line defined by the two points p1,p2; and a plane defined by the
// normal n and point p0, compute an intersection. The parametric
// coordinate along the line is returned in u, and the coordinates of
// intersection are returned in p2. A zero is returned if the plane and line
// do not intersect between (0<=u<=1). If the plane and line are parallel,
// zero is returned and u is set to FV_LARGE_DOUBlE.
int Plane::IntersectWithLine(const double pn[4], double p2[3], double p1[3], double & u)
{
	double dot, val, v[3];
	double fabsden, fabstolerance;

	// Compute line vector 
	v[0] = p2[0] - p1[0];
	v[1] = p2[1] - p1[1];
	v[2] = p2[2] - p1[2];

	// Compute denominator.  If ~0, line and plane are parallel.
	// 
	dot = pn[0]*v[0]+pn[1]*v[1]+pn[2]*v[2];

	if (dot == 0.0)
	{
		u = FV_LARGE_DOUBlE;
		return 0;
	}


	// If denominator with respect to numerator is "zero", then the line and
	// plane are considered parallel. 
	//

	// trying to avoid an expensive call to fabs()
	if (dot < 0.0)
    {
		fabsden = -dot;
    }
	else
    {
		fabsden = dot;
    }
	val = Evaluate(pn,p1);
	if (val < 0.0)
    {
		fabstolerance = -val*FV_PLANE_TOL;
    }
	else
    {
		fabstolerance = val*FV_PLANE_TOL;
    }
	if (fabsden <= fabstolerance)
    {
		u = FV_LARGE_DOUBlE;
		return 0;
    }

	// valid intersection
	u = -val / dot;
	//std::cout<<"u= " << u << std::endl;
	p2[0] = p1[0] + u*v[0];
	p2[1] = p1[1] + u*v[1];
	p2[2] = p1[2] + u*v[2];
	//std::cout<<"wektor p0 = ("<< p0[0] <<", "<< p0[1] <<", "<< p0[2] << ")\n";
	if ( u >= 0.0 && u <= 1.0 )
    {
		return 1;
    }
	else
    {
		return 0;
    }
}

int  Plane::GetCuttedLine(double pn[4], double p2[3], double p1[3], double p0[3])
{
	//int res1,res2,res;
	//double u;

	int res2 = Check(pn,p2);
	int res1 = Check(pn,p1);

	if (res2 > 0 && res1 > 0) return(-1);

	if ((res1 == 0 && res2 == 0) || (res1 < 0 && res2 < 0)) 
	{
		return ( 2);
	}

	double u;
	int res = IntersectWithLine(pn,p2,p1,u);
	if (res == 1)
	{
		if (res2 < 0) return( 0);
		else return( 1);
	}
	else return (-1);
}

#define V(b,e)	el_coords[(e)] - el_coords[(b)];
int Plane::CreatePlanesForPrism(const vec3 el_coords[6],vec4 el_planes[5])
{
	// This is for prism
	cvec3 v01 = V(0,1);
	cvec3 v12 = V(1,2);
	//printf("thr = %d v01 = {%f %f %f}\n",omp_get_thread_num(),v01.x,v01.y,v01.z);
	//printf("thr = %d v12 = {%f %f %f}\n",omp_get_thread_num(),v12.x,v12.y,v12.z);
	cvec3 n;
	int index = 0;
	fvmath::Cross(&v12,&v01,&n);
	//printf("thr = %d n = {%f %f %f}\n",omp_get_thread_num(),n.x,n.y,n.z);
	Extract(&el_coords[0],&n,&el_planes[index]);
	index++;
	vec3 v34 = V(3,4);
	vec3 v35 = V(3,5);
	fvmath::Cross(&v34,&v35,&n);
	Extract(&el_coords[3],&n,&el_planes[index]);
	index++;
	vec3 v14 = V(1,4);
	fvmath::Cross(&v01,&v14,&n);
	Extract(&el_coords[4],&n,&el_planes[index]);
	index++;
	fvmath::Cross(&v12,&v14,&n);
	Extract(&el_coords[4],&n,&el_planes[index]);
	index++;
	vec3 v25 = V(2,5);
	fvmath::Cross(&v25,&v35,&n);
	Extract(&el_coords[5],&n,&el_planes[index]);

	return index;

}


int Plane::CreatePlanesForTetra(const vec3 el_coords[4],vec4 el_planes[4])
{
	// This is for tetrahedron
	cvec3 v01 = V(0,1);
	cvec3 v12 = V(1,2);
	cvec3 n;
	int index = 0;
	fvmath::Cross(&v01,&v12,&n);
	Extract(&el_coords[0],&n,&el_planes[index]);
	index++;
	vec3 v03 = V(0,3);
	fvmath::Cross(&v03,&v01,&n);
	Extract(&el_coords[0],&n,&el_planes[index]);
	index++;
	vec3 v23 = V(2,3);
	fvmath::Cross(&v03,&v23,&n);
	Extract(&el_coords[0],&n,&el_planes[index]);
	index++;
	fvmath::Cross(&v23,&v12,&n);
	Extract(&el_coords[1],&n,&el_planes[index]);

	return index;
}
#undef V

void Plane::Extract(const vec3* pt,const vec3* n,vec4* params)
{
	params->x = (CoordType)n->x;
	params->y = (CoordType)n->y;
	params->z = (CoordType)n->z;
	params->w = (CoordType)-(pt->x*n->x + pt->y*n->y + pt->z*n->z);
}

void Plane::SetPlane(double n[3], double p[3])
{
	try {
		
		fvmath::CVec3d nv(n);

		if( 0.0 == fvmath::Normalize(nv) )
			throw "Can't normalize null vector for cut plane!";
		
		this->p.a = nv[0];
		this->p.b = nv[1];
		this->p.c = nv[2];
		this->p.d = -(n[0]*p[0]+n[1]*p[1]+n[2]*p[2]);

	}
	catch(const char* ex) {
		std::cerr << ex << std::endl;
	}
	return;

}

} // end namespace Femviewer
