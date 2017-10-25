#ifndef _FIRE_H_
#define _FIRE_H_

#include <vector>
#include <sstream>
#include <string>

//#include "GraphicElem.hpp"

namespace FemViewer {

class Mesh;
class Legend;

template <typename T> class Point3D;
template <typename T> class GraphElement2;

class Fire
{
	static const int indexCase[6][3];

public:

	const Mesh* 	pMesh;
	const Legend& 	legend;
public:

	Fire(const Mesh* meshp, const Legend& leg)
			: pMesh(meshp), legend(leg) {}


	//! Nowa wersja dla przeciezc
	void Exe3(GraphElement2<double>& elIn, std::vector<GraphElement2<double> >& elOut);

	//! Tworzy linie kontu
	void CreateLine(GraphElement2<double>& elIn, std::vector<GraphElement2<double> >& elOut);

	
	// sortuje vierzcho�ki w tr�jk�cie od min do max warto�ci
	int GetSorted(const double* values, const int*& index);

protected:
	void ExtractPentagram2(GraphElement2<double>& el, const std::vector<double>& values, std::vector<GraphElement2<double> >& elOut,
			std::vector<GraphElement2<double> >& subTriangles, int& line1, int& line2);

	// procedura do kolorowania sub-triangles
	void CutTriangle(GraphElement2<double>& elIn, const std::vector<double>& values, int line, std::vector<GraphElement2<double> >&elOut);

	// fill triangle
	bool FullTriangleFallage(GraphElement2<double>& elIn, const std::vector<double>& values, std::vector<GraphElement2<double> >& elOut);

	// check if value is between contour bounds
	bool ValueWithin(const double& val, const double& bnd1, const double& bnd2)
	{ return (val >= bnd1 && val <= bnd2) || (val >= bnd2 && val <= bnd1); }

	// get contour point
	template<typename T>
	static Point3D<T> GetContourPoint(Point3D<T>& p1, Point3D<T>& p2, T& v, T& v1, T& v2)
	{
		T t(0);
		T diff = v2 - v1;

		if(diff != 0) {
			t = (v - v1) / diff;
			Point3D<T> p;
			p.SetInterp(p1,p2,t);
			
			return(p);
		}
	}
};


} // end namespace FemViewer


#endif /* _FIRE_H_
*/
