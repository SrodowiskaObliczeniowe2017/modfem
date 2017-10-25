#ifndef _ELEM_TABLES_H_
#define _ELEM_TABLES_H_

#include "MathHelper.h"
#include "CutPlane.h"
#include "GraphicElem.hpp"


extern int iQuadEdgesIds[4][2];
extern int iQuadFlags[16];
extern int iQuadEdges[16][7];
extern const int iPrizmEdges[9][2];
///* Square Edge Flags Table*/


typedef FemViewer::fvmath::CVec3d myVec3d;
typedef FemViewer::Node_t myVertex;
typedef FemViewer::Vertex Vertex;
typedef FemViewer::CutPlane myPlane;

struct comp {
inline bool operator()(const Vertex& a,const Vertex& b){
	if(!FemViewer::fvmath::is_near( a.position.x , b.position.x, 0.001f )) return a.position.x < b.position.x;
	if(!FemViewer::fvmath::is_near( a.position.y , b.position.y, 0.001f )) return a.position.y < b.position.y;
	if(!FemViewer::fvmath::is_near( a.position.z , b.position.z, 0.001f )) return a.position.z < b.position.z;
	return false; // Vertices are equal
}
};


//extern int  DoSliceTerahedron(const FemViewer::CutPlane* plane,const int faces[5],double tetra[16], std::vector<FemViewer::GraphElement2<double> >& grEls,FILE* fp);

//extern int  DoSlicePrizm(const FemViewer::CutPlane* plane,const int faces[5],double tetra[16], std::vector<FemViewer::GraphElement2<double> >& grEls,FILE* fp);

extern int  IsPrismSelected(const double coords[6*4],const double* plane);

extern int  IsTetraSelected(const double coords[4*4],const double* plane);

extern unsigned int  CutElement(const double vertices[18], const int nodes[7],const myPlane*  cut_plane,std::vector<myVertex>& out_vertices);

extern void indexVBO(std::vector< Vertex >& inout_vertices,std::vector< unsigned int >& inout_indices);

extern unsigned TessellateTriangle(const myVertex* a,const myVertex* b,const myVertex* c,
		int ndiv,
		std::vector<myVertex>& vertices, std::vector<unsigned>& indices);

#endif
