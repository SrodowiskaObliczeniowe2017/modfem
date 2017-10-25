#include "Enums.h"
#include "CutPlane.h"
#include "MathHelper.h"
#include "elem_tables.h"
#include <cassert>
#include <algorithm>


namespace FemViewer {

int CutPlane::Defaults(const BBox3D& bbox, std::vector<CutPlane>& out_planes)
{
	out_planes.clear();
	if (bbox.isInitialized()) {
		fvmath::CVec3d center(bbox.getCenter());
		//std::cout << center << "\n";
		//out_planes.push_back(CutPlane::SetPlane(fvmath::CVec3d(1.0,0.0,0.0),center));
		//out_planes.push_back(CutPlane::SetPlane(fvmath::CVec3d(0.0,1.0,0.0),center));
		out_planes.push_back(CutPlane::SetPlane(fvmath::CVec3d(0.0,0.0,1.0),center));
	}

	return (int)out_planes.size();
}


CutPlane CutPlane::SetPlane(fvmath::CVec3d normal,fvmath::CVec3d point)
{
	CutPlane plane;
	plane.Plane::SetPlane(normal.v,point.v);
	return plane;
}

	unsigned int CutPlane::CutElement( /* return the number of elems */
			const CutPlane* pl,
			const double vts[18],
			const int nodes[7],
			const int index,
			int div,
			std::vector<Node_t>& vertices,
			std::vector<unsigned>& triangles
			)
	{
		assert(div >= 1);
		std::vector<myVertex> slice_vertices(5);
		int type = int(nodes[0] == 5);


		if (!pl->IsElement(elem_idx(type,index))) return 0;

		// Cut the element
		int num_tris = ::CutElement(vts,nodes,pl,slice_vertices);
		assert(num_tris >= 3);
		num_tris -= 2;
		//mfp_debug("::CutPlane produce: %d-triangular cross section\n",num_tris);
		unsigned int size = 0;
		myVertex *pa, *pb, *pc;
		switch(num_tris) {

		case 1: // All edges are drawn - 1 triangle
		{
			// One triangle
			pa = &slice_vertices[0]; pa->type = VertexType::EDGE_SHAPE;
			pb = &slice_vertices[1]; pb->type = VertexType::EDGE_SHAPE;
			pc = &slice_vertices[2]; pc->type = VertexType::EDGE_SHAPE;
			size = TessellateTriangle(pa,pb,pc,div,vertices,triangles);

		} break;

		case 2: // Leave one middle edge
		{
			// First triangle
			pa = &slice_vertices[0]; pa->type = VertexType::EDGE_SHAPE;
			pb = &slice_vertices[1]; pb->type = VertexType::EDGE_SHAPE;
			pc = &slice_vertices[2]; pc->type = VertexType::INSIDE_SHAPE;
			size = TessellateTriangle(pa,pb,pc,div,vertices,triangles);

			// Second triangle
			pa->type = VertexType::INSIDE_SHAPE;
			pb = &slice_vertices[3]; pb->type = VertexType::EDGE_SHAPE;
			pc->type = VertexType::EDGE_SHAPE;
			size += TessellateTriangle(pa,pc,pb,div,vertices,triangles);

		} break;

		case 3: // Leave two middle edges
		{
			// First triangle
			pa = &slice_vertices[0]; pa->type = VertexType::EDGE_SHAPE;
			pb = &slice_vertices[1]; pb->type = VertexType::EDGE_SHAPE;
			pc = &slice_vertices[2]; pc->type = VertexType::INSIDE_SHAPE;
			size = TessellateTriangle(pa,pb,pc,div,vertices,triangles);

			// Second triangle
			pa->type = VertexType::INSIDE_SHAPE;
			pb = &slice_vertices[3]; pb->type = VertexType::EDGE_SHAPE;
			pc->type = VertexType::EDGE_SHAPE;
			size += TessellateTriangle(pa,pc,pb,div,vertices,triangles);

			// Third triangle
			pb->type = VertexType::EDGE_SHAPE;
			pc = &slice_vertices[4]; pc->type = VertexType::EDGE_SHAPE;
			size += TessellateTriangle(pa,pc,pb,div,vertices,triangles);

		} break;

		default:
			assert("!Unexpected number of tris"); exit(-1);
		}
		// tessellate patches
		return size;
	}
    void CutPlane::AddIndex(const elem_idx&& idx) {
    	_elementIndices.insert(idx);
    }

    bool CutPlane::IsElement(const elem_idx&& idx) const {
    	return _elementIndices.find(idx) != _elementIndices.end();
    }

    unsigned CutPlane::GetMaxNUmVertices(int max_div)
    {
    	unsigned size = 0;
    	unsigned npts = (max_div+2)*(max_div+1)/2;
    	unsigned ntris;
    	for (const auto& it : _elementIndices) {
    		ntris = it.first ? 2U : 3U;
    		size += ntris*npts;
		}
    	return size;
    }

    unsigned CutPlane::GetMaxNumIndices(int max_div)
    {
    	unsigned size = 0;
    	unsigned nitems = (max_div+1)*max_div/2;
    	unsigned ntris;
    	for (const auto& it : _elementIndices) {
    		ntris = it.first ? 2U : 3U;
    		size += ntris*nitems;
		}
		return size;
    }



	bool CutPlane::operator ==(const CutPlane& rhs)
	{
		double diva = rhs.p.a != 0.0 ? p.a / rhs.p.a : FV_LARGE;
		double divb = rhs.p.b != 0.0 ? p.b / rhs.p.b : FV_LARGE;
		double divc = rhs.p.c != 0.0 ? p.c / rhs.p.c : FV_LARGE;
		double divd = rhs.p.d != 0.0 ? p.d / rhs.p.d : FV_LARGE;

		return ((diva == divb) && (diva == divc) && (diva == divd));
	}

//	int CutPlane::CheckElementLocation(const double x[], std::vector<FemViewer::GraphElement2<double> >& vels, int size)
//	{
//		static /*const*/ double tetra[16] = { x[0],x[1],  x[2],-1, x[3], x[4], x[5],-1, x[6], x[7], x[8],-1,x[9],x[10],x[11],-1};
//		static /*const*/ double prizm[24] = { x[0], x[1], x[2],-1, x[3], x[4], x[5],-1, x[6], x[7], x[8],-1,
//											  x[9],x[10],x[11],-1,x[12],x[13],x[14],-1,x[15],x[16],x[17],-1 };
//		static std::vector<FemViewer::GraphElement2<double> > _vels;
//		static GraphElement2<double> el;
//		if (!x) return (-1);
//
//		int result(0);
//		const double* table;
//		_vels.clear();
//
//		int nel(-1);
//		if (size == 4)
//		{
//			nel = DoSliceTerahedron(this,NULL,tetra,_vels,NULL);
//		}
//		else
//		{
//			nel = DoSlicePrizm(this,NULL,prizm,_vels,NULL);
//		}
//
//		return nel;
//
//	}

	int CutPlane::CheckOrientation(const double v[])
	{
		fvmath::CVec3d nrl(p.n);
		fvmath::CVec3d pos(v);
		fvmath::Normalize(pos);
		
		// calculate dot product
		double dot = fvmath::Dot(nrl,pos);

		if(fvmath::Compare(dot,0.0)) return(0); // perpendicular

		if(dot < 0.0f) return(-1); // the angle is > 90

		return(1); // the angle < 90

	}

	int CutPlane::InversPlane(const double v[])
	{
		fvmath::CVec3d nrl(p.n);
		fvmath::CVec3d pos(v);
		fvmath::Normalize(pos);

		// calculate dot product
		double dot = fvmath::Dot(nrl,pos);

		if(fvmath::Compare(dot,1.0)) 
			return(0); // do nothing - the same oritntation

		if(fvmath::Compare(dot,-1.0)) 
		{
			// Invers normal vector
			p.a = -p.a; p.b = -p.b; p.c = -p.c;
			return(1);
		}

		// Error
		return(-1);

	}


	bool CutPlane::InitOutlines(const BBox3D& bbox)
	{
		bool result = false;
		if (bbox.isInitialized())
		{
			mfp_debug("initOutlines\n");
			CVec3d vertices[8], center, beginVertex;
			double u;
			std::vector<CVec3d> plane_vertices;
			uint_t num = bbox.InitVertices(&bbox,vertices);
			for (uint_t ed = 0; ed < num; ++ed)
			{
				uint begin = edges[ed].begin();
				uint end   = edges[ed].end();
				beginVertex = vertices[begin];
				int res = IntersectWithLine(beginVertex.v,vertices[end].v,u);
				if (res && u >= 0.0 && u <= 1.0) {
					center += beginVertex;
					plane_vertices.push_back(beginVertex);
				}
			}

			if (!plane_vertices.empty()) {
				center *= 1.0/static_cast<double>(plane_vertices.size());
				//std::cout << "CENTRUM: " << center << std::endl;
				assert(fvmath::Compare(Plane::EvaluatePlane(center.x,center.y,center.z),0.0));
				// Make them local
				for_each(plane_vertices.begin(),plane_vertices.end(),
						[&](CVec3d& lhs) { lhs -= center; });
				// Sort in CCW order
				const CVec3d planeNormal = CVec3d(p.a,p.b,p.c);
				std::sort(plane_vertices.begin(),plane_vertices.end(),
						[&](const CVec3d& lhs,const CVec3d& rhs) -> bool {
							CVec3d v = (lhs* rhs);
							return fvmath::Dot(v,planeNormal) > 0;
						}
						);
				// Stroe in array of float coordinates
				//_vplaneCorrners.resize(plane_vertices.size()*2 + 1);
				// Store the center
				_vplaneCorrners.clear();
				_vplaneCorrners.push_back(center);
				for (size_t i=0;i<plane_vertices.size();++i) {
					_vplaneCorrners.push_back(CVec3f(plane_vertices[i]));
					_vplaneCorrners.push_back(CVec3f(planeNormal));
				}
				//for(auto & item : _vplaneCorrners)
				//	std::cout << item << std::endl;
				// Put into vector the center at the begining
				// Init scale
				_scale = CVec3f(1.2f,1.2f,1.2f);
				result = true;
			}

		}

		return result;
	}


	std::ostream& operator << (std::ostream& os, const CutPlane& rhs)
	{
		os << "Plane status: "
		   << (rhs.IsActive() ? "active |" : "inactive |")
		   << (rhs.IsVisible() ? " visible |" : " invisible |" )
		   << (rhs.IsNormalized() ? " normalized |" : " unnormalized |")
		   << (rhs.IsPlaneChanged() ? " changed" : " unchanged")
		   << std::endl
		   << "Plane parameters: "
		   << "A: " << rhs.p.a << " B: " << rhs.p.b << " C: " << rhs.p.c << " D: " << rhs.p.d
		   << std::endl;
		return os;

	}

}// end namespace
