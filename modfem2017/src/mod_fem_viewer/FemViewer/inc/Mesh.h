#ifndef _MESH_H_
#define _MESH_H_

#include "Log.h"
#include "fv_float.h"
#include "../../utils/fv_assert.h"
#include "../../include/fv_compiler.h"
#include "../../include/BaseMesh.h"
#include "mmh_intf.h"
#include "aph_intf.h"
#include "elem_tables.h"
#include "BBox3D.h"
#include "Color.h"
#include "ElemId.hpp"
#include "ArrayT.hpp"
#include "MathHelper.h"
#include "Enums.h"



#include<string>
#include <iostream>
#include<map>
#include<string.h>

#define MAX_FIELDS

namespace FemViewer {
	class ModelControler;
	class Object;
	class VtxAccumulator;
	class RenderParams;
	class CutPlane;
	class Accelerator;

	template<class T> class Ray;
	//struct IsectData;



	class  Mesh : public BaseMesh
	{
		friend class ModelControler;

	public:

		//typedef std::vector<int>			nodes_t;
		typedef ArrayT<id_t>				arBndElems;
		typedef std::vector<ElemInfo>		arElems;
		//typedef ArrayT<id_t>	arElems;
		typedef arBndElems::iterator		BndElemsIter;
		typedef arBndElems::const_iterator	BndElemsConstIter;
		typedef arBndElems::reference		BndElemsRef;
		typedef arBndElems::const_reference BndElemsConstRef;
		typedef arElems::iterator			arElemsIter;
		typedef arElems::const_iterator		arElConstIter;

		static const int tetra[4][3];// = { {0,1,2},{0,1,3},{0,2,3},{1,2,3} };
		static const int prizm[5][4];// = { {0,1,2,-1},{3,4,5,-1},{0,1,4,3},{1,2,5,4},{0,2,5,3} };
	public:
		// Constructors & destructor
		explicit Mesh();
		explicit Mesh(const std::string& name_);
		virtual ~Mesh();

		bool IsInit() const { return (!_arAllElements.empty()); }

	public: // Accessors
		      /*BBox3D&     GetMeshBBox3D() {
		    	  return (_bbox); }*/
		const BBox3D&     GetMeshBBox3D()    const { return( _bbox); }
		//const arBndElems& GetBoundaryElems() const { return _arElems; }

		template<typename T = double>
		int GetReferenceCoordinates(int nel,int type,int pdeg,int base,
			                        std::vector<fvmath::CVec3<T> >& coords,
									std::vector<short>& indices);
		template<typename T>
		int GetReferenceCoordinates(int nel,int type,int pdeg,int base,int side,
			                        std::vector<fvmath::CVec3<T> >& coords,
									std::vector<short>& indices);
		template<typename T = double>
		inline int GetElementCoordinates(int nel,int* Nodes,T* Coords) const {
			return BaseMesh::Get_el_node_coor(this->idx(),nel,Nodes,Coords);
		}

		template<typename T = double>
		int GetElemenCoordinatesFromNodes(const int* nodes, T* Coord) const {
			int num_nno = *nodes++;
			int index = 0;
			while (num_nno > 0) {
				GetNodeCoor(*nodes++,&Coord[index]);
				//printf("Get: %f %f %f\n",Coord[index],Coord[index+1],Coord[index+2]);
				--num_nno;
				index += 3;
			}
			return num_nno;
		}
		inline int GetElementFaces(int nel,int* Faces,int* Orient) {
			return BaseMesh::Get_el_faces(this->idx(),nel,Faces,Orient);
		}
		inline int GetNextActiveElement(int nel) {
			return BaseMesh::Get_next_act_el(this->idx(),nel);
		}
		inline int GetMaxNodeId(void) const {
			return BaseMesh::Get_nmno(this->idx());
		}


	public: // Other methods

		int   Init(const char* name);
		int   Create(const int type = eSelectionCategory::All);
		//int   Free();
		int	  Reload();
		int   Update();
		void  Reset();

		BBox3D   GetElemBBox3D (const int nodes[]) const {
			mfp_check_debug((nodes[0] == 6 || nodes[0] == 4), "Incorrect number of nodes: %d\n",nodes[0]);
			BBox3D bbox;
			Vec3d vrtx;
			for (int i = 1; i <= nodes[0]; ++i) {
				GetNodeCoor(nodes[i],vrtx.v);
				bbox += vrtx;
			}

			return bbox;
		}
        unsigned GetNumVerticesSelected(int div);
        unsigned GetNumIndicesSelected(int div);
		bool Render(unsigned int options,Object* hObj,CutPlane* plane=NULL,arBndElems* vbnds=NULL);

	private:

		//Object*		_pobject;
		eSelectionCategory m_selectType;
		//arBndElems  _arElems;
		arElems		_arAllElements;
		BBox3D 	    _bbox, _bbox_cut;
		size_type   _nno, _ned, _nfa, _nel;
		size_type	_lnel; // total number of nodes
		size_type   _ianno; // number of inactive nodes, info for VBO Creator
		Accelerator* _paccel;

	public:

		//Object*&  GetGLHandle() { return _pobject; }
		//int&      GetGLid() 	{ return _GLid; }
		size_type GetNumNodes() const { return _nno; }
		size_type GetNumEdges() const { return _ned; }
		size_type GetNumFaces() const { return _nfa; }
		size_type GetNumElems() const { return _nel; }

		arElems& GetElems() { return _arAllElements; }
		const arElems& GetElems() const { return _arAllElements; }

		size_type GetNumActiveNodes() const { return _nno - _ianno; }
		//nodes_t& GetNodes() { return _nodes; }

		BBox3D ExtractMeshBBoxFromNodes(bool do_parallel=false);

		int Select(std::vector<CutPlane>& vec_pcpls);

		int Select(const bool boundaryOnly = false);

		static uint32_t SetFaces(Mesh &rmesh);

		unsigned int RenderWireframe(
				std::vector<BaseVertex>& out_vertices,
				std::vector<unsigned>& out_edge_indices) const;

		unsigned int RenderWireframeCuts(const std::vector<CutPlane>& cut_planes,
				std::vector<Vertex>& outline_vertices,
				std::vector<Origin>& out_origins,
				std::vector<unsigned>& out_edge_indices,
				std::vector<int>& out_counts) const;

		//bool RenderWireframe2() const;

		bool RenderWireframeCutted(Object*  nObj,CutPlane* plane) const;

		int GetBoundaryVertices(std::vector<BaseVertex>& nodes);

		size_t GetElementOrigins(std::vector<Origin>& origins);



		int GetOriginsOfBoundaryFaces(std::vector<BaseVertex>& origins,
				std::vector<BaseVertex>& nodes);

		size_t GetVerticesOfCuts(const std::vector<CutPlane>& cut_planes,std::vector<BaseVertex>& vertices);
		template<typename T>
		inline int GetNodeCoor(int node,T* coords) const
		{
			return BaseMesh::Get_node_coor(this->idx(),node,coords);
	    }

		inline int GetEdgeNodes(int edge,int* Edge) const
		{
			return BaseMesh::Get_edge_nodes(this->idx(),edge,Edge);
		}

		inline int GetFaceoriginf(int ned,int nodes[],float Coords[]) const
		{
			return 0;
		}

		inline int GetElementCoordsf(int nel,int nodes[],float Coords[]) const
		{
			double lCoords[3*MMC_MAXELVNO];
			int nno = Get_el_node_coor(this->idx(),nel,nodes,lCoords);
			for (int i(0); i < nno; ++i) {
				Coords[3*i  ] = static_cast<float>(lCoords[3*i]);
				Coords[3*i+1] = static_cast<float>(lCoords[3*i+1]);
				Coords[3*i+2] = static_cast<float>(lCoords[3*i+2]);
			}
			return nno;
		}

		inline int GetElementOriginf(int nel,int nodes[],float Coords[]) const
		{
			double lCoords[3*MMC_MAXELVNO];
			double lOrigin[] = { 0.0, 0.0 ,0.0};
			//int nodes1[7];
			//printf("Mesh id: %d Element id: %d\n",this->id(),nel);
			int nno = Get_el_node_coor(this->idx(),nel,nodes,lCoords);
			//for(int i=0;i<3*MMC_MAXELVNO;i++) printf("%f\n",lCoords[i]);
			double inv_nno = 1.0 / (double)nno;
			for (int i(0); i < nno; ++i) {
				lOrigin[0] += lCoords[3*i];
				lOrigin[1] += lCoords[3*i+1];
				lOrigin[2] += lCoords[3*i+2];
				//printf("LCords: {%f, %f, %f}\n",lCoords[3*i],Coords[3*i+1],Coords[3*i+2]);
			}
			//printf*"Cords: {%f, %f, %f}\n",Coords[0],Coords[1],Coords[2]);
			Coords[0] = static_cast<float>(lOrigin[0]*inv_nno);
			Coords[1] = static_cast<float>(lOrigin[1]*inv_nno);
			Coords[2] = static_cast<float>(lOrigin[2]*inv_nno);
			//printf("Cords: {%f, %f, %f}\n",lOrigin[0],lOrigin[1],lOrigin[2]);
			return nno;
		}

		inline int GetElementCoordsd(int nel,int nodes [],double Coords[]) const
		{
			int nno = Get_el_node_coor(this->idx(),nel,nodes,Coords);
			return nno;
		}
		inline void GetFaceNeigInfo(int Fa,int Fa_neig[]) const {
			Get_face_neig(this->idx(),Fa,Fa_neig,NULL,NULL,NULL,NULL,NULL);
		}
		inline int GetElementFaces(int El,int Faces[],int Orient[]) const {
			return Get_el_faces(this->idx(),El,Faces,Orient);
		}

		inline int GetNodeStatus(int node) const
		{
			return BaseMesh::Get_node_status(this->idx(),node);
		}

		inline int GetEdgeStatus(int edge) const
		{
			return BaseMesh::Get_edge_status(this->idx(),edge);
		}
		inline int GetElemStatus(int nel) const {
			return BaseMesh::Get_el_status(this->idx(),nel);
		}

		template<typename T>
		inline int GetNextActiveElement(T nel = T(0)) const {
			return BaseMesh::Get_next_act_el(this->idx(),nel);
		}

		size_t PackElementCoord(std::vector<fvmath::Vec4<CoordType> >& ElemPrismCoords,
				bool onlyboundary = false);

		void   PackElementPlanes(std::vector<fvmath::Vec4<CoordType> >& ElemPlanes,
				bool onlyboundary = false);

		static void FillArrayWithElemCoords(const Mesh* pMesh,CoordType ElemCoords[]);

		static uint32_t FillArrayWithIndicesOfElemNodes(const Mesh *pMesh,int Indices[]);

		static void ConvertNodes(const Mesh *pMesh, float NodeCoords[]);
	private:
		// Not implemented
		Mesh(const Mesh&);
		Mesh& operator=(const Mesh&);
	};


template<>
inline int Mesh::GetNodeCoor<float>(int node,float* coords) const
{
	double lcoords[3];
	if(BaseMesh::Get_node_coor(this->idx(),node,lcoords) <0) return(-1);
	coords[0] = static_cast<float>(lcoords[0]);
	coords[1] = static_cast<float>(lcoords[1]);
	coords[2] = static_cast<float>(lcoords[2]);
	return(1);
}

template<typename T>
int Mesh::GetReferenceCoordinates(int nel,int type,int pdeg,int base_q,std::vector<fvmath::CVec3<T> >& coords,std::vector<short>& indices)
{
	// Specify a reference element
	const T * refCoords = (type == 4 ? &XlocTetra[0] : &XlocPrizm[0]);
	// Clear output containers
	coords.clear();
	indices.clear();
	// For liner or triangle elemenets
	if (base_q > APC_BASE_COMPLETE_DG || pdeg == 101) 
	{
		for (int i=0;i<type;i++) {
			int kk=3*i;
			//mfp_debug("type = %d refcords[%d] = { %lf, %lf, %lf}\n",type,i,refCoords[kk],refCoords[kk+1],refCoords[kk+2]);
			coords.push_back(fvmath::CVec3<T>((refCoords[kk]),
				                       (refCoords[kk+1]),
									   (refCoords[kk+2])));
			indices.push_back(i);
		}
		/*for (int i=0;i <type;++i) {
			//std::cout << coords[i] << std::endl;
			mfp_debug("type = %d refcords[%d] = { %lf, %lf, %lf}\n",type,i,coords[i].x,coords[i].y,coords[i].z);
		}*/
	} 
	else // Prism and high order
	{
		// Get the number of points and specify length
		int ptsx, ptsz;
		if (base_q == APC_BASE_TENSOR_DG) 
		{
			ptsx = pdeg % 100  + 1;
			ptsz = pdeg / 100  + 1;
		} 
		else // APC_BASE_COMPLETE_DG
		{
			ptsx = pdeg + 1;
	        ptsz = pdeg + 1;
		}

		T dl1 = T(1.0) / (T)ptsx;
		T dl2 = dl1;
		T dl3 = T(1.0) / (T)ptsz;

		// Loop over points
		//fvmath::Vec3<T> v;
		int ino = 0;
		for (int i=0; i < ptsx; i++) 
		{
			int ptsy = ptsx - i;
			for (int j=0; j < ptsy; j++)
			{
			 for (int k=0; k < ptsz; k++)
			 {
				 const T x = i*dl1;
			     const T y = j*dl2;
				 const T z = 2*k*dl3 - (T)1.0;
				 const T x1 = x + dl1;
				 const T y1 = y + dl2;
				 const T z1 = z + 2*dl3;
				 for (int itr=0; itr < 2; itr++)
				 {
					 // First triangle
					 if (itr==0 || j<(ptsy-1)) {
						if (itr==0) {
							// Store
							coords.push_back(fvmath::CVec3<T>(x,y,z));
							coords.push_back(fvmath::CVec3<T>(x1,y,z));
					 		coords.push_back(fvmath::CVec3<T>(x,y1,z));
					 		coords.push_back(fvmath::CVec3<T>(x,y,z1));
							coords.push_back(fvmath::CVec3<T>(x1,y,z1));
							coords.push_back(fvmath::CVec3<T>(x,y1,z1));
							for (int ii=0;ii<6;ii++) indices.push_back(ino+ii);
							ino+=6;
						} 
						else {
							// Second triangle
							coords.push_back(fvmath::CVec3<T>(x1,y,z));
							coords.push_back(fvmath::CVec3<T>(x1,y1,z));
							coords.push_back(fvmath::CVec3<T>(x,y1,z));
							coords.push_back(fvmath::CVec3<T>(x1,y,z1));
							coords.push_back(fvmath::CVec3<T>(x1,y1,z1));
							coords.push_back(fvmath::CVec3<T>(x,y1,z1));
							for (int ii=0;ii<6;ii++) indices.push_back(ino+ii);
							ino+=6;
						}// itr ==0
					 }
				 }// end itr
			 }// end k
			}// end j
		 }// end i
	}

	return (int)(indices.size());

}


template<typename T>
int Mesh::GetReferenceCoordinates(int nel,int type,int pdeg,int base_q,int side,std::vector<fvmath::CVec3<T> >& coords,std::vector<short>& indices)
{
	// Clear output containers
	coords.clear();
	indices.clear();
	// For liner or triangle elemenets
	if (base_q > APC_BASE_COMPLETE_DG || pdeg == 101) 
	{
		// Specify a reference element
		const T * refCoords = (type == 4 ? &XlocTetra[0] : &XlocPrizm[0]);
		for (int i=0;i<type;i++) {
			int kk=3*i;
			coords.push_back(fvmath::CVec3<T>(static_cast<T>(refCoords[kk++]),
				                       static_cast<T>(refCoords[kk++]),
									   static_cast<T>(refCoords[kk  ])));
			indices.push_back(i);
		}
	} 
	else // Prism and high order
	{
		// Get the number of points and specify length
		int ptsx, ptsz;
		if (base_q == APC_BASE_TENSOR_DG) 
		{
			ptsx = pdeg % 100  + 1;
			ptsz = pdeg / 100  + 1;
		} 
		else // APC_BASE_COMPLETE_DG
		{
			ptsx = pdeg + 1;
	        ptsz = pdeg + 1;
		}

		T dl1 = T(1.0) / (T)ptsx;
		T dl2 = dl1;
		T dl3 = T(1.0) / (T)ptsz;

		// Loop over points
		fvmath::Vec3<T> v;
		int ino = 0, k = 0;
		switch(side)
		{
			// Upper and lower base of a prizm
		case 1: k = ptsz - 1;
		case 0:
			{
				for (int i=0; i < ptsx; i++)
				{
					int ptsy = ptsx - i;
			        for (int j=0; j < ptsy; j++)
					{
						for (int itr=0; itr < 2; itr++)
						{
							const T x = i*dl1;
							const T y = j*dl2;
							const T z = 2*k*dl3 - (T)1.0;
							const T x1 = x + dl1;
							const T y1 = y + dl2;
							const T z1 = z + 2*dl3;
							// First triangle
							if (itr==0 || j<(ptsy-1)) {
								if (itr==0) {
									// Store
									coords.push_back(fvmath::CVec3<T>(x,y,z));
									coords.push_back(fvmath::CVec3<T>(x1,y,z));
					 				coords.push_back(fvmath::CVec3<T>(x,y1,z));
					 				coords.push_back(fvmath::CVec3<T>(x,y,z1));
									coords.push_back(fvmath::CVec3<T>(x1,y,z1));
									coords.push_back(fvmath::CVec3<T>(x,y1,z1));
									for (int ii=0;ii<6;ii++) indices.push_back(ino+ii);
									ino+=6;
								} 
								else {
									// Second triangle
									coords.push_back(fvmath::CVec3<T>(x1,y,z));
									coords.push_back(fvmath::CVec3<T>(x1,y1,z));
									coords.push_back(fvmath::CVec3<T>(x,y1,z));
									coords.push_back(fvmath::CVec3<T>(x1,y,z1));
									coords.push_back(fvmath::CVec3<T>(x1,y1,z1));
									coords.push_back(fvmath::CVec3<T>(x,y1,z1));
									for (int ii=0;ii<6;ii++) indices.push_back(ino+ii);
									ino+=6;
								}// itr ==0
							}
						}// end itr
					}// end for j
				}// end for i
			} break;
		case 2:
			{
				int j = 0;
				for (int i=0; i < ptsx; i++)
				{
			        for (int k=0; k < ptsz; k++)
					{
						const T x = i*dl1;
						const T y = j*dl2;
						const T z = 2*k*dl3 - (T)1.0;
						const T x1 = x + dl1;
						const T y1 = y + dl2;
						const T z1 = z + 2*dl3;
						coords.push_back(fvmath::CVec3<T>(x,y,z));
						coords.push_back(fvmath::CVec3<T>(x1,y,z));
					 	coords.push_back(fvmath::CVec3<T>(x,y1,z));
					 	coords.push_back(fvmath::CVec3<T>(x,y,z1));
						coords.push_back(fvmath::CVec3<T>(x1,y,z1));
						coords.push_back(fvmath::CVec3<T>(x,y1,z1));
						for (int ii=0;ii<6;ii++) indices.push_back(ino+ii);
						ino+=6;
					}// end for k
				}// end for i
			} break;
			case 5:
			{
				int i = 0;
				for (int j=0; i < ptsx; j++)
				{
			        for (int k=0; k < ptsz; k++)
					{
						const T x = i*dl1;
						const T y = j*dl2;
						const T z = 2*k*dl3 - (T)1.0;
						const T x1 = x + dl1;
						const T y1 = y + dl2;
						const T z1 = z + 2*dl3;
						coords.push_back(fvmath::CVec3<T>(x,y,z));
						coords.push_back(fvmath::CVec3<T>(x1,y,z));
					 	coords.push_back(fvmath::CVec3<T>(x,y1,z));
					 	coords.push_back(fvmath::CVec3<T>(x,y,z1));
						coords.push_back(fvmath::CVec3<T>(x1,y,z1));
						coords.push_back(fvmath::CVec3<T>(x,y1,z1));
						for (int ii=0;ii<6;ii++) indices.push_back(ino+ii);
						ino+=6;
					}// end for k
				}// end for i
			} break;
			case 4:
			{
				for (int i=0,j=ptsx; i < ptsx; i++,j-=i)
				{
			        for (int k=0; k < ptsz; k++)
					{
						const T x = i*dl1;
						const T y = j*dl2;
						const T z = 2*k*dl3 - (T)1.0;
						const T x1 = x + dl1;
						const T y1 = y + dl2;
						const T z1 = z + 2*dl3;
						coords.push_back(fvmath::CVec3<T>(x,y,z));
						coords.push_back(fvmath::CVec3<T>(x1,y,z));
					 	coords.push_back(fvmath::CVec3<T>(x,y1,z));
					 	coords.push_back(fvmath::CVec3<T>(x,y,z1));
						coords.push_back(fvmath::CVec3<T>(x1,y,z1));
						coords.push_back(fvmath::CVec3<T>(x,y1,z1));
						for (int ii=0;ii<6;ii++) indices.push_back(ino+ii);
						ino+=6;
					}// end for k
				}// end for i
			} break;
		}// switch		
	}// end else

	return (int)(indices.size());

}

template<>
inline int Mesh::GetElementCoordinates<float>(int nel,int* Nodes,float* Coords) const {
	const unsigned dim = 3*MMC_MAXELVNO;
	double lcoords[dim];
	int result = BaseMesh::Get_el_node_coor(this->idx(),nel,Nodes,lcoords);
	if (Coords) {
		for (unsigned i(0);i<dim;++i) {
			Coords[i] = static_cast<float>(lcoords[i]);
		}
	}
	return result;
}

} // end namespace

#endif /* _MESH_H_
*/
