#include "Mesh.h"
#include "Vec3D.h"
#include "Ray.h"
#include "Geometry.h"
#include "RenderParams.h"
//#include "Object.h"
#include "ViewManager.h"
#include "ModelControler.h"
#include "CutPlane.h"
#include "elem_tables.h"
#include "../../include/Enums.h"
#include "../../utils/fv_exception.h"
//#include "mmh_intf.h"
#include "defs.h"
#include "fv_config.h"
#include "fv_timer.h"
#include "Log.h"
#ifdef _USE_FV_EXT_MOD
#include "../../include/ApproxModule.h"
#endif
//#include "elem_tables.h"
#include "../../include/fv_compiler.h"
#include "../../utils/fv_assert.h"
#include "MathHelper.h"
//#include "VtxAccumulator.h"
//#include "VtxPrimitiveAccumulator.h"

// include system
#include <iostream>
#include <algorithm>
#include <set>
#include <omp.h>

#ifdef FV_DEBUG
#include<fstream>
#include<stdio.h>
#include"fv_conv_utils.h"
#endif


namespace FemViewer {

const int Mesh::tetra[4][3] = { {0,2,1},{0,1,3},{0,3,2},{1,2,3} };
const int Mesh::prizm[5][4] = { {0,2,1,-1},{3,4,5,-1},{0,1,4,3},{1,2,5,4},{2,0,3,5} };


static
fvmath::CVec3d baryCenter(myVertex vertices[],int num)
{
	fvmath::CVec3d center;
	double inv_num = 1.0 / (double)num;

	for (int iv=0; iv<num; ++iv)
		center += vertices[iv].position;

	center *= inv_num;

	return center;
}



Mesh::Mesh()
: BaseMesh("Unknown")
, m_selectType(All)
, _arAllElements()
, _bbox()
, _bbox_cut()
, _nno(0)
, _ned(0)
, _nfa(0)
, _nel(0)
, _lnel(0)
, _ianno(0)
, _paccel(nullptr)
{
	//mfp_log_debug("Mesh ctr with unknown name\n");
}


Mesh::Mesh(const std::string& name_)
: BaseMesh(name_)
, m_selectType(All)
, _arAllElements()
, _bbox()
, _bbox_cut()
, _nno(0)
, _ned(0)
, _nfa(0)
, _nel(0)
, _lnel(0)
, _ianno(0)
, _paccel(nullptr)
{
	//mfp_log_debug("Mesh ctr\n");
	if (Init(name_.c_str()) < 0) throw fv_exception("Can't init mesh data!");
}

Mesh::~Mesh()
{
	//mfp_log_debug("Mesh dtr\n");
	Reset();
}

int Mesh::Init(const char* fname_path)
{
	if (BaseMesh::Init(fname_path) < 0) return(-1);

	// Initialize basic parameters
	_nno = BaseMesh::Get_nmno(this->idx());
	_ned = BaseMesh::Get_nmed(this->idx());
	_nfa = BaseMesh::Get_nmfa(this->idx()); 
	_nel = BaseMesh::Get_nmel(this->idx());

	// Initialize bounding box info
	_bbox = ExtractMeshBBoxFromNodes();

	return 1;
}

int Mesh::Create(const int type)
{
	if (m_selectType != type && All <= type && type <= Cutted) {
		_arAllElements.clear();
		m_selectType = static_cast<eSelectionCategory>(type);
	}

	int num_selected = static_cast<int>(_arAllElements.size());
	if (_arAllElements.empty())
	{
		bool state = false;
		// Select appropriate elements according to given (type) policy
		switch (m_selectType)
		{
		case eSelectionCategory::Boundary: state = true;
		case eSelectionCategory::All:
			num_selected = Select(state);
			break;

		case eSelectionCategory::Cutted:
		{
			mfp_debug("Init cut elements\n");
			if (ModelCtrlInst().GetCutPlanes().empty())
				throw fv_exception("Empty set of cut-planes!");
			else
			{
				num_selected = Select(ModelCtrlInst().GetCutPlanes());
				for (unsigned int i = 0; i<ModelCtrlInst().GetCutPlanes().size(); ++i) {
					// Init outlines for drawing
					assert(_bbox.isInitialized());
					ModelCtrlInst().GetCutPlanes()[i].InitOutlines(_bbox);
				}
			}
		} break;

		default:
			assert(!"Unknow selection type");
			throw fv_exception("Unknown selcetion type in mesh module!");
		}

	}

	return num_selected;

}



//int Mesh::Free()
//{
//	Reset();
//	return this->BaseMesh::Free();
//}

int Mesh::Reload()
{
	Reset();
	return Init(this->_name.c_str());
}

int Mesh::Update()
{
	return(-1);
}

void Mesh::Reset()
{
	//mfp_log_debug("Reset\n");
	BaseMesh::Reset();
	m_selectType = All;
	_arAllElements.clear();
	_bbox.Reset();
	_bbox_cut.Reset();
	_nno = _ned = _nfa = _nel = 0;

	_lnel = 0;
}


BBox3D Mesh::ExtractMeshBBoxFromNodes(bool do_parallel)
{
	using fvmath::Vec3d;

	BBox3D mesh_bbox;
	int num_inact = 0;

	int num_nodes = GetNumNodes();
	fv_timer t;
	t.start();

	#pragma omp parallel if(do_parallel) default(none) \
		firstprivate(num_nodes) shared(mesh_bbox,num_inact)
	{
		BBox3D lbbox;
		Vec3d lcoord;
		int lianno = 0;

		#pragma omp for schedule(dynamic,4)
		for (int i=1; i <= num_nodes; ++i)
		{
			if (GetNodeCoor(i,lcoord.v)) {
				lbbox += lcoord;
			}
			else {
				++lianno;
			}
		}

		#pragma omp critical(bbox)
		{
			mesh_bbox += lbbox;
			num_inact += lianno;
		}
	}
	t.stop();
	_ianno = num_inact;
	assert(mesh_bbox.isInitialized());
	mfp_log_info("Building%bounding box info in %lf sec",(do_parallel ? " in parallel " : " "), t.get_time_sec() );
	mfp_log_info("Found %d number of inactive nodes", num_inact);
	return mesh_bbox;
}

// Fill in a vector of graphic elements
// from boundary faces of elements

int Mesh::Select(std::vector<CutPlane>& vec_planes)
{
	mfp_log_debug("Select elemensts with %u plane(s)\n",vec_planes.size());
	typedef std::vector<CutPlane>::iterator cpIter;
	const cpIter ed(vec_planes.end());

	int curr_index = 0;
	int non_selected = 0;
	int faces[MMC_MAXELFAC+1];
	int Fa_neig[2];
	Vec3d v[MMC_MAXELVNO];

	ElemInfo el_info;
	const int id = this->Id();

	assert(_nel > 0);
	_arAllElements.reserve(_nel);
    _arAllElements.clear();
	fv_timer t;
	t.start();
	int nel = 0;
	while ((nel=Get_next_act_el(id, nel)) > 0)
	{
		// Get element physical coordinates
		int nno = Get_el_node_coor(id,nel,el_info.nodes,v[0].v);
		int selected = 0;
#if 0
		// I don't know but this produce not exacly elements
		// Calculate bounding box
		AAbbd bb;
		if (el_info.nodes[0] == 5) {
			// Prism
			bb.mn.x = fv_min(fv_min(fv_min(fv_min(fv_min(v[0].x,v[1].x),v[2].x),v[3].x),v[4].x),v[5].x);
			bb.mn.y = fv_min(fv_min(fv_min(fv_min(fv_min(v[0].y,v[1].y),v[2].y),v[3].y),v[4].y),v[5].y);
			bb.mn.z = fv_min(fv_min(fv_min(fv_min(fv_min(v[0].z,v[1].z),v[2].z),v[3].z),v[4].z),v[5].z);
			bb.mx.x = fv_max(fv_max(fv_max(fv_max(fv_max(v[0].x,v[1].x),v[2].x),v[3].x),v[4].x),v[5].x);
			bb.mx.y = fv_max(fv_max(fv_max(fv_max(fv_max(v[0].y,v[1].y),v[2].y),v[3].y),v[4].y),v[5].y);
			bb.mx.z = fv_max(fv_max(fv_max(fv_max(fv_max(v[0].z,v[1].z),v[2].z),v[3].z),v[4].z),v[5].z);
		}
		else {
			// Terahedron
			bb.mn.x = fv_min(fv_min(fv_min(fv_min(v[0].x,v[1].x),v[2].x),v[3].x),v[4].x);
			bb.mn.y = fv_min(fv_min(fv_min(fv_min(v[0].y,v[1].y),v[2].y),v[3].y),v[4].y);
			bb.mn.z = fv_min(fv_min(fv_min(fv_min(v[0].z,v[1].z),v[2].z),v[3].z),v[4].z);
			bb.mx.x = fv_max(fv_max(fv_max(fv_max(v[0].x,v[1].x),v[2].x),v[3].x),v[4].x);
			bb.mx.y = fv_max(fv_max(fv_max(fv_max(v[0].y,v[1].y),v[2].y),v[3].y),v[4].y);
			bb.mx.z = fv_max(fv_max(fv_max(fv_max(v[0].z,v[1].z),v[2].z),v[3].z),v[4].z);
		}
		// Little extend
		bb.mn.x -= 0.01;
		bb.mn.y -= 0.01;
		bb.mn.z -= 0.01;
		bb.mx.x += 0.01;
		bb.mx.y += 0.01;
		bb.mn.z += 0.01;

		// First check for bbox intersection with plane
		for (size_t i = 0; i < vec_planes.size(); ++i) {
			bool val = vec_planes[i].IntersectBBox3D(bb) != 0;
			selected += static_cast<int>(val);
		}

		// Skip non-selected elements
		if (selected) { non_selected++; continue; }
#endif
		// Select by marching-prism/tetra algorithm
		for (cpIter bg(vec_planes.begin()); bg != ed; ++bg)
		{
			int type;
			bool res;
			if (el_info.nodes[0] == 6) {
				type = 1;
				res = (IsPrismSelected(v[0].v,bg->GetParams()) >= 0);
			}
			else {
				type = 0;
				res = (IsTetraSelected(v[0].v,bg->GetParams()) >= 0);
			}

			if (res) {
				bg->AddIndex(CutPlane::elem_idx(type,curr_index));
				selected++;
			}
		}
		// Reset id
		el_info.id = 0;
		// If the element was selected by any plane
		// denote it as active
		if (selected) {
			if (*faces == 5) { SET_PRISM(el_info.id); };
			el_info.id |= fvmath::abs(nel);
			_arAllElements.push_back(el_info);
			curr_index++;
			//mfp_debug("Element%u has ",_arAllElements.size());
			//for (int i=0,j=1;i<el_info.nodes[0];i++,j++)
			//	mfp_debug("%d ",el_info.nodes[j]);
		}
		// insert element into vector
	}
	t.stop();
	double stop = t.get_duration_us();//t.test_time(1);
	mfp_log_info("Selected %u elementss from %u total\nTime: %lf us",_arAllElements.size(),_nel,stop);
	mfp_log_info("Unselected elements: %u",non_selected);
	return static_cast<int>(_arAllElements.size());
}


int Mesh::Select(const bool boundaryOnly)
{
	//mfp_log_debug("SetElems\n");
	int nel=0,cnt,Fa;
	int faces[MMC_MAXELFAC+1];
	int Fa_neig[2];
	ElemInfo el_info;

	_arAllElements.reserve(_nel);
	_arAllElements.clear();
	uint32_t totalNumOfNodes = 0;
	//_arAllElements.clear();
	fv_timer t;
	t.start();

	while ((nel=Get_next_act_el(this->idx(), nel)) > 0)
	{
		// Get element physical coordinates
		int nno = Get_el_node_coor(this->idx(),nel,el_info.nodes,NULL);

		// Get the ids of faces
		Get_el_faces(this->idx(),nel,faces,NULL);
		assert((*faces == 5) || (*faces == 4));

		// Resets the counter of boundary faces for the elment
		cnt = 0;
		el_info.id = 0;
		for (int j = 1; j <= *faces; ++j) {
			Fa = fvmath::abs(faces[j]);
			Get_face_neig(this->idx(),Fa,Fa_neig,NULL,NULL,NULL,NULL,NULL);
			if (Fa_neig[1] == MMC_BOUNDARY) {
				//printf("mesh_select: nr_i = %d, nr_Fa = %d, bound = %d\n",j,Fa,Fa_neig[1]);
				SET_STFACEN(el_info.id,(j-1));
				++cnt;
			}
		}// end for

		// Skip element if not boundary
		if (!cnt && boundaryOnly) continue;

		// Set boundary flag
		if (cnt) SET_BOUND(el_info.id);

		// Set type of an element
		if (*faces == 5) SET_PRISM(el_info.id);

		// Insert element into vector
		el_info.id |= fvmath::abs(nel);

		_arAllElements.push_back(el_info);
		//totalNumOfNodes += nno;
	}//end while

	//mfp_log_debug("Toatal number of nodes = %d\n",rmesh._arAllElements.size());
	double stop = t.get_duration_sec();
	mfp_log_info("Selected %u elementss form %d total\nTime: %lf sec",_arAllElements.size(),_nel,stop);
	return static_cast<int>(_arAllElements.size());
}

unsigned Mesh::GetNumVerticesSelected(int div)
{
	const int ntris_prism = 3*(div + 1)*(div + 2)/2;
	const int ntris_tetra = 2*(div + 1)*(div + 2)/2;
	unsigned num_vrts = 0;

	for (unsigned int iel = 0; iel < _arAllElements.size(); ++iel) {
//		if (_arAllElements[iel].is_prism()) num_vrts  += ntris_prism;
//		else if (_arAllElements[iel].is_tetra()) num_vrts  += ntris_tetra;
		num_vrts += _arAllElements[iel].is_prism() ? ntris_prism : ntris_tetra;
	}

	return num_vrts;
}

unsigned Mesh::GetNumIndicesSelected(int div)
{
	const int ntris_prism = 3*(div)*(div + 1)/2;
	const int ntris_tetra = 2*(div)*(div + 1)/2;
	unsigned num_idxs = 0;

	for (unsigned int iel = 0; iel < _arAllElements.size(); ++iel) {
//		if (_arAllElements[iel].is_prism()) num_idxs  += ntris_prism;
//		else if (_arAllElements[iel].is_tetra()) num_idxs  += ntris_tetra;
		num_idxs += _arAllElements[iel].is_prism() ? ntris_prism : ntris_tetra;
	}

	return num_idxs;
}


bool Mesh::Render(const uint_t options,Object* hObj,CutPlane* plane,arBndElems* vbnd)
{
//	FV_ASSERT(hObj != NULL);
//	FV_ASSERT(_nno > 0);
//	std::vector<CVec3f> lVertices(_nno);
//	CVec3f lcoord;
//	_ianno = 0;
//	// Convert vertices to its floats coordinates
//	for (int i=1; i <= _nno; ++i) {
//		if ( MMC_ACTIVE == Get_node_status(_idx,i)) {
//			GetNodeCoor(this->idx(),lcoord.v);
//			lVertices.push_back(lcoord);
//		} else {
//			++_ianno;
//		}
//	}
//	hObj->GetCurrentVtxAccumulator();
//
//	return true;
}

unsigned int Mesh::RenderWireframe(
	std::vector<BaseVertex>& out_vertices, // original mesh vertices
	std::vector<unsigned>& out_edge_indices // indices of edges
	) const
{

	// Loop over nodes
	out_vertices.reserve(_nno);
	out_vertices.clear();
	BaseVertex vertex;

	std::vector<int> local_index(_nno);
	// Store vertices of the mesh
	for (int inode=1,lindex=0; inode<=_nno; ++inode)
	{
		if (MMC_ACTIVE == GetNodeStatus(inode)) {
			GetNodeCoor(inode,vertex.position.v);
			vertex.info = inode;
			out_vertices.push_back(vertex);
			local_index[lindex] = inode;
			++lindex;
		}
	}

	// Loop over edges
	out_edge_indices.reserve(2*_ned);
	out_edge_indices.clear();

	auto const p = local_index.data();
	unsigned int ii[2];
	int edge[2], next;
	for (int iedge=1; iedge<=_ned;++iedge)
	{
		if (MMC_ACTIVE == GetEdgeStatus(iedge) ) {
			GetEdgeNodes(iedge,edge);
			assert(edge[0] > 0 || edge[0] <= out_vertices.size());
			assert(edge[1] > 0 || edge[1] <= out_vertices.size());
			// Swap for increase of performance
			if (edge[0] > edge[1]) {
				int tmp = edge[0];
				edge[0] = edge[1];
				edge[1] = tmp;
			}
			// We have to eliminate inactive nodes if there are any
			ii[1]=ii[0]=0;
			next = 0;
			do {
				if (p[ii[next] ] != edge[next]) {
					ii[next]++;
					next = 1 - next;
			    }
				else break;
			} while(1);
			next = 1 - next;
			do {
				if (p[ii[next]] == edge[next]) break;
			    ii[next]++;
			} while(1);
			// Store edge vertices according to local indices
			out_edge_indices.push_back(ii[0]);
			out_edge_indices.push_back(ii[1]);
		}
	}


	return out_edge_indices.size();//vertex_counter;
}

#define MAX_VERT_PER_ELEM 8
//
unsigned int Mesh::RenderWireframeCuts(const std::vector<CutPlane>& cut_planes,
		std::vector<Vertex>& outline_vertices,
		std::vector<Origin>& out_origins,
		std::vector<unsigned>& out_edge_indices,
		std::vector<int>& out_counts) const
{
  double Coords[3*MMC_MAXELVNO];
  std::vector<myVertex> cutted_vertices(MAX_VERT_PER_ELEM);
  const unsigned int num_corners = _arAllElements.size()*MAX_VERT_PER_ELEM*cut_planes.size() / 2;
  outline_vertices.reserve(num_corners);
  out_edge_indices.reserve(num_corners);
  out_counts.reserve(_arAllElements.size());
  Origin center;
  int nel = 0;
  std::pair<int,int> curr_el;

  assert(!_arAllElements.empty());

  for (const auto & el : _arAllElements)
  {
    // Get the index of an element in mesh
    center.info = el.eid;
      int type = (el.nodes[0] == 6) ? 1 : 0;
      for (const auto & pl : cut_planes)
	{
    	 // mfp_debug("At the begining of loop\n");
	  // find element in plane's elements
	  if (pl.IsElement(CutPlane::elem_idx(type,nel)))
	    {
	      // Get coordinates of an element
		  //for (int i=0,j=1;i<el.nodes[0];i++,j++)
		  //{
		//	  printf("%d ",el.nodes[j]);
		 // }
		  //printf("\n");
	      GetElemenCoordinatesFromNodes(el.nodes,Coords);

	      // Cut an element
	      int num_vertices = static_cast<int>(::CutElement(Coords,el.nodes,&pl,cutted_vertices));
	      //printf("Cutted vertices:\n");
	      //for (const auto& vt : cutted_vertices)
	      //{
	    	//  std::cout << vt.position << std::endl;
	      //}

	      //out_edge_indices.push_back(num_vertices);
	      out_counts.push_back(num_vertices);

	      // Convert and store vertices in out containers
	      const unsigned int base = out_edge_indices.size();
	      unsigned int count  = 0;
	      for (const auto& v : cutted_vertices)
	      {
	    	  outline_vertices.push_back(Vertex(v));
	    	  out_edge_indices.push_back(base + count++);
	      }

	      // Get centrum of a slice
	      center.position = baryCenter(cutted_vertices.data(),num_vertices);
	      //std::cout << "Edge out size: " << out_edge_indices.size() << std::endl;
	      //std::cout << "Calculated barycenter: " << center.position << std::endl;
	      out_origins.push_back(center);
	      
	      // Clear vector
	      cutted_vertices.clear();

	      //mfp_debug("here is OK\n");
	    }
	 }
  }

  // We have to remove duplicates in outline_vertices
  //mfp_debug("Before removing duplicatesz: outlinevertices %u edge_indices %u\n",outline_vertices.size(),out_edge_indices.size());
  //for (const auto& vt : outline_vertices)
//	  std::cout << vt.position << std::endl;
  //for (const auto& vt : out_edge_indices)
//	  std::cout << "edge = " << vt << std::endl;
  //exit(-1);
  ::indexVBO(outline_vertices,out_edge_indices);
  //std::cout << "After\n";
  //for (const auto& vt : outline_vertices)
  //	  std::cout << vt.position << std::endl;
   // for (const auto& vt : out_edge_indices)
  	 // std::cout << "edge = " << vt << std::endl;

  return out_counts.size();
}

size_t Mesh::GetElementOrigins(std::vector<Origin>& origins)
{
	origins.reserve(_nel);
	origins.clear();

	Origin center;
	int nel = 0;
	while ((nel=Get_next_act_el(this->idx(), nel)) > 0)
	{
		int nno = GetElementOriginf(nel,NULL,center.position.v);
		center.info = nel;
		//std::cout << "The centerof an element " << nel << " with " << nno <<" vertices is: " << center.position << std::endl;
		origins.push_back(center);
	}

	return origins.size();
}

#undef MAX_VERT_PER_ELEM



int Mesh::GetBoundaryVertices(std::vector<BaseVertex>& out_nodes)
{
	BaseVertex node;
	arElConstIter it = _arAllElements.begin();
	const arElConstIter it_e = _arAllElements.end();
	const int * elem_p = nullptr;
	std::set<int> local_indices;

	if (Create(Boundary)==0) {
		mfp_log_err("No element selected!");
		return 0;
	}

	while (it != it_e) {

		if (IS_BOUND(it->id))
		{
			int nr_faces(4);
			if (it->is_tetra()) {
				elem_p = Mesh::tetra[0];
			}
			else {
				nr_faces++;
				elem_p = Mesh::prizm[0];
			}

			//local_indices.clear();
			for (int nfa(0); nfa < nr_faces; ++nfa)
			{
				// Skip a non boundary face
				if (!IS_SET_FACEN(it->id, nfa)) continue;
				//printf("id = %u nf_faces = %d nfa = %d\n",it->id,nr_faces,nfa);
				int of = nfa*(nr_faces-1);
				// Get offsets in table of nodes
				int ino = elem_p[of++]+1;
				local_indices.insert(it->nodes[ino]);
				ino = elem_p[of++]+1;
				local_indices.insert(it->nodes[ino]);
				ino = elem_p[of++]+1;
				local_indices.insert(it->nodes[ino]);
				if (nr_faces == 5 && nfa >= 2) {
					ino = elem_p[of]+1;
					local_indices.insert(it->nodes[ino]);
				}
			}
		}
		++it;
	}

	assert(local_indices.empty()==false);
	out_nodes.reserve(local_indices.size());
	out_nodes.clear();

	for (auto it = local_indices.begin(), ite = local_indices.end(); it != ite; ++it)
	{
		GetNodeCoor(*it, node.position.v);
		node.info = *it;
		out_nodes.push_back(node);
	}

	return static_cast<int>(out_nodes.size());
}

int Mesh::GetOriginsOfBoundaryFaces(
		std::vector<BaseVertex>& origins,
		std::vector<BaseVertex>& nodes
		)
{
	// Select bounday elements
	Create(Boundary);
	// Reserve necessary memory space
	origins.reserve(_arAllElements.size());
	origins.clear();

	BaseVertex origin;
	arElConstIter it = _arAllElements.begin();
	const arElConstIter it_e = _arAllElements.end();
	const int * elem_p = nullptr;

	while (it != it_e)
	{
		if (IS_BOUND(it->id))
		{
			int nr_faces(4);
			if (it->is_tetra()) {
				elem_p = Mesh::tetra[0];
			}
			else {
				nr_faces++;
				elem_p = Mesh::prizm[0];
			}

			for (int nfa(0); nfa < nr_faces; ++nfa)
			{
				// Skip a non boundary face
				if (!IS_SET_FACEN(it->id, nfa)) continue;
				//printf("id = %u nf_faces = %d nfa = %d\n",it->id,nr_faces,nfa);
				int of = nfa*(nr_faces-1);
				// Get offsets in table of nodes
				int ino = elem_p[of++]+1;
				fvmath::CVec3d orig;
				GetNodeCoor(it->nodes[ino],orig.v);
				//std::cout <<"first: " << origin;
				ino = elem_p[of++]+1;
				fvmath::CVec3d node;
				GetNodeCoor(it->nodes[ino],node.v);
				//std::cout <<"second: " << node;
				orig += node;
				ino = elem_p[of++]+1;
				GetNodeCoor(it->nodes[ino],node.v);
				//std::cout <<"third: " << node;
				orig += node;
				double inv_div = 0.333333;
				if (nr_faces == 5 && nfa >= 2) {
					ino = elem_p[of]+1;
					GetNodeCoor(it->nodes[ino],node.v);
					//std::cout <<"forth: " << origin;
					orig += node;
					inv_div = 0.25;
				}
				orig *= inv_div;
				origin.position = orig;
				origin.info = int(it->eid);
				origins.push_back(origin);
			}
		}
		++it;
	}

	return GetBoundaryVertices(nodes);
}
#define MAX_VERT_PER_ELEM 8
size_t Mesh::GetVerticesOfCuts(const std::vector<CutPlane>& cut_planes,std::vector<BaseVertex>& vertices)
{
	double Coords[3*MMC_MAXELVNO];
	std::vector<myVertex> cutted_vertices(MAX_VERT_PER_ELEM);
	assert(!_arAllElements.empty());
	const size_t num_corners = _arAllElements.size()*MAX_VERT_PER_ELEM*cut_planes.size() / 2;
	vertices.reserve(num_corners + _nel);
	vertices.clear();

	Origin cut_node;
	int nel = 0;
	std::pair<int,int> curr_el;


	std::set<int> unique_vindices;
	for (const auto & el : _arAllElements)
	{
		// Get the index of an element in mesh
		cut_node.info = el.eid;
		int type = (el.nodes[0] == 6) ? 1 : 0;
		// Loop over each slice in element
		for (const auto & pl : cut_planes)
		{
			// Find element in plane's elements
			if (pl.IsElement(CutPlane::elem_idx(type,nel)))
		    {
		      // Get coordinates of an element
		      GetElemenCoordinatesFromNodes(el.nodes,Coords);

		      // Cut an element
		      int num_vertices = static_cast<int>(::CutElement(Coords,el.nodes,&pl,cutted_vertices));

		      // Store only unique global ids
		      for (const auto& v : cutted_vertices)
		      {
		    	  if (v.info) unique_vindices.insert(v.info);
		      }

		      // Get centrum of a slice
		      cut_node.position = baryCenter(cutted_vertices.data(),num_vertices);
		      vertices.push_back(cut_node);

		      // Clear vector
		      cutted_vertices.clear();
		    }
		}
	}

	// Store the number of roigins
	size_t num_elem_origs = vertices.size();

	// Get vertes coordinates
	for (const auto & id: unique_vindices)
	{
		GetNodeCoor(id, cut_node.position.v);
		cut_node.info = id;
		vertices.push_back(cut_node);
	}

	return num_elem_origs;
}
#undef MAX_VERT_PER_ELEM

uint32_t Mesh::FillArrayWithIndicesOfElemNodes(const Mesh *pMesh,int Indices[])
{
	assert(Indices != NULL);
	Mesh::arElConstIter itb(pMesh->GetElems().begin());
	const Mesh::arElConstIter ite(pMesh->GetElems().end());
	uint32_t index(0);
	for (uint32_t numNodes(0);itb != ite; ++itb, index += (numNodes+1))
	{
		numNodes = static_cast<uint32_t>(itb->nodes[0]);
		Indices[0] = -numNodes;
		for (uint32_t i(1);i <= numNodes; ++i) {
			Indices[index + i] = static_cast<uint32_t>(itb->nodes[i]);
		}
	}

	return index;
}


size_t Mesh::PackElementCoord(std::vector<fvmath::Vec4<CoordType> >& ElemPrismCoords,
		bool onlyboundary)
{
	// Check if mesh data is present
	if (_arAllElements.empty()) Select(onlyboundary);

	// Adjust size of a container to the number of elements
	fvmath::Vec4<CoordType> node_coords = {0, 0 , 0, 0};
	size_t num_elems = _arAllElements.size();
	assert(num_elems > 0);
	ElemPrismCoords.resize(6 * num_elems);

	// Start proceeding
	for (size_t nel = 0; nel < num_elems; ++nel) {
		ElemInfo* info_ptr = &_arAllElements[nel];
		assert(info_ptr->nodes[0] == 6);
		// Pack nodes one by one up to 6 of prism
		for (int i = 0; i < info_ptr->nodes[0]; ++i) {
			if (GetNodeCoor<CoordType>(info_ptr->nodes[i+1],node_coords.v) > 0)
				ElemPrismCoords.push_back(node_coords);
		}
	}

	return num_elems;
}

void Mesh::PackElementPlanes(std::vector<fvmath::Vec4<CoordType> >& ElemPlanes,
		bool onlyboundary)
{
	// Check if mesh data is present
	if (_arAllElements.empty()) Select(onlyboundary);

	// Adjust size of a container to the number of elements
	double node_coords[MMC_MAXELVNO];
	size_t num_elems = _arAllElements.size();
	assert(num_elems > 0);
	ElemPlanes.resize(num_elems);

	// Start proceeding
	for (size_t nel = 0; nel < num_elems; ++nel) {
		ElemInfo* info_ptr = &_arAllElements[nel];
		assert(info_ptr->nodes[0] == 6);
		// Pack nodes one by one up to 6 of prism in local buffer
		for (int i = 0; i < info_ptr->nodes[0]; ++i) {
			GetNodeCoor<double>(info_ptr->nodes[i+1],&node_coords[3*i]);
		}
		// Extract planes for faces
		if (info_ptr->nodes[0] == 6) Plane::CreatePlanesForPrism((Plane::vec3 *)&node_coords[0],&ElemPlanes[nel]);
		else if (info_ptr->nodes[0] == 4) Plane::CreatePlanesForTetra((Plane::vec3 *)&node_coords[0],&ElemPlanes[nel]);
	}
}

void Mesh::FillArrayWithElemCoords(const Mesh* pMesh,CoordType ElemCoords[])
{
	mfp_log_debug("Fill array with coords of element nodes");
	assert(ElemCoords != NULL);
	Mesh::arElConstIter itb(pMesh->GetElems().begin());
	Mesh::arElConstIter ite(pMesh->GetElems().begin());
	for (uint32_t index(0);itb != ite; ++itb, ++index)
	{
		for (int i(0);i < itb->nodes[0]; ++i) {
			pMesh->GetNodeCoor(itb->nodes[i],&ElemCoords[3*index]);
		}
	}
}

void Mesh::ConvertNodes(const Mesh *pMesh, float Nodes[])
{
	for (size_type node(1); node < pMesh->GetNumNodes(); ++node)
		pMesh->GetNodeCoor((int)node, &Nodes[node]);
}



}// end namespace FemViewer

