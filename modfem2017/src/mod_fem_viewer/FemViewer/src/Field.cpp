#include "fv_config.h"
#include "Geometry.h"
#include "Field.h"
#include "Object.h"
#include "Fire.h"

#include "aph_intf.h"
#include "mmh_intf.h"
#include "pdh_intf.h"
#include "defs.h"
#include "Enums.h"
#include "Log.h"
#include "Matrix.h"
#include "fv_compiler.h"
#include "fv_float.h"
#include "../../include/Enums.h"
#include "../../utils/fv_exception.h"
#include "../../utils/fv_mathparser.h"
#include "elem_tables.h"
#include "ViewManager.h"
#include "VtxAccumulator.h"
#include<cstdio>
#include<omp.h>
#include<iterator>
#include<algorithm>

#ifdef FV_DEBUG
#include <iostream>
using namespace std;
#endif

namespace FemViewer {

/*typedef struct _Node_t {
  fvmath::CVec3d position;
  int  info;
} Node_t;*/

//struct AprElemCalc3D {
//
//};

static inline
int maximal_pdeg(int pdeg)
{
	int zdeg  = abs(pdeg / 100);
	int xydeg = abs(pdeg % 100);
	int ydeg  = xydeg / 10;
	int xdeg  = xydeg % 10;

	return fv_max(fv_max(xdeg,ydeg),zdeg);
}

static inline
int minimal_pdeg(int pdeg)
{
	int zdeg  = abs(pdeg / 100);
	int xydeg = abs(pdeg % 100);
	int ydeg  = xydeg / 10;
	int xdeg  = xydeg % 10;

	return fv_min(fv_min(xdeg,ydeg),zdeg);
}

// returns the number of verices in unit triangle
// for given cout of division
static int num_pts_in_triangle(
		int ndiv	// The number of uniform division
		)
{
	return (ndiv+1)*(ndiv+2)/2;
}

// It returns the number of smsll triangles
// in big one
static int num_triples_intriangle(int ndiv)
{
	return ndiv * ndiv;
}

// Check if a given index of node is on the list of given face list
static inline
int is_on_boundary_face(const int** Elem, const int type, const int nr,
		const int itr) {
	int result = 0;
	for (int i = 0; i < (type - 1); ++i) {
		if (Elem[nr][i] == itr) {
			result++;
			break;
		}
	}
	return result;
}

template<typename T>
T mark_edges_in_triangle(int s,int u,const int dim)
{
	// Mark corners
	if ((u + s == 0) || (s == dim) || (u == dim)) return T(1);
	// Mark top edge
	else if (u == 0) return T(2);
	// Mark skew edge
	else if (u + s == dim) return(3);
	// Mark left edge
	else if (s == 0) return T(4);
	// Mark internal vertices
	return T(0);
}

template<typename T>
T mark_edges_in_quad(int s,int u,const int s_dim,const int u_dim)
{
	// Mark corners
	if ((u + s == 0) || (u + s == u_dim + s_dim) || (s== 0 && u == u_dim) || (u == 0 && s == s_dim))
		return T(1);
	// Mark top edge
	else if (u == 0) return T(2);
	// Mark right edge
	else if (s == s_dim) return T(3);
	// Mark bottom edge
	else if (u == u_dim) return T(4);
	// Mark left edge
	else if (s == 0) return T(5);
	// Mark internal vertices
	return T(0);
}

static //unsignedhttp://wiadomosci.onet.pl/kraj/bulgaria-rsk-mig-nie-udzielal-polsce-licencji-na-remont-bulgarskich-mysliwcow/bjcmsj
int make_plane_XY0(int nptsx, int nptsy, double dlx, double dly,
		std::vector<Node_t>& vertices, std::vector<ubyte_t>& indices)
{
	assert(nptsx == nptsy);
	Node_t node;
	// Generate point coords in columns
	for (int y(0); y < nptsy; y++) {
		int nx = nptsy - y;
		for (int x(0); x < nx; x++) {
			node.position = fvmath::CVec3d(x * dlx, y * dly, -1.0);
			//node.edge_flag =  ((x > 0) && (y > 0) && ((x + y) < (nptsy -1))) ? 0.0f : 1.0f;
			node.info = mark_edges_in_triangle<Vertex::info_t>(x,y,nptsy-1);
			vertices.push_back(node);
			//std::cout << "p["<<y<<"]["<<x<<"] = " << node.position << " edf = " << node.info << std::endl;
		}
	}

	int base(0);
	nptsy--;
	//i = 0;
	for (int y(0); y < nptsy; y++) {
		int n = nptsx - y;
		bool cond = y < nptsy - 1;
		for (int x = 0; x < n; x++) {
			if ((x == (n - 1)) && cond) {
				// Repeat last index
				/*indices[i++] =*///indices.push_back( (ubyte_t)(base + 2*(n - 1)) );
				/*indices[i++] = */indices.push_back((ubyte_t) (base + x));
				indices.push_back((ubyte_t) (base + x));
			} else {
				indices.push_back((ubyte_t) (base + x));
				if (x < n - 1) /*indices[i++] = */
					indices.push_back((ubyte_t) (base + n + x));
			}
		}
		// add a degenerate triangle (except in a last row)
		if (cond) {
			/*indices[i++]*/indices.push_back((ubyte_t) (base + n - 1));
			/*indices[i++]*/indices.push_back((ubyte_t) (base + n));
		}

		base += n;
	}

	return static_cast<int>(indices.size());
}

static //unsigned
int make_plane_XY1(int nptsx, int nptsy, double dlx, double dly,
		std::vector<Node_t>& vertices, std::vector<ubyte_t>& indices)
{
	assert(nptsx == nptsy);
	Node_t node;
	// Generate point coords in columns
	for (int x(0); x < nptsx; x++) {
		int ny = nptsx - x;
		for (int y(0); y < ny; y++) {
			node.position = fvmath::CVec3d(x * dlx, y * dly, 1.0);
			//node.edge_flag =  ((x > 0) && (y > 0) && ((x + y + 1) < nptsx)) ? 0.0f : 1.0f;
			node.info = mark_edges_in_triangle<Vertex::info_t>(y,x,nptsx-1);
			vertices.push_back(node);
		}
	}

	int base(0);
	nptsx--;
	//i = 0;
	for (int x(0); x < nptsx; x++) {
		int n = nptsy - x;
		bool cond = x < nptsx - 1;
		for (int y = 0; y < n; y++) {
			if ((y == (n - 1)) && cond) {
				// Repeat last index
				indices.push_back((ubyte_t) (base + y));
				indices.push_back((ubyte_t) (base + y));
			} else {
				indices.push_back((ubyte_t) (base + y));
				if (y < n - 1)
					indices.push_back((ubyte_t) (base + n + y));
			}
		}
		// add a degenerate triangle (except in a last row)
		if (cond) {
			indices.push_back((ubyte_t) (base + n - 1));
			indices.push_back((ubyte_t) (base + n));
		}

		base += n;
	}

	return static_cast<int>(indices.size());
}

static //unsigned
int make_plane_XZ(int nptsx, int nptsz, double dlx, double dlz,
		std::vector<Node_t>& vertices, std::vector<ubyte_t>& indices)
{
	Node_t node;
	// Generate point coords in columns
	for (int x(0); x < nptsx; x++) {

		//int base = x * nptsz;

		for (int z(0); z < nptsz; z++) {
			node.position = fvmath::CVec3d(x * dlx, 0.0, 2 * z * dlz - 1.0);
			//node.edge_flag = ((x > 0) && ((x+1) < nptsx) && (z > 0) && ((z+1) < nptsz)) ? 0.0f : 1.0f;
			node.info = mark_edges_in_quad<Vertex::info_t>(z,x,nptsz-1,nptsx-1);
			vertices.push_back(node);
		}
	}

	//unsigned int i(0);
	nptsx--;

	for (int x(0); x < nptsx; x++) {

		int base = x * nptsz;

		for (int z(0); z < nptsz; z++) {
			indices.push_back((ubyte_t) (base + z));
			indices.push_back((ubyte_t) (base + nptsz + z));
		}
		// Add a degenerate triangle (except in a last column)
		if (x < nptsx - 1) {
			indices.push_back((ubyte_t) ((x + 1) * nptsz + (nptsz - 1)));
			indices.push_back((ubyte_t) ((x + 1) * nptsz));
		}
	}

	return (indices.size());
}

static //unsigned
int make_plane_XYZ(int nptsxy, int nptsz, double dlxy,
		double dlz, std::vector<Node_t>& vertices,
		std::vector<ubyte_t>& indices)
{
	Node_t node;
	// Generate point coords in columns
	for (int x(0); x < nptsxy; x++) {

		//int base = x * nptsz;

		for (int z(0); z < nptsz; z++) {
			//int index = base + z;
			//fvmath::Vec3d * pv = vertices + index;
			const double val = x * dlxy;
			node.position = fvmath::CVec3d(1.0 - val, val, 2 * z * dlz - 1.0);
			//node.edge_flag = ((x > 0) && ((x + 1) < nptsxy) && (z > 0) && ((z + 1) < nptsz)) ? 0.0f : 1.0f;
			node.info = mark_edges_in_quad<Vertex::info_t>(z,x,nptsz-1,nptsxy-1);
			vertices.push_back(node);
		}
	}

	//unsigned int i(0);
	nptsxy--;

	for (int x(0); x < nptsxy; x++) {

		int base = x * nptsz;

		for (int z(0); z < nptsz; z++) {
			indices.push_back((ubyte_t) (base + z));
			indices.push_back((ubyte_t) (base + nptsz + z));
		}
		// add a degenerate triangle (except in a last row)
		if (x < nptsxy - 1) {
			indices.push_back((ubyte_t) ((x + 1) * nptsz + (nptsz - 1)));
			indices.push_back((ubyte_t) ((x + 1) * nptsz));
		}
	}

	return static_cast<int>(indices.size());
}

static //unsigned
int make_plane_YZ(int nptsy, int nptsz, double dly, double dlz,
		std::vector<Node_t>& vertices, std::vector<ubyte_t>& indices)
{
	Node_t node;
	// Generate point coords in columns
	for (int y(0); y < nptsy; y++) {

		//int base = y * nptsz;

		for (int z(0); z < nptsz; z++) {
			//int index = base + z;
			//fvmath::Vec3d * pv = vertices + index;
			const double val = y * dly;
			node.position = fvmath::CVec3d(0.0, 1.0 - val, 2 * z * dlz - 1.0);
			//node.edge_flag = ((y > 0) && (y < nptsy) && (z > 0) && (z < nptsz)) ? 0.0f : 1.0f;
			node.info = mark_edges_in_quad<Vertex::info_t>(z,y,nptsz-1,nptsy-1);
			vertices.push_back(node);
		}
	}

	//int i(0);
	nptsy--;

	for (int y(0); y < nptsy; y++) {

		int base = y * nptsz;

		for (int z(0); z < nptsz; z++) {
			indices.push_back((ubyte_t) (base + z));
			indices.push_back((ubyte_t) (base + nptsz + z));
		}
		// add a degenerate triangle (except in a last row)
		if (y < nptsy - 1) {
			indices.push_back((ubyte_t) ((y + 1) * nptsz + (nptsz - 1)));
			indices.push_back((ubyte_t) ((y + 1) * nptsz));
		}
	}

	return static_cast<int>(indices.size());
}



//static unsigned int (* make_plane[5])(int,int,double,double,fvmath::Vec3d*,ubyte_t*) = {
//  make_plane_XY0, make_plane_XY1, make_plane_XZ, make_plane_XYZ, make_plane_YZ
//};

Field::Field(Mesh& pmesh_, const std::string& name_)
: BaseField()
, m_mesh(pmesh_)
, m_name(name_)
, m_min_deg(0)
, m_max_deg(0)
, m_min_val(FV_LARGE)
, m_max_val(-FV_LARGE)
{
	//mfp_log_debug("Field ctr\n");

	const std::size_t num = m_mesh.GetNumElems();
	m_elem_degs.reserve(num);
	m_elem_values.reserve(2*num);
}

Field::~Field() {
	mfp_log_debug("Dtr\n");
	this->Free();
}

int Field::Init(const int parent_id_, const char* fname_)
{
	//mfp_log_debug("Init\n");
	int result;
	if ((result = BaseField::Init(fname_, parent_id_)) < 0)
		return (-1);
	//mfp_debug("before of\n");
	m_name = fname_ ? std::string(fname_) : std::string("Unknown");
	if (this->type() == FieldSTD) {
		//mfp_debug("in std\n");
		// set min/max pdeg to default value
		m_min_deg = FV_DEF_MIN_P_STD;
		m_max_deg = FV_DEF_MIN_P_STD;
		// set render routine
		//_renderer = &Field::RenderSTD;

	}
	else if (this->type() == FieldDG) {
		//mfp_debug("In dg cectuion of init\n");
		// set render routine
		//_renderer = &Field::RenderDG;
		m_min_deg = FV_DEF_MAX_P_DG;
		m_max_deg = FV_DEF_MIN_P_DG;
		// get min/max value of pdeg
		int pdeg[PDC_MAXEQ];
		// here is the best way to do omp stuff
		Mesh::arElems::const_iterator it(m_mesh.GetElems().begin());
		const Mesh::arElems::const_iterator it_e(m_mesh.GetElems().end());
		for (; it != it_e; ++it) {
			(void) get_el_pdeg(this->idx(), it->eid, pdeg);
			m_min_deg = (m_min_deg < pdeg[0]) ? m_min_deg : pdeg[0];
			m_max_deg = (m_max_deg > pdeg[0]) ? m_max_deg : pdeg[0];
		}
		//
	} else {
		// Unknown field type
		assert(!"Unknown approximation type");
		return (-1);
	}

	CalculateMinMaxValues(m_min_val,m_max_val,m_elem_values);

	return (result);
}


// Calculate the minimal and maximal value
void Field::CalculateMinMaxValues(
		double & min_value,
		double & max_value,
		std::vector<double> & vec_values,
		const int ctrl) const
{
	double el_dofs[APC_MAXELSD], el_coords[3*MMC_MAXELVNO], sol[PDC_MAXEQ];
	int nodes[MMC_MAXELVNO+1];
	int pdeg;
	size_type nel(m_mesh.GetNumElems()), iel;
	min_value = FV_LARGE;
	max_value = -FV_LARGE;
	vec_values.reserve(2 *m_mesh.GetNumElems());
	vec_values.clear();
	std::vector<double> vSol;
	std::vector<fvmathParser::MathElement> mathElems;
	fvmathParser::MathCalculator::ONP(_solution.sFormula, _solution.nrEquetions,
				mathElems);
	fv_timer t;
	t.start();
	// Loop over elements in parent mesh
	for (iel = 1; iel <= nel; iel++)
	{
		// Skip inactive elements
		if (m_mesh.GetElemStatus(iel) != MMC_ACTIVE) continue;

		// Collect parameters
		int nr_dofs   = GetElementDofs(iel,el_dofs);
		int nr_shup   = GetNumberOfShapeFunc(iel,&pdeg);
		int base_type = GetElementBaseType(iel);

		// Get element coordinates
		m_mesh.GetElementCoordinates(iel,nodes,el_coords);

		// Loops over 10 point in each directions
		// Start at the first vertex in prism
		double elmin = FV_LARGE, elmax = -FV_LARGE;
		CVec3d pt;
		if (nodes[0] == 6) { // prism
			//mfp_debug("Prism%d\n",iel);
			for (pt.z = -1.0; pt.z <= 1.001; pt.z += 0.2) {
				for (pt.y = 0.0; pt.y <= 1.001; pt.y+= 0.1) {
					for (pt.x = 0.0; pt.x <= 1.001 - pt.y; pt.x += 0.1) {
						CalculateSolution(ctrl,pdeg,base_type,el_coords, pt.v, el_dofs, sol);
#if 1
						vSol.clear();
						for (int nreq = 0; nreq < _solution.nrEquetions; ++nreq)
								vSol.push_back(sol[nreq]);

						fvmathParser::MathCalculator::ONPCalculate(mathElems, vSol, sol[0]);
#endif
						// Element compare
						elmin = fv_min(elmin,sol[0]);
						elmax = fv_max(elmax,sol[0]);
						//printf("pt: {%lf %lf %lf}\n",pt.x,pt.y,pt.z);
						//printf("elmin: %lf elmax: %lf\n",elmin,elmax);
					}
				}
			}
		}
		else { // tetrahedron
			for (pt.z = 0.0; pt.z <= 1.001; pt.z += 0.1) {
				for (pt.y = 0.0; pt.y <= 1.001 - pt.z; pt.y += 0.1) {
					for (pt.x = 0.0; pt.x <= 1.001 - pt.y - pt.z; pt.x += 0.1) {
						CalculateSolution(ctrl,pdeg,base_type,el_coords, pt.v, el_dofs, sol);
#if 1
						vSol.clear();
						for (int nreq = 0; nreq < _solution.nrEquetions; ++nreq)
								vSol.push_back(sol[nreq]);

						fvmathParser::MathCalculator::ONPCalculate(mathElems, vSol, sol[0]);
#endif
						// Element compare
						elmin = fv_min(elmin,sol[0]);
						elmax = fv_max(elmax,sol[0]);
					}
				}
			}
		}
		// Global compare
		min_value = fv_min(min_value,elmin);
		max_value = fv_max(max_value,elmax);

		// Store element values in vector
		vec_values.push_back(elmin);
		vec_values.push_back(elmax);
	}

	t.stop();
	//mfp_log_debug("Min: %lf\tMax: %lf\n",min_value,max_value);
	//mfp_log_debug("Time for min/max value: %lf\n",t.get_time_sec() );
}

size_t Field::CountCoeffs(std::vector<int>* degree_ptr) const
{
	int pdeg;
	const int nreq = _solution.nrEquetions;
	if (degree_ptr) degree_ptr->resize(m_mesh.GetElems().size());

	size_t num_coeffs = 0;
	for (size_t i(0); i < m_mesh.GetElems().size(); ++i) {
		int el_id = m_mesh.GetElems()[i].eid;
		int num_shup = GetNumberOfShapeFunc(el_id, &pdeg);
		num_coeffs += get_el_pdeg_numshap(this->idx(), el_id, &pdeg) * nreq;
		if (degree_ptr) degree_ptr->push_back(pdeg);
	}
	return num_coeffs;
}

int Field::Free() {
	return BaseField::Free();
}

void Field::Reset() {
	this->BaseField::Reset();
}

void Field::Clear() {
}


int Field::RenderDG (
		std::vector<Vertex>& out_vertices,
		std::vector<unsigned>& out_indices,
		std::vector<int>& out_counts,
		std::vector<int64_t>& out_offsets
		)
{
	int pdeg[PDC_MAXEQ];
	int nr_nodes, nr_dofs;
	const int iaux = 2;
	double Node_coor[3 * MMC_MAXELVNO]; /* coord of nodes of El */
	double Dofs_loc[APC_MAXELSD]; /* element solution dofs */
	double Xcoor[3]; /* global coord of gauss point */
	double U_val[PDC_MAXEQ]; /* computed solution */
	double U_x[PDC_MAXEQ]; /* gradient of computed solution */
	double U_y[PDC_MAXEQ]; /* gradient of computed solution */
	double U_z[PDC_MAXEQ]; /* gradient of computed solution */
	double Base_phi[APC_MAXELVD]; /* basis functions */
	double Base_dphix[APC_MAXELVD]; /* x-derivatives of basis function */
	double Base_dphiy[APC_MAXELVD]; /* y-derivatives of basis function */
	double Base_dphiz[APC_MAXELVD]; /* y-derivatives of basis function */
	// Vector of solutions
	std::vector<double> vSol;
	std::vector<fvmathParser::MathElement> mathElems;

	//std::cout<<"Evaluating formula is: " << _solution.sFormula << std::endl;
	fvmathParser::MathCalculator::ONP(_solution.sFormula, _solution.nrEquetions,
			mathElems);

	Mesh::arElConstIter it = m_mesh.GetElems().begin();
	const Mesh::arElConstIter ite = m_mesh.GetElems().end();

	assert(m_mesh.GetElems().size()!=0);

	// Resize containers to the larges evaliabe number according pdeg = 707
	std::vector<Node_t>  nodeCoords(121);	// vector container for coords (9+2) * (9+2) points
	std::vector<ubyte_t> indexNode(207);// vecot container for indexes of nodes 9*21 + 16


	int vertex_counter[2] = {0,0};

	//vEvaluated.reserve(FV_MAX_POINTS);
	int my_couter = 0;

	std::vector<MeshVertex> *vertices_ptr, tvertices, qvertices;
	std::vector<int> *indices_ptr,tindices, qindices, *counts_ptr, tcounts, qcounts;
	std::vector<int64_t> offsets;

	MeshVertex vt;
	while (it != ite)
	{
		// Get element id form it
		int nel(it->eid);
		// Specify type of an alement
		int nr_faces = it->is_prism() ? 5 : 4;
		assert(nr_faces == 5);
		//mfp_debug("nr_faces = %d\n",nr_faces);
		// Set startup handles
		vertices_ptr = &tvertices;
		indices_ptr = &tindices;
		counts_ptr = &tcounts;
		// Get element real coordinates
		nr_nodes = m_mesh.GetElementCoordinates(nel, NULL, Node_coor);
		// Get element dofs
		nr_dofs = get_element_dofs(this->idx(), nel, _solution.currSolution,
				Dofs_loc);
		// Element base function type
		const int base = get_base_type(this->idx(), nel);
		// Get its corresponding pdeg
		(void) get_el_pdeg(this->idx(), nel, pdeg);
		// Specify a pdeg and number of discretizations
		int nptsx, nptsy;
		if (base == APC_BASE_TENSOR_DG) {
			nptsx = pdeg[0] % 100 + 1;
			nptsy = pdeg[0] / 100 + 1;
		} else {	// APC_BASE_COMPLETE_DG
			nptsx = pdeg[0] + 1;
			nptsy = pdeg[0] + 1;
		}
		nptsx+= dg_density;
		nptsy+= dg_density;
		//printf("Pdeg od pdeg[0] = %d\n",pdeg[0]);
		// Lenghts of segments
		double dlx = 1.0 / nptsx; // in x and y directions
		double dly = 1.0 / nptsy; // in z directions

		// increase to real number of points, 
		// because counter of points is +1
		nptsx++;
		nptsy++;

		// Loop over boundary faces
		int flag = 0;
		for (int nfa(0); nfa < nr_faces; nfa++) {

			// Skip a non boundary face
			if (!IS_SET_FACEN(it->id, nfa))
				continue;

			bool cond = nfa < 2;
			const int nno  = cond ? (nptsx + 1) * nptsx / 2 : nptsx * nptsy;
			const int dest = cond ? 0 : 1;

			// Change handles for quadratic faces
			if (nr_faces == 5 && !cond && !flag) {
				//mfp_debug("This is happedn only one\n");
				vertices_ptr = &qvertices;
				indices_ptr = &qindices;
				counts_ptr = &qcounts;
				flag++;
			}

			// Clear temporary containers
			nodeCoords.clear();
			indexNode.clear();

			unsigned int size;
			//CVec3f Normal;
			switch (nfa)
			{
			case 0:
				// Tesselate face
				size = make_plane_XY0(nptsx, nptsx, dlx, dlx, nodeCoords,
						indexNode);
				//Normal = fvmath::GetNormal2Plane<float, double>(&Node_coor[0],
				//		&Node_coor[3], &Node_coor[6]);
				//mfp_debug("Traingle base0: inidces: %u nodes: %u\n",size,nodeCoords.size());
				break;
			case 1:
				// Tesselate face
				size = make_plane_XY1(nptsx, nptsx, dlx, dlx, nodeCoords,
						indexNode);
				//Normal = fvmath::GetNormal2Plane<float, double>(&Node_coor[9],
				//		&Node_coor[15], &Node_coor[12]);
				//mfp_debug("Traingle base1: inidces: %u nodes: %u\n",size,nodeCoords.size());
				break;
			case 2:
				// Tesselate face
				size = make_plane_XZ(nptsx, nptsy, dlx, dly, nodeCoords,
						indexNode);
				//Normal = fvmath::GetNormal2Plane<float, double>(&Node_coor[0],
				//		&Node_coor[9], &Node_coor[3]);
				//mfp_debug("Quad face1: inidces: %u nodes: %u\n",size,nodeCoords.size());
				break;
			case 3:
				// Tesselate face
				size = make_plane_XYZ(nptsx, nptsy, dlx, dlx, nodeCoords,
						indexNode);
				//Normal = fvmath::GetNormal2Plane<float, double>(&Node_coor[3],
				//		&Node_coor[12], &Node_coor[6]);
				//mfp_debug("Quad face1: inidces: %u nodes: %u\n",size,nodeCoords.size());
				break;
			case 4:
				// Tesselate face
				size = make_plane_YZ(nptsx, nptsy, dlx, dlx, nodeCoords,
						indexNode);
				//Normal = fvmath::GetNormal2Plane<float, double>(&Node_coor[6],
				//		&Node_coor[15], &Node_coor[0]);
				//mfp_debug("Quad face2: inidces: %u nodes: %u\n",size,nodeCoords.size());
				break;
			}					// switch
			// Add count of primitives
			assert(size == indexNode.size());

			for (size_t i(0); i < indexNode.size(); ++i) {
				int index = vertex_counter[dest] + indexNode[i];
				indices_ptr->push_back(index);
				//std::cout << "dest" << dest << " index" << i << " = " << index << std::endl;
			}


			counts_ptr->push_back(size);

			Node_t * pXloc(nodeCoords.data());
			// Loop over discretizated local points
			for (int ino(0); ino < nno; ino++) {
				if (iaux == 2) {
					(void) apr_elem_calc_3D_mod(iaux, _solution.nrEquetions,
							pdeg, base, pXloc[ino].position.v, Node_coor, Dofs_loc,
							Base_phi, Base_dphix, Base_dphiy, Base_dphiz, vt.position.v,
							U_val, U_x, U_y, U_z, NULL);
				} else {
					(void) apr_elem_calc_3D_mod(iaux, _solution.nrEquetions,
							pdeg, base, pXloc[ino].position.v, Node_coor, Dofs_loc,
							Base_phi,
							NULL, NULL, NULL, vt.position.v, U_val,
							NULL, NULL, NULL, NULL);
				}

				// Evaluate solution formula
				vSol.clear();

				for (int sol = 0; sol < _solution.nrEquetions; ++sol) {
					vSol.push_back(U_val[sol]);
				}

				fvmathParser::MathCalculator::ONPCalculate(mathElems, vSol,	vt.value);
				// Add vertex
				vt.info = pXloc[ino].info;
				//std::cout << "info: " << vt.info << "\n";
				//vt.value = static_cast<float>(dTmp);
				vertices_ptr->push_back(vt);
				//std::cout << "vertex: " << vt.position << " info: " << vt.info << " value: " << vt.value << std::endl;
				//accum_ptr->addVertex(vt,dest);

				//FV_ASSERT(dTmp>=_values.first && dTmp<=_values.second);
			}						// fo ino

			//accum_ptr->addBaseIndex(vertex_counter[dest],dest);
			vertex_counter[dest] += nno;
		}// for nfa
		//}// else pdeg
		++it;
	}			//while



	// Update the legend
	Legend& rleg = ViewManagerInst().GetLegend();

	// getter vertices into single container
	out_vertices.resize(tvertices.size()+qvertices.size());
	std::vector<Vertex>::iterator vt_p = out_vertices.begin();

	for (const auto& v : tvertices) {
		vt_p->position = v.position;
		vt_p->info = static_cast<Vertex::info_t>(v.info);
		rleg.GetColor(v.value);
		vt_p->value = static_cast<float>(v.value);
		//vt_p->color.x = rleg.GetColor3b(v.value);
		vt_p->color.x = rleg.color.rgb.r;
		vt_p->color.y = rleg.color.rgb.g;
		vt_p->color.z = rleg.color.rgb.b;
		//std::cout << "color: " << vt_p->color << std::endl;
		++vt_p;
	}
	//mfp_debug("Number of tri vertices: %u and qua: %u\n",tvertices.size(),qvertices.size());
	for (const auto& v : qvertices) {
		vt_p->position = v.position;
		vt_p->info = static_cast<Vertex::info_t>(v.info);
		rleg.GetColor(v.value);
		vt_p->value = static_cast<float>(v.value);
		//vt_p->color.x = rleg.GetColor3b(v.value);
		vt_p->color.x = rleg.color.rgb.r;
		vt_p->color.y = rleg.color.rgb.g;
		vt_p->color.z = rleg.color.rgb.b;
		//std::cout << "color: " << vt_p->color << std::endl;
		++vt_p;
	}

	// Update indices for quad faces
	// and gether all indices into single container
	const unsigned int lsize = tvertices.size();
	out_indices.reserve(lsize + qindices.size());
	out_indices.insert(out_indices.end(), tindices.begin(), tindices.end());
	for (auto index : qindices)
		out_indices.push_back(index + lsize);

	// Merge counters
	out_counts.reserve(tcounts.size()+qcounts.size());
	out_counts.insert(out_counts.end(),tcounts.begin(),tcounts.end());
	out_counts.insert(out_counts.end(),qcounts.begin(),qcounts.end());
	//for (const auto& cnt : out_counts)
	//	std::cout << "count: " << cnt << std::endl;
	// Calculate addreses
	size_t p(0);
	std::vector<int>::const_iterator id_p = out_counts.begin();
	//std::cout << "Adresses:\n";
	while (id_p != out_counts.end()) {
		out_offsets.push_back(p);
		//std::cout << "\t" << p << std::endl;
		p += sizeof(unsigned int)*(*id_p);
		++id_p;
	}
	//mfp_debug("Number of qu vertices: %u\n",ntri_vertices);
	//accum_ptr->UseVBO() = ViewManagerInst().IsVBOUsed();
	//accum_ptr->aRenderTraingleStrips = true;
	//FV_CHECK_ERROR_GL();

	return static_cast<int>(tindices.size());
}

int Field::RenderDGCutted(
		const std::vector<CutPlane>& vPlanes,
		std::vector<Vertex>& out_vertices,
		std::vector<unsigned>& out_indices
		)
{
	//static const int FV_MAX_POINTS = 1000000;
	//std::vector<double> vEvaluated;

	int pdeg[PDC_MAXEQ];
	int nr_nodes, nr_dofs;
	const int iaux = 2;
	double Node_coor[3 * MMC_MAXELVNO]; /* coord of nodes of El */
	double Dofs_loc[APC_MAXELSD]; /* element solution dofs */
	double Xcoor[3]; /* global coord of gauss point */
	double U_val[PDC_MAXEQ]; /* computed solution */
	double U_x[PDC_MAXEQ]; /* gradient of computed solution */
	double U_y[PDC_MAXEQ]; /* gradient of computed solution */
	double U_z[PDC_MAXEQ]; /* gradient of computed solution */
	double Base_phi[APC_MAXELVD]; /* basis functions */
	double Base_dphix[APC_MAXELVD]; /* x-derivatives of basis function */
	double Base_dphiy[APC_MAXELVD]; /* y-derivatives of basis function */
	double Base_dphiz[APC_MAXELVD]; /* y-derivatives of basis function */
	// Vector of solutions
	std::vector<double> vSol;
	std::vector<fvmathParser::MathElement> mathElems;

	//std::cout<<"Evaluating formula is: " << _solution.sFormula << std::endl;
	fvmathParser::MathCalculator::ONP(_solution.sFormula, _solution.nrEquetions,
			mathElems);

	Mesh::arElConstIter it = m_mesh.GetElems().begin();
	const Mesh::arElConstIter ite = m_mesh.GetElems().end();

	// Resize containers to the larges evaliabe number according pdeg = 707
	std::vector<Node_t>  nodeCoords(121);	// vector container for coords (9+2) * (9+2) points
	std::vector<ubyte_t> indexNode(207);// vecot container for indexes of nodes 9*21 + 16


	int vertex_counter[2] = {0,0};

	//vEvaluated.reserve(FV_MAX_POINTS);
	int my_couter = 0;
	//double minVal = LARGE_F;
	const unsigned maxdiv = 5;
	const unsigned max_num_verts = (maxdiv+1)*(maxdiv+2)/2*vPlanes.size()*m_mesh.GetElems().size();
	const unsigned max_num_idxs = maxdiv*maxdiv*vPlanes.size()*m_mesh.GetElems().size();
	mfp_log_info("Max division: %u\tmax vertices size: %u\tmax indices size: %u\n",maxdiv,max_num_verts,max_num_idxs);

	std::vector<myVertex> mesh_vertices(max_num_verts);
	mesh_vertices.clear();
	out_indices.reserve(max_num_idxs);
	out_indices.clear(); // ?

	std::vector<int64_t> offsets;
	Vertex vt;
	Legend& rleg = ViewManagerInst().GetLegend();
	int counter = 0;

	mfp_debug("In DG Cutted");
	while (it != ite)
	{
		//mfp_debug("Fisrt loop\n");
		// Get element id form it-erator
		int nel(it->eid);
		// Set type
		int type(it->is_prism());
		// Specify type of an alement
		const int nr_faces = type ? 5 : 4;
		// Get element real coordinates
		nr_nodes = m_mesh.GetElementCoordinates(nel, NULL, Node_coor);
		// Get element dofs
		nr_dofs = get_element_dofs(this->idx(), nel, this->_solution.currSolution,Dofs_loc);
		// Element base function type
		const int base = get_base_type(this->idx(), nel);
		// Get its corresponding pdeg
		(void) get_el_pdeg(this->idx(), nel, pdeg);
		int ndiv = maximal_pdeg(pdeg[0]);

		// Search for element in plane
		for (int ip=0, npls=(int)vPlanes.size();ip<npls;++ip) {
			// Check if element is in plane's elements
			// Slice element
			unsigned new_size = CutPlane::CutElement(&vPlanes[ip],Node_coor,it->nodes,counter,ndiv,mesh_vertices,out_indices);
			//for (const auto & v : mesh_vertices) std::cout << v.position << std::endl;
			if (0 == new_size) continue;
			//mfp_debug("just after cutt in cutplane\n");
			// Get reference coordinates
			unsigned offset = mesh_vertices.size() - new_size;

			for (unsigned iv = 0; iv < new_size; ++iv)
			{
				CVec3d resultPt = Prizm::TransformToReference(Node_coor,mesh_vertices[offset+iv].position);
				// Correct point into FEM prism
				CVec3d refPt(resultPt.x*0.5+0.5,resultPt.y*0.5+0.5,resultPt.z);
				assert(0.0 <= refPt.x <=1.0);
				assert(0.0 <= refPt.y <=1.0);
				assert(-1.0 <= refPt.z <=1.0);
				//std::cout << "World point: " << mesh_vertices[iv].position << " reference point: " << refPt << std::endl;
				if (iaux == 2) {
					(void) apr_elem_calc_3D_mod(iaux, _solution.nrEquetions,
							pdeg, base, refPt.v, Node_coor, Dofs_loc,
							Base_phi, Base_dphix, Base_dphiy, Base_dphiz, Xcoor,
							U_val, U_x, U_y, U_z, NULL);
				} else {
					(void) apr_elem_calc_3D_mod(iaux, _solution.nrEquetions,
							pdeg, base, refPt.v, Node_coor, Dofs_loc,
							Base_phi,
							NULL, NULL, NULL, /*Xcoor*/NULL, U_val,
							NULL, NULL, NULL, NULL);
				}

				// Evaluate solution formula
				vSol.clear();

				for (int sol = 0; sol < _solution.nrEquetions; ++sol) {
					vSol.push_back(U_val[sol]);
				}

				double dTmp;
				fvmathParser::MathCalculator::ONPCalculate(mathElems, vSol,	dTmp);

				// Add vertex
				vt.position = mesh_vertices[offset+iv].position;
				vt.info = static_cast<Vertex::info_t>(mesh_vertices[offset+iv].type);
				vt.value = static_cast<float>(dTmp);
				rleg.GetColor(dTmp);
				vt.color.x = rleg.color.rgb.r;
				vt.color.y = rleg.color.rgb.g;
				vt.color.z = rleg.color.rgb.b;
				//std::cout << "vertex: " << vt.position << " info: " << vt.info << std::endl;
				out_vertices.push_back(vt);
				// Cut element
			}// end over slice points
		}// end over cut planes
		++it;
		++counter;
	}//while


	//FV_CHECK_ERROR_GL();

	return static_cast<int>(out_indices.size());
}

int Field::RenderSTD(
		std::vector<Vertex>& out_vertices,
		std::vector<unsigned>& out_indices
		)
{
	//mfp_log_debug("RenderSTD\n");
	int pdeg[PDC_MAXEQ];
	int nr_nodes, nr_dofs;
	const int iaux = 2;
	double Node_coor[3 * MMC_MAXELVNO]; /* coord of nodes of El */
	double Dofs_loc[APC_MAXELSD]; /* element solution dofs */
	double Xcoor[3]; /* global coord of gauss point */
	double U_val[PDC_MAXEQ]; /* computed solution */
	double U_x[PDC_MAXEQ]; /* gradient of computed solution */
	double U_y[PDC_MAXEQ]; /* gradient of computed solution */
	double U_z[PDC_MAXEQ]; /* gradient of computed solution */
	double Base_phi[APC_MAXELVD]; /* basis functions */
	double Base_dphix[APC_MAXELVD]; /* x-derivatives of basis function */
	double Base_dphiy[APC_MAXELVD]; /* y-derivatives of basis function */
	double Base_dphiz[APC_MAXELVD]; /* y-derivatives of basis function */
	// Vector of solutions
	std::vector<double> vSol;
	std::vector<fvmathParser::MathElement> mathElems;

	//std::cout<<"Evaluating formula is: " << _solution.sFormula << std::endl;
	fvmathParser::MathCalculator::ONP(_solution.sFormula, _solution.nrEquetions,
			mathElems);

	Mesh::arElConstIter it = m_mesh.GetElems().begin();
	const Mesh::arElConstIter ite = m_mesh.GetElems().end();

	std::vector<unsigned int> quad_indices;

	std::vector<unsigned>* ptr[2] = { &out_indices, &quad_indices };
	MeshVertex vts[MMC_MAXELVNO];

	Legend& rleg = ViewManagerInst().GetLegend();
	std::vector<MeshVertex> mesh_vertices(m_mesh.GetNumNodes());
	std::vector<int> lindices;
	mesh_vertices.clear();
	//_pobject->SetPointSize(2.0);
	fvmath::CVec3f normal;
	unsigned int iv = 0;
	int target;
	//mfp_debug("Before while loop\n");
	double* Xloc;
	const int* elem_p;
	int nr_faces, Nodes[MMC_MAXELVNO+1];

	while (it != ite)
	{
		// Retrieve element index
		int nel(it->eid);
		// Element base function type
		const int base = get_base_type(this->idx(), nel);
		// Set reference element
		if (IS_PRISM(it->id)) {
			Xloc = (double*) &XlocPrizm[0];
			elem_p = Mesh::prizm[0];
			nr_faces = 5;

		}
		else { // Tetrahedorn
			Xloc = (double*) &XlocTetra[0];
			elem_p = Mesh::tetra[0];
			nr_faces = 4;
		}
		// Get element dofs
		int nr_dofs = get_element_dofs(this->idx(), nel, _solution.currSolution,
				Dofs_loc);
		// Get its corresponding pdeg
		(void) get_el_pdeg(this->idx(), nel, pdeg);
		// Get element real coordinates
		nr_nodes = m_mesh.GetElementCoordinates(nel, Nodes, Node_coor);

		// Loop over local indices
		lindices.clear();
		const unsigned int offset = mesh_vertices.size();
		unsigned int no_status = 0;
		for (int nf = 0; nf < nr_faces; ++nf) {
			// Check if current face is boundary
			if (!IS_SET_FACEN(it->id, nf)) continue;
			// Mark which node is calculated
			int of = nf * (nr_faces - 1);

			int v0 = elem_p[of];
			no_status |= (1 << v0);
			//lindices.push_back(v0);

			int v1 = elem_p[of+1];
			no_status |= (1 << v1);
			//lindices.push_back(v1);

			int v2 = elem_p[of+2];
			no_status |= (1 << v2);
			//lindices.push_back(v2);

			// For prism
			if (nr_faces == 5 && nf >= 2)
			{
				int v3 = elem_p[of+3];
				no_status |= (1 << v3);
				//lindices.push_back(v0);
				//lindices.push_back(v2);
				//lindices.push_back(v3);
			}
		}

		// Calculate values
		for (unsigned int ino=0, cnt = 0; ino < MMC_MAXELVNO; ++ino)
		{
			if (!(no_status & (1 << ino))) continue;

			vts[ino].position = fvmath::CVec3d(&Node_coor[3*ino]);

			if (iaux == 2) {
				apr_elem_calc_3D_mod(iaux, _solution.nrEquetions, pdeg, base,
						&Xloc[3 * ino], Node_coor, Dofs_loc, Base_phi,
						Base_dphix, Base_dphiy, Base_dphiz, Xcoor, U_val, U_x,
						U_y, U_z, NULL);
			} else {
				apr_elem_calc_3D_mod(iaux, _solution.nrEquetions, pdeg, base,
						&Xloc[3 * ino], Node_coor, Dofs_loc, Base_phi,
						NULL, NULL, NULL, Xcoor, U_val,
						NULL, NULL, NULL, NULL);
			}
			//mfp_debug("U_val[%d] = %lf\n",ino,U_val[0]);
			vSol.clear();

			for (int sol = 0; sol < _solution.nrEquetions; ++sol)
				vSol.push_back(U_val[sol]);

			fvmathParser::MathCalculator::ONPCalculate(mathElems, vSol, vts[ino].value);

			vts[ino].info = Nodes[ino + 1];
			vts[ino].type = offset + cnt;
			cnt++;

			mesh_vertices.push_back(vts[ino]);
		}

		// Now we have to add indices of vertces
		for (int nf = 0; nf < nr_faces; ++nf) {

			// Check if current face is boundary
			if (!IS_SET_FACEN(it->id, nf)) continue;

			int target = (nr_faces == 5 && nf >= 2);
			// for each face add triangle; for quads there will be 2 triangles
			int of = nf * (nr_faces - 1);

			//mfp_debug("of = %d %d\n",of,nr_faces);
			int v0 = elem_p[of++];
			int v1 = elem_p[of++];
			int v2 = elem_p[of];

			ptr[target]->push_back( vts[v0].type );
			ptr[target]->push_back( vts[v1].type );
			ptr[target]->push_back( vts[v2].type );

			if (target) {
				// Quad face = 2 triangles to add
				int v3 = elem_p[++of];
				ptr[target]->push_back( vts[v2].type );
				ptr[target]->push_back( vts[v3].type );
				ptr[target]->push_back( vts[v0].type );
			}
		}
		++it;
	}

	int num_tri_fa = static_cast<int>(out_indices.size()) / 3;
	mfp_log_debug("Number of tri faces is %d",num_tri_fa);

	std::copy(quad_indices.begin(),quad_indices.end(),std::back_inserter(out_indices));

	mfp_log_debug("vertices size: %u",mesh_vertices.size());
	fvmath::removeDuplicates(mesh_vertices,out_indices);

	out_vertices.resize(mesh_vertices.size());
	std::vector<Vertex>::iterator vt_p = out_vertices.begin();

	for (const auto& v : mesh_vertices)
	{
		vt_p->position = v.position;
		vt_p->info = static_cast<Vertex::info_t>(v.info);
		rleg.GetColor(v.value);
		vt_p->value = static_cast<float>(v.value);
		//vt_p->color.x = rleg.GetColor3b(v.value);
		vt_p->color.x = rleg.color.rgb.r;
		vt_p->color.y = rleg.color.rgb.g;
		vt_p->color.z = rleg.color.rgb.b;
		//std::cout << "color: " << vt_p->color << std::endl;
		++vt_p;
	}



	// Remove dupliactes
	//fvmath::removeDuplicates(out_vertices,out_indices);
	mfp_log_debug("Total number of triangles is %u",out_indices.size());
	return num_tri_fa;
}

size_t Field::PackCoeffs(std::vector<mfvFloat_t>& Coeffs,std::vector<int>& Degree)
{
	size_t num_elems = m_mesh.GetElems().size();
	assert(num_elems > 0);

	// Get information about amount of all coefficients
	// and resize the buffer
	size_t total_ndofs = CountCoeffs(&Degree);
	Coeffs.resize(total_ndofs);

	// Accumulate coefficients
	size_t nr_dofs = 0;
	for (size_t nel = 0; nel < num_elems; ++nel) {
		const ElemInfo *info_ptr = &m_mesh.GetElems()[nel];
		int el_id = info_ptr->eid;
		nr_dofs += GetElementDofs(el_id,&Coeffs[nr_dofs]);
	}

	return nr_dofs;
}

uint32_t Field::FillArrayWithCoeffs(const Field *pField, CoordType Coeffs[]) {
	Mesh::arElConstIter itb(pField->m_mesh.GetElems().begin());
	const Mesh::arElConstIter ite(pField->m_mesh.GetElems().end());
	uint32_t totalNumCoeffs(0);
	for (; itb != ite; ++itb) {
		int elId = itb->eid;
		totalNumCoeffs += pField->GetElementDofs(elId, &Coeffs[totalNumCoeffs]);
	}
	return totalNumCoeffs;
}

int Field::GetDegreeVector(int nel, int Pdeg[]) const {
	int pdeg = GetElementDegree(nel);
	int base = GetElementBaseType(nel);
	switch (base) {
	case APC_BASE_TENSOR_DG:
		Pdeg[0] = Pdeg[1] = pdeg % 100;
		Pdeg[2] = pdeg / 100;
		break;
	case APC_BASE_COMPLETE_DG:
	case APC_BASE_PRISM_STD:
	case APC_BASE_TETRA_STD:
		Pdeg[0] = Pdeg[1] = Pdeg[2] = pdeg;
		break;
	default:
		exit(-1);
	}

	return base;
}

int Field::CalculateSolution(int control, int pdeg, int base, double* ElCoords,
		double* Ploc, double* Coeffs, double* Sol, double* DSolx, double* DSoly,
		double* DSolz) const {
	double BasePhi[APC_MAXELSD];
	double BaseDPhix[APC_MAXELSD];
	double BaseDPhiy[APC_MAXELSD];
	double BaseDPhiz[APC_MAXELSD];

	if (control == 2) {
		return apr_elem_calc_3D_mod(control, _solution.nrEquetions, &pdeg, base,
				Ploc, ElCoords, Coeffs, BasePhi, BaseDPhix, BaseDPhiy,
				BaseDPhiz, NULL, Sol, DSolx, DSoly, DSolz, NULL);
	}
	return apr_elem_calc_3D_mod(control, _solution.nrEquetions, &pdeg, base,
			Ploc, ElCoords, Coeffs, BasePhi,
			NULL, NULL, NULL, NULL, Sol,
			NULL, NULL, NULL, NULL);

}

}				// end namespace
