#ifndef _APPROX_MANAGER_H_
#define _APPROX_MANAGER_H_

#include"Plugin.h"
#include"ApproxServer.h"

#include<string>
#include<map>

namespace Approx {

	class IApproximator {
	public:
		virtual ~IApproximator() {}
		virtual int  init_mesh(int Control, char * Filename) = 0;
		virtual int  free_mesh(int Mesh_id) = 0;
		virtual int  get_nmno(int Mesh_id) = 0;
		virtual int  get_nmed(int Mesh_id) = 0;
		virtual int  get_nmfa(int Mesh_id) = 0;
		virtual int  get_nmel(int Mesh_id) = 0;
		virtual int  get_next_act_elem(int Mesh_id,int El) = 0;
		virtual int  get_next_elem(int Mesh_id,int El) = 0;
		virtual int  get_next_face(int Mesh_id,int Fa) = 0;
		virtual int  get_next_edge(int Mesh_id,int Ed) = 0;
		virtual int  get_next_node(int Mesh_id,int Node) = 0;
		virtual int  get_el_status(int Mesh_id,int el) = 0;
		virtual int  get_el_type(int Mesh_id,int El) = 0;
		virtual int  get_el_mate(int Mesh_id,int El) = 0;
		virtual int  get_el_faces(int Mesh_id,int El,int* Faces, int* Orient) = 0;
		virtual int  get_el_fam(int Mesh_id,int El,int* Elsons,int* Type) = 0;
		virtual int  get_el_struct(int Mesh_id,int El_id,int *El_struct) = 0;
		virtual int  get_el_node_coor(int Mesh_id,int El,int* Nodes, double *Xcoor) = 0;
		virtual int  get_fa_status(int Mesh_id,int Fa) = 0;
		virtual int  get_fa_type(int Mesh_id,int Fa) = 0;
		virtual int  get_fa_sub_bnd(int Mesh_id,int Fa) = 0;
		virtual int  get_fa_node_coor(int Mesh_id,int Fa,int* Nodes,double* Coords) = 0;
		virtual void get_fa_neig(int Mesh_id, int Fa, int* Fa_neig,
			int* Neig_sides, int* Node_shift, int* Diff_gen,
			double* Acceff, double* Bcceff) = 0;
		virtual int  get_face_edges(int Mesh_id,int Fa,int *Edges,int *Orient) = 0;
		virtual int  get_face_struct(int Mesh_id,int Fa,int *Fa_struct) = 0;
		virtual int  get_edge_status(int Mesh_id,int Ed) = 0;
		virtual int  get_edge_struct(int Mesh_id,int Ed,int *Ed_struct) = 0;
		virtual int  get_edge_nodes(int Mesh_id,int Ed,int *Ed_nodes) = 0;
		virtual int  get_node_status(int Mesh_id,int Node) = 0;
		virtual int  get_node_coor(int Mesh_id,int Node,double *Coor) = 0;
		// Aproximation interface
		virtual int  init_field(char Field_type,int Mesh_id,const char * Filename) = 0;
		virtual int  free_field(int Field_id) = 0;
		virtual int  get_nreq(int Field_id) = 0;
		virtual int  get_nr_sol(int Field_id) = 0;
		virtual int  get_base_type(int Field_id) = 0;
		virtual int  get_el_pdeg(int Field_id,int El,int* Pdeg_vec) = 0;
		virtual int  get_element_dofs(int Field_id,int El_id,int Vec_id,double *Dofs) = 0;
		virtual int  create_constr_data(int Field_id) = 0;
		virtual double apr_elem_calc(
		/* returns: Jacobian determinant at a point, either for */
		/* 	volume integration if Vec_norm==NULL,  */
		/* 	or for surface integration otherwise */
		int Control,	    /* in: control parameter (what to compute): */
					/*	1  - shape functions and values */
					/*	2  - derivatives and jacobian */
					/* 	>2 - computations on the (Control-2)-th */
					/*	     element's face */
		int Nreq,	    /* in: number of equations */
		int *Pdeg_vec,	    /* in: element degree of polynomial */
		int Base_type,	    /* in: type of basis functions: */
					/* 	1 (APC_TENSOR) - tensor product */
					/* 	2 (APC_COMPLETE) - complete polynomials */
		double *Eta,	    /* in: local coordinates of the input point */
		double *Node_coor,  /* in: array of coordinates of vertices of element */
		double *Sol_dofs,   /* in: array of element' dofs */
		double *Base_phi,   /* out: basis functions */
		double *Base_dphix, /* out: x-derivatives of basis functions */
		double *Base_dphiy, /* out: y-derivatives of basis functions */
		double *Base_dphiz, /* out: z-derivatives of basis functions */
		double *Xcoor,	    /* out: global coordinates of the point*/
		double *Sol,        /* out: solution at the point */
		double *Dsolx,      /* out: derivatives of solution at the point */
		double *Dsoly,      /* out: derivatives of solution at the point */
		double *Dsolz,      /* out: derivatives of solution at the point */
		double *Vec_nor     /* out: outward unit vector normal to the face */
		) = 0;
	};

	class ApproxManager {
	private:
		static ApproxManager *pInstance;

		typedef std::map<std::string,Plugin> ApproxModuleMap;
		ApproxModuleMap LoadedModules;
		ApproxServer    ServerModules;

		// Private constructor
		ApproxManager();
	public:
		static ApproxManager & GetInstance();

		ApproxServer & GetApproxServer() { return ServerModules; }

		void LoadApproxModule(const std::string & sModuleName);
	};




} // end namespace Approx

#endif /* _APPROX_MANAGER_H_
*/
