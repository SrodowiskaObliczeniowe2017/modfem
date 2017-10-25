#include "../include/ApproxModule.h"
#include "../include/fv_compiler.h"
#include "../include/fv_config.h"
#include "ApproxManager.h"
#include "Plugin.h"

#include <stdexcept>
#include <iostream>
#include <sstream>

using namespace Approx;

#define FV_SIZEOF_TABLE(A) (sizeof(A)/sizeof(A[0]))

static const char* moduleNames[] = {dgDllName, stdDllName, stdhDllName};

static IApproximator* pDLLApprox = NULL;

static inline bool is_init()
{
	if( ! pDLLApprox ) {
		std::cerr << "Error! Approximation module isn't set!" << std::endl;
		return(false);
	}
	return(true);
}

int init()
{
	int res(1);
	Plugin::hlibModule hPlug;

	for(uint_t i(0);i<FV_SIZEOF_TABLE(moduleNames);i++)
	{
		if(!(hPlug = Plugin::loadPlugin(moduleNames[i])))
		{
			std::cerr << "Warning! Can't load module: " << moduleNames[i] << std::endl;
			std::cerr << "Check if module is available in an executable directory.\n";
			res = -1;
		}
		else {
			std::cout << "Module: " << moduleNames[i] << " is present\n";
			Plugin::unloadPlugin(hPlug);
		}
	}

	return(res);
}

void destroy()
{
	delete & ApproxManager::GetInstance();
}

int set_approx_module(apr_type type)
{
	try 
	{
		uint_t key;
		if(type == APR_DG_PRSIM) 		{ key = 0; }
		else if(type == APR_STD_PRISM) 	{ key = 1; }
		else if(type == APR_STD_HYBRID) { key = 2; }
		else
		{
			std::ostringstream os("Unknown filed type: ");
			os << int(type) << " \n";
			throw std::runtime_error(os.str());
		}

		std::cout << "Loading " << moduleNames[key] << std::endl;
		std::string module(moduleNames[key]);
		ApproxManager::GetInstance().LoadApproxModule(module);
		std::cout << "Getting " << moduleNames[key] << std::endl;
		pDLLApprox = ApproxManager::GetInstance().GetApproxServer().GetModule(module).GetApproximator();
		
	} 
	catch(std::runtime_error & Ex)
	{
		const char * descr;
		if(type == APR_DG_PRSIM) 
			descr = " dg field type for prism mesh!";
		else if(type == APR_STD_PRISM)
			descr = " std field type for prism mesh!";
		else if(type == APR_STD_HYBRID)
			descr = " std filed type for hybrid mesh!";
		else
			descr = " unknown field type!";
		
		std::cerr << Ex.what() << std::endl;
		std::cerr << "Can't set active module for" << descr << std::endl;

		return(-1);
	}
	catch(...)
	{
		std::cerr << "Something is wrong!\n";
		return(-1);
	}
	return(1);
}

int init_mesh(int Control, char * Filename)
{
	if(is_init()) {
		return pDLLApprox->init_mesh(Control,Filename);
	}
	return(-1);
}

int free_mesh(int Mesh_id)
{
	if(is_init()) {
		return pDLLApprox->free_mesh(Mesh_id);
	}
	return(-1);	
}

int fun_call get_nmno(int Mesh_id)
{
	return pDLLApprox->get_nmno(Mesh_id);
}

int fun_call get_nmed(int Mesh_id)
{
	return pDLLApprox->get_nmed(Mesh_id);
}

int fun_call get_nmfa(int Mesh_id)
{
	return pDLLApprox->get_nmfa(Mesh_id);
}

int fun_call get_nmel(int Mesh_id)
{
	return pDLLApprox->get_nmel(Mesh_id);
}

int fun_call get_el_status(int Mesh_id,int El)
{
	return pDLLApprox->get_el_status(Mesh_id,El);
}

int fun_call get_next_act_elem(int Mesh_id,int Fa)
{
	return pDLLApprox->get_next_act_elem(Mesh_id,Fa);
}

int fun_call get_next_elem(int Mesh_id,int Fa)
{
	return pDLLApprox->get_next_elem(Mesh_id,Fa);
}

int fun_call get_next_face(int Mesh_id,int Fa)
{
	return pDLLApprox->get_next_face(Mesh_id,Fa);
}

int fun_call get_next_edge(int Mesh_id,int Ed)
{
	return pDLLApprox->get_next_edge(Mesh_id,Ed);
}

int fun_call get_next_node(int Mesh_id,int Node)
{
	return pDLLApprox->get_next_node(Mesh_id,Node);
}



int fun_call get_el_type(int Mesh_id,int El)
{
	return pDLLApprox->get_el_type(Mesh_id,El);
}

int fun_call get_el_mate(int Mesh_id,int El)
{
	return pDLLApprox->get_el_mate(Mesh_id,El);
}

int fun_call get_el_faces(int Mesh_id,int El,int* Faces, int* Orient)
{
	return pDLLApprox->get_el_faces(Mesh_id,El,Faces,Orient);
}

int fun_call get_el_fam(int Mesh_id,int El,int* Elsons,int* Type)
{
	return pDLLApprox->get_el_fam(Mesh_id,El,Elsons,Type);
}

int  fun_call get_el_struct(int Mesh_id,int El_id,int *El_struct)
{
	return pDLLApprox->get_el_struct(Mesh_id,El_id,El_struct);
}

int  fun_call get_el_node_coor(int Mesh_id,int El_id,int *Nodes,double *Node_coor)
{
	return pDLLApprox->get_el_node_coor(Mesh_id,El_id,Nodes,Node_coor);
}
	
int fun_call get_fa_status(int Mesh_id,int Fa)
{
	return pDLLApprox->get_fa_status(Mesh_id,Fa);
}

int fun_call get_fa_type(int Mesh_id,int Fa)
{
	return pDLLApprox->get_fa_type(Mesh_id,Fa);
}

int fun_call get_fa_sub_bnd(int Mesh_id,int Fa)
{
	return pDLLApprox->get_fa_sub_bnd(Mesh_id,Fa);
}

int fun_call get_fa_node_coor(int Mesh_id,int Fa,int* Nodes,double* Coords)
{
	return pDLLApprox->get_fa_node_coor(Mesh_id,Fa,Nodes,Coords);
}


void fun_call get_fa_neig(int Mesh_id, int Fa, int* Fa_neig,
			int* Neig_sides, int* Node_shift, int* Diff_gen,
			double* Acceff, double* Bcceff)
{
	pDLLApprox->get_fa_neig(Mesh_id,Fa,Fa_neig,Neig_sides,
		Node_shift,Diff_gen,Acceff,Bcceff);
}

int fun_call get_face_edges(int Mesh_id,int Fa,int *Edges,int *Orient)
{
	return pDLLApprox->get_face_edges(Mesh_id,Fa,Edges,Orient);
}

int  fun_call get_face_struct(int Mesh_id,int Fa,int *Fa_struct)
{
	return pDLLApprox->get_face_struct(Mesh_id,Fa,Fa_struct);
}

int fun_call get_edge_status(int Mesh_id,int Ed)
{
	return pDLLApprox->get_edge_status(Mesh_id,Ed);
}

int fun_call get_edge_nodes(int Mesh_id,int Ed,int *Ed_nodes)
{
	return pDLLApprox->get_edge_nodes(Mesh_id,Ed,Ed_nodes);
}

int fun_call get_edge_struct(int Mesh_id,int Ed,int *Ed_struct)
{
	return pDLLApprox->get_edge_struct(Mesh_id,Ed,Ed_struct);
}

int fun_call get_node_status(int Mesh_id,int Node)
{
	return pDLLApprox->get_node_status(Mesh_id,Node);
}

int fun_call get_node_coor(int Mesh_id,int Node,double *Coor)
{
	return pDLLApprox->get_node_coor(Mesh_id,Node,Coor);
}

// Approximation interface
int std_call init_field(char Field_type,int Mesh_id,const char *Filename)
{
	return pDLLApprox->init_field(Field_type,Mesh_id,Filename);
}

int free_field(int Field_id)
{
	return pDLLApprox->free_field(Field_id);
}

int fun_call get_nreq(int Field_id)
{
	return pDLLApprox->get_nreq(Field_id);
}

int fun_call get_nr_sol(int Field_id)
{
	return pDLLApprox->get_nr_sol(Field_id);
}

int  fun_call get_base_type(int Field_id)
{
	return pDLLApprox->get_base_type(Field_id);
}

int  fun_call get_el_pdeg(int Field_id,int El,int* Pdeg_vec)
{
	return pDLLApprox->get_el_pdeg(Field_id,El,Pdeg_vec);
}

int fun_call get_element_dofs(int Field_id,int El_id,int Vec_id,double* Dofs)
{
	return pDLLApprox->get_element_dofs(Field_id,El_id,Vec_id,Dofs);
}

int fun_call create_constr_data(int Field_id)
{
	return pDLLApprox->create_constr_data(Field_id);
}

double apr_elem_calc_3D_mod(
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
		) 
{
	return pDLLApprox->apr_elem_calc(Control,Nreq,Pdeg_vec,Base_type,Eta,Node_coor,
			Sol_dofs,Base_phi,Base_dphix,Base_dphiy,Base_dphiz,Xcoor,Sol,
			Dsolx,Dsoly,Dsolz,Vec_nor);
}


