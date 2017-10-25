#ifndef _APPROX_MODULE_H_
#define _APPROX_MODULE_H_

#include "../../include/mod_fem_viewer.h"

#if defined(__cplusplus)
extern "C" {
#endif

#define APC_MAX_NUM_FIELD 10

#ifdef WIN32
# ifndef fun_call
#   define fun_call _fastcall
# endif
# ifndef std_call
#   define std_call _stdcall
# endif
#else
# define fun_call
# define std_call
#endif


typedef enum {
	APR_NONE		= 0x0,
	APR_MESH	    = 0x1,
	APR_DG_PRSIM	= 0x2,
	APR_STD_PRISM	= 0x4,
	APR_STD_HYBRID  = 0x8,
} apr_type;

enum mod_approx
{
	approx_dg_std_int,
	approx_dg_ext,
	approx_std_ext,
	approx_std_hibrid_ext,
};


/*** FUNCTIONS DECLARATIONS - headers for internal functions ***/

/* init - function for checking module's presentation */
int init_approx_modules();

/* destroy - function to destroying manager and do cleanup */
void destroy();

/* set_approx - set current approximation */
int set_approx_module(apr_type type);

/* init_mesh - Function to init mesh */
int init_mesh(int Control, char * Filename);

/* free_mesh - release mesh */
int free_mesh(int Mesh_id);

int  fun_call get_nmno(int Mesh_id);

int  fun_call get_nmed(int Mesh_id);

int  fun_call get_nmfa(int Mesh_id);

int  fun_call get_nmel(int Mesh_id);

int  fun_call get_next_act_elem(int Mesh_id,int El);

int  fun_call get_next_elem(int Mesh_id,int El);

int  fun_call get_next_face(int Mesh_id,int Fa);

int  fun_call get_next_edge(int Mesh_id,int Ed);
int  fun_call get_next_node(int Mesh_id,int Node);
int  fun_call get_el_status(int Mesh_id,int el);
int  fun_call get_el_type(int Mesh_id,int El);
int  fun_call get_el_mate(int Mesh_id,int El);
int  fun_call get_el_faces(int Mesh_id,int El,int* Faces, int* Orient);
int  fun_call get_el_fam(int Mesh_id,int El,int* Elsons,int* Type);
int  fun_call get_el_struct(int Mesh_id,int El_id,int *El_struct);

int  fun_call get_el_node_coor(int Mesh_id,int El,int *Nodes,double *Xcoor);

int  fun_call get_fa_status(int Mesh_id,int Fa);
int  fun_call get_fa_type(int Mesh_id,int Fa);
int  fun_call get_fa_sub_bnd(int Mesh_id,int Fa);
int  fun_call get_fa_node_coor(int Mesh_id,int Fa,int* Nodes,double* Coords);
void fun_call get_fa_neig(int Mesh_id, int Fa, int* Fa_neig,
			int* Neig_sides, int* Node_shift, int* Diff_gen,
			double* Acceff, double* Bcceff);

int fun_call get_face_edges(int Mesh_id,int Fa,int *Edges,int *Orient);

int  fun_call get_face_struct(int Mesh_id,int fa,int *Fa_struct);
int  fun_call get_edge_status(int Mesh_id,int Ed);
int  fun_call get_edge_nodes(int Mesh_id,int Ed,int *Ed_nodes);
int  fun_call get_edge_struct(int Mesh_id,int Ed,int *Ed_struct);
int  fun_call get_node_status(int Mesh_id,int Node);
int  fun_call get_node_coor(int Mesh_id,int Node,double *Coor);

int  std_call init_field(char Field_type,int Mesh_id,const char *Filename);
int  free_field(int Field_id);

int  fun_call get_nreq(int Field_id);

int fun_call get_nr_sol(int Field_id);

int  fun_call get_base_type(int Field_id);

int  fun_call get_el_pdeg(int Field_id,int El,int* Pdeg_vec);

/*------------------------------------------------------------
  get_element_dofs - to return the list of standard dofs for a given element
                        each dof corresponds to a standard shape function
------------------------------------------------------------*/
int  fun_call get_element_dofs( /* returns: >0 - the number of dofs
                                         <0 - error code */
  int Field_id,     /* in: approximation field ID  */
  int El_id,        /* in: element ID */
  int Vect_id,       /* in: vector ID in case of multiple solution vectors */
  double *El_dofs_std     /* out: the list of values of element dofs */
  );

/*---------------------------------------------------------
  apr_select_field - to sel ect the proper field
---------------------------------------------------------*/
//apt_field_dg* apr_select_field_dg( /* returns: pointer to the selected field */
//  int Field_id  /* in: field ID */
//  );

/*---------------------------------------------------------
  apr_select_field - to select the proper field
---------------------------------------------------------*/
//apt_field_std* apr_select_field_std( /* returns: pointer to the selected field */
//  int Field_id  /* in: field ID */
//  );

/*------------------------------------------------------------
  create_constr_data - to create and fill constraint arrays
                   (assuming approximation and mesh data agree)
------------------------------------------------------------*/
int fun_call create_constr_data( /* returns: >=0 - success code, <0 - error code */
  int Field_id			/* in: approximation field ID  */
  );



/*------------------------------------------------------------
  apr_get_constr_data - to return constraints data for a node (dof entity) 
------------------------------------------------------------*/
//int apr_get_constr_data(
//                  /* returns: >=0 - success code, <0 - error code */
//  int Field_id,   /* in: approximation field ID  */
//  int Nvert,      /* in: vertex (i.e. node, i.e. dof entity) ID  */
//  int *Constr     /* out: table with constraints data; */
//  /* constr[0] - number of constraints (2 - mid-edge node, 4 - mid-side node */
//  );

/*------------------------------------------------------------
  apr_get_el_constr_data - to return the number and the list of element's real 
                           nodes, with the corresponding constraint coefficients
------------------------------------------------------------*/
//int apr_get_el_constr_data(
//                  /* returns: >=0 - success code, <0 - error code */
//  int Field_id,   /* in: approximation field ID  */
//  int El_id,        /* in: element ID */
//  int* Nodes,  	  /* out: list of vertex node IDs */
//       		  /*	(Nodes[0] - number of nodes) */
//  int* Nr_constr, /* out: list with the numbers of constraints for each vertex*/ 
//  int* Constr_id, /* out: ID's of parent (constraining) nodes */
//  double* Constr_val  /* out: the corresponding constraint coefficients */
//  );

/*------------------------------------------------------------------
apr_elem_calc_3D - to perform element calculations (to provide data on
	coordinates, solution, shape functions, etc. for a given point
	inside element (given local coordinates Eta[i]);
	for geometrically multi-linear or linear 3D elements
-------------------------------------------------------------------*/
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
	);


   

#if defined(__cplusplus)
}
#endif

#endif /* _APPROX_MANAGER_H_
*/
