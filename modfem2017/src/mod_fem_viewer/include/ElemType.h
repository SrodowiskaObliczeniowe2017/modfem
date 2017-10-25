#ifndef _ELEM_TYPE_H_
#define _ELEM_TYPE_H_

namespace FemViewer {
namespace ELEMS {
	typedef unsigned int uint;

	typedef enum /*EL_TYPE*/{
		VERTEX = 0, 
		EDGE,
		FACE,
		PRIZM, 
		FIELD,
		INVALID_ELM = 999,
	} ElemType;

	typedef struct {
			float x, y, z;
		} Vertex;

	typedef struct {
			int edge_type;
			uint indexes[2];
		} Edge;

	typedef enum {
			TRI_FACE	= 3,
			QUAD_FACE	= 4,
		} FaceType;

	typedef struct {
			FaceType face_type;
			int bound_cond;
			int shift_pos;
			int neighbour;
			int index_oriented_edges[4];
		} Face;

	typedef struct {
			uint nr_of_faces;	// numer of faces constituing the element
			uint mat_id;		// material id
			uint parent_id;		// parent id
			uint refine_type;	// type of refinement 0 - no refinement
			int face_idx[5];
		} Element;

	typedef struct {
			int	   pdeg;     /* indicator of polynomial degrees */
							/* for different types of elements: */
							/* TETRA - pdeg = pdegxyz one for all directions */
							/* PRISM - pdeg = 100*pdegz + pdegxy */
							/* BRICK - pdeg = 100*pdegz + pdegy*10 + pdegx */
							/* (pdeg = -1 - sol_* vectors not allocated) */
			double *vec_dof_1;   /* pointer to array of field dofs */
			double *vec_dof_2;   /* pointer to array of field dofs */
			double *vec_dof_3;   /* pointer to array of field dofs */
		} apt_dg_dof_ent;

	typedef enum {
		APC_TENSOR		= 0,
		APC_COMPLETE	= 1,	
	} BaseFunType;


	/* discretization structure, including time and space discretization */
	typedef struct {
			int     mesh_id;      /* identifier of the mesh associated with field */
			int		n_els;		  /* number of elemsnts in asociated mesh */
			int     nreq;         /* number of components in solution vector */
			int     nr_sol;       /* number of solution vectors for each element */
			int     pdeg_coarse;  /* degree of approximation index for coarse grid */
			BaseFunType     base;		/* type of shape functions:  */
			/* 	1 (APC_TENSOR) - tensor product */
			/* 	2 (APC_COMPLETE) - complete polynomials */
			/* ! for both types the order of the first four dofs  */
			/* ! as returned from apr_shape_fun_3D is: 1,x,y,z ! */
			//  int     nr_dof_ent;   /* number of dof entities */
			apt_dg_dof_ent  *dof_ents;   /* pointer to structure with the vectors of dofs */
		} ElField;

} // end of ELEMS
} // end of FemViewer


#endif /* _ELEM_TYPE_H_
*/