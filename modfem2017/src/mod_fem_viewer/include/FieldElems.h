#ifndef _FIELD_ELEMS_H_
#define _FIELD_ELEMS_H_

/*** CONSTANTS ***/

#define APC_MAX_NUM_FIELD   10   /* maximal number of fields */

/*** DATA TYPES ***/

/* solution degrees of freedom structure (parallel to element structure) */
typedef struct {
		int pdeg;   /* indicator of polynomial degrees */
					/* for different types of elements: */
					/* TETRA - pdeg = pdegxyz one for all directions */
					/* PRISM - pdeg = 100*pdegz + pdegxy */
					/* BRICK - pdeg = 100*pdegz + pdegy*10 + pdegx */
					/* (pdeg = -1 - sol_* vectors not allocated) */
		double *vec_dof_1;   /* pointer to array of field dofs */
		double *vec_dof_2;   /* pointer to array of field dofs */
		double *vec_dof_3;   /* pointer to array of field dofs */
} apt_dg_dof_ent;


/* discretization structure, including time and space discretization */
typedef struct {
		int     mesh_id;      /* identifier of the mesh associated with field */
		int     nreq;         /* number of components in solution vector */
		int     nr_sol;       /* number of solution vectors for each element */
		int     uniform;      /* whether the pdeg is the same for all elements */
		int     base;	/* type of shape functions:  */
						/* 	1 (APC_TENSOR) - tensor product */
						/* 	2 (APC_COMPLETE) - complete polynomials */
						/* ! for both types the order of the first four dofs  */
						/* ! as returned from apr_shape_fun_3D is: 1,x,y,z ! */
		//  int     nr_dof_ent;   /* number of dof entities */
		apt_dg_dof_ent  *dof_ents;   /* pointer to structure with the vectors of dofs */
} apt_dg_field;

/* new structure */
typedef struct {
	//int node_id;
	double* vec_dof_1/*[3]*/; //Uwaga!
	double* vec_dof_2;
	double* vec_dof_3;
	int* constr;
} apt_std_dof_ent;
/* new discretization structure, including time and space discretization */
typedef struct {
	int mesh_id;
	int nreq;
	int nr_sol;
	int base;
	int pdeg;
	int constr;
	int uniform;
	int nr_nodes;
	int nmno_nodes;
	apt_std_dof_ent *dof_ents;
} apt_std_field;

#endif /* _FIELD_ELEMS_H_
*/