/************************************************************************
File aph_dg_prism.h - internal information for discontinuous Galerkin
                          approximation module for prismatic elements       

Contains:
  - constants
  - data types 
  - global variables (for the whole module)
  - headers for internal functions:                     
  apr_select_field - to select the proper field   
  apr_loc_loc - to compute local coordinates within an element,
	given local coordinates within an element of the same family
  apr_sol_xglob - to return the solution at a point with global
	coordinates specified

------------------------------  			
History:        
        02.2002 - Krzysztof Banas, initial version
*************************************************************************/

#ifndef _aph_dg_prism_
#define _aph_dg_prism_


/*** CONSTANTS ***/
/* maximal number of equations */
#define APC_MAXEQ 1

#define APC_MAX_NUM_FIELD   10   /* maximal number of fields */

/*** DATA TYPES ***/

/* solution degrees of freedom structure (parallel to element structure) */
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
} apt_dof_ent;


/* discretization structure, including time and space discretization */
typedef struct {
  int     mesh_id;      /* identifier of the mesh associated with field */
  int     nreq;         /* number of components in solution vector */
  int     nr_sol;       /* number of solution vectors for each element */
  int     base;		/* type of shape functions:  */
			/* 	1 (APC_TENSOR) - tensor product */
			/* 	2 (APC_COMPLETE) - complete polynomials */
			/* ! for both types the order of the first four dofs  */
			/* ! as returned from apr_shape_fun_3D is: 1,x,y,z ! */
  int uniform;          /* 0 - non-uniform field */
                        /* 1 - uniform field (all pdeg equal) */
  //  int     nr_dof_ent;   /* number of dof entities */
  apt_dof_ent  *dof_ents;   /* pointer to structure with the vectors of dofs */
} apt_field;

/*** GLOBAL VARIABLES for the whole module ***/

extern int       apv_nr_fields;     /* the number of fields in the problem */
extern int       apv_cur_field_id;              /* ID of the current field */
extern apt_field  apv_fields[APC_MAX_NUM_FIELD];        /* array of fields */


/*** FUNCTIONS DECLARATIONS - headers for internal functions ***/
/**--------------------------------------------------------
  apr_select_field - to select the proper field   
---------------------------------------------------------*/
extern apt_field* apr_select_field( /* returns: pointer to the selected field */
  int Field_id  /* in: field ID */
  );





/*
In file aps_gauss_util.c - utility routines for Gauss-Legendre integration 
*/

/**-------------------------------------------------------------------
  apr_gauss_init - to initialize Gauss-Legendre quadratures
--------------------------------------------------------------------*/
extern void apr_gauss_init();

/**-------------------------------------------------------------------
  apr_gauss_select - to select the proper quadrature
--------------------------------------------------------------------*/
extern void apr_gauss_select(
	int Type,	/* in: type of integration considered */
			/* 	1 - LINE INTEGRATION over [-1,1] */
			/* 	2 - GAUSSIAN QUADRATURES FOR TRIANGLES */
	int Order, 	/* in: order of integration == degree */
			/*     of polynomials integrated exactly */
        int *Ng,	/* out: number of gauss points */
	double **Xg_p, 	/* out: gauss points */
	double **Wg_p 	/* out: gauss weights */
	);



#endif
