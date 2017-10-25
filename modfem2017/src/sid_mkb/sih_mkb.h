/************************************************************************
File sih_mkb.h - internal information for the interface module
   between the iterative block based Krylow solver and the finite element
   code (the module forms part of the fem code): definition of parameters,
   data types, global variables and external functions)


Contains declarations of data types, constants and global variables 

------------------------------  			
History:        
	02.2002 - Krzysztof Banas, initial version		
*************************************************************************/

#ifndef _sih_mkb_
#define _sih_mkb_

#include "sih_intf.h"

#ifdef __cplusplus
extern "C" {
#endif

/*** CONSTANTS ***/
#define SIC_MAX_NUM_LEV  20

#define SIC_PDEG_COARSE_HIGH        -4
#define SIC_PDEG_COARSE_LOW         -3
#define SIC_PDEG_COARSE_HIGH_ALL    -2
#define SIC_PDEG_COARSE_LOW_ALL     -1
#define SIC_PDEG_FINEST            -10

#define SIC_STD_APPROX 1
#define SIC_DG_APPROX  2

// other constants in include/sih_intf.h 

/*** DATA TYPES ***/

/* dof structure with data useful for creating flexible interfaces*/
/*  between FEM code and different solvers */
typedef struct{
  int dof_ent_type;  /* type of the associated FEM code (mesh) entity */
  int dof_ent_id;    /* ID of the associated FEM code (mesh) entity */
  int nr_int_ent; /* number of  integration entities providing SMs and LVs*/
  int l_int_ent_index[SIC_MAX_INT_PER_DOF]; 
                 /* list of integration entities providing SMs and LVs*/
  int block_id; /* ID for solver - used for renumbering */
  int nrdofs;      /* number of DOFs */
  int posglob;    /* position in a global load vector (and stiffness matrix) */
  int nrneig;     /* number of neighboring DOF structures */
  //* Two lists of neighbours - the order on lists may be different !!! */
  int l_neig[SIC_MAX_DOF_STR_NGB];  /* list of neighboring DOF structures */
  int l_neig_bl[SIC_MAX_DOF_STR_NGB]; /* list of IDs for solver (block IDs) */
}  sit_dof_struct;


/* solver data structure for a single level */
typedef struct {

  int nr_int_ent;   /* number of integration entities - entities that */
                    /* provide solver with stiffness matrices and load vectors*/
  int nr_dof_ent;   /* number of dof entities - mesh entities with which */
                    /* degrees of freedom are associated */
  int nr_dof_ent_internal ; /* number of internal (not ghost i.e. overlap dof entities) */
                /* for which stiffness matrix rows are assembled and processed */
  int pdeg_coarse;  /* degree of approximation for coarse levels */
  int nrdofs_glob;     /* the global number of degrees of freedom */
  int nrdofs_internal; /* the number of degrees of freedom for internal DOF entities */
  int block_size;   /* for stiffness matrices with constant block_size */
                    /* !!! makes block.nrdofs and block.posglob obsolete !!! */
  int max_dofs_int_ent; /* maximal number of dofs per integration entity, i.e. */
                       /* maximal size of the local stiffness matrix */
  int max_dofs_dof_ent; /* maximal number of dofs per dof entity */

/* arrays for assembling local stiffness matrices into global stiffness matrix*/
  int* l_int_ent_type; /*list of types of entities providing local SMs and LVs */
  int* l_int_ent_id; /* list of ID's of entities providing local SMs and LVs */
  /* data for element coloring - to enable concurrent assembly */
  int nr_colors_elems; /* elements and faces have the same l_int_ent_... arrays but ... */
  int* l_color_index_elems; /* indices in permuted l_int_ent_.. for subsequent colors */
  int nr_colors_faces; /* ... elements and faces have different l_colors_... arrays!!! */
  int* l_color_index_faces; /* indices in permuted l_int_ent_.. for subsequent colors */

  // data for assembly tables: for each local SM (stiffness matrix)
  // the positions in the global SM of first entries associated with all matrix blocks 
  // in local SM are stored (each matrix block corresponds to a pair of DOF entities,
  // i.e. a pair of DOF structures, i.e. - in the current implementation - a pair of 
  // DOF blocks) 
  // 1. assembly data is stored in a consecutive way for all integration entities, 
  //    according to the order in l_int_ent_... lists)
  // 2. positions in the global SM depend on the storage format, for each format
  //    different tables can be used (as the first implementation the coordinate (COO)
  //    format is used - it can be further used for other formats) 
  int nr_asse_blocks_all_int_ent; // the sum of the numbers of DOF entities for all 
                          // integration entities i.e. the length of the tables below
  int* asse_pos_first_dof_int_ent;
  // integration routines return stiffness matrices with the associated lists of 
  // DOF_entities - for each pair of DOF entities a single SM block is created
  // for each SM block its position in sparse storage arrays of the solver  
  // (like CRS or BCRS) is stored
  // (since DOF entities can have more than one associated DOF the whole storage scheme
  // is based on the notions of entities-structures-blocks, not single DOFs!!!)
  int* assembly_table; //

  sit_dof_struct *l_dof_struct; /* list of dof structures with data useful for */
                               /* creating flexible interfaces between FEM code*/
                               /* and different solvers */

  /* for each possible type of dof entity - its corresponding dof structure */
  int* l_dof_vert_to_struct;  
  int* l_dof_edge_to_struct;  
  int* l_dof_face_to_struct;  
  int* l_dof_elem_to_struct;
  /* dimensions of the above arrays */
  int max_dof_vert_id;
  int max_dof_edge_id;
  int max_dof_face_id;
  int max_dof_elem_id; 

/* renumbering arrays */
//  int* l_bl_to_struct;    /* block to structure correspondence */

} sit_levels;

/* definition of sit_solvers - data type for multi-level iterative solver */
typedef struct {
  int solver_id;           /* solver_id */
  int problem_id;          /* ID of the problem associated with the solver */
  int solver_type;  /* linear equations solver: */
		    /*	0 - direct solver */
		    /*	1 - GMRES */
		    /*	2 - multi-level GMRES */
                    /* 99 - MKB_CORE_NS_SUPG_SOLVER */
  int storage_type;
  int parallel;            /* parameter specifying sequential (SIC_SEQUENTIAL) */
                           /* or parallel (SIC_PARALLEL) execution */
  int nr_levels;	   /* number of levels in multi-level GMRES */
  int cur_level;           /* current level number in multi-level GMRES */
  sit_levels level[SIC_MAX_NUM_LEV];    /* array of solver data structures */
                           /* corresponding to different levels */
} sit_solvers;		    

/*** GLOBAL VARIABLES (for the solver module only) ***/

extern int   siv_nr_solvers;     /* the number of solvers in the problem */
extern int   siv_cur_solver_id;              /* ID of the current solver */
extern sit_solvers siv_solver[SIC_MAX_NUM_SOLV];     /* array of solvers */

#ifdef __cplusplus
}
#endif


#endif

