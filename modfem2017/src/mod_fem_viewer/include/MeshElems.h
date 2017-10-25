#ifndef _FV_MSH_ELEMS_H_
#define _FV_MSH_ELEMS_H_


/* structure with mesh parameters */
typedef struct {
	int    nrno;        /* number of active nodes */
	int    nmno;        /* maximal index of initiated node */
	int    mxno;        /* maximal available index of node */
	int    pfno;        /* pointer to first free node space*/
	int    nred;        /* number of active edges */
	int    nmed;        /* maximal index of initiated edge */
	int    mxed;        /* maximal available index of edge */
	int    pfed;        /* pointer to first free edge space*/
	int    nrfa;        /* number of active faces */
	int    nmfa;        /* maximal index of initiated face */
	int    mxfa;        /* maximal available index of face */
	int    pffa;        /* pointer to first free face space */
	int    nrel;        /* number of active elements */
	int    nmel;        /* maximal index of initiated element */
	int    mxel;        /* maximal available index of element */
	int    pfel;        /* pointer to first free element space */
	int    maxgen;      /* maximum generation level for elements */
	int    maxgendiff;  /* maximum difference of generation levels */
	     			  /* for neighboring elements */
} mmt_meshp;

/* node structure */
typedef struct {
/*||begin||*/
	//int   ipid;        
/*||end||*/
	double   x;        /* x-coor of the node (or -1e11 for free spaces ) */
	double   y;        /* y-coor of the node (or next free space ID) */
	double   z;        /* z-coor of the node */
	int		 ipid;	   /* global identifier that include ownership info */ 
}	mmt_nodes;

/* edge structure */
typedef struct {
	int type;		/*  1 (MMC_EDGE) - active edge */
		            /* <0 - inactive edge, number of attempted divisions */
					/*  0 (MMC_FREE) - free space */
/*||begin||*/
	//int ipid;		/* global identifier that include ownership info */ 
/*||end||*/
	int node[2];	/* for active edges: node's IDs */
					/* for free spaces: node[0] - the next free space */
					/* for inactive edges: son's IDs */
} mmt_edges;

/* face structure */
typedef struct {
	int type;       /* type and status of the face: */
					/*   1 (MMC_EDGE) - 1D face for 2D elements, */
					/*   2 (MMC_QUAD) - quadrilateral */
					/*   3 (MMC_TRIA) - triangle, */
					/*   0 (MMC_FREE) - free space */
					/*  >0 - active, not refined */
					/*  <0 - inactive, refined */
/*||begin||*/
	//int ipid;       /* global identifier that include ownership info */ 
/*||end||*/
	int bc;         /* bc flag for the face (obtained from mesh generator and */
					/*   interpreted by the problem dependent module) */
					/* for free spaces: next free space's ID */
	int edge[4];	/* edges' IDs - in proper order */
					/*  >0 - edge with the same orientation as face */
					/*  <0 - edge with opposite orientation */
	int neig[2];	/* neighboring elements' global IDs */
					/* (first neighbor has the same orientation */
					/*  as the face and ordering of its nodes */
					/*  defines ordering of nodes on the face) */
					/*  >0 - ID of equal size neighbor */
					/*  -1 (MMC_BIG_NGB) - indicates big neighbor */
					/*   0 (MMC_BOUNDARY) - boundary (always second neighbor) */
					/*  -2 (MMC_SUB_BND) - inter-subdomain boundary */
	int *sons;	/* sons */
} mmt_faces;

/* element structure */
typedef struct {
	int    type;    /* 7 (MMC_TETRA) - tetrahedron */
		             /* 5 (MMC_PRISM) - prism */
					/* 6 (MMC_BRICK) - brick */
					/* 0 (MMC_FREE) - free space (face and sons not allocated) */
					/* >0 - active, not refined */
					/* <0 - inactive, refined */
/*||begin||*/
	//int    ipid;    /* global identifier that include ownership info */ 
/*||end||*/
	int    mate;    /* material ID */
					/* for free space: ID of the next free structure */
	int    fath;    /* parent's global ID */
					/* 0 (MMC_NO_FATH) - if the element has no parent */
	int    refi;    /* type of the last refinement */
					/* 0 (MMC_NOT_REF) - sons not allocated */
					/* 1 (MMC_REF_TYP_ISO) - isotropic refinement */
	int    *face;   /* faces' global IDs (sign=orientation)*/
	int    *sons;   /* children */
} mmt_elems;

/* structure with mesh data */
typedef struct {
	int     space_dim;        /* number of space dimensions */
	int       monitor;        /* printing flag */
	mmt_meshp    parm;        /* structure with mesh parameters */ 
	mmt_nodes   *node;        /* pointer to array of nodes */
	mmt_edges   *edge;        /* pointer to array of edges */
	mmt_faces   *face;        /* pointer to array of faces */
	mmt_elems   *elem;        /* pointer to array of elements */
} mmt_mesh;


/*||begin||*/
/* info about refinements - data exchanged between subdoamins */
typedef struct{
	  /* for each refined element (prism or hexahedron):
		 element_ID, element_type, new_node_ID (if any), 
		 new_edges_IDs, new_faces_IDs (horizontal+vertical), new_elems_IDs (sons)
		 1 + 1 + [1] + 3..6 + 10..12 + 8 = 23..29
	  */
	  int* elem;
	  /* for each refined face:
		 face_ID, face_type,  
		 new_node_ID (if any), new_edges_IDs, new_faces_IDs (sons) 
		 1 + 1 + [1] + 3..4 + 4 = 9..11
	  */
	  int* face;
	  /* for each refined edge:
		 edge_ID, new_node_ID, new_edges_IDs (sons) 
		 1 + 1 + 2 = 4
	  */
	  int* edge;
} mmt_ref_info;


#endif /* _FV_MSH_ELEMS_H_
*/