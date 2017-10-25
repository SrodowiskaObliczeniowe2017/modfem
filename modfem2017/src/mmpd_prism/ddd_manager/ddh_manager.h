/************************************************************************
File

Contains definitions of routines:   

------------------------------  			
History:        
	02.2002 - Krzysztof Banas, initial version		
*************************************************************************/

#ifndef _ddh_manager_
#define _ddh_manager_


/**--------------------------------------------------------
  ddr_vert_elems - to prepare lists of all elements (active and
                   inactive) to which belong subsequent nodes
---------------------------------------------------------*/
extern int ddr_vert_elems(
  int Mesh_id,      /* in: mesh ID */
  int Gen_lev,       /* in: generation level as basis for decomposition */
  int **Vert_elems  /* out: list of lists of vertices' elements */
  );

/**--------------------------------------------------------
  ddr_create_overlap
---------------------------------------------------------*/
extern int ddr_create_overlap(/* returns: >=0 - success code, <0 - error code */
  int Mesh_id, /* in: mesh ID */
  int Ovl_size,
  int Control  /* in: variable to select different kinds of overlap: */
	       /*        1 - overlap for dofs within elements */
               /*        2 - overlap for dofs at nodes */
  );

/**--------------------------------------------------------
  ddr_add_sons_ovl - to add ancestors to the overlap in a recursive manner
---------------------------------------------------------*/
extern int ddr_add_sons_ovl(
  int Mesh_id, /* in: mesh ID */
  int El       /* in: element ID */
  );

/**--------------------------------------------------------
ddr_chk_list - to check whether Num is on the list List
	with length Ll
---------------------------------------------------------*/
extern int ddr_chk_list(	/* returns: */
			/* >0 - position on the list */
            		/* 0 - not found on the list */
	int Num, 	/* number to be checked */
	int* List, 	/* list of numbers */
	int Ll		/* length of the list */
	);

/**--------------------------------------------------------
ddr_put_list - to put Num on the list List with length Ll 
	(filled with numbers and zeros at the end)
---------------------------------------------------------*/
extern int ddr_put_list( /* returns*/
		/*  >0 - position already occupied on the list */
             	/*  <0 - position at which put on the list */
            	/*   0 - list full, not found on the list */
	int Num, 	/* in: number to put on the list */
	int* List, 	/* in: list */
	int Ll		/* in: total list's lengths */
	);


/**--------------------------------------------------------
  ddr_create_patch - to create a patch surrounding element
---------------------------------------------------------*/
int ddr_create_patch(
  int Mesh_id, 
  int El_id, 
  int Gen_lev,       /* in: generation level as basis for decomposition */
  int List_length, 
  int* List_el, 
  int* List_face_int, 
  int* List_face_bnd,
  int* List_edge_int, 
  int* List_edge_bnd, 
  int* List_vert_int, 
  int* List_vert_bnd
  );

#endif
