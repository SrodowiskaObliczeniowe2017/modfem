/************************************************************************
File

Contains definitions of routines:   

------------------------------  			
History:        
	02.2002 - Krzysztof Banas, initial version		
*************************************************************************/

#ifndef _ddh_front_
#define _ddh_front_

#define DDC_MAX_NGB_SUBD 10


/**--------------------------------------------------------
  dds_prepare_nodes
---------------------------------------------------------*/
extern int dds_prepare_nodes(
  int Mesh_id,      /* in: mesh ID */
  int Level_id      /* in: generation level as basis for decomposition */
  );

/**--------------------------------------------------------
  ddr_create_subdomain - to create a subdomain
---------------------------------------------------------*/
extern int ddr_create_subdomain (
  int Mesh_id,	     /* in: mesh ID */
  int Level_id,      /* in: generation level as basis for decomposition */
  int Nr_sub,
  int Sub_id,
  int* Nr_sub_elem,  /* in/out: requested and realized number of elements */
                     /*         in the subdoamin */
  int** L_sub_elem    /* out: list of subdomain elements */
  );

/**--------------------------------------------------------
choose_front - to pick a node from the front according
		to prescribed weights
---------------------------------------------------------*/
extern int choose_front(
  int Mesh_id	     /* in: mesh ID */
  );

/**--------------------------------------------------------
  upd_front - to update the front after adding an element to
	the subdoamin
---------------------------------------------------------*/
extern void upd_front(
  int Mesh_id,	     /* in: mesh ID */
  int Node,
  int Sub_id,
  int Ilev
  );


/**--------------------------------------------------------
put_sons - to put sons of element Nel together with their nodes
	on the lists of subdomain elements and nodes
---------------------------------------------------------*/
extern void put_sons(
  int Mesh_id,	     /* in: mesh ID */
  int El, 
  int Sub_id,
  int *Sub_size_p ,     /* in/out: subdomain size to be updated */
  int *Sub_size_tot_p   /* in/out: subdomain size to be updated */
  );

#endif
