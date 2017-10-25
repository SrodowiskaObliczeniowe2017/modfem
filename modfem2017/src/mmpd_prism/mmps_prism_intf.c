/************************************************************************
File mmps_prism_intf.c - implementation of parallel interface routines for 
                        meshes of prismatic elements 

Contains definitions of interface routines:   
  mmpr_init_mesh - to initialize the parallel mesh data structure 

  mmpr_el_owner - to get an owning processor identifier of an element 
  mmpr_el_set_owner - to set an owning processor identifier of an element
  mmpr_fa_owner - to get an owning processor identifier of a face
  mmpr_fa_set_owner - to set an owning processor identifier of a face
  mmpr_ed_owner - to get an owning processor identifier of an edge 
  mmpr_ed_set_owner - to set an owning processor identifier of an edge
  mmpr_ve_owner - to get an owning processor identifier of a vertex
  mmpr_ve_set_owner - to set an owning processor identifier of a vertex
  mmpr_el_id_at_owner - to get an local identifier of an element 
  mmpr_el_set_id_at_owner - to set an local identifier of an element
  mmpr_fa_id_at_owner - to get an local identifier of a face
  mmpr_fa_set_id_at_owner - to set an local identifier of a face
  mmpr_ed_id_at_owner - to get an local identifier of an edge 
  mmpr_ed_set_id_at_owner - to set an local identifier of an edge
  mmr_ve_id_at_owner - to get an local identifier of a vertex
  mmpr_ve_set_id_at_owner - to set an local identifier of a vertex

  mmpr_init_ref - to initialize the parallel process of refinement
  mmpr_refine - to refine an element or the whole mesh in parallel
  mmpr_derefine - to derefine an element or the whole mesh in parallel
  mmpr_final_ref - to finalize the process of parallel refinement 
  mmpr_free_mesh - to free space allocated for parallel mesh data structure

------------------------------  			
History:        
	02.2002 - Krzysztof Banas, initial version		
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

/* interface for the mesh manipulation module */
#include "mmh_intf.h"	

/* interface for the mesh manipulation module */
#include "mmph_intf.h"	

/* mesh manipulation data structure and headers for internal routines */
#include "./mmph_prism.h"

/* primitive utilities - for the time being */
#include "ddd_manager/ddh_manager.h"

int       mmpv_nr_meshes;   /* the number of meshes in the problem */
int       mmpv_cur_mesh_id;              /* ID of the current mesh */
mmpt_mesh  mmpv_meshes[MMC_MAX_NUM_MESH];        /* array of meshes */


int mmpv_nr_sub;
int mmpv_my_proc_id;

/*------------------------------------------------------------
  mmpr_init_mesh - to initialize the parallel mesh data structure
------------------------------------------------------------*/
int mmpr_init_mesh(  /* returns: >0 - Mesh ID, <0 - error code */
  int Control,
  int Mesh_id,
  int Nr_proc,
  int My_proc_id
  )
{

  mmpt_mesh* mesh;
  /* auxiliary variables */
  int iaux, ned;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  mmpv_nr_sub = Nr_proc;

  mmpv_my_proc_id = My_proc_id;

  /* increase the counter for meshes */
  mmpv_nr_meshes++;

  /* set the current mesh number */
  mmpv_cur_mesh_id = mmpv_nr_meshes;

  if(mmpv_cur_mesh_id != Mesh_id){
    printf("Wrong mesh_ids in initialization of parallel mesh!!! Exiting.\n");
    exit(-1);
  } 

/* select the proper mesh data structure */
  mesh = mmpr_select_mesh(Mesh_id);

  iaux = mmr_get_max_elem_max(Mesh_id);
/* allocate space for nodes' structures */
  mesh->elem =
       (mmpt_elems*)malloc((iaux+1)*sizeof(mmpt_elems));

/* report */
#ifdef DEBUG_MMM
    printf("Mesh %d:\n",Mesh_id);
    printf("Elems   : allocated %d structures for parallel management\n", iaux);
#endif


  iaux = mmr_get_max_face_max(Mesh_id);
/* allocate space for nodes' structures */
  mesh->face =
       (mmpt_faces*)malloc((iaux+1)*sizeof(mmpt_faces));

/* report */
#ifdef DEBUG_MMM
    printf("Mesh %d:\n",Mesh_id);
    printf("Faces   : allocated %d structures for parallel management\n", iaux);
#endif


  iaux = mmr_get_max_edge_max(Mesh_id);
/* allocate space for nodes' structures */
  mesh->edge =
       (mmpt_edges*)malloc((iaux+1)*sizeof(mmpt_edges));

/* report */
#ifdef DEBUG_MMM
    printf("Mesh %d:\n",Mesh_id);
    printf("Edges   : allocated %d structures for parallel management\n", iaux);
#endif


  iaux = mmr_get_max_node_max(Mesh_id);
/* allocate space for nodes' structures */
  mesh->node =
       (mmpt_nodes*)malloc((iaux+1)*sizeof(mmpt_nodes));

/* report */
#ifdef DEBUG_MMM
  printf("Mesh %d:\n",Mesh_id);
  printf("Nodes   : allocated %d structures for parallel management\n", iaux);
#endif

  mesh->ref_loc = NULL;
  mesh->ref_ali = NULL;
  mesh->del_loc = NULL;
  mesh->mesh_ovls = NULL;

  mmpr_create_subdomains(Mesh_id, Control); 

  return(mmpv_cur_mesh_id);
}


/*---------------------------------------------------------
  mmpr_select_mesh - to select the proper mesh   
---------------------------------------------------------*/
mmpt_mesh* mmpr_select_mesh( /* returns pointer to the chosen mesh */
			   /* to avoid errors if input is not valid */
			   /* it returns the pointer to the current mesh */
  int Mesh_id    /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  )
{

  /* select the proper mesh from the array of meshes */
  if( Mesh_id == MMC_CUR_MESH_ID ) {
    return(&mmpv_meshes[mmpv_cur_mesh_id-1]);
  }
  else if( Mesh_id>0 && Mesh_id<=mmpv_nr_meshes ) {
    return(&mmpv_meshes[Mesh_id-1]);
  }
  else {
    return(&mmpv_meshes[mmpv_cur_mesh_id-1]);
    /* alternative:  return(NULL);   */
  }
}


/*--------------------------------------------------------------------------
  mmpr_el_owner - to get an owning processor identifier of an element 
---------------------------------------------------------------------------*/
int mmpr_el_owner(/* returns: >0 -an owning processor identifier of an element */
	         /*          <0 - error code */ 
  int Mesh_id,	 /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El         /* in: global (within the subdomain) element ID */
  )
{
/* auxiliary variables */
  mmpt_mesh* mesh;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh = mmpr_select_mesh(Mesh_id);

  //if(mesh->elem[El].owner==0) return(mmpv_my_proc_id);
  return(mesh->elem[El].owner);
}

/*--------------------------------------------------------------------------
  mmpr_el_set_owner - to set an owning processor identifier of an element
---------------------------------------------------------------------------*/
int mmpr_el_set_owner( /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El,	/* in: global (within the subdomain) element ID */
  int Owner /* in: owning processor identifier of the element */
  )
{
/* auxiliary variables */
  mmpt_mesh* mesh;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh = mmpr_select_mesh(Mesh_id);
  
  mesh->elem[El].owner = Owner;

  return(1);
}

/*--------------------------------------------------------------------------
  mmpr_fa_owner - to get an owning processor identifier of a face
---------------------------------------------------------------------------*/
int mmpr_fa_owner(  /* returns: >0 - an owning processor identifier of a face */
	          /*          <0 - error code */ 
  int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa        /* in: global (within the subdomain)  ID */
  )
{
/* auxiliary variables */
  mmpt_mesh* mesh;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh = mmpr_select_mesh(Mesh_id);
  
  //if(mesh->face[Fa].owner==0) return(mmpv_my_proc_id);
  return(mesh->face[Fa].owner);
}

/*--------------------------------------------------------------------------
  mmpr_fa_set_owner - to set an owning processor identifier of a face
---------------------------------------------------------------------------*/
int mmpr_fa_set_owner(  /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa,       /* in: global (within the subdomain) face ID */
  int Owner      /* in: owning processor identifier of the face */
  )
{
/* auxiliary variables */
  mmpt_mesh* mesh;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh = mmpr_select_mesh(Mesh_id);
  
  mesh->face[Fa].owner = Owner;

  return(1);
}

/*--------------------------------------------------------------------------
  mmpr_ed_owner - to get an owning processor identifier of an edge 
---------------------------------------------------------------------------*/
int mmpr_ed_owner(  /* returns: >0 an owning processor identifier of an edge */
	          /*          <0 - error code */ 
  int Mesh_id,    /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ed          /* in: global (within the subdomain) edge ID */
  )
{
/* auxiliary variables */
  mmpt_mesh* mesh;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh = mmpr_select_mesh(Mesh_id);
  
  //if(mesh->edge[Ed].owner==0) return(mmpv_my_proc_id);
  return(mesh->edge[Ed].owner);

}

/*--------------------------------------------------------------------------
  mmpr_ed_set_owner - to set an owning processor identifier of an edge
---------------------------------------------------------------------------*/
int mmpr_ed_set_owner(  /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,  /* in: mesh ID or 0 (MMPC_CUR_MESH_ID) for the current mesh */
  int Ed,       /* in: global (within the subdomain) edge ID */
  int Owner      /* in: owning processor identifier of the edge */
  )
{
/* auxiliary variables */
  mmpt_mesh* mesh;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh = mmpr_select_mesh(Mesh_id);
  
  mesh->edge[Ed].owner = Owner;

  return(1);
}

/*--------------------------------------------------------------------------
  mmpr_ve_owner - to get an owning processor identifier of a vertex
---------------------------------------------------------------------------*/
int mmpr_ve_owner( /* returns: >0 - an owning processor identifier of a vertex */
	         /*          <0 - error code */ 
  int Mesh_id,  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ve        /* in: global (within the subdomain) vertex ID */
  )
{
/* auxiliary variables */
  mmpt_mesh* mesh;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh = mmpr_select_mesh(Mesh_id);
  
  //if(mesh->node[Ve].owner==0) return(mmpv_my_proc_id);
  return(mesh->node[Ve].owner);

}

/*--------------------------------------------------------------------------
  mmpr_ve_set_owner - to set an owning processor identifier of a vertex
---------------------------------------------------------------------------*/
int mmpr_ve_set_owner(  /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ve,       /* in: global (within the subdomain) vertex ID */
  int Owner      /* in: owning processor identifier of the vertex */
  )
{
/* auxiliary variables */
  mmpt_mesh* mesh;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh = mmpr_select_mesh(Mesh_id);
  
  mesh->node[Ve].owner = Owner;

  return(1);
}



/*--------------------------------------------------------------------------
  mmpr_el_id_at_owner - to get an local identifier of an element 
---------------------------------------------------------------------------*/
int mmpr_el_id_at_owner( /* returns: >0 - an local identifier of an element */
	         /*          <0 - error code */ 
  int Mesh_id,	 /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El         /* in: global (within the subdomain) element ID */
  )
{
/* auxiliary variables */
  mmpt_mesh* mesh;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh = mmpr_select_mesh(Mesh_id);
  
  return(mesh->elem[El].id_at_owner);
}

/*--------------------------------------------------------------------------
  mmpr_el_set_id_at_owner - to set an local identifier of an element
---------------------------------------------------------------------------*/
int mmpr_el_set_id_at_owner( /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El,	/* in: global (within the subdomain) element ID */
  int Id_At_Owner /* in: local identifier of the element */
  )
{
/* auxiliary variables */
  mmpt_mesh* mesh;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh = mmpr_select_mesh(Mesh_id);
  
  mesh->elem[El].id_at_owner = Id_At_Owner;

  return(1);
}

/*--------------------------------------------------------------------------
  mmpr_fa_id_at_owner - to get an local identifier of a face
---------------------------------------------------------------------------*/
int mmpr_fa_id_at_owner(  /* returns: >0 - an local identifier of a face */
	          /*          <0 - error code */ 
  int Mesh_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa        /* in: global (within the subdomain)  ID */
  )
{
/* auxiliary variables */
  mmpt_mesh* mesh;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh = mmpr_select_mesh(Mesh_id);
  
  return(mesh->face[Fa].id_at_owner);
}

/*--------------------------------------------------------------------------
  mmpr_fa_set_id_at_owner - to set an local identifier of a face
---------------------------------------------------------------------------*/
int mmpr_fa_set_id_at_owner(  /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Fa,       /* in: global (within the subdomain) face ID */
  int Id_At_Owner      /* in: local identifier of the face */
  )
{
/* auxiliary variables */
  mmpt_mesh* mesh;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh = mmpr_select_mesh(Mesh_id);
  
  mesh->face[Fa].id_at_owner = Id_At_Owner;

  return(1);
}

/*--------------------------------------------------------------------------
  mmpr_ed_id_at_owner - to get an local identifier of an edge 
---------------------------------------------------------------------------*/
int mmpr_ed_id_at_owner(  /* returns: >0 an local identifier of an edge */
	          /*          <0 - error code */ 
  int Mesh_id,    /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ed          /* in: global (within the subdomain) edge ID */
  )
{
/* auxiliary variables */
  mmpt_mesh* mesh;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh = mmpr_select_mesh(Mesh_id);
  
  return(mesh->edge[Ed].id_at_owner);

}

/*--------------------------------------------------------------------------
  mmpr_ed_set_id_at_owner - to set an local identifier of an edge
---------------------------------------------------------------------------*/
int mmpr_ed_set_id_at_owner(  /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ed,       /* in: global (within the subdomain) edge ID */
  int Id_At_Owner      /* in: local identifier of the edge */
  )
{
/* auxiliary variables */
  mmpt_mesh* mesh;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh = mmpr_select_mesh(Mesh_id);
  
  mesh->edge[Ed].id_at_owner = Id_At_Owner;

  return(1);
}

/*--------------------------------------------------------------------------
  mmr_ve_id_at_owner - to get an local identifier of a vertex
---------------------------------------------------------------------------*/
int mmpr_ve_id_at_owner( /* returns: >0 - an local identifier of a vertex */
	         /*          <0 - error code */ 
  int Mesh_id,  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ve        /* in: global (within the subdomain) vertex ID */
  )
{
/* auxiliary variables */
  mmpt_mesh* mesh;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh = mmpr_select_mesh(Mesh_id);
  
  return(mesh->node[Ve].id_at_owner);

}

/*--------------------------------------------------------------------------
  mmpr_ve_set_id_at_owner - to set an local identifier of a vertex
---------------------------------------------------------------------------*/
int mmpr_ve_set_id_at_owner(  /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,  /* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int Ve,       /* in: global (within the subdomain) vertex ID */
  int Id_At_Owner      /* in: local identifier of the vertex */
  )
{
/* auxiliary variables */
  mmpt_mesh* mesh;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh = mmpr_select_mesh(Mesh_id);
  
  mesh->node[Ve].id_at_owner = Id_At_Owner;

  return(1);
}


/*---------------------------------------------------------
  mmpr_create_subdomains - to decompose the mesh and create subdomains
---------------------------------------------------------*/
int mmpr_create_subdomains( /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,      /* in: mesh ID  */
  int Control        /* in: generation level as basis for decomposition */
  )
{


/*++++++++++++++++ executable statements ++++++++++++++++*/
  

  ddr_create_subdomains(Mesh_id, Control);

  return(0);
}


/*------------------------------------------------------------
  mmpr_init_ref - to initialize the process of refinement
------------------------------------------------------------*/
int mmpr_init_ref(  /* returns: >=0 - success code, <0 - error code */
  int Mesh_id	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  )
{

/* local variables */
  mmpt_mesh* par_mesh;
  int iaux;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper mesh data structure */
  par_mesh = mmpr_select_mesh( Mesh_id);

  iaux=mmr_get_nr_elem(Mesh_id);
  par_mesh->ref_loc->elem = (int*) malloc(29*iaux*sizeof(int));
 if(!par_mesh->ref_loc->elem){
   printf("Not enough space for allocating mmpv_ref.... ! Exiting\n");
   exit(-1);
 } 
  par_mesh->ref_loc->elem[0]=0;
  par_mesh->ref_ali->elem = (int*) malloc(29*iaux*sizeof(int));
 if(!par_mesh->ref_ali->elem){
   printf("Not enough space for allocating mmpv_ref.... ! Exiting\n");
   exit(-1);
 } 
  par_mesh->ref_ali->elem[0]=0;

  iaux=mmr_get_nr_face(Mesh_id);
  par_mesh->ref_loc->face = (int*) malloc(11*iaux*sizeof(int));
 if(!par_mesh->ref_loc->face){
   printf("Not enough space for allocating mmpv_ref.... ! Exiting\n");
   exit(-1);
 } 
  par_mesh->ref_loc->face[0]=0;
  par_mesh->ref_ali->face = (int*) malloc(11*iaux*sizeof(int));
 if(!par_mesh->ref_ali->face){
   printf("Not enough space for allocating mmpv_ref.... ! Exiting\n");
   exit(-1);
 } 
  par_mesh->ref_ali->face[0]=0;

  par_mesh->del_loc->face = (int*) malloc(iaux*sizeof(int));
 if(!par_mesh->del_loc->face){
   printf("Not enough space for allocating mmpv_ref.... ! Exiting\n");
   exit(-1);
 } 
  par_mesh->del_loc->face[0]=0;

  iaux=mmr_get_nr_edge(Mesh_id);
  par_mesh->ref_loc->edge = (int*) malloc(4*iaux*sizeof(int));
 if(!par_mesh->ref_loc->edge){
   printf("Not enough space for allocating mmpv_ref.... ! Exiting\n");
   exit(-1);
 } 
  par_mesh->ref_loc->edge[0]=0;
  par_mesh->ref_ali->edge = (int*) malloc(4*iaux*sizeof(int));
 if(!par_mesh->ref_ali->edge){
   printf("Not enough space for allocating mmpv_ref.... ! Exiting\n");
   exit(-1);
 } 
  par_mesh->ref_ali->edge[0]=0;

  par_mesh->del_loc->edge = (int*) malloc(iaux*sizeof(int));
 if(!par_mesh->del_loc->edge){
   printf("Not enough space for allocating mmpv_ref.... ! Exiting\n");
   exit(-1);
 } 
  par_mesh->del_loc->edge[0]=0;

  return(0);
}

/**-----------------------------------------------------------
  mmr_is_ready_for_proj_dof_ref - to check if mesh module is ready
                                  for dofs projection
------------------------------------------------------------*/
int mmpr_is_ready_for_proj_dof_ref(
  int Mesh_id
)
{
    return 1;
}


/*------------------------------------------------------------
  mmpr_final_ref - to finalize the process of parallel refinement 
------------------------------------------------------------*/
int mmpr_final_ref(  /* returns: >=0 - success code, <0 - error code */
  int Mesh_id 	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  )
{
/* local variables */
  mmpt_mesh* par_mesh;
  int i, iel, nmel, num_sons, ison, son;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  par_mesh = mmpr_select_mesh(Mesh_id);
  
  
  /* free lists with the last refinement info */
  if(par_mesh->ref_loc->elem!=NULL){
    free(par_mesh->ref_loc->elem);
  }
  if(par_mesh->ref_loc->face!=NULL){
    free(par_mesh->ref_loc->face);
  }
  if(par_mesh->ref_loc->edge!=NULL){
    free(par_mesh->ref_loc->edge);
  }
  if(par_mesh->ref_ali->elem!=NULL){
    free(par_mesh->ref_ali->elem);
  }
  if(par_mesh->ref_ali->face!=NULL){
    free(par_mesh->ref_ali->face);
  }
  if(par_mesh->ref_ali->edge!=NULL){
    free(par_mesh->ref_ali->edge);
  }
  if(par_mesh->del_loc->face!=NULL){
    free(par_mesh->del_loc->face);
  }
  if(par_mesh->del_loc->edge!=NULL){
    free(par_mesh->del_loc->edge);
  }

  return(0);
}

/*------------------------------------------------------------
  mmpr_free_mesh - to free space allocated for parallel mesh data structure
------------------------------------------------------------*/
int mmpr_free_mesh(  /* returns: >=0 - success code, <0 - error code */
  int Mesh_id 	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  )
{

/* local variables */
  mmpt_mesh* par_mesh;
  int iaux;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper mesh data structure */
  par_mesh = mmpr_select_mesh( Mesh_id);

  if(par_mesh->node != NULL) free(par_mesh->node);
  if(par_mesh->edge != NULL) free(par_mesh->edge);
  if(par_mesh->face != NULL) free(par_mesh->face);
  if(par_mesh->elem != NULL) free(par_mesh->elem);
  //if(par_mesh->ref_loc != NULL) free(par_mesh->ref_loc);
  //if(par_mesh->ref_ali != NULL) free(par_mesh->ref_ali);
  //if(par_mesh->del_loc != NULL) free(par_mesh->del_loc);
  //if(par_mesh->mesh_ovls != NULL) free(par_mesh->mesh_ovls);

  return(0);
}


/*---------------------------------------------------------
  mmpr_update_ref_list (old ddr_exchange_list_ref) - to update the list 
                       of refined elements using inter-processor communication
---------------------------------------------------------*/
int mmpr_update_ref_list( /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,      /* in: mesh ID */
  int *Nr_ref_p,    /* in/out: number of refined elements */
  int **List_ref_p     /* in/out: list of refined elements */
  //int *List_ref     /* in/out: list of refined elements */ - old interface
)
{

  // OLD VERSION - NOT NECESSARILY WORKING

  /* ddt_mesh_ovl* mesh_ovl_p; /\* lists of overlap entities *\/ */

  /* int i, iaux, isub, iproc, buffer_id; */
  /* int nr_ref, iel, ient, nel; */

/*++++++++++++++++ executable statements ++++++++++++++++*/


  
/*   nr_ref = *Nr_ref_p; */
/*   List_ref = *List_ref_p; */

/*   /\* proceed only if necessary *\/ */
/*   if(ddv_nr_proc==1) return(0); */

/*   /\* select overlap for the proper mesh *\/ */
/*   mesh_ovl_p = &ddv_mesh_ovls[Mesh_id-1]; */

/*   /\* check whether refined entities are shared with other processors *\/ */
/*   /\* if yes - inform other processors on ipids of refined entities *\/ */
/*   /\* pack info into buffers and send *\/ */

/*   //  buffer_id = PCC_UPD_REF_1; */
/*   buffer_id = 7821; */
/*   for(isub=0;isub<ddv_nr_proc;isub++){ */
/* /\*kbw */
/*       printf("processor %d, nr owned elems %d\n", */
/* 	     isub+1,mesh_ovl_p->nr_elem_ovl[isub]); */
/* /\*kew*\/ */

/*     if(mesh_ovl_p->nr_elem_ovl[isub]>0){ */

/*       pcr_send_buffer_open(buffer_id); */
/*       pcr_buffer_pack_int(buffer_id, 1, &ddv_my_proc_id); */

/*       for(ient=0;ient<nr_ref;ient++){ */
/* 	nel=List_ref[ient]; */
/* 	for(iel=0;iel<mesh_ovl_p->nr_elem_ovl[isub];iel++){ */
/* 	  if(nel==mesh_ovl_p->l_elem_ovl_loc[isub][iel]){ */
/* 	    pcr_buffer_pack_int(buffer_id, 1,  */
/* 				&mesh_ovl_p->l_elem_ovl_ali[isub][iel]); */
/* /\*kbw */
/*       printf("SENT element %d (%d) to processor %d\n", */
/* 	     nel,mesh_ovl_p->l_elem_ovl_ali[isub][iel],isub+1); */
/* /\*kew*\/ */
/* 	    break; */
/* 	  } /\* end if shared element has to be refined *\/ */
/* 	} /\* end loop over shared (overlap) elements *\/ */
/*       } /\* end loop over refined elements *\/ */
/*       /\* pack end marker *\/ */
/*       iaux = -1; */
/*       pcr_buffer_pack_int(buffer_id, 1, &iaux); */
/*       pcr_buffer_send(buffer_id, isub+1); */

/* /\*kbw */
/*       printf("SENT elements to processor %d\n",isub+1); */
/* /\*kew*\/ */

/*     } /\* end if there are overlap elements *\/ */
/*   } /\* end loop over subdomains-processors *\/ */

/*   /\*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*\/ */

/*   //  buffer_id = PCC_UPD_REF_1; */
/*   buffer_id = 7821; */

/*   for(iproc=0;iproc<ddv_nr_proc;iproc++){ */

/*     if(mesh_ovl_p->nr_elem_alien[iproc]>0){ */

/*       pcr_buffer_receive_any(buffer_id); */

/*       /\* sender ID *\/ */
/*       pcr_buffer_unpack_int(buffer_id, 1, &i); */
/*       isub=i-1; */

/* /\*kbw */
/*       printf("DATA from processor %d\n",i); */
/* /\*kew*\/ */

/*       /\* refined element ID *\/ */
/*       pcr_buffer_unpack_int(buffer_id, 1, &ient); */
/*       while(ient>0){ */

/* /\*kbw */
/* 	printf("FROM processor %d, refined element %d\n",isub+1,ient); */
/* /\*kew*\/ */


/* 	/\* check element is not on the original list *\/   */
/* 	i=ddr_chk_list(ient, List_ref, nr_ref); */
/* 	if(i==0){ */
/* 	  List_ref[nr_ref]=ient; */
/* 	  nr_ref++; */
/* 	} */

/* 	/\* unpack next element ID or end mark *\/ */
/* 	pcr_buffer_unpack_int(buffer_id, 1, &ient); */

/*       } */
	  
/*       pcr_recv_buffer_close(buffer_id); */

/*     } /\* end if subdomain with shared entities *\/ */
    
/*   } /\* end loop over all subdoamins *\/ */
	
/*   *Nr_ref_p = nr_ref; */

  return(1);
}

