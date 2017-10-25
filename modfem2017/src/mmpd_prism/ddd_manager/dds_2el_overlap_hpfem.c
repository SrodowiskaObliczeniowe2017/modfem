/************************************************************************
File dds_2el_overlap_hpfem.c - overlap management routines, support exchange 
                     of data between processors 

Contains definitions of routines:   
  ddr_node_elems - to prepare lists of all elements (active and
                   inactive) to which belong subsequent nodes
  ddr_create_overlap

------------------------------  			
History:        
	02.2002 - Krzysztof Banas, initial version		
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

/* interface of the mesh manipulation module */
#include "mmh_intf.h"	

/* interface of the mesh manipulation module */
#include "mmph_intf.h"	

/* interface of the mesh manipulation module */
#include "../mmph_prism.h"	

/* domain decomposition interface specification */
#include "./ddh_front.h"

/* internal info for domain decomposition manager */
#include "./ddh_manager.h"

#define MAXELNO 100
static int **elems=NULL;

static int max_node=0;

/*---------------------------------------------------------
  ddr_vert_elems - to prepare lists of all elements (active and
                   inactive) to which belong subsequent nodes
---------------------------------------------------------*/
int ddr_vert_elems(
  int Mesh_id,      /* in: mesh ID */
  int Gen_lev,       /* in: generation level as basis for decomposition */
  int **Vert_elems  /* out: list of lists of vertices' elements */
)
{

  const int zero =  0;
  int max_elem_id, max_node_id;
  int i, iaux, nel, nno, ino, ipr;
  int el_nodes[MMC_MAXELVNO+1];

/*++++++++++++++++ executable statements ++++++++++++++++*/
  
  max_node_id = mmr_get_max_node_id(Mesh_id);

  if(Vert_elems==NULL){
    printf("Old entry to vert-elems with NULL array - change\n");
    exit(-1);
  } /* end if Vert_elems==NULL */
  else{

    /* divide nodes into external and other */
    nel=0;
    while((nel=mmr_get_next_elem_all(Mesh_id, nel))!=0){

      //      if(mmr_el_gen(Mesh_id,nel)<=Gen_lev &&
      //	 mmpr_el_owner(Mesh_id,nel)==mmpv_my_proc_id){
      if(mmr_el_gen(Mesh_id,nel)==Gen_lev &&
	 mmpr_el_owner(Mesh_id,nel)==mmpv_my_proc_id){

	/* for all owned element's nodes */ 
	mmr_el_node_coor(Mesh_id, nel, el_nodes, NULL);

	for(ino=0;ino<el_nodes[0];ino++){

	  nno = el_nodes[ino+1];

/*kbw
	  printf("proc %d, nel %d, nno %d\n",
		 mmpv_my_proc_id, nel, nno);
/*kew*/

	  /* allocate Vert_elems for internal or boundary node */
	  if(Vert_elems[nno]==NULL){

	    Vert_elems[nno] = malloc(MAXELNO*sizeof(int));
	    Vert_elems[nno][0]=0;
#ifdef DEBUG_MMM
	    for(i=1;i<MAXELNO;i++) Vert_elems[nno][i]=0;
#endif
	  }

	} /* end for all element's nodes: ino */
      } /* end if owned element of proper generation*/
    } /* end loop over all elements: nel */


    /* loop over elements */
    nel=0;
    while((nel=mmr_get_next_elem_all(Mesh_id, nel))!=0){

      if(mmr_el_gen(Mesh_id,nel)<=Gen_lev){

	mmr_el_node_coor(Mesh_id,nel,el_nodes,NULL);
	
	for(nno=1;nno<=el_nodes[0];nno++){
	  ino=el_nodes[nno];
	  
	  if(Vert_elems[ino]!=NULL){
	    
#ifdef DEBUG_MMM
	    i=ddr_chk_list(nel,&Vert_elems[ino][1],MAXELNO-1);
	    if(i>0){
	      printf("Something wrong in filling NODES->Vert_elems for node %d\n",ino);
	      for(ipr=1;ipr<MAXELNO;ipr++) 
		printf("%d\n",Vert_elems[ino][ipr]);
	      printf("\n");
	      exit(-1);;
	    }
#endif
	    
	    Vert_elems[ino][0]++;
	    Vert_elems[ino][Vert_elems[ino][0]] = nel;

/*kbw
	printf("node %d (%d), added element %d, gen %d, list (length %d): ",
	       ino, el_nodes[nno], nel, mmr_el_gen(Mesh_id, nel), Vert_elems[ino][0]);
	for(ipr=1;ipr<=Vert_elems[ino][0];ipr++){
	  printf("%d ",Vert_elems[ino][ipr]);
	}
	printf("\n");
/*kew*/
	  } /* end if node internal or boundary */
	} /* end if element of proper generation */
      } /* end for all element's nodes: ino */
    } /* end loop over all elements: nel */
  }

#ifdef DEBUG_MMM
  iaux=0;
  for(ino=1;ino<=max_node_id;ino++){
    if(Vert_elems[ino]!=NULL&&Vert_elems[ino][0]>iaux) iaux = Vert_elems[ino][0];
  }
  printf("Maximal number of vert_elems - %d\n", iaux);
#endif

  return(0);
}

/*---------------------------------------------------------
  ddr_add_elem_ovl - to add an element to the overlap with all
                     its antecedents in a recursive manner
---------------------------------------------------------*/
int ddr_add_elem_ovl(
  int Mesh_id, /* in: mesh ID */
  int El       /* in: element ID */
  );

/*---------------------------------------------------------
  ddr_add_face_ovl - to add a face to the overlap with all
                     its antecedents in a recursive manner
---------------------------------------------------------*/
int ddr_add_face_ovl(
  int Mesh_id, /* in: mesh ID */
  int Fa       /* in: face ID */
  );

/*---------------------------------------------------------
  ddr_add_edge_ovl - to add an edge to the overlap with all
                     its antecedents in a recursive manner
---------------------------------------------------------*/
int ddr_add_edge_ovl(
  int Mesh_id, /* in: mesh ID */
  int Ed       /* in: edge ID */
  );

/*---------------------------------------------------------
  ddr_create_overlap
---------------------------------------------------------*/
int ddr_create_overlap(  /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,      /* in: mesh ID */
  int Ovl_size,
  int Gen_lev       /* in: generation level as basis for decomposition */
  )
{

  int nno, nel, ino, iel, ied, ifa, nfa, owner, iaux, jaux, i, ient;
  int max_node_id, nr_alien;
  int el_nodes[MMC_MAXELVNO+1];
  const int zero =  0;
  const int one  =  1;
  const int mone = -1;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  max_node_id = mmr_get_max_node_id(Mesh_id);

  /* if elems arrays do not exist */
  if(elems==NULL || max_node_id!=max_node){

    /* if number of nodes changed */ 
    if(max_node_id!=max_node){

      for(ino=1;ino<=max_node;ino++) {
	if(elems[ino]!=NULL) free(elems[ino]);
      }
      free(elems);

    }
    max_node = max_node_id;

    elems = malloc((max_node_id+1)*sizeof(int *));
    for(ino=1;ino<=max_node_id;ino++){
      elems[ino]=NULL;  
    }
    /* prepare array of element connectivities - all vertices belonging to */
    /* owned elements get elems array with lists of all adjacent elements */
    ddr_vert_elems(Mesh_id, Gen_lev, elems); 

  }

// reset temporarily ownership
  i=0; iaux = 0;
  while((i=mmr_get_next_elem_all(Mesh_id, i))!=0){
    if(mmpr_el_owner(Mesh_id,i)==mmpv_my_proc_id)
      mmpr_el_set_owner(Mesh_id,i,iaux);
  }
  i=0;
  while((i=mmr_get_next_face_all(Mesh_id, i))!=0){
    if(mmpr_fa_owner(Mesh_id,i)==mmpv_my_proc_id)
      mmpr_fa_set_owner(Mesh_id,i,iaux);
  }
  i=0;
  while((i=mmr_get_next_edge_all(Mesh_id, i))!=0){
    if(mmpr_ed_owner(Mesh_id,i)==mmpv_my_proc_id)
      mmpr_ed_set_owner(Mesh_id,i,iaux);
  }
  i=0;
  while((i=mmr_get_next_node_all(Mesh_id, i))!=0){
    if(mmpr_ve_owner(Mesh_id,i)==mmpv_my_proc_id)
      mmpr_ve_set_owner(Mesh_id,i,iaux);
  }


#ifdef DEBUG_MMM
  i=0; 
  while((i=mmr_get_next_elem_all(Mesh_id, i))!=0){
    if(mmpr_el_owner(Mesh_id,i)<0){
      printf("Error before creating overlap: element %d, owner %d ???\n",
	   i, mmpr_el_owner(Mesh_id,i));
      exit(-1);
    }
  }
#endif

  /* create two-element overlap (each subdomain extended by one element) */
  /* - suitable for DG-FEM approximation */
  ino=0;
  while((ino=mmr_get_next_node_all(Mesh_id, ino))!=0){

#ifdef DEBUG_MMM
    if(mmpr_ve_owner(Mesh_id, ino)>0 &&
       mmpr_ve_owner(Mesh_id, ino)==mmpv_my_proc_id){
      printf("Node %d, owned by %d - not zeroed in %d!!!\n", ino, 
	     mmpr_ve_owner(Mesh_id, ino), mmpv_my_proc_id);
      exit(-1);
    }
#endif

    /* for each old boundary and internal node */
    if(elems[ino]!=NULL){

      nr_alien=0;
      for(iel=1;iel<=elems[ino][0];iel++){
	if( mmpr_el_owner(Mesh_id,elems[ino][iel]) != mmpv_my_proc_id 
	    && mmpr_el_owner(Mesh_id,elems[ino][iel]) != 0 ){
	  nr_alien++;
	}
      }

#ifdef DEBUG_MMM
      if(nr_alien>=elems[ino][0]){
	printf("Node %d, internal to %d - in overlap of %d!!!\n", ino, 
	       mmpr_ve_owner(Mesh_id, ino), mmpv_my_proc_id);
	exit(-1);
      }
#endif

/*kbw
  if(nel==323){
    printf("adding element %d (owner %d, nr_alien %d) to overlap, father %d\n",
	   nel, mmpr_el_owner(Mesh_id, nel), nr_alien,
	   mmr_el_fam(Mesh_id, nel, NULL, NULL));
  }
/*kew*/

      if(nr_alien>0){

	/* old boundary node - will become internal */
/*kbw
	printf("Node %d, owned by %d - old boundary of %d\n", ino, 
	       mmpr_ve_owner(Mesh_id, ino), mmpv_my_proc_id);
/*kew*/

	/* for all elements adjacent to old boundary node */
	for(iel=1;iel<=elems[ino][0];iel++){ 

	  nel=elems[ino][iel];


#ifdef DEBUG_MMM
	  if(mmpr_el_owner(Mesh_id,nel)==mmpv_my_proc_id){
	    printf("Error : internal element %d, subdoamin %d - not zeroed\n",
		   nel, mmpv_my_proc_id);
	    exit(-1);
	  }
#endif

	  /* to create two-element overlap between subdomains */
	  if(mmr_el_gen(Mesh_id,nel)==Gen_lev && 
	     mmpr_el_owner(Mesh_id,nel)>0){
//for one element overlap ? mmpv_my_proc_id>mmpr_el_owner(Mesh_id,nel)){

	    //printf("Warning: 2 INITIAL elements overlap !!!\n");

/*kbw
	    printf("Overlap element %d, owned by %d\n",
		   nel,mmpr_el_owner(Mesh_id,nel));
/*kew*/

	    int el_fath, f_neig[MMC_MAXELFAC+1];
	  
	    /* find father and check whether its neighbors are active */
	    el_fath=mmr_el_fam(Mesh_id, nel, NULL, NULL);
	    
	    if(el_fath>0){
#ifdef DEBUG_MMM
	      printf("Overlap element %d, generation %d > 0\n",
		     nel,mmr_el_gen(Mesh_id,nel));
	      exit(-1);
#endif
	      mmr_el_eq_neig(Mesh_id, el_fath, f_neig, NULL);
	      for(i=1;i<=f_neig[0];i++){
#ifdef DEBUG_MMM
		if(f_neig[i]==-1){
		  printf("broken 1-irregularity - cannot create overlap\n");
		  printf("element %d, father %d, fath_neig%d %d\n",
			 nel, el_fath, i, f_neig[i]);
		  exit(-1);
		}
#endif
	      
		if(f_neig[i]>0 && mmr_el_status(Mesh_id,f_neig[i])==MMC_ACTIVE){
		  /* neighbor active, include father to overlap */
		  nel=el_fath;
		  continue;
		}
	      }
	    }
	    
	    ddr_add_elem_ovl(Mesh_id, nel);
	    
	  } /* end if element not owned and of proper generation */
	} /* end for all node's elements */
    
      } /* end if old boundary node */ 
    } /* end if old boundary or internal node */
  } /* end loop over all nodes */


  /* check element ownership */
  iel=0;
  while((iel=mmr_get_next_elem_all(Mesh_id, iel))!=0){

#ifdef DEBUG_MMM
    if(mmpr_el_owner(Mesh_id,iel)==mmpv_my_proc_id){
      printf("Error : internal element %d, subdoamin %d - not zeroed\n",
	     iel, mmpv_my_proc_id);
      exit(-1);
    }
#endif

    owner=mmpr_el_owner(Mesh_id,iel);
    if(owner<0) {
      mmpr_el_set_owner(Mesh_id,iel,-owner);
/*kbw
      printf("element %d, in overlap of subdoamin %d (owner %d)\n",
	     iel, mmpv_my_proc_id, (-owner));
/*kew*/
    }
    else if(owner>0) {
#ifdef DEBUG_MMM
      if(mmpr_el_owner(Mesh_id,iel)==mmpv_my_proc_id){
	printf("element %d, external to %d owned by %d\n",
	       iel, mmpv_my_proc_id, (owner));
	exit(-1);
      }
#endif
      mmpr_el_set_owner(Mesh_id,iel,mone);
/*kbw
      printf("element %d, external to %d owned by %d\n",
	     iel, mmpv_my_proc_id, (owner));
/*kew*/
    }
/*kbw
    else{
      printf("element %d, internal to %d owned by %d\n",
	     iel, mmpv_my_proc_id, (owner));
	exit(-1);
    }
/*kew*/
  }

  /* check again for internal and boundary nodes */ 
  ino=0;
  while((ino=mmr_get_next_node_all(Mesh_id, ino))!=0){

#ifdef DEBUG_MMM
    if(abs(mmpr_ve_owner(Mesh_id,ino))==mmpv_my_proc_id){
      printf("Error : owned node %d, subdoamin %d - not zeroed\n",
	     ino, mmpv_my_proc_id);
      mmpr_ve_set_owner(Mesh_id, ino, zero);
      exit(-1);
    }
#endif

    owner=mmpr_ve_owner(Mesh_id,ino);

    if(owner==0){
/*kbw
      printf("Node %d, owned by %d - internal to %d\n", ino, 
	     mmpr_ve_owner(Mesh_id,  ino), mmpv_my_proc_id);
/*kew*/
    }
    else if(owner<0){
/*kbw
      printf("Node %d, owned by %d - close to boundary of %d\n", ino, 
	     mmpr_ve_owner(Mesh_id,  ino), mmpv_my_proc_id);
/*kew*/
      mmpr_ve_set_owner(Mesh_id, ino, -owner);
    }
    else{
/*kbw
      printf("Node %d, owned by %d - external to %d\n", ino, 
	     mmpr_ve_owner(Mesh_id,  ino), mmpv_my_proc_id);
      getchar();
/*kew*/
      mmpr_ve_set_owner(Mesh_id, ino, mone);
    }
  }

  /* check */

  /* indicate ownership for edges */
  ied=0;
  while((ied=mmr_get_next_edge_all(Mesh_id, ied))!=0){

    iaux=mmpr_ed_owner(Mesh_id, ied);
    if(iaux<0){
      mmpr_ed_set_owner(Mesh_id, ied, -iaux);
    }
    else if(iaux>0){
      mmpr_ed_set_owner(Mesh_id, ied, mone);
/*kbw
      if(ied==15443||ied==15444){
      printf("Edge %d, owned by %d - external to %d\n", ied, 
	     mmpr_ed_owner(Mesh_id,  ied), mmpv_my_proc_id);
      getchar();
      }
/*kew*/
    }

/*kbw
    if(ied==15443||ied==15444){
      printf("Edge %d, owner %d\n",ied,iaux);
      getchar();
    }
/*kew*/

  }

  /* indicate ownership for faces */
  ifa=0;
  while((ifa=mmr_get_next_face_all(Mesh_id, ifa))!=0){

    iaux=mmpr_fa_owner(Mesh_id, ifa);

/*kbw
	if(ifa>=1976&&ifa<=1979){
	  printf("face %d, iaux %d\n",ifa,iaux);
	}
/*kew*/

    if(iaux<0){
      mmpr_fa_set_owner(Mesh_id, ifa, -iaux);
    }
    else if(iaux>0){
      mmpr_fa_set_owner(Mesh_id, ifa, mone);
/*kbw
      printf("Face %d, owned by %d - external to %d\n", ifa, 
	     mmpr_fa_owner(Mesh_id,  ifa), mmpv_my_proc_id);
      getchar();
/*kew*/
    }

  }


  max_node_id = mmr_get_max_node_id(Mesh_id);
  for(ino=1;ino<=max_node_id;ino++){

    if(elems[ino]!=NULL){
      free(elems[ino]);
    }

  }
  free(elems);
  elems=NULL;


  /* indicate boundary faces */
  iel=0;
  while((iel=mmr_get_next_elem_all(Mesh_id, iel))!=0){

    if(mmpr_el_owner(Mesh_id,iel)>=0){

      int el_faces[MMC_MAXELFAC+1];

      mmr_el_faces(Mesh_id, iel, el_faces, NULL);

      for(nfa=0;nfa<el_faces[0];nfa++) {

	int fa_neig[2], fa_sons[5], neig, ifa, ineig, side;

	ifa = el_faces[nfa+1];

/*kbw
	if(ifa==10325||ifa==96){
	  //if(ifa>=1976&&ifa<=1979){
	  printf("element %d, face %d\n",iel, ifa);
	}
/*kew*/

	mmr_fa_eq_neig(Mesh_id, ifa, fa_neig, NULL, NULL);
/*kbw
	if(ifa==10325||ifa==96){
	  //if(ifa>=1976&&ifa<=1979){
	  printf("equal size neighbors %d, %d\n",fa_neig[0], fa_neig[1]);
	}
/*kew*/
	for(ineig=0;ineig<2;ineig++){
	  if(fa_neig[ineig]!=iel){
	    side=ifa;
	    neig=fa_neig[ineig];
	    while(neig==MMC_BIG_NGB){
	      side=mmr_fa_fam(Mesh_id, side, NULL, NULL);
	      if(side==0){
		neig=MMC_SUB_BND;
	      }
	      else{
		mmr_fa_eq_neig(Mesh_id, side, fa_neig, NULL, NULL);
		neig=fa_neig[ineig];
	      }
	    }
	    
	    if(neig==MMC_SUB_BND){
/*kbw
		if(ifa==10325||ifa==96){
		  printf("sub bnd %d, side %d\n",neig, ineig);
		}
/*kew*/
	      mmr_fa_set_sub_bnd(Mesh_id,ifa,ineig);
	      if(mmr_fa_status(Mesh_id,ifa)==MMC_INACTIVE){
		/* add children */
		mmr_fa_fam(Mesh_id, ifa, fa_sons, NULL);
		for(i=1;i<=fa_sons[0];i++){
/*kbw
		if(fa_sons[i]==10325){
		  printf("face %d, as child of %d, sub bnd %d, side %d\n",
			 fa_sons[i], ifa, neig, ineig);
		}
/*kew*/
		  mmr_fa_set_sub_bnd(Mesh_id,fa_sons[i],ineig);
		}
	      }

	      continue;
	    }
	    else if(neig>0){
	      if(mmpr_el_owner(Mesh_id,neig)<0){
/*kbw
		if(ifa==10325||ifa==96){
		  printf("external neighbor %d, side %d\n",neig, ineig);
		}
/*kew*/
		mmr_fa_set_sub_bnd(Mesh_id,ifa,ineig);
		if(mmr_fa_status(Mesh_id,ifa)==MMC_INACTIVE){
		  /* add children */
		  mmr_fa_fam(Mesh_id, ifa, fa_sons, NULL);
		  for(i=1;i<=fa_sons[0];i++){
/*kbw
		    if(fa_sons[i]==10325){
		      printf("face %d, as child of %d, sub bnd %d, side %d\n",
			     fa_sons[i], ifa, neig, ineig);
		    }
/*kew*/
		    mmr_fa_set_sub_bnd(Mesh_id,fa_sons[i],ineig);
		  }
		}
		continue;
	      }
	    }
	  }
	}
      }
    }
  }
  

  /* delete external entities */
  /* since we delete structures we cannot use standard iterators */

  /* elements */
  iaux = mmr_get_max_elem_id(Mesh_id);
  for(ient=1;ient<=iaux;ient++){
    if(mmr_el_status(Mesh_id,ient)!=MMC_FREE&&mmpr_el_owner(Mesh_id,ient)<0){
/*kbw
      printf("Subdomain %d, deleting element %d (owner %d)\n", 
	     mmpv_my_proc_id, ient, mmpr_el_owner(Mesh_id,ient));
/*kew*/
      mmr_del_elem(Mesh_id,ient);
    }
  }

  /* faces */  
  iaux = mmr_get_max_face_id(Mesh_id);
  for(ient=1;ient<=iaux;ient++){
    if(mmr_fa_status(Mesh_id,ient)!=MMC_FREE&&mmpr_fa_owner(Mesh_id,ient)<0){
/*kbw
      printf("Subdomain %d, deleting face %d (owner %d)\n", 
	     mmpv_my_proc_id, ient, mmpr_fa_owner(Mesh_id,ient));
/*kew*/
      mmr_del_face(Mesh_id,ient);
    }
  }

  /* edges */ 
  iaux = mmr_get_max_edge_id(Mesh_id);
  for(ient=1;ient<=iaux;ient++){
    if(mmr_edge_status(Mesh_id,ient)!=MMC_FREE&&mmpr_ed_owner(Mesh_id,ient)<0){
/*kbw
      printf("Subdomain %d, deleting edge %d (owner %d)\n", 
	     mmpv_my_proc_id, ient, mmpr_ed_owner(Mesh_id,ient));
/*kew*/
/*kbw
      if(ient==15443||ient==15444){
      printf("Edge %d, owned by %d - external to %d - DELETING\n", ient, 
	     mmpr_ed_owner(Mesh_id,  ient), mmpv_my_proc_id);
      getchar();
      }
/*kew*/
      mmr_del_edge(Mesh_id,ient);
    }
  }

  /* vertices */ 
  iaux = mmr_get_max_node_id(Mesh_id);
  for(ient=1;ient<=iaux;ient++){
    if(mmr_node_status(Mesh_id,ient)!=MMC_FREE&&mmpr_ve_owner(Mesh_id,ient)<0){
/*kbw
      printf("Subdomain %d, deleting vertex %d (owner %d)\n", 
	     mmpv_my_proc_id, ient, mmpr_ve_owner(Mesh_id,ient));
/*kew*/
      mmr_del_node(Mesh_id,ient);
    }
  }



  /* check consistency */ 
#ifdef DEBUG_MMM
  iel=0;
  while((iel=mmr_get_next_elem_all(Mesh_id, iel))!=0){

    int el_faces[MMC_MAXELFAC+1];
    int j, ifa, num_edges, fa_edges[5];

    if(mmpr_el_owner(Mesh_id,iel)<0) {
      printf("element %d, external to %d present in data structure\n",
	     iel, mmpv_my_proc_id);
      exit(-1);
    }

    mmr_el_faces(Mesh_id, iel, el_faces, NULL);

    for(i=0;i<el_faces[0];i++) {

      ifa = el_faces[i+1];

      if(mmr_fa_status(Mesh_id, ifa)==MMC_FREE || mmpr_fa_owner(Mesh_id, ifa)<0){

	printf("element %d, face %d, external to %d present in data structure\n",
	       iel, ifa, mmpv_my_proc_id);
	exit(-1);
      }

      num_edges = mmr_fa_edges(Mesh_id, ifa, fa_edges, NULL);

      for(j=0;j<num_edges;j++){

	ied = fa_edges[j];

	if(mmr_edge_status(Mesh_id,ied)==MMC_FREE||mmpr_ed_owner(Mesh_id,ied)<0){
	  
	  printf("element %d (father %d), face %d, edge %d, external to %d, status %d, owner %d - error\n",
		 iel,mmr_el_fam(Mesh_id,iel,NULL,NULL),ifa,ied,mmpv_my_proc_id, 
		 mmr_edge_status(Mesh_id,ied),
		 mmpr_ed_owner(Mesh_id,ied) );
	  exit(-1);
	}

      } 

    } 

    mmr_el_node_coor(Mesh_id, iel, el_nodes, NULL);

    for(i=0;i<el_nodes[0];i++){

      ino = el_nodes[i+1];

      if(mmr_node_status(Mesh_id,ino)==MMC_FREE||mmpr_ve_owner(Mesh_id, ino)<0){
	printf("node %d, external to %d present in data structure\n",
	       ino, mmpv_my_proc_id);
	exit(-1);
      }
	
    }

  } 
#endif
/*kew*/

// restore back proper ownership info
  i=0; iaux = mmpv_my_proc_id;
  while((i=mmr_get_next_elem_all(Mesh_id, i))!=0){
    if(mmpr_el_owner(Mesh_id,i)==0)
      mmpr_el_set_owner(Mesh_id,i,iaux);
  }
  i=0;
  while((i=mmr_get_next_face_all(Mesh_id, i))!=0){
    if(mmpr_fa_owner(Mesh_id,i)==0)
      mmpr_fa_set_owner(Mesh_id,i,iaux);
  }
  i=0;
  while((i=mmr_get_next_edge_all(Mesh_id, i))!=0){
    if(mmpr_ed_owner(Mesh_id,i)==0)
/*kbw
      if(i==15443||i==15444){
      printf("Edge %d, owned by %d - external to %d - SETTING OWNERSHIP\n", i, 
	     mmpr_ed_owner(Mesh_id,  i), mmpv_my_proc_id);
      getchar();
      }
/*kew*/
      mmpr_ed_set_owner(Mesh_id,i,iaux);
  }
  i=0;
  while((i=mmr_get_next_node_all(Mesh_id, i))!=0){
    if(mmpr_ve_owner(Mesh_id,i)==0)
      mmpr_ve_set_owner(Mesh_id,i,iaux);
  }



  return(0);
}


/*---------------------------------------------------------
  ddr_add_elem_ovl - to add an element to the overlap with all
                     its antecedents in a recursive manner
---------------------------------------------------------*/
int ddr_add_elem_ovl(
  int Mesh_id, /* in: mesh ID */
  int El       /* in: element ID */
  )
{

  int i,j,k, ient, ifa, jfa, kfa, nfa, ison, owner;
  int iaux, nr_sub, isub, iel, nel, ned, ied, ino, nr_elem, gen_lev, mesh_id;
  int elsons[MMC_MAXELSONS+1];
  int el_faces[MMC_MAXELFAC+1];

/*++++++++++++++++ executable statements ++++++++++++++++*/

/*kbw
  if(El==4105||El==323){
    printf("adding element %d to overlap, father %d\n",El, 
	   mmr_el_fam(Mesh_id, El, NULL, NULL));
  }
 /*kew*/

  /* if element external - indicate inclusion into overlap */
  owner=mmpr_el_owner(Mesh_id,El);
  if(owner>0) {
    mmpr_el_set_owner(Mesh_id,El,-owner);
  }	
  else{
    return(0);
  }

  /* add faces to overlap */
  mmr_el_faces(mesh_id, El, el_faces, NULL);
  
  for(nfa=0;nfa<el_faces[0];nfa++) {
	
    ifa = el_faces[nfa+1];

    owner=mmpr_fa_owner(Mesh_id,ifa);

/*kbw
	if((ifa>=1976&&ifa<=1979)||ifa==416){
	  printf("elem %d, face %d, owner %d\n",El,ifa,owner);
	}
/*kew*/


    if(owner>0) {
      ddr_add_face_ovl(Mesh_id,ifa);
    }

  }

  if(mmr_el_status(mesh_id, El)==MMC_INACTIVE){

    /* add sons to overlap */
    mmr_el_fam(Mesh_id,El,elsons,NULL);
    
    for(ison=1;ison<=elsons[0];ison++){
      ddr_add_elem_ovl(Mesh_id,elsons[ison]);
    }

  }

  return(1);
}


/*---------------------------------------------------------
  ddr_add_face_ovl - to add a face to the overlap with all
                     its antecedents in a recursive manner
---------------------------------------------------------*/
int ddr_add_face_ovl(
  int Mesh_id, /* in: mesh ID */
  int Fa       /* in: face ID */
  )
{

  int i,j,k, jfa, kfa, ison, owner;
  int ned, ied;
  int num_edges, fa_sons[5], fa_edges[5];

/*++++++++++++++++ executable statements ++++++++++++++++*/

  //printf("\tadding face %d to overlap\n",Fa);

  owner=mmpr_fa_owner(Mesh_id,Fa);


/*kbw
	if(Fa>=1976&&Fa<=1979){
	  printf("adding to overlap face %d, owner %d\n",Fa,owner);
	}
/*kew*/

  if(owner>0) {
    mmpr_fa_set_owner(Mesh_id,Fa,-owner);
  }	
  else{
    return(0);
  }

  /* add edges to overlap */
  num_edges = mmr_fa_edges(Mesh_id, Fa, fa_edges, NULL);
	  
  for(ned=0;ned<num_edges;ned++){
	    
    ied = fa_edges[ned];
	
    owner=mmpr_ed_owner(Mesh_id,ied);
    if(owner>0) {
/*kbw
      if(ied==15443||ied==15444){
      printf("Edge %d, owned by %d - external to %d - ADDING TO OVERLAP\n", ied, 
	     mmpr_ed_owner(Mesh_id,  ied), mmpv_my_proc_id);
      getchar();
      }
/*kew*/
      ddr_add_edge_ovl(Mesh_id,ied);
      }

  }

  /* add children to overlap */
  if(mmr_fa_status(Mesh_id, Fa)==MMC_INACTIVE){

    mmr_fa_fam(Mesh_id, Fa, fa_sons, NULL);

    for(jfa=1; jfa<=fa_sons[0]; jfa++){

      kfa=fa_sons[jfa];

/*kbw
	if(kfa>=1976&&kfa<=1979){
	  printf("father face %d, owner %d\n",Fa,owner);
	}
/*kew*/

      owner = mmpr_fa_owner(Mesh_id,kfa);

      if(owner>0){
	ddr_add_face_ovl(Mesh_id,kfa);
      }
      
    }

  }

  return(1);
}


/*---------------------------------------------------------
  ddr_add_edge_ovl - to add an edge to the overlap with all
                     its antecedents in a recursive manner
---------------------------------------------------------*/
int ddr_add_edge_ovl(
  int Mesh_id, /* in: mesh ID */
  int Ed       /* in: edge ID */
  )
{

  int i,j,k, ison, owner;
  int num_edges, edge_sons[2], edge_nodes[2];

/*++++++++++++++++ executable statements ++++++++++++++++*/

  //printf("\t\tadding edge %d to overlap\n",Ed);
/*kbw
  if(Ed==15443||Ed==15444){
      printf("Edge %d, owned by %d - external to %d - adding to overlap\n", Ed, 
	     mmpr_ed_owner(Mesh_id,  Ed), mmpv_my_proc_id);
      getchar();
      }
/*kew*/

  owner=mmpr_ed_owner(Mesh_id,Ed);
  if(owner>0) {
    mmpr_ed_set_owner(Mesh_id,Ed,-owner);
  }	
  else{
    return(0);
  }

  /* add nodes to overlap */	      
  mmr_edge_nodes(Mesh_id, Ed, edge_nodes);
	      
  if((i=mmpr_ve_owner(Mesh_id, edge_nodes[0])) > 0){
    mmpr_ve_set_owner(Mesh_id, edge_nodes[0], -i);
  }
  if((i=mmpr_ve_owner(Mesh_id, edge_nodes[1])) > 0){
    mmpr_ve_set_owner(Mesh_id, edge_nodes[1], -i);
  }
  
  /* add children to overlap */
  if(mmr_edge_status(Mesh_id, Ed)==MMC_INACTIVE){	      

    mmr_edge_sons(Mesh_id,Ed,edge_sons,NULL);
    ddr_add_edge_ovl(Mesh_id,edge_sons[0]);
    ddr_add_edge_ovl(Mesh_id,edge_sons[1]);

  }

  //printf("\t\tfinished adding edge %d to overlap\n",Ed);

  return(1);
}

/*---------------------------------------------------------
  ddr_add_sons_patch - to add antecedents to patch in a recursive manner
---------------------------------------------------------*/
int ddr_add_sons_patch(
  int Mesh_id,  /* in: mesh ID */
  int El,       /* in: element ID */
  int Ll,       /* in: length of array of elements */
  int* List_el   /* in/out: list of elements in a patch to be updated */
  );

/*---------------------------------------------------------
  ddr_create_patch - to create a patch surrounding element
---------------------------------------------------------*/
int ddr_create_patch(
  int Mesh_id, 
  int El_id,
  int Gen_lev, 
  int List_length, 
  int* List_el, 
  int* List_face_int, 
  int* List_face_bnd,
  int* List_edge_int, 
  int* List_edge_bnd, 
  int* List_vert_int, 
  int* List_vert_bnd
  )
{

  int el_faces[MMC_MAXELFAC+1]; /* element faces */
  int el_nodes[MMC_MAXELVNO+1];
  int i, iel, nel, iaux, max_node_id;
  int ifa, nfa, face_ngb[2];
  int ied, ned, num_edges, fa_edges[4];
  int ino, nno, edge_vert[2];

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* initiate the lists - assume all have the same length */
  for(iel=0;iel<=List_length;iel++) List_el[iel]=0; 
  for(i=0;i<=List_length;i++) List_face_int[i]=0; 
  for(i=0;i<=List_length;i++) List_face_bnd[i]=0; 
  for(i=0;i<=List_length;i++) List_edge_int[i]=0; 
  for(i=0;i<=List_length;i++) List_edge_bnd[i]=0; 
  for(i=0;i<=List_length;i++) List_vert_int[i]=0; 
  for(i=0;i<=List_length;i++) List_vert_bnd[i]=0; 

  /* start list of elements with the central one */
  List_el[0]=1;
  List_el[1]=El_id;

  max_node_id = mmr_get_max_node_id(Mesh_id);

  /* if elems arrays do not exist */
  if(elems==NULL || max_node_id!=max_node){

    /* if number of nodes changed */ 
    if(max_node_id!=max_node){

      for(ino=1;ino<=max_node;ino++) {
	if(elems[ino]!=NULL) free(elems[ino]);
      }
      free(elems);

    }
    max_node = max_node_id;

    elems = malloc((max_node_id+1)*sizeof(int *));
    for(ino=1;ino<=max_node_id;ino++){
      elems[ino]=NULL;  
    }
    /* prepare array of element connectivities - all vertices belonging to */
    /* owned elements get elems array with lists of all adjacent elements */
    ddr_vert_elems(Mesh_id, Gen_lev, elems); 

  }

  /* for all nodes of the central element */
  mmr_el_node_coor(Mesh_id, El_id, el_nodes, NULL);

  for(ino=0;ino<el_nodes[0];ino++){

    nno = el_nodes[ino+1];

    /* for all elements adjacent to the node */
    for(iel=1;iel<=elems[nno][0];iel++){ 
      
      nel=elems[nno][iel];

      iaux=ddr_put_list(nel, &List_el[1], List_length);
      if(iaux<0) List_el[0]++;
      if(iaux==0 || List_el[0]>List_length){
	printf("too big patch in ddr_create_patch - increase size of arrays\n");
	exit(-1);
      }

    }

  }

  iaux=List_el[0];
  for(iel=1;iel<=iaux;iel++) {
    ddr_add_sons_patch(Mesh_id, List_el[iel], List_length, List_el);
  }

/*kbw
  printf("for element %d\n list of %d elements in a patch:\n",
	 El_id, List_el[0]);
  for(iel=1;iel<=List_el[0];iel++) printf("%d  ",List_el[iel]);
  printf("\n");
/*kew*/

  /* consider faces and divide them into internal and boundary */
  for(iel=1;iel<=List_el[0];iel++) {

    nel = List_el[iel];

    mmr_el_faces(Mesh_id, nel, el_faces, NULL);

    for(ifa=0; ifa<el_faces[0]; ifa++){
      nfa=el_faces[ifa+1];

      iaux = ddr_chk_list(nfa, &List_face_bnd[1], List_length);

      if(iaux==0){

	iaux = ddr_chk_list(nfa, &List_face_int[1], List_length);

	if(iaux==0){

	  mmr_fa_eq_neig(Mesh_id, nfa, face_ngb, NULL, NULL); 

	  i=ddr_chk_list(face_ngb[0], &List_el[1], List_length);
	  if(i==0){
	    List_face_bnd[0]++;
	    List_face_bnd[List_face_bnd[0]]=nfa;
	  }
	  else{
	    i=ddr_chk_list(face_ngb[1], &List_el[1], List_length);
	    if(i==0){
	      List_face_bnd[0]++;
	      List_face_bnd[List_face_bnd[0]]=nfa;
	    }
	    else{

	      List_face_int[0]++;
	      List_face_int[List_face_int[0]]=nfa;

	    }

	  }

	  if(List_face_int[0]>=List_length || List_face_bnd[0]>=List_length){
	    printf("too much faces in ddr_create_patch - increase arrays\n");
	    exit(-1);
	  }

	} /* end if face not yet internal */

      } /* end if face not yet boundary */

    } /* end loop over element faces */
  
  } /* end loop over elements */

/*kbw
  printf("list of %d internal faces in a patch:\n",
	 List_face_int[0]);
  for(i=1;i<=List_face_int[0];i++) printf("%d  ", List_face_int[i]);
  printf("\n");
  printf("list of %d boundary faces in a patch:\n",
	 List_face_bnd[0]);
  for(i=1;i<=List_face_bnd[0];i++) printf("%d  ", List_face_bnd[i]);
  printf("\n");
/*kew*/

  /* consider edges and divide them into internal and boundary */

  /* for all boundary faces and their edges */
  for(ifa=1;ifa<=List_face_bnd[0];ifa++) {
    nfa = List_face_bnd[ifa];

    num_edges = mmr_fa_edges(Mesh_id, nfa, fa_edges, NULL);

    for(ied=0;ied<num_edges;ied++){

      ned = fa_edges[ied];

      /* edges of boundary faces belong to the boundary */
      iaux = ddr_put_list(ned, &List_edge_bnd[1], List_length);
      if(iaux<0) List_edge_bnd[0]++;
      if(iaux==0 || List_edge_bnd[0]>List_length){
	printf("too big patch in ddr_create_patch - increase size of arrays\n");
	exit(-1);
      }

   } /* end loop over face edges */

  } /* end loop over boundary faces */

  /* for all internal faces and their edges */
  for(ifa=1;ifa<=List_face_int[0];ifa++) {
    nfa = List_face_int[ifa];

    num_edges = mmr_fa_edges(Mesh_id, nfa, fa_edges, NULL);

    for(ied=0;ied<num_edges;ied++){

      ned = fa_edges[ied];

      /* if edge is not on the boundary */
      iaux = ddr_chk_list(ned, &List_edge_bnd[1], List_length);
      if(iaux==0){

	/* edges of boundary faces belong to the boundary */
	iaux = ddr_put_list(ned, &List_edge_int[1], List_length);
	if(iaux<0) List_edge_int[0]++;
	if(iaux==0 || List_edge_int[0]>List_length){
	  printf("too big patch in ddr_create_patch - increase size of arrays\n");
	  exit(-1);
	}

      } /* end if edge is not on the boundary */

    } /* end loop over face edges */

  } /* end loop over boundary faces */

/*kbw
  printf("list of %d internal edges in a patch:\n",
	 List_edge_int[0]);
  for(i=1;i<=List_edge_int[0];i++) printf("%d  ", List_edge_int[i]);
  printf("\n");
  printf("list of %d boundary edges in a patch:\n",
	 List_edge_bnd[0]);
  for(i=1;i<=List_edge_bnd[0];i++) printf("%d  ", List_edge_bnd[i]);
  printf("\n");
/*kew*/

  /* consider vertices and divide them into internal and boundary */

  for(ied=1;ied<=List_edge_bnd[0];ied++) {
    ned = List_edge_bnd[ied];

    mmr_edge_nodes( Mesh_id, ned, edge_vert);

    nno=edge_vert[0];

    /* vertices of boundary edges belong to the boundary */
    iaux = ddr_put_list(nno, &List_vert_bnd[1], List_length);
      if(iaux<0) List_vert_bnd[0]++;
      if(iaux==0 || List_vert_bnd[0]>List_length){
	printf("too big patch in ddr_create_patch - increase size of arrays\n");
	exit(-1);
      }

    nno=edge_vert[1];

    /* vertices of boundary edges belong to the boundary */
    iaux = ddr_put_list(nno, &List_vert_bnd[1], List_length);
      if(iaux<0) List_vert_bnd[0]++;
      if(iaux==0 || List_vert_bnd[0]>List_length){
	printf("too big patch in ddr_create_patch - increase size of arrays\n");
	exit(-1);
      }

  }

  for(ied=1;ied<=List_edge_int[0];ied++) {
    ned = List_edge_int[ied];

    nno=edge_vert[0];

    /* if vertex is not on the boundary */
    iaux = ddr_chk_list(ned, &List_vert_bnd[1], List_length);
    if(iaux==0){

      /* vertices of boundary edges belong to the boundary */
      iaux = ddr_put_list(ned, &List_vert_int[1], List_length);
      if(iaux<0) List_vert_int[0]++;
      if(iaux==0 || List_vert_int[0]>List_length){
	printf("too big patch in ddr_create_patch - increase size of arrays\n");
	exit(-1);
      }

    } /* end if vertex is not on the boundary */

    nno=edge_vert[1];

    /* if vertex is not on the boundary */
    iaux = ddr_chk_list(ned, &List_vert_bnd[1], List_length);
    if(iaux==0){

      /* vertices of boundary edges belong to the boundary */
      iaux = ddr_put_list(ned, &List_vert_int[1], List_length);
      if(iaux<0) List_vert_int[0]++;
      if(iaux==0 || List_vert_int[0]>List_length){
	printf("too big patch in ddr_create_patch - increase size of arrays\n");
	exit(-1);
      }

    } /* end if vertex is not on the boundary */

  } /* end loop over internal edges */

/*kbw
  printf("list of %d internal vertices in a patch:\n",
	 List_vert_int[0]);
  for(i=1;i<=List_vert_int[0];i++) printf("%d  ", List_vert_int[i]);
  printf("\n");
  printf("list of %d boundary vertices in a patch:\n",
	 List_vert_bnd[0]);
  for(i=1;i<=List_vert_bnd[0];i++) printf("%d  ", List_vert_bnd[i]);
  printf("\n");
/*kew*/

  return(1);
}

/*---------------------------------------------------------
  ddr_add_sons_patch - to add antecedents to patch in a recursive manner
---------------------------------------------------------*/
int ddr_add_sons_patch(
  int Mesh_id,  /* in: mesh ID */
  int El,       /* in: element ID */
  int Ll,       /* in: length of array of elements */
  int* List_el   /* in/out: list of elements in a patch to be updated */
  )
{

  int ison, elson, iaux;
  int elsons[MMC_MAXELSONS+1];

/*++++++++++++++++ executable statements ++++++++++++++++*/

  mmr_el_fam(Mesh_id,El,elsons,NULL);

  for(ison=1;ison<=elsons[0];ison++){

    elson=elsons[ison];
	  
/*kbw
      printf("son %d, element %d", ison, elson);
/*kew*/

    iaux=ddr_put_list(elson, &List_el[1], Ll);
    if(iaux<0) List_el[0]++;
    if(List_el[0]>Ll){
      printf("too big patch in ddr_create_patch - increase size of arrays\n");
      exit(-1);
    }

    /* add grandsons to overlap */
    ddr_add_sons_patch(Mesh_id, elson, Ll, List_el);

  }	

  return(1);
}

/*** utilities ***/

/*---------------------------------------------------------
ddr_put_list - to put Num on the list List with length Ll 
	(filled with numbers and zeros at the end)
---------------------------------------------------------*/
int ddr_put_list( /* returns*/
		/*  >0 - position already occupied on the list */
             	/*  <0 - position at which put on the list */
            	/*   0 - list full, not found on the list */
	int Num, 	/* in: number to put on the list */
	int* List, 	/* in: list */
	int Ll		/* in: total list's lengths */
	)
{

  int i, il;
  
  for(i=0;i<Ll;i++){
    if((il=List[i])==0) break;
    /* found on the list on (i+1) position */
    if(Num==il) return(i+1);
  }
  /* if list is full return error message */
  if(i==Ll) return(0);
  /* update the list and return*/
  List[i]=Num;
  return(-(i+1));
}

/*---------------------------------------------------------
ddr_chk_list - to check whether Num is on the list List
	with length Ll
---------------------------------------------------------*/
int ddr_chk_list(	/* returns: */
			/* >0 - position on the list */
            		/* 0 - not found on the list */
	int Num, 	/* number to be checked */
	int* List, 	/* list of numbers */
	int Ll		/* length of the list */
	)
{

  int i, il;

  for(i=0;i<Ll;i++){
    if((il=List[i])==0) break;
    /* found on the list on (i+1) position */
    if(Num==il) return(i+1);
  }
  /* not found on the list */
  return(0);
}

