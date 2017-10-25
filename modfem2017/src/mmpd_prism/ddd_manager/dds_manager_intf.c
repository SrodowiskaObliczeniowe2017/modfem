/************************************************************************
File dds_manager_intf - to manage all tasks 

Contains definitions of routines:   
  ddr_create_subdomains - to decompose the mesh and create subdomains
  ddr_mesh_ent_ovl - to create tables (for each mesh entity type) of local 
                     numbers and storing processors for each mesh entity owned 
                     by a given processor but stored also on other processors
  ddr_free_mesh_ovl - to free data structure with overlap info
  ddr_update_ipid - to update interprocessor ID for all mesh entities
                    using interprocessor communication
  ddr_exchange_list_ref - to update the list of refined elements
                          using inter-processor communication
  ddr_update_ref_list - to update list of refined elements due to irregularity 
                      constraint using inter-processor communication

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

/* interface of the approximation module */
#include "aph_intf.h"

/* parallel communication interface specification */
#include "pch_intf.h"

/* parallel mesh manipulation interface specification */
#include "mmph_intf.h"

/* parallel mesh manipulation interface for prismatic elements */
#include "../mmph_prism.h"

/* internal info for domain decomposition manager */
#include "./ddh_manager.h"

/* internal info for front based domain decomposition */
#include "./ddh_front.h"


/* CONSTANTS */
#define DDC_MAX_NUM_PROC MMPC_MAX_NUM_SUB;  /* maximal number of processors */


/* arrays for domain decomposition */
static int *nr_sub_elem, **l_sub_elem;

/*---------------------------------------------------------
  ddr_set_edge_sons_owner - local utility
---------------------------------------------------------*/
int ddr_set_edge_sons_owner(
  int Mesh_id, /* in: as is */
  int Ed,      /* in: edge whose sons will be checked */
  int Sub      /* in: subdomain ID to be set for edges */
  );

/*---------------------------------------------------------
  ddr_create_subdomains - to decompose the mesh and create subdomains
---------------------------------------------------------*/
int ddr_create_subdomains( /* returns: >=0 - success code, <0 - error code */
  int Mesh_id,      /* in: mesh ID  */
  int Control        /* in: generation level as basis for decomposition */
  )
{

  int i,j,k, ient, ifa, jfa, kfa, nfa, ison, owner;
  int iaux, nr_sub, isub, iel, nel, ned, ied, ino, nr_elem, gen_lev, mesh_id;
  int el_faces[MMC_MAXELFAC+1], el_nodes[MMC_MAXELVNO+1], el_sons[9];
  int num_edges, fa_sons[5], fa_edges[5], edge_sons[2], edge_nodes[2];

  double daux;
  const int zero =  0;

/*++++++++++++++++ executable statements ++++++++++++++++*/
  
  /* associate the suitable mesh with approximation field */ 
  mesh_id = Mesh_id;

  /* simplifying assumption */
  nr_sub = mmpv_nr_sub;
  nr_sub_elem = (int *)malloc(nr_sub*sizeof(int));
  l_sub_elem = (int **)malloc(nr_sub*sizeof(int *));
  l_sub_elem[0] = NULL;

/* check weights - if existent 
#ifdef DEBUG
  if(L_weight!=NULL){
    daux=0.0;
    for(isub=0; isub<Nr_sub; isub++){
      daux+=L_weight[isub];
    }
    daux /= Nr_sub;
    if(fabs(daux-1.0)>1e-3){
      printf("Wrong weights in the input to create subdomains!\n");
      exit(-1);
    }
  }
#endif
end checking weights */

  /* temporary setting for testing performance */
/*kbw
  if(nr_sub==1){
    printf("single subdomain: changing for testing (1-no, >1 - new number)\n");
    scanf("%d",&nr_sub); getchar();
    if(nr_sub>1){
      free(nr_sub_elem);
      free(l_sub_elem);
      nr_sub_elem = (int *)malloc(nr_sub*sizeof(int));
      l_sub_elem = (int **)malloc(nr_sub*sizeof(int *));
      l_sub_elem[0] = NULL;
    }
  }
/*kew*/

  /* reset ownership info */
  if(nr_sub==1) iaux=mmpv_my_proc_id;
  else iaux=0;

/*kbw
  printf("starting mesh partitioning: nr_sub = %d, my_proc_id = %d, iaux = %d\n",
	 nr_sub,mmpv_my_proc_id,iaux); 

/*kew*/

  i=0;
  while((i=mmr_get_next_elem_all(mesh_id, i))!=0){

    mmpr_el_set_owner(mesh_id,i,iaux);
    mmpr_el_set_id_at_owner(mesh_id,i,i);
  }
/*kbw
    printf("Initialized ownership and id_at_owner for %d elems \n", i);
/*kew*/
  i=0;
  while((i=mmr_get_next_face_all(mesh_id, i))!=0){
    mmpr_fa_set_owner(mesh_id,i,iaux);
    mmpr_fa_set_id_at_owner(mesh_id,i,i);
  }
/*kbw
    printf("Initialized ownership and id_at_owner for %d faces \n", i);
/*kew*/
  i=0;
  while((i=mmr_get_next_edge_all(mesh_id, i))!=0){
    mmpr_ed_set_owner(mesh_id,i,iaux);
    mmpr_ed_set_id_at_owner(mesh_id,i,i);
  }
/*kbw
    printf("Initialized ownership and id_at_owner for %d edges \n", i);
/*kew*/
  i=0;
  while((i=mmr_get_next_node_all(mesh_id, i))!=0){
    mmpr_ve_set_owner(mesh_id,i,iaux);
    mmpr_ve_set_id_at_owner(mesh_id,i,i);
  }
/*kbw
    printf("Initialized ownership and id_at_owner for %d nodes \n", i);
/*kew*/

  if(nr_sub==1) return(1);


  /*kb!!! begin test usage of front partitioner !!! */
  /* for example loop over subdoamins */
  nr_elem = mmr_get_nr_elem(mesh_id);
  
  nr_sub_elem[0] = nr_elem/nr_sub;
  if(nr_sub_elem[0]*nr_sub<nr_elem){
    nr_sub_elem[0] = nr_elem/nr_sub+1;
  }

  for(i=1;i<nr_sub-1;i++) nr_sub_elem[i] = nr_sub_elem[0];

  nr_sub_elem[nr_sub-1] = nr_elem - (nr_sub-1)*nr_sub_elem[0];

/*kb!!!
  printf("set generation level as basis for domain decomposition: ");
  scanf("%d",&gen_lev);getchar(); 
kbw!!!*/
  gen_lev = Control;
      
  if( l_sub_elem[0] != NULL) {
    for(isub=0;isub<nr_sub;isub++) free(l_sub_elem[isub]);
  }

  for(isub=0;isub<nr_sub;isub++){

/*kbw
    printf("\nmesh %d, gen_lev %d, nr_sub %d, isub %d, nr_sub_elem %d\n",
	   mesh_id, gen_lev, nr_sub, isub, nr_sub_elem[isub]);

/*kew*/

    ddr_create_subdomain(mesh_id, gen_lev, nr_sub, isub+1,  
			 &nr_sub_elem[isub], &l_sub_elem[isub]);
/*kbw
#ifdef DEBUG
    printf("subdomain %d, nr_elem %d, list_elem:\n",
	   isub, nr_sub_elem[isub]);
    //    for(i=0;i<nr_sub_elem[isub];i++){
    //      printf("%d ", l_sub_elem[isub][i]);
    //    }
    printf("\n");
#endif
/*kew*/


  }
  /*kb!!! end test usage of front partitioner !!! */

  for(isub=1;isub<=nr_sub;isub++){

    for(iel=0;iel<nr_sub_elem[isub-1];iel++){

      nel = l_sub_elem[isub-1][iel];

      if(mmr_el_gen(mesh_id, nel)==gen_lev){

	/*kb!!! when setting the ownership for boundary entities the rule is: */
	/*      the owner is the subdoamin with smaller ID. */
	
	/* set ownership for elements */
	ddr_set_owner_elem(mesh_id, nel, isub);
	

/*kbw
#ifdef DEBUG
      if(ddr_get_owner(mesh_id, MMC_ELEMENT, l_sub_elem[isub-1][iel]) != isub){
	printf("subdoamin %d, setting owner for element %d to %d (%d)\n",
	       isub, l_sub_elem[isub-1][iel], 
	       l_sub_elem[isub-1][iel]*MMPC_MAX_NUM_SUB+isub,
	       mmpr_el_owner(mesh_id, l_sub_elem[isub-1][iel]));
      }
#endif
/*kew*/

      }

    } /* end loop over all subdomain elements */

  }


  iel=0;
  while((iel=mmr_get_next_elem_all(mesh_id, iel))!=0){

    if(mmpr_el_owner(mesh_id,iel)==0){
      printf("Error : element %d - without owner \n",iel);
      exit(-1);
    }

  }
  ifa=0;
  while((ifa=mmr_get_next_face_all(mesh_id, ifa))!=0){
    if(mmpr_fa_owner(mesh_id,ifa)==0){
      printf("Error : face %d - without owner \n",ifa);
      exit(-1);
    }
  }
  ied=0;
  while((ied=mmr_get_next_edge_all(mesh_id, ied))!=0){
    if(mmpr_ed_owner(mesh_id,ied)==0){
      printf("Error : edge %d - without owner \n",ied);
      exit(-1);
    }
  }
  i=0;
  while((i=mmr_get_next_node_all(mesh_id, i))!=0){
    if(mmpr_ve_owner(mesh_id,i)==0){
      printf("Error : vertex %d - without owner \n",i);
      exit(-1);
    }
  }


/*kbw
  printf("Subdomain %d\n", mmpv_my_proc_id);
  printf("Active elements (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_elem_all(mesh_id, ient))!=0){
    if(mmr_el_status(mesh_id,ient)==MMC_ACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_el_owner(mesh_id,ient),
	     mmpr_el_id_at_owner(mesh_id,ient));
    }
  }
  printf("\n");
  printf("Inactive elements (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_elem_all(mesh_id, ient))!=0){
    if(mmr_el_status(mesh_id,ient)==MMC_INACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_el_owner(mesh_id,ient),
	     mmpr_el_id_at_owner(mesh_id,ient));
    }
  }
  printf("\n");
  printf("Active faces (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_face_all(mesh_id, ient))!=0){
    if(mmr_fa_status(mesh_id,ient)==MMC_ACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_fa_owner(mesh_id,ient),
	     mmpr_fa_id_at_owner(mesh_id,ient));
    }
  }
  printf("\n");
  printf("Inactive faces (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_face_all(mesh_id, ient))!=0){
    if(mmr_fa_status(mesh_id,ient)==MMC_INACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_fa_owner(mesh_id,ient),
	     mmpr_fa_id_at_owner(mesh_id,ient));
    }
  }
  printf("\n");
  printf("Active edges (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_edge_all(mesh_id, ient))!=0){
    if(mmr_edge_status(mesh_id,ient)==MMC_ACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_ed_owner(mesh_id,ient),
	     mmpr_ed_id_at_owner(mesh_id,ient));
    }
  }
  printf("\n");
  printf("Inactive edges (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_edge_all(mesh_id, ient))!=0){
    if(mmr_edge_status(mesh_id,ient)==MMC_INACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_ed_owner(mesh_id,ient),
	     mmpr_ed_id_at_owner(mesh_id,ient));
    }
  }
  printf("\n");
  printf("Vertices (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_node_all(mesh_id, ient))!=0){
    printf("%d (%d %d), ",
	   ient, mmpr_ve_owner(mesh_id,ient),
	   mmpr_ve_id_at_owner(mesh_id,ient));
  }
  printf("\n");
/*kew*/


  /* specify internal (owned) and ghost (alien) entities */ 
  iaux = 2; // two element overlap
  ddr_create_overlap(mesh_id, iaux, gen_lev);

/*kbw
  sleep(mmpv_my_proc_id);
  printf("Subdomain %d\n", mmpv_my_proc_id);
  printf("Owned Active elements (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_elem_all(mesh_id, ient))!=0){
    if(mmpr_el_owner(mesh_id,ient)==mmpv_my_proc_id){
    if(mmr_el_status(mesh_id,ient)==MMC_ACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_el_owner(mesh_id,ient),
	     mmpr_el_id_at_owner(mesh_id,ient));
    }
    }
  }
  printf("\n");
  printf("Owned Inactive elements (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_elem_all(mesh_id, ient))!=0){
    if(mmpr_el_owner(mesh_id,ient)==mmpv_my_proc_id){
    if(mmr_el_status(mesh_id,ient)==MMC_INACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_el_owner(mesh_id,ient),
	     mmpr_el_id_at_owner(mesh_id,ient));
    }
    }
  }
  printf("\n");
  printf("Owned Active faces (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_face_all(mesh_id, ient))!=0){
    if(mmpr_fa_owner(mesh_id,ient)==mmpv_my_proc_id){
    if(mmr_fa_status(mesh_id,ient)==MMC_ACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_fa_owner(mesh_id,ient),
	     mmpr_fa_id_at_owner(mesh_id,ient));
    }
    }
  }
  printf("\n");
  printf("Owned Inactive faces (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_face_all(mesh_id, ient))!=0){
    if(mmpr_fa_owner(mesh_id,ient)==mmpv_my_proc_id){
    if(mmr_fa_status(mesh_id,ient)==MMC_INACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_fa_owner(mesh_id,ient),
	     mmpr_fa_id_at_owner(mesh_id,ient));
    }
    }
  }
  printf("\n");
  printf("Owned Active edges (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_edge_all(mesh_id, ient))!=0){
    if(mmpr_ed_owner(mesh_id,ient)==mmpv_my_proc_id){
    if(mmr_edge_status(mesh_id,ient)==MMC_ACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_ed_owner(mesh_id,ient),
	     mmpr_ed_id_at_owner(mesh_id,ient));
    }
    }
  }
  printf("\n");
  printf("Owned Inactive edges (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_edge_all(mesh_id, ient))!=0){
    if(mmpr_ed_owner(mesh_id,ient)==mmpv_my_proc_id){
    if(mmr_edge_status(mesh_id,ient)==MMC_INACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_ed_owner(mesh_id,ient),
	     mmpr_ed_id_at_owner(mesh_id,ient));
    }
    }
  }
  printf("\n");
  printf("Owned Vertices (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_node_all(mesh_id, ient))!=0){
    if(mmpr_ve_owner(mesh_id,ient)==mmpv_my_proc_id){
    printf("%d (%d %d), ",
	   ient, mmpr_ve_owner(mesh_id,ient),
	   mmpr_ve_id_at_owner(mesh_id,ient));
    }
  }
  printf("\n");
/*kew*/
/*kbw
  printf("Subdomain %d\n", mmpv_my_proc_id);
  printf("Alien Active elements (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_elem_all(mesh_id, ient))!=0){
    if(mmpr_el_owner(mesh_id,ient) != mmpv_my_proc_id){
    if(mmr_el_status(mesh_id,ient)==MMC_ACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_el_owner(mesh_id,ient),
	     mmpr_el_id_at_owner(mesh_id,ient));
    }
    }
  }
  printf("\n");
  printf("Alien Inactive elements (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_elem_all(mesh_id, ient))!=0){
    if(mmpr_el_owner(mesh_id,ient) != mmpv_my_proc_id){
    if(mmr_el_status(mesh_id,ient)==MMC_INACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_el_owner(mesh_id,ient),
	     mmpr_el_id_at_owner(mesh_id,ient));
    }
    }
  }
  printf("\n");
  printf("Alien Active faces (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_face_all(mesh_id, ient))!=0){
    if(mmpr_fa_owner(mesh_id,ient) != mmpv_my_proc_id){
    if(mmr_fa_status(mesh_id,ient)==MMC_ACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_fa_owner(mesh_id,ient),
	     mmpr_fa_id_at_owner(mesh_id,ient));
    }
    }
  }
  printf("\n");
  printf("Alien Inactive faces (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_face_all(mesh_id, ient))!=0){
    if(mmpr_fa_owner(mesh_id,ient) != mmpv_my_proc_id){
    if(mmr_fa_status(mesh_id,ient)==MMC_INACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_fa_owner(mesh_id,ient),
	     mmpr_fa_id_at_owner(mesh_id,ient));
    }
    }
  }
  printf("\n");
  printf("Alien Active edges (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_edge_all(mesh_id, ient))!=0){
    if(mmpr_ed_owner(mesh_id,ient) != mmpv_my_proc_id){
    if(mmr_edge_status(mesh_id,ient)==MMC_ACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_ed_owner(mesh_id,ient),
	     mmpr_ed_id_at_owner(mesh_id,ient));
    }
    }
  }
  printf("\n");
  printf("Alien Inactive edges (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_edge_all(mesh_id, ient))!=0){
    if(mmpr_ed_owner(mesh_id,ient) != mmpv_my_proc_id){
    if(mmr_edge_status(mesh_id,ient)==MMC_INACTIVE){
      printf("%d (%d %d), ",
	     ient, mmpr_ed_owner(mesh_id,ient),
	     mmpr_ed_id_at_owner(mesh_id,ient));
    }
    }
  }
  printf("\n");
  printf("Alien Vertices (owner, local ID): ");
  ient=0;
  while((ient=mmr_get_next_node_all(mesh_id, ient))!=0){
    if(mmpr_ve_owner(mesh_id,ient) != mmpv_my_proc_id){
    printf("%d (%d %d), ",
	   ient, mmpr_ve_owner(mesh_id,ient),
	   mmpr_ve_id_at_owner(mesh_id,ient));
    }
  }
  printf("\n");
/*kew*/

   return(0);
}



/*---------------------------------------------------------
  ddr_set_owner_elem - to add an element to the Sub with all
                     its antecedents in a recursive manner
---------------------------------------------------------*/
int ddr_set_owner_elem(
  int Mesh_id, /* in: mesh ID */
  int El,       /* in: element ID */
  int Sub
  )
{

  int i,j,k, ient, ifa, jfa, kfa, nfa, ison, owner;
  int iaux, nr_sub, isub, iel, nel, ned, ied, ino, nr_elem, gen_lev, mesh_id;
  int elsons[MMC_MAXELSONS+1];
  int el_faces[MMC_MAXELFAC+1];

/*++++++++++++++++ executable statements ++++++++++++++++*/

 /*kbw
  if(El==4105||El==323){
    printf("adding element %d to SUB %d , father %d\n",El, Sub,
	   mmr_el_fam(Mesh_id, El, NULL, NULL));
  }
/*kew*/

  owner=mmpr_el_owner(Mesh_id,El);
  if(owner==0){

    mmpr_el_set_owner(Mesh_id,El,Sub);
 
    /* add faces to Sub */
    mmr_el_faces(mesh_id, El, el_faces, NULL);
    
    for(nfa=0;nfa<el_faces[0];nfa++) {
      
      ifa = el_faces[nfa+1];
      
      owner=mmpr_fa_owner(Mesh_id,ifa);
      
      /*kbw
	if((ifa>=1976&&ifa<=1979)||ifa==416){
	printf("elem %d, face %d, owner %d\n",El,ifa,owner);
	}
	/*kew*/
      
      
      if(owner==0) {
	ddr_set_owner_face(Mesh_id,ifa,Sub);
      }
      
    }
    
    if(mmr_el_status(mesh_id, El)==MMC_INACTIVE){
      
      /* add sons to Sub */
      mmr_el_fam(Mesh_id,El,elsons,NULL);
      
      for(ison=1;ison<=elsons[0];ison++){
	ddr_set_owner_elem(Mesh_id,elsons[ison],Sub);
      }
      
    }
  }
  return(1);
}


/*---------------------------------------------------------
  ddr_set_owner_face - to add a face to the Sub with all
                     its antecedents in a recursive manner
---------------------------------------------------------*/
int ddr_set_owner_face(
  int Mesh_id, /* in: mesh ID */
  int Fa,       /* in: face ID */
  int Sub
  )
{

  int i,j,k, jfa, kfa, ison, owner;
  int ned, ied;
  int num_edges, fa_sons[5], fa_edges[5];

/*++++++++++++++++ executable statements ++++++++++++++++*/

  //printf("\tadding face %d to Sub\n",Fa);

  owner=mmpr_fa_owner(Mesh_id,Fa);


/*kbw
	if(Fa>=1976&&Fa<=1979){
	  printf("adding to Sub face %d, owner %d\n",Fa,owner);
	}
/*kew*/

  if(owner==0) {
    
    mmpr_fa_set_owner(Mesh_id,Fa,Sub);
 
    /* add edges to Sub */
    num_edges = mmr_fa_edges(Mesh_id, Fa, fa_edges, NULL);
    
    for(ned=0;ned<num_edges;ned++){
      
      ied = fa_edges[ned];
      
      owner=mmpr_ed_owner(Mesh_id,ied);
      if(owner==0) {

/*kbw
	if(ied==15443||ied==15444){
	  printf("Edge %d, owned by %d - external to %d - ADDING TO SUB %d\n",  
		 ied, mmpr_ed_owner(Mesh_id,  ied), mmpv_my_proc_id, Sub);
	  getchar();
	}
/*kew*/
	ddr_set_owner_edge(Mesh_id,ied,Sub);
      }
      
    }
    
    /* add children to Sub */
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
	
	if(owner==0){
	  ddr_set_owner_face(Mesh_id,kfa,Sub);
	}
	
      }
    }
  }

  return(1);
}


/*---------------------------------------------------------
  ddr_set_owner_edge - to add an edge to the Sub with all
                     its antecedents in a recursive manner
---------------------------------------------------------*/
int ddr_set_owner_edge(
  int Mesh_id, /* in: mesh ID */
  int Ed,       /* in: edge ID */
  int Sub
  )
{

  int i,j,k, ison, owner;
  int num_edges, edge_sons[2], edge_nodes[2];

/*++++++++++++++++ executable statements ++++++++++++++++*/

  //printf("\t\tadding edge %d to Sub\n",Ed);
/*kbw
  if(Ed==15443||Ed==15444){
      printf("Edge %d, owned by %d - external to %d - adding to Sub %d\n", Ed, 
	     mmpr_ed_owner(Mesh_id,  Ed), mmpv_my_proc_id, Sub);
      getchar();
      }
/*kew*/

  owner=mmpr_ed_owner(Mesh_id,Ed);
  if(owner==0) {

    mmpr_ed_set_owner(Mesh_id,Ed,Sub);

    /* add nodes to Sub */	      
    mmr_edge_nodes(Mesh_id, Ed, edge_nodes);
	      
    if(mmpr_ve_owner(Mesh_id, edge_nodes[0]) == 0){
      mmpr_ve_set_owner(Mesh_id, edge_nodes[0], Sub);
    }
    if(mmpr_ve_owner(Mesh_id, edge_nodes[1]) == 0){
      mmpr_ve_set_owner(Mesh_id, edge_nodes[1], Sub);
    }
  
    /* add children to Sub */
    if(mmr_edge_status(Mesh_id, Ed)==MMC_INACTIVE){	      

      mmr_edge_sons(Mesh_id,Ed,edge_sons,NULL);
      ddr_set_owner_edge(Mesh_id,edge_sons[0],Sub);
      ddr_set_owner_edge(Mesh_id,edge_sons[1],Sub);

    }

  }

  //printf("\t\tfinished adding edge %d to Sub\n",Ed);

  return(1);
}
