/************************************************************************
File

Contains definitions of routines:   

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

/* domain decomposition interface specification */
#include "mmph_intf.h"

/* internal info for front based domain decomposition */
#include "./ddh_front.h"

#define MAXELNO 20
static int nr_node_front;
static int **elems=NULL, *sub_ind_el=NULL, *front_ind=NULL;

static int put_list(int Num, int* List, int Ll);


/*---------------------------------------------------------
  dds_prepare_nodes
---------------------------------------------------------*/
int dds_prepare_nodes(
		      int Mesh_id,      /* in: mesh ID */
		      int Gen_lev /* generation to consider */
)
{

  int max_elem_id, max_node_id;
  int i, nel, nno, ino, ipr;
  int el_nodes[MMC_MAXELVNO+1];

/*++++++++++++++++ executable statements ++++++++++++++++*/
  
  max_elem_id = mmr_get_max_elem_id(Mesh_id);
  sub_ind_el = malloc((max_elem_id+1)*sizeof(int));
  for(i=0;i<=max_elem_id;i++) sub_ind_el[i]=0;
  
  max_node_id = mmr_get_max_node_id(Mesh_id);
  front_ind = malloc((max_node_id+1)*sizeof(int));
  for(i=0;i<=max_node_id;i++) front_ind[i]=0;
  elems = malloc((max_node_id+1)*sizeof(int *));
  for(ino=0;ino<=max_node_id;ino++){
    elems[ino] = malloc(MAXELNO*sizeof(int));
    for(i=0;i<MAXELNO;i++) elems[ino][i]=0;
  }
  
  /* loop over elements */
  nel=0;
  while((nel=mmr_get_next_elem_all(Mesh_id, nel))!=0){

    if(mmr_el_gen(Mesh_id,nel)==Gen_lev){
    
      mmr_el_node_coor(Mesh_id,nel,el_nodes,NULL);
    
      for(nno=1;nno<=el_nodes[0];nno++){
	ino=el_nodes[nno];
	i=put_list(nel,elems[ino],MAXELNO);
	if(i<0){
	  printf("Something wrong with filling NODES->elems for node %d\n",ino);
	  for(ipr=0;ipr<MAXELNO;ipr++) 
	    printf("%d\n",elems[ino][ipr]);
	  printf("\n");
	  getchar();
	}

/*kbw
#ifdef DEBUG
      printf("node %d (%d), added element %d, list:",
	     ino, el_nodes[nno], nel);
      ipr=0;
      while(elems[ino][ipr]!=0) {
	printf("%d ",elems[ino][ipr]);
	ipr++;
      }
      printf("\n");
#endif
/*kew*/

      }
    }
  }
  return(0);
}


/*---------------------------------------------------------
  ddr_create_subdomain - to create a subdomain
---------------------------------------------------------*/
int ddr_create_subdomain (
  int Mesh_id,	     /* in: mesh ID */
  int Gen_lev,       /* in: generation level as basis for decomposition */
  int Nr_sub,
  int Sub_id,
  int* Nr_sub_elem,
  int** L_sub_elem
)
{

  int max_elem_id, max_node_id;
  int size, size_tot, size_max, iell, iel, ielaux, nno, nodaux, ino;
  int el_nodes[MMC_MAXELVNO+1];
  int i, iaux;

/***************** beginning of the great loop *****************/

  if(elems==NULL) dds_prepare_nodes(Mesh_id,Gen_lev);

  size=0; size_tot=0;
  size_max = *Nr_sub_elem;

/*kbw
#ifdef DEBUG
  printf("starting creation for subdoamin %d, gen_lev %d, size_max %d\n",
	 Sub_id, Gen_lev, size_max);
#endif
/*kew*/

  /*initialize front*/ 
  if(Sub_id==1){
    /*kb!!! - should ensure initial node is on the boundary */ 
    front_ind[1]=1;
    nr_node_front=1;
  }

  /* while the size of a subdomain
     does not exceed the optimal and the front is not empty */
  while((size<size_max||Sub_id==Nr_sub)&&nr_node_front>0){

    /* pick the next node from the front*/
    ino=choose_front(Mesh_id);

/*kbw
#ifdef DEBUG
    printf("chosen node %d, nr_node_front %d, size %d\n",
	   ino, nr_node_front, size);
#endif
/*kew*/
    
    /* option: add number of DOF to the size */
    
    /* indicate that node has been added ??? */

    /* for each element sharing the node */
    for(iell=0;iell<MAXELNO&&(ielaux=elems[ino][iell])!=0;iell++){
      
      iel = ielaux;

      if(mmr_el_gen(Mesh_id,iel)!=Gen_lev){
	printf("Wrong generation. Exiting!\n");
	exit(-1);
      }
      
/*kbw
      if(iel==34){
	printf("considering element %d, mesh_gen %d (active %d) - status %d\n",
	       iel, mmr_el_gen(Mesh_id, iel), ielaux, sub_ind_el[iel]);
      }
/*kew*/

      /* if element not yet considered */
      if(sub_ind_el[iel]==0) {
	
	/* put it on the list of subdomain elements */
	sub_ind_el[iel]=Sub_id;
	size_tot++;
	if(mmr_el_status(Mesh_id,iel)==MMC_ACTIVE) size++;

	
/*kbw
      if(iel==34){
      printf("element %d, gen %d  added - status %d, size %d, size_tot %d\n",
	     iel, mmr_el_gen(Mesh_id, iel), sub_ind_el[iel], size, size_tot);
      }
/*kew*/

	/* get element's nodes */
	mmr_el_node_coor(Mesh_id,iel,el_nodes,NULL);
	
	/* check whether element has constrained nodes */
	for(nno=1;nno<=el_nodes[0];nno++){
	  /* for each of element's nodes */
	  
	  nodaux=el_nodes[nno];
	  
	  if(nodaux!=ino){

	    upd_front(Mesh_id, nodaux, Sub_id, Gen_lev);
	    
	  }
        }
	/* ^finished with all nodes of a given element sharing a front node */
	
	
	/* for inactive elements put all its antecedents and their nodes */
        if(mmr_el_status(Mesh_id,iel)==MMC_INACTIVE) {
	  put_sons(Mesh_id, iel, Sub_id, &size, &size_tot);
	}

	/* break when size reaches maximum */
	if(size>=size_max) break;

      } /* end if element not yet considered */
    }
    /* ^ finished with all elements sharing a front node */
    
    /* check whether the node should be deleted from the front */
    upd_front(Mesh_id, ino, Sub_id, Gen_lev);
    
    /* check for number of elements */
    if(size>size_max){
/*kbw
      printf("size for subdomain %d exceeds maximum ( %d > %d )\n",
	     Sub_id, size, size_max);
/*kew*/
    }
    
  } 
  /* ^ the end of while loop */
  
  *Nr_sub_elem = 0; iaux=0;
  *L_sub_elem = (int *)malloc(size_tot*sizeof(int));
  max_elem_id = mmr_get_max_elem_id(Mesh_id);
  for(i=1;i<=max_elem_id;i++){
    if(sub_ind_el[i]==Sub_id){
      (*L_sub_elem)[*Nr_sub_elem] = i;
      (*Nr_sub_elem)++;
      if(mmr_el_status(Mesh_id,i)==MMC_ACTIVE) iaux++;
    }
  }

#ifdef DEBUG
  if( *Nr_sub_elem != size_tot || size != iaux ){
    printf("error 5723 in create subdoamin! %d %d %d %d\n",
	   *Nr_sub_elem, size_tot,  size, iaux);
    exit(-1);
  }
#endif

  
  if(Sub_id==Nr_sub){

#ifdef DEBUG
    for(i=1;i<=max_elem_id;i++){
      if(sub_ind_el[i]<=0 && mmr_el_status(Mesh_id,i) != MMC_FREE
	 && mmr_el_gen(Mesh_id,i) >= Gen_lev ){
	printf("element %d, gen %d left after decomposition\n",
	       i, mmr_el_gen(Mesh_id,i) );
	exit(-1);
      }      
    }
#endif

    max_node_id = mmr_get_max_node_id(Mesh_id);
    free(sub_ind_el);
    free(front_ind);
    for(ino=0;ino<=max_node_id;ino++){
      free(elems[ino]);
    }
    free(elems);
    elems=NULL;
  }

  return(0);
}


/*---------------------------------------------------------
choose_front - to pick a node from the front according
		to prescribed weights
---------------------------------------------------------*/
int choose_front(
  int Mesh_id	     /* in: mesh ID */
  )
{

  int max_node_id, ino, max_val=10000000, max_ind=0;
  double coor[3], x_min=1.e10;

/***************** beginning of the great loop *****************/

  max_node_id = mmr_get_max_node_id(Mesh_id);
  for(ino=1;ino<=max_node_id;ino++){

/*kbw
#ifdef DEBUG
    printf("testing node %d, front_ind %d, max_val %d, max_ind %d\n",
           ino,front_ind[ino],max_val,max_ind);
#endif
/*kew*/

/*kb!!! - old simple
    if(front_ind[ino]>0 && front_ind[ino]<max_val){
      max_val=front_ind[ino];
      max_ind = ino;
    }
  }

/*kb!!!*/

    if(front_ind[ino]>0){

      mmr_node_coor(Mesh_id,ino,coor);
      if(coor[0]<x_min){
	x_min=coor[0];
	max_ind = ino;
	
      }
    }

  }
  if(max_ind==0){
    max_ind = 1;
  }

/*kbw
#ifdef DEBUG
    printf("choose front on exit: max_val %d, max_ind %d\n",
           max_val,max_ind);
#endif
/*kew*/

  return(max_ind);
}


/*---------------------------------------------------------
  upd_front - to update the front after adding an element to
	the subdoamin
---------------------------------------------------------*/
void upd_front(
  int Mesh_id,	     /* in: mesh ID */
  int Node,
  int Sub_id,
  int Gen_lev
)
{

  int m, iell, ielaux, iel; 

/***************** beginning of the great loop *****************/

/*kbw
#ifdef DEBUG
  printf("\nupdating front for node %d (mesh %d, subdoamin %d, gen_lev %d)\n",
	 Node, Mesh_id, Sub_id, Gen_lev);
#endif
/*kew*/

  /* for each element sharing the node */
  m=0;
  for(iell=0;iell<MAXELNO&&(ielaux=elems[Node][iell])!=0;iell++){
    
    iel=ielaux;

/*kbw
#ifdef DEBUG
  printf("in element %d (active %d) - indyk %d\n",
	 iel, ielaux, sub_ind_el[iel]);
#endif
/*kew*/

    
    if(sub_ind_el[iel]==0){
      /* if there is at least one element sharing the node
	 and not yet considered then add node to the front and break*/
      if(front_ind[Node]==0){
	front_ind[Node]=Sub_id;
	nr_node_front++;

/*kbw
#ifdef DEBUG
  printf("adding node %d (indyk %d) to front - length %d\n",
	 Node,front_ind[Node],nr_node_front);
#endif
/*kew*/

      }

/*kbw
#ifdef DEBUG
      else{
	printf("node %d (indyk %d) already inside front - length %d\n",
	       Node,front_ind[Node],nr_node_front);
      }
#endif
/*kew*/

      m=1;
      break;
    }
  }
  
  /* if all elements belonging already to the subdomain delete the node
     from the front */
  if(m==0){
    if(front_ind[Node]>0){
      front_ind[Node] *= -1;
      nr_node_front--;
/*kbw
#ifdef DEBUG
  printf("deleting node %d (indyk %d) from front - length %d\n",
	 Node,front_ind[Node],nr_node_front);
#endif
/*kew*/
    }
/*kbw
#ifdef DEBUG
      else{
	printf("node %d (indyk %d) already outside front - length %d\n",
	       Node,front_ind[Node],nr_node_front);
      }
#endif
/*kew*/
  }
   
  return;
}

/*---------------------------------------------------------
put_sons - to put sons of element Nel together with their nodes
	on the lists of subdomain elements and nodes
---------------------------------------------------------*/
void put_sons(
  int Mesh_id,	        /* in: mesh ID */
  int El,               /* in: element ID */
  int Sub_id,           /* in: subdomain ID */
  int *Sub_size_p,      /* in/out: subdomain size to be updated */
  int *Sub_size_tot_p   /* in/out: subdomain size to be updated */
  ) 
{

  int elsons[MMC_MAXELSONS+1];
  int ison, iel;

/***************** beginning of the great loop *****************/

  mmr_el_fam(Mesh_id,El,elsons,NULL);
  for(ison=1;ison<=elsons[0];ison++){
    iel=elsons[ison];

/*kbw
#ifdef DEBUG
      printf("considering element %d, mesh_gen %d - status %d\n",
	     iel, mmr_el_gen(Mesh_id, iel), sub_ind_el[iel]);
#endif
/*kew*/

    /* if element not yet considered */
    if(sub_ind_el[iel]==0) {
      
      /* put it on the list of subdomain elements */
      sub_ind_el[iel]=Sub_id;
      (*Sub_size_tot_p)++;
      if(mmr_el_status(Mesh_id,iel)==MMC_ACTIVE) (*Sub_size_p)++;
      
/*kbw
#ifdef DEBUG
      printf("element %d added - status %d, size %d, size_tot %d\n",
	     iel, sub_ind_el[iel], (*Sub_size_p), (*Sub_size_tot_p));
#endif
/*kew*/

      /* for inactive elements put all its antecedents and their nodes */
      if(mmr_el_status(Mesh_id,iel)==MMC_INACTIVE) {
	put_sons(Mesh_id, iel, Sub_id, Sub_size_p, Sub_size_tot_p);
      }
    }
  }
}

/************************************************************
put_list - to put Num on the list List
with length Ll (filled with numbers and zeros at the end)
*************************************************************/
int put_list(int Num, int* List, int Ll)
/*
in: see above
out:
	* - >0 - position already occupied on the list
             0 - put on the list
            -1 - list full, not found on the list
*/
{

  int i, il;
  
  for(i=0;i<Ll;i++){
    if((il=List[i])==0) break;
    /* found on the list on (i+1) position */
    if(Num==il) return(i+1);
  }
  /* if list is full return error message */
  if(i==Ll) return(-1);
  /* update the list and return*/
  List[i]=Num;
  return(0);
}


