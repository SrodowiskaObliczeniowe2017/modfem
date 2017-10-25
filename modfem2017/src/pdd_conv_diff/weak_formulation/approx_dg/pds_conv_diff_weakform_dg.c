/************************************************************************
File pds_conv_diff_weakform_dg.c - weakform functions for conv_diff and DG approximation

Contains definitions of routines:

  pdr_conv_diff_comp_el_stiff_mat_app - to construct UNCONSTRAINED stiffness matrix and
                          a load vector for an element (DG approximation specific version)

  pdr_conv_diff_comp_fa_stiff_mat_app - to construct UNCONSTRAINED stiffness matrix and
                          a load vector for a face (DG approximation specific version)

------------------------------
History:
	2012    - Krzysztof Banas (pobanas@cyf-kr.edu.pl)
*************************************************************************/


#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

/* USED DATA STRUCTURES AND INTERFACES FROM OTHER MODULES */

/* utilities - including simple time measurement library */
#include "uth_intf.h"

/* interface for all mesh manipulation modules */
#include "mmh_intf.h"	

/* interface for all approximation modules */
#include "aph_intf.h"	


/* USED AND IMPLEMENTED PROBLEM DEPENDENT DATA STRUCTURES AND INTERFACES */

/* problem dependent module interface */
#include "pdh_intf.h"	

/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include "pdh_control_intf.h"

/* header files for the problem dependent module for Conv_Diff equations */
#include "../../include/pdh_conv_diff.h"	


/*------------------------------------------------------------
  pdr_conv_diff_comp_el_stiff_mat_app - to construct a stiffness matrix and  
                          a load vector for an element
------------------------------------------------------------*/
int pdr_conv_diff_comp_el_stiff_mat_app(/*returns: >=0 -success code, <0 -error code*/
  int Problem_id,     /* in: approximation field ID  */
  int El_id,        /* in: unique identifier of the element */ 
  int Comp_sm,       /* in: indicator for the scope of computations: */
                     /*   PDC_NO_COMP  - do not compute anything */
                     /*   PDC_COMP_SM - compute entries to stiff matrix only */
                     /*   PDC_COMP_RHS - compute entries to rhs vector only */
                     /*   PDC_COMP_BOTH - compute entries for sm and rhsv */
  int Pdeg_in,        /* in: enforced degree of polynomial (if > 0 ) */
  int* Nr_dof_ent,   /* in: size of arrays, */
                /* out: no of filled entries, i.e. number of mesh entities*/
                /* with which dofs and stiffness matrix blocks are associated */
  int* List_dof_ent_type, /* out: list of no of dofs for 'dof' entity */
  int* List_dof_ent_id,   /* out: list of no of dofs for 'dof' entity */
  int* List_dof_ent_nrdof,/* out: list of no of dofs for 'dof' entity */
  int* Nrdofs_loc,        /* in(optional): size of Stiff_mat and Rhs_vect */
                /* out(optional): actual number of dofs per integration entity*/
  double* Stiff_mat,      /* out(optional): stiffness matrix stored columnwise*/
  double* Rhs_vect,       /* out(optional): rhs vector */
  char* Rewr_dofs        /* out(optional): flag to rewrite or sum up entries */
                         /*   'T' - true, rewrite entries when assembling */
                         /*   'F' - false, sum up entries when assembling */
  )
{

  int pdeg;		/* degree of polynomial */
  int el_mate;		/* element material */
  int sol_vec_id;       /* indicator for the solution dofs */
  int num_shap;         /* number of element shape functions */
  double sol_dofs[APC_MAXELSD]; /* solution dofs */
  double implicit;

/*   static int ndofs_max=0;            /\* local dimension of the problem *\/ */
/* #pragma omp threadprivate(ndofs_max) */

  // to make OpenMP working
  int ndofs_max=0;

  int field_id, mesh_id;  
  int kk, idofs, jdofs;
  int i,j,k;
  int max_dof_ent, max_nrdofs;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper mesh */
  pdv_conv_diff_current_problem_id = Problem_id; 
  field_id = Problem_id;
  mesh_id = apr_get_mesh_id(field_id);

/* get material number */
  el_mate =  mmr_el_groupID(mesh_id, El_id);

  /* find degree of polynomial and number of element scalar dofs */
  if(Pdeg_in>0) pdeg = Pdeg_in;
  else{
#ifdef DEBUG
    /* check input for DG approximation */
    if(mmr_el_status(mesh_id,El_id)!=MMC_ACTIVE){
      printf("Asking for pdeg of inactive element %d in comp_el_stiff_mat !\n",
	     El_id);
      exit(-1);
    }
#endif
      pdeg = apr_get_ent_pdeg(field_id, APC_ELEMENT, El_id);
  }
  num_shap = apr_get_el_pdeg_numshap(field_id, El_id, &pdeg);

  max_dof_ent = *Nr_dof_ent;
  *Nr_dof_ent = 1;

#ifdef DEBUG
  if(max_dof_ent < 1){
    printf("Too small arrays List_dof_... passed to comp_el_stiff_mat\n");
    printf("%d < %d. Exiting !!!",max_dof_ent, *Nr_dof_ent);
    exit(-1);
  }
#endif

  List_dof_ent_type[0] = PDC_ELEMENT;
  List_dof_ent_id[0] = El_id;
  /* scalar problem nrdofs = numshap */
  List_dof_ent_nrdof[0] = num_shap;
  int num_dofs = num_shap;

/*kbw
    printf("In pdr_com_el_stiff_mat: field_id %d, El_id %d, Comp_sm %d, Nr_dof_ent %d, pdeg %d\n",
	   field_id, El_id, Comp_sm, *Nr_dof_ent, pdeg);
    printf("For each block: \ttype, \tid, \tnrdof\n");
    for(i=0;i<*Nr_dof_ent;i++){
      printf("\t\t\t%d\t%d\t%d\n",
	     List_dof_ent_type[i],List_dof_ent_id[i],List_dof_ent_nrdof[i]);
    }
    getchar();getchar();
/*kew*/

  if(Comp_sm!=PDC_NO_COMP){

#ifdef DEBUG
    if(Nrdofs_loc == NULL || Stiff_mat == NULL || Rhs_vect == NULL){
printf("NULL arrays Stiff_mat and Rhs_vect in pdr_comp_stiff_el_mat. Exiting!");
      exit(-1);
    }
    if(*Nrdofs_loc<num_shap){
printf("Too small arrays Stiff_mat and Rhs_vect passed to comp_el_stiff_mat\n");
      printf("%d < %d. Exiting !!!", *Nrdofs_loc, num_shap);
      exit(-1);
    }

#endif


    *Nrdofs_loc = num_shap;

    /* get the most recent solution degrees of freedom */
    if(mmr_el_status(mesh_id,El_id)==MMC_ACTIVE){

      /* select solution vector to be used */
      sol_vec_id = 1;
      apr_read_ent_dofs(field_id,APC_ELEMENT,El_id,num_shap,sol_vec_id,sol_dofs);
    }
    else{

      /*!!! coarse element dofs should be supplied by calling routine !!!*/
      for(i=0;i<num_shap;i++) sol_dofs[i]=0;

    }


    if(Comp_sm==PDC_COMP_SM||Comp_sm==PDC_COMP_BOTH){

      /* initialize the matrices to zero */
      for(i=0;i<num_dofs*num_dofs;i++) Stiff_mat[i]=0.0;

    }

    if(Comp_sm==PDC_COMP_RHS||Comp_sm==PDC_COMP_BOTH){
    
      /* initialize the vector to zero */
      for(i=0;i<num_dofs;i++) Rhs_vect[i]=0.0; 

    }

    int diagonal[5]={0,0,0,0,0}; // diagonality of: M, A_ij, B_j, T_i and C  
    // coefficient matrices returned to apr_num_int_el by pdr_el_coeff
    /* perform numerical integration of terms from the weak formualation */
    apr_num_int_el(Problem_id,field_id,El_id,Comp_sm,&pdeg,sol_dofs,sol_dofs,
		   diagonal, Stiff_mat, Rhs_vect);


/*kbw
    if(El_id>0){
      printf("Element %d: Stiffness matrix:\n",El_id);
      for (idofs=0;idofs<num_shap*PDC_CONV_DIFF_NREQ;idofs++) {
	for (jdofs=0;jdofs<num_shap*PDC_CONV_DIFF_NREQ;jdofs++) {
	  printf("%20.12lf",Stiff_mat[idofs+jdofs*num_shap*PDC_CONV_DIFF_NREQ]);        
	} 
	printf("\n");
      }
      printf("Element %d: Rhs_vect:\n",El_id);
      for (idofs=0;idofs<num_shap*PDC_CONV_DIFF_NREQ;idofs++) {
	printf("%20.12lf",Rhs_vect[idofs]);
      }
      printf("\n");
      //getchar(); 
    } 
/*kew*/

    if(Rewr_dofs != NULL) *Rewr_dofs = 'F';

  } /* end if computing SM and/or RHSV */

  return(1);
}

/*------------------------------------------------------------
  pdr_conv_diff_comp_fa_stiff_mat_app - to construct a stiffness matrix and  
                          a load vector for a face - special routine for discontinuous 
                          Galerkin approximation fields
------------------------------------------------------------*/
int pdr_conv_diff_comp_fa_stiff_mat_app(/*returns: >=0 -success code, <0 -error code*/
  int Problem_id,     /* in: approximation field ID  */
  int Fa_id,        /* in: unique identifier of the face */ 
  int Comp_sm,       /* in: indicator for the scope of computations: */
                     /*   PDC_NO_COMP  - do not compute anything */
                     /*   PDC_COMP_SM - compute entries to stiff matrix only */
                     /*   PDC_COMP_RHS - compute entries to rhs vector only */
                     /*   PDC_COMP_BOTH - compute entries for sm and rhsv */
  int Pdeg_in,     /* in: enforced degree of polynomial (if > 0 ) */
  int* Nr_dof_ent,   /* in: size of arrays List_dof_ent_... */
                     /* out: number of mesh entities with which dofs and */
                     /*      stiffness matrix blocks are associated */
  int* List_dof_ent_type, /* out: list of no of dofs for 'dof' entity */
  int* List_dof_ent_id,   /* out: list of no of dofs for 'dof' entity */
  int* List_dof_ent_nrdof,/* out: list of no of dofs for 'dof' entity */
  int* Nrdofs_loc,        /* in(optional): size of Stiff_mat and Rhs_vect */
                /* out(optional): actual number of dofs per integration entity*/
  double* Stiff_mat,      /* out(optional): stiffness matrix stored columnwise*/
  double* Rhs_vect,       /* out(optional): rhs vector */
  char* Rewr_dofs        /* out(optional): flag to rewrite or sum up entries */
                         /*   'T' - true, rewrite entries when assembling */
                         /*   'F' - false, sum up entries when assembling */
  )
{

/* local variables */
  int fa_type;		/* type of face: quad or triangle */
  int fa_bc;	        /* bc flag for a face */
  int bc_type;	        /* type of bc*/
  int base_1, base_2;	/* type of basis functions */
  int nreq;		/* number of equations */
  int pdeg1=0, num_shap1=0;/* local dimension for neigbhor 1 */
  int pdeg2=0, num_shap2=0;/* local dimension for neigbhor 2 */
  int el_mate1=0;	/* material number for neigbhor 1 */
  int el_mate2=0;	/* material number for neigbhor 2 */
  int ngauss;		/* number of gauss points in 1 dimension */
  int face_neig[2];	/* list of neighbors */
  int neig_sides[2]={0,0};/* sides of face wrt neigs */
  int node_shift;
  double acoeff[4], bcoeff[2];/* to transform coordinates between faces */
  double loc_xg[6];	/* local coord of gauss points for neighbors */
  double xg[450];   	/* gauss points in 2D*/
  double wg[150];      	/* gauss weights in 2D*/
  double determ; 	/* determinant of jacobi matrix */
  double hsize1, hsize2; /* size of elements */
  double hf_size;        /* size of face */
  double area; 		/* area for integration rule */
  double time;          /* current time instant */
  double xcoor[3],xcoor2[3];/* global coord of gauss point */
  double vec_norm[3];	/* normal vector to the considered face */
  double u_val1[PDC_CONV_DIFF_NREQ]; /* computed solution */
  double u_x1[PDC_CONV_DIFF_NREQ];   /* gradient of computed solution */
  double u_y1[PDC_CONV_DIFF_NREQ];   /* gradient of computed solution */
  double u_z1[PDC_CONV_DIFF_NREQ];   /* gradient of computed solution */
  double u_val2[PDC_CONV_DIFF_NREQ]; /* computed solution */
  double u_x2[PDC_CONV_DIFF_NREQ];   /* gradient of computed solution */
  double u_y2[PDC_CONV_DIFF_NREQ];   /* gradient of computed solution */
  double u_z2[PDC_CONV_DIFF_NREQ];   /* gradient of computed solution */
  double base_phi1[APC_MAXELVD];    /* basis functions */
  double base_dphix1[APC_MAXELVD];  /* x-derivatives of basis function */
  double base_dphiy1[APC_MAXELVD];  /* y-derivatives of basis function */
  double base_dphiz1[APC_MAXELVD];  /* y-derivatives of basis function */
  double base_phi2[APC_MAXELVD];    /* basis functions */
  double base_dphix2[APC_MAXELVD];  /* x-derivatives of basis function */
  double base_dphiy2[APC_MAXELVD];  /* y-derivatives of basis function */
  double base_dphiz2[APC_MAXELVD];  /* y-derivatives of basis function */
  int el_nodes1[MMC_MAXELVNO+1];        /* list of nodes of El */
  int el_nodes2[MMC_MAXELVNO+1];        /* list of nodes of El */
  double node_coor1[3*MMC_MAXELVNO];  /* coord of nodes of El */
  double node_coor2[3*MMC_MAXELVNO];  /* coord of nodes of El */
  double dofs_loc1[APC_MAXELSD]; /* element solution dofs */
  double dofs_loc2[APC_MAXELSD]; /* element solution dofs */
  double anx1[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];  /* coeff of PDE */
  double any1[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];  /* coeff of PDE */
  double anz1[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];  /* coeff of PDE */
  double bn1[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];   /* coeff of PDE */
  double fval1[PDC_CONV_DIFF_NREQ]; /* coeff of PDE */
  double gval1[PDC_CONV_DIFF_NREQ]; /* coeff of PDE */
  double qn1[PDC_CONV_DIFF_NREQ];      /* source for the first neigbor */
  double anx2[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];  /* coeff of PDE */
  double any2[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];  /* coeff of PDE */
  double anz2[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];  /* coeff of PDE */
  double bn2[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];   /* coeff of PDE */
  double fval2[PDC_CONV_DIFF_NREQ]; /* coeff of PDE */
  double gval2[PDC_CONV_DIFF_NREQ]; /* coeff of PDE */
  double qn2[PDC_CONV_DIFF_NREQ];      /* source for the second neigbor */
  double norm_vel[PDC_CONV_DIFF_NREQ]; /* normal velocity for each component of u */
  double norm_vel1[PDC_CONV_DIFF_NREQ]; /* normal velocity for each component of u */
  double norm_vel2[PDC_CONV_DIFF_NREQ]; /* normal velocity for each component of u */

  int i, ki, kk, ieq1, ieq2, time_dis;
  int iaux, jaux;
  int idofs, jdofs, ndofs_fa, ndofs1, ndofs2, neig1_id, neig2_id;
  double daux, faux, eaux, gaux, haux, penalty, implicit;


  int field_id, mesh_id;  
  int sol_vec_id;       /* indicator for the solution dofs */


/*++++++++++++++++ executable statements ++++++++++++++++*/

  time=pdr_time_d_params(PDC_USE_CURRENT_PROBLEM_ID, 4);
  //int problem_name = pdr_ctrl_i_params(PDC_USE_CURRENT_PROBLEM_ID, 1); changed for string
  char problem_name[300];
  pdr_problem_name(PDC_USE_CURRENT_PROBLEM_ID, problem_name);

#ifdef DEBUG
  if(*Nr_dof_ent<2){ // two elements = dof entities per face
    printf("Too small arrays List_dof_... passed to comp_fa_stiff_mat\n");
    printf("%d < 2. Exiting !!!", *Nr_dof_ent);
    exit(-1);
  }
#endif

  /* select the proper mesh */
  pdv_conv_diff_current_problem_id = Problem_id; 
  field_id = Problem_id;
  mesh_id = apr_get_mesh_id(field_id);

  /* get type of face and bc flag */
  fa_type=mmr_fa_type(mesh_id,Fa_id);
  fa_bc=mmr_fa_bc(mesh_id,Fa_id);
  mmr_fa_area(mesh_id,Fa_id,&daux,NULL);
  hf_size = sqrt(daux);

  bc_type = pdr_get_bc_type(fa_bc);

  /* get approximation parameters */
  nreq=PDC_CONV_DIFF_NREQ;
  if(nreq!=pdr_ctrl_i_params(PDC_USE_CURRENT_PROBLEM_ID, 5)){
    printf("nreq(%d)!=PDC_CONV_DIFF_NREQ - performance penalty. Exiting.\n",nreq);
    exit(-1);
  }

/* get neigbors list with corresponding neighbors' sides numbers*/
  mmr_fa_neig(mesh_id,Fa_id,face_neig,neig_sides,&node_shift,
					NULL,acoeff,bcoeff);

  /* get material number */
  neig1_id = abs(face_neig[0]);
  el_mate1 =  mmr_el_groupID(mesh_id, neig1_id);
  hsize1 = mmr_el_hsize(mesh_id,abs(face_neig[0]),NULL,NULL,NULL);
  base_1=apr_get_base_type(field_id, neig1_id);

  if(Pdeg_in>0) pdeg1 = Pdeg_in;
  else{
    /* find degree of polynomial and number of element scalar dofs */
    /*!!! need some trick for coarse solve with different Pdeg !!!*/
#ifdef DEBUG
    /* check input for DG approximation */
    if(mmr_el_status(mesh_id,neig1_id)!=MMC_ACTIVE){
      printf("Asking for pdeg of inactive element in comp_fa_stiff_mat !\n");
      exit(-1);
    }
#endif
    pdeg1 = apr_get_ent_pdeg(field_id, PDC_ELEMENT, neig1_id);
  }
  num_shap1 = apr_get_el_pdeg_numshap(field_id, neig1_id, &pdeg1);
  ndofs1 = num_shap1*nreq;

  *Nr_dof_ent = 1;
  List_dof_ent_type[0] = PDC_ELEMENT;
  List_dof_ent_id[0] = neig1_id;
  List_dof_ent_nrdof[0] = ndofs1;

  if(Comp_sm!=PDC_NO_COMP){
    mmr_el_node_coor(mesh_id,neig1_id,el_nodes1,node_coor1);

    // for Poisson / Laplace
    //if(problem_name==1){ 
    if(strcmp(problem_name, "DIFF_IN_CUBE") == 0 || strcmp(problem_name, "DIFF_IN_BOX") == 0 ){

      for(i=0;i<ndofs1; i++) dofs_loc1[i]=0.0;

    }
    else{

    /* get the most recent solution degrees of freedom */
    /*!!! coarse element dofs should be supplied by calling routine !!!*/
    sol_vec_id = 1;
    apr_read_ent_dofs(field_id,PDC_ELEMENT,neig1_id,ndofs1,sol_vec_id,
    	      dofs_loc1);

    }
  }


/* if not on the wall, i.e. if(bc_type==INTERNAL) */
  if(face_neig[1]!=0){

    neig2_id = abs(face_neig[1]);

    el_mate2 =  mmr_el_groupID(mesh_id, neig2_id);
    hsize2 = mmr_el_hsize(mesh_id,abs(face_neig[1]),NULL,NULL,NULL);
    base_2=apr_get_base_type(field_id, neig2_id);

    if(Pdeg_in>0) pdeg2 = Pdeg_in;
    else{
  /* find degree of polynomial and number of element scalar dofs */
  /*!!! need some trick for coarse solve with different Pdeg !!!*/
#ifdef DEBUG
      /* check input for DG approximation */
      if(mmr_el_status(mesh_id,neig2_id)!=MMC_ACTIVE){
	printf("Asking for pdeg of inactive element in comp_fa_stiff_mat !\n");
	exit(-1);
      }
#endif
      pdeg2 = apr_get_ent_pdeg(field_id, PDC_ELEMENT, neig2_id);
    }
    num_shap2 = apr_get_el_pdeg_numshap(field_id, neig2_id, &pdeg2);
    ndofs2 = num_shap2*nreq;

    *Nr_dof_ent = 2;
    List_dof_ent_type[1] = PDC_ELEMENT;
    List_dof_ent_id[1] = neig2_id;
    List_dof_ent_nrdof[1] = ndofs2;

    if(Comp_sm!=PDC_NO_COMP){
      mmr_el_node_coor(mesh_id,neig2_id,el_nodes2,node_coor2);
      
      // for Poisson / Laplace
      //if(problem_name==1){ 
      if(strcmp(problem_name, "DIFF_IN_CUBE") == 0 || strcmp(problem_name, "DIFF_IN_BOX") == 0 ){
	
	for(i=0;i<ndofs1; i++) dofs_loc1[i]=0.0;
	
      }
      else{
	
	/* get the most recent solution degrees of freedom */
	/* !!! coarse element dofs should be supplied by calling routine !!!*/
	sol_vec_id = 1;
	apr_read_ent_dofs(field_id,PDC_ELEMENT,neig2_id,ndofs2,sol_vec_id,
			  dofs_loc2);
      }
      
    }
  } /* end if second neighbor not boundary */ 
  else{
    pdeg2 = 0;
    num_shap2 = 0;
    ndofs2 = 0;
  }

  /*kbw 
#ifdef DEBUG
    printf("In pdr_comp_fa_stiff_mat: field_id %d, mesh_id %d, Fa_id %d, Comp_sm %d\n",
	   field_id, mesh_id, Fa_id, Comp_sm);
    printf("For each block: \ttype, \tid, \tnrdof\n");
    for(i=0;i<*Nr_dof_ent;i++){
      printf("\t\t\t%d\t%d\t%d\n",
	     List_dof_ent_type[i],List_dof_ent_id[i],List_dof_ent_nrdof[i]);
    }
#endif
  /*kew*/

  if(Comp_sm!=PDC_NO_COMP){

   /* initialize the matrices to zero */
   ndofs_fa = ndofs1+ndofs2;

#ifdef DEBUG
   if(*Nrdofs_loc < ndofs_fa ){
    printf("Too small arrays Stiff_mat and Rhs_vect passed to comp_el_stiff_mat\n");
    printf("%d < %d. Exiting !!!", *Nrdofs_loc, ndofs_fa);
    exit(-1);
  }
#endif

   *Nrdofs_loc = ndofs_fa;

   for(i=0;i<ndofs_fa*ndofs_fa;i++) Stiff_mat[i]=0.0;
   for(i=0;i<ndofs_fa;i++) Rhs_vect[i]=0.0; 
    
    /* prepare data for gaussian integration */
   if(face_neig[1]!=0){
     apr_set_quadr_2D(fa_type,base_1,&pdeg1, &iaux,NULL,NULL);
     apr_set_quadr_2D(fa_type,base_2,&pdeg2, &jaux,NULL,NULL);
     if(iaux>jaux) apr_set_quadr_2D(fa_type,base_1,&pdeg1,&ngauss,xg,wg); 
     else apr_set_quadr_2D(fa_type,base_2,&pdeg2,&ngauss,xg,wg);
   }
   else apr_set_quadr_2D(fa_type,base_1,&pdeg1,&ngauss,xg,wg);
   
   /* set penalty and implicitness parameters */ 
   if(base_1==APC_BASE_COMPLETE_DG){
     if(face_neig[1]!=0) penalty= utm_min(pdeg1,pdeg2);
     else penalty = pdeg1;
   }
   else{
     if(face_neig[1]!=0){
       daux=utm_min(pdeg1/100,pdeg1%100);
       eaux=utm_min(pdeg2/100,pdeg2%100);
       penalty=utm_min(daux,eaux);
     }
     else penalty=utm_min(pdeg1/100,pdeg1%100);
   }
   
/*kb!!!!!!!!!!!!!!*/
   penalty *= pdr_ctrl_d_params(PDC_USE_CURRENT_PROBLEM_ID, 12);
/*kb!!!!!!!!!!!!!!*/

   implicit = pdr_time_d_params(PDC_USE_CURRENT_PROBLEM_ID, 2);

   

    /*kbw 
    if(Fa_id>0){
      printf("in stiff_face: face %d, type %d, bc_type %d, NREQ %d, base-type %d, node_shift %d\n",
	     Fa_id,fa_type, bc_type,nreq,base_1, node_shift);
      printf("neigh1: %d, side %d, pdeg %d, num_shap %d\n",
	     face_neig[0],neig_sides[0],pdeg1,num_shap1);
      for(i=0;i<el_nodes1[0];i++){
	printf("node %d (global - %d): x - %f, y - %f, z - %f\n",
	       i,el_nodes1[i+1],node_coor1[3*i],
	       node_coor1[3*i+1],node_coor1[3*i+2]);
      }
      if(face_neig[1]!=0){
	printf("neigh2: %d, side %d, pdeg %d, num_shap %d\n",
	       face_neig[1],neig_sides[1],pdeg2,num_shap2);
	for(i=0;i<el_nodes2[0];i++){
	  printf("node %d (global - %d): x - %f, y - %f, z - %f\n",
		 i,el_nodes2[i+1],node_coor2[3*i],
		 node_coor2[3*i+1],node_coor2[3*i+2]);
	}
      }
      printf("transformation coefficients: a %lf, b %lf\n",acoeff,bcoeff);
      printf("Gauss points:\n");
      for (ki=0;ki<ngauss;ki++) {
	printf("Xg: %lf, %lf, Weight: %lf\n", xg[2*ki], xg[2*ki+1], wg[ki]);
      }
      getchar();
    }
  /*kew*/

 
/* loop over integration points */
   for (ki=0;ki<ngauss;ki++) {

/* find coordinates within neighboring elements for a point on face */
    mmr_fa_elem_coor(mesh_id,&xg[2*ki],face_neig,neig_sides,node_shift,
			acoeff,bcoeff,loc_xg);

/* at the gauss point for neig 1, compute basis functions */
    iaux = 3+neig_sides[0]; /* boundary data for face neig_sides[0] */

    determ = apr_elem_calc_3D(iaux, nreq, &pdeg1, base_1,
			      loc_xg,node_coor1,dofs_loc1,
			      base_phi1,base_dphix1,base_dphiy1,base_dphiz1,
			      xcoor,u_val1,u_x1,u_y1,u_z1,vec_norm);

/* get coefficients of convection-diffusion-reaction equations */
    pdr_fa_coeff(Fa_id,bc_type,neig1_id,el_mate1,
		   hsize1,time,pdeg1,
		   xcoor,vec_norm,u_val1,u_x1,u_y1,u_z1,
		   anx1,any1,anz1,bn1,fval1,gval1,qn1,norm_vel1);

/* make the normal velocity returned by the first neighbor default */
    for(i=0;i<PDC_CONV_DIFF_NREQ;i++) norm_vel[i]=norm_vel1[i];

    if(face_neig[1]!=0) {

/* if first neighbor is big, we need proper determinant of Jacobian */
      if(face_neig[0]<0){
        iaux = 3+neig_sides[1]; /* boundary data for face neig_sides[1] */

        determ = apr_elem_calc_3D(iaux, nreq, &pdeg2, base_2,
				  &loc_xg[3],node_coor2,dofs_loc2,
				  base_phi2,base_dphix2,base_dphiy2,base_dphiz2,
				  xcoor2,u_val2,u_x2,u_y2,u_z2,NULL);
      }
      else{
        iaux = 2; /* no boundary data */
        apr_elem_calc_3D( iaux, nreq, &pdeg2, base_2,
			   &loc_xg[3],node_coor2,dofs_loc2,
			   base_phi2,base_dphix2,base_dphiy2,base_dphiz2,
			   xcoor2,u_val2,u_x2,u_y2,u_z2,NULL);
      }

/* get coefficients from the second neighbor */
      pdr_fa_coeff(Fa_id,bc_type,abs(face_neig[1]),
		    el_mate2,hsize2,time,pdeg2,
		    xcoor,vec_norm,u_val2,u_x2,u_y2,u_z2,
		    anx2,any2,anz2,bn2,fval2,gval2,qn2,norm_vel2);


#ifdef DEBUG
      for(i=0;i<PDC_CONV_DIFF_NREQ;i++){
	if(bc_type==PDC_BC_MIXED && norm_vel1[i] * norm_vel2[i] < -1e-6){
	printf("opposite signs for normal velocities on face %d:\n",Fa_id);
	printf("%lf != %lf\n",norm_vel1[i],norm_vel2[i]);
	}
      }
#endif

      for(i=0;i<PDC_CONV_DIFF_NREQ;i++) norm_vel[i] = 0.5 * ( norm_vel1[i] + norm_vel2[i] );

#ifdef DEBUG
      for(i=0;i<3;i++){
        if(fabs(xcoor2[i]-xcoor[i])>1e-9){
          printf("%lf != %lf\n",xcoor[i],xcoor2[i]);
          printf("Coord in neigh1 != coord in neigh2 for face %d!\n",Fa_id);
          return(-1);
        }
      }
#endif

    } /* end if exist second neighbor */

/* coefficient for 2D numerical integration */
    area = determ*wg[ki];

    /*kbw
if(Fa_id>0){
printf("at integration point %d, coor: local %lf, %lf, weight %lf\n",
ki,xg[2*ki],xg[2*ki+1],wg[ki]);
printf("coor global %lf, %lf, %lf, determ %lf, area %lf\n",
xcoor[0],xcoor[1],xcoor[2],determ,area);
printf("normal %lf, %lf, %lf\n", 
	vec_norm[0], vec_norm[1], vec_norm[2]);
printf("NORMAL velocities: \neig1= %lf, neig2= %lf, average= %lf\n",
	 norm_vel1[0],norm_vel2[0],norm_vel[0]);
printf("coord for neigh1: %lf, %lf, %lf, hsize %lf\n",
loc_xg[0],loc_xg[1],loc_xg[2],hsize1);
if(face_neig[1]!=0)
  printf("coord for neigh2: %lf, %lf, %lf, hsize %lf\n",
loc_xg[3],loc_xg[4],loc_xg[5],hsize2);
printf("weight %lf, integration coeff. %lf\n", wg[ki], area);
printf("%d shape functions and derivatives for neigh1: \n",num_shap1);
for(i=0;i<num_shap1;i++){
  printf("fun - %lf, der: x - %lf, y - %lf, z - %lf\n",
  base_phi1[i],base_dphix1[i],base_dphiy1[i],base_dphiz1[i]);
}
if(face_neig[1]!=0){
 printf("%d shape functions and derivatives for neigh2: \n",num_shap2);
 for(i=0;i<num_shap2;i++){
  printf("fun - %lf, der: x - %lf, y - %lf, z - %lf\n",
  base_phi2[i],base_dphix2[i],base_dphiy2[i],base_dphiz2[i]);
 }
}
  printf("solution %lf, der %lf, %lf, %lf\n", 
	 u_val1[0],u_x1[0],u_y1[0],u_z1[0]);
  printf("diffusion coeff : %lf %lf %lf\n",anx1[0],any1[0],anz1[0]);
  printf("convection coeff: %lf\n",bn1[0]);
  printf("rhs coeff       : %lf %lf %lf\n",fval1[0],gval1[0],qn1[0]);
  if(face_neig[1]!=0){
  printf("solution for neigh2: %lf, der %lf, %lf, %lf\n", 
	 u_val2[0],u_x2[0],u_y2[0],u_z2[0]);
    printf("diffusion coeff : %lf %lf %lf\n",anx2[0],any2[0],anz2[0]);
    printf("convection coeff: %lf\n",bn2[0]);
    printf("rhs coeff       : %lf %lf %lf\n",fval2[0],gval2[0],qn2[0]);
  }
  printf("penalty %lf : fval1[0] %lf, hf_size %lf\n",
       penalty,fval1[0],hf_size);
 
getchar();
}
/*kew*/ 


    if(bc_type==PDC_INTERIOR){

      for(ieq2=0;ieq2<PDC_CONV_DIFF_NREQ;ieq2++){
        for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){

          kk=(ndofs_fa*ieq2+ieq1)*num_shap1;

          for (jdofs=0;jdofs<num_shap1;jdofs++) {
            for (idofs=0;idofs<num_shap1;idofs++) {

/* stiffness matrices */
              Stiff_mat[kk+idofs] += (
			        (
	 base_phi1[jdofs] * (
            anx1[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphix1[idofs] +
            any1[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphiy1[idofs] +
            anz1[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphiz1[idofs]
                          )	  
	-base_phi1[idofs] * (
            anx1[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphix1[jdofs] +
            any1[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphiy1[jdofs] +
            anz1[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphiz1[jdofs]
                          ) 
				) *  0.5 * implicit
/*penalty*/
	 + penalty*base_phi1[jdofs]*base_phi1[idofs]/hf_size
/*penalty*/
				       ) * area;

            }/* idofs */

            kk+=ndofs_fa;

          } /* jdofs */

          kk=ndofs_fa*(ndofs1+ieq2*num_shap2)+ieq1*num_shap1;

          for (jdofs=0;jdofs<num_shap2;jdofs++) {
            for (idofs=0;idofs<num_shap1;idofs++) {

              Stiff_mat[kk+idofs] += (
			        (
	-base_phi2[jdofs] * (
            anx1[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphix1[idofs] +
            any1[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphiy1[idofs] +
            anz1[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphiz1[idofs]
                          )	  
	-base_phi1[idofs] * (
            anx2[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphix2[jdofs] +
            any2[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphiy2[jdofs] +
            anz2[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphiz2[jdofs]
                          ) 
				) *  0.5 * implicit
/*penalty*/
	 - penalty*base_phi2[jdofs]*base_phi1[idofs]/hf_size
/*penalty*/
				       ) * area;

            }/* idofs */

            kk+=ndofs_fa;

          } /* jdofs */

          kk=ndofs_fa*ieq2*num_shap1+ndofs1+ieq1*num_shap2;

          for (jdofs=0;jdofs<num_shap1;jdofs++) {
            for (idofs=0;idofs<num_shap2;idofs++) {

              Stiff_mat[kk+idofs] += (
			        (
	 base_phi1[jdofs] * (
            anx2[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphix2[idofs] +
            any2[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphiy2[idofs] +
            anz2[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphiz2[idofs]
                          )	  
	 +base_phi2[idofs] * (
            anx1[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphix1[jdofs] +
            any1[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphiy1[jdofs] +
            anz1[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphiz1[jdofs]
                          ) 
				) *  0.5 * implicit 
/*penalty*/
	 - penalty*base_phi1[jdofs]*base_phi2[idofs]/hf_size
/*penalty*/
				       ) * area;


            }/* idofs */

            kk+=ndofs_fa;

          } /* jdofs */

          kk=ndofs_fa*(ndofs1+ieq2*num_shap2)+ ndofs1+ieq1*num_shap2;

          for (jdofs=0;jdofs<num_shap2;jdofs++) {
            for (idofs=0;idofs<num_shap2;idofs++) {

              Stiff_mat[kk+idofs] += (
			        (
	-base_phi2[jdofs] * (
            anx2[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphix2[idofs] +
            any2[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphiy2[idofs] +
            anz2[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphiz2[idofs]
                          )	  
	+base_phi2[idofs] * (
            anx2[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphix2[jdofs] +
            any2[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphiy2[jdofs] +
            anz2[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * base_dphiz2[jdofs]
                          ) 
				) *  0.5* implicit
/*penalty*/
	 + penalty*base_phi2[jdofs]*base_phi2[idofs]/hf_size
/*penalty*/
				       ) * area;

            }/* idofs */

            kk+=ndofs_fa;

          } /* jdofs */

        }/* ieq1 */
      } /* ieq2 */
      for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){

/* upwind terms - norm_vel[ieq1] gives normal velocity  !!! */
	if(norm_vel[ieq1]>0) {

          kk=(ndofs_fa*ieq1+ieq1)*num_shap1;
	  
	  for (jdofs=0;jdofs<num_shap1;jdofs++) {
	    for (idofs=0;idofs<num_shap1;idofs++) {
	      
/* inflow-outflow condition - we assume b_i*n_i is diagonal !!! */
	      Stiff_mat[kk+idofs] += implicit * area *
	bn1[ieq1*PDC_CONV_DIFF_NREQ+ieq1] * base_phi1[idofs] * base_phi1[jdofs];

	    }/* idofs */

	    kk+=ndofs_fa;

	  } /* jdofs */

          kk=ndofs_fa*ieq1*num_shap1+ndofs1+ieq1*num_shap2;
	  
	  for (jdofs=0;jdofs<num_shap1;jdofs++) {
	    for (idofs=0;idofs<num_shap2;idofs++) {

	      Stiff_mat[kk+idofs] -= implicit * area *
	bn2[ieq1*PDC_CONV_DIFF_NREQ+ieq1] * base_phi2[idofs] * base_phi1[jdofs];

	    }/* idofs */

	    kk+=ndofs_fa;

	  } /* jdofs */
	}
	else if(norm_vel[ieq1]<0){

          kk=ndofs_fa*(ndofs1+ieq1*num_shap2)+ieq1*num_shap1;

	  for (jdofs=0;jdofs<num_shap2;jdofs++) {
	    for (idofs=0;idofs<num_shap1;idofs++) {
		
/* inflow-outflow condition - we assume b_i*n_i is diagonal !!! */
	      Stiff_mat[kk+idofs] += implicit * area *
	bn1[ieq1*PDC_CONV_DIFF_NREQ+ieq1] * base_phi1[idofs] * base_phi2[jdofs];

	    }/* idofs */
	    
	    kk+=ndofs_fa;

	  } /* jdofs */

          kk=ndofs_fa*(ndofs1+ieq1*num_shap2)+ ndofs1+ieq1*num_shap2;

	  for (jdofs=0;jdofs<num_shap2;jdofs++) {
	    for (idofs=0;idofs<num_shap2;idofs++) {
	      Stiff_mat[kk+idofs] -= implicit * area *
	bn2[ieq1*PDC_CONV_DIFF_NREQ+ieq1] * base_phi2[idofs] * base_phi2[jdofs];

	    }/* idofs */
	    
	    kk+=ndofs_fa;

	  } /* jdofs */

	} /* if neig2 upwind */

      } /* ieq1 */

/* rhs due to explicit fluxes */
      for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){

	if(norm_vel[ieq1]>0) {
	  haux = -qn1[ieq1];
	}
	else {
	  haux = -qn2[ieq1];
	}
	    
	kk = ieq1*num_shap1;

	for (idofs=0;idofs<num_shap1;idofs++) {
	  
	  /* upwind terms - jump in test functions, upwinded fluxes */
	  Rhs_vect[kk+idofs] +=
	    haux * base_phi1[idofs] * area;

	}/* idofs */

	kk = PDC_CONV_DIFF_NREQ*num_shap1+ieq1*num_shap2;

	for (idofs=0;idofs<num_shap2;idofs++) {
	  
	  Rhs_vect[kk+idofs] += 
	    - haux * base_phi2[idofs] * area;

	}/* idofs */

      }/* ieq1 */

/* right hand side correction due to an and bn for time dependent problems */
      if(implicit<1.0){

	for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){

	  daux=0.0;faux=0.0;eaux=0.0;gaux=0.0;
	  for(ieq2=0;ieq2<PDC_CONV_DIFF_NREQ;ieq2++){
	    daux +=
           	  anx1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_val1[ieq2]
           	- anx1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_val2[ieq2];
	    faux +=
           	  any1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_val1[ieq2]
           	- any1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_val2[ieq2];
	    gaux +=
           	  anz1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_val1[ieq2]
           	- anz1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_val2[ieq2];
	    eaux += 
           	- anx1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_x1[ieq2]
           	- any1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_y1[ieq2]
           	- anz1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_z1[ieq2];
	  }     

	  kk = ieq1*num_shap1;

	  for (idofs=0;idofs<num_shap1;idofs++) {

	    Rhs_vect[kk+idofs] += (
		daux*base_dphix1[idofs]  +
		faux*base_dphiy1[idofs]  +
		gaux*base_dphiz1[idofs]  +
                eaux*base_phi1[idofs]
					) * (implicit-1.0) * 0.5*area;

	  }/* idofs */

	  kk = PDC_CONV_DIFF_NREQ*num_shap1+ieq1*num_shap2;

	  for (idofs=0;idofs<num_shap2;idofs++) {

	    Rhs_vect[kk+idofs] += (
                -eaux*base_phi2[idofs]
					) * (implicit-1.0) * 0.5*area;


	  }/* idofs */


	  daux=0.0;faux=0.0;eaux=0.0;gaux=0.0;
	  for(ieq2=0;ieq2<PDC_CONV_DIFF_NREQ;ieq2++){
	    daux +=
           	  anx2[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_val1[ieq2]
           	- anx2[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_val2[ieq2];
	    faux +=
           	  any2[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_val1[ieq2]
           	- any2[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_val2[ieq2];
	    gaux +=
           	  anz2[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_val1[ieq2]
           	- anz2[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_val2[ieq2];
	    eaux += 
           	  anx2[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_x2[ieq2]
           	+ any2[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_y2[ieq2]
           	+ anz2[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_z2[ieq2];
	  }     

	  kk = ieq1*num_shap1;

	  for (idofs=0;idofs<num_shap1;idofs++) {

	    Rhs_vect[kk+idofs] += (
                -eaux*base_phi1[idofs]
					) * (implicit-1.0) * 0.5*area;

	  }/* idofs */

	  kk = PDC_CONV_DIFF_NREQ*num_shap1+ieq1*num_shap2;

	  for (idofs=0;idofs<num_shap2;idofs++) {

	    Rhs_vect[kk+idofs] += (
		daux*base_dphix2[idofs]  +
		faux*base_dphiy2[idofs]  +
		gaux*base_dphiz2[idofs]  +
                eaux*base_phi2[idofs]
					) * (implicit-1.0) * 0.5*area;

	  }/* idofs */
	  
/* upwind terms - norm_vel[ieq1] gives normal velocity  !!! */
	  if(norm_vel[ieq1]>0) {
/* upwind terms - we assume b_i*n_i is diagonal !!! */
	    haux = bn1[ieq1*PDC_CONV_DIFF_NREQ+ieq1] * u_val1[ieq1];
	  }
	  else if(norm_vel[ieq1]<0) {
	    haux = bn1[ieq1*PDC_CONV_DIFF_NREQ+ieq1] * u_val2[ieq1];
	  }
	    
	  kk = ieq1*num_shap1;

	  for (idofs=0;idofs<num_shap1;idofs++) {

	    Rhs_vect[kk+idofs] +=
	      haux * base_phi1[idofs] * (implicit-1.0) * area;

	  }/* idofs */

 /* upwind terms */
	  if(norm_vel[ieq1]>0) {
	    haux = bn2[ieq1*PDC_CONV_DIFF_NREQ+ieq1] * u_val1[ieq1];
	  }
	  else if(norm_vel[ieq1]<0) {
	    haux = bn2[ieq1*PDC_CONV_DIFF_NREQ+ieq1] * u_val2[ieq1];
	  }
	    
	  kk = PDC_CONV_DIFF_NREQ*num_shap1+ieq1*num_shap2;

	  for (idofs=0;idofs<num_shap2;idofs++) {

	    Rhs_vect[kk+idofs] += 
	      - haux * base_phi2[idofs] * (implicit-1.0) * area;

	  }/* idofs */

	}/* ieq1 */
      } /* if time dependent and partially explicit */

    }
    else if(bc_type==PDC_BC_DIRI){

      for(ieq2=0;ieq2<PDC_CONV_DIFF_NREQ;ieq2++){
        for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){

/* fortran-like style - columnwise storage of matrices */
          kk=(ieq2*PDC_CONV_DIFF_NREQ*num_shap1+ieq1)*num_shap1;

          for (jdofs=0;jdofs<num_shap1;jdofs++) {
            for (idofs=0;idofs<num_shap1;idofs++) {

/* stiffness matrix */ 
              Stiff_mat[kk+idofs] += (
			(
		(
	    anx1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]*base_dphix1[idofs]
	   +any1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]*base_dphiy1[idofs]
	   +anz1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]*base_dphiz1[idofs]
		)*base_phi1[jdofs] +
		(
	   -anx1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]*base_dphix1[jdofs]
	   -any1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]*base_dphiy1[jdofs]
	   -anz1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]*base_dphiz1[jdofs]
		)*base_phi1[idofs]
			)  * implicit
/*penalty*/
	 + penalty*base_phi1[jdofs]*base_phi1[idofs]/hf_size
/*penalty*/
                                         ) * area;

            }/* idofs */

            kk+=num_shap1*PDC_CONV_DIFF_NREQ;

          } /* jdofs */
        }/* ieq1 */
      } /* ieq2 */

      for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){
        for (idofs=0;idofs<num_shap1;idofs++) {

/* right hand side vector */
          daux=0.0;faux=0.0;eaux=0.0;gaux=0.0;
          for(ieq2=0;ieq2<PDC_CONV_DIFF_NREQ;ieq2++){
            daux +=
           	  anx1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * fval1[ieq2];
            faux +=
           	  any1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * fval1[ieq2];
            gaux +=
           	  anz1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * fval1[ieq2];
	  }

          Rhs_vect[ieq1*num_shap1+idofs] += (
		daux*base_dphix1[idofs] +
		faux*base_dphiy1[idofs] +
		gaux*base_dphiz1[idofs] 
/*penalty*/
	 + penalty*fval1[ieq1]*base_phi1[idofs]/hf_size
/*penalty*/
					) * area;

        }/* ieq1 */
      }/* idofs */

/* explicit fluxes (no upwinding on Dirichlet boundary !!! - mixed BC better) */
      for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){
	for (idofs=0;idofs<num_shap1;idofs++) {
 /* right hand side vector */
	  Rhs_vect[ieq1*num_shap1+idofs] += 
	    -qn1[ieq1] * base_phi1[idofs] * area;

	}/* idofs */
      }/* ieq1 */

/* inflow-outflow condition - we assume b_i*n_i is diagonal !!! */
      for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){

/* Dirichlet inflow boundary */
        if(norm_vel[ieq1]<0) {

          for (idofs=0;idofs<num_shap1;idofs++) {

/* general case
	    daux=0.0;
	    for(ieq2=0;ieq2<PDC_CONV_DIFF_NREQ;ieq2++){
	      daux += bn1[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * fval1[ieq2];
	    }
*/
/*  we assume b_i*n_i is diagonal !!! */
	    daux = bn1[ieq1*PDC_CONV_DIFF_NREQ+ieq1] * fval1[ieq1];

 /* right hand side vector */
            Rhs_vect[ieq1*num_shap1+idofs] += 
		-daux * base_phi1[idofs] * area;

          }/* idofs */
        }

/* Dirichlet outflow boundary */
        else {

/* fortran-like style - columnwise storage of matrices */
          kk=(ieq1*PDC_CONV_DIFF_NREQ*num_shap1+ieq1)*num_shap1;

          for (jdofs=0;jdofs<num_shap1;jdofs++) {
            for (idofs=0;idofs<num_shap1;idofs++) {

              Stiff_mat[kk+idofs] += implicit * area *
	bn1[ieq1*PDC_CONV_DIFF_NREQ+ieq1] * base_phi1[idofs] * base_phi1[jdofs];

            }/* idofs */

            kk+=num_shap1*PDC_CONV_DIFF_NREQ;

          } /* jdofs */
        } /* if outflow boundary */
      }/* ieq1 */



/* right hand side correction due to an and bn for time dependent problems */
      if(implicit<1.0){
	for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){
	  for (idofs=0;idofs<num_shap1;idofs++) {


	    daux=0.0;faux=0.0;eaux=0.0;gaux=0.0;
	    for(ieq2=0;ieq2<PDC_CONV_DIFF_NREQ;ieq2++){
	      daux +=
		  anx1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_val1[ieq2];
	      faux +=
           	  any1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_val1[ieq2];
	      gaux +=
           	  anz1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_val1[ieq2];
	      eaux += 
           	- anx1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_x1[ieq2]
           	- any1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_y1[ieq2]
           	- anz1[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  * u_z1[ieq2];
	    }     


/* upwind terms */
	    if(norm_vel[ieq1]>0) {
	      eaux += bn1[ieq1*PDC_CONV_DIFF_NREQ+ieq1] * u_val1[ieq1];
	    }

	    Rhs_vect[ieq1*num_shap1+idofs] += (
		daux*base_dphix1[idofs]  +
		faux*base_dphiy1[idofs]  +
		gaux*base_dphiz1[idofs]  +
                eaux*base_phi1[idofs]
					) * (implicit-1.0) * area;

	  }/* idofs */
	}/* ieq1 */
      } /* if time dependent */

    }
    else if(bc_type==PDC_BC_NEUM){

      for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){
	if(norm_vel[ieq1]>0) {

          kk=(ieq1*PDC_CONV_DIFF_NREQ*num_shap1+ieq1)*num_shap1;

          for (jdofs=0;jdofs<num_shap1;jdofs++) {
            for (idofs=0;idofs<num_shap1;idofs++) {

              Stiff_mat[kk+idofs] += implicit * area *
	bn1[ieq1*PDC_CONV_DIFF_NREQ+ieq1] * base_phi1[idofs] * base_phi1[jdofs];

            }/* idofs */

            kk+=num_shap1*PDC_CONV_DIFF_NREQ;

          } /* jdofs */

        }/* end if outflow */

      } /* ieq1 */
      for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){
        kk=ieq1*num_shap1;

        for (idofs=0;idofs<num_shap1;idofs++) {

/* right hand side vector */
          Rhs_vect[kk] += gval1[ieq1]*base_phi1[idofs]*area;
          kk++;

        }/* idofs */
      }/* ieq1 */

/* explicit fluxes (no upwinding) */
      for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){
	for (idofs=0;idofs<num_shap1;idofs++) {
 /* right hand side vector */
	  Rhs_vect[ieq1*num_shap1+idofs] += 
	    -qn1[ieq1] * base_phi1[idofs] * area;
	}/* idofs */
      }/* ieq1 */

/* right hand side correction due to bn for time dependent problems */
      if(implicit<1.0){

	for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){

/* inflow-outflow condition - we assume b_i*n_i is diagonal !!! */
          if(norm_vel[ieq1]>0) {
	    for (idofs=0;idofs<num_shap1;idofs++) {

	      Rhs_vect[ieq1*num_shap1+idofs] += (implicit-1.0) * area *
                bn1[ieq1*PDC_CONV_DIFF_NREQ+ieq1] * u_val1[ieq1] * base_phi1[idofs];

	    }/* idofs */
	  }/* end if outflow */
	}/* ieq1 */
      } /* if time dependent */


    } /* if Neumann BC */
    else if(bc_type==PDC_BC_MIXED){

      for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){

/* fortran-like style - columnwise storage of matrices */
	kk=ieq1*num_shap1;

	for (jdofs=0;jdofs<num_shap1;jdofs++) {
	  for (idofs=0;idofs<num_shap1;idofs++) {

/* stiffness matrix */ 
	    Stiff_mat[kk+idofs] += (
/* normal velocity norm_vel gives proportionality constant between */
/* the difference in values and the flux */
	 - norm_vel[ieq1] * base_phi1[jdofs]*base_phi1[idofs] 
/*kb!!! penalty?
	 + penalty/hf_size * base_phi1[jdofs]*base_phi1[idofs]
/*kb!!! penalty*/
	                              ) * area;

	  }/* idofs */

	  kk+=num_shap1*PDC_CONV_DIFF_NREQ;

	} /* jdofs */
      }/* ieq1 */

      for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){
        for (idofs=0;idofs<num_shap1;idofs++) {

/* right hand side vector */

          Rhs_vect[ieq1*num_shap1+idofs] += (
	 - norm_vel[ieq1] *fval1[ieq1]*base_phi1[idofs]
/*kb!!! penalty?
	 + penalty/hf_size *fval1[ieq1]*base_phi1[idofs]
/*kb!!! penalty*/
					) * area;

        }/* ieq1 */
      }/* idofs */

/* explicit fluxes (no upwinding !!!)  */
      for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){
	for (idofs=0;idofs<num_shap1;idofs++) {
 /* right hand side vector */
	  Rhs_vect[ieq1*num_shap1+idofs] += 
	    -qn1[ieq1] * base_phi1[idofs] * area;

	}/* idofs */
      }/* ieq1 */

/* inflow-outflow condition - we assume b_i*n_i is diagonal !!! */
      for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){

/* mixed inflow boundary */
        if(norm_vel[ieq1]<0) {

          for (idofs=0;idofs<num_shap1;idofs++) {

	    daux = bn1[ieq1*PDC_CONV_DIFF_NREQ+ieq1] * fval1[ieq1];

 /* right hand side vector */
            Rhs_vect[ieq1*num_shap1+idofs] += 
		-daux * base_phi1[idofs] * area;

          }/* idofs */
        }

/* mixed outflow boundary */
        else {

/* fortran-like style - columnwise storage of matrices */
          kk=(ndofs_fa*ieq1+ieq1)*num_shap1;

          for (jdofs=0;jdofs<num_shap1;jdofs++) {
            for (idofs=0;idofs<num_shap1;idofs++) {

              Stiff_mat[kk+idofs] += implicit * area *
	bn1[ieq1*PDC_CONV_DIFF_NREQ+ieq1] * base_phi1[idofs] * base_phi1[jdofs];

            }/* idofs */

            kk+=ndofs_fa;

          } /* jdofs */
        } /* if outflow boundary */
      }/* ieq1 */

/* right hand side correction due to an and bn for time dependent problems */
      if(implicit<1.0){

	for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){
/* upwind terms */
	  if(norm_vel[ieq1]>0) {

	    for (idofs=0;idofs<num_shap1;idofs++) {
	      Rhs_vect[ieq1*num_shap1+idofs] +=
	        bn1[ieq1*PDC_CONV_DIFF_NREQ+ieq1] * u_val1[ieq1] * (implicit-1.0) * area;
	    }/* idofs */

	  } /* if outflow */
	}/* ieq1 */
      } /* if time dependent */
    }/* if mixed type BC */
 
  } /* ki */


/*kbw
   if(Fa_id>0){
     printf("Face %d, neig1 %d, neig2 %d\n",
	    Fa_id,neig1_id,neig2_id);
     printf("node_coor %lf, %lf, %lf, %lf\n",
	    node_coor1[0],node_coor1[1],node_coor1[2],node_coor1[3]);
      printf("Stiffness matrix:\n");
     for (idofs=0;idofs<ndofs_fa;idofs++) {
       for (jdofs=0;jdofs<ndofs_fa;jdofs++) {
	 printf("%20.12lf",Stiff_mat[idofs+jdofs*ndofs_fa]);        
       } 
       printf("\n");
     }
     printf("Rhs_vect:\n");
     for (idofs=0;idofs<ndofs_fa;idofs++) {
       printf("%20.12lf",Rhs_vect[idofs]);
     }
     printf("\n");
     getchar();  
   }
/*kew*/

   if(Rewr_dofs != NULL) *Rewr_dofs = 'F';

  } /* end if computing entries to the stiffness matrix */

/*kbw
#ifdef DEBUG
    {
      int ibl;
      printf("leaving pdr_fa_stiff_mat:\n");
      printf("face %d, nrdofbl %d\n", Fa_id, *Nr_dof_ent);
      for(ibl=0;ibl<*Nr_dof_ent; ibl++){
	printf("bl_id %d, bl_nrdof %d\n",
	  List_dof_ent_id[ibl],List_dof_ent_nrdof[ibl]);
      }
    }
#endif
/*kew*/
  
  return(1);
}
