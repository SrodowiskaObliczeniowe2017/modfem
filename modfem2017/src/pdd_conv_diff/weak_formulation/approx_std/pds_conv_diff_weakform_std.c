/************************************************************************
File pds_conv_diff_weakform_std.c - weakform functions for conv_diff and std approximation

Contains definitions of routines:

  pdr_conv_diff_comp_el_stiff_mat_app - to construct UNCONSTRAINED stiffness matrix and
                          a load vector for an element

  pdr_conv_diff_comp_fa_stiff_mat_app - to construct UNCONSTRAINED stiffness matrix and
                          a load vector for a face

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
                          and a load vector for an element - special routine for
                          standard approximation with constrained nodes
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
  int* List_dof_ent_nrdofs,/* out: list of no of dofs for 'dof' entity */
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


  // to make OpenMP working
  int ndofs_max=0;

  int field_id, mesh_id, nreq;  
  int kk, idofs, jdofs, num_dofs;
  int i,j,k;
  int max_dof_ent, max_nrdofs;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper field */
  field_id = Problem_id;

  if(Comp_sm!=PDC_NO_COMP){

    /* select the proper mesh */
    mesh_id = apr_get_mesh_id(field_id);
    nreq = PDC_CONV_DIFF_NREQ;

    if(nreq!=1){
      printf("conv_diff module for scalar problems only !!! nreq!=1. Exiting\n");
      exit(-1);
    }

    /* get material number */
    el_mate =  mmr_el_groupID(mesh_id, El_id);
    
    /* find degree of polynomial and number of element scalar dofs */
    if(Pdeg_in>0) pdeg = Pdeg_in;
    else{
      apr_get_el_pdeg(field_id, El_id, &pdeg);
    }
    num_shap = apr_get_el_pdeg_numshap(field_id, El_id, &pdeg);
    num_dofs = num_shap*nreq;
    
    //printf("pdeg %d, num_shap=%d\n",pdeg,num_shap);
    
    max_dof_ent = *Nr_dof_ent;
    
#ifdef DEBUG
    if(max_dof_ent < 1){
      printf("Too small arrays List_dof_... passed to comp_el_stiff_mat\n");
      printf("%d < %d. Exiting !!!",max_dof_ent, *Nr_dof_ent);
      exit(-1);
    }
#endif
    
    
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


    /* get the most recent solution degrees of freedom */
    if(mmr_el_status(mesh_id,El_id)==MMC_ACTIVE){

      sol_vec_id = 1;
      apr_get_el_dofs(field_id, El_id, sol_vec_id, sol_dofs);
      
    }
    else{

     /*!!! coarse element dofs should be supplied by calling routine !!!*/

      for(i=0;i<num_dofs;i++) sol_dofs[i]=0;

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
    apr_num_int_el(Problem_id, field_id, El_id, Comp_sm, &pdeg, sol_dofs, sol_dofs,
		   diagonal, Stiff_mat, Rhs_vect);


/*kbw
    if(El_id==13753||El_id==13754){
      printf("Element %d: Stiffness matrix:\n",El_id);
      for (idofs=0;idofs<num_shap*nreq;idofs++) {
	for (jdofs=0;jdofs<num_shap*nreq;jdofs++) {
	  printf("%20.12lf",Stiff_mat[idofs+jdofs*num_shap*nreq]);        
	} 
	printf("\n");
      }
      printf("Element %d: Rhs_vect:\n",El_id);
      for (idofs=0;idofs<num_shap*nreq;idofs++) {
	printf("%20.12lf",Rhs_vect[idofs]);
      }
      printf("\n");
      //getchar(); 
    } 
/*kew*/

    if(Rewr_dofs != NULL) *Rewr_dofs = 'F';

  } /* end if computing SM and/or RHSV */

  /* change the option compute SM and RHSV to rewrite SM and RHSV */
  if(Comp_sm!=PDC_NO_COMP) Comp_sm += 3;

  /* obligatory procedure to fill Lists of dof_ents and rewite SM and RHSV */
  /* the reason is to take into account POSSIBLE CONSTRAINTS (HANGING NODES) */
  apr_get_stiff_mat_data(field_id, El_id, Comp_sm, 'N', Pdeg_in, 0, Nr_dof_ent, 
	         List_dof_ent_type, List_dof_ent_id, List_dof_ent_nrdofs,
			   Nrdofs_loc, Stiff_mat, Rhs_vect);

/* matrix displayed by rows, altghough stored by columns !!!!!!!!!*/
/*kbw
if(Comp_sm!=PDC_NO_COMP && El_id>0){
  printf("Element %d: Modified stiffness matrix:\n",El_id);
  for (idofs=0;idofs<*Nrdofs_loc;idofs++) { //for each row!!!!
    for (jdofs=0;jdofs<*Nrdofs_loc;jdofs++) { // for each element in row !!!
      printf("%7.3lf",Stiff_mat[idofs+jdofs*(*Nrdofs_loc)]);
    }
    printf("\n");
  }
  printf("Element %d: Rhs_vect:\n",El_id);
  for (idofs=0;idofs<*Nrdofs_loc;idofs++) {
    printf("%7.3lf",Rhs_vect[idofs]);
  }
  printf("\n");
  getchar();
 }
/*kew*/

/*kbw
    printf("In pdr_comp_el_stiff_mat: field_id %d, El_id %d, Comp_sm %d, Nr_dof_ent %d\n",
	   field_id, El_id, Comp_sm, *Nr_dof_ent);
    printf("For each block: \ttype, \tid, \tnrdofs\n");
    for(i=0;i<*Nr_dof_ent;i++){
      printf("\t\t\t%d\t%d\t%d\n",
	     List_dof_ent_type[i],List_dof_ent_id[i],List_dof_ent_nrdofs[i]);
    }
    getchar();getchar();
/*kew*/

  return(1);
}

/*------------------------------------------------------------
  pdr_conv_diff_comp_fa_stiff_mat_app - to construct a stiffness matrix  
                          and a load vector for a face - special routine for  
                          standard approximation with constrained nodes
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
  int* List_dof_ent_nrdofs,/* out: list of no of dofs for 'dof' entity */
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
  int base_1;		/* type of basis functions */
  int nreq;		/* number of equations */
  int pdeg1=0, num_shap1=0;/* local dimension for neigbhor 1 */
  int el_mate1=0;	/* material number for neigbhor 1 */
  int ngauss;		/* number of gauss points in 1 dimension */
  int face_neig[2];	/* list of neighbors */
  int neig_sides[2]={0,0};/* sides of face wrt neigs */
  int node_shift;
  double acoeff[4], bcoeff[2];/* to transform coordinates between faces */
  double loc_xg[6];	/* local coord of gauss points for neighbors */
  double xg[450];   	/* gauss points in 2D*/
  double wg[150];      	/* gauss weights in 2D*/
  double determ; 	/* determinant of jacobi matrix */
  double hsize1=0.0;        /* size of elements */
  double hf_size;        /* size of face */
  double area; 		/* area for integration rule */
  double time;          /* current time instant */
  double xcoor[3],xcoor2[3];/* global coord of gauss point */
  double vec_norm[3];	/* normal vector to the considered face */
  double u_val1[PDC_CONV_DIFF_NREQ]={0.0}; /* computed solution */
  double u_x1[PDC_CONV_DIFF_NREQ]={0.0};   /* gradient of computed solution */
  double u_y1[PDC_CONV_DIFF_NREQ]={0.0};   /* gradient of computed solution */
  double u_z1[PDC_CONV_DIFF_NREQ]={0.0};   /* gradient of computed solution */
  double base_phi1[APC_MAXELVD]={0.0};    /* basis functions */
  double base_dphix1[APC_MAXELVD]={0.0};  /* x-derivatives of basis function */
  double base_dphiy1[APC_MAXELVD]={0.0};  /* y-derivatives of basis function */
  double base_dphiz1[APC_MAXELVD]={0.0};  /* y-derivatives of basis function */
  int el_nodes1[MMC_MAXELVNO+1]={0.0};        /* list of nodes of El */
  double node_coor1[3*MMC_MAXELVNO]={0.0};  /* coord of nodes of El */
  double dofs_loc1[APC_MAXELSD]={0.0}; /* element solution dofs */
  double anx1[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ]={0.0};  /* coeff of PDE */
  double any1[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ]={0.0};  /* coeff of PDE */
  double anz1[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ]={0.0};  /* coeff of PDE */
  double bn1[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ]={0.0};   /* coeff of PDE */
  double fval1[PDC_CONV_DIFF_NREQ]={0.0}; /* coeff of PDE */
  double gval1[PDC_CONV_DIFF_NREQ]={0.0}; /* coeff of PDE */
  double qn1[PDC_CONV_DIFF_NREQ]={0.0};      /* source for the first neigbor */
  double norm_vel[PDC_CONV_DIFF_NREQ]={0.0}; /* normal velocity for each component of u */

  int i, ki, kk, ieq1, ieq2, time_dis;
  int iaux, jaux;
  int idofs, jdofs, ndofs_el, neig1_id;
  double daux, faux, eaux, gaux, haux, penalty=0.0, implicit;


  int field_id, mesh_id;  
  int sol_vec_id;       /* indicator for the solution dofs */


/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper mesh */
  pdv_conv_diff_current_problem_id = Problem_id; 
  field_id = Problem_id;
  mesh_id = apr_get_mesh_id(field_id);
  
  time=pdr_time_d_params(PDC_USE_CURRENT_PROBLEM_ID, 4);

  /* get type of face and bc flag */
  fa_type=mmr_fa_type(mesh_id,Fa_id);
  fa_bc=mmr_fa_bc(mesh_id,Fa_id);
  mmr_fa_area(mesh_id,Fa_id,&daux,NULL);
  hf_size = sqrt(daux);

  bc_type = pdr_get_bc_type(fa_bc);

/*kbw
  printf("in pdr_conv_diff_comp_fa_stiff_mat_app for STD_LIN and face %d\n",
	 Fa_id);
  printf("Problem_id %d, field_id %d, mesh_id %d, fa_type %d, fa_bc %d, hf_size %lf, bc_type %d\n",
	 Problem_id, field_id, mesh_id, fa_type, fa_bc, hf_size, bc_type);
/*kew*/

  /* get approximation parameters */
  nreq=PDC_CONV_DIFF_NREQ;

  if(nreq!=1){
    printf("conv_diff module for scalar problems only !!! nreq!=1. Exiting\n");
    exit(-1);
  }

  /* get neigbors list with corresponding neighbors' sides numbers*/
  /* matrices acoeff, bcoeff are necessary for mmr_fa_elem_coor to compute */
  /* coordinates within element for a point with given coordinates on a face */
  mmr_fa_neig(mesh_id,Fa_id,face_neig,neig_sides,&node_shift,NULL,acoeff,bcoeff);

  /* get material number */
  neig1_id = abs(face_neig[0]);
  el_mate1 =  mmr_el_groupID(mesh_id, neig1_id);
  base_1=apr_get_base_type(field_id, neig1_id);

  if(Pdeg_in>0) pdeg1 = Pdeg_in;
  else{
    /* find degree of polynomial and number of element scalar dofs */
    /*!!! need some trick for coarse solve with different Pdeg !!!*/
    apr_get_el_pdeg(field_id, neig1_id, &pdeg1);
  }

  num_shap1 = apr_get_el_pdeg_numshap(field_id, neig1_id, &pdeg1);
  ndofs_el = num_shap1*nreq;
  
  if(Comp_sm!=PDC_NO_COMP){

    /* find vertices */
    mmr_el_node_coor(mesh_id,neig1_id,el_nodes1,node_coor1);

    /* get the most recent solution degrees of freedom */
    /*!!! coarse element dofs should be supplied by calling routine !!!*/
    sol_vec_id = 1;
    apr_get_el_dofs(field_id,neig1_id,sol_vec_id,dofs_loc1);



#ifdef DEBUG
    if(*Nrdofs_loc < ndofs_el ){
      printf("Too small arrays Stiff_mat and Rhs_vect passed to comp_el_stiff_mat\n");
      printf("%d < %d. Exiting !!!", *Nrdofs_loc, ndofs_el);
      exit(-1);
    }
#endif

    /* initialize the matrices to zero - number of DOFs from element !!!!!! */
    for(i=0;i<ndofs_el*ndofs_el;i++) Stiff_mat[i]=0.0;
    for(i=0;i<ndofs_el;i++) Rhs_vect[i]=0.0;

   /* prepare data for gaussian integration */
   /* !!!!! for triangular faces coordinates are standard [0,1][0,1] */
   /* !!!!! for quadrilateral faces coordinates are [-1,1][-1,1] - */
   /* !!!!! which means that they do not conform to element coordinates */
   /* !!!!! proper care is taken in mmr_fa_elem_coor !!!!!!!!!!!!!!!!!! */
   //   if(bc_type==PDC_BC_DIRI){
   //apr_set_quadr_2D_penalty(fa_type,base_1,&pdeg1,&ngauss,xg,wg);
   //}
   //else{
    apr_set_quadr_2D(fa_type,base_1,&pdeg1,&ngauss,xg,wg);
     //}

/* loop over integration points */
   for (ki=0;ki<ngauss;ki++) {

/* find coordinates within neighboring elements for a point on face */
   /* !!!!! for triangular faces coordinates are standard [0,1][0,1] */
   /* !!!!! for quadrilateral faces coordinates are [-1,1][-1,1] - */
   /* !!!!! which means that they do not conform to element coordinates */
   /* !!!!! proper care is taken in mmr_fa_elem_coor !!!!!!!!!!!!!!!!!! */
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
		   anx1,any1,anz1,bn1,fval1,gval1,qn1,norm_vel);


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
printf("NORMAL velocity: %lf\n", norm_vel[0]);
printf("coord for neigh1: %lf, %lf, %lf, hsize %lf\n",
loc_xg[0],loc_xg[1],loc_xg[2],hsize1);
printf("weight %lf, integration coeff. %lf\n", wg[ki], area);
printf("%d shape functions and derivatives for neigh1: \n",num_shap1);
for(i=0;i<num_shap1;i++){
  printf("fun - %lf, der: x - %lf, y - %lf, z - %lf\n",
  base_phi1[i],base_dphix1[i],base_dphiy1[i],base_dphiz1[i]);
}
  printf("solution %lf, der %lf, %lf, %lf\n", 
	 u_val1[0],u_x1[0],u_y1[0],u_z1[0]);
  printf("diffusion coeff : %lf %lf %lf\n",anx1[0],any1[0],anz1[0]);
  printf("convection coeff: %lf\n",bn1[0]);
  printf("rhs coeff       : %lf %lf %lf\n",fval1[0],gval1[0],qn1[0]);
  printf("penalty %lf : fval1[0] %lf, hf_size %lf\n",
       penalty,fval1[0],hf_size);
 
getchar();
}
/*kew*/ 



    penalty = 1.0e7;
    if(bc_type==PDC_BC_DIRI){
      
      kk=0;
      
      for (jdofs=0;jdofs<num_shap1;jdofs++) {
	for (idofs=0;idofs<num_shap1;idofs++) {
	  
	  /* stiffness matrix */
	  Stiff_mat[kk+idofs] += 
	    /*penalty*/
	    penalty*base_phi1[jdofs]*base_phi1[idofs]
	    /*penalty*/
	    * area;
	  
	}/* idofs */
	
	kk+=num_shap1;

      } /* jdofs */
      
      for (idofs=0;idofs<num_shap1;idofs++) {
	
	
	Rhs_vect[idofs] +=
	  /*penalty*/
	  penalty*fval1[0]*base_phi1[idofs]
	  /*penalty*/
	  * area;
	
      }/* idofs */

    }
    else if(bc_type==PDC_BC_NEUM){
      
      for (idofs=0;idofs<num_shap1;idofs++) {
	
	/* right hand side vector */
	Rhs_vect[idofs] += gval1[0]*base_phi1[idofs]*area;
	
      }/* idofs */

    } /* if Neumann BC */
    else if(bc_type==PDC_BC_MIXED){
      
      

    }/* if mixed type BC */


  } /* end loop over integration points ki */

/*kbw
   if(Fa_id>0){
     printf("Stiffness matrix face=%d:\n",Fa_id);
     for (idofs=0;idofs<ndofs_el;idofs++) {
       for (jdofs=0;jdofs<ndofs_el;jdofs++) {
	 printf("%20.12lf",Stiff_mat[idofs+jdofs*ndofs_el]);
       }
       printf("\n");
     }
     printf("Rhs_vect:\n");
     for (idofs=0;idofs<ndofs_el;idofs++) {
       printf("%20.12lf",Rhs_vect[idofs]);
     }
     printf("\n");
     getchar();
   }
/*kew*/

   if(Rewr_dofs != NULL) *Rewr_dofs = 'F';

/* matrix displayed by rows, altghough stored by columns !!!!!!!!!*/
/*kbw
if(Fa_id>0){
  int num_dofs=ndofs_el;
  printf("Face %d, Element %d: Standard stiffness matrix:\n",Fa_id,neig1_id);
  for (idofs=0;idofs<num_dofs;idofs++) { //for each row!!!!
    for (jdofs=0;jdofs<num_dofs;jdofs++) { // for each element in row !!!
      printf("%11.3le",Stiff_mat[idofs+jdofs*num_dofs]);
    }
    printf("\n");
  }
  printf("Element %d: Rhs_vect:\n",neig1_id);
  for (idofs=0;idofs<num_dofs;idofs++) {
    printf("%11.3le",Rhs_vect[idofs]);
  }
  printf("\n");
  getchar();
 }
/*kew*/

  } /* end if computing entries to the stiffness matrix */

/*kbw
#ifdef DEBUG
    {
      int ibl;
      printf("leaving pdr_fa_stiff_mat:\n");
      printf("face %d, nrdofsbl %d\n", Fa_id, *Nr_dof_ent);
      for(ibl=0;ibl<*Nr_dof_ent; ibl++){
	printf("bl_id %d, bl_nrdofs %d\n",
	  List_dof_ent_id[ibl],List_dof_ent_nrdofs[ibl]);
      }
    }
#endif
/*kew*/

  /* change the option compute SM and RHSV to rewrite SM and RHSV */
  if(Comp_sm!=PDC_NO_COMP) Comp_sm += 3;

  /* obligatory procedure to fill Lists of dof_ents and rewite SM and RHSV */
  /* the reason is to take into account POSSIBLE CONSTRAINTS (HANGING NODES) */
  apr_get_stiff_mat_data(field_id,neig1_id,Comp_sm,'N',Pdeg_in,0,Nr_dof_ent, 
	         List_dof_ent_type, List_dof_ent_id, List_dof_ent_nrdofs,
			   Nrdofs_loc, Stiff_mat, Rhs_vect);

/* matrix displayed by rows, altghough stored by columns !!!!!!!!!*/
/*kbw
if(Comp_sm!=PDC_NO_COMP && Fa_id>0){
  int num_dofs=ndofs_el;
  printf("Face %d, Element %d: Modified stiffness matrix:\n",Fa_id,neig1_id);
  for (idofs=0;idofs<*Nrdofs_loc;idofs++) { //for each row!!!!
    for (jdofs=0;jdofs<*Nrdofs_loc;jdofs++) { // for each element in row !!!
      printf("%11.3le",Stiff_mat[idofs+jdofs*(*Nrdofs_loc)]);
    }
    printf("\n");
  }
  printf("Element %d: Rhs_vect:\n",neig1_id);
  for (idofs=0;idofs<*Nrdofs_loc;idofs++) {
    printf("%11.3le",Rhs_vect[idofs]);
  }
  printf("\n");
  getchar();
 }
/*kew*/

/*kbw
    printf("In pdr_comp_el_stiff_mat: field_id %d, El_id %d, Comp_sm %d, Nr_dof_ent %d\n",
	   field_id, El_id, Comp_sm, *Nr_dof_ent);
    printf("For each block: \ttype, \tid, \tnrdofs\n");
    for(i=0;i<*Nr_dof_ent;i++){
      printf("\t\t\t%d\t%d\t%d\n",
	     List_dof_ent_type[i],List_dof_ent_id[i],List_dof_ent_nrdofs[i]);
    }
    getchar();getchar();
/*kew*/

  return(1);
}

