/************************************************************************
File pds_conv_diff_weakform.c - weakform functions for conv_diff 

Contains definitions of routines:

MODULE PROVIDES IMPLEMENTATION (in pds_conv_diff_weakform.c) 
PROCEDURES CAN BE CALLED BY ALL OTHER PROBLEM MODULES 
  pdr_conv_diff_select_el_coeff_vect - to select coefficients 
                        returned to approximation routines for element integrals
  OBSOLETE pdr_conv_diff_select_el_coeff - to select coefficients 
                        returned to approximationroutines for element integrals

  pdr_conv_diff_el_coeff - to return coefficients for element integrals

  pdr_conv_diff_comp_el_stiff_mat - to construct UNCONSTRAINED stiffness matrix and
                          a load vector for an element


  pdr_select_fa_coeff - to select coefficients returned to approximation
                        routines for face integrals
  pdr_fa_coeff - to return coefficients for face integrals
  pdr_pde_coeff - to return coefficients of the original pdes
  pdr_bc_diri_coeff - to return DIRICHLET boundary coeficients
  pdr_bc_neum_coeff - to return NEUMANN boundary coeficients
  pdr_bc_mixed_coeff - to return ROBIN (mixed) boundary coeficients
  pdr_bc_val - to get the boundary condition values
  pdr_ic_val - to get the initial condition values

  pdr_conv_diff_comp_fa_stiff_mat - to construct UNCONSTRAINED stiffness matrix and
                          a load vector for a face
  pdr_conv_diff_compute_CFL - to compute global CFL numbers (for a subdomain)


------------------------------
History:
	2012    - Krzysztof Banas (pobanas@cyf-kr.edu.pl)
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<string.h>


/* USED DATA STRUCTURES AND INTERFACES FROM OTHER MODULES */

/* utilities - including simple time measurement library */
#include "uth_intf.h"

/* interface for all mesh manipulation modules */
#include "mmh_intf.h"	

/* interface for all approximation modules */
#include "aph_intf.h"	

/* interface for thread management modules */
#include "tmh_intf.h"

/* interface for all solver modules */
#include "sih_intf.h"	


/* USED AND IMPLEMENTED PROBLEM DEPENDENT DATA STRUCTURES AND INTERFACES */

/* problem dependent module interface */
#include "pdh_intf.h"	

/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include "pdh_control_intf.h"

/* header files for the problem dependent module for Conv_Diff equations */
#include "../include/pdh_conv_diff.h"	
#include "../include/pdh_conv_diff_problem.h"	


/***************************************/
/* DEFINITIONS OF PROCEDURES */
/***************************************/
#define SMALL      1.0e-10

/*------------------------------------------------------------
  pdr_conv_diff_select_el_coeff_vect - to select coefficients 
                     returned to approximation routines for element integrals
------------------------------------------------------------*/
int pdr_conv_diff_select_el_coeff_vect( // returns success indicator
  int Problem_id,
  int *Coeff_vect_ind	/* out: coefficient indicator */
				   )
{

  // initialize Coeff_vect_ind
  int i;
  for(i=0; i<=22; i++) Coeff_vect_ind[i]=0;

  // indicate the vector has been filled
  Coeff_vect_ind[0] = 1;


  //int name=pdr_ctrl_i_params(Problem_id,1); changed for string
  char name[300];
  pdr_problem_name(Problem_id, name);

  if(strcmp(name, "DIFF_IN_CUBE") == 0 || strcmp(name, "DIFF_IN_BOX") == 0){
    
    //if(name>=1&&name<100){

    // Axx
    Coeff_vect_ind[2] = 1;
    // Ayy
    Coeff_vect_ind[6] = 1;
    // Azz
    Coeff_vect_ind[10] = 1;
    // Sval
    Coeff_vect_ind[22] = 1;
  
  }
  else if(strcmp(name, "TEST_NUM_INT") == 0){
    
    // Mval - possibly mass matrix in time integration (may be included in Cval)
    Coeff_vect_ind[1] = 0;
    // Axx
    Coeff_vect_ind[2] = 1;
    // Axy
    Coeff_vect_ind[3] = 1;
    // Axz
    Coeff_vect_ind[4] = 1;
    // Ayx
    Coeff_vect_ind[5] = 1;
    // Ayy
    Coeff_vect_ind[6] = 1;
    // Ayz
    Coeff_vect_ind[7] = 1;
    // Azx
    Coeff_vect_ind[8] = 1;
    // Azy
    Coeff_vect_ind[9] = 1;
    // Azz
    Coeff_vect_ind[10] = 1;
    // Bx
    Coeff_vect_ind[11] = 1;
    // By
    Coeff_vect_ind[12] = 1;
    // Bz
    Coeff_vect_ind[13] = 1;
    // Tx
    Coeff_vect_ind[14] = 1;
    // Ty
    Coeff_vect_ind[15] = 1;
    // Tz
    Coeff_vect_ind[16] = 1;
    // Cval
    Coeff_vect_ind[17] = 1;
    // Lval - possibly for time integration (may be included in Sval)
    Coeff_vect_ind[18] = 0;
    // Qx
    Coeff_vect_ind[19] = 1;
    // Qy
    Coeff_vect_ind[20] = 1;
    // Qz
    Coeff_vect_ind[21] = 1;
    // Sval
    Coeff_vect_ind[22] = 1;

  }
  // pure convection problems not implemented yet...
  /* else if(name>=100&&name<200){ */

  /*   // Axx */
  /*   Coeff_vect_ind[2] = 1; */
  /*   // Ayy */
  /*   Coeff_vect_ind[6] = 1; */
  /*   // Azz */
  /*   Coeff_vect_ind[10] = 1; */
  /*   // Bx */
  /*   Coeff_vect_ind[11] = 1; */
  /*   // By */
  /*   Coeff_vect_ind[12] = 1; */
  /*   // Bz */
  /*   Coeff_vect_ind[13] = 1; */
  /*   // Cval */
  /*   Coeff_vect_ind[17] = 1; */
  /*   // Sval */
  /*   Coeff_vect_ind[22] = 1; */

  /* } */
  /* else if(name>=200&&name<1000){ */

  /*   // Axx */
  /*   Coeff_vect_ind[2] = 1; */
  /*   // Axy */
  /*   Coeff_vect_ind[3] = 1; */
  /*   // Axz */
  /*   Coeff_vect_ind[4] = 1; */
  /*   // Ayx */
  /*   Coeff_vect_ind[5] = 1; */
  /*   // Ayy */
  /*   Coeff_vect_ind[6] = 1; */
  /*   // Ayz */
  /*   Coeff_vect_ind[7] = 1; */
  /*   // Azx */
  /*   Coeff_vect_ind[8] = 1; */
  /*   // Azy */
  /*   Coeff_vect_ind[9] = 1; */
  /*   // Azz */
  /*   Coeff_vect_ind[10] = 1; */
  /*   // Bx */
  /*   Coeff_vect_ind[11] = 1; */
  /*   // By */
  /*   Coeff_vect_ind[12] = 1; */
  /*   // Bz */
  /*   Coeff_vect_ind[13] = 1; */
  /*   // Cval */
  /*   Coeff_vect_ind[17] = 1; */
  /*   // Sval */
  /*   Coeff_vect_ind[22] = 1; */
    
  /* } */
  else{
    
    // Mval - possibly mass matrix in time integration (may be included in Cval)
    Coeff_vect_ind[1] = 0;
    // Axx
    Coeff_vect_ind[2] = 1;
    // Axy
    Coeff_vect_ind[3] = 1;
    // Axz
    Coeff_vect_ind[4] = 1;
    // Ayx
    Coeff_vect_ind[5] = 1;
    // Ayy
    Coeff_vect_ind[6] = 1;
    // Ayz
    Coeff_vect_ind[7] = 1;
    // Azx
    Coeff_vect_ind[8] = 1;
    // Azy
    Coeff_vect_ind[9] = 1;
    // Azz
    Coeff_vect_ind[10] = 1;
    // Bx
    Coeff_vect_ind[11] = 1;
    // By
    Coeff_vect_ind[12] = 1;
    // Bz
    Coeff_vect_ind[13] = 1;
    // Tx
    Coeff_vect_ind[14] = 1;
    // Ty
    Coeff_vect_ind[15] = 1;
    // Tz
    Coeff_vect_ind[16] = 1;
    // Cval
    Coeff_vect_ind[17] = 1;
    // Lval - possibly for time integration (may be included in Sval)
    Coeff_vect_ind[18] = 0;
    // Qx
    Coeff_vect_ind[19] = 1;
    // Qy
    Coeff_vect_ind[20] = 1;
    // Qz
    Coeff_vect_ind[21] = 1;
    // Sval
    Coeff_vect_ind[22] = 1;

  }

  return(1);
}

/*!!!!!! OLD OBSOLETE VERSION !!!!!!*/
/*------------------------------------------------------------
  pdr_conv_diff_select_el_coeff - to select coefficients returned to approximation
                        routines for element integrals
------------------------------------------------------------*/
double* pdr_conv_diff_select_el_coeff( 
			  /* returns: pointer !=NULL to indicate selection */
  int Problem_id,
  double **Mval,	/* out: mass matrix coefficient */
  double **Axx,double **Axy,double **Axz, /* out:diffusion coefficients, e.g.*/
  double **Ayx,double **Ayy,double **Ayz, /* Axy denotes scalar or matrix */
  double **Azx,double **Azy,double **Azz, /* related to terms with dv/dx*du/dy */
  /* second order derivatives in weak formulation (scalar for scalar problems */
  /* matrix for vector problems) */
  /* WARNING: if axy==NULL only diagonal (axx, ayy, azz) terms are considered */
  /* in apr_num_int_el */
  double **Bx,double **By,double **Bz,	/* out: convection coefficients */
  /* Bx denotes scalar or matrix related to terms with du/dx*v in weak form */
  double **Tx,double **Ty,double **Tz,	/* out: convection coefficients */
  /* Tx denotes scalar or matrix related to terms with u*dv/dx in weak form */
  double **Cval,/* out: reaction coefficients - for terms without derivatives */
  /*  in weak form (as usual: scalar for scalar problems, matrix for vectors) */
  double **Lval,/* out: rhs coefficient for time term, Lval denotes scalar */
  /* or matrix corresponding to time derivative - similar as mass matrix but  */
  /* with known solution at the previous time step (usually denoted by u_n) */
  double **Qx,/* out: rhs coefficients for terms with derivatives */
  double **Qy,/* Qy denotes scalar or matrix corresponding to terms with dv/dy */
  double **Qz,/* derivatives in weak formulation */
  double **Sval	/* out: rhs coefficients without derivatives (source terms) */
  )
{

/*++++++++++++++++ executable statements ++++++++++++++++*/

  int nreq = pdr_ctrl_i_params(Problem_id,5);

  // We allocate space for all possible coefficients!
  // (in versions for specific problems only necessary terms should be allocated)
  if(*Mval!=NULL) free(*Mval);
  *Mval = (double *) malloc(nreq*nreq*sizeof(double));
  if(*Axx!=NULL) free(*Axx);
  *Axx = (double *) malloc(nreq*nreq*sizeof(double));
  if(*Axy!=NULL) free(*Axy);
  *Axy = (double *) malloc(nreq*nreq*sizeof(double));
  if(*Axz!=NULL) free(*Axz);
  *Axz = (double *) malloc(nreq*nreq*sizeof(double));
  if(*Ayx!=NULL) free(*Ayx);
  *Ayx = (double *) malloc(nreq*nreq*sizeof(double));
  if(*Ayy!=NULL) free(*Ayy);
  *Ayy = (double *) malloc(nreq*nreq*sizeof(double));
  if(*Ayz!=NULL) free(*Ayz);
  *Ayz = (double *) malloc(nreq*nreq*sizeof(double));
  if(*Azx!=NULL) free(*Azx);
  *Azx = (double *) malloc(nreq*nreq*sizeof(double));
  if(*Azy!=NULL) free(*Azy);
  *Azy = (double *) malloc(nreq*nreq*sizeof(double));
  if(*Azz!=NULL) free(*Azz);
  *Azz = (double *) malloc(nreq*nreq*sizeof(double));
  if(*Bx!=NULL) free(*Bx);
  *Bx = (double *) malloc(nreq*nreq*sizeof(double));
  if(*By!=NULL) free(*By);
  *By = (double *) malloc(nreq*nreq*sizeof(double));
  if(*Bz!=NULL) free(*Bz);
  *Bz = (double *) malloc(nreq*nreq*sizeof(double));
  if(*Tx!=NULL) free(*Tx);
  *Tx = (double *) malloc(nreq*nreq*sizeof(double));
  if(*Ty!=NULL) free(*Ty);
  *Ty = (double *) malloc(nreq*nreq*sizeof(double));
  if(*Tz!=NULL) free(*Tz);
  *Tz = (double *) malloc(nreq*nreq*sizeof(double));
  if(*Cval!=NULL) free(*Cval); 
  *Cval = (double *) malloc(nreq*nreq*sizeof(double));
  if(*Lval!=NULL) free(*Lval);
  *Lval = (double *) malloc(nreq*sizeof(double));
  if(*Qx!=NULL) free(*Qx);
  *Qx = (double *) malloc(nreq*sizeof(double));
  if(*Qy!=NULL) free(*Qy);
  *Qy = (double *) malloc(nreq*sizeof(double));
  if(*Qz!=NULL) free(*Qz);
  *Qz = (double *) malloc(nreq*sizeof(double));
  if(*Sval!=NULL) free(*Sval);
  *Sval = (double *) malloc(nreq*sizeof(double));

  return (*Axx);
}


/*---------------------------------------------------------
  pdr_conv_diff_el_coeff - to return coefficients for element integrals
----------------------------------------------------------*/
int pdr_conv_diff_el_coeff(
  /* GENERIC arguments as in pdr_el_coeff */
  int Problem_id,
  int Elem,	/* in: element number */
  int Mat_num,	/* in: material number */
  double Hsize,	/* in: size of an element */
  int Pdeg,	/* in: local degree of polynomial */
  double *X_loc,      /* in: local coordinates of point within element */
  double *Base_phi,   /* in: basis functions */
  double *Base_dphix, /* in: x-derivatives of basis functions */
  double *Base_dphiy, /* in: y-derivatives of basis functions */
  double *Base_dphiz, /* in: z-derivatives of basis functions */
  double *Xcoor,	/* in: global coordinates of a point */
  double *Uk_val, 	/* in: computed solution from previous iteration */
  double *Uk_x, 	/* in: gradient of computed solution Uk_val */
  double *Uk_y,   	/* in: gradient of computed solution Uk_val */
  double *Uk_z,   	/* in: gradient of computed solution Uk_val */
  double *Un_val, 	/* in: computed solution from previous time step */
  double *Un_x, 	/* in: gradient of computed solution Un_val */
  double *Un_y,   	/* in: gradient of computed solution Un_val */
  double *Un_z,   	/* in: gradient of computed solution Un_val */
  /* OBSOLETE!!!!!! */double *Mval,	/* out: mass matrix coefficient OBSOLETE!!!!!! */
  double *Axx, double *Axy, double *Axz,  /* out:diffusion coefficients */
  double *Ayx, double *Ayy, double *Ayz,  /* e.g. Axy denotes scalar or matrix */
  double *Azx, double *Azy, double *Azz,  /* related to terms with dv/dx*du/dy */
  /* second order derivatives in weak formulation (scalar for scalar problems */
  /* matrix for vector problems) */
  /* WARNING: if axy==NULL only diagonal (axx, ayy, azz) terms are considered */
  /* in apr_num_int_el */
  double *Bx, double *By, double *Bz,	/* out: convection coefficients */
  /* Bx denotes scalar or matrix related to terms with du/dx*v in weak form */
  double *Tx, double *Ty, double *Tz,	/* out: convection coefficients */
  /* Tx denotes scalar or matrix related to terms with u*dv/dx in weak form */
  double *Cval,	/* out: reaction coefficients - for terms without derivatives */
  /*  in weak form (as usual: scalar for scalar problems, matrix for vectors) */
  /* OBSOLETE!!!!!! */double *Lval, /* out: rhs coefficient for time term, Lval denotes scalar*/
  /* or matrix corresponding to time derivative - similar as mass matrix but  */
  /* with known solution at the previous time step (usually denoted by u_n) */
  double *Qx, /* out: rhs coefficients for terms with derivatives */
  double *Qy, /* Qy denotes scalar or matrix corresponding to terms with dv/dy */
  double *Qz, /* derivatives in weak formulation */
  double *Sval	/* out: rhs coefficients without derivatives (source terms) */
)
{

/*++++++++++++++++ executable statements ++++++++++++++++*/

  double time=pdr_time_d_params(PDC_USE_CURRENT_PROBLEM_ID, 4);

/* when not using any additional modifications (e.g. artificial viscosity) */
/* get directly coefficients of the corresponding pde */
  pdr_pde_coeff(Mat_num,Xcoor,time,Un_val,Un_x,Un_y,Un_z,
			Axx,Axy,Axz,Ayx,Ayy,Ayz,Azx,Azy,Azz,Bx,By,Bz,
			Tx,Ty,Tz,Cval,Mval,Lval,Sval,Qx,Qy,Qz);

  //printf("\nelem %d, xcoor %lf, %lf, %lf sval %lf\n", Elem, Xcoor[0], Xcoor[1], Xcoor[2], Sval[0]);
 
  //// pure convection problems not implemented yet....
  //double implicit = pdr_time_d_params(PDC_USE_CURRENT_PROBLEM_ID, 2);
  ////int name = pdr_ctrl_i_params(PDC_USE_CURRENT_PROBLEM_ID, 1); - changed for string
  //char name[300];
  //pdr_problem_name(Problem_id, name);

  /* if(name>=100 && name<1000){ */

  /*   //printf("Implicitness parameter %lf\n",implicit); */
  /*   if(implicit < 1.0){ */
      
  /*     int nreq = pdr_ctrl_i_params(Problem_id, 5); */
  /*     int weq, ueq, i; */
      
  /*     for(weq=0;weq<nreq;weq++){ */
	
  /* 	for(ueq=0;ueq<nreq;ueq++){ */
  /* 	  Qx[weq] += (implicit-1.0) * ( */
  /* 				       Axx[weq+nreq*ueq]  * Un_x[ueq] */
  /* 				       + Axy[weq+nreq*ueq]  * Un_y[ueq] */
  /* 				       + Axz[weq+nreq*ueq]  * Un_z[ueq] */
  /* 				       + Tx[weq+nreq*ueq]   * Un_val[ueq] */
  /* 				       ); */
  /* 	  Qy[weq] += (implicit-1.0) * ( */
  /* 				       Ayx[weq+nreq*ueq]  * Un_x[ueq] */
  /* 				       + Ayy[weq+nreq*ueq]  * Un_y[ueq] */
  /* 				       + Ayz[weq+nreq*ueq]  * Un_z[ueq] */
  /* 				       + Ty[weq+nreq*ueq]   * Un_val[ueq] */
  /* 				       ); */
  /* 	  Qz[weq] += (implicit-1.0) * ( */
  /* 				       Azx[weq+nreq*ueq]  * Un_x[ueq] */
  /* 				       + Azy[weq+nreq*ueq]  * Un_y[ueq] */
  /* 				       + Azz[weq+nreq*ueq]  * Un_z[ueq] */
  /* 				       + Tz[weq+nreq*ueq]   * Un_val[ueq] */
  /* 				       ); */
  /* 	  Sval[weq] += (implicit-1.0) * ( */
  /* 	                        	 Cval[weq+nreq*ueq] * Un_val[ueq] */
  /* 					 + Bx[weq+nreq*ueq] * Un_x[ueq] */
  /* 					 + By[weq+nreq*ueq] * Un_y[ueq] */
  /* 					 + Bz[weq+nreq*ueq] * Un_z[ueq] */
  /* 					 ); */
  /* 	} */
	
  /*     } */

  /*     for(weq=0;weq<nreq;weq++){ */
	
  /* 	for(ueq=0;ueq<nreq;ueq++){ */
	  
  /* 	  i = weq+nreq*ueq; */
	  
  /* 	  Axx[i] *= implicit; */
  /* 	  Axy[i] *= implicit; */
  /* 	  Axz[i] *= implicit; */
  /* 	  Ayx[i] *= implicit; */
  /* 	  Ayy[i] *= implicit; */
  /* 	  Ayz[i] *= implicit; */
  /* 	  Azx[i] *= implicit; */
  /* 	  Azy[i] *= implicit; */
  /* 	  Azz[i] *= implicit; */
  /* 	  Bx[i] *= implicit; */
  /* 	  By[i] *= implicit; */
  /* 	  Bz[i] *= implicit; */
  /* 	  Tx[i] *= implicit; */
  /* 	  Ty[i] *= implicit; */
  /* 	  Tz[i] *= implicit; */
  /* 	Cval[i] *= implicit; */
  /* 	} */
  /*     } */
      
  /*   } */

  /* } */

  return(1);
}


/*------------------------------------------------------------
  pdr_conv_diff_comp_el_stiff_mat_app - to construct a stiffness matrix and  
                          a load vector for an element (approximation specific version)
------------------------------------------------------------*/
extern int pdr_conv_diff_comp_el_stiff_mat_app(/*returns: >=0 -success code, <0 -error code*/
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
				       );


/*------------------------------------------------------------
  pdr_conv_diff_comp_el_stiff_mat - to construct a stiffness matrix and  
                          a load vector for an element
------------------------------------------------------------*/
int pdr_conv_diff_comp_el_stiff_mat(/*returns: >=0 -success code, <0 -error code*/
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


    pdr_conv_diff_comp_el_stiff_mat_app(Problem_id, El_id, Comp_sm, Pdeg_in, 
				  Nr_dof_ent, List_dof_ent_type, 
				  List_dof_ent_id, List_dof_ent_nrdofs, 
				  Nrdofs_loc, Stiff_mat, Rhs_vect, Rewr_dofs);

    /* if(Rhs_vect != NULL){ */
    /*   printf("ielem %d, rhs_vect: %lf, %lf, %lf, %lf, %lf, %lf\n", El_id, */
    /* 	     Rhs_vect[0],Rhs_vect[1],Rhs_vect[2],Rhs_vect[3],Rhs_vect[4],Rhs_vect[5]); */
    /* } */

  return (1);

}




/*---------------------------------------------------------
pdr_fa_coeff - to return coefficients for face integrals
	in DGFEM for convection-diffusion-reaction-etc equations
----------------------------------------------------------*/
void pdr_fa_coeff(
	int Face,	/* in: face number */
        int BC_flag,	/* in: boundary condition flag */
	int Elem,	/* in: element number */
	int Mat_num,	/* in: material number */
	double Hsize,	/* in: size of an element */
	double Time,    /* in: time instant */
	int Pdeg,	/* in: local degree of polynomial */
	double* Xcoor,	/* in: global coordinates of a point */
	double* Vec_norm, /* in: unit normal to the boundary */
	double* U_val, 	/* in: computed solution */
	double* U_x, 	/* in: gradient of computed solution */
	double* U_y,   	/* in: gradient of computed solution */
	double* U_z,   	/* in: gradient of computed solution */
	double* Anx,	/* out: diffusion coefficients */
	double* Any, 	
	double* Anz, 	
	double* Bn,	/* out: convection coefficients */
	double* Fval,	/* out: rhs  Dirichlet data f */
	double* Gval,	/* out: rhs  Neumann data g */
	double* Qn,	/* out: rhs normal fluxes */
	double* Vel_norm /* out: velocity normal to the boundary */
	)
{

  int i,iaux, ieq1, ieq2;
  double x, y, z, f, f_x, f_y, f_z, daux;
  double axx[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];  /* coeff of PDE */
  double axy[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];  /* coeff of PDE */
  double axz[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];  /* coeff of PDE */
  double ayx[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];  /* coeff of PDE */
  double ayy[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];  /* coeff of PDE */
  double ayz[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];  /* coeff of PDE */
  double azx[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];  /* coeff of PDE */
  double azy[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];  /* coeff of PDE */
  double azz[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];  /* coeff of PDE */
  double bx[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];   /* coeff of PDE */
  double by[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];   /* coeff of PDE */
  double bz[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];   /* coeff of PDE */
  double tx[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];   /* coeff of PDE */
  double ty[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];   /* coeff of PDE */
  double tz[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ];   /* coeff of PDE */
  double cval[PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ]; /* coeff of PDE */
  double mval[PDC_CONV_DIFF_NREQ];            /* coeff of PDE */
  double lval[PDC_CONV_DIFF_NREQ];           /* coeff of PDE */
  double sval[PDC_CONV_DIFF_NREQ];            /* coeff of PDE */
  double qx[PDC_CONV_DIFF_NREQ];              /* coeff of PDE */
  double qy[PDC_CONV_DIFF_NREQ];              /* coeff of PDE */
  double qz[PDC_CONV_DIFF_NREQ];              /* coeff of PDE */

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* get coefficients of the corresponding pde */
  pdr_pde_coeff(Mat_num,Xcoor,Time,U_val,U_x,U_y,U_z,
			axx,axy,axz,ayx,ayy,ayz,azx,azy,azz,bx,by,bz,
			tx,ty,tz,cval,mval,lval,sval,qx,qy,qz);

/*kbw
  if(Face>0){
    printf("diffusion coeff : %lf %lf %lf\n",axx[0],axy[0],axz[0]);
    printf("diffusion coeff : %lf %lf %lf\n",ayx[0],ayy[0],ayz[0]);
    printf("diffusion coeff : %lf %lf %lf\n",azx[0],azy[0],azz[0]);
    printf("convection coeff: %lf %lf %lf\n",bx[0],by[0],bz[0]);
    printf("convection coeff: %lf %lf %lf\n",tx[0],ty[0],tz[0]);
    printf("reaction coeff  : %lf\n",cval[0]);
    printf("time coeff LHS  : %lf RHS %lf, source %lf\n",
	   mval[0],lval[0], sval[0]);
    printf("s_xn coeff: %lf %lf %lf\n",qx[0],qy[0],qz[0]);
  }
/*kew*/

/* initialize coefficients when necessary */
  utr_d_zero(Qn,PDC_CONV_DIFF_NREQ);
  utr_d_zero(Gval,PDC_CONV_DIFF_NREQ);
  utr_d_zero(Fval,PDC_CONV_DIFF_NREQ);

/* internal faces */
  if(BC_flag==PDC_INTERIOR){

/* compute normal fluxes */
    for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){
      for(ieq2=0;ieq2<PDC_CONV_DIFF_NREQ;ieq2++){
	Anx[ieq1*PDC_CONV_DIFF_NREQ+ieq2] =   axx[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[0] + 
	                        ayx[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[1] +
	                        azx[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[2] ;
	Any[ieq1*PDC_CONV_DIFF_NREQ+ieq2] =   axy[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[0] + 
	                        ayy[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[1] +
	                        azy[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[2] ;
	Anz[ieq1*PDC_CONV_DIFF_NREQ+ieq2] =   axz[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[0] + 
	                        ayz[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[1] +
	                        azz[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[2] ;
	Bn[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  =   bx[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[0] + 
	                        by[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[1] +
	                        bz[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[2] ;
      }
    }
  }
  else if(BC_flag==PDC_BC_DIRI){
/* Dirichlet faces */

    //printf("Dirichlet");

/* compute normal fluxes */
    for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){
      for(ieq2=0;ieq2<PDC_CONV_DIFF_NREQ;ieq2++){
	Anx[ieq1*PDC_CONV_DIFF_NREQ+ieq2] =   axx[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[0] + 
	                        ayx[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[1] +
	                        azx[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[2] ;
	Any[ieq1*PDC_CONV_DIFF_NREQ+ieq2] =   axy[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[0] + 
	                        ayy[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[1] +
	                        azy[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[2] ;
	Anz[ieq1*PDC_CONV_DIFF_NREQ+ieq2] =   axz[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[0] + 
	                        ayz[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[1] +
	                        azz[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[2] ;
	Bn[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  =   bx[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[0] + 
	                        by[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[1] +
	                        bz[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[2] ;
      }
    }

/* get boundary coefficient of the corresponding pde */
    pdr_bc_diri_coeff(Face,Mat_num,Xcoor,Vec_norm,Time,
		       U_val,U_x,U_y,U_z, 
		       Fval);


  }
  else if(BC_flag==PDC_BC_NEUM){
/* Neumann faces */

    //printf("Neumann");

/* compute normal fluxes */
    for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++){
      for(ieq2=0;ieq2<PDC_CONV_DIFF_NREQ;ieq2++){
	Bn[ieq1*PDC_CONV_DIFF_NREQ+ieq2]  =   bx[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[0] + 
	                        by[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[1] +
	                        bz[ieq1*PDC_CONV_DIFF_NREQ+ieq2] * Vec_norm[2] ;
      }
    }

/* get boundary coefficient of the corresponding pde */
    pdr_bc_neum_coeff(Face,Mat_num,Xcoor,Vec_norm,Time,
		       U_val,U_x,U_y,U_z,
		       Gval);

  }
  else if(BC_flag==PDC_BC_MIXED){
/* Robin faces */

/* get boundary coefficient of the corresponding pde */
    pdr_bc_mixed_coeff(Face,Mat_num,Xcoor,Vec_norm,Time,
		       U_val,U_x,U_y,U_z, 
		       Fval, Vel_norm);

  }

  for(ieq1=0;ieq1<PDC_CONV_DIFF_NREQ;ieq1++) Vel_norm[ieq1]=Bn[ieq1*PDC_CONV_DIFF_NREQ+ieq1];

/*kbw
    printf("returning from pdr_fa_coeff for face %d (type %d, neig %d, mate %d)\n",
	   Face, BC_flag, Elem, Mat_num);
    printf("coor global %lf, %lf, %lf, area %lf\n",
	   Xcoor[0],Xcoor[1],Xcoor[2], Hsize);
    printf("normal %lf, %lf, %lf\n", 
	   Vec_norm[0], Vec_norm[1], Vec_norm[2]);
    printf("solution %lf, der %lf, %lf, %lf\n", 
	   U_val[0],U_x[0],U_y[0],U_z[0]);
    printf("diffusion coeff : %lf %lf %lf\n",Anx[0],Any[0],Anz[0]);
    printf("convection coeff: %lf\n",Bn[0]);
    printf("rhs coeff       : %lf %lf %lf\n",Fval[0],Gval[0],Qn[0]);
/*kew*/ 

  return;
}

/*---------------------------------------------------------
pdr_pde_coeff - to return coefficients for internal integrals
	in DGFEM for convection-diffusion-reaction-etc equations
----------------------------------------------------------*/
void pdr_pde_coeff(
	int Mat_num,	/* in: material number */
	double* Xcoor,	/* in: global coordinates of a point */
	double  Time,	/* in: time instant */
	double* U_val, 	/* in: computed solution */
	double* U_x, 	/* in: gradient of computed solution */
	double* U_y,   	/* in: gradient of computed solution */
	double* U_z,   	/* in: gradient of computed solution */
	double* Axx, 	/* out: diffusion coefficients */		
	double* Axy,		
	double* Axz,		
	double* Ayx,
	double* Ayy,
	double* Ayz,
	double* Azx,
	double* Azy,
	double* Azz,
	double* Bx,	/* out: convection coefficients */	
	double* By,	
	double* Bz,	
	double* Tx,	/* out: convection coefficients */	
	double* Ty,	
	double* Tz,	
	double* Cval,	/* out: reaction coefficients */
	double* Mval,	/* out: mass matrix coefficients */
	double* Lval,	/* out: rhs time coefficients */
	double* Sval,	/* out: rhs coefficients without derivatives */	
	double* Qx,	/* out: rhs coefficients with derivatives */	
	double* Qy,	/* out: rhs coefficients with derivatives */	
	double* Qz	/* out: rhs coefficients with derivatives */	
	)
{

  int i, iaux;
  double lapl,daux,cd1d_coeff,time, dtime;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* initialize coefficients */
  if(Mval!=NULL) utr_d_zero(Mval,PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ);
  if(Axx!=NULL) utr_d_zero(Axx,PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ);
  if(Axy!=NULL) utr_d_zero(Axy,PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ);
  if(Axz!=NULL) utr_d_zero(Axz,PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ);
  if(Ayx!=NULL) utr_d_zero(Ayx,PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ);
  if(Ayy!=NULL) utr_d_zero(Ayy,PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ);
  if(Ayz!=NULL) utr_d_zero(Ayz,PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ);
  if(Azx!=NULL) utr_d_zero(Azx,PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ);
  if(Azy!=NULL) utr_d_zero(Azy,PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ);
  if(Azz!=NULL) utr_d_zero(Azz,PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ);
  if(Bx!=NULL) utr_d_zero(Bx,PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ);
  if(By!=NULL) utr_d_zero(By,PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ);
  if(Bz!=NULL) utr_d_zero(Bz,PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ);
  if(Tx!=NULL) utr_d_zero(Tx,PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ);
  if(Ty!=NULL) utr_d_zero(Ty,PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ);
  if(Tz!=NULL) utr_d_zero(Tz,PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ);
  if(Cval!=NULL) utr_d_zero(Cval,PDC_CONV_DIFF_NREQ*PDC_CONV_DIFF_NREQ);
  if(Lval!=NULL) utr_d_zero(Lval,PDC_CONV_DIFF_NREQ);
  if(Qx!=NULL) utr_d_zero(Qx,PDC_CONV_DIFF_NREQ);
  if(Qy!=NULL) utr_d_zero(Qy,PDC_CONV_DIFF_NREQ);
  if(Qz!=NULL) utr_d_zero(Qz,PDC_CONV_DIFF_NREQ);
  if(Sval!=NULL) utr_d_zero(Sval,PDC_CONV_DIFF_NREQ);

/* rewrite and check parameters*/
  //name = pdr_ctrl_i_params(PDC_USE_CURRENT_PROBLEM_ID, 1);
  char name[300];
  pdr_problem_name(PDC_USE_CURRENT_PROBLEM_ID, name);

  //printf("Problem name: %s\n", name);

  //if(name==1||name==2){
  if(strcmp(name, "DIFF_IN_CUBE") == 0 || strcmp(name, "DIFF_IN_BOX") == 0){

    /* simple scalar elliptic problem with known exact solution */
    Axx[0]=1;	
    Ayy[0]=1;	
    Azz[0]=1;	

    if(Xcoor==NULL) return;

    double x = Xcoor[0];
    double y = Xcoor[1];
    double z = Xcoor[2];

    /* right hand side corresponding to exact solution of test problem */
    /* u = exp(-(x^2+y^2+z^2)) */
    pdr_exact_sol(Mat_num,x,y,z,Time,&daux,&daux,&daux,&daux,&lapl);

    Sval[0]= -lapl;

    //printf("Problem name: %s, lapl %lf\n", name, lapl);


  }
  // pure convection problems not implemented yet....
  /* else if(name==101){ */

  /*   /\* simple 1D scalar convection-diffusion problem *\/ */
  /*   dtime = pdr_time_d_params(PDC_USE_CURRENT_PROBLEM_ID, 6); */

  /*   //!!!!!!!!!! OLD OBSOLETE INTERFACE !!!!!!!!!!!! */
  /*   //Mval[0]= 1.0/dtime ; */
  /*   //Lval[0]= U_val[0]/dtime ; */
  /*   // new interface     */
  /*   Cval[0]= 1.0/dtime ; */
  /*   Sval[0]= U_val[0]/dtime ; */

  /*   // hard coded PDE coefficient */
  /*   cd1d_coeff=1.0; */

  /*   //printf("diffusion coeff : %lf\n",cd1d_coeff); */

  /*   Axx[0] = cd1d_coeff; */
  /*   Ayy[0] = cd1d_coeff; */
  /*   Azz[0] = cd1d_coeff; */

  /*   Bx[0]=1.0; */
    
  /* } */
  /* else if(name==105){ */

  /*   /\* simple 1D scalar convection-diffusion problem *\/ */
  /*   dtime = pdr_time_d_params(PDC_USE_CURRENT_PROBLEM_ID, 6); */

  /*   //!!!!!!!!!! OLD OBSOLETE INTERFACE !!!!!!!!!!!! */
  /*   //Mval[0]= 1.0/dtime ; */
  /*   //Lval[0]= U_val[0]/dtime ; */
  /*   // new interface     */
  /*   Cval[0]= 1.0/dtime ; */
  /*   Sval[0]= U_val[0]/dtime ; */

  /*   // hard coded PDE coefficient */
  /*   cd1d_coeff=1.0; */

  /*   //printf("diffusion coeff : %lf\n",cd1d_coeff); */
    
  /*   Axx[0] = cd1d_coeff; */
  /*   Ayy[0] = cd1d_coeff; */
  /*   Azz[0] = cd1d_coeff; */
    
  /*   Bx[0]=1.0/3.0; */
  /*   By[0]=1.0/3.0; */
  /*   Bz[0]=1.0/3.0; */

    
  /* } */
  /* else if(name==201){ */

  /*   /\* simple 1D scalar convection-diffusion problem *\/ */
  /*   dtime = pdr_time_d_params(PDC_USE_CURRENT_PROBLEM_ID, 6); */

  /*   //!!!!!!!!!! OLD OBSOLETE INTERFACE !!!!!!!!!!!! */
  /*   //Mval[0]= 1.0/dtime ; */
  /*   //Lval[0]= U_val[0]/dtime ; */
  /*   // new interface     */
  /*   Cval[0]= 1.0/dtime ; */
  /*   Sval[0]= U_val[0]/dtime ; */

  /*   // hard coded PDE coefficient */
  /*   cd1d_coeff=1.0; */

  /*   //printf("diffusion coeff before TG: %lf\n",cd1d_coeff); */

  /*   // Taylor-Galerkin ? */
  /*   cd1d_coeff += dtime/2.0; */
    
  /*   //printf("diffusion coeff : %lf\n",cd1d_coeff); */

  /*   Bx[0]=1.0; */
    
  /*   Axx[0] = cd1d_coeff + dtime/2.0*Bx[0]*Bx[0]; */
  /*   Axy[0] = dtime/2.0*Bx[0]*By[0]; */
  /*   Axz[0] = dtime/2.0*Bx[0]*Bz[0]; */
  /*   Ayx[0] = dtime/2.0*By[0]*Bx[0]; */
  /*   Ayy[0] = cd1d_coeff + dtime/2.0*By[0]*By[0]; */
  /*   Ayz[0] = dtime/2.0*By[0]*Bz[0]; */
  /*   Azx[0] = dtime/2.0*Bz[0]*Bx[0]; */
  /*   Azy[0] = dtime/2.0*Bz[0]*By[0]; */
  /*   Azz[0] = cd1d_coeff + dtime/2.0*Bz[0]*Bz[0]; */

  /* } */
  /* else if(name==205){ */

  /*   /\* simple 1D scalar convection-diffusion problem *\/ */
  /*   dtime = pdr_time_d_params(PDC_USE_CURRENT_PROBLEM_ID, 6); */

  /*   //!!!!!!!!!! OLD OBSOLETE INTERFACE !!!!!!!!!!!! */
  /*   //Mval[0]= 1.0/dtime ; */
  /*   //Lval[0]= U_val[0]/dtime ; */
  /*   // new interface     */
  /*   Cval[0]= 1.0/dtime ; */
  /*   Sval[0]= U_val[0]/dtime ; */

  /*   // hard coded PDE coefficient */
  /*   cd1d_coeff=1.0; */

  /*   //printf("diffusion coeff : %lf\n",cd1d_coeff); */
    
  /*   Bx[0]=1.0/3.0; */
  /*   By[0]=1.0/3.0; */
  /*   Bz[0]=1.0/3.0; */

  /*   Axx[0] = cd1d_coeff + dtime/2.0*Bx[0]*Bx[0]; */
  /*   Axy[0] = dtime/2.0*Bx[0]*By[0]; */
  /*   Axz[0] = dtime/2.0*Bx[0]*Bz[0]; */
  /*   Ayx[0] = dtime/2.0*By[0]*Bx[0]; */
  /*   Ayy[0] = cd1d_coeff + dtime/2.0*By[0]*By[0]; */
  /*   Ayz[0] = dtime/2.0*By[0]*Bz[0]; */
  /*   Azx[0] = dtime/2.0*Bz[0]*Bx[0]; */
  /*   Azy[0] = dtime/2.0*Bz[0]*By[0]; */
  /*   Azz[0] = cd1d_coeff + dtime/2.0*Bz[0]*Bz[0]; */
    
    
  /* } */
  //else if(name==1234){
  else if(strcmp(name, "TEST_NUM_INT") == 0){

    /* simple 1D scalar convection-diffusion problem with full set of coeff matrices*/
    //dtime = pdr_time_d_params(PDC_USE_CURRENT_PROBLEM_ID, 6);
    dtime = 0.1;

    //!!!!!!!!!! OLD OBSOLETE INTERFACE !!!!!!!!!!!!
    //Mval[0]= 1.0/dtime ;
    //Lval[0]= U_val[0]/dtime ;
    // new interface    
    Cval[0]= 1.0/dtime ;
    //Sval[0]= U_val[0]/dtime ;
    Sval[0]= 2.0/dtime ;

    Qx[0] = 1.0;
    Qy[0] = 2.0;
    Qz[0] = 3.0;

    // hard coded PDE coefficient
    cd1d_coeff=1.0;
  
    Bx[0]=1.0/3.0;
    By[0]=2.0/3.0;
    Bz[0]=3.0/3.0;

    Tx[0]=3.0/3.0;
    Ty[0]=2.0/3.0;
    Tz[0]=1.0/3.0;

    Axx[0] = cd1d_coeff + dtime/2.0*Bx[0]*Bx[0];
    Axy[0] = dtime/2.0*Bx[0]*By[0];
    Axz[0] = dtime/2.0*Bx[0]*Bz[0];
    Ayx[0] = dtime/2.0*By[0]*Bx[0];
    Ayy[0] = cd1d_coeff + dtime/2.0*By[0]*By[0];
    Ayz[0] = dtime/2.0*By[0]*Bz[0];
    Azx[0] = dtime/2.0*Bz[0]*Bx[0];
    Azy[0] = dtime/2.0*Bz[0]*By[0];
    Azz[0] = cd1d_coeff + dtime/2.0*Bz[0]*Bz[0];
    
    
  }


return;
}

/*---------------------------------------------------------
pdr_bc_diri_coeff - to return DIRICHLET boundary coeficients
                     (routine substitutes only non-zero values)
----------------------------------------------------------*/
void pdr_bc_diri_coeff(
	int Face,       /* in: index in an array of BC coefficients */
	int Mat_num,	/* in: material number */
	double* Xcoor,	/* in: global coordinates of a point */
	double* Vec_norm,/* in: unit normal to the boundary */
	double  Time,	/* in: time instant */
	double* U_val, 	/* in: computed solution */
	double* U_x, 	/* in: gradient of computed solution */
	double* U_y,   	/* in: gradient of computed solution */
	double* U_z,   	/* in: gradient of computed solution */
	double* Fval	/* out: BC DIRICHLET coefficients */	
	)
{

  int i, bc_num;
  double x, y, z, f, f_x, f_y, f_z, daux;

/*++++++++++++++++ executable statements ++++++++++++++++*/


/* rewrite and check parameters*/
  x = Xcoor[0];
  y = Xcoor[1];
  z = Xcoor[2];

/* initialize coefficients */
  utr_d_zero(Fval,PDC_CONV_DIFF_NREQ);

  // name = pdr_ctrl_i_params(PDC_USE_CURRENT_PROBLEM_ID, 1); changed for string
  char name[300];
  pdr_problem_name(PDC_USE_CURRENT_PROBLEM_ID, name);

  if(strcmp(name, "DIFF_IN_CUBE") == 0 || strcmp(name, "DIFF_IN_BOX") == 0 ){

    //if(name==1||name==2){

/* simple scalar elliptic problem with known exact solution */
/* right hand side corresponding to exact solution of test problem */
/* u = exp(-(x^2+y^2+z^2)) */
    pdr_exact_sol(Mat_num,x,y,z,Time,&f,&f_x,&f_y,&f_z,&daux);

    Fval[0] = f;

/*kbw
    printf("Problem name %s, Dirichlet condition value %lf\n",
	   name, Fval[0]); 
/*kew*/

  }
/*   else if(name==101||name==201){ */

/*     Fval[0] = 1; */
/* /\* travelling square *\/ */
/*     if(Time>200) Fval[0]=0; */

/*   } */
/*   else if(name==105||name==205){ */

/* /\* simple 3D scalar convection problem with known exact solution *\/ */
/*     if((Xcoor[0]+Xcoor[1]+Xcoor[2]-Time)<SMALL&& */
/*        (Xcoor[0]+Xcoor[1]+Xcoor[2]-Time)>-0.4) Fval[0]=1; */

/*   } */

  return;
}

/*---------------------------------------------------------
pdr_bc_neum_coeff - to return NEUMANN boundary coeficients
                     (routine substitutes only non-zero values)
----------------------------------------------------------*/
void pdr_bc_neum_coeff(
	int Face,	/* in: face number */
	int Mat_num,	/* in: material number */
	double* Xcoor,	/* in: global coordinates of a point */
	double* Vec_norm,/* in: unit normal to the boundary */
	double  Time,	/* in: time instant */
	double* U_val, 	/* in: computed solution */
	double* U_x, 	/* in: gradient of computed solution */
	double* U_y,   	/* in: gradient of computed solution */
	double* U_z,   	/* in: gradient of computed solution */
	double* Gval	/* out: BC NEUMANN coefficients */	
	)
{

  int i, iaux;
  double x, y, z, f, f_x, f_y, f_z, daux;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  //printf("Neumann condition!\n");

/* rewrite and check parameters*/
  x = Xcoor[0];
  y = Xcoor[1];
  z = Xcoor[2];

/* initialize coefficients */
  utr_d_zero(Gval,PDC_CONV_DIFF_NREQ);

  // name = pdr_ctrl_i_params(PDC_USE_CURRENT_PROBLEM_ID, 1); changed for string
  char name[300];
  pdr_problem_name(PDC_USE_CURRENT_PROBLEM_ID, name);

  if(strcmp(name, "DIFF_IN_CUBE") == 0 || strcmp(name, "DIFF_IN_BOX") == 0 ){

    //if(name==1||name==2||name==101||name==105||name==201||name==205){

/* simple scalar elliptic problem with known exact solution */
/* right hand side corresponding to exact solution of test problem */
/* u = exp(-(x^2+y^2+z^2)) */
    pdr_exact_sol(Mat_num,x,y,z,Time,&f,&f_x,&f_y,&f_z,&daux);
    Gval[0]=f_x*Vec_norm[0]+f_y*Vec_norm[1]+f_z*Vec_norm[2];

/*kbw
    printf("Problem name %s, Neumann condition value %lf\n",
	   name, Gval[0]); 
/*kew*/

  }


  return;
}


/*---------------------------------------------------------
pdr_bc_mixed_coeff - to return ROBIN (mixed) boundary coeficients
----------------------------------------------------------*/
void pdr_bc_mixed_coeff(
	int Face,	/* in: face number */
	int Mat_num,	/* in: material number */
	double* Xcoor,	/* in: global coordinates of a point */
	double* Vec_norm,/* in: unit normal to the boundary */
	double  Time,	/* in: time instant */
	double* U_val, 	/* in: computed solution */
	double* U_x, 	/* in: gradient of computed solution */
	double* U_y,   	/* in: gradient of computed solution */
	double* U_z,   	/* in: gradient of computed solution */
	double* Fval,	/* out: BC mixed coefficients */	
	double* Kr	/* out: proportionality coefficient */	
	)
{

  int i;
  double x,y,z,daux;
/*++++++++++++++++ executable statements ++++++++++++++++*/

/* rewrite and check parameters*/
  x = Xcoor[0];
  y = Xcoor[1];
  z = Xcoor[2];

  return;
}



/*--------------------------------------------------------------------------
 pdr_bc_val - to get the boundary condition values
---------------------------------------------------------------------------*/
int pdr_bc_val( /* returns: 1 - success, <=0 - failure */
	int Problem_id,	/* in: data structure to be used */
	int Num,       	/* in: index in the array of boundary values */
	double *Val     /* out: vector of boundary values */
	)
{
/* auxiliary variables */
  pdt_conv_diff_problem* problem;
  int i,nreq;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  problem = pdr_select_problem(Problem_id);

  if(problem->ctrl.bc_val==NULL || problem->ctrl.bc_val[Num]==NULL){
    return(-1);
  }

  i=5; nreq=pdr_ctrl_i_params(Problem_id,i);
  for(i=0;i<nreq;i++) Val[i]=problem->ctrl.bc_val[Num][i];

  return(1);
}

/*--------------------------------------------------------------------------
 pdr_ic_val - to get the initial condition values
---------------------------------------------------------------------------*/
int pdr_ic_val( /* returns: 1 - success, <=0 - failure */
	int Problem_id,	/* in: data structure to be used */
	int Num,       	/* in: index in the array of initial values */
	double *Val     /* out: vector of initial values */
	)
{
/* auxiliary variables */
  pdt_conv_diff_problem* problem;
  int i,nreq;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  problem = pdr_select_problem(Problem_id);

  if(problem->ctrl.ic_val==NULL || problem->ctrl.ic_val[Num]==NULL){
    return(-1);
  }

  i=5; nreq=pdr_ctrl_i_params(Problem_id,i);
  for(i=0;i<nreq;i++) Val[i]=problem->ctrl.ic_val[Num][i];

  return(1);
}


/*------------------------------------------------------------
  pdr_get_bc_type - to get BC type given BC flag from mesh data structure
                    !!! according to some adopted convention !!!
------------------------------------------------------------*/
int pdr_get_bc_type( /* returns: bc type for a face: */
		/*  	1 (PDC_INTERIOR) - interior face */
		/*  	2 (PDC_BC_DIRI) - Dirichlet boundary face */
		/*  	3 (PDC_BC_NEUM) - Neumann boundary face */
		/*  	4 (PDC_BC_MIXED) - Robin boundary face */
   int Fa_bc	/* in: BC flag (returned by mesh data structure) */
   )
{

/*++++++++++++++++ executable statements ++++++++++++++++*/

  if(Fa_bc<=0) return(PDC_INTERIOR);
  else if(Fa_bc<PDC_MAXBCVAL) return(PDC_BC_DIRI);
  else if(Fa_bc<2*PDC_MAXBCVAL) return(PDC_BC_NEUM);
  else if(Fa_bc<3*PDC_MAXBCVAL) return(PDC_BC_MIXED);
  else {
    printf("Wrong Fa_bc flag %d in pdr_get_bc_type \n", Fa_bc);
    exit(1);
  }

}


/*------------------------------------------------------------
  pdr_conv_diff_comp_fa_stiff_mat_app - to construct a stiffness matrix and  
        a load vector for a face - special routine for different approximation fields
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
				       );



/*------------------------------------------------------------
  pdr_conv_diff_comp_fa_stiff_mat - to construct a stiffness matrix and  
                          a load vector for a face - special routine for discontinuous 
                          Galerkin approximation fields
------------------------------------------------------------*/
int pdr_conv_diff_comp_fa_stiff_mat(/*returns: >=0 -success code, <0 -error code*/
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

    pdr_conv_diff_comp_fa_stiff_mat_app(Problem_id, Fa_id, Comp_sm, Pdeg_in, 
				  Nr_dof_ent, List_dof_ent_type, 
				  List_dof_ent_id, List_dof_ent_nrdofs, 
				  Nrdofs_loc, Stiff_mat, Rhs_vect, Rewr_dofs);

  return (1);

}


/*------------------------------------------------------------
pdr_conv_diff_compute_CFL - to compute global CFL numbers (for a subdomain)
------------------------------------------------------------*/
int pdr_conv_diff_compute_CFL(
  int Problem_id,
  double *CFL_min_p,
  double *CFL_max_p,
  double *CFL_ave_p
)
{


}


