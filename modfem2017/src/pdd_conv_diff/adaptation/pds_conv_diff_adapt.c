/************************************************************************
File contains routines:
  pdr_adapt - to enforce adaptation strategy for a given problem 
  pdr_err_indi - to return error indicator for an element (internal procedure)

------------------------------  			
History:     
	02.2002 - Krzysztof Banas, initial version		
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<string.h>


/* USED DATA STRUCTURES AND INTERFACES FROM OTHER MODULES */

/* utilities - including simple time measurement library */
#include "uth_intf.h"
#include "uth_system.h"

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
#include "../include/pdh_conv_diff.h"	
#include "../include/pdh_conv_diff_problem.h"	


/***************************************/
/* DEFINITIONS OF PROCEDURES */
/***************************************/

/*---------------------------------------------------------
pdr_err_indi_exact - to return error indicator for an element,
	based on the knowledge of the exact solution 
----------------------------------------------------------*/
double pdr_err_indi_exact (	/* returns error indicator for an element */
        int Problem_id,	/* in: data structure to be used  */
	int El		/* in: element number */
);

/*---------------------------------------------------------
pdr_err_indi_ZZ - to return error indicator for an element,
	based on ZZ first derivative recovery 
----------------------------------------------------------*/
double pdr_err_indi_ZZ (	/* returns error indicator for an element */
        int Problem_id,	/* in: data structure to be used  */
	int El		/* in: element number */
);

/*---------------------------------------------------------
pdr_conv_diff_adapt - to enforce adaptation strategy for a given problem 
---------------------------------------------------------*/
int pdr_conv_diff_adapt( /* returns: >0 - success, <=0 - failure */
  int Problem_id,       /* in: problem data structure to be used */
  char* Work_dir,
  FILE *Interactive_input, 
  FILE *Interactive_output
			 );

/*---------------------------------------------------------
pdr_conv_diff_adapt - to enforce adaptation strategy for a given problem 
---------------------------------------------------------*/
int pdr_conv_diff_adapt( /* returns: >0 - success, <=0 - failure */
  int Problem_id,       /* in: problem data structure to be used */
  char* Work_dir,
  FILE *Interactive_input, 
  FILE *Interactive_output
  )
{

  double daux, clock_time;
  int i, iaux, ino, nno, nr_patches, mesh_id, field_id, type;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  i=2; mesh_id=pdr_ctrl_i_params(Problem_id,i);
  i=3; field_id = pdr_ctrl_i_params(Problem_id,i);  

  i=1; type=pdr_adapt_i_params(Problem_id,i);
  // default adaptation strategy:
  // depending on type of adaptation (type=pdr_adapt_i_params(Problem_id,i);)
  // type is specified as the first ADAPTATION_PARAMETERS in problem input file
  // type == -2 - uniform derefinement
  // type == -1 - uniform refinement
  // type == 0 - no adaptations
  // type > 0  - adaptive refinement: element error indicators are computed
  //             by pdr_err_indi, elements with error greater then 
  //             eps*average_error are broken, families of elements (sons
  //             of a single father) with total error less than 
  //             ratio*eps*average_error are clustered back
  // particular values of type may be problem dependent
  // most popular types are encoded in include/pdh_intf.h as:
  //  PDC_ADAPT_EXACT = 1 - adaptations based on the knowledge of exact solution
  //  PDC_ADAPT_ZZ    = 2 - adaptations based on Zienkiewicz-Zhu error estimate
  //  PDC_ADAPT_EXPL  = 3 - adaptations based on explicit residual error estimate
  if(type==PDC_ADAPT_ZZ){

    // create patches for nno_old nodes
    //clock_time = time_clock();
    //nr_patches = utr_create_patches(field_id, &pdv_patches);
    nr_patches = utr_create_patches_small(field_id, &pdv_patches);
    //clock_time = time_clock()-clock_time;
    //printf("New patch creation - time %lf\n", clock_time);
    
    
    // !!! it is assumed that the most recent dofs are in the last vector
    //int sol_vec_id = apr_get_nr_sol(Field_id);
    // !!! it is assumed that the most recent dofs are in the first vector
    int sol_vec_id = 1;
    // recover all first derivatives of all solution components
    //utr_recover_derivatives(field_id, sol_vec_id, nr_patches, pdv_patches);
    utr_recover_derivatives_small(field_id, sol_vec_id, nr_patches, pdv_patches);

/*kbw
{
  nreq, ider, nr_deriv;
  nreq = apr_get_nreq(field_id);
  nr_deriv = 3*nreq;
  nno=0;
  while((nno=mmr_get_next_node_all(mesh_id, nno))!=0){
    double f,f_x,f_y,f_z,xcoor[3],eaux;
    printf("\nReal derivatives for node %d :\n",nno);
    mmr_node_coor(mesh_id, nno, xcoor);
    pdr_exact_sol(iaux,xcoor[0],xcoor[1],xcoor[2],daux,
		    &f,&f_x,&f_y,&f_z,&eaux);
    printf("%8.4lf%8.4lf%8.4lf", f_x, f_y, f_z);
    printf("\nRecovered derivatives for node %d :\n",nno);
    for(ider=0;ider<nr_deriv;ider++)
      {
	printf("%8.4lf", pdv_patches[nno].deriv[ider]);
      }
    printf("\n");
  }
}
/*kew*/

  }

  // enforce default adaptation strategy (utr_adapt):
  // eps is specified as the first of ADAPT_TOLERANCE_REF_UNREF
  // ratio is specified as the second of ADAPT_TOLERANCE_REF_UNREF
  // (there is one more hook - if eps is less than 0.1 it is considered
  // as global parameter and elements with error greater then eps are broken)
  utr_adapt(Problem_id, Work_dir, Interactive_input, Interactive_output);

  // option to consider: utr_adapt accepting the following parameters:
  //char name[300]; pdr_problem_name(Problem_id, name);
  //i=2; mesh_id=pdr_ctrl_i_params(Problem_id,i);
  //i=3; field_id=pdr_ctrl_i_params(Problem_id,i);
  //i=1; type=pdr_adapt_i_params(Problem_id,i);
  //i=5; eps=pdr_adapt_d_params(Problem_id,i);
  //i=6; ratio=pdr_adapt_d_params(Problem_id,i);
  //i=7; iprint=pdr_adapt_i_params(Problem_id,i);


  if(type==PDC_ADAPT_ZZ){

    // free the space - only nno_old patches; no patches for new nodes !!! */
    for(ino=1; ino<=nr_patches; ino++){ 
      if(pdv_patches[ino].deriv!=NULL) free(pdv_patches[ino].deriv);
    }
    free(pdv_patches);

  }

  return(1);

}


/*---------------------------------------------------------
pdr_err_indi - to return error indicator for an element
----------------------------------------------------------*/
double pdr_err_indi (	/* returns error indicator for an element */
        int Problem_id,	/* in: data structure to be used  */
	int Mode,	/* in: mode of operation */
	int El		/* in: element number */
)
{

  if(Mode==PDC_ADAPT_EXACT){

    return pdr_err_indi_exact(Problem_id, El);

  } else if(Mode==PDC_ADAPT_ZZ){

    return pdr_err_indi_ZZ(Problem_id, El);

  } else {

    printf("Use PDC_ADAPT_EXACT = 1 (exact) or 2 (ZZ) \n");
    printf("in problem_conv_diff.dat as the first of ADAPTATION_PARAMETERS)\n");

  }

  return(0.0);
}


/*---------------------------------------------------------
pdr_err_indi_exact - to return error indicator for an element,
	based on the knowledge of the exact solution 
----------------------------------------------------------*/
double pdr_err_indi_exact (	/* returns error indicator for an element */
        int Problem_id,	/* in: data structure to be used  */
	int El		/* in: element number */
)
{

  double l2_err = 0.0, h1_err = 0.0; /* error norms */
  int el_type;		/* element type */
  int nreq;		/* number of equations */
  int num_shap;         /* number of element shape functions */
  int ndofs;            /* local dimension of the problem */
  double determ;        /* determinant of jacobi matrix */
  double vol;           /* volume for integration rule */
  double xcoor[3];      /* global coord of gauss point */
  double u_val[PDC_CONV_DIFF_NREQ]; /* computed solution */
  double u_x[PDC_CONV_DIFF_NREQ];   /* gradient of computed solution */
  double u_y[PDC_CONV_DIFF_NREQ];   /* gradient of computed solution */
  double u_z[PDC_CONV_DIFF_NREQ];   /* gradient of computed solution */
  double base_phi[APC_MAXELVD];    /* basis functions */
  double base_dphix[APC_MAXELVD];  /* x-derivatives of basis function */
  double base_dphiy[APC_MAXELVD];  /* y-derivatives of basis function */
  double base_dphiz[APC_MAXELVD];  /* y-derivatives of basis function */
  int el_nodes[MMC_MAXELVNO+1];        /* list of nodes of El */
  double node_coor[3*MMC_MAXELVNO];  /* coord of nodes of El */
  double dofs_loc[APC_MAXELSD]; /* element solution dofs */


/* for scalar test cases - exact solution and its derivatives */
  double f, f_x, f_z, f_y;
/* auxiliary variables */
  int field_id, mesh_id;
  int i, ki, iaux, mat_num, sol_vec_id;
  double daux, eaux,volume, average, time;

  /* quadrature rules */
/*   int pdeg;		/\* degree of polynomial *\/ */
/*   static int pdeg_old=-1; /\* indicator for recomputing quadrature data *\/ */
/*   static int base_q;		/\* type of basis functions for quadrilaterals *\/ */
/*   static int ngauss;            /\* number of gauss points *\/ */
/*   static double xg[3000];   	 /\* coordinates of gauss points in 3D *\/ */
/*   static double wg[1000];       /\* gauss weights *\/ */

/* #pragma omp threadprivate (pdeg_old) */
/* #pragma omp threadprivate (base_q) */
/* #pragma omp threadprivate (ngauss) */
/* #pragma omp threadprivate (xg) */
/* #pragma omp threadprivate (wg) */

  // to make OpenMP working
  int pdeg;		/* degree of polynomial */
  int pdeg_old=-1; /* indicator for recomputing quadrature data */
  int base_q;		/* type of basis functions for quadrilaterals */
  int ngauss;            /* number of gauss points */
  double xg[3000];   	 /* coordinates of gauss points in 3D */
  double wg[1000];       /* gauss weights */

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* get formulation parameters */
  //i=1; name=pdr_ctrl_i_params(Problem_id,i); - changed for string
  char name[300];
  pdr_problem_name(Problem_id, name);
  i=2; mesh_id=pdr_ctrl_i_params(Problem_id,i);
  i=3; field_id=pdr_ctrl_i_params(Problem_id,i);
  i=5; nreq=pdr_ctrl_i_params(Problem_id,i);
  base_q=apr_get_base_type(field_id, El);

  time=pdr_time_d_params(Problem_id, 4);

/* only if exact solution is known */
  //if(name==1||name==101||name==105||name==201||name==205){
  if(strcmp(name, "DIFF_IN_BOX") == 0 || strcmp(name, "DIFF_IN_CUBE") == 0){

/* find element type */
    el_type = mmr_el_type(mesh_id,El);
    mat_num = mmr_el_groupID(mesh_id,El);

/* find degree of polynomial and number of element scalar dofs */
    apr_get_el_pdeg(field_id, El, &pdeg);
    num_shap = apr_get_el_pdeg_numshap(field_id,El,&pdeg);
    ndofs = nreq*num_shap;

/* get the coordinates of element nodes in the right order */
    mmr_el_node_coor(mesh_id,El,el_nodes,node_coor);

    /* get the most recent solution degrees of freedom */
    sol_vec_id = 1;
    apr_get_el_dofs(field_id,El,sol_vec_id,dofs_loc);

/* prepare data for gaussian integration */
    if(pdeg!=pdeg_old){
      base_q=apr_get_base_type(field_id, El);
      apr_set_quadr_3D(base_q, &pdeg, &ngauss, xg, wg);
      pdeg_old = pdeg;
    }

/*kbw
      if(El==269){
printf("In err_indi for element %d, type %d\n",El,el_type);
printf("pdeg %d, ngauss %d\n",pdeg,ngauss);
printf("nreq %d, ndof %d, local_dim %d\n",nreq,num_shap,ndofs);
printf("%d nodes with coordinates:\n",el_nodes[0]);
for(i=0;i<el_nodes[0];i++){
  printf("node %d: x - %f, y - %f, z - %f\n",
el_nodes[i+1],node_coor[3*i],node_coor[3*i+1],node_coor[3*i+2]);
}
printf("solution dofs:\n");
for(i=0;i<ndofs;i++){
  printf("%20.12lf",dofs_loc[i]);
}
printf("\n");
getchar();
      }
/*kew*/

    volume=0.0; average=0.0; double err_elem_L2 = 0.0; double err_elem_H1 = 0.0;
    for (ki=0;ki<ngauss;ki++) {
      
      /* at the gauss point, compute basis functions, determinant etc*/
      iaux = 2; /* calculations with jacobian but not on the boundary */
      determ = apr_elem_calc_3D(iaux, nreq, &pdeg, base_q,
				&xg[3*ki],node_coor,dofs_loc,
				base_phi,base_dphix,base_dphiy,base_dphiz,
				xcoor,u_val,u_x,u_y,u_z,NULL);
      
      vol = determ * wg[ki];
      
      pdr_exact_sol(mat_num,xcoor[0],xcoor[1],xcoor[2],time,
		    &f,&f_x,&f_y,&f_z,&eaux);
      

/*kbw
      if(El==269){
printf("at gauss point %d, local coor %lf, %lf, %lf\n", 
ki,xg[3*ki],xg[3*ki+1],xg[3*ki+2]);
printf("global coor %lf %lf %lf\n",xcoor[0],xcoor[1],xcoor[2]);
printf("weight %lf, determ %lf, coeff %lf\n",
	wg[ki],determ,vol);
printf("%d shape functions and derivatives: \n", num_shap);
for(i=0;i<num_shap;i++){
  printf("fun - %lf, der: x - %lf, y - %lf, z - %lf\n",
  base_phi[i],base_dphix[i],base_dphiy[i],base_dphiz[i]);
}
printf("u_val = %f, exactv = %f\n", u_val[0], f); 
printf("u_x = %f, exactv_x = %f\n", u_x[0], f_x); 
printf("u_y = %f, exactv_y = %f\n", u_y[0], f_y); 
printf("u_z = %f, exactv_z = %f\n", u_z[0], f_z); 
getchar();
      }
/*kew*/

      volume += vol;
      average+= u_val[0]*vol;
      err_elem_L2 +=  (f-u_val[0])*(f-u_val[0])*vol;
      
      err_elem_H1 += ((f_x-u_x[0])*(f_x-u_x[0])+(f_y-u_y[0])*(f_y-u_y[0])
		      +(f_z-u_z[0])*(f_z-u_z[0]))*vol;
      
    } /* ki */
    
    l2_err += err_elem_L2;
    h1_err += err_elem_H1;
/*kbw
      if(El==269){
printf("Average in element %d = %lf (at center = %lf)\n",
El,average/volume,dofs_loc[0]+dofs_loc[2]/3.0+dofs_loc[3]/3.0);
      }
/*kew*/



  } /* if exact solution known */
  else{

    printf("Cannot compute error for unknown exact solution...\n");
    return(-1.0);

  }

  return((l2_err + h1_err));

}

/*---------------------------------------------------------
pdr_zzhu_error - to compute estimated norm of error based on recovered 
		 first derivatives - the notorious ZZ error estimate
---------------------------------------------------------*/
double pdr_zzhu_error( /* returns  - Zienkiewicz-Zhu error for the whole mesh */
  int Field_id,    /* in: approximation field ID  */
  FILE *interactive_output
		)
{

/* local variables */
  double err_ZZ = 0.0; /* error norms */
  int i,j,k, mesh_id, nreq, nr_deriv, ino, nno, nr_patches, nel;
  double clock_time;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the corresponding mesh */
  mesh_id = apr_get_mesh_id(Field_id);

  nreq = apr_get_nreq(Field_id);
  nr_deriv = 3*nreq;

  // create patches for nno_old nodes
  clock_time = time_clock();
  nr_patches = utr_create_patches(Field_id, &pdv_patches);
  //nr_patches = utr_create_patches_small(Field_id, &pdv_patches);
  clock_time = time_clock()-clock_time;
  fprintf(interactive_output,"Patch creation - time %lf\n", clock_time);


  // recover all first derivatives of all solution components
  clock_time = time_clock();
  
  // !!! it is assumed that the most recent dofs are in the last vector
  //int sol_vec_id = apr_get_nr_sol(Field_id);
  // !!! it is assumed that the most recent dofs are in the first vector
  int sol_vec_id = 1;
  // recover all first derivatives of all solution components
  utr_recover_derivatives(Field_id, sol_vec_id, nr_patches, pdv_patches);
  //utr_recover_derivatives_small(Field_id, sol_vec_id, nr_patches, pdv_patches);
  
  clock_time = time_clock()-clock_time;
  fprintf(interactive_output,"Derivative recovery - time %lf\n", clock_time);

/*kbw
  nno=0;
  while((nno=mmr_get_next_node_all(mesh_id, nno))!=0){
    double f,f_x,f_y,f_z,xcoor[3],eaux;
    printf("\nReal derivatives for node %d :\n",nno);
    mmr_node_coor(mesh_id, nno, xcoor);
    int iaux = 0; double daux = 0.0;
    pdr_exact_sol(iaux,xcoor[0],xcoor[1],xcoor[2],daux,
		    &f,&f_x,&f_y,&f_z,&eaux);
    printf("%8.4lf%8.4lf%8.4lf", f_x, f_y, f_z);
    printf("\nRecovered derivatives for node %d :\n",nno);
    int ider;
    for(ider=0;ider<nr_deriv;ider++)
      {
	printf("%8.4lf", pdv_patches[nno].deriv[ider]);
      }
    printf("\n");
  }
/*kew*/

  nel=0;
  err_ZZ=0.0;
  while((nel=mmr_get_next_act_elem(mesh_id, nel))!=0)
    {
      err_ZZ += pdr_err_indi_ZZ(Field_id,nel);
/*kbw
printf("in active element %d - error_indi %20.15lf \n",
nel,pdr_err_indi_ZZ(Field_id,nel));
/*kew*/
    }
  err_ZZ = sqrt(err_ZZ);
  fprintf(interactive_output,"\nZienkiewicz-Zhu error estimator =  %lf\n",err_ZZ);
  
  // free the space - only nno_old patches; no patches for new nodes !!! */
  for(ino=1; ino<=nr_patches; ino++){ 
    if(pdv_patches[ino].deriv!=NULL) free(pdv_patches[ino].deriv);
  }
  free(pdv_patches);

  return(err_ZZ);

}

/*---------------------------------------------------------
pdr_zzhu_error_small - to compute estimated norm of error based on recovered 
		 first derivatives - the notorious ZZ error estimate
(version with small patches - for hanging nodes as well)
---------------------------------------------------------*/
double pdr_zzhu_error_small( /* returns  - Zienkiewicz-Zhu error for the whole mesh */
  int Field_id,    /* in: approximation field ID  */
  FILE *interactive_output
		)
{

/* local variables */
  double err_ZZ = 0.0; /* error norms */
  int i,j,k, mesh_id, nreq, nr_deriv, ino, nno, nr_patches, nel;
  double clock_time;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the corresponding mesh */
  mesh_id = apr_get_mesh_id(Field_id);

  nreq = apr_get_nreq(Field_id);
  nr_deriv = 3*nreq;

  // create patches for nno_old nodes
  clock_time = time_clock();
  //nr_patches = utr_create_patches(Field_id, &pdv_patches);
  nr_patches = utr_create_patches_small(Field_id, &pdv_patches);
  clock_time = time_clock()-clock_time;
  fprintf(interactive_output,"Patch creation - time %lf\n", clock_time);


  // recover all first derivatives of all solution components
  clock_time = time_clock();
  
  // !!! it is assumed that the most recent dofs are in the last vector
  //int sol_vec_id = apr_get_nr_sol(Field_id);
  // !!! it is assumed that the most recent dofs are in the first vector
  int sol_vec_id = 1;
  // recover all first derivatives of all solution components
  //utr_recover_derivatives(Field_id, sol_vec_id, nr_patches, pdv_patches);
  utr_recover_derivatives_small(Field_id, sol_vec_id, nr_patches, pdv_patches);
  
  clock_time = time_clock()-clock_time;
  fprintf(interactive_output,"Derivative recovery - time %lf\n", clock_time);

/*kbw
  nno=0;
  while((nno=mmr_get_next_node_all(mesh_id, nno))!=0){
    double f,f_x,f_y,f_z,xcoor[3],eaux;
    printf("\nReal derivatives for node %d :\n",nno);
    mmr_node_coor(mesh_id, nno, xcoor);
    int iaux = 0; double daux = 0.0;
    pdr_exact_sol(iaux,xcoor[0],xcoor[1],xcoor[2],daux,
		    &f,&f_x,&f_y,&f_z,&eaux);
    printf("%8.4lf%8.4lf%8.4lf", f_x, f_y, f_z);
    printf("\nRecovered derivatives for node %d :\n",nno);
    int ider;
    for(ider=0;ider<nr_deriv;ider++)
      {
	printf("%8.4lf", pdv_patches[nno].deriv[ider]);
      }
    printf("\n");
  }
/*kew*/

  nel=0;
  err_ZZ=0.0;
  while((nel=mmr_get_next_act_elem(mesh_id, nel))!=0)
    {
      err_ZZ += pdr_err_indi_ZZ(Field_id,nel);
/*kbw
printf("in active element %d - error_indi %20.15lf \n",
nel,pdr_err_indi_ZZ(Field_id,nel));
/*kew*/
    }
  err_ZZ = sqrt(err_ZZ);
  fprintf(interactive_output,"\nZienkiewicz-Zhu error estimator =  %lf\n",err_ZZ);
  
  // free the space - only nno_old patches; no patches for new nodes !!! */
  for(ino=1; ino<=nr_patches; ino++){ 
    if(pdv_patches[ino].deriv!=NULL) free(pdv_patches[ino].deriv);
  }
  free(pdv_patches);

  return(err_ZZ);

}


/*---------------------------------------------------------
pdr_err_indi_ZZ - to return error indicator for an element
----------------------------------------------------------*/
double pdr_err_indi_ZZ (	/* returns error indicator for an element */
        int Problem_id,	/* in: data structure to be used  */
	int El		/* in: element number */
)
{
  double err=0;
  int i,ki,j,k;
  int field_id,mesh_id,pdeg,num_shap,sol_vec_id,iaux,nreq;
  int base;		/* type of basis functions  */
  int el_nodes[MMC_MAXELVNO+1];        /* list of nodes of El */
  double node_coor[3*MMC_MAXELVNO];  /* coord of nodes of El */
  double dofs_loc[APC_MAXELSD]; /* element solution dofs */
  double xcoor[3];      /* global coord of gauss point */
  double u_val[PDC_CONV_DIFF_NREQ]; /* computed solution */
  double u_x[PDC_CONV_DIFF_NREQ];   /* gradient of computed solution */
  double u_y[PDC_CONV_DIFF_NREQ];   /* gradient of computed solution */
  double u_z[PDC_CONV_DIFF_NREQ];   /* gradient of computed solution */
  double base_phi[APC_MAXELVD];    /* basis functions */
  double base_dphix[APC_MAXELVD];  /* x-derivatives of basis function */
  double base_dphiy[APC_MAXELVD];  /* y-derivatives of basis function */
  double base_dphiz[APC_MAXELVD];  /* y-derivatives of basis function */
  double determ;
  int ngauss;           /* number of gauss points */
  double xg[3000];   	/* coordinates of gauss points in 3D */
  double wg[1000];      /* gauss weights */
  double vol;           /* volume for integration rule */
  double deriv_el_nodes[3][9]; /* derivatives at element nodes */
  double deriv_gauss_points[3];
  
/*++++++++++++++++ executable statements ++++++++++++++++*/

  i=3; field_id=pdr_ctrl_i_params(Problem_id,i);
  mesh_id = apr_get_mesh_id(field_id);
  base=apr_get_base_type(field_id, El);
  nreq = apr_get_nreq(field_id);

  //printf("\nErr indi Element %d:\n",El);

  /* find degree of polynomial and number of element scalar dofs */
  apr_get_el_pdeg(field_id, El, &pdeg);

  /* get the coordinates of the nodes of El in the right order */
  mmr_el_node_coor(mesh_id,El,el_nodes,node_coor);
  
  /* prepare data for gaussian integration */
  apr_set_quadr_3D(base, &pdeg, &ngauss, xg, wg);
  
  /* get the most recent solution degrees of freedom */
  sol_vec_id = 1;
  apr_get_el_dofs(field_id, El, sol_vec_id, dofs_loc);

  // get 3 derivatives for each solution component
  for(j=0;j<3;j++)
    {
      k=0;
      for(i=1;i<=el_nodes[0];i++)
	{
	  deriv_el_nodes[j][k]=pdv_patches[el_nodes[i]].deriv[j];
	  k++;
	}
    }
  for(j=0;j<3;j++) deriv_gauss_points[j]=0.0;

  
  for (ki=0;ki<ngauss;ki++)
    {
      /* at the gauss point, compute basis functions, determinant etc*/
      iaux = 2; /* calculations with jacobian but not on the boundary */
      determ = apr_elem_calc_3D(iaux, nreq, &pdeg, base, 
				&xg[3*ki],node_coor,dofs_loc, 
				base_phi,base_dphix,base_dphiy,base_dphiz, 
				xcoor,u_val,u_x,u_y,u_z,NULL);
      vol = determ * wg[ki];
      
      // get derivatives at gauss points
      for(i=0;i<3;i++){
	iaux = 1;
	apr_elem_calc_3D(iaux, 1, &pdeg, base, &xg[3*ki], 
			 node_coor,&deriv_el_nodes[i][0],
			 base_phi,base_dphix,base_dphiy,base_dphiz, 
			 xcoor,&deriv_gauss_points[i],NULL,NULL,NULL,NULL);
      }

      err += vol * (
		   (deriv_gauss_points[0]-u_x[0])*(deriv_gauss_points[0]-u_x[0])+
		   (deriv_gauss_points[1]-u_y[0])*(deriv_gauss_points[1]-u_y[0])+
		   (deriv_gauss_points[2]-u_z[0])*(deriv_gauss_points[2]-u_z[0])
		   );
            
    }
    
  return(err);

}

