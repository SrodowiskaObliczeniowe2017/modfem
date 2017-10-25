/************************************************************************
File aps_dg_util.c - utilities functions necessary for discontinuous
                         approximation on prismatic elements

Contains implementation of routines:   
  apr_select_field - to select the proper field   
  apr_get_base_type - to return the type of basis functions
  apr_elem_calc_3D - to perform element calculations (to provide data on
	coordinates, solution, shape functions, etc. for a given point
	inside element (given local coordinates Eta[i]);
	for geometrically multi-linear or linear 3D elements
  apr_set_quadr_3D - to prepare quadrature data for a given element
  apr_set_quadr_2D - to prepare quadrature data for a given face
  apr_L2_proj - to L2 project a function onto an element
  apr_sol_xglob - to return the solution at a point with global
	coordinates specified
  apr_spec_ini_con - to specify initial condition

------------------------------  			
History:     
	02.2002 - Krzysztof Banas, initial version		
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

/* interface for all approximation modules */
#include "aph_intf.h"	

/* internal header file for the dg approximation module */
#include "./aph_dg_prism.h"	

/* interface of the mesh manipulation module */
#include "mmh_intf.h"	

/* LAPACK procedures */
#include "lin_alg_intf.h"

#include "uth_log.h"

#define SMALL 1e-10

/* internal utility procedures */
double ut_mat3_inv(
	/* returns: determinant of matrix to invert */
	double *mat,	/* matrix to invert */
	double *mat_inv	/* inverted matrix */
	);
void ut_vec3_prod(
	double* vec_a, 	/* in: vector a */
	double* vec_b, 	/* in: vector b */
	double* vec_c	/* out: vector product axb */
	);
double ut_vec3_mxpr( /* returns: mixed product [a,b,c] */
	double* vec_a, 	/* in: vector a */
	double* vec_b, 	/* in: vector b */
	double* vec_c	/* in: vector c */
	);
double ut_vec3_length(	/* returns: vector length */
	double* vec	/* in: vector */
	);
void ut_mat3vec(
	double* m1, 	/* in: matrix (stored by rows as a vector!) */
	double* v1, 	/* in: vector */
	double* v2	/* out: resulting vector */
	);
void ut_mat3mat(
	double* m1,	/* in: matrix */
	double* m2,	/* in: matrix */ 
	double* m3	/* out: matrix m1*m2 */
	);
int ut_chk_list(	/* returns: */
			/* >0 - position on the list */
            		/* 0 - not found on the list */
	int Num, 	/* number to be checked */
	int* List, 	/* list of numbers */
	int Ll		/* length of the list */
	);

/* standard macro for max and min and abs */
#define ut_max(x,y) ((x)>(y)?(x):(y))
#define ut_min(x,y) ((x)<(y)?(x):(y))
#define ut_abs(x)   ((x)<0?-(x):(x))

/*---------------------------------------------------------
  apr_select_field - to select the proper field   
---------------------------------------------------------*/
apt_field* apr_select_field( /* returns: pointer to the selected field */
  int Field_id  /* in: field ID */
  )
{
  if( Field_id == APC_CUR_FIELD_ID ) {
    return(&apv_fields[apv_cur_field_id-1]);
  }
  else if( Field_id>0 && Field_id<=apv_nr_fields ) {
    return(&apv_fields[Field_id-1]);
  }
  else {
    return(&apv_fields[apv_cur_field_id-1]);
    /* alternative:  return(NULL);   */
  }
}


/*------------------------------------------------------------------
apr_elem_calc_3D - to perform element calculations (to provide data on
	coordinates, solution, shape functions, etc. for a given point
	inside element (given local coordinates Eta[i]);
	for geometrically multi-linear or linear 3D elements
-------------------------------------------------------------------*/
double apr_elem_calc_3D( 
	/* returns: Jacobian determinant at a point, either for */
	/* 	volume integration if Vec_norm==NULL,  */
	/* 	or for surface integration otherwise */
	int Control,	    /* in: control parameter (what to compute): */
			    /*	1  - shape functions and values */
			    /*	2  - derivatives and jacobian */
			    /* 	>2 - computations on the (Control-2)-th */
			    /*	     element's face */
	int Nreq,	    /* in: number of equations */
	int *Pdeg_vec,	    /* in: element degree of polynomial */
	int Base_type,	    /* in: type of basis functions: */
			    /* 	(APC_BASE_TENSOR_DG) - tensor product */
			    /* 	(APC_BASE_COMPLETE_DG) - complete polynomials */
        double *Eta,	    /* in: local coordinates of the input point */
	double *Node_coor,  /* in: array of coordinates of vertices of element */
	double *Sol_dofs,   /* in: array of element' dofs */
	double *Base_phi,   /* out: basis functions */
	double *Base_dphix, /* out: x-derivatives of basis functions */
	double *Base_dphiy, /* out: y-derivatives of basis functions */
	double *Base_dphiz, /* out: z-derivatives of basis functions */
	double *Xcoor,	    /* out: global coordinates of the point*/
	double *Sol,        /* out: solution at the point */
	double *Dsolx,      /* out: derivatives of solution at the point */
	double *Dsoly,      /* out: derivatives of solution at the point */
	double *Dsolz,      /* out: derivatives of solution at the point */
	double *Vec_nor     /* out: outward unit vector normal to the face */
	)
{
/* local variables */
  int nrgeo;		/* number of element geometry dofs */
  int nrdofs;		/* number of element solution dofs */
  int pdeg;
  double det,dxdeta[9],detadx[9]; /* jacobian, jacobian matrix 
	and its inverse; for geometrical transformation */

  double ds[3]; /* coordinates of vector normal to a face */

  double geo_phi[8];  	/* geometry shape functions */
  double geo_dphix[8]; 	/* derivatives of geometry shape functions */
  double geo_dphiy[8];  /* derivatives of geometry shape functions */
  double geo_dphiz[8];  /* derivatives of geometry shape functions */

/* auxiliary variables */
  int i, ieq;
  double daux, faux;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* check parameter Nreq */
#ifdef DEBUG_APM
  if(Nreq>APC_MAXEQ){
    printf("Number of equations %d greater than the limit %d.\n",Nreq,APC_MAXEQ);
    printf("Change APC_MAXEQ in aph_dg_prism.h and recompile the code.\n");
    exit(-1);
  }
#endif

  pdeg = Pdeg_vec[0];

/* set the number of geometry dofs */
  nrgeo = 6;

/* get the values of shape functions and their LOCAL derivatives */
  if(Base_phi!=NULL){
    if(Control==1)
      nrdofs=apr_shape_fun_3D(Base_type, pdeg, Eta,
	Base_phi, NULL, NULL, NULL);
    else
      nrdofs=apr_shape_fun_3D(Base_type, pdeg, Eta,
	Base_phi, Base_dphix, Base_dphiy, Base_dphiz);
  }

  geo_phi[0]=(1.0-Eta[0]-Eta[1])*(1.0-Eta[2])/2.0;
  geo_phi[1]=Eta[0]*(1.0-Eta[2])/2.0;
  geo_phi[2]=Eta[1]*(1.0-Eta[2])/2.0;
  geo_phi[3]=(1.0-Eta[0]-Eta[1])*(1.0+Eta[2])/2.0;
  geo_phi[4]=Eta[0]*(1.0+Eta[2])/2.0;
  geo_phi[5]=Eta[1]*(1.0+Eta[2])/2.0;

/*physical coordinates*/
  if(Xcoor!=NULL){

#ifdef DEBUG_APM
      if(Node_coor==NULL){
	printf("Error 23769483 in shape_fun. Exiting!\n");
	exit(-1);
      }
#endif

    Xcoor[0] = 0;
    Xcoor[1] = 0;
    Xcoor[2] = 0;
    for(i=0;i<nrgeo;i++){
      Xcoor[0] += Node_coor[3*i]*geo_phi[i];
      Xcoor[1] += Node_coor[3*i+1]*geo_phi[i];
      Xcoor[2] += Node_coor[3*i+2]*geo_phi[i];
    }
  }

/* function value */
  if(Sol!=NULL){

#ifdef DEBUG_APM
      if(Sol_dofs==NULL){
	printf("Error 23945483 in shape_fun. Exiting!\n");
	exit(-1);
      }
#endif

    for(ieq=0;ieq<Nreq;ieq++){
      Sol[ieq]=0.0;
      for(i=0;i<nrdofs;i++){
        Sol[ieq] += Sol_dofs[nrdofs*ieq+i]*Base_phi[i];
      }
    }
  }
 
/*if no computations involving derivatives of shape functions*/
  if(Control==1) return(0.0);

/* local derivatives of geometrical shape functions*/
  geo_dphix[0] = -(1.0-Eta[2])/2.0;  
  geo_dphix[1] =  (1.0-Eta[2])/2.0;  
  geo_dphix[2] =  0.0; 
  geo_dphix[3] = -(1.0+Eta[2])/2.0;  
  geo_dphix[4] =  (1.0+Eta[2])/2.0;  
  geo_dphix[5] =  0.0; 
  geo_dphiy[0] = -(1.0-Eta[2])/2.0;
  geo_dphiy[1] =  0.0;
  geo_dphiy[2] =  (1.0-Eta[2])/2.0;
  geo_dphiy[3] = -(1.0+Eta[2])/2.0;
  geo_dphiy[4] =  0.0;
  geo_dphiy[5] =  (1.0+Eta[2])/2.0;
  geo_dphiz[0] = -(1.0-Eta[0]-Eta[1])/2.0;
  geo_dphiz[1] = -Eta[0]/2.0;
  geo_dphiz[2] = -Eta[1]/2.0;
  geo_dphiz[3] =  (1.0-Eta[0]-Eta[1])/2.0;
  geo_dphiz[4] =  Eta[0]/2.0;
  geo_dphiz[5] =  Eta[1]/2.0;

/* Jacobian matrix J */
  dxdeta[0] = 0.0; dxdeta[1] = 0.0; dxdeta[2] = 0.0;
  dxdeta[3] = 0.0; dxdeta[4] = 0.0; dxdeta[5] = 0.0;
  dxdeta[6] = 0.0; dxdeta[7] = 0.0; dxdeta[8] = 0.0;
  for(i=0;i<nrgeo;i++){
    dxdeta[0] += Node_coor[3*i]  *geo_dphix[i];
    dxdeta[1] += Node_coor[3*i]  *geo_dphiy[i];
    dxdeta[2] += Node_coor[3*i]  *geo_dphiz[i];
    dxdeta[3] += Node_coor[3*i+1]*geo_dphix[i];
    dxdeta[4] += Node_coor[3*i+1]*geo_dphiy[i];
    dxdeta[5] += Node_coor[3*i+1]*geo_dphiz[i];
    dxdeta[6] += Node_coor[3*i+2]*geo_dphix[i];
    dxdeta[7] += Node_coor[3*i+2]*geo_dphiy[i];
    dxdeta[8] += Node_coor[3*i+2]*geo_dphiz[i];
  }

/* Jacobian |J| and inverse of the Jacobian matrix*/
  det = ut_mat3_inv(dxdeta,detadx);

/* global derivatives of geometrical shape functions */
  for(i=0;i<nrgeo;i++){
    daux = geo_dphix[i]*detadx[0] 
         + geo_dphiy[i]*detadx[3]
         + geo_dphiz[i]*detadx[6];
    faux = geo_dphix[i]*detadx[1]
         + geo_dphiy[i]*detadx[4]
         + geo_dphiz[i]*detadx[7];
    geo_dphiz[i] = geo_dphix[i]*detadx[2]
         + geo_dphiy[i]*detadx[5]
         + geo_dphiz[i]*detadx[8];
    geo_dphix[i] = daux;
    geo_dphiy[i] = faux;
  }

/* global derivatives of solution shape functions */
  if(Base_dphix!=NULL){

#ifdef DEBUG_APM
      if(Base_dphiy==NULL&&Base_dphiz==NULL){
	printf("Error 23983483 in shape_fun. Exiting!\n");
	exit(-1);
      }
#endif

    for(i=0;i<nrdofs;i++){
      daux = Base_dphix[i]*detadx[0] 
         +   Base_dphiy[i]*detadx[3]
         +   Base_dphiz[i]*detadx[6];
      faux = Base_dphix[i]*detadx[1]
         +   Base_dphiy[i]*detadx[4]
         +   Base_dphiz[i]*detadx[7];
      Base_dphiz[i] =   Base_dphix[i]*detadx[2]
         +   Base_dphiy[i]*detadx[5]
         +   Base_dphiz[i]*detadx[8];
      Base_dphix[i] = daux;
      Base_dphiy[i] = faux;
    }
  }

/* global derivatives of solution */
  if(Dsolx!=NULL){

#ifdef DEBUG_APM
      if(Dsoly==NULL&&Dsolz==NULL){
	printf("Error 239828483 in shape_fun. Exiting!\n");
	exit(-1);
      }
#endif

    for(ieq=0;ieq<Nreq;ieq++){
      Dsolx[ieq]=0.0;
      Dsoly[ieq]=0.0;
      Dsolz[ieq]=0.0;
      for(i=0;i<nrdofs;i++){
        Dsolx[ieq] += Sol_dofs[nrdofs*ieq+i]*Base_dphix[i];
        Dsoly[ieq] += Sol_dofs[nrdofs*ieq+i]*Base_dphiy[i];
        Dsolz[ieq] += Sol_dofs[nrdofs*ieq+i]*Base_dphiz[i];
      }
    }
  }

  if(Control==2) return(fabs(det));


/* area element dS = vector normal */
  if(Control==3){
    ds[0] = - dxdeta[3]*dxdeta[7] + dxdeta[6]*dxdeta[4];
    ds[1] = - dxdeta[6]*dxdeta[1] + dxdeta[0]*dxdeta[7];
    ds[2] = - dxdeta[0]*dxdeta[4] + dxdeta[3]*dxdeta[1];
  }
  if(Control==4){
    ds[0] = dxdeta[3]*dxdeta[7] - dxdeta[6]*dxdeta[4];
    ds[1] = dxdeta[6]*dxdeta[1] - dxdeta[0]*dxdeta[7];
    ds[2] = dxdeta[0]*dxdeta[4] - dxdeta[3]*dxdeta[1];
  }
  else if(Control==5){
    ds[0] = dxdeta[3]*dxdeta[8] - dxdeta[6]*dxdeta[5];
    ds[1] = dxdeta[6]*dxdeta[2] - dxdeta[0]*dxdeta[8];
    ds[2] = dxdeta[0]*dxdeta[5] - dxdeta[3]*dxdeta[2];
  }
  else if(Control==6){
    ds[0] = (dxdeta[4]-dxdeta[3])*dxdeta[8] 
      - (dxdeta[7]-dxdeta[6])*dxdeta[5];
    ds[1] = (dxdeta[7]-dxdeta[6])*dxdeta[2] 
      - (dxdeta[1]-dxdeta[0])*dxdeta[8];
    ds[2] = (dxdeta[1]-dxdeta[0])*dxdeta[5] 
      - (dxdeta[4]-dxdeta[3])*dxdeta[2];
  }
  else if(Control==7){
    ds[0] = dxdeta[5]*dxdeta[7] - dxdeta[8]*dxdeta[4];
    ds[1] = dxdeta[8]*dxdeta[1] - dxdeta[2]*dxdeta[7];
    ds[2] = dxdeta[2]*dxdeta[4] - dxdeta[5]*dxdeta[1];
  }
  

  det = ut_vec3_length(ds);

  if (Vec_nor != NULL) {

/* normalize vector normal */
    Vec_nor[0] = ds[0]/det;
    Vec_nor[1] = ds[1]/det;
    Vec_nor[2] = ds[2]/det;

  }

  return(det);

}

/* local functions for apr_shape_fun_3D and apr_shape_fun_2D */
int apr_shape_fun_2D( /* returns: the number of shape functions */
	int Base_type,	   /* in: type of basis functions: */
			   /*   (APC_BASE_TENSOR_DG) - tensor product */
			   /* 	(APC_BASE_COMPLETE_DG) - complete polynomials */
	int Pdeg, 	   /* in: degree of polynomial - can be either */
			   /*	a single number, for isotropic p, */
			   /*	or a combination pdegy*10+pdegx */
	double *Eta,	   /* in: local coord of the considered point */
	double *Base_phi,  /* out: basis functions */
	double *Base_dphix,/* out: x derivative of basis functions */
	double *Base_dphiy /* out: y derivative of basis functions */
	);

int apr_shape_fun_1D( /* returns: the number of shape functions */
	int Pdeg, 	   /* in: degree of polynomial */
	double Eta,	   /* in: local coord of the considered point */
	double *Base_phi,  /* out: basis functions */
	double *Base_dphix/* out: x derivative of basis functions */
	);

/*---------------------------------------------------------
apr_shape_fun_3D - to compute values of shape functions and their
	local derivatives at a point within the master 3D element
----------------------------------------------------------*/
int apr_shape_fun_3D( /* returns: the number of shape functions (<=0 - failure) */
	int Base_type,	   /* in: type of basis functions: */
			   /*   (APC_BASE_TENSOR_DG) - tensor product */
			   /* 	(APC_BASE_COMPLETE_DG) - complete polynomials */
	int Pdeg, 	   /* in: degree of polynomial - can be either */
			   /*	a single number, for isotropic p, */
			   /*	or a combination pdegy*10+pdegx */
	double *Eta,	   /* in: local coord of the considered point */
	double *Base_phi,  /* out: basis functions */
	double *Base_dphix,/* out: x derivative of basis functions */
	double *Base_dphiy,/* out: y derivative of basis functions */
	double *Base_dphiz /* out: z derivative of basis functions */
	)
{

/* auxiliary variables */
  int i, j, k, ii, kk, pdegx, pdegz, num_shap, num_shap_2D, base_type_2D;
  double  ssn[APC_MAXELVD];
  double  uun[APC_MAXELVD];
  double dssn[APC_MAXELVD];
  double dttn[APC_MAXELVD];
  double duun[APC_MAXELVD];

  int num_shap_old;
  double Base_phi_old[APC_MAXELVD];
  double Base_dphix_old[APC_MAXELVD];
  double Base_dphiy_old[APC_MAXELVD];
  double Base_dphiz_old[APC_MAXELVD];

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef DEBUG_APM
  if(Base_phi==NULL) {
    printf("No space provided for shape functions in apr_shape_fun_3D!\n");
    exit(-1);
  }
#endif

  if(Pdeg==0){
    /* for piecewise constant approximation */
    
    Base_phi[0] = 1.0;
    
    if(Base_dphix!=NULL){

#ifdef DEBUG_APM
      if(Base_dphiy==NULL&&Base_dphiz==NULL){
	printf("Error 239483 in shape_fun. Exiting!\n");
	exit(-1);
      }
#endif

      Base_dphix[0] = 0.0;
      Base_dphiy[0] = 0.0;
      Base_dphiz[0] = 0.0;
    }
    
    num_shap=1;
    
    return(num_shap);
    
  }
  else if(Pdeg==100){
    /* for linear in z and piecewise constant in xy approximation */
    
    Base_phi[0] = 1.0;
    Base_phi[1] = Eta[2];
    
    if(Base_dphix!=NULL){
      
#ifdef DEBUG_APM
      if(Base_dphiy==NULL&&Base_dphiz==NULL){
	printf("Error 23759483 in shape_fun. Exiting!\n");
	exit(-1);
      }
#endif

      Base_dphix[0] = 0.0;
      Base_dphix[1] = 0.0;
      Base_dphiy[0] = 0.0;
      Base_dphiy[1] = 0.0;
      Base_dphiz[0] = 0.0;
      Base_dphiz[1] = 1.0;
      
    }
    
    num_shap=2;
    
    return(num_shap);
    
  }
  else if(Pdeg==1&&Base_type == APC_BASE_TENSOR_DG){
    /* for piecewise constant in z and linear in xy approximation */
    
    Base_phi[0] = 1.0;
    Base_phi[1] = Eta[0];
    Base_phi[2] = Eta[1];
    
    if(Base_dphix!=NULL){
      
#ifdef DEBUG_APM
      if(Base_dphiy==NULL&&Base_dphiz==NULL){
	printf("Error 237339483 in shape_fun. Exiting!\n");
	exit(-1);
      }
#endif

      Base_dphix[0] = 0.0;
      Base_dphix[1] = 1.0;
      Base_dphix[2] = 0.0;
      
      Base_dphiy[0] = 0.0;
      Base_dphiy[1] = 0.0;
      Base_dphiy[2] = 1.0;
      
      Base_dphiz[0] = 0.0;
      Base_dphiz[1] = 0.0;
      Base_dphiz[2] = 0.0;
      
    }
    
    num_shap=3;
    
    return(num_shap);
    
  }
  else if(Pdeg==2&&Base_type == APC_BASE_TENSOR_DG){
/* for piecewise constant in z and quadratic in xy approximation */
    
    Base_phi[0] = 1.0;
    Base_phi[1] = Eta[0];
    Base_phi[2] = Eta[1];
    Base_phi[3] = Eta[0]*Eta[0];
    Base_phi[4] = Eta[0]*Eta[1];
    Base_phi[5] = Eta[1]*Eta[1];
    
    if(Base_dphix!=NULL){
      
#ifdef DEBUG_APM
      if(Base_dphiy==NULL&&Base_dphiz==NULL){
	printf("Error 23986483 in shape_fun. Exiting!\n");
	exit(-1);
      }
#endif

      Base_dphix[0] = 0.0;
      Base_dphix[1] = 1.0;
      Base_dphix[2] = 0.0;
      Base_dphix[3] = 2*Eta[0];
      Base_dphix[4] = Eta[1];
      Base_dphix[5] = 0.0;
      
      Base_dphiy[0] = 0.0;
      Base_dphiy[1] = 0.0;
      Base_dphiy[2] = 1.0;
      Base_dphiy[3] = 0.0;
      Base_dphiy[4] = Eta[0];
      Base_dphiy[5] = 2*Eta[1];
      
      Base_dphiz[0] = 0.0;
      Base_dphiz[1] = 0.0;
      Base_dphiz[2] = 0.0;
      Base_dphiz[3] = 0.0;
      Base_dphiz[4] = 0.0;
      Base_dphiz[5] = 0.0;
      
    }
    
    num_shap=6;
    
    return(num_shap);
    
  }
  else{
    /* for linear and higher order approximation */
    
    Base_phi[0] = 1.0;
    Base_phi[1] = Eta[0];
    Base_phi[2] = Eta[1];
    Base_phi[3] = Eta[2];
    
    if(Base_dphix!=NULL){
      
#ifdef DEBUG_APM
      if(Base_dphiy==NULL&&Base_dphiz==NULL){
	printf("Error 239474483 in shape_fun. Exiting!\n");
	exit(-1);
      }
#endif

      Base_dphix[0] = 0.0;
      Base_dphix[1] = 1.0;
      Base_dphix[2] = 0.0;
      Base_dphix[3] = 0.0;
      
      Base_dphiy[0] = 0.0;
      Base_dphiy[1] = 0.0;
      Base_dphiy[2] = 1.0;
      Base_dphiy[3] = 0.0;
      
      Base_dphiz[0] = 0.0;
      Base_dphiz[1] = 0.0;
      Base_dphiz[2] = 0.0;
      Base_dphiz[3] = 1.0;
      
    }
    
    if(Pdeg==1){
      /* to speed up computations with linear shape functions */
      
      num_shap=4;
      return(num_shap);
      
    }
    
    if(Base_type==APC_BASE_TENSOR_DG&&(Pdeg<100||Pdeg%100==0)){
      printf("Sorry %d, at least linears required both for z and xy polynomials!\n", Pdeg);
      exit(-1);
    }
    
    
    Base_phi[4] = Eta[2]*Eta[0];
    Base_phi[5] = Eta[2]*Eta[1];
    
    if(Base_dphix!=NULL){

#ifdef DEBUG_APM
      if(Base_dphiy==NULL&&Base_dphiz==NULL){
	printf("Error 23948283 in shape_fun. Exiting!\n");
	exit(-1);
      }
#endif
      
      Base_dphix[4] = Eta[2];
      Base_dphix[5] = 0.0;
      
      Base_dphiy[4] = 0.0;
      Base_dphiy[5] = Eta[2];
      
      Base_dphiz[4] = Eta[0];
      Base_dphiz[5] = Eta[1];
      
    }
    
    if(Pdeg==101){
      /* to speed up computations with bi/linear shape functions */
      
      num_shap=6;
      return(num_shap);
      
    }
    
    Base_phi[6] = Eta[1]*Eta[0];
    if(Base_dphix!=NULL){

#ifdef DEBUG_APM
      if(Base_dphiy==NULL&&Base_dphiz==NULL){
	printf("Error 23964483 in shape_fun. Exiting!\n");
	exit(-1);
      }
#endif

      Base_dphix[6] = Eta[1];
      Base_dphiy[6] = Eta[0];
      Base_dphiz[6] = 0.0;
    }
    
    if(Pdeg==2&&Base_type==APC_BASE_COMPLETE_DG){
      /* to speed up computations with 2nd order complete polynomial shape functions*/
      
      Base_phi[7] = Eta[0]*Eta[0];
      Base_phi[8] = Eta[1]*Eta[1];
      Base_phi[9] = Eta[2]*Eta[2];
      
      if(Base_dphix!=NULL){
	
#ifdef DEBUG_APM
      if(Base_dphiy==NULL&&Base_dphiz==NULL){
	printf("Error 239483 in shape_fun. Exiting!\n");
	exit(-1);
      }
#endif
	Base_dphix[7] = 2.0*Eta[0];
	Base_dphix[8] = 0.0;
	Base_dphix[9] = 0.0;
	
	Base_dphiy[7] = 0.0;
	Base_dphiy[8] = 2.0*Eta[1];
	Base_dphiy[9] = 0.0;
	
	Base_dphiz[7] = 0.0;
	Base_dphiz[8] = 0.0;
	Base_dphiz[9] = 2.0*Eta[2];
	
      }
      
      num_shap=10;
      return(num_shap);
      
    }
    
/* for second (TENSOR) and higher order of approximation */
    Base_phi[7] = Eta[2]*Eta[1]*Eta[0];
    if(Base_dphix!=NULL){
#ifdef DEBUG_APM
      if(Base_dphiy==NULL&&Base_dphiz==NULL){
	printf("Error 239483 in shape_fun. Exiting!\n");
	exit(-1);
      }
#endif
      Base_dphix[7] = Eta[2]*Eta[1];
      Base_dphiy[7] = Eta[2]*Eta[0];
      Base_dphiz[7] = Eta[1]*Eta[0];
    }
    
    base_type_2D=APC_BASE_COMPLETE_DG;
    
    if(Base_type==APC_BASE_TENSOR_DG){
      
      pdegx = Pdeg%100;
      pdegz = Pdeg/100;
      
    }
    else if(Base_type==APC_BASE_COMPLETE_DG){
      
      if(Pdeg<100){
	pdegx = Pdeg;
	pdegz = Pdeg;
      }
      else{
	pdegx = Pdeg%100;
	pdegz = Pdeg/100;
      }
      
    }
    else {
      
      printf("Type of base not valid for prisms!\n");
      exit(-1);
      
    }
    
    if(Base_dphix!=NULL){
      
#ifdef DEBUG_APM
      if(Base_dphiy==NULL&&Base_dphiz==NULL){
	printf("Error 239483 in shape_fun. Exiting!\n");
	exit(-1);
      }
#endif

      /* get 2D shape functions for bases */
      num_shap_2D = apr_shape_fun_2D(base_type_2D,pdegx,Eta,ssn,dssn,dttn);
      
      /* get 1D shape functions for vertical direction */
      apr_shape_fun_1D(pdegz,Eta[2],uun,duun);
      
    }
    else{
      
      /* get 2D shape functions for bases */
      num_shap_2D = apr_shape_fun_2D(base_type_2D,pdegx,Eta,ssn,NULL,NULL);
      
      /* get 1D shape functions for vertical direction */
      apr_shape_fun_1D(pdegz,Eta[2],uun,NULL);
      
    }
    
#ifdef DEBUG_APM
    if(num_shap_2D!=(pdegx+1)*(pdegx+2)/2){
      printf("Wrong number of shape functions %d != %d in shape_fun_3D!\n"
	     ,num_shap_2D,(pdegx+1)*(pdegx+2)/2);
      exit(1);
    } 
#endif
    
/*kbw
printf("In 3D shape: base_type %d, pdeg %d (pdegxy %d, pdegz %d\n", 
Base_type, Pdeg, pdegx, pdegz);
printf("Eta: %lf %lf %lf\n",Eta[0],Eta[1],Eta[2]);
printf("2D: \n");
for(j=0;j<num_shap_2D;j++) {
printf("fun: %lf, fun,x: %lf, fun,y: %lf \n",ssn[j],dssn[j],dttn[j]);
}
printf("1D: \n");
for(j=0;j<pdegz+1;j++) printf("fun: %lf, fun,x: %lf \n",uun[j],duun[j]);
getchar();
/*kew*/


    if (Base_type == APC_BASE_TENSOR_DG) {
	
/* kb-old version 
	k=0;
	for(i=0;i<num_shap_2D;i++){
	  for(j=0;j<pdegz+1;j++){
	    Base_phi_old[k++]=ssn[i]*uun[j];
	  }
	}
	num_shap_old=k;
/* kb-old version */
      
      k=6;
      for(j=0;j<2;j++){
	for(i=3;i<num_shap_2D;i++){
	  Base_phi[k++]=ssn[i]*uun[j];
	}
      }
      for(j=2;j<pdegz+1;j++){
	for(i=0;i<num_shap_2D;i++){
	  Base_phi[k++]=ssn[i]*uun[j];
	}
      }
      num_shap=k;
      
    }
    else if(Base_type == APC_BASE_COMPLETE_DG ) {
      
/* kb-old version 
	ii=0;kk=0;
	for(i=0;i<=Pdeg;i++){
	  for(j=0;j<=Pdeg-i;j++){
	    for(k=0;k<=Pdeg-i-j;k++){
	      Base_phi_old[kk++]=ssn[ii]*uun[k];
	    }
	    ii++; 
	  }
	}
	num_shap_old=kk;
/* kb-old version */
	
	
      kk=8; ii=0;
      for (j=0;j<2;j++) {
	for (i=0;i<2;i++) {
	  for(k=2;k<=Pdeg-i-j;k++){
	    Base_phi[kk++]=ssn[ii]*uun[k];
	  } 
	  ii++;
	} 
      } 
      for (j=0;j<2;j++) {
	for (i=2;i<=Pdeg-j;i++) {
	  for(k=0;k<=Pdeg-i-j;k++){
	    Base_phi[kk++]=ssn[ii]*uun[k];
	  } 
	  ii++;
	}
      } 
      for (j=2;j<=Pdeg;j++) {
	for (i=0;i<=Pdeg-j;i++) {
	  for(k=0;k<=Pdeg-i-j;k++){
	    Base_phi[kk++]=ssn[ii]*uun[k];
	  } 
	  ii++;
	} 
      }
      num_shap=kk;
      
    }
    
    if(Base_dphix!=NULL){
      
#ifdef DEBUG_APM
      if(Base_dphiy==NULL&&Base_dphiz==NULL){
	printf("Error 239483 in shape_fun. Exiting!\n");
	exit(-1);
      }
#endif

      if (Base_type == APC_BASE_TENSOR_DG) {
	  
/* kb-old version 
	  k=0;
	  for(i=0;i<num_shap_2D;i++){
	    for(j=0;j<pdegz+1;j++){
	      Base_dphix_old[k]=dssn[i]*uun[j];
	      Base_dphiy_old[k]=dttn[i]*uun[j];
	      Base_dphiz_old[k]=ssn[i]*duun[j];
	      k++;
	    }
	  }
/* kb-old version */
	  
	k=6;
	for(j=0;j<2;j++){
	  for(i=3;i<num_shap_2D;i++){
	    Base_dphix[k]=dssn[i]*uun[j];
	    Base_dphiy[k]=dttn[i]*uun[j];
	    Base_dphiz[k]=ssn[i]*duun[j];
	    k++;
	  }
	}
	for(j=2;j<pdegz+1;j++){
	  for(i=0;i<num_shap_2D;i++){
	    Base_dphix[k]=dssn[i]*uun[j];
	    Base_dphiy[k]=dttn[i]*uun[j];
	    Base_dphiz[k]=ssn[i]*duun[j];
	    k++;
	  }
	}
	
      }
      else if(Base_type == APC_BASE_COMPLETE_DG ) {
	
/* kb-old version 
	  ii=0;kk=0;
	  for(i=0;i<=Pdeg;i++){
	    for(j=0;j<=Pdeg-i;j++){
	      for(k=0;k<=Pdeg-i-j;k++){
		Base_dphix_old[kk]=dssn[ii]*uun[k];
		Base_dphiy_old[kk]=dttn[ii]*uun[k];
		Base_dphiz_old[kk]=ssn[ii]*duun[k];
		kk++;
	      }
	      ii++;
	    }
	  }
/* kb-old version */
	  
	kk=8; ii=0;
	for (j=0;j<2;j++) {
	  for (i=0;i<2;i++) {
	    for(k=2;k<=Pdeg-i-j;k++){
	      Base_dphix[kk]=dssn[ii]*uun[k];
	      Base_dphiy[kk]=dttn[ii]*uun[k];
	      Base_dphiz[kk]=ssn[ii]*duun[k];
	      kk++;
	    } 
	    ii++;
	  } 
	} 
	for (j=0;j<2;j++) {
	  for (i=2;i<=Pdeg-j;i++) {
	    for(k=0;k<=Pdeg-i-j;k++){
	      Base_dphix[kk]=dssn[ii]*uun[k];
	      Base_dphiy[kk]=dttn[ii]*uun[k];
	      Base_dphiz[kk]=ssn[ii]*duun[k];
	      kk++;
	    } 
	    ii++;
	  }
	} 
	for (j=2;j<=Pdeg;j++) {
	  for (i=0;i<=Pdeg-j;i++) {
	    for(k=0;k<=Pdeg-i-j;k++){
	      Base_dphix[kk]=dssn[ii]*uun[k];
	      Base_dphiy[kk]=dttn[ii]*uun[k];
	      Base_dphiz[kk]=ssn[ii]*duun[k];
	      kk++;
	    } 
	    ii++;
	  } 
	}

      } /* end if set of complete polynomials as shape functions */


    } /* end if derivatives of shape functions required */

/*kbw
printf("3D old:\n");
for(i=0;i<num_shap_old;i++){
  printf("fun - %lf, der: x - %lf, y - %lf, z - %lf\n",
  Base_phi_old[i],Base_dphix_old[i],Base_dphiy_old[i],Base_dphiz_old[i]);
}
printf("3D new:\n");
for(i=0;i<num_shap;i++){
  printf("fun - %lf, der: x - %lf, y - %lf, z - %lf\n",
  Base_phi[i],Base_dphix[i],Base_dphiy[i],Base_dphiz[i]);
}
getchar();
/*kew*/

    return(num_shap);

  } /* end if second or higher order of approximation */


  printf("Wrong parameters for apr_shape_fun_3D !\n");
  exit(1);
 
}


/*---------------------------------------------------------
apr_shape_fun_2D - to compute values of shape functions and their
	local derivatives at a point within the master 2D element
----------------------------------------------------------*/
int apr_shape_fun_2D( /* returns: the number of shape functions */
	int Base_type,	   /* in: type of basis functions: */
			   /*   (APC_BASE_TENSOR_DG) - tensor product */
			   /* 	(APC_BASE_COMPLETE_DG) - complete polynomials */
	int Pdeg, 	   /* in: degree of polynomial - can be either */
			   /*	a single number, for isotropic p, */
			   /*	or a combination pdegy*10+pdegx */
	double *Eta,	   /* in: local coord of the considered point */
	double *Base_phi,  /* out: basis functions */
	double *Base_dphix,/* out: x derivative of basis functions */
	double *Base_dphiy /* out: y derivative of basis functions */
	)
{

/* auxiliary variables */
  int i, j, ii, pdegx, pdegy, pdego, num_shap;
  double  ssn[APC_MAXELVD];
  double  ttn[APC_MAXELVD];
  double dssn[APC_MAXELVD];
  double dttn[APC_MAXELVD];

  int num_shap_old;
  double Base_phi_old[APC_MAXELVD];
  double Base_dphix_old[APC_MAXELVD];
  double Base_dphiy_old[APC_MAXELVD];

/*++++++++++++++++ executable statements ++++++++++++++++*/

  if(Base_phi==NULL) {
    printf("No space provided for shape functions in apr_shape_fun_2D!\n");
    exit(1);
  }
  
  if(Pdeg==0){
/* for piecewise constant approximation */

    Base_phi[0]=1.0;
    
    if(Base_dphix!=NULL){
      Base_dphix[0]=0.0;
      Base_dphiy[0]=0.0;
    }
    
    num_shap=1;
    
    return(num_shap);
    
  }
  else if(Pdeg==10){
/* for linear in y and piecewise constant in x approximation */

    Base_phi[0]=1.0;
    Base_phi[1]=Eta[1];
    
    if(Base_dphix!=NULL){
      
      Base_dphix[0]=0.0;
      Base_dphix[1]=0.0;
      
      Base_dphiy[0]=0.0;
      Base_dphiy[1]=1.0;
      
    }
    
    num_shap=2;
    
    return(num_shap);
    
  }
  else if(Pdeg==1 && Base_type==APC_BASE_TENSOR_DG){
/* for linear in x and piecewise constant in y approximation */

    Base_phi[0]=1.0;
    Base_phi[1]=Eta[0];
    
    if(Base_dphix!=NULL){
      
      Base_dphix[0]=0.0;
      Base_dphix[1]=1.0;
      
      Base_dphiy[0]=0.0;
      Base_dphiy[1]=0.0;
      
    }
    
    num_shap=2;
    
    return(num_shap);
    
  }
  else {
/* for linear and higher order approximation */
    
    Base_phi[0]=1.0;
    Base_phi[1]=Eta[0];
    Base_phi[2]=Eta[1];
  
    if(Base_dphix!=NULL){

      Base_dphix[0]=0.0;
      Base_dphix[1]=1.0;
      Base_dphix[2]=0.0;

      Base_dphiy[0]=0.0;
      Base_dphiy[1]=0.0;
      Base_dphiy[2]=1.0;

    }

    if(Pdeg==1){
/* for linear approximation */
      num_shap=3;
      return(num_shap);
    }

    Base_phi[3]=Eta[0]*Eta[1];
    if(Base_dphix!=NULL){
      Base_dphix[3]=Eta[1];
      Base_dphiy[3]=Eta[0];
    }

    if(Pdeg==11){
/* for bi-linear approximation */
      num_shap=4;
      return(num_shap);
    }

/* take into account different pdeg for anisotropic quads */
    pdego=Pdeg;
    if(Pdeg>10){
      pdegx=Pdeg%10;
      pdegy=Pdeg/10;
      Pdeg=ut_max(pdegx,pdegy);
    }
    else{
      pdegx=Pdeg;
      pdegy=Pdeg;
    }
  
/* get 1D shape functions in directions x and y */
    if(Base_dphix!=NULL){
      apr_shape_fun_1D(pdegx,Eta[0],ssn,dssn);
      apr_shape_fun_1D(pdegy,Eta[1],ttn,dttn);
    }
    else{
      apr_shape_fun_1D(pdegx,Eta[0],ssn,NULL);
      apr_shape_fun_1D(pdegy,Eta[1],ttn,NULL);
    }

    if(pdego<10){

      if(Base_type == APC_BASE_COMPLETE_DG ) {
	
/* kb-old version 
	ii=0;
	for (i=0;i<=Pdeg;i++) {
	  for (j=0;j<=Pdeg-i;j++) {
	    Base_phi_old[ii] = ssn[i]*ttn[j];
	    ii++;
	  } 
	} 
	num_shap_old=ii;
/* kb-old version */

	ii=4;
	for (j=0;j<2;j++) {
	  for (i=2;i<=Pdeg-j;i++) {
	    Base_phi[ii++] = ssn[i]*ttn[j];
	  }
	} 
	for (j=2;j<=Pdeg;j++) {
	  for (i=0;i<=Pdeg-j;i++) {
	    Base_phi[ii++] = ssn[i]*ttn[j];
	  } 
	} 

	num_shap=ii;
	
      }
      else if (Base_type == APC_BASE_TENSOR_DG) {

/* kb-old version 
	ii=0;
	for (i=0;i<=Pdeg;i++) {
	  for (j=0;j<=Pdeg;j++) {
	    Base_phi_old[ii] = ssn[i]*ttn[j];
	    ii++;
	  } 
	} 
	num_shap_old=ii;
/* kb-old version */

	ii=4;
	for (j=0;j<2;j++) {
	  for (i=2;i<=Pdeg;i++) {
	    Base_phi[ii++] = ssn[i]*ttn[j];
	  }
	}
	for (j=2;j<=Pdeg;j++) {
	  for (i=0;i<=Pdeg;i++) {
	    Base_phi[ii++] = ssn[i]*ttn[j];
	  } 
	} 
		
	num_shap=ii;

      }
    }
    else{

/* kb-old version 
      ii=0;
      for (i=0;i<=pdegx;i++) {
	for (j=0;j<=pdegy;j++) {
	  Base_phi_old[ii] = ssn[i]*ttn[j];
	  ii++;
	} 
      }
      num_shap_old=ii;
/* kb-old version */

/* for anisotropic pdeg */
      ii=4;
      for (j=0;j<2;j++) {
	for (i=2;i<=pdegx;i++) {
	  Base_phi[ii++] = ssn[i]*ttn[j];
	}
      }
      for (j=2;j<=pdegy;j++) {
	for (i=0;i<=pdegx;i++) {
	  Base_phi[ii++] = ssn[i]*ttn[j];
	} 
      } 

      num_shap=ii;
      
    }


    if(Base_dphix!=NULL){
      
      if(pdego<10){
	
	if(Base_type == APC_BASE_COMPLETE_DG ) {
	  
/* kb-old version 
	  ii=0;
	  for (i=0;i<=Pdeg;i++) {
	    for (j=0;j<=Pdeg-i;j++) {
	      Base_dphix_old[ii] = dssn[i]*ttn[j];
	      Base_dphiy_old[ii] = ssn[i]*dttn[j];
	      ii++;
	    } 
	  }
/* kb-old version */

	  ii=4;
	  for (j=0;j<2;j++) {
	    for (i=2;i<=Pdeg-j;i++) {
/* local derivatives */
	      Base_dphix[ii] = dssn[i]*ttn[j];
	      Base_dphiy[ii] = ssn[i]*dttn[j];
	      ii++;
	    } /* j */
	  } /* i */
	  for (j=2;j<=Pdeg;j++) {
	    for (i=0;i<=Pdeg-j;i++) {
/* local derivatives */
	      Base_dphix[ii] = dssn[i]*ttn[j];
	      Base_dphiy[ii] = ssn[i]*dttn[j];
	      ii++;
	    } /* j */
	  } /* i */
	  
	}
	else if (Base_type == APC_BASE_TENSOR_DG) {
	  
/* kb-old version 
	  ii=0;
	  for (i=0;i<=Pdeg;i++) {
	    for (j=0;j<=Pdeg;j++) {
	      Base_dphix_old[ii] = dssn[i]*ttn[j];
	      Base_dphiy_old[ii] = ssn[i]*dttn[j];
	      ii++;
	    } 
	  } 
/* kb-old version */

	  ii=4;
	  for (j=0;j<2;j++) {
	    for (i=2;i<=Pdeg;i++) {
/* local derivatives */
	      Base_dphix[ii] = dssn[i]*ttn[j];
	      Base_dphiy[ii] = ssn[i]*dttn[j];
	      ii++;
	    } /* j */
	  } /* i */
	  for (j=2;j<=Pdeg;j++) {
	    for (i=0;i<=Pdeg;i++) {
/* local derivatives */
	      Base_dphix[ii] = dssn[i]*ttn[j];
	      Base_dphiy[ii] = ssn[i]*dttn[j];
	      ii++;
	    } /* j */
	  } /* i */

	}
      }
      else{

/* kb-old version 
	ii=0;
	for (i=0;i<=pdegx;i++) {
	  for (j=0;j<=pdegy;j++) {
	    Base_dphix_old[ii] = dssn[i]*ttn[j];
	    Base_dphiy_old[ii] = ssn[i]*dttn[j];
	    ii++;
	  }
	}
/* kb-old version */

/* for anisotropic pdeg */
	ii=4;
	for (j=0;j<2;j++) {
	  for (i=2;i<=pdegx;i++) {
/* local derivatives */
	    Base_dphix[ii] = dssn[i]*ttn[j];
	    Base_dphiy[ii] = ssn[i]*dttn[j];
	    ii++;
	  }
	}
	for (j=2;j<=pdegy;j++) {
	  for (i=0;i<=pdegx;i++) {
/* local derivatives */
	    Base_dphix[ii] = dssn[i]*ttn[j];
	    Base_dphiy[ii] = ssn[i]*dttn[j];
	    ii++;
	  } 
	} 


      } /* end if anisotropic pdeg */
    
    } /* endif derivatives of shape functions required */

/*kbw
#ifdef DEBUG_APM
printf("3D new:\n");
for(i=0;i<num_shap;i++){
  printf("fun - %lf, der: x - %lf, y - %lf\n",
  Base_phi[i],Base_dphix[i],Base_dphiy[i]);
}
printf("3D old:\n");
for(i=0;i<num_shap_old;i++){
  printf("fun - %lf, der: x - %lf, y - %lf\n",
  Base_phi_old[i],Base_dphix_old[i],Base_dphiy_old[i]);
}
getchar();
#endif
/*kew*/

    return(num_shap);

  } /* end if linear and higher order approximation */

  printf("Wrong parameters for apr_shape_fun_2D !\n");
  exit(1);
 
}

/*---------------------------------------------------------
apr_shape_fun_1D - to compute values of 1D shape functions and their
	local derivatives at a point within the master element
----------------------------------------------------------*/
int apr_shape_fun_1D( /* returns: the number of shape functions */
	int Pdeg, 	   /* in: degree of polynomial */
	double Eta,	   /* in: local coord of the considered point */
	double *Base_phi,  /* out: basis functions */
	double *Base_dphix/* out: x derivative of basis functions */
	)
{

/* auxiliary variables */
  int i;

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef DEBUG_APM
  if(Base_phi==NULL) {
    printf("No space provided for shape functions in apr_shape_fun_1D!\n");
    exit(1);
  }
#endif
  
/* monomial basis */
  Base_phi[0] = 1.0;
  if(Pdeg>0) Base_phi[1] = Eta;
  for (i=2;i<=Pdeg; i++) {
    Base_phi[i]  = Eta*Base_phi[i-1];
  }   

  if(Base_dphix!=NULL){
    Base_dphix[0]  = 0.0;
    if(Pdeg>0) Base_dphix[1]  = 1.0;
    for (i=2;i<=Pdeg; i++) {
      Base_dphix[i] = Eta*Base_dphix[i-1];
    }   
    for (i=2;i<=Pdeg; i++) {
      Base_dphix[i] *= i;
    }   
  } /* end if derivatives of 1D shape functions required */

  return(Pdeg+1);
}

/*---------------------------------------------------------
apr_set_quadr_3D - to prepare quadrature data for a given element
---------------------------------------------------------*/
int apr_set_quadr_3D(
	int Base_type,	   /* in: type of basis functions: */
			   /* 	(APC_BASE_TENSOR_DG) - tensor product */
			   /* 	(APC_BASE_COMPLETE_DG) - complete polynomials */
	int *Pdeg_vec,	/* in: element degree of polynomial */
	int *Ngauss,	/* out: number of gaussian points */
	double *Xg,	/* out: coordinates of gaussian points */
	double *Wg	/* out: weights associated with points */
	)
{

/* local variables */
  int pdeg,orderx,orderz; 	/* orders of approximation in x,y,z */
  int ngauss2,ngaussz;   /* numbers of gaussian points */
  double *xg2,*xgz;/* gauss points in 2D and 1D*/
  double *wg2,*wgz; /* gauss weights in 2D and 1D*/
  int iaux,ki,kj;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  pdeg = Pdeg_vec[0];

  if(Base_type==APC_BASE_TENSOR_DG){
  
/* decipher horizontal and vertical orders of approximation */
    orderz=2*(pdeg/100);
    orderx=2*(pdeg%100);
    if(orderz==0) orderz=1;
    if(orderx==0) orderx=1;

  }
  else if(Base_type==APC_BASE_COMPLETE_DG){

    orderz=2*pdeg;
    orderx=2*pdeg;
    if(orderz==0) orderz=1;
    if(orderx==0) orderx=1;

  }
   
  iaux=1; /* 1D integration for vertical direction */
  apr_gauss_select(iaux,orderz,&ngaussz,&xgz,&wgz);
  iaux=2; /* 2D integration for triangular bases */
  apr_gauss_select(iaux,orderx,&ngauss2,&xg2,&wg2);

/*kbw
#ifdef DEBUG_APM
  printf("selecting gauss for orderz %d and orderx %d\n",
	 orderz, orderx);
  for(ki=0;ki<ngaussz;ki++){
    printf("point %d, weight %lf, coor %lf\n",
	   ki, wgz[ki], xgz[ki]);
  }
  for(ki=0;ki<ngauss2;ki++){
    printf("point %d, weight %lf, coor %lf %lf\n",
	   ki, wg2[ki], xg2[3*ki+1], xg2[3*ki+2]);
  }
#endif
/*kew*/


/* set the total number of Gauss points */
  *Ngauss=ngaussz*ngauss2;

/* fill the arrays of Gauss points and Gauss coefficients */
  for(ki=0;ki<ngaussz;ki++){
    for(kj=0;kj<ngauss2;kj++){

      Xg[3*(ki*ngauss2+kj)] = xg2[3*kj+1];
      Xg[3*(ki*ngauss2+kj)+1] = xg2[3*kj+2];
      Xg[3*(ki*ngauss2+kj)+2] = xgz[ki];
/* we correct weights for 2D integration on triangles */
      Wg[ki*ngauss2+kj] = 0.5*wg2[kj]*wgz[ki];

    }
  }


  return(1);
}

/*---------------------------------------------------------
apr_set_quadr_2D - to prepare quadrature data for a given face
---------------------------------------------------------*/
int apr_set_quadr_2D(
	int Fa_type,	/* in: type of a face */
	int Base_type,	   /* in: type of basis functions: */
			   /* 	(APC_BASE_TENSOR_DG) - tensor product */
			   /* 	(APC_BASE_COMPLETE_DG) - complete polynomials */
	int *Pdeg_vec,	/* in: element degree of polynomial */
	int *Ngauss,	/* out: number of gaussian points */
	double *Xg,	/* out: coordinates of gaussian points */
	double *Wg	/* out: weights associated with points */
	)
{

/* local variables */
  int pdeg,orderx,orderz; 	/* orders of approximation in x,y,z */
  int ngauss2,ngaussx,ngaussy;/* numbers of gaussian points */
  double *xg2,*xgx,*xgy;/* gauss points in 2D and 1D*/
  double *wg2,*wgx,*wgy; /* gauss weights in 2D and 1D*/
  int iaux,ki,kj;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  pdeg = Pdeg_vec[0];
  
  if(Base_type==APC_BASE_TENSOR_DG){
  
/* decipher horizontal and vertical orders of approximation */
    orderz=2*(pdeg/100);
    orderx=2*(pdeg%100);
    if(orderz==0) orderz=1;
    if(orderx==0) orderx=1;

  }
  else if(Base_type==APC_BASE_COMPLETE_DG){

    orderz=2*pdeg;
    orderx=2*pdeg;
    if(orderz==0) orderz=1;
    if(orderx==0) orderx=1;

  }
   
/* fill the arrays of Gauss points and Gauss coefficients */

  if (Fa_type == MMC_QUAD) {

    iaux=1; /* 1D line integration */
    apr_gauss_select(iaux,orderx,&ngaussx,&xgx,&wgx);
    apr_gauss_select(iaux,orderz,&ngaussy,&xgy,&wgy);
    *Ngauss = ngaussx*ngaussy;

    if(Xg!=NULL){

      for (ki=0;ki<ngaussx;ki++) {
        for (kj=0;kj<ngaussy;kj++) {

          Xg[2*(ki*ngaussy+kj)] = xgx[ki];
          Xg[2*(ki*ngaussy+kj)+1] = xgy[kj];
/* we correct weights because master element side is from 0 to 1 */
/* but Xg is unchanged so apr_fa_elem_coor must take care of this */
          Wg[ki*ngaussy+kj] = 0.5*wgx[ki]*wgy[kj];

        }
      }

    }
  }
  else if (Fa_type == MMC_TRIA) {

    iaux=2; /* 2D integration for triangles */
    apr_gauss_select(iaux,orderx,&ngauss2,&xg2,&wg2);
    *Ngauss = ngauss2;

    if(Xg!=NULL){

      for (ki=0;ki<ngauss2;ki++) {
        Xg[2*ki] = xg2[3*ki+1]; 
        Xg[2*ki+1] = xg2[3*ki+2];
/* we correct weights for 2D integration on triangles */
        Wg[ki] = 0.5*wg2[ki];
      }

    }
  }


  return(1);
}

/*---------------------------------------------------------
apr_L2_proj - to L2 project a function onto an element
---------------------------------------------------------*/
int apr_L2_proj(	/* returns: >0 - success, <=0 - failure */
        int Field_id,   /* in: field ID */
	int Mode,	/* in: mode of operation */
			/*    <-1 - projection from ancestor to father */
			/*          the value is the number of ancestors */
			/*     -1 - projection from father to son */
                        /*     >0 - projection of function, the routine */
                        /*          returning value at point is specified */
			/*          by the pointer in the last argument */
	int El,		/* in: element number */
	int *Pdeg_vec,	/* in: element degree of approximation */
	double* Dofs,	/* out: workspace for degress of freedom of El */
			/* 	NULL - write to  data structure */
	int* El_from,	/* in: list of elements to provide function */
	int *Pdeg_vec_from,	/* in: degree of polynomial for each El_from */
	double* Dofs_from, /* in: Dofs of El_from or...*/
        double (*Fun_p)(double*,double*,double*,double*)   /* in: pointer to */
	                   /* function with field values and its derivatives */
	)
{

  /* pointer to field structure */
  apt_field *field_p;
  int mesh_id;

/* local variables */
  int pdeg;             /* element degree of approximation - one number in DG */
  int base_q;		/* type of basis functions for quadrilaterals */
  int nreq;		/* number of equations (solution components)  */
  int num_shap;         /* number of element shape functions */
  int ndofs;            /* local dimension of the problem */
  int ngauss;           /* number of gauss points */
  double xg[3000];   	/* coordinates of gauss points in 3D */
  double wg[1000];      /* gauss weights */
  double determ;        /* determinant of jacobi matrix */
  double vol;           /* volume for integration rule */
  double xcoor[3];      /* global coord of gauss point */
  double xloc[3];       /* local point coordinates */
  double base_phi[APC_MAXELVD];    /* basis functions */
  double base_dphix[APC_MAXELVD];  /* x-derivatives of basis function */
  double base_dphiy[APC_MAXELVD];  /* y-derivatives of basis function */
  double base_dphiz[APC_MAXELVD];  /* y-derivatives of basis function */
  double base_aux[APC_MAXELVD];    /* basis functions */
  int el_nodes[MMC_MAXELVNO+1];      /* list of nodes of El */
  double node_coor[3*MMC_MAXELVNO];  /* coord of nodes of El */
  double stiff_loc[APC_MAXELSD*APC_MAXELSD];
				/* stiffness matrix for local problem */
  double f_loc[APC_MAXELSD]; /* rhs vector for local problem */
  int ipiv[APC_MAXELSD];
  double value[APC_MAXEQ]; /* solution vector to be projected */
  double u_val[APC_MAXEQ]; /* projected solution vector */

/* auxiliary variables */
  int i, ki, kk, iaux, iel, iel_from, ifound;
  int idofs, jdofs, ieq;
  double daux, time;

/* constants */
  int ione = 1;	
/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);

  /* select the corresponding mesh */
  mesh_id = apr_get_mesh_id(Field_id);

  pdeg = Pdeg_vec[0];

/* get formulation parameters */
  base_q=apr_get_base_type(Field_id, El);
  nreq = apr_get_nreq(Field_id);

/* find number of element shape functions and scalar dofs */
  num_shap = apr_get_el_pdeg_numshap(Field_id, El, &pdeg);
  ndofs = num_shap*nreq;

/* initialize the matrices to zero */
  for(i=0;i<ndofs*ndofs;i++) stiff_loc[i]=0.0;
  for(i=0;i<ndofs;i++) f_loc[i]=0.0; 

/* get the coordinates of the nodes of El in the right order */
  mmr_el_node_coor(mesh_id,El,el_nodes,node_coor);

/* prepare data for gaussian integration */
  apr_set_quadr_3D(base_q, &pdeg, &ngauss, xg, wg);

/*kbw
//if(Mode<-1){
printf("In l2_proj for element %d\n",El);
printf("pdeg %d, ngauss %d\n",pdeg,ngauss);
printf("nreq %d, ndof %d, local_dim %d\n",nreq,num_shap,ndofs);
printf("%d nodes with coordinates:\n",el_nodes[0]);
for(i=0;i<el_nodes[0];i++){
  printf("node %d (global - %d): x - %f, y - %f, y - %f\n", i, el_nodes[i+1],
	node_coor[3*i],node_coor[3*i+1],node_coor[3*i+2]);
}
getchar();
//}
/*kew*/

  for (ki=0;ki<ngauss;ki++) {

/* at the gauss point, compute basis functions, determinant etc*/
      iaux = 2; /* calculations with jacobian but not on the boundary */
      determ = apr_elem_calc_3D(iaux, nreq, &pdeg, base_q,
		&xg[3*ki],node_coor,NULL,
  		base_phi,base_dphix,base_dphiy,base_dphiz,
		xcoor,NULL,NULL,NULL,NULL,NULL);

      vol = determ * wg[ki];

/*kbw
if(Mode<-1){
printf("at gauss point %d, local coor %lf, %lf, %lf\n", 
ki,xg[3*ki],xg[3*ki+1],xg[3*ki+2]);
printf("global coor %lf %lf %lf\n",xcoor[0],xcoor[1],xcoor[2]);
printf("weight %lf, determ %lf, coeff %lf\n",
	wg[ki],determ,vol);
getchar();
}
/*kew*/

/* get value of projected function at integration point */
      if(Mode>0){
	//	(*Fun_p)(xcoor, fun_val, fun_dx, fun_dy, fun_dz);
      }
      else if(Mode<=-1){

/* we assume all elements are of the same type and have the same
basis functions */

/* initiate counter for Dofs */
        idofs=0;

/* for each element on a list */
        for(iel=0;iel<-Mode;iel++){
          iel_from=El_from[iel];

/* find local coordinates within the element */
          ifound = mmr_loc_loc(Field_id,El,&xg[3*ki],iel_from,xloc);

/*kbw
if(Mode<-1){
if(ifound==1) printf("in element %d, pdeg %d, idofs %d, at local coor %lf, %lf, %lf\n", 
iel_from,Pdeg_vec_from[iel][0], idofs,xloc[0],xloc[1],xloc[2]);
else printf("in element %d, pdeg %d, idofs %d - not found\n",
iel_from,Pdeg_vec_from[iel][0], idofs);
printf("Dofs_from:\n");
for (i=0;i<ndofs;i++) {
  printf("%20.15lf",Dofs_from[idofs+i]);
}
printf("\n");
}
/*kew*/


          if(ifound==1){

/* find value of function at integration point */
            iaux = 1; /* calculations without jacobian */
            apr_elem_calc_3D(iaux, nreq, &Pdeg_vec_from[iel], 
                base_q,xloc,NULL,&Dofs_from[idofs],
		base_aux,NULL,NULL,NULL,
		NULL,value,NULL,NULL,NULL,NULL);

            break;
          }
          
/* update counter for dofs */
          idofs += nreq*apr_get_el_pdeg_numshap(Field_id, iel, 
						&Pdeg_vec_from[iel]);

        }

      }
      else {
        printf("Unknown mode for L2 projection!\n");
        exit(-1);
      }
    
/*kbw
if(Mode<-1){
  printf("found value(s) for projection: ");
  for(ieq=0;ieq<nreq;ieq++){
    printf("%20.15lf", value[ieq]);
  }
  printf("\n");
}
/*kew*/

      for(ieq=0;ieq<nreq;ieq++){

        kk=(ieq*ndofs+ieq)*num_shap;

        for (jdofs=0;jdofs<num_shap;jdofs++) {
          for (idofs=0;idofs<num_shap;idofs++) {

         
/* mass matrix */
            stiff_loc[kk+idofs] += 
                     base_phi[jdofs] * base_phi[idofs] * vol;

          }/* idofs */
          kk+=ndofs;

        } /* jdofs */
      }/* ieq */

      kk=0;
      for(ieq=0;ieq<nreq;ieq++){
        for (idofs=0;idofs<num_shap;idofs++) {

/* right hand side vector */
        f_loc[kk] += value[ieq] * base_phi[idofs] * vol;
        kk++;

        }/* idofs */
      }/* ieq */

  } /* ki */

/*kbw
if(El>0){
printf("Stiffness matrix:\n");
for (idofs=0;idofs<ndofs;idofs++) {
  for (jdofs=0;jdofs<ndofs;jdofs++) {
    printf("%20.15lf",stiff_loc[idofs+jdofs*ndofs]);        
  } 
  printf("\n");
}
printf("F_loc:\n");
for (idofs=0;idofs<ndofs;idofs++) {
  printf("%20.15lf",f_loc[idofs]);
}
printf("\n");
getchar(); 
} 
/*kew*/

/* solve the local problem */
  dgetrf_(&ndofs,&ndofs,stiff_loc,&ndofs,ipiv,&iaux);
   mf_check(iaux == 0, "DGETRF exited with error %d",iaux);
  dgetrs_("N",&ndofs,&ione,stiff_loc,&ndofs,ipiv,f_loc,&ndofs,&iaux);

/*kbw
if(El>0){
printf("In l2_proj for element %d, ",El);
printf("Solution:\n");
for (idofs=0;idofs<ndofs;idofs++) {
  printf("%20.15lf",f_loc[idofs]);
}
printf("\n");
getchar(); 
} 
kew*/

/* rewrite solution */
  if(Dofs==NULL){
/* to global data structure */
    ndofs = apr_get_ent_nrdofs(Field_id,APC_ELEMENT,El);
    iaux=1; /* current solution dofs_1 */
    apr_write_ent_dofs(Field_id,APC_ELEMENT,El,ndofs,iaux,f_loc);
  }
  else {
/* to vector Dofs */
    for(i=0;i<ndofs; i++) Dofs[i]=f_loc[i];
  }

/* check for Mode = -1 */
#ifdef DEBUG_APM
if(Mode == -1 && Dofs!=NULL){

/* in a loop over integration points */
  for (ki=0;ki<ngauss;ki++) {

/* at a gauss point compute solution etc*/
    iaux = 1; /* calculations without jacobian */
    apr_elem_calc_3D(iaux, nreq, &pdeg, base_q,
		&xg[3*ki],NULL,Dofs,
		base_phi,NULL,NULL,NULL,
		NULL,u_val,NULL,NULL,NULL,NULL);


/* get solution in element_from */
    iel_from=*El_from;

/* find local coordinates within the element */
    ifound = mmr_loc_loc(Field_id,El,&xg[3*ki],iel_from,xloc);

    if(ifound==1){

/* find value of function at integration point */
      iaux = 1; /* calculations without jacobian */
      apr_elem_calc_3D(iaux, nreq, &Pdeg_vec_from[0], 
                base_q,xloc,NULL,Dofs_from,
		base_phi,NULL,NULL,NULL,
		NULL,value,NULL,NULL,NULL,NULL);

    }
    else {
      printf("Something wrong in checking L2 proj!\n");
      exit(-1);
    }

/*kbw
if(El>0){
printf("Checking l2_proj from element %d (pdeg %d to element %d, (pdeg %d)",
*El_from,Pdeg_vec_from[0],El,pdeg);
printf("Dofs_from:\n");
for (idofs=0;idofs<ndofs;idofs++) {
  printf("%20.15lf",Dofs_from[idofs]);
}
printf("\n");
printf("Dofs_to:\n");
for (idofs=0;idofs<ndofs;idofs++) {
  printf("%20.15lf",Dofs[idofs]);
}
printf("\n");
printf("Data:\n");
for(ieq=0;ieq<nreq;ieq++){
  printf("%20.15lf",value[ieq]);
}
printf("\n");
printf("Solution:\n");
for(ieq=0;ieq<nreq;ieq++){
  printf("%20.15lf",u_val[ieq]);
}
printf("\n");
getchar(); 
} 
/*kew*/

    for(ieq=0;ieq<nreq;ieq++){

      if(fabs(value[ieq]-u_val[ieq])>SMALL){
        printf("Something wrong in checking L2 proj! %lf != %lf\n",
		value[ieq],u_val[ieq]);
        exit(-1);
      }
    }

  }
}
#endif

  return(1);
}



/*---------------------------------------------------------
apr_sol_xglob - to return the solution at a point with global
	coordinates specified. The procedure finds,
	for a given global point, the respective local coordinates
	within the proper initial mesh element, than computes 
	corresponding local coordinates within an active ancestor of 
	the initial element and finally finds the value of solution. 
	There may be several initial mesh elements, several ancestors
	and several values due to the discontinuity of approximate solution.
---------------------------------------------------------*/
int apr_sol_xglob( /* returns: >0 - success, <=0 - failure */
        int Field_id,   /* in: field ID */
        double *Xglob,	/* in: global coordinates of a point */
    int Nb_sol,    /* in: which solution to take: 1 - sol_1, 2 - sol_2 */
    int* El,	/* out: list of element numbers,  */
                    /*      El[0] - number of elements on the list */
    double* Xloc,	/* out: list of local coordinates within elements */
        double *Sol,	/* out: list of solutions at the point */
        double *Dxsol,  /* out: list of derivatives wrt x of solution */
        double *Dysol,  /* out: list of derivatives wrt y of solution */
        double *Dzsol,   /* out: list of derivatives wrt z of solution */
        double Close_proximity_precision, //< default to APC_CLOSE; how far from element point is considered to be within an element
        int Sol_xglob_flags //< see ape_sol_xlgob
    )
{


  /* pointer to field structure */
  apt_field *field_p;
  int mesh_id;

  int i, iel, iaux, jaux, row, col, isys, pdeg, nrel, el_aux;
  int nmel, ndofs, base_q, ifound, nreq;

  int ipiv[4], el_son, ison, elsons[MMC_MAXELSONS+1];
  int nrel_check, el_check[10*MMC_MAXELSONS+1];
  double amat[16], bvec[4], atrans[9], btrans[3], daux, faux;

  double dofs_loc[APC_MAXELSD]; 
  double base_phi[APC_MAXELVD];   
  double base_dphix[APC_MAXELVD];  /* x-derivatives of basis function */
  double base_dphiy[APC_MAXELVD];  /* y-derivatives of basis function */
  double base_dphiz[APC_MAXELVD];  /* y-derivatives of basis function */
  double node_coor[3*MMC_MAXELVNO];
  double xloc[10*MMC_MAXELSONS*3], xloc_aux[3], xcoor[3];
  double ref_points_prism[6][3] = 
/* data for reference points for three dimensional prismatic element */
	{ { 0.0, 0.0, -1.0 }, { 1.0, 0.0, -1.0 }, { 0.0, 1.0, -1.0 }, 
	  { 0.0, 0.0,  1.0 }, { 1.0, 0.0,  1.0 }, { 0.0, 1.0,  1.0 } };

  int ione=1;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);

  /* get the corresponding mesh ID */
  mesh_id = field_p->mesh_id;

/* get necessary mesh parameters */
  nmel=mmr_get_max_elem_id(mesh_id);

/* get formulation parameters */
  nreq = apr_get_nreq(Field_id);

/* prepare the matrix for finding coefficients of linear transformation */

  nrel=0;

  int first_id = 1;
  if(Sol_xglob_flags & APE_SOL_XGLOB_CHECK_ONLY_GIVEN_ELEMENT) {
      first_id = El[0];
      nmel = El[0];
  }

/* loop over all initial mesh elements */
  for (iel=first_id;iel<=nmel;iel++) {
    if (mmr_el_status(mesh_id,iel)!=MMC_FREE &&
	mmr_el_gen(mesh_id,iel)==0) {


/* get the coordinates of element nodes in the right order */
      mmr_el_node_coor(mesh_id,iel,NULL,node_coor);

/*kbw
if(iel>0){
printf("In sol_xglob for element %d\n",iel);
printf("nodes with coordinates:\n");
for(i=0;i<6;i++){
  printf("node %d : x - %f, y - %f, y - %f\n", i, 
	node_coor[3*i],node_coor[3*i+1],node_coor[3*i+2]);
}
getchar();
}
kew*/


/*kb!!!
we assume that the geometrical transformation from the master element 
is linear
kb!!!*/

/* find coefficients of the transformation by solving a system of
linear equations */

        ndofs = 4; /* four unknown coefficients of the transformation */

        i=0;
        for(col=0;col<3;col++){
          for(row=0;row<4;row++){
            amat[i++] = node_coor[3*row+col];
          }
        }

        for(row=0;row<4;row++){
          amat[i++] = 1.0;
        }

/*kbw
printf("before decomposition amat: \n");
for(i=0;i<(ndofs)*(ndofs);i++) printf("%20.15lf",amat[i]) ; 
printf("\n");
kew*/

        dgetrf_(&ndofs,&ndofs,amat,&ndofs,ipiv,&iaux);
         mf_check(iaux == 0, "DGETRF exited with error %d",iaux);

/*kbw
printf("after decomposition amat: \n");
for(i=0;i<(ndofs)*(ndofs);i++) printf("%20.15lf",amat[i]) ; 
printf("\n");
kew*/

        for(isys=0;isys<3;isys++){

          for(row=0;row<4;row++){
            bvec[row] = ref_points_prism[row][isys];
          }

/*kbw
printf("bvec: \n");
for(i=0;i<4;i++) printf("%20.15lf",bvec[i]) ; 
printf("\n");
/*kew*/

/* solve the local problem */
          dgetrs_("N",&ndofs,&ione,amat,&ndofs,ipiv,bvec,&ndofs,&iaux);
          for(i=0;i<3;i++) atrans[3*isys+i]=bvec[i] ; 
          btrans[isys] = bvec[3];

/*kbw
printf("solution: \n");
for(i=0;i<ndofs;i++) printf("%20.15lf",bvec[i]) ; 
printf("\n");
/*kew*/

    
/*kbw
          for(row=0;row<3;row++){

            daux=0;
            for(col=0;col<3;col++){

              daux += atrans[3*isys+col]*node_coor[3*row+col];

            }

            daux += btrans[isys];

            printf("Checking: %20.15lf ?= %20.15lf\n", 
		daux,ref_points_prism[row][isys]);

          }
/*kew*/

        } /* end of loop over systems */

/* checking remaining element vertices to confirm the transformation
is linear */
/*kbw
printf("Vertex 5: \n");
        for(row=0;row<3;row++){

          daux=0;
          for(col=0;col<3;col++){

            daux += atrans[3*row+col]*node_coor[12+col];

          }

          daux += btrans[row];

          printf("Checking coordinate %d : %20.15lf ?= %20.15lf\n", 
		row, daux,ref_points_prism[4][row]);

        }
printf("Vertex 6: \n");
        for(row=0;row<3;row++){

          daux=0;
          for(col=0;col<3;col++){

            daux += atrans[3*row+col]*node_coor[15+col];

          }

          daux += btrans[row];

          printf("Checking coordinate %d : %20.15lf ?= %20.15lf\n", 
		row, daux,ref_points_prism[5][row]);

        }
/*kew*/


        ut_mat3vec(atrans,Xglob,xloc);
        for(i=0;i<3;i++) xloc[i] += btrans[i];

/*checking*/
#ifdef DEBUG_APM
  iaux = 1; /* calculations without jacobian */
  apr_elem_calc_3D(iaux, nreq, &iaux, iaux,
		xloc,node_coor,NULL,
		NULL,NULL,NULL,NULL,
		xcoor,NULL,NULL,NULL,NULL,NULL);
  for(i=0;i<3;i++) {
    if(fabs(Xglob[i]-xcoor[i])>APC_SMALL){
      printf("Error in finding local coordinates for global point\n");
      exit(-1);
    }
  }
/*kbw
printf("Checking :\n");
for(i=0;i<3;i++) 
    printf("coord %d, local %15.12lf, global %15.12lf ?= %15.12lf\n", 
		i,xloc[i],Xglob[i],xcoor[i]);
/*kew*/
#endif

/* if global point within initial mesh element */
        if( xloc[0] > -Close_proximity_precision && xloc[1] > -Close_proximity_precision &&
          xloc[0]+xloc[1] < 1.0+Close_proximity_precision &&
          xloc[2] > -1.0-Close_proximity_precision && xloc[2] < 1.0+Close_proximity_precision ) {

	  if(xloc[0] < 0.0) xloc[0]=0.0;
	  if(xloc[1] < 0.0) xloc[1]=0.0;
	  if(xloc[0]+xloc[1] > 1.0) xloc[0]=1.0-xloc[1];
	  if(xloc[2] <-1.0) xloc[2]=-1.0;
	  if(xloc[2] > 1.0) xloc[2]=1.0;

/*kbw
printf("Found respective point in initial element %d\n local coordinates:",iel);
for(i=0;i<3;i++) printf("%15.12lf",xloc[i]);printf("\n");
/*kew*/
	  nrel_check=1; el_check[0]=iel;

	  while(nrel_check>0){

	    nrel_check--;
	    el_aux = el_check[nrel_check];

/*kbw
printf("checking element %d, waiting list:",el_aux);
for(i=0;i<nrel_check;i++) printf(" %d",el_check[i]);
printf("\n");
/*kew*/

	    while( mmr_el_status(mesh_id,el_aux)<=0){

	      for(i=0;i<3;i++) xloc_aux[i] = xloc[3*nrel_check+i];

	      mmr_el_fam(mesh_id,el_aux,elsons,NULL);

/*find proper son and its local coordinates*/
	      ifound=0;
          if(xloc_aux[2]<Close_proximity_precision){
           if(xloc_aux[0]>0.5-Close_proximity_precision) {
		el_check[nrel_check]=elsons[2];
		ifound += mmr_loc_loc(Field_id,el_aux,xloc_aux,
				  el_check[nrel_check],&xloc[3*nrel_check]);
		nrel_check++;
	       }
           if(xloc_aux[1]>0.5-Close_proximity_precision) {
		el_check[nrel_check]=elsons[3];
		ifound += mmr_loc_loc(Field_id,el_aux,xloc_aux,
				  el_check[nrel_check],&xloc[3*nrel_check]);
		nrel_check++;
	       }
           if(xloc_aux[0]+xloc_aux[1]<0.5+Close_proximity_precision)  {
		el_check[nrel_check]=elsons[1];
		ifound += mmr_loc_loc(Field_id,el_aux,xloc_aux,
				  el_check[nrel_check],&xloc[3*nrel_check]);
		nrel_check++;
	       }
           if(xloc_aux[0]+xloc_aux[1]>0.5-Close_proximity_precision &&
         xloc_aux[0]<0.5+Close_proximity_precision &&
         xloc_aux[1]<0.5+Close_proximity_precision )  {
		el_check[nrel_check]=elsons[4];
		ifound += mmr_loc_loc(Field_id,el_aux,xloc_aux,
				  el_check[nrel_check],&xloc[3*nrel_check]);
		nrel_check++;
	       }
	      }
          if(xloc_aux[2]>-Close_proximity_precision){
           if(xloc_aux[0]>0.5-Close_proximity_precision) {
		el_check[nrel_check]=elsons[6];
		ifound += mmr_loc_loc(Field_id,el_aux,xloc_aux,
				  el_check[nrel_check],&xloc[3*nrel_check]);
		nrel_check++;
	       }
           if(xloc_aux[1]>0.5-Close_proximity_precision) {
		el_check[nrel_check]=elsons[7];
		ifound += mmr_loc_loc(Field_id,el_aux,xloc_aux,
				  el_check[nrel_check],&xloc[3*nrel_check]);
		nrel_check++;
	       }
           if(xloc_aux[0]+xloc_aux[1]<0.5+Close_proximity_precision)  {
		el_check[nrel_check]=elsons[5];
		ifound += mmr_loc_loc(Field_id,el_aux,xloc_aux,
				  el_check[nrel_check],&xloc[3*nrel_check]);
		nrel_check++;
	       }
           if(xloc_aux[0]+xloc_aux[1]>0.5-Close_proximity_precision &&
         xloc_aux[0]<0.5+Close_proximity_precision &&
         xloc_aux[1]<0.5+Close_proximity_precision )  {
		el_check[nrel_check]=elsons[8];
		ifound += mmr_loc_loc(Field_id,el_aux,xloc_aux,
				  el_check[nrel_check],&xloc[3*nrel_check]);
		nrel_check++;
	       }
	      }

	      if(ifound<=0){
		printf("Local coordinates not found within a family\n");
		exit(-1);
	      }

/*kbw
printf("after checking family for element %d, waiting list:\n",el_aux);
for(i=0;i<nrel_check;i++) {
printf("element %d, xcoor %lf %lf %lf\n",
el_check[i],xloc[3*i],xloc[3*i+1],xloc[3*i+2]);
}
getchar();
/*kew*/

	      nrel_check--;
	      el_aux=el_check[nrel_check];
	  
	    } /* end while elements inactive */

/* found active element with global coordinates Xglob */
	    nrel++;
	    El[nrel]=el_aux;
	    for(i=0;i<3;i++) Xloc[(nrel-1)*3+i] = xloc[3*nrel_check+i];

/*kbw
printf("found next (%d) element %d, coor %lf %lf %lf\n",
       nrel,El[nrel],Xloc[(nrel-1)*3],Xloc[(nrel-1)*3+1],Xloc[(nrel-1)*3+2]);

getchar();
/*kew*/

/* if solution required */
	    if(Sol!=NULL){

	      base_q=apr_get_base_type(Field_id, el_aux);

/* find degree of polynomial and number of element scalar dofs */
	      pdeg = apr_get_ent_pdeg(Field_id, APC_ELEMENT, el_aux);

/* get the coordinates of the nodes of El in the right order */
	      mmr_el_node_coor(mesh_id,el_aux,NULL,node_coor);

/* get the most recent solution degrees of freedom */
	      ndofs = apr_get_ent_nrdofs(Field_id,APC_ELEMENT,el_aux);
	      apr_read_ent_dofs(Field_id,APC_ELEMENT,el_aux,ndofs,Nb_sol,
				dofs_loc);

	      if(Dxsol!=NULL){

/* find value of solution at a point */
		iaux = 2; /* calculations with jacobian */
		jaux = (nrel-1)*nreq;
		apr_elem_calc_3D(iaux, nreq, &pdeg, base_q,
				  &xloc[3*nrel_check],node_coor,dofs_loc,
				  base_phi,base_dphix,base_dphiy,base_dphiz,
				  xcoor,&Sol[jaux],
				  &Dxsol[jaux],&Dysol[jaux],&Dzsol[jaux],NULL);

	      }
	      else{

/* find value of solution at a point */
		iaux = 1; /* calculations without jacobian */
		jaux = (nrel-1)*nreq;
		apr_elem_calc_3D(iaux, nreq, &pdeg, base_q,
				  &xloc[3*nrel_check],node_coor,dofs_loc,
				  base_phi,NULL,NULL,NULL,
				  xcoor,&Sol[jaux],NULL,NULL,NULL,NULL);

	      }

/*kbw
{
double f, fx,fy,fz, lalp;

printf("\nChecking :\n");
for(i=0;i<3;i++) 
    printf("coord %d, local %15.12lf, global %15.12lf ?= %15.12lf\n", 
		i,xloc[3*nrel_check+i],Xglob[i],xcoor[i]);

dgpr_exact_sol(name,xcoor[0],xcoor[1],xcoor[2],&f,&fx,&fy,&fz,&lalp);
printf("Solution (%d) : exact %20.15lf, approximate %20.15lf\n",
	jaux,f,Sol[jaux]);
if(Dxsol!=NULL){
printf("x derivative  : exact %20.15lf, approximate %20.15lf\n",
	fx,Dxsol[jaux]);
printf("y derivative  : exact %20.15lf, approximate %20.15lf\n",
	fy,Dysol[jaux]);
printf("z derivative  : exact %20.15lf, approximate %20.15lf\n",
	fz,Dzsol[jaux]);
}
getchar();
}

/*kew*/

/*checking*/
#ifdef DEBUG_APM
for(i=0;i<3;i++) 
  if(fabs(Xglob[i]-xcoor[i])>APC_SMALL){
    printf("Error in finding local coordinates for global point\n");
    exit(-1);
  }
#endif

	    } /* end rewriting solution */

	  } /* end different branches from initial mesh element */

	} /* end if point found within initial mesh lement */


    } /* end if element in initial mesh */
  } /* end of loop over elements */

  El[0]=nrel;

  return(1);
}


int apr_get_profile( // returns: >0 - success, <=0 - failure
  FILE* filePtr,//in: pointer to file, where profile will be printed
  int fieldId,		//in: field ID
  int solNr,			//in: solution number
  int nSol,			//in: solution length (no. of components)
  double * pt1,		//in: begin point coords
  double * pt2,		//in: end point coords
  int nPoints			//in: number of points (including end points)
		     )
{
  if(filePtr != NULL
     && fieldId>0
     && solNr>0
     && nSol>0
     && pt1!=NULL
     && pt2!=NULL
     && nPoints>1)
    {
      double nDp=nPoints-1;
      int p=0,s=0;
      int el[30]={0};
      double dp[3]={(pt2[0]-pt1[0])/nDp,(pt2[1]-pt1[1])/nDp,(pt2[2]-pt1[2])/nDp};
      double coord[3]={pt1[0],pt1[1],pt1[2]};
      double solution[100]={0.0},xloc[100]={0.0};
      
      fprintf(filePtr,"x\ty\tz\t");
      for(s=0; s < nSol; ++s)
	{
	  fprintf(filePtr,"sol%d\t",s);
	}
      
      for(p=0; p < nPoints;++p)
	{
      apr_sol_xglob(fieldId,coord,solNr,el,xloc,solution,NULL,NULL,NULL,APC_CLOSE,
                    APE_SOL_XGLOB_DEFAULT | APE_SOL_XGLOB_MATCH_WITH_ADAPTATION);
	  fprintf(filePtr,"\n%lf\t%lf\t%lf\t",coord[0],coord[1],coord[2]);
	  // coords print
	  
	  for(s=0;s<nSol; ++s)
	    {
	      fprintf(filePtr,"%lf\t",solution[s]);	// val print
	    }
	  coord[0]+=dp[0];
	  coord[1]+=dp[1];
	  coord[2]+=dp[2];
	}
    }
  return 1;
}


/* Internal utilities */

/*---------------------------------------------------------
ut_mat3_inv - to invert a 3x3 matrix (stored as a vector!)
---------------------------------------------------------*/
double ut_mat3_inv(	/* returns: determinant of matrix to invert */
	double *mat,	/* matrix to invert */
	double *mat_inv	/* inverted matrix */
	)
{
double s0,s1,s2,rjac,rjac_inv;

  rjac = mat[0]*mat[4]*mat[8] + mat[3]*mat[7]*mat[2]
       + mat[6]*mat[1]*mat[5] - mat[6]*mat[4]*mat[2]
       - mat[0]*mat[7]*mat[5] - mat[3]*mat[1]*mat[8];
  rjac_inv = 1/rjac;

  s0 =   mat[4]*mat[8] - mat[7]*mat[5];
  s1 = - mat[3]*mat[8] + mat[6]*mat[5];
  s2 =   mat[3]*mat[7] - mat[6]*mat[4];

  mat_inv[0] = s0*rjac_inv;
  mat_inv[3] = s1*rjac_inv;
  mat_inv[6] = s2*rjac_inv;

  s0 =   mat[7]*mat[2] - mat[1]*mat[8];
  s1 =   mat[0]*mat[8] - mat[6]*mat[2];
  s2 = - mat[0]*mat[7] + mat[6]*mat[1];

  mat_inv[1] = s0*rjac_inv;
  mat_inv[4] = s1*rjac_inv;
  mat_inv[7] = s2*rjac_inv;

  s0 =   mat[1]*mat[5] - mat[4]*mat[2];
  s1 = - mat[0]*mat[5] + mat[3]*mat[2];
  s2 =   mat[0]*mat[4] - mat[3]*mat[1];

  mat_inv[2] = s0*rjac_inv;
  mat_inv[5] = s1*rjac_inv;
  mat_inv[8] = s2*rjac_inv;

  return(rjac);
}

/*---------------------------------------------------------
ut_vec3_prod - to compute vector product of 3D vectors
---------------------------------------------------------*/
void ut_vec3_prod(
	double* vec_a, 	/* in: vector a */
	double* vec_b, 	/* in: vector b */
	double* vec_c	/* out: vector product axb */
	)
{

vec_c[0]=vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1];
vec_c[1]=vec_a[2]*vec_b[0]-vec_a[0]*vec_b[2];
vec_c[2]=vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0];

return;
}

/*---------------------------------------------------------
ut_vec3_mxpr - to compute mixed vector product of 3D vectors
---------------------------------------------------------*/
double ut_vec3_mxpr( /* returns: mixed product [a,b,c] */
	double* vec_a, 	/* in: vector a */
	double* vec_b, 	/* in: vector b */
	double* vec_c	/* in: vector c */
	)
{
double daux;

daux  = vec_c[0]*(vec_a[1]*vec_b[2]-vec_a[2]*vec_b[1]);
daux += vec_c[1]*(vec_a[2]*vec_b[0]-vec_a[0]*vec_b[2]);
daux += vec_c[2]*(vec_a[0]*vec_b[1]-vec_a[1]*vec_b[0]);

return(daux);
}

/*---------------------------------------------------------
ut_vec3_length - to compute length of a 3D vector
---------------------------------------------------------*/
double ut_vec3_length(	/* returns: vector length */
	double* vec	/* in: vector */
	)
{
return(sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]));
}

/*---------------------------------------------------------
ut_mat3vec - to compute matrix vector product in 3D space
---------------------------------------------------------*/
void ut_mat3vec(
	double* m1, 	/* in: matrix (stored by rows as a vector!) */
	double* v1, 	/* in: vector */
	double* v2	/* out: resulting vector */
	)
{

v2[0] = m1[0]*v1[0] + m1[1]*v1[1] + m1[2]*v1[2] ;
v2[1] = m1[3]*v1[0] + m1[4]*v1[1] + m1[5]*v1[2] ;
v2[2] = m1[6]*v1[0] + m1[7]*v1[1] + m1[8]*v1[2] ;
return;
}

/*---------------------------------------------------------
ut_mat3mat - to compute matrix matrix product in 3D space
	(all matrices are stored by rows as vectors!)
---------------------------------------------------------*/
void ut_mat3mat(
	double* m1,	/* in: matrix */
	double* m2,	/* in: matrix */ 
	double* m3	/* out: matrix m1*m2 */
	)
{

m3[0] = m1[0]*m2[0] + m1[1]*m2[3] + m1[2]*m2[6] ;
m3[1] = m1[0]*m2[1] + m1[1]*m2[4] + m1[2]*m2[7] ;
m3[2] = m1[0]*m2[2] + m1[1]*m2[5] + m1[2]*m2[8] ;

m3[3] = m1[3]*m2[0] + m1[4]*m2[3] + m1[5]*m2[6] ;
m3[4] = m1[3]*m2[1] + m1[4]*m2[4] + m1[5]*m2[7] ;
m3[5] = m1[3]*m2[2] + m1[4]*m2[5] + m1[5]*m2[8] ;

m3[6] = m1[6]*m2[0] + m1[7]*m2[3] + m1[8]*m2[6] ;
m3[7] = m1[6]*m2[1] + m1[7]*m2[4] + m1[8]*m2[7] ;
m3[8] = m1[6]*m2[2] + m1[7]*m2[5] + m1[8]*m2[8] ;

}

/*---------------------------------------------------------
ut_chk_list - to check whether Num is on the list List
	with length Ll
---------------------------------------------------------*/
int ut_chk_list(	/* returns: */
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

