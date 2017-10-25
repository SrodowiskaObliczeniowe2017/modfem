#include<cassert>
#include"fv_mathparser.h"
#include"fv_str_utils.h"
#include"fv_conv_utils.h"
#include"fv_exception.h"

#include "aph_intf.h"

#include<iostream>
#include<sstream>
#include<string>
#include<cmath>
#include<string>

namespace fvmathParser {
#define MY_PRISM 0
#if MY_PRISM
inline void geom_func(double geo_phi[6],double Eta[3]);
inline void deriv_geom_func(double geo_dphix[6],double geo_dphiy[6],double geo_dphiz[6],double Eta[3]);
#endif

#define ut_max(x,y) ((x)>(y)?(x):(y))
#define ut_min(x,y) ((x)<(y)?(x):(y))
#define ut_abs(x)   ((x)<0?-(x):(x))


double ut_vec3_length(	/* returns: vector length */
	double* vec	/* in: vector */
	)
{
   return(sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]));
}

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
	);

double my_apr_elem_calc_3D_mod(
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
	#if MY_PRISM
	  geom_func(geo_phi,Eta);
	#else
	  geo_phi[0]=(1.0-Eta[0]-Eta[1])*(1.0-Eta[2])/2.0;
	  geo_phi[1]=Eta[0]*(1.0-Eta[2])/2.0;
	  geo_phi[2]=Eta[1]*(1.0-Eta[2])/2.0;
	  geo_phi[3]=(1.0-Eta[0]-Eta[1])*(1.0+Eta[2])/2.0;
	  geo_phi[4]=Eta[0]*(1.0+Eta[2])/2.0;
	  geo_phi[5]=Eta[1]*(1.0+Eta[2])/2.0;
	#endif
	/*physical coordinates*/
	  if(Xcoor!=NULL){


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
	#if MY_PRISM
	  deriv_geom_func(geo_dphix,geo_dphiy,geo_dphiz,Eta);
	#else
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
	#endif
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

#if MY_PRISM
/** Reference prizm reorganization */
inline
void geom_func(double geo_phi[6],double Eta[3])
{
	geo_phi[0]=-(Eta[0] + Eta[2])*(1.0 - Eta[1])*0.25;
	geo_phi[1]= (Eta[0] + 1.0)*(1.0 - Eta[1])*0.25;
	geo_phi[2]= (Eta[0] + 1.0)*(1.0 + Eta[1])*0.25;
	geo_phi[3]=-(Eta[0] + Eta[2])*(1.0 + Eta[1])*0.25;
	geo_phi[4]= (1.0 - Eta[1])*(1.0 + Eta[2])*0.25;
	geo_phi[5]= (1.0 + Eta[1])*(1.0 + Eta[2])*0.25;
}

inline
void deriv_geom_func(double geo_dphix[6],double geo_dphiy[6],double geo_dphiz[6],double Eta[3])
{
	const double a = (1.0-Eta[1])*0.25;
	const double b = (1.0+Eta[1])*0.25;
	const double c = (Eta[0]+Eta[2])*0.25;
	const double d = (1.0+Eta[0])*0.25;
	const double e = (1.0+Eta[2])*0.25;
	geo_dphix[0] = -a;
	geo_dphix[1] = a;
	geo_dphix[2] = b;
	geo_dphix[3] = -b;
	geo_dphix[4] = 0.0;
	geo_dphix[5] = 0.0;

	geo_dphiy[0] = c;
	geo_dphiy[1] = -d;
	geo_dphiy[2] = d;
	geo_dphiy[3] = -c;
	geo_dphiy[4] = -e;
	geo_dphiy[5] = e;

	geo_dphiz[0] = -a;
	geo_dphiz[1] = 0.0;
	geo_dphiz[2] = 0.0;
	geo_dphiz[3] = -b;
	geo_dphiz[4] = a;
	geo_dphiz[5] = b;
}
#endif

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



//namespace fvmathParser {

	int MathCalculator::GetElement(
		ElemType &type, 
		double &value, 
		std::string &str)
	{
		// if only one character

		// Check sign/operation
		if (str.length() == 1)
		{
			if (str[0] == '^')
			{
				type = Power;
				return 1;
			}
			else if (str[0] == '*')
			{
				type = Multiply;
				return 1;
			}
			else if (str[0] == '/')
			{
				type = Divide;
				return 1;
			}
			else if (str[0] == '+')
			{
				type = Plus;
				return 1;
			}
			else if (str[0] == '-')
			{
				type = Minus;
				return 1;
			}
			else if (str[0] == '(')
			{
				type = Open;
				return 1;
			}
			else if (str[0] == ')')
			{
				type = Close;
				return 1;
			}
		}//if (str.length() == 1)

		// Check if a vector
		if ( (str[0] == 'v') || (str[0] == 'V'))
		{
			try
			{
				std::string str_nr(str.substr(1));
				// Get suffix ex. v0 -> 0
				value = Str2T<double>( str_nr );
				type = Vector;
				return 1;
			}
			catch(...)
			{
				return(-1);
			}
		}

		// Check if a number
		try
		{
			value = (double)Str2T<double>(str);
			type = Number;
			return(1);
		}
		catch(...)
		{
			return(-1);
		}

		// Unrconizable element
		return(-1);
		
	}

	void MathCalculator::GetElement(MathElement& el, std::string &str)
	{
		if (str.length() == 1)
		{
			switch(str[0])
			{
				case '^':
						el.type = Power;
						return;
					break;

				case '*':
						el.type = Multiply;
						return;
					break;

				case '/':
						el.type = Divide;
						return;
					break;
				case '+':
						el.type = Plus;
						return;
					break;
				case '-':
						el.type = Minus;
						return;
					break;
				case '(':
						el.type = Open;
						return;
					break;
				case ')':
						el.type = Close;
						return;
					break;
				//default:
				//    throw CException(-1, "znak nieobslugiwany");
			}//switch(str[0])
		}//if (str.length() == 1)

		if ( (str[0] == 'v') || (str[0] == 'V'))
		{
			try
			{
				//sprawdzenie czy skalda sie tylko z cyfr
				if (str.substr(1).find_first_not_of("0123456789") != std::string::npos)
					throw fv_exception("Error: unable to specify vector's component");
				std::string str_nr(str.substr(1));
				el.value = Str2T<double>(str_nr);
				el.type = Vector;
				return;
				//return 1;
			}
			catch(fv_exception& ex)
			{
				throw ex;
			}
			catch(...)
			{
				//return -1;
				//wyjatek: blad odczytu numeru wektora
				throw fv_exception("Error: unable to specify vector's component");
			}
		}//if ( (str[0] == 'v') || (str[0] == 'V'))

		//14.50

		//jesli pierwszy znakjest + lub -
		//to wykonac dlasza analize
		if ( (str[0] == '+') || (str[0] == '-'))
		{
			//sprawdzenie czy po za pierwszym znakiem dalsze to cyfry pozostala czesc to cyfry
			if (str.substr(1).find_first_not_of("0123456789.") != std::string::npos)
				throw fv_exception("Blad danych elementu - nie mozliwa interpretacja.");
			//sprawdznie czy jest tylko jedna kropka
			try
			{
				el.value = Str2T<double>(str);
				el.type = Number;
				return;
			}
			catch(fv_exception& ex)
			{
				throw ex;
			}
			catch(...)
			{
				throw fv_exception("Error while analizing data - str2double.");
			}

			return;
		}

		if (str.find_first_not_of("0123456789.") != std::string::npos)
		//if (str.find_first_not_of("0123456789.-") != string::npos)
			throw fv_exception("Blad danych elementu - nie mozliwa interpretacja(nie cyfra i nie znak przecinkowy).");

		//sprawdznie czy jest tylko jedna kropka

		try
		{
			el.value = Str2T<double>(str);
			el.type = Number;
			return;
		}
		catch(...)
		{
			throw fv_exception( "Blad analizy danych elementu(strToDouble).");
		}

		//nieznany blad analizy elementu
		throw fv_exception("Nieznany blad analizy elementu.");

	}



	int MathCalculator::ONP(std::string func, unsigned int vectorSize, std::vector<MathElement> &wyjscie)
	{
		std::stack<MathElement> stos;
		//funkcja po dodaniu ewentualnej spacji
		std::string strFunc;

		//stos.clear();
		wyjscie.clear();

		//jesli string jest pusty
		if( func.empty())
			return -1;

		//jesli ostatni element nie jest spacja to ja dodajemy na koncu

		if( func[func.length()-1] != ' ')
		{
			strFunc = func + " ";
		}
		else
			strFunc = func;

		int cntSpaces = countSubString(strFunc, " ");


		int start, stop;
		std::string str;

		//je�li liczba spacji jest mniejsza lub r�wna  zero
		//to sprawdzic czy to jest wektor lub liczba (konewrsja do double) je�li co innego to b�ad

		if ( (cntSpaces < 0) || (cntSpaces == 2) )
			return -1;

		MathElement elTmp;

		// Case - only one expression
		if(cntSpaces == 1)
		{
			MathCalculator::GetElement(elTmp, func);

			// If given element is a number or vector
			if((elTmp.type == Number) || (elTmp.type == Vector))
			{
				wyjscie.push_back(elTmp);
				return(1);
			}
			else
			{
				return(-1);
			}
		}

		// Case - number of operation >= 3

		// Loop over all elements seperated by spaces
		for(int i = 0; i < cntSpaces; ++i)
		{

			start = getSubStrPosition(func, " ", i);
			stop = getSubStrPosition(func, " ", i+1);
			str = func.substr(start+1, stop-start-1);

			//*************************************************************************
			//*************************************************************************
			//  odczytanie elementow
			//*************************************************************************

			//2
			try
			{
				MathCalculator::GetElement(elTmp, str);
				//jesli element jest wektorem to sprawdzic czy skladowa miesci sie w zasiegu
				if (elTmp.type == Vector)
				{
					if ( (elTmp.value < 0) || (  elTmp.value >  vectorSize-1  ) )
						throw fv_exception("Skladowa wektora z po za zakresu.");
				}
			}
			catch(fv_exception& ex)
			{
				//int i = 42;
				std::ostringstream ss;
				ss << "Blad analizy elementu ";
				ss << i;
				ss << " : " << ex;
				throw fv_exception(ss.str() );
			}

			//*************************************************************************
			//*************************************************************************
			//  analiza elementow
			//*************************************************************************

			//test polega na analizie i sprawdzaniu czy na zmiane za podawane
			//wartosci z znakami i koncza sie wartosciami

			//jesli wartosc to wypisa� na wyjscie
			if( (elTmp.type == Number) || (elTmp.type == Vector) )
			{

				//zapisanie na wyjscie
				wyjscie.push_back(elTmp);

			}

			//jesli jeden ze znakow
			//PLUS, MIN, Multiply, Divide, POW
			//to
			if (
				(elTmp.type == Plus) ||
				(elTmp.type == Minus) ||
				(elTmp.type == Multiply) ||
				(elTmp.type == Divide) ||
				(elTmp.type == Power)
				)
			{

				//jesli wektor jest pusty to dodac
				if (stos.size() == 0)
				{
					stos.push(elTmp);
				}
				else//jesli stos nie jest pusty
				{

					//jesli nowy operator jest wyrzszy od aktualnego
					//to zapisujemy element na stos
					if ( MathCalculator::GetPriority(elTmp.type) > MathCalculator::GetPriority(stos.top().type) )
					{
						stos.push(elTmp);
					}
					else
					{
						//jesli nowy operator jest mniejszy lub rwony to
						//przenosimy na wyjscie wszytskie o priorytecie r�wnym lub wyrzszym
						//az do uzyskania priorytetu nizszego
						//badany element zapisac na stos

						//petla z warunkiem
						////while (mathCalculator::GetPriority(elTmp.type) >= mathCalculator::GetPriority(stos.top().type))
						//while (mathCalculator::GetPriority(elTmp.type) < mathCalculator::GetPriority(stos.top().type))
						while (MathCalculator::GetPriority(elTmp.type) <= MathCalculator::GetPriority(stos.top().type))
						{
							wyjscie.push_back(stos.top());
							stos.pop();
							if (stos.size() == 0)
								break;
						}
						stos.push(elTmp);
					}



				}

				//jesli nawias otwierajacy
				//to zapisac go na stos

				//jesli nawias zamykajacy to zdejmowec ze stosu
				//i zapisywac na wyjscie az do napotkania na otwierajacy
			}


			//jesli nawias otwieraj�cy to
			if ( elTmp.type == Open )
			{
				stos.push(elTmp);
			}

			//jesli nawias zamykajacy
			if ( elTmp.type == Close )
			{

				while (stos.top().type != Open)
				{
					wyjscie.push_back(stos.top());
					stos.pop();
				}
				stos.pop();

			}

			#ifdef __WXDEBUG__

			std::cout << std::endl;

			std::cout << "wyjscie: " << mathCalculator::GetString(wyjscie, stos) << std::endl;

			std::cout << std::endl;

			#endif



			//sprawdzenie elementu czym jest
			//i stworzeniem odpowiedniego elementu
			//
			//znakiem: ^, *, /, +, -
			//nawiazem: (, )
			//liczba
			//

			#ifdef __WXDEBUG__
			std::cout << std::endl;
			#endif

		}//for (int i = 0; i < liczbaSpacji; ++i)

		//po zakonczeniu analizy
		//przepisanie wszystkich elementow z stosu na wyjscie
		while (!stos.empty())
		{
			wyjscie.push_back(stos.top());
			stos.pop();
		}

		//#ifdef __WXDEBUG__
		//    std::cout << "Ostateczne rozwiazanie:" << std::endl;
		//    std::cout << "wyjscie: " << mathCalculator::GetString(wyjscie, stos) << std::endl;
		//    std::cout << std::endl;
		//#endif

		if ( (wyjscie.size() % 2) != 1 )
			throw fv_exception( "Niepoprawna liczba elementow funkcji." );

		return 1;

	}

	std::string MathCalculator::GetString(std::vector<MathElement> &vElement, std::stack<MathElement> &stos)
	{

		std::stringstream oss;

		for(size_t i = 0; i < vElement.size(); ++i)
		{

		//if(vElement[i].type

		//POW
		//Multiply, Divide
		//PLUS, Minus
		//OPEN, CLOSE
		//NUM, VEC

			switch(vElement[i].type)
			{
				case Power:
					oss << "^ ";
					break;

				case Multiply:
					oss << "* ";
					break;

				case Divide:
					oss << "/ ";
					break;

				case Plus:
					oss << "+ ";
					break;

				case Minus:
					oss << "- ";
					break;

				case Number:
					oss << vElement[i].value << " ";
					break;

				case Vector:
					oss << "v" << vElement[i].value << " ";
					break;
			}

		}//for (unsigned int i = 0; i < vElement.size(); ++i)

		oss << "        stos top(" << stos.size() << "): ";

		if (stos.size() == 0)
		{
			oss << "non";
		}
		else
		//for (unsigned int i = 0; i < stos.size(); ++i)
		{

			switch(stos.top().type)
			{
				case Power:
					oss << "^ ";
					break;

				case Multiply:
					oss << "* ";
					break;

				case Divide:
					oss << "/ ";
					break;

				case Plus:
					oss << "+ ";
					break;

				case Minus:
					oss << "- ";
					break;

				case Open:
					oss << "( ";
					break;

				case Close:
					oss << ") ";
					break;
			}
		}//if (stos.size() == 0)

		return oss.str();
	}

	std::string MathCalculator::GetString(std::vector<MathElement> &vElement)
	{

		std::stringstream oss;
		for (size_t i = 0; i < vElement.size(); ++i)
		{
			switch(vElement[i].type)
			{
				case Power:
					oss << "^ ";
				break;

				case Multiply:
					oss << "* ";
				break;

				case Divide:
					oss << "/ ";
				break;

				case Plus:
					oss << "+ ";
				break;

				case Minus:
					oss << "- ";
				break;

				case Number:
					oss << vElement[i].value << " ";
				break;

				case Vector:
					oss << "v" << vElement[i].value << " ";
				break;
			}

		}//for (unsigned int i = 0; i < vElement.size(); ++i)

		return oss.str();

	}

	int MathCalculator::GetPriority(ElemType ElType)
	{

		//NUM, VEC, PLUS, MIN, Multiply, Divide, POW, OPEN, CLOSE

		switch(ElType)
		{
			case Power :
				return 3;
				break;

			case Multiply :
				return 2;
				break;
			case Divide :
				return 2;
				break;

			case Plus :
				return 1;
				break;
			case Minus :
				return 1;
				break;

			default:
				return -1;
				break;
		}

	}

	void MathCalculator::Calculate(MathElement &a, MathElement &b, MathElement &sign, MathElement &result, std::vector<double> v)
	{

		double aa, bb;
		int i;

		//jesli wektor lub wartosc numeryczna

		if ( (a.type == Vector) || (a.type == Number) )
		{
		//czy pozycja wektora jest
		if (a.type == Number)
		{
			aa = a.value;
		}
		//jesli jest wektorem to sprawdzic czy miesci sie w zakresie
		else if (a.type == Vector)
		{
			i = (int)a.value;
			if ( (i < 0) || ( i >  (int)v.size() ) )
			{
				//blad zakresu wektora
				//return -2;
				throw fv_exception(  "Wspolzedna wektora z po za zakresu.");
			}
			else
				aa = v[i];
		}
		}
		else
		//return -4;
		throw fv_exception( "Niepoprawna wspolrzedna dzialania");

		if ( (b.type == Vector) || (b.type == Number) )
		{
		//czy pozycja wektora jest
		if (b.type == Number)
		{
			bb = b.value;
		}
		else if (b.type == Vector)
		{
			i = (int)b.value;
			if ( (i < 0) || ( i > (int)v.size() ) )
			{
				//blad zakresu wektora
				//return -2;
				throw fv_exception( "Wspolzedna wektora z po za zakresu.");
			}
			else
				bb = v[i];
		}
		}
		else
		throw fv_exception( "Niepoprawna wspolrzedna dzialania");
		//return -5;

		//czy znak jest znakiem

		//POW
		//Multiply, Divide
		//PLUS, MIN

		double r;

		if (
		(sign.type == Power) ||
		(sign.type == Multiply) ||
		(sign.type == Divide) ||
		(sign.type == Plus) ||
		(sign.type == Minus)
		)
		{
		//wykonanie obliczenia

		if (sign.type == Power)
		{
			r = pow(aa, bb);

			#ifdef __WXDEBUG__
				std::cout << aa << " ^ " << bb << " = " << r << std::endl;
			#endif
		}
		else if (sign.type == Multiply)
		{
			r = aa * bb;

			#ifdef __WXDEBUG__
				std::cout << aa << " * " << bb << " = " << r << std::endl;
			#endif
		}
		else if (sign.type == Divide)
		{
			r = aa / bb;
			#ifdef __WXDEBUG__
				std::cout << aa << " / " << bb << " = " << r << std::endl;
			#endif
		}
		else if (sign.type == Plus)
		{
			r = aa + bb;
			#ifdef __WXDEBUG__
				std::cout << aa << " + " << bb << " = " << r << std::endl;
			#endif
		}
		else if (sign.type == Minus)
		{
			r = aa - bb;
			#ifdef __WXDEBUG__
				std::cout << aa << " - " << bb << " = " << r << std::endl;
			#endif
		}
		}
		else
		//return -3;
		throw fv_exception( "Niepoprawna operacja dzialania");

		result.type = Number;
		result.value = r;

		return;

	}




	int MathCalculator::ONPCalculate(std::vector<MathElement> &vElement, std::vector<double> v, double &result)

	{

		std::stack<MathElement> stos;

		MathElement a, b, r;

		if(vElement.size() == 0 || vElement.size() == 2)
		{
			throw fv_exception( "Blad wykonania elmentow wejsciowych - 0");
		}

		if(vElement.size() == 1)
		{
			//jesli vec lib liczba to ok
			//w przeciwnym wypaku blad
			if(vElement[0].type == Number)
			{
				result = vElement[0].value;
				return 1;
			}
			else if (vElement[0].type == Vector)
			{
				assert(vElement[0].value >= 0);
				result = v[static_cast<int>(vElement[0].value)];
				return 1;
			}
			else
				throw fv_exception( "Blad wykonania elmentow wejsciowych 2");
		}

		//petla po wszytskich elementach
		for(size_t i = 0; i < vElement.size(); ++i)
		{

			//jesli stala lub vec to zapisac na stos
			if((vElement[i].type == Number) || (vElement[i].type == Vector) )
			{
				#ifdef __WXDEBUG__
					std::cout << "Dodanie elementu na stos" << std::endl;
				#endif

				stos.push(vElement[i]);
			}

			if(
				(vElement[i].type == Power) ||
				(vElement[i].type == Multiply) ||
				(vElement[i].type == Divide) ||
				(vElement[i].type == Plus) ||
				(vElement[i].type == Minus)
				)
			{

				#ifdef __WXDEBUG__
					std::cout << "Wykonanie dzialania: " << std::endl;
				#endif

				b = stos.top();
				stos.pop();

				a = stos.top();
				stos.pop();

				#ifdef __WXDEBUG__
					std::cout << "Wykonanie dzialania 2: " << std::endl;
				#endif

				try
				{
					MathCalculator::Calculate(a, b, vElement[i], r, v);
				}
				catch(fv_exception& ex)
				{
					std::ostringstream os;
					os << ex;
					throw fv_exception("Blad wykonania dialana: " + os.str());
				}

				stos.push(r);

			}

		}//for (unsigned int i = 0; i < vElement.size(); ++i)

		result = stos.top().value;

		return 1;
	}
	

}// end namespace fvmatchParser
