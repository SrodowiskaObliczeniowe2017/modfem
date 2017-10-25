#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"mod_fem.h"
#include"../include/aph_intf.h"

#define APC_MAXEQ 1
/* Some forward declarations */

__inline double ut_vec3_length(	/* returns: vector length */
	double* vec	/* in: vector */
	)
{
	return(sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]));
}

static 
double ut_mat3_inv(	/* returns: determinant of matrix to invert */
	double *mat,	/* matrix to invert */
	double *mat_inv	/* inverted matrix */
	);


int apr_shape_fun_1D( /* returns: the number of shape functions */
	int Pdeg, 	   /* in: degree of polynomial */
	double Eta,	   /* in: local coord of the considered point */
	double *Base_phi,  /* out: basis functions */
	double *Base_dphix/* out: x derivative of basis functions */
	);


int apr_shape_fun_2D( /* returns: the number of shape functions */
	int Base_type,	   /* in: type of basis functions: */
			   /* 	1 (APC_TENSOR) - tensor product */
			   /* 	2 (APC_COMPLETE) - complete polynomials */
	int Pdeg, 	   /* in: degree of polynomial - can be either */
			   /*	a single number, for isotropic p, */
			   /*	or a combination pdegy*10+pdegx */
	double *Eta,	   /* in: local coord of the considered point */
	double *Base_phi,  /* out: basis functions */
	double *Base_dphix,/* out: x derivative of basis functions */
	double *Base_dphiy /* out: y derivative of basis functions */
	);


int apr_shape_fun_3D( /* returns: the number of shape functions (<=0 - failure) */
	int Base_type,	   /* in: type of basis functions: */
			   /* 	1 (APC_TENSOR) - tensor product */
			   /* 	2 (APC_COMPLETE) - complete polynomials */
	int Pdeg, 	   /* in: degree of polynomial - can be either */
			   /*	a single number, for isotropic p, */
			   /*	or a combination pdegy*10+pdegx */
	double *Eta,	   /* in: local coord of the considered point */
	double *Base_phi,  /* out: basis functions */
	double *Base_dphix,/* out: x derivative of basis functions */
	double *Base_dphiy,/* out: y derivative of basis functions */
	double *Base_dphiz /* out: z derivative of basis functions */
	);


/*------------------------------------------------------------------
apr_elem_calc_3D - to perform element calculations (to provide data on
	coordinates, solution, shape functions, etc. for a given point
	inside element (given local coordinates Eta[i]);
	for geometrically multi-linear or linear 3D elements
-------------------------------------------------------------------*/
double apr_elem_calc_3D_dg(
	/* returns: Jacobian determinant at a point, either for */
	/* 	volume integration if Vec_norm==NULL,  */
	/* 	or for surface integration otherwise */
	int Control,	    /* in: control parameter (what to compute): */
			    /*	1  - shape functions and values */
			    /*	2  - derivatives and jacobian */
			    /* 	>2 - computations on the (Control-2)-th */
			    /*	     element's face */
	int Nreq,	    /* in: number of equations */
	int Pdeg,	    /* in: element degree of polynomial */
	int Base_type,	    /* in: type of basis functions: */
			    /* 	1 (APC_TENSOR) - tensor product */
			    /* 	2 (APC_COMPLETE) - complete polynomials */
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
#ifdef FV_DEBUG
  if(Nreq>APC_MAXEQ){
    printf("Number of equations %d greater than the limit %d.\n",Nreq,APC_MAXEQ);
    printf("Change APC_MAXEQ in aph_dg_prism.h and recompile the code.\n");
    exit(-1);
  }
#endif

/* set the number of geometry dofs */
  nrgeo = 6;

/* get the values of shape functions and their LOCAL derivatives */
  if(Base_phi!=NULL){
    if(Control==1)
      nrdofs=apr_shape_fun_3D(Base_type, Pdeg, Eta,
	Base_phi, NULL, NULL, NULL);
    else
      nrdofs=apr_shape_fun_3D(Base_type, Pdeg, Eta,
	Base_phi, Base_dphix, Base_dphiy, Base_dphiz);
  }

  geo_phi[0]=(1.0-Eta[0]-Eta[1])*(1.0-Eta[2])/2.0;
  geo_phi[1]=Eta[0]*(1.0-Eta[2])/2.0;
  geo_phi[2]=Eta[1]*(1.0-Eta[2])/2.0;
  geo_phi[3]=(1.0-Eta[0]-Eta[1])*(1.0+Eta[2])/2.0;
  geo_phi[4]=Eta[0]*(1.0+Eta[2])/2.0;
  geo_phi[5]=Eta[1]*(1.0+Eta[2])/2.0;

/*physical coordinates*/
  if(Xcoor!=NULL) {
#ifdef FV_DEBUG
	  if(Node_coor==NULL) {
		  printf("Error: null coordinate table. Exiting!\n");
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
  if(Sol!=NULL) {
#ifdef FV_DEBUG
	  if(Sol_dofs==NULL) {
		  printf("Error: null dofs table. Exiting!\n");
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
  if(Base_dphix!=NULL) {
#ifdef FV_DEBUG
	  if(Base_dphiy==NULL&&Base_dphiz==NULL) {
		  printf("Error: null base phi table. Exiting!\n");
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
  if(Dsolx!=NULL) {
#ifdef FV_DEBUG
	  if(Dsoly==NULL&&Dsolz==NULL) {
		  printf("Error: null gradient table. Exiting!\n");
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


int apr_shape_fun_3D( /* returns: the number of shape functions (<=0 - failure) */
	int Base_type,	   /* in: type of basis functions: */
			   /* 	1 (APC_TENSOR) - tensor product */
			   /* 	2 (APC_COMPLETE) - complete polynomials */
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

  //int num_shap_old;
//  double Base_phi_old[APC_MAXELVD];
//  double Base_dphix_old[APC_MAXELVD];
 // double Base_dphiy_old[APC_MAXELVD];
 // double Base_dphiz_old[APC_MAXELVD];

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef FV_DEBUG
  if(Base_phi==NULL) {
    printf("No space provided for shape functions in apr_shape_fun_3D!\n");
    exit(1);
  }
#endif

  if(Pdeg==0){
    /* for piecewise constant approximation */

    Base_phi[0] = 1.0;

	if(Base_dphix!=NULL) {

#ifdef FV_DEBUG
	  if(Base_dphiy==NULL&&Base_dphiz==NULL) {
		  printf("Error: null base phi table. Exiting!\n");
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

	if(Base_dphix!=NULL) {

#ifdef FV_DEBUG
	  if(Base_dphiy==NULL&&Base_dphiz==NULL) {
		  printf("Error: null base phi table. Exiting!\n");
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
  else if(Pdeg==1&&Base_type == APC_TENSOR_DG){
    /* for piecewise constant in z and linear in xy approximation */

    Base_phi[0] = 1.0;
    Base_phi[1] = Eta[0];
    Base_phi[2] = Eta[1];

	if(Base_dphix!=NULL) {

#ifdef FV_DEBUG
	  if(Base_dphiy==NULL&&Base_dphiz==NULL) {
		  printf("Error: null base phi table. Exiting!\n");
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
  else if(Pdeg==2&&Base_type == APC_TENSOR_DG){
/* for piecewise constant in z and quadratic in xy approximation */

    Base_phi[0] = 1.0;
    Base_phi[1] = Eta[0];
    Base_phi[2] = Eta[1];
    Base_phi[3] = Eta[0]*Eta[0];
    Base_phi[4] = Eta[0]*Eta[1];
    Base_phi[5] = Eta[1]*Eta[1];

	if(Base_dphix!=NULL) {
		
#ifdef FV_DEBUG
	  if(Base_dphiy==NULL&&Base_dphiz==NULL) {
		  printf("Error: null base phi table. Exiting!\n");
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

	if(Base_dphix!=NULL) {
		
#ifdef FV_DEBUG
	  if(Base_dphiy==NULL&&Base_dphiz==NULL) {
		  printf("Error: null base phi table. Exiting!\n");
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

    if(Base_type==APC_TENSOR_DG&&(Pdeg<100||Pdeg%100==0)){
      printf("Sorry %d, at least linears required both for z and xy polynomials!\n", Pdeg);
      exit(1);
    }


    Base_phi[4] = Eta[2]*Eta[0];
    Base_phi[5] = Eta[2]*Eta[1];

	if(Base_dphix!=NULL) {
		
#ifdef FV_DEBUG
	  if(Base_dphiy==NULL&&Base_dphiz==NULL) {
		  printf("Error: null base phi table. Exiting!\n");
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
	if(Base_dphix!=NULL) {
		
#ifdef FV_DEBUG
	  if(Base_dphiy==NULL&&Base_dphiz==NULL) {
		  printf("Error: null base phi table. Exiting!\n");
		  exit(-1);
	  }
#endif

      Base_dphix[6] = Eta[1];
      Base_dphiy[6] = Eta[0];
      Base_dphiz[6] = 0.0;
    }

    if(Pdeg==2&&Base_type==APC_COMPLETE_DG){
      /* to speed up computations with 2nd order complete polynomial shape functions*/

      Base_phi[7] = Eta[0]*Eta[0];
      Base_phi[8] = Eta[1]*Eta[1];
      Base_phi[9] = Eta[2]*Eta[2];

	  if(Base_dphix!=NULL) {
		  
#ifdef FV_DEBUG
	  if(Base_dphiy==NULL&&Base_dphiz==NULL) {
		  printf("Error: null base phi table. Exiting!\n");
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
	if(Base_dphix!=NULL) {
#ifdef FV_DEBUG
	  if(Base_dphiy==NULL&&Base_dphiz==NULL) {
		  printf("Error: null base phi table. Exiting!\n");
		  exit(-1);
	  }
#endif
      Base_dphix[7] = Eta[2]*Eta[1];
      Base_dphiy[7] = Eta[2]*Eta[0];
      Base_dphiz[7] = Eta[1]*Eta[0];
    }

    base_type_2D=APC_COMPLETE_DG;

    if(Base_type==APC_TENSOR_DG){

      pdegx = Pdeg%100;
      pdegz = Pdeg/100;

    }
    else if(Base_type==APC_COMPLETE_DG){

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

      printf("Type of base (BASE in control.dat) not valid for prisms!\n");
      return(-1);

    }

	if(Base_dphix!=NULL) {
		
#ifdef FV_DEBUG
	  if(Base_dphiy==NULL&&Base_dphiz==NULL) {
		  printf("Error: null base phi table. Exiting!\n");
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

#ifdef FV_DEBUG
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


    if (Base_type == APC_TENSOR_DG) {

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
    else if(Base_type == APC_COMPLETE_DG ) {

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

	if(Base_dphix!=NULL) {
		
#ifdef FV_DEBUG
	  if(Base_dphiy==NULL&&Base_dphiz==NULL) {
		  printf("Error: null base phi table. Exiting!\n");
		  exit(-1);
	  }
#endif

      if (Base_type == APC_TENSOR_DG) {

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
      else if(Base_type == APC_COMPLETE_DG ) {

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
			   /* 	1 (APC_TENSOR) - tensor product */
			   /* 	2 (APC_COMPLETE) - complete polynomials */
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

//  int num_shap_old;
  //double Base_phi_old[APC_MAXELVD];
  //double Base_dphix_old[APC_MAXELVD];
  //double Base_dphiy_old[APC_MAXELVD];

/*++++++++++++++++ executable statements ++++++++++++++++*/

  if(Base_phi==NULL) {
    printf("No space provided for shape functions in apr_shape_fun_2D!\n");
    exit(1);
  }

  if(Pdeg==0){
/* for piecewise constant approximation */

    Base_phi[0]=1.0;

    if(Base_dphix!=NULL && Base_dphiy!=NULL){
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

    if(Base_dphix!=NULL && Base_dphiy!=NULL){

      Base_dphix[0]=0.0;
      Base_dphix[1]=0.0;

      Base_dphiy[0]=0.0;
      Base_dphiy[1]=1.0;

    }

    num_shap=2;

    return(num_shap);

  }
  else if(Pdeg==1 && Base_type==APC_TENSOR_DG){
/* for linear in x and piecewise constant in y approximation */

    Base_phi[0]=1.0;
    Base_phi[1]=Eta[0];

    if(Base_dphix!=NULL && Base_dphiy!=NULL){

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

    if(Base_dphix!=NULL && Base_dphiy!=NULL){

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
    if(Base_dphix!=NULL && Base_dphiy!=NULL){
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

      if(Base_type == APC_COMPLETE_DG ) {

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
      else if (Base_type == APC_TENSOR_DG) {

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

	if(Base_type == APC_COMPLETE_DG ) {

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
	else if (Base_type == APC_TENSOR_DG) {

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
#ifdef DEBUG
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

#ifdef FV_DEBUG
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

#undef APC_MAXEQ


