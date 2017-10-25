/* This is necessary for Approx(DG|STD|STDhybrid) modules
 * because of linking external adequate modules
 */

#ifndef _EXTERNAL_STUFF_H_
#define _EXTERNAL_STUFF_H_

/* Here are imported functions */

/* For mesh modules */
#include "mmh_intf.h"

/* For approximation module */
#include "aph_intf.h"

/* Some additonal functions */
extern "C" int apr_create_constr_data( /* returns: >=0 - success code, <0 - error code */
  int Field_id    /* in: approximation field ID  */
  );

/* This are unused functions but have to be implemented simply */

/* LINPACK */
extern "C" void    dgetrf_(int *m,int *n,double *a,int *lda,int *ipiv,int *info){}
extern "C" void    dgetrs_(char *trans,int *n,int *nrhs,double *a,int *lda,int *ipiv,double *b,int *ldb,int *info){}

/* Provblem dependend */
extern "C" double* pdr_select_el_coeff( /* returns: pointer !=NULL to indicate selection */
  double** Mval, /* out: mass matrix coefficient */
  double** Axx,  /* out: diffusion coefficients */		
  double** Axy,		
  double** Axz,		
  double** Ayx,
  double** Ayy,
  double** Ayz,
  double** Azx,
  double** Azy,
  double** Azz,
  double** Bx,	 /* out: convection coefficients */	
  double** By,	
  double** Bz,	
  double** Cval, /* out: reaction coefficients */
  double** Lval, /* out: rhs coefficient for time term */
  double** Qx,	 /* out: rhs coefficients with derivatives */
  double** Qy,	 /* out: rhs coefficients with derivatives */
  double** Qz,	 /* out: rhs coefficients with derivatives */
  double** Sval	 /* out: rhs coefficients without derivatives */	
  ){return(NULL);}

extern "C" int pdr_el_coeff(
	int Elem,	/* in: element number */
	int Mat_num,	/* in: material number */
	double Hsize,	/* in: size of an element */
	int Pdeg,	/* in: local degree of polynomial */
	double* Xcoor,	/* in: global coordinates of a point */
	double* U_val, 	/* in: computed solution */
	double* U_x, 	/* in: gradient of computed solution */
	double* U_y,   	/* in: gradient of computed solution */
	double* U_z,   	/* in: gradient of computed solution */
	double* Mval,	/* out: mass matrix coefficient */
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
	double* Cval,	/* out: reaction coefficients */
	double* Lval,	/* out: rhs coefficient for time term */
	double* Qx,	/* out: rhs coefficients with derivatives */
	double* Qy,	/* out: rhs coefficients with derivatives */
	double* Qz,	/* out: rhs coefficients with derivatives */
	double* Sval	/* out: rhs coefficients without derivatives */	
	) {return(0);}


#endif /* _EXTERNAL_STUFF_H_
*/