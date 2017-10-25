#ifndef _ALLCALC3D_H_
#define _ALLCALC3D_H_

#include<cstdlib>
#include<cstdio>
#include<cmath>

#include "aph_intf.h"



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

double ut_mat3_inv(
	/* returns: determinant of matrix to invert */
	double *mat,	/* matrix to invert */
	double *mat_inv	/* inverted matrix */
	);

double ut_vec3_length(	/* returns: vector length */
	double* vec	/* in: vector */
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

int apr_shape_fun_1D( /* returns: the number of shape functions */
	int Pdeg, 	   /* in: degree of polynomial */
	double Eta,	   /* in: local coord of the considered point */
	double *Base_phi,  /* out: basis functions */
	double *Base_dphix/* out: x derivative of basis functions */
	);

/**--------------------------------------------------------
ut_mat3_inv - to invert a 3x3 matrix (stored as a vector!)
---------------------------------------------------------*/
double ut_mat3_inv(	/* returns: determinant of matrix to invert */
	double *mat,	/* matrix to invert */
	double *mat_inv	/* inverted matrix */
	);

double ut_vec3_length(	/* returns: vector length */
	double* vec	/* in: vector */
	);

#endif /* _ALLCAL3D_H_ 
*/
