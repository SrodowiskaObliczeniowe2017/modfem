#ifndef _MESH_PHP_H_
#define _MESH_PHP_H_

#if defined(__cplusplus)
extern "C" {
#endif

#ifndef FV_DEBUG
	#if defined(DEBUG) || defined(_DEBUG)
		#define FV_DEBUG
	#endif
#endif


int apr_get_pdeg_nrdofs(int base, int pdeg);  

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
			    /*  (APC_TENSOR) - tensor product */
			    /*  (APC_COMPLETE) - complete polynomials */
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

double apr_elem_calc_3D_std_lin(
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
			    /*  (APC_TENSOR) - tensor product */
			    /*  (APC_COMPLETE) - complete polynomials */
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

#if defined(__cplusplus)
}
#endif

#endif /*_MESH_PHP_H_
*/
