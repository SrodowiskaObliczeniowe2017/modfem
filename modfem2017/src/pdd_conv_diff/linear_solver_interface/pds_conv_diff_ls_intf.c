/************************************************************************
File pds_conv_diff_ls_intf.c - interface between the problem dependent 
         module for convection-diffusion time dependent equations 
         and linear solver modules (direct and iterative)


Contains definitions of routines:   

  pdr_get_list_ent - to return to the solver module :
                          1. the list of integration entities - entities
                             for which stiffness matrices and load vectors are
                             provided by the FEM code
                          2. the list of DOF entities - entities with which  
                             there are dofs associated by the given approximation
  pdr_get_list_ent_coarse - the same as above but for COARSE level and
                            given the corresponding lists from the fine level 
  pdr_get_max_num_grid_levels - for limiting nr_levels in multigrid 
                              based on mesh and field data
  pdr_create_assemble_stiff_mat - to create element stiffness matrices
                                   and assemble them to the global SM
  pdr_assemble_local_stiff_mat_with_table - to assemble an element stiffness matrix
                                   to the global SM using assembly table
  pdr_assemble_local_stiff_mat - to assemble an element stiffness matrix
                                   to the global SM
  pdr_comp_stiff_mat - to provide a solver with a stiffness matrix 
                      and a load vector corresponding to the specified
                      mesh entity, together with information on how to
                      assemble entries into the global stiffness matrix
                      and the global load vector
  OBSOLETE pdr_select_el_coeff - to select coefficients for internal integrals
  pdr_select_el_coeff_vect - to select coefficients for internal integrals
  pdr_el_coeff - to return coefficients for internal integrals


  pdr_read_sol_dofs - to read a vector of dofs associated with a given 
                   mesh entity from approximation field data structure
  pdr_write_sol_dofs - to write a vector of dofs associated with a given 
                   mesh entity to approximation field data structure

  pdr_L2_proj_sol - to project solution between elements of different generations
  pdr_renum_coeff - to return a coefficient being a basis for renumbering

  pdr_get_ent_pdeg - to return the degree of approximation index 
                      associated with a given mesh entity
  pdr_dof_ent_sons - to return info whether the entity is owned
                     and a list of dof entity sons for owned entity

  pdr_proj_sol_lev - to project solution between mesh levels
  pdr_vec_norm - to compute a norm of global vector in parallel
  pdr_sc_prod - to compute a scalar product of two global vectors
  pdr_create_exchange_tables - to create tables to exchange dofs
  pdr_exchange_dofs - to exchange dofs between processors

------------------------------  			
History:        
	02.2002 - Krzysztof Banas, initial version		
	12.2012 - Krzysztof Banas, rewritten with utr_... calls		
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

#ifdef PARALLEL
/* interface of parallel mesh manipulation modules */
#include "mmph_intf.h"
/* interface for all parallel approximation modules */
#include "apph_intf.h"
/* interface for parallel communication modules */
#include "pch_intf.h"
#endif

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
#include "../include/pdh_conv_diff_weakform.h"	



/*** CONSTANTS ***/

const int PDC_ELEMENT = APC_ELEMENT;
const int PDC_FACE    = APC_FACE   ;
const int PDC_EDGE    = APC_EDGE   ;
const int PDC_VERTEX  = APC_VERTEX ;

const int PDC_NO_COMP  = APC_NO_COMP  ;/* do not compute stiff mat and rhs vect*/
const int PDC_COMP_SM  = APC_COMP_SM  ;/* compute entries to stiff matrix only */
const int PDC_COMP_RHS = APC_COMP_RHS ;/* compute entries to rhs vector only */
const int PDC_COMP_BOTH = APC_COMP_BOTH; /* compute entries for sm and rhsv */


/**************************************/
/* PDH_INTF.H PROCEDURES              */
/**************************************/
/* Rules:
/* - name always begins with pdr_ */
/* - argument names start uppercase */


/*------------------------------------------------------------
  pdr_get_list_ent - to return the list of integration entities - entities
                         for which stiffness matrices and load vectors are
                         provided by the FEM code to the solver module,
                         and DOF entities - entities with which there are dofs 
                         associated by the given approximation 
------------------------------------------------------------*/
int pdr_get_list_ent( /* returns: >=0 - success code, <0 - error code */
  int Problem_id,     /* in:  problem (and solver) identification */
  int* Nr_int_ent,    /* out: number of integration entitites */
  int** List_int_ent_type,/* out: list of types of integration entitites */
  int** List_int_ent_id,  /* out: list of IDs of integration entitites */
  int* Nr_dof_ent,    /* out: number of dof entities (entities with which there
		              are dofs associated by the given approximation) */
  int** List_dof_ent_type,/* out: list of types of integration entitites */
  int** List_dof_ent_id,  /* out: list of IDs of integration entitites */
  int** List_dof_ent_nrdofs,/* out: list of no of dofs for 'dof' entity */
  int* Nrdofs_glob,    /* out: global number of degrees of freedom (unknowns) */
  int* Max_dofs_per_dof_ent/* out: maximal number of dofs per dof entity */
  )
{

/*++++++++++++++++ executable statements ++++++++++++++++*/

// default procedure for returning lists of integration and dof entities
// for standard and discontinuous Galerkin approximations
  utr_get_list_ent(Problem_id, Nr_int_ent, 
		       List_int_ent_type, List_int_ent_id, 
                       Nr_dof_ent, List_dof_ent_type, 
		       List_dof_ent_id, List_dof_ent_nrdofs, 
		       Nrdofs_glob, Max_dofs_per_dof_ent);
  
  return(1);
}


/*------------------------------------------------------------
  pdr_get_list_ent_coarse - to return the list of integration entities - entities
                         for which stiffness matrices and load vectors are
                         provided by the FEM code to the solver module,
                         and DOF entities - entities with which there are dofs 
                         associated by the given approximation for COARSE level
                         given the corresponding lists from the fine level 
------------------------------------------------------------*/
int pdr_get_list_ent_coarse( /* returns: >=0 - success code, <0 - error code */
  int Problem_id,       /* in:  problem (and solver) identification */
  int Nr_int_ent_fine, /* in: number of integration entitites */
  int* List_int_ent_type_fine,/* in: list of types of integration entitites */
  int* List_int_ent_id_fine,  /* in: list of IDs of integration entitites */
  int Nr_dof_ent_fine, /* in: number of dof entities (entities with which there
		             are dofs associated by the given approximation) */
  int* List_dof_ent_type_fine,/* in: list of types of integration entitites */
  int* List_dof_ent_id_fine,  /* in: list of IDs of integration entitites */
  int* List_dof_ent_nrdofs_fine,/* in: list of no of dofs for 'dof' entity */
  int Nrdofs_glob_fine,  /* in: global number of degrees of freedom (unknowns) */
  int Max_dofs_per_dof_ent_fine, /* in: maximal number of dofs per dof entity */
  int* Pdeg_coarse_p,  /* in/out: degree of approximation for coarse space */
  int* Nr_int_ent_p,    /* out: number of integration entitites */
  int** List_int_ent_type,/* out: list of types of integration entitites */
  int** List_int_ent_id,  /* out: list of IDs of integration entitites */
  int* Nr_dof_ent_p,    /* out: number of dof entities (entities with which there
		              are dofs associated by the given approximation) */
  int** List_dof_ent_type,/* out: list of types of integration entitites */
  int** List_dof_ent_id,  /* out: list of IDs of integration entitites */
  int** List_dof_ent_nrdofs,/* out: list of no of dofs for 'dof' entity */
  int* Nrdofs_glob,    /* out: global number of degrees of freedom (unknowns) */
  int* Max_dofs_per_dof_ent/* out: maximal number of dofs per dof entity */
  )
{

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* associate the suitable field with the problem and the solver */
  int i=2; int mesh_id=pdr_ctrl_i_params(Problem_id,i);
  i=3; int field_id = pdr_ctrl_i_params(Problem_id,i);
  //i=6; solver_id = pdr_ctrl_i_params(Problem_id,i);
 
  
  if(*Pdeg_coarse_p<=0){
    int nel=0; int pdeg;
    while((nel=mmr_get_next_act_elem(mesh_id, nel))!=0){

      pdeg = apr_get_el_pdeg(field_id, nel, &pdeg);
      break;

    }
    *Pdeg_coarse_p = pdeg;
  }



// default procedure for returning lists of integration and dof entities
// for coarse mesh and standard and discontinuous Galerkin approximations
  utr_get_list_ent_coarse(Problem_id, Nr_int_ent_fine, 
			  List_int_ent_type_fine, List_int_ent_id_fine, 
			  Nr_dof_ent_fine, List_dof_ent_type_fine, 
			  List_dof_ent_id_fine, List_dof_ent_nrdofs_fine, 
			  Nrdofs_glob_fine, Max_dofs_per_dof_ent_fine,
			  Pdeg_coarse_p, Nr_int_ent_p,
			  List_int_ent_type, List_int_ent_id, 
			  Nr_dof_ent_p, List_dof_ent_type, 
			  List_dof_ent_id, List_dof_ent_nrdofs, 
			  Nrdofs_glob, Max_dofs_per_dof_ent);

  
  return(1);

}

/*------------------------------------------------------------
pdr_get_max_num_grid_levels - for limiting nr_levels in multigrid 
                              based on mesh and field data
------------------------------------------------------------*/
int pdr_get_max_num_grid_levels(
  int Problem_id
)
{
  return(utr_get_max_num_grid_levels(Problem_id));
}


/*------------------------------------------------------------
 pdr_create_assemble_stiff_mat - to create element stiffness matrices
                                 and assemble them to the global SM
------------------------------------------------------------*/
int pdr_create_assemble_stiff_mat(
  int Problem_id, 
  int Level_id, 
  int Comp_type,         /* in: indicator for the scope of computations: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int* Pdeg_coarse_p, // indicator or value for pdeg on coarse meshes
  int Nr_int_ent,
  int* L_int_ent_type,
  int* L_int_ent_id,
  int Nr_colors_elems, 
  int* L_color_index_elems,
  int Nr_colors_faces,
  int* L_color_index_faces,
  int* Asse_pos_first_dof_int_ent,
  int* Assembly_table,
  int Max_dofs_int_ent
)
{

  utr_create_assemble_stiff_mat(Problem_id, Level_id, Comp_type, Pdeg_coarse_p,
				Nr_int_ent, L_int_ent_type, L_int_ent_id, 
				Nr_colors_elems, L_color_index_elems,
				Nr_colors_faces, L_color_index_faces,
				Asse_pos_first_dof_int_ent, Assembly_table, 
				Max_dofs_int_ent);

  return(1);

}
 

/*------------------------------------------------------------
  pdr_assemble_local_stiff_mat - to assemble an element stiffness matrix
                                   to the global SM
------------------------------------------------------------*/
int pdr_assemble_local_stiff_mat( 
                         /* returns: >=0 - success code, <0 - error code */
  int Problem_id,        /* in: solver ID (used to identify the subproblem) */
  int Level_id,          /* in: level ID */
  int Comp_type,         /* in: indicator what was computed: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int Nr_dof_ent,        /* in: number of global dof blocks */
                         /*     associated with the local stiffness matrix */
  int* L_dof_ent_type,   /* in: list of dof blocks' IDs */
  int* L_dof_ent_id,     /* in: list of dof blocks' IDs */
  int* L_dof_ent_nrdofs, /* in: list of blocks' numbers of dof */
  double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /* in: rhs vector */
  char* Rewr_dofs        /* in: flag to rewrite or sum up entries */
                         /*   'T' - true, rewrite entries when assembling */
                         /*   'F' - false, sum up entries when assembling */
)
{
  int i=6;
  int solver_id = pdr_ctrl_i_params(Problem_id,i);

  sir_assemble_local_stiff_mat(solver_id, Level_id, Comp_type,
			       Nr_dof_ent, L_dof_ent_type,
			       L_dof_ent_id,L_dof_ent_nrdofs, 
			       Stiff_mat, Rhs_vect, Rewr_dofs);

  return(1);

}

/*------------------------------------------------------------
  pdr_assemble_local_stiff_mat_with_table - to assemble an element stiffness matrix
                                   to the global SM using assembly table
------------------------------------------------------------*/
int pdr_assemble_local_stiff_mat_with_table( 
                         /* returns: >=0 - success code, <0 - error code */
  int Problem_id,        /* in: solver ID (used to identify the subproblem) */
  int Level_id,          /* in: level ID */
  int Comp_type,         /* in: indicator what was computed: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int Nr_dof_ent,        /* in: number of global dof blocks */
                         /*     associated with the local stiffness matrix */
  int* Assembly_table_int_ent, /* part of the global assembly table */
  double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /* in: rhs vector */
  char* Rewr_dofs        /* in: flag to rewrite or sum up entries */
                         /*   'T' - true, rewrite entries when assembling */
                         /*   'F' - false, sum up entries when assembling */
)
{
  int i=6;
  int solver_id = pdr_ctrl_i_params(Problem_id,i);

  sir_assemble_local_stiff_mat_with_table(solver_id, Level_id, Comp_type,
			       Nr_dof_ent, Assembly_table_int_ent, 
			       Stiff_mat, Rhs_vect, Rewr_dofs);

  return(1);

}



/*------------------------------------------------------------
  pdr_comp_stiff_mat - to provide a solver with a stiffness matrix
                      and a load vector corresponding to the specified
                      mesh entity, together with information on how to
                      assemble entries into the global stiffness matrix
                      and the global load vector
------------------------------------------------------------*/
int pdr_comp_stiff_mat( /* returns: >=0 - success code, <0 - error code */
  int Problem_id,      /* in: approximation field ID  */
  int Int_ent_type,    /* in: unique identifier of the integration entity */ 
  int Int_ent_id,    /* in: unique identifier of the integration entity */ 
  int Comp_sm,       /* in: indicator for the scope of computations: */
                     /*   PDC_NO_COMP  - do not compute anything */
                     /*   PDC_COMP_SM - compute entries to stiff matrix only */
                     /*   PDC_COMP_RHS - compute entries to rhs vector only */
                     /*   PDC_COMP_BOTH - compute entries for sm and rhsv */
  int *Pdeg_vec,        /* in: enforced degree of polynomial (if != NULL ) */
  int* Nr_dof_ent,   /* in: size of arrays, */
                          /* out: number of mesh entities with which dofs and */
                          /*      stiffness matrix blocks are associated */
  int* List_dof_ent_type, /* out: list of types for 'dof' entities */
  int* List_dof_ent_id,   /* out: list of ids for 'dof' entities */
  int* List_dof_ent_nrdof,/* out: list of no of dofs for 'dof' entity */
  int* Nrdofs_loc,        /* in(optional): size of Stiff_mat and Rhs_vect */
                /* out(optional): actual number of dofs per integration entity*/
  double* Stiff_mat,      /* out(optional): stiffness matrix stored columnwise*/
  double* Rhs_vect,       /* out(optional): rhs vector */
  char* Rewr_dofs         /* out(optional): flag to rewrite or sum up entries */
                          /*   'T' - true, rewrite entries when assembling */
                          /*   'F' - false, sum up entries when assembling */
  )
{

  int pdeg;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* in higher order approximation pdeg may be a vector - in DG single number */
  if(Pdeg_vec==NULL) pdeg = -10; // use pdeg from approximation on finest grid
  else pdeg = Pdeg_vec[0];

/*kbw
    printf("in pdr_comp_stiff_mat: problem_id %d, int_ent: type %d, id %d, pdeg %d\n", 
	   Problem_id, Int_ent_type, Int_ent_id, pdeg);
/*kew*/

  if( Int_ent_type == PDC_ELEMENT ){

    //#pragma omp critical (pdr_comp_el_stiff_mat)
    pdr_conv_diff_comp_el_stiff_mat(Problem_id, Int_ent_id, Comp_sm, pdeg,
			  Nr_dof_ent, List_dof_ent_type, 
			  List_dof_ent_id, List_dof_ent_nrdof,
			  Nrdofs_loc, Stiff_mat, Rhs_vect, Rewr_dofs); 


  }
  else if(Int_ent_type==PDC_FACE){

    pdr_conv_diff_comp_fa_stiff_mat(Problem_id, Int_ent_id, Comp_sm, pdeg,
			  Nr_dof_ent, List_dof_ent_type, 
			  List_dof_ent_id, List_dof_ent_nrdof,
			  Nrdofs_loc, Stiff_mat, Rhs_vect, Rewr_dofs); 

  }
  else{
    printf("Wrong integration entity type in pdr_comp_stiff_mat!\n");
    exit(-1);
  }

  /*kbw
  //if(Comp_sm!=PDC_NO_COMP)
	  {
		int i;
    printf("In pdr_com_el_stiff_mat: field_id %d, El_id %d, Comp_sm %d, Nr_dof_ent %d\n",
     field_id, Int_ent_id, Comp_sm, *Nr_dof_ent);
     printf("For each block: \ttype, \tid, \tnrdof\n");
     for(i=0;i<*Nr_dof_ent;i++){
     printf("\t\t\t%d\t%d\t%d\n",
     List_dof_ent_type[i],List_dof_ent_id[i],List_dof_ent_nrdofs[i]);
     }
     printf("\n\n");
     }
     //getchar();getchar();
     /*kew */

  /* matrix displayed by rows, altghough stored by columns !!!!!!!!! */
  /*kbw
     if(Comp_sm!=PDC_NO_COMP)
     {
     printf("\nElement %d: Modified stiffness matrix:\n",El_id);
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
     printf("\n\n");
     //getchar();
     }
     /* */


  return(1);
}


/**-----------------------------------------------------------
  pdr_comp_stiff_mat_uncon - to construct UNCONSTRAINED stiffness matrix and 
                      a load vector for some given mesh entity
------------------------------------------------------------*/
int pdr_comp_stiff_mat_uncon( /** returns: >=0 - success code, <0 - error code */
  int Problem_id,     /** in: approximation field ID  */
  int Int_ent_type,    /** in: type of the integration entity */ 
  int Int_ent_id,    /** in: unique identifier of the integration entity */ 
  int Comp_sm,       /** in: indicator for the scope of computations: */
                     /**   PDC_NO_COMP  - do not compute anything */
                     /**   PDC_COMP_SM - compute entries to stiff matrix only */
                     /**   PDC_COMP_RHS - compute entries to rhs vector only */
                     /**   PDC_COMP_BOTH - compute entries for sm and rhsv */
  int *Pdeg_vec,        /** in: enforced degree of polynomial (if != NULL ) */
  int* Nrdofs_loc,        /** in(optional): size of Stiff_mat and Rhs_vect */
                /** out(optional): actual number of dofs per integration entity*/
  double* Stiff_mat,      /** out(optional): stiffness matrix stored columnwise*/
  double* Rhs_vect,       /** out(optional): rhs vector */
  char* Rewr_dofs        /** out(optional): flag to rewrite or sum up entries */
                         /**   'T' - true, rewrite entries when assembling */
                         /**   'F' - false, sum up entries when assembling */
  )
{

    printf("Due to DG/STD splitting pdr_comp_stiff_mat_uncon not implemented!\n");
    exit(-1);


}

/*------------------------------------------------------------
  pdr_select_el_coeff_vect - to select coefficients returned to approximation
                        routines for element integrals in weak formulation
           (the procedure indicates which terms are non-zero in weak form)
------------------------------------------------------------*/
int pdr_select_el_coeff_vect(  // returns success indicator
  int Problem_id,
  int *Coeff_vect_ind	/* out: coefficient indicator */

  )
{

  pdr_conv_diff_select_el_coeff_vect(Problem_id, Coeff_vect_ind);

  return(1);

}


/*!!!!!! OLD OBSOLETE VERSION !!!!!!*/
/*------------------------------------------------------------
  pdr_select_el_coeff - to select coefficients returned to approximation
                        routines for element integrals in weak formulation
           (the procedure indicates which terms are non-zero in weak form)
------------------------------------------------------------*/
double *pdr_select_el_coeff( /* returns: pointer !=NULL to indicate selection */
  int Problem_id,
  double **Mval,	/* out: mass matrix coefficient */
  double **Axx,double **Axy,double **Axz, /* out:diffusion coefficients, e.g.*/
  double **Ayx,double **Ayy,double **Ayz, /* Axy denotes scalar or matrix */
  double **Azx,double **Azy,double **Azz, /* related to terms with dv/dx*du/dy */
  /* second order derivatives in weak formulation (scalar for scalar problems */
  /* matrix for vector problems) */
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


  double *select_coeff=NULL;

  /* allocate storage for coefficients and select the needed ones */
  select_coeff=pdr_conv_diff_select_el_coeff(Problem_id, Mval,
	 Axx,Axy,Axz,Ayx,Ayy,Ayz,Azx,Azy,Azz,
	 Bx,By,Bz,Tx,Ty,Tz,Cval,Lval,Qx,Qy,Qz,Sval);

  return(select_coeff);

}


/*------------------------------------------------------------
pdr_el_coeff - to return coefficients for internal integrals
------------------------------------------------------------*/
int pdr_el_coeff(
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
  double* Uk_val, 	/* in: computed solution from previous iteration */
  double* Uk_x, 	/* in: gradient of computed solution Uk_val */
  double* Uk_y,   	/* in: gradient of computed solution Uk_val */
  double* Uk_z,   	/* in: gradient of computed solution Uk_val */
  double* Un_val, 	/* in: computed solution from previous time step */
  double* Un_x, 	/* in: gradient of computed solution Un_val */
  double* Un_y,   	/* in: gradient of computed solution Un_val */
  double* Un_z,   	/* in: gradient of computed solution Un_val */
  double* Mval,	/* out: mass matrix coefficient */
  double *Axx, double *Axy, double *Axz,  /* out:diffusion coefficients */
  double *Ayx, double *Ayy, double *Ayz,  /* e.g. Axy denotes scalar or matrix */
  double *Azx, double *Azy, double *Azz,  /* related to terms with dv/dx*du/dy */
  /* second order derivatives in weak formulation (scalar for scalar problems */
  /* matrix for vector problems) */
  double *Bx, double *By, double *Bz,	/* out: convection coefficients */
  /* Bx denotes scalar or matrix related to terms with du/dx*v in weak form */
  double *Tx, double *Ty, double *Tz,	/* out: convection coefficients */
  /* Tx denotes scalar or matrix related to terms with u*dv/dx in weak form */
  double *Cval,	/* out: reaction coefficients - for terms without derivatives */
  /*  in weak form (as usual: scalar for scalar problems, matrix for vectors) */
  double *Lval,	/* out: rhs coefficient for time term, Lval denotes scalar */
  /* or matrix corresponding to time derivative - similar as mass matrix but  */
  /* with known solution at the previous time step (usually denoted by u_n) */
  double *Qx, /* out: rhs coefficients for terms with derivatives */
  double *Qy, /* Qy denotes scalar or matrix corresponding to terms with dv/dy */
  double *Qz, /* derivatives in weak formulation */
  double *Sval	/* out: rhs coefficients without derivatives (source terms) */
  )
{

  /* get coefficients for conv_diff weak formulation */
  pdr_conv_diff_el_coeff(Problem_id, Elem, Mat_num, Hsize, Pdeg, X_loc,
		 Base_phi, Base_dphix, Base_dphiy, Base_dphiz,
		 Xcoor, Uk_val, Uk_x, Uk_y, Uk_z, Un_val, Un_x, Un_y, Un_z,
		 Mval, Axx, Axy, Axz, Ayx, Ayy, Ayz, Azx, Azy, Azz,
		 Bx, By, Bz, Tx, Ty, Tz, Cval, Lval, Qx, Qy, Qz, Sval);

  return (0);
}



/*------------------------------------------------------------
  pdr_read_sol_dofs - to read a vector of dofs associated with a given 
                   mesh entity from approximation field data structure
------------------------------------------------------------*/
int pdr_read_sol_dofs(/* returns: >=0 - success code, <0 - error code */
  int Problem_id,     /* in: solver ID (used to identify the subproblem) */
  int Dof_ent_type,
  int Dof_ent_id,
  int Nrdof,
  double* Vect_dofs  /* in: dofs to be written */
  )
{

  int i;
  int vect_id = 0;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /*kbw
  printf("in pdr_read_sol_dofs before apr_read_ent_dofs\n");
  printf("problem_id %d, Dof_ent_type %d, Dof_ent_id %d, nrdof %d\n",
	 Problem_id, Dof_ent_type, Dof_ent_id,  Nrdof);
  /*kew*/

  i=3; int field_id=pdr_ctrl_i_params(Problem_id,i);
  apr_read_ent_dofs(field_id, Dof_ent_type, Dof_ent_id, 
		    Nrdof, vect_id, Vect_dofs);

  /*kbw
  printf("in pdr_read_sol_dofs after apr_read_ent_dofs\n");
  for(i=0;i<Nrdof;i++){
    printf("%10.6lf",Vect_dofs[i]);
  }
  printf("\n");
  /*kew*/


  return(1);
}

/*------------------------------------------------------------
  pdr_write_sol_dofs - to write a vector of dofs associated with a given 
                   mesh entity to approximation field data structure
------------------------------------------------------------*/
int pdr_write_sol_dofs(/* returns: >=0 - success code, <0 - error code */
  int Problem_id,     /* in: solver ID (used to identify the subproblem) */
  int Dof_ent_type,
  int Dof_ent_id,
  int Nrdof,
  double* Vect_dofs  /* in: dofs to be written */
  )
{

  int vect_id = 0;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  int field_id=pdr_ctrl_i_params(Problem_id,3);
  apr_write_ent_dofs(field_id, Dof_ent_type, Dof_ent_id, 
		     Nrdof, vect_id, Vect_dofs);

  return(1);
}


/*---------------------------------------------------------
  pdr_L2_proj_sol - to project solution between elements of different generations
----------------------------------------------------------*/
int pdr_L2_proj_sol(
  int Problem_id, /* in: problem ID */
  int El,		/* in: element number */
  int *Pdeg_vec,	/* in: element degree of approximation */
  double* Dofs,	/* out: workspace for degress of freedom of El */
  /* 	NULL - write to  data structure */
  int* El_from,	/* in: list of elements to provide function */
  int *Pdeg_vec_from,	/* in: degree of polynomial for each El_from */
  double* Dofs_from /* in: Dofs of El_from or...*/
  )
{

  int field_id = pdr_ctrl_i_params(Problem_id,3);
  int i= -1; /* mode: -1 - projection from father to son */

  apr_L2_proj(field_id, i, El, Pdeg_vec, Dofs, 
	                   El_from, Pdeg_vec_from, Dofs_from, NULL);

  return(1);

}


/*---------------------------------------------------------
pdr_renum_coeff - to return a coefficient being a basis for renumbering
----------------------------------------------------------*/
int pdr_renum_coeff(
  int Problem_id, /* in: problem ID */
  int Ent_type,	/* in: type of mesh entity */
  int Ent_id,	/* in: mesh entity ID */
  double* Ren_coeff  /* out: renumbering coefficient */
	)
{


/*++++++++++++++++ executable statements ++++++++++++++++*/

// no problem dependent renumbering
  *Ren_coeff = 1.0;
  return(1);

}


/*------------------------------------------------------------
  pdr_get_ent_pdeg - to return the degree of approximation index 
                      associated with a given mesh entity
------------------------------------------------------------*/
extern int pdr_get_ent_pdeg( /* returns: >0 - approximation index,
                                          <0 - error code */
  int Problem_id,     /* in: approximation field ID  */
  int Ent_type,      /* in: type of mesh entity */
  int Ent_id         /* in: mesh entity ID */
  )
{

  int field_id = pdr_ctrl_i_params(Problem_id,3);
  return(apr_get_ent_pdeg(field_id, Ent_type, Ent_id));

}


/*---------------------------------------------------------
  pdr_dof_ent_sons - to return info whether the entity is owned
                     and a list of dof entity sons for owned entity
---------------------------------------------------------*/
int pdr_dof_ent_sons( /* returns: info whether the entity is owned */
  int Problem_id,     /* in: problem ID  */
  int Ent_type,      /* in: type of mesh entity */
  int Ent_id,        /* in: mesh entity ID */
  int *Ent_sons     /* out: list of dof entity sons (for owned entities) */
               	     /* 	Ent_sons[0] - number of sons */
  )
{
  return( utr_dof_ent_sons( Problem_id, Ent_type, Ent_id, Ent_sons) );
}


/*---------------------------------------------------------
  pdr_proj_sol_lev - to L2 project solution dofs between mesh levels
---------------------------------------------------------*/
int pdr_proj_sol_lev( /* returns: >=0 - success; <0 - error code*/
  int Problem_id, /* in: problem ID */
  int Solver_id,        /* in: solver data structure to be used */
  int Ilev_from,    /* in: level number to project from */
  double* Vec_from, /* in: vector of values to project */
  int Ilev_to,      /* in: level number to project to */
  double* Vec_to    /* out: vector of projected values */
  )
{

  // currently not used
  return(0);
}

/*---------------------------------------------------------
  pdr_vec_norm - to compute a norm of global vector (in parallel)
---------------------------------------------------------*/
double pdr_vec_norm( /* returns: L2 norm of global Vector */
  int Problem_id, /* in: problem ID */
  int Solver_id,        /* in: solver data structure to be used */
  int Level_id,         /* in: level number */
  int Nrdof,            /* in: number of vector components */
  double* Vector        /* in: local part of global Vector */
  )
{

  const int IONE=1;
  double vec_norm = 0.0;
  int i, field_id;

/*++++++++++++++++ executable statements ++++++++++++++++*/
  
#ifdef PARALLEL
  // simplified setting, only one problem and one field
  vec_norm = appr_sol_vec_norm(pdv_exchange_table_index, 
			       Level_id, Nrdof, Vector);
#endif

  return(vec_norm);
}


/*---------------------------------------------------------
  pdr_sc_prod - to compute a scalar product of two global vectors 
---------------------------------------------------------*/
double pdr_sc_prod( /* retruns: scalar product of Vector1 and Vector2 */
  int Problem_id, /* in: problem ID */
  int Solver_id,        /* in: solver data structure to be used */
  int Level_id,         /* in: level number */
  int Nrdof,           /* in: number of vector components */
  double* Vector1,     /* in: local part of global Vector */
  double* Vector2      /* in: local part of global Vector */
  )
{

  const int IONE=1;
  double sc_prod = 0.0;
  int i, field_id;

/*++++++++++++++++ executable statements ++++++++++++++++*/
  
#ifdef PARALLEL
  // simplified setting, only one problem and one field
  sc_prod = appr_sol_sc_prod(pdv_exchange_table_index, 
			     Level_id, Nrdof, Vector1, Vector2);
#endif

  return(sc_prod);
}

/*---------------------------------------------------------
  pdr_create_exchange_tables - to create tables to exchange dofs 
---------------------------------------------------------*/
int pdr_create_exchange_tables( 
                      /* returns: >=0 -success code, <0 -error code */
  int Problem_id, /* in: problem ID */
  int Solver_id,        /* in: solver data structure to be used */
  int Level_id,       /* in: level ID */
  int Nr_dof_ent,     /* in: number of DOF entities in the level */
  /* all four subsequent arrays are indexed by block IDs with 1(!!!) offset */
  int* L_dof_ent_type,/* in: list of DOF entities associated with DOF blocks */
  int* L_dof_ent_id,  /* in: list of DOF entities associated with DOF blocks */
  int* L_bl_nrdof,    /* in: list of nrdofs for each dof block */
  int* L_bl_posg,     /* in: list of positions within the global */
                      /*     vector of dofs for each dof block */
  int* L_elem_to_bl,  /* in: list of DOF blocks associated with DOF entities */
  int* L_face_to_bl,  /* in: list of DOF blocks associated with DOF entities */
  int* L_edge_to_bl,  /* in: list of DOF blocks associated with DOF entities */
  int* L_vert_to_bl  /* in: list of DOF blocks associated with DOF entities */
  )
{

  int i, field_id;

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef PARALLEL
  // simplified setting, only one problem and one field
  appr_create_exchange_tables(pdv_exchange_table_index, 
			      Level_id, Nr_dof_ent, L_dof_ent_type,
			      L_dof_ent_id, L_bl_nrdof, L_bl_posg, 
			      L_elem_to_bl, L_face_to_bl,
			      L_edge_to_bl, L_vert_to_bl);
#endif

  return(0);

}

/*---------------------------------------------------------
  pdr_exchange_dofs - to exchange dofs between processors
---------------------------------------------------------*/
int pdr_exchange_dofs(
  int Problem_id, /* in: problem ID */
  int Solver_id,        /* in: solver data structure to be used */
  int Level_id,         /* in: level number */
  double* Vec_dofs  /* in: vector of dofs to be exchanged */
)
{

  int i, field_id;

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef PARALLEL
  // simplified setting, only one problem and one field
  appr_exchange_dofs(pdv_exchange_table_index, Level_id, Vec_dofs);
#endif

  return(0);
}


