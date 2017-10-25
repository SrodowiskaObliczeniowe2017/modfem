/************************************************************************
File aps_std_quad_intf.c - implementation of the approximation module
                  interface for the standard quadratic discretization
                  of PDs on 40D meshes with prismatic elements

Contains definitions of routines:
  apr_module_introduce - to return the approximation method's name
  apr_init_field - to initiate new approximation field and read its data
  apr_write_field - to dump-out field data in the standard MOD_FEM format
  apr_check_field - to check approximation field data structure
  apr_get_mesh_id - to return the ID of the associated mesh
  apr_get_nreq - to return the number of components in solution vector
  apr_get_nr_sol - to return the number of solution vectors stored
  apr_get_base_type - to return the type of basis functions
  apr_get_ent_pdeg - to return the degree of approximation symbol
                      associated with a given mesh entity
  apr_set_ent_pdeg - to set the degree of approximation index 
                      associated with a given mesh entity
  apr_get_ent_numshap - to return the number of shape functions (vector DOFs)
                        associated with a given mesh entity
  apr_get_ent_nrdofs - to return the number of dofs associated with 
                      a given mesh entity
  apr_get_el_pdeg - to return the degree of approximation vector 
                      associated with a given element
  apr_set_el_pdeg - to set the degree of approximation vector 
                      associated with a given element
  apr_get_el_pdeg_numshap - to return the number of shape functions
                            (scalar DOFs) for an element given its 
                            degree of approximation symbol or vector pdeg
  apr_get_el_dofs - to return the number and the list of element's degrees
                        of freedom (corresponding to standard shape functions)
  apr_get_nrdofs_glob - to return a global dimension of the problem
  apr_read_ent_dofs - to read a vector of dofs associated with a given
                  mesh entity from approximation field data structure
  apr_write_ent_dofs - to write a vector of dofs associated with a given
                   mesh entity to approximation field data structure
  apr_create_ent_dofs - to write a vector of dofs associated with a given
                   mesh entity to approximation field data structure
  apr_set_ini_con - to set an initial condition
  apr_num_int_el - to perform numerical integration for an element
  apr_num_int_fa - to perform numerical integration for a face
  apr_get_stiff_mat_data - to return data on dof entities for an element and 
                      to compute or rewrite element's stiffness matrix and RHSV
  apr_proj_dofref - to rewrite dofs after modifying the mesh
  apr_rewr_sol - to rewrite solution from one vector to another
  apr_free_field - to free approximation field data structure

  apr_limit_deref - to return whether derefinement is allowed or not
  apr_limit_ref - to return whether refinement is allowed or not
    the last two routines does not use approximation data structures but 
    the result depends on the approximation method
  apr_refine - to refine an element or the whole mesh checking mesh irregularity
  apr_derefine - to derefine an element or the whole mesh with irregularity check

------------------------------
History:
	02.2002 - Krzysztof Banas, initial version
	10.2013 - Jan Bielański, implementation quadratic approximation for tetrahedron and prisms without constrained nodes
	11.2013 - Jan Bielański, implementation quadratic approximation for tetrahedron and prisms with constrained nodes
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<assert.h>
/* interface of the mesh manipulation module */
#include "mmh_intf.h"

/* interface for all approximation modules */
// moved to local header where we redifne some constants
//#include "aph_intf.h"

#include "uth_intf.h"

/* internal header file for the std lin approximation module */
#include "./aph_std_quad.h"

/* debugging/logging/checking macros */
#include "uth_log.h"

#define ut_min(x,y) ((x)<(y)?(x):(y))

#ifndef NULL
#define NULL NULL
#endif

#ifndef SAFE_FREE
#define SAFE_FREE(ptr) if(ptr!=NULL) free(ptr); ptr=NULL;
#endif

/* coefficients of constraint nodes/edges */
extern const struct {
  int id[7];
  double val[15];
} constr_coeff;

/* two functions needed from problem dependent module */
/*------------------------------------------------------------
  pdr_select_el_coeff_vect - to select coefficients returned to approximation
                        routines for element integrals in weak formulation
           (the procedure indicates which terms are non-zero in weak form)
------------------------------------------------------------*/
int pdr_select_el_coeff_vect( // returns success indicator
  int Problem_id,
  int *Coeff_vect_ind	/* out: coefficient indicator */
  // input to pdr_select_el_coeff_vect:
  // 0 - perform substitutions in coeff_vect_ind only if coeff_vect_ind[0]==0
  // output from pdr_select_el_coeff_vect:
  // 0 - coeff_vect_ind[0]==1 - always
  // 1 - mval, 2 - axx, 3 - axy, 4 - axz, 5 - ayx, 6 - ayy, 7 - ayz, 
  // 8 - azx, 9 - azy, 10 - azz, 11 - bx, 12 - by, 13 - bz
  // 14 - tx, 15 - ty, 16 - tz, 17 - cval
  // 18 - lval, 19 - qx, 20 - qy, 21 - qz, 22 - sval 
			      );


/*!!!!!! OLD OBSOLETE VERSION !!!!!!*/
/*------------------------------------------------------------
  pdr_select_el_coeff - to select coefficients returned to approximation
                        routines for element integrals
------------------------------------------------------------*/
double* pdr_select_el_coeff( /* returns: pointer !=NULL to indicate selection */
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
  );

/*---------------------------------------------------------
  pdr_el_coeff - to return coefficients for element integrals
----------------------------------------------------------*/
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
  double *Uk_val, 	/* in: computed solution from previous iteration */
  double *Uk_x, 	/* in: gradient of computed solution Uk_val */
  double *Uk_y,   	/* in: gradient of computed solution Uk_val */
  double *Uk_z,   	/* in: gradient of computed solution Uk_val */
  double *Un_val, 	/* in: computed solution from previous time step */
  double *Un_x, 	/* in: gradient of computed solution Un_val */
  double *Un_y,   	/* in: gradient of computed solution Un_val */
  double *Un_z,   	/* in: gradient of computed solution Un_val */
  double *Mval,	/* out: mass matrix coefficient */
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
		 );
/* end of forward declaration of functions from problem dependent module */


/*** GLOBAL VARIABLES for the whole module ***/

int       apv_nr_fields=0;     /* the number of fields in the problem */
int       apv_cur_field_id;              /* ID of the current field */
apt_field  apv_fields[APC_MAX_NUM_FIELD];        /* array of fields */


/*------------------------------------------------------------
 apr_set_and_test_pdeg - set and test pdeg value exit if failed
------------------------------------------------------------*/
extern void apr_set_and_test_pdeg( /* exit if failed */
				  int Field_id,
				  int *Pdeg_in
);

/*------------------------------------------------------------
  apr_module_introduce - to return the approximation method's name
------------------------------------------------------------*/
int apr_module_introduce(
                  /* returns: >=0 - success code, <0 - error code */
  char* Approx_name /* out: the name of the approximation method */
  )
{
  const char string[] = "STANDARD_QUADRATIC";
  strncpy(Approx_name,string,sizeof(string));
  return(1);
}


/*------------------------------------------------------------
  apr_init_field - to initiate new approximation field and
                   read its control parameters
------------------------------------------------------------*/
int apr_init_field(  /* returns: >0 - field ID, <0 - error code */
  char Field_type,
  int Control,	 /* in: control variable: */
                 /*      APC_ZERO - to initialize the field to zero */
                 /*      APC_READ - to read field values from the file */
                 /*      APC_INIT - to initialize the field using function */
                 /*                 provided by the problem dependent module */
  int Mesh_id,	 /* in: ID of the corresponding mesh */
  int Nreq,	 /* in: number of equations - solution vector components */
  int Nr_sol,	 /* in: number of solution vectors for each dof entity */
  int Pdeg_in,   /* in: degree of approximating polynomial */
  char *Filename, /* in: name of the file to read approximation data */
  double (*Fun_p)(int, double*, int) /* pointer to function that provides */
		 /* problem dependent initial condition data */
  )
{
  /* pointer to field structure */
  apt_field *field_p;
  int nreq;
  int mxvert, nvert, ivert, nrvert_glob, nr_dof_vert, nr_constr_dof_vert, nrvert_constr_glob; // Variable for VERTEX data
  int mxedge, nedge, iedge, nredge_glob, nr_dof_edge, nr_constr_dof_edge, nredge_constr_glob; // Variable for EDGE data 
  int mxelem, nelem, ielem, nrelem_glob, nr_dof_elem; // Variable for ELEM data [future approximation]

  int nr_dof_ents; // nr_dof_vert + nr_dof_edge + nr_dof_face + nr_dof_elem

  /* auxiliary variables */
  int i; int error = 0; //iaux,,j,k,jaux
  FILE *fp;
  char mesh_name[100];

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* first check that the mesh is compatible */
  mmr_module_introduce(mesh_name);

#ifdef DEBUG_APM
  if(strcmp("3D_PRISM",mesh_name)!=0 && strcmp("3D_Hybrid",mesh_name)!=0 ){
    printf("Wrong mesh %s for standard quadratic approximation with prisms!\n",mesh_name);
    exit(-1);
  }
#endif
  /*  */

  /* initiate quadrature data for Gauss-Legendre integration */
  apr_gauss_init();

  /* increase the counter for fields */
  apv_nr_fields++;

  /* set the current field ID */
  apv_cur_field_id = apv_nr_fields;

  field_p = &apv_fields[apv_cur_field_id-1];
  field_p->mesh_id = Mesh_id;

  /* at the beginning we assume uniform mesh without constrained nodes */
  field_p->uniform = APC_TRUE;
  field_p->constr = APC_FALSE;
  
  field_p->nreq = Nreq;
  field_p->nr_sol = Nr_sol;

  /* set and test the current field PDEG */
  apr_set_and_test_pdeg(apv_cur_field_id,&Pdeg_in);

  /* set solution vector size */
  nreq = field_p->nreq;

  /* Check NREQ value */
  if(field_p->nreq>APC_MAXEQ) {
    printf("Requested number of equations %d greater than APC_MAXEQ=%d!!!\n",field_p->nreq, APC_MAXEQ);
    printf("Change APC_MAXEQ in aph_std_lin_prism.h and recompile the code.\n");
    printf("Exiting.\n");
    exit(-1);
  }
  
  /* depending on the control variable */
  if(Control==APC_READ) {

    /* open the input file with approximation parameters */
    if((fp=fopen(Filename,"r")) == NULL) {
      printf("Not found file '%s' with field data!!! Exiting.\n",Filename);
      exit(-1);
    }
	
    /* read the number of solution vectors for each dof entity and */
    /* the number of components in each vector = number of equations */
    fscanf(fp,"%d %d",&field_p->nreq, &field_p->nr_sol);

#ifdef DEBUG_APM
    printf("(apr_init_field) field_p->nreq = %d ; field_p->nr_sol = %d\n",field_p->nreq, field_p->nr_sol);
#endif
	
#ifdef DEBUG_APM
    if(field_p->nreq!=Nreq) {      
      printf("Requested number of equations %d different than read from file %d!!!\n",Nreq, field_p->nreq);
      printf("Requested number of components %d different than read from file %d!!!\n",Nr_sol, field_p->nr_sol);
      printf("Exiting.\n");
      exit(-1);
    }
#endif
    
  }

  /* Initialize data structures */
  field_p->capacity_dof_vert = 0;
  field_p->capacity_dof_edge = 0;
  field_p->capacity_dof_elem = 0;
  field_p->dof_vert = NULL;
  field_p->dof_edge = NULL;
  field_p->dof_elem = NULL;

  /* get the maximal admissible number of vertices */
  mxvert=mmr_get_max_node_max(Mesh_id);
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    mxedge=mmr_get_max_edge_max(Mesh_id);

    // Future changes for higher order approximation or lagrange shape fuctions
    //mxelem = ???; 
  }

  /* Create apt_dof_ent structure */
  field_p->capacity_dof_vert=mxvert;
  field_p->dof_vert = (apt_dof_ent *) malloc((mxvert+1)*sizeof(apt_dof_ent)); // allocation space for vertexes
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    field_p->capacity_dof_edge=mxedge;
    field_p->dof_edge = (apt_dof_ent *) malloc((mxedge+1)*sizeof(apt_dof_ent)); // allocation space for edges

    // Future changes for higher order approximation or lagrange shape fuctions
    //field_p->capacity_dof_elem=mxelem;
    //field_p->dof_elem = (apt_dof_ent *) malloc((mxelem+1)*sizeof(apt_dof_ent)); // allocation space for elements
  }

  /* Check structures */
  if(field_p->dof_vert==NULL) {
    SAFE_FREE(field_p->dof_edge);
    SAFE_FREE(field_p->dof_elem);
    printf("Dofs structures for vertex and edges not allocated\n");
    exit(-1);
  }
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG && field_p->dof_edge==NULL) {
    SAFE_FREE(field_p->dof_vert);
    SAFE_FREE(field_p->dof_elem);
    printf("Dofs structures for vertex and edges not allocated\n");
    exit(-1);
  }

  /* Initialize apt_dof_ents structures for field */
  for(nvert=1;nvert<=mxvert;nvert++){
    field_p->dof_vert[nvert].vec_dof_1 = NULL;
    field_p->dof_vert[nvert].vec_dof_2 = NULL;
    field_p->dof_vert[nvert].vec_dof_3 = NULL;
    field_p->dof_vert[nvert].constr = NULL;
    field_p->dof_vert[nvert].constr_type = NULL;
    field_p->dof_vert[nvert].active = 0;
  }

  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    for(nedge=1;nedge<=mxedge;nedge++){
      field_p->dof_edge[nedge].vec_dof_1 = NULL;
      field_p->dof_edge[nedge].vec_dof_2 = NULL;
      field_p->dof_edge[nedge].vec_dof_3 = NULL;
      field_p->dof_edge[nedge].constr = NULL;
      field_p->dof_edge[nedge].constr_type = NULL;
      field_p->dof_edge[nedge].active = 0;
    }

    // Future changes for higher order approximation or lagrange shape fuctions

    //for(nelem=1;nelem<=mxelem;nelem++){
    //  field_p->dof_elem[nelem].vec_dof_1 = NULL;
    //  field_p->dof_elem[nelem].vec_dof_2 = NULL;
    //  field_p->dof_elem[nelem].vec_dof_3 = NULL;
    //  field_p->dof_elem[nelem].constr = NULL;
    //  field_p->dof_elem[nelem].constr_type = NULL;
    //  field_p->dof_elem[nelem].inactive_used_with_constr = 0;
    //}
  }

  #ifdef DEBUG_APM
  printf("\n(apr_init_field) -> Approximation type: ");
  switch(field_p->pdeg)
    {
    case APC_LINEAR_APPROXIMATION_PDEG:
      printf("LINEAR\n"); break;
    case APC_QUADRATIC_APPROXIMATION_PDEG:
      printf("QUADRATIC\n"); break;
    default:
      printf("-\n"); break;

    };
  #endif

  /* depending on the control variable */
  if(field_p->pdeg == APC_LINEAR_APPROXIMATION_PDEG) {
    if(Control==APC_READ){
      fscanf(fp,"%d ", &nrvert_glob); nr_dof_ents = nrvert_glob;
      fscanf(fp,"%d ", &nrvert_constr_glob);

#ifdef DEBUG_APM
      printf("\n(apr_init_field) nrvert_glob = %d ; nrvert_constr_glob = %d ;; mmr_get_max_node_id = %d\n",nrvert_glob,nrvert_constr_glob,mmr_get_max_node_id(Mesh_id));
#endif

      assert((nrvert_glob+nrvert_constr_glob) <= mmr_get_max_node_id(Mesh_id));

    }
    else{
      nrvert_glob=mmr_get_max_node_id(Mesh_id);
      nrvert_constr_glob=0;
    }
  }
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    if(Control==APC_READ){
      fscanf(fp,"%d %d %d", &nr_dof_ents, &nrvert_glob, &nredge_glob);
      fscanf(fp,"%d %d", &nrvert_constr_glob, &nredge_constr_glob);

#ifdef DEBUG_APM
      printf("\n(apr_init_field) nr_dof_ents = %d ; nrvert_glob = %d ; nredge_glob = %d\n",nr_dof_ents, nrvert_glob, nredge_glob);
      printf("\n(apr_init_field) mmr_get_max_node_id = %d ; mmr_get_max_edge_id = %d\n",mmr_get_max_node_id(Mesh_id),mmr_get_max_edge_id(Mesh_id));
      printf("\n(apr_init_field) nrvert_constr_glob = %d ; nredge_constr_glob = %d\n", nrvert_constr_glob, nredge_constr_glob);
#endif

      assert((nrvert_glob+nrvert_constr_glob) <= mmr_get_max_node_id(Mesh_id));
      assert((nredge_glob+nredge_constr_glob) <= mmr_get_max_edge_id(Mesh_id));
      assert(nr_dof_ents == (nrvert_glob+nredge_glob));
    }
    else{
      nrvert_glob=mmr_get_max_node_id(Mesh_id);
      nredge_glob=mmr_get_max_edge_id(Mesh_id);
      nrvert_constr_glob=0;
      nredge_constr_glob=0;
    }
  }

  /* read data from file */
  if(Control==APC_READ){ 

    nr_dof_vert = 0;
    nr_constr_dof_vert = 0;

    /* create vertex data structure */
    for(ivert=1;ivert<=(nrvert_glob+nrvert_constr_glob);ivert++){
      fscanf(fp,"%d",&nvert);

      /*jbw
      if(nvert>0) printf("Read REAL VERT: %d\n",abs(nvert));
      else printf("Read CONSTR VERT: %d\n",abs(nvert));
      /*jbw*/


      field_p->dof_vert[abs(nvert)].active = 1;

      if(field_p->nr_sol>1){
	field_p->dof_vert[nvert].vec_dof_2 = (double *) malloc(nreq*sizeof(double));
      }
      if(field_p->nr_sol>2){
	field_p->dof_vert[nvert].vec_dof_3 = (double *) malloc(nreq*sizeof(double));
      }
      
      for(i=0;i<nreq;i++) {
	fscanf(fp,"%lg",&field_p->dof_vert[nvert].vec_dof_1[i]);
	//assert(field_p->dof_vert[nvert].vec_dof_1[i] > -10e200 && field_p->dof_vert[nvert].vec_dof_1[i] < 10e200);
      }
      if(field_p->nr_sol>1) {
	for(i=0;i<nreq;i++){
	  fscanf(fp,"%lg",&field_p->dof_vert[nvert].vec_dof_2[i]);
	  //assert(field_p->dof_vert[nvert].vec_dof_2[i] > -10e200 && field_p->dof_vert[nvert].vec_dof_2[i] < 10e200);
	}
      }
      if(field_p->nr_sol>2) {
	for(i=0;i<nreq;i++){
	  fscanf(fp,"%lg",&field_p->dof_vert[nvert].vec_dof_3[i]);
	  //assert(field_p->dof_vert[nvert].vec_dof_3[i] > -10e200 && field_p->dof_vert[nvert].vec_dof_3[i] < 10e200);
	}
      }
    }

    /* for active constrained vertices */
    if(nvert < 0) {
      int num_related_nodes; // Number of related nodes

      nr_constr_dof_vert++;
	
      nvert = abs(nvert); // Vertex id

      // Read number of related nodes
      fscanf(fp,"%d",&num_related_nodes);

      // Alloc memory
      field_p->dof_vert[nvert].constr = (int *) malloc((num_related_nodes+1)*sizeof(int));
      field_p->dof_vert[nvert].constr_type = (int *) malloc((num_related_nodes+1)*sizeof(int));

      //Read data
      field_p->dof_vert[nvert].constr[0] = num_related_nodes;
      field_p->dof_vert[nvert].constr_type[0] = 1;

      for(i=1; i<=num_related_nodes; i++) {
	fscanf(fp,"%d %d",&field_p->dof_vert[nvert].constr[i],&field_p->dof_vert[nvert].constr_type[i]);
      }
    }

    /* create edge data structure */
    if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {

      nr_dof_edge = 0;
      nr_constr_dof_edge = 0;
 
      for(iedge=1;iedge<=(nredge_glob+nredge_constr_glob);iedge++){

	fscanf(fp,"%d",&nedge);

	field_p->dof_edge[abs(nedge)].active = 1;

	/*jbw
	  if(nedge>0) printf("Read REAL EDGE: %d\n",abs(nedge));
	  else printf("Read CONSTR EDGE: %d\n",abs(nedge));
	/*jbw*/
	
	/* for active unconstrained edges */
	if(nedge > 0) {

	  nr_dof_edge++;

	  field_p->dof_edge[nedge].vec_dof_1 = (double *) malloc(nreq*sizeof(double));
	  
	  if(field_p->nr_sol>1){
	    field_p->dof_edge[nedge].vec_dof_2 = (double *) malloc(nreq*sizeof(double));
	  }
	  if(field_p->nr_sol>2){
	    field_p->dof_edge[nedge].vec_dof_3 = (double *) malloc(nreq*sizeof(double));
	  }
 	
	  for(i=0;i<nreq;i++) {
	    fscanf(fp,"%lg",&field_p->dof_edge[nedge].vec_dof_1[i]);
	    //assert(field_p->dof_edge[nedge].vec_dof_1[i] > -10e200 && field_p->dof_edge[nedge].vec_dof_1[i] < 10e200);
	  }
	  if(field_p->nr_sol>1) {
	    for(i=0;i<nreq;i++){
	      fscanf(fp,"%lg",&field_p->dof_edge[nedge].vec_dof_2[i]);
	      //assert(field_p->dof_edge[nedge].vec_dof_2[i] > -10e200 && field_p->dof_edge[nedge].vec_dof_2[i] < 10e200);
	    }
	  }
	  if(field_p->nr_sol>2){
	    for(i=0;i<nreq;i++){
	      fscanf(fp,"%lg",&field_p->dof_edge[nedge].vec_dof_3[i]);
	      //assert(field_p->dof_vert[nedge].vec_dof_3[i] > -10e200 && field_p->dof_edge[nedge].vec_dof_3[i] < 10e200);
	    }
	  }
	}
      
	/* for active constrained edges */
	if(nedge< 0) {
	  int num_related_nodes; // Number of related nodes
	  
	  nr_constr_dof_edge++;
	  nedge = abs(nedge); // Edge id

	  // Read number of related nodes
	  fscanf(fp,"%d",&num_related_nodes);

	  // Alloc memory
	  field_p->dof_edge[nedge].constr = (int *) malloc((num_related_nodes+1)*sizeof(int));
	  field_p->dof_edge[nedge].constr_type = (int *) malloc((num_related_nodes+1)*sizeof(int));

	  //Read data
	  field_p->dof_edge[nedge].constr[0] = num_related_nodes;
	  field_p->dof_edge[nedge].constr_type[0] = 1;

	  for(i=1; i<=num_related_nodes; i++) {
	    fscanf(fp,"%d %d",&field_p->dof_edge[nedge].constr[i],&field_p->dof_edge[nedge].constr_type[i]);
	  }
	}
      }
    }
  } else {

    /* iniiate the field to zero */
    nr_dof_vert = 0;
    nr_dof_ents = 0;
    nr_constr_dof_vert = 0;
    
    /* create vertex data structure */
    for(nvert=1;nvert<=nrvert_glob;nvert++){
    
      /* for all active vertices */
      if(mmr_node_status(Mesh_id,nvert)==MMC_ACTIVE){

	nr_dof_ents++; nr_dof_vert++;

	field_p->dof_vert[nvert].active = 1;

	field_p->dof_vert[nvert].vec_dof_1 = (double *) malloc(nreq*sizeof(double));

	if(field_p->nr_sol>1){
	  field_p->dof_vert[nvert].vec_dof_2 = (double *) malloc(nreq*sizeof(double));
	}
	if(field_p->nr_sol>2){
	  field_p->dof_vert[nvert].vec_dof_3 = (double *) malloc(nreq*sizeof(double));
	}
      
	/*kbw
	  printf("setting ivert %d, nrvert_glob %d\n", nvert, nrvert_glob);
	/*kew*/
	
	for(i=0;i<nreq;i++){
	  field_p->dof_vert[nvert].vec_dof_1[i]=0.0;
	}
	if(field_p->nr_sol>1){
	  for(i=0;i<nreq;i++){
	    field_p->dof_vert[nvert].vec_dof_2[i]=0.0;
	  }
	}
	if(field_p->nr_sol>2){
	  for(i=0;i<nreq;i++){
	    field_p->dof_vert[nvert].vec_dof_3[i]=0.0;
	  }
	}
      }
    }
    
    /* create edges data structure */
    if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
      nr_dof_edge = 0;
      nr_constr_dof_edge = 0;
      
      for(nedge=1;nedge<=nredge_glob;nedge++){
    
	/* for all active edges */
	if(mmr_edge_status(Mesh_id,nedge)==MMC_ACTIVE){

	  nr_dof_ents++; nr_dof_edge++;

	  field_p->dof_edge[nedge].active = 1;
	  field_p->dof_edge[nedge].vec_dof_1 = (double *) malloc(nreq*sizeof(double));

	  if(field_p->nr_sol>1){
	    field_p->dof_edge[nedge].vec_dof_2 = (double *) malloc(nreq*sizeof(double));
	  }
	  if(field_p->nr_sol>2){
	    field_p->dof_edge[nedge].vec_dof_3 = (double *) malloc(nreq*sizeof(double));
	  }
	  
	  for(i=0;i<nreq;i++){
	    field_p->dof_edge[nedge].vec_dof_1[i]=0.0;
	  }
	  if(field_p->nr_sol>1){
	    for(i=0;i<nreq;i++){
	      field_p->dof_edge[nedge].vec_dof_2[i]=0.0;
	    }
	  }
	  if(field_p->nr_sol>2){
	    for(i=0;i<nreq;i++){
	      field_p->dof_edge[nedge].vec_dof_3[i]=0.0;
	    }
	  }
	} 
      }
    }
  } /* end if not read from file */

  /* number of dof entities does not include constrained nodes ! */
  field_p->nr_dof_ents = nr_dof_ents;
  field_p->nr_dof_vert = nr_dof_vert;
  field_p->nr_constr_dof_vert = nr_constr_dof_vert;
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    field_p->nr_dof_edge = nr_dof_edge;
    field_p->nr_constr_dof_edge = nr_constr_dof_edge;
  }

#ifdef DEBUG_APM
  if(field_p->pdeg == APC_LINEAR_APPROXIMATION_PDEG) {
    printf("\n\n LINEAR => field_id = %d field_p->nr_dof_ents = %d ; field_p->nr_dof_vert = %d",apv_cur_field_id,field_p->nr_dof_ents,field_p->nr_dof_vert);
  }
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    printf("\n\n QUADRATIC => field_id = %d field_p->nr_dof_ents = %d ; field_p->nr_dof_vert = %d ; field_p->nr_dof_edge = %d",apv_cur_field_id,field_p->nr_dof_ents,field_p->nr_dof_vert,field_p->nr_dof_edge);
  }
#endif
  
  /* close file with approximation data */
  if(Control==APC_READ) { fclose(fp); }

  /* initiate values using provided function */
  if(Control==APC_INIT) {
    apr_set_ini_con(apv_cur_field_id,Fun_p);
  }

  /*initiate constraints data*/
  //apr_create_constr_data(apv_cur_field_id);
  //if(Control==APC_READ) apr_create_constr_data(apv_cur_field_id);

#ifdef DEBUG_APM
  apr_check_field(apv_cur_field_id);
#endif

  return(apv_cur_field_id);
}

/*---------------------------------------------------------
  apr_write_field - to dump-out field data in the standard HP_FEM format
---------------------------------------------------------*/
int apr_write_field( /* returns: >=0 - success code, <0 - error code */
  int Field_id,    /* in: field ID */
  int Nreq,        /* in: number of equations (scalar dofs) */
  int Select,      /* in: parameter to select written vectors */
  double Accuracy, /* in: parameter specyfying accuracy - significant digits */
		   /* (put 0 for full accuracy in "%g" format) */  
  char *Filename   /* in: name of the file to write field data */
  )
{

  /* pointer to field structure */
  apt_field *field_p;
  
  int mesh_id;
  int nreq, nr_sol, selected[3];
  int mxvert, nvert, nr_dof_vert, nr_constr_dof_vert;
  int mxedge, nedge, nr_dof_edge, nr_constr_dof_edge;
  int nr_dof_ents;
  int i;
  FILE *fp;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);

  /* select the corresponding mesh */
  mesh_id = apr_get_mesh_id(Field_id);

  /* open the input file with approximation parameters */
  fp=fopen(Filename,"w");

  /* check which vectors to write to file */
  if(Select<=0 || Select>=7) {
    selected[0] = 1;
    selected[1] = 1;
    selected[2] = 1;
    nr_sol = field_p->nr_sol;
  }
  else {
    // scheme: 1 - 1, 2 - 2, 3 - 4, 1+2 - 3, 1+3 - 5, 2+3 - 6, 1+2+3 - 7
    selected[0] = Select%2;
    selected[1] = (Select/2)%2;
    selected[2] = (Select/4)%2;
    nr_sol=selected[0]+selected[1]+selected[2];
  }

  /* write the number of solution vectors for each dof entity */
  printf("writing field %d %d into file %s \n",field_p->nreq, nr_sol, Filename);
  fprintf(fp,"%d %d\n",field_p->nreq, nr_sol);
  nreq=apr_get_nreq(Field_id);

  /* loop over dofs data structure including free spaces */
  mxvert=mmr_get_max_node_id(mesh_id);
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    mxedge=mmr_get_max_edge_id(mesh_id);
  }

#ifdef DEBUG_APM
  nr_dof_ents=0;
  nr_dof_vert=0;
  nr_constr_dof_vert=0;
  for(nvert=1;nvert<=mxvert;nvert++){
    /* for active unconstrained vertices */
    if(mmr_node_status(mesh_id,nvert)==MMC_ACTIVE && field_p->dof_vert[nvert].vec_dof_1 != NULL) { 
      nr_dof_ents++; 
      nr_dof_vert++; 
    }

    /* for active constrained vertices */
    if(mmr_node_status(mesh_id,nvert)==MMC_ACTIVE && field_p->dof_vert[nvert].constr != NULL) {
      nr_constr_dof_vert++;
    }
  }
  if(nr_dof_vert!=field_p->nr_dof_vert){
    printf("wrong number of dof_vert in apr_write_field: %d %d !!! Exiting!\n",nr_dof_vert, field_p->nr_dof_vert);
    exit(-1);
  }
  if(nr_constr_dof_vert!=field_p->nr_constr_dof_vert) {
    printf("wrong number of constr_dof_vert in apr_write_field: %d %d !!! Exiting!\n",nr_constr_dof_vert, field_p->nr_constr_dof_vert);
    exit(-1);
  }

  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    nr_dof_edge=0;
    nr_constr_dof_edge=0;
    for(nedge=1;nedge<=mxedge;nedge++){
      /* for active unconstrained edges */
      if((mmr_edge_status(mesh_id,nedge)==MMC_ACTIVE || mmr_edge_status(mesh_id,nedge)==MMC_INACTIVE) && field_p->dof_edge[nedge].vec_dof_1 != NULL) { 
	nr_dof_ents++; 
	nr_dof_edge++; 
      }

      /* for active constrained edges */
      if(mmr_edge_status(mesh_id,nedge)==MMC_ACTIVE && field_p->dof_edge[nedge].constr != NULL) {
	nr_constr_dof_edge++;
      }
    }
    if(nr_dof_edge!=field_p->nr_dof_edge){
      printf("wrong number of dof_edge in apr_write_field: %d %d !!! Exiting!\n",nr_dof_edge, field_p->nr_dof_edge);
      exit(-1);
    }
    if(nr_constr_dof_edge!=field_p->nr_constr_dof_edge){
      printf("wrong number of constr_dof_edge in apr_write_field: %d %d !!! Exiting!\n",nr_constr_dof_edge, field_p->nr_constr_dof_edge);
      exit(-1);
    }
  }
  
  if(nr_dof_ents!=field_p->nr_dof_ents){
    printf("wrong number of dof_entities in apr_write_field! %d %d Exiting!\n",nr_dof_ents,field_p->nr_dof_ents);
    exit(-1);
  }
#endif

  nr_dof_ents=field_p->nr_dof_ents;
  nr_dof_vert=field_p->nr_dof_vert;
  nr_constr_dof_vert=field_p->nr_constr_dof_vert;
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    nr_dof_edge=field_p->nr_dof_edge;
    nr_constr_dof_edge=field_p->nr_constr_dof_edge;
  }

  /* write the number of dof entities */
  if(field_p->pdeg == APC_LINEAR_APPROXIMATION_PDEG) {
    fprintf(fp,"%d\n",nr_dof_ents);
    fprintf(fp,"%d\n",nr_constr_dof_vert);
  }
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    fprintf(fp,"%d %d %d\n",nr_dof_ents,nr_dof_vert,nr_dof_edge);
    fprintf(fp,"%d %d\n",nr_constr_dof_vert,nr_constr_dof_edge);
  }

  /* Save data for active vertex vertices */
  for(nvert=1;nvert<=mxvert;nvert++){

    /* for active unconstrained vertices */
    if(mmr_node_status(mesh_id,nvert)==MMC_ACTIVE && field_p->dof_vert[nvert].vec_dof_1 != NULL){

      fprintf(fp,"%d ",nvert);

      if(selected[0]==1){
	for(i=0;i<nreq;i++){
	  
	  if(fabs(field_p->dof_vert[nvert].vec_dof_1[i])<Accuracy){
	    double temp=0.0;
	    fprintf(fp,"%.12lg ",temp);
	  }
	  else{
            fprintf(fp,"%.12lg ",field_p->dof_vert[nvert].vec_dof_1[i]);
          }
          /*
	    utr_fprintf_double(fp, Accuracy, field_p->dof_vert[nvert].vec_dof_1[i]);
	    fprintf(fp," ");
          */     
	}
	fprintf(fp,"\n");
      }

      if(field_p->nr_sol>1 && selected[1]==1){

	for(i=0;i<nreq;i++){

	  if(fabs(field_p->dof_vert[nvert].vec_dof_2[i])<Accuracy){
	    double temp=0.0;
	    fprintf(fp,"%.12lg ",temp);
	  }
	  else{
            fprintf(fp,"%.12lg ",field_p->dof_vert[nvert].vec_dof_2[i]);
          }
	  /*
	    utr_fprintf_double(fp, Accuracy, field_p->dof_vert[nvert].vec_dof_2[i]);
	    fprintf(fp," ");
          */
	}
	fprintf(fp,"\n");
      }
      if(field_p->nr_sol>2 && selected[2]==1){

	for(i=0;i<nreq;i++){
	  
	  if(fabs(field_p->dof_vert[nvert].vec_dof_3[i])<Accuracy){
	    double temp=0.0;
	    fprintf(fp,"%.12lg ",temp);
	  }
	  else{
            fprintf(fp,"%.12lg ",field_p->dof_vert[nvert].vec_dof_3[i]);
          }
          /*
	    utr_fprintf_double(fp, Accuracy, field_p->dof_vert[nvert].vec_dof_3[i]);
	    fprintf(fp," ");
          */         
	}
	fprintf(fp,"\n");
      }
    } /* end if vertex active */

    /* for active constrained vertices */
    if(mmr_node_status(mesh_id,nvert)==MMC_ACTIVE && field_p->dof_vert[nvert].constr != NULL){
      fprintf(fp,"%d ",(-nvert));
      fprintf(fp,"%d ",field_p->dof_vert[nvert].constr[0]);
      for(i=1; i<=field_p->dof_vert[nvert].constr[0]; i++) {
	fprintf(fp,"%d %d ",field_p->dof_vert[nvert].constr[i],field_p->dof_vert[nvert].constr_type[i]);
      }
      fprintf(fp,"\n");
    }
    
  } /* end loop over all dofs vertexes structures */
  
  /* Save data for active edges vertices */
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    for(nedge=1;nedge<=mxedge;nedge++){
    
    /* for active edges */
      if((mmr_edge_status(mesh_id,nedge)==MMC_ACTIVE && field_p->dof_edge[nedge].vec_dof_1 != NULL) || 
	 (mmr_edge_status(mesh_id,nedge)==MMC_INACTIVE && field_p->dof_edge[nedge].vec_dof_1 != NULL)) {
       
	fprintf(fp,"%d ",nedge);

	if(selected[0]==1){
	  for(i=0;i<nreq;i++){
	  
	    if(fabs(field_p->dof_edge[nedge].vec_dof_1[i])<Accuracy){
	      double temp=0.0;
	      fprintf(fp,"%.12lg ",temp);
	    }
	    else{
	      fprintf(fp,"%.12lg ",field_p->dof_edge[nedge].vec_dof_1[i]);
	    }
	    /*
	      utr_fprintf_double(fp, Accuracy, field_p->dof_edge[nedge].vec_dof_1[i]);
	      fprintf(fp," ");
	    */
	}
	  fprintf(fp,"\n");
	}

	if(field_p->nr_sol>1 && selected[1]==1){
	  
	  for(i=0;i<nreq;i++){

	    if(fabs(field_p->dof_edge[nedge].vec_dof_2[i])<Accuracy){
	      double temp=0.0;
	      fprintf(fp,"%.12lg ",temp);
	    }
	    else{
	      fprintf(fp,"%.12lg ",field_p->dof_edge[nedge].vec_dof_2[i]);
	    }
	    /*
	      utr_fprintf_double(fp, Accuracy, field_p->dof_edge[nedge].vec_dof_2[i]);
	      fprintf(fp," ");
	    */
	  }
	  fprintf(fp,"\n");
	}
	if(field_p->nr_sol>2 && selected[2]==1){
	  
	  for(i=0;i<nreq;i++){
	    
	    if(fabs(field_p->dof_edge[nedge].vec_dof_3[i])<Accuracy){
	      double temp=0.0;
	      fprintf(fp,"%.12lg ",temp);
	    }
	    else{
	      fprintf(fp,"%.12lg ",field_p->dof_edge[nedge].vec_dof_3[i]);
	    }
	    /*
	      utr_fprintf_double(fp, Accuracy, field_p->dof_edge[nedge].vec_dof_3[i]);
	      fprintf(fp," ");
	    */         
	  }
	  fprintf(fp,"\n");
	}
      } /* end if edge active */

      /* for active constrained edges */
      if(mmr_edge_status(mesh_id,nedge)==MMC_ACTIVE && field_p->dof_edge[nedge].constr != NULL){
	fprintf(fp,"%d ",(-nedge));
	fprintf(fp,"%d ",field_p->dof_edge[nedge].constr[0]);
	for(i=1; i<=field_p->dof_edge[nedge].constr[0]; i++) {
	  fprintf(fp,"%d %d ",field_p->dof_edge[nedge].constr[i],field_p->dof_edge[nedge].constr_type[i]);
	}
	fprintf(fp,"\n");
      }
    } /* end loop over all dofs edges structures */
  }

  /* close file with control data */
  fclose(fp);

  return(1);
}


/*------------------------------------------------------------
  apr_check_field - to free approximation field data structure
------------------------------------------------------------*/
int apr_check_field(
  int Field_id    /* in: approximation field ID  */
  )
{
  /* pointer to field structure */
  apt_field *field_p;
  int mesh_id;
  int nvert, mxvert, nr_dof_vert, nr_constr_dof_vert;
  int nedge, mxedge, nr_dof_edge, nr_constr_dof_edge;
  int nr_dof_ents;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);

  /* check field PDEG */
  apr_set_and_test_pdeg(Field_id,NULL);

#ifdef DEBUG_APM
  printf("\n\t(apr_check_field) -> Approximation type: ");
  switch(field_p->pdeg)
    {
    case APC_LINEAR_APPROXIMATION_PDEG:
      printf("LINEAR\n"); break;
    case APC_QUADRATIC_APPROXIMATION_PDEG:
      printf("QUADRATIC\n"); break;
    default:
      printf("-\n"); break;
    };
#endif

  /* select the corresponding mesh */
  mesh_id = apr_get_mesh_id(Field_id);

  /* loop over dofs data structure including free spaces */
  mxvert=mmr_get_max_node_id(mesh_id);
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    mxedge=mmr_get_max_edge_id(mesh_id);
  }

  /* Loop over dof vertexes */
  assert(mxvert >= field_p->nr_dof_vert);
  for(nvert=1;nvert<=mxvert;nvert++){	
    if(mmr_node_status(mesh_id, nvert)==MMC_FREE){
      if(field_p->dof_vert[nvert].vec_dof_1 != NULL){
	printf("[\e[1;31mERROR\e[0;0m] Dofs VERTEX associated with free space %d in check_field, exiting\n",nvert);
	exit(-1);
      }
      if(field_p->dof_vert[nvert].constr != NULL){
	printf("[\e[1;31mERROR\e[0;0m] Constr VERTEX associated with free space %d in check_field, exiting\n",nvert);
	exit(-1);
      }
      if(apr_get_ent_pdeg(Field_id, APC_VERTEX, nvert)!=-1){
	printf("[\e[1;31mERROR\e[0;0m] Error 3843 in approximation data structures for %d. Exiting!\n",nvert);
	assert(!"Holding execution before exit!");
	exit(-1);
      }
    }
    else{
      if(apr_get_ent_pdeg(Field_id, APC_VERTEX, nvert)==1){
	if(field_p->dof_vert[nvert].vec_dof_1 == NULL){
	  printf("[\e[1;32mERROR\e[0;0m] Error 38 in approximation data structures for %d. Exiting!\n",nvert);
	  exit(-1);
	}

	if(field_p->dof_vert[nvert].constr != NULL){
	  printf("[\e[1;31mERROR\e[0;0m] Constr VERTEX associated with real node %d in check_field, exiting\n",nvert);
	  exit(-1);
	}
      }
      else if(apr_get_ent_pdeg(Field_id, APC_VERTEX, nvert)==0){		
	if(field_p->dof_vert[nvert].constr == NULL){
	  assert(field_p->dof_vert[nvert].constr != NULL);
	  printf("[\e[1;31mERROR\e[0;0m] Error 43 in approximation data structures for %d in check field. Exiting!\n",nvert);
	  exit(-1);
	}

	if(field_p->dof_vert[nvert].vec_dof_1 != NULL){
	  printf("[\e[1;31mERROR\e[0;0m] Dofs VERTEX associated with constr node %d in check_field, exiting\n",nvert);
	  exit(-1);
	}
      }
    }
  }

  /* Loop over dof edges */
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    assert(mxedge >= field_p->nr_dof_edge);
    for(nedge=1;nedge<=mxedge;nedge++){	
      if(mmr_edge_status(mesh_id, nedge)==MMC_FREE){
	if(field_p->dof_edge[nedge].vec_dof_1 != NULL){
	  printf("[\e[1;31mERROR\e[0;0m] Dofs EDGE associated with free space %d in check_field, exiting\n",nedge);
	  exit(-1);
	}
	if(field_p->dof_edge[nedge].constr != NULL){
	  printf("[\e[1;31mERROR\e[0;0m] Constr EDGE associated with free space %d in check_field, exiting\n",nedge);
	  exit(-1);
	}
	if(apr_get_ent_pdeg(Field_id, APC_EDGE, nedge)!=-1){
	  printf("[\e[1;31mERROR\e[0;0m] Error 3944 in approximation data structures for %d. Exiting!\n",nedge);
	  assert(!"Holding execution before exit!");
	  exit(-1);
	}
      }
      else{
	
	if(apr_get_ent_pdeg(Field_id, APC_EDGE, nedge)>0){
	  
	  if(field_p->dof_edge[nedge].vec_dof_1 == NULL){
	    printf("[\e[1;31mERROR\e[0;0m] Error 39 in approximation data structures for %d. Exiting!\n", nedge);
	    exit(-1);
	  }
	  
	  if(field_p->dof_edge[nedge].constr != NULL){
	    printf("[\e[1;31mERROR\e[0;0m] Constr EDGE associated with real node %d in check_field, exiting\n", nedge);
	    exit(-1);
	  }
	}
	else if(apr_get_ent_pdeg(Field_id, APC_EDGE, nedge)==0){		
	  if(field_p->dof_edge[nedge].constr == NULL){
	    assert(field_p->dof_edge[nedge].constr != NULL);
	    printf("[\e[1;31mERROR\e[0;0m] Error 44 in approximation data structures for %d in check field. Exiting!\n", nedge);
	    exit(-1);
	  }

	  if(field_p->dof_edge[nedge].vec_dof_1 != NULL){
	    printf("[\e[1;31mERROR\e[0;0m] Dofs EDGE associated with constr node %d in check_field, exiting\n", nedge);
	    exit(-1);
	  }
	}
      }
    }
  }

  /* Check number of dofs (edges+vertex) */
  nr_dof_vert=0; nr_constr_dof_vert=0;
  nr_dof_edge=0; nr_constr_dof_edge=0;
  nr_dof_ents=0;
  
  /* Loop over vertexes */
  for(nvert=1;nvert<=mxvert;nvert++){
    /* for active unconstrained vertices */
    if(mmr_node_status(mesh_id,nvert)==MMC_ACTIVE && field_p->dof_vert[nvert].vec_dof_1 != NULL) { 
      nr_dof_ents++; 
      nr_dof_vert++; 
    }

    /* for active constrained vertices */
    if(mmr_node_status(mesh_id,nvert)==MMC_ACTIVE && field_p->dof_vert[nvert].constr != NULL) {
      nr_constr_dof_vert++;
    }
  }
  
  /* Loop over edges */
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    for(nedge=1;nedge<=mxedge;nedge++){
      /* for active unconstrained edges */
      if((mmr_edge_status(mesh_id,nedge)==MMC_ACTIVE || mmr_edge_status(mesh_id,nedge)==MMC_INACTIVE) && field_p->dof_edge[nedge].vec_dof_1 != NULL) { 
	nr_dof_ents++; 
	nr_dof_edge++; 
      }

      /* for active constrained edges */
      if((mmr_edge_status(mesh_id,nedge)==MMC_ACTIVE || mmr_edge_status(mesh_id,nedge)==MMC_INACTIVE) && field_p->dof_edge[nedge].constr != NULL) { 
	nr_constr_dof_edge++; 
      }
    }
  }

  if(nr_dof_ents!=field_p->nr_dof_ents ||	\
     nr_dof_vert!=field_p->nr_dof_vert ||	\
     nr_dof_edge!=field_p->nr_dof_edge) {    
    printf("[\e[1;31mERROR\e[0;0m] Wrong number of dof_entities in apr_check_field! %d %d Exiting!\n",nr_dof_ents,field_p->nr_dof_ents);
    printf("[\e[1;31mERROR\e[0;0m] Number of vertexes %d %d\n",nr_dof_vert,field_p->nr_dof_vert);
    if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
      printf("[\e[1;31mERROR\e[0;0m] Number of edges %d %d\n",nr_dof_edge,field_p->nr_dof_edge);
    }
    exit(-1);
  }
  if(nr_constr_dof_vert!=field_p->nr_constr_dof_vert ||
     nr_constr_dof_edge!=field_p->nr_constr_dof_edge) {
    printf("[\e[1;31mERROR\e[0;0m] Wrong number of dof_entities in apr_check_field! %d %d Exiting!\n",nr_dof_ents,field_p->nr_dof_ents);
    printf("[\e[1;31mERROR\e[0;0m] Number of constrained vertexes %d %d\n",nr_constr_dof_vert,field_p->nr_constr_dof_vert);
    if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
      printf("[\e[1;31mERROR\e[0;0m] Number of constrained edges %d %d\n",nr_constr_dof_edge,field_p->nr_constr_dof_edge);
    }
    exit(-1);
  }

  return(1);
}

/*------------------------------------------------------------
  apr_get_mesh_id - to return the ID of the associated mesh
------------------------------------------------------------*/
int apr_get_mesh_id( /* returns: >0 - ID of the associated mesh,
                                        <0 - error code */
  int Field_id     /* in: approximation field ID  */
  )
{
  /* pointer to field structure */
  apt_field *field_p;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);

  return(field_p->mesh_id);
}

/*------------------------------------------------------------
  apr_get_nreq - to return the number of components in solution vector
------------------------------------------------------------*/
int apr_get_nreq( /* returns: >0 - number of solution components,
                                        <0 - error code */
  int Field_id     /* in: approximation field ID  */
  )
{
  /* pointer to field structure */
  apt_field *field_p;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);

  return(field_p->nreq);
}

/*------------------------------------------------------------
  apr_get_nr_sol - to return the number of solution vectors stored
------------------------------------------------------------*/
int apr_get_nr_sol( /* returns: >0 - number of solution vectors,
                                        <0 - error code */
  int Field_id     /* in: approximation field ID  */
  )
{
  /* pointer to field structure */
  apt_field *field_p;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);

  return(field_p->nr_sol);
}



/*---------------------------------------------------------
  apr_get_base_type - to return the type of basis functions
  REMARK: type of basis functions differentiates element types as well 
  examples for standard linear approximation (from include/aph_intf.h): 
  #define APC_BASE_PRISM_STD  3   // for linear prismatic elements 
  #define APC_BASE_TETRA_STD  4  // for linear tetrahedral elements 
---------------------------------------------------------*/
int apr_get_base_type(/* returns: >0 - type of basis functions,
                                  <0 - error code */
  int Field_id,  /* in: field ID */
  int El_id      /* in: element ID */
  )
{

  /* pointer to field structure */
  apt_field *field_p;
  int mesh_id, el_type;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);

  /* get the corresponding mesh ID */
  mesh_id = field_p->mesh_id;

// base type is used to differentiate element types
// get element type
  el_type = mmr_el_type(mesh_id, El_id);
  
  if(el_type == MMC_PRISM) {
    return(APC_BASE_PRISM_STD);
  } else
  if(el_type == MMC_TETRA) {
    return(APC_BASE_TETRA_STD);
  } else {
    printf("Wrong element type in apr_get_base_type. Exiting!\n");
    exit(-1);
  }
}


/*------------------------------------------------------------
  apr_get_ent_pdeg - to return the degree of approximation index
                      associated with a given mesh entity
------------------------------------------------------------*/
int apr_get_ent_pdeg( /* returns: >0 - approximation index,
                                   0 - dof entity inactive (constrained)
                                  <0 - error code */
  int Field_id,     /* in: approximation field ID  */
  int Ent_type,      /* in: type of mesh entity */
  int Ent_id         /* in: mesh entity ID */
  )
{
  /* pointer to field structure */
  apt_field *field_p;
  int mesh_id;

/* auxiliary variables */
  int pdeg;

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef DEBUG_APM
  /* check input for standard approximation */
  if(Ent_type!=APC_VERTEX && Ent_type!=APC_EDGE){
    assert(!"Holding execution.");
    printf("Wrong dof entity type in get_ent_pdeg !\n");
    exit(-1);
  }
#endif

  field_p = apr_select_field(Field_id); //get the field ID
  mesh_id = field_p->mesh_id; //get the corresponding mesh ID

  switch(Ent_type) {
    case APC_VERTEX:

      if(mmr_node_status(mesh_id,Ent_id) == MMC_FREE || mmr_node_status(mesh_id,Ent_id) == MMC_INACTIVE || field_p->dof_vert[Ent_id].active == 0) {
	pdeg = -1;
      } else {
	if(field_p->dof_vert[Ent_id].vec_dof_1==NULL) { pdeg = 0; }
	else { pdeg = 1; }
      } 
      break;
      
    case APC_EDGE:

      if(mmr_edge_status(mesh_id,Ent_id)==MMC_FREE) {
	pdeg = -1; //-2
      } else if(mmr_edge_status(mesh_id,Ent_id)==MMC_INACTIVE) {
	if(field_p->dof_edge[Ent_id].active == 1) {
	  pdeg = 1; //2
	} else {
	  pdeg = -1; //-2
	}
      } else {
	if(field_p->dof_edge[Ent_id].vec_dof_1==NULL) { pdeg = 0; }
	else { pdeg = 1; } //2
      }
      break;
      
    default:
      exit(-1);
  };

  return(pdeg);
}


/*------------------------------------------------------------
  apr_set_ent_pdeg - to set the degree of approximation index
                      associated with a given mesh entity
------------------------------------------------------------*/
int apr_set_ent_pdeg( /* returns: >0 - success code,
                                          <0 - error code */
  int Field_id,     /* in: approximation field ID  */
  int Ent_type,      /* in: type of mesh entity */
  int Ent_id,         /* in: mesh entity ID */
  int Pdeg          /* in: degree of approximation */
  )
{
  /* pointer to field structure */
//  apt_field *field_p;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  return(1);

}

/*------------------------------------------------------------
  apr_get_ent_numshap - to return the number of shape functions (vector
                        dofs) associated with a given mesh entity
------------------------------------------------------------*/
int apr_get_ent_numshap( /* returns: >0 - the number of shape functions,
                                          <0 - error code */
  int Field_id,     /* in: approximation field ID  */
  int Ent_type,      /* in: type of mesh entity */
  int Ent_id         /* in: mesh entity ID */
  )
{
  /* auxiliary variables */
//  int nrdof;

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef DEBUG_APM
  /* check input for standard quadratic approximation */
  if(abs(Ent_type!=APC_VERTEX) && abs(Ent_type!=APC_EDGE)){
    printf("Wrong dof entity type in get_ent_pdeg for STD quadratic approximation !\n");
    assert(!"Holding execution!");
    exit(-1);
  }
#endif

  return(1);
}


/*------------------------------------------------------------
  apr_get_ent_nrdofs - to return the number of dofs associated with
                      a given mesh entity
------------------------------------------------------------*/
int apr_get_ent_nrdofs( /* returns: >0 - the number of dofs,
                                          <0 - error code */
  int Field_id,     /* in: approximation field ID  */
  int Ent_type,      /* in: type of mesh entity */
  int Ent_id         /* in: mesh entity ID */
  )
{
  /* auxiliary variables */
//  int nrdof;

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef DEBUG_APM
  /* check input for standard quadratic approximation */
  if(abs(Ent_type!=APC_VERTEX) && abs(Ent_type!=APC_EDGE)){
    printf("Wrong dof entity type in get_ent_pdeg for STD quadratic approximation !\n");
    assert(!"Holding execution!");
    exit(-1);
  }
#endif

  return(apr_get_nreq(Field_id));
}

/*------------------------------------------------------------
  apr_get_el_pdeg - to return the degree of approximation vector 
                      associated with a given element
------------------------------------------------------------*/
int apr_get_el_pdeg( /* returns:  >0 - success code - scalar pdeg
                                         <0 - error code */
  int Field_id,     /* in: approximation field ID  */
  int Ent_id,         /* in: element ID */
  int *Pdeg_vec       /* out: degree of approximation symbol or vector */
  )
{

  /* pointer to field structure */
  apt_field *field_p;
  int mesh_id;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  field_p = apr_select_field(Field_id);

  if(Pdeg_vec != NULL) *Pdeg_vec = field_p->pdeg;
    
  return(field_p->pdeg);

}


/*------------------------------------------------------------
  apr_set_el_pdeg - to set the degree of approximation vector 
                      associated with a given element
------------------------------------------------------------*/
int apr_set_el_pdeg( /* returns: >0 - success code,
                                          <0 - error code */
  int Field_id,     /* in: approximation field ID  */
  int Ent_id,         /* in: element ID */
  int *Pdeg_vec       /* in: degree of approximation symbol or vector */
  )
{

/*++++++++++++++++ executable statements ++++++++++++++++*/

  return(1);

}

/*---------------------------------------------------------
  apr_get_el_pdeg_numshap - to return the number of shape functions
                            (scalar DOFs) for an element given its 
                            degree of approximation symbol or vector pdeg
---------------------------------------------------------*/
int apr_get_el_pdeg_numshap(
		 /* returns: >=0 - success code, <0 - error code*/
  int Field_id,  /* in: field ID */
  int El_id,     /* in: element ID */
  int *Pdeg      /* in: degree of approximation symbol or vector */
  )
{

  int base = apr_get_base_type(Field_id, El_id);
  apt_field *field_p = apr_select_field(Field_id);

  /* REMARK: type of basis functions differentiates element types as well 
       examples for standard linear/quadratic approximation (from include/aph_intf.h): 
          #define APC_BASE_PRISM_STD  3   // for linear prismatic elements 
          #define APC_BASE_TETRA_STD  4  // for linear tetrahedral elements 
  */

  switch(base) {
  case APC_BASE_PRISM_STD:
    switch(field_p->pdeg) {
    case APC_LINEAR_APPROXIMATION_PDEG:
      return(6);
    case APC_QUADRATIC_APPROXIMATION_PDEG:
      return(15);
    }
   
  case APC_BASE_TETRA_STD:
    switch(field_p->pdeg) {
    case APC_LINEAR_APPROXIMATION_PDEG:
      return(4);
    case APC_QUADRATIC_APPROXIMATION_PDEG:
      return(10);
    }
      
  default:
    printf("Wrong element type in apr_get_base_type. Exiting!\n");
    exit(-1);
  }
}

/*---------------------------------------------------------
  apr_get_el_dofs - to return the number and the list of element's degrees
                        of freedom (corresponding to standard shape functions)
!!! The order is the same as for mmr_el_node_coor (loop over ieq inside):
(node1_1, node1_2, ..., node1_NREQ, node2_1, node2_2, ..., node2_NREQ, etc. !!!
---------------------------------------------------------*/
int apr_get_el_dofs( /* returns: >0 - the number of dofs
                                         <0 - error code */
  int Field_id,     /* in: approximation field ID  */
  int El_id,        /* in: element ID */
  int Vect_id,       /* in: vector ID in case of multiple solution vectors */
  double *El_dofs_std     /* out: the list of values of element dofs */
  /* the size of the table MUST BE >= APC_MAXELSD */
  )
{
  /* pointer to field structure */
  apt_field *field_p;
  int mesh_id, nreq, idof, idofent, ieq;
  int icon, ive, ied, ido, nr_dofs;
  int end_vert, end_edge;

  int *el_nodes = NULL; // list of nodes of El 
  int *nr_constr = NULL; // list of numbers of constraining nodes
  int *constr_id = NULL; // list with IDs of constraining nodes
  int *constr_type = NULL; // list with TYPEs of constraining nodes
  double *constr_value = NULL; // list with constraint coefficients
    
  double dofs_loc[APC_MAXEQ];

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);

  /* check PDEG value */
  apr_set_and_test_pdeg(Field_id,NULL);

  /* Alloc memory */
  if(field_p->pdeg == APC_LINEAR_APPROXIMATION_PDEG) {
    el_nodes = (int*) malloc ((MMC_MAXELVNO+1)*sizeof(int)); // [vertex + 1] (num/..vertices..)
    nr_constr = (int*) malloc ((MMC_MAXELVNO+1)*sizeof(int)); // [vertex + 1] (num/..nr_constr..)
    constr_id = (int*) malloc ((4*MMC_MAXELVNO)*sizeof(int)); // list with IDs of constraining nodes
    constr_type = (int*) malloc ((4*MMC_MAXELVNO)*sizeof(int)); // list with TYPEs of constraining nodes
    constr_value = (double*) malloc ((4*MMC_MAXELVNO)*sizeof(double)); // list with constraint coefficients
  }
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    el_nodes = (int*) malloc ((MMC_MAXELVNO+APC_MAXELVNO_ON_EDGES+2)*sizeof(int)); // [vertex + edges + 2] (num/..vertices../num/..edges..)
    nr_constr = (int*) malloc ((MMC_MAXELVNO+APC_MAXELVNO_ON_EDGES+2)*sizeof(int)); // [vertex + edges + 2] (num/..nr_constr../num/..nr_constr..)
    constr_id = (int*) malloc ((8*MMC_MAXELVNO+8*APC_MAXELVNO_ON_EDGES)*sizeof(int));  // list with IDs of constraining nodes [vertex + edges]
    constr_type = (int*) malloc ((8*MMC_MAXELVNO+8*APC_MAXELVNO_ON_EDGES)*sizeof(int)); // list with TYPEs of constraining nodes [vertex + edges]
    constr_value = (double*) malloc ((8*MMC_MAXELVNO+8*APC_MAXELVNO_ON_EDGES)*sizeof(double)); // list with constraint coefficients [vertex + edges]
  }

  /* get mesh ID */
  mesh_id = apr_get_mesh_id(Field_id);

  nreq = apr_get_nreq(Field_id);

  for(idof=0;idof<APC_MAXELSD;idof++) El_dofs_std[idof]=0.0;

  /* if there are no constraints in the mesh */
  if(field_p->constr==APC_FALSE){
    
    /* get vertices */
    ive = 0; //Number of VERTICES [0] / First value [1] 
    mmr_el_node_coor(mesh_id,El_id,&el_nodes[ive],NULL);
    
    /* get edges */
    if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
      ied = el_nodes[ive]+1; //Number of EDGES [number_of_vertices+1] / First value [number_of_vertices+2] 
      mmr_el_edges(mesh_id,El_id,&el_nodes[ied]);
    }
    
    /* simply rewrite vertex dofs */
    for(ido=0; ido<el_nodes[ive]; ido++) {
      apr_read_ent_dofs(Field_id,APC_VERTEX,el_nodes[ido+1],nreq, Vect_id, &El_dofs_std[ido*nreq]);
    }

    if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
      for(ido=el_nodes[ive]; ido<(el_nodes[ive]+el_nodes[ied]); ido++) {
	apr_read_ent_dofs(Field_id,APC_EDGE,el_nodes[ido+2],nreq, Vect_id, &El_dofs_std[ido*nreq]);
      }
    }
  }
  else{
    /* if there are constraints */

    /* get the list of real nodes with constraint coefficients */
    apr_get_el_constr_data(Field_id,El_id,el_nodes,nr_constr,constr_id,constr_type,constr_value);

    ive = 1; end_vert = el_nodes[0]+ive;

    if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
      ied = end_vert+1; end_edge = ied+el_nodes[end_vert];
    }
    
    /*jbw
    idofent=0;
    printf("(apr_get_el_dofs - constraints data print\n");
    printf("Num vertices: %d\n",el_nodes[ive-1]);
    for(ido=ive; ido<end_vert; ido++) {
      printf("Constraint VERTEX: %d; number of constraint: %d\n",el_nodes[ido],nr_constr[ido]);
      for(icon=0; icon<nr_constr[ido]; icon++,idofent++) {
	printf("\tConstr ID: %d;\t",constr_id[idofent]);
	if(constr_type[idofent] == APC_VERTEX) {
	  printf("type: APC_VERTEX;");
	}
	if(constr_type[idofent] == APC_EDGE) {
	  printf("type: APC_EDGE;");
	}
	//printf("\tvalue: %.2lf\n",constr_value[idofent]);
	printf("\n");
      }
    }
    if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
      printf("Num edges: %d\n",el_nodes[ied-1]);
      for(ido=ied; ido<end_edge; ido++) {
	printf("Constraint EDGE: %d; number of constraint: %d\n",el_nodes[ido],nr_constr[ido]);
      
	for(icon=0; icon<nr_constr[ido]; icon++,idofent++) {
	  printf("\tConstr ID: %d;\t",constr_id[idofent]);
	  if(constr_type[idofent] == APC_VERTEX) {
	    printf("type: APC_VERTEX;\t");
	  }
	  if(constr_type[idofent] == APC_EDGE) {
	    printf("type: APC_EDGE;\n");
	  }
	  //printf("\tvalue: %.2lf\n",constr_value[idofent]);
	  printf("\n");
	}
      }
    }
    /*jbw*/    

    idofent=0;
    for(ido=ive; ido<end_vert; ido++) {
      for(icon=0; icon<nr_constr[ido]; icon++,idofent++) {
	apr_read_ent_dofs(Field_id,constr_type[idofent],constr_id[idofent],nreq,Vect_id,dofs_loc);
	for(ieq=0;ieq<nreq;ieq++) {
	  El_dofs_std[(ido-1)*nreq+ieq] += dofs_loc[ieq]*constr_value[idofent];

#ifdef DEBUG_APM
	  //printf("[\e[1;30mDEBUG INFO\e[0;0m] (line: %d) IDO: %d\tdof[%4d] coeff: %lf\tvalue: %lf\tEL_DOFS_STD[%4d] = %lf\n",__LINE__,(ido-1),constr_id[idofent],constr_value[idofent],dofs_loc[ieq],((ido-1)*nreq+ieq),El_dofs_std[(ido-1)*nreq+ieq]);
#endif
	}
      }
    }
    if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
      for(ido=ied; ido<end_edge; ido++) {
	for(icon=0; icon<nr_constr[ido]; icon++,idofent++) {
	  apr_read_ent_dofs(Field_id,constr_type[idofent],constr_id[idofent],nreq,Vect_id,dofs_loc);
	  for(ieq=0;ieq<nreq;ieq++) {
	    El_dofs_std[(ido-2)*nreq+ieq] += dofs_loc[ieq]*constr_value[idofent];

#ifdef DEBUG_APM
	    //printf("[\e[1;30mDEBUG INFO\e[0;0m] (line: %d) IDO: %d\tdof[%4d] coeff: %lf\tvalue: %lf\tEL_DOFS_STD[%4d] = %lf\n",__LINE__,(ido-2),constr_id[idofent],constr_value[idofent],dofs_loc[ieq],((ido-2)*nreq+ieq),El_dofs_std[(ido-2)*nreq+ieq]);
#endif	    
	  }
	}
      }
    }
  }

  if(field_p->pdeg == APC_LINEAR_APPROXIMATION_PDEG) {
    nr_dofs = nreq*(el_nodes[ive]);
  }
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    nr_dofs = nreq*(el_nodes[ive]+el_nodes[ied]);
  }

  assert(El_dofs_std == NULL || (El_dofs_std[0] > -10e200 && El_dofs_std[0] < 10e200));

  SAFE_FREE(el_nodes);
  SAFE_FREE(nr_constr);
  SAFE_FREE(constr_id);
  SAFE_FREE(constr_type);
  SAFE_FREE(constr_value);
  
  return(nr_dofs);
} 


/*---------------------------------------------------------
  apr_get_nrdofs_glob - to return a global dimension of the problem
---------------------------------------------------------*/
int apr_get_nrdofs_glob(	/* returns: global dimension of the problem */
  int Field_id    /* in: field ID */
	)
{

  /* pointer to field structure */
  apt_field *field_p;
  int mesh_id;

/* auxiliary variables */
  int glob_dim=0, nno; 

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);

  /* check PDEG value */
  apr_set_and_test_pdeg(Field_id,NULL);

#ifdef DEBUG_APM
  /* get the corresponding mesh ID */
  mesh_id = field_p->mesh_id;

  /* loop over vertices-nodes */
  nno=0;
  while((nno=mmr_get_next_node_all(mesh_id,nno))!=0) {
    if(mmr_node_status(mesh_id,nno)==MMC_ACTIVE && field_p->dof_vert[nno].vec_dof_1!=NULL) {
      glob_dim ++;
    }
  }
  
  /* loop over points on edges */
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    nno=0;
    while((nno=mmr_get_next_edge_all(mesh_id,nno))!=0) {
      if(mmr_edge_status(mesh_id,nno)==MMC_ACTIVE && field_p->dof_edge[nno].vec_dof_1!=NULL) {
	glob_dim ++;
      }
    }
  }

  if(glob_dim!=field_p->nr_dof_ents){
    printf("[\e[1;31mERROR\e[0;0m] Wrong number of dof_entities in get_nrdofs_glob! %d %d Exiting!\n",glob_dim,field_p->nr_dof_ents);
    exit(-1);
  }
#endif

  return(field_p->nr_dof_ents*apr_get_nreq(Field_id));
}

/*------------------------------------------------------------
  apr_read_ent_dofs - to read a vector of dofs associated with a given
                  mesh entity from approximation field data structure
------------------------------------------------------------*/
int apr_read_ent_dofs(/* returns: >=0 - success code, <0 - error code */
  int Field_id,     /* in: approximation field ID  */
  int Ent_type,      /* in: type of mesh entity */
  int Ent_id,        /* in: mesh entity ID */
  int Ent_nrdof,     /* in: number of dofs associated with the entity */
  int Vect_id,       /* in: vector ID in case of multiple solution vectors */
  double* Vect_dofs  /* out: dofs read from data structure */
  )
{
  /* pointer to field structure */
  apt_field *field_p;

  /* auxiliary variables */
  double* glob_dofs;
  int i; //mesh_id;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);

  /* get the corresponding mesh ID */
  //mesh_id = field_p->mesh_id;

  /* Select dofs type */
  switch(Ent_type) {
    case APC_VERTEX:
      #ifdef DEBUG_APM
      {
	/* check input */
	if(field_p->dof_vert[Ent_id].vec_dof_1 == NULL){
	  if(field_p->dof_vert[Ent_id].constr != NULL)
	    printf("[\e[1;31mERROR\e[0;0m] Wrong dof VERTEX entity \e[1;37m%d\e[0m in read_ent_dofs it is \e[1;33mconstrained\e[0m !\n",Ent_id);
	  else
	    printf("[\e[1;31mERROR\e[0;0m] Wrong dof VERTEX entity \e[1;37m%d\e[0m in read_ent_dofs it is \e[1;37mfree space\e[0m !\n",Ent_id);
	  getchar();getchar();
	  exit(-1);
	}
      }
      #endif

      if(Vect_id<=1) glob_dofs = field_p->dof_vert[Ent_id].vec_dof_1;
      else if(Vect_id==2) glob_dofs = field_p->dof_vert[Ent_id].vec_dof_2;
      else if(Vect_id==3) glob_dofs = field_p->dof_vert[Ent_id].vec_dof_3;
      else glob_dofs = field_p->dof_vert[Ent_id].vec_dof_1;

      assert(Vect_dofs != NULL);

      #ifdef DEBUG_APM
	if(glob_dofs==NULL){	  
	  printf("[\e[1;31mERROR\e[0;0m] Field %d, ent_type %d, ent_nrdof %d, vect_id %d\n",Field_id, Ent_type, Ent_nrdof, Vect_id);
	  printf("\t\tdof VERTEX entity %d in read_ent_dofs !\n",Ent_id);
	  getchar();getchar();
	  exit(-1);
	}
      #endif
      assert(glob_dofs != NULL);

      for(i=0;i<Ent_nrdof;i++)
	  Vect_dofs[i] = glob_dofs[i];
      
      break;
      
    case APC_EDGE:
      #ifdef DEBUG_APM
      {
	/* check input */
	if(field_p->dof_edge[Ent_id].vec_dof_1 == NULL) {
	  if(field_p->dof_edge[Ent_id].constr != NULL)
	    printf("[\e[1;31mERROR\e[0;0m] Wrong dof EDGE entity \e[1;37m%d\e[0m in read_ent_dofs it is \e[1;33mconstrained\e[0m !\n",Ent_id);
	  else
	    printf("[\e[1;31mERROR\e[0;0m] Wrong dof EDGE entity \e[1;37m%d\e[0m in read_ent_dofs it is \e[1;30mfree space\e[0m !\n",Ent_id);
	  getchar();getchar();
	  exit(-1);
	}
      }
      #endif

      if(Vect_id<=1) glob_dofs = field_p->dof_edge[Ent_id].vec_dof_1;
      else if(Vect_id==2) glob_dofs = field_p->dof_edge[Ent_id].vec_dof_2;
      else if(Vect_id==3) glob_dofs = field_p->dof_edge[Ent_id].vec_dof_3;
      else glob_dofs = field_p->dof_edge[Ent_id].vec_dof_1;

      assert(Vect_dofs != NULL);

      #ifdef DEBUG_APM
	if(glob_dofs==NULL){
	  printf("[\e[1;31mERROR\e[0;0m] Field %d, ent_type %d, ent_nrdof %d, vect_id %d\n",Field_id, Ent_type, Ent_nrdof, Vect_id); 
	  printf("\t\tdof EDGE entity %d in read_ent_dofs !\n",Ent_id);
	  getchar();getchar();
	  exit(-1);
	}
      #endif
      assert(glob_dofs != NULL);

      for(i=0;i<Ent_nrdof;i++)
	  Vect_dofs[i] = glob_dofs[i];
      
      break;
    
    default:
      exit(-1);
  };

  assert(Vect_dofs == NULL || (Vect_dofs[0] > -10e200 && Vect_dofs[0] < 10e200));
  return(1);
}

/*------------------------------------------------------------
  apr_write_ent_dofs - to write a vector of dofs associated with a given
                   mesh entity to approximation field data structure
------------------------------------------------------------*/
int apr_write_ent_dofs(/* returns: >=0 - success code, <0 - error code */
  int Field_id,     /* in: approximation field ID  */
  int Ent_type,      /* in: type of mesh entity */
  int Ent_id,        /* in: mesh entity ID */
  int Ent_nrdof,     /* in: number of dofs associated with the entity */
  int Vect_id,       /* in: vector ID in case of multiple solution vectors */
  double* Vect_dofs  /* in: dofs to be written */
  )
{
  /* pointer to field structure */
  apt_field *field_p;

  /* auxiliary variables */
  double* glob_dofs;
  int i;//,j

/*++++++++++++++++ executable statements ++++++++++++++++*/

  //assert(Vect_dofs[0] > -10e200 && Vect_dofs[0] < 10e200);
  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);

  /* Select dofs type */
  switch(Ent_type) {
    case APC_VERTEX:
      #ifdef DEBUG_APM
      {
	/* check input */
	if(field_p->dof_vert[Ent_id].vec_dof_1==NULL){
	  printf("Wrong dof VERTEX entity in write_ent_dofs !\n");
	  exit(-1);
	}
      }
      #endif

      if(Vect_id<=1) glob_dofs = field_p->dof_vert[Ent_id].vec_dof_1;
      else if(Vect_id==2) glob_dofs = field_p->dof_vert[Ent_id].vec_dof_2;
      else if(Vect_id==3) glob_dofs = field_p->dof_vert[Ent_id].vec_dof_3;
      else glob_dofs = field_p->dof_vert[Ent_id].vec_dof_1;

      for(i=0;i<Ent_nrdof;i++)  glob_dofs[i] = Vect_dofs[i];
      
      break;
      
    case APC_EDGE:
      #ifdef DEBUG_APM
      {
	/* check input */
	if(field_p->dof_edge[Ent_id].vec_dof_1==NULL){
	  printf("Wrong dof EDGE entity in write_ent_dofs !\n");
	  exit(-1);
	}
      }
      #endif

      if(Vect_id<=1) glob_dofs = field_p->dof_edge[Ent_id].vec_dof_1;
      else if(Vect_id==2) glob_dofs = field_p->dof_edge[Ent_id].vec_dof_2;
      else if(Vect_id==3) glob_dofs = field_p->dof_edge[Ent_id].vec_dof_3;
      else glob_dofs = field_p->dof_vert[Ent_id].vec_dof_1;

      for(i=0;i<Ent_nrdof;i++)  glob_dofs[i] = Vect_dofs[i];
      break;
    
    default:
      exit(-1);
  };

  return(1);
}


/*------------------------------------------------------------
  apr_create_ent_dofs - to create a vector of dofs associated with a given
                   mesh entity within approximation field data structure
------------------------------------------------------------*/
int apr_create_ent_dofs(/* returns: >=0 - success code, <0 - error code */
  int Field_id,     /* in: approximation field ID  */
  int Ent_type,      /* in: type of mesh entity */
  int Ent_id,        /* in: mesh entity ID */
  int Ent_nrdof,     /* in: number of dofs associated with the entity */
  int Vect_id,       /* in: vector ID in case of multiple solution vectors */
  double* Vect_dofs  /* in: dofs to be written */
  )
{
  /* pointer to field structure */
  apt_field *field_p;

  /* auxiliary variables */
  //double* glob_dofs;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);

  /* Select dofs type */
  switch(Ent_type) {
    case APC_VERTEX:
      #ifdef DEBUG_APM
      {
	/* check input */
	if(Ent_nrdof != apr_get_ent_nrdofs(Field_id, Ent_type, Ent_id)){
	  printf("Wrong number of dofs VERTEX in write_ent_dofs !\n");
	  exit(-1);
	}
      }
      #endif

      if(Vect_id<=1) {
	field_p->dof_vert[Ent_id].vec_dof_1 = (double *) malloc(Ent_nrdof*sizeof(double));
	//glob_dofs = field_p->dof_vert[Ent_id].vec_dof_1;
      }
      else if(Vect_id==2) {
	field_p->dof_vert[Ent_id].vec_dof_2 = (double *) malloc(Ent_nrdof*sizeof(double));
	//glob_dofs = field_p->dof_vert[Ent_id].vec_dof_2;
      }
      else if(Vect_id==3) {
	field_p->dof_vert[Ent_id].vec_dof_3 = (double *) malloc(Ent_nrdof*sizeof(double));
	//glob_dofs = field_p->dof_vert[Ent_id].vec_dof_3;
      }
      else  {
	field_p->dof_vert[Ent_id].vec_dof_1 = (double *) malloc(Ent_nrdof*sizeof(double));
	//glob_dofs = field_p->dof_vert[Ent_id].vec_dof_1;
      }
      
      break;
      
    case APC_EDGE:
      #ifdef DEBUG_APM
      {
	/* check input */
	if(Ent_nrdof != apr_get_ent_nrdofs(Field_id, Ent_type, Ent_id)){
	  printf("Wrong number of dofs VERTEX in write_ent_dofs !\n");
	  exit(-1);
	}
      }
      #endif

      if(Vect_id<=1) {
	field_p->dof_edge[Ent_id].vec_dof_1 = (double *) malloc(Ent_nrdof*sizeof(double));
	//glob_dofs = field_p->dof_edge[Ent_id].vec_dof_1;
      }
      else if(Vect_id==2) {
	field_p->dof_edge[Ent_id].vec_dof_2 = (double *) malloc(Ent_nrdof*sizeof(double));
	//glob_dofs = field_p->dof_edge[Ent_id].vec_dof_2;
      }
      else if(Vect_id==3) {
	field_p->dof_edge[Ent_id].vec_dof_3 = (double *) malloc(Ent_nrdof*sizeof(double));
	//glob_dofs = field_p->dof_edge[Ent_id].vec_dof_3;
      }
      else  {
	field_p->dof_edge[Ent_id].vec_dof_1 = (double *) malloc(Ent_nrdof*sizeof(double));
	//glob_dofs = field_p->dof_edge[Ent_id].vec_dof_1;
      }
      
      break;
      
    default:
      exit(-1);
  };

  return(1);
}

/*------------------------------------------------------------
  apr_set_ini_con - to set an initial condition
------------------------------------------------------------*/
int apr_set_ini_con(/* returns: >=0 - success code, <0 - error code */
  int Field_id,     /* in: approximation field ID  */
  double (*Fun_p)(int, double*, int) /* pointer to function that provides */
		 /* problem dependent initial condition data */
  )
{
  double x[3];
  int sol_comp;
  double Uval[APC_MAXEQ];
  
  int nreq = apr_get_nreq(Field_id);
  int mesh_id = apr_get_mesh_id(Field_id);
  
  // in a loop over all vertices (nodes)
  int node_id=0;
  while ((node_id = mmr_get_next_node_all(mesh_id, node_id)) != 0) {
    if (apr_get_ent_pdeg(Field_id, APC_VERTEX, node_id) > 0) {
      
      // get node coordinates
      mmr_node_coor(mesh_id, node_id, x);   
      
      // for each field vector component
      for(sol_comp = 0; sol_comp < nreq; sol_comp++){
	
	// get initial condition value
	Uval[sol_comp] = (*Fun_p)(Field_id, x, sol_comp);
	
/*kbw
	debug("in set_ini_con - get value %lf from pd function\n",(*Fun_p)(Field_id, x, sol_comp));
	printf("in set_ini_con - get value %lf from pd function\n",(*Fun_p)(Field_id, x, sol_comp));
/*kew*/
	
	// set value as initial condition for a given node and solution component
	apr_write_ent_dofs(Field_id, APC_VERTEX, node_id, nreq, 0, Uval); 
	
      } // end loop over solution components
    } // end if active node
  } // end loop over vertices (nodes)
    
  return(1);
}


/*------------------------------------------------------------
  apr_prepare_integration_parameters - used e.g. by apr_num_int_el
------------------------------------------------------------*/
int apr_prepare_integration_parameters( 
  int Field_id,     /* in: approximation field ID  */
  int El_id,        /* in: element ID */
  int *Geo_order,   /* out: geometrical order of approximation */
  int *Num_geo_dofs,/* out: number of geometrical degrees of freedom */
  double* Geo_dofs, /* out: geometrical degrees of freedom */
  int *El_mate,     /* out: material index (in materials database) */
  int *Base,        /* out: type of basis functions */
  int *Pdeg_vec,    /* out: degree of approximation symbol or vector */
  int *Num_shap,    /* out: number of shape functions */
  int *Nreq,        /* out: number of components for unknowns */
  int *Num_dofs     /* out: number of scalar dofs (num_dofs = nreq*num_shap) */
)
{
  int mesh_id;
  int pdeg;		/* degree of polynomial */
  int el_nodes[MMC_MAXELVNO+1]; /* list of nodes of El */
  
/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the corresponding mesh */
  mesh_id = apr_get_mesh_id(Field_id);

  *Geo_order = 1; // geometrically (multi)linear elements
  
  /* get the coordinates of the nodes of El in the right order */
  mmr_el_node_coor(mesh_id,El_id,el_nodes,Geo_dofs);
  
  /* for geometrically (multi)linear elements number of degrees of freedom */
  /* is equal to the number of vertices - classical FEM nodes */
  *Num_geo_dofs = el_nodes[0];

#ifdef DEBUG_APM
  if(mmr_el_status(mesh_id,El_id)!=MMC_ACTIVE){
    printf("Asking for pdeg of inactive element in apr_num_int_el !\n");
    exit(-1);
  }
#endif

  *El_mate =  mmr_el_groupID(mesh_id, El_id);
  *Base = apr_get_base_type(Field_id, El_id);

  apr_get_el_pdeg(Field_id, El_id, &pdeg);
  Pdeg_vec[0]=pdeg;

  *Num_shap = apr_get_el_pdeg_numshap(Field_id,El_id,&pdeg);
  *Nreq=apr_get_nreq(Field_id);
  *Num_dofs = (*Nreq)*(*Num_shap);

  return(1);
}

/*------------------------------------------------------------
  apr_num_int_el - to perform numerical integration for an element
------------------------------------------------------------*/
int apr_num_int_el(
  int Problem_id,
  int Field_id,    /* in: approximation field ID  */
  int El_id,       /* in: unique identifier of the element */ 
  int Comp_sm,     /* in: indicator for the scope of computations: */
                   /*   APC_NO_COMP  - do not compute anything */
                   /*   APC_COMP_SM - compute entries to stiff matrix only */
                   /*   APC_COMP_RHS - compute entries to rhs vector only */
                   /*   APC_COMP_BOTH - compute entries for sm and rhsv */
  int *Pdeg_vec,        /* in: enforced degree of polynomial (if !=NULL ) */
  double *Sol_dofs_k,   /* in: solution dofs from previous iteration */
                        /*     (for nonlinear problems) */
  double *Sol_dofs_n,   /* in: solution dofs from previous time step */ 
			/*     (for nonlinear problems) */
  /* out: various matrices and vectors containing integrals that appear in */
  /*      the weak form */
  /* 1. matrices for the left hand side with the size num_dofs*num_dofs */
  /*      (num_dofs=nreq*num_shap) */
  int *Diagonal,   /* array of indicators whether matrices are diagonal */
  double *Stiff_mat,	/* out: stiffness matrix stored columnwise */
  double *Rhs_vect	/* out: rhs vector */
  // REMARK:
  //   1. udofs, ueq - correspond to the solution u and go from left to right
  //   2. wdofs, weq - correspond to test functions w and go up and down
  //   3. matrices are stored columnwise in vectors
  //   4. solution indices change in rows (go from left to right)
  //   5. test functions indices change in columns (go up and down)
  //   6. when matrices are stored in a vector the index is computed as:
  //      vector[udofs*nreq*num_dofs+wdofs*nreq+ueq*num_dofs+weq]
  //      (num_dofs=num_shap*nreq)
  //   7. for each pair (wdofs,udofs) there is a small submatrix nreq x nreq
  //      with indices ueq (from left to right) and weq (from top to bottom)
  //   8. when stiffness matrix entries are computed for vector problems with 
  //      the same shape functions for each component, then for each pair 
  //      (wdofs,udofs) a matrix of coefficients with the size nreq x nreq 
  //      is provided by the problem dependent module; its index is computed as 
  //      vector[ueq*nreq+weq]
  //   9. when load vector entries are computed for vector problems with 
  //      the same shape functions for each component, then for each wdofs index 
  //      a vector of coefficients with the size nreq is provided by the problem 
  //      dependent module
		   )
{

  /* pde coefficients */
  static int coeff_ind = 0;
  int coeff_vect_ind[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  // input to pdr_select_el_coeff_vect:
  // 0 - perform substitutions in coeff_vect_ind only if coeff_vect_ind[0]==0
  // output from pdr_select_el_coeff_vect:
  // 0 - coeff_vect_ind[0]==1 - always
  // 1 - mval, 2 - axx, 3 - axy, 4 - axz, 5 - ayx, 6 - ayy, 7 - ayz, 
  // 8 - azx, 9 - azy, 10 - azz, 11 - bx, 12 - by, 13 - bz
  // 14 - tx, 15 - ty, 16 - tz, 17 - cval
  // 18 - lval, 19 - qx, 20 - qy, 21 - qz, 22 - sval 

  double axx[APC_MAXEQ*APC_MAXEQ];
  double axy[APC_MAXEQ*APC_MAXEQ];
  double axz[APC_MAXEQ*APC_MAXEQ];
  double ayx[APC_MAXEQ*APC_MAXEQ];
  double ayy[APC_MAXEQ*APC_MAXEQ];
  double ayz[APC_MAXEQ*APC_MAXEQ];
  double azx[APC_MAXEQ*APC_MAXEQ];
  double azy[APC_MAXEQ*APC_MAXEQ];
  double azz[APC_MAXEQ*APC_MAXEQ];
  double bx[APC_MAXEQ*APC_MAXEQ];
  double by[APC_MAXEQ*APC_MAXEQ];
  double bz[APC_MAXEQ*APC_MAXEQ];
  double tx[APC_MAXEQ*APC_MAXEQ];
  double ty[APC_MAXEQ*APC_MAXEQ];
  double tz[APC_MAXEQ*APC_MAXEQ];
  double cval[APC_MAXEQ*APC_MAXEQ];
  double mval[APC_MAXEQ*APC_MAXEQ];
  double qx[APC_MAXEQ];
  double qy[APC_MAXEQ];
  double qz[APC_MAXEQ];
  double sval[APC_MAXEQ];
  double lval[APC_MAXEQ];

/*!!!!!! OLD OBSOLETE VERSION !!!!!!*/
  /* pde coefficients */
/*   static int ndofs_max=0;            /\* local dimension of the problem *\/ */
/*   static double *select_coeff=NULL; */
/*   static double  *axx=NULL, *axy=NULL, *axz=NULL, *ayx=NULL, *ayy=NULL; */
/*   static double  *ayz=NULL, *azx=NULL, *azy=NULL, *azz=NULL; */
/*   static double  *bx=NULL, *by=NULL, *bz=NULL, *tx=NULL, *ty=NULL, *tz=NULL;  */
/*   static double  *cval=NULL, *mval=NULL, *lval=NULL, *sval=NULL; */
/*   static double  *qx=NULL, *qy=NULL, *qz=NULL; */

/* #pragma omp threadprivate (ndofs_max) */
/* #pragma omp threadprivate (select_coeff) */
/* #pragma omp threadprivate (axx) */
/* #pragma omp threadprivate (axy) */
/* #pragma omp threadprivate (ayx) */
/* #pragma omp threadprivate (axz) */
/* #pragma omp threadprivate (ayz) */
/* #pragma omp threadprivate (ayy) */
/* #pragma omp threadprivate (ayz) */
/* #pragma omp threadprivate (azy) */
/* #pragma omp threadprivate (azz) */
/* #pragma omp threadprivate (bx) */
/* #pragma omp threadprivate (by) */
/* #pragma omp threadprivate (bz) */
/* #pragma omp threadprivate (tx) */
/* #pragma omp threadprivate (ty) */
/* #pragma omp threadprivate (tz) */
/* #pragma omp threadprivate (cval) */
/* #pragma omp threadprivate (mval) */
/* #pragma omp threadprivate (lval) */
/* #pragma omp threadprivate (sval) */
/* #pragma omp threadprivate (qx) */
/* #pragma omp threadprivate (qy) */
/* #pragma omp threadprivate (qz) */


  /* quadrature rules */
/*   static int problem_id_old=-1; /\* indicator for recomputing data *\/ */
/*   int pdeg;		/\* degree of polynomial *\/ */
/*   static int pdeg_old=-1; /\* indicator for recomputing quadrature data *\/ */
/*   int base;		/\* type of basis functions *\/ */
/*   static int base_old=-1; /\* indicator for recomputing quadrature data *\/ */
/*   static int ngauss;            /\* number of gauss points *\/ */
/*   static double xg[3000];   	 /\* coordinates of gauss points in 3D *\/ */
/*   static double wg[1000];       /\* gauss weights *\/ */

/* #pragma omp threadprivate (problem_id_old) */
/* #pragma omp threadprivate (pdeg_old) */
/* #pragma omp threadprivate (base_old) */
/* #pragma omp threadprivate (ngauss) */
/* #pragma omp threadprivate (xg) */
/* #pragma omp threadprivate (wg) */

  // to make old OpenMP compilers working
  int problem_id_old=-1; /* indicator for recomputing data */
  int pdeg;		/* degree of polynomial */
  int pdeg_old=-1; /* indicator for recomputing quadrature data */
  int base;		/* type of basis functions */
  int base_old=-1; /* indicator for recomputing quadrature data */
  int ngauss;            /* number of gauss points */
  double xg[3000];   	 /* coordinates of gauss points in 3D */
  double wg[1000];       /* gauss weights */

  int geo_order, num_geo_dofs;
  double geo_dofs[3*MMC_MAXELVNO];  /* coord of nodes of El */

  int nreq;	   /* number of equations */
  int el_mate;     /* material index in materials database */
  int num_shap;    /* number of element shape functions */
  int num_dofs;    /* number of element DOFs; usually num_dofs=nreq*num_shap */
  double determ;   /* determinant of jacobi matrix */
  double hsize=0;  /* size of an element */
  double vol;      /* volume for integration rule */
  double xcoor[3]={0.0};   /* global coord of gauss point */
  double uk_val[APC_MAXEQ]; /* computed solution from previous iteration */
  double uk_x[APC_MAXEQ];   /* x-derivatives of components of uk_val */
  double uk_y[APC_MAXEQ];   /* y-derivatives of components of uk_val */
  double uk_z[APC_MAXEQ];   /* z-derivatives of components of uk_val */
  double un_val[APC_MAXEQ]; /* computed solution from previous time step */
  double un_x[APC_MAXEQ];   /* x-derivatives of components of un_val */
  double un_y[APC_MAXEQ];   /* y-derivatives of components of un_val */
  double un_z[APC_MAXEQ];   /* z-derivatives of components of un_val */
  double base_phi[APC_MAXELVD];    /* basis functions */
  double base_dphix[APC_MAXELVD];  /* x-derivatives of basis function */
  double base_dphiy[APC_MAXELVD];  /* y-derivatives of basis function */
  double base_dphiz[APC_MAXELVD];  /* y-derivatives of basis function */

  int i, ki, kk, udofs, wdofs, iaux, weq, ueq , j, k;


/*++++++++++++++++ executable statements ++++++++++++++++*/

//#define COUNT_OPERATIONS
#ifdef COUNT_OPERATIONS
  int nr_oper=0;
#endif

/* prepare data necessary for numerical integration: */
/* geo_dofs - geometry degrees of freedom - used by apr_elem_calc_3D  */
/*            (for linear elements these are vertices coordinates) */
/* pdeg - indicator of element degree of approximation (this may be a single */
/*        number for simple (linear, quadratic, etc.) elements or a vector */
/*        with separate degrees of approximation for each mesh entity */
/*        (vertices, edges, faces, interior) forming an element */
/* num_shap - number of shape functions */
/* nreq - number of components for unknowns being vector fields */
/* num_dofs - number of degrees of freedom (often num_dofs = nreq*num_shap) */
/* base - type of basis functions (differentiates element types as well)*/
/*          examples (from include/aph_intf.h): */
/*          #define APC_BASE_PRISM_STD  3   // for linear prismatic elements */
/*          #define APC_BASE_TETRA_STD  4  // for linear tetrahedral elements */
  apr_prepare_integration_parameters(Field_id, El_id, 
				     &geo_order, &num_geo_dofs, geo_dofs, 
				     &el_mate, &base, &pdeg, 
				     &num_shap, &nreq, &num_dofs);
 

  /* degree of polynomial passed as argument overrides default element values */
  if(Pdeg_vec != NULL) {
    pdeg = Pdeg_vec[0]; // for simple one number Pdeg_vec
    num_shap = apr_get_el_pdeg_numshap(Field_id, El_id, &pdeg);
    nreq=apr_get_nreq(Field_id);
    num_dofs = nreq*num_shap;
  }

#ifdef DEBUG_APM
  if(nreq>APC_MAXEQ){
    printf("Number of equations %d greater than the limit %d.\n",nreq,APC_MAXEQ);
    printf("Change APC_MAXEQ in aph_std_lin_prism.h and recompile the code.\n");
    exit(-1);
  }
  if(num_shap>APC_MAXELVD){
    printf("Number of shape functions %d greater than the limit %d.\n",num_shap, APC_MAXELVD);
    printf("Change APC_MAXELVD in include/aph_intf.h and recompile the code.\n");
    exit(-1);
  }
  if(num_dofs>APC_MAXELSD){
    printf("Number of element degrees of freedom %d greater than the limit %d.\n", num_dofs, APC_MAXELSD);
    printf("Change APC_MAXELSD in include/aph_intf.h and recompile the code.\n");
    exit(-1);
  }
#endif


  // if selection concerned a different problem or has not been done yet
  if(Problem_id!=problem_id_old || coeff_ind == 0){
    pdr_select_el_coeff_vect(Problem_id, coeff_vect_ind); 
  // input to pdr_select_el_coeff_vect:
  // 0 - perform substitutions in coeff_vect_ind only if coeff_vect_ind[0]==0
  // output from pdr_select_el_coeff_vect:
  // 0 - coeff_vect_ind[0]==1 - always
  // 1 - mval, 2 - axx, 3 - axy, 4 - axz, 5 - ayx, 6 - ayy, 7 - ayz, 
  // 8 - azx, 9 - azy, 10 - azz, 11 - bx, 12 - by, 13 - bz
  // 14 - tx, 15 - ty, 16 - tz, 17 - cval
  // 18 - lval, 19 - qx, 20 - qy, 21 - qz, 22 - sval 
  }



/*!!!!!! OLD OBSOLETE VERSION !!!!!!*/
  /* if(Problem_id!=problem_id_old || select_coeff==NULL){ */
  /*   /\* allocate storage for coefficients and select the neeeded ones *\/ */
  /*   select_coeff=pdr_select_el_coeff(Problem_id, &mval, */
  /* 	 &axx,&axy,&axz,&ayx,&ayy,&ayz,&azx,&azy,&azz, */
  /* 	 &bx,&by,&bz,&tx,&ty,&tz,&cval,&lval,&qx,&qy,&qz,&sval); */
  /* } */

  /* prepare data for gaussian integration */
  if(pdeg!=pdeg_old || base != base_old){
    apr_set_quadr_3D(base, &pdeg, &ngauss, xg, wg);
    pdeg_old = pdeg;
    base_old = base;
  }

/*kbw
  //if(Field_id==2 && El_id==13753){
  if(Field_id>=0 && El_id>0){
    printf("In num_int_el: Field_id %d, element %d\n",Field_id, El_id);
    printf("pdeg %d, ngauss %d\n",pdeg, ngauss);
    printf("NREQ %d, ndof %d, local_dim %d\n",nreq, num_shap, num_dofs);
    printf("%d geo_dofs (nodes) with coordinates:\n", num_geo_dofs);
    for(i=0;i<num_geo_dofs;i++){
      printf("geo_dofs (node) %d: x - %f, y - %f, z - %f\n", 
	     i,
	     geo_dofs[3*i],geo_dofs[3*i+1],geo_dofs[3*i+2]);
    }
    printf("DOFS k:\n");
    for(i=0;i<num_dofs;i++) printf("%20.15lf",Sol_dofs_k[i]);
    printf("\n");
    //printf("DOFS n:\n");
    //for(i=0;i<num_dofs;i++) printf("%20.15lf",Sol_dofs_n[i]);
    printf("\n");
    //getchar();
  } 
/*kew*/

//#define TIME_TEST_2
#ifdef TIME_TEST_2
    double t00=0.0, t01=0.0, t02=0.0, t03=0.0, t04=0.0, t05=0.0;
#endif

  /*! ----------------------------------------------------------------------! */
  /*! --------------------- LOOP OVER GAUSS POINTS -------------------------! */
  /*! ----------------------------------------------------------------------! */
  for (ki=0;ki<ngauss;ki++) {
    
#ifdef TIME_TEST_2
      t00 = time_clock();
#endif

    /* at the gauss point, compute basis functions, determinant etc*/
    iaux = 2; /* calculations with jacobian but not on the boundary */
    determ = apr_elem_calc_3D(iaux, nreq, &pdeg, base, 
			      &xg[3*ki], geo_dofs, Sol_dofs_k,
			      base_phi,base_dphix,base_dphiy,base_dphiz,
			      xcoor,uk_val,uk_x,uk_y,uk_z,NULL);
 
#ifdef DEBUG_APM
/*jbw
    printf("ki=%d,x=%.4lf,y=%.4lf,z=%.4lf\n",ki,xg[3*ki],xg[3*ki+1],xg[3*ki+2]);
    for(i=0; i<15; i++) {
      printf("[%d] basep=%.4lf,basex=%.4lf,basey=%.4lf,basez=%.4lf\n",(i+1), base_phi[i],base_dphix[i],base_dphiy[i],base_dphiz[i]);
    }
    printf("u=%lf,ux=%lf,uy=%lf,uz=%lf\n",uk_val[0],uk_x[0],uk_y[0],uk_z[0]);
    getchar();
/*jbw*/
#endif
    
    
    vol = determ * wg[ki];

#ifdef TIME_TEST_2
      t01 += time_clock()-t00;
#endif

#ifdef COUNT_OPERATIONS
    // for DG
    //nr_oper += ((4+15)*num_shap + 3*(pdeg+1)*(pdeg+2)/2 + 6*(pdeg+1) + 220);
    // for STD - circa 500; definitely < 1000
    nr_oper += 21 + 24 + 15*(6+num_shap) + 12*6 + 21 + 8*num_shap*nreq;
#endif


/*kbw
    //if(norm_u>0.00010){
	if(El_id==49031){
	  int idofs;
	  double hsize = mmr_el_hsize(apr_get_mesh_id(Field_id),El_id,NULL,NULL,NULL);
	  printf("element %d, num_shap %d, u_x %lf, u_y %lf, u_z %lf\n",El_id, num_shap, uk_val[0], uk_val[1], uk_val[2]);
	  for (idofs = 0; idofs < num_shap; idofs++) {
	    printf("size: vol %.12lf, standard %.12lf, dphix %lf, dphiy %lf, dphiz %lf\n",vol, hsize,base_dphix[idofs], base_dphiy[idofs], base_dphiz[idofs]);
	  }
	  getchar();  getchar();  getchar();
	}
    //}
/*kew*/

    // if necessary compute the values of solution at the previous time step
    if(Sol_dofs_n != NULL){

      int ieq;
      for(ieq=0;ieq<nreq;ieq++) un_val[ieq]=0.0;
      for(i=0;i<num_shap;i++){
	for(ieq=0;ieq<nreq;ieq++){
	  un_val[ieq] += Sol_dofs_n[i*nreq+ieq]*base_phi[i];
	}
      }
      for(ieq=0;ieq<nreq;ieq++){
	un_x[ieq]=0.0;
	un_y[ieq]=0.0;
	un_z[ieq]=0.0;
      }
      for(i=0;i<num_shap;i++){
	for(ieq=0;ieq<nreq;ieq++){
	  un_x[ieq] += Sol_dofs_n[i*nreq+ieq]*base_dphix[i];
	  un_y[ieq] += Sol_dofs_n[i*nreq+ieq]*base_dphiy[i];
	  un_z[ieq] += Sol_dofs_n[i*nreq+ieq]*base_dphiz[i];
	}
      }

#ifdef COUNT_OPERATIONS
      // < 333 for std
      nr_oper += nreq*(num_shap*8+4);
      printf("\nNumber of operations in apr_num_int_el for single element %d - before\n", nr_oper);
#endif

    }

    // REMARK:
  //   1. udofs, ueq - correspond to the solution u and go from left to right
  //   2. wdofs, weq - correspond to test functions w and go up and down
  //   3. matrices are stored columnwise in vectors
  //   4. solution indices change in rows (go from left to right)
  //   5. test functions indices change in columns (go up and down)
  //   6. when matrices are stored in a vector the index is computed as:
  //      vector[udofs*nreq*num_dofs+wdofs*nreq+ueq*num_dofs+weq]
  //      (num_dofs=num_shap*nreq)
  //   7. for each pair (wdofs,udofs) there is a small submatrix nreq x nreq
  //      with indices ueq (from left to right) and weq (from top to bottom)
  //   8. when stiffness matrix entries are computed for vector problems with 
  //      the same shape functions for each component, then for each pair 
  //      (wdofs,udofs) a matrix of coefficients with the size nreq x nreq 
  //      is provided by the problem dependent module - the matrix is stored
  //      columnwise in the vector and the index is computed as [ueq*nreq+weq] 
  //   9. when load vector entries are computed for vector problems with 
  //      the same shape functions for each component, then for each wdofs index 
  //      a vector of coefficients with the size nreq is provided by the problem 
  //      dependent module
    
#ifdef TIME_TEST_2
      t02 += time_clock()-t00;
#endif

    /* get coefficients of convection-diffusion-reaction equations */
    pdr_el_coeff(Problem_id, El_id, el_mate, hsize, pdeg, &xg[3*ki],
		 base_phi, base_dphix, base_dphiy, base_dphiz,
		 xcoor, uk_val, uk_x, uk_y, uk_z, un_val, un_x, un_y, un_z,
		 mval, axx, axy, axz, ayx, ayy, ayz, azx, azy, azz,
		 bx, by, bz, tx, ty, tz, cval, lval, qx, qy, qz, sval);
    
#ifdef TIME_TEST_2
      t03 += time_clock()-t00;
#endif

/*kbw
    //static int juz_bylo=0;
    //if(juz_bylo==0){
    //juz_bylo = 1;
    if(Field_id==1 && El_id==13753) {

*//*
  if(Field_id>0 && El_id >= 0) {

      int i,ieq;
      
      printf("In num_int_el: Field_id %d, element %d\n",
	   Field_id, El_id);
      printf("pdeg %d, ngauss %d\n",pdeg, ngauss);
      printf("NREQ %d, ndof %d, local_dim %d\n",nreq, num_shap, num_dofs);
      printf("%d geo_dofs (nodes) with coordinates:\n", num_geo_dofs);
      for(i=0;i<num_geo_dofs;i++){
	printf("geo_dofs (node) %d: x - %f, y - %f, z - %f\n", 
	      i,
	      geo_dofs[3*i],geo_dofs[3*i+1],geo_dofs[3*i+2]);
      }

      printf("at gauss point %d, local coor %lf, %lf, %lf\n", 
	     ki,xg[3*ki],xg[3*ki+1],xg[3*ki+2]);
      printf("global coor %lf %lf %lf\n",xcoor[0],xcoor[1],xcoor[2]);
      printf("weight %lf, determ %lf, vol %lf\n",
	     wg[ki],determ,vol);
      printf("%d shape functions and derivatives: \n", num_shap);
      for(i=0;i<num_shap;i++){
	printf("fun - %lf, der: x - %lf, y - %lf, z - %lf\n",
	       base_phi[i],base_dphix[i],base_dphiy[i],base_dphiz[i]);
      }
      
      printf("previous iteration solution and derivatives: \n");
      for(ieq=0;ieq<nreq;ieq++){
	printf("component %d: ukx - %lf, der: x - %lf, y - %lf, z - %lf\n",
	       ieq, uk_val[ieq],uk_x[ieq],uk_y[ieq],uk_z[ieq]);
      }
      
      printf("previous timestep solution and derivatives: \n");
      for(ieq=0;ieq<nreq;ieq++){
      printf("component %d: unx - %lf, der: x - %lf, y - %lf, z - %lf\n",
	     ieq, un_val[0],un_x[0],un_y[0],un_z[0]);
      }

      if(mval!=NULL){
	printf("\ntime coeff LHS:\n");
	for(ueq=0;ueq<nreq*nreq;ueq++) printf("%20.12lf",mval[ueq]);
      }
      if(axx!=NULL){
	printf("\ndiffusion coeff axx:\n");
	for(ueq=0;ueq<nreq*nreq;ueq++) printf("%20.12lf",axx[ueq]);
	printf("\ndiffusion coeff axy:\n");
	for(ueq=0;ueq<nreq*nreq;ueq++) printf("%20.12lf",axy[ueq]);
	printf("\ndiffusion coeff axz:\n");
	for(ueq=0;ueq<nreq*nreq;ueq++) printf("%20.12lf",axz[ueq]);
      }
      if(ayy!=NULL){
	printf("\ndiffusion coeff ayx:\n");
	for(ueq=0;ueq<nreq*nreq;ueq++) printf("%20.12lf",ayx[ueq]);
	printf("\ndiffusion coeff ayy:\n");
	for(ueq=0;ueq<nreq*nreq;ueq++) printf("%20.12lf",ayy[ueq]);
	printf("\ndiffusion coeff ayz:\n");
	for(ueq=0;ueq<nreq*nreq;ueq++) printf("%20.12lf",ayz[ueq]);
      }
      if(azz!=NULL){
	printf("\ndiffusion coeff azx:\n");
	for(ueq=0;ueq<nreq*nreq;ueq++) printf("%20.12lf",azx[ueq]);
	printf("\ndiffusion coeff azy:\n");
	for(ueq=0;ueq<nreq*nreq;ueq++) printf("%20.12lf",azy[ueq]);
	printf("\ndiffusion coeff azz:\n");
	for(ueq=0;ueq<nreq*nreq;ueq++) printf("%20.12lf",azz[ueq]);
      }
      if(bx!=NULL){
	printf("\nconvection coeff bx:\n");
	for(ueq=0;ueq<nreq*nreq;ueq++) printf("%20.12lf",bx[ueq]);
	printf("\nconvection coeff by:\n");
	for(ueq=0;ueq<nreq*nreq;ueq++) printf("%20.12lf",by[ueq]);
	printf("\nconvection coeff bz:\n");
	for(ueq=0;ueq<nreq*nreq;ueq++) printf("%20.12lf",bz[ueq]);
      }
      if(tx!=NULL){
	printf("\nconvection coeff tx:\n");
	for(ueq=0;ueq<nreq*nreq;ueq++) printf("%20.12lf",tx[ueq]);
	printf("\nconvection coeff ty:\n");
	for(ueq=0;ueq<nreq*nreq;ueq++) printf("%20.12lf",ty[ueq]);
	printf("\nconvection coeff tz:\n");
	for(ueq=0;ueq<nreq*nreq;ueq++) printf("%20.12lf",tz[ueq]);
      }
      if(cval!=NULL){
	printf("\nreaction coeff c:\n");
	for(ueq=0;ueq<nreq*nreq;ueq++) printf("%20.12lf",cval[ueq]);
      }
      if(lval!=NULL){
	printf("\ntime coeff RHS:\n");
	for(ueq=0;ueq<nreq;ueq++) printf("%20.12lf",lval[ueq]);
      }
      if(qx!=NULL){
	printf("\nqx coeff RHS:\n");
	for(ueq=0;ueq<nreq;ueq++) printf("%20.12lf",qx[ueq]);
	printf("\nqy coeff RHS:\n");
	for(ueq=0;ueq<nreq;ueq++) printf("%20.12lf",qy[ueq]);
	printf("\nqz coeff RHS:\n");
	for(ueq=0;ueq<nreq;ueq++) printf("%20.12lf",qz[ueq]);
      }
      if(sval!=NULL){
	printf("\nsource:\n");
	for(ueq=0;ueq<nreq;ueq++) printf("%20.12lf",sval[ueq]);
      }
      getchar();
    }
/*kew*/

    if(Comp_sm==APC_COMP_SM||Comp_sm==APC_COMP_BOTH){
     
    if(coeff_vect_ind[1]==1){
       
	if(nreq==1){
	  
	  kk=0;
	  for (udofs=0;udofs<num_shap;udofs++) {
	    for (wdofs=0;wdofs<num_shap;wdofs++) {	
	      
	      Stiff_mat[kk+wdofs] += 
		mval[0] * base_phi[udofs] * base_phi[wdofs] * vol;
	      
	    }/* wdofs */
	    kk+=num_shap;      
	  } /* udofs */
	  
	}
	else if (nreq>1){  
		
	  kk=0;
	  for (udofs=0;udofs<num_shap;udofs++) {
	    for (wdofs=0;wdofs<num_shap;wdofs++) {

	      if(Diagonal[0]==0){
	
		for(ueq=0;ueq<nreq;ueq++){
		  for(weq=0;weq<nreq;weq++){

		    Stiff_mat[kk+wdofs*nreq+ueq*num_dofs+weq] +=
		      mval[ueq*nreq+weq] * base_phi[udofs] * base_phi[wdofs] * vol;

		  }
		}

	      } else {

		for(weq=0;weq<nreq;weq++){

		  Stiff_mat[kk+wdofs*nreq+weq*num_dofs+weq] +=
		    mval[weq*nreq+weq] * base_phi[udofs] * base_phi[wdofs] * vol;

		}

	      }

	    } /* wdofs */
	    kk += nreq*num_dofs;
	
	  } /* udofs */

	}

#ifdef COUNT_OPERATIONS
	if(Diagonal[0]==0) nr_oper += (2*nreq*nreq+2)*num_shap*num_shap;
	else nr_oper += (2*nreq+2)*num_shap*num_shap;
	printf("\nNumber of operations in apr_num_int_el for single element %d - m\n", nr_oper);
#endif

      }
      
      if(coeff_vect_ind[10]==1){
	if(coeff_vect_ind[2]==1){
	  if(nreq==1){
	    
	    kk=0;
	    for (udofs=0;udofs<num_shap;udofs++) {
	      for (wdofs=0;wdofs<num_shap;wdofs++) {	
		
		Stiff_mat[kk+wdofs] += ( 
			axx[0] *base_dphix[udofs] *base_dphix[wdofs] +
			axy[0] *base_dphiy[udofs] *base_dphix[wdofs] +
			axz[0] *base_dphiz[udofs] *base_dphix[wdofs] +
			ayx[0] *base_dphix[udofs] *base_dphiy[wdofs] +
			ayy[0] *base_dphiy[udofs] *base_dphiy[wdofs] +
			ayz[0] *base_dphiz[udofs] *base_dphiy[wdofs] +
			azx[0] *base_dphix[udofs] *base_dphiz[wdofs] +
			azy[0] *base_dphiy[udofs] *base_dphiz[wdofs] +
			azz[0] *base_dphiz[udofs] *base_dphiz[wdofs]  
				      ) * vol;
	      
	      }/* wdofs */
	      kk+=num_shap;
	    } /* udofs */
	  }
	  else if (nreq>1){  
	    
	    kk=0;
	    for (udofs=0;udofs<num_shap;udofs++) {
	      for (wdofs=0;wdofs<num_shap;wdofs++) {
		
		if(Diagonal[1]==0){
		  
		  for(ueq=0;ueq<nreq;ueq++){
		    for(weq=0;weq<nreq;weq++){
		      
		      
		      Stiff_mat[kk+wdofs*nreq+ueq*num_dofs+weq] += ( 
		    // axx[ueq*nreq+weq] * dw[weq]/dx * du[ueq]/dx 
		        axx[ueq*nreq+weq] *base_dphix[udofs] *base_dphix[wdofs] +
		    // axy[ueq*nreq+weq] * dw[weq]/dx * du[ueq]/dy
			axy[ueq*nreq+weq] *base_dphiy[udofs] *base_dphix[wdofs] +
		    // axz[ueq*nreq+weq] * dw[weq]/dx * du[ueq]/dz
			axz[ueq*nreq+weq] *base_dphiz[udofs] *base_dphix[wdofs] +
		    // ayx[ueq*nreq+weq] * dw[weq]/dy * du[ueq]/dx
			ayx[ueq*nreq+weq] *base_dphix[udofs] *base_dphiy[wdofs] +
		    // ayy[ueq*nreq+weq] * dw[weq]/dy * du[ueq]/dy
			ayy[ueq*nreq+weq] *base_dphiy[udofs] *base_dphiy[wdofs] +
		    // ayz[ueq*nreq+weq] * dw[weq]/dy * du[ueq]/dz
			ayz[ueq*nreq+weq] *base_dphiz[udofs] *base_dphiy[wdofs] +
		    // azx[ueq*nreq+weq] * dw[weq]/dz * du[ueq]/dx
			azx[ueq*nreq+weq] *base_dphix[udofs] *base_dphiz[wdofs] +
		    // azy[ueq*nreq+weq] * dw[weq]/dz * du[ueq]/dy
			azy[ueq*nreq+weq] *base_dphiy[udofs] *base_dphiz[wdofs] +
		    // azz[ueq*nreq+weq] * dw[weq]/dz * du[ueq]/dz
			azz[ueq*nreq+weq] *base_dphiz[udofs] *base_dphiz[wdofs]  
								    ) * vol;

		    }
		  }
		  
		} else {
		  
		  for(weq=0;weq<nreq;weq++){
		    
		    Stiff_mat[kk+wdofs*nreq+weq*num_dofs+weq] += ( 
		        axx[weq*nreq+weq] *base_dphix[udofs] *base_dphix[wdofs] +
			axy[weq*nreq+weq] *base_dphiy[udofs] *base_dphix[wdofs] +
			axz[weq*nreq+weq] *base_dphiz[udofs] *base_dphix[wdofs] +
			ayx[weq*nreq+weq] *base_dphix[udofs] *base_dphiy[wdofs] +
			ayy[weq*nreq+weq] *base_dphiy[udofs] *base_dphiy[wdofs] +
			ayz[weq*nreq+weq] *base_dphiz[udofs] *base_dphiy[wdofs] +
			azx[weq*nreq+weq] *base_dphix[udofs] *base_dphiz[wdofs] +
			azy[weq*nreq+weq] *base_dphiy[udofs] *base_dphiz[wdofs] +
			azz[weq*nreq+weq] *base_dphiz[udofs] *base_dphiz[wdofs]  
								  ) * vol;

		  }
		}
		
	      } /* wdofs */
	      kk += nreq*num_dofs;
	      
	    } /* udofs */
	    
	  }

#ifdef COUNT_OPERATIONS
	if(nreq==1) nr_oper += 28*num_shap*num_shap;
	else{
	  if(Diagonal[1]==0) nr_oper += 9*(2*nreq*nreq+2)*num_shap*num_shap;
	  else nr_oper += 9*(2*nreq+2)*num_shap*num_shap;
	}
	printf("\nNumber of operations in apr_num_int_el for single element %d - a\n", nr_oper);
#endif
	
	} else { // if axy == NULL 
	  if(nreq==1){
	  
	    kk=0;
	    for (udofs=0;udofs<num_shap;udofs++) {
	      for (wdofs=0;wdofs<num_shap;wdofs++) {	
		
		Stiff_mat[kk+wdofs] += ( 
			axx[0] *base_dphix[udofs] *base_dphix[wdofs] +
			ayy[0] *base_dphiy[udofs] *base_dphiy[wdofs] +
			azz[0] *base_dphiz[udofs] *base_dphiz[wdofs]  
				      ) * vol;
	      
	      }/* wdofs */
	      kk+=num_shap;
	    } /* udofs */
	  }
	  else if (nreq>1){  
	    
	    kk=0;
	    for (udofs=0;udofs<num_shap;udofs++) {
	      for (wdofs=0;wdofs<num_shap;wdofs++) {
		
		if(Diagonal[1]==0){
		  
		  for(ueq=0;ueq<nreq;ueq++){
		    for(weq=0;weq<nreq;weq++){

		      Stiff_mat[kk+wdofs*nreq+ueq*num_dofs+weq] += ( 
		        axx[ueq*nreq+weq] *base_dphix[udofs] *base_dphix[wdofs] +
			ayy[ueq*nreq+weq] *base_dphiy[udofs] *base_dphiy[wdofs] +
			azz[ueq*nreq+weq] *base_dphiz[udofs] *base_dphiz[wdofs]  
								    ) * vol;

		    }
		  }
		  
		} else {
		  
		  for(weq=0;weq<nreq;weq++){
		    
		    Stiff_mat[kk+wdofs*nreq+weq*num_dofs+weq] += ( 
		        axx[weq*nreq+weq] *base_dphix[udofs] *base_dphix[wdofs] +
			ayy[weq*nreq+weq] *base_dphiy[udofs] *base_dphiy[wdofs] +
			azz[weq*nreq+weq] *base_dphiz[udofs] *base_dphiz[wdofs]  
								  ) * vol;

		  }
		}
		
	      } /* wdofs */
	      kk += nreq*num_dofs;
	      
	    } /* udofs */
	    
	  }
 
#ifdef COUNT_OPERATIONS
	if(nreq==1) nr_oper += 10*num_shap*num_shap;
	else{
	  if(Diagonal[1]==0) nr_oper += 3*(2*nreq*nreq+2)*num_shap*num_shap;
	  else nr_oper += 3*(2*nreq+2)*num_shap*num_shap;
	}
#endif

	} // endif if azz!= NULL and axy == NULL 

      }
      else if(coeff_vect_ind[1]==1){
	
	if(nreq==1){
	  
	  kk=0;
	  for (udofs=0;udofs<num_shap;udofs++) {
	    for (wdofs=0;wdofs<num_shap;wdofs++) {	
	      
	      Stiff_mat[kk+wdofs] += (
			   axx[0] *base_dphix[udofs] *base_dphix[wdofs] +
			   axy[0] *base_dphiy[udofs] *base_dphix[wdofs] +
			   ayx[0] *base_dphix[udofs] *base_dphiy[wdofs] +
			   ayy[0] *base_dphiy[udofs] *base_dphiy[wdofs] 
				     ) * vol;
	      
	    }/* wdofs */
	    kk+=num_shap;
	  } /* udofs */
	}
	else if (nreq>1){  
		
	  kk=0;
	  for (udofs=0;udofs<num_shap;udofs++) {
	    for (wdofs=0;wdofs<num_shap;wdofs++) {
	
	      if(Diagonal[1]==0){
	
		for(ueq=0;ueq<nreq;ueq++){
		  for(weq=0;weq<nreq;weq++){

		    Stiff_mat[kk+wdofs*nreq+ueq*num_dofs+weq] += ( 
		        axx[ueq*nreq+weq] *base_dphix[udofs] *base_dphix[wdofs] +
			axy[ueq*nreq+weq] *base_dphiy[udofs] *base_dphix[wdofs] +
			ayx[ueq*nreq+weq] *base_dphix[udofs] *base_dphiy[wdofs] +
			ayy[ueq*nreq+weq] *base_dphiy[udofs] *base_dphiy[wdofs]
				      ) * vol;

		  }
		}

	      } else {

		for(weq=0;weq<nreq;weq++){

		  Stiff_mat[kk+wdofs*nreq+weq*num_dofs+weq] += ( 
		        axx[weq*nreq+weq] *base_dphix[udofs] *base_dphix[wdofs] +
			axy[weq*nreq+weq] *base_dphiy[udofs] *base_dphix[wdofs] +
			ayx[weq*nreq+weq] *base_dphix[udofs] *base_dphiy[wdofs] +
			ayy[weq*nreq+weq] *base_dphiy[udofs] *base_dphiy[wdofs]
				      ) * vol;

		}
	      }
	    } /* wdofs */
	    kk += nreq*num_dofs;
	
	  } /* udofs */
	
	} 

#ifdef COUNT_OPERATIONS
	if(nreq==1) nr_oper += 13*num_shap*num_shap;
	else{
	  if(Diagonal[1]==0) nr_oper += 4*(2*nreq*nreq+2)*num_shap*num_shap;
	  else nr_oper += 4*(2*nreq+2)*num_shap*num_shap;
	}
#endif

      }
      
      if(coeff_vect_ind[13]==1){
	
	if(nreq==1){
	  
	  kk=0;
	  for (udofs=0;udofs<num_shap;udofs++) {
	    for (wdofs=0;wdofs<num_shap;wdofs++) {	
	      
	      Stiff_mat[kk+wdofs] += ( bx[0] *base_dphix[udofs] *base_phi[wdofs] +
				      by[0] *base_dphiy[udofs] *base_phi[wdofs] +
				      bz[0] *base_dphiz[udofs] *base_phi[wdofs] 
				      ) *vol;
	      
	    }/* wdofs */
	    kk+=num_shap;      
	  } /* udofs */
	}
	else if (nreq>1){  
		
	  kk=0;
	  for (udofs=0;udofs<num_shap;udofs++) {
	    for (wdofs=0;wdofs<num_shap;wdofs++) {
	
	      if(Diagonal[2]==0){
	
		for(ueq=0;ueq<nreq;ueq++){
		  for(weq=0;weq<nreq;weq++){
		    
		    Stiff_mat[kk+wdofs*nreq+ueq*num_dofs+weq] += (
			   //  bx[ueq*nreq+weq] * w[weq] * du[ueq]/dx
			   bx[ueq*nreq+weq] *base_dphix[udofs] *base_phi[wdofs] +
			   //  by[ueq*nreq+weq] * w[weq] * du[ueq]/dy
			   by[ueq*nreq+weq] *base_dphiy[udofs] *base_phi[wdofs] +
			   //  bz[ueq*nreq+weq] * w[weq] * du[ueq]/dz
			   bz[ueq*nreq+weq] *base_dphiz[udofs] *base_phi[wdofs] 
				      ) *vol;
	      
		  }
		}

	      } else {

		for(weq=0;weq<nreq;weq++){
		  
		  Stiff_mat[kk+wdofs*nreq+weq*num_dofs+weq] += ( 
			   bx[weq*nreq+weq] *base_dphix[udofs] *base_phi[wdofs] +
			   by[weq*nreq+weq] *base_dphiy[udofs] *base_phi[wdofs] +
			   bz[weq*nreq+weq] *base_dphiz[udofs] *base_phi[wdofs] 
				      ) *vol;
	      
		}
	      }
	    } /* wdofs */
	    kk += nreq*num_dofs;
	
	  } /* udofs */
	
	} 
	
#ifdef COUNT_OPERATIONS
	if(nreq==1) nr_oper += 10*num_shap*num_shap;
	else{
	  if(Diagonal[2]==0) nr_oper += 3*(2*nreq*nreq+2)*num_shap*num_shap;
	  else nr_oper += 3*(2*nreq+2)*num_shap*num_shap;
	}
	printf("\nNumber of operations in apr_num_int_el for single element %d - b\n", nr_oper);
#endif

      }
      else if(coeff_vect_ind[11]==1){
	
	if(nreq==1){
	  
	  kk=0;
	  for (udofs=0;udofs<num_shap;udofs++) {
	    for (wdofs=0;wdofs<num_shap;wdofs++) {	
	      
	      Stiff_mat[kk+wdofs] += ( bx[0] *base_dphix[udofs] *base_phi[wdofs] +
				       by[0] *base_dphiy[udofs] *base_phi[wdofs] 
				      ) *vol;
	      
	    }/* wdofs */
	    kk+=num_shap;      
	  } /* udofs */
	}
	else if (nreq>1){  
		
	  kk=0;
	  for (udofs=0;udofs<num_shap;udofs++) {
	    for (wdofs=0;wdofs<num_shap;wdofs++) {
	
	      if(Diagonal[2]==0){
	
		for(ueq=0;ueq<nreq;ueq++){
		  for(weq=0;weq<nreq;weq++){
		    
		    Stiff_mat[kk+wdofs*nreq+ueq*num_dofs+weq] += ( 
			   bx[ueq*nreq+weq] *base_dphix[udofs] *base_phi[wdofs] +
			   by[ueq*nreq+weq] *base_dphiy[udofs] *base_phi[wdofs] 
				      ) *vol;
	      
		  }
		}

	      } else {

		for(weq=0;weq<nreq;weq++){
		  
		  Stiff_mat[kk+wdofs*nreq+weq*num_dofs+weq] += ( 
			   bx[weq*nreq+weq] *base_dphix[udofs] *base_phi[wdofs] +
			   by[weq*nreq+weq] *base_dphiy[udofs] *base_phi[wdofs] 
				      ) *vol;
	      
		}
	      }
	    } /* wdofs */
	    kk += nreq*num_dofs;
	
	  } /* udofs */
	
	} 

#ifdef COUNT_OPERATIONS
	if(nreq==1) nr_oper += 7*num_shap*num_shap;
	else{
	  if(Diagonal[2]==0) nr_oper += 2*(2*nreq*nreq+2)*num_shap*num_shap;
	  else nr_oper += 2*(2*nreq+2)*num_shap*num_shap;
	}
#endif

      }
      
      if(coeff_vect_ind[16]==1){
	
	if(nreq==1){
	  
	  kk=0;
	  for (udofs=0;udofs<num_shap;udofs++) {
	    for (wdofs=0;wdofs<num_shap;wdofs++) {	
	      
	      Stiff_mat[kk+wdofs] += ( 
				     tx[0] *base_phi[udofs] *base_dphix[wdofs] +
				     ty[0] *base_phi[udofs] *base_dphiy[wdofs] +
				     tz[0] *base_phi[udofs] *base_dphiz[wdofs] 
				      )	*vol;
	      
	    }/* wdofs */
	    kk+=num_shap;      
	  } /* udofs */
	}
	else if (nreq>1){  
		
	  kk=0;
	  for (udofs=0;udofs<num_shap;udofs++) {
	    for (wdofs=0;wdofs<num_shap;wdofs++) {
	
	      if(Diagonal[3]==0){
	
		for(ueq=0;ueq<nreq;ueq++){
		  for(weq=0;weq<nreq;weq++){

		    Stiff_mat[kk+wdofs*nreq+ueq*num_dofs+weq] += ( 
			   // tx[ueq*nreq+weq] * dw[weq]/dx * u[ueq] 
		           tx[ueq*nreq+weq] *base_phi[udofs] *base_dphix[wdofs] +
			   // ty[ueq*nreq+weq] * dw[weq]/dz * u[ueq] 
		           ty[ueq*nreq+weq] *base_phi[udofs] *base_dphiy[wdofs] +
			   // tz[ueq*nreq+weq] * dw[weq]/dz * u[ueq] 
			   tz[ueq*nreq+weq] *base_phi[udofs] *base_dphiz[wdofs] 
				      )	*vol;
	      
		  }
		}

	      } else {

		for(weq=0;weq<nreq;weq++){

		  Stiff_mat[kk+wdofs*nreq+weq*num_dofs+weq] += ( 
		           tx[weq*nreq+weq] *base_phi[udofs] *base_dphix[wdofs] +
		           ty[weq*nreq+weq] *base_phi[udofs] *base_dphiy[wdofs] +
			   tz[weq*nreq+weq] *base_phi[udofs] *base_dphiz[wdofs] 
				      )	*vol;
	      
		}
	      }
	    } /* wdofs */
	    kk += nreq*num_dofs;
	
	  } /* udofs */
	
	} 
	
#ifdef COUNT_OPERATIONS
	if(nreq==1) nr_oper += 10*num_shap*num_shap;
	else{
	  if(Diagonal[3]==0) nr_oper += 3*(2*nreq*nreq+2)*num_shap*num_shap;
	  else nr_oper += 3*(2*nreq+2)*num_shap*num_shap;
	}
	printf("\nNumber of operations in apr_num_int_el for single element %d - t\n", nr_oper);
#endif

      }
      else if(coeff_vect_ind[14]==1){
	
	if(nreq==1){
	  
	  kk=0;
	  for (udofs=0;udofs<num_shap;udofs++) {
	    for (wdofs=0;wdofs<num_shap;wdofs++) {	
	      
	      Stiff_mat[kk+wdofs] += ( tx[0] *base_phi[udofs] *base_dphix[wdofs] +
				       ty[0] *base_phi[udofs] *base_dphiy[wdofs] 
				     ) *vol;
	      
	    }/* wdofs */
	    kk+=num_shap;      
	  } /* udofs */
	}
	else if (nreq>1){  
		
	  kk=0;
	  for (udofs=0;udofs<num_shap;udofs++) {
	    for (wdofs=0;wdofs<num_shap;wdofs++) {
	
	      if(Diagonal[3]==0){
	
		for(ueq=0;ueq<nreq;ueq++){
		  for(weq=0;weq<nreq;weq++){
		    
		    Stiff_mat[kk+wdofs*nreq+ueq*num_dofs+weq] += ( 
		           tx[ueq*nreq+weq] *base_phi[udofs] *base_dphix[wdofs] +
		           ty[ueq*nreq+weq] *base_phi[udofs] *base_dphiy[wdofs] 
		      )	*vol;
		  }
		}
	      } else {

		for(weq=0;weq<nreq;weq++){
		  
		  Stiff_mat[kk+wdofs*nreq+weq*num_dofs+weq] += ( 
		           tx[weq*nreq+weq] *base_phi[udofs] *base_dphix[wdofs] +
		           ty[weq*nreq+weq] *base_phi[udofs] *base_dphiy[wdofs] 
				      )	*vol;
	      
		}
	      }
	    } /* wdofs */
	    kk += nreq*num_dofs;
	
	  } /* udofs */
	
	} 
		
#ifdef COUNT_OPERATIONS
	if(nreq==1) nr_oper += 7*num_shap*num_shap;
	else{
	  if(Diagonal[3]==0) nr_oper += 2*(2*nreq*nreq+2)*num_shap*num_shap;
	  else nr_oper += 2*(2*nreq+2)*num_shap*num_shap;
	}
#endif

      }
      
      if(coeff_vect_ind[17]==1){
	if(nreq==1){
	  kk=0;
	  for (udofs=0;udofs<num_shap;udofs++) {
	    for (wdofs=0;wdofs<num_shap;wdofs++) {	
	      
	      Stiff_mat[kk+wdofs] += cval[0]*base_phi[udofs]*base_phi[wdofs]*vol;
	      
	    }/* wdofs */
	    kk+=num_shap;      
	  } /* udofs */
	}
	else if (nreq>1){  
	  kk=0;
	  for (udofs=0;udofs<num_shap;udofs++) {
	    for (wdofs=0;wdofs<num_shap;wdofs++) {
	      if(Diagonal[4]==0){
		for(ueq=0;ueq<nreq;ueq++){
		  for(weq=0;weq<nreq;weq++){  
		    Stiff_mat[kk+wdofs*nreq+ueq*num_dofs+weq] += cval[ueq*nreq+weq]*base_phi[udofs]*base_phi[wdofs]*vol;
		  }
		}
	      } else {
		for(weq=0;weq<nreq;weq++){
		  Stiff_mat[kk+wdofs*nreq+weq*num_dofs+weq] += cval[weq*nreq+weq]*base_phi[udofs]*base_phi[wdofs]*vol;
		}
	      }
	    } /* wdofs */
	    kk += nreq*num_dofs;
	  } /* udofs */
	} 
	
#ifdef COUNT_OPERATIONS
	if(nreq==1) nr_oper += 4*num_shap*num_shap;
	else{
	  if(Diagonal[4]==0) nr_oper += (2*nreq*nreq+2)*num_shap*num_shap;
	  else nr_oper += (2*nreq+2)*num_shap*num_shap;
	}
	printf("\nNumber of operations in apr_num_int_el for single element %d - c\n", nr_oper);
#endif

      }
      
    } /* end if computing SM */
    
    if(Comp_sm==APC_COMP_RHS||Comp_sm==APC_COMP_BOTH){
      if(coeff_vect_ind[18]==1){
	if(nreq==1){
	  kk=0;
	  for (wdofs=0;wdofs<num_shap;wdofs++) {  
	    Rhs_vect[kk] += lval[0] * base_phi[wdofs] * vol; 
	    kk++;      
	  }/* wdofs */    
	}
	else if (nreq>1){  	
	  for (wdofs=0;wdofs<num_shap;wdofs++) {
	    for(weq=0;weq<nreq;weq++){
	      Rhs_vect[wdofs*nreq+weq] += lval[weq] * base_phi[wdofs] * vol; 
	    }
	  } /* wdofs */
	} 
	
#ifdef COUNT_OPERATIONS
	if(nreq==1) nr_oper += 3*num_shap;
	else nr_oper += (2*nreq+1)*num_shap;
	printf("\nNumber of operations in apr_num_int_el for single element %d - l\n", nr_oper);
#endif

      }

      if(coeff_vect_ind[21]==1){
	if(nreq==1){
	  kk=0;
	  for (wdofs=0;wdofs<num_shap;wdofs++) {
	    Rhs_vect[kk] += ( qx[0] * base_dphix[wdofs] +
			      qy[0] * base_dphiy[wdofs] +
			      qz[0] * base_dphiz[wdofs] ) * vol;
	    kk++;      
	  }/* wdofs */    
	}
	else if (nreq>1){  	
	  for (wdofs=0;wdofs<num_shap;wdofs++) {
	    for(weq=0;weq<nreq;weq++){
	      Rhs_vect[wdofs*nreq+weq] += ( 
			      qx[weq] * base_dphix[wdofs] +
			      qy[weq] * base_dphiy[wdofs] +
			      qz[weq] * base_dphiz[wdofs] ) * vol;
	    }
	  } /* wdofs */
	} 
	
#ifdef COUNT_OPERATIONS
	if(nreq==1) nr_oper += 7*num_shap;
	else nr_oper += 3*(2*nreq+1)*num_shap;
	printf("\nNumber of operations in apr_num_int_el for single element %d - q\n", nr_oper);
#endif

      } 
      else if(coeff_vect_ind[19]==1){
	if(nreq==1){ 
	  kk=0;
	  for (wdofs=0;wdofs<num_shap;wdofs++) { 
	    Rhs_vect[kk] += ( qx[0] * base_dphix[wdofs] + qy[0] * base_dphiy[wdofs] ) * vol; 
	    kk++;      
	  }/* wdofs */
	}
	else if (nreq>1){  	
	  for (wdofs=0;wdofs<num_shap;wdofs++) {
	    for(weq=0;weq<nreq;weq++){
	      Rhs_vect[wdofs*nreq+weq] += (qx[weq] * base_dphix[wdofs] + qy[weq] * base_dphiy[wdofs] ) * vol;
	    }
	  } /* wdofs */
	} 
	
#ifdef COUNT_OPERATIONS
	if(nreq==1) nr_oper += 5*num_shap;
	else nr_oper += 2*(2*nreq+1)*num_shap;
#endif

      }
      
      if(coeff_vect_ind[22]==1){
	
	if(nreq==1) { 
	  kk=0;
	  for (wdofs=0;wdofs<num_shap;wdofs++) {
	    
	    Rhs_vect[kk] += sval[0] * base_phi[wdofs] * vol;
	    
	    kk++;      
	  }/* wdofs */
	}
	
	else if (nreq>1){  
	  for (wdofs=0;wdofs<num_shap;wdofs++) {
	    for(weq=0;weq<nreq;weq++){
	      Rhs_vect[wdofs*nreq+weq] += sval[weq] * base_phi[wdofs] * vol; 
	    }
	  } /* wdofs */ 
	} 
	
#ifdef COUNT_OPERATIONS
	if(nreq==1) nr_oper += 3*num_shap;
	else nr_oper += (2*nreq+1)*num_shap;
	printf("\nNumber of operations in apr_num_int_el for single element %d - s\n", nr_oper);
#endif

      } 
    } /* end if computing RHSV */
    
#ifdef TIME_TEST_2
    t04 += time_clock()-t00;
#endif

  } /* end loop over integration points: ki */


#ifdef TIME_TEST_2
  printf("EXECUTION TIME: t01 %lf, t02 %lf, t03 %lf, t04 %lf, t05 %lf\n",t01, t02-t01, t03-t02, t04-t03, t05-t04);
#endif

#ifdef COUNT_OPERATIONS
  printf("\nNumber of operations in apr_num_int_el for single element %d\n", nr_oper);
#endif
  
  return(1);
}
  


/*------------------------------------------------------------
  apr_get_stiff_mat_data - to return data on dof entities for an element and 
                      to compute or rewrite element's stiffness matrix and RHSV
------------------------------------------------------------*/
int apr_get_stiff_mat_data(
  int Field_id,   /* in: approximation field ID  */
  int El_id,      /* in: unique identifier of the element */
  int Comp_sm,    /* in: indicator for the scope of computations: */
                  /*   APC_NO_COMP  - do not compute anything */
                  /*   APC_COMP_SM - compute entries to stiff matrix only */
                  /*   APC_COMP_RHS - compute entries to rhs vector only */
                  /*   APC_COMP_BOTH - compute entries for sm and rhsv */
                  /*   APC_REWR_SM - rewrite only entries to stiff matrix only */
                  /*   APC_REWR_RHS - rewrite only entries to rhs vector only */
                  /*   APC_REWR_BOTH - rewrite only entries for sm and rhsv */
  char  Transpose,/* in: perform transposition while rewriting */
                  /*     'y' or 'Y' - yes, otherwise - no */
  int Pdeg_in,    /* in: enforced degree of polynomial (if > 0 ) */
  int Nreq_in,    /* in: enforced nreq (if > 0 ) */
  int* Nr_dof_ent,/* in: size of arrays, */
                  /* out: no of filled entries, i.e. number of mesh entities*/
                  /* with which dofs and stiffness matrix blocks are associated */
  int* List_dof_ent_type,   /* out: list of no of dofs for 'dof' entity */
  int* List_dof_ent_id,     /* out: list of no of dofs for 'dof' entity */
  int* List_dof_ent_nrdofs, /* out: list of no of dofs for 'dof' entity */
  /* below: optional means when Stiff_mat and Rhs_vect are computed */
  int* Nrdofs_loc,/* in(optional): size of Stiff_mat and Rhs_vect */
                  /* out(optional): actual number of dofs per integration entity*/
  /* for matrices passing NULL is possible if suitable Comp_sm is passed */
  double* Stiff_mat,      /* out(optional): stiffness matrix stored columnwise*/
  double* Rhs_vect        /* out(optional): rhs vector */
  // REMARK:
  //   1. udofs, ueq - correspond to the solution and go from left to right
  //   2. wdofs, weq - correspond to test functions and go up and down
  //   3. matrices are stored columnwise in vectors
  //   4. solution indices change in rows (go from left to right)
  //   5. test functions indices change in columns (go up and down)
  //   6. when matrices are stored in a vector the index is computed as:
  //      vector[udofs*nreq*num_dofs+wdofs*nreq+ueq*num_dofs+weq]
  //      (num_dofs=num_shap*nreq)
  //   7. for each pair (wdofs,udofs) there is a small submatrix nreq x nreq
  //      with indices ueq (from left to right) and weq (from top to bottom)
			     )
{
  apt_field *field_p;
  int mesh_id, nreq, is_constrained;
  int num_shap, num_dofs, pdeg;
    
  int i,j,ino,ive,ied;
  int constr[9], constr_type[9];
  int el_verts[MMC_MAXELVNO+1];  /* list of nodes of El */
  
  int *el_edges; /* list of edges of El */
  int *el_nodes_constr; /* constraint vert & edges list */
  
  int num_shap_constr;  /* number of element shape functions with constraints */
  int num_dofs_constr;  /* number of element degrees of freedom */
                        /* with constraints taken into account */
  double constr_val;
  double **constr_el_loc;

  /* in general APC_MAXELVD should be used instead of 15 */
  double *stiff_loc;
  double *rhs_loc;
  
  /* auxiliary constans for TETRA & PRISMS */
  const int num_shap_tetra[3] = { 10, 4, 6 };
  const int num_shap_prism[3] = { 15, 6, 9 };

  /* auxiliary variable */
  int num_shap_vert; int nr_vert;
  int num_shap_edge; int nr_edge;
  int constr_el_loc_last_index;


/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select field and allocate data structures */
  field_p = apr_select_field(Field_id);
  
  /* check PDEG value */
  apr_set_and_test_pdeg(Field_id,NULL);

  /* checking */
  if(Transpose=='y'||Transpose=='Y'){
    printf("\e[1;31mERROR\e[0;0m] Row-major stiffnes matrices not yet supported in apr_get_stiff_mata_data!\n");
    exit(-1);
  }

  /* alloc memory for data structures */
  if(field_p->pdeg == APC_LINEAR_APPROXIMATION_PDEG) {
    el_edges = NULL;
    el_nodes_constr = (int*) malloc ((MMC_MAXELVNO+1)*sizeof(int));
    
    constr_el_loc = (double**) malloc ((MMC_MAXELVNO+1)*sizeof(double*));
    for(i=0; i<(MMC_MAXELVNO+1); i++) {
      constr_el_loc[i] = (double*) malloc (4*sizeof(double));
    }

    stiff_loc = (double*) malloc ((APC_MAXEQ*6*APC_MAXEQ*6)*sizeof(double));
    rhs_loc = (double*) malloc ((APC_MAXEQ*6)*sizeof(double));
    
  }
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    el_edges = (int*) malloc ((APC_MAXELVNO_ON_EDGES+1)*sizeof(int));
    el_nodes_constr = (int*) malloc ((MMC_MAXELVNO+APC_MAXELVNO_ON_EDGES+1)*sizeof(int));

    constr_el_loc = (double**) malloc ((MMC_MAXELVNO+APC_MAXELVNO_ON_EDGES+1)*sizeof(double*));
    for(i=0; i<(MMC_MAXELVNO+APC_MAXELVNO_ON_EDGES+1); i++) {
      constr_el_loc[i] = (double*) malloc (8*sizeof(double));
    }

    stiff_loc = (double*) malloc ((APC_MAXEQ*15*APC_MAXEQ*15)*sizeof(double));
    rhs_loc = (double*) malloc ((APC_MAXEQ*15)*sizeof(double));
  }

  /* select the corresponding mesh */
  mesh_id = apr_get_mesh_id(Field_id);
  
  if(Nreq_in>0) nreq = Nreq_in;
  else{
    nreq=apr_get_nreq(Field_id);
  }
  
  /* auxiliary number of active vertex and edges */
  nr_vert = mmr_get_nr_node(mesh_id);
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    nr_edge = mmr_get_nr_edge(mesh_id);
  }

  /* find degree of polynomial and number of element scalar dofs */
  //if(Pdeg_in>0) pdeg = Pdeg_in;
  //else{
  //  apr_get_el_pdeg(Field_id, El_id, &pdeg);
  //}
  pdeg = field_p->pdeg;
  num_shap = apr_get_el_pdeg_numshap(Field_id, El_id, &pdeg);

  if(field_p->pdeg == APC_LINEAR_APPROXIMATION_PDEG) {
    if(num_shap == num_shap_prism[1]) {
      num_shap_vert=num_shap_prism[1];
    } else if(num_shap == num_shap_tetra[1]) {
      num_shap_vert=num_shap_tetra[1];
    }
  }
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    if(num_shap == num_shap_prism[0]) {
      num_shap_vert=num_shap_prism[1];
      num_shap_edge=num_shap_prism[2];
    } else if(num_shap == num_shap_tetra[0]) {
      num_shap_vert=num_shap_tetra[1];
      num_shap_edge=num_shap_tetra[2];
    } else {
      printf("\e[1;31mERROR\e[0;0m] Object with %d shape function not definited!!!\n",num_shap);
      exit(-1);
    }
  }
  num_dofs = num_shap*nreq;

  mmr_el_node_coor(mesh_id,El_id,el_verts,NULL);
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    mmr_el_edges(mesh_id,El_id,el_edges);
  }

  is_constrained = 0;
  num_shap_constr = 0;
  constr_el_loc_last_index = 1;

  /* Loop over vertices */
  for(ive=1; ive<=num_shap_vert; ive++,constr_el_loc_last_index++) {

    if(apr_get_ent_pdeg(Field_id,APC_VERTEX,el_verts[ive])>0) { //real node

      List_dof_ent_type[num_shap_constr] = APC_VERTEX;
      List_dof_ent_id[num_shap_constr] = el_verts[ive];
      List_dof_ent_nrdofs[num_shap_constr] = nreq;
      constr_el_loc[constr_el_loc_last_index][0]=1.0;
      el_nodes_constr[constr_el_loc_last_index]=1;

      /*jbw
      printf("List_dof_ent_id[%d] = %d ; constr_el_loc[%d][0] = %lf\n",num_shap_constr,List_dof_ent_id[num_shap_constr],constr_el_loc_last_index,constr_el_loc[constr_el_loc_last_index][0]);
      /*jbw*/

      num_shap_constr++;

    } else
    if(apr_get_ent_pdeg(Field_id,APC_VERTEX,el_verts[ive])==0) { //constraint node
      
      is_constrained = 1;

      //Get constraint data
      apr_get_constr_data(Field_id,el_verts[ive],APC_VERTEX,constr,constr_type);

      //Save data
      el_nodes_constr[constr_el_loc_last_index]=constr[0];

      /*jbw
      printf("el_verts_constr[%d] = %d\n",ive,constr[0]);
      /*jbw*/

      for(i=0,j=1; i<constr[0]; i++,j++) {

	/*jbw
	printf("\tnr: %d ; constr: %d ; constr_type: ",j,constr[j]);
	if(constr_type[j] == APC_VERTEX)
	  printf("APC_VERTEX ; ");
	else if(constr_type[j] == APC_EDGE)
	  printf("APC_EDGE ; ");
	/*jbw*/

	List_dof_ent_id[num_shap_constr] = constr[j];
	List_dof_ent_type[num_shap_constr] = constr_type[j];
	List_dof_ent_nrdofs[num_shap_constr] = nreq;
	constr_el_loc[constr_el_loc_last_index][i]=constr_coeff.val[constr_coeff.id[constr_type[0]]+i];

	/*jbw
	printf("List_dof_ent_id[%d] = %d ; constr_el_loc[%d][%d] = %lf\n",num_shap_constr,List_dof_ent_id[num_shap_constr],constr_el_loc_last_index,i,constr_el_loc[constr_el_loc_last_index][i]);
	/*jbw*/

	num_shap_constr++;
      }
    }
  }
  
  /* Loop over edges */
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    for(ied=1; ied<=num_shap_edge; ied++, constr_el_loc_last_index++) {
      if(apr_get_ent_pdeg(Field_id,APC_EDGE,el_edges[ied])>0) { //real node

	List_dof_ent_type[num_shap_constr] = APC_EDGE;
	List_dof_ent_id[num_shap_constr] = el_edges[ied];
	List_dof_ent_nrdofs[num_shap_constr] = nreq;
	constr_el_loc[constr_el_loc_last_index][0]=1.0;
	el_nodes_constr[constr_el_loc_last_index]=1;

	/*jbw
	  printf("List_dof_ent_id[%d] = %d ; constr_el_loc[%d][0] = %lf\n",num_shap_constr,List_dof_ent_id[num_shap_constr],constr_el_loc_last_index,constr_el_loc[constr_el_loc_last_index][0]);
	/*jbw*/

	num_shap_constr++;
	
      } else {
	if(apr_get_ent_pdeg(Field_id,APC_EDGE, el_edges[ied])==0) { //constraint node
      
	  is_constrained = 1;

	  //Get constraint data
	  apr_get_constr_data(Field_id,el_edges[ied],APC_EDGE,constr,constr_type);

	  //Save data
	  el_nodes_constr[constr_el_loc_last_index]=constr[0];

	  /*jbw
	    printf("el_edge_constr[%d] = %d\n",ied,constr[0]);
	  /*jbw*/

	  for(i=0,j=1; i<constr[0]; i++,j++) {

	    /*jbw
	      printf("\tnr: %d ; constr: %d ; constr_type: ",j,constr[j]);
	      if(constr_type[j] == APC_VERTEX)
	      printf("APC_VERTEX ; ");
	      else if(constr_type[j] == APC_EDGE)
	      printf("APC_EDGE ; ");
	    /*jbw*/

	    List_dof_ent_id[num_shap_constr] = constr[j];
	    List_dof_ent_type[num_shap_constr] = constr_type[j];
	    List_dof_ent_nrdofs[num_shap_constr] = nreq;
	    constr_el_loc[constr_el_loc_last_index][i]=constr_coeff.val[constr_coeff.id[constr_type[0]]+i];

	    /*jbw
	      printf("List_dof_ent_id[%d] = %d ; constr_el_loc[%d][%d] = %lf\n",num_shap_constr,List_dof_ent_id[num_shap_constr],constr_el_loc_last_index,i,constr_el_loc[constr_el_loc_last_index][i]);
	    /*jbw*/

	    num_shap_constr++;
	  }
	}
      }
    }
  }

  /*jbw
    printf("Press ENTER to continue..."); getchar(); getchar();
  /*jbw*/

  num_dofs_constr = num_shap_constr*nreq;
  *Nr_dof_ent = num_shap_constr;

  if(Comp_sm != APC_NO_COMP){
  
    if(*Nrdofs_loc<nreq*num_shap_constr){
      printf("\e[1;31mERROR\e[0;0m] Too small arrays Stiff_mat and Rhs_vect passed to apr_get_stiff_mat_data\n");
      printf("\e[1;31mERROR\e[0;0m] %d < %d. Exiting !!!", *Nrdofs_loc, num_shap_constr*nreq);
      exit(-1);
    }

    *Nrdofs_loc = num_dofs_constr;

    if(Comp_sm==APC_COMP_SM || Comp_sm==APC_COMP_RHS || Comp_sm==APC_COMP_BOTH || Comp_sm==APC_COMP_MM){
      printf("\e[1;33mWARNING\e[0;0m] May be some day... (APC_COMP_... not implemented yet\n");
      exit(-1);
    }

    if(is_constrained != 0) {

      int i_constr, j_constr, kk, kk_constr, idofs, jdofs;
      int nci, ncj, ieq, jeq;

      if(Comp_sm==APC_REWR_SM || Comp_sm==APC_REWR_BOTH || Comp_sm==APC_REWR_MM){

	/* rewrite Stiff_mat to local array */
	for(i=0; i<num_dofs*num_dofs; i++) { stiff_loc[i]=Stiff_mat[i]; }
	for(i=0; i<num_dofs_constr*num_dofs_constr; i++) { Stiff_mat[i]=0.0; }

	/*jbw
	printf("\t\e[1;33m LOC_STIFF MAT DATA\e[0m\n");
	printf("Matrix size: %d, Number of elements %d\n",num_dofs,(num_dofs*num_dofs));
	for(i=0; i<num_dofs; i++) {
	  printf("[%d] ",(i*num_dofs));
	  for(j=0; j<num_dofs; j++) {
	    printf("%4.4lf, ",stiff_loc[i*num_dofs+j]);
	  }
	  printf("\n");
	}
	printf("Press ENTER to continue..."); getchar(); getchar();
	/*jbw*/

	j_constr = 0; i_constr = 0;
	kk=0; kk_constr = 0;
	for (jdofs=0;jdofs<num_shap;jdofs++) {
	  //printf("jdofs = %d\n",jdofs);

	  for(ncj=0; ncj<el_nodes_constr[jdofs+1]; ncj++) {
	    i_constr = 0;
	    
	    for (idofs=0; idofs<num_shap; idofs++) {	      
	      for(nci=0; nci<el_nodes_constr[idofs+1]; nci++) {		
		for(jeq=0;jeq<nreq;jeq++){		  
		  for(ieq=0;ieq<nreq;ieq++){
		    Stiff_mat[kk_constr+i_constr*nreq+jeq*num_dofs_constr+ieq] += 
		      constr_el_loc[jdofs+1][ncj]*
		      constr_el_loc[idofs+1][nci]*
		      stiff_loc[kk+idofs*nreq+jeq*num_dofs+ieq];

		    /*jbw
		    printf("constr_el_loc[%d][%d] = %lf ; ",(jdofs+1),ncj,constr_el_loc[jdofs+1][ncj]);
		    printf("constr_el_loc[%d][%d] = %lf ; ",(idofs+1),nci,constr_el_loc[idofs+1][nci]);
		    printf("Stiff_mat id: %d ; stiff_loc[%d] = %lf\n",(kk_constr+i_constr*nreq+jeq*num_dofs_constr+ieq),(kk+idofs*nreq+jeq*num_dofs+ieq),stiff_loc[kk+idofs*nreq+jeq*num_dofs+ieq]);
		    /*jbw*/
	
		    /*kbw
		      printf("VERT: rewriting:\nold: \tjdofs %d, idofs %d, kk %d, index %d value %lf\n",jdofs, idofs, kk, kk+idofs*nreq, stiff_loc[kk+idofs*nreq]);
		      printf("new: \tncj %d, nci %d, kk_constr %d, i_constr %d, index %d, constr %lf, value %lf\n",
		      ncj, nci, kk_constr, i_constr, 
		      kk_constr+i_constr*nreq+jeq*num_dofs_constr+ieq,
		      constr_el_loc[jdofs+1][ncj]*constr_el_loc[idofs+1][nci],
		      Stiff_mat[kk_constr+i_constr*nreq+jeq*num_dofs_constr+ieq]
		      );
		    /*kew*/

		    /*jbw
		      printf("\e[1;37mStiff_mat[%3d] = %lf\e[0m\n",(kk_constr+i_constr*nreq+jeq*num_dofs_constr+ieq),Stiff_mat[kk_constr+i_constr*nreq+jeq*num_dofs_constr+ieq]);
		    /*jbw*/	
		  }
		}
		i_constr++;
	      }  
	    } /* idofs */
	    j_constr++;
	    kk_constr += nreq*num_dofs_constr;
	  }

	  /*jbw
	  printf("Press ENTER to continue..."); getchar(); getchar();
	  /*jbw*/
	  
	  kk += nreq*num_dofs;
	} /* jdofs */

	/*jbw
	printf("\t\e[1;33m LOC_STIFF MAT DATA\e[0m\n");
	printf("Matrix size: %d, Number of elements %d\n",num_dofs,(num_dofs*num_dofs));
	for(i=0; i<num_dofs; i++) {
	  for(j=0; j<num_dofs; j++) {
	    printf("%4.4lf, ",stiff_loc[i*num_dofs+j]);
	    
	  }
	  printf("\n");
	}

	printf("\t\e[1;33m STIFF MAT DATA\e[0m\n");
	printf("Matrix size: %d, Number of elements %d\n",num_dofs_constr,(num_dofs_constr*num_dofs_constr));
	for(i=0; i<num_dofs_constr; i++) {
	  for(j=0; j<num_dofs_constr; j++) {
	    printf("%2.2lf, ",Stiff_mat[i*num_dofs_constr+j]);
	    
	  }
	  printf("\n");
	}
	printf("\n");
	printf("Press ENTER to continue..."); getchar(); getchar();
	/*jbw*/

      } // end if rewriting SM

      if(Comp_sm==APC_REWR_RHS || Comp_sm==APC_REWR_BOTH){

	//printf("APC_REWR_RGS or APC_REWR_BOTH- Enabled!!!\n");

	/* rewrite Rhs_vect to local array */
	for(i=0;i<num_dofs;i++) rhs_loc[i]=Rhs_vect[i];
	for(i=0;i<num_dofs_constr;i++) Rhs_vect[i]=0;

	i_constr = 0;
	for (idofs=0;idofs<num_shap;idofs++) {
	  for(nci=0; nci<el_nodes_constr[idofs+1]; nci++) {
	    for(ieq=0;ieq<nreq;ieq++){
	      Rhs_vect[i_constr*nreq+ieq] += constr_el_loc[idofs+1][nci] * rhs_loc[idofs*nreq+ieq];	
	    }
	    i_constr++;
	  }
	}/* idofs */     
      } // end if rewriting RHS
    } // end if constrained nodes exist    
  } // end if considering SM and RHS

  /* free memory */
  SAFE_FREE(el_edges);
  if(field_p->pdeg == APC_LINEAR_APPROXIMATION_PDEG) {
    for(i=0; i<(MMC_MAXELVNO+1); i++) {
      SAFE_FREE(constr_el_loc[i]);
    }
  }
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) {
    for(i=0; i<(MMC_MAXELVNO+APC_MAXELVNO_ON_EDGES+1); i++) {
      SAFE_FREE(constr_el_loc[i]);
    }
  }
  SAFE_FREE(el_nodes_constr); SAFE_FREE(constr_el_loc); SAFE_FREE(stiff_loc); SAFE_FREE(rhs_loc);

  return(1);
}


/*------------------------------------------------------------
  apr_proj_dof_ref - to rewrite dofs after modifying the mesh
              (the procedure also recreates constraints data) 
------------------------------------------------------------*/
int apr_proj_dof_ref(
  int Field_id,    /* in: approximation field ID  */
  int El,	   /* in:  >0 - rewrite after one [de]refinement of el */
                   /*     <=0 - rewrite after massive [de]refinements */
  int Max_elem_id_before, /* in: maximal element (face, etc.) id */
  int Max_face_id_before, /*     before !!!!!!! refinements */
  int Max_edge_id_before,
  int Max_vert_id_before
  )
{
  /* pointer to field structure */
  apt_field *field_p;
  int mesh_id;
  int nr_dof_ents, nr_dof_vert, nr_dof_edge;
  int nr_constr_dof_vert, nr_constr_dof_edge, nr_inactive_dof_edge;
  int nreq, nvert, nedge;

  int mid_node; //Central face node
  int mid_tria; //Central traingle face in divided face
  int coeff_id; //Coefficient id (auxiliary variable)

  int edge_nodes[3][3]; //Edge nodes
  int edge_sons[3]; //Edge sons

  int quad_fa_nodes[2][5]; //Quadangle face nodes
  int quad_fa_edges[2][4]; //Quadangle face edges
  int tria_fa_edges[4][3]; //Triangle face edges

  int fa_sons[5]; //Face sons

  int ned, nfa, ido, ino, ied, ifa; 
  int num_nodes[4], num_edges[4];

  int i,j,k,iaux;
  int *edge_elems;

  int node_shift;
  int face_neig[2];	/* list of neighbors */
  int neig_sides[2]={0,0};/* sides of face wrt neigs */

  int max_elem_id, max_face_id, max_edge_id, max_vert_id; 
	
/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* check PDEG value */
  apr_set_and_test_pdeg(Field_id,NULL);
  
  /* select the corresponding mesh */
  mesh_id = apr_get_mesh_id(Field_id);
  
  /* select the pointer to the approximation field */
  field_p = apr_select_field(Field_id);

  // check if mesh module is ready
  if(mmr_is_ready_for_proj_dof_ref(mesh_id) < 0) {
    printf("[\e[1;31mERROR\e[0;0m] Mesh module not ready for dofs projection!");
    exit(-1);
  }

  edge_elems = (int*) malloc((MMC_MAX_EDGE_ELEMS+1)*sizeof(int));

  max_elem_id = mmr_get_max_elem_id(mesh_id); 
  if(Max_elem_id_before>max_elem_id) max_elem_id=Max_elem_id_before;
  max_face_id = mmr_get_max_face_id(mesh_id); 
  if(Max_face_id_before>max_face_id) max_face_id=Max_face_id_before;
  max_edge_id = mmr_get_max_edge_id(mesh_id); 
  if(Max_edge_id_before>max_edge_id) max_edge_id=Max_edge_id_before;
  max_vert_id = mmr_get_max_node_id(mesh_id); 
  if(Max_vert_id_before>max_vert_id) max_vert_id=Max_vert_id_before;
  
  printf("Sizing dof_vert from %d to %d \n",Max_vert_id_before,max_vert_id);
  // check if there is a need of resizing dofs structure
  if(max_vert_id > field_p->capacity_dof_vert) {
    apt_dof_ent* new_dof_vert = (apt_dof_ent*) malloc((2*max_vert_id+1)*sizeof(apt_dof_ent));
    if(new_dof_vert==NULL) {
      printf("\e[1;31mERROR\e[0;0m] Dofs VERTEX structures not reallocated (new size:%d)!",2*max_vert_id+1);
      exit(-1);
    }
	
    // fast rewriting existing solution
    memcpy(new_dof_vert,field_p->dof_vert,(Max_vert_id_before+1)*sizeof(apt_dof_ent));

    // fast NULL-ing new space for dof_ents
    memset(new_dof_vert+(Max_vert_id_before+1)*sizeof(apt_dof_ent),0,2*max_vert_id-Max_vert_id_before);
    // KB: shouldn't it be multiplied by sizeof(apt_dof_ent)
    printf("[\e[1;33mWARNING\e[0;0m] Usage of memset in apr_proj_dof_ref\n");
    getchar();getchar();getchar();

    SAFE_FREE(field_p->dof_vert);
    field_p->dof_vert = new_dof_vert;
    field_p->capacity_dof_vert = 2*max_vert_id;
  }

  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
    printf("Sizing dof_edge from %d to %d \n",Max_edge_id_before,max_edge_id);
    // check if there is a need of resizing dofs structure
    if(max_edge_id > field_p->capacity_dof_edge) {
      apt_dof_ent* new_dof_edge = (apt_dof_ent*) malloc((2*max_edge_id+1)*sizeof(apt_dof_ent));
      if(new_dof_edge==NULL) {
	printf("\e[1;31mERROR\e[0;0m] Dofs EDGE structures not reallocated (new size:%d)!",2*max_edge_id+1);
	exit(-1);
      }
	
      // fast rewriting existing solution
      memcpy(new_dof_edge,field_p->dof_edge,(Max_edge_id_before+1)*sizeof(apt_dof_ent));

      // fast NULL-ing new space for dof_ents
      memset(new_dof_edge+(Max_edge_id_before+1)*sizeof(apt_dof_ent),0,2*max_edge_id-Max_edge_id_before);
      // KB: shouldn't it be multiplied by sizeof(apt_dof_ent)
      printf("[\e[1;33mWARNING\e[0;0m] Usage of memset in apr_proj_dof_ref\n");
      getchar();getchar();getchar();

      SAFE_FREE(field_p->dof_edge);
      field_p->dof_edge = new_dof_edge;
      field_p->capacity_dof_edge = 2*max_edge_id;
    }
  }
  
  nreq=apr_get_nreq(Field_id);

  /* to begin indicate that there are NO constrained nodes in the mesh */
  field_p->constr = APC_FALSE;

  /* allocate the space for elems structure for active and inactive edges */
  mmr_create_edge_elems(mesh_id, max_edge_id);

  /* REMOVE CONSTRAINED NODES FROM PREVIOUS MESH */
  for(nvert=1; nvert<=Max_vert_id_before; nvert++) {

    //Check vertices
    if(mmr_node_status(mesh_id,nvert) == MMC_ACTIVE) {
      //Only for constrained vertices
      if(field_p->dof_vert[nvert].vec_dof_1 == NULL) {

	//Alloc memory
	field_p->dof_vert[nvert].vec_dof_1 = (double *) malloc(nreq*sizeof(double));
	for(i=0;i<nreq;i++) { field_p->dof_vert[nvert].vec_dof_1[i]=0.0; }
	if(field_p->nr_sol>1){
	  field_p->dof_vert[nvert].vec_dof_2 = (double *) malloc(nreq*sizeof(double));
	  for(i=0;i<nreq;i++) { field_p->dof_vert[nvert].vec_dof_2[i]=0.0; }
	}
	if(field_p->nr_sol>2){
	  field_p->dof_vert[nvert].vec_dof_3 = (double *) malloc(nreq*sizeof(double));
	  for(i=0;i<nreq;i++) { field_p->dof_vert[nvert].vec_dof_3[i]=0.0; }
	}

	if(field_p->dof_vert[nvert].constr != NULL) {

	  //Calculate values
	  for(ino=1,j=0; ino<=field_p->dof_vert[nvert].constr[0]; ino++,j++) {	  
	    switch(field_p->dof_vert[nvert].constr_type[ino]) {
	    
	    case APC_VERTEX:
	      // -- DEBUG: constrained node points to an empty pointer 
	      assert(field_p->dof_vert[field_p->dof_vert[nvert].constr[ino]].vec_dof_1 != NULL);

	      for(i=0;i<nreq;i++) {
		field_p->dof_vert[nvert].vec_dof_1[i] +=
		  constr_coeff.val[constr_coeff.id[field_p->dof_vert[nvert].constr_type[0]]+j]*
		  field_p->dof_vert[field_p->dof_vert[nvert].constr[ino]].vec_dof_1[i];
	      }
	      if(field_p->nr_sol>1){
		for(i=0;i<nreq;i++) { 
		  field_p->dof_vert[nvert].vec_dof_2[i] +=
		    constr_coeff.val[constr_coeff.id[field_p->dof_vert[nvert].constr_type[0]]+j]*
		    field_p->dof_vert[field_p->dof_vert[nvert].constr[ino]].vec_dof_2[i];
		}
	      }
	      if(field_p->nr_sol>2){
		for(i=0;i<nreq;i++) { 
		  field_p->dof_vert[nvert].vec_dof_3[i] +=
		    constr_coeff.val[constr_coeff.id[field_p->dof_vert[nvert].constr_type[0]]+j]*
		    field_p->dof_vert[field_p->dof_vert[nvert].constr[ino]].vec_dof_3[i];
		}
	      }
	      break;

	    
	    case APC_EDGE:
	      // -- DEBUG: constrained node points to an empty pointer 
	      assert(field_p->dof_edge[field_p->dof_vert[nvert].constr[ino]].vec_dof_1 != NULL);
	    
	      for(i=0;i<nreq;i++) {
		field_p->dof_vert[nvert].vec_dof_1[i] +=
		  constr_coeff.val[constr_coeff.id[field_p->dof_vert[nvert].constr_type[0]]+j]*
		  field_p->dof_edge[field_p->dof_vert[nvert].constr[ino]].vec_dof_1[i];

	      }
	      if(field_p->nr_sol>1){
		for(i=0;i<nreq;i++) { 
		  field_p->dof_vert[nvert].vec_dof_2[i] +=
		    constr_coeff.val[constr_coeff.id[field_p->dof_vert[nvert].constr_type[0]]+j]*
		    field_p->dof_edge[field_p->dof_vert[nvert].constr[ino]].vec_dof_2[i];
		}
	      }
	      if(field_p->nr_sol>2){
		for(i=0;i<nreq;i++) { 
		  field_p->dof_vert[nvert].vec_dof_3[i] +=
		    constr_coeff.val[constr_coeff.id[field_p->dof_vert[nvert].constr_type[0]]+j]*
		    field_p->dof_edge[field_p->dof_vert[nvert].constr[ino]].vec_dof_3[i];
		}
	      }
	      break;

	    default:
	      exit(-1);
	    };
	  }

	  //Free constr value
	  SAFE_FREE(field_p->dof_vert[nvert].constr);
	  SAFE_FREE(field_p->dof_vert[nvert].constr_type);
	}
      }
      field_p->dof_vert[nvert].active = 1; // Set activity flag to true
    } else {
      field_p->dof_vert[nvert].active = 0; // Set activity flag to false
    }
  }

  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation

    for(ned=1; ned<=Max_edge_id_before; ned++) {

      //Check edges
      if(mmr_edge_status(mesh_id,ned)==MMC_INACTIVE || mmr_edge_status(mesh_id,ned)==MMC_ACTIVE) {

	//Only for constrained vertices
	if(field_p->dof_edge[ned].vec_dof_1 == NULL) {

	  //Alloc memory
	  field_p->dof_edge[ned].vec_dof_1 = (double *) malloc(nreq*sizeof(double));
	  for(i=0;i<nreq;i++) { field_p->dof_edge[ned].vec_dof_1[i]=0.0; }
	  if(field_p->nr_sol>1){
	    field_p->dof_edge[ned].vec_dof_2 = (double *) malloc(nreq*sizeof(double));
	    for(i=0;i<nreq;i++) { field_p->dof_edge[ned].vec_dof_2[i]=0.0; }
	  }
	  if(field_p->nr_sol>2){
	    field_p->dof_edge[ned].vec_dof_3 = (double *) malloc(nreq*sizeof(double));
	    for(i=0;i<nreq;i++) { field_p->dof_edge[ned].vec_dof_3[i]=0.0; }
	  }

	  if(field_p->dof_edge[ned].constr != NULL) {
	    
	    //Calculate values
	    for(ino=1,j=0; ino<=field_p->dof_edge[ned].constr[0]; ino++,j++) {
	      
	      // -- DEBUG: constrained node points to an empty pointer 
	      assert(field_p->dof_edge[field_p->dof_edge[ned].constr[ino]].vec_dof_1 != NULL);
	      
	      for(i=0;i<nreq;i++) {
		field_p->dof_edge[ned].vec_dof_1[i] +=
		  constr_coeff.val[constr_coeff.id[field_p->dof_edge[ned].constr_type[0]]+j]*
		  field_p->dof_edge[field_p->dof_edge[ned].constr[ino]].vec_dof_1[i];
	      }
	      if(field_p->nr_sol>1){
		for(i=0;i<nreq;i++) { 
		  field_p->dof_edge[ned].vec_dof_2[i] +=
		    constr_coeff.val[constr_coeff.id[field_p->dof_edge[ned].constr_type[0]]+j]*
		    field_p->dof_edge[field_p->dof_edge[ned].constr[ino]].vec_dof_2[i];
		}
	      }

	      if(field_p->nr_sol>2){
		for(i=0;i<nreq;i++) { 
		  field_p->dof_edge[ned].vec_dof_3[i] +=
		    constr_coeff.val[constr_coeff.id[field_p->dof_edge[ned].constr_type[0]]+j]*
		    field_p->dof_edge[field_p->dof_edge[ned].constr[ino]].vec_dof_3[i];
		}
	      }
	    }

	    //Free constr value
	    SAFE_FREE(field_p->dof_edge[ned].constr);
	    SAFE_FREE(field_p->dof_edge[ned].constr_type);
	  }
	}
	//Change activity flag
	if(field_p->dof_edge[ned].vec_dof_1 != NULL) {
	  field_p->dof_edge[ned].active = 1; // Set activity flag to true
	} else {
	  field_p->dof_edge[ned].active = 0; // Set activity flag to false
	}
      } else {
	field_p->dof_edge[ned].active = 0; // Set activity flag to false
      }
    }
  }
  
  /* PROJECTION OF SOLUTION IN NEW MESH */

  /* For all edges in the current and the previous mesh */
  /*
   *   o        o
   *   *        *    Note: If new vertex in the middle of edge is CONSTRAINT 
   *   *       (c)   then new edges are CONSTRAINT too.
   *   *        *
   *   x  <==   C
   *   *        *
   *   *       (c)
   *   *        *
   *   o        o
   *
   */  
  for(ned=1; ned<=max_edge_id; ned++) {

    //Check inactive edge
    if(mmr_edge_status(mesh_id,ned)==MMC_INACTIVE) {

	//Get edge sons
    mmr_edge_sons(mesh_id, ned, edge_sons, NULL);

      /* the node in the middle of a divided edge */
      mmr_edge_nodes(mesh_id, ned, edge_nodes[0]);
      mmr_edge_nodes(mesh_id, edge_sons[0], edge_nodes[1]);
      mmr_edge_nodes(mesh_id, edge_sons[1], edge_nodes[2]);
      
      //Set active middle node 
      nvert = edge_nodes[1][1];

      // -- DEBUG: vertex outside the space
      assert(nvert <= mmr_get_max_node_id(mesh_id));
      assert(nvert <= mmr_get_max_node_max(mesh_id));
      assert(nvert <= field_p->capacity_dof_vert);

      /* for inactive edges check whether all surrounding elements are active */
      mmr_edge_elems(mesh_id,ned,edge_elems);
      
      
      iaux = 0;
      for(i=1; i<=edge_elems[0]; i++) {
	if(mmr_el_status(mesh_id,edge_elems[i]) == MMC_INACTIVE) iaux++;
      }

      
      if(iaux == edge_elems[0]) {
	
	// -- DEBUG: ancestors are constrained
	assert(field_p->dof_vert[edge_nodes[0][0]].vec_dof_1 != NULL);
	assert(field_p->dof_vert[edge_nodes[0][1]].vec_dof_1 != NULL);
	if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
	  assert(field_p->dof_edge[ned].vec_dof_1 != NULL);
	}
	
	/* all edge elements are divided - a real node */
	
	// -------------------- CREATE REAL DOF_VERT -------------------- //
	
	//Remove constrained data if it needed
	SAFE_FREE(field_p->dof_vert[nvert].constr);
	SAFE_FREE(field_p->dof_vert[nvert].constr_type);
	
	//Alloc memory and calculate values for new node if it needed
	if(field_p->dof_vert[nvert].vec_dof_1 == NULL) {
	  
	  //Alloc memory
	  field_p->dof_vert[nvert].vec_dof_1 = (double *) malloc(nreq*sizeof(double));
	  for(i=0;i<nreq;i++) { field_p->dof_vert[nvert].vec_dof_1[i]=0.0; }
	  if(field_p->nr_sol>1){
	    field_p->dof_vert[nvert].vec_dof_2 = (double *) malloc(nreq*sizeof(double));
	    for(i=0;i<nreq;i++) { field_p->dof_vert[nvert].vec_dof_2[i]=0.0; }
	  }
	  if(field_p->nr_sol>2){
	    field_p->dof_vert[nvert].vec_dof_3 = (double *) malloc(nreq*sizeof(double));
	    for(i=0;i<nreq;i++) { field_p->dof_vert[nvert].vec_dof_3[i]=0.0; }
	  }

	  field_p->dof_vert[nvert].active = 1; //Set activity glag

	  //Calculate new values
	  if(field_p->pdeg == APC_LINEAR_APPROXIMATION_PDEG) { // Do it only for linear approximation
	    for(i=0;i<nreq;i++) {
	      field_p->dof_vert[nvert].vec_dof_1[i] += 
		constr_coeff.val[constr_coeff.id[1]+0]*field_p->dof_vert[edge_nodes[0][0]].vec_dof_1[i] + 
		constr_coeff.val[constr_coeff.id[1]+1]*field_p->dof_vert[edge_nodes[0][1]].vec_dof_1[i];

	      if(field_p->nr_sol>1) {
		field_p->dof_vert[nvert].vec_dof_2[i] += 
		  constr_coeff.val[constr_coeff.id[1]+0]*field_p->dof_vert[edge_nodes[0][0]].vec_dof_2[i] +
		  constr_coeff.val[constr_coeff.id[1]+1]*field_p->dof_vert[edge_nodes[0][1]].vec_dof_2[i];
	      }
	      
	      if(field_p->nr_sol>2) {
		field_p->dof_vert[nvert].vec_dof_3[i] +=
		  constr_coeff.val[constr_coeff.id[1]+0]*field_p->dof_vert[edge_nodes[0][0]].vec_dof_3[i] +
		  constr_coeff.val[constr_coeff.id[1]+1]*field_p->dof_vert[edge_nodes[0][1]].vec_dof_3[i];
	      }
	    }
	  }
	  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
	    for(i=0;i<nreq;i++) {
	      field_p->dof_vert[nvert].vec_dof_1[i] += 
		constr_coeff.val[constr_coeff.id[1]+0]*field_p->dof_vert[edge_nodes[0][0]].vec_dof_1[i] + 
		constr_coeff.val[constr_coeff.id[1]+1]*field_p->dof_vert[edge_nodes[0][1]].vec_dof_1[i] +
		constr_coeff.val[constr_coeff.id[1]+2]*field_p->dof_edge[ned].vec_dof_1[i];
	      
	      if(field_p->nr_sol>1) {
		field_p->dof_vert[nvert].vec_dof_2[i] += 
		  constr_coeff.val[constr_coeff.id[1]+0]*field_p->dof_vert[edge_nodes[0][0]].vec_dof_2[i] +
		  constr_coeff.val[constr_coeff.id[1]+1]*field_p->dof_vert[edge_nodes[0][1]].vec_dof_2[i] +
		  constr_coeff.val[constr_coeff.id[1]+2]*field_p->dof_edge[ned].vec_dof_2[i];
	      }
	      
	      if(field_p->nr_sol>2) {
		field_p->dof_vert[nvert].vec_dof_3[i] +=
		  constr_coeff.val[constr_coeff.id[1]+0]*field_p->dof_vert[edge_nodes[0][0]].vec_dof_3[i] +
		  constr_coeff.val[constr_coeff.id[1]+1]*field_p->dof_vert[edge_nodes[0][1]].vec_dof_3[i] +
		  constr_coeff.val[constr_coeff.id[1]+2]*field_p->dof_edge[ned].vec_dof_3[i];
	      }
	    }
	  }
	}

	// -------------------- CREATE REAL DOF_EDGE -------------------- //
	if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
	  //Remove constrained data if it needed
	  SAFE_FREE(field_p->dof_edge[edge_sons[0]].constr);
	  SAFE_FREE(field_p->dof_edge[edge_sons[0]].constr_type);

	  SAFE_FREE(field_p->dof_edge[edge_sons[1]].constr);
	  SAFE_FREE(field_p->dof_edge[edge_sons[1]].constr_type);
	
	  //Alloc memory and calculate values for new node if it needed
	  
	  // - First edge:
	  if(field_p->dof_edge[edge_sons[0]].vec_dof_1 == NULL) {
	    //Alloc memory
	    field_p->dof_edge[edge_sons[0]].vec_dof_1 = (double *) malloc(nreq*sizeof(double));
	    for(i=0;i<nreq;i++) { field_p->dof_edge[edge_sons[0]].vec_dof_1[i]=0.0; }
	    if(field_p->nr_sol>1){
	      field_p->dof_edge[edge_sons[0]].vec_dof_2 = (double *) malloc(nreq*sizeof(double));
	      for(i=0;i<nreq;i++) { field_p->dof_edge[edge_sons[0]].vec_dof_2[i]=0.0; }
	    }
	    if(field_p->nr_sol>2){
	      field_p->dof_edge[edge_sons[0]].vec_dof_3 = (double *) malloc(nreq*sizeof(double));
	      for(i=0;i<nreq;i++) { field_p->dof_edge[edge_sons[0]].vec_dof_3[i]=0.0; }
	    }

	    field_p->dof_edge[edge_sons[0]].active = 1; //Set child status to active
	    
	    //Calculate new values
	    for(i=0;i<nreq;i++) {
	      field_p->dof_edge[edge_sons[0]].vec_dof_1[i] += 
		constr_coeff.val[constr_coeff.id[2]]*field_p->dof_edge[ned].vec_dof_1[i];
	      
	      if(field_p->nr_sol>1) {
		field_p->dof_edge[edge_sons[0]].vec_dof_2[i] += 
		  constr_coeff.val[constr_coeff.id[2]]*field_p->dof_edge[ned].vec_dof_2[i];
	      }
	      
	      if(field_p->nr_sol>2) {
		field_p->dof_edge[edge_sons[0]].vec_dof_3[i] += 
		  constr_coeff.val[constr_coeff.id[2]]*field_p->dof_edge[ned].vec_dof_3[i];
	      }
	    }
	  }

	  // - Second edge:
	  if(field_p->dof_edge[edge_sons[1]].vec_dof_1 == NULL) {
	    //Alloc memory
	    field_p->dof_edge[edge_sons[1]].vec_dof_1 = (double *) malloc(nreq*sizeof(double));
	    for(i=0;i<nreq;i++) { field_p->dof_edge[edge_sons[1]].vec_dof_1[i]=0.0; }
	    if(field_p->nr_sol>1){
	      field_p->dof_edge[edge_sons[1]].vec_dof_2 = (double *) malloc(nreq*sizeof(double));
	      for(i=0;i<nreq;i++) { field_p->dof_edge[edge_sons[1]].vec_dof_2[i]=0.0; }
	    }
	    if(field_p->nr_sol>2){
	      field_p->dof_edge[edge_sons[1]].vec_dof_3 = (double *) malloc(nreq*sizeof(double));
	      for(i=0;i<nreq;i++) { field_p->dof_edge[edge_sons[1]].vec_dof_3[i]=0.0; }
	    }

	    field_p->dof_edge[edge_sons[1]].active = 1; //Set child status to active

	    //Calculate new values
	    for(i=0;i<nreq;i++) {
	      field_p->dof_edge[edge_sons[1]].vec_dof_1[i] += 
		constr_coeff.val[constr_coeff.id[2]]*field_p->dof_edge[ned].vec_dof_1[i];
	      
	      if(field_p->nr_sol>1) {
		field_p->dof_edge[edge_sons[1]].vec_dof_2[i] += 
		  constr_coeff.val[constr_coeff.id[2]]*field_p->dof_edge[ned].vec_dof_2[i];
	      }
	      
	      if(field_p->nr_sol>2) {
		field_p->dof_edge[edge_sons[1]].vec_dof_3[i] += 
		  constr_coeff.val[constr_coeff.id[2]]*field_p->dof_edge[ned].vec_dof_3[i];
	      }
	    }
	  }
 
	  field_p->dof_edge[ned].active = 0; //Set parent status to inactive
	}
      } /* end if REAL node */
      else {
	
	/* CHECK PARENTS NODES */
	for(i=0; i<2; i++) {
	  if(field_p->dof_vert[edge_nodes[0][i]].vec_dof_1 == NULL && field_p->dof_vert[edge_nodes[0][i]].constr != NULL) {
	    ido = edge_nodes[0][i];
	    
	    //Alloc memory
	    field_p->dof_vert[ido].vec_dof_1 = (double *) malloc(nreq*sizeof(double));
	    for(i=0;i<nreq;i++) { field_p->dof_vert[ido].vec_dof_1[i]=0.0; }
	    if(field_p->nr_sol>1){
	      field_p->dof_vert[ido].vec_dof_2 = (double *) malloc(nreq*sizeof(double));
	      for(i=0;i<nreq;i++) { field_p->dof_vert[ido].vec_dof_2[i]=0.0; }
	    }
	    if(field_p->nr_sol>2){
	      field_p->dof_vert[ido].vec_dof_3 = (double *) malloc(nreq*sizeof(double));
	      for(i=0;i<nreq;i++) { field_p->dof_vert[ido].vec_dof_3[i]=0.0; }
	    }

	    field_p->dof_vert[ido].active = 1;
	    
	    //Calculate values
	    for(ino=1,j=0; ino<=field_p->dof_vert[ido].constr[0]; ino++,j++) {	  
	      switch(field_p->dof_vert[ido].constr_type[ino]) {
		
	      case APC_VERTEX:
		// -- DEBUG: constrained node points to an empty pointer 
		assert(field_p->dof_vert[field_p->dof_vert[ido].constr[ino]].vec_dof_1 != NULL);
		
		for(i=0;i<nreq;i++) {
		  field_p->dof_vert[ido].vec_dof_1[i] +=
		    constr_coeff.val[constr_coeff.id[field_p->dof_vert[ido].constr_type[0]]+j]*
		    field_p->dof_vert[field_p->dof_vert[ido].constr[ino]].vec_dof_1[i];
		}
		if(field_p->nr_sol>1){
		  for(i=0;i<nreq;i++) { 
		    field_p->dof_vert[ido].vec_dof_2[i] +=
		      constr_coeff.val[constr_coeff.id[field_p->dof_vert[ido].constr_type[0]]+j]*
		      field_p->dof_vert[field_p->dof_vert[ido].constr[ino]].vec_dof_2[i];
		  }
		}
		if(field_p->nr_sol>2){
		  for(i=0;i<nreq;i++) { 
		    field_p->dof_vert[ido].vec_dof_3[i] +=
		      constr_coeff.val[constr_coeff.id[field_p->dof_vert[ido].constr_type[0]]+j]*
		      field_p->dof_vert[field_p->dof_vert[ido].constr[ino]].vec_dof_3[i];
		  }
		}
		break;
		
		
	      case APC_EDGE:
		// -- DEBUG: constrained node points to an empty pointer 
		assert(field_p->dof_edge[field_p->dof_vert[ido].constr[ino]].vec_dof_1 != NULL);
		
		for(i=0;i<nreq;i++) {
		  field_p->dof_vert[ido].vec_dof_1[i] +=
		    constr_coeff.val[constr_coeff.id[field_p->dof_vert[ido].constr_type[0]]+j]*
		    field_p->dof_edge[field_p->dof_vert[ido].constr[ino]].vec_dof_1[i];
		  
		}
		if(field_p->nr_sol>1){
		  for(i=0;i<nreq;i++) { 
		    field_p->dof_vert[ido].vec_dof_2[i] +=
		      constr_coeff.val[constr_coeff.id[field_p->dof_vert[ido].constr_type[0]]+j]*
		      field_p->dof_edge[field_p->dof_vert[ido].constr[ino]].vec_dof_2[i];
		  }
		}
		if(field_p->nr_sol>2){
		  for(i=0;i<nreq;i++) { 
		    field_p->dof_vert[ido].vec_dof_3[i] +=
		      constr_coeff.val[constr_coeff.id[field_p->dof_vert[ido].constr_type[0]]+j]*
		      field_p->dof_edge[field_p->dof_vert[ido].constr[ino]].vec_dof_3[i];
		  }
		}
		break;
		
	      default:
		exit(-1);
	      };
	    }

	    //Free constr value
	    SAFE_FREE(field_p->dof_vert[ido].constr);
	    SAFE_FREE(field_p->dof_vert[ido].constr_type);	
	  }
	}

	if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
	  if(field_p->dof_edge[ned].vec_dof_1 == NULL && field_p->dof_edge[ned].constr != NULL ) {
	    //Alloc memory
	    field_p->dof_edge[ned].vec_dof_1 = (double *) malloc(nreq*sizeof(double));
	    for(i=0;i<nreq;i++) { field_p->dof_edge[ned].vec_dof_1[i]=0.0; }
	    if(field_p->nr_sol>1){
	      field_p->dof_edge[ned].vec_dof_2 = (double *) malloc(nreq*sizeof(double));
	      for(i=0;i<nreq;i++) { field_p->dof_edge[ned].vec_dof_2[i]=0.0; }
	    }
	    if(field_p->nr_sol>2){
	      field_p->dof_edge[ned].vec_dof_3 = (double *) malloc(nreq*sizeof(double));
	      for(i=0;i<nreq;i++) { field_p->dof_edge[ned].vec_dof_3[i]=0.0; }
	    }
	    
	    field_p->dof_edge[ned].active = 1;

	    //Calculate values
	    for(ino=1,j=0; ino<=field_p->dof_edge[ned].constr[0]; ino++,j++) {
	      
	      // -- DEBUG: constrained node points to an empty pointer 
	      assert(field_p->dof_edge[field_p->dof_edge[ned].constr[ino]].vec_dof_1 != NULL);
	      
	      for(i=0;i<nreq;i++) {
		field_p->dof_edge[ned].vec_dof_1[i] +=
		  constr_coeff.val[constr_coeff.id[field_p->dof_edge[ned].constr_type[0]]+j]*
		  field_p->dof_edge[field_p->dof_edge[ned].constr[ino]].vec_dof_1[i];
	      }
	      if(field_p->nr_sol>1){
		for(i=0;i<nreq;i++) { 
		  field_p->dof_edge[ned].vec_dof_2[i] +=
		    constr_coeff.val[constr_coeff.id[field_p->dof_edge[ned].constr_type[0]]+j]*
		    field_p->dof_edge[field_p->dof_edge[ned].constr[ino]].vec_dof_2[i];
		}
	      }
	      if(field_p->nr_sol>2){
		for(i=0;i<nreq;i++) { 
		  field_p->dof_edge[ned].vec_dof_3[i] +=
		    constr_coeff.val[constr_coeff.id[field_p->dof_edge[ned].constr_type[0]]+j]*
		    field_p->dof_edge[field_p->dof_edge[ned].constr[ino]].vec_dof_3[i];
		}
	      }
	    }
	    
	    //Free constr value
	    SAFE_FREE(field_p->dof_edge[ned].constr);
	    SAFE_FREE(field_p->dof_edge[ned].constr_type);
	  }
	}

	// -- DEBUG: ancestors are constrained
	assert(field_p->dof_vert[edge_nodes[0][0]].vec_dof_1 != NULL);
	assert(field_p->dof_vert[edge_nodes[0][1]].vec_dof_1 != NULL);
	if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
	  assert(field_p->dof_edge[ned].vec_dof_1 != NULL);
	}

	/* indicate that there are constrained nodes in the mesh */
	field_p->constr = APC_TRUE;
	
	
	// -------------------- CREATE CONSTR DOF_VERT -------------------- //
	
	//Remove real data if it needed
	SAFE_FREE(field_p->dof_vert[nvert].vec_dof_1);	  
	SAFE_FREE(field_p->dof_vert[nvert].vec_dof_2);
	SAFE_FREE(field_p->dof_vert[nvert].vec_dof_3);
	
	//Alloc memory and create constr new node if it needed
	if(field_p->dof_vert[nvert].constr == NULL) {
	  
	  if(field_p->pdeg == APC_LINEAR_APPROXIMATION_PDEG) { // Do it only for linear approximation
	    
	    //Alloc memory
	    field_p->dof_vert[nvert].constr = (int *) malloc(3*sizeof(int));
	    field_p->dof_vert[nvert].constr_type = (int *) malloc(3*sizeof(int));
	    field_p->dof_vert[nvert].constr[0] = 2;
	    field_p->dof_vert[nvert].constr_type[0] = 1;
	    
	    //Create constr data
	    field_p->dof_vert[nvert].constr[1] = edge_nodes[0][0];
	    field_p->dof_vert[nvert].constr[2] = edge_nodes[0][1];
	    field_p->dof_vert[nvert].constr_type[1] = APC_VERTEX;
	    field_p->dof_vert[nvert].constr_type[2] = APC_VERTEX;
	    
	  }

	  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
	    
	    //Alloc memory
	    field_p->dof_vert[nvert].constr = (int *) malloc(4*sizeof(int));
	    field_p->dof_vert[nvert].constr_type = (int *) malloc(4*sizeof(int));
	    field_p->dof_vert[nvert].constr[0] = 3;
	    field_p->dof_vert[nvert].constr_type[0] = 1;
	    
	    //Create constr data
	    field_p->dof_vert[nvert].constr[1] = edge_nodes[0][0];
	    field_p->dof_vert[nvert].constr[2] = edge_nodes[0][1];
	    field_p->dof_vert[nvert].constr[3] = ned;
	    field_p->dof_vert[nvert].constr_type[1] = APC_VERTEX;
	    field_p->dof_vert[nvert].constr_type[2] = APC_VERTEX;
	    field_p->dof_vert[nvert].constr_type[3] = APC_EDGE;
	  }
	}
	field_p->dof_vert[nvert].active = 1;
	
	// -------------------- CREATE CONSTR DOF_EDGE -------------------- //
	
	if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
	  
	  //Remove real data if it needed
	  SAFE_FREE(field_p->dof_edge[edge_sons[0]].vec_dof_1);
	  SAFE_FREE(field_p->dof_edge[edge_sons[0]].vec_dof_2);
	  SAFE_FREE(field_p->dof_edge[edge_sons[0]].vec_dof_3);

	  SAFE_FREE(field_p->dof_edge[edge_sons[1]].vec_dof_1);
	  SAFE_FREE(field_p->dof_edge[edge_sons[1]].vec_dof_2);
	  SAFE_FREE(field_p->dof_edge[edge_sons[1]].vec_dof_3);
	  
	  //Alloc memory and create constr new node if it needed
	  // - First edge:
	  if(field_p->dof_edge[edge_sons[0]].constr == NULL) {
	    
	    //Alloc memory
	    field_p->dof_edge[edge_sons[0]].constr = (int *) malloc(2*sizeof(int));
	    field_p->dof_edge[edge_sons[0]].constr_type = (int *) malloc(2*sizeof(int));

	    field_p->dof_edge[edge_sons[0]].active = 1;
	    
	    //Create constr data
	    field_p->dof_edge[edge_sons[0]].constr[0] = 1;
	    field_p->dof_edge[edge_sons[0]].constr_type[0] = 2;
	    
	    field_p->dof_edge[edge_sons[0]].constr[1] = ned;
	    field_p->dof_edge[edge_sons[0]].constr_type[1] = APC_EDGE;
	  }
	  
	  // - Second edge:
	  if(field_p->dof_edge[edge_sons[1]].constr == NULL) {
	    
	    //Alloc memory
	    field_p->dof_edge[edge_sons[1]].constr = (int *) malloc(2*sizeof(int));
	    field_p->dof_edge[edge_sons[1]].constr_type = (int *) malloc(2*sizeof(int));

	    field_p->dof_edge[edge_sons[1]].active = 1;
	    
	    //Create constr data
	    field_p->dof_edge[edge_sons[1]].constr[0] = 1;
	    field_p->dof_edge[edge_sons[1]].constr_type[0] = 2;
	    
	    field_p->dof_edge[edge_sons[1]].constr[1] = ned;
	    field_p->dof_edge[edge_sons[1]].constr_type[1] = APC_EDGE;
	  }

	  field_p->dof_edge[ned].active = 1;
	}
      } /* end if CONSTR node, after checking the current status */
    } /* end if edge inactive */
  } /* end loop over edges */

  /* for new mid-face nodes */
  for(nfa=1; nfa<=max_face_id; nfa++) {

    //Detect inactive faces
    if(mmr_fa_type(mesh_id,nfa) == MMC_QUAD && mmr_fa_status(mesh_id,nfa) == MMC_INACTIVE) {

      // QUADANGLE FACE //
      /*
       *    o* * * * x * * * *o        o* (c) * C * (c) *o    Note: If new vertices in divided quadangle face are 
       *    *        +        +        *        +        +    CONSTRAINT then new edges in this face are CONSTRAINTS too.
       *    *        +        +       (c)      (c)      (c)
       *    *        +        +        *        +        + 
       *    *        +        +        *        +        +  
       *    x+ + + + x + + + +x  <==   C+ (c) + C + (c) +C
       *    *        +        +        *        +        + 
       *    *        +        +        *        +        +
       *    *        +        +       (c)      (c)      (c) 
       *    *        +        +        *        +        +     
       *    o* * * * x * * * *o        o* (c) * C * (c) *o  
       *
       */
      
      //Get face information
      mmr_fa_fam(mesh_id, nfa, fa_sons, &mid_node);
      
      /* for divided quadrilateral faces */
      if(mid_node>0 && mmr_node_status(mesh_id, mid_node)==MMC_ACTIVE) {
	  
	/* active node in the middle */
	nvert = mid_node;

	/* check the status of the node */
	mmr_fa_neig(mesh_id,nfa,face_neig,neig_sides,&node_shift,NULL,NULL,NULL);

	/* if face between two elements and only one divided (inactive) */      
	if(face_neig[0]!=0 && face_neig[1]!=0 &&
	   ((mmr_el_status(mesh_id, face_neig[0])==MMC_ACTIVE && mmr_el_status(mesh_id, face_neig[1])==MMC_INACTIVE) ||
	    (mmr_el_status(mesh_id, face_neig[0])==MMC_INACTIVE && mmr_el_status(mesh_id, face_neig[1])==MMC_ACTIVE))
	   ) {
	  
	  /* indicate that there are constrained nodes in the mesh */
	  field_p->constr = APC_TRUE;

	  //Get face nodes and edges
	  num_nodes[0] = mmr_fa_node_coor(mesh_id, nfa, quad_fa_nodes[0], NULL);
	  
	  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
	    num_edges[0] = mmr_fa_edges(mesh_id, nfa, quad_fa_edges[0], NULL);
	  }
	    
	  // -- DEBUG: ancestors are constrained
	  assert(field_p->dof_vert[quad_fa_nodes[0][1]].vec_dof_1 != NULL);
	  assert(field_p->dof_vert[quad_fa_nodes[0][2]].vec_dof_1 != NULL);
	  assert(field_p->dof_vert[quad_fa_nodes[0][3]].vec_dof_1 != NULL);
	  assert(field_p->dof_vert[quad_fa_nodes[0][4]].vec_dof_1 != NULL);
	  
	  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
	    assert(field_p->dof_edge[quad_fa_edges[0][0]].vec_dof_1 != NULL);
	    assert(field_p->dof_edge[quad_fa_edges[0][1]].vec_dof_1 != NULL);
	    assert(field_p->dof_edge[quad_fa_edges[0][2]].vec_dof_1 != NULL);
	    assert(field_p->dof_edge[quad_fa_edges[0][3]].vec_dof_1 != NULL);
	  }

	  if(field_p->dof_vert[nvert].constr != NULL){
	    /* currently CONSTR node, previously CONSTR node - do nothing */
	  }
	  else {
	    /* currently CONSTR node, previously REAL node  */
	    SAFE_FREE(field_p->dof_vert[nvert].vec_dof_1);
	    SAFE_FREE(field_p->dof_vert[nvert].vec_dof_2);
	    SAFE_FREE(field_p->dof_vert[nvert].vec_dof_3);

	    field_p->dof_vert[nvert].active = 1; //Set vertex status to active

	    /* Alloc memory for new contrained vertex node */
	    if(field_p->pdeg == APC_LINEAR_APPROXIMATION_PDEG) { // Do it only for linear approximation
	      field_p->dof_vert[nvert].constr = (int*) malloc(5*sizeof(int));
	      field_p->dof_vert[nvert].constr_type = (int*) malloc(5*sizeof(int));	    
	      field_p->dof_vert[nvert].constr[0] = num_nodes[0];
	      field_p->dof_vert[nvert].constr_type[0] = 3;
	    }
	    
	    if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
	      field_p->dof_vert[nvert].constr = (int*) malloc(9*sizeof(int));
	      field_p->dof_vert[nvert].constr_type = (int*) malloc(9*sizeof(int));	    
	      field_p->dof_vert[nvert].constr[0] = num_nodes[0] + num_edges[0];
	      field_p->dof_vert[nvert].constr_type[0] = 3;
	    }

	    ido = 1;
	    //Loop over vertices
	    for(ino=1; ino<=num_nodes[0]; ino++,ido++) {
	      field_p->dof_vert[nvert].constr[ido] = quad_fa_nodes[0][ino];
	      field_p->dof_vert[nvert].constr_type[ido] = APC_VERTEX;
	    }

	    if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
	      //Loop over edges
	      for(ied=0; ied<num_edges[0]; ied++,ido++) {
		field_p->dof_vert[nvert].constr[ido] = quad_fa_edges[0][ied];
		field_p->dof_vert[nvert].constr_type[ido] = APC_EDGE;
	      }
	    }

	    /* Add CONSTRAINT EDGES */
	    /*
	     *    o* * * * x * * * *o
	     *    *        +        +
	     *    * F[3]   + F[4]   +
	     *    *       (2)check  + 
	     *    *        +        +  
	     *    x+ +(2)+ C +(1)+ +x
	     *    *        +        +
	     *    * F[1]   + F[2]   +
	     *    * check (1)       + 
	     *    *        +        +    
	     *    o* * * * x * * * *o  
	     *
	     */

	    if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
	      
	      //List of edges to modify and list of constraints coeff types
	      int list_of_edges[5] = { 0, -1, -1, -1, -1 };
	    
	      // Add constraint EDGE to subface 1
	      {
		//Get face nodes and edges
		num_edges[1] = mmr_fa_edges(mesh_id, fa_sons[1], quad_fa_edges[1], NULL);
		
		//First edge
		nedge = quad_fa_edges[1][1];
		list_of_edges[++list_of_edges[0]] = nedge;

		//Second edge
		nedge = quad_fa_edges[1][2];
		list_of_edges[++list_of_edges[0]] = nedge;
	      }

	      // Add constraint EDGE to subface 3
	      {
		//Get face nodes and edges
		num_edges[1] = mmr_fa_edges(mesh_id, fa_sons[3], quad_fa_edges[1], NULL);
	      
		//First edge
		nedge = quad_fa_edges[1][0];
		list_of_edges[++list_of_edges[0]] = nedge;

		//Second edge
		nedge = quad_fa_edges[1][3];
		list_of_edges[++list_of_edges[0]] = nedge;
	      }

	      //Create constraint EDGE
	      for(i=1; i<=list_of_edges[0]; i++) {
		nedge = list_of_edges[i];

		field_p->dof_edge[nedge].active = 1; //Set edge status to active
	      
		//Alloc memory for new contrained edge node
		field_p->dof_edge[nedge].constr = (int*) malloc(3*sizeof(int));
		field_p->dof_edge[nedge].constr_type = (int*) malloc(3*sizeof(int));
	      
		field_p->dof_edge[nedge].constr[0] = 2;
		field_p->dof_edge[nedge].constr_type[0] = 4;

		//Dealloc memory for real edge
		SAFE_FREE(field_p->dof_edge[nedge].vec_dof_1);
		SAFE_FREE(field_p->dof_edge[nedge].vec_dof_2);
		SAFE_FREE(field_p->dof_edge[nedge].vec_dof_3);	
	      }

	      //First edge
	      field_p->dof_edge[list_of_edges[1]].constr[1] = quad_fa_edges[0][1]; field_p->dof_edge[list_of_edges[1]].constr_type[1] = APC_EDGE;
	      field_p->dof_edge[list_of_edges[1]].constr[2] = quad_fa_edges[0][3]; field_p->dof_edge[list_of_edges[1]].constr_type[2] = APC_EDGE;

	      //Second edge
	      field_p->dof_edge[list_of_edges[2]].constr[1] = quad_fa_edges[0][0]; field_p->dof_edge[list_of_edges[2]].constr_type[1] = APC_EDGE;
	      field_p->dof_edge[list_of_edges[2]].constr[2] = quad_fa_edges[0][2]; field_p->dof_edge[list_of_edges[2]].constr_type[2] = APC_EDGE;

	      //Third edge
	      field_p->dof_edge[list_of_edges[3]].constr[1] = quad_fa_edges[0][0]; field_p->dof_edge[list_of_edges[3]].constr_type[1] = APC_EDGE;
	      field_p->dof_edge[list_of_edges[3]].constr[2] = quad_fa_edges[0][2]; field_p->dof_edge[list_of_edges[3]].constr_type[2] = APC_EDGE;
	      
	      //Fourth edge
	      field_p->dof_edge[list_of_edges[4]].constr[1] = quad_fa_edges[0][1]; field_p->dof_edge[list_of_edges[4]].constr_type[1] = APC_EDGE;
	      field_p->dof_edge[list_of_edges[4]].constr[2] = quad_fa_edges[0][3]; field_p->dof_edge[list_of_edges[4]].constr_type[2] = APC_EDGE;
	    }
	  } /* end if CONSTR node, previously not CONSTR node */
	} /* end if currently constrained active node */
	else {

	  //Get face nodes and edges
	  num_nodes[0] = mmr_fa_node_coor(mesh_id, nfa, quad_fa_nodes[0], NULL);

	  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
	    num_edges[0] = mmr_fa_edges(mesh_id, nfa, quad_fa_edges[0], NULL);
	  }
	  
	  // -- DEBUG: ancestors are constrained
	  assert(field_p->dof_vert[quad_fa_nodes[0][1]].vec_dof_1 != NULL);
	  assert(field_p->dof_vert[quad_fa_nodes[0][2]].vec_dof_1 != NULL);
	  assert(field_p->dof_vert[quad_fa_nodes[0][3]].vec_dof_1 != NULL);
	  assert(field_p->dof_vert[quad_fa_nodes[0][4]].vec_dof_1 != NULL);

	  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
	    assert(field_p->dof_edge[quad_fa_edges[0][0]].vec_dof_1 != NULL);
	    assert(field_p->dof_edge[quad_fa_edges[0][1]].vec_dof_1 != NULL);
	    assert(field_p->dof_edge[quad_fa_edges[0][2]].vec_dof_1 != NULL);
	    assert(field_p->dof_edge[quad_fa_edges[0][3]].vec_dof_1 != NULL);
	  }
	  
	  /* free constr table */
	  SAFE_FREE(field_p->dof_vert[nvert].constr);
	  SAFE_FREE(field_p->dof_vert[nvert].constr_type);

	  /* create real node */
	  if(field_p->dof_vert[nvert].vec_dof_1 != NULL) {
	    /* currently REAL node, previously REAL node  - do nothing */
	  }
	  else {

	    /* Create list of edges
	     *    o* * * * x * * * *o
	     *    *        +        +
	     *    * F[3]   + F[4]   +
	     *    *       (2)check  + 
	     *    *        +        +  
	     *    x+ +(2)+ R +(1)+ +x
	     *    *        +        +
	     *    * F[1]   + F[2]   +
	     *    * check (1)       + 
	     *    *        +        +    
	     *    o* * * * x * * * *o  
	     *
	     */
	    
	    //List of edges to modify
	    int list_of_edges[5] = { 0, -1, -1, -1, -1 };

	    if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
	      // Add real EDGE to subface 1
	      {
		//Get face nodes and edges
		num_edges[1] = mmr_fa_edges(mesh_id, fa_sons[1], quad_fa_edges[1], NULL);

		//First edge
		nedge = quad_fa_edges[1][1];
		list_of_edges[++list_of_edges[0]] = nedge;

		//Second edge
		nedge = quad_fa_edges[1][2];
		list_of_edges[++list_of_edges[0]] = nedge;
	      }

	      // Add real EDGE to subface 3
	      {
		//Get face nodes and edges
		num_edges[1] = mmr_fa_edges(mesh_id, fa_sons[3], quad_fa_edges[1], NULL);
	      
		//First edge
		nedge = quad_fa_edges[1][0];
		list_of_edges[++list_of_edges[0]] = nedge;

		//Second edge
		nedge = quad_fa_edges[1][3];
		list_of_edges[++list_of_edges[0]] = nedge;
	      }
	    }

	    /* Create real VERTICES and EDGES */
	    
	    //Operation on VERTICES
	    field_p->dof_vert[nvert].vec_dof_1 = (double *) malloc(nreq*sizeof(double));
	    for(i=0;i<nreq;i++) { field_p->dof_vert[nvert].vec_dof_1[i] = 0.0; }
	    if(field_p->nr_sol>1){
	      field_p->dof_vert[nvert].vec_dof_2 = (double *) malloc(nreq*sizeof(double));
	      for(i=0;i<nreq;i++) { field_p->dof_vert[nvert].vec_dof_2[i] = 0.0; }
	    }
	    if(field_p->nr_sol>2){
	      field_p->dof_vert[nvert].vec_dof_3 = (double *) malloc(nreq*sizeof(double));
	      for(i=0;i<nreq;i++) { field_p->dof_vert[nvert].vec_dof_3[i] = 0.0; }
	    }
	    field_p->dof_vert[nvert].active = 1; //Set vertex status to active

	    //Operation on Edges
	    if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
	      for(i=1; i<=list_of_edges[0]; i++) {
		nedge = list_of_edges[i];

		field_p->dof_edge[nedge].active = 1; //Set edge status to active

		//Alloc memory for real edge
		if(field_p->dof_edge[nedge].vec_dof_1 == NULL) {
		  field_p->dof_edge[nedge].vec_dof_1 = (double *) malloc(nreq*sizeof(double));
		  for(j=0;j<nreq;j++) { field_p->dof_edge[nedge].vec_dof_1[j]=0.0; }

		  if(field_p->nr_sol>1){
		    field_p->dof_edge[nedge].vec_dof_2 = (double *) malloc(nreq*sizeof(double));
		    for(j=0;j<nreq;j++) { field_p->dof_edge[nedge].vec_dof_2[j]=0.0; }
		  }
		  if(field_p->nr_sol>2){
		    field_p->dof_edge[nedge].vec_dof_3 = (double *) malloc(nreq*sizeof(double));
		    for(j=0;j<nreq;j++) { field_p->dof_edge[nedge].vec_dof_3[j]=0.0; }
		  }		  
		}

		//Dealloc memory for constrained edge
		SAFE_FREE(field_p->dof_edge[nedge].constr);
		SAFE_FREE(field_p->dof_edge[nedge].constr_type);
	      }
	    }
	 
	    //Calculate value for new vertices
	    for(i=0;i<nreq;i++){
	      for(ino=1,j=constr_coeff.id[3]; ino<=num_nodes[0]; ino++,j++){
		field_p->dof_vert[nvert].vec_dof_1[i] += constr_coeff.val[j]*field_p->dof_vert[quad_fa_nodes[0][ino]].vec_dof_1[i];
		if(field_p->nr_sol>1){
		  field_p->dof_vert[nvert].vec_dof_2[i] += constr_coeff.val[j]*field_p->dof_vert[quad_fa_nodes[0][ino]].vec_dof_2[i];
		}
		if(field_p->nr_sol>2){
		  field_p->dof_vert[nvert].vec_dof_3[i] += constr_coeff.val[j]*field_p->dof_vert[quad_fa_nodes[0][ino]].vec_dof_3[i];
		}
	      }
	      
	      if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
		for(ied=0,j=(constr_coeff.id[3]+4); ied<num_edges[0]; ied++,j++) {
		  if(field_p->dof_edge[quad_fa_edges[0][ied]].vec_dof_1 != NULL) {
		    field_p->dof_vert[nvert].vec_dof_1[i] += constr_coeff.val[j]*field_p->dof_edge[quad_fa_edges[0][ied]].vec_dof_1[i];
		    if(field_p->nr_sol>1) {
		      field_p->dof_vert[nvert].vec_dof_2[i] += constr_coeff.val[j]*field_p->dof_edge[quad_fa_edges[0][ied]].vec_dof_2[i];
		    }
		    if(field_p->nr_sol>2) {
		      field_p->dof_vert[nvert].vec_dof_3[i] += constr_coeff.val[j]*field_p->dof_edge[quad_fa_edges[0][ied]].vec_dof_3[i];
		    }
		  }
		}
	      }
	    }

	    if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation

	      // -- DEBUG: element does not existed
	      assert(list_of_edges[1] != -1);
	      assert(list_of_edges[2] != -1);
	      assert(list_of_edges[3] != -1);
	      assert(list_of_edges[4] != -1);

	      //Calculate value for new edges
	      for(i=0;i<nreq;i++) {
		//Calcualte values of first edge
		field_p->dof_edge[list_of_edges[1]].vec_dof_1[i] +=
		  constr_coeff.val[constr_coeff.id[4]+0]*field_p->dof_edge[quad_fa_edges[0][1]].vec_dof_1[i] +
		  constr_coeff.val[constr_coeff.id[4]+1]*field_p->dof_edge[quad_fa_edges[0][3]].vec_dof_1[i];
		
		if(field_p->nr_sol>1) {
		  field_p->dof_edge[list_of_edges[1]].vec_dof_2[i] +=
		    constr_coeff.val[constr_coeff.id[4]+0]*field_p->dof_edge[quad_fa_edges[0][1]].vec_dof_2[i] +
		    constr_coeff.val[constr_coeff.id[4]+1]*field_p->dof_edge[quad_fa_edges[0][3]].vec_dof_2[i];
		}
	      
		if(field_p->nr_sol>2) {
		  field_p->dof_edge[list_of_edges[1]].vec_dof_3[i] +=
		    constr_coeff.val[constr_coeff.id[4]+0]*field_p->dof_edge[quad_fa_edges[0][1]].vec_dof_3[i] +
		    constr_coeff.val[constr_coeff.id[4]+1]*field_p->dof_edge[quad_fa_edges[0][3]].vec_dof_3[i];
		}
		
		
		//Calcualte values of second edge
		field_p->dof_edge[list_of_edges[2]].vec_dof_1[i] +=
		  constr_coeff.val[constr_coeff.id[4]+0]*field_p->dof_edge[quad_fa_edges[0][0]].vec_dof_1[i] +
		  constr_coeff.val[constr_coeff.id[4]+1]*field_p->dof_edge[quad_fa_edges[0][2]].vec_dof_1[i];
		
		if(field_p->nr_sol>1) {
		  field_p->dof_edge[list_of_edges[2]].vec_dof_2[i] +=
		    constr_coeff.val[constr_coeff.id[4]+0]*field_p->dof_edge[quad_fa_edges[0][0]].vec_dof_2[i] +
		    constr_coeff.val[constr_coeff.id[4]+1]*field_p->dof_edge[quad_fa_edges[0][2]].vec_dof_2[i];
		}
		
		if(field_p->nr_sol>2) {
		  field_p->dof_edge[list_of_edges[2]].vec_dof_3[i] +=
		    constr_coeff.val[constr_coeff.id[4]+0]*field_p->dof_edge[quad_fa_edges[0][0]].vec_dof_3[i] +
		    constr_coeff.val[constr_coeff.id[4]+1]*field_p->dof_edge[quad_fa_edges[0][2]].vec_dof_3[i];
		}

		
		//Calcualte values of third edge
		field_p->dof_edge[list_of_edges[3]].vec_dof_1[i] +=
		  constr_coeff.val[constr_coeff.id[4]+0]*field_p->dof_edge[quad_fa_edges[0][0]].vec_dof_1[i] +
		  constr_coeff.val[constr_coeff.id[4]+1]*field_p->dof_edge[quad_fa_edges[0][2]].vec_dof_1[i];
		
		if(field_p->nr_sol>1) {
		  field_p->dof_edge[list_of_edges[3]].vec_dof_2[i] +=
		    constr_coeff.val[constr_coeff.id[4]+0]*field_p->dof_edge[quad_fa_edges[0][0]].vec_dof_2[i] +
		    constr_coeff.val[constr_coeff.id[4]+1]*field_p->dof_edge[quad_fa_edges[0][2]].vec_dof_2[i];
		}
		
		if(field_p->nr_sol>2) {
		  field_p->dof_edge[list_of_edges[3]].vec_dof_3[i] +=
		    constr_coeff.val[constr_coeff.id[4]+0]*field_p->dof_edge[quad_fa_edges[0][0]].vec_dof_3[i] +
		    constr_coeff.val[constr_coeff.id[4]+1]*field_p->dof_edge[quad_fa_edges[0][2]].vec_dof_3[i];
		}
		
		//Calculate values of fourth edge
		field_p->dof_edge[list_of_edges[4]].vec_dof_1[i] +=
		  constr_coeff.val[constr_coeff.id[4]+0]*field_p->dof_edge[quad_fa_edges[0][1]].vec_dof_1[i] +
		  constr_coeff.val[constr_coeff.id[4]+1]*field_p->dof_edge[quad_fa_edges[0][3]].vec_dof_1[i];
		
		if(field_p->nr_sol>1) {
		  field_p->dof_edge[list_of_edges[4]].vec_dof_2[i] +=
		    constr_coeff.val[constr_coeff.id[4]+0]*field_p->dof_edge[quad_fa_edges[0][1]].vec_dof_2[i] +
		    constr_coeff.val[constr_coeff.id[4]+1]*field_p->dof_edge[quad_fa_edges[0][3]].vec_dof_2[i];
		}

		if(field_p->nr_sol>2) {
		  field_p->dof_edge[list_of_edges[4]].vec_dof_3[i] +=
		    constr_coeff.val[constr_coeff.id[4]+0]*field_p->dof_edge[quad_fa_edges[0][1]].vec_dof_3[i] +
		    constr_coeff.val[constr_coeff.id[4]+1]*field_p->dof_edge[quad_fa_edges[0][3]].vec_dof_3[i];
		}
	      }
	    }
	  } /* end new REAL node previously not REAL node */
	} /* end REAL node */
      } /* end active mid-edge node */
    } /* end faces */
  } /* end loop over faces */

  /* check mid triangle faces */
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
    for(nfa=1; nfa<=max_face_id; nfa++) {

      //Detect inactive faces
      if((mmr_fa_status(mesh_id,nfa)) == MMC_INACTIVE && mmr_fa_type(mesh_id,nfa) == MMC_TRIA) {

      // TRIANGLE FACE //
      /*
       *    o*                            C +(c)+ C      R +(c)+ C       R +(c)+ C
       *    *  *                            +     +        +     +         +     +
       *    *    *                           (c) (c)        (c) (c)         (r) (c)
       *    *      *                            + +            + +             + +
       *    x+ + + + x            <==             C              C               R
       *    *  +     +  *                               
       *    *    +   +    *               R +(c)+ R        Note: Edge in middle triangle face become CONSTRAINTS when one
       *    *      + +      *               +     +        or all vertices on the ends of edge are CONSTRAINT. In the
       *    o* * * * x * * * *o              (r) (r)       otherwise edge is REAL.
       *                                        + +
       *                                          R
       */
      
        //Is constraint
	int is_constraint[3] = { 0,0,0 };

	//Get face information
	mmr_fa_fam(mesh_id, nfa, fa_sons, NULL);

	//Get parent vertices and edges
	num_edges[0] = mmr_fa_edges(mesh_id, nfa, tria_fa_edges[0], NULL);

	//Set center triagle
	mid_tria = fa_sons[4];

	if((mmr_fa_status(mesh_id,mid_tria)) == MMC_ACTIVE) {
	  
	  //Get triangle edges
	  num_edges[1] = mmr_fa_edges(mesh_id, mid_tria, tria_fa_edges[1], NULL);
	
	  // Look for constraint coefficient
	  for(ned=0; ned<num_edges[1]; ned++) {
	    nedge = tria_fa_edges[1][ned];
	    
	    //Get edge nodes
	    mmr_edge_nodes(mesh_id,nedge,edge_nodes[0]);
	    
	    //Check edge ends node status
	    if(field_p->dof_vert[edge_nodes[0][0]].constr != NULL && field_p->dof_vert[edge_nodes[0][1]].constr != NULL) {

	      // indicate that there are constrained nodes in the mesh
	      field_p->constr = APC_TRUE;

	      //Create contstraint
	      SAFE_FREE(field_p->dof_edge[nedge].vec_dof_1);
	      SAFE_FREE(field_p->dof_edge[nedge].vec_dof_2);
	      SAFE_FREE(field_p->dof_edge[nedge].vec_dof_3);

	      field_p->dof_edge[nedge].active = 1; //Set edge status to active
	      
	      if(field_p->dof_edge[nedge].constr == NULL) {
		field_p->dof_edge[nedge].constr = (int*) malloc(2*sizeof(int));
		field_p->dof_edge[nedge].constr_type = (int*) malloc(2*sizeof(int));
		field_p->dof_edge[nedge].constr[0] = 1;
		field_p->dof_edge[nedge].constr_type[0] = 2;
	      }
	      
	      field_p->dof_edge[nedge].constr[1] = tria_fa_edges[0][ned];
	      field_p->dof_edge[nedge].constr_type[1] = APC_EDGE;	      

	      // Is constraint face
	      is_constraint[ned] = 1;
	    } else {

	      //Create real edge
	      if(field_p->dof_edge[nedge].vec_dof_1 == NULL) {
		
		//Dealloc memory for constrained edge
		SAFE_FREE(field_p->dof_edge[nedge].constr);
		SAFE_FREE(field_p->dof_edge[nedge].constr_type);

		//Alloc memory for new edges
		field_p->dof_edge[nedge].vec_dof_1 = (double *) malloc(nreq*sizeof(double));
		for(i=0;i<nreq;i++) { field_p->dof_edge[nedge].vec_dof_1[i]=0.0; }
		if(field_p->nr_sol>1){
		  field_p->dof_edge[nedge].vec_dof_2 = (double *) malloc(nreq*sizeof(double));
		  for(i=0;i<nreq;i++) { field_p->dof_edge[nedge].vec_dof_2[i]=0.0; }
		}
		if(field_p->nr_sol>2){
		  field_p->dof_edge[nedge].vec_dof_3 = (double *) malloc(nreq*sizeof(double));
		  for(i=0;i<nreq;i++) { field_p->dof_edge[nedge].vec_dof_3[i]=0.0; }
		}
	      }
	      
	      field_p->dof_edge[nedge].constr[1] = tria_fa_edges[0][ned];
	      field_p->dof_edge[nedge].constr_type[1] = APC_EDGE; 

	      field_p->dof_edge[nedge].active = 1; //Set edge status to active
	    }
	  }
	
	  // Refine od derefine edges in mid triangle
	  for(ned=0; ned<num_edges[1]; ned++) {
	    nedge = tria_fa_edges[1][ned];

	    if(is_constraint[ned] == 0) {

	      // -- DEBUG: null pointer
	      assert(field_p->dof_edge[tria_fa_edges[0][ned]].vec_dof_1 != NULL);

	      //Calculate value for new edge
	      for(i=0;i<nreq;i++) {
		field_p->dof_edge[nedge].vec_dof_1[i] += constr_coeff.val[constr_coeff.id[2]]*field_p->dof_edge[tria_fa_edges[0][ned]].vec_dof_1[i];
		
		if(field_p->nr_sol>1) {
		  field_p->dof_edge[nedge].vec_dof_2[i] += constr_coeff.val[constr_coeff.id[2]]*field_p->dof_edge[tria_fa_edges[0][ned]].vec_dof_2[i];
		}
		if(field_p->nr_sol>2) {
		  field_p->dof_edge[nedge].vec_dof_3[i] += constr_coeff.val[constr_coeff.id[2]]*field_p->dof_edge[tria_fa_edges[0][ned]].vec_dof_3[i];
		}
	      }
	    }
	  }
	}
      }
    }

    /* looking mid triangle faces with no data */
    for(nfa=1; nfa<=max_face_id; nfa++) {

      //Detect inactive faces
      if(mmr_fa_status(mesh_id,nfa) == MMC_ACTIVE || mmr_fa_status(mesh_id,nfa) == MMC_INACTIVE) {
	int face_neig[2] = {0,0};
	int elem_neig[6] = {0,0,0,0,0,0};
	int elem_faces[4][6] = {{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};
	int elem_face_id = 0;

	//Get edges
	num_edges[0] = mmr_fa_edges(mesh_id, nfa, tria_fa_edges[0], NULL);
      
	//Face with no data
	if(field_p->dof_edge[tria_fa_edges[0][0]].vec_dof_1 == NULL && 
	   field_p->dof_edge[tria_fa_edges[0][1]].vec_dof_1 == NULL &&
	   field_p->dof_edge[tria_fa_edges[0][2]].vec_dof_1 == NULL &&
	   field_p->dof_edge[tria_fa_edges[0][0]].constr == NULL &&
	   field_p->dof_edge[tria_fa_edges[0][1]].constr == NULL &&
	   field_p->dof_edge[tria_fa_edges[0][2]].constr == NULL) {

	  //Get face  neighbours
	  mmr_fa_eq_neig(mesh_id,nfa,face_neig,NULL,NULL);
	
	  //Find face id in element
	  mmr_el_faces(mesh_id,face_neig[0],elem_faces[0],NULL);
	  for(i=1; i<elem_faces[0][0]; i++) {
	    if(elem_faces[0][i] == nfa) {
	      elem_face_id = i;
	      break;
	    }
	  }

	  //Get element neighbour
	  mmr_el_eq_neig(mesh_id,face_neig[0],elem_neig,NULL);

	  //Get faces of neighbour elements 
	  mmr_el_faces(mesh_id,elem_neig[3],elem_faces[1],NULL);
	  mmr_el_faces(mesh_id,elem_neig[4],elem_faces[2],NULL);
	  mmr_el_faces(mesh_id,elem_neig[5],elem_faces[3],NULL);
	  
	  //Get face edges
	  num_edges[1] = mmr_fa_edges(mesh_id, elem_faces[1][elem_face_id], tria_fa_edges[1], NULL);
	  num_edges[2] = mmr_fa_edges(mesh_id, elem_faces[2][elem_face_id], tria_fa_edges[2], NULL);
	  num_edges[3] = mmr_fa_edges(mesh_id, elem_faces[3][elem_face_id], tria_fa_edges[3], NULL);
	  
	  //Dependence inside tringle face
	  // Mid triangle face F4
	  // F4->E1 = 0.5 * ( F2->E1 + F3->E1 )
	  // F4->E2 = 0.5 * ( F1->E2 + F3->E2 )
	  // F4->E3 = 0.5 * ( F1->E3 + F2->E3 )
	  int used_faces[3][2] = {{2,3},{1,3},{1,2}};
	  double values[2];

	  //Calculate values of edges
	  for(ned=0; ned<num_edges[0]; ned++) {
	    nedge = tria_fa_edges[0][ned];
	    field_p->dof_edge[nedge].active = 1;

	    //Alloc memory for real edge and calculate values
	    field_p->dof_edge[nedge].vec_dof_1 = (double *) malloc(nreq*sizeof(double));
	    for(i=0;i<nreq;i++) { 

	      //Edge 1
	      if(field_p->dof_edge[tria_fa_edges[used_faces[ned][0]][ned]].constr != NULL) {
		values[0] = 0.0;
		for(j=1; j<field_p->dof_edge[tria_fa_edges[used_faces[ned][0]][ned]].constr[0]; j++) {
		  values[0] += 
		    constr_coeff.val[constr_coeff.id[field_p->dof_edge[tria_fa_edges[used_faces[ned][0]][ned]].constr_type[0]]+j]*
		    field_p->dof_edge[field_p->dof_edge[tria_fa_edges[used_faces[ned][0]][ned]].constr[j]].vec_dof_1[i];
		}
	      } else {
		values[0] = field_p->dof_edge[tria_fa_edges[used_faces[ned][0]][ned]].vec_dof_1[i];
	      }

	      //Edge 2
	      if(field_p->dof_edge[tria_fa_edges[used_faces[ned][1]][ned]].constr != NULL) {
		values[1] = 0.0;
		for(j=1; j<field_p->dof_edge[tria_fa_edges[used_faces[ned][1]][ned]].constr[0]; j++) {
		  values[1] +=
		    constr_coeff.val[constr_coeff.id[field_p->dof_edge[tria_fa_edges[used_faces[ned][1]][ned]].constr_type[0]]+j]*
		    field_p->dof_edge[field_p->dof_edge[tria_fa_edges[used_faces[ned][1]][ned]].constr[j]].vec_dof_1[i];
		}
	      } else {
		values[1] = field_p->dof_edge[tria_fa_edges[used_faces[ned][1]][ned]].vec_dof_1[i];
	      }
	      
	      field_p->dof_edge[nedge].vec_dof_1[i] = constr_coeff.val[constr_coeff.id[5]]*(values[0]+values[1]);
	    }

	    if(field_p->nr_sol>1){
	      field_p->dof_edge[nedge].vec_dof_2 = (double *) malloc(nreq*sizeof(double));

	      for(i=0;i<nreq;i++) { 

		//Edge 1
		if(field_p->dof_edge[tria_fa_edges[used_faces[ned][0]][ned]].constr != NULL) {
		  values[0] = 0.0;
		  for(j=1; j<field_p->dof_edge[tria_fa_edges[used_faces[ned][0]][ned]].constr[0]; j++) {
		    values[0] += 
		      constr_coeff.val[constr_coeff.id[field_p->dof_edge[tria_fa_edges[used_faces[ned][0]][ned]].constr_type[0]]+j]*
		      field_p->dof_edge[field_p->dof_edge[tria_fa_edges[used_faces[ned][0]][ned]].constr[j]].vec_dof_2[i];
		  }
		} else {
		  values[0] = field_p->dof_edge[tria_fa_edges[used_faces[ned][0]][ned]].vec_dof_2[i];
		}

		//Edge 2
		if(field_p->dof_edge[tria_fa_edges[used_faces[ned][1]][ned]].constr != NULL) {
		  values[1] = 0.0;
		  for(j=1; j<field_p->dof_edge[tria_fa_edges[used_faces[ned][1]][ned]].constr[0]; j++) {
		    values[1] +=
		      constr_coeff.val[constr_coeff.id[field_p->dof_edge[tria_fa_edges[used_faces[ned][1]][ned]].constr_type[0]]+j]*
		      field_p->dof_edge[field_p->dof_edge[tria_fa_edges[used_faces[ned][1]][ned]].constr[j]].vec_dof_2[i];
		  }
		} else {
		  values[1] = field_p->dof_edge[tria_fa_edges[used_faces[ned][1]][ned]].vec_dof_2[i];
		}
		
		field_p->dof_edge[nedge].vec_dof_2[i] = constr_coeff.val[constr_coeff.id[5]]*(values[0]+values[1]);	  
	      }
	    }
	    if(field_p->nr_sol>2){
	      field_p->dof_edge[nedge].vec_dof_3 = (double *) malloc(nreq*sizeof(double));
	      
	      for(i=0;i<nreq;i++) { 
		
		//Edge 1
		if(field_p->dof_edge[tria_fa_edges[used_faces[ned][0]][ned]].constr != NULL) {
		  values[0] = 0.0;
		  for(j=1; j<field_p->dof_edge[tria_fa_edges[used_faces[ned][0]][ned]].constr[0]; j++) {
		    values[0] += 
		      constr_coeff.val[constr_coeff.id[field_p->dof_edge[tria_fa_edges[used_faces[ned][0]][ned]].constr_type[0]]+j]*
		      field_p->dof_edge[field_p->dof_edge[tria_fa_edges[used_faces[ned][0]][ned]].constr[j]].vec_dof_3[i];
		  }
		} else {
		  values[0] = field_p->dof_edge[tria_fa_edges[used_faces[ned][0]][ned]].vec_dof_3[i];
		}

		//Edge 2
		if(field_p->dof_edge[tria_fa_edges[used_faces[ned][1]][ned]].constr != NULL) {
		  values[1] = 0.0;
		  for(j=1; j<field_p->dof_edge[tria_fa_edges[used_faces[ned][1]][ned]].constr[0]; j++) {
		    values[1] +=
		      constr_coeff.val[constr_coeff.id[field_p->dof_edge[tria_fa_edges[used_faces[ned][1]][ned]].constr_type[0]]+j]*
		      field_p->dof_edge[field_p->dof_edge[tria_fa_edges[used_faces[ned][1]][ned]].constr[j]].vec_dof_3[i];
		  }
		} else {
		  values[1] = field_p->dof_edge[tria_fa_edges[used_faces[ned][1]][ned]].vec_dof_3[i];
		}

		field_p->dof_edge[nedge].vec_dof_2[i] = constr_coeff.val[constr_coeff.id[5]]*(values[0]+values[1]);	  
	      }
	    }
	  }
	}
      }
    }
  }

  // Vertices nodes
  nr_dof_vert=0; nr_constr_dof_vert=0;
  for(nvert=1; nvert<=max_vert_id; nvert++) {
    if(mmr_node_status(mesh_id,nvert) == MMC_FREE) {
      SAFE_FREE(field_p->dof_vert[nvert].vec_dof_1);
      SAFE_FREE(field_p->dof_vert[nvert].vec_dof_2);
      SAFE_FREE(field_p->dof_vert[nvert].vec_dof_3);
      SAFE_FREE(field_p->dof_vert[nvert].constr);
      SAFE_FREE(field_p->dof_vert[nvert].constr_type);
      field_p->dof_vert[nvert].active = 0; //Set vertex status to inactive
    } else {
      if(mmr_node_status(mesh_id,nvert) == MMC_ACTIVE) {
	if(field_p->dof_vert[nvert].constr != NULL) {
	  nr_constr_dof_vert++;
	  continue;
	}
	nr_dof_vert++;
	field_p->dof_vert[nvert].active = 1; //Set vertex status to active
      }
    }
  }
  field_p->nr_dof_vert=nr_dof_vert;
  field_p->nr_constr_dof_vert=nr_constr_dof_vert;

  // Edge nodes
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
    nr_dof_edge=0; nr_constr_dof_edge=0;

    // Count active real and contrained edges & find used inactive edges
    for(ned=1; ned<=max_edge_id; ned++) {
      if(mmr_edge_status(mesh_id,ned) == MMC_ACTIVE) {
	if(field_p->dof_edge[ned].vec_dof_1 != NULL) {
	  nr_dof_edge++;
	} else if(field_p->dof_edge[ned].constr != NULL) {
	  nr_constr_dof_edge++;
	}
      }

      if(mmr_edge_status(mesh_id,ned) == MMC_INACTIVE) {
	if(field_p->dof_edge[ned].vec_dof_1 != NULL) {
	  if(field_p->dof_edge[ned].active == 1) {
	    nr_dof_edge++;
	  }
	  nr_constr_dof_edge++;
	  continue;
	}
	nr_dof_edge++;
      }
      /*
	if(field_p->dof_edge[ned].constr != NULL) {
	  for(i=1; i<=field_p->dof_edge[ned].constr[0]; i++) {
	    if(mmr_edge_status(mesh_id,field_p->dof_edge[ned].constr[i]) == MMC_INACTIVE) {
	      if(field_p->dof_edge[field_p->dof_edge[ned].constr[i]].active != 1) {
		field_p->dof_edge[field_p->dof_edge[ned].constr[i]].active = 1;
		nr_dof_edge++;
	      }
	    }
	  }
	  nr_constr_dof_edge++;
	  continue;
	}
	nr_dof_edge++;
      */
    }

    // Free data if it needed
    int free_constr_data; int free_real_data;

    for(ned=1; ned<=max_edge_id; ned++) {

      free_constr_data = 0; free_real_data = 0;
      
      if(mmr_edge_status(mesh_id,ned)==MMC_FREE) {
	free_constr_data = 1; free_real_data = 1;
      } else
      if(mmr_edge_status(mesh_id,ned) == MMC_INACTIVE) {
	free_constr_data = 1;
	if(field_p->dof_edge[ned].active == 0) {
	  free_real_data = 1;
	}
      }

      if(free_constr_data == 1) {
	  SAFE_FREE(field_p->dof_edge[ned].constr);
	  SAFE_FREE(field_p->dof_edge[ned].constr_type);
      }

      if(free_real_data == 1) {
	SAFE_FREE(field_p->dof_edge[ned].vec_dof_1);
	SAFE_FREE(field_p->dof_edge[ned].vec_dof_2);
	SAFE_FREE(field_p->dof_edge[ned].vec_dof_3);
	field_p->dof_edge[ned].active = 0;
      }
    }

    field_p->nr_dof_edge=nr_dof_edge;
    field_p->nr_constr_dof_edge=nr_constr_dof_edge;
  }

  if(field_p->pdeg == APC_LINEAR_APPROXIMATION_PDEG) { // Do it only for linear approximation
    nr_dof_ents = nr_dof_vert;
  }
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
    nr_dof_ents = nr_dof_vert + nr_dof_edge;
  }
  field_p->nr_dof_ents=nr_dof_ents;


  #ifdef DEBUG_APM
  /*
  int aV, iV, cV; int aE, iE, cE;
  aV = 0; iV = 0; cV = 0; aE = 0; iE = 0; cE = 0;
  
  printf("\e[1;31m--- (begin) apr_proj_dof_ref - debug print (begin) ---\e[0;0m\n");
  printf("\n\t\e[1;35mVertices nodes:\e[0;0m\n\n");
  
  for(nvert=1; nvert<=max_vert_id; nvert++) {
    if(mmr_node_status(mesh_id,nvert) == MMC_ACTIVE) {
      if(field_p->dof_vert[nvert].vec_dof_1 != NULL) {
	aV++;
        printf("[\e[1;32mREAL\e[0;0m]\t\t dof VERT id: %d,\tdata: ",nvert);
        for(i=0;i<nreq;i++) { printf("%lf\t",field_p->dof_vert[nvert].vec_dof_1[i]); }
        printf("\n");
      } else
      if(field_p->dof_vert[nvert].constr != NULL) {
	cV++;
        printf("[\e[1;33mCONSTRAINT\e[0;0m]\t dof VERT id: %d\n",nvert);
	for(i=0; i<=field_p->dof_vert[nvert].constr[0]; i++) {
          if(i == 0) {
	    printf("\tNumber of constraints: %d\n",field_p->dof_vert[nvert].constr[i]);
	    printf("\tCoefficient vector id: %d\n",field_p->dof_vert[nvert].constr_type[i]);
          } else {
	    printf("\tConstr = %d\tConstr_type = ",field_p->dof_vert[nvert].constr[i]);
	    if(field_p->dof_vert[nvert].constr_type[i] == APC_VERTEX) {
	      printf("APC_VERTEX\n");
	    } else {
	      if(field_p->dof_vert[nvert].constr_type[i] == APC_EDGE) {
	        printf("APC_EDGE\n");
	      }
	    }
          }
	}
      }
    } else
    if(mmr_node_status(mesh_id,nvert) == MMC_INACTIVE) {
      if(field_p->dof_vert[nvert].vec_dof_1 != NULL) {
	iV++;
	printf("[\e[1;30mINACTIVE ");
	if(field_p->dof_vert[nvert].active == 1)
	  printf("\e[1;32m(U)\e[0m");
        printf("\e[0;0m]\t dof VERT id: %d,\tdata: ",nvert);
	for(i=0;i<nreq;i++) { printf("%lf\t",field_p->dof_vert[nvert].vec_dof_1[i]); }
        printf("\n");
      } else
      if(field_p->dof_vert[nvert].constr != NULL) {
        printf("[\e[1;33mINACTIVE CONSTRAINT\e[0;0m] dof VERT id: %d\n",nvert);
	for(i=0; i<=field_p->dof_vert[nvert].constr[0]; i++) {
          if(i == 0) {
	    printf("\tNumber of constraints: %d\n",field_p->dof_vert[nvert].constr[i]);
	    printf("\tCoefficient vector id: %d\n",field_p->dof_vert[nvert].constr_type[i]);
          } else {
	    printf("\tConstr = %d\tConstr_type = ",field_p->dof_vert[nvert].constr[i]);
	    if(field_p->dof_vert[nvert].constr_type[i] == APC_VERTEX) {
	      printf("APC_VERTEX\n");
	    } else {
	      if(field_p->dof_vert[nvert].constr_type[i] == APC_EDGE) {
	        printf("APC_EDGE\n");
	      }
	    }
          }
	}
      }
    }
  }
  printf("VERTEX -> ACTIVE: %d ; INACTIVE: %d ; CONSTRAINED: %d\n",aV,iV,cV);

  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
    printf("\n\t\e[1;35mEdges nodes:\e[0;0m\n\n");

    for(nedge=1; nedge<=max_edge_id; nedge++) {
      if(mmr_edge_status(mesh_id,nedge) == MMC_ACTIVE) {
	if(field_p->dof_edge[nedge].vec_dof_1 != NULL) {
	  aE++;
	  printf("[\e[1;32mREAL\e[0;0m]\t\t dof EDGE id: %d,\tdata: ",nedge);
	  for(i=0;i<nreq;i++) { printf("%lf\t",field_p->dof_edge[nedge].vec_dof_1[i]); }
	  printf("\n");
	} else
	  if(field_p->dof_edge[nedge].constr != NULL) {
	    cE++;
	    printf("[\e[1;33mCONSTRAINT\e[0;0m]\t dof EDGE id: %d\n",nedge);
	    for(i=0; i<=field_p->dof_edge[nedge].constr[0]; i++) {
	      if(i == 0) {
		printf("\tNumber of constraints: %d\n",field_p->dof_edge[nedge].constr[i]);
		printf("\tCoefficient vector id: %d\n",field_p->dof_edge[nedge].constr_type[i]);
	      } else {
		printf("\tConstr = %d\tConstr_type = ",field_p->dof_edge[nedge].constr[i]);
		if(field_p->dof_edge[nedge].constr_type[i] == APC_VERTEX) {
		  printf("APC_VERTEX\n");
		} else {
		  if(field_p->dof_edge[nedge].constr_type[i] == APC_EDGE) {
		    printf("APC_EDGE\n");
		  }
		}
	      }
	    }
	  }
      } else
	if(mmr_edge_status(mesh_id,nedge) == MMC_INACTIVE) {
	  if(field_p->dof_edge[nedge].vec_dof_1 != NULL) {
	    iE++;
	    printf("[\e[1;30mINACTIVE ");
	    if(field_p->dof_edge[nedge].active == 1)
	      printf("\e[1;32m(U)\e[0m");
	    printf("\e[0;0m]\t dof EDGE id: %d,\tdata: ",nedge);
	    for(i=0;i<nreq;i++) { printf("%lf\t",field_p->dof_edge[nedge].vec_dof_1[i]); }
	    printf("\n");
	  } else
	    if(field_p->dof_edge[nedge].constr != NULL) {
	      iE++;
	      printf("[\e[1;30mINACTIVE CONSTRAINT\e[0;0m] dof EDGE id: %d\n",nedge);
	      for(i=0; i<=field_p->dof_edge[nedge].constr[0]; i++) {
		if(i == 0) {
		  printf("\tNumber of constraints: %d\n",field_p->dof_edge[nedge].constr[i]);
		  printf("\tCoefficient vector id: %d\n",field_p->dof_edge[nedge].constr_type[i]);
		} else {
		  printf("\tConstr = %d\tConstr_type = ",field_p->dof_edge[nedge].constr[i]);
		  if(field_p->dof_edge[nedge].constr_type[i] == APC_VERTEX) {
		    printf("APC_VERTEX\n");
		  } else {
		    if(field_p->dof_edge[nedge].constr_type[i] == APC_EDGE) {
		      printf("APC_EDGE\n");
		    }
		  }
		}
	      }
	    }
	}
    }
    printf("EDGES -> ACTIVE: %d ; INACTIVE: %d ; CONSTRAINED: %d\n",aE,iE,cE);
  }

  getchar();
  printf("Press ENTER to continue...\n"); getchar();
  printf("\e[1;31m--- (end) apr_proj_dof_ref - debug print (end) ---\e[0;0m\n");
  /**/

  apr_check_field(Field_id);
  #endif

  SAFE_FREE(edge_elems);
  
  return(1);
}


/*------------------------------------------------------------
  apr_rewr_sol - to rewrite solution from one vector to another
------------------------------------------------------------*/
int apr_rewr_sol(
  int Field_id,      /* in: data structure to be used  */
  int Sol_from,      /* in: ID of vector to read solution from */
  int Sol_to         /* in: ID of vector to write solution to */
  )
{

  int nno, mesh_id, num_shap, nreq; // Changed on 12.2010 - nreq added
  double dofs_loc[APC_MAXELSD];

/*++++++++++++++++ executable statements ++++++++++++++++*/

  mesh_id = apr_get_mesh_id(Field_id);
  nreq = apr_get_nreq(Field_id); // Changed on 12.2010 - line added

  nno=0;
  while((nno=mmr_get_next_node_all(mesh_id,nno))!=0){
    
    if(apr_get_ent_pdeg(Field_id, APC_VERTEX, nno)>0){
      // faster: select field_p and next
      // if(field_p->dof_ents[Ent_id].vec_dof_1==NULL){

      num_shap = apr_get_ent_numshap(Field_id, APC_VERTEX, nno);
      apr_read_ent_dofs(Field_id, APC_VERTEX, nno, // Changed on 12.2010
			num_shap*nreq, Sol_from, dofs_loc); // - num_shap*nreq
      apr_write_ent_dofs(Field_id, APC_VERTEX, nno, // Changed on 12.2010
			 num_shap*nreq, Sol_to, dofs_loc); // - num_shap*nreq

    }
  }

  return 0;
}


/*------------------------------------------------------------
  apr_free_field - to free approximation field data structure
------------------------------------------------------------*/
int apr_free_field(
  int Field_id    /* in: approximation field ID  */
  )
{
  /* pointer to field structure */
  apt_field *field_p;
  int i, nmvert, nmedge, mesh_id;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);

  /* select the corresponding mesh */
  mesh_id = apr_get_mesh_id(Field_id);

  /* loop over dof vertexes data structure including free spaces */
  nmvert=mmr_get_max_node_id(mesh_id);
  for(i=1;i<=nmvert;i++){
    /* clean the space if it has dofs */
    SAFE_FREE(field_p->dof_vert[i].vec_dof_1);
    SAFE_FREE(field_p->dof_vert[i].vec_dof_2);
    SAFE_FREE(field_p->dof_vert[i].vec_dof_3);

    SAFE_FREE(field_p->dof_vert[i].constr);
    SAFE_FREE(field_p->dof_vert[i].constr_type);
  }
  SAFE_FREE(field_p->dof_vert);
  
  /* loop over dof edges data structure including free spaces */
  if(field_p->pdeg == APC_QUADRATIC_APPROXIMATION_PDEG) { // Do it only for quadratic approximation
    nmedge=mmr_get_max_edge_id(mesh_id);
    for(i=1;i<=nmedge;i++){
      /* clean the space if it has dofs */
	SAFE_FREE(field_p->dof_edge[i].vec_dof_1);
	SAFE_FREE(field_p->dof_edge[i].vec_dof_2);
	SAFE_FREE(field_p->dof_edge[i].vec_dof_3);

	SAFE_FREE(field_p->dof_edge[i].constr);
	SAFE_FREE(field_p->dof_edge[i].constr_type);
    }
    SAFE_FREE(field_p->dof_edge);
  }
  return(1);
}


/*------------------------------------------------------------
  apr_limit_deref - to return whether derefinement is allowed or not
    the routine does not use approximation data structures but the result
    depends on the approximation method
------------------------------------------------------------*/
int apr_limit_deref(
  int Field_id,    /* in: approximation field ID  */
  int El_id      /* in: unique identifier of the element */
  )
{
  int ied, ned, el_edges[13], el_edges_2[13];
  int i, j; 
  int edge_elems[20+1];
  int elsons[MMC_MAXELSONS+1]; /* family information */
  int mesh_id, father, son, ison, max_gen_diff, ignore_edges, ref_type;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  mesh_id = apr_get_mesh_id(Field_id);
  max_gen_diff=mmr_get_max_gen_diff(mesh_id);

  if(max_gen_diff== 1){
	      
    father=mmr_el_fam(mesh_id,El_id,NULL,NULL);
    mmr_el_fam(mesh_id,father,elsons,NULL);

    /* we can cluster elements only if all their edges do not contain */
    /* constrained nodes, i.e. are active */
    /* and also edges of all their neighbours across edges... */
    ref_type = mmr_el_type_ref(mesh_id,father);

    if(ref_type==MMC_REF_ANI) ignore_edges=3;
    else ignore_edges=0;
    
    for(ison=1;ison<=elsons[0];ison++){
      
      son = elsons[ison];

      mmr_el_edges(mesh_id, son, el_edges);
      
      /* for all edges */
      for(ied=1;ied<=el_edges[0]-ignore_edges;ied++){
	
	ned=el_edges[ied];
	
	/* find elements containing the edge */
	mmr_edge_elems(mesh_id, ned, edge_elems);
	
	for(i=1;i<=edge_elems[0];i++){
	  
	  mmr_el_edges(mesh_id, edge_elems[i], el_edges_2);
	  
	  for(j=1;j<el_edges_2[0]-ignore_edges;j++){
	    
	    if(mmr_edge_status(mesh_id, el_edges_2[j])
	       ==MMC_INACTIVE) {
	      
	      return(APC_DEREF_DENIED);
	      
	    }
	  } /* end loop over edges */
	} /* end loop over elements sharing edges with sons */
      } /* end loop over sons' edgse */
    } /* end loop over sons */
  } /* end if(max_gen_diff==1) */
  return(APC_DEREF_ALLOWED);
}

/*------------------------------------------------------------
  apr_limit_ref - to return whether refinement is allowed or not
    the routine does not use approximation data structures but the result
    depends on the approximation method
------------------------------------------------------------*/
int apr_limit_ref(
  int Field_id,    /* in: approximation field ID  */
  int El_id      /* in: unique identifier of the element */
  )
{

  int mesh_id, gen_el, max_gen;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  mesh_id = apr_get_mesh_id(Field_id);
  max_gen=mmr_get_max_gen(mesh_id);

  if(max_gen==0) return(APC_REF_DENIED);
  else if(max_gen>0){
    
    /* check whether max generation is not exceeded */
    gen_el=mmr_el_gen(mesh_id,El_id);

    if(gen_el<max_gen) return(APC_REF_ALLOWED);
    else return(APC_REF_DENIED);

  }

  return(APC_REF_ALLOWED);
}

/*------------------------------------------------------------
  apr_refine - to refine an element or the whole mesh checking mesh irregularity
------------------------------------------------------------*/
int apr_refine( /* returns: >=0 - success code, <0 - error code */
  int Field_id,  /* in: field ID */
  int El   /* in: element ID or -1 (MMC_DO_UNI_REF) for uniform refinement */
)
{

/* local variables */
  int i,iaux,iel,nmel_old,nrel_old,nrel_div,nelref,nrwait,listwait[100];
  int el_type, gen_el, max_gen, max_gen_diff; 
  int ifa, num_face, face, face_sons[5], ifaneig;
  int iprint=MMC_PRINT_INFO+1;
  int mesh_id;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh_id = apr_get_mesh_id(Field_id);

  max_gen=mmr_get_max_gen(mesh_id);
  max_gen_diff=mmr_get_max_gen_diff(mesh_id);
  
  if(El!=MMC_DO_UNI_REF&&El>0) {

/*set number of elements waiting for refinement (they will be stored in listwait)*/
    nrwait=0;

/*element to refine*/
    nelref=El;

    beginning:{}

/* check whether element is active */
    if(mmr_el_status(mesh_id, nelref)!=MMC_ACTIVE){
#ifdef DEBUG_APM
      if(iprint>MMC_PRINT_INFO){ 
	printf("not active element or free space %d in refine!\n", nelref);
      }
#endif
      return(-2);

      /*kb!!!
      if(nrwait==0) return(-2);
      nrwait--;
      nelref=listwait[nrwait];
      goto beginning;
      !!!kb*/
    }

/* barrier */
    gen_el=mmr_el_gen(mesh_id,nelref);

    if(gen_el>=max_gen){

#ifdef DEBUG_APM
      if(iprint>MMC_PRINT_INFO) {
	printf("element %d generation (%d) > max_gen (%d) in divide8_p!\n",nelref,gen_el,max_gen);
      }
#endif
      return(-3);
      /*kb!!!
      if(nrwait==0) return(-3);
      nrwait--;
      nelref=listwait[nrwait];
      goto beginning;
      !!!kb*/
    }


/* barrier */
    if(max_gen_diff == 0) {
      El = MMC_DO_UNI_REF;
    }
    else if(max_gen_diff == 1){

      /* 1-irregular meshes - no double constrained nodes */

      /* proper version: */
/* 1. find the parent */
/* 2. check the position within parent */
/* 3. identify parent's edges concerned  */
/* 4. check the respective parent's edges do not contain constrained nodes */

      /* simplified version with over-refinements: */
      /* all parent's edges must not contain constrained nodes */

      int parent;
      int ned, ied, nr_edges, el_edges[13];
      int edge_elems[20+1];

      parent = mmr_el_fam(mesh_id, nelref, NULL, NULL);

      if(parent!=MMC_NO_FATH){
	
	mmr_el_edges(mesh_id, parent, el_edges);
	
	/* for all parent's edges */
	nr_edges = el_edges[0];
	/* for 2D refinements we check only horizontal edges */
	if(mmr_el_type_ref(mesh_id,parent)==MMC_REF_ANI) {
	  nr_edges = 6;
	}

	for(ied=1;ied<=nr_edges;ied++){
	  
	  ned = el_edges[ied];
	  
	  mmr_edge_elems(mesh_id, ned, edge_elems);
	  
	  for(i=1;i<=edge_elems[0];i++){
	    
/*kbw
if(nelref==724||parent==724){
    printf("Edge %d edge_elems: ", ned);
    for(i=1; i<=edge_elems[0]; i++){
      printf("%d  ", edge_elems[i]);
    }
    printf("\n");
  }
/*kew*/
	    if( edge_elems[i]!=nelref &&
		mmr_el_status(mesh_id,edge_elems[i])==MMC_ACTIVE) {
	      
	      listwait[nrwait]=nelref;
	      nrwait++;
	      nelref= edge_elems[i];
	      
	      
/*kbw
#ifdef DEBUG_APM
  if(iprint>MMC_PRINT_INFO){
    printf("adding element %d to refine because of DIFF_GEN\n",nelref);
    printf("nelref %d, parent %d, edge %d, elem_new %d\n",
	   listwait[nrwait-1], parent, ned, edge_elems[i]);

  }
#endif
/*kew*/
	      
	      goto beginning;
	      
	    }
	  }
	}
      } /* end if not initial mesh element */

    } /* end if max_gen_diff == -1 i.e. 1-irregular meshes */
    else{
      printf("\e[1;31mERROR\e[0;0m] Wrong parameter max_gen_diff %d in apr_refine! Exiting.\n", max_gen_diff);
      exit(-1);
    }
    

/* divide an element */
    iaux=mmr_refine_el(mesh_id,nelref);

    if(iaux<0) return(iaux-10);

/*get the next element from the list (if any)*/      
    if(nrwait==0) return(iaux);
    nrwait--;
    nelref=listwait[nrwait];
    goto beginning;

  }
  else if(El==MMC_DO_UNI_REF){
    return(mmr_refine_mesh(mesh_id));
  }
  else{
    printf("\e[1;31mERROR\e[0;0m] Wrong parameter %d in apr_refine! Exiting.\n", El);
    exit(-1);
  }


#ifdef DEBUG_APM
  printf("Error 230956 in apr_refine!\n");
#endif

/* error condition - that point should not be reached */
  return(-1);
}

/*------------------------------------------------------------
  apr_derefine - to derefine an element or the whole mesh with irregularity check
------------------------------------------------------------*/
int apr_derefine( /* returns: >=0 - success code, <0 - error code */
  int Field_id,	/* in: mesh ID or 0 (MMC_CUR_MESH_ID) for the current mesh */
  int El  /* in: element ID or -2 (MMC_DO_UNI_DEREF) for uniform derefinement */
  )
{

/* local variables */
  int i, iaux,iel,nmel_old,nrel_old,nrel_div;
  int iprint=5;
  int mesh_id;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  mesh_id = apr_get_mesh_id(Field_id);

  if(El!=MMC_DO_UNI_DEREF && El>0) {
    
    if(MMC_ACTIVE == mmr_el_status(mesh_id,El)
       && APC_DEREF_ALLOWED == apr_limit_deref(Field_id,El) )
      {
	iaux=mmr_derefine_el(mesh_id,El);
      }
    
  }
  else if(El==MMC_DO_UNI_DEREF){
    iaux=mmr_derefine_mesh(mesh_id);
  }
  else{
    printf("\e[1;31mERROR\e[0;0m] Wrong parameter %d in apr_derefine! Exiting.\n", El);
    exit(-1);
  }

  return(iaux);
}

