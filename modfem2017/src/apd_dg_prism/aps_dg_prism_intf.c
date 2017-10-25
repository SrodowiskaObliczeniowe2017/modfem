/************************************************************************
File aps_dg_prism_intf.c - implementation of the approximation module 
                  interface for the discontinuous Galerkin discretization 
                  of PDEs on 3D meshes with prismatic elements 
 
Contains definitions of routines:   
  apr_module_introduce - to return the approximation method's name
  apr_init_field - to initiate new approximation field and read its data
  apr_write_field - to dump-out field data in the standard MOD_FEM format
  apr_check_field - to check approximation field data structure
  apr_get_mesh_id - to return the ID of the associated mesh
  apr_get_nreq - to return the number of components in solution vector
  apr_get_nr_sol - to return the number of solution vectors stored
  apr_get_el_pdeg_numshap - to return the number of shape functions
                            (scalar DOFs) for an element given its 
                            degree of approximation symbol or vector pdeg
  apr_get_base_type - to return the type of basis functions
  apr_get_ent_pdeg - to return the degree of approximation symbol
                      associated with a given mesh entity
  apr_set_ent_pdeg - to set the degree of approximation index 
                      associated with a given mesh entity
  apr_get_el_pdeg - to return the degree of approximation vector 
                      associated with a given element
  apr_set_el_pdeg - to set the degree of approximation vector 
                      associated with a given element
  apr_get_ent_numshap - to return the number of shape functions (vector DOFs)
                        associated with a given mesh entity
  apr_get_ent_nrdofs - to return the number of dofs associated with 
                      a given mesh entity
  apr_get_nrdofs_glob - to return a global dimension of the problem
  apr_read_ent_dofs - to read a vector of dofs associated with a given
                  mesh entity from approximation field data structure
  apr_write_ent_dofs - to write a vector of dofs associated with a given
                   mesh entity to approximation field data structure
  apr_create_ent_dofs - to create a vector of dofs associated with a given
                   mesh entity within approximation field data structure
  apr_set_ini_con - to set an initial condition
  apr_num_int_el - to perform numerical integration for an element
  apr_proj_dof_ref - to rewrite dofs after modifying the mesh
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
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

/* interface of the mesh manipulation module */
#include "mmh_intf.h"	

/* interface for all approximation modules */
#include "aph_intf.h"	

#include "uth_intf.h"

/* internal header file for the dg approximation module */
#include "./aph_dg_prism.h"	

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
  apr_module_introduce - to return the approximation method's name
------------------------------------------------------------*/
int apr_module_introduce( 
                  /* returns: >=0 - success code, <0 - error code */
  char* Approx_name /* out: the name of the approximation method */
  )
{
  char* string = "DG_SCALAR_PRISM";

  strcpy(Approx_name,string);

  return(1);
}

/*------------------------------------------------------------
  apr_init_field - to initiate new approximation field and 
                   read its control parameters
------------------------------------------------------------*/
int apr_init_field(  /* returns: >0 - field ID, <0 - error code */
  char Field_type, /* options: s - standard */ 
		   /*          c - discontinuous with complete basis */
                   /*          t - discontinuous with tensor product basis */
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
  int pdeg, mxel, nel, iel, num_dof;
  
  /* auxiliary variables */
  int i,j,k,iaux,jaux;
  FILE *fp;
  char mesh_name[100];

/*++++++++++++++++ executable statements ++++++++++++++++*/


#ifdef DEBUG_APM
  if(Field_type != 'c' && Field_type != 't'){
    printf("Wrong Field_type %c in apr_init_field for DG approximation!\n"
           , Field_type);
    exit(-1);
  }
#endif

/*kbw
  printf("Initializing field: type %c, control %c\n", Field_type, Control);
  printf("Mesh_id %d, Nreq %d, Nr_sol %d, Pdeg %d\n",
	 Mesh_id, Nreq, Nr_sol, Pdeg_in);
/*kew*/


  /* first check that the mesh is compatible */
  mmr_module_introduce(mesh_name);

#ifdef DEBUG_APM
  if(strcmp("3D_PRISM",mesh_name)){
    // if(monitor>APC_NO_PRINT){
    printf("Wrong mesh %s for discontinuous Galerkin approximation with prisms!\n",
	   mesh_name);
    // {
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

  /* set type of basis functions */
  if(Field_type == 'c'){
    field_p->base=APC_BASE_COMPLETE_DG;
  }
  else {
    field_p->base=APC_BASE_TENSOR_DG;
  }

  /* depending on the control variable */
  if(Control==APC_READ){

    /* open the input file with approximation parameters */
    fp=fopen(Filename,"r");
#ifdef DEBUG_APM
    if(fp==NULL) {
      printf("Not found file '%s' with field data!!! Exiting.\n",Filename);
      exit(-1);
    } 
#endif

    /* read the number of solution vectors for each dof entity and */
    /* the type of shape functions */
    fscanf(fp,"%d %d",
	   &field_p->nreq, &field_p->nr_sol);

#ifdef DEBUG_APM
    if(field_p->nreq!=Nreq) {
      printf("Requested number of equations %d different than read from file %d!!!\n",
	     Nreq, field_p->nreq);
      printf("Requested number of components %d different than read from file %d!!!\n",
	     Nr_sol, field_p->nr_sol);
      printf("Exiting.\n");
      exit(-1);
    }
#endif

  }
  else{

    field_p->nreq = Nreq;
    field_p->nr_sol = Nr_sol;

  }

  if(field_p->nreq>APC_MAXEQ) {
    printf("Requested number of equations %d greater than APC_MAXEQ=%d!!!\n",
	   field_p->nreq, APC_MAXEQ);
    printf("Change APC_MAXEQ in aph_dg_prism.h and recompile the code.\n");
    printf("Exiting.\n");
    exit(-1);
  }
 
 /* get the maximal admisible number of elements */
  mxel=mmr_get_max_elem_max(Mesh_id);

  /* allocate the space for degrees of freedom */
  field_p->dof_ents = (apt_dof_ent *) malloc((mxel+1)*sizeof(apt_dof_ent));
  if(field_p->dof_ents==NULL) {
    printf("Dofs structures not allocated\n");
    exit(-1);
  }

  /* depending on the control variable */
  if(Control==APC_READ){
    /* read degree of polynomial for each element (or uniform p if pdeg <0 ) */
    fscanf(fp,"%d",&pdeg);

  }
  else{

    pdeg = -Pdeg_in;

  }

  if(pdeg<0) {
    pdeg = -pdeg;

    field_p->uniform = 1;

    if((field_p->base==APC_BASE_COMPLETE_DG && pdeg>APC_MAXELP_COMP) || 
       (field_p->base==APC_BASE_TENSOR_DG && pdeg>APC_MAXELP_TENS) ){
      printf("Too big degree of polynomial: %d\n",pdeg);
      exit(-1);
    }

    iaux=0;
    for(nel=1;nel<=mxel;nel++){
      if(mmr_el_status(Mesh_id,nel)==MMC_ACTIVE){
	iaux++;
	field_p->dof_ents[nel].pdeg=pdeg;
      }
      else {
	field_p->dof_ents[nel].pdeg = APC_NO_DOFS;
      }
    }
  }
  else{
    field_p->uniform = 0;
    iaux=0;
    for(nel=1;nel<=mxel;nel++){
      if(mmr_el_status(Mesh_id,nel)==MMC_ACTIVE){
	iaux++;
	if(iaux==1)  field_p->dof_ents[nel].pdeg=pdeg;
	else fscanf(fp,"%d\n", &field_p->dof_ents[nel].pdeg);
      }
      else {
	field_p->dof_ents[nel].pdeg = APC_NO_DOFS;
      }

    }
  }

  /* create dofs data structure */
  for(nel=1;nel<=mxel;nel++){
    if(mmr_el_status(Mesh_id,nel)==MMC_ACTIVE){
	
      jaux=field_p->dof_ents[nel].pdeg;
      num_dof=apr_get_ent_nrdofs(APC_CUR_FIELD_ID,APC_ELEMENT,nel);

#ifdef DEBUG_APM
      if(num_dof!=
         field_p->nreq*apr_get_el_pdeg_numshap(APC_CUR_FIELD_ID,nel,&jaux)){
	printf("Something wrong in init_field!!!\n");
	exit(-1);
      }
#endif

      field_p->dof_ents[nel].vec_dof_1 = 
	(double *) malloc(num_dof*sizeof(double));
      if(field_p->nr_sol>1){
	field_p->dof_ents[nel].vec_dof_2 = 
	  (double *) malloc(num_dof*sizeof(double));
      }	  
      if(field_p->nr_sol>2){
	field_p->dof_ents[nel].vec_dof_3 = 
	  (double *) malloc(num_dof*sizeof(double));
      }	  

      /* depending on the control variable */
      if(Control==APC_READ){
	
	/* read values from a dump file */
	for(i=0;i<num_dof;i++){
	  fscanf(fp,"%lg",&field_p->dof_ents[nel].vec_dof_1[i]);
	} 
	if(field_p->nr_sol>1){
	  for(i=0;i<num_dof;i++){
	    fscanf(fp,"%lg",&field_p->dof_ents[nel].vec_dof_2[i]);
	  }
	}
	if(field_p->nr_sol>2){
	  for(i=0;i<num_dof;i++){
	    fscanf(fp,"%lg",&field_p->dof_ents[nel].vec_dof_3[i]);
	  }
	}
      }
      else{
	/* initiate the field to zero */
	for(i=0;i<num_dof;i++){
	  field_p->dof_ents[nel].vec_dof_1[i]=0.0;
	}
	if(field_p->nr_sol>1){
	  for(i=0;i<num_dof;i++){
	    field_p->dof_ents[nel].vec_dof_2[i]=0.0;
	  }
	}
	if(field_p->nr_sol>2){
	  for(i=0;i<num_dof;i++){
	    field_p->dof_ents[nel].vec_dof_3[i]=0.0;
	  }
	}
      }
    }
    else {
      field_p->dof_ents[nel].vec_dof_1 = NULL;
    }
  }

/* close file with approximation data */
  if(Control==APC_READ) fclose(fp);

  /* initiate values using provided function */
  if(Control==APC_INIT) {
    apr_set_ini_con(apv_cur_field_id,Fun_p);
#ifdef DEBUG_APM
    printf("\nSpecified initial data for approximation field %d!\n",
	   apv_cur_field_id);
#endif
  }

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
  int mesh_id, nmel, nel, numshap, i;
  FILE *fp;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);
  
  /* select the corresponding mesh */
  mesh_id = apr_get_mesh_id(Field_id);

  /* open the input file with approximation parameters */
  fp=fopen(Filename,"w");

  /* read the number of solution vectors for each dof entity and */
  /* the type of shape functions */
  fprintf(fp,"%d %d\n",
	  field_p->nreq, field_p->nr_sol);

  /* loop over dofs data structure including free spaces */
  nmel=mmr_get_max_elem_id(mesh_id);


  if(field_p->uniform==1){
    
    for(nel=1;nel<=nmel;nel++){
      
      /* for the first active element */
      if(mmr_el_status(mesh_id,nel)==MMC_ACTIVE){
	fprintf(fp,"%d\n", -field_p->dof_ents[nel].pdeg);
	break;
      }
    }
  }
  else{
    for(nel=1;nel<=nmel;nel++){
      
      /* for all active elements */
      if(mmr_el_status(mesh_id,nel)==MMC_ACTIVE){
	fprintf(fp,"%d\n", -field_p->dof_ents[nel].pdeg);
      }
    }
  }


  for(nel=1;nel<=nmel;nel++){

/* for active elements */
    if(mmr_el_status(mesh_id,nel)==MMC_ACTIVE){

      numshap=apr_get_ent_numshap(Field_id,APC_ELEMENT,nel);

      for(i=0;i<Nreq*numshap;i++){

	if(fabs(field_p->dof_ents[nel].vec_dof_1[i])<1e-15){
	  field_p->dof_ents[nel].vec_dof_1[i]=0.0;
	}
	fprintf(fp,"%.12lg ",field_p->dof_ents[nel].vec_dof_1[i]);
        /*
          utr_fprintf_double(fp, Accuracy, field_p->dof_ents[nel].vec_dof_1[i]);
          fprintf(fp," ");        
        */
      }
      fprintf(fp,"\n");

      if(field_p->nr_sol>1){

	for(i=0;i<Nreq*numshap;i++){

          
	  if(fabs(field_p->dof_ents[nel].vec_dof_2[i])<1e-15){
	    field_p->dof_ents[nel].vec_dof_2[i]=0.0;
	  }
	  fprintf(fp,"%.12lg ",field_p->dof_ents[nel].vec_dof_2[i]);
          /*
          utr_fprintf_double(fp, Accuracy, field_p->dof_ents[nel].vec_dof_2[i]);
          fprintf(fp," ");    
	  */             
	}
	fprintf(fp,"\n");

      }
    } /* end if element active */
  } /* end loop over all dofs structures */


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
  int i, nel, nmel, mesh_id;

/*++++++++++++++++ executable statements ++++++++++++++++*/
  
  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);
  
  /* select the corresponding mesh */
  mesh_id = apr_get_mesh_id(Field_id);

  /* loop over dofs data structure including free spaces */
  nmel=mmr_get_max_elem_id(mesh_id);

/*kbw
  printf("Checking field %d associated with mesh %d, nmel=%d\n",
	 Field_id, mesh_id, nmel);
/*kew*/

  for(nel=1;nel<=nmel;nel++){

/*kbw
    printf("Structure %d, el_status %d\n",
	   nel, mmr_el_status(mesh_id, nel));
/*kew*/

    if(mmr_el_status(mesh_id, nel)==MMC_FREE){

      if(field_p->dof_ents[nel].pdeg != APC_NO_DOFS){
	printf("1 dofs associated with free space %d, exiting\n", nel);
	exit(-1);
      }
      if(field_p->dof_ents[nel].vec_dof_1 != NULL){
	printf("2 dofs associated with free space %d, exiting\n", nel);
	exit(-1);
      }

    }
    else if(mmr_el_status(mesh_id, nel)==MMC_INACTIVE){

/*kbw
      int el_sons[MMC_MAXELSONS+1], ison, el_fath, ref_type;
      el_fath=mmr_el_fam(mesh_id, nel, el_sons, &ref_type);
      printf("Inactive el %d, type %d, father %d, ref_type %d, sons :",
	     nel, mmr_el_type(mesh_id,nel), el_fath, ref_type);
      for(ison=1;ison<=el_sons[0];ison++){
	printf("%d ",el_sons[ison]);
      }
      printf("\n");
/*kew*/

      if(field_p->dof_ents[nel].pdeg != APC_NO_DOFS){
	printf("dofs associated with inactive element %d, exiting\n", nel);
	exit(-1);
      }
      if(field_p->dof_ents[nel].vec_dof_1 != NULL){
	printf("dofs associated with inactive element %d, exiting\n", nel);
	exit(-1);
      }

    }
    else{

      int el_sons[MMC_MAXELSONS+1], ison, el_fath, ref_type;
      el_fath=mmr_el_fam(mesh_id, nel, el_sons, &ref_type);
/*kbw
      printf("Active el %d, type %d, father %d\n",
	     nel, mmr_el_type(mesh_id,nel), el_fath);
/*kew*/
      if(el_sons[0]!=0 || ref_type!=MMC_NOT_REF){
	printf("Inconsistency in data structure for active el %d, exiting\n",
	       nel);
	exit(-1);
      }

      if(field_p->dof_ents[nel].pdeg == APC_NO_DOFS){
	printf("no dofs associated with active element %d, exiting\n", nel);
	exit(-1);
      }
      if(field_p->dof_ents[nel].vec_dof_1 == NULL){
	printf("no dofs associated with active element %d, exiting\n", nel);
	exit(-1);
      }

    }

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
  apr_get_el_pdeg_numshap - to return the number of shape functions
                            for an element given its 
                            degree of approximation symbol or vector pdeg
---------------------------------------------------------*/
int apr_get_el_pdeg_numshap(
		 /* returns: >=0 - success code, <0 - error code*/
  int Field_id,  /* in: field ID */
  int El_id,     /* in: element ID (not used - all elements have the same base */
  int *Pdeg_vec       /* in: degree of approximation symbol or vector */
  )
{
  /* pointer to field structure */
  apt_field *field_p;

  /* auxiliary variables */
  int base, porder, numshap, porderz, pdeg;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper field */  
  field_p = apr_select_field(Field_id);

  base = field_p->base;

  /* in DG approximation element pdeg concerns only interior - single number */
  pdeg = Pdeg_vec[0];

/* check */
#ifdef DEBUG_APM
  if(pdeg<0||(base==APC_BASE_COMPLETE_DG&&pdeg>100)){
    printf("Wrong pdeg %d in apr_get_el_pdeg_numshap\n",pdeg);
    exit(-1);
  }
#endif

  if(base==APC_BASE_TENSOR_DG){

    /* decipher horizontal and vertical orders of approximation */
    porderz=pdeg/100;
    porder=pdeg%100;
    numshap = (porderz+1)*(porder+1)*(porder+2)/2;

  }
  else if(base==APC_BASE_COMPLETE_DG){

    porder=pdeg;
    numshap = (porder+1)*(porder+2)*(porder+3)/6;

  }
  else {

    printf("Type of base in get_el_pdeg_numshap not valid for prisms!\n");
    exit(-1);

  }

  return(numshap);
}


/*---------------------------------------------------------
  apr_get_el_dofs - to return the number and the list of element's degrees
                        of freedom (corresponding to standard shape functions)
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

  int nreq, numshap, numdofs;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  nreq = apr_get_nreq(Field_id);

  numshap = apr_get_ent_numshap(Field_id, APC_ELEMENT, El_id);

  numdofs = nreq*numshap;

  apr_read_ent_dofs(Field_id,APC_ELEMENT,El_id,numdofs,Vect_id,El_dofs_std);

  return(numdofs);

}



/*---------------------------------------------------------
  apr_get_base_type - to return the type of basis functions
---------------------------------------------------------*/
int apr_get_base_type(/* returns: >0 - type of basis functions, 
                                  <0 - error code */
  int Field_id, /* in: field ID */
  int El_id     /* in: element ID (not used - all elements have the same base */

  )
{
  /* pointer to field structure */
  apt_field *field_p;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper field */  
  field_p = apr_select_field(Field_id);

  return(field_p->base);
}

/*------------------------------------------------------------
  apr_get_ent_pdeg - to return the degree of approximation index 
                      associated with a given mesh entity
------------------------------------------------------------*/
int apr_get_ent_pdeg( /* returns: >0 - approximation index,
                                          <0 - error code */
  int Field_id,     /* in: approximation field ID  */
  int Ent_type,      /* in: type of mesh entity */
  int Ent_id         /* in: mesh entity ID */
  )
{
  /* pointer to field structure */
  apt_field *field_p;
  
/*++++++++++++++++ executable statements ++++++++++++++++*/
  
#ifdef DEBUG_APM
  /* check input for DG approximation */
  if(abs(Ent_type)!=APC_ELEMENT){
    printf("Wrong dof entity type in get_ent_pdeg for DG approximation !\n");
    exit(-1);
  }
#endif
  
  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);
  
#ifdef DEBUG_APM
  /* check whether element is active */
  if(Ent_type>0){
    int mesh_id;
    mesh_id = field_p->mesh_id;
    if(mmr_el_status(mesh_id,Ent_id)!=MMC_ACTIVE){
      printf("Wrong element %d status in get_ent_pdeg for DG approximation !\n",
	     Ent_id);
      getchar();
      exit(-1);
    }
  }
#endif

  return(field_p->dof_ents[Ent_id].pdeg);  

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
  apt_field *field_p;
  
/*++++++++++++++++ executable statements ++++++++++++++++*/
  
#ifdef DEBUG_APM
  /* check input for DG approximation */
  if(abs(Ent_type)!=APC_ELEMENT){
    printf("Wrong dof entity type in get_ent_pdeg for DG approximation !\n");
    exit(-1);
  }
#endif
  
  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);
  
  field_p->dof_ents[Ent_id].pdeg=Pdeg;

  return(1);  

}


/*------------------------------------------------------------
  apr_get_el_pdeg - to return the degree of approximation vector 
                      associated with a given element
------------------------------------------------------------*/
int apr_get_el_pdeg( /* returns: >0 - approximation index,
				          0 - dof entity inactive (constrained)
                                         <0 - error code */
  int Field_id,     /* in: approximation field ID  */
  int Ent_id,         /* in: element ID */
  int *Pdeg_vec       /* out: degree of approximation symbol or vector */
  )
{

  /* pointer to field structure */
  apt_field *field_p;
  
/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);
  
  /* in DG approximation element pdeg concerns only interior - single number */
  if(Pdeg_vec != NULL) *Pdeg_vec = field_p->dof_ents[Ent_id].pdeg;

  return(field_p->dof_ents[Ent_id].pdeg);
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

  /* pointer to field structure */
  apt_field *field_p;
  
/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);
  
  /* in DG approximation element pdeg concerns only interior - single number */
  field_p->dof_ents[Ent_id].pdeg=Pdeg_vec[0];

  return(1);  
}


/*------------------------------------------------------------
  apr_get_ent_numshap - to return the number of dofs associated with 
                      a given mesh entity
------------------------------------------------------------*/
int apr_get_ent_numshap( /* returns: >0 - the number of dofs, 
                                          <0 - error code */
  int Field_id,     /* in: approximation field ID  */
  int Ent_type,      /* in: type of mesh entity */
  int Ent_id         /* in: mesh entity ID */
  )
{
  /* auxiliary variables */
  int pdeg, mesh_id;  

/*++++++++++++++++ executable statements ++++++++++++++++*/

#ifdef DEBUG_APM
  /* check input for DG approximation */
  if(abs(Ent_type)!=APC_ELEMENT){
    printf("Wrong dof entity type in get_ent_pdeg for DG approximation !\n");
    exit(-1);
  }
#endif

#ifdef DEBUG_APM
  mesh_id = apr_get_mesh_id(Field_id);
  /* check input for DG approximation */
  if(mmr_el_status(mesh_id,Ent_id)!=MMC_ACTIVE){
    printf("Asking for pdeg of inactive element in get_ent_pdeg !\n");
    exit(-1);
  }
#endif

  pdeg = apr_get_ent_pdeg(Field_id, APC_ELEMENT, Ent_id);

  return(apr_get_el_pdeg_numshap(Field_id, Ent_id, &pdeg));
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

 return(apr_get_nreq(Field_id)*apr_get_ent_numshap(Field_id, Ent_type, Ent_id));

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
  int nel, glob_dim=0;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);

  /* get the corresponding mesh ID */
  mesh_id = field_p->mesh_id;

  /* get number of equations */
  
  /* loop over elements */
  nel=0;
  while((nel=mmr_get_next_act_elem(mesh_id, nel))!=0){
    glob_dim += apr_get_ent_nrdofs(Field_id, APC_ELEMENT, nel);
  }

  return(glob_dim);
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
  int i, mesh_id;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);
  
  /* get the corresponding mesh ID */
  mesh_id = field_p->mesh_id;

#ifdef DEBUG_APM
  {
  /* check input */
    //    if(mmr_el_status(mesh_id,Ent_id)==MMC_ACTIVE &&
    if(Ent_nrdof != apr_get_ent_nrdofs(Field_id, Ent_type, Ent_id)){
      printf("Wrong number of dofs in read_ent_dofs !\n");
      exit(-1);
    }
  }
#endif
  
  if(Vect_id<=1) glob_dofs = field_p->dof_ents[Ent_id].vec_dof_1;
  else if(Vect_id==2) glob_dofs = field_p->dof_ents[Ent_id].vec_dof_2;
  else if(Vect_id==3) glob_dofs = field_p->dof_ents[Ent_id].vec_dof_3;
  else field_p->dof_ents[Ent_id].vec_dof_1;

  for(i=0;i<Ent_nrdof;i++) Vect_dofs[i] = glob_dofs[i];

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
  double* Vect_dofs  /* in: dofs to be written */)
{
  /* pointer to field structure */
  apt_field *field_p;
  
  /* auxiliary variables */
  double* glob_dofs;
  int i;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);
  
#ifdef DEBUG_APM
  {
  /* check input */
    if(Ent_nrdof != apr_get_ent_nrdofs(Field_id, Ent_type, Ent_id)){
      printf("Wrong number of dofs in write_ent_dofs !\n");
      exit(-1);
    }
  }
#endif
  
  if(Vect_id<=1) glob_dofs = field_p->dof_ents[Ent_id].vec_dof_1;
  else if(Vect_id==2) glob_dofs = field_p->dof_ents[Ent_id].vec_dof_2;
  else if(Vect_id==3) glob_dofs = field_p->dof_ents[Ent_id].vec_dof_3;
  else field_p->dof_ents[Ent_id].vec_dof_1;

  for(i=0;i<Ent_nrdof;i++) {
      glob_dofs[i] = Vect_dofs[i];
  }

  return(1);
}

/*------------------------------------------------------------
  apr_create_ent_dofs - to write a vector of dofs associated with a given 
                   mesh entity to approximation field data structure
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
  double* glob_dofs;
  int i;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);
  
#ifdef DEBUG_APM
  {
  /* check input */
    if(Ent_nrdof != apr_get_ent_nrdofs(Field_id, Ent_type, Ent_id)){
      printf("Wrong number of dofs in write_ent_dofs !\n");
      exit(-1);
    }
  }
#endif
  
  if(Vect_id<=1) {
    field_p->dof_ents[Ent_id].vec_dof_1 = 
      (double *) malloc(Ent_nrdof*sizeof(double));
    glob_dofs = field_p->dof_ents[Ent_id].vec_dof_1;
  }
  else if(Vect_id==2) {
    field_p->dof_ents[Ent_id].vec_dof_2 = 
      (double *) malloc(Ent_nrdof*sizeof(double));
    glob_dofs = field_p->dof_ents[Ent_id].vec_dof_2;
  }
  else if(Vect_id==3) {
    field_p->dof_ents[Ent_id].vec_dof_3 = 
      (double *) malloc(Ent_nrdof*sizeof(double));
      glob_dofs = field_p->dof_ents[Ent_id].vec_dof_3;
  }
  else  {
    field_p->dof_ents[Ent_id].vec_dof_1 = 
      (double *) malloc(Ent_nrdof*sizeof(double));
      glob_dofs = field_p->dof_ents[Ent_id].vec_dof_1;
  }

  for(i=0;i<Ent_nrdof;i++)  glob_dofs[i] = Vect_dofs[i];

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
  int sol_comp = 0;

  x[0]=1;x[1]=1;x[2]=1;

#ifdef DEBUG_APM
  printf("in set_ini_con - get value %lf from pd function\n",
	 (*Fun_p)(Field_id,x,sol_comp)); 
#endif

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
  int *diagonal,   /* array of indicators whether coefficient matrices are diagonal */
  double *Stiff_mat,	/* out: stiffness matrix stored columnwise */
  double *Rhs_vect	/* out: rhs vector */
  // REMARK:
  //   1. jdofs, jeq - correspond to the solution and go from left to right
  //   2. idofs, ieq - correspond to test functions and go up and down
  //   3. matrices are stored columnwise in vectors
  //   4. solution indices change in rows (go from left to right)
  //   5. test functions indices change in columns (go up and down)
  //   6. when matrices are stored in a vector the index is computed as:
  //      vector[jdofs*nreq*num_dofs+idofs*nreq+jeq*num_dofs+ieq]
  //      (num_dofs=num_shap*nreq)
  //   7. for each pair (idofs,jdofs) there is a small submatrix nreq x nreq
  //      with indices jeq (from left to right) and ieq (from top to bottom)
  //   8. when stiffness matrix entries are computed for vector problems with 
  //      the same shape functions for each component, then for each pair 
  //      (idofs,jdofs) a matrix of coefficients with the size nreq x nreq 
  //      is provided by the problem dependent module
  //   9. when load vector entries are computed for vector problems with 
  //      the same shape functions for each component, then for each idofs index 
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
  int pdeg;		/* degree of polynomial */

/*   static int pdeg_old; /\* indicator for recomputing quadrature data *\/ */
/*   static int base_q;		/\* type of basis functions for quadrilaterals *\/ */
/*   static int ngauss;            /\* number of gauss points *\/ */
/*   static double xg[3000];   	 /\* coordinates of gauss points in 3D *\/ */
/*   static double wg[1000];       /\* gauss weights *\/ */

/* #pragma omp threadprivate (pdeg_old) */
/* #pragma omp threadprivate (base_q) */
/* #pragma omp threadprivate (ngauss) */
/* #pragma omp threadprivate (xg) */
/* #pragma omp threadprivate (wg) */

  // to make old OpenMP compilers working
  int pdeg_old=-1; /* indicator for recomputing quadrature data */
  int base_q;		/* type of basis functions for quadrilaterals */
  int ngauss;            /* number of gauss points */
  double xg[3000];   	 /* coordinates of gauss points in 3D */
  double wg[1000];       /* gauss weights */

  int el_nodes[MMC_MAXELVNO+1];        /* list of nodes of El */
  double node_coor[3*MMC_MAXELVNO];  /* coord of nodes of El */

  int el_mate;		/* element material */
  int nreq;		/* number of equations */
  int num_shap;         /* number of element shape functions */
  double determ;        /* determinant of jacobi matrix */
  double hsize;         /* size of an element */
  double vol;           /* volume for integration rule */
  double xcoor[3];      /* global coord of gauss point */
  double uk_val[APC_MAXEQ]; /* computed solution */
  double uk_x[APC_MAXEQ];   /* gradient of computed solution */
  double uk_y[APC_MAXEQ];   /* gradient of computed solution */
  double uk_z[APC_MAXEQ];   /* gradient of computed solution */
  double un_val[APC_MAXEQ]; /* computed solution */
  double un_x[APC_MAXEQ];   /* gradient of computed solution */
  double un_y[APC_MAXEQ];   /* gradient of computed solution */
  double un_z[APC_MAXEQ];   /* gradient of computed solution */
  double base_phi[APC_MAXELVD];    /* basis functions */
  double base_dphix[APC_MAXELVD];  /* x-derivatives of basis function */
  double base_dphiy[APC_MAXELVD];  /* y-derivatives of basis function */
  double base_dphiz[APC_MAXELVD];  /* y-derivatives of basis function */

  int mesh_id, ndofs;  
  int i, j, k, ki, kk, jdofs, idofs, iaux;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the corresponding mesh */
  mesh_id = apr_get_mesh_id(Field_id);

  // should be
  //if(coeff_ind == 0 || Problem_id != problem_id_old)
  // but because of problems of some OpenMP compilers with static variables is
  {
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
//
//  if(select_coeff==NULL){
//    select_coeff=pdr_select_el_coeff(Problem_id, &mval,
//	 &axx,&axy,&axz,&ayx,&ayy,&ayz,&azx,&azy,&azz,
//	 &bx,&by,&bz,&tx,&ty,&tz,&cval,&lval,&qx,&qy,&qz,&sval);
//  }

  el_mate =  mmr_el_groupID(mesh_id, El_id);
  //printf("mat num %d for el %d\n", el_mate, El_id);

  /* get the coordinates of the nodes of El in the right order */
  mmr_el_node_coor(mesh_id,El_id,el_nodes,node_coor);
  
  /* find degree of polynomial and number of element scalar dofs */
  if(Pdeg_vec!=NULL) pdeg = Pdeg_vec[0];
  else{
#ifdef DEBUG_APM
    /* check input for DG approximation */
    if(mmr_el_status(mesh_id,El_id)!=MMC_ACTIVE){
      printf("Asking for pdeg of inactive element in apr_num_int_el !\n");
      exit(-1);
    }
#endif
      pdeg = apr_get_ent_pdeg(Field_id, APC_ELEMENT, El_id);
  }

  num_shap = apr_get_el_pdeg_numshap(Field_id,El_id,&pdeg);
  nreq=apr_get_nreq(Field_id);
  ndofs = nreq*num_shap;

#ifdef DEBUG_APM
  if(nreq>APC_MAXEQ){
    printf("Number of equations %d greater than the limit %d.\n",nreq,APC_MAXEQ);
    printf("Change APC_MAXEQ in aph_dg_prism.h and recompile the code.\n");
    exit(-1);
  }
  if(num_shap>APC_MAXELVD){
    printf("Number of element shape functions %d greater than the limit %d.\n",
	   num_shap, APC_MAXELVD);
    printf("Change APC_MAXELVD in include/aph_intf.h and recompile the code.\n");
    exit(-1);
  }
  if(ndofs>APC_MAXELSD){
    printf("Number of element degrees of freedom %d greater than the limit %d.\n"
	   , ndofs, APC_MAXELSD);
    printf("Change APC_MAXELSD in include/aph_intf.h and recompile the code.\n");
    exit(-1);
  }
#endif


  /* prepare data for gaussian integration */
  if(pdeg!=pdeg_old){
    base_q=apr_get_base_type(Field_id, El_id);
    apr_set_quadr_3D(base_q, &pdeg, &ngauss, xg, wg);
    pdeg_old = pdeg;
  }

/*kbw
  if(El_id>=0){
    printf("In num_int_el: Field_id %d, mesh_id %d,  element %d\n",
	   Field_id, mesh_id, El_id);
    printf("pdeg %d, ngauss %d\n",pdeg, ngauss);
    printf("NREQ %d, ndof %d, local_dim %d\n",nreq, num_shap, ndofs);
    printf("%d nodes with coordinates:\n",el_nodes[0]);
    for(i=0;i<el_nodes[0];i++){
      printf("node %d (global - %d): x - %f, y - %f, y - %f\n", 
	     i, el_nodes[i+1],
	     node_coor[3*i],node_coor[3*i+1],node_coor[3*i+2]);
    }
    printf("DOFS:\n");
    for(i=0;i<ndofs;i++) printf("%20.15lf",Sol_dofs_n[i]);
    printf("\n");
    getchar();
  } 
/*kew*/


  for (ki=0;ki<ngauss;ki++) {
    
    /* at the gauss point, compute basis functions, determinant etc*/
    iaux = 2; /* calculations with jacobian but not on the boundary */
    determ = apr_elem_calc_3D(iaux, nreq, &pdeg, base_q,
			       &xg[3*ki],node_coor,Sol_dofs_k,
			       base_phi,base_dphix,base_dphiy,base_dphiz,
			       xcoor,uk_val,uk_x,uk_y,uk_z,NULL);
    
    vol = determ * wg[ki];
    
    if(Sol_dofs_n != NULL){

      iaux = 1; /* calculations without jacobian */
      apr_elem_calc_3D(iaux, nreq, &pdeg, base_q,
			       &xg[3*ki],node_coor,Sol_dofs_n,
			       base_phi,base_dphix,base_dphiy,base_dphiz,
			       xcoor,un_val,un_x,un_y,un_z,NULL);
    }

    /* get coefficients of convection-diffusion-reaction equations */
    pdr_el_coeff(Problem_id,El_id,el_mate,hsize,pdeg,&xg[3*ki],
		 base_phi,base_dphix,base_dphiy,base_dphiz,
		 xcoor,uk_val,uk_x,uk_y,uk_z,un_val,un_x,un_y,un_z,
		 mval,axx,axy,axz,ayx,ayy,ayz,azx,azy,azz,
		 bx,by,bz,tx,ty,tz,cval,lval,qx,qy,qz,sval);
    
    /*kbw
    if(El_id>=0){
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
      printf("solution and derivatives: \n");
      printf("un - %lf, der: x - %lf, y - %lf, z - %lf\n",
	     *un_val,*un_x,*un_y,*un_z);
      if(mval!=NULL)
	printf("time coeff : LHS %lf, RHS %lf\n",mval[0],lval[0]);
      if(axx!=NULL)
	printf("diffusion coeff : %lf %lf %lf\n",axx[0],axy[0],axz[0]);
      if(ayy!=NULL)
	printf("diffusion coeff : %lf %lf %lf\n",ayx[0],ayy[0],ayz[0]);
      if(azz!=NULL)
	printf("diffusion coeff : %lf %lf %lf\n",azx[0],azy[0],azz[0]);
      if(bx!=NULL)
	printf("convection coeff: %lf %lf %lf\n",bx[0],by[0],bz[0]);
      if(cval!=NULL)
	printf("reaction coeff  : %lf\n",cval[0]);
      if(qx!=NULL)
	printf("s_xn coeff: %lf %lf %lf\n",qx[0],qy[0],qz[0]);
      if(sval!=NULL)
	printf("source: %lf\n",sval[0]);
      getchar();
    }
      /*kew*/
    

   if(Comp_sm==APC_COMP_SM||Comp_sm==APC_COMP_BOTH){

    if(coeff_vect_ind[1]==1){

      kk=0;
      for (jdofs=0;jdofs<num_shap;jdofs++) {
	for (idofs=0;idofs<num_shap;idofs++) {	
       
	  Stiff_mat[kk+idofs] += 
	                 mval[0] * base_phi[jdofs] * base_phi[idofs] * vol;

	}/* idofs */
	kk+=num_shap;      
      } /* jdofs */
    }


    if(coeff_vect_ind[10]==1){

      kk=0;
      for (jdofs=0;jdofs<num_shap;jdofs++) {
	for (idofs=0;idofs<num_shap;idofs++) {	
       
	  Stiff_mat[kk+idofs] += ( 
                               axx[0] *base_dphix[jdofs] *base_dphix[idofs] +
	                       axy[0] *base_dphiy[jdofs] *base_dphix[idofs] +
	                       axz[0] *base_dphiz[jdofs] *base_dphix[idofs] +
                               ayx[0] *base_dphix[jdofs] *base_dphiy[idofs] +
                               ayy[0] *base_dphiy[jdofs] *base_dphiy[idofs] +
                               ayz[0] *base_dphiz[jdofs] *base_dphiy[idofs] +
                               azx[0] *base_dphix[jdofs] *base_dphiz[idofs] +
                               azy[0] *base_dphiy[jdofs] *base_dphiz[idofs] +
                               azz[0] *base_dphiz[jdofs] *base_dphiz[idofs]  
	                     ) * vol;

	}/* idofs */
	kk+=num_shap;
      } /* jdofs */
    }
    else if(coeff_vect_ind[2]==1){

      kk=0;
      for (jdofs=0;jdofs<num_shap;jdofs++) {
	for (idofs=0;idofs<num_shap;idofs++) {	
	  
	  Stiff_mat[kk+idofs] += (
                               axx[0] *base_dphix[jdofs] *base_dphix[idofs] +
	                       axy[0] *base_dphiy[jdofs] *base_dphix[idofs] +
                               ayx[0] *base_dphix[jdofs] *base_dphiy[idofs] +
                               ayy[0] *base_dphiy[jdofs] *base_dphiy[idofs] 
	                     ) * vol;

	}/* idofs */
	kk+=num_shap;
      } /* jdofs */
    }


    if(coeff_vect_ind[13]==1){

      kk=0;
      for (jdofs=0;jdofs<num_shap;jdofs++) {
	for (idofs=0;idofs<num_shap;idofs++) {	
       
	  Stiff_mat[kk+idofs] += ( bx[0] *base_phi[jdofs] *base_dphix[idofs] +
                               by[0] *base_phi[jdofs] *base_dphiy[idofs] +
                               bz[0] *base_phi[jdofs] *base_dphiz[idofs] ) *vol;
	
	}/* idofs */
	kk+=num_shap;      
      } /* jdofs */
    }
    else if(coeff_vect_ind[11]==1){

      kk=0;
      for (jdofs=0;jdofs<num_shap;jdofs++) {
	for (idofs=0;idofs<num_shap;idofs++) {	
       
	  Stiff_mat[kk+idofs] += ( bx[0] *base_phi[jdofs] *base_dphix[idofs] +
                               by[0] *base_phi[jdofs] *base_dphiy[idofs] ) *vol;

	}/* idofs */
	kk+=num_shap;      
      } /* jdofs */
    }


    // !!!!!!!!!!! Tx, Ty, Tz MISSING !!!!!!!!!!!!!!!!!!!!!!!!

    if(coeff_vect_ind[17]==1){

      kk=0;
      for (jdofs=0;jdofs<num_shap;jdofs++) {
	for (idofs=0;idofs<num_shap;idofs++) {	
       
	  Stiff_mat[kk+idofs] += cval[0] *base_phi[jdofs] *base_phi[idofs] *vol;
	
	}/* idofs */
	kk+=num_shap;      
      } /* jdofs */
    }
    
   } /* end if computing SM */

   if(Comp_sm==APC_COMP_RHS||Comp_sm==APC_COMP_BOTH){

    if(coeff_vect_ind[18]==1){

      kk=0;
      for (idofs=0;idofs<num_shap;idofs++) {
      
	Rhs_vect[kk] += lval[0] * base_phi[idofs] * vol; 
	
	kk++;      
      }/* idofs */    
    }

  
    if(coeff_vect_ind[21]==1){

      kk=0;
      for (idofs=0;idofs<num_shap;idofs++) {
      
	Rhs_vect[kk] += ( qx[0] * base_dphix[idofs] +
		          qy[0] * base_dphiy[idofs] +
		          qz[0] * base_dphiz[idofs] ) * vol;
	    
	kk++;      
      }/* idofs */    
    } 
    else if(coeff_vect_ind[19]==1){

      kk=0;
      for (idofs=0;idofs<num_shap;idofs++) {
      
	Rhs_vect[kk] += ( qx[0] * base_dphix[idofs] +
		          qy[0] * base_dphiy[idofs] ) * vol;

	kk++;      
      }/* idofs */
    }

  
    if(coeff_vect_ind[22]==1){

      kk=0;
      for (idofs=0;idofs<num_shap;idofs++) {
      
	Rhs_vect[kk] += sval[0] * base_phi[idofs] * vol;

	kk++;      
      }/* idofs */
    } 
  
   } /* end if computing RHSV */
    
  } /* end loop over integration points: ki */
  
  return(1);
}



/*------------------------------------------------------------
  apr_num_int_fa - to perform numerical integration for a face
------------------------------------------------------------*/
int apr_num_int_fa(
  )
{


  return(1);
}

/*------------------------------------------------------------
  apr_proj_dof_ref - to rewrite dofs after modifying the mesh 
------------------------------------------------------------*/
int apr_proj_dof_ref(
  int Field_id,    /* in: approximation field ID  */
  int El,	   /* in:  >0 - rewrite after one [de]refinement of el */
                   /*     <=0 - rewrite after massive [de]refinements */
  int Max_elem_id, /* in: maximal element (face, etc.) id */
  int Max_face_id, /*     before and after refinements */
  int Max_edge_id,
  int Max_vert_id
  )
{

  /* pointer to field structure */
  apt_field *field_p;
  int i, iaux, idofs, mesh_id, nel, nmel, ref_type, nreq;
  int el_fath, pdeg_fath, pmax, ison, son, numshap, num_dof_scal;
  double dofs_fath[APC_MAXELSD]; /* father solution dofs */
  double dofs_son[MMC_MAXELSONS*APC_MAXELSD]; /* sons' solution dofs */
  int el_sons[MMC_MAXELSONS+1], pdeg_from[MMC_MAXELSONS+1];
/*++++++++++++++++ executable statements ++++++++++++++++*/
  
  /* select the pointer to the approximation field */
  field_p = apr_select_field(Field_id);
  
  /* select the corresponding mesh */
  mesh_id = apr_get_mesh_id(Field_id);

  /* get the number of equations - solution vector components */
  nreq = apr_get_nreq(Field_id);

  /* in a loop over active elements */
  nel=0;
  while((nel=mmr_get_next_act_elem(mesh_id,nel))!=0){


    /* check whether element results from refinement */
    el_fath = mmr_el_fam(mesh_id, nel, NULL, NULL);
    if(el_fath!=MMC_NO_FATH) pdeg_fath = field_p->dof_ents[el_fath].pdeg;

/* 
 * refinement algorithm leaves pdeg indicator unchanged for the father
 * element - this inconsistency (pdeg should be set to APC_NO_DOFS) indicates
 * that dofs should be allocated and filled for the sons
*/

    /* if father exists and still has dofs */
    if(el_fath!=MMC_NO_FATH && pdeg_fath!=APC_NO_DOFS){

#ifdef DEBUG_APM
      if(mmr_el_status(mesh_id, el_fath) != MMC_INACTIVE){
	printf("active father %d for son %d in rewr_dofs\n", el_fath, nel);
	exit(-1);
      }
#endif

      /* project solution from father element to its son */

      /* get the number of dofs */
      numshap=apr_get_el_pdeg_numshap(Field_id,el_fath,&pdeg_fath);
      num_dof_scal=nreq*numshap;

      /* set element's pdeg */
      field_p->dof_ents[nel].pdeg = pdeg_fath;

      /* allocate space for dofs */
      field_p->dof_ents[nel].vec_dof_1 = 
	(double *) malloc(num_dof_scal*sizeof(double));

      /* rewrite dofs from global data structure */
      for(i=0;i<num_dof_scal;i++)
	dofs_fath[i]=field_p->dof_ents[el_fath].vec_dof_1[i];

/*kbw
      {
	int ino;
	printf("found element %d with no dofs, el_fath %d, pdeg %d, numshap %d (%d, nr_sol %d)\n",
	       nel,el_fath,pdeg_fath,numshap,num_dof_scal,field_p->nr_sol);
	printf("father dofs_1 (address %x): \n", dofs_fath);
	for(ino=0;ino<num_dof_scal;ino++) printf("%20.15lf",dofs_fath[ino]);
	printf("\n");
	iaux=1;
	apr_read_ent_dofs(Field_id,APC_ELEMENT,nel,num_dof_scal,iaux,
			   dofs_son);
	printf("element's dofs_1 (address %x): \n", dofs_son);
	for(ino=0;ino<num_dof_scal;ino++) printf("%20.15lf",dofs_son[ino]);
	printf("\n");
	printf("APC_MAXELVD %d, MMC_MAXELVNO %d, APC_MAXELSD %d, APC_MAXEQ %d\n",
	       APC_MAXELVD, MMC_MAXELVNO, APC_MAXELSD, APC_MAXEQ);
	
      }
/*kew*/
      
      /* perform L2 projection to set dofs */
      iaux = -1;
      apr_L2_proj(Field_id,iaux,nel,&pdeg_fath,dofs_son,
		  &el_fath,&pdeg_fath,dofs_fath,NULL);
      
      
      iaux=1;
      apr_write_ent_dofs(Field_id,APC_ELEMENT,nel,num_dof_scal,iaux,
			 dofs_son);
      

/*kbw
      {
	int ino;
	printf("found element %d with no dofs, el_fath %d, pdeg %d, numshap %d (%d, nr_sol %d)\n",
	       nel,el_fath,pdeg_fath,numshap,num_dof_scal,field_p->nr_sol);
	printf("father dofs_1: \n");
	for(ino=0;ino<num_dof_scal;ino++) printf("%20.15lf",dofs_fath[ino]);
	printf("\n");
	iaux=1;
	apr_read_ent_dofs(Field_id,APC_ELEMENT,nel,num_dof_scal,iaux,
			   dofs_fath);
	printf("element's dofs_1: \n");
	for(ino=0;ino<num_dof_scal;ino++) printf("%20.15lf",dofs_fath[ino]);
	printf("\n");
	printf("APC_MAXELVD %d, MMC_MAXELVNO %d, APC_MAXELSD %d, APC_MAXEQ %d\n",
	       APC_MAXELVD, MMC_MAXELVNO, APC_MAXELSD, APC_MAXEQ);

      }
/*kew*/

      if(field_p->nr_sol>1){
	field_p->dof_ents[nel].vec_dof_2 = 
	  (double *) malloc(num_dof_scal*sizeof(double));
	/* rewrite dofs from global data structure */
	for(i=0;i<num_dof_scal;i++)
	  dofs_fath[i]=field_p->dof_ents[el_fath].vec_dof_2[i];
	/* perform L2 projection to set dofs */
	iaux = -1;
	apr_L2_proj(Field_id,iaux,nel,&pdeg_fath,dofs_son,
		    &el_fath,&pdeg_fath,dofs_fath,NULL);	
	/* rewrite dofs to global data structure */
	iaux=2;
	apr_write_ent_dofs(Field_id,APC_ELEMENT,nel,num_dof_scal,iaux,
			   dofs_son);
      }	  
      
      
      if(field_p->nr_sol>2){
	field_p->dof_ents[nel].vec_dof_3 = 
	  (double *) malloc(num_dof_scal*sizeof(double));
	/* rewrite dofs from global data structure */
	for(i=0;i<num_dof_scal;i++)
	  dofs_fath[i]=field_p->dof_ents[el_fath].vec_dof_3[i];
	iaux = -1;
	apr_L2_proj(Field_id,iaux,nel,&pdeg_fath,dofs_son,
		    &el_fath,&pdeg_fath,dofs_fath,NULL);	
	/* rewrite dofs to global data structure */
	iaux=3;
	apr_write_ent_dofs(Field_id,APC_ELEMENT,nel,num_dof_scal,iaux,
			   dofs_son);
      }	  

      /* free father dofs - postponed to the final check due to other sons
      free(field_p->dof_ents[el_fath].vec_dof_1);
      field_p->dof_ents[el_fath].vec_dof_1 = NULL;
      if(field_p->nr_sol>1) {
	free(field_p->dof_ents[el_fath].vec_dof_2);
	field_p->dof_ents[el_fath].vec_dof_2 = NULL;
      }
      if(field_p->nr_sol>2) {
	free(field_p->dof_ents[el_fath].vec_dof_3);
	field_p->dof_ents[el_fath].vec_dof_3 = NULL;
      }
      field_p->dof_ents[el_fath].pdeg = APC_NO_DOFS;
      */

    } /* if element results from refinement and still his father has dofs */

    /* check whether element results from derefinement */
    ref_type = mmr_el_type_ref(mesh_id,nel);

/*
 * derefinement algorithm leaves ref_type unchanged for the father element
 * but indicates that its active again - this inconsistency indicates that
 * dofs should be allocated and filled for the father (ref_type tells also
 * about the number of sons, which is important for rewriting dofs)
 * (the sons' pdeg indicators can also be used here - active element has sons
 * that have dofs allocated - inconsistency -> rewrite dofs, but still
 * ref_type is necessary to indicate the number of sons, so why not using it)
*/
    /* if element is active but ref_type indicates it was refined */
    if(ref_type != MMC_NOT_REF){

      mmr_el_fam(mesh_id, nel, el_sons, NULL);
      
      pmax=0;num_dof_scal=0;	
      /* for all element's sons */ 
      for(ison=1;ison<=el_sons[0];ison++){
	son=el_sons[ison];
	pdeg_from[ison] = field_p->dof_ents[son].pdeg;

#ifdef DEBUG_APM
	if(pdeg_from[ison]<0){
	  printf("Wrong pdeg %d in rewriting dofs for element %d\n",
		 pdeg_from[ison],son);
	  exit(-1);
	}
#endif

	/* find the greatest degree of polynomial - simplified version!*/
	numshap=apr_get_el_pdeg_numshap(Field_id,son,&pdeg_from[ison]);
	if(num_dof_scal<nreq*numshap) num_dof_scal=nreq*numshap;
	if(pmax<pdeg_from[ison]) pmax=pdeg_from[ison]; 
	
/*kbw
	printf("son %d (%d), pdeg %d, pdeg_max %d, num_dof_scal %d\n",
	       ison,el_sons[ison],pdeg_from[ison],pmax,num_dof_scal); 
/*kew*/

      }

      field_p->dof_ents[nel].pdeg=pmax;

      idofs=0;
      /* for each derefined element */
      for(ison=1;ison<=el_sons[0];ison++){
	son=el_sons[ison];
	
	  
	num_dof_scal=nreq*apr_get_el_pdeg_numshap(Field_id,son,&pdeg_from[ison]);
	
#ifdef DEBUG_APM
	if(mmr_el_status(mesh_id,son)!=MMC_FREE){
	  printf("Inconsistency in derefinement: not free space %d after clustering!\n",
		 son);
	  exit(-1);
	}
#endif

	/* rewrite son's dofs to local array */
	for(i=0;i<num_dof_scal;i++){
	  dofs_son[idofs+i]=field_p->dof_ents[son].vec_dof_1[i];
	}
	
/*kbw
	  printf("son %d (%d), status %d, pdeg %d, dofs:\n", son, ison, 
		 mmr_el_status(mesh_id,son), field_p->dof_ents[son].pdeg);
	  for(i=0;i<num_dof_scal;i++){
	    printf("%20.15lf",dofs_son[idofs+i]);
	  }
	  printf("\n");
/*kew*/

/*kbw
	printf("son %d (%d) dofs:\n",ison,el_sons[ison]);
	for(i=0;i< num_dof_scal;i++) printf("%20.15lf",dofs_son[idofs+i]);
 	printf("\n");
/*kew*/

	/* update counter for dofs */
	idofs += num_dof_scal;
	
      }
      
      /* perform L2 projection to set dofs for the father */
      iaux= -el_sons[0];

/*kbw
      if(nel==255||nel==145||nel==178){
	printf("Input to l2_proj: iaux %d, iel %d, pmax %d\n",
	       iaux, nel, pmax);
	for(ison=1;ison<=el_sons[0];ison++){
	  son=el_sons[ison];
	  printf("son %d (%d), pdeg %d\n", son, ison, pdeg_from[ison]);
	}
	for(i=0;i<idofs;i++){
	  printf("%20.15lf", dofs_son[i]);
	}
	printf("\n");
      }
/*kew*/

      apr_L2_proj(Field_id, iaux, nel, &pmax, dofs_fath,
		  &el_sons[1], &pdeg_from[1], dofs_son, NULL);
      
      
      num_dof_scal= nreq*apr_get_el_pdeg_numshap(Field_id, el_sons[1], &pmax);
      
/*kbw
      if(nel==255||nel==145||nel==178){
	printf("Output from l2_proj: iaux %d, iel %d, pmax %d, dofs_fath:\n",
	       iaux, nel, pmax);
	for(i=0;i<num_dof_scal;i++){
	  printf("%20.15lf", dofs_fath[i]);
	}
	printf("\n");
      }
/*kew*/

/*kbw
	printf("projected dofs:\n");
	for(i=0;i< num_dof_scal;i++) printf("%20.15lf",dofs_fath[i]);
 	printf("\n");
/*kew*/

      field_p->dof_ents[nel].vec_dof_1=
	(double *) malloc(num_dof_scal*sizeof(double));
      
      /* write projected dofs to data structure */
      for(i=0;i<num_dof_scal;i++) 
	field_p->dof_ents[nel].vec_dof_1[i]=dofs_fath[i];
      
      
/*kbw
	printf("dofs read from data structure:\n");
	iaux=1;
	apr_read_ent_dofs(Field_id,APC_ELEMENT,nel,num_dof_scal,iaux,
			    dofs_fath);
	for(i=0;i< num_dof_scal;i++) printf("%20.15lf",dofs_fath[i]);
 	printf("\n");
/*kew*/

      if(field_p->nr_sol>1){
	
	idofs=0;
	/* for each derefined element */
	for(ison=1;ison<=el_sons[0];ison++){
	  son=el_sons[ison];

	  num_dof_scal=nreq*apr_get_el_pdeg_numshap(Field_id, son,
						     &pdeg_from[ison]);
	  
	  /* rewrite son's dofs to local array */
	  for(i=0;i<num_dof_scal;i++){
	    dofs_son[idofs+i]=field_p->dof_ents[son].vec_dof_2[i];
	  }
	
	  
	  /* update counter for dofs */
	  idofs += num_dof_scal;
	  
	}
	
	/* perform L2 projection to set dofs */
	iaux= -el_sons[0];
	apr_L2_proj(Field_id, iaux, nel, &pmax, dofs_fath,
		    &el_sons[1], &pdeg_from[1], dofs_son, NULL);
	
	
	num_dof_scal= nreq*apr_get_el_pdeg_numshap(Field_id, nel, &pmax);
	
	field_p->dof_ents[nel].vec_dof_2=
	  (double *) malloc(num_dof_scal*sizeof(double));
	
	/* write projected dofs to data structure */
	for(i=0;i<num_dof_scal;i++) 
	  field_p->dof_ents[nel].vec_dof_2[i]=dofs_fath[i];
	
      }
      
      if(field_p->nr_sol>2){
	
	idofs=0;
	/* for each derefined element */
	for(ison=1;ison<=el_sons[0];ison++){
	  son=el_sons[ison];
	  
	  num_dof_scal=nreq*apr_get_el_pdeg_numshap(Field_id,son,&pdeg_from[ison]);
	  
	  /* rewrite son's dofs to local array */
	  for(i=0;i<num_dof_scal;i++){
	    dofs_son[idofs+i]=field_p->dof_ents[son].vec_dof_3[i];
	  }
		  
	  /* update counter for dofs */
	  idofs += num_dof_scal;

	}
	
	/* perform L2 projection to set dofs */
	iaux= -el_sons[0];
	apr_L2_proj(Field_id, iaux, nel, &pmax, dofs_fath,
		    &el_sons[1], &pdeg_from[1], dofs_son, NULL);
	
	
	num_dof_scal= nreq*apr_get_el_pdeg_numshap(Field_id, nel, &pmax);
	
	field_p->dof_ents[nel].vec_dof_3=
	  (double *) malloc(num_dof_scal*sizeof(double));
	
	/* write projected dofs to data structure */
	for(i=0;i<num_dof_scal;i++) 
	  field_p->dof_ents[nel].vec_dof_3[i]=dofs_fath[i];
	
      }

    } /* end if element results from refinement */


  } /* end loop over active elements: nel */

  
  /* loop over dofs data structure including free spaces */
  nmel=mmr_get_max_elem_id(mesh_id);
  for(nel=1;nel<=nmel;nel++){

    /* if element inactive or free space */
    if(mmr_el_status(mesh_id, nel) != MMC_ACTIVE){

      /* clean the space if it has dofs */
      if(field_p->dof_ents[nel].pdeg != APC_NO_DOFS){

/*kbw
	printf("found inactive element %d with dofs\n",nel);
	getchar();
/*kew*/

	free(field_p->dof_ents[nel].vec_dof_1);
	field_p->dof_ents[nel].vec_dof_1 = NULL;
	if(field_p->nr_sol>1) {
	  free(field_p->dof_ents[nel].vec_dof_2);
	  field_p->dof_ents[nel].vec_dof_2 = NULL;
	}
	if(field_p->nr_sol>2) {
	  free(field_p->dof_ents[nel].vec_dof_3);
	  field_p->dof_ents[nel].vec_dof_3 = NULL;
	}
      }
      /* set pdeg */
      field_p->dof_ents[nel].pdeg = APC_NO_DOFS;

    } /* end if element inactive */

  }

  return(1);
}


/*------------------------------------------------------------
  apr_rewr_sol - to rewrite solution from one vector to another
------------------------------------------------------------*/
extern int apr_rewr_sol(
  int Field_id,      /* in: data structure to be used  */
  int Sol_from,      /* in: ID of vector to read solution from */
  int Sol_to         /* in: ID of vector to write solution to */
  )
{

  int nel, mesh_id, num_shap;
  double dofs_loc[APC_MAXELSD];

/*++++++++++++++++ executable statements ++++++++++++++++*/

  mesh_id = apr_get_mesh_id(Field_id);

  nel=0;
  while((nel=mmr_get_next_act_elem(mesh_id,nel))!=0){
    
    num_shap = apr_get_ent_numshap(Field_id, APC_ELEMENT, nel);
    apr_read_ent_dofs(Field_id, APC_ELEMENT, nel, num_shap, Sol_from, dofs_loc);
    apr_write_ent_dofs(Field_id, APC_ELEMENT, nel, num_shap, Sol_to, dofs_loc);

  }

  return(0);
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
  int i, nmel, mesh_id;

/*++++++++++++++++ executable statements ++++++++++++++++*/
  
  /* select the proper approximation field */
  field_p = apr_select_field(Field_id);
  
  /* select the corresponding mesh */
  mesh_id = apr_get_mesh_id(Field_id);

  /* loop over dofs data structure including free spaces */
  nmel=mmr_get_max_elem_id(mesh_id);
  for(i=1;i<=nmel;i++){
    /* clean the space if it has dofs */
    if(field_p->dof_ents[i].pdeg != APC_NO_DOFS){
      free(field_p->dof_ents[i].vec_dof_1);
      if(field_p->nr_sol>1) free(field_p->dof_ents[i].vec_dof_2); 
      if(field_p->nr_sol>2) free(field_p->dof_ents[i].vec_dof_3); 
    }
  }
  free(field_p->dof_ents);

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

  int mesh_id, father, son, ison, max_gen_diff;
  int elsons[MMC_MAXELSONS+1]; /* family information */
  int ifa, face, son_faces[MMC_MAXELFAC+1];


/*++++++++++++++++ executable statements ++++++++++++++++*/

  mesh_id = apr_get_mesh_id(Field_id);
  max_gen_diff=mmr_get_max_gen_diff(mesh_id);

  if(max_gen_diff== 1){
	      
    father=mmr_el_fam(mesh_id,El_id,NULL,NULL);
    mmr_el_fam(mesh_id,father,elsons,NULL);

    /* check generation difference is not exceeded */
    /* check whether derefinement will not cause too big generation difference */
    /* simple test for max_gen_diff==1 */

    for(ison=1;ison<=elsons[0];ison++){
      
      son = elsons[ison];
      
      mmr_el_faces(mesh_id, son, son_faces, NULL);
      
      
      for(ifa=1;ifa<=son_faces[0];ifa++){
	
	face=son_faces[ifa];
	
/*kbw
	if(!abs(face)==669||abs(face)==739){
	  printf("father %d son %d, face %d, status %d, type %d, bc %d\n",
	  father, son, face, mmr_fa_status(mesh_id,face),
	  mmr_fa_type(mesh_id,face), mmr_fa_bc(mesh_id,face));
	}
/*kew*/
	
	
	if(mmr_fa_status(mesh_id,face)==MMC_INACTIVE){
	  
/*kbw
  printf("In derefinements: too big gen diff for element %d across the face %d!\n", son,face);
/*kew*/

	  return(APC_DEREF_DENIED);
	      
	}
      }
    }
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

/*kbw
printf("apr_limit_ref: Field_id %d, mesh_id %d, max_gen %d (checked %d)\n",
       Field_id, mesh_id, max_gen, mmr_get_max_gen(mesh_id));
/*kew*/

  if(max_gen==0) return(APC_REF_DENIED);
  else if(max_gen>0){
    
    /* check whether max generation is not exceeded */
    gen_el=mmr_el_gen(mesh_id,El_id);

/*kbw
printf("apr_limit_ref: element %d to refine: gen %d <> maxgen %d\n",
       El_id, gen_el,max_gen);
/*kew*/

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

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper mesh data structure */
  int mesh_id = apr_get_mesh_id(Field_id);

  //i=21; max_gen=mmr_get_mesh_i_params(mesh,i);
  max_gen=mmr_get_max_gen(mesh_id);
  //i=22; max_gen_diff=mmr_get_mesh_i_params(mesh,i);
  max_gen_diff=mmr_get_max_gen_diff(mesh_id);
  
  if(El!=MMC_DO_UNI_REF&&El>0) {

/*set number of elements waiting for refinement (they will be stored in listwait)*/
    nrwait=0;

/*element to refine*/
    nelref=El;

    beginning:{}

/* check whether element is active */
    if(mmr_el_status(mesh_id, nelref)!=MMC_ACTIVE){
#ifdef DEBUG_MMM
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

#ifdef DEBUG_MMM
      if(iprint>MMC_PRINT_INFO) {
	printf("element %d generation (%d) > max_gen (%d) in divide8_p!\n", 
		nelref,gen_el,max_gen);
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
    else if(max_gen_diff > 0){

 /* max_gen_diff >0 - regularity based on generation difference */

      int face_eq_neig[2];
      int el_faces[7];

      el_type=mmr_el_type(mesh_id,nelref);
      mmr_el_faces(mesh_id, nelref, el_faces, NULL);

      num_face = el_faces[0];
      for(ifa=0;ifa<num_face;ifa++){
	
	face=abs(el_faces[ifa+1]);
	mmr_fa_eq_neig(mesh_id, face, face_eq_neig, NULL, NULL); 
	
	if(abs(face_eq_neig[0])==nelref) ifaneig=1;
	else ifaneig=0;
	
	if(face_eq_neig[ifaneig]==MMC_SUB_BND){
	  
#ifdef DEBUG_MMM
	  if(iprint>MMC_PRINT_INFO){
	    printf("element %d: refinements may go outside because of DIFF_GEN\n",
		   nelref);
	  }
#endif
	  
	}
	else if(face_eq_neig[ifaneig]==-1){
	  
	  int face_neig[2];
	  
	  mmr_fa_neig(mesh_id, face, face_neig, NULL, NULL, NULL, NULL, NULL); 
	  
	  iaux = gen_el - mmr_el_gen(mesh_id,abs(face_neig[ifaneig]));
	  if(iaux>=max_gen_diff){ 
	    
	    listwait[nrwait]=nelref;
	    nrwait++;
	    nelref= abs(face_neig[ifaneig]);
	    
	    
#ifdef DEBUG_MMM
	    if(iprint>MMC_PRINT_INFO){
	      printf("adding element %d to refine because of DIFF_GEN\n",nelref);
	    }
#endif
	    
	    goto beginning;
	    
	  }
	}
      }
    } /* end if max_gen_diff >0 - regularity based on generation difference */


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
    iaux=mmr_refine_mesh(mesh_id);
    return(iaux);
  }
  else{
    printf("Wrong parameter %d in apr_refine! Exiting.\n"
           , El);
    exit(-1);
  }


#ifdef DEBUG_MMM
  printf("Error 230956 in mmr_refine!\n");
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
    iaux=mmr_derefine_el(mesh_id,El);
  }
  else if(El==MMC_DO_UNI_DEREF){
    iaux=mmr_derefine_mesh(mesh_id);
  }
  else{
    printf("Wrong parameter %d in apr_derefine! Exiting.\n"
           , El);
    exit(-1);
  }

  return(iaux);
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
  int* List_dof_ent_type, /* out: list of no of dofs for 'dof' entity */
  int* List_dof_ent_id,   /* out: list of no of dofs for 'dof' entity */
  int* List_dof_ent_nrdofs,/* out: list of no of dofs for 'dof' entity */
  /* below: optional means when Stiff_mat and Rhs_vect are computed */
  int* Nrdofs_loc,/* in(optional): size of Stiff_mat and Rhs_vect */
                /* out(optional): actual number of dofs per integration entity*/
  /* for matrices passing NULL is possible if suitable Comp_sm is passed */
  double* Stiff_mat,      /* out(optional): stiffness matrix stored columnwise*/
  double* Rhs_vect        /* out(optional): rhs vector */
			     )
{

  printf("apr_get_stiff_mat_data not yet implemented for DG!\n");
  exit(-1);
  
}

/*------------------------------------------------------------
  apr_get_constr_data - to return constraints data for a node (dof entity) 
------------------------------------------------------------*/
int apr_get_constr_data( /* returns: >=0 - success code, <0 - error code */
  int Field_id,     /* in: approximation field ID  */
  int Node_id,      /* in: id of node (can be vertex or edge) */
  int Node_type,    /* in: type of input node APC_VERTEX or APC_EDGE */

  int *Constr,      /* out: table of constraints data, 
		     *      Constr[0] - number of constraints 
		     *                  (in linear approximation
		     *                  2 - mid-edge node, 4 - mid-side node)
		     */
  int *Constr_type  /* out: table of constraints element type,
		     *      Constr_type[0] - id of coefficient vector
		     */

  )
{
  printf("There are no constrained nodes for DG!\n");
  exit(-1);
  
}

