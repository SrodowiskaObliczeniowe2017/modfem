/****************************************************************
File sis_mkb_intf.c - implementation of the interface module
   between iterative solvers and the finite element
   code (the module forms part of the fem code)
file contains: definition of parameters,
   data types, global variables and external functions)


Contains definitions of interface routines:
  sir_module_introduce - to return the solver name
  sir_solve_lin_sys - to perform the five steps below in one call
  sir_init - to create a new solver instance and read its control parameters
  sir_create - to create and initialize solver data structure
  sir_solve - to solve the system for a given data
  sir_free - to free memory for a given solver instance
  sir_destroy - to destroy a solver instance
  sir_assemble_local_stiff_mat_with_table - to assemble an element stiffness matrix
                                   to the global SM using assembly table
  sir_assemble_local_stiff_mat - to assemble an element stiffness matrix
                                   to the global SM

and auxiliary procedures:
  sir_mkb_show_blocks_matrix - save sm_matrix to file

  sir_reverse_cuthill_mckee - bandwidth reduction algorithm (renumbering)


The module uses several notions associated with degrees of freedom (DOF)
1. DOF entity - mesh entity (such as vertex, edge, face, element interior)
                with which DOFs are associated in FEM approximation module
2. DOF structure - internal structure of sid_mkb for renumbering of DOFs
                   and managing the structure of global stiffness matrix
3. DOF block - block of DOFs corresponding usually to a single DOF structure
               and a single DOF entity, used in the interface with the solver
               module
DOF entities are identified by their type and ID (obtained from FEM code)
DOF structures are stored in an array of structures (offset 0)
The interface with the solver module is based on simple arrays of integers
or doubles and DOF blocks' IDs (numbered from 0) are used as indices for these arrays
 
------------------------------  			
History:        
	02.2002 - Krzysztof Banas, initial version
        2015 - Krzysztof Banas, Kazimierz Chlon - extensions for renumbering 
               and assembly tables
**************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>

//#include "../amg_mkb/lad_petsc/GhostBlockExchanger.h"
/* problem dependent interface between approximation and solver modules */
#include "pdh_intf.h"

/* interface for general purpose linear solver modules */
#include "sih_intf.h" 

/* provided interface of the mkb solver adapter module */
#include "../lsd_mkb/lsh_mkb_intf.h"

#ifdef TIME_TEST
#define MKB_SOLVE_TIME_TEST
#endif

#define MKB_SOLVE_TIME_TEST

//#define PRINT_MATRIX
#include "uth_system.h"

#ifdef PARALLEL
  #undef RENUMBERING
  #undef INTERNAL_RENUMBERING
  #define RENUMBERING_FOR_PARALLEL
#else
  #define RENUMBERING
  //#define INTERNAL_RENUMBERING
#endif

  #undef RENUMBERING
  #undef INTERNAL_RENUMBERING

/* internal information for the mkb solver interface module */
#include "./sih_mkb.h"

/*** CONSTANTS ***/
/* Solver (preconditioner) types */
const int SIC_SINGLE_LEVEL = 1;
const int SIC_TWO_LEVEL    = 2;
const int SIC_MULTI_LEVEL  = 0;


/*** GLOBAL VARIABLES (for the solver module) ***/
int   siv_nr_solvers = 0;    /* the number of solvers in the problem */
int   siv_cur_solver_id = -1;                /* ID of the current problem */
sit_solvers siv_solver[SIC_MAX_NUM_SOLV];        /* array of solvers */


const int FULL_MULTIGRID = 0;
const int GALERKIN_PROJ = 0; 

/* declarations of local utility procedures */
int sir_put_list( /* returns*/
		/*  >0 - position already occupied on the list */
             	/*  <0 - position at which put on the list */
            	/*   0 - list full, not found on the list */
	int Num, 	/* in: number to put on the list */
	int* List, 	/* in: list */
	int Ll		/* in: total list's lengths */
	);
int sir_sort_short(
		   int* A, // array
		   int p,  //first index
		   int k   // last index
		   );

#ifdef INTERNAL_RENUMBERING
/*------------------------------------------------------------
  sir_reverse_cuthill_mckee - bandwidth reduction algorithm (renumbering)
    - we keep sir_... version in case utr_renumber takes too much time
------------------------------------------------------------*/
void sir_reverse_cuthill_mckee(sit_dof_struct * L_dof_struct, int Nr_dof_ent);
#endif

/*------------------------------------------------------------
  sir_mkb_show_blocks_matrix - save sm_matrix to file
------------------------------------------------------------*/
int sir_mkb_show_blocks_matrix(
                         /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id          /* in: level ID */
			       );


/********************* INTERFACE PROCEDURES *******************/

/*------------------------------------------------------------
  sir_module_introduce - to return the solver name
------------------------------------------------------------*/
int sir_module_introduce( /* returns: >=0 - success code, <0 - error code */
  char* Solver_name  /* out: the name of the solver */
  )
{
  char* string = "MKB";

  strcpy(Solver_name,string);

  return(1);
}

/*------------------------------------------------------------
  sir_solve_lin_sys - to solve the system of linear equations for the current
          setting of parameters (with initialization and finalization phases)
------------------------------------------------------------*/
int sir_solve_lin_sys( /* returns: >=0 - success code, <0 - error code */
  int Problem_id,    /* ID of the problem associated with the solver */
  int Solver_type, /* type of solver: 0 - direct, >0 - iterative (passed to mkb_core) */
  int Parallel,      /* parameter specifying sequential (SIC_SEQUENTIAL) */
                     /* or parallel (SIC_PARALLEL) execution */
  int Nr_levels, // maximal number of levels for multigrid
  char* Filename,  /* in: name of the file with control parameters */
  int Max_iter, /* maximal number of iterations, -1 for values from Filename */
  int Error_type, /* type of error norm (stopping criterion), -1 for Filename*/
  double Error_tolerance, /* value for stopping criterion, -1.0 for Filename */ 
  int Monitoring_level /* Level of output, -1 for Filename */
  )
{


  int comp_type, solver_id, monitor, ini_guess;
  char solver_name[100];
  double conv_rate;

/*++++++++++++++++ executable statements ++++++++++++++++*/


  /*kbw
  printf("In sir_solve_lin_sys before init: problem_id %d\n",
	 problem_id);
  /*kew*/

  /* initialize the solver */
  solver_id = sir_init(Solver_type, Parallel, Nr_levels,
		       Filename, Max_iter, Error_type, Error_tolerance,
		       Monitoring_level);
  siv_cur_solver_id = solver_id;
  /*kbw
  printf("In sir_solve_lin_sys before create: solver_id %d, problem_id %d\n",
	 solver_id, problem_id);
  /*kew*/

  /* create the solver data structure and associate it with a given problem */
  sir_create(solver_id, Problem_id);

  /*kbw
  printf("In sir_solve_lin_sys before solve: solver_id %d, problem_id %d\n",
	 solver_id, problem_id);
  /*kew*/

  /* launch the solver */
  comp_type = SIC_SOLVE;
  //monitor = SIC_PRINT_INFO;
  //monitor = SIC_PRINT_NOT;
  ini_guess = 0; // do not get initial guess
  sir_solve(solver_id, comp_type, ini_guess, Monitoring_level, 
	    &Max_iter, &Error_tolerance, &conv_rate);

  /*kbw
  printf("In sir_solve_lin_sys before free: solver_id %d, problem_id %d\n",
	 solver_id, problem_id);
  /*kew*/

  /* free the solver data structure - together with renumbering data */
  sir_free(solver_id);
  
  /*kbw
  printf("In sir_solve_lin_sys before destroy: solver_id %d, problem_id %d\n",
	 solver_id, problem_id);
  /*kew*/

  /* destroy the solver */
  sir_destroy(solver_id);

  return(0);
}

/*------------------------------------------------------------
  sir_init - to create a new solver, read its control parameters
             and initialize its data structure
------------------------------------------------------------*/
int sir_init( /* returns: >0 - solver ID, <0 - error code */
  int Solver_type, /* type of solver: 0 - direct, >0 - iterative (passed to mkb_core) */
  int Parallel,      /* parameter specifying sequential (SIC_SEQUENTIAL) */
                     /* or parallel (SIC_PARALLEL) execution */
  int Nr_levels, // maximal number of levels for multigrid
  char* Filename,  /* in: name of the file with control parameters */
  int Max_iter, /* maximal number of iterations, -1 for values from Filename */
  int Error_type, /* type of error norm (stopping criterion), -1 for Filename*/
  double Error_tolerance, /* value for stopping criterion, -1.0 for Filename */ 
  int Monitoring_level /* Level of output, -1 for Filename */
  )
{

/*++++++++++++++++ executable statements ++++++++++++++++*/

  int info;

/*ok_kbw*/
  printf("\nEntering sir_init: Solver_type %d, Parallel %d, Nr_levels %d\n",
	 Solver_type, Parallel, Nr_levels);
  printf("Filename %s\n", Filename);
  printf("Max_iter %d, Error_type %d, Error_tolerance %lf, Monitor %d\n",
	 Max_iter, Error_type, Error_tolerance, Monitoring_level);
//*kew*/

  /* check constants compatibility */
  if(SIC_SEQUENTIAL != LSC_SEQUENTIAL || SIC_PARALLEL != LSC_PARALLEL){
    printf("Wrong constants in sir_init\n"); exit(0);
  }    
  if(SIC_SOLVE != LSC_SOLVE || SIC_RESOLVE != LSC_RESOLVE){
    printf("Wrong constants in sir_init\n"); exit(0);
  }    

  /* increase the counter for solvers */
  siv_nr_solvers++;

  /* set the current solver ID */
  siv_cur_solver_id = siv_nr_solvers-1;
  siv_solver[siv_cur_solver_id].solver_id = siv_cur_solver_id;
  siv_solver[siv_cur_solver_id].solver_type = Solver_type;
  siv_solver[siv_cur_solver_id].parallel = Parallel;

  // SIC_SINGLE_LEVEL == 1
  siv_solver[siv_cur_solver_id].nr_levels = Nr_levels;
  info = lsr_mkb_init(Solver_type, Parallel,  &siv_solver[siv_cur_solver_id].nr_levels, 
		      &siv_solver[siv_cur_solver_id].storage_type, 
		      Filename, Max_iter, Error_type, Error_tolerance,
		      Monitoring_level);


  /*kbw
  printf("\nIn sir_init: solver_id %d, Parallel %d, Nr_levels %d\n",
	 siv_cur_solver_id, siv_solver[siv_cur_solver_id].parallel, 
	 siv_solver[siv_cur_solver_id].nr_levels);
  /*kew*/

  // numbering in lsr_mkb MUST be the same as in sis_mkb!!!
  if(info!=siv_cur_solver_id){
    printf("Error 213758 in sir_init!!! Exiting.\n");
    exit(-1);
  }

  return(siv_cur_solver_id);
}


/*------------------------------------------------------------
  sir_mkb_show_blocks_matrix - save sm_matrix to file
------------------------------------------------------------*/
int sir_mkb_show_blocks_matrix(
                         /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id          /* in: level ID */


){

  int nr_levels, ineg, nrdofgl, max_nrdof;
  sit_levels *level_p; /* mesh levels */
  sit_dof_struct *dof_struct_p, *dof_struct_neig;

  int i,j,k, iaux, ient, iblock,dim,idof,neig;
  double *wiersz;

  char* sm_trace;
  unsigned char tmp = 1;
  unsigned char res = 0;
  unsigned char* sm_trace_bit;
  
  
/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* set the current solver ID */

  level_p = &(siv_solver[Solver_id].level[Level_id]);


  nrdofgl=level_p->nrdofs_glob;
  dim=nrdofgl;
  if(nrdofgl%8==0)dim=nrdofgl;
  else dim=nrdofgl+(8-nrdofgl%8);
  sm_trace= (char *)malloc((dim*nrdofgl)*sizeof(char));
  sm_trace_bit= (unsigned char *)malloc((dim*nrdofgl/8)*sizeof(char));
  for(i=0;i<(dim*nrdofgl);++i)sm_trace[i]='0';
  
  for(idof = 0; idof < level_p->nr_dof_ent; idof++){
    dof_struct_p = &level_p->l_dof_struct[idof];
    iaux=dim*dof_struct_p->posglob+dof_struct_p->posglob;
    //printf("\nW %d block=%d dim=%d posglob%d W\n",iaux,dof_struct_p->block_id,dim,dof_struct_p->posglob);
    for(i=0;i<dof_struct_p->nrdofs;i++,iaux++){
      for(j=0;j<dof_struct_p->nrdofs;j++){
        sm_trace[iaux+j*dim]='1';
      }
    }
    for(neig=0;neig<dof_struct_p->nrneig;neig++){
      for(i=0;level_p->l_dof_struct[i].block_id!=dof_struct_p->l_neig_bl[neig];i++);
      dof_struct_neig= &level_p->l_dof_struct[i];
      iaux=dim*dof_struct_p->posglob+dof_struct_neig->posglob;
      //printf("\nS %d block=%d dim=%d posglob%d %d S\n",iaux,dof_struct_neig->block_id,dim,dof_struct_neig->posglob,dof_struct_p->l_neig_bl[neig]);
      for(i=0;i<dof_struct_neig->nrdofs;i++,iaux++){
        for(j=0;j<dof_struct_p->nrdofs;j++){
          sm_trace[iaux+j*dim]='1';
        }
      }
    }
  }
  
  /*display sm_matrix 
  for(i=0;i<nrdofgl;i++){
    for(j=0;j<dim;j++)printf("%c",sm_trace[i*dim+j]);
    printf("\n");
  }
  /**/
  
  /* rewrite to binary format */
  for(j=0;j<(dim*nrdofgl/8);j++){
    res=0;
    for(i=0; i<8;i++) {
      //printf("%c",liczba[i]);
      if(sm_trace[8*j+i] == '1') {
        tmp <<= (7-i);
        res += tmp;
      }
    tmp = 1;
    }sm_trace_bit[j]=res;
  }



  //utr_io_write_img_to_pam(".","tescik",dim,nrdofgl,1,1,"BLACKANDWHITE",sm_trace,NULL);
  utr_io_write_img_to_pbm(".","sm_matrix","#sm_matrix",dim,nrdofgl,1,4,sm_trace_bit,NULL);

  free(sm_trace_bit);
  free(sm_trace);

return(0);

}



/*------------------------------------------------------------
  sir_create - to create solver's data structure
------------------------------------------------------------*/
int sir_create( /* returns: >0 - success code, <0 - error code */
  int Solver_id,    /* in: solver identification */
  int Problem_id    /* ID of the problem associated with the solver */
  )
{


/* local variables */
  sit_levels *level_p; /* mesh levels */

  /* pointer to dofs structure */
  sit_dof_struct *dof_struct_p;
  
  /* the number of (different) mesh entities for which entries to the global */
  /* stiffness matrix and load vector will be provided */
  int nr_int_ent, nr_dof_ent, max_dofs_per_dof_ent, max_dofs_int_ent, max_dof_ent_id;
  /* the global number of degrees of freedom */
  int nrdofs_glob, int_ent_id, nrdofs_ent_loc, dof_ent_id, pos_glob;
  int nr_dof_struct, dof_struct_id, nr_levels, nr_dof_ent_loc;
  int* temp_list_dof_type;
  int* temp_list_dof_id;
  int* temp_list_dof_nrdofs;
  int l_dof_ent_types[SIC_MAX_DOF_PER_INT], l_dof_ent_ids[SIC_MAX_DOF_PER_INT];
  int l_dof_ent_nrdofs[SIC_MAX_DOF_PER_INT];
  int *l_bl_nrdofs, *l_bl_posglob, *l_bl_nrneig, **l_bl_l_neig;
  int *l_bl_dof_ent_type, *l_bl_dof_ent_id;

  /* auxiliary variables */
  int i,j,k, idof, index, iaux, ient, ineig, level_id, pdeg_coarse, approx_type;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  siv_cur_solver_id=Solver_id;
  siv_solver[Solver_id].problem_id = Problem_id;

  nr_levels = siv_solver[Solver_id].nr_levels;

  /*ok_kbw*/
  printf("\nEntering sir_create (mkb): solver_id %d, Problem_id %d, nr_levels %d\n",
	 Solver_id, Problem_id, nr_levels);
  /*kew*/

  /* in a loop over levels */
  for(level_id=nr_levels-1;level_id>=0;level_id--){
    
    siv_solver[Solver_id].cur_level = level_id;
    level_p = &(siv_solver[Solver_id].level[level_id]);
    level_id = level_id;
    
    /*kbw
    printf("\nIn sir_create (mkb): solver_id %d, Problem_id %d, level %d\n",
	   Solver_id, Problem_id, level_id);
    /*kew*/
    
    /* get lists of integration and dof entities (types and IDs and nrdofs) */
    // the order on the lists for DOF entities should take into account 
    // parallel (message passing) execution (internal entities first!)
    /* the order on the lists of DOF entities determine the order of DOF structures !!! */
    // (but may not be the same as order passed to the solver due to renumbering)
    if(level_id==nr_levels-1){
      
      pdr_get_list_ent(Problem_id, &nr_int_ent, 
		       &level_p->l_int_ent_type, &level_p->l_int_ent_id, 
		       /* GHOST DOF ENTITIES HAVE NEGATIVE TYPE !!! */
		       &nr_dof_ent, &temp_list_dof_type, // &nr_dof_ent_internal,
		       &temp_list_dof_id, &temp_list_dof_nrdofs, 
		       &nrdofs_glob, &max_dofs_per_dof_ent);
      
      // pdeg_coarse specified for the finest level is used as input for coarse meshes
      pdeg_coarse = SIC_PDEG_FINEST;
      level_p->pdeg_coarse = pdeg_coarse;
      
/*kbw
    printf("In sir_create after pdr_get_list_ent\n");
    printf("nr_int_ent %d, pdeg_coarse %d\n", nr_int_ent, level_p->pdeg_coarse);
    for(i=0;i<nr_int_ent;i++)  printf("type %d, id %d\n",
	    level_p->l_int_ent_type[i],level_p->l_int_ent_id[i]);
    printf("\nNr_dof_ent %d, Nrdofs_glob %d, max_dofs_per_dof_ent %d\n",
	   nr_dof_ent, nrdofs_glob, max_dofs_per_dof_ent); 
    for(i=0;i<nr_dof_ent;i++)  printf("type %d, id %d, nrdof %d\n",
	    temp_list_dof_type[i], temp_list_dof_id[i],
				      temp_list_dof_nrdofs[i]);
    printf("level_p %u, level_p->l_dof_struct %u, level_p->l_int_ent_type %u, id %u\n",
	   level_p, level_p->l_dof_struct, level_p->l_int_ent_type, level_p->l_int_ent_id);
/*kew*/
      
    }
    else{
      
      int nr_int_ent_fine = siv_solver[Solver_id].level[level_id+1].nr_int_ent;
      int nr_dof_ent_fine = siv_solver[Solver_id].level[level_id+1].nr_dof_ent;
      int nrdof_glob_fine = siv_solver[Solver_id].level[level_id+1].nrdofs_glob;
      int max_dof_per_ent_fine = max_dofs_per_dof_ent;
      
      int * temp_list_dof_type_fine = 
	(int *) malloc( (nr_dof_ent_fine+1)*sizeof(int) );
      int *temp_list_dof_id_fine = 
	(int *) malloc( (nr_dof_ent_fine+1)*sizeof(int) );
      int *temp_list_dof_nrdof_fine = 
	(int *) malloc( (nr_dof_ent_fine+1)*sizeof(int) );
      
      for(i=0;i<nr_dof_ent_fine;i++){
	temp_list_dof_type_fine[i] = temp_list_dof_type[i];
	temp_list_dof_id_fine[i] = temp_list_dof_id[i];
	temp_list_dof_nrdof_fine[i] = temp_list_dof_nrdofs[i];
	
      }
      
      free(temp_list_dof_type);
      free(temp_list_dof_id);
      free(temp_list_dof_nrdofs);
      
      //works only for uniform approximation !!!!
      //pdeg_coarse = lsr_get_pdeg_coarse(Solver_id, level_id); 
      pdr_get_list_ent_coarse(Problem_id, nr_int_ent_fine, 
			      siv_solver[Solver_id].level[level_id+1].l_int_ent_type,
			      siv_solver[Solver_id].level[level_id+1].l_int_ent_id,
			      nr_dof_ent_fine, temp_list_dof_type_fine, 
			      temp_list_dof_id_fine, temp_list_dof_nrdof_fine, 
			      nrdof_glob_fine, max_dof_per_ent_fine,
			      &pdeg_coarse, &nr_int_ent, 
			      &level_p->l_int_ent_type, &level_p->l_int_ent_id, 
			      /* GHOST DOF ENTITIES HAVE NEGATIVE TYPE !!! */
			      &nr_dof_ent, &temp_list_dof_type,  // &nr_dof_ent_internal,
			      &temp_list_dof_id, &temp_list_dof_nrdofs, 
			      &nrdofs_glob, &max_dofs_per_dof_ent);
      

      // in the old setting the sons of all coarse elements were included in the fine level...
      // in the new setting projections between levels are always followed by exchange of DOFs

#ifdef DEBUG_SIM
      if(pdeg_coarse<=0){
	printf("Error in pdeg_coarse obtained from pdr_get_list_ent_coarse, %d\n"
	       , pdeg_coarse);
	exit(-1);
      }
#endif
      
      level_p->pdeg_coarse = pdeg_coarse;
      
      free(temp_list_dof_type_fine);
      free(temp_list_dof_id_fine);
      free(temp_list_dof_nrdof_fine);
      
    }
    
    
    level_p->nr_int_ent = nr_int_ent;
    level_p->nr_dof_ent = nr_dof_ent;
    level_p->nrdofs_glob = nrdofs_glob;
    level_p->block_size = -1; // initially indicate non-constant block_size
    

    /* array of structures storing DOF data */
    level_p->l_dof_struct = 
      (sit_dof_struct *) malloc( nr_dof_ent*sizeof(sit_dof_struct) );
    
/*kbw
    printf("In sir_create after pdr_get_list_ent_coarse\n");
    printf("nr_int_ent %d, pdeg_coarse %d\n", nr_int_ent, level_p->pdeg_coarse);
    for(i=0;i<nr_int_ent;i++)  printf("type %d, id %d\n",
	    level_p->l_int_ent_type[i],level_p->l_int_ent_id[i]);
    printf("\nNr_dof_ent %d, Nrdofs_glob %d, max_dofs_per_dof_ent %d\n",
	   nr_dof_ent, nrdofs_glob, max_dofs_per_dof_ent); 
    for(i=0;i<nr_dof_ent;i++)  printf("type %d, id %d, nrdof %d\n",
	    temp_list_dof_type[i], temp_list_dof_id[i],
				      temp_list_dof_nrdofs[i]);
    printf("level_p %u, level_p->l_dof_struct %u, level_p->l_int_ent_type %u, id %u\n",
	   level_p, level_p->l_dof_struct, level_p->l_int_ent_type, level_p->l_int_ent_id);
/*kew*/

    
    /* dof_ent_index to dof_struct_index (based on which dof_ent_id and */
    /* dof_ent_type can be find */
    level_p->max_dof_vert_id = -1;
    level_p->max_dof_edge_id = -1;
    level_p->max_dof_face_id = -1;
    level_p->max_dof_elem_id = -1;
    for(i=0; i< nr_dof_ent; i++){
      
      if(abs(temp_list_dof_type[i]) == PDC_ELEMENT ){
	if ( temp_list_dof_id[i] > level_p->max_dof_elem_id )
	  level_p->max_dof_elem_id = temp_list_dof_id[i];
      } else if(abs(temp_list_dof_type[i]) == PDC_FACE ){
	if ( temp_list_dof_id[i] > level_p->max_dof_face_id )
	  level_p->max_dof_face_id = temp_list_dof_id[i];
      } else if(abs(temp_list_dof_type[i]) == PDC_EDGE ){
	if ( temp_list_dof_id[i] > level_p->max_dof_edge_id )
	  level_p->max_dof_edge_id = temp_list_dof_id[i];
      } else if(abs(temp_list_dof_type[i]) == PDC_VERTEX ){
	if ( temp_list_dof_id[i] > level_p->max_dof_vert_id )
	  level_p->max_dof_vert_id = temp_list_dof_id[i];
      } else {
	printf("DOF_ENT: level %d, index %d, type %d, id %d, nrdof %d\n",
	       level_id, i, temp_list_dof_type[i], temp_list_dof_id[i],
				      temp_list_dof_nrdofs[i]);
	printf("Error 87373732 in sir_create!!! Exiting\n");
	exit(-1);
      }
      
      
    }
    
    if(level_p->max_dof_elem_id>=0){
      level_p->l_dof_elem_to_struct=
	(int*)malloc((level_p->max_dof_elem_id+1)*sizeof(int)); 
      for(i=0;i<=level_p->max_dof_elem_id;i++) 
	level_p->l_dof_elem_to_struct[i] = -1;
    }
    if(level_p->max_dof_face_id>=0){
      level_p->l_dof_face_to_struct=
	(int*)malloc((level_p->max_dof_face_id+1)*sizeof(int)); 
      for(i=0;i<=level_p->max_dof_face_id;i++) 
	level_p->l_dof_face_to_struct[i] = -1;
    }
    if(level_p->max_dof_edge_id>=0){
      level_p->l_dof_edge_to_struct=
	(int*)malloc((level_p->max_dof_edge_id+1)*sizeof(int)); 
      for(i=0;i<=level_p->max_dof_edge_id;i++) 
	level_p->l_dof_edge_to_struct[i] = -1;
    }
    if(level_p->max_dof_vert_id>=0){
      level_p->l_dof_vert_to_struct=
	(int*)malloc((level_p->max_dof_vert_id+1)*sizeof(int)); 
      for(i=0;i<=level_p->max_dof_vert_id;i++) 
	level_p->l_dof_vert_to_struct[i] = -1;
    }
    
    
    // PLACE FOR POSSIBLE CALL TO RENUMBERING PROCEDURE
    //!!! If renumbering not used, the order in pdr_ and sir_ routines the same and 
    //!!! the order in dof_struct table the same as returned by pdr_get_list_ent!!! 

    int nr_dofs_glob = 0;
    int nr_dofs_internal = 0;
    int nr_dof_ent_internal = 0;
    int const_block_size = temp_list_dof_nrdofs[0];
    for(index = 0; index < nr_dof_ent; index++){
      
      /* in the case of no renumbering index on temp_list_dof_... arrays */
      /* is the same as index=ID of DOF structure in the l_dof_struct array */
      idof = index;
      /* otherwise idof should be given by some output array */
      /* of the renumbering procedure idof = renumber_array[index]; */
      
      dof_struct_p = &level_p->l_dof_struct[idof];   
      
      nrdofs_ent_loc = temp_list_dof_nrdofs[index];
      nr_dofs_glob += nrdofs_ent_loc;
      if(nrdofs_ent_loc > level_p->max_dofs_dof_ent) {
	level_p->max_dofs_dof_ent = nrdofs_ent_loc;
      }
      if(nrdofs_ent_loc != const_block_size) const_block_size = -1;

      // alien (ghost) DOF entities (e.g. from overlap in parallel version) have negative type
      if(temp_list_dof_type[index]<0){
	// in mkb solver DOFs of ghost entities are used for vector product
	// but stiffness matrix rows are not stored for them
	// because of that they should be placed as the last DOFs!!!!!!!!
	dof_struct_p->nrneig = -1;
	temp_list_dof_type[index] = abs(temp_list_dof_type[index]);
      } else {
	dof_struct_p->nrneig = 0;
	nr_dofs_internal += nrdofs_ent_loc;
	nr_dof_ent_internal++;
#ifdef DEBUG_SIM
	if(nr_dof_ent_internal != index+1){
	  printf("Internal DOFs not returned as first from pdr_get_list_ent!\n");
	  printf("PARALLEL version will not work in this case - unless suitable\n");
	  printf("renumbering is introduced in sir_create!!!!!!!!!! Exiting.\n");
	  exit(-1);
	}
#endif
      }
      dof_struct_p->dof_ent_type = temp_list_dof_type[index];
      dof_struct_p->dof_ent_id = temp_list_dof_id[index];
      dof_struct_p->nrdofs = nrdofs_ent_loc;
      
      /* initialize lists of integration entities and neighbouring dof_ent */
      dof_struct_p->nr_int_ent = 0;
      for(i=0;i<SIC_MAX_INT_PER_DOF;i++) 
	dof_struct_p->l_int_ent_index[i]=SIC_LIST_END_MARK;
      for(i=0;i<SIC_MAX_DOF_STR_NGB;i++) 
	dof_struct_p->l_neig[i]=SIC_LIST_END_MARK;
      
      /* arrays dof_ent to dof_struct */
      if(dof_struct_p->dof_ent_type == PDC_ELEMENT ) {
	level_p->l_dof_elem_to_struct[dof_struct_p->dof_ent_id] = idof;
	
#ifdef DEBUG_SIM
	if(dof_struct_p->dof_ent_id > level_p->max_dof_elem_id){
	  printf("Error 84543732 in sir_create!!! Exiting\n");
	  printf("%d > %d\n", dof_struct_p->dof_ent_id, level_p->max_dof_elem_id);
	  exit(-1);
	}
#endif
	
      } else if(dof_struct_p->dof_ent_type == PDC_FACE ) {
	level_p->l_dof_face_to_struct[dof_struct_p->dof_ent_id] = idof;
	
#ifdef DEBUG_SIM
	if(dof_struct_p->dof_ent_id > level_p->max_dof_face_id){
	  printf("Error 84543733 in sir_create!!! Exiting\n");
	  printf("%d > %d\n", dof_struct_p->dof_ent_id, level_p->max_dof_face_id);
	  exit(-1);
	}
#endif
	
      } else if(dof_struct_p->dof_ent_type == PDC_EDGE ) {
	level_p->l_dof_edge_to_struct[dof_struct_p->dof_ent_id] = idof;
	
#ifdef DEBUG_SIM
	if(dof_struct_p->dof_ent_id > level_p->max_dof_edge_id){
	  printf("Error 84543734 in sir_create!!! Exiting\n");
	  printf("%d > %d\n", dof_struct_p->dof_ent_id, level_p->max_dof_edge_id);
	  exit(-1);
	}
#endif
	
      } else if(dof_struct_p->dof_ent_type == PDC_VERTEX ) {
	level_p->l_dof_vert_to_struct[dof_struct_p->dof_ent_id] = idof;
#ifdef DEBUG_SIM
	if(dof_struct_p->dof_ent_id > level_p->max_dof_vert_id){
	  printf("Error 84543735 in sir_create!!! Exiting\n");
	  printf("%d > %d\n", dof_struct_p->dof_ent_id, level_p->max_dof_vert_id);
	  exit(-1);
	}
#endif
      } else {
	printf("Error 8963732 in sir_create!!! Exiting\n");
	exit(-1);
      }
      
      
/*kbw
    printf("In sir_create after filling dof_struct no %d\n", idof);
    printf("dof_ent_type %d, dof_ent_id %d, nrdof %d\n",
	   dof_struct_p->dof_ent_type , dof_struct_p->dof_ent_id, 
	   dof_struct_p->nrdofs);
    printf("Initialized lists of int_ent %d, neig %d\n",
	   dof_struct_p->l_int_ent_index[0], dof_struct_p->l_neig[0]); 
    getchar();
/*kew*/
        
    }
    
    if(nr_dofs_glob != level_p->nrdofs_glob){
      printf("Error in counting nr_dofs_glob in sir_create (%d != %d). Exiting!\n",
	     nr_dofs_glob, level_p->nrdofs_glob);
    }
    level_p->nrdofs_internal = nr_dofs_internal;
    level_p->nr_dof_ent_internal = nr_dof_ent_internal;

    // initialized number of DOF entities summed over all integration entities
    level_p->nr_asse_blocks_all_int_ent = 0;

    // test whether block size is constant for all DOFs
    if(const_block_size > 0) level_p->block_size =  const_block_size;   

    pdeg_coarse = SIC_PDEG_FINEST;
    if(level_id<nr_levels-1) pdeg_coarse = level_p->pdeg_coarse;

    /* getting information on the structure of the global stiffness matrix */
    nr_dof_struct = 0;
    max_dofs_int_ent = 0;
    for(ient=0; ient<level_p->nr_int_ent;ient++){
      
      int idofent;     
      int nrdofs_int_ent = 0;
      int_ent_id = level_p->l_int_ent_id[ient];
      nr_dof_ent_loc = SIC_MAX_DOF_PER_INT;
      
      pdr_comp_stiff_mat(Problem_id, level_p->l_int_ent_type[ient], 
			 level_p->l_int_ent_id[ient], PDC_NO_COMP, &pdeg_coarse,
			 &nr_dof_ent_loc, l_dof_ent_types, 
			 l_dof_ent_ids, l_dof_ent_nrdofs,
			 NULL, NULL, NULL, NULL);
      
      
      level_p->nr_asse_blocks_all_int_ent += (nr_dof_ent_loc+1)*nr_dof_ent_loc;

/*kbw
   printf("in sir_create after pdr_comp_stiff_mat for int_ent no %d (id %d):\n", 
	   ient, int_ent_id);
    printf("nr_dof_ent_loc %d (nr_asse_blocks_all_int_ent %d), types, ids, nrdofs:\n", 
    nr_dof_ent_loc, level_p->nr_asse_blocks_all_int_ent);
    for(idofent=0; idofent<nr_dof_ent_loc; idofent++){
      printf("%d %d %d\n", l_dof_ent_types[idofent],
    	     l_dof_ent_ids[idofent], l_dof_ent_nrdofs[idofent]);
    }
/*kew*/

      for(idofent=0; idofent<nr_dof_ent_loc; idofent++){
	
	int dof_ent_id = l_dof_ent_ids[idofent];
	int dof_ent_type = l_dof_ent_types[idofent];
	
	if(dof_ent_type == PDC_ELEMENT ) {
	  
	  dof_struct_id = level_p->l_dof_elem_to_struct[dof_ent_id];
	  
#ifdef DEBUG_SIM
	  if( level_p->l_dof_elem_to_struct[dof_ent_id] == -1){
	    printf("Error 3472941 in sir_create!!! Exiting\n");
	    exit(-1);
	  }
#endif
	  
	} else if(dof_ent_type == PDC_FACE ) {
	  
	  dof_struct_id = level_p->l_dof_face_to_struct[dof_ent_id];
	  
#ifdef DEBUG_SIM
	  if( level_p->l_dof_face_to_struct[dof_ent_id] == -1){
	    printf("Error 3472942 in sir_create!!! Exiting\n");
	    exit(-1);
	  }
#endif
	  
	} else if(dof_ent_type == PDC_EDGE ) {
	  
	  dof_struct_id = level_p->l_dof_edge_to_struct[dof_ent_id];
	  
#ifdef DEBUG_SIM
	  if( level_p->l_dof_edge_to_struct[dof_ent_id] == -1){
	    printf("Error 3472943 in sir_create!!! Exiting\n");
	    exit(-1);
	  }
#endif
	  
	} else if(dof_ent_type == PDC_VERTEX ) {
	  
	  dof_struct_id = level_p->l_dof_vert_to_struct[dof_ent_id];
	  
#ifdef DEBUG_SIM
	  if( level_p->l_dof_vert_to_struct[dof_ent_id] == -1){
	    printf("Error 3472944 for %d in sir_create!!! Exiting\n",
		   dof_ent_id);
	    exit(-1);
	  }
#endif
	  
	} else {
	  printf("Error 34fsf7294 in sir_create!!! Exiting\n");
	  exit(-1);
	}
	
	dof_struct_p = &level_p->l_dof_struct[dof_struct_id];
	
/*kbw
	printf("for int_type %d, int_id %d, dof_type %d, dof_id %d, struct %d\n",
	       level_p->l_int_ent_type[ient], level_p->l_int_ent_id[ient],
	       l_dof_ent_types[idofent], l_dof_ent_ids[idofent],dof_struct_id );
/*kew*/

#ifdef DEBUG_SIM
	if((dof_struct_p->dof_ent_id != dof_ent_id) || 
	   (dof_struct_p->dof_ent_type != dof_ent_type) || 
	   (dof_struct_p->nrdofs != l_dof_ent_nrdofs[idofent]) ){
	  printf("Error 3827 in sir_create!!! Exiting");
	  exit(-1);
	}
#endif      
	
	nrdofs_int_ent += l_dof_ent_nrdofs[idofent];
	
/*kbw
      printf("putting int_ent no %d on the list of int_ent, nr_int_ent %d\n", 
	     ient, dof_struct_p->nr_int_ent);
      printf("before:");
      for(i=0;i<SIC_MAX_INT_PER_DOF;i++) {
	printf("%d",dof_struct_p->l_int_ent_index[i]) ;
      }
      printf("\n");
/*kew*/

	iaux=sir_put_list(ient, 
			  dof_struct_p->l_int_ent_index, SIC_MAX_INT_PER_DOF);
	if(iaux<0) dof_struct_p->nr_int_ent++;
	
/*kbw
      if(dof_struct_p->nr_int_ent>77){
	printf("For dof_ent %d, list of int_ent:\n", dof_ent_id);
	for(i=0;i<dof_struct_p->nr_int_ent;i++){
	  printf("%d  ", dof_struct_p->l_int_ent_index[i]);
	}
	printf("\n");
      }
/*kew*/

#ifdef DEBUG_SIM
	if(iaux == 0){ // list full - increase SIC_MAX_INT_PER_DOF
	  printf("Error 383627 (nr_int_ent %d) in sir_create!!! Exiting", 
		 dof_struct_p->nr_int_ent);
	  exit(-1);
	}
#endif      
/*kbw
      printf("putting int_ent no %d on the list of int_ent, nr_int_ent %d\n", 
	     ient,dof_struct_p->nr_int_ent);
      printf("after:");
      for(i=0;i<SIC_MAX_INT_PER_DOF;i++) {
	printf("%d",dof_struct_p->l_int_ent_index[i]) ;
      }
      printf("\n");
/*kew*/

	for(ineig = 0; ineig<nr_dof_ent_loc; ineig++){ 
	  
	  //	if(ineig != idofent){
	  /* change for constrained approximation */
	  /* ineig may be != idofent but this is the same dof_ent */
	  if(dof_ent_type != l_dof_ent_types[ineig] ||
	     dof_ent_id != l_dof_ent_ids[ineig]){	
	    
	    int neig_id = l_dof_ent_ids[ineig];
	    int neig_type = l_dof_ent_types[ineig];

/*kbw
      printf("dof_ent %d: putting ineig no %d (id %d) on the list of neig, nrneig %d\n", 
	     idofent, ineig, neig_id, dof_struct_p->nrneig);
      printf("before:");
      for(i=0;i<SIC_MAX_DOF_STR_NGB;i++) {
	printf("%d",dof_struct_p->l_neig[i]) ;
      }
      printf("\n");
/*kew*/

	    int neig_index = 0;
	    if(neig_type==PDC_ELEMENT){
	      neig_index = level_p->l_dof_elem_to_struct[neig_id];
	    } else if(neig_type==PDC_FACE){
	      neig_index = level_p->l_dof_face_to_struct[neig_id];
	    } else if(neig_type==PDC_EDGE){
	      neig_index = level_p->l_dof_edge_to_struct[neig_id];
	    } else if(neig_type==PDC_VERTEX){
	      neig_index = level_p->l_dof_vert_to_struct[neig_id];
	    } 
	    
/*kbw
      printf("dof_ent %d: putting ineig no %d (id %d, index %d) on the list of neig, nrneig %d\n", 
	     idofent, ineig, neig_id, neig_index, dof_struct_p->nrneig);
      printf("before:");
      for(i=0;i<SIC_MAX_DOF_STR_NGB;i++) {
	printf("%d",dof_struct_p->l_neig[i]) ;
      }
      printf("\n");
/*kew*/

	    /* only active elements=DOF entities store neighbors */
	    if(dof_struct_p->nrneig>=0){
	      
	      iaux=sir_put_list(neig_index, 
				dof_struct_p->l_neig, SIC_MAX_DOF_STR_NGB);
	      if(iaux<0) {
		
		dof_struct_p->nrneig++;
		
	      }

#ifdef DEBUG_SIM
	      if(iaux == 0){ // list full - increase SIC_MAX_DOF_STR_NGB
		printf("Error 385627 in sir_create!!! Exiting");
		exit(-1);
	      }
#endif      

/*kbw
      printf("putting ineig no %d (id %d) on the list of neig, nrneig %d\n", 
	     ineig, neig_id, dof_struct_p->nrneig);
      printf("after:");
      for(i=0;i<SIC_MAX_DOF_STR_NGB;i++) {
	printf("%d",dof_struct_p->l_neig[i]) ;
      }
      printf("\n");
/*kew*/

	    }
	    
	  }
	  
	}
	
      } /* end loop over dof_ents of a given int_ent (idofent) */
      
      if(nrdofs_int_ent > max_dofs_int_ent) max_dofs_int_ent = nrdofs_int_ent; 
      
    } /* end loop over int_ent */

    // ANOTHER OPTION FOR RENUMBERING:
    // ordering in si module is the same as given by FEM part in get_list_ent
    // but each struct has different block_id and posglob passed to solver module

    // IT IS ASSUMED THAT IN PARALLEL VERSION INTERNAL DOFS ARE GROUPED BEFORE OVERLAP
    // DOFS - HENCE RENUMBERING SHOULD BE DONE IN SUCH A WAY, THAT FIRST
    // INTERNAL DOFS ARE PLACED FIRST (IF NOT RETURNED IN SUCH A WAY FROM FEM PART)
    // AND THEN ONLY INTERNAL DOFS ARE RENUMBERED

    // MOREOVER: MULTITHREADING DIVIDES AND RENUMBERS ONLY INTERNAL DOFS
    // HENCE: WE NEED ANOTHER PARAMETER:
    // nr_dof_ent_internal - the number of internal DOF entities

    
    // HERE INTERNAL DOFS SHOULD BE PLACED AT THE BEGINNING
    // ......................


    // if NO RENUMBERING:
    // pos_glob is an index of the first DOF associated with a given dof entity
    // (i.e. DOF block in solver nomenclature) in a vector of unknowns
    // (it also indicates the row in global SM associated with this DOF entity).
    // For systems with the constant number of DOFs per DOF entity, 
    // i.e. level_p->l_dof_struct[dof_ent_id].nrdofs, for which this constant
    // number is stored as level_p->block_size
    // (currently all, except DG with different p in different elements),
    // pos_glob can be computed as level_p->block_size * block_id !!!
    // (this information can be used by solvers internally, the situation
    // of constant block_size is assumed when  nrdofgl is a multiple of 
    // nr_dof_ent, i.e. level_p->nrdofs_glob % level_p->nr_dof_ent == 0 )
    // 
    pos_glob = 0; int idofent;
    for (idofent = 0; idofent< nr_dof_ent; idofent++){
      level_p->l_dof_struct[idofent].block_id = idofent;
      level_p->l_dof_struct[idofent].posglob = pos_glob;
      pos_glob += level_p->l_dof_struct[idofent].nrdofs;
      
#ifdef DEBUG_SIM
      if(level_p->block_size > 0){
    	if(level_p->l_dof_struct[idofent].posglob != 
	   level_p->l_dof_struct[idofent].block_id * level_p->block_size){
	  printf("posglob not in agreement with block_id in sir_create!!! Exiting\n");
	  exit(-1);
	}
      }    
#endif

    }
    
#ifdef DEBUG_SIM
    if ( level_p->nrdofs_glob != pos_glob ){
      printf("Error 843732135 in sir_create!!! Exiting\n");
      exit(-1);
    }
#endif
  
    // PREPARE INPUT FOR RENUMBERING FOR INTERNAL DOFS AND THEN
    //call renumbering procedures:

    // old internal procedure operating on sid_mkb data structure
#ifdef INTERNAL_RENUMBERING
    //printf("Internal renumbering !!!!\n"); exit(-1);
    sir_reverse_cuthill_mckee(level_p->l_dof_struct, level_p->nr_dof_ent);
#endif


    // generic procedure for graph vertex permutation (renumbering)
#ifdef RENUMBERING
    //printf("Renumbering !!!!\n"); exit(-1);
    int *permutation_array =(int *)malloc( nr_dof_ent*sizeof(int) );
    int *nrneig=(int *) malloc( nr_dof_ent*sizeof(int) );
    int **neig=(int **) malloc( nr_dof_ent*sizeof(int*) );
    for(idofent = 0; idofent< nr_dof_ent; idofent++){
      nrneig[idofent]=level_p->l_dof_struct[idofent].nrneig;
      neig[idofent]=level_p->l_dof_struct[idofent].l_neig;
    }
    utr_renumber(nr_dof_ent, nrneig, neig, permutation_array);
    free(nrneig);
    free(neig);

    // HERE POS_GLOB IS ASSIGNED AFTER RENUMBERING !!!
    // pos_glob is an index of the first DOF associated with a given dof entity
    // (i.e. DOF block in solver nomenclature) in a vector of unknowns
    // (it also indicates the row in global SM associated with this DOF entity).
    // For systems with the constant number of DOFs per DOF entity, 
    // i.e. level_p->l_dof_struct[dof_ent_id].nrdofs, for which this constant
    // number is stored as level_p->block_size
    // (currently all, except DG with different p in different elements),
    // pos_glob can be computed as level_p->block_size * block_id !!!
    // (this information can be used by solvers internally, the situation
    // of constant block_size is assumed when  nrdofgl is a multiple of 
    // nr_dof_ent, i.e. level_p->nrdofs_glob % level_p->nr_dof_ent == 0 )
    // 

    int start_dof=0; 
    for(idofent = 0; idofent< nr_dof_ent; idofent++){ 
      level_p->l_dof_struct[permutation_array[idofent]].block_id = idofent; 
      level_p->l_dof_struct[permutation_array[idofent]].posglob=start_dof; 
      start_dof+=level_p->l_dof_struct[permutation_array[idofent]].nrdofs; 

#ifdef DEBUG_SIM 
      if(level_p->block_size > 0){ 
	if(level_p->l_dof_struct[idofent].posglob !=  
	   level_p->l_dof_struct[idofent].block_id * level_p->block_size){ 
	  printf("posglob not in agreement with block_id in sir_create!!! Exiting\n"); 
	  exit(-1); 
	} 
      }     
#endif
      
    }
    free(permutation_array);	 

#endif  //RENUMBERING

#ifdef RENUMBERING_FOR_PARALLEL
    printf("RENUMBERING_FOR_PARALLEL not implemented yet!!!\n");
#endif
  
    int istruct;
    // use renumbering data to prepare data for assembling
    for(istruct = 0; istruct< nr_dof_ent; istruct++){
      
      //* set the list of neighboring blocks (NOT DOF ENTITIES!!!)*/
      for(ineig=0;ineig<level_p->l_dof_struct[istruct].nrneig;ineig++){
	iaux=level_p->l_dof_struct[istruct].l_neig[ineig];
	level_p->l_dof_struct[istruct].l_neig_bl[ineig]=level_p->l_dof_struct[iaux].block_id;
      }
      
/*kbw
    printf("dof_struct %d, dof_ent_type %d, dof_ent_id %d, block_id %d, posglob %d, nrneig %d\n",
	   istruct, level_p->l_dof_struct[istruct].dof_ent_type , level_p->l_dof_struct[istruct].dof_ent_id , 
	   level_p->l_dof_struct[istruct].block_id , level_p->l_dof_struct[istruct].posglob , 
	   level_p->l_dof_struct[istruct].nrneig);
    printf("neighbors before sorting: \n");
    for(ineig=0;ineig<level_p->l_dof_struct[istruct].nrneig;ineig++){
      printf("neig_id %d, struct_id %d, block_id %d\n",
	     ineig, level_p->l_dof_struct[istruct].l_neig[ineig],
	     level_p->l_dof_struct[istruct].l_neig_bl[ineig]);
    }
/*kew*/

      sir_sort_short(level_p->l_dof_struct[istruct].l_neig_bl, 0, 
		     level_p->l_dof_struct[istruct].nrneig-1);

/*kbw
    printf("neighbors after sorting: \n");
    for(ineig=0;ineig<level_p->l_dof_struct[istruct].nrneig;ineig++){
      printf("neig_id %d, struct_id %d, block_id %d\n",
	     ineig, level_p->l_dof_struct[istruct].l_neig[ineig],
	     level_p->l_dof_struct[istruct].l_neig_bl[ineig]);
    }
    getchar();
/*kew*/
      

    }

    
    level_p->max_dofs_int_ent = max_dofs_int_ent;
    
#ifdef PRINT_MATRIX
    sir_mkb_show_blocks_matrix(Solver_id,level_id);
#endif

  } // end loop over levels


  // temporary lists created by pdr_get_list_ent... (used by all levels) are free here
  free(temp_list_dof_type);
  free(temp_list_dof_id);
  free(temp_list_dof_nrdofs);

  // COLORING!!!
  //
  // DIVIDE INTEGRATION ENTITIES INTO GROUPS WITH THE SAME TYPE, PDEG AND THE SAME COLOUR
  // (by permuting L_int_ent_id list and creating CRS like l_color_index arrays)
  //


#ifdef MKB_SOLVE_TIME_TEST
  double time_coloring=time_clock();
#endif

  int active_coloring = 1; // set coloring on as default


  nr_levels = siv_solver[Solver_id].nr_levels;

  // simple test for type of approximation
  approx_type = SIC_STD_APPROX;
  if(siv_solver[Solver_id].level[nr_levels-1].max_dof_elem_id>=0) {
    approx_type = SIC_DG_APPROX;
  }
  if(approx_type==SIC_DG_APPROX) active_coloring = 0;

  // switch off colouring for small meshes (the same limit as for assembly tables!!!!!)
  if(siv_solver[Solver_id].level[nr_levels-1].nrdofs_glob<1000){

#ifdef OPENCL_GPU
    printf("Problem too small for GPU computing! Use CPU-only version!\n");
    exit(0);
#endif

    active_coloring = 0; 
  }
  
  /* in a loop over levels */
  for(level_id=nr_levels-1;level_id>=0;level_id--){
  
   siv_solver[Solver_id].cur_level = level_id;
   level_p = &(siv_solver[Solver_id].level[level_id]);
    
    // divide list into elements and faces
    int nr_elems=0;  
    int int_entity; 
    int previous_type = PDC_ELEMENT;
    int ok = 1;
    for(int_entity=0;int_entity<level_p->nr_int_ent;int_entity++){
      if(level_p->l_int_ent_type[int_entity]==PDC_ELEMENT){
	nr_elems++;
	if(previous_type != PDC_ELEMENT) ok = 0;
      }
      else previous_type = PDC_FACE;
    }
    int nr_faces =  level_p->nr_int_ent - nr_elems;
    
    if(ok!=1){
      printf("Elements are not first on the list to integrate. \n");
      printf("OpenMP is not optimized for that case.\n");
      exit(-1);
    }

    if(active_coloring==0){ // if coloring deactivated
      
      level_p->nr_colors_elems = 1;
      level_p->l_color_index_elems=(int*)malloc((level_p->nr_colors_elems+1)*sizeof(int));
      level_p->l_color_index_elems[0] = 0;
      level_p->l_color_index_elems[1] = nr_elems;
      
      level_p->nr_colors_faces = 1;
      level_p->l_color_index_faces =(int*)malloc((level_p->nr_colors_faces+1)*sizeof(int));
      level_p->l_color_index_faces[0] = nr_elems;
      level_p->l_color_index_faces[1] = nr_elems+nr_faces;
      
    }
    else{ // if coloring activated (default)
      
      level_p->l_color_index_elems = NULL;
      level_p->l_color_index_faces = NULL;  
        
    int** l_int_ent_index = (int **)malloc((level_p->nr_dof_ent)*sizeof(int *));
    int* nr_int_ent = (int *)malloc((level_p->nr_dof_ent)*sizeof(int));
    int i_dof_ent;
    for(i_dof_ent=0;i_dof_ent<level_p->nr_dof_ent;i_dof_ent++){

/*kbw
      printf("dof_ent %d, number of int_ents %d, list: \n",
	     i_dof_ent, level_p->l_dof_struct[i_dof_ent].nr_int_ent);
      for(i=0; i < level_p->l_dof_struct[i_dof_ent].nr_int_ent; i++){
	printf("%10d  ", level_p->l_dof_struct[i_dof_ent].l_int_ent_index[i]);
      }
      printf("\n");
/*kew*/

      l_int_ent_index[i_dof_ent]=&level_p->l_dof_struct[i_dof_ent].l_int_ent_index[0];
      nr_int_ent[i_dof_ent]=level_p->l_dof_struct[i_dof_ent].nr_int_ent;
    }

      utr_color_int_ent_for_assembly(Problem_id, level_id, nr_elems, nr_faces,
                                   &level_p->l_int_ent_type[0],
                                   &level_p->l_int_ent_id[0],
                                   level_p->nr_dof_ent,
                                   &nr_int_ent[0],
                                   &l_int_ent_index[0],
                                   &level_p->nr_colors_elems,
                                   &level_p->l_color_index_elems,
                                   &level_p->nr_colors_faces,
                                   &level_p->l_color_index_faces);
                                   
    free(l_int_ent_index);
    free(nr_int_ent);
    
    
    /* //elements colouring */
    /* utr_color_int_ent_for_assembly_with_graph_creation(Problem_id, */
    /*                                level_id, nr_elems, 0, */
    /*                                &level_p->l_int_ent_type[0],  */
    /*                                &level_p->l_int_ent_id[0],  */
    /*                                &level_p->nr_colors_elems,  */
    /*                                &level_p->l_color_index_elems); */
    
    /* //faces colouring */
    /* utr_color_int_ent_for_assembly_with_graph_creation(Problem_id, */ 
    /*                                level_id, nr_faces, nr_elems, */
    /*                                &level_p->l_int_ent_type[nr_elems],  */
    /*                                &level_p->l_int_ent_id[nr_elems],  */
    /*                                &level_p->nr_colors_faces,  */
    /*                                &level_p->l_color_index_faces); */
    
    /* //or all integration entities together(elems+faces) colouring */
    /* // utr_color_int_ent_for_assembly_with_graph_creation(Problem_id, */ 
    /* // Level_id, Nr_int_ent, 0, */
    /* // &L_int_ent_type[0], &L_int_ent_id[0],  */
    /* //&nr_colors_elems, &l_color_index_elems); */
    
    
/*ok_kbw*/
#ifdef MKB_SOLVE_TIME_TEST
    printf("Time for coloring: %lf\n", time_clock() - time_coloring);
#endif
/*kew*/

   
   
    } // end if coloring active (not DG yet)
  
  } // end loop over levels

  // At that moment we assume that there is a set of fixed data used for assembly
  // and system solution (including parallel solution with message passing)
  //
  // 1. The list of integration entities for which SMs and LVs are passed to the solver
  //    (the order is fixed as well, so possible colouring must be done before!!!)
  //    the lists: l_int_ent_type, l_int_ent_id for all nr_int_ent entities.
  //    As an option additional lists: elem_pos_for_colours and face_pos_for_colours
  //    storing in elem_pos_for_colours[icolour] the position in int_ent tables
  //    of the first element of a given colour (and in the same way for faces)
  //    (the scheme is the same as for rows in CRS)
  //
  // 2. The data passed to the solver for DOF blocks (at the current implementation
  //    identified also with DOF entities (the difference is that DOF entities have
  //    IDs and types associated with the FEM mesh). DOF blocks are characterized by:
  //    - id: simply the position on the lists passed to the solver
  //    - posglob: the position in global LV of the first DOF, i.e. the position 
  //      of the associated row in SM (inducing also the position in SM of the other
  //      associated entries in other rows and the associated column)
  //      (for CRS storage posglob==id, for BCRS posglob==id*block_size)
  //    - nrdofs: the number of individual scalar DOFs
  //      (for CRS storage nrdofs==1, for BCRS nrdofs==block_size)
  //    - l_neig_bl: the list of neighbouring blocks, i.e. the blocks that have
  //      non-zero entries in the row of SM associated with the block
  //      IMPORTANT assumption: the list is sorted
  //      (there exist also l_neig lists where neighbouring DOF structures are listed
  //      i.e. instead of block IDs, DOF structures indexes are used)
  //
  // 3. Additional data structure l_int_ent_index stores for each DOF structure the indices
  //    (in tables l_int_ent_type, l_int_ent_id) of the associated integration entities



  nr_levels = siv_solver[Solver_id].nr_levels;

  int storage_type;
  if(approx_type==SIC_DG_APPROX && siv_solver[Solver_id].storage_type != LSC_STORAGE_PETSC) {
    // because of special preconditioners we do not want (may be temporarily) BCRS to work for DG
    storage_type = LSC_STORAGE_BLOCK;
  }
  else{
    storage_type = siv_solver[Solver_id].storage_type; // returned from lsd_...
  }

  /* in a loop over levels */
  for(level_id=nr_levels-1;level_id>=0;level_id--){
    siv_solver[Solver_id].cur_level = level_id;
    level_p = &(siv_solver[Solver_id].level[level_id]);
    
/*kbw
    printf("filling temp lists for mkb_create_matrix: solver %d, level %d (%lu)\n",
	   Solver_id, level_id, level_p);
    printf("nrblocks %d, max_sm_size %d, nrdofs_glob %d\n",
	   level_p->nr_dof_ent, max_dofs_int_ent, level_p->nrdofs_glob);
/*kew*/

    // !!!!!!!!!!!!! important test: !!!!!!!!!!!!!!!
    // if nrdofgl is a multiple of nr_dof_ent - it is assumed that
    // block size is constant !!!!!!!!!!!!!!!!!!!!!!
    // lsr_mkb_create_matrix uses this information internally !!!!!!!
    if(level_p->nrdofs_glob % level_p->nr_dof_ent == 0 ){

      int induced_block_size = level_p->nrdofs_glob / level_p->nr_dof_ent;
/*kbw
      printf("test induced block_size %d, stored in data structure %d\n",
	     induced_block_size, level_p->block_size);
/*kew*/

      if(induced_block_size != level_p->block_size){
	printf("Inconsistency in passing block_size to lsd_mkb_create_matrix:\n");
	printf("test induced block_size %d != stored in data structure %d !!!\n",
	       induced_block_size, level_p->block_size);
	printf("Exiting.\n");
	exit(-1);
      }
    } // end block_size test

    /* allocate memory for temporary lists */
    if(siv_solver[Solver_id].parallel){
      l_bl_dof_ent_id =  (int *)malloc((level_p->nr_dof_ent+1)*sizeof(int));
      l_bl_dof_ent_type =  (int *)malloc((level_p->nr_dof_ent+1)*sizeof(int));
    }
    l_bl_nrdofs =  (int *)malloc((level_p->nr_dof_ent+1)*sizeof(int));
    l_bl_posglob = (int *)malloc((level_p->nr_dof_ent+1)*sizeof(int));
    l_bl_nrneig = (int *)malloc((level_p->nr_dof_ent+1)*sizeof(int));
    l_bl_l_neig = (int **)malloc((level_p->nr_dof_ent+1)*sizeof(int *));
    
    int istruct;
    for(istruct = 0; istruct< level_p->nr_dof_ent; istruct++){
      
      dof_struct_id = istruct;
      
      dof_struct_p = &level_p->l_dof_struct[dof_struct_id];
      
      /* !!! offset 1 numbering of blocks moved to lad_block !!! */
      //int ibl = dof_struct_p->block_id +1;
      int ibl = dof_struct_p->block_id;

/*kbw
      printf("dof_struct_p %lu, id %d, dof_ent %d\n", dof_struct_p,
	     dof_struct_id, dof_struct_p->dof_ent_id);
      printf("nrdofs %d, posglob %d, nrneig %d\nneig: ",
	     dof_struct_p->nrdofs, dof_struct_p->posglob, dof_struct_p->nrneig);
      for(ineig=0;ineig<dof_struct_p->nrneig;ineig++){
	int iaux=dof_struct_p->l_neig[ineig];
	printf("%6d",level_p->l_dof_struct[iaux].dof_ent_id);
      }
      printf("\nint_ent: ");
	for(i=0;i<dof_struct_p->nr_int_ent;i++){
	  printf("%d  ", dof_struct_p->l_int_ent_index[i]);
	}
	printf("\n");
	getchar();
/*kew*/

      l_bl_nrdofs[ibl] = dof_struct_p->nrdofs;
      l_bl_posglob[ibl] = dof_struct_p->posglob;
      l_bl_nrneig[ibl] = dof_struct_p->nrneig;

      if(dof_struct_p->nrneig>0){
	l_bl_l_neig[ibl]= (int *)malloc(l_bl_nrneig[ibl]*sizeof(int));
	for(ineig=0;ineig<l_bl_nrneig[ibl];ineig++){
	  // int iaux=dof_struct_p->l_neig[ineig]; - old without renumbering
	  int iaux=dof_struct_p->l_neig_bl[ineig];
	  /* !!! offset 1 numbering of blocks moved to lad_block !!! */
	  //l_bl_l_neig[ibl][ineig]=iaux+1;
	  l_bl_l_neig[ibl][ineig]=iaux;
	}

      }
    }
  
  
/*kbw
  printf("\nbefore calling mkb_create_matrix: solver %d, level %d, storage_type %d (approx %d)\n",
	 Solver_id, level_id, storage_type, approx_type );
  printf("nrblocks %d, max_sm_size %d, nrdofs_glob %d, block_size %d\n",
	 level_p->nr_dof_ent, max_dofs_int_ent, level_p->nrdofs_glob, level_p->block_size);
  printf("block_id,\tnr_dof_ent,\tposglob,\tnroffbl\tneighbors:\n");
// !!! offset 1 numbering of blocks moved to lad_block !!! 
  int ibl;
//  for(ibl = 1; ibl<= level_p->nr_dof_ent; ibl++){
  for(ibl = 0; ibl< level_p->nr_dof_ent; ibl++){
    printf("%d\t\t%d\t\t%d\t\t%d\t",ibl, l_bl_nrdofs[ibl], l_bl_posglob[ibl], 
	   l_bl_nrneig[ibl]);
    for(j=0;j<l_bl_nrneig[ibl]; j++){
      printf("%10d",l_bl_l_neig[ibl][j]);
    }
    printf("\n");
  }
/*kew*/

  
    lsr_mkb_create_matrix(Solver_id, level_id, 
			  storage_type, level_p->nr_dof_ent, 
			  level_p->nrdofs_glob, level_p->block_size,
			  max_dofs_int_ent, l_bl_nrdofs, l_bl_posglob, 
			  l_bl_nrneig, l_bl_l_neig);
    
 
    /* create tables for exchanging data between processors */ 
    if(siv_solver[Solver_id].parallel){

      // TABLES FOR BLOCKS ARE REUSED BUT THE CONTENT IS FOR DOF STRUCTURES
      
      /* !!! offset 0 numbering of DOF structs  !!! */
      for(dof_struct_id = 0; dof_struct_id< level_p->nr_dof_ent; dof_struct_id++){
	
	dof_struct_p = &level_p->l_dof_struct[dof_struct_id];

/*kbw
      printf("dof_struct_p %lu, id %d, dof_ent %d\n", dof_struct_p,
	     dof_struct_id, dof_struct_p->dof_ent_id);
      printf("nrdofs %d, posglob %d, nrneig %d\nneig: ",
	     dof_struct_p->nrdofs, dof_struct_p->posglob, dof_struct_p->nrneig);
      for(ineig=0;ineig<dof_struct_p->nrneig;ineig++){
	int iaux=dof_struct_p->l_neig[ineig];
	printf("%6d",level_p->l_dof_struct[iaux].dof_ent_id);
      }
      printf("\nint_ent: ");
	for(i=0;i<dof_struct_p->nr_int_ent;i++){
	  printf("%d  ", dof_struct_p->l_int_ent_index[i]);
	}
	printf("\n");
	getchar();
/*kew*/

	l_bl_nrdofs[dof_struct_id] = dof_struct_p->nrdofs;
	l_bl_posglob[dof_struct_id] = dof_struct_p->posglob;
	l_bl_dof_ent_id[dof_struct_id] = dof_struct_p->dof_ent_id;
	l_bl_dof_ent_type[dof_struct_id] = dof_struct_p->dof_ent_type;
	
      }


/*kbw
    printf("solver: sending to exchange_tables - nrblocks %d\n",
	   level_p->nr_dof_ent);
    //for(ibl=0;ibl<it_level->Nrblocks;ibl++){
    for(ibl=0;ibl<10;ibl++){
      printf("ibl %d, nrdof %d, posglob %d\n",
	     ibl, l_bl_nrdofs[ibl], l_bl_posglob[ibl]);
    }
/*kew*/

      pdr_create_exchange_tables(Problem_id, Solver_id, level_id, 
				 level_p->nr_dof_ent, 
				 l_bl_dof_ent_type, l_bl_dof_ent_id,
				 l_bl_nrdofs, l_bl_posglob, 
				 level_p->l_dof_elem_to_struct,
				 level_p->l_dof_face_to_struct,
				 level_p->l_dof_edge_to_struct,
				 level_p->l_dof_vert_to_struct
				 );

      //exchange_global_row_indices(level_p->nr_dof_ent,l_bl_dof_ent_id,l_bl_nrdofs,l_bl_posglob,level_p->l_dof_vert_to_struct,0);
    }
    
    
    for(istruct = 0; istruct < level_p->nr_dof_ent; istruct++){
      dof_struct_id = istruct;      
      dof_struct_p = &level_p->l_dof_struct[dof_struct_id];
      if(dof_struct_p->nrneig>0) {
	/* !!! offset 1 numbering of blocks moved to lad_block !!! */
	//ibl = dof_struct_p->block_id+1;
	int ibl = dof_struct_p->block_id;
	free(l_bl_l_neig[ibl]);
      }
    }
    if(siv_solver[Solver_id].parallel){
      free(l_bl_dof_ent_id);
      free(l_bl_dof_ent_type);
    }
    free(l_bl_nrdofs);
    free(l_bl_posglob);
    free(l_bl_nrneig);
    free(l_bl_l_neig);

    /* create preconditioner data structure */
    lsr_mkb_create_precon(Solver_id, level_id);


  } // end loop over levels


  // We make a loop over levels to initialize assembly tables
  nr_levels = siv_solver[Solver_id].nr_levels;
  /* in a loop over levels */
  for(level_id=nr_levels-1;level_id>=0;level_id--){
    
    level_p = &(siv_solver[Solver_id].level[level_id]);
    
  level_p->asse_pos_first_dof_int_ent = NULL;
  level_p->assembly_table = NULL;
    
  }

  // PREPARING ASSEMBLY TABLES!!!
  int assembly_table_active = 1;

  if(siv_solver[Solver_id].level[nr_levels-1].block_size<=0) assembly_table_active = 0;
  if(approx_type==SIC_DG_APPROX) assembly_table_active = 0;
  if(storage_type==LSC_STORAGE_BLOCK) assembly_table_active = 0;
  if(storage_type==LSC_STORAGE_PETSC) assembly_table_active = 0;
  //if(storage_type==LSC_STORAGE_UNDEFINED) assembly_table_active = 0;

  // switch off assembly tables for small meshes (the same limit as for coloring!!!!!)
  if(siv_solver[Solver_id].level[nr_levels-1].nrdofs_glob<1000) assembly_table_active = 0;

  printf("Storage type %d, assembly_table_active %d\n",
	 storage_type, assembly_table_active);
  
  //!!!!!!!!!!*********?????????? FOR DEBUGGING ??????????**********!!!!!!!!!!
  //assembly_table_active = 1;

  if(assembly_table_active){

    int l_bl_id[SIC_MAX_DOF_PER_INT];
    
    // We make a loop over levels to create assembly tables
    nr_levels = siv_solver[Solver_id].nr_levels;
    
    /* in a loop over levels */
    for(level_id=nr_levels-1;level_id>=0;level_id--){
      siv_solver[Solver_id].cur_level = level_id;
      
      level_p = &(siv_solver[Solver_id].level[level_id]);
      
      // Assembly tables can now be allocated
      level_p->asse_pos_first_dof_int_ent = 
	(int *) malloc( level_p->nr_int_ent*sizeof(int) );
      level_p->assembly_table = 
	(int *) malloc( level_p->nr_asse_blocks_all_int_ent*sizeof(int) );
      
      
      int current_block = 0; // blocks are numbered consecutively in the whole assembly tables
      for(ient=0; ient<level_p->nr_int_ent;ient++){
	
	level_p->asse_pos_first_dof_int_ent[ient] = current_block;
	
	/*kbw
	  printf("filling part of assembly table for integration entity %d (type %d, id %d)\n", 
	  ient, level_p->l_int_ent_type[ient], level_p->l_int_ent_id[ient]);
	  printf("position of first entry: %d\n", level_p->asse_pos_first_dof_int_ent[ient]);
	  /*kew*/
	
	int idofent, jdofent, i_dof_struct_id, j_dof_struct_id;     
	
	int_ent_id = level_p->l_int_ent_id[ient];
	nr_dof_ent_loc = SIC_MAX_DOF_PER_INT;
	pdeg_coarse = SIC_PDEG_FINEST;
	if(level_id<nr_levels-1) pdeg_coarse = level_p->pdeg_coarse;
	
	pdr_comp_stiff_mat(Problem_id, level_p->l_int_ent_type[ient], 
			   level_p->l_int_ent_id[ient], PDC_NO_COMP, &pdeg_coarse,
			   &nr_dof_ent_loc, l_dof_ent_types, 
			   l_dof_ent_ids, l_dof_ent_nrdofs,
			   NULL, NULL, NULL, NULL);
	
	// should be returned by pdr_comp_stiff_mat :(
	int nrdofs_int_ent = 0;
	for(idofent=0; idofent<nr_dof_ent_loc; idofent++){
	  nrdofs_int_ent += l_dof_ent_nrdofs[idofent];
	}
	
	
	for(idofent=0; idofent<nr_dof_ent_loc; idofent++){
	  
	  int i_dof_ent_id = l_dof_ent_ids[idofent];
	  int i_dof_ent_type = l_dof_ent_types[idofent];
	  
	  if(i_dof_ent_type == PDC_ELEMENT ) {  
	    i_dof_struct_id = level_p->l_dof_elem_to_struct[i_dof_ent_id];
	  } else if(i_dof_ent_type == PDC_FACE ) {
	    i_dof_struct_id = level_p->l_dof_face_to_struct[i_dof_ent_id];
	  } else if(i_dof_ent_type == PDC_EDGE ) {
	    i_dof_struct_id = level_p->l_dof_edge_to_struct[i_dof_ent_id];
	  } else if(i_dof_ent_type == PDC_VERTEX ) {
	    i_dof_struct_id = level_p->l_dof_vert_to_struct[i_dof_ent_id];
	  }
	  
	  /* !!! offset 1 numbering of blocks moved to lad_block !!! */
	  //l_bl_id[idofent] = level_p->l_dof_struct[dof_struct_id].block_id + 1;
	  l_bl_id[idofent] = level_p->l_dof_struct[i_dof_struct_id].block_id;
	  // for the time being assembly tables are for
	  // storage formats with constant block size
	  //l_bl_nrdofs[idofent] = level_p->l_dof_struct[i_dof_struct_id].nrdofs;
	  
	}
	
	
/*kbw
      //if(level_p->l_int_ent_id[intent]==1 || level_p->l_int_ent_id[intent]==4) {
      printf("In sir_solve before filling assembly table: Solver_id %d, level_id %d\n", 
	     Solver_id, level_id);
      printf("ient %d, int_ent_id %d, int_ent_type %d, nr_dof_ent_loc %d, position %d\n", 
	     ient, level_p->l_int_ent_id[ient],  
	     level_p->l_int_ent_type[ient], nr_dof_ent_loc,
	     level_p->asse_pos_first_dof_int_ent[ient]);
      int ibl;
      for(ibl=0;ibl<nr_dof_ent_loc; ibl++) printf("bl_id %d ", l_bl_id[ibl]);
      printf("\n");
/*kew*/

	int position = level_p->asse_pos_first_dof_int_ent[ient];
  //if(ient<10){printf("\n\n ient=%d, position=%d %d\n\n",ient,position, (nr_dof_ent_loc+1)*nr_dof_ent_loc);getchar();}
	lsr_mkb_fill_assembly_table_int_ent(Solver_id, level_id,
					    nr_dof_ent_loc, l_bl_id, NULL, // l_bl_nrdofs,
					    // for the time being assembly tables are for
					    // storage formats with constant block size
					    &level_p->assembly_table[position]);
	
	current_block += (nr_dof_ent_loc+1)*nr_dof_ent_loc;
	
      } // the end of loop over integration entities 
      
/*kbw
      printf("After filling assembly tables at level %d\n", level_id);
      printf("total number of blocks in tables %d (%d)\n",
	     level_p->nr_asse_blocks_all_int_ent, current_block);
/*kew*/

    } // end if assembly tables created for level (block_size>0 and approx_type != DG)
	
  } // the end of loop over levels 


  return(0);
}


/*------------------------------------------------------------
sir_solve - to solve the system for a given data
------------------------------------------------------------*/
int sir_solve(/* returns: >=0 - success code, <0 - error code */
  int Solver_id,     /* in: solver identification */
  int Comp_type,     /* in: indicator for the scope of computations: */
                     /*   SIC_SOLVE - solve the system */
                     /*   SIC_RESOLVE - resolve for the new right hand side */
  int Ini_guess,     /* in: indicator on whether to set initial guess (>0), */
                     /*     or to initialize it to zero (0) */
         /* if >0 then it indicates from which solution vector to take data */
  int Monitor,       /* in: monitoring flag with options (-1 for defaults): */
                     /*   SIC_PRINT_NOT - do not print anything */ 
                     /*   SIC_PRINT_ERRORS - print error messages only */
                     /*   SIC_PRINT_INFO - print most important information */
                     /*   SIC_PRINT_ALLINFO - print all available information */
  int* Nr_iter,      /* in:	the maximum iterations to be performed */
                     /* out:	actual number of iterations performed */
  double* Conv_meas, /* in:	tolerance level for chosen measure */
		     /* out:	the final value of convergence measure */
  double *Conv_rate  /* out (optional): the total convergence rate */ 
  )
{

  /* pointer to solver structure */
  sit_solvers *solver_p;
  sit_levels *level_p; /* mesh levels */

  /* pointer to dofs structure */
  sit_dof_struct *dof_struct_p;
  
  /* auxiliary variables */
  int nrdofs_glob, max_nrdofs, nr_levels, posglob, nrdofs_int_ent;
  int l_dof_ent_id[SIC_MAX_DOF_PER_INT], l_dof_ent_nrdof[SIC_MAX_DOF_PER_INT];
  int l_dof_ent_posglob[SIC_MAX_DOF_PER_INT];
  int l_dof_ent_type[SIC_MAX_DOF_PER_INT];
  int pdr_comp_type;
  double *stiff_mat, *rhs_vect, *x_ini, normb;
  int i,j,k, iaux, kaux, intent, ient;
  int level_id=0;
  char rewrite;
  int pdeg_coarse;

/* variables to store timing results */
  double t_int_el=0.0;
  //double t_int_fa=0.0;
  double t_fac_dia=0.0;
  double t_iter=0.0;
  double t_temp=0.0;
  double t_total=0.0;
  //double time_clock();

/*++++++++++++++++ executable statements ++++++++++++++++*/

  siv_cur_solver_id=Solver_id;
  solver_p = &siv_solver[siv_cur_solver_id];

  nr_levels = siv_solver[Solver_id].nr_levels;
  level_id = nr_levels-1;
  level_p = &(siv_solver[Solver_id].level[level_id]);

  /*ok_kbw*/
  printf("Entering sir_solve: solver %d, nr_levels %d, level_id %d, nrdofs_glob %d\n",
     Solver_id, nr_levels, level_id, level_p->nrdofs_glob);
  printf("Monitor %d, Nr_iter %d, Conv_meas %15.12lf\n",
	 Monitor, *Nr_iter, *Conv_meas);
  /*kew*/


#ifdef MKB_SOLVE_TIME_TEST
  t_total=time_clock();
#endif

  /* allocate memory for the initial guess and solution vectors */
  nrdofs_glob = siv_solver[Solver_id].level[nr_levels-1].nrdofs_glob;
  x_ini = (double *)malloc(nrdofs_glob*sizeof(double));
  int ini_zero=1; // default: initialize initial guess with zeros
  for(i=0;i<nrdofs_glob;i++) x_ini[i]=0.0;
  
/*------------------------------------------------------------*/
/* GET INITIAL GUESS WHEN POSSIBLE AND NECESSARY              */
/*------------------------------------------------------------*/
  if(Ini_guess>0){
    
    ini_zero=0;
    
    int idofent;
    /* get initial guess */
    for(idofent=0;idofent<level_p->nr_dof_ent;idofent++){
      int dof_struct_id = idofent;
      sit_dof_struct dof_struct = level_p->l_dof_struct[dof_struct_id];

  /*kbw
      printf("in sir_solve before pdr_read_sol_dofs: idofent %d, struct_id %d \n",
	     idofent, dof_struct_id);
      printf("dof_struct_p %lu\n", dof_struct_p);
  printf("problem_id %d, Dof_ent_type %d, Dof_ent_id %d, nrdof %d\n",
	solver_p->problem_id, dof_struct.dof_ent_type, dof_struct.dof_ent_id,  
	 dof_struct.nrdofs);
  /*kew*/

      
      pdr_read_sol_dofs(solver_p->problem_id, 
			dof_struct.dof_ent_type,
			dof_struct.dof_ent_id,
			dof_struct.nrdofs,
			&x_ini[dof_struct.posglob]);
      
/*kbw
      if(idofent<5){
printf("Initial guess in dof_ent %d (local %d)\n", dof_struct.dof_ent_id, idofent);
for (i=0;i<dof_struct.nrdofs;i++) printf("%20.10lf",x_ini[dof_struct.posglob+i]);
printf("\n");
    getchar();
      }
/*kew*/

    }
    
/*kbw
    //if(siv_solver[Solver_id].parallel){
      printf("x_ini before exchange dofs\n");
      for(i=0;i<nrdofs_glob;i++) printf("%5d%15.10lf",i/4+1,x_ini[i]);
      printf("\n");
      getchar();
      //}
/*kew*/

/*||begin||
    if(siv_solver[Solver_id].parallel){
      pdr_exchange_dofs(solver_p->problem_id, Solver_id, nr_levels-1, x_ini);
    }
/*||end||*/



/*kbw
    //if(siv_solver[Solver_id].parallel){
      printf("x_ini after exchange dofs\n");
      for(i=0;i<nrdofs_glob;i++) printf("%5d%15.10lf",i/4+1,x_ini[i]);
      printf("\n");
      getchar();
      //}
/*kew*/
  
  }


#ifdef MULTITHREADED
  // Problem module controls assembling - one call for many elements - OpenMP possible
  int Assemble_control = SIC_PROBLEM_ASSEMBLE;
  //int Assemble_control = SIC_SOLVER_ASSEMBLE;
#else
  // Solver controls assembling - one solver thread per one problem thread - no OpenMP
  int Assemble_control = SIC_SOLVER_ASSEMBLE;
  //int Assemble_control = SIC_PROBLEM_ASSEMBLE;
#endif
  
#ifdef MKB_SOLVE_TIME_TEST
  printf("\nbeginning integration of %d elements and faces\n",
	 level_p->nr_int_ent);
  t_temp=time_clock();
#endif    

    // simple test for type of approximation
  int approx_type = SIC_STD_APPROX;
  if(siv_solver[Solver_id].level[nr_levels-1].max_dof_elem_id>=0) {
    approx_type = SIC_DG_APPROX;
  }


  if(Assemble_control==SIC_SOLVER_ASSEMBLE || approx_type==SIC_DG_APPROX){ // SEQUENTIAL!!!
  
    /* allocate memory for the stiffness matrices and RHS */
    max_nrdofs = level_p->max_dofs_int_ent;
    stiff_mat = (double *)malloc(max_nrdofs*max_nrdofs*sizeof(double));
    rhs_vect = (double *)malloc(max_nrdofs*sizeof(double));
    
    /* for each level */
    for(level_id=nr_levels-1;level_id>=0;level_id--){

      level_p = &(siv_solver[Solver_id].level[level_id]);

      lsr_mkb_clear_matrix(Solver_id, level_id, Comp_type);   
	
      pdeg_coarse = SIC_PDEG_FINEST;
      if(level_id<nr_levels-1) pdeg_coarse = level_p->pdeg_coarse;

      // some day we can try coloring and OpenMP for loop here...
      for(intent=0;intent<level_p->nr_int_ent;intent++){
	
	int idofent;
	int nr_dof_ent_loc = SIC_MAX_DOF_PER_INT;
	int nrdofs_int_ent = max_nrdofs;
	
	/* compute local stiffness matrices */
	if(Comp_type==SIC_SOLVE) pdr_comp_type = PDC_COMP_BOTH;
	else pdr_comp_type = PDC_COMP_RHS;
	
	pdr_comp_stiff_mat(solver_p->problem_id, level_p->l_int_ent_type[intent], 
			   level_p->l_int_ent_id[intent], pdr_comp_type, &pdeg_coarse,
			   &nr_dof_ent_loc,l_dof_ent_type,l_dof_ent_id,l_dof_ent_nrdof, 
			   &nrdofs_int_ent, stiff_mat, rhs_vect, &rewrite);
	
	
#ifdef DEBUG_SIM
	if(nrdofs_int_ent>max_nrdofs){
	  printf("Too small arrays stiff_mat and rhs_vect passed to comp_el_stiff_mat\n");
	  printf("from sir_create. %d < %d. Exiting !!!",
		 max_nrdofs, nrdofs_int_ent);
	  exit(-1);
	}
#endif
	
    
	
/*kbw
    //if(level_p->l_int_ent_id[intent]==1 || level_p->l_int_ent_id[intent]==4) {
	int position = level_p->asse_pos_first_dof_int_ent[intent];
      printf("In sir_solve before assemble: Solver_id %d, level_id %d, nr_dof_ent %d, nr_dofs_glob %d\n", 
	     Solver_id, level_id, level_p->nr_dof_ent, level_p->nrdofs_glob);
      int ibl,jbl,pli,plj,nri,nrj,nrdof,jaux;
      printf("ient %d, int_ent_id %d, int_ent_type %d, nr_dof_ent_loc %d\n", 
	     intent, level_p->l_int_ent_id[intent],  
	     level_p->l_int_ent_type[intent], nr_dof_ent_loc);
      pli = 0; nrdof=0;
      for(ibl=0;ibl<nr_dof_ent_loc; ibl++) nrdof+=l_dof_ent_nrdof[ibl];
      for(ibl=0;ibl<nr_dof_ent_loc; ibl++){
	printf("bl_id %d, bl_nrdof %d\n",
	  l_dof_ent_id[ibl],l_dof_ent_nrdof[ibl]);
	nri = l_dof_ent_nrdof[ibl];
	plj=0;
	for(jbl=0;jbl<nr_dof_ent_loc;jbl++){
	  //printf("Stiff_mat transposed! (blocks %d:%d -> group of rows %d, group of columns %d)\n",
	  //	 jbl,ibl,jbl,ibl);
	  //nrj = l_dof_ent_nrdof[jbl];
	  //for(i=0;i<nri;i++){
   	  //  jaux = plj+(pli+i)*nrdof;
	  //  for(j=0;j<nrj;j++){
	  //    printf("%20.15lf",stiff_mat[jaux+j]);
	  //  }
	  //  printf("\n");
	  //}
	  printf("Stiff_mat  (blocks %d:%d -> group of rows %d, group of columns %d)\n",
		 jbl,ibl,jbl,ibl);
	  nrj = l_dof_ent_nrdof[jbl];
	  for(j=0;j<nrj;j++){
	    for(i=0;i<nri;i++){
	      jaux = plj+(pli+i)*nrdof;
	      printf("%20.15lf",stiff_mat[jaux+j]);
	    }
	    printf("\n");
	  }
	  plj += nrj;
	}
	nri = l_dof_ent_nrdof[ibl];
	plj=0;
	for(jbl=0;jbl<nr_dof_ent_loc;jbl++){
	  printf("Assembly_table  (blocks %d:%d -> group of rows %d, group of columns %d)\n",
		 jbl,ibl,jbl,ibl);
	  nrj = l_dof_ent_nrdof[jbl];
	  for(j=0;j<nrj;j++){
	    for(i=0;i<nri;i++){
	      jaux = plj+(pli+i)*nrdof;
	      printf("%20d",level_p->assembly_table[position+jaux+j]);
	    }
	    printf("\n");
	  }
	  plj += nrj;
	}
/*kbw
    //if(level_p->l_int_ent_id[intent]==1 || level_p->l_int_ent_id[intent]==4) {

	int position = level_p->asse_pos_first_dof_int_ent[intent];
	if(level_p->assembly_table[position+nr_dof_ent_loc*nr_dof_ent_loc+4]==36){
      printf("In sir_solve before assemble: Solver_id %d, level_id %d, nr_dof_ent %d, nr_dofs_glob %d\n", 
	     Solver_id, level_id, level_p->nr_dof_ent, level_p->nrdofs_glob);
      int ibl,jbl,pli,plj,nri,nrj,nrdof,jaux;
      printf("ient %d, int_ent_id %d, int_ent_type %d, nr_dof_ent_loc %d\n", 
	     intent, level_p->l_int_ent_id[intent],  
	     level_p->l_int_ent_type[intent], nr_dof_ent_loc);
      pli = 0;
      for(ibl=0;ibl<nr_dof_ent_loc; ibl++){
	printf("bl_id %d, bl_nrdof %d\n",
	  l_dof_ent_id[ibl],l_dof_ent_nrdof[ibl]);
	nri = l_dof_ent_nrdof[ibl];
	printf("Rhs_vect (block %d -> group of rows %d):\n", ibl, ibl);
	for(i=0;i<nri;i++){
	  printf("%20.15lf\n",rhs_vect[pli+i]);
	}
	printf("\n");
	pli += nri;    
      }
      getchar();
      }
/*kew*/
    
      /* printf("In sir_solve before assemble: Solver_id %d, level_id %d, nr_dof_ent %d, nr_dofs_glob %d\n",  */
      /* 	     Solver_id, level_id, level_p->nr_dof_ent, level_p->nrdofs_glob); */
      /* printf("ient %d, int_ent_id %d, int_ent_type %d, nr_dof_ent_loc %d\n",  */
      /* 	     intent, level_p->l_int_ent_id[intent],   */
      /* 	     level_p->l_int_ent_type[intent], nr_dof_ent_loc); */
 
	//level_p->assembly_table = NULL; // for testing block solver !!!!!!!!!!!

	if(level_p->assembly_table != NULL){
	  //if(level_p->l_int_ent_type[intent]==3){

	  int position = level_p->asse_pos_first_dof_int_ent[intent];
	  sir_assemble_local_stiff_mat_with_table(Solver_id, level_id, pdr_comp_type,
						      nr_dof_ent_loc,
						      &level_p->assembly_table[position],
						      stiff_mat, rhs_vect, &rewrite);
	  
	}
	else{
	  sir_assemble_local_stiff_mat(Solver_id, level_id, pdr_comp_type,
				       nr_dof_ent_loc, l_dof_ent_type,
				       l_dof_ent_id, l_dof_ent_nrdof, 
				       stiff_mat, rhs_vect, &rewrite);

//	  if(level_p->l_int_ent_type[intent] == PDC_ELEMENT){
//			pdr_comp_stiff_mat(solver_p->problem_id, level_p->l_int_ent_type[intent],
//					   level_p->l_int_ent_id[intent], PDC_COMP_MM, &pdeg_coarse,
//					   &nr_dof_ent_loc,l_dof_ent_type,l_dof_ent_id,l_dof_ent_nrdof,
//					   &nrdofs_int_ent, stiff_mat, rhs_vect, &rewrite);
//
//			  sir_assemble_local_stiff_mat(Solver_id, level_id, PDC_COMP_MM,
//							   nr_dof_ent_loc, l_dof_ent_type,
//							   l_dof_ent_id, l_dof_ent_nrdof,
//							   stiff_mat, rhs_vect, &rewrite);
//	  }
	}

      } /* end loop over integration entities: ient */
      
    } /* end loop over levels: level_id */

  } // end if solver controls assembly

  // if problem module controls assembly
  else{

    /* for each level */
    for(level_id=nr_levels-1;level_id>=0;level_id--){

      level_p = &(siv_solver[Solver_id].level[level_id]);

      lsr_mkb_clear_matrix(Solver_id, level_id, Comp_type);
    
      int pdeg_coarse = SIC_PDEG_FINEST;
      if(level_id<nr_levels-1) pdeg_coarse = level_p->pdeg_coarse;

      if(Comp_type==SIC_SOLVE) pdr_comp_type = PDC_COMP_BOTH;
      else pdr_comp_type = PDC_COMP_RHS;

/*kbw
  printf("before  pdr_create_assemble_stiff_mat: Problem_id %d, Level_id %d\n",
	 solver_p->problem_id, level_id);
  printf("Comp_type %d, Nr_int_ent %d, Max_dofs_int_ent %d\n",
	 pdr_comp_type, level_p->nr_int_ent, level_p->max_dofs_int_ent);
/*kew*/

      pdr_create_assemble_stiff_mat(solver_p->problem_id, level_id, 
				    pdr_comp_type, &pdeg_coarse,
				    level_p->nr_int_ent,
				    level_p->l_int_ent_type, 
				    level_p->l_int_ent_id,
				    level_p->nr_colors_elems, level_p->l_color_index_elems,
				    level_p->nr_colors_faces, level_p->l_color_index_faces,
				    level_p->asse_pos_first_dof_int_ent, 
				    level_p->assembly_table, 
				    level_p->max_dofs_int_ent);
    

      
    }

  }
  

    
#ifdef MKB_SOLVE_TIME_TEST
  t_int_el+=time_clock()-t_temp;
  printf("beginning factorization of diagonal blocks\n");
  t_temp=time_clock();
#endif

  for(level_id=nr_levels-1;level_id>=0;level_id--){
    lsr_mkb_fill_precon(Solver_id, level_id);
    
  }

#ifdef MKB_SOLVE_TIME_TEST
  t_fac_dia+=time_clock()-t_temp;
#endif
  



  /*kbw
    printf("Initial guess:\n");
    for(i=0;i<nrdofs_glob;i++)printf("%20.15lf",x_ini[i]);
    printf("\n");
    getchar();getchar();
/*kew*/


#ifdef MKB_SOLVE_TIME_TEST
  printf("beginning iterations\n");
  t_temp=time_clock();
#endif

/*------------------------------------------------------------*/
/* SOLVE THE PROBLEM                                          */
/*------------------------------------------------------------*/



  lsr_mkb_solve(Solver_id, Comp_type, ini_zero, x_ini, rhs_vect,
		Nr_iter, Conv_meas, Monitor, Conv_rate);

#ifdef MKB_SOLVE_TIME_TEST
  t_iter+=time_clock()-t_temp;
#endif



  /* rewrite the solution */
  level_id = nr_levels-1;
  level_p = &(siv_solver[Solver_id].level[level_id]);
  int istruct;
  for(istruct=0;istruct<level_p->nr_dof_ent;istruct++){
    int dof_struct_id = istruct;
    sit_dof_struct dof_struct = level_p->l_dof_struct[dof_struct_id];
    
    /* print out solution
    //if(istruct<10){
    //if(dof_struct.dof_ent_id<10){
      printf("Solution in struct %d, block %d, dof_ent %d, posglob=%d\n", 
	     istruct, dof_struct.block_id, dof_struct.dof_ent_id, dof_struct.posglob);
      for (i=0;i<dof_struct.nrdofs;i++) {
	printf("%20.10lf",x_ini[dof_struct.posglob+i]);
	//printf("%20.10lf",x_ini[dof_struct_id+i]);
      }
      printf("\n");
//    getchar();
//    }
/**/
    
    pdr_write_sol_dofs(solver_p->problem_id, 
		       dof_struct.dof_ent_type,
		       dof_struct.dof_ent_id,
		       dof_struct.nrdofs,
		       &x_ini[dof_struct.posglob]);



  }

  if(Assemble_control == SIC_SOLVER_ASSEMBLE){
    free(stiff_mat);
    free(rhs_vect);
  }
  free(x_ini);

/*ok_kbw*/
#ifdef MKB_SOLVE_TIME_TEST
  t_total=time_clock()-t_total;
  printf("\nTotal solver times:\n");
  printf("\tintegration of elements and faces \t%lf\n", t_int_el);
  printf("\tfactorization of diagonal blocks  \t%lf\n", t_fac_dia);
  printf("\titerations / direct solution      \t%lf\n", t_iter);
  t_temp=t_int_el+t_fac_dia+t_iter;
  printf("\tsuma                        \t%lf\n", t_temp);
  printf("\ttotal                       \t%lf\n", t_total);
  //printf("%lf  %lf  %lf  %lf  %lf  %lf\n",
  //	 t_int_el,t_fac_dia,t_iter,t_temp,t_total);
#endif
/*kew*/


  return(1);
}


/*------------------------------------------------------------
  sir_assemble_local_stiff_mat_with_table - to assemble an element stiffness matrix
                                   to the global SM using assembly table
------------------------------------------------------------*/
int sir_assemble_local_stiff_mat_with_table( 
                         /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id,          /* in: level ID */
  int Comp_type,         /* in: indicator what was computed: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int Nr_dof_ent_loc,        /* in: number of global dof blocks */
                         /*     associated with the local stiffness matrix */
  int* Assembly_table_int_ent, /* part of the global assembly table */
  double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /* in: rhs vector */
  char* Rewr_dofs         /* in: flag to rewrite or sum up entries */
                         /*   'T' - true, rewrite entries when assembling */
                         /*   'F' - false, sum up entries when assembling */
)
{

  int solver_comp_type;
  if(Comp_type==PDC_COMP_BOTH) solver_comp_type = SIC_SOLVE;
  else if(Comp_type==PDC_COMP_RHS) solver_comp_type = SIC_RESOLVE;
  else{
    printf("How the hell you want to solve something without SM and RHS?\n");
    exit(-1);
  }

  
  return lsr_mkb_assemble_local_stiff_mat_with_table(Solver_id, Level_id, solver_comp_type,
					      Nr_dof_ent_loc, Assembly_table_int_ent,
					      Stiff_mat, Rhs_vect, Rewr_dofs);

  
}


/*------------------------------------------------------------
  sir_assemble_local_stiff_mat - to assemble an element stiffness matrix
                                   to the global SM
------------------------------------------------------------*/
int sir_assemble_local_stiff_mat( 
                         /* returns: >=0 - success code, <0 - error code */
  int Solver_id,         /* in: solver ID (used to identify the subproblem) */
  int Level_id,          /* in: level ID */
  int Comp_type,         /* in: indicator what was computed: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int Nr_dof_ent_loc,        /* in: number of global dof blocks */
                         /*     associated with the local stiffness matrix */
  int* L_dof_ent_type,   /* in: list of dof blocks' IDs */
  int* L_dof_ent_id,     /* in: list of dof blocks' IDs */
  int* L_dof_ent_nrdofs, /* in: list of blocks' numbers of dof */
  double* Stiff_mat,     /* in: stiffness matrix stored columnwise */
  double* Rhs_vect,      /* in: rhs vector */
  char* Rewr_dofs         /* in: flag to rewrite or sum up entries */
                         /*   'T' - true, rewrite entries when assembling */
                         /*   'F' - false, sum up entries when assembling */
)
{
  int l_bl_id[SIC_MAX_DOF_PER_INT], l_bl_nrdofs[SIC_MAX_DOF_PER_INT];

  sit_levels *level_p;       /* mesh levels */
  level_p = &siv_solver[Solver_id].level[Level_id];
  
  int solver_comp_type;
  if(Comp_type==PDC_COMP_BOTH) solver_comp_type = SIC_SOLVE;
  else if(Comp_type==PDC_COMP_RHS) solver_comp_type = SIC_RESOLVE;
//  else if(Comp_type == PDC_COMP_MM){
//	  solver_comp_type = PDC_COMP_MM;
//  }
  else {
    printf("How the hell you want to solve something without SM and RHS?\n");
    exit(-1);
  }
  
  int idofent;
  int nrdofbl = Nr_dof_ent_loc;
  for(idofent=0; idofent<nrdofbl; idofent++){
    int dof_ent_id = L_dof_ent_id[idofent];
    int dof_ent_type = L_dof_ent_type[idofent];
    int dof_struct_id;
    
    if(dof_ent_type == PDC_ELEMENT ) {   
      dof_struct_id = level_p->l_dof_elem_to_struct[dof_ent_id];
    } else if(dof_ent_type == PDC_FACE ) {
      dof_struct_id = level_p->l_dof_face_to_struct[dof_ent_id];
    } else if(dof_ent_type == PDC_EDGE ) {
      dof_struct_id = level_p->l_dof_edge_to_struct[dof_ent_id];
    } else if(dof_ent_type == PDC_VERTEX ) {
      dof_struct_id = level_p->l_dof_vert_to_struct[dof_ent_id];
    }
    else{
      printf("Unknown dof_ent_type %d in sir_assemble!!! Exiting\n",
	     dof_ent_type);
	exit(-1);
    }
    
    /* !!! offset 1 numbering of blocks moved to lad_block !!! */
    //l_bl_id[idofent] = level_p->l_dof_struct[dof_struct_id].block_id + 1;
    l_bl_id[idofent] = level_p->l_dof_struct[dof_struct_id].block_id;
    l_bl_nrdofs[idofent] = level_p->l_dof_struct[dof_struct_id].nrdofs;
    
  }
  
  
/*kbw
  //if(level_p->l_int_ent_id[intent]==1 || level_p->l_int_ent_id[intent]==4) {
#pragma omp critical(printing)
  {
      printf("In sir_assemble_local_stiff_mat: Solver_id %d, Level_id %d, sol_typ %d\n", Solver_id, Level_id, solver_comp_type);
      int i,j,ibl,jbl,pli,plj,nri,nrj,nrdof,jaux;
      printf("nrdofbl %d\n", nrdofbl);
      pli = 0; nrdof=0;
      for(ibl=0;ibl<nrdofbl; ibl++) nrdof+=l_bl_nrdofs[ibl];
      for(ibl=0;ibl<nrdofbl; ibl++){
	printf("bl_id %d, bl_nrdof %d\n",  l_bl_id[ibl],l_bl_nrdofs[ibl]);
      }
      for(ibl=0;ibl<nrdofbl; ibl++){
	printf("bl_id %d, bl_nrdof %d\n",
	  l_bl_id[ibl],l_bl_nrdofs[ibl]);
	nri = l_bl_nrdofs[ibl];
	plj=0;
	for(jbl=0;jbl<nrdofbl;jbl++){
	  printf("Stiff_mat (blocks %d:%d)\n",jbl,ibl);
	  nrj = l_bl_nrdofs[jbl];
	  for(i=0;i<nri;i++){
   	    jaux = plj+(pli+i)*nrdof;
	    for(j=0;j<nrj;j++){
	      printf("%20.15lf",Stiff_mat[jaux+j]);
	    }
	    printf("\n");
	  }
	  plj += nrj;
	}
	printf("Rhs_vect:\n");
	for(i=0;i<nri;i++){
	  printf("%20.15lf",Rhs_vect[pli+i]);
	}
	printf("\n");
	pli += nri;    
      }
      getchar();
      }
/*kew*/
    

  lsr_mkb_assemble_local_sm(Solver_id, Level_id, solver_comp_type,
			    nrdofbl, l_bl_id, l_bl_nrdofs, 
			    Stiff_mat, Rhs_vect, Rewr_dofs);
  
  return(1);
}

/*------------------------------------------------------------
  sir_free - to free memory for stiffness and preconditioner matrices
             and make room for next solvers
------------------------------------------------------------*/
int sir_free(/* returns: >=0 - success code, <0 - error code */
  int Solver_id   /* in: solver identification */
  )
{

  int nr_levels, level_id, ilev;
  sit_levels *level_p;       /* mesh levels */

/*++++++++++++++++ executable statements ++++++++++++++++*/

  siv_cur_solver_id=Solver_id;

/* free solver data structures */
  //lsr_mkb_free_precon(Solver_id); - may be should be separated...
  lsr_mkb_free_matrix(Solver_id);

  nr_levels = siv_solver[Solver_id].nr_levels;

  /* in a big loop over mesh levels free renumbering data structures */
  for(ilev=0;ilev<nr_levels;ilev++){
    
    level_id = nr_levels - 1 - ilev;
    level_p = &siv_solver[Solver_id].level[level_id];

    // tables created by pdr_get_list_ent...    
    free(level_p->l_int_ent_type);
    free(level_p->l_int_ent_id);

    if(level_p->max_dof_elem_id>=0) free(level_p->l_dof_elem_to_struct);
    if(level_p->max_dof_face_id>=0) free(level_p->l_dof_face_to_struct);
    if(level_p->max_dof_edge_id>=0) free(level_p->l_dof_edge_to_struct);
    if(level_p->max_dof_vert_id>=0) free(level_p->l_dof_vert_to_struct);

    // Assembly tables
    if(level_p->asse_pos_first_dof_int_ent!=NULL){
      free(level_p->asse_pos_first_dof_int_ent);
      free(level_p->assembly_table);
    }

    // Coloring tables
    if(level_p->l_color_index_elems != NULL) free(level_p->l_color_index_elems);
    if(level_p->l_color_index_faces != NULL) free(level_p->l_color_index_faces);

    free(level_p->l_dof_struct);
    
  }

  return(0);
  
}

/*------------------------------------------------------------
  sir_destroy - to make room for next solvers - in LIFO manner !!!!!
------------------------------------------------------------*/
int sir_destroy(/* returns: >=0 - success code, <0 - error code */
  int Solver_id   /* in: solver identification */
  )
{

  siv_cur_solver_id=Solver_id;

  /* destroy the solver instance */
  lsr_mkb_destroy(Solver_id);

  /* decrease the counter for solvers */
  siv_nr_solvers--;

  /* set the current solver ID */
  if(siv_cur_solver_id == siv_nr_solvers) siv_cur_solver_id = siv_nr_solvers-1;
  else{
    printf("Solver destroyed is not the last solver in sir_destroy!!!\n");
    exit(0);
  }

  return(1);
}


/*---------------------------------------------------------
sir_put_list - to put Num on the list List with length Ll 
	(filled with numbers and SIC_LIST_END_MARK at the end)
---------------------------------------------------------*/
int sir_put_list( /* returns*/
		/*  >0 - position already occupied on the list */
             	/*  <0 - position at which put on the list */
            	/*   0 - list full, not found on the list */
	int Num, 	/* in: number to put on the list */
	int* List, 	/* in: list */
	int Ll		/* in: total list's lengths */
	)
{

  int i, il;
  
  for(i=0;i<Ll;i++){
    if((il=List[i])==SIC_LIST_END_MARK) break;
    /* found on the list on (i+1) position */
    if(Num==il) return(i+1);
  }
  /* if list is full return error message */
  if(i==Ll) return(0);
  /* update the list and return*/
  List[i]=Num;
  return(-(i+1));
}

/*---------------------------------------------------------
  sir_sort_short - to sort an array by insertion (for small integer arrays)
---------------------------------------------------------*/
int sir_sort_short(
		   int* A, // array
		   int p,  //first index
		   int k   // last index
		   )
{
  int i,j;
  int t;
  
  for(i=p+1;i<=k;i++) {
    t=A[i];
    j=i-1;
    while( j>=p && A[j]>t ) {
      A[j+1]=A[j];
      j--;
    }
    A[j+1]=t;
  }
  return 0;
}
 

#ifdef INTERNAL_RENUMBERING

/*------------------------------------------------------------
  sir_reverse_cuthill_mckee - bandwidth reduction algorithm (renumbering)
------------------------------------------------------------*/
typedef struct list {
	int row;
	struct list *next;
}list_row;


void sir_reverse_cuthill_mckee(sit_dof_struct * L_dof_struct, int Nr_dof_ent){
  
  
  list_row *head=NULL, *tail=NULL, *lhead=NULL, *ltail=NULL, *elem=NULL, *lelem=NULL;
  int sum_nreig, icrs=0, jcrs=0;
  
  /* array of visited DOF - 1 visited, 0 not visited */
  int *vt=NULL, *permutation=NULL, *rpermutation=NULL;
  
  permutation =(int *)malloc( Nr_dof_ent*sizeof(int) );
  rpermutation =(int *)malloc( Nr_dof_ent*sizeof(int) );
  if(permutation==NULL){
    printf("no memory in sir_reverse_cuthill_mckee permutation"); 
    exit(0);
  }
  if(rpermutation==NULL){
    printf("no memory in sir_reverse_cuthill_mckee rpermutation"); 
    exit(0);
  }
  
  vt =(int *)malloc( Nr_dof_ent*sizeof(int));
  if(vt==NULL){
    printf("no memory in sir_reverse_cuthill_mckee vt"); 
    exit(0);
  }
  
  
  sit_dof_struct *dof_struct_p;
  int idofent,ineig,i;
  int iaux;
  
  //printf("Nr_dof_ent=%d",Nr_dof_ent);

  for(idofent = 0; idofent< Nr_dof_ent; idofent++){vt[idofent]=0;}
  
  int start_dof=0;
  sum_nreig=L_dof_struct[start_dof].nrneig;
  
  for(idofent = 1; idofent< Nr_dof_ent; idofent++){ 
    sum_nreig+=L_dof_struct[idofent].nrneig;
    if(L_dof_struct[idofent].nrneig <= L_dof_struct[start_dof].nrneig) {start_dof=idofent;}   
  }
  
  for(idofent = 0; idofent< Nr_dof_ent; idofent++)rpermutation[idofent]=-1;
  
  
  // printf("start_dof %d, min %d\n", start_dof, L_dof_struct[start_dof].nrneig);
  
  
  /* start compute permutation vector */
  vt[start_dof]=1;
  int iper=Nr_dof_ent-1;
  permutation[iper]=start_dof;
  rpermutation[permutation[iper]]=iper;
  iper=iper-1;
  while(iper>=0){
    for(ineig=0;ineig<L_dof_struct[start_dof].nrneig;ineig++){
      iaux=L_dof_struct[start_dof].l_neig[ineig];
            
      if(vt[iaux]==0){
	vt[iaux]=1;
	elem =(list_row*)malloc(sizeof(list_row));
	elem->row=iaux;
	elem->next=NULL;
	if(lhead==NULL){
	  lhead = elem;
	  ltail = lhead;
	}
	else{
	  if(L_dof_struct[iaux].nrneig < L_dof_struct[lhead->row].nrneig){
	    elem->next = lhead;
	    lhead = elem;
	  }else{
	    lelem = lhead;
	    while((lelem->next != NULL) && 
		  (L_dof_struct[iaux].nrneig >= L_dof_struct[lelem->row].nrneig)){
	      lelem = lelem->next;
	    }
	    elem->next = lelem->next;
	    lelem->next = elem;
	    if(lelem == ltail)ltail=elem;
	  }
	}
      }
      
      //printf("%6d",L_dof_struct[iaux].dof_ent_id);
    }
    
    if(head == NULL) { head=lhead; tail=ltail; lhead=NULL; ltail=NULL; }
    else if(lhead!=NULL){ tail->next = lhead; tail=ltail; lhead=NULL; ltail=NULL;}
    
    start_dof=head->row;  
    elem=head;
    head=head->next;
    free(elem);
    permutation[iper]=start_dof;
    rpermutation[permutation[iper]]=iper;
    iper=iper-1;
    
    
  }	
  
  
  
  /*
    printf("permutation:%d \n",Nr_dof_ent);
    for(idofent = 0; idofent< Nr_dof_ent; idofent++)printf("%d ",permutation[idofent]);
    printf("\nrpermutation:\n");
    for(idofent = 0; idofent< Nr_dof_ent; idofent++)printf("%d ",rpermutation[idofent]);
    /**/
  
  // setting block_id, pos_glob, l_neig_bl 
  
  start_dof=0;
  for(idofent = 0; idofent< Nr_dof_ent; idofent++){
    // block_id is offset 0 !!!
    L_dof_struct[rpermutation[idofent]].block_id = idofent;
    L_dof_struct[rpermutation[idofent]].posglob=start_dof;
   start_dof+=L_dof_struct[rpermutation[idofent]].nrdofs;
   
  }
  
  
  
  free(vt);
  free(permutation);
  free(rpermutation);
  
}

#endif
