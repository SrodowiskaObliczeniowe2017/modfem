/************************************************************************
File pds_conv_diff_time.c - time integration schemes for convection-
                            diffusion problems

Contains routines:
  pdr_conv_diff_time_integration - time integration driver.

  pdr_conv_diff_sol_diff_norm - returns max difference between current and old
			  solutions vectors

------------------------------  			
History:    
	02.2002 - Krzysztof Banas, initial version		
	2012    - Krzysztof Banas (pobanas@cyf-kr.edu.pl)
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<string.h>


/* USED DATA STRUCTURES AND INTERFACES FROM OTHER MODULES */

/* interface for all mesh manipulation modules */
#include "mmh_intf.h"		/* USES */
/* interface for all approximation modules */
#include "aph_intf.h"		/* USES */
/* interface for all solver modules */
#include "sih_intf.h"		/* USES */

#ifdef PARALLEL
/* interface of parallel mesh manipulation modules */
#include "mmph_intf.h"		/* USES */
/* interface for all parallel approximation modules */
#include "apph_intf.h"		/* USES */
/* interface for parallel communication modules */
#include "pch_intf.h"		/* USES */
#endif

/* utilities - including simple time measurement library */
#include "uth_intf.h"

/* interface for thread management modules */
#include "tmh_intf.h"


/* USED AND IMPLEMENTED PROBLEM DEPENDENT DATA STRUCTURES AND INTERFACES */

/* problem dependent module interface */
#include "pdh_intf.h"

/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include "pdh_control_intf.h"

/* header files for the problem dependent module for Conv_Diff's equation */
#include "../include/pdh_conv_diff.h"		/* USES & IMPLEMENTS */


/**************************************/
/* INTERNAL PROCEDURES                */
/**************************************/
/* Rules:
/* - name always begins with pdr_ */
/* - argument names start uppercase */

/*------------------------------------------------------------
pdr_conv_diff_sol_diff_norm - returns max difference in current and old
			  solutions vectors
------------------------------------------------------------*/
void pdr_conv_diff_sol_diff_norm(
  int Current, // in: current vector id 
  int Old,     // in: old vector id
  double* sol_diff_norm_conv_diff_p // norm of difference current-old for conv_diff
  )
{

  double sol_dofs_current[APC_MAXELSD];	/* solution dofs */
  double sol_dofs_old[APC_MAXELSD];	/* solution dofs */
  int field_id, mesh_id;
  int node_id = 0;
  double temp, norm_conv_diff = 0.0;
  int i;

/*++++++++++++++++ executable statements ++++++++++++++++*/
  
  i=3; field_id = pdr_ctrl_i_params(pdv_conv_diff_current_problem_id,i);
  mesh_id = apr_get_mesh_id(field_id);
  
  while ((node_id = mmr_get_next_node_all(mesh_id, node_id)) != 0) {
    if (apr_get_ent_pdeg(field_id, APC_VERTEX, node_id) > 0) {
      int numdofs = apr_get_ent_nrdofs(field_id, APC_VERTEX, node_id);
      apr_read_ent_dofs(field_id, APC_VERTEX, node_id, numdofs, 
			Current, sol_dofs_current);
      apr_read_ent_dofs(field_id, APC_VERTEX, node_id, numdofs, 
			Old, sol_dofs_old);
      i = 0; 
      for (i = 0; i < numdofs; ++i) {

/*kbw
  double norm=0.0;
  printf("node_id %d\ti %d\tcurrent %lf\told %lf\tnorm %lf\n",
	 node_id,i,sol_dofs_current[i],sol_dofs_old[i],norm);
  norm += (sol_dofs_current[i] - sol_dofs_old[i]) * (sol_dofs_current[i] - sol_dofs_old[i]);
/*kew*/

	temp = fabs(sol_dofs_current[i] - sol_dofs_old[i]);
	if (norm_conv_diff < temp)
	  norm_conv_diff = temp;
      }
    }
  }
  
  
  
#ifdef PARALLEL
  pcr_allreduce_max_double(1, &norm_conv_diff, &temp);
  norm_conv_diff = temp;
#endif

  *sol_diff_norm_conv_diff_p = norm_conv_diff;
  
  return;
}


/*------------------------------------------------------------
pdr_conv_diff_time_integration - time integration driver. 
------------------------------------------------------------*/  
void pdr_conv_diff_time_integration(
  int Problem_id,       /* in: problem data structure to be used */
  char* Work_dir,
  FILE *Interactive_input, 
  FILE *Interactive_output
)
{
/* time integration parameters */
  int time_discr, cur_step, final_step, time_monitor;
  double cur_time, cur_dtime;
  double sol_norm_un_conv_diff;

/* linear solver parameters */
  int solver_id;
  int control, comp_type, monitor, nr_levels, ilev, ini_guess;
  char solver_filename[100];

/* auxiliary variables */ 
  int my_proc_id, problem_id, field_id, mesh_id, iaux, slope;
  int i, iadapt, interval, adapt;
  double daux, faux;

/* constants */
  int ione   =  1;
  int itwo   =  2;
  int ithree =  3;
  int izero  =  0;

  /*  conv_diff problem parameters */
  pdt_conv_diff_problem *problem_conv_diff = 
    &pdv_conv_diff_problems[pdv_conv_diff_current_problem_id];
  pdt_conv_diff_ctrls *ctrl_conv_diff = &problem_conv_diff->ctrl;
  pdt_conv_diff_times *time_conv_diff = &problem_conv_diff->time;
  pdt_conv_diff_nonls *nonl_conv_diff = &problem_conv_diff->nonl;
  pdt_conv_diff_linss *lins_conv_diff = &problem_conv_diff->lins;
  pdt_conv_diff_adpts *adpt_conv_diff = &problem_conv_diff->adpt;


/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* get suitable parameters for first problem */
  i=1; time_discr=pdr_time_i_params(PDC_USE_CURRENT_PROBLEM_ID, i);
  
  fprintf(Interactive_output,"time_scheme %d\n",time_discr);
  
  /* in case of time integration */ 
  if(time_discr==1){
    
    problem_id=pdv_conv_diff_current_problem_id;
    i=2; mesh_id=pdr_ctrl_i_params(problem_id,i);
    i=3; field_id=pdr_ctrl_i_params(problem_id,i);
    
    i=11; slope=pdr_ctrl_i_params(problem_id,i);
      
    /* get adaptation parameters */
    i=1; adapt = pdr_adapt_i_params(problem_id,i);
    i=2; interval = pdr_adapt_i_params(problem_id,i);
    
    /* set parameter indicating new mesh */
    iadapt=1;
    
    i=2; cur_step=pdr_time_i_params(problem_id,i);
    i=3; final_step=pdr_time_i_params(problem_id,i);
    i=10; time_monitor=pdr_time_i_params(problem_id,i);
    
    if(final_step<=cur_step){
      
      if(Interactive_input == stdin){
	
#ifdef PARALLEL
	if (pcr_my_proc_id() == pcr_print_master()) {
#endif
	  printf("How many steps to perform?\n");
	  scanf("%d",&iaux);
#ifdef PARALLEL
	}
	pcr_bcast_int(pcr_print_master(), 1, &iaux);
#endif
	
      } else{
	
	// read control data from interactive_input file
	fscanf(Interactive_input,"%d\n",&iaux);
	
      }
      
      final_step += iaux;
      i=3; pdr_set_time_i_params(problem_id,i,final_step);
      
    } 
    
    
    
    if(pdr_lins_i_params(problem_id, 1) >= 0){ // mkb interface for solvers

      int solver_type = pdr_lins_i_params(problem_id, 1);
      int max_iter = -1;
      int error_type = -1;
      double error_tolerance = -1;
      int monitoring_level = -1;

      // when no parameter file passed - take control parameters from problem input file
      if(0==strlen(ctrl_conv_diff->solver_filename)){

	strcpy(solver_filename, ctrl_conv_diff->solver_filename);
	max_iter = pdr_lins_i_params(problem_id, 2); // max_iter
	error_type = pdr_lins_i_params(problem_id, 3); // error_type
	error_tolerance = pdr_lins_d_params(problem_id, 4); // error_tolerance
	monitoring_level = pdr_lins_i_params(problem_id, 5);  // monitoring level

      }
      else{

	/* solver file */
	sprintf(solver_filename, "%s/%s", 
		ctrl_conv_diff->work_dir, ctrl_conv_diff->solver_filename); 

      }

      int parallel = SIC_SEQUENTIAL;
#ifdef PARALLEL
      parallel = SIC_PARALLEL;
#endif

/*kbw
    fprintf(interactive_output, "initializing conv_diff solver (file %s)\n", solver_filename);
    fprintf(interactive_output, "parameters: parallel %d, maxiter %d, error_type %d, error_meas %.15lf, monitor %d\n", 
	    parallel, max_iter,  error_type,  error_tolerance, monitoring_level);
/*kbw*/

      solver_id = sir_init(solver_type, parallel, pdr_get_max_num_grid_levels(problem_id), 
			   solver_filename, max_iter,
			   error_type, error_tolerance, monitoring_level);
      ctrl_conv_diff->solver_id = solver_id;
      fprintf(Interactive_output, "Assigned solver ID: %d for CONV_DIFF\n", 
	      solver_id);

   }
    
    
    /* set parameter indicating new mesh */
    iadapt=1;
    
    /* start loop over time steps */
    while (utv_SIGINT_not_caught) {
      cur_step++;
      i=2; pdr_set_time_i_params(problem_id,i,cur_step);
      i=4; cur_time=pdr_time_d_params(problem_id,i);
      i=6; cur_dtime=pdr_time_d_params(problem_id,i);
      cur_time+=cur_dtime;
      i=4; pdr_set_time_d_params(problem_id,i,cur_time);
      
      if(time_monitor){
	fprintf(Interactive_output,"\nSolving time step %d (<=%d), dtime = %lf, time = %lf\n",
		cur_step,final_step,cur_dtime,cur_time);
      }
      
      apr_rewr_sol(field_id, 1,2 );

/*kbw
	printf("After rewriting:\n");
	{
	  int nel, mesh_id, num_shap, el_mate;
	  double dofs_loc[APC_MAXELSD];

	nel=0;
	while((nel=mmr_get_next_act_elem(mesh_id,nel))!=0){
	  if(nel==110){
	  el_mate =  mmr_el_groupID(mesh_id, nel);
	  printf("Element %d, material %d\n",nel, el_mate);
	  num_shap = apr_get_ent_numshap(field_id, APC_ELEMENT, nel);
	  apr_read_ent_dofs(field_id, APC_ELEMENT, nel, num_shap, 1, dofs_loc);
	  printf("DOFS_1: ");
	  for(i=0;i<num_shap;i++) printf("%15.12lf",dofs_loc[i]);
	  printf("\n");
	  apr_read_ent_dofs(field_id, APC_ELEMENT, nel, num_shap, 2, dofs_loc);
	  printf("DOFS_2: ");
	  for(i=0;i<num_shap;i++) printf("%15.12lf",dofs_loc[i]);
	  printf("\n");
	  }
	}

	}
/*kew*/
      
      
      /* create data structures for stiffness matrix and RHS vector */
      if(iadapt==1) {
	
#ifdef PARALLEL
	/* initiate exchange tables for DOFs - for one field, max_num_levels */
	// one DOF - offset = 0, nreq = 1;
	int nr_levels = pdr_get_max_num_grid_levels(problem_id);
	pdv_exchange_table_index = appr_init_exchange_tables(pcr_nr_proc(), 
							     pcr_my_proc_id(), 
							     field_id, 0, 1, 
							     nr_levels);
#endif
	if(pdr_lins_i_params(problem_id, 1) >= 0){ // mkb interface for solvers
	  sir_create(solver_id, problem_id);
	}
	iadapt = 0;
	
	
      }
      
      if(pdr_lins_i_params(problem_id, 1) >= 0){ // mkb interface for solvers
	
	int nr_iter = -1;
	double conv_meas = -1.0;
	double conv_rate = -1.0;
	int monitor = -1;
	
	// when no parameter file passed - take control parameters from problem input file
	if(0==strlen(ctrl_conv_diff->solver_filename)){
	  nr_iter = pdr_lins_i_params(problem_id, 2); // max_iter
	  conv_meas = pdr_lins_d_params(problem_id, 4); // error_tolerance
	  monitor = pdr_lins_i_params(problem_id, 5);  // monitoring level
	}

	int ini_guess = 1; // get initial guess from data structure

	/*---------- CALLING DIRECT OR ITERATIVE SOLVER THROUGH MKB INTERFACE ---------*/
	sir_solve(solver_id, SIC_SOLVE, ini_guess, monitor, 
		  &nr_iter, &conv_meas, &conv_rate);
	if(pdr_lins_i_params(problem_id, 1) > 0){ // iterative solver
	  fprintf(Interactive_output, 
		  "\nAfter %d iterations of linear solver for conv_diff problem\n", 
		  nr_iter); 
	  fprintf(Interactive_output, 
		  "Convergence measure: %lf, convergence rate %lf\n", 
		  conv_meas, conv_rate); 
	}
	
      } else {
	
#ifdef PARALLEL
	fprintf(Interactive_output,
		"Do not use direct solver for parallel execution!\n");
	exit(-1); 
#endif
	
#ifndef PARALLEL
	sir_direct_solve_lin_sys(problem_id, SIC_SEQUENTIAL, solver_filename);	
#endif
	
      }
      
      
      /* post-process solution using slope limiter */
      if(slope) {
	iaux=pdr_slope_limit(problem_id); 
	if(iaux<0) { printf("\nError in slope!\n");getchar();}
      }
      
      
/*kbw
	printf("Before adapt and balance\n");
	if(mmr_el_status(mesh_id,3039)==MMC_INACTIVE){
	  int elem_struct[100];
	  int face_struct[100];
	  int num_faces=5;
	  int ied,ned,ifa,nfa,num_edges;
	  int fath = mmr_el_fam(mesh_id, 3039, NULL, NULL);
	  int edge_struct[100];
	  mmr_elem_structure(mesh_id, fath, elem_struct);
	  printf("elem 3039: fath %d, grandpa %d, refi %d\n",
		 fath, elem_struct[3], elem_struct[4]);
	  for(ifa=0;ifa<num_faces;ifa++) {	    
	    nfa = abs(elem_struct[5+ifa]);
	    printf("\tface%d, sign %d, id %d, ",
		   ifa, nfa/elem_struct[5+ifa], nfa);
	    mmr_face_structure(mesh_id, nfa, face_struct);
	    printf("type %d, bc %d\n",
		   face_struct[0], face_struct[2]);
	    if(abs(face_struct[0])==MMC_QUAD) num_edges=4;
	    else if(abs(face_struct[0])==MMC_TRIA) num_edges=3;
	    for(ied=0;ied<num_edges;ied++){
	      ned=abs(face_struct[3+ied]);
	      printf("\t\tedge%d sign %d, id %d\n",
		     ied, ned/abs(face_struct[3+ied]), ned);
	      mmr_edge_structure(mesh_id, ned, edge_struct);
	      printf("\t\ttype %d, n/s_1 %d, n/s_2 %d, IPID %d\n",
		     edge_struct[0],edge_struct[1],
		     edge_struct[2],edge_struct[3]);
	    }
	  }
	  printf("clustering 3039\n");
	  apr_derefine(mesh_id, 3039);
	  apr_refine(mesh_id, 3039);
	}
	else if(mmr_el_status(mesh_id,3039)==MMC_ACTIVE){
	  int elem_struct[100];
	  int face_struct[100];
	  int num_faces=5;
	  int ied,ned,ifa,nfa,num_edges;
	  int fath = mmr_el_fam(mesh_id, 3039, NULL, NULL);
	  int edge_struct[100];
	  mmr_elem_structure(mesh_id, fath, elem_struct);
	  printf("elem 3039: fath %d, grandpa %d, refi %d\n",
		 fath, elem_struct[3], elem_struct[4]);
	  for(ifa=0;ifa<num_faces;ifa++) {	    
	    nfa = abs(elem_struct[5+ifa]);
	    printf("\tface%d, sign %d, id %d, ",
		   ifa, nfa/elem_struct[5+ifa], nfa);
	    mmr_face_structure(mesh_id, nfa, face_struct);
	    printf("type %d, bc %d\n",
		   face_struct[0], face_struct[2]);
	    if(abs(face_struct[0])==MMC_QUAD) num_edges=4;
	    else if(abs(face_struct[0])==MMC_TRIA) num_edges=3;
	    for(ied=0;ied<num_edges;ied++){
	      ned=abs(face_struct[3+ied]);
	      printf("\t\tedge%d sign %d, id %d\n",
		     ied, ned/abs(face_struct[3+ied]), ned);
	      mmr_edge_structure(mesh_id, ned, edge_struct);
	      printf("\t\ttype %d, n/s_1 %d, n/s_2 %d, IPID %d\n",
		     edge_struct[0],edge_struct[1],
		     edge_struct[2],edge_struct[3]);
	    }
	  }
	}
	else{
	  printf("FREE SPACE 3039!!!\n");
	}
/*kew*/

      /* graphics data dumps */
      if (time_conv_diff->intv_graph > 0 && 
	  time_conv_diff->cur_step % time_conv_diff->intv_graph == 0) {
	
	fprintf(Interactive_output, "(Writing field to disk (graphics)...)\n");
	
	/* sprintf(autodump_filename, "%s/%s_g%d.dmp", ctrl_conv_diff->work_dir, ctrl_conv_diff->field_dmp_filepattern, time_conv_diff->cur_step); */
	/* if (apr_write_field(field_id, PDC_NREQ, 1, time_conv_diff->intv_graph_accu, autodump_filename) < 0) */
	/* 	fprintf(Interactive_output, "Error in writing field data!\n"); */
	
	//pdr_conv_diff_write_paraview(Work_dir, 
	//			       Interactive_input, Interactive_output);
      }
      
      /* full data dumps (for restarting) */
      if (time_conv_diff->intv_dumpout > 0 
	  && time_conv_diff->cur_step % time_conv_diff->intv_dumpout == 0) {
	
	fprintf(Interactive_output, "\n(Writing field to disk - not implemented)\n");
	
	/* sprintf(autodump_filename, "%s/%s_f.dmp", ctrl_conv_diff->work_dir, ctrl_conv_diff->field_dmp_filepattern); */
	/* if (apr_write_field(field_id, PDC_NREQ, 0, 0, autodump_filename) < 0) */
	/* 	fprintf(Interactive_output, "Error in writing field data!\n"); */
	
      }
      
      /* check for stop conditions */
      
      /* stop when final time reached */
      if (time_conv_diff->cur_time >= time_conv_diff->final_time) {
	fprintf(Interactive_output, "\nFinal time reached. Stopping.\n");
	break;
      }    
      
      /* cur_dtime could change through simulation (time step adaptation) */
      /* thus we need to check also if cur_step not bigger than final_step */
      if (time_conv_diff->cur_step >= time_conv_diff->final_step) {
	fprintf(Interactive_output, "\nFinal step reached. Stopping.\n");
	break;
      }
      
      /* stop if convergence reached (for stationary problems) */
      if (sol_norm_un_conv_diff < time_conv_diff->conv_meas) {
	fprintf(Interactive_output, 
		"\nNorm conv_diff (u_n, prev. u_n) below n_epsilon: %lf. Stopping.\n", 
		time_conv_diff->conv_meas);
	break;
      }
      
      // when time for adaptation
      if( adpt_conv_diff->type>0 && adpt_conv_diff->interval>0 &&
	  (time_conv_diff->cur_step+1)%adpt_conv_diff->interval==0 ) {
	
#ifdef PARALLEL
	/* free exchange tables for DOFs */
	appr_free_exchange_tables(pdv_exchange_table_index);
	pdv_exchange_table_index = -1;
#endif
	
	if(pdr_lins_i_params(problem_id, 1) >= 0){ // mkb interface for solvers
	  /* free solver data structures */
	  sir_free(solver_id);
	}
	pdr_conv_diff_adapt(Work_dir,
			    Interactive_input, Interactive_output);
	/* indicate to recreate block structure */
	iadapt=1;
      }
      

/*kbw
	printf("After adapt\n");
	if(mmr_el_status(mesh_id,3039)==MMC_INACTIVE){
	  int elem_struct[100];
	  int face_struct[100];
	  int num_faces=5;
	  int ied,ned,ifa,nfa,num_edges;
	  int fath = mmr_el_fam(mesh_id, 3039, NULL, NULL);
	  int edge_struct[100];
	  mmr_elem_structure(mesh_id, fath, elem_struct);
	  printf("elem 3039: fath %d, grandpa %d, refi %d\n",
		 fath, elem_struct[3], elem_struct[4]);
	  for(ifa=0;ifa<num_faces;ifa++) {	    
	    nfa = abs(elem_struct[5+ifa]);
	    printf("\tface%d, sign %d, id %d, ",
		   ifa, nfa/elem_struct[5+ifa], nfa);
	    mmr_face_structure(mesh_id, nfa, face_struct);
	    printf("type %d, bc %d\n",
		   face_struct[0], face_struct[2]);
	    if(abs(face_struct[0])==MMC_QUAD) num_edges=4;
	    else if(abs(face_struct[0])==MMC_TRIA) num_edges=3;
	    for(ied=0;ied<num_edges;ied++){
	      ned=abs(face_struct[3+ied]);
	      printf("\t\tedge%d sign %d, id %d\n",
		     ied, ned/abs(face_struct[3+ied]), ned);
	      mmr_edge_structure(mesh_id, ned, edge_struct);
	      printf("\t\ttype %d, n/s_1 %d, n/s_2 %d, IPID %d\n",
		     edge_struct[0],edge_struct[1],
		     edge_struct[2],edge_struct[3]);
	    }
	  }
	  printf("clustering 3039\n");
	  apr_derefine(mesh_id, 3039);
	  apr_refine(mesh_id, 3039);
	}
	else if(mmr_el_status(mesh_id,3039)==MMC_ACTIVE){
	  int elem_struct[100];
	  int face_struct[100];
	  int num_faces=5;
	  int ied,ned,ifa,nfa,num_edges;
	  int fath = mmr_el_fam(mesh_id, 3039, NULL, NULL);
	  int edge_struct[100];
	  mmr_elem_structure(mesh_id, fath, elem_struct);
	  printf("elem 3039: fath %d, grandpa %d, refi %d\n",
		 fath, elem_struct[3], elem_struct[4]);
	  for(ifa=0;ifa<num_faces;ifa++) {	    
	    nfa = abs(elem_struct[5+ifa]);
	    printf("\tface%d, sign %d, id %d, ",
		   ifa, nfa/elem_struct[5+ifa], nfa);
	    mmr_face_structure(mesh_id, nfa, face_struct);
	    printf("type %d, bc %d\n",
		   face_struct[0], face_struct[2]);
	    if(abs(face_struct[0])==MMC_QUAD) num_edges=4;
	    else if(abs(face_struct[0])==MMC_TRIA) num_edges=3;
	    for(ied=0;ied<num_edges;ied++){
	      ned=abs(face_struct[3+ied]);
	      printf("\t\tedge%d sign %d, id %d\n",
		     ied, ned/abs(face_struct[3+ied]), ned);
	      mmr_edge_structure(mesh_id, ned, edge_struct);
	      printf("\t\ttype %d, n/s_1 %d, n/s_2 %d, IPID %d\n",
		     edge_struct[0],edge_struct[1],
		     edge_struct[2],edge_struct[3]);
	    }
	  }
	}
	else{
	  printf("FREE SPACE 3039!!!\n");
	}
/*kew*/
      

/*kbw
	printf("After adapt and balance\n");
	if(mmr_el_status(mesh_id,3039)==MMC_INACTIVE){
	  int elem_struct[100];
	  int face_struct[100];
	  int num_faces=5;
	  int ied,ned,ifa,nfa,num_edges;
	  int fath = mmr_el_fam(mesh_id, 3039, NULL, NULL);
	  int edge_struct[100];
	  mmr_elem_structure(mesh_id, fath, elem_struct);
	  printf("elem 3039: fath %d, grandpa %d, refi %d\n",
		 fath, elem_struct[3], elem_struct[4]);
	  for(ifa=0;ifa<num_faces;ifa++) {	    
	    nfa = abs(elem_struct[5+ifa]);
	    printf("\tface%d, sign %d, id %d, ",
		   ifa, nfa/elem_struct[5+ifa], nfa);
	    mmr_face_structure(mesh_id, nfa, face_struct);
	    printf("type %d, bc %d\n",
		   face_struct[0], face_struct[2]);
	    if(abs(face_struct[0])==MMC_QUAD) num_edges=4;
	    else if(abs(face_struct[0])==MMC_TRIA) num_edges=3;
	    for(ied=0;ied<num_edges;ied++){
	      ned=abs(face_struct[3+ied]);
	      printf("\t\tedge%d sign %d, id %d\n",
		     ied, ned/abs(face_struct[3+ied]), ned);
	      mmr_edge_structure(mesh_id, ned, edge_struct);
	      printf("\t\ttype %d, n/s_1 %d, n/s_2 %d, IPID %d\n",
		     edge_struct[0],edge_struct[1],
		     edge_struct[2],edge_struct[3]);
	    }
	  }
	  printf("clustering 3039\n");
	  apr_derefine(mesh_id, 3039);
	  apr_refine(mesh_id, 3039);
	}
	else if(mmr_el_status(mesh_id,3039)==MMC_ACTIVE){
	  int elem_struct[100];
	  int face_struct[100];
	  int num_faces=5;
	  int ied,ned,ifa,nfa,num_edges;
	  int fath = mmr_el_fam(mesh_id, 3039, NULL, NULL);
	  int edge_struct[100];
	  mmr_elem_structure(mesh_id, fath, elem_struct);
	  printf("elem 3039: fath %d, grandpa %d, refi %d\n",
		 fath, elem_struct[3], elem_struct[4]);
	  for(ifa=0;ifa<num_faces;ifa++) {	    
	    nfa = abs(elem_struct[5+ifa]);
	    printf("\tface%d, sign %d, id %d, ",
		   ifa, nfa/elem_struct[5+ifa], nfa);
	    mmr_face_structure(mesh_id, nfa, face_struct);
	    printf("type %d, bc %d\n",
		   face_struct[0], face_struct[2]);
	    if(abs(face_struct[0])==MMC_QUAD) num_edges=4;
	    else if(abs(face_struct[0])==MMC_TRIA) num_edges=3;
	    for(ied=0;ied<num_edges;ied++){
	      ned=abs(face_struct[3+ied]);
	      printf("\t\tedge%d sign %d, id %d\n",
		     ied, ned/abs(face_struct[3+ied]), ned);
	      mmr_edge_structure(mesh_id, ned, edge_struct);
	      printf("\t\ttype %d, n/s_1 %d, n/s_2 %d, IPID %d\n",
		     edge_struct[0],edge_struct[1],
		     edge_struct[2],edge_struct[3]);
	    }
	  }
	}
	else{
	  printf("FREE SPACE 3039!!!\n");
	}
/*kew*/

      
    }  //end loop over timesteps
    
    if (iadapt == 0) {

#ifdef PARALLEL
      /* free exchange tables for DOFs */
      appr_free_exchange_tables(pdv_exchange_table_index);
      pdv_exchange_table_index = -1;
#endif

      if(pdr_lins_i_params(problem_id, 1) >= 0){ // mkb interface for solvers
	sir_free(solver_id);
      }
      
    }

    if(pdr_lins_i_params(problem_id, 1) >= 0){ // mkb interface for solvers
      sir_destroy(solver_id);
    }

  }
  
  return;
}
