/************************************************************************
File pds_conv_diff_main - main control routine for interactive solution of
                          convection-diffusion equations 

Contains definitions of routines:   

  main() - just for C
  pdr_conv_diff_init - initialize problem data

 
REMARK:
The code uses solution_1 for character reading problem in C - namely
scanf("%c", &var);getchar();
Possible solution_2 is:
scanf(" %c", &var);

------------------------------  			
History:        
	08.2008 - Krzysztof Banas, pobanas@cyf-kr.edu.pl, initial version	
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>
#include<string.h>


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

/* graphics module */
#include "mod_fem_viewer.h"

#include "uth_log.h"
#include "uth_io_results.h"
#include "uth_system.h"

/**************************************/
/* GLOBAL CONSTANTS                   */
/**************************************/
/* Rules:
/* - constants always uppercase and start with PDC_ */


/**************************************/
/* GLOBAL VARIABLES                   */
/**************************************/
/* Rules:
/* - name always begins with pdv_conv_diff_ */
int  pdv_conv_diff_nr_problems;   /* the number of problems in the simulation */
int  pdv_conv_diff_current_problem_id;     /* ID of the current problem */
pdt_conv_diff_problem pdv_conv_diff_problems[PDC_CONV_DIFF_MAX_NUM_PROB];      
                                                         /* array of problems */
#ifdef PARALLEL
int pdv_exchange_table_index = -1;
#endif

// important generic data structure for several phases in calculations
// in multi-physics problems patches may be shared by several problems
// or there may be several patches structures 
// data structure is made global, similarly to mesh 
// it is not included in mesh, since it may also depend on approximation
utt_patches *pdv_patches=NULL;
// other procedures may create their own structures of the same type

/* global variable */
char work_dir[200]; /* the name of the working directory */

const int Current_solution_ID       = 1;

/***************************************/
/* DECLARATIONS OF INTERNAL PROCEDURES */
/***************************************/
/* Rules:
/* - name always begins with pdr_conv_diff_ */
/* - argument names start uppercase */

/* to test setting initial approximation field */
double Initial_condition(int Field_id, double* Coor, int Sol_comp_id);

/*---------------------------------------------------------
  pdr_conv_diff_init - to initialize problem dependent data
----------------------------------------------------------*/
extern int pdr_conv_diff_init(
  char* Work_dir,
  FILE *Interactive_input,
  FILE *Interactive_output
  );


/***************************************/
/* DEFINITIONS OF PROCEDURES */
/***************************************/
/*------------------------------------------------------------
main function
------------------------------------------------------------*/
int main(int argc, char **argv)
{

/* local variables */
  char interactive_input_name[300];
  char interactive_output_name[300], tmp[100];
  char field_module_name[100];
  FILE *interactive_input, *interactive_output;
  char c, d, arg[300]; /* string variable to read menu */
  int i, iaux, jaux, kaux;
  int nreq, nr_sol, iel, nr_mat, base;
  int problem_id, field_id, mesh_id;
  double daux, t_wall;
  int nrno, nmno, nred, nmed, nrfa, nmfa, nrel, nmel;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /****************************************/
  /* input/output settings                */
  /****************************************/

  /* very simple - for the beginning instead of GUI */
  if(argv[1]==NULL){
    strcpy(work_dir,".");
  } else {
    sprintf(work_dir,"%s",argv[1]);
  }

  utr_set_interactive(work_dir, argc, argv,
		      &interactive_input, &interactive_output);

  /****************************************/
  /* input/output settings - done         */
  /****************************************/

#ifdef DEBUG
  fprintf(interactive_output,"\nStarting program in debug mode.\n");
#endif
#ifdef DEBUG_MMM
  fprintf(interactive_output,"Starting mesh module in debug mode.\n");
#endif
#ifdef DEBUG_APM
  fprintf(interactive_output,"Starting approximation module in debug mode.\n");
#endif
#ifdef DEBUG_SIM
  fprintf(interactive_output,"Starting solver interface in debug mode.\n");
#endif
#ifdef DEBUG_LSM
  fprintf(interactive_output,"Starting linear solver (adapter) module in debug mode.\n");
#endif

  /******************************************************************/
  /* initialization of problem data (including mesh and field) */
  /******************************************************************/

  /*++++ initialization of problem, mesh and field data ++++*/
  pdr_conv_diff_init(work_dir, interactive_input, interactive_output);

  problem_id = pdv_conv_diff_current_problem_id;
  i=2; mesh_id=pdr_ctrl_i_params(problem_id,i);
  i=3; field_id = pdr_ctrl_i_params(problem_id,i);  
  i=4; nr_sol = pdr_ctrl_i_params(problem_id,i);
  i=5; nreq = pdr_ctrl_i_params(problem_id,i);
  // currently supported fields:
  // STANDARD_LINEAR - for continuous, vector, linear approximations
  // DG_SCALAR_PRISM - for discontinuous, scalar, high order approximations
  apr_module_introduce(field_module_name);

  /****************************************/
  /* main menu loop                       */
  /****************************************/
  do {

      fflush(interactive_output);

#ifdef PARALLEL
    if(pcr_my_proc_id()==pcr_print_master()){
#endif

     if(interactive_input == stdin){
    
      do {
	/* define a menu */
	printf("\nChoose a command from the menu:\n");
	printf("\tt - integrate in time \n");
	printf("\ts - solve single convection-diffusion problem \n");
	printf("\tu - solve single convection-diffusion problem using direct solver \n");
	printf("\tp - postprocessing\n");
	printf("\tk - launch graphic module\n");
	printf("\te - compute error for test problem\n");
	printf("\tz - compute ZZ error estimator for test problem\n");
	printf("\tx - compute better ZZ error estimator for test problem\n");
	printf("\ta - automatic mesh adaptation based on ZZ error estimate\n");
	printf("\tm - perform manual mesh adaptation \n");
	printf("\td - dump out data \n");
	if( strncmp(field_module_name, "STANDARD_LINEAR", 15) == 0){
	printf("\tv - dump out ParaView graphics data \n");
	}
	printf("\tg - dump out ModFemViewer graphics data \n");
	printf("\tc - change control data \n");
	printf("\tr - perform random test of h-adaptivity\n");
	printf("\tq - exit the program\n");
	
	scanf(" %c",&c); 
      } while ( c != 'm' && c != 'd' && c != 'z'  && c != 'a' 
		&& c != 'g' && c != 'r' && c != 'v' && c != 'u'
		&& c != 's' && c != 'e' && c != 'p' && c != 'x'
		&& c != 't' && c != 'c' && c != 'q' && c != 'k');
      
     } else{

      // read control data from interactive_input file
      fscanf(interactive_input,"%c\n",&c);

     }

#ifdef PARALLEL
    }

    pcr_bcast_char(pcr_print_master(),1,&c);
    //printf("After BCAST %c\n",c);

#endif

    // dumping data
    if(c=='d'){
      
/* different mesh files can be taken into account using MESHFIL_TYP parameter */ 
      sprintf(arg,"%s/mesh_prism.dmp", work_dir);
      iaux = utr_io_export_mesh(interactive_output,mesh_id, MMC_MOD_FEM_MESH_DATA, arg);
      if(iaux<0) {
	/* if error in writing data */
	fprintf(interactive_output,"Error in writing mesh data!\n");
      }

      nreq = pdr_ctrl_i_params(problem_id, 5);
      nr_sol =  pdr_ctrl_i_params(problem_id, 4);

      if(nr_sol==1) iaux=1; // parameter to select vectors written to file
      if(nr_sol==2) iaux=3; // both vectors written to file
      // scheme: 1 - 1, 2 - 2, 3 - 4, 1+2 - 3, 1+3 - 5, 2+3 - 6, 1+2+3 - 7

      sprintf(arg,"%s/field.dmp", work_dir);


      // write field to file arg with accuracy parameter = 0 - full accuracy 
      iaux = utr_export_field(field_id, nreq, iaux, 0, arg);
      if(iaux<0) {
	/* if error in writing data */
	fprintf(interactive_output,"Error in writing field data!\n");
      }

      sprintf(arg,"%s/problem_conv_diff.dmp", work_dir);
      iaux = pdr_conv_diff_write_problem(problem_id,  arg); 
      if(iaux<0) {
	/* if error in writing data */
	fprintf(interactive_output,"Error in writing problem data!\n");
      }
      
    }

    // storing ParaView graphics data
    else if(c=='v'){

      if( strncmp(field_module_name, "STANDARD_LINEAR", 15) == 0){  
#ifdef PARALLEL
      sprintf(arg, "%s/dump_paraview_proc%d.vtk", work_dir,pcr_my_proc_id());
#else
    sprintf(arg, "dump_paraview");
#endif

    utr_write_paraview_std_lin(mesh_id, field_id, work_dir,arg, UTE_VTK_XML);
      }
      else{
	fprintf(interactive_output,"ParaView dump only for std_lin approx!\n");
      }

    }
    else if(c=='k'){
      	c='\0';
#ifdef MOD_FEM_VIEWER
      	init_mod_fem_viewer(argc,argv,interactive_output);
#else
	fprintf(interactive_output,"The ModFEM Viewer module is unavailable\n");
#endif
    }

    // solving the problem with direct solver
    else if(c=='u'){
      
#ifndef PARALLEL

      /* very simple time measurement */
      t_wall=time_clock();

      /* problem dependent interface with linear solvers */
      sprintf(arg,"%s/pardiso.dat", work_dir);

      sir_direct_solve_lin_sys(problem_id, SIC_SEQUENTIAL, arg);

      /* very simple time measurement */
      t_wall=time_clock()-t_wall;
      
      fprintf(interactive_output,
	      "\nTime for solving a system of linear equations %lf\n\n",t_wall);
#endif

    }

    // solving the problem
    else if(c=='s'){

        utr_io_result_write_values_and_proceed();

        nrno = mmr_get_nr_node(pdr_select_problem(pdv_conv_diff_current_problem_id)->ctrl.mesh_id);
        utr_io_result_add_value_int(RESULT_N_DOFS, nrno);

      
#ifdef PARALLEL
    if(pcr_my_proc_id()==pcr_print_master()){
#endif
      fprintf(interactive_output,
"\nBeginning solution of a single convection-diffusion problem\n\n");
#ifdef PARALLEL
    }
#endif

      /* very simple time measurement */
      t_wall=time_clock();
      
#ifdef PARALLEL
      {
/* initiate exchange tables for DOFs - for one field,  max_num_levels */
// one DOF - offset = 0, nreq = 1;
	int nr_levels = pdr_get_max_num_grid_levels(problem_id);
	pdv_exchange_table_index = appr_init_exchange_tables(pcr_nr_proc(), 
							     pcr_my_proc_id(), 
							     field_id, 0, 1, 
							     nr_levels); 
      }
#endif
      /* problem dependent interface with linear solvers */
      sprintf(arg,"%s/mkb.dat", work_dir);

#ifndef PARALLEL
      pdt_conv_diff_problem* problem = pdr_select_problem(problem_id);
      sir_solve_lin_sys(problem_id, problem->lins.type, SIC_SEQUENTIAL,
			pdr_get_max_num_grid_levels(problem_id), arg, 
			// below: defaults
			-1, //pdr_lins_i_params(problem_id, 2), // max_iter
			-1, //pdr_lins_i_params(problem_id, 3), // error_type
			-1.0, //pdr_lins_d_params(problem_id, 4), // error_tol
			-1 //pdr_lins_i_params(problem_id, 5) // monitor level
			);
#endif
      
#ifdef PARALLEL
      sir_solve_lin_sys(problem_id, 1, SIC_PARALLEL,
			pdr_get_max_num_grid_levels(problem_id), arg, 
			// below: defaults
			-1, //pdr_lins_i_params(problem_id, 2), // max_iter
			-1, //pdr_lins_i_params(problem_id, 3), // error_type
			-1.0, //pdr_lins_d_params(problem_id, 4), // error_tol
			-1 //pdr_lins_i_params(problem_id, 5) // monitor level
			);
/* free exchange tables for DOFs - for one field = one solver */
      pdr_error_test_exact(field_id,interactive_output);

      appr_free_exchange_tables(pdv_exchange_table_index);
#endif

      /* very simple time measurement */
      t_wall=time_clock()-t_wall;
       utr_io_result_add_value_double(RESULT_TIME_TO_SOL, t_wall);
#ifdef PARALLEL
    if(pcr_my_proc_id()==pcr_print_master()){
#endif
      fprintf(interactive_output,
	      //printf(
"\nTime for solving the problem %lf\n\n",t_wall);


#ifdef PARALLEL
    }
#endif
     utr_ctrl_pts_add_values(field_id);
     utr_io_result_write_values_and_proceed();
      
    }

    // computing error
    else if(c=='e'){
#ifndef PARALLEL
      double error_h1;
      /* simple L2 and H1 error for the test problem (with known exact solution) */
      error_h1 = pdr_error_test_exact(field_id,interactive_output);
#else
        mf_log_err("Unable to do so, because of appr_free_exchange_tables(pdv_exchange_table_index) called earlier!");
#endif

    }
    else if(c=='z'){
      
      double error_ZZ;
      /* Zienkiewicz-Zhu error estimator */
      error_ZZ = pdr_zzhu_error_small(field_id,interactive_output);
      
    }
    else if(c=='x'){
      
      double error_ZZ;
      /* Zienkiewicz-Zhu error estimator */
      error_ZZ = pdr_zzhu_error(field_id,interactive_output);
      
    }

    // performing simple post-processing
    else if(c=='p'){

      pdr_post_process(field_id,interactive_input,interactive_output);

    }

    else if(c=='a'){
      
      /* adaptation based on Zienkiewicz-Zhu error estimator */
      pdr_conv_diff_adapt(problem_id, work_dir, 
		interactive_input, interactive_output );
      
    }

    // manual mesh adaptation
    else if(c=='m'){
      
      utr_manual_refinement( problem_id, work_dir, 
			     interactive_input, interactive_output );

    }

    // testing refinement/derefinement procedures
    else if(c=='r'){

      utr_test_refinements(problem_id, work_dir, 
			   interactive_input, interactive_output );

    } /* end test of mesh refinements */
  
    // performing time integration
    else if(c=='t'){

#ifdef PARALLEL
    if(pcr_my_proc_id()==pcr_print_master()){
#endif
      fprintf(interactive_output, "\nBeginning time integration\n\n");
#ifdef PARALLEL
    }
#endif

      /* very simple time measurement */
      t_wall=time_clock();
      
      /* problem dependent interface with linear solvers */
      pdr_conv_diff_time_integration(problem_id, work_dir, 
			   interactive_input, interactive_output); 

      /* very simple time measurement */
      t_wall=time_clock()-t_wall;

#ifdef PARALLEL
    if(pcr_my_proc_id()==pcr_print_master()){
#endif
      fprintf(interactive_output,
	      //printf(
"\nTotal time of integration in time %lf\n\n",t_wall);
#ifdef PARALLEL
    }
#endif
      
    }

    // changing control data
    else if(c=='c'){
      
      pdr_change_data(problem_id);
      
    }

  } while(c != 'q');
 
  /* free allocated space */
  apr_free_field(field_id);
  mmr_free_mesh(mesh_id); 

#ifdef PARALLEL
  if(pcv_my_proc_id == PCC_MASTER_PROC_ID) {
    utr_io_result_clear_rows();
    utr_io_result_gather_and_avg_all_files(pcv_nr_proc,pcv_my_proc_id,3,0); // setup
    utr_io_result_gather_and_avg_all_files(pcv_nr_proc,pcv_my_proc_id,4,1); // setup
  }
  pcr_exit_parallel();
#endif

  fclose(interactive_input);
  fclose(interactive_output);

  return(0);
  
}


/* to test setting initial approximation field */
double Initial_condition(int Field_id, double* Coor, int Sol_comp_id)
{

  return(Coor[0]+Coor[1]+Coor[2]);
}




/*---------------------------------------------------------
  pdr_conv_diff_init_problem - to initialize problem (including mesh and field) data
----------------------------------------------------------*/
int
pdr_conv_diff_init(
  char* Work_dir,
  FILE *interactive_input,
  FILE *interactive_output
  )
{

  /* pointer to problem structure */
  pdt_conv_diff_problem *problem;

/* auxiliary variables */
  FILE *problem_input;
  char c, field_type, keyword[300], arg[100];
  char problem_input_name[300];
  int num_prob, iprob, meshfil_typ, num_bc, ibc, mesh_id, field_id, pdeg;
  int nrno, nmno, nred, nmed, nrfa, nmfa, nrel, nmel;
  int nreq, nr_sol, imat, base;
  int max_gen, max_gen_diff;
  int i, j;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /*+++++++++++   opening file with problem data   ++++++++++++*/
  sprintf(problem_input_name,"%s/problem_conv_diff.dat", Work_dir);
  problem_input = fopen(problem_input_name,"r");
  
  if(problem_input == NULL){
    fprintf(interactive_output,
	    "\nCannot find problem input file %s. \n", problem_input_name);
    exit(-1); 
  }
  
  /* read control data */
  fscanf(problem_input,"%s\n",keyword);
  if(strncmp(keyword,"NUMBER_OF_PROBLEMS", 12) != 0) {
    printf("ERROR in reading input file for NUMBER_OF_PROBLEMS\n");
    exit(1);
  }
  fscanf(problem_input,"%d",&num_prob);
  utr_skip_rest_of_line(problem_input);
  
  pdv_conv_diff_nr_problems=num_prob;


  char csv_name[128]={0};
  strcpy(csv_name, pdr_select_problem(pdv_conv_diff_current_problem_id)->ctrl.name);
  if(strlen(csv_name) == 0 ) {
      strcpy(csv_name, "conv_diff");
  }
#ifdef PARALLEL
  char proc_id_string[8]={0};
  sprintf(proc_id_string,"%d",pcv_my_proc_id);
  strcat(csv_name,proc_id_string);
#endif
  strcat(csv_name,".csv");
  utr_io_result_set_filename(work_dir,csv_name);
  utr_io_result_clear_columns();
  utr_io_result_add_column(RESULT_STEP);
  utr_io_result_add_column(RESULT_N_DOFS);
  utr_io_result_add_column(RESULT_TIME_TO_SOL);
  utr_io_result_add_column(RESULT_TIME_PCM);
  utr_io_result_add_column(RESULT_TIME_DDM);
  utr_io_result_add_column(RESULT_TIME_MMPM);
  utr_io_result_add_column(RESULT_TIME_APPM);
  utr_io_result_add_column(RESULT_TIME_SIM);
  utr_io_result_add_column(RESULT_OUT_FILENAME);
  utr_ctrl_pts_init(Work_dir,"control_points.dat",1,"control_potins.csv");
  utr_io_result_write_columns(0);
  
  for(iprob=0;iprob<num_prob;iprob++){
    
    pdv_conv_diff_current_problem_id = iprob+1;
    problem = &pdv_conv_diff_problems[iprob];
    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"PROBLEM_TYPE", 10) != 0) {
      printf("ERROR in reading input file for PROBLEM_TYPE\n");
      exit(1);
    }
    fscanf(problem_input,"%s",&problem->ctrl.name );
    utr_skip_rest_of_line(problem_input);
    



    /*+++++++++++   initialization of mesh data   ++++++++++++*/
    
    // mesh/meshfile types:
    //   j - gradmesh 2D mesh
    //   p - prismatic mesh in standard (text) format
    //   t - tetrahedral mesh in standard format
    //   h - hybrid mesh in standard format
    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"MESH_TYPE",8) != 0) {
      printf("ERROR in reading input file for MESH_TYPE\n");
      exit(-1);
    }
    fscanf(problem_input,"%c",&c);
    utr_skip_rest_of_line(problem_input);
    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"MESH_FILE",8) != 0) {
      printf("ERROR in reading input file for MESH_FILE\n");
      exit(-1);
    }
    fscanf(problem_input,"%s",&arg);
    utr_skip_rest_of_line(problem_input);
    
    mesh_id = utr_initialize_mesh(interactive_output, work_dir, c, arg);
    problem->ctrl.mesh_id = mesh_id;
    
    
    /*+++++++++++   initialization of approximation field data   ++++++++++++*/
    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"FIELD_SOURCE",8) != 0) {
      printf("ERROR in reading input file for FIELD_SOURCE\n");
      exit(-1);
    }
    
    // options for field initialization
    // r - read from file
    // i - initialize to specified function values
    // z - initialize to zero
    fscanf(problem_input,"%c",&c);
    utr_skip_rest_of_line(problem_input);
    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"NREQ",4) != 0) {
      printf("ERROR in reading input file for NREQ\n");
      exit(-1);
    }
    fscanf(problem_input,"%d", &problem->ctrl.nreq );
    utr_skip_rest_of_line(problem_input);
    if(problem->ctrl.nreq > PDC_CONV_DIFF_NREQ){
      printf("Parameter NREQ read from input file > PDC_CONV_DIFF_NREQ\n");
      printf("Type any key (twice) to continue or CTRL C to quit\n"); 
      getchar();getchar();
    }
    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"NRSOL",4) != 0) {
      printf("ERROR in reading input file for NRSOL\n");
      exit(-1);
    }
    fscanf(problem_input,"%d", &problem->ctrl.nr_sol );
    utr_skip_rest_of_line(problem_input);
    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"BASIS_FUNCTIONS",4) != 0) {
      fprintf(interactive_output,"ERROR in reading input file for BASIS_FUNCTIONS\n");
      exit(1);
    }
    fscanf(problem_input,"%c",&field_type);
    if(field_type=='c'){
      // set of discontinuous complete basis functions
      problem->ctrl.base=APC_BASE_COMPLETE_DG;
    }
    else if(field_type=='t'){
      // set of discontinuous tensor product basis functions
      problem->ctrl.base=APC_BASE_TENSOR_DG;
    }
    else{
      // for continuous basis functions problem->ctrl.base is not used
      problem->ctrl.base=-1;
    }

    utr_skip_rest_of_line(problem_input);

    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"PDEG",4) != 0) {
      fprintf(interactive_output,"ERROR in reading input file for PDEG\n");
      exit(-1);
    }
    fscanf(problem_input,"%d", &pdeg );
    utr_skip_rest_of_line(problem_input);

    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"FIELD_FILE",8) != 0) {
      printf("ERROR in reading input file for FIELD_FILE\n");
      exit(-1);
    }
    fscanf(problem_input,"%s",&arg);
    utr_skip_rest_of_line(problem_input);
    

    /* reading field values saved by previous runs */
    if(c=='r'){
      
      field_id = utr_initialize_field(interactive_output,field_type,c,mesh_id, 
        problem->ctrl.nreq, problem->ctrl.nr_sol, pdeg, arg, NULL);
      
    }
    /* initializing field values according to the specified function */
    else if(c=='i'){
      
      field_id = utr_initialize_field(interactive_output,field_type,c,mesh_id, 
	problem->ctrl.nreq,problem->ctrl.nr_sol,pdeg,NULL,&Initial_condition);
      
    }
    /* initializing field values to zero */
    else if(c=='z'){
      
      field_id = utr_initialize_field(interactive_output,field_type,c,mesh_id, 
      problem->ctrl.nreq, problem->ctrl.nr_sol, pdeg, NULL, NULL);
      
    }
    
    problem->ctrl.field_id = field_id;
    
    nreq = pdr_ctrl_i_params(pdv_conv_diff_current_problem_id, 5);
    nr_sol =  pdr_ctrl_i_params(pdv_conv_diff_current_problem_id, 6);
    
    
    
    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"NUMBER_OF_BOUNDARY_CONDITIONS", 25) != 0) {
      printf("ERROR in reading input file for NUMBER_OF_BOUNDARY_CONDITIONS\n");
      exit(1);
    }
    fscanf(problem_input,"%d",&num_bc);
    utr_skip_rest_of_line(problem_input);
    if(num_bc>0) {
      fscanf(problem_input,"%s\n",keyword);
      if(strncmp(keyword,"BC_VALUES", 8) != 0) {
	printf("ERROR in reading input file for BC_VALUES\n");
	exit(1);
      }
      problem->ctrl.bc_val=(double **)malloc(3*PDC_CONV_DIFF_MAX_BC_VAL*sizeof(double *));
      for(i=0;i<3*PDC_CONV_DIFF_MAX_BC_VAL;i++) problem->ctrl.bc_val[i]=NULL;
      for(i=0;i<num_bc;i++){
	fscanf(problem_input,"%d",&ibc);
	utr_skip_rest_of_line(problem_input);
	problem->ctrl.bc_val[ibc] =
	  utr_dvector(problem->ctrl.nreq,"bc_val in read");
	for(j=0;j<problem->ctrl.nreq;j++){
	  fscanf(problem_input,"%lg",&problem->ctrl.bc_val[ibc][j]);
	}
	utr_skip_rest_of_line(problem_input);
      }
    }
    else problem->ctrl.bc_val=NULL;
    
    if(field_type!='s'){
      fscanf(problem_input,"%s\n",keyword);
      if(strncmp(keyword,"PENALTY_CONSTANT",12) != 0) {
	fprintf(interactive_output,"ERROR in reading input file for PENA\n");
	exit(1);
      }
      fscanf(problem_input,"%lg",&problem->ctrl.penalty);
      utr_skip_rest_of_line(problem_input);

      fscanf(problem_input,"%s\n",keyword);
      if(strncmp(keyword,"SLOPE_LIMITER_SWITCH",25) != 0) {
	fprintf(interactive_output,"ERROR in reading input file for SLOP\n");
	exit(1);
    }
      fscanf(problem_input,"%d",&problem->ctrl.slope);
      utr_skip_rest_of_line(problem_input);
    }    

    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"TIME_INTEGRATION_PARAMETERS",25) != 0) {
      printf("ERROR in reading input file for TIME_INTEGRATION_PARAMETERS\n");
      exit(1);
    }
    fscanf(problem_input,"%d %d", &problem->time.type, &problem->time.monitor);
    utr_skip_rest_of_line(problem_input);
    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"IMPLICITNESS_PARAMETER",15) != 0) {
      printf("ERROR in reading input file for IMPL\n");
      exit(1);
    }
    fscanf(problem_input,"%lg",&problem->time.alpha);
    utr_skip_rest_of_line(problem_input);
    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"TIME_STEP_CURRENT_AND_FINAL",25) != 0) {
      printf("ERROR in reading input file for TIME_STEP\n");
      exit(1);
    }
    fscanf(problem_input,"%d %d",  &problem->time.cur_step ,
	   &problem->time.final_step );
    utr_skip_rest_of_line(problem_input);
    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"TIME_INSTANT_CURRENT_AND_FINAL",25) != 0) {
      printf("ERROR in reading input file for TIME\n");
      exit(1);
    }
    fscanf(problem_input,"%lg %lg", 
	   &problem->time.cur_time, &problem->time.final_time );
    utr_skip_rest_of_line(problem_input);
    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"DELTA_TIME_CURRENT_AND_PREVIOUS",25) != 0) {
      printf("ERROR in reading input file for DELTA_TIME\n");
      exit(1);
    }
    fscanf(problem_input,"%lg %lg", &problem->time.cur_dtime, 
	   &problem->time.prev_dtime );
    utr_skip_rest_of_line(problem_input);
    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"CONVERGENCE_TYPE_AND_LIMIT",25) != 0) {
      printf("ERROR in reading input file for CONVERGENCE\n");
      exit(1);
    }
    fscanf(problem_input,"%d %lg", &problem->time.conv_type,
	   &problem->time.conv_meas );
    utr_skip_rest_of_line(problem_input);
    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"INTERVAL_FOR_DUMPS_STANDARD_AND_GRAPHICS",25) != 0) {
      printf("ERROR in reading input file for INTERVAL_FOR_DUMPS\n");
      exit(1);
    }
    fscanf(problem_input,"%d %d", &problem->time.intv_dumpout,
	   &problem->time.intv_graph );
    utr_skip_rest_of_line(problem_input);
    
    /* fscanf(problem_input,"%s\n",keyword); */
    /* if(strncmp(keyword,"NONL_PARAMS") != 0) { */
    /*   printf("ERROR in reading input file for NONL_PARAMS\n"); */
    /*   exit(1); */
    /* } */
    /* fscanf(problem_input,"%d %d", &problem->nonl.type,  */
    /* 	&problem->nonl.monitor); */
    /* utr_skip_rest_of_line(problem_input); */
    
    /* fscanf(problem_input,"%s\n",keyword); */
    /* if(strncmp(keyword,"CONV") != 0) { */
    /*   printf("ERROR in reading input file for CONV\n"); */
    /*   exit(1); */
    /* } */
    /* fscanf(problem_input,"%d %d %lg", &problem->nonl.max_iter,  */
    /*     &problem->nonl.conv_type, &problem->nonl.conv_meas); */
    /* utr_skip_rest_of_line(problem_input); */
    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"LINEAR_SOLVER_PARAMETERS",20) != 0) {
      printf("ERROR in reading input file for LINEAR_SOLVER_PARAMETERS\n");
      exit(1);
    }
    fscanf(problem_input,"%d %d", &problem->lins.type,
	   &problem->lins.monitor  ); 
    utr_skip_rest_of_line(problem_input);
    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"SOLVER_CONVERGENCE_PARAMETERS",15) != 0) {
      printf("ERROR in reading input file for CONVERGENCE_PARAMETERS\n");
      exit(1);
    }
    fscanf(problem_input,"%d %d %lg", 
	   &problem->lins.max_iter , &problem->lins.conv_type ,
	   &problem->lins.conv_meas  ); 
    utr_skip_rest_of_line(problem_input);
    
    
    fscanf(problem_input,"%s\n",keyword);
    if(strncmp(keyword,"ADAPTATION_PARAMETERS",10) != 0) {
      printf("ERROR in reading input file for ADAPTATION__PARAMETERS\n");
      exit(1);
    }
    fscanf(problem_input,"%d %d",&problem->adpt.type,
	   &problem->adpt.monitor );
    utr_skip_rest_of_line(problem_input);
    
    if(problem->adpt.type==0){ /* adaptivity switched off */
      mmr_set_max_gen(mesh_id, 10);
      mmr_set_max_gen_diff(mesh_id, 1);
    }
    else{
      
      fscanf(problem_input,"%s\n",keyword);
      if(strncmp(keyword,"ADAPT_INTERVAL_MAXGEN_ETC",15) != 0) {
	printf("ERROR in reading input file for ADAPT_INTERVAL_MAXGEN\n");
	exit(1);
      }
      fscanf(problem_input,"%d %d %d",
	     &problem->adpt.interval, &problem->adpt.maxgen,
	     &problem->adpt.maxgendiff);
      utr_skip_rest_of_line(problem_input);
      
      // maximal generation for elements
      mmr_set_max_gen(mesh_id, problem->adpt.maxgen); 
      
      // irregularity constraints switched on/off 
      // (maxgendiff==0 - no constraints) 
      // (for constrained meshes only 1-irregular meshes allowed)
      if(problem->adpt.maxgendiff != 0) mmr_set_max_gen_diff(mesh_id, 1);
      else mmr_set_max_gen_diff(mesh_id, 0); 
      
      fscanf(problem_input,"%s\n",keyword);
      if(strncmp(keyword,"ADAPT_TOLERANCE_REF_UNREF",25) != 0) {
	printf("ERROR in reading input file for ADAPT_TOLERANCE_REF_UNREF\n");
	exit(1);
      }
      fscanf(problem_input,"%lg %lg",
	     &problem->adpt.eps, &problem->adpt.ratio);
      utr_skip_rest_of_line(problem_input);
      
      
      fprintf(interactive_output,"\nmaximal generation level set to %d, maximal generation difference set to %d\n", 
	      problem->adpt.maxgen, problem->adpt.maxgendiff);
      
    }
  }

  apr_check_field(field_id);

  //i=1; nrno=mmr_mesh_i_params(mesh_id,i);
  nrno = mmr_get_nr_node(mesh_id);
  //i=2; nmno=mmr_mesh_i_params(mesh_id,i);
  nmno = mmr_get_max_node_id(mesh_id);
  //i=5; nred=mmr_mesh_i_params(mesh_id,i);
  nred = mmr_get_nr_edge(mesh_id);
  //i=6; nmed=mmr_mesh_i_params(mesh_id,i);
  nmed = mmr_get_max_edge_id(mesh_id);
  //i=9; nrfa=mmr_mesh_i_params(mesh_id,i);
  nrfa = mmr_get_nr_face(mesh_id);
  //i=10; nmfa=mmr_mesh_i_params(mesh_id,i);
  nmfa = mmr_get_max_face_id(mesh_id);
  //i=13; nrel=mmr_mesh_i_params(mesh_id,i);
  nrel = mmr_get_nr_elem(mesh_id);
  //i=14; nmel=mmr_mesh_i_params(mesh_id,i);
  nmel = mmr_get_max_elem_id(mesh_id);
  //i=21; max_gen=mmr_get_mesh_i_params(mesh,i);
  max_gen=mmr_get_max_gen(mesh_id);
  //i=22; max_gen_diff=mmr_get_mesh_i_params(mesh,i);
  max_gen_diff=mmr_get_max_gen_diff(mesh_id);

  fprintf(interactive_output,"\nAfter reading initial data.\n");
  fprintf(interactive_output,"Parameters (number of active, maximal index):\n");
  fprintf(interactive_output,"Elements: nrel %d, nmel %d\n", nrel, nmel);
  fprintf(interactive_output,"Faces:    nrfa %d, nmfa %d\n", nrfa, nmfa); 
  fprintf(interactive_output,"Edges:    nred %d, nmed %d\n", nred, nmed); 
  fprintf(interactive_output,"Nodes:    nrno %d, nmno %d\n", nrno, nmno); 

  fclose(problem_input);

  return(1);
}


/*---------------------------------------------------------
  pdr_conv_diff_write_problem - to write problem dependent data to a disk file
----------------------------------------------------------*/
int pdr_conv_diff_write_problem(
  int Problem_id, 
  char* Filename
  )
{
  /* pointer to problem structure */
  pdt_conv_diff_problem *problem;

/* auxiliary variables */
  FILE *fp;
  int num_prob, iprob, num_bc, ibc, base;
  int i, j;

  int meshfil_typ=MMC_MOD_FEM_MESH_DATA;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* open the input file with control variables */
  fp=fopen(Filename,"w");
  if(fp==NULL) {
    printf("Cannot open file %s to write problem dependent data\n",Filename);
    return(-1);
  } 

/* write control data */
  fprintf(fp,"%s\n","NUMBER_OF_PROBLEMS");
  fprintf(fp,"%d\n",pdv_conv_diff_nr_problems);

  for(iprob=0;iprob<pdv_conv_diff_nr_problems;iprob++){

    problem = &pdv_conv_diff_problems[iprob];

    fprintf(fp,"%s\n","PROBLEM_ID");
    fprintf(fp,"%s\n",problem->ctrl.name );

    fprintf(fp,"%s\n","MESH_TYPE");
    fprintf(fp,"%c\n",'p' );

    fprintf(fp,"%s\n","MESH_FILE");
    fprintf(fp,"%s\n","./mesh_prism.dmp" );

    fprintf(fp,"%s\n","FIELD_SOURCE");
    fprintf(fp,"%c\n",'r' );

    fprintf(fp,"%s\n","NREQ");
    fprintf(fp,"%d\n",problem->ctrl.nreq);

    fprintf(fp,"%s\n","NRSOL");
    fprintf(fp,"%d\n",problem->ctrl.nr_sol);

    base =  pdr_ctrl_i_params(Problem_id, 13);
    fprintf(fp,"%s\n","BASIS_FUNCTIONS");
    if(base==APC_BASE_TENSOR_DG){
      fprintf(fp,"%c\n",'t');
    }
    else if(base==APC_BASE_COMPLETE_DG){
      fprintf(fp,"%c\n",'c');
    }
    else{
      fprintf(fp,"%c\n",'s');
    }

    fprintf(fp,"%s\n","FIELD_FILE");
    if(base==APC_BASE_TENSOR_DG || APC_BASE_COMPLETE_DG){
      fprintf(fp,"%s\n","./field_dg.dmp" );
    }
    else{
      fprintf(fp,"%s\n","./field_std.dmp" );
    }


    fprintf(fp,"%s\n","NUMBER_OF_BOUNDARY_CONDITIONS");
    if(problem->ctrl.bc_val==NULL){
      i=0; fprintf(fp,"%d\n", i);
    }
    else{
      num_bc=0;
      for(i=0;i<3*PDC_CONV_DIFF_MAX_BC_VAL;i++){
	if(problem->ctrl.bc_val[i]!=NULL) num_bc++;
      }
      fprintf(fp,"%d\n", num_bc);
      if(num_bc>0) {
	fprintf(fp,"%s\n","BC_VALUES");
	for(i=0;i<3*PDC_CONV_DIFF_MAX_BC_VAL;i++){
	  if(problem->ctrl.bc_val[i]!=NULL) {
	    fprintf(fp,"%d\n", i);
	    for(j=0;j<problem->ctrl.nreq;j++){
	      fprintf(fp,"%lg",problem->ctrl.bc_val[i][j]);
	    }
	    fprintf(fp,"\n");
	  }
	}
      }
    }

    if(base==APC_BASE_TENSOR_DG || APC_BASE_COMPLETE_DG){

      fprintf(fp,"%s\n","PENALTY_CONSTANT");
      fprintf(fp,"%lg\n",problem->ctrl.penalty);
      
      fprintf(fp,"%s\n","SLOPE_LIMITER_SWITCH");
      fprintf(fp,"%d\n",problem->ctrl.slope);

    }

    /* fprintf(fp,"%s\n","ICON"); */
    /* if(problem->ctrl.ic_val==NULL){ */
    /*   i=0; fprintf(fp,"%d\n", i); */
    /* } */
    /* else{ */
    /*   num_bc=0; */
    /*   for(i=0;;i++){ */
    /* 	if(problem->ctrl.ic_val[i]!=NULL) num_bc++; */
    /* 	else break; */
    /*   } */
    /*   fprintf(fp,"%d\n", num_bc); */
    /*   if(num_bc>0) { */
    /* 	for(i=0;i<num_bc;i++){ */
    /* 	  if(problem->ctrl.ic_val[i]!=NULL) { */
    /* 	    for(j=0;j<problem->ctrl.nreq;j++){ */
    /* 	      fprintf(fp,"%lg",problem->ctrl.ic_val[i][j]); */
    /* 	    } */
    /* 	    fprintf(fp,"\n"); */
    /* 	  } */
    /* 	} */
    /*   } */
    /* } */


    fprintf(fp,"%s\n","TIME_INTEGRATION_PARAMETERS");
    fprintf(fp,"%d %d\n", problem->time.type, problem->time.monitor);

    fprintf(fp,"%s\n","IMPLICITNESS_PARAMETER");
    fprintf(fp,"%lg\n",problem->time.alpha);

    fprintf(fp,"%s\n","TIME_STEP_CURRENT_AND_FINAL");
    fprintf(fp,"%d %d\n",  problem->time.cur_step ,
	problem->time.final_step );

    fprintf(fp,"%s\n","TIME_INSTANT_CURRENT_AND_FINAL");
    fprintf(fp,"%lg %lg\n", 
	problem->time.cur_time, problem->time.final_time );

    fprintf(fp,"%s\n","DELTA_TIME_CURRENT_AND_PREVIOUS");
    fprintf(fp,"%lg %lg\n", problem->time.cur_dtime, 
	problem->time.prev_dtime );

    fprintf(fp,"%s\n","CONVERGENCE_TYPE_AND_LIMIT");
    fprintf(fp,"%d %lg\n", problem->time.conv_type,
	problem->time.conv_meas );

    fprintf(fp,"%s\n","INTERVAL_FOR_DUMPS_STANDARD_AND_GRAPHICS");
    fprintf(fp,"%d %d\n", problem->time.intv_dumpout,
	problem->time.intv_graph );

    /* fprintf(fp,"%s\n","NONL_PARAMS"); */
    /* fprintf(fp,"%d %d\n", problem->nonl.type,  */
    /* 	problem->nonl.monitor ); */
    /* fprintf(fp,"%s\n","CONV"); */
    /* fprintf(fp,"%d %d %lg\n",  */
    /* 	problem->nonl.max_iter, problem->nonl.conv_type, */
    /* 	problem->nonl.conv_meas ); */

    fprintf(fp,"%s\n","LINEAR_SOLVER_PARAMETERS");
    fprintf(fp,"%d %d\n", problem->lins.type,
	problem->lins.monitor ); 
    fprintf(fp,"%s\n","SOLVER_CONVERGENCE_PARAMETERS");
    fprintf(fp,"%d %d %lg\n", 
	problem->lins.max_iter , problem->lins.conv_type ,
	problem->lins.conv_meas ); 

    fprintf(fp,"%s\n","ADAPTATION_PARAMETERS");
    fprintf(fp,"%d %d\n",problem->adpt.type,
	problem->adpt.monitor );
    fprintf(fp,"%s\n","ADAPT_INTERVAL_MAXGEN_ETC");
    fprintf(fp,"%d %d %d\n",
	problem->adpt.interval, problem->adpt.maxgen,
	problem->adpt.maxgendiff );
    fprintf(fp,"%s\n","ADAPT_TOLERANCE_REF_UNREF");
    fprintf(fp,"%lg %lg\n",
	problem->adpt.eps, problem->adpt.ratio );

  }

  fclose(fp);

  return(1);
}



