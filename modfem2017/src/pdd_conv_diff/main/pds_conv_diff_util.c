/************************************************************************
File pds_conv_diff_util - utility routines for interactive solution of
                          convection-diffusion equations 

Contains definitions of routines:   

  pdr_module_introduce - to return the problem module's name
  pdr_problem_name - to return the problem's name
  pdr_select_problem - to select the proper problem   
  pdr_ctrl_i_params - to return one of control parameters
  pdr_ctrl_d_params - to return one of control parameters
  pdr_time_i_params - to return parameters of time discretization
  pdr_time_d_params - to return parameters of time discretization
  pdr_set_time_i_params - to change parameters of time discretization
  pdr_set_time_d_params - to change parameters of time discretization
  pdr_lins_i_params - to return parameters of linear equations solver
  pdr_lins_d_params - to return parameters of linear equations solver
  pdr_adapt_i_params - to return parameters of adaptation
  pdr_adapt_d_params - to return parameters of adaptation
  pdr_change_data - to change some of control data 

  pdr_error_test_exact - to compute error norm for test examples
  pdr_post_process - to privide simple interactive post-processing
  pdr_average_sol_el - to compute the average of solution over element
  pdr_exact_sol - to return values and derivatives at a point
	for functions used as exact solutions for test problems


 
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

#include "uth_log.h"

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
#include "../include/pdh_conv_diff_problem.h"	


/***************************************/
/* DEFINITIONS OF PROCEDURES */
/***************************************/
#define SMALL      1.0e-10

/*------------------------------------------------------------
  pdr_module_introduce - to return the problem module's name
------------------------------------------------------------*/
int pdr_module_introduce(
                  /* returns: >=0 - success code, <0 - error code */
  char* Problem_module_name /* out: the name of the problem module */
  )
{

  char* string = "CONV_DIFF";

  strcpy(Problem_module_name,string);

  return(1);
}

/*------------------------------------------------------------
  pdr_problem_name - to return the problem's name
------------------------------------------------------------*/
int pdr_problem_name(
                  /* returns: >=0 - success code, <0 - error code */
  int Problem_id,  /* in: problem ID or 
		         PDC_USE_CURRENT_PROBLEM_ID for the current problem */
  char* Problem_name /* out: the name of the problem solved */
  )
{

/* auxiliary variables */
  pdt_conv_diff_problem* problem;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  problem = pdr_select_problem(Problem_id);

  strcpy(Problem_name, problem->ctrl.name);

  return(1);
}


/*---------------------------------------------------------
  pdr_select_problem - to select the proper problem   
---------------------------------------------------------*/
pdt_conv_diff_problem* pdr_select_problem( /* returns pointer to the chosen problem */
			         /* (if input is not valid it returns */
		                 /* the pointer to the current problem) */
  int Problem_id  /* in: problem ID or 
		         PDC_USE_CURRENT_PROBLEM_ID for the current problem */
  )
{

  /* select the proper problem from the array of problems */
  if( Problem_id == PDC_USE_CURRENT_PROBLEM_ID ) {
    return(&pdv_conv_diff_problems[pdv_conv_diff_current_problem_id-1]);
  }
  else if( Problem_id>0 && Problem_id<=pdv_conv_diff_nr_problems ) {
    return(&pdv_conv_diff_problems[Problem_id-1]);
  }
  else {
    return(&pdv_conv_diff_problems[pdv_conv_diff_current_problem_id-1]);
    /* alternative:  return(NULL);   */
  }
}

/*---------------------------------------------------------
pdr_ctrl_i_params - to return one of control parameters
---------------------------------------------------------*/
int pdr_ctrl_i_params( /* returns: integer problem parameter */
	int Problem_id,	/* in: problem ID  */
	int Num         /* in: parameter number in control structure */
	)
{
/* auxiliary variables */
  pdt_conv_diff_problem* problem;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  problem = pdr_select_problem(Problem_id);

  //if(Num==1) return(problem->ctrl.name); // changed for string
  if(Num==2) return(problem->ctrl.mesh_id);
  else if(Num==3) return(problem->ctrl.field_id);
  else if(Num==4) return(problem->ctrl.nr_sol);
  else if(Num==5) return(problem->ctrl.nreq);
  else if(Num==6) return(problem->ctrl.solver_id); 
  // else if(Num==10) return(problem->ctrl.solved_problem_type); removed
  else if(Num==11) return(problem->ctrl.slope);
  else if(Num==13) return(problem->ctrl.base);
  else {
      mf_log_err("Wrong parameter number(%d) in ctrl_i_params!",Num);
  }

/* error condition - that point should not be reached */
  return(-1);
}

/*---------------------------------------------------------
pdr_ctrl_d_params - to return one of control parameters
---------------------------------------------------------*/
double pdr_ctrl_d_params( /* returns: real problem parameter */
	int Problem_id,	/* in: problem ID  */
	int Num         /* in: parameter number in control structure */
	)
{
/* auxiliary variables */
  pdt_conv_diff_problem* problem;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  problem = pdr_select_problem(Problem_id);

  if(Num==12) return(problem->ctrl.penalty);
  //if(Num==8) return(problem->ctrl.coeff1);
  //else if(Num==9) return(problem->ctrl.coeff2);
  //else if(Num==10) return(problem->ctrl.coeff3);
  else 
  {
    printf("Wrong parameter number in ctrl_d_params!");
    exit(1);
  }

/* error condition - that point should not be reached */
  return(-1);
}


/*---------------------------------------------------------
pdr_time_i_params - to return parameters of time discretization
---------------------------------------------------------*/
int pdr_time_i_params( /* returns: integer time integration parameter */
	int Problem_id,	/* in: data structure to be used  */
	int Num         /* in: parameter number in control structure */
	)
{
/* auxiliary variables */
  pdt_conv_diff_problem* problem;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  problem = pdr_select_problem(Problem_id);

  if(Num==1) return(problem->time.type);
  else if(Num==2) return(problem->time.cur_step);
  else if(Num==3) return(problem->time.final_step);
  else if(Num==8) return(problem->time.conv_type);
  else if(Num==10) return(problem->time.monitor);
  else if(Num==11) return(problem->time.intv_dumpout);
  else if(Num==12) return(problem->time.intv_graph);
  else {
    printf("Wrong parameter number in time_i_params!");
    exit(1);
  }

/* error condition - that point should not be reached */
  return(-1);
}

/*---------------------------------------------------------
pdr_time_d_params - to return parameters of time discretization
---------------------------------------------------------*/
double pdr_time_d_params( /* returns: real time integration parameter */
	int Problem_id,	/* in: data structure to be used  */
	int Num         /* in: parameter number in control structure */
	)
{
/* auxiliary variables */
  pdt_conv_diff_problem* problem;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  problem = pdr_select_problem(Problem_id);

  if (Num == 2) return (problem->time.alpha);
  else if(Num==4) return(problem->time.cur_time);
  else if(Num==5) return(problem->time.final_time);
  else if(Num==6) return(problem->time.cur_dtime);
  else if(Num==7) return(problem->time.prev_dtime);
  else if(Num==9) return(problem->time.conv_meas);
  else {
    printf("Wrong parameter number in time_d_params!");
    exit(1);
  }

/* error condition - that point should not be reached */
  return(-1);
}

/*---------------------------------------------------------
pdr_set_time_i_params - to change parameters of time discretization
---------------------------------------------------------*/
void pdr_set_time_i_params( 
        int Problem_id,	     /* in: data structure to be used  */
	int Num,             /* in: parameter number in control structure */
	int Value            /* in: parameter value */
	)
{
/* auxiliary variables */
  pdt_conv_diff_problem* problem;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  problem = pdr_select_problem(Problem_id);

  if(Num==1) problem->time.type=Value;
  else if(Num==2) problem->time.cur_step=Value;
  else if(Num==3) problem->time.final_step=Value;
  else if(Num==8) problem->time.conv_type=Value;
  else if(Num==10) problem->time.monitor=Value;
  else if(Num==11) problem->time.intv_dumpout=Value;
  else if(Num==12) problem->time.intv_graph=Value;
  else {
    printf("Wrong parameter number in set_time_i_params!");
    exit(1);
  }

  return;
}

/*---------------------------------------------------------
pdr_set_time_d_params - to change parameters of time discretization
---------------------------------------------------------*/
void pdr_set_time_d_params( 
        int Problem_id,	     /* in: data structure to be used  */
	int Num,             /* in: parameter number in control structure */
	double Value         /* in: parameter value */
	)
{
/* auxiliary variables */
  pdt_conv_diff_problem* problem;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  problem = pdr_select_problem(Problem_id);

  if(Num==2) problem->time.alpha=Value;
  else if(Num==4) problem->time.cur_time=Value;
  else if(Num==5) problem->time.final_time=Value;
  else if(Num==6) problem->time.cur_dtime=Value;
  else if(Num==7) problem->time.prev_dtime=Value;
  else if(Num==9) problem->time.conv_meas=Value;
  else {
    printf("Wrong parameter number in set_time_d_params!");
    exit(1);
  }

  return;
}


/*---------------------------------------------------------
pdr_lins_i_params - to return parameters of linear equations solver
---------------------------------------------------------*/
int pdr_lins_i_params( /* returns: integer linear solver parameter */
	int Problem_id,	/* in: data structure to be used  */
	int Num         /* in: parameter number in control structure */
	)
{
/* auxiliary variables */
  pdt_conv_diff_problem* problem;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  problem = pdr_select_problem(Problem_id);

  if(Num==1) return(problem->lins.type);
  else if(Num==2) return(problem->lins.max_iter);
  else if(Num==3) return(problem->lins.conv_type);
  else if(Num==5) return(problem->lins.monitor);
  else {
    printf("Wrong parameter number in lins_i_params!");
    exit(-1);
  }

/* error condition - that point should not be reached */
  return(-1);
}

/*---------------------------------------------------------
pdr_lins_d_params - to return parameters of linear equations solver
---------------------------------------------------------*/
double pdr_lins_d_params( /* returns: real linear solver parameter */
	int Problem_id,	/* in: data structure to be used  */
	int Num         /* in: parameter number in control structure */
	)
{
/* auxiliary variables */
  pdt_conv_diff_problem* problem;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  problem = pdr_select_problem(Problem_id);

  if(Num==4) return(problem->lins.conv_meas);
  else {
    printf("Wrong parameter number in lins_d_params!");
    exit(-1);
  }

/* error condition - that point should not be reached */
  return(-1);
}

/*---------------------------------------------------------
pdr_adapt_i_params - to return parameters of adaptation
---------------------------------------------------------*/
int pdr_adapt_i_params( /* returns: integer adaptation parameter */
	int Problem_id,	/* in: data structure to be used  */
	int Num         /* in: parameter number in control structure */
	)
{
/* auxiliary variables */
  pdt_conv_diff_problem* problem;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  problem = pdr_select_problem(Problem_id);

  if(Num==1) return(problem->adpt.type);
  else if(Num==2) return(problem->adpt.interval);
  else if(Num==3) return(problem->adpt.maxgen);
  else if(Num==4) return(problem->adpt.maxgendiff);
  else if(Num==7) return(problem->adpt.monitor);
  else {
    printf("Wrong parameter number in adapt_i_params!");
    exit(1);
  }

/* error condition - that point should not be reached */
  return(-1);
}

/*---------------------------------------------------------
pdr_adapt_d_params - to return parameters of adaptation
---------------------------------------------------------*/
double pdr_adapt_d_params( /* returns: real adaptation parameter */
	int Problem_id,	/* in: data structure to be used  */
	int Num         /* in: parameter number in control structure */
	)
{
/* auxiliary variables */
  pdt_conv_diff_problem* problem;

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  problem = pdr_select_problem(Problem_id);

  if(Num==5) return(problem->adpt.eps);
  else if(Num==6) return(problem->adpt.ratio);
  else {
    printf("Wrong parameter number in adapt_d_params!");
    exit(1);
  }

/* error condition - that point should not be reached */
  return(-1);
}

/*---------------------------------------------------------
pdr_change_data - to change some of control data 
---------------------------------------------------------*/
void pdr_change_data(
	int Problem_id	/* in: data structure to be used  */
        )
{

/* local variables */
  pdt_conv_diff_problem* problem;
  char c, d, pans[100]; /* string variable to read menu */

/*++++++++++++++++ executable statements ++++++++++++++++*/

/* select the proper problem data structure */
  problem = pdr_select_problem(Problem_id);

    do {

      do {
/* define a menu */
        printf("\nChoose a group of data:\n");
        printf("\tc - general control data \n"); 
        printf("\tt - time integration parameters \n"); 
        printf("\ta - adaptation parameters \n"); 
        printf("\tl - linear solver parameters \n"); 
        printf("\tq - quit changing data for problem %d\n",Problem_id);

        scanf("%s",pans);getchar();
      } while ( *pans != 'c' && *pans != 't' && *pans != 'a' 
             && *pans != 'l' && *pans != 'q' && *pans != 'q' );

      c = *pans;

      if(c=='c'){

        do {

          do {
/* define a menu */
            printf("\nChoose variable to change:\n");
            printf("\tl - activation of slope limiter\n"); 
            //printf("\tc - values of three method coefficients\n"); 
            printf("\tq - quit changing general control data\n");

            scanf("%s",pans);getchar();
          } while ( *pans != 's' && *pans != 's' && *pans != 's' 
               && *pans != 'l' && *pans != 'c' && *pans != 'q' );

          d = *pans;

	  if(d=='l'){

            printf("Old value: %d, new value: ",problem->ctrl.slope);
            scanf("%d",&problem->ctrl.slope); getchar();

          }
          else if(d=='c'){

            /* printf("Old values: %lf %lf %lf, new values: ",  */
	    /* 	   problem->ctrl.coeff1,  */
	    /* 	   problem->ctrl.coeff2, problem->ctrl.coeff3); */
            /* scanf("%lg %lg %lg",&problem->ctrl.coeff1,&problem->ctrl.coeff2, */
	    /* 	&problem->ctrl.coeff3); getchar(); */

          }

        } while(d != 'q');

      }
      else if(c=='t'){

        do {

          do {
/* define a menu */
            printf("\nChoose variable to change:\n");
            printf("\ta - method identifier\n"); 
            printf("\tc - current time-step number\n"); 
            printf("\td - final time-step number\n"); 
            printf("\te - current time-step length\n"); 
            printf("\tf - previous time-step length\n"); 
            printf("\tg - current time\n"); 
            printf("\th - final time\n"); 
            printf("\ti - convergence in time criterion number\n"); 
            printf("\tj - convergence in time treshold value\n"); 
            //printf("\tk - values of three method coefficients\n"); 
            printf("\tq - quit changing time integration data\n");

            scanf("%s",pans);getchar();
          } while ( *pans != 'a' && *pans != 'c' && *pans != 'c' 
                 && *pans != 'd' && *pans != 'e' && *pans != 'f' 
                 && *pans != 'g' && *pans != 'h' && *pans != 'i' 
                 && *pans != 'j' && *pans != 'k' && *pans != 'q' );

          d = *pans;

          if(d=='a'){

            printf("Old value: %d, new value: ",problem->time.type);
            scanf("%d",&problem->time.type); getchar();

          }
          else if(d=='c'){

            printf("Old value: %d, new value: ",problem->time.cur_step);
            scanf("%d",&problem->time.cur_step); getchar();

          }
          else if(d=='d'){

            printf("Old value: %d, new value: ",problem->time.final_step);
            scanf("%d",&problem->time.final_step); getchar();

          }
          else if(d=='e'){

            printf("Old value: %lg, new value: ",problem->time.cur_dtime);
            scanf("%lg",&problem->time.cur_dtime); getchar();

          }
          else if(d=='f'){

            printf("Old value: %lg, new value: ",problem->time.prev_dtime);
            scanf("%lg",&problem->time.prev_dtime); getchar();

          }
          else if(d=='g'){

            printf("Old value: %lg, new value: ",problem->time.cur_time);
            scanf("%lg",&problem->time.cur_time); getchar();

          }
          else if(d=='h'){

            printf("Old value: %lg, new value: ",problem->time.final_time);
            scanf("%lg",&problem->time.final_time); getchar();

          }
          else if(d=='i'){

            printf("Old value: %d, new value: ",problem->time.conv_type);
            scanf("%d",&problem->time.conv_type); getchar();

          }
          else if(d=='j'){

            printf("Old value: %lg, new value: ",problem->time.conv_meas);
            scanf("%lg",&problem->time.conv_meas); getchar();

          }
          else if(d=='k'){

            /* printf("Old values: %lg %lg %lg, new values: ",  */
	    /* 	   problem->ctrl.coeff1,  */
	    /* 	   problem->ctrl.coeff2, problem->ctrl.coeff3); */
            /* scanf("%lg %lg %lg",&problem->ctrl.coeff1,&problem->ctrl.coeff2, */
	    /* 	  &problem->ctrl.coeff3); getchar(); */

          }

        } while(d != 'q');

      }
      else if(c=='a'){

        do {

          do {
/* define a menu */
            printf("\nChoose variable to change:\n");
            printf("\tt - strategy number\n"); 
            printf("\ti - time interval between adaptations\n"); 
            printf("\tm - maximal generation level for elements\n"); 
            printf("\td - maximal generation difference for neighboring elements\n"); 
            printf("\te - global treshold value for adaptation\n"); 
            printf("\tr - ratio for indicating derefinements\n"); 
            printf("\tq - quit changing adaptation data\n");

            scanf("%s",pans);getchar();
          } while ( *pans != 't' && *pans != 'i' && *pans != 'm' && *pans != 'd'
               && *pans != 'e' && *pans != 'r' && *pans != 'q' );

          d = *pans;

          if(d=='t'){

            printf("Old value: %d, new value: ",problem->adpt.type);
            scanf("%d",&problem->adpt.type); getchar();

          }
          else if(d=='i'){

            printf("Old value: %d, new value: ",problem->adpt.interval);
            scanf("%d",&problem->adpt.interval); getchar();

          }
          else if(d=='m'){

            printf("Old value: %d, new value: ",problem->adpt.maxgen);
            scanf("%d",&problem->adpt.maxgen); getchar();

          }
          else if(d=='d'){

            printf("Old value: %d, new value: ",problem->adpt.maxgendiff);
            scanf("%d",&problem->adpt.maxgendiff); getchar();

          }
          else if(d=='e'){

            printf("Old value: %lg, new value: ",problem->adpt.eps);
            scanf("%lg",&problem->adpt.eps); getchar();

          }
          else if(d=='r'){

            printf("Old value: %lg, new value: ",problem->adpt.ratio);
            scanf("%lg",&problem->adpt.ratio); getchar();

          }

        } while(d != 'q');

      }

      else if(c=='l'){

        do {

          do {
/* define a menu */
            printf("\nChoose variable to change:\n");
            printf("\tt - solver type\n"); 
            printf("\tm - maximal number of iterations\n"); 
            printf("\tc - convergence criterion number\n"); 
            printf("\te - convergence treshold value\n"); 
            printf("\tq - quit changing linear solver data\n");

            scanf("%s",pans);getchar();
          } while ( *pans != 't' && *pans != 'm' && *pans != 'c' 
                 && *pans != 'e' && *pans != 'p' && *pans != 'k' 
                 && *pans != 'b' && *pans != 'q' && *pans != 'q' );

          d = *pans;

          if(d=='t'){

            printf("Old value: %d, new value: ",problem->lins.type);
            scanf("%d",&problem->lins.type); getchar();

          }
          else if(d=='m'){

            printf("Old value: %d, new value: ",problem->lins.max_iter);
            scanf("%d",&problem->lins.max_iter); getchar();

          }
          else if(d=='c'){

            printf("Old value: %d, new value: ",problem->lins.conv_type);
            scanf("%d",&problem->lins.conv_type); getchar();

          }
          else if(d=='e'){

            printf("Old value: %lg, new value: ",problem->lins.conv_meas);
            scanf("%lg",&problem->lins.conv_meas); getchar();

          }

        } while(d != 'q');

      }

    } while(c != 'q');

    return;
}



/*---------------------------------------------------------
  pdr_error_test_exact - to compute error norm for test examples
----------------------------------------------------------*/
double pdr_error_test_exact( /* returns H1 norm of error for */
	         	/* several known exact solutions */
  int Problem_id,    /* in: problem ID  */
  FILE *Interactive_output
        )
{

/* local variables */
  double l2_err = 0.0, h1_err = 0.0; /* error norms */

  int pdeg;		/* degree of polynomial */
  int pdeg_old=-1; /* indicator for recomputing quadrature data */
  int base;		/* type of basis functions */
  int ngauss;            /* number of gauss points */
  double xg[3000];   	 /* coordinates of gauss points in 3D */
  double wg[1000];       /* gauss weights */

  int num_shap;         /* number of element shape functions */
  int ndofs;            /* local dimension of the problem */
  double determ;        /* determinant of jacobi matrix */
  double vol;           /* volume for integration rule */
  double xcoor[3];      /* global coord of gauss point */
  double u_val[PDC_MAXEQ]; /* computed solution */
  double u_x[PDC_MAXEQ];   /* gradient of computed solution */
  double u_y[PDC_MAXEQ];   /* gradient of computed solution */
  double u_z[PDC_MAXEQ];   /* gradient of computed solution */
  double base_phi[APC_MAXELVD];    /* basis functions */
  double base_dphix[APC_MAXELVD];  /* x-derivatives of basis function */
  double base_dphiy[APC_MAXELVD];  /* y-derivatives of basis function */
  double base_dphiz[APC_MAXELVD];  /* y-derivatives of basis function */
  int el_nodes[MMC_MAXELVNO+1];        /* list of nodes of El */
  double node_coor[3*MMC_MAXELVNO];  /* coord of nodes of El */
  double dofs_loc[APC_MAXELSD]={0}; /* element solution dofs */

/* for scalar test cases - exact solution and its derivatives */
  double f, f_x, f_y, f_z;
/* auxiliary variables */
  int i, iel, ki, iaux, mat_num, mesh_id, field_id, nel, sol_vec_id, nreq, counter;
  double daux=0.0, eaux=0.0, volume=0.0, average=0.0, time=0.0;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  /* select the corresponding mesh */
  i=3; field_id=pdr_ctrl_i_params(Problem_id,i);
  mesh_id = apr_get_mesh_id(field_id);
  time=pdr_time_d_params(Problem_id, 4);

/*kbw
  // checking mesh
  int nr_mine_active=0;
  int nr_mine_inactive=0;
  int nr_overlap_active=0;
  int nr_overlap_inactive=0;
  int nr_free=0;
  int nr_total=0;
  int ient=0;
  for(ient=1;ient<=mmr_get_max_elem_id(mesh_id);ient++){
    if(mmr_el_status(mesh_id,ient)==MMC_FREE) nr_free++;
    else if(mmr_el_status(mesh_id,ient)==MMC_ACTIVE){
#ifdef PARALLEL
      if(mmpr_el_owner(mesh_id,ient)==pcr_my_proc_id()) nr_mine_active++;
      else if(mmpr_el_owner(mesh_id,ient)>0) nr_overlap_active++;
      else {
	printf("Active element %d not mine (%d) nor in overlap (owner %d)\n",
	       ient, pcr_my_proc_id(), mmpr_el_owner(mesh_id,ient));
	exit(-1);
      }
#else
      nr_mine_active++;
#endif
    }
    else if(mmr_el_status(mesh_id,ient)==MMC_INACTIVE){
#ifdef PARALLEL
      if(mmpr_el_owner(mesh_id,ient)==pcr_my_proc_id()) nr_mine_inactive++;
      else if(mmpr_el_owner(mesh_id,ient)>0) nr_overlap_inactive++;
      else {
	printf("Inactive element %d not mine (%d) nor in overlap (owner %d)\n",
	       ient, pcr_my_proc_id(), mmpr_el_owner(mesh_id,ient));
	exit(-1);
      }
#else
      nr_mine_inactive++;
#endif
    }
    else{
      printf("Element %d not FREE, ACTIVE or INACTIVE \n",
	     ient);
      exit(-1);
    }
  }
  nr_total = nr_mine_active + nr_mine_inactive + nr_overlap_active + nr_overlap_inactive + nr_free;
  printf("element statistics (owned active %d, inactive %d)\n", 
	 nr_mine_active, nr_mine_inactive);
  printf("(overlap active %d, inactive %d), free %d, total %d\n", 
	 nr_overlap_active, nr_overlap_inactive, nr_free, nr_total);
/*kew*/


  l2_err=0.0; h1_err=0.0;
  nel=0; counter=0;
  while((nel=mmr_get_next_act_elem(mesh_id,nel))!=0){

#ifdef PARALLEL
   if(mmpr_el_owner(mesh_id,nel)==pcr_my_proc_id()){

     //printf("element %d, owner %d, my_id %d\n", 
     //	    nel, mmpr_el_owner(mesh_id,nel), pcr_my_proc_id());
#endif

    counter++;

    /* find degree of polynomial and number of element scalar dofs */
    apr_get_el_pdeg(field_id, nel, &pdeg);
    num_shap = apr_get_el_pdeg_numshap(field_id,nel,&pdeg);
    nreq=apr_get_nreq(field_id);
    ndofs = nreq*num_shap;
    mat_num = mmr_el_groupID(mesh_id,nel);
    
    /* get the coordinates of the nodes of El in the right order */
    mmr_el_node_coor(mesh_id,nel,el_nodes,node_coor);

    /* get the most recent solution degrees of freedom */
    sol_vec_id = 0;
    apr_get_el_dofs(field_id, nel, sol_vec_id, dofs_loc);

    /* prepare data for gaussian integration */
    if(pdeg!=pdeg_old){
      base=apr_get_base_type(field_id, nel);
      apr_set_quadr_3D(base, &pdeg, &ngauss, xg, wg);
      pdeg_old = pdeg;
    }

/*kbw
//  if(Field_id==2 && El_id==13753)
    if(counter>=0 && nel==15425)
      {
	printf("Field_id %d, element %d, nodes: \n",
	       field_id, nel);
	for(i=0;i<ndofs;i++) printf("%d  ", el_nodes[i+1]);
	printf("\n");
	printf("DOFS :\n");
	for(i=0;i<ndofs;i++) printf("%20.15lf",dofs_loc[i]);
	printf("\n");
	//getchar();
      } 
/*kew*/

    volume=0.0; average=0.0; double err_elem_L2 = 0.0; double err_elem_H1 = 0.0;
    for (ki=0;ki<ngauss;ki++) {
    
      /* at the gauss point, compute basis functions, determinant etc*/
      iaux = 2; /* calculations with jacobian but not on the boundary */
      determ = apr_elem_calc_3D(iaux, nreq, &pdeg, base,
				&xg[3*ki],node_coor,dofs_loc,
				base_phi,base_dphix,base_dphiy,base_dphiz,
				xcoor,u_val,u_x,u_y,u_z,NULL);
    
      vol = determ * wg[ki];

      pdr_exact_sol(mat_num,xcoor[0],xcoor[1],xcoor[2],time,
      		    &f,&f_x,&f_y,&f_z,&eaux);


      volume += vol;
      average+= u_val[0]*vol;
      err_elem_L2 +=  (f-u_val[0])*(f-u_val[0])*vol;
      
      err_elem_H1 += ((f_x-u_x[0])*(f_x-u_x[0])+(f_y-u_y[0])*(f_y-u_y[0])
		      +(f_z-u_z[0])*(f_z-u_z[0]))*vol;


/*kbw
      if(counter>=0 && nel==15425){
      printf("at integration point %d, vol %lf, f %lf, sol %lf, l2_err %lf\n",
	     ki, vol, f, u_val[0], sqrt(l2_err));
      }
/*kew*/

    } /* ki */

    l2_err += err_elem_L2;
    h1_err += err_elem_H1;
/*kbw
    if(counter<10){
printf("In element %d: local L2 %lf (global %lf), local H1 %lf (global %lf)\n",
       nel, sqrt(err_elem_L2), sqrt(l2_err), sqrt(err_elem_H1), sqrt(h1_err));
}
/*kew*/

/*kbw
printf("Average in element %d = %lf (at center = %lf, element volume %lf)\n",
nel,average/volume,dofs_loc[0]+dofs_loc[1]/3.0+dofs_loc[2]/3.0,volume);
/*kew*/

#ifdef PARALLEL
   } /* end if owned element */
#endif
    
/*kbw
   printf("in active element %d - error_exact %20.15lf \n",
	  nel,err_elem_H1);
/*kew*/
  } /* end of loop over active elements */

  printf("\nnumber of owned active elements: %d\n",	 counter);

  counter = apr_get_nrdofs_glob(field_id);

#ifndef PARALLEL
  fprintf(Interactive_output,
"\nLocal number of degrees of freedom: %d\n",	 counter);
  fprintf(Interactive_output,
"Local L2 norm of error      = %20.15lf\n",	 sqrt(l2_err));
  fprintf(Interactive_output,
"Local H1 seminorm of error  = %20.15lf\n\n",	 sqrt(h1_err));
#endif

#ifdef PARALLEL
  pcr_allreduce_sum_int(1, &counter, &iaux);
  int glob_all = iaux;
  // At this point we are after freeing exchange table wiht pdv_exchange_table_index!!!
  counter = appr_get_nr_dofs_owned(pdv_exchange_table_index);
  pcr_allreduce_sum_int(1, &counter, &iaux);
  pcr_allreduce_sum_double(1, &l2_err, &daux);
  l2_err = daux;
  pcr_allreduce_sum_double(1, &h1_err, &daux);
  h1_err = daux;
 
  if(pcr_my_proc_id()==pcr_print_master()){
    fprintf(Interactive_output,
	    "\nGlobal number of degrees of freedom: %d (all %d)\n",
	    iaux, glob_all);
    fprintf(Interactive_output,
	    "Global L2 norm of error      = %20.15lf\n",	 
	    sqrt(l2_err));
    fprintf(Interactive_output,
	    "Global H1 seminorm of error  = %20.15lf\n\n",	 
	    sqrt(h1_err));
  }
  counter = iaux;
#endif

  if(Interactive_output != stdout){
    printf("\nGlobal number of degrees of freedom: %d\n",	 counter);
    printf("Local L2 norm of error      = %20.15lf\n",	 sqrt(l2_err));
    printf("Local H1 seminorm of error  = %20.15lf\n\n",	 sqrt(h1_err));
  }

  return(l2_err+h1_err);
}




/*---------------------------------------------------------
  pdr_post_process - to privide simple Interactive post-processing
----------------------------------------------------------*/
double pdr_post_process(
  int Field_id,    /* in: approximation field ID  */
  FILE *Interactive_input,
  FILE *Interactive_output
        )
{

/* local variables */

  int pdeg;		/* degree of polynomial */
  int base;		/* type of basis functions */
  int max_gen, max_gen_diff;

  int num_shap;         /* number of element shape functions */
  int ndofs;            /* local dimension of the problem */
  double determ;        /* determinant of jacobi matrix */
  double vol;           /* volume for integration rule */
  double u_val[PDC_MAXEQ]; /* computed solution */
  double u_x[PDC_MAXEQ];   /* gradient of computed solution */
  double u_y[PDC_MAXEQ];   /* gradient of computed solution */
  double u_z[PDC_MAXEQ];   /* gradient of computed solution */
  double base_phi[APC_MAXELVD];    /* basis functions */
  double base_dphix[APC_MAXELVD];  /* x-derivatives of basis function */
  double base_dphiy[APC_MAXELVD];  /* y-derivatives of basis function */
  double base_dphiz[APC_MAXELVD];  /* y-derivatives of basis function */
  int el_nodes[MMC_MAXELVNO+1];        /* list of nodes of El */
  double node_coor[3*MMC_MAXELVNO];  /* coord of nodes of El */
  double dofs_loc[APC_MAXELSD]; /* element solution dofs */
  double dofs_loc2[APC_MAXELSD]; /* element solution dofs */
  double xcoor[3];      /* global coord of gauss point */

/* for scalar test cases - exact solution and its derivatives */
  double f, f_x, f_y, f_z;
/* auxiliary variables */
  int i, j, iel, ki, iaux, mat_num, mesh_id, nel, sol_vec_id, nreq;
  double daux, eaux, volume, average, time;

  double x[3],xg[600],u_val_list[200],u_x_list[200],u_y_list[200],u_z_list[200];
  int list_el[200];

/*++++++++++++++++ executable statements ++++++++++++++++*/

  if(Interactive_input!=stdin) return(-1);
  if(Interactive_output!=stdout) return(-1);

  /* select the corresponding mesh */
  mesh_id = apr_get_mesh_id(Field_id);
  nreq=apr_get_nreq(Field_id);;

  printf("Give global coordinates of a point (x,y,z):\n");
  scanf("%lf",&x[0]);
  scanf("%lf",&x[1]);
  scanf("%lf",&x[2]);

  apr_sol_xglob(Field_id,x,1,list_el,xg,u_val_list,u_x_list,u_y_list,u_z_list,APC_CLOSE,
                APE_SOL_XGLOB_DEFAULT
                | APE_SOL_XGLOB_MATCH_ALL_ELEMENTS
                | APE_SOL_XGLOB_MATCH_WITH_ADAPTATION);

  for(iel=1; iel<=list_el[0]; iel++){

    printf("\nSolution at local point %.2lf %.2lf %.2lf in element %d:\n\n",
	   xg[(iel-1)*3+0],xg[(iel-1)*3+1],xg[(iel-1)*3+2],list_el[iel]);

    for(j=0;j<nreq;j++)  {
      printf("u_val[%d]=%lf, u_x[%d]=%lf, u_y[%d]=%lf, u_z[%d]=%lf\n",
	     j, u_val_list[(iel-1)*nreq+j], j, u_x_list[(iel-1)*nreq+j],
	     j, u_y_list[(iel-1)*nreq+j], j, u_z_list[(iel-1)*nreq+j]);
    }
  }

  printf("\nExact solution:\n");
  
  iaux = 0; daux = 0.0;
  pdr_exact_sol(iaux,x[0],x[1],x[2],daux,&f,&f_x,&f_y,&f_z,&eaux);
  
  for(j=0;j<nreq;j++) {
    printf("f[%d]=%lf, u_x[%d]=%lf, u_y[%d]=%lf, u_z[%d]=%lf\n"
	   ,j,f,j,f_x,j,f_y,j,f_z);
  }
  printf("\n");

  //  return(0);

  printf("Give element number:\n");
  scanf("%d",&nel);

  printf("Give local coordinates of a point (x,y,z):\n");
  scanf("%lf",&x[0]);
  scanf("%lf",&x[1]);
  scanf("%lf",&x[2]);
  //printf("x=%lf,y=%lf,z=%lf\n",x[0],x[1],x[2]);

  if(mmr_el_status(mesh_id, nel)==MMC_ACTIVE){

    base=apr_get_base_type(Field_id, nel);
     apr_get_el_pdeg(Field_id, nel, &pdeg);
    num_shap = apr_get_el_pdeg_numshap(Field_id, nel, &pdeg);
    ndofs = nreq*num_shap;

    /* get the coordinates of the nodes of El in the right order */
    mmr_el_node_coor(mesh_id,nel,el_nodes,node_coor);

    /* get the most recent element solution degrees of freedom */
    sol_vec_id = 1;
    apr_get_el_dofs(Field_id, nel, sol_vec_id, dofs_loc2);

    iaux = 2; /* calculations with jacobian but not on the boundary */

    apr_elem_calc_3D(iaux, nreq, &pdeg, base,
		     x,node_coor,dofs_loc2,
		     base_phi,base_dphix,base_dphiy,base_dphiz,
		     xcoor,u_val,u_x,u_y,u_z,NULL);

    printf("\nSolution at point %.2lf %.2lf %.2lf in element %d:\n\n",
	   x[0],x[1],x[2],nel);

    for(j=0;j<nreq;j++) {
      printf("u_val[%d]=%lf, u_x[%d]=%lf, u_y[%d]=%lf, u_z[%d]=%lf\n",
	     j,u_val[j],j,u_x[j],j,u_y[j],j,u_z[j]);
    }
    printf("\nGlobal coordinates of the point:  %.2lf %.2lf %.2lf:\n\n",
	   xcoor[0],xcoor[1],xcoor[2]);

    printf("\nExact solution:\n");

    iaux = 0; daux = 0.0;
    pdr_exact_sol(iaux,xcoor[0],xcoor[1],xcoor[2],daux,
		  &f,&f_x,&f_y,&f_z,&eaux);

    for(j=0;j<nreq;j++)  {
      printf("f[%d]=%lf, u_x[%d]=%lf, u_y[%d]=%lf, u_z[%d]=%lf\n"
	     ,j,f,j,f_x,j,f_y,j,f_z);
    }

    printf("\n");

  } /* end of loop over active elements */


  return(1);
}


/*---------------------------------------------------------
pdr_average_sol_el - to compute the average of solution over element
----------------------------------------------------------*/
double pdr_average_sol_el( /* returns: the average of solution over element */
        int Field_id,      /* in: data structure to be used  */
        int El          /* in: element number */
)
{

  utr_average_sol_el(Field_id, El);
  return(1);

}

/*---------------------------------------------------------
pdr_exact_sol - to return values and derivatives at a point
	for functions used as exact solutions for test problems
----------------------------------------------------------*/
int pdr_exact_sol(   /* returns: 1-found exact solution, 0-not found */
	int Flag,    /* in: flag, e.g. material number */
	double X,	/* in: coordinates of a point */
	double Y,
	double Z,
	double Time,	/* in: time instant */
	double *Exact,	/* out: values of solution and derivatives */
	double *Exact_x,
	double *Exact_y,
	double *Exact_z,
	double *Lapl	/* out: Laplacian of function */
	)
{

/* auxiliary variables */
  int i;
  double a1,daux,cd1d_coeff,x,y,z;

/*++++++++++++++++ executable statements ++++++++++++++++*/

  x=X; y=Y; z=Z;

  *Exact   = 0.0;
  *Exact_x = 0.0;
  *Exact_y = 0.0;
  *Exact_z = 0.0;
  *Lapl    = 0.0;

  // name = pdr_ctrl_i_params(PDC_USE_CURRENT_PROBLEM_ID, 1); changed for string
  char name[300];
  pdr_problem_name(PDC_USE_CURRENT_PROBLEM_ID, name);

  if(strcmp(name, "DIFF_IN_CUBE") == 0){
    //if(name==1){
/* u = exp(-(x^2+y^2)/a1) 
a1=1.0;
*Exact   =  exp((-x*x-y*y)/a1);
*Exact_x = -2.0*x*exp((-x*x-y*y)/a1)/a1;
*Exact_y = -2.0*y*exp((-x*x-y*y)/a1)/a1;
*Exact_z = 0.0;
*Lapl    = -(2.0/a1)*((1-2*x*x/a1)+(1-2*y*y/a1))
			*exp((-x*x-y*y)/a1);
/**/

/* u = exp(-(x^2+y^2+z^2)/a1) */ 
    a1=1.0;
    *Exact   =  exp((-x*x-y*y-z*z)/a1);
    *Exact_x = -2.0*x*exp((-x*x-y*y-z*z)/a1)/a1;
    *Exact_y = -2.0*y*exp((-x*x-y*y-z*z)/a1)/a1;
    *Exact_z = -2.0*z*exp((-x*x-y*y-z*z)/a1)/a1;
    *Lapl    = -(2.0/a1)*((1-2*x*x/a1)+(1-2*y*y/a1)+(1-2*z*z/a1))
			*exp((-x*x-y*y-z*z)/a1);


/* u = -(x^2+y^2+z^2) 
*Exact   =  -x*x-y*y-z*z;
*Exact_x = -2.0*x;
*Exact_y = -2.0*y;
*Exact_z = -2.0*z;
*Lapl    = -6.0;
/**/

  }
  if(strcmp(name, "DIFF_IN_BOX") == 0){
    //else if(name==2){

/* u = exp(-(x^2+y^2)/a1) 
a1=1.0;
*Exact   =  exp((-x*x-y*y)/a1);
*Exact_x = -2.0*x*exp((-x*x-y*y)/a1)/a1;
*Exact_y = -2.0*y*exp((-x*x-y*y)/a1)/a1;
*Exact_z = 0.0;
*Lapl    = -(2.0/a1)*((1-2*x*x/a1)+(1-2*y*y/a1))
			*exp((-x*x-y*y)/a1);
/**/

/* u = exp(-(x^2+y^2+z^2)/a1) */ 
    a1=100.0;
    *Exact   =  exp((-x*x-y*y-z*z)/a1);
    *Exact_x = -2.0*x*exp((-x*x-y*y-z*z)/a1)/a1;
    *Exact_y = -2.0*y*exp((-x*x-y*y-z*z)/a1)/a1;
    *Exact_z = -2.0*z*exp((-x*x-y*y-z*z)/a1)/a1;
    *Lapl    = -(2.0/a1)*((1-2*x*x/a1)+(1-2*y*y/a1)+(1-2*z*z/a1))
			*exp((-x*x-y*y-z*z)/a1);


/* u = -(x^2+y^2+z^2) 
*Exact   =  -x*x-y*y-z*z;
*Exact_x = -2.0*x;
*Exact_y = -2.0*y;
*Exact_z = -2.0*z;
*Lapl    = -6.0;
/**/

  }
  // pure convection problems not implemented yet....
  /* else if(name==101||name==201){ */
  /*   if(X<Time-200||X>Time) *Exact=0; */
  /*   else *Exact=1.0; */
  /*   *Exact_x = 0; */
  /*   *Exact_y = 0; */
  /*   *Exact_z = 0; */
  /*   *Lapl    = 0; */
  /* } */
  /* else if(name==105||name==205){ */
  /*   if(X+Y+Z<Time-0.4||X+Y+Z>Time) *Exact=0; */
  /*   else *Exact=1.0; */
  /*   *Exact_x = 0; */
  /*   *Exact_y = 0; */
  /*   *Exact_z = 0; */
  /*   *Lapl    = 0; */
  /* } */

  return(1);
}


