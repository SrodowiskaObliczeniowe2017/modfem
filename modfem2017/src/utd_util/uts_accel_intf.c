/***********************************************************************
File uts_accel_intf - utility routines supporting streaming processing 
                      on accelerators (common to all problem modules)

Contains definitions of routines:

local:
  utr_create_assemble_stiff_mat_elem_accel - to create element stiffness matrices
                           and assemble them to the global SM using ACCELERATOR
------------------------------
History:
	08.2015 - Krzysztof Banas, pobanas@cyf-kr.edu.pl, initial version
        10.2016 - Jan Bielanski, Kazimierz Chlon - assembling into crs on gpu
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<signal.h>

#ifdef _OPENMP
#include<omp.h>
#endif

#include <CL/cl.h> 

/* interface for all mesh manipulation modules */
#include "mmh_intf.h"

/* interface for all approximation modules */
#include "aph_intf.h"

/* interface for general purpose utilities - for all problem dependent modules*/
#include "uth_intf.h"
// interface for accelerator utilities - for all problem dependent modules*/
#include "./uth_accel_intf.h"

/* interface for linear algebra packages */
#include "lin_alg_intf.h"

/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include "pdh_control_intf.h"

#include "pdh_intf.h"

#include"tmh_intf.h"

// interface of opencl implementation
#include"../tmd_opencl/tmh_ocl.h"

#define TIME_TEST
#ifdef TIME_TEST
double t_begin;
double t_end;
double total_time;
#endif

// Two important switches for different variants of the algorithm:
// SWITCH 1: float versus double (MUST BE COMPATIBLE WITH KERNEL SWITCH!!!!)
// data type for integration
//#define SCALAR float
#define SCALAR double

#define SIC_MAX_DOF_PER_INT 40 /* maximal number of dof_ent (blocks) per int_ent */

// used thread management dependent acelerator routines
int tmr_prepare_create_assemble_accel(
  int Problem_id,
  int Monitor
				      );

int tmr_perform_create_assemble_accel(
  int Problem_id, 
  int Level_id, 
  int Comp_type,         /* in: indicator for the scope of computations: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int Nr_elems,
  int First_elem_in_color_index,
  int* L_elem_id,
  int* Asse_pos_first_dof_int_ent,
  int* Assembly_table,
  int Ngauss,
  SCALAR* Gauss_dat_host,
  int Num_shap,
  SCALAR* Shape_fun_dat_host,
  int Num_geo_dofs,
  SCALAR* El_geo_dat_host,
  int Num_dofs, 
  int Size_global_pde_data,  // size for all elements
  int Size_per_element_pde_data, // size for one element and all integration point
  int Size_per_int_point_pde_data,  // size for one element and one integration point
  SCALAR* El_pde_dat_host,
  int Max_dofs_int_ent
				      );

/*------------------------------------------------------------
 utr_prepare_geo_dat_for_all_color_streaming - to create table with geo data for all colors
------------------------------------------------------------*/
void  utr_prepare_geo_dat_for_all_color_streaming (
			       int Problem_id, 
			       int Level_id, 
			       int Nr_colors_accel,
			       int* L_color_index_elems,
			       int* L_int_ent_id,

			       int* Ngauss_p,
			       SCALAR** Gauss_dat_host_p,

			       int* Nshape_p,
			       SCALAR** Shape_fun_dat_host_p,

			       int* Ngeo_p,
			       int** Ngeo_color_p,
			       SCALAR** El_geo_dat_host_p						   
			      );

/*------------------------------------------------------------
 utr_create_assemble_stiff_mat_elem_accel - to create element stiffness matrices
                                 and assemble them to the global SM using ACCELERATOR
------------------------------------------------------------*/
int utr_create_assemble_stiff_mat_elem_accel(
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
  int Nr_colors_accel,
  int* Asse_pos_first_dof_int_ent,
  int* Assembly_table,
  int Max_dofs_int_ent
					 )
{

#if (defined OPENCL_CPU || defined OPENCL_GPU || defined OPENCL_PHI)

  int iselect=1;
  //printf("select: old style - 0, new style - 1:"); scanf("%d", &iselect);

  if(iselect==0){


    int icolor;
    
    
    // old style - problem dependent routine does all
    for(icolor = 0; icolor<Nr_colors_accel; icolor++){
      
      int nr_elems_this_color = L_color_index_elems[icolor+1] - L_color_index_elems[icolor];
    

/*kbw
  printf("In utr_create_assemble_stiff_mat, thread %d before opencl\n", my_thread_id);
  printf("Problem_id %d, Level_id %d, Comp_type %d, Pdeg_coarse_p %d\n",
      Problem_id, Level_id, Comp_type, *Pdeg_coarse_p );
  printf("icolor %d, nr_elems_this_color %d, L_color_index_elems[icolor] %d\n", 
      icolor, nr_elems_this_color, L_color_index_elems[icolor]);

/*kew*/

      pdr_create_assemble_stiff_mat_opencl_elem(Problem_id, Level_id, Comp_type,  
						Pdeg_coarse_p, nr_elems_this_color,
						&L_int_ent_type[L_color_index_elems[icolor]], 
						&L_int_ent_id[L_color_index_elems[icolor]],
						Max_dofs_int_ent);
      
    }

  } // end old style
  else{

    
  // new style...
#ifdef TIME_TEST
  t_begin = time_clock();
#endif

  int monitor = 3; // the level of message outputs
  // prepare integration and assembly for a given environment
  //  tmr_prepare_create_assemble_accel(Problem_id, monitor);


  // ----------- PREPARE GAUSS/SHAPE_FUN/GEO DATA FOR ALL COLORS -----------

  int ngauss; 
  SCALAR* gauss_dat_all_color_host;

  int nshape;
  SCALAR* shape_fun_dat_all_color_host;

  int ngeo;
  int* ngeo_color_p;
  SCALAR* el_geo_dat_all_color_host;
  

  // Read GAUSS/SHAPE_FUN/GEO data
  utr_prepare_geo_dat_for_all_color_streaming (
			       Problem_id, 
			       Level_id, 
			       Nr_colors_accel,
			       L_color_index_elems,
			       L_int_ent_id,
			       &ngauss,
			       &gauss_dat_all_color_host,
			       &nshape,
			       &shape_fun_dat_all_color_host,
			       &ngeo,
			       &ngeo_color_p,
			       &el_geo_dat_all_color_host
				);

  // we split processing into colors - although it is possible to prepare data for
  // all elements and then launch several kernels
  // but splitting into colors we can eventually get software pipelining
  // one kernel is executed on accelerator, while the data is prepared for another kernel
  int icolor;
  for(icolor = 0; icolor<Nr_colors_accel; icolor++){
    

    int nr_elems_this_color = L_color_index_elems[icolor+1] - L_color_index_elems[icolor];
    //int ngauss; 
    SCALAR* gauss_dat_host = gauss_dat_all_color_host;
    //int nshape;
    SCALAR* shape_fun_dat_host = shape_fun_dat_all_color_host;
    //int ngeo;
    SCALAR* el_geo_dat_host = &el_geo_dat_all_color_host[ngeo_color_p[icolor]];
    int num_dofs;
    int all_el_pde_coeff_size;
    int one_el_pde_coeff_size;
    int one_int_p_pde_coeff_size;
    SCALAR* el_pde_dat_host;
    
    // prepare geometrical data for streaming assembly (if necessary?)
    // e.g. set new constants for Comp_type - solve on the old mesh, SOLVE_OLD_MESH ?  

    /* ---> move outside the color loop
    utr_prepare_geo_dat_streaming (Problem_id, Level_id, 
				   nr_elems_this_color, 
				   &L_int_ent_id[L_color_index_elems[icolor]],
				   &ngauss, &gauss_dat_host,
				   &nshape, &shape_fun_dat_host,
				   &ngeo, &el_geo_dat_host);
    /* */

/*kbw*/
    printf("\nReceived from utr_prepare_geo_dat_streaming (Problem %d, Level %d, Nr_elems %d:\n",
	   Problem_id, Level_id,  nr_elems_this_color);
    printf("ngauss %d, gauss_dat_host %lu\n",
	   ngauss, gauss_dat_host);
    printf("nshape %d, shape_fun_host %lu\n",
	   nshape, shape_fun_dat_host);
    printf("ngeo (for each element) %d, el_geo_dat_host %lu\n",
	   ngeo, el_geo_dat_host);

/*kew*/


    // it is possible to get software pipelining also here:
    // send geometry input data to accelerator
    // launch (asynchronously) kernel for geometry calculations
    // prepare (concurrently) PDE data
    // send PDE data to accelerator

/*kbw
    printf("in parallel %d\n",omp_in_parallel());
    printf("\nBefore pdr_prepare_pde_coeff_streaming (Problem %d, Level %d, Nr_elems %d:\n",
	   Problem_id, Level_id,  nr_elems_this_color);
    printf("num_dofs %d, all_el_pde_coeff_size %d\n",
	   num_dofs, all_el_pde_coeff_size);
    printf("one_el_pde_coeff_size %d, one_int_p_pde_coeff_size %d\n",
	   one_el_pde_coeff_size, one_int_p_pde_coeff_size);
    printf("el_pde_dat_host %lu\n",
	   el_pde_dat_host);

/*kew*/

    // prepare PDE dependent data for streaming assembly (if necessary?)
    // (depending on problem name we decide whether new coefficients are necessary)
    pdr_prepare_pde_coeff_streaming(Problem_id, Level_id, Comp_type,
				    nr_elems_this_color, 
				    &L_int_ent_id[L_color_index_elems[icolor]],
				    &num_dofs, &all_el_pde_coeff_size,
				    &one_el_pde_coeff_size, &one_int_p_pde_coeff_size,
				    &el_pde_dat_host);

/*kbw*/
    printf("\nReceived from pdr_prepare_pde_coeff_streaming (Problem %d, Level %d, Nr_elems %d:\n",
	   Problem_id, Level_id,  nr_elems_this_color);
    printf("num_dofs %d, all_el_pde_coeff_size %d\n",
	   num_dofs, all_el_pde_coeff_size);
    printf("one_el_pde_coeff_size %d, one_int_p_pde_coeff_size %d\n",
	   one_el_pde_coeff_size, one_int_p_pde_coeff_size);
    printf("el_pde_dat_host %lu\n",
	   el_pde_dat_host);

/*kew*/


/*kbw
  printf("In utr_create_assemble_stiff_mat, thread %d before opencl\n", my_thread_id);
  printf("Problem_id %d, Level_id %d, Comp_type %d, Pdeg_coarse_p %d\n",
      Problem_id, Level_id, Comp_type, *Pdeg_coarse_p );
  printf("icolor %d, nr_elems_this_color %d, L_color_index_elems[icolor] %d\n", 
      icolor, nr_elems_this_color, L_color_index_elems[icolor]);

/*kew*/

    /*jbw
    printf("\n --> color: %d ; L_color_index_elems[%d]=%d ; L_int_ent_id[%d]=%d\n\n",icolor,icolor,L_color_index_elems[icolor],L_color_index_elems[icolor],L_int_ent_id[L_color_index_elems[icolor]]);
    /*jbw*/

    tmr_perform_create_assemble_accel(Problem_id, Level_id, Comp_type,  
				      nr_elems_this_color,
				      L_color_index_elems[icolor],
				      &L_int_ent_id[L_color_index_elems[icolor]],
				      Asse_pos_first_dof_int_ent, Assembly_table, 
				      ngauss, gauss_dat_host,
				      nshape, shape_fun_dat_host,
				      ngeo, el_geo_dat_host,
				      num_dofs, all_el_pde_coeff_size,
				      one_el_pde_coeff_size, one_int_p_pde_coeff_size,
				      el_pde_dat_host,
				      Max_dofs_int_ent);  

   

    /*
    free(gauss_dat_host);
    free(shape_fun_dat_host);
    free(el_geo_dat_host);
    */
    free(el_pde_dat_host);

  }

  free(gauss_dat_all_color_host);
  free(shape_fun_dat_all_color_host);
  free(el_geo_dat_all_color_host);
  
#else
  printf("Accelerator integration and assembly called without linked code!!!\n");
  exit(-1);
#endif



/* ------------------------------------------------------------------------- GPU SOLVER ------------------------------------------------------------------------------ */
#ifdef GPU_ASSEMBLING

#ifdef GPU_DEBUG

  int i;

  // choose the platform 
  int platform_index = tmr_ocl_get_current_platform_index();
  int device_tmc_type = tmr_ocl_get_current_device_type();
  int device_index = tmr_ocl_get_current_device(); //tmr_ocl_select_device(platform_index, device_tmc_type);
  cl_context context = tmr_ocl_select_context(platform_index, device_index);
  cl_command_queue command_queue = tmr_ocl_select_command_queue(platform_index, device_index);

  FILE *test_out;

  /*jbw*/
  test_out = fopen("gpu_data_asse.txt","w+");
  
  // Test ALLOCATE MEMEORY FOR ASSEMBLY TABLE 1
  int size_asse = tmv_ocl_crs_struct.ocl_asse_pos_first_dof_int_ent_bytes;
  int *test_asse_pos_first_dof_int_ent = (int*) malloc (size_asse);

  clEnqueueReadBuffer(command_queue, tmv_ocl_crs_struct.ocl_asse_pos_first_dof_int_ent, CL_TRUE, 0, size_asse, test_asse_pos_first_dof_int_ent, 0, NULL, NULL);

  fprintf(test_out,"tmv_ocl_crs_struct.ocl_asse_pos_first_dof_int_ent: \n");
  for(i=0; i<(size_asse/sizeof(int)); i++) {
    fprintf(test_out,"%d ",test_asse_pos_first_dof_int_ent[i]);
  }
  fprintf(test_out,"\n\n");

  int size_asse_t = tmv_ocl_crs_struct.ocl_assembly_table_bytes;
  int *test_ocl_assembly_table = (int*) malloc (size_asse_t);

  clEnqueueReadBuffer(command_queue, tmv_ocl_crs_struct.ocl_assembly_table, CL_TRUE, 0, size_asse_t, test_ocl_assembly_table, 0, NULL, NULL);

  fprintf(test_out,"tmv_ocl_crs_struct.ocl_assembly_table: \n");
  for(i=0; i<(size_asse_t/sizeof(int)); i++) {
    fprintf(test_out,"%d ",test_ocl_assembly_table[i]);
  }
  fprintf(test_out,"\n\n");

  fclose(test_out);


  // Test ALLOCATE MEMEORY CRS STRUCTURE 2
  test_out = fopen("gpu_data_crs.txt","w+");
  
  int size_col = tmv_ocl_crs_struct.ocl_crs_col_ind_bytes;
  int *test_ocl_crs_col = (int*) malloc (size_col);

  clEnqueueReadBuffer(command_queue, tmv_ocl_crs_struct.ocl_crs_col_ind, CL_TRUE, 0, size_col, test_ocl_crs_col, 0, NULL, NULL);
  fprintf(test_out,"tmv_ocl_crs_struct.ocl_crs_col_ind: \n");
  for(i=0; i<(size_col/sizeof(int)); i++) {
    fprintf(test_out,"%d ",test_ocl_crs_col[i]);
  }
  fprintf(test_out,"\n\n");
  
  int size_row = tmv_ocl_crs_struct.ocl_crs_row_ptr_bytes;
  int *test_ocl_crs_row = (int*) malloc (size_row);
  clEnqueueReadBuffer(command_queue, tmv_ocl_crs_struct.ocl_crs_row_ptr, CL_TRUE, 0, size_row, test_ocl_crs_row, 0, NULL, NULL);
  fprintf(test_out,"tmv_ocl_crs_struct.ocl_crs_row_ptr: \n");
  for(i=0; i<(size_row/sizeof(int)); i++) {
    fprintf(test_out,"%d ",test_ocl_crs_row[i]);
  }
  fprintf(test_out,"\n\n");
  fclose(test_out );

  // Test ALLOCATE MEMEORY RHS
  test_out = fopen("gpu_data_rhs.txt","w+");
  
  int size_rhs = tmv_ocl_crs_struct.ocl_rhs_bytes;

  clEnqueueReadBuffer(command_queue, tmv_ocl_crs_struct.ocl_rhs_val, CL_TRUE, 0, size_rhs, tmv_ocl_crs_struct.rhs_val_gpu, 0, NULL, NULL);
  fprintf(test_out,"tmv_ocl_crs_struct.ocl_rhs: \n");
  for(i=0; i<(size_rhs/sizeof(double)); i++) {
    fprintf(test_out,"(%d: %lf); \n ",i,tmv_ocl_crs_struct.rhs_val_gpu[i]);
  }
  fprintf(test_out,"\n\n");
  
  fprintf(test_out,"\n\n");
  fclose(test_out );

  free(test_asse_pos_first_dof_int_ent);
  free(test_ocl_assembly_table);
  free(test_ocl_crs_col);
  free(test_ocl_crs_row);
  
  /**/
#endif

  // Copy the rest of CRS data from CPU to GPU
  // tmr_perform_copy_crs(...)
  // !!! TODO !!!
  //

#ifdef GPU_DEBUG
  
  clEnqueueReadBuffer(command_queue, tmv_ocl_crs_struct.ocl_crs_val, CL_TRUE, 0, tmv_ocl_crs_struct.ocl_crs_val_bytes, tmv_ocl_crs_struct.crs_val_gpu, 0, NULL, NULL);

  test_out = fopen("gpu_data_crs.txt","w+");
  
  fprintf(test_out,"tmv_ocl_crs_struct.ocl_crs_val: \n");
  for(i=0; i<tmv_ocl_crs_struct.Nnz; i++) {
    fprintf(test_out,"(%d: %lf) \n; ",i, tmv_ocl_crs_struct.crs_val_gpu[i]);  
  }

  fclose(test_out);

#endif


  // Prepare data for solver
  //tmr_prepare_final_crs();



#ifdef GPU_DEBUG

  test_out = fopen("gpu_data_crs_finalized.txt","w+");
  
  fprintf(test_out,"tmv_ocl_crs_struct.ocl_crs_val: \n");
  for(i=0; i<tmv_ocl_crs_struct.Nnz; i++) {
    fprintf(test_out,"(%d: %lf) \n; ",i, tmv_ocl_crs_struct.crs_val_cpu[i]);  
  }

  fclose(test_out);

#endif


  // Initialize and create SOLVER kernel
  // tmr_perform_solve(...)
  // !!! TODO !!!
  //


  // Read the result of calculations from GPU
  // tmr_read_solution
  // !!! TODO !!!
  //


  // CLEANUP CRS STRUCTURE
  //tmr_cleanup_crs();
  
#endif
  /* -------------------------------------------------------------------- !!END!! GPU SOLVER !!END!!  ------------------------------------------------------------------------- */

  } // end new style

  return(0);

}

/*------------------------------------------------------------
 utr_prepare_geo_dat_for_all_color_streaming - to create table with geo data for all colors
------------------------------------------------------------*/
void  utr_prepare_geo_dat_for_all_color_streaming (
			       int Problem_id, 
			       int Level_id, 
			       int Nr_colors_accel,
			       int* L_color_index_elems,
			       int* L_int_ent_id,

			       int* Ngauss_p,
			       SCALAR** Gauss_dat_host_p,

			       int* Nshape_p,
			       SCALAR** Shape_fun_dat_host_p,

			       int* Ngeo_p,
			       int** Ngeo_color_p,
			       SCALAR** El_geo_dat_host_p
			      )
{
  int i,j;
  int icolor, ielem;
  
  int mesh_id=pdr_ctrl_i_params(Problem_id,2);
  int field_id=pdr_ctrl_i_params(Problem_id,3);

  int nr_elems_all_gpu_colors; // Number of elements for all colors
  int nr_elems_this_color; // Number of elements per color
  int* L_elem_id; // Element id

  SCALAR *gauss_dat_host; int ref_el_quadr_dat_size;

  SCALAR *shape_fun_host;
  int ref_el_shape_fun_size; int shape_fun_host_bytes;

  int num_geo_dofs; int one_el_geo_dat_size; int position_geo_dofs;
  SCALAR* el_geo_dat; int *el_geo_pos_color;


  // Number of elements for all colors
  nr_elems_all_gpu_colors = L_color_index_elems[Nr_colors_accel] - L_color_index_elems[0];
  
  // 1. ------------- QUADRATURE DATA AND JACOBIAN TERMS -------------

  L_elem_id = &L_int_ent_id[L_color_index_elems[0]];
  
  // get the size of quadrature data (we assume all elements are of the same type)
  int base = apr_get_base_type(field_id, L_elem_id[0]);
  int ngauss;            /* number of gauss points */
  double xg[300];   	 /* coordinates of gauss points in 3D */
  double wg[100];       /* gauss weights */
  
  // choose an example element for a given pdeg and color
  int el_id = L_elem_id[0];
  int pdeg = apr_get_el_pdeg(field_id, el_id, NULL);  
  int num_shap = apr_get_el_pdeg_numshap(field_id, el_id, &pdeg); 
  apr_set_quadr_3D(base, &pdeg, &ngauss, xg, wg);
  
  // we may need quadrature data for the reference element
  // for each gauss point its coordinates and weight are sent for reference element
  ref_el_quadr_dat_size = ngauss*4;
  
  gauss_dat_host = (SCALAR*) malloc(ngauss*4*sizeof(SCALAR));

  /*kbw*/
  printf("Allocated %d bytes for gauss_dat_host %lu\n",
	 ngauss*4*sizeof(SCALAR), gauss_dat_host);
  /*kew*/


  // This is the version when we send only necessary variables but we do not care about sending
  for(i=0; i<ngauss; i++){
      gauss_dat_host[4*i] = xg[3*i];
      gauss_dat_host[4*i+1] = xg[3*i+1];
      gauss_dat_host[4*i+2] = xg[3*i+2];
      gauss_dat_host[4*i+3] = wg[i];
  }
  
  /*kbw
    for(i=0;i<5;i++){
  
    printf("gauss_dat_host[%d] = %f\n", i, gauss_dat_host[i]);
    //printf("gauss_dat_host[%d] = %lf\n", i, gauss_dat_host[i]);
  
    }
    for(i=ngauss-5;i<ngauss;i++){
  
    printf("gauss_dat_host[%d] = %f\n", i, gauss_dat_host[i]);
    //printf("gauss_dat_host[%d] = %lf\n", i, gauss_dat_host[i]);
  
    }
  //kew*/
  
  // 2. ------------- SHAPE FUNCTIONS ---------------

  // space for element shape functions' values and derivatives in global memory
  // we need all shape functions and their derivatives
  // at all integration points for the reference element 
  ref_el_shape_fun_size = 4*num_shap*ngauss; // for the reference element
  
  // !!!!!!!!!!!!!!!***************!!!!!!!!!!!!!!!!!
  // fill and send buffers with shape function values for the reference element
  shape_fun_host_bytes = ref_el_shape_fun_size * sizeof(SCALAR);
  shape_fun_host = (SCALAR*)malloc(shape_fun_host_bytes);
  
/*kbw*/
  printf("Allocated %d bytes for shape_fun_host %lu\n",
	 shape_fun_host_bytes, shape_fun_host);
/*kew*/

  // shape functions and derivatives are computed here but used also later on
  double base_phi_ref[APC_MAXELVD];    /* basis functions */
  double base_dphix_ref[APC_MAXELVD];  /* x-derivatives of basis function */
  double base_dphiy_ref[APC_MAXELVD];  /* y-derivatives of basis function */
  double base_dphiz_ref[APC_MAXELVD];  /* z-derivatives of basis function */
  
  // we need all shape functions and their derivatives
  // at all integration points for the reference element
  int ki;
  for (ki=0;ki<ngauss;ki++) {
    
    int temp = apr_shape_fun_3D(base, pdeg, &xg[3*ki],
				base_phi_ref, base_dphix_ref,base_dphiy_ref,base_dphiz_ref);
    assert(temp==num_shap);
    
    for(i=0;i<num_shap;i++){
      
      shape_fun_host[ki*4*num_shap+4*i] = base_phi_ref[i];
      shape_fun_host[ki*4*num_shap+4*i+1] = base_dphix_ref[i];
      shape_fun_host[ki*4*num_shap+4*i+2] = base_dphiy_ref[i];
      shape_fun_host[ki*4*num_shap+4*i+3] = base_dphiz_ref[i];   
    } 
  }
  
/*kbw
	for(i=0;i<5;i++){
	
	printf("shape_fun_host[%d] = %f\n", i, shape_fun_host[i]);
	//printf("shape_fun_host[%d] = %lf\n", i, shape_fun_host[i]);
	
	}
	for(i=num_shap-5;i<num_shap;i++){
	
	printf("shape_fun_host[%d] = %f\n", i, shape_fun_host[i]);
	//printf("shape_fun_host[%d] = %lf\n", i, shape_fun_host[i]);
	
	}
//kew*/
  
  // 3. ------------- GEO_DOFS ---------------
  
  // Create data structure for geo data
  L_elem_id = &L_int_ent_id[L_color_index_elems[0]];
  num_geo_dofs = mmr_el_node_coor(mesh_id, L_elem_id[0], NULL, NULL);
  one_el_geo_dat_size = 3*num_geo_dofs;

  el_geo_dat = (SCALAR *)malloc(3*num_geo_dofs*nr_elems_all_gpu_colors*sizeof(SCALAR));
  el_geo_pos_color = (int*)malloc(Nr_colors_accel*sizeof(int));

  /*kbw*/
  printf("Allocated %d bytes for el_geo_dat_host %lu\n",
	 3*num_geo_dofs*nr_elems_all_gpu_colors*sizeof(SCALAR), el_geo_dat);
  /*kew*/

  position_geo_dofs = 0;
    
  // Loop over colors for accelerator
  for(icolor = 0; icolor<Nr_colors_accel; icolor++){

    nr_elems_this_color = L_color_index_elems[icolor+1] - L_color_index_elems[icolor]; // Elements per color
    L_elem_id = &L_int_ent_id[L_color_index_elems[icolor]]; // Element id

    el_geo_pos_color[icolor] = position_geo_dofs;
    
    for(ielem=0; ielem<nr_elems_this_color; ielem++) {
    
      // element ID
      int el_id = L_elem_id[ielem];
    
      // checking whether this element has the same data as assumed for this color
      assert( pdeg == apr_get_el_pdeg(field_id, el_id, NULL) );
      assert( num_shap == apr_get_el_pdeg_numshap(field_id, el_id, &pdeg) );
      assert( base == apr_get_base_type(field_id, el_id) );
    
      // element geo_dofs for computing Jacobian terms on GPU
    
      // IDs of element vertices (nodes) and their coordinates as geo_dofs
      double geo_dofs[3*MMC_MAXELVNO];  /* coord of nodes of El */
      mmr_el_node_coor(mesh_id, el_id, NULL, geo_dofs);
    
    
      for(i=0;i<num_geo_dofs;i++){
	el_geo_dat[position_geo_dofs] = geo_dofs[3*i];
	el_geo_dat[position_geo_dofs+1] = geo_dofs[3*i+1];
	el_geo_dat[position_geo_dofs+2] = geo_dofs[3*i+2];
	position_geo_dofs += 3;
      }
      
    }
  }

  /*kbw*/
  printf("\nData structure sizes:\n");
  printf("\tQuadrature data for reference element %d\n", 
	 ref_el_quadr_dat_size);
  printf("\tShape functions and derivatives for reference element %d\n",
	 ref_el_shape_fun_size);
  printf("\tGeo dofs for each element: %d\n", one_el_geo_dat_size);
  /*kew*/
  
  *Ngauss_p = ngauss;
  *Gauss_dat_host_p = gauss_dat_host;
  *Nshape_p = num_shap;
  *Shape_fun_dat_host_p = shape_fun_host;
  *Ngeo_p = num_geo_dofs;
  *Ngeo_color_p = el_geo_pos_color;
  *El_geo_dat_host_p = el_geo_dat;

}


utr_prepare_geo_dat_streaming (
			       int Problem_id, 
			       int Level_id, 
			       int Nr_elems,
			       int* L_elem_id,
			       int *Ngauss_p,
			       SCALAR** Gauss_dat_host_p,
			       int* Nshape_p,
			       SCALAR** Shape_fun_dat_host_p,
			       int* Ngeo_p,
			       SCALAR** El_geo_dat_host_p
			       )
{
  int i;
  int mesh_id=pdr_ctrl_i_params(Problem_id,2);
  int field_id=pdr_ctrl_i_params(Problem_id,3);
  
  // 1. QUADRATURE DATA AND JACOBIAN TERMS
  
  // get the size of quadrature data (we assume all elements are of the same type)
  int base = apr_get_base_type(field_id, L_elem_id[0]);
  int ngauss;            /* number of gauss points */
  double xg[300];   	 /* coordinates of gauss points in 3D */
  double wg[100];       /* gauss weights */
  
  // choose an example element for a given pdeg and color
  int el_id = L_elem_id[0];
  int pdeg = apr_get_el_pdeg(field_id, el_id, NULL);  
  int num_shap = apr_get_el_pdeg_numshap(field_id, el_id, &pdeg); 
  apr_set_quadr_3D(base, &pdeg, &ngauss, xg, wg);
  
  // we may need quadrature data for the reference element
  // for each gauss point its coordinates and weight are sent for reference element
  int ref_el_quadr_dat_size = ngauss*4;
  
  SCALAR *gauss_dat_host = (SCALAR*) malloc(ngauss*4*sizeof(SCALAR));

/*kbw*/
  printf("Allocated %d bytes for gauss_dat_host %lu\n",
	 ngauss*4*sizeof(SCALAR), gauss_dat_host);
/*kew*/


  // This is the version when we send only necessary variables but we do not care about sending
  for(i=0; i<ngauss; i++){
      gauss_dat_host[4*i] = xg[3*i];
      gauss_dat_host[4*i+1] = xg[3*i+1];
      gauss_dat_host[4*i+2] = xg[3*i+2];
      gauss_dat_host[4*i+3] = wg[i];
  }
  
/*kbw
  for(i=0;i<5;i++){
  
  printf("gauss_dat_host[%d] = %f\n", i, gauss_dat_host[i]);
  //printf("gauss_dat_host[%d] = %lf\n", i, gauss_dat_host[i]);
  
  }
  for(i=ngauss-5;i<ngauss;i++){
  
  printf("gauss_dat_host[%d] = %f\n", i, gauss_dat_host[i]);
  //printf("gauss_dat_host[%d] = %lf\n", i, gauss_dat_host[i]);
  
  }
//kew*/

  // 2. SHAPE FUNCTIONS
  
  
  // space for element shape functions' values and derivatives in global memory
  // we need all shape functions and their derivatives
  // at all integration points for the reference element 
  int ref_el_shape_fun_size = 4*num_shap*ngauss; // for the reference element
  
  
  // !!!!!!!!!!!!!!!***************!!!!!!!!!!!!!!!!!
  // fill and send buffers with shape function values for the reference element
  int shape_fun_host_bytes = ref_el_shape_fun_size * sizeof(SCALAR);
  SCALAR *shape_fun_host = (SCALAR*)malloc(shape_fun_host_bytes);
  

/*kbw*/
  printf("Allocated %d bytes for shape_fun_host %lu\n",
	 shape_fun_host_bytes, shape_fun_host);
/*kew*/

  // shape functions and derivatives are computed here but used also later on
  double base_phi_ref[APC_MAXELVD];    /* basis functions */
  double base_dphix_ref[APC_MAXELVD];  /* x-derivatives of basis function */
  double base_dphiy_ref[APC_MAXELVD];  /* y-derivatives of basis function */
  double base_dphiz_ref[APC_MAXELVD];  /* z-derivatives of basis function */
  
  // we need all shape functions and their derivatives
  // at all integration points for the reference element
  int ki;
  for (ki=0;ki<ngauss;ki++) {
    
    int temp = apr_shape_fun_3D(base, pdeg, &xg[3*ki],
				base_phi_ref, base_dphix_ref,base_dphiy_ref,base_dphiz_ref);
    
    assert(temp==num_shap);
    
    for(i=0;i<num_shap;i++){
      
      shape_fun_host[ki*4*num_shap+4*i] = base_phi_ref[i];
      shape_fun_host[ki*4*num_shap+4*i+1] = base_dphix_ref[i];
      shape_fun_host[ki*4*num_shap+4*i+2] = base_dphiy_ref[i];
      shape_fun_host[ki*4*num_shap+4*i+3] = base_dphiz_ref[i];
      
    }
    
  }
  
  
  
/*kbw
	for(i=0;i<5;i++){
	
	printf("shape_fun_host[%d] = %f\n", i, shape_fun_host[i]);
	//printf("shape_fun_host[%d] = %lf\n", i, shape_fun_host[i]);
	
	}
	for(i=num_shap-5;i<num_shap;i++){
	
	printf("shape_fun_host[%d] = %f\n", i, shape_fun_host[i]);
	//printf("shape_fun_host[%d] = %lf\n", i, shape_fun_host[i]);
	
	}
//kew*/

  
  // 3. GEO_DOFS    
    
  // get the size of geometry data for one element - we assume multi-linear elements
  /* for geometrically (multi)linear elements number of geometrical  */
  /* degrees of freedom is equal to the number of vertices - classical FEM nodes */

  // choose an example element for a given pdeg and color
  el_id = L_elem_id[0];
  int num_geo_dofs = mmr_el_node_coor(mesh_id, el_id, NULL, NULL);
  
  int one_el_geo_dat_size = 3*num_geo_dofs;
  
  SCALAR* el_geo_dat = (SCALAR *)malloc(3*num_geo_dofs*Nr_elems*sizeof(SCALAR));

/*kbw*/
  printf("Allocated %d bytes for el_geo_dat_host %lu\n",
	 3*num_geo_dofs*Nr_elems*sizeof(SCALAR), el_geo_dat);
/*kew*/
  

  int ielem;
  int position_geo_dofs = 0;
  for(ielem=0; ielem<Nr_elems; ielem++){
    
    // element ID
    int el_id = L_elem_id[ielem];
    
    // checking whether this element has the same data as assumed for this color
    assert( pdeg == apr_get_el_pdeg(field_id, el_id, NULL) );
    assert( num_shap == apr_get_el_pdeg_numshap(field_id, el_id, &pdeg) );
    assert( base == apr_get_base_type(field_id, el_id) );
    
    // element geo_dofs for computing Jacobian terms on GPU
    
    // IDs of element vertices (nodes) and their coordinates as geo_dofs
    double geo_dofs[3*MMC_MAXELVNO];  /* coord of nodes of El */
    mmr_el_node_coor(mesh_id, el_id, NULL, geo_dofs);
    
    
    for(i=0;i<num_geo_dofs;i++){
      el_geo_dat[position_geo_dofs] = geo_dofs[3*i];
      el_geo_dat[position_geo_dofs+1] = geo_dofs[3*i+1];
      el_geo_dat[position_geo_dofs+2] = geo_dofs[3*i+2];
      position_geo_dofs += 3;
    }
    
  }
    
/*kbw*/
  printf("\nData structure sizes:\n");
  printf("\tQuadrature data for reference element %d\n", 
	 ref_el_quadr_dat_size);
  printf("\tShape functions and derivatives for reference element %d\n",
	 ref_el_shape_fun_size);
  printf("\tGeo dofs for each element: %d\n", one_el_geo_dat_size);
/*kew*/


  *Ngauss_p = ngauss;
  *Gauss_dat_host_p = gauss_dat_host;
  *Nshape_p = num_shap;
  *Shape_fun_dat_host_p = shape_fun_host;
  *Ngeo_p = num_geo_dofs;
  *El_geo_dat_host_p = el_geo_dat;

  return(0);
}
  

