/************************************************************************
File pds_conv_diff_ls_intf_opencl - interactions with linear solver
                          for conv_diff module implemented for opencl

Contains definitions of routines:


  pdr_create_assemble_stiff_mat_opencl_elem - 
                                 to create element stiffness matrices
                                 and assemble them to the global SM

------------------------------
History:
	08.2008 - Krzysztof Banas, pobanas@cyf-kr.edu.pl, initial version
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<limits.h>

#ifdef _OPENMP
#include<omp.h>
#endif

#include"tmh_intf.h"

#include <CL/cl.h> 

// interface of opencl implementation
#include"../../../../tmd_opencl/tmh_ocl.h"

/* interface for all mesh manipulation modules */
#include "mmh_intf.h"

/* interface for all approximation modules */
#include "aph_intf.h"


/* interface for general purpose utilities - for all problem dependent modules*/
#include "uth_intf.h"

/* interface for linear algebra packages */
#include "lin_alg_intf.h"

/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include "pdh_control_intf.h"

#include "pdh_intf.h"


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

// Less important switches - hacks for specific versions of kernels
// SWITCH 3: generic conv-diff (with plenty of coeffcients) versus Laplace
//#define GENERIC_CONV_DIFF
#define LAPLACE
// artificial example - coefficients constant for all integration points
//#define TEST_SCALAR


int pdr_prepare_pde_coeff_streaming(
  int Problem_id,
  int Level_id, 
  int Comp_type,         /* in: indicator for the scope of computations: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int Nr_elems,
  int* L_elem_id,
  int* Size_SM_LV_p, 
  int* Size_global_pde_data_p,
  int* Size_per_element_pde_data_p,
  int* Size_per_int_point_pde_data_p,
  SCALAR** El_pde_dat_host_p
				    )
{
  int i,j,k;
  
/*kbw
    printf("\nEntering pdr_prepare_pde_coeff_streaming (Problem %d, Level %d, Nr_elems %d:\n",
	   Problem_id, Level_id,  Nr_elems);
    printf("num_dofs %d, all_el_pde_coeff_size %d\n",
	   * Size_SM_LV_p, * Size_global_pde_data_p);
    printf("one_el_pde_coeff_size %d, one_int_p_pde_coeff_size %d\n",
	   * Size_per_element_pde_data_p, * Size_per_int_point_pde_data_p);
    printf("el_pde_dat_host %lu\n",
	   * El_pde_dat_host_p);

/*kew*/

  // find out the module name
  char module_name[100];
  pdr_module_introduce(module_name);
  
  // get problem name
  char problem_name[100];
  pdr_problem_name(Problem_id, problem_name);
  
  //
  // PREPARE PROBLEM DEPENDENT DATA FOR INTEGRATION
  //
  
  int field_id = pdr_ctrl_i_params(Problem_id, 3);
  int mesh_id = apr_get_mesh_id(field_id);
  int nreq = apr_get_nreq(field_id);
  
  if(nreq>PDC_MAXEQ){
    printf("nreq (%d) > PDC_MAXEQ (%d) in pdr_create_assemble_stiff_mat_opencl_elem\n",
	   nreq, PDC_MAXEQ);
    printf("Exiting!\n"); exit(-1);
  }
  
  // get the active PDE coefficient matrices
  /* pde coefficients */
  static int coeff_ind = 0;
  int coeff_ind_vect[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  // 1 - mval, 2 - axx, 3 - axy, 4 - axz, 5 - ayx, 6 - ayy, 7 - ayz, 
  // 8 - azx, 9 - azy, 10 - azz, 11 - bx, 12 - by, 13 - bz
  // 14 - tx, 15 - ty, 16 - tz, 17 - cval
  // 18 - lval, 19 - qx, 20 - qy, 21 - qz, 22 - sval 
  
  double axx[PDC_MAXEQ*PDC_MAXEQ];
  double axy[PDC_MAXEQ*PDC_MAXEQ];
  double axz[PDC_MAXEQ*PDC_MAXEQ];
  double ayx[PDC_MAXEQ*PDC_MAXEQ];
  double ayy[PDC_MAXEQ*PDC_MAXEQ];
  double ayz[PDC_MAXEQ*PDC_MAXEQ];
  double azx[PDC_MAXEQ*PDC_MAXEQ];
  double azy[PDC_MAXEQ*PDC_MAXEQ];
  double azz[PDC_MAXEQ*PDC_MAXEQ];
  double bx[PDC_MAXEQ*PDC_MAXEQ];
  double by[PDC_MAXEQ*PDC_MAXEQ];
  double bz[PDC_MAXEQ*PDC_MAXEQ];
  double tx[PDC_MAXEQ*PDC_MAXEQ];
  double ty[PDC_MAXEQ*PDC_MAXEQ];
  double tz[PDC_MAXEQ*PDC_MAXEQ];
  double cval[PDC_MAXEQ*PDC_MAXEQ];
  double mval[PDC_MAXEQ*PDC_MAXEQ];
  double qx[PDC_MAXEQ];
  double qy[PDC_MAXEQ];
  double qz[PDC_MAXEQ];
  double sval[PDC_MAXEQ];
  double lval[PDC_MAXEQ];
  
  pdr_select_el_coeff_vect(Problem_id, coeff_ind_vect); 
  
  // there are two choices:
  // 1. consider all terms (16 arrays and 4 vectors)
  // 2. use coeff_vector_indicator to select which terms are non-zero
  int coeff_array_indicator[16];
  int nr_coeff_arrays = 0;
  for(i=0; i<16; i++) {
    coeff_array_indicator[i]=coeff_ind_vect[i+2];
    if(coeff_ind_vect[i+2]==1) nr_coeff_arrays++;
  }
  // in create_assemble we combine mval (1) with cval (17)
  // we use this for vectorization of kernels
  if(coeff_ind_vect[1]==1) {
    if(coeff_array_indicator[15]!=1) nr_coeff_arrays++;
    coeff_array_indicator[15]=1;
  }
  int coeff_vector_indicator[4];
  int nr_coeff_vectors = 0; 
  for(i=0; i<4; i++) {
    coeff_vector_indicator[i]=coeff_ind_vect[i+19];
    if(coeff_ind_vect[i+19]==1) nr_coeff_vectors++;
  }
  // in create_assemble we combine lval (18) with sval (22)
  // we use this for vectorization of kernels
  if(coeff_ind_vect[18]==1) {
    if(coeff_vector_indicator[3]!=1) nr_coeff_vectors++;
    coeff_vector_indicator[3]=1;
  }
  

  
/*kbw*/
  printf("Problem ID: %d, problem_name %s\n", Problem_id, problem_name);
  printf("\nthe number of coefficients arrays %d, indicator = \n",nr_coeff_arrays);
  for(i=0;i<16;i++){
    printf("%2d",coeff_array_indicator[i]);
  }
  printf("\n");
  printf("the number of coefficients vectors %d, indicator = \n",nr_coeff_vectors);
  for(i=0;i<4;i++){
    printf("%2d",coeff_vector_indicator[i]);
  }
  printf("\n");
/*kew*/

  //
  // BASED ON PDEG COMPUTE EXECUTION CHARACTERISTICS
  //
  
  // choose an example element for a given pdeg and color
  int el_id = L_elem_id[0];
  
  int pdeg = apr_get_el_pdeg(field_id, el_id, NULL);  
  int num_shap = apr_get_el_pdeg_numshap(field_id, el_id, &pdeg); 
  int num_dofs = num_shap * nreq;
  
  // get the size of quadrature data (we assume all elements are of the same type)
  int base = apr_get_base_type(field_id, L_elem_id[0]);
  int ngauss;            /* number of gauss points */
  double xg[300];   	 /* coordinates of gauss points in 3D */
  double wg[100];       /* gauss weights */
  apr_set_quadr_3D(base, &pdeg, &ngauss, xg, wg);
  
  int el_nodes[MMC_MAXELVNO+1];        /* list of nodes of El */
  mmr_el_node_coor(mesh_id, el_id, el_nodes, NULL);
  /* for geometrically (multi)linear elements number of geometrical  */
  /* degrees of freedom is equal to the number of vertices - classical FEM nodes */
  int num_geo_dofs = el_nodes[0];
  int pdeg_single = pdeg;
  if(pdeg>100) {
    pdeg_single = pdeg/100;
    if(pdeg != pdeg_single*100 + pdeg_single){
      printf("wrong pdeg %d ( > 100 but not x0x )\n", pdeg);
      exit(-1); 
    }
  }
  
  /*kbw*/
  printf("problem and element characteristics: nreq %d, pdeg %d, num_shap %d, num_dofs %d!\n",
	 nreq, pdeg, num_shap, num_dofs);
  //kew*/
  
  
  // 4. PDE COEFFICIENTS
  
  // there are two choices:
  // 1. consider all terms (16 arrays and 4 vectors)
  // 2. use coeff_vector_indicator to select which terms are non-zero
  // this options have to be taken into account when rewriting coefficients returned 
  // by problem dependent module to coeff array 
  // HERE coeff_vector_indicator is used to calculate pde_coeff_size
  int pde_coeff_size = nr_coeff_arrays*nreq*nreq + nr_coeff_vectors*nreq;
  
  
  // for different types of problems solved, there is a multitude of options:
  // 1. no coefficients at all (e.g. Laplace)
  // 2. coefficients constant for all elements (elasticity and uniform material)
  // 3. coefficients constant for all integration points (elasticity)
  // 4. coefficients different for all elements and integration points
  //    a. for nonlinear problems with coefficients sent to GPU
  //    b. for multi-scale problems??? (data different for each integration point)
  // also possible
  // coefficients are sent to GPU but SM calculations involve not only coefficients
  // but some other data as well (previous solution, etc.)
  // then additional factors have to be taken into account (time for sending
  // solution or degrees of freedom for computing solution, time for computing
  // solution from degrees of freedom, time for computing SM entries given
  // coefficients and the other data, etc.)
  
  // FOR SIMPLE PROBLEMS WE SUGGEST TO USE:
  //      REG_JAC, REG_NOJAC, SHM_JAC or SHM_NOJAC - THESE DO NOT REQUIRE
  //      DIFFERENT COEFFICIENTS FOR EACH ELEMENT AND INTEGRATION POINT
  // FOR COMPLEX PROBLEMS WE SUGGEST TO USE:
  //      REG_GL_DER or SHM_GL_DER WITH READY MADE COEFFICIENTS SENT
  //      TO GPU (A LOT OF DATA!!!!!!)
  
  // different parameters to differentiate between different cases
  int all_el_pde_coeff_size = 0; // size for all elements
  int one_el_pde_coeff_size = 0; // size for one element and all integration point
  int one_int_p_pde_coeff_size = 0; // size for one element and one integration point
  
  // default - not practical: all coefficients at all integration points sent
  one_int_p_pde_coeff_size  = pde_coeff_size;
  
  // special versions
#ifdef LAPLACE
  all_el_pde_coeff_size = 0;
  one_el_pde_coeff_size = 0;
  one_int_p_pde_coeff_size = 1; // one RHS coefficient per integration point
#endif
#ifdef TEST_SCALAR
  all_el_pde_coeff_size = 0;
  one_el_pde_coeff_size = pde_coeff_size; // all coeff are constant over element
  one_int_p_pde_coeff_size  = 0;
#endif
  
  int el_pde_data_size = all_el_pde_coeff_size + 
              Nr_elems*(one_el_pde_coeff_size + ngauss*one_int_p_pde_coeff_size);
  
  // 5. COMPUTED STIFFNESS MATRIX AND LOAD VECTOR
  int one_el_stiff_mat_size = num_dofs*num_dofs;
  int one_el_load_vec_size = num_dofs;
  
/*kbw*/
  printf("\nData structure sizes:\n");
  printf("\tPDE coefficients for each element %d or each integration point %d\n", 
	 one_el_pde_coeff_size, one_int_p_pde_coeff_size);
  printf("\tPDE coefficients for all elements %d and total %d\n", 
	 all_el_pde_coeff_size, el_pde_data_size);
  printf("\tSM for each element %d\n", one_el_stiff_mat_size);
  printf("\tLoad vector for each element %d\n", one_el_load_vec_size);
/*kew*/

  
  
  // !!!!!!!!!!!!!!!***************!!!!!!!!!!!!!!!!!
  // finally fill element input data
  int ielem;
  int packed_bytes = 0;
  int final_position = 0; // for testing packing procedure 
  
  int el_pde_data_bytes = el_pde_data_size * sizeof(SCALAR);
  SCALAR* el_pde_data = (SCALAR *)malloc(el_pde_data_bytes);
  
  //#ifndef OPENCL_HSA
  memset(el_pde_data, 0, el_pde_data_bytes);
  //#endif
  
  int position_coeff = 0;
  for(ielem=0; ielem<Nr_elems; ielem++){
    
    // element ID
    el_id = L_elem_id[ielem];
    
    int el_mate = mmr_el_groupID(mesh_id, el_id);
    double hsize = mmr_el_hsize(mesh_id, el_id, NULL,NULL,NULL);
    
    // checking whether this element has the same data as assumed for this color
    assert( pdeg == apr_get_el_pdeg(field_id, el_id, NULL) );  
    assert( num_shap == apr_get_el_pdeg_numshap(field_id, el_id, &pdeg) ); 
    assert( base == apr_get_base_type(field_id, el_id) );
    
    
    // PDE coefficients
    
#ifdef LAPLACE
    
    // get the size of quadrature data
    int base = apr_get_base_type(field_id, el_id);
    int ngauss;            /* number of gauss points */
    double xg[3000];   	 /* coordinates of gauss points in 3D */
    double wg[1000];       /* gauss weights */
    apr_set_quadr_3D(base, &pdeg, &ngauss, xg, wg);
    
    double geo_dofs[3*MMC_MAXELVNO];  /* coord of nodes of El */
    mmr_el_node_coor(mesh_id, el_id, NULL, geo_dofs);
    
    int ki;
    for(ki=0;ki<ngauss;ki++)  {
      double xcoor[3] = {0,0,0};
      
      apr_elem_calc_3D(2, nreq, &pdeg, base,
		       &xg[3*ki], geo_dofs, NULL,
		       NULL,NULL,NULL,NULL,
		       xcoor,NULL,NULL,NULL,NULL,NULL);
      
      pdr_el_coeff(Problem_id, el_id, el_mate, hsize, pdeg, NULL,
		   NULL, NULL, NULL, NULL,
		   xcoor, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
		   //base_phi_ref, base_dphix, base_dphiy, base_dphiz,
		   //xcoor, uk_val, uk_x, uk_y, uk_z, un_val, un_x, un_y, un_z,
		   NULL, axx, axy, axz, ayx, ayy, ayz, azx, azy, azz,
		   bx, by, bz, tx, ty, tz, cval, NULL, qx, qy, qz, sval);
      
      // for Laplace (Poisson) problem we send RHS coefficients only
      // (stiffness matrix coefficients are all 1.0)
      el_pde_data[position_coeff+ki] = sval[0];
      
/*kbw
	    if(field_id==1 && el_id==13753)
	    {
	    int i,ueq;

		  printf("Problem_id %d, el_id %d, el_mate %d, hsize %lf, pdeg %d\n",
		   Problem_id, el_id, el_mate, hsize, pdeg);

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
    }
    
    position_coeff += ngauss;
    packed_bytes += ngauss*sizeof(SCALAR);
    
#endif // LAPLACE
    
#ifdef TEST_SCALAR
    
    double xcoor_middle[3];
    xcoor_middle[0] = (geo_dofs[0+0]+geo_dofs[3+0]+geo_dofs[6+0]
		       +geo_dofs[9+0]+geo_dofs[12+0]+geo_dofs[15+0])/6.0;
    xcoor_middle[1] = (geo_dofs[0+1]+geo_dofs[3+1]+geo_dofs[6+1]
		       +geo_dofs[9+1]+geo_dofs[12+1]+geo_dofs[15+1])/6.0;
    xcoor_middle[2] = (geo_dofs[0+2]+geo_dofs[3+2]+geo_dofs[6+2]
		       +geo_dofs[9+2]+geo_dofs[12+2]+geo_dofs[15+2])/6.0;
    pdr_el_coeff(Problem_id, el_id, el_mate, hsize, pdeg, NULL, 
		 NULL, NULL, NULL, NULL, 
		 xcoor_middle, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
		 //base_phi_ref, base_dphix, base_dphiy, base_dphiz, 
		 //xcoor, uk_val, uk_x, uk_y, uk_z, un_val, un_x, un_y, un_z, 
		 NULL, axx, axy, axz, ayx, ayy, ayz, azx, azy, azz, 
		 bx, by, bz, tx, ty, tz, cval, NULL, qx, qy, qz, sval); 
    
      
      
/*kbw
    if(field_id==1 && el_id==13753)
      {  
      int i,ueq;

      printf("Problem_id %d, el_id %d, el_mate %d, hsize %lf, pdeg %d\n",
	   Problem_id, el_id, el_mate, hsize, pdeg);

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
	
    
    el_pde_data[position_coeff+0] = axx[0];
    el_pde_data[position_coeff+1] = axy[0];
    el_pde_data[position_coeff+2] = axz[0];
    el_pde_data[position_coeff+3] = ayx[0];
    el_pde_data[position_coeff+4] = ayy[0];
    el_pde_data[position_coeff+5] = ayz[0];
    el_pde_data[position_coeff+6] = azx[0];
    el_pde_data[position_coeff+7] = azy[0];
    el_pde_data[position_coeff+8] = azz[0];
    
    el_pde_data[position_coeff+9] = bx[0];
    el_pde_data[position_coeff+10] = by[0];
    el_pde_data[position_coeff+11] = bz[0];
    
    el_pde_data[position_coeff+12] = tx[0];
    el_pde_data[position_coeff+13] = ty[0];
    el_pde_data[position_coeff+14] = tz[0];
    
    el_pde_data[position_coeff+15] = cval[0];
    
    el_pde_data[position_coeff+16] = qx[0];
    el_pde_data[position_coeff+17] = qy[0];
    el_pde_data[position_coeff+18] = qz[0];
    
    el_pde_data[position_coeff+19] = sval[0];
    position_coeff += 20;
    packed_bytes += 20*sizeof(SCALAR);
    
    assert(pde_coeff_size==20);
    
    
#endif // TEST_SCALAR
  
  } // end loop over elements
  
  *Size_SM_LV_p = num_dofs;
  *Size_global_pde_data_p = all_el_pde_coeff_size;
  *Size_per_element_pde_data_p = one_el_pde_coeff_size;
  *Size_per_int_point_pde_data_p = one_int_p_pde_coeff_size;
  // total size = Size_global_pde_data + Nr_elems * 
  //              (Size_per_element_pde_data + ngauss*Size_per_int_point_pde_data) 
  *El_pde_dat_host_p = el_pde_data;
  
  return(0);
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////        OLD OBSOLETE    ////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

// Master switch: GPU versus PHI - may be controlled by compilation options !!!!!
// MUST BE COMPATIBLE WITH: tmd_ocl/tms_ocl_intf.c and tmd_ocl/tmh_ocl.h 

//Opencl_HSA is for Heterogenous System Architecture with Shared Virtual Memory eg. APU

#ifdef OPENCL_HSA
	#define OPENCL_GPU
#endif

//#define OPENCL_CPU
//#define OPENCL_GPU
//#define OPENCL_PHI


// SWITCH 2: Jacobian on GPU or CPU (CPU version obsolete - but may require two kernel 
//           passes - one for computing Jacobian terms and the second for integration)
//#define JACOBIAN_CALCULATIONS_CPU
#define JACOBIAN_CALCULATIONS_GPU

// SWITCH 2: one_el_one_thread strategy versus one_el_one_workgroup strategy
#define ONE_EL_ONE_THREAD
//#define ONE_EL_ONE_WORKGROUP
//#define ONE_EL_TWO_THREADS



// SWITCH 4: size for  work-group
#ifdef OPENCL_CPU
#define WORK_GROUP_SIZE 8
#elif defined OPENCL_GPU
#define WORK_GROUP_SIZE 64
#elif defined OPENCL_PHI
#define WORK_GROUP_SIZE 16
#endif

// this should always be defined for GPUs
#define COAL_WRITE


/*------------------------------------------------------------
 pdr_create_assemble_stiff_mat_opencl_elem - to create element stiffness matrices
                                 and assemble them to the global SM
// GENERIC PROCEDURE - MANY PROBLEM DEPENDENT OPTIMIZATIONS ARE POSSIBLE
// SEVERAL SUITABLE PLACES ARE INDICATED BY:
// !!!OPT_PDT!!!
------------------------------------------------------------*/
int pdr_create_assemble_stiff_mat_opencl_elem(
  int Problem_id, 
  int Level_id, 
  int Comp_type,         /* in: indicator for the scope of computations: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int* Pdeg_coarse_p, // indicator or value for pdeg on coarse meshes
  int Nr_elems_opencl,
  int* L_int_ent_type,
  int* L_int_ent_id,
  int Max_dofs_int_ent
					      )
{
  int i,j,k;

#ifdef TIME_TEST
    t_begin = time_clock();
#endif

  // choose the platform 
  int platform_index = tmr_ocl_get_current_platform_index();
  
  int device_tmc_type = tmr_ocl_get_current_device_type();

  // choose device_index
  int device_index = tmr_ocl_get_current_device_type(); //tmr_ocl_select_device(platform_index, device_tmc_type);
  
  // OpenCL device characteristics stored in data structure
  tmt_ocl_device_struct device_struct = 
    tmv_ocl_struct.list_of_platforms[platform_index].list_of_devices[device_index];
  double global_mem_bytes = device_struct.global_mem_bytes;
  double global_max_alloc = device_struct.global_max_alloc;
  double shared_mem_bytes = device_struct.shared_mem_bytes;
  double constant_mem_bytes = device_struct.constant_mem_bytes;
  int max_num_comp_units = device_struct.max_num_comp_units;
  int max_work_group_size = device_struct.max_work_group_size;

/*kbw*/
  printf("\n\nStarting OpenCL GPU computations: platform %d, device %d, nr_elems %d\n", 
	 platform_index, device_index, Nr_elems_opencl);
  printf("OpenCL device characteristics: max_num_comp_units %d, max_work_group_size %d\n"
	 , max_num_comp_units, max_work_group_size);
  printf("\tglobal_mem_bytes %.0lf, global_max_alloc %.0lf, shared_mem_bytes %.0lf\n",
	 global_mem_bytes, global_max_alloc, shared_mem_bytes);
/*kew*/
  
  // choose the context
  cl_context context = tmr_ocl_select_context(platform_index, device_index);  
  
  // choose the command queue
  cl_command_queue command_queue = 
    tmr_ocl_select_command_queue(platform_index, device_index);  
  
  // THERE MAY BE SEVERAL KERNELS FOR THE CODE, THE INDEX IN OPENCL DATA 
  // STRUCTURE FOR THE NUMERICAL INTEGRATION KERNEL IS ASSIGNED IN (tmd_ocl/tmh_ocl.h)
  int kernel_index=TMC_OCL_KERNEL_NUM_INT_INDEX; 
  
  // !!!!! FOR THE TIME BEING WE ALWAYS READ AND COMPILE tmr_ocl_num_int_el.cl
  // FROM THE CURRENT DIRECTORY AND ASSIGN DEFAULT VERSION = 0 

  cl_kernel kernel = tmr_ocl_select_kernel(platform_index, device_index, kernel_index); 

  int kernel_version_alg = 0 ; 
  
#define ONE_EL_ONE_THREAD_KERNEL 1
#define ONE_EL_TWO_THREADS_KERNEL 2
#define ONE_EL_ONE_WORKGROUP_KERNEL 3
  
#ifdef ONE_EL_ONE_THREAD
  kernel_version_alg = ONE_EL_ONE_THREAD_KERNEL;
#elif defined(ONE_EL_TWO_THREADS)
  kernel_version_alg = ONE_EL_TWO_THREADS_KERNEL;
#elif defined ONE_EL_ONE_WORKGROUP
  kernel_version_alg = ONE_EL_ONE_WORKGROUP_KERNELP;
#else
  printf("wrong kernel version specified!!!). Exiting.\n");
  exit(-1);
#endif
  
/*kbw*/
  printf("\nSelected: kernel_index %d, version %d (", kernel_index, kernel_version_alg);
#ifdef ONE_EL_ONE_THREAD
  printf("ONE_EL_ONE_THREAD)\n");
#elif defined(ONE_EL_TWO_THREADS)
  printf("ONE_EL_TWO_THREADS)\n");
#elif defined ONE_EL_ONE_WORKGROUP
  printf("ONE_EL_ONE_WORKGROUP)\n");
#else
  printf("wrong kernel version specified!!!). Exiting.\n");
  exit(-1);
#endif
/*kew*/


  if(context == NULL || command_queue == NULL || kernel == NULL){ 
    
    printf("failed to restore kernel for platform %d, device %d, kernel %d\n", 
	   platform_index, device_index, kernel_index);
    printf("context %lu, command queue %lu, kernel %lu\n", 
	   context, command_queue, kernel);
    exit(-1);
  }

  //
  // PREPARE PROBLEM DEPENDENT DATA FOR INTEGRATION
  //

  //int name=pdr_ctrl_i_params(Problem_id,1); changed for string
  char name[300];
  pdr_problem_name(Problem_id, name);
  int field_id = pdr_ctrl_i_params(Problem_id, 3);
  int mesh_id = apr_get_mesh_id(field_id);
  int nreq = apr_get_nreq(field_id);

  if(nreq>PDC_MAXEQ){
    printf("nreq (%d) > PDC_MAXEQ (%d) in pdr_create_assemble_stiff_mat_opencl_elem\n",
	   nreq, PDC_MAXEQ);
    printf("Exiting!\n"); exit(-1);
  }

  // get the active PDE coefficient matrices
  /* pde coefficients */
  static int coeff_ind = 0;
  int coeff_ind_vect[23] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  // 1 - mval, 2 - axx, 3 - axy, 4 - axz, 5 - ayx, 6 - ayy, 7 - ayz, 
  // 8 - azx, 9 - azy, 10 - azz, 11 - bx, 12 - by, 13 - bz
  // 14 - tx, 15 - ty, 16 - tz, 17 - cval
  // 18 - lval, 19 - qx, 20 - qy, 21 - qz, 22 - sval 

  double axx[PDC_MAXEQ*PDC_MAXEQ];
  double axy[PDC_MAXEQ*PDC_MAXEQ];
  double axz[PDC_MAXEQ*PDC_MAXEQ];
  double ayx[PDC_MAXEQ*PDC_MAXEQ];
  double ayy[PDC_MAXEQ*PDC_MAXEQ];
  double ayz[PDC_MAXEQ*PDC_MAXEQ];
  double azx[PDC_MAXEQ*PDC_MAXEQ];
  double azy[PDC_MAXEQ*PDC_MAXEQ];
  double azz[PDC_MAXEQ*PDC_MAXEQ];
  double bx[PDC_MAXEQ*PDC_MAXEQ];
  double by[PDC_MAXEQ*PDC_MAXEQ];
  double bz[PDC_MAXEQ*PDC_MAXEQ];
  double tx[PDC_MAXEQ*PDC_MAXEQ];
  double ty[PDC_MAXEQ*PDC_MAXEQ];
  double tz[PDC_MAXEQ*PDC_MAXEQ];
  double cval[PDC_MAXEQ*PDC_MAXEQ];
  double mval[PDC_MAXEQ*PDC_MAXEQ];
  double qx[PDC_MAXEQ];
  double qy[PDC_MAXEQ];
  double qz[PDC_MAXEQ];
  double sval[PDC_MAXEQ];
  double lval[PDC_MAXEQ];

  pdr_select_el_coeff_vect(Problem_id, coeff_ind_vect); 

  // there are two choices:
  // 1. consider all terms (16 arrays and 4 vectors)
  // 2. use coeff_vector_indicator to select which terms are non-zero
  int coeff_array_indicator[16];
  int nr_coeff_arrays = 0;
  for(i=0; i<16; i++) {
    coeff_array_indicator[i]=coeff_ind_vect[i+2];
    if(coeff_ind_vect[i+2]==1) nr_coeff_arrays++;
  }
  // in create_assemble we combine mval (1) with cval (17)
  // we use this for vectorization of kernels
  if(coeff_ind_vect[1]==1) {
    if(coeff_array_indicator[15]!=1) nr_coeff_arrays++;
    coeff_array_indicator[15]=1;
  }
  int coeff_vector_indicator[4];
  int nr_coeff_vectors = 0; 
  for(i=0; i<4; i++) {
    coeff_vector_indicator[i]=coeff_ind_vect[i+19];
    if(coeff_ind_vect[i+19]==1) nr_coeff_vectors++;
  }
  // in create_assemble we combine lval (18) with sval (22)
  // we use this for vectorization of kernels
  if(coeff_ind_vect[18]==1) {
    if(coeff_vector_indicator[3]!=1) nr_coeff_vectors++;
    coeff_vector_indicator[3]=1;
  }


  
/*kbw*/
  printf("Problem ID: %d, name %s\n", Problem_id, name);
  printf("\nthe number of coefficients arrays %d, indicator = \n",nr_coeff_arrays);
  for(i=0;i<16;i++){
    printf("%2d",coeff_array_indicator[i]);
  }
  printf("\n");
  printf("the number of coefficients vectors %d, indicator = \n",nr_coeff_vectors);
  for(i=0;i<4;i++){
    printf("%2d",coeff_vector_indicator[i]);
  }
  printf("\n");
/*kew*/


  //
  // DIVIDE ELEMENTS INTO GROUPS WITH THE SAME TYPE, PDEG AND THE SAME COLOUR
  //

  // ARTIFICIAL ASSUMPTIONS - COLORING TO BE DONE !!!!!!!!!!!!!!!!!!!
  int nr_colors = 1;


  int icolor;
  for(icolor = 0; icolor<nr_colors; icolor++){

    // nr_elems_per_color is used to allocate memory
    int nr_elems_per_color =  ceil(((double) Nr_elems_opencl)/nr_colors); 
    // nr_elems_per_color may be too large so that for the last color there are less elements
    // nr_elems_this_color is used for actual calculations
    int nr_elems_this_color = nr_elems_per_color;
    if(icolor==nr_colors-1) 
      nr_elems_this_color = Nr_elems_opencl-nr_elems_per_color*(nr_colors-1);

    if(nr_elems_this_color > nr_elems_per_color){
      printf("nr_elems_this_color %d > nr_elems_per_color %d !!!\n",
	     nr_elems_this_color, nr_elems_per_color);
      exit(-1);
    }

    //
    // BASED ON PDEG COMPUTE EXECUTION CHARACTERISTICS
    //

    // choose an example element for a given pdeg and color
    int el_id = L_int_ent_id[nr_elems_per_color*icolor];

    int pdeg = apr_get_el_pdeg(field_id, el_id, NULL);  
    int num_shap = apr_get_el_pdeg_numshap(field_id, el_id, &pdeg); 
    int num_dofs = num_shap * nreq;
    int el_nodes[MMC_MAXELVNO+1];        /* list of nodes of El */
    mmr_el_node_coor(mesh_id, el_id, el_nodes, NULL);
    /* for geometrically (multi)linear elements number of geometrical  */
    /* degrees of freedom is equal to the number of vertices - classical FEM nodes */
    int num_geo_dofs = el_nodes[0];
    int pdeg_single = pdeg;
    if(pdeg>100) {
      pdeg_single = pdeg/100;
      if(pdeg != pdeg_single*100 + pdeg_single){
	printf("wrong pdeg %d ( > 100 but not x0x )\n", pdeg);
	exit(-1); 
      }
    }

/*kbw*/
      printf("problem and element characteristics: nreq %d, pdeg %d, num_shap %d, num_dofs %d!\n",
	     nreq, pdeg, num_shap, num_dofs);
//kew*/


    // SIZES OF DATA STRUCTURES 
    

    // 1. QUADRATURE DATA AND JACOBIAN TERMS

    // get the size of quadrature data
    int base = apr_get_base_type(field_id, el_id);
    int ngauss;            /* number of gauss points */
    double xg[3000];   	 /* coordinates of gauss points in 3D */
    double wg[1000];       /* gauss weights */
    apr_set_quadr_3D(base, &pdeg, &ngauss, xg, wg);

    // we may need quadrature data for the reference element
    int ref_el_quadr_dat_size;
    // but we may need also/instead quadrature related Jacobian data for each element
    int one_el_quadr_dat_size;
    
#ifdef JACOBIAN_CALCULATIONS_GPU
    if(num_dofs==4)
    {
        // for tetrahedrons only weight are sent for reference element
        ref_el_quadr_dat_size = ngauss;
        // we do not need Jacobian terms, but we need geometry DOFs
        one_el_quadr_dat_size = 0;
    }else{
		// for each gauss point its coordinates and weight are sent for reference element
		ref_el_quadr_dat_size = ngauss*4;
		// we do not need Jacobian terms, but we need geometry DOFs
		one_el_quadr_dat_size = 0;
    }
#elif defined JACOBIAN_CALCULATIONS_CPU
    // we do not need base quadrature data (coordinates and weights)
    ref_el_quadr_dat_size = 0;
    // for each gauss point and element - 10 parameters are sent (wg*det, detadx[9])
    one_el_quadr_dat_size = ngauss*10;
#else
    printf("Wrong switch for Jacobian terms in opencl integration. Exiting.");
    exit(-1);
#endif
    
    
    // 2. GEO_DOFS    
    
    // get the size of geometry data for one element - we assume multi-linear elements
    double geo_dofs[3*MMC_MAXELVNO];  /* coord of nodes of El */
    
    int one_el_geo_dat_size;
#ifdef JACOBIAN_CALCULATIONS_GPU
    one_el_geo_dat_size = 3*num_geo_dofs;
#elif defined JACOBIAN_CALCULATIONS_CPU
    one_el_geo_dat_size = 0;
#else
    printf("Wrong switch for Jacobian terms in opencl integration. Exiting.");
    exit(-1);
#endif
    
    
    // 3. SHAPE FUNCTIONS
    
    // space for element shape functions' values and derivatives in global memory
    int ref_el_shape_fun_size; // for the reference element
    int one_el_shape_fun_size; // for each processed element
    
    // we need all shape functions and their derivatives
    // at all integration points for the reference element 
    ref_el_shape_fun_size = 4*num_shap*ngauss;
    // we do not store any data for particular elements
    one_el_shape_fun_size = 0;
    
    
    
    // 4. PDE COEFFICIENTS
    
    // there are two choices:
    // 1. consider all terms (16 arrays and 4 vectors)
    // 2. use coeff_vector_indicator to select which terms are non-zero
    // this options have to be taken into account when rewriting coefficients returned 
    // by problem dependent module to coeff array 
    // HERE coeff_vector_indicator is used to calculate pde_coeff_size
    int pde_coeff_size = nr_coeff_arrays*nreq*nreq + nr_coeff_vectors*nreq;
    
    
    // for different types of problems solved, there is a multitude of options:
    // 1. no coefficients at all (e.g. Laplace)
    // 2. coefficients constant for all elements (elasticity and uniform material)
    // 3. coefficients constant for all integration points (elasticity)
    // 4. coefficients different for all elements and integration points
    //    a. for nonlinear problems with coefficients sent to GPU
    //    b. for multi-scale problems??? (data different for each integration point)
    // also possible
    // coefficients are sent to GPU but SM calculations involve not only coefficients
    // but some other data as well (previous solution, etc.)
    // then additional factors have to be taken into account (time for sending
    // solution or degrees of freedom for computing solution, time for computing
    // solution from degrees of freedom, time for computing SM entries given
    // coefficients and the other data, etc.)

    // FOR SIMPLE PROBLEMS WE SUGGEST TO USE:
    //      REG_JAC, REG_NOJAC, SHM_JAC or SHM_NOJAC - THESE DO NOT REQUIRE
    //      DIFFERENT COEFFICIENTS FOR EACH ELEMENT AND INTEGRATION POINT
    // FOR COMPLEX PROBLEMS WE SUGGEST TO USE:
    //      REG_GL_DER or SHM_GL_DER WITH READY MADE COEFFICIENTS SENT
    //      TO GPU (A LOT OF DATA!!!!!!)

    // different parameters to differentiate between different cases
    int all_el_pde_coeff_size = 0; // size for all elements
    int one_el_pde_coeff_size = 0; // size for one element and all integration point
    int one_int_p_pde_coeff_size = 0; // size for one element and one integration point

    // default - not practical: all coefficients at all integration points sent
    one_int_p_pde_coeff_size  = pde_coeff_size;

      // special versions
#ifdef LAPLACE
    all_el_pde_coeff_size = 0;
    one_el_pde_coeff_size = 0;
    one_int_p_pde_coeff_size = 1; // one RHS coefficient per integration point
#endif
#ifdef TEST_SCALAR
    all_el_pde_coeff_size = 0;
    one_el_pde_coeff_size = pde_coeff_size; // all coeff are constant over element
    one_int_p_pde_coeff_size  = 0;
#endif


    // 5. COMPUTED STIFFNESS MATRIX AND LOAD VECTOR
    int one_el_stiff_mat_size = num_dofs*num_dofs;
    int one_el_load_vec_size = num_dofs;

/*kbw*/
    printf("\nAssumed data structure sizes:\n");
    printf("\tQuadrature data: global %d, for each element %d\n", 
	   ref_el_quadr_dat_size, one_el_quadr_dat_size);
    printf("\tGeo dofs for each element: %d\n", one_el_geo_dat_size);
    printf("\tShape functions and derivatives: reference el. %d, each el. %d\n",
	   ref_el_shape_fun_size, one_el_shape_fun_size);
    printf("\tPDE coefficients for each element %d or each integration point %d\n", 
	   one_el_pde_coeff_size, one_int_p_pde_coeff_size);
    printf("\tSM for each element %d\n", one_el_stiff_mat_size);
    printf("\tLoad vector for each element %d\n", one_el_load_vec_size);
/*kew*/

#define NR_EXEC_PARAMS 16  // size of array with execution parameters
    // here: the smallest work-group for reading data is selected
    // exec_params are read from global to shared memory and used when needed
    // if shared memory resources are scarce this can be reduced

    // ACTUAL GLOBAL!!! MEMORY CALCULATIONS
    int global_memory_req_in = NR_EXEC_PARAMS; // size of array with execution parameters
    int global_memory_req_one_el_in = 0;
    int global_memory_req_out = 0; // size of array with execution parameters
    int global_memory_req_one_el_out = 0;

    // 1. QUADRATURE DATA AND JACOBIAN TERMS
    //  - coordinates and weights for JAC - for the reference element
    //  - Jacobian terms for NOJAC - for all elements 
    global_memory_req_in += ref_el_quadr_dat_size;
    global_memory_req_one_el_in += one_el_quadr_dat_size;

    // 2. GEO_DOFS - geo_dofs for JAC, nothing for NOJAC 
    global_memory_req_one_el_in += one_el_geo_dat_size;

    // 3. SHAPE FUNCTIONS
    // for the reference element
    //  - all functions and derivatives at all integration points - JAC and NOJAC
    global_memory_req_in += ref_el_shape_fun_size;
    // for each considered element - nothing for the time being
    global_memory_req_one_el_in += one_el_shape_fun_size;

    // 4. PDE COEFFICIENTS
    // global coefficients - the same for all elements
    global_memory_req_in += all_el_pde_coeff_size;
    // COEFFICIENTS DIFFERENT FOR EACH ELEMENT BUT THE SAME FOR ALL INTEGRATION POINTS
    global_memory_req_one_el_in += one_el_pde_coeff_size;
    // COEFFICIENTS DIFFERENT FOR EACH ELEMENT AND EACH INTEGRATION POINT ?????????????
    global_memory_req_one_el_in += ngauss*one_int_p_pde_coeff_size;


    // 5. COMPUTED STIFFNESS MATRIX AND LOAD VECTOR
    global_memory_req_one_el_out += one_el_stiff_mat_size + one_el_load_vec_size;

    int global_memory_req = global_memory_req_in + global_memory_req_out;
    int global_memory_req_one_el = global_memory_req_one_el_in + global_memory_req_one_el_out;

/*kbw*/
  printf("\nGlobal memory requirements (before computing number of elements considered):\n");
  printf("\treference element %d, for each element %d\n",
	 global_memory_req, global_memory_req_one_el);
/*kew*/

 
// ASSUMED WORK_GROUP_SIZE
  int work_group_size = WORK_GROUP_SIZE;


// GLOBAL GPU MEMORY CONSIDERATIONS FOR SETTING THE NUMBER OF ELEMENTS PER KERNEL

    // nr_elems_per_kernel is related to the size of input and output data for kernel
    // therefore:
    // 1. nr_elems_per_kernel is optimized to make the best use of global GPU memory
    // 2. nr_elems_per_kernel is optimized to speed up transfer times for input and output
    // 3. nr_elems_per_kernel may be optimized for overlapping communication with computations


    // we take global_max_alloc as the first indicator for nr_elems_per_kernel
    int global_limiter;
    if(global_memory_req_one_el_in > global_memory_req_one_el_out) {
      global_limiter = global_memory_req_one_el_in;
    }
    else{
      global_limiter = global_memory_req_one_el_out;
    }
    int nr_elems_per_kernel = global_max_alloc / (global_limiter*sizeof(SCALAR));

    // CHECK 1
    if( (global_memory_req + global_memory_req_one_el*nr_elems_per_kernel)*sizeof(SCALAR)
	> global_mem_bytes) {
      nr_elems_per_kernel = 
	(global_mem_bytes/sizeof(SCALAR) - global_memory_req)/global_memory_req_one_el;
    }


    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // AT THAT MOMENT WE CAN CHANGE nr_elems_per_kernel SO THAT THERE ARE 
    // TODO: SEVERAL KERNELS WORKING CONCURRENTLY
    // moreover we can assume that not only integration kernels are working concurrently
    // but also solver kernels are working concurrently
    // even more: we can limit nr_elems_per_kernel so we can use many times buffers
    // allocated once in global GPU memory
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


// THE SETTING OF NR_WORK_GROUPS, NR_THREADS AND NR_ELEMS_PER_WORK_GROUP
// THIS IS BASED MORE ON EXPERIMENTATION THAN ON THEORETICAL GROUNDS

// how many workgroups will be enough to make GPU cores (comp_units) busy?
#ifdef OPENCL_CPU
#define NR_WORK_GROUPS_PER_COMP_UNIT  1 // ????
#elif defined OPENCL_GPU
#define NR_WORK_GROUPS_PER_COMP_UNIT  8 // ????
#elif defined OPENCL_PHI
#define NR_WORK_GROUPS_PER_COMP_UNIT  1 // ????
#endif

    // usually for GPUs it is good to maximize the number of work_groups and threads
    int nr_work_groups = NR_WORK_GROUPS_PER_COMP_UNIT * max_num_comp_units;
    int nr_elems_per_work_group;

    // CHECK 2 nr_elems_per_kernel can be larger than nr_elems_per_color but not too much
    // i.e. we can allocate more memory for kernel execution but not much more
    // than necessary
      
    if( kernel_version_alg==ONE_EL_ONE_THREAD_KERNEL) {
      
      if(nr_elems_per_kernel > nr_elems_per_color + nr_work_groups*work_group_size) {
	// if nr_elems_per_kernel > nr_elems_per_color + nr_work_groups*work_group_size
	// we set nr_elems_per_work_group in such a way that it is the multiple of work_group_size
	// but nr_elems_per_color < nr_elems_per_work_group*nr_work_groups 
	
	nr_elems_per_work_group = 
	  ceil(((double)nr_elems_per_color) / nr_work_groups / work_group_size) * work_group_size;
	
      }
      else{
	
	// if nr_elems_per_kernel <= nr_elems_per_color + nr_work_groups*work_group_size
	// we set nr_elems_per_work_group in such a way that it is the multiple of work_group_size
	// but nr_elems_per_kernel >= nr_elems_per_work_group*nr_work_groups 
	
	int temp = nr_elems_per_kernel / nr_work_groups;
	temp = temp / work_group_size;
	nr_elems_per_work_group = temp * work_group_size;
	
      }
      
    }
    else{
      // the strategy 1 work_group per 1 element (but many elements per work_group)
      
      if(nr_elems_per_kernel > nr_elems_per_color + nr_work_groups) {
	
	// if nr_elems_per_kernel > nr_elems_per_color + nr_work_groups
	// nr_elems_per_color < nr_elems_per_work_group*nr_work_groups 
	nr_elems_per_work_group =  ceil(((double)nr_elems_per_color) / nr_work_groups);
	
      }
      else{
	// if nr_elems_per_kernel <= nr_elems_per_color + nr_work_groups
	// nr_elems_per_kernel >= nr_elems_per_work_group*nr_work_groups 
	
	nr_elems_per_work_group = nr_elems_per_kernel / nr_work_groups;
	
      }
      
      
    }

    // we adjust nr_elems_per_kernel (we always decrease it !!!)
    // nr_elems_per_kernel = nr_elems_per_work_group*nr_work_groups 
    nr_elems_per_kernel = nr_work_groups * nr_elems_per_work_group;


/*kbw
    printf("\nInitial settings: nr_elems_per_work_group %d = nr_elems_per_kernel %d / nr_work_groups %d\n",
	   nr_elems_per_work_group, nr_elems_per_kernel, nr_work_groups);
//kew*/

#ifdef ONE_EL_ONE_THREAD

#define MIN_NR_ELEMS_PER_WORK_GROUP WORK_GROUP_SIZE // ??? 

#endif

#ifdef ONE_EL_ONE_WORKGROUP

#define MIN_NR_ELEMS_PER_WORK_GROUP 1 // ??? 

#endif

    // CHECK 3
    // we adjust the number of work_groups so that either each work_group operates
    // 1. on one element
    // 2. or on several elements ( >MIN_NR_ELEMS_PER_WORK_GROUP )
    if(nr_elems_per_work_group < MIN_NR_ELEMS_PER_WORK_GROUP) {
      if( kernel_version_alg==ONE_EL_ONE_WORKGROUP_KERNEL) {
	nr_work_groups = nr_elems_per_kernel;
	nr_elems_per_work_group = 1;
      }
      else{
	printf("Too low number of elements per workgroup %d \n", nr_elems_per_work_group );
	printf("Check parameters of execution. Exiting!!!\n");
	exit(-1);
      }
    }

/*kbw
      printf("\tADJUSTED: nr_work_groups %d, nr_elems_per_work_group %d\n", 
	     nr_work_groups, nr_elems_per_work_group);
//kew*/
    
    // we check whether index within stiff_mat in global memory do not exceed MAXINT
    if((double)one_el_stiff_mat_size*nr_elems_per_kernel > (double)INT_MAX){
      printf("Too big size of stiffness matrices in global device memory\n");
      printf("Index exceeds MAX_INT. Exiting. Decrease nr_elems_per_kernel = %d.\n",
	     nr_elems_per_kernel);
      exit(-1);
    }


    // finally we decide how many threads will perform calculations
    int nr_threads = nr_work_groups * work_group_size;
    int nr_elems_per_thread = 0;

    if( kernel_version_alg==ONE_EL_ONE_THREAD_KERNEL) {
      nr_elems_per_thread = nr_elems_per_work_group / work_group_size;
    }
    else{
      nr_elems_per_thread = nr_elems_per_work_group;
    }

/*kbw*/
    printf("\nAFTER GLOBAL MEMORY COSIDERATIONS (with size %.3lf[MB] and max alloc %.3lf[MB])\n",
	   global_mem_bytes*1.0e-6, global_max_alloc*1.0e-6);
    printf("\tnr_elems_per_kernel %d, nr_elems_per_work_group %d, nr_elems_per_thread %d\n",
	   nr_elems_per_kernel,  nr_elems_per_work_group, nr_elems_per_thread);
    printf("\tnr_work_groups %d, work_group_size %d, nr_threads %d\n",
	   nr_work_groups, work_group_size, nr_threads );
	   //printf("\tnumber of elements for kernel %d, el_data_in size %.3lf[MB], el_data_out size %.3lf [MB]\n",
	   //	   nr_elems_per_kernel, (double)el_data_in_bytes_max*1.0e-6, 
	   //	   (double)nr_elems_per_kernel*one_el_stiff_mat_size*sizeof(SCALAR)*1.0e-6);
//kew*/

    if( kernel_version_alg==ONE_EL_ONE_THREAD_KERNEL) {
      if(nr_elems_per_work_group != nr_elems_per_thread*work_group_size){
	printf("nr_elems_per_work_group %d != %d nr_elems_per_thread*work_group_size\n",
	       nr_elems_per_work_group, nr_elems_per_thread*work_group_size );
	printf("Exiting!\n");
	exit(-1);
      }
    }



    
// WE ALLOCATE OPENCL BUFFERS FOR nr_elems_per_kernel ELEMENTS
// this should be good for all kernel invocations and all colors !!!
// (actual kernel calculations will be for nr_elems)

    cl_int retval = CL_SUCCESS;

// there are general data, integration points data, data for reference elements etc.
// with size independent of the number of elements
// element data are in el_data_in buffer with size el_data_in_bytes
// output data are in el_data_out buffer with size el_data_out_bytes
// the final set of buffers are workspaces on the device (in global or shared mameory)


// arrays for kernel parameters (execution parameters, integration points data, shape
// functions and their derivatives for reference element)

    // there are NR_EXEC_PARAMS execution parameters (ALLWAYS CHECK!!!)

    const size_t execution_parameters_dev_bytes = NR_EXEC_PARAMS * sizeof(int);

    cl_mem execution_parameters_dev;

    execution_parameters_dev = clCreateBuffer(context, CL_MEM_READ_ONLY,
    					execution_parameters_dev_bytes, NULL, NULL);
    retval |= clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &execution_parameters_dev);

/*
    int* execution_parameters_host=(int*)clSVMAlloc(context,CL_MEM_READ_ONLY|CL_MEM_SVM_FINE_GRAIN_BUFFER,execution_parameters_dev_bytes,NULL);

    retval |= clSetKernelArgSVMPointer(kernel, 0, execution_parameters_host);
*/
/*kbw
    printf("\nAllocated buffer on device for exec. params %u with size: %u\n",
	   execution_parameters_dev, execution_parameters_dev_bytes);
//kew*/
 
    if (retval != CL_SUCCESS) {
      printf("Failed to Set the kernel argument %d.\n", 0);
      exit(0);
    }

    const size_t gauss_dat_dev_bytes = ref_el_quadr_dat_size * sizeof(SCALAR);
    cl_mem gauss_dat_dev = NULL;

    // if integration points data are not necessary
    if(gauss_dat_dev_bytes == 0){

      gauss_dat_dev = clCreateBuffer(context, CL_MEM_READ_ONLY, 1, NULL, NULL);
      retval |= clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &gauss_dat_dev);

    }
    else{

      gauss_dat_dev = clCreateBuffer(context, CL_MEM_READ_ONLY,
      				   gauss_dat_dev_bytes, NULL, NULL);
      retval |= clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &gauss_dat_dev);

    }

/*kbw
    printf("\nAllocated buffer on device for Gauss data with size: %u\n",
	   gauss_dat_dev_bytes);
//kew*/

    if (retval != CL_SUCCESS) {
      printf("Failed to Set the kernel argument %d.\n", 1);
      exit(0);
    }


    const size_t shape_fun_dev_bytes  = ref_el_shape_fun_size * sizeof(SCALAR);
    cl_mem shape_fun_dev;

    // if reference element shape functions data are not necessary
    if(shape_fun_dev_bytes == 0){

      shape_fun_dev = clCreateBuffer(context, CL_MEM_READ_ONLY, 1, NULL, NULL);
      retval |= clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) NULL);

    }
    else{

      shape_fun_dev = clCreateBuffer(context, CL_MEM_READ_ONLY,
      				     shape_fun_dev_bytes, NULL, NULL);
      retval |= clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &shape_fun_dev);

    }

    if (retval != CL_SUCCESS) {
      printf("Failed to Set the kernel argument %d.\n", 2);
      exit(0);
    }

/*kbw
    printf("\nAllocated buffer on device for shape functions with size: %u\n",
	   shape_fun_dev_bytes);
//kew*/

    
// el_data_in and el_data_out large arrays with element data
// since we aim at NVIDIA graphics cards we follow rules from
// OpenCL Best Practices Guide (NVIDIA, 14.02.2011) for using
// the page-locked memory for faster transfers

    // size of input buffer for element data
    const size_t el_data_in_bytes =
      global_memory_req_one_el_in * nr_elems_per_kernel * sizeof(SCALAR);

// For CPU better result is when we don't use pinned buffer
#ifndef OPENCL_CPU
	#define PINNED
#endif

#ifdef OPENCL_HSA

    SCALAR* el_data_in = (SCALAR *) clSVMAlloc(context,CL_MEM_READ_ONLY|CL_MEM_SVM_FINE_GRAIN_BUFFER,el_data_in_bytes,NULL);
    if(el_data_in==NULL)
        	printf("Error 536536356 in allocating buffer");

    retval |= clSetKernelArgSVMPointer(kernel, 3, el_data_in);

#else



	#ifdef PINNED
		cl_mem cmPinnedBufIn = NULL; // input buffer mapped to host input buffer
		cmPinnedBufIn = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_ALLOC_HOST_PTR,
					   el_data_in_bytes, NULL, NULL);
	#endif

		cl_mem cmDevBufIn = NULL; // device input buffer allocated in card's memory
		cmDevBufIn = clCreateBuffer(context, CL_MEM_READ_ONLY,
					el_data_in_bytes, NULL, NULL);

		SCALAR* el_data_in = NULL; // host input buffer (standard array)

	#ifdef PINNED
		el_data_in = (SCALAR*)clEnqueueMapBuffer(command_queue, cmPinnedBufIn, CL_TRUE,
								 CL_MAP_WRITE, 0, el_data_in_bytes, 0,
								 NULL, NULL, NULL);
	#else
		el_data_in = (SCALAR *)malloc(el_data_in_bytes);
	#endif

		retval |= clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &cmDevBufIn);

#endif //ifdef HSA

		// writing input data to global GPU memory
		//clEnqueueWriteBuffer(command_queue, cmDevBufIn, CL_FALSE, 0,
		//			 el_data_in_bytes, el_data_in, 0, NULL, NULL);


	const size_t el_data_out_bytes =
	global_memory_req_one_el_out * nr_elems_per_kernel * sizeof(SCALAR);


#ifdef OPENCL_HSA

    SCALAR* el_data_out = (SCALAR *) clSVMAlloc(context,CL_MEM_WRITE_ONLY|CL_MEM_SVM_FINE_GRAIN_BUFFER,el_data_out_bytes,NULL);
    if(el_data_out==NULL)
    	printf("Error 536536356 in allocating buffer");

    retval |= clSetKernelArgSVMPointer(kernel, 4, el_data_out);

#else

	#ifdef PINNED
		cl_mem cmPinnedBufOut = NULL; // output buffer mapped to host output buffer
		cmPinnedBufOut = clCreateBuffer(context, CL_MEM_WRITE_ONLY | CL_MEM_ALLOC_HOST_PTR,
							el_data_out_bytes, NULL, NULL);
	#endif

		cl_mem cmDevBufOut = NULL; // device output buffer allocated in card's memory
		cmDevBufOut = clCreateBuffer(context, CL_MEM_WRITE_ONLY,
						 el_data_out_bytes, NULL, NULL);

		SCALAR* el_data_out = NULL; // host output buffer (standard array)

	#ifdef PINNED
		el_data_out = (SCALAR*)clEnqueueMapBuffer(command_queue, cmPinnedBufOut, CL_TRUE,
							   CL_MAP_READ, 0, el_data_out_bytes, 0,
							   NULL, NULL, NULL);
	#else
		el_data_out = (SCALAR *)malloc(el_data_out_bytes);
	#endif

		retval |= clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &cmDevBufOut);

#endif

    // reading output data from global GPU memory
    // clEnqueueReadBuffer(command_queue, cmDevBufOut, CL_TRUE, 0,
    //		  	   el_data_out_bytes, el_data_out, 0, NULL, NULL);

    
    if (retval != CL_SUCCESS) {
      printf("Failed to Set the kernel argument %d.\n", 4);
      exit(0);
    }
  
/*kbw*/
    printf("\nAllocated buffer on device for %d elements data:\n", nr_elems_per_kernel );
    printf("\tel_data_in bytes %.3lf[MB], el_data_out bytes %.3lf [MB]\n",
	   (double)el_data_in_bytes*1.0e-6, (double)el_data_out_bytes*1.0e-6);
//kew*/

#ifdef TUNING

    //system("more header.csv |wc -l >num.txt");
    unsigned int line_count=1;
    FILE *numf;
    numf=fopen("num","r");
    if(numf)
    {
    	line_count = 0;
    	fclose(numf);
    	system("rm num");
    }

    printf("line_count=%d\n",line_count);
	if(line_count==0)
	{
		fprintf(headf,"el_data_in [MB],el_data_out [MB],");
	}
    fprintf(resf,"%.3lf,%.3lf,",
	   (double)el_data_in_bytes*1.0e-6, (double)el_data_out_bytes*1.0e-6);
#endif


#ifdef TIME_TEST
    //clFinish(command_queue);
    t_end = time_clock();
    printf("\nEXECUTION TIME: initial settings on CPU and creating buffers on GPU %lf\n",
	   t_end-t_begin);
    total_time += t_end-t_begin; 
#endif


// AFTER ALLOCATING DATA BUFFERS ON GPU CARD IN A LOOP OVER KERNEL CALLS
// BUFFERS ARE FILLED WITH DATA, DATA ARE SEND TO GPU, CALCULATIONS PERFORMED
// AND OUTPUT DATA TRANSFERRED BACK TO HOST COMPUTER MEMORY

    // the number of kernel calls per single color depends on the number of elements
    int nr_kernel_calls=1;
    if(nr_elems_per_kernel < nr_elems_this_color){

      nr_kernel_calls = ceil(((double) nr_elems_this_color) / nr_elems_per_kernel);

    }

    int nr_elems_per_kercall = nr_elems_per_kernel;


    if( kernel_version_alg==ONE_EL_ONE_THREAD_KERNEL) {
      // if the last kernel call gets too low number of elements
      if(nr_elems_this_color - (nr_kernel_calls-1)*nr_elems_per_kernel < 
	 nr_work_groups*work_group_size){
	
	// we decrease the number of elements per kernel call
	nr_elems_per_kercall = nr_elems_per_kernel-nr_work_groups*work_group_size;
	
      }
    }

/*kbw
      printf("\nnr_colors %d, nr_elems_per_color %d, nr_elems_this_color %d\n",
	     nr_colors, nr_elems_per_color, nr_elems_this_color);
//kew*/


// for each color we set the ID for the first and the last element
    int first_elem_color=nr_elems_per_color*icolor;
    int last_elem_color=nr_elems_per_color*(icolor+1) - 1;
    if(icolor==nr_colors-1) last_elem_color=Nr_elems_opencl-1; // IDs from 0 to Nr_elems_opencl-1

/*kbw
      printf("\nColor %d, first_elem_color %d, last_elem_color %d\n",
	     icolor, first_elem_color, last_elem_color);
//kew*/

    // in a loop over several kernel calls for elements of the same color
    int ikercall;
    for(ikercall=0; ikercall<nr_kernel_calls; ikercall++){
      
      /* int first_elem_kercall = first_elem_color + ikercall*nr_elems_per_kernel; */
      /* int last_elem_kercall = first_elem_kercall + nr_elems_per_kernel - 1; */
      /* if(last_elem_kercall > last_elem_color) last_elem_kercall = last_elem_color; */
      /* int nr_elems_this_kercall = last_elem_kercall - first_elem_kercall + 1; */
      
      int first_elem_kercall = first_elem_color + ikercall*nr_elems_per_kercall;
      int last_elem_kercall = first_elem_kercall + nr_elems_per_kercall - 1;
      if(last_elem_kercall > last_elem_color) last_elem_kercall = last_elem_color;
      int nr_elems_this_kercall = last_elem_kercall - first_elem_kercall + 1;

      /*kbw*/
      printf("\nKernel call %d, first_elem_kercall %d, last_elem_kercall %d\n",
	     ikercall, first_elem_kercall, last_elem_kercall);
      printf("nr_elems_this_kercall %d, nr_elems_per_kernel %d\n",
	     nr_elems_this_kercall, nr_elems_per_kernel);
      //kew*/
      
      
      // for the last time we adjust the number of elements per work-group
      // to take into account the number of elements for this kernel call
      if( kernel_version_alg==ONE_EL_ONE_THREAD_KERNEL) {
	
	// we set nr_elems_per_work_group in such a way that it is the multiple of work_group_size
	// but nr_elems_this_kercall < nr_elems_per_work_group*nr_work_groups 
	nr_elems_per_work_group = 
	  ceil(((double)nr_elems_this_kercall) / nr_work_groups / work_group_size) * work_group_size;
	
	
	
      }
      else{
	// the strategy 1 work_group per 1 element (but many elements per work_group)
	
	
	// if nr_elems_per_kernel > Nr_elems + nr_work_groups
	// Nr_elems < nr_elems_per_work_group*nr_work_groups 
	nr_elems_per_work_group =  ceil(((double)nr_elems_this_kercall) / nr_work_groups);
	
	
      
      }
      
      
      if(nr_elems_per_work_group < MIN_NR_ELEMS_PER_WORK_GROUP) {
	if( kernel_version_alg==ONE_EL_ONE_WORKGROUP_KERNEL) {
	  nr_work_groups = nr_elems_this_kercall;
	  nr_elems_per_work_group = 1;
	}
	else{
	  printf("Too low number of elements per workgroup %d \n", nr_elems_per_work_group );
	  printf("Check parameters of execution. Exiting!!!\n");
	  exit(-1);
	}
	
      }
      

      if( kernel_version_alg==ONE_EL_ONE_THREAD_KERNEL) {
	nr_elems_per_thread = nr_elems_per_work_group / work_group_size;
      }
      else{
	nr_elems_per_thread = nr_elems_per_work_group;
      }
      
/*kbw*/
      printf("ADJUSTED TO THIS KERNEL CALL:\n");
      printf("\tnr_elems: nr_elems_per_work_group %d (%d), nr_elems_per_thread %d (%d)(%d)\n",
	     nr_elems_per_work_group, nr_elems_per_work_group*nr_work_groups,
	     nr_elems_per_thread, nr_elems_per_thread*nr_threads,
	     nr_elems_per_thread*work_group_size*nr_work_groups);
      printf("\tnr_work_groups %d, work_group_size %d, nr_threads %d\n",
	     nr_work_groups, work_group_size, nr_threads );
//kew*/
#ifdef TUNING
      if(line_count==0)
      {
      	fprintf(headf,"nr_elems_per_work_group,nr_elems,nr_elems_per_thread,nr_work_groups,work_group_size,nr_threads,");
      }
      fprintf(resf,"%d,%d,%d,",
    		     nr_elems_per_work_group, nr_elems_per_work_group*nr_work_groups, nr_elems_per_thread);
      fprintf(resf,"%d,%d,%d,",
    		     nr_work_groups, work_group_size, nr_threads);
#endif

      if( kernel_version_alg==ONE_EL_ONE_THREAD_KERNEL) {
	if(nr_elems_per_work_group != nr_elems_per_thread*work_group_size){
	  printf("nr_elems_per_work_group %d != %d nr_elems_per_thread*work_group_size\n",
		 nr_elems_per_work_group, nr_elems_per_thread*work_group_size );
	  printf("Exiting!\n");
	  exit(-1);
	}
	if(nr_elems_per_work_group*nr_work_groups != nr_elems_per_thread*nr_threads){
	  printf("nr_elems_per_work_group*nr_work_groups %d != %d nr_elems_per_thread*nr_threads\n",
		 nr_elems_per_work_group*nr_work_groups, nr_elems_per_thread*nr_threads );
	  printf("Exiting!\n");
	  exit(-1);
	}
	if(nr_elems_per_work_group*nr_work_groups < nr_elems_this_kercall){
	  printf("nr_elems_per_work_group*nr_work_groups %d < %d nr_elems_this_kercall\n",
		 nr_elems_per_work_group*nr_work_groups, nr_elems_this_kercall );
	  printf("Exiting!\n");
	  exit(-1);
	}

      }
      
      
      if(nr_elems_per_work_group*nr_work_groups > nr_elems_per_kernel) {
	printf("nr_elems_per_work_group*nr_work_groups %d > %d nr_elems_per_kernel\n", 
	       nr_elems_per_work_group*nr_work_groups, nr_elems_per_kernel );
	printf("Check parameters of execution. Exiting!!!\n");
	exit(-1);
      }
      
   
#ifdef TIME_TEST
    clFinish(command_queue);
    t_begin = time_clock();
#endif

// FILL MEMORY OBJECTS WITH DATA AND SEND INPUT DATA TO GPU MEMORY

      // !!!!!!!!!!!!!!!***************!!!!!!!!!!!!!!!!!
      // first  fill buffers for assumed and computed execution parameters
      // there are NR_EXEC_PARAMS execution parameters (ALLWAYS CHECK!!!)
      int execution_parameters_host[NR_EXEC_PARAMS];

      execution_parameters_host[0] = nr_elems_per_thread;
      execution_parameters_host[1] = nr_elems_this_kercall;
      //execution_parameters_host[1] = nr_parts_of_stiff_mat;
      //execution_parameters_host[2] = nr_blocks_per_thread;
      //execution_parameters_host[3] = nreq;
      //execution_parameters_host[4] = pdeg_single;
      //execution_parameters_host[5] = num_shap;
      //execution_parameters_host[6] = base;
      //execution_parameters_host[7] = ngauss;
      //execution_parameters_host[8] = nr_work_groups;
      //execution_parameters_host[9] = nr_elems_per_kernel; 
      //execution_parameters_host[9] = work_group_size; - taken from OpenCL in kernels
 
      // send data to device - non-blocking call (CL_FALSE)
      clEnqueueWriteBuffer(command_queue, execution_parameters_dev, CL_TRUE, 0,
      			   execution_parameters_dev_bytes, execution_parameters_host,
      			   0, NULL, NULL);

      // !!!!!!!!!!!!!!!***************!!!!!!!!!!!!!!!!!
      // fill and send buffers with integration data when necessary
      if(gauss_dat_dev_bytes != 0){
	
#define MAX_SIZE_ARRAY_GAUSS 1344
	// 1344 - maximal value of arrGaussSize (for p=707), (for 909 is 2920)
	SCALAR gauss_dat_host[MAX_SIZE_ARRAY_GAUSS];
// This is the version when we send only necessary variables but we do not care about sending
	if(num_dofs==4)
	{
		for(i=0; i<ngauss; i++){
		  gauss_dat_host[i] = wg[i];
		}
	}
	else
	{
		for(i=0; i<ngauss; i++){
		  gauss_dat_host[4*i] = xg[3*i];
		  gauss_dat_host[4*i+1] = xg[3*i+1];
		  gauss_dat_host[4*i+2] = xg[3*i+2];
		  gauss_dat_host[4*i+3] = wg[i];
		}
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

	// send data to device - blocking version to ensure gauss_dat_host is still valid
	clEnqueueWriteBuffer(command_queue, gauss_dat_dev, CL_TRUE, 0,
			     gauss_dat_dev_bytes, gauss_dat_host, 0, NULL, NULL);
	
      }
      
      
      // shape functions and derivatives are computed here but used also later on
      double base_phi_ref[APC_MAXELVD];    /* basis functions */
      double base_dphix_ref[APC_MAXELVD];  /* x-derivatives of basis function */
      double base_dphiy_ref[APC_MAXELVD];  /* y-derivatives of basis function */
      double base_dphiz_ref[APC_MAXELVD];  /* z-derivatives of basis function */
      
      // !!!!!!!!!!!!!!!***************!!!!!!!!!!!!!!!!!
      // fill and send buffers with shape function values for the reference element
      if(shape_fun_dev_bytes != 0){
	
	int shape_fun_host_bytes = shape_fun_dev_bytes;
	SCALAR *shape_fun_host = (SCALAR*)malloc(shape_fun_host_bytes);
	
	
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
	for(i=ngauss*num_shap-5;i<ngauss*num_shap;i++){
	
	printf("shape_fun_host[%d] = %f\n", i, shape_fun_host[i]);
	//printf("shape_fun_host[%d] = %lf\n", i, shape_fun_host[i]);
	
	}
//kew*/

	// send data to device - CL_TRUE i.e. blocking to free shape_fun_host
	clEnqueueWriteBuffer(command_queue, shape_fun_dev, CL_TRUE, 0,
			     shape_fun_dev_bytes, shape_fun_host, 0, NULL, NULL);
	
	free(shape_fun_host);
    
      }
      
      
      
      // !!!!!!!!!!!!!!!***************!!!!!!!!!!!!!!!!!
      // finally fill element input data
      int ielem;
      int packed_bytes = 0;
      int final_position = 0; // for testing packing procedure 

//#ifndef OPENCL_HSA
      memset(el_data_in, 0, el_data_in_bytes);
//#endif
      
#ifdef ONE_EL_ONE_WORKGROUP
      // for one_el_one_workgroup strategy we keep geo_dofs and coeffs for each element together
      int position_all = 0; // actual position in el_data_in and also the number of packed items
#endif
      
#ifdef ONE_EL_ONE_THREAD
      // for one_el_one_thread strategy we keep geo_dofs and coeffs for each element separately
      // i.e. geo_dofs for all elements and next coeffs for all elements
      int position_geo_dofs = 0;
      int position_coeff = nr_elems_this_kercall*num_geo_dofs*3;
#endif
      
      for(ielem=first_elem_kercall; ielem<=last_elem_kercall; ielem++){
	
	// element ID
	el_id = L_int_ent_id[ielem];
	
	int el_mate = mmr_el_groupID(mesh_id, el_id);
	double hsize = mmr_el_hsize(mesh_id, el_id, NULL,NULL,NULL);
	
	// checking whether this element has the same data as assumed for this color
	assert( pdeg == apr_get_el_pdeg(field_id, el_id, NULL) );  
	assert( num_shap == apr_get_el_pdeg_numshap(field_id, el_id, &pdeg) ); 
	assert( base == apr_get_base_type(field_id, el_id) );
	
	// element geo_dofs for computing Jacobian terms on GPU
#ifdef  JACOBIAN_CALCULATIONS_GPU
	
	// IDs of element vertices (nodes) and their coordinates as geo_dofs
	mmr_el_node_coor(mesh_id, el_id, NULL, geo_dofs);
	
	
#ifdef ONE_EL_ONE_THREAD
	for(i=0;i<num_geo_dofs;i++){
	  el_data_in[position_geo_dofs] = geo_dofs[3*i];
	  el_data_in[position_geo_dofs+1] = geo_dofs[3*i+1];
	  el_data_in[position_geo_dofs+2] = geo_dofs[3*i+2];
	  position_geo_dofs += 3;
	}
#endif
	
#ifdef ONE_EL_ONE_WORKGROUP
	for(i=0;i<num_geo_dofs;i++){
	  el_data_in[position_all] = geo_dofs[3*i];
	  el_data_in[position_all+1] = geo_dofs[3*i+1];
	  el_data_in[position_all+2] = geo_dofs[3*i+2];
	  position_all += 3;
	}
#endif
	
	packed_bytes += num_geo_dofs*3*sizeof(SCALAR);
	
#endif

	
	// PDE coefficients
	
#ifdef LAPLACE

	int ki;
	for(ki=0;ki<ngauss;ki++)
	  {
	    double xcoor[3] = {0,0,0};
	    
	    apr_elem_calc_3D(2, nreq, &pdeg, base,
			     &xg[3*ki], geo_dofs, NULL,
			     NULL,NULL,NULL,NULL,
			     xcoor,NULL,NULL,NULL,NULL,NULL);
	    
	    pdr_el_coeff(Problem_id, el_id, el_mate, hsize, pdeg, NULL,
			 NULL, NULL, NULL, NULL,
			 xcoor, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
			 //base_phi_ref, base_dphix, base_dphiy, base_dphiz,
			 //xcoor, uk_val, uk_x, uk_y, uk_z, un_val, un_x, un_y, un_z,
			 NULL, axx, axy, axz, ayx, ayy, ayz, azx, azy, azz,
			 bx, by, bz, tx, ty, tz, cval, NULL, qx, qy, qz, sval);
	    
	    // for Laplace (Poisson) problem we send RHS coefficients only
	    // (stiffness matrix coefficients are all 1.0)
	    el_data_in[position_coeff+ki] = sval[0];
	    
/*kbw
	    if(field_id==1 && el_id==13753)
	    {
	    int i,ueq;

		  printf("Problem_id %d, el_id %d, el_mate %d, hsize %lf, pdeg %d\n",
		   Problem_id, el_id, el_mate, hsize, pdeg);

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
	  }
	
	position_coeff += ngauss;
	packed_bytes += ngauss*sizeof(SCALAR);
	
#endif // LAPLACE
	
#ifdef TEST_SCALAR
	
	double xcoor_middle[3];
	xcoor_middle[0] = (geo_dofs[0+0]+geo_dofs[3+0]+geo_dofs[6+0]
			   +geo_dofs[9+0]+geo_dofs[12+0]+geo_dofs[15+0])/6.0;
	xcoor_middle[1] = (geo_dofs[0+1]+geo_dofs[3+1]+geo_dofs[6+1]
			   +geo_dofs[9+1]+geo_dofs[12+1]+geo_dofs[15+1])/6.0;
	xcoor_middle[2] = (geo_dofs[0+2]+geo_dofs[3+2]+geo_dofs[6+2]
			   +geo_dofs[9+2]+geo_dofs[12+2]+geo_dofs[15+2])/6.0;
	pdr_el_coeff(Problem_id, el_id, el_mate, hsize, pdeg, NULL, 
		     NULL, NULL, NULL, NULL, 
		     xcoor_middle, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
		     //base_phi_ref, base_dphix, base_dphiy, base_dphiz, 
		     //xcoor, uk_val, uk_x, uk_y, uk_z, un_val, un_x, un_y, un_z, 
		     NULL, axx, axy, axz, ayx, ayy, ayz, azx, azy, azz, 
		     bx, by, bz, tx, ty, tz, cval, NULL, qx, qy, qz, sval); 
	
      
      
/*kbw
    if(field_id==1 && el_id==13753)
      {  
      int i,ueq;

      printf("Problem_id %d, el_id %d, el_mate %d, hsize %lf, pdeg %d\n",
	   Problem_id, el_id, el_mate, hsize, pdeg);

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
	
	    
	el_data_in[position_coeff+0] = axx[0];
	el_data_in[position_coeff+1] = axy[0];
	el_data_in[position_coeff+2] = axz[0];
	el_data_in[position_coeff+3] = ayx[0];
	el_data_in[position_coeff+4] = ayy[0];
	el_data_in[position_coeff+5] = ayz[0];
	el_data_in[position_coeff+6] = azx[0];
	el_data_in[position_coeff+7] = azy[0];
	el_data_in[position_coeff+8] = azz[0];
	
	el_data_in[position_coeff+9] = bx[0];
	el_data_in[position_coeff+10] = by[0];
	el_data_in[position_coeff+11] = bz[0];
	
	el_data_in[position_coeff+12] = tx[0];
	el_data_in[position_coeff+13] = ty[0];
	el_data_in[position_coeff+14] = tz[0];
	
	el_data_in[position_coeff+15] = cval[0];
	
	el_data_in[position_coeff+16] = qx[0];
	el_data_in[position_coeff+17] = qy[0];
	el_data_in[position_coeff+18] = qz[0];
	
	el_data_in[position_coeff+19] = sval[0];
	position_coeff += 20;
	packed_bytes += 20*sizeof(SCALAR);

	assert(pde_coeff_size==20);
	
	
#endif // TEST_SCALAR
	
	
#ifdef ONE_EL_ONE_THREAD
	int position_jac_data = position_coeff; // position for jac_data
#endif
      
#ifdef ONE_EL_ONE_WORKGROUP
	int position_jac_data = position_all; // position for jac_data
#endif
	
	
#ifdef JACOBIAN_CALCULATIONS_CPU
	
	// IDs of element vertices (nodes) and their coordinates as geo_dofs
	mmr_el_node_coor(mesh_id, el_id, NULL, geo_dofs);
	
	double det, dxdeta[9], detadx[9]; /* jacobian, jacobian matrix */
	/*  and its inverse; for geometrical transformation */
	
	if(num_geo_dofs==4){ // geometrically linear tetrahedral element
	  
	  const int X=0;
	  const int Y=1;
	  const int Z=2;
	  const double one_sixth = 0.16666666666666666666666666666667;
	  
	  dxdeta[0]=geo_dofs[3+X]-geo_dofs[0+X];
	  dxdeta[3]=geo_dofs[3+Y]-geo_dofs[0+Y];
	  dxdeta[6]=geo_dofs[3+Z]-geo_dofs[0+Z];
	  
	  dxdeta[1]=geo_dofs[6+X]-geo_dofs[0+X];
	  dxdeta[4]=geo_dofs[6+Y]-geo_dofs[0+Y];
	  dxdeta[7]=geo_dofs[6+Z]-geo_dofs[0+Z];
	  
	  dxdeta[2]=geo_dofs[9+X]-geo_dofs[0+X];
	  dxdeta[5]=geo_dofs[9+Y]-geo_dofs[0+Y];
	  dxdeta[8]=geo_dofs[9+Z]-geo_dofs[0+Z];
	  
	  /* Jacobian |J| and inverse of the Jacobian matrix*/
	  det = utr_mat3_inv(dxdeta, detadx) * one_sixth; // *1/6 fot tetrahedrons
	  
	}
	
	for (ki=0;ki<ngauss;ki++) {
	  
	  // local coordinates of integration point
	  double Eta[3] = { xg[3*ki], xg[3*ki+1], xg[3*ki+2]};
	  
	  if(num_geo_dofs==6){ // geometrically bi-linear prismatic element
	    
	    double geo_dphix[8]; 	/* derivatives of geometry shape functions */
	    double geo_dphiy[8];  /* derivatives of geometry shape functions */
	    double geo_dphiz[8];  /* derivatives of geometry shape functions */
	    
	    /* local derivatives of geometrical shape functions*/
	    geo_dphix[0] = -(1.0-Eta[2])/2.0;
	    geo_dphix[1] =  (1.0-Eta[2])/2.0;
	    geo_dphix[2] =  0.0;
	    geo_dphix[3] = -(1.0+Eta[2])/2.0;
	    geo_dphix[4] =  (1.0+Eta[2])/2.0;
	    geo_dphix[5] =  0.0;
	    geo_dphiy[0] = -(1.0-Eta[2])/2.0;
	    geo_dphiy[1] =  0.0;
	    geo_dphiy[2] =  (1.0-Eta[2])/2.0;
	    geo_dphiy[3] = -(1.0+Eta[2])/2.0;
	    geo_dphiy[4] =  0.0;
	    geo_dphiy[5] =  (1.0+Eta[2])/2.0;
	    geo_dphiz[0] = -(1.0-Eta[0]-Eta[1])/2.0;
	    geo_dphiz[1] = -Eta[0]/2.0;
	    geo_dphiz[2] = -Eta[1]/2.0;
	    geo_dphiz[3] =  (1.0-Eta[0]-Eta[1])/2.0;
	    geo_dphiz[4] =  Eta[0]/2.0;
	    geo_dphiz[5] =  Eta[1]/2.0;
	    /* Jacobian matrix J */
	    dxdeta[0] = 0.0; dxdeta[1] = 0.0; dxdeta[2] = 0.0;
	    dxdeta[3] = 0.0; dxdeta[4] = 0.0; dxdeta[5] = 0.0;
	    dxdeta[6] = 0.0; dxdeta[7] = 0.0; dxdeta[8] = 0.0;
	    for(i=0;i<num_geo_dofs;i++){
	      dxdeta[0] += geo_dofs[3*i]  *geo_dphix[i];
	      dxdeta[1] += geo_dofs[3*i]  *geo_dphiy[i];
	      dxdeta[2] += geo_dofs[3*i]  *geo_dphiz[i];
	      dxdeta[3] += geo_dofs[3*i+1]*geo_dphix[i];
	      dxdeta[4] += geo_dofs[3*i+1]*geo_dphiy[i];
	      dxdeta[5] += geo_dofs[3*i+1]*geo_dphiz[i];
	      dxdeta[6] += geo_dofs[3*i+2]*geo_dphix[i];
	      dxdeta[7] += geo_dofs[3*i+2]*geo_dphiy[i];
	      dxdeta[8] += geo_dofs[3*i+2]*geo_dphiz[i];
	    }
	    
	    /* Jacobian |J| and inverse of the Jacobian matrix*/
	    det = utr_mat3_inv(dxdeta,detadx);
	    
	  }
	  
	  // 10 values for each integration point (wg*det, detadx[9])
	  for(i=0;i<9;i++){
	    el_data_in[position_jac_data] = detadx[i];
	    position_jac_data++;
	  }
	  el_data_in[position_jac_data] = det * wg[ki];
	  position_jac_data++;
	  packed_bytes += 10*sizeof(SCALAR);
	  
/*kbw
	    if(ielem==0){
	      //if(ielem==last_elem_kercall){

	      printf("Eta %lf, %lf, %lf, wg %lf, vol %lf\n", 
		     Eta[0],Eta[1],Eta[2], wg[ki], det * wg[ki]);
	      for(i=0;i<10;i++){
		printf("detadx[%d] = %lf\n", i, detadx[i]);
		printf("dxdeta[%d] = %lf\n", i, dxdeta[i]);
	      }
	      for(i=position-10;i<position-10+5;i++){
		
		printf("el_data_in[%d] = %lf\n", i, el_data_in[i]);
		
	      }
	    for(i=position-5;i<position;i++){
	      printf("el_data_in[%d] = %lf\n", i, el_data_in[i]);
	      //printf("el_data_in[%d] = %lf\n", i, el_data_in[i]);
	      
	    }
	  }
//kew*/
      	
	}
      
      
#endif // JACOBIAN_CALCULATIONS_CPU
     
	final_position = position_jac_data;  
	       
      } // end for all elements in input (nr_elems_this_kercall)

      assert(packed_bytes == final_position*sizeof(SCALAR));
      
      
      // in create_assemble we have to send to GPU coefficients for all elements
      // processed by kernel and for all integration points
      // WE REWRITE COEFFICIENT MATRICESS TO THE FORM THAT IS SUITABLE FOR THREADS
      // EXECUTING KERNELS AND ACCESSING GLOBAL AND LOCAL MEMORY
      
// !!!OPT_PDT!!!

      // for many problems the coefficients can be sent not as full matrices but
      // in the form of problem specific parameters for which matrices are calculated
      // by kernels
      
      
/*kbw
      for(i=0;i<50;i++){
	
	printf("el_data_in[%d] = %f\n", i, el_data_in[i]);
	//printf("el_data_in[%d] = %lf\n", i, el_data_in[i]);
	
      }
      for(i=position_coeff-5;i<position_coeff;i++){
	printf("el_data_in[%d] = %f\n", i, el_data_in[i]);
	//printf("el_data_in[%d] = %lf\n", i, el_data_in[i]);
	
      }
//kew*/

#ifdef TIME_TEST
      clFinish(command_queue);
      t_end = time_clock();
      printf("\nEXECUTION TIME: filling buffers and sending to GPU memory %lf\n", 
	     t_end-t_begin);
      total_time += t_end-t_begin; 
      t_begin = time_clock();
#endif

      
      assert(packed_bytes<=el_data_in_bytes);
      // writing input data to global GPU memory

#ifndef OPENCL_HSA
      clEnqueueWriteBuffer(command_queue, cmDevBufIn, CL_FALSE, 0,
      			   el_data_in_bytes, el_data_in, 0, NULL, NULL);
#endif
      
#ifdef TIME_TEST
      clFinish(command_queue);
      t_end = time_clock();
      printf("\nEXECUTION TIME: sending el_data_in to GPU memory %lf (speed %lf [GB/s])\n", 
	     t_end-t_begin, el_data_in_bytes*1.0e-9/(t_end-t_begin));
#ifdef TUNING
      if(line_count==0)
	{
	fprintf(headf,"sending el_data_in to GPU memory,[GB/s],");
	}
      fprintf(resf,"%lf,%lf,",
	     t_end-t_begin, el_data_in_bytes*1.0e-9/(t_end-t_begin));
#endif
      total_time += t_end-t_begin; 
      t_begin = time_clock();
#endif
      
      // FINALLY EXECUTE THE KERNEL
      // set kernel invocation parameters 
      size_t globalWorkSize[1] = { 0 };
      size_t localWorkSize[1] = { 0 };
      globalWorkSize[0] = nr_threads ;
      localWorkSize[0] = work_group_size ;
      
      cl_event ndrEvt;
      
      // Queue the kernel up for execution
      retval = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL,
      				      globalWorkSize, localWorkSize,
      				      0, NULL,  &ndrEvt);
      
      
#ifdef TIME_TEST
      clWaitForEvents(1, &ndrEvt);
      clFinish(command_queue);

      t_end = time_clock();
      
      // Get kernel profiling info 
      cl_ulong startTime;
      cl_ulong endTime;
      clGetEventProfilingInfo(ndrEvt,
			      CL_PROFILING_COMMAND_START,
			      sizeof(cl_ulong),
			      &startTime,
			      0);
      clGetEventProfilingInfo(ndrEvt,
			      CL_PROFILING_COMMAND_END,
			      sizeof(cl_ulong),
			      &endTime,
			      0);
      double time_internal = ((double)endTime - (double)startTime)*1.0e-9;
      
      double kernel_execution_time=0.0;
      kernel_execution_time = t_end-t_begin;
      printf("EXECUTION TIME: executing kernel: %lf (internal %lf)\n",
	     kernel_execution_time, time_internal);
#ifdef TUNING
      if(line_count==0)
    	{
    	fprintf(headf,"executing kernel,internal,");
    	}
      fprintf(resf,"%lf,%lf,",
      	     kernel_execution_time, time_internal);
#endif

      total_time += t_end-t_begin; 
      t_begin = time_clock();
#endif
      
    
     //!!!!!!!!!!!!!******************TESTING INPUT/OUTPUT to GPU MEMORY

//#define TESTING_GPU_TRANSFER

#ifdef TESTING_GPU_TRANSFER 

      memset(execution_parameters_host, 0, NR_EXEC_PARAMS*sizeof(int));
      clEnqueueReadBuffer(command_queue, execution_parameters_dev, CL_TRUE, 0,
      			  execution_parameters_dev_bytes, execution_parameters_host,
      			  0, NULL, NULL);
      for(i=0;i<NR_EXEC_PARAMS;i++){
	
	printf("execution_parameters_host[%d] = %d\n", i, execution_parameters_host[i]);
	
      }
      
      if(gauss_dat_dev_bytes != 0){
      
	SCALAR gauss_dat_host[MAX_SIZE_ARRAY_GAUSS];
	memset(gauss_dat_host, 0, gauss_dat_dev_bytes);
	clEnqueueReadBuffer(command_queue, gauss_dat_dev, CL_TRUE, 0,
			    gauss_dat_dev_bytes, gauss_dat_host, 0, NULL, NULL);
	
	for(i=0;i<5;i++){
	  
	  //printf("gauss_dat_host[%d] = %f\n", i, gauss_dat_host[i]);
	  printf("gauss_dat_host[%d] = %lf\n", i, gauss_dat_host[i]);
	  
	}
	for(i=ngauss-5;i<ngauss;i++){
	  
	  //printf("gauss_dat_host[%d] = %f\n", i, gauss_dat_host[i]);
	  printf("gauss_dat_host[%d] = %lf\n", i, gauss_dat_host[i]);
	  
	}
      }     

      if(shape_fun_dev_bytes!=0){

	int shape_fun_host_bytes = shape_fun_dev_bytes;
	SCALAR *shape_fun_host = (SCALAR*)malloc(shape_fun_host_bytes);
	memset(shape_fun_host, 0, shape_fun_dev_bytes);
	clEnqueueReadBuffer(command_queue, shape_fun_dev, CL_TRUE, 0,
			    shape_fun_dev_bytes, shape_fun_host, 0, NULL, NULL);
	
	for(i=0;i<5;i++){
	  
	  //printf("shape_fun_host[%d] = %f\n", i, shape_fun_host[i]);
	  printf("shape_fun_host[%d] = %lf\n", i, shape_fun_host[i]);
	  
	}
	for(i=num_shap-5;i<num_shap;i++){
	  //printf("shape_fun_host[%d] = %f\n", i, shape_fun_host[i]);
	  printf("shape_fun_host[%d] = %lf\n", i, shape_fun_host[i]);
	  
	}
	free(  shape_fun_host);
	
      }

      memset(el_data_in, 0, el_data_in_bytes);
      // writing input data to global GPU memory
      clEnqueueReadBuffer(command_queue, cmDevBufIn, CL_TRUE, 0,
      			   el_data_in_bytes, el_data_in, 0, NULL, NULL);

	for(i=0;i<5;i++){
	  
	  //printf("el_data_in[%d] = %f\n", i, el_data_in[i]);
	  printf("el_data_in[%d] = %lf\n", i, el_data_in[i]);
	  
	}
	for(i=final_position-5;i<final_position;i++){
	  //printf("el_data_in[%d] = %f\n", i, el_data_in[i]);
	  printf("el_data_in[%d] = %lf\n", i, el_data_in[i]);
	  
	}

	clFinish(command_queue);
#endif
     //!!!!!!!!!!!!!******************THE END OF TESTING INPUT/OUTPUT to GPU MEMORY


#ifdef TIME_TEST
      t_begin = time_clock();
#endif

#ifndef OPENCL_HSA

      memset(el_data_out, 0, el_data_out_bytes);

      // transfer output data
      clEnqueueReadBuffer(command_queue, cmDevBufOut, CL_TRUE, 0,
      			  el_data_out_bytes, el_data_out,
      			  0, NULL, NULL);
#endif

/*kbw
      for(i=0;i<100;i++){

      	printf("el_data_out[%d] = %lf\n", i, el_data_out[i]);

      }
      getchar();
//kew*/
    
#ifdef TIME_TEST
      clFinish(command_queue);
      t_end = time_clock();
      printf("EXECUTION TIME: copying output buffer: %lf (total with kernel %lf)\n",
	     t_end-t_begin, t_end-t_begin+kernel_execution_time);
      
#ifdef TUNING
      if(line_count==0)
	{
    	  fprintf(headf,"copying output buffer,[GB/s],");
	}
      fprintf(resf,"%lf,%lf,",
	      t_end-t_begin, el_data_out_bytes*1.0e-9/(t_end-t_begin));
#endif
      
      total_time += t_end-t_begin; 
      t_begin = time_clock();
#endif
      
#ifdef TUNING
#ifdef COUNT_OPER
      fprintf(resf,"%.0lf, %.0lf, %.0lf\n",el_data_out[0],el_data_out[1],el_data_out[2]);
      fclose(resf);
      if(line_count==0)
	{
	  fprintf(headf,"Nr Oper,Nr acces,Nr_global_access\n");
	  fclose(headf);
	  //system("cat header.csv result.csv >result.csv");
	}
      //exit(0);
#else
      fprintf(resf,"\n");
      fclose(resf);
      if(line_count==0)
	{
	  fprintf(headf,"\n");
	  fclose(headf);
	  //system("cat header.csv result.csv >result.csv");
	}
      //exit(0);
#endif
#endif
      
      //#define NO_TESTING_CORRECTNESS
      //#define NO_COPYING
#define TESTING_CORRECTNESS
#define COPYING
      
      int nr_parts_of_stiff_mat = 1;
      int nr_blocks_per_thread = 1;
      int group_id;
      
#ifdef TUNING
      FILE *cw;
      cw=fopen("coal_write","r");
      if(cw)
	{
	  SCALAR* el_data_out_tmp =(SCALAR*) malloc(el_data_out_bytes);
	  memcpy(el_data_out_tmp,el_data_out,el_data_out_bytes);
	  int wg,thr,elem;
	  int offset_gpu,offset_cpu, elem_index;
	  
	  for(wg=0;wg<nr_work_groups;wg++)
	    {
	      for(ielem=0;ielem<(nr_elems_per_work_group/work_group_size);ielem++)
		{
		  
		  for(thr=0;thr<work_group_size;thr++)
		    {
		      elem_index=wg*nr_elems_per_work_group+ielem*work_group_size+thr;
		      offset_cpu=elem_index*(num_shap*num_shap+num_shap);
		      offset_gpu=(elem_index-thr)*(num_shap*num_shap+num_shap);
		      for(i=0;i<num_shap*num_shap;i++)
			{
			  el_data_out[offset_cpu+i]=el_data_out_tmp[offset_gpu+i*work_group_size+thr];
			}
		      
		      for(i=0;i<num_shap;i++)
			{
			  el_data_out[offset_cpu+num_shap*num_shap+i]=el_data_out_tmp[offset_gpu+(num_shap*num_shap+i)*work_group_size+thr];
			}
		    }
		}
	    }
	}
      system("rm coal_write");
#else
#ifdef COAL_WRITE
      SCALAR* el_data_out_tmp =(SCALAR*) malloc(el_data_out_bytes);
      memcpy(el_data_out_tmp,el_data_out,el_data_out_bytes);
      int wg,thr,elem;
      int offset_gpu,offset_cpu, elem_index;
      
      for(wg=0;wg<nr_work_groups;wg++)
	{
	  for(ielem=0;ielem<(nr_elems_per_work_group/work_group_size);ielem++)
	    {
	      
	      for(thr=0;thr<work_group_size;thr++)
		{
		  elem_index=wg*nr_elems_per_work_group+ielem*work_group_size+thr;
		  offset_cpu=elem_index*(num_shap*num_shap+num_shap);
		  offset_gpu=(elem_index-thr)*(num_shap*num_shap+num_shap);
		  for(i=0;i<num_shap*num_shap;i++)
		    {
		      el_data_out[offset_cpu+i]=el_data_out_tmp[offset_gpu+i*work_group_size+thr];
		    }
		  
		  for(i=0;i<num_shap;i++)
		    {
		      el_data_out[offset_cpu+num_shap*num_shap+i]=el_data_out_tmp[offset_gpu+(num_shap*num_shap+i)*work_group_size+thr];
		    }
		}
	    }
	}

 /*kbw
      for(i=0;i<100;i++){

      	printf("el_data_out[%d] = %lf\n", i, el_data_out[i]);

      }
      getchar();
//kew*/
   

#endif
#endif
      

#ifdef TESTING_CORRECTNESS


#pragma omp parallel for default(none) firstprivate(nr_work_groups, nr_elems_per_work_group, nr_parts_of_stiff_mat) firstprivate(num_shap, work_group_size, kernel_index, nreq, num_dofs, el_data_out, first_elem_kercall, last_elem_kercall, Max_dofs_int_ent, Problem_id, L_int_ent_type, L_int_ent_id, Comp_type, field_id, pdeg, Level_id, nr_blocks_per_thread, kernel_version_alg  )
/* #pragma omp parallel for default(none) firstprivate(nr_work_groups, nr_elems_per_work_group, nr_parts_of_stiff_mat, nr_iter_within_part, num_shap, work_group_size, kernel_index, nreq, num_dofs, stiff_mat, el_data_out) shared(first_elem_kercall) */
    for(group_id=0;group_id<nr_work_groups;group_id++){
      int ielem_CPU = first_elem_kercall+nr_elems_per_work_group*group_id;
      double *stiff_mat = (double *)malloc(num_dofs*num_dofs*sizeof(double));
      double *rhs_vect = (double *)malloc(num_dofs*sizeof(double));
      int ielem; int i;
      for(ielem=0;ielem<nr_elems_per_work_group;ielem++){
       if(ielem_CPU<=last_elem_kercall){

	int intent = ielem_CPU;
	int l_dof_ent_id[PDC_MAX_DOF_PER_INT], l_dof_ent_nrdof[PDC_MAX_DOF_PER_INT];
	int l_dof_ent_posglob[PDC_MAX_DOF_PER_INT];
	int l_dof_ent_type[PDC_MAX_DOF_PER_INT];
	int nr_dof_ent = PDC_MAX_DOF_PER_INT;
	int nrdofs_int_ent = Max_dofs_int_ent;
	char rewrite;
	// FOR REWRITING STIFF_MAT AND LOAD_VEC THEY ARE NOT COPUTED
	// only data necessary for asembling are obtained
	/* pdr_comp_stiff_mat(Problem_id, L_int_ent_type[intent], */
	/* 		   L_int_ent_id[intent], PDC_NO_COMP, NULL, */
	/* 		   &nr_dof_ent,l_dof_ent_type,l_dof_ent_id,l_dof_ent_nrdof, */
	/* 		   &nrdofs_int_ent, NULL, NULL, &rewrite); */
	

	// FOR TESTING CORRECTNESS STIFF_MAT AND LOAD_VEC ARE COMPUTED
      /* initialize the matrices to zero */
      for(i=0;i<num_dofs*num_dofs;i++) stiff_mat[i]=0.0;
  
      /* initialize the vector to zero */
      for(i=0;i<num_dofs;i++) rhs_vect[i]=0.0;

      
      	pdr_comp_stiff_mat(Problem_id, L_int_ent_type[intent],
      			   L_int_ent_id[intent], PDC_COMP_BOTH, NULL,
      			   &nr_dof_ent,l_dof_ent_type,l_dof_ent_id,l_dof_ent_nrdof,
      			   &nrdofs_int_ent, stiff_mat, rhs_vect, &rewrite);
      
	

	int ipart;
	int nr_entries=0;
        for(ipart=0;ipart<nr_parts_of_stiff_mat;ipart++){  
	  int iiter;
	  int nr_iter_within_part = nr_blocks_per_thread;
	  int last_iter = nr_iter_within_part;
	  if(kernel_version_alg == -1 ) {
	    if(ipart == nr_parts_of_stiff_mat-1){
	      last_iter = ceil(((double)(num_shap*num_shap-ipart*nr_iter_within_part*work_group_size))/work_group_size);
	      //printf("ipart %d, nr_parts_of_stiff_mat %d, last_iter %d, num_shap2 %d\n",
	      //	     ipart , nr_parts_of_stiff_mat, last_iter, num_shap*num_shap);
	    }
	    else{
	      last_iter = nr_iter_within_part;
	    }
	  }
	  for(iiter = 0; iiter < last_iter; iiter++ ){
	    //int aux_offset = nreq*nreq*work_group_size*(iiter+nr_iter_within_part*(ipart+nr_parts_of_stiff_mat*(ielem+nr_elems_per_work_group*group_id)));
	    int local_offset = iiter*work_group_size*(nreq*nreq);
	    int aux_offset = (group_id*nr_elems_per_work_group+ielem)
	                                *(num_shap*num_shap+num_shap)*nreq*nreq+
	      + ipart*nr_iter_within_part*work_group_size*nreq*nreq + local_offset;
	    int j;
	    int stride;
	    int last_entry;
	    if(kernel_version_alg==ONE_EL_ONE_WORKGROUP_KERNEL ) {
	      if(ipart==nr_parts_of_stiff_mat-1){
		j = num_shap*num_shap -  ipart*work_group_size;
	      }
	      else{
		j = work_group_size;
	      }
	      stride = j;
	      last_entry = j;
	    }
	    else{
	      if(ipart == nr_parts_of_stiff_mat-1 && iiter==last_iter-1){
		j = num_shap*num_shap - work_group_size*(iiter+nr_iter_within_part*ipart);
		//printf("ipart %d, nr_parts_of_stiff_mat %d, iiter %d, last_iter %d, num_shap2 %d, j %d\n",
		//      ipart , nr_parts_of_stiff_mat, iiter, last_iter, num_shap*num_shap, j);
	      }
	      else{
		j = work_group_size;
	      }
	      stride=j;
	      last_entry = j;
	    }
	    int ientry;
	    for(ientry=0;ientry<last_entry;ientry++){
	      int block_offset = work_group_size*(iiter+nr_iter_within_part*ipart)+ientry;
	      int jdofs = block_offset/num_shap;
	      int idofs = (block_offset-jdofs*num_shap); // jdofs = (block_offset%num_shap);
	      int ieq2;
	      for(ieq2=0;ieq2<nreq*nreq;ieq2++){ // block_size = nreq*nreq
		int jeq = ieq2/nreq;
		int ieq = ieq2 - jeq*nreq; // ieq = ieq2%nreq;
		int index_GPU = aux_offset + ieq2*stride + ientry;
		int index_CPU = jdofs*nreq*num_dofs+idofs*nreq+jeq*num_dofs+ieq;
/*kbw*/
		double tol = 1.e-6;
		if((fabs(stiff_mat[index_CPU]) <  tol &&
		     fabs(el_data_out[index_GPU])> tol)
		    || (fabs(stiff_mat[index_CPU]) >  tol &&
		        fabs(stiff_mat[index_CPU]-el_data_out[index_GPU])/
		        fabs(stiff_mat[index_CPU]) >  tol))
{
#ifdef TUNING
			printf("ERROR!!\n");
			exit(-1);
#endif
		  printf("aux_offset %d, stride %d, block_offset %d\n",
			 aux_offset, stride, block_offset);
		  printf("group_id %d, ielem %d, ipart %d, iiter %d, ientry %d, jdofs %d, idofs %d, jeq %d, ieq %d\n",
			 group_id, ielem, ipart, iiter, ientry, jdofs, idofs, jeq, ieq);
		  printf("index GPU %d,\tindex CPU %d,\tvalue GPU %12.6lf,\tvalue CPU %12.6lf\n",
			 index_GPU, index_CPU, el_data_out[index_GPU],
			 stiff_mat[index_CPU]);
		  getchar();
		}
//kew*/
		stiff_mat[index_CPU] = el_data_out[index_GPU];
		nr_entries++;
	      }
	    }
	    int idofs;
	    for(idofs=0;idofs<num_shap;idofs++){
	      int ieq2;
	      for(ieq2=0;ieq2<nreq;ieq2++){
		int index_GPU_v = aux_offset + num_dofs*num_dofs + idofs*nreq + ieq2;
		int index_CPU_v = idofs*nreq + ieq2;
/*kbw*/
		double tol = 1.e-6;
		if((fabs(rhs_vect[index_CPU_v]) <  tol && 
		    fabs(el_data_out[index_GPU_v])> tol) 
		   || (fabs(rhs_vect[index_CPU_v]) >  tol && 
		       fabs(rhs_vect[index_CPU_v]-el_data_out[index_GPU_v])/
		       fabs(rhs_vect[index_CPU_v]) >  tol)){
#ifdef TUNING
			printf("ERROR!!\n");
			exit(-1);
#endif
		  printf("rhs: aux_offset %d, group_id %d, ielem %d (%d), ipart %d, iiter %d, idofs %d, ieq %d\n",
			 aux_offset, group_id, ielem, L_int_ent_id[intent], ipart, iiter, idofs, ieq2);
		  printf("index GPU_v %d,\tindex CPU_v %d,\tvalue GPU_v %12.6lf,\tvalue CPU_v %12.6lf\n",
			 index_GPU_v, index_CPU_v, el_data_out[index_GPU_v],
			 rhs_vect[index_CPU_v]);
		  getchar();
		}
//kew*/
		rhs_vect[index_CPU_v] = el_data_out[index_GPU_v];
	      }
	    }

	  
	  }
	}
	    
/*kbw
#pragma omp critical(printing)
{
    if(L_int_ent_id[intent]!=-1) {
      printf("pdr_create_assemble: before assemble: Solver_id %d, level_id %d, sol_typ %d\n", 
	     Problem_id, level_id, Comp_type);
      int ibl,jbl,pli,plj,nri,nrj,nrdof,jaux;
      printf("ient %d, int_ent_id %d, nr_dof_ent %d\n", 
	     intent, L_int_ent_id[intent], nrdofbl);
      pli = 0; nrdof=0;
      for(ibl=0;ibl<nrdofbl; ibl++) nrdof+=l_bl_nrdofs[ibl];
      for(ibl=0;ibl<nrdofbl; ibl++){
	printf("bl_id %d, bl_nrdof %d\n",
	  l_bl_id[ibl],l_bl_nrdofs[ibl]);
	nri = l_bl_nrdofs[ibl];
	plj=0;
	for(jbl=0;jbl<nrdofbl;jbl++){
	  printf("stiff_mat (blocks %d:%d)\n",jbl,ibl);
	  nrj = l_bl_nrdofs[jbl];
	  for(i=0;i<nri;i++){
   	    jaux = plj+(pli+i)*nrdof;
	    for(j=0;j<nrj;j++){
	      printf("%20.15lf",stiff_mat[jaux+j]);
	    }
	    printf("\n");
	  }
	  plj += nrj;
	}
	printf("Rhs_vect:\n");
	for(i=0;i<nri;i++){
	  printf("%20.15lf",rhs_vect[pli+i]);
	}
	printf("\n");
	pli += nri;    
      }
      getchar();
    }
      }
/*kew*/

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


	/* change the option compute SM and RHSV to rewrite SM and RHSV */
	int Comp_sm = Comp_type;
	if(Comp_sm!=PDC_NO_COMP) Comp_sm += 3;
	apr_get_stiff_mat_data(field_id, L_int_ent_id[intent], Comp_sm, 'N',
			       pdeg, 0, &nr_dof_ent, l_dof_ent_type,
			       l_dof_ent_id, l_dof_ent_nrdof,
			       &nrdofs_int_ent, stiff_mat, rhs_vect);

	
#pragma omp critical(assembling)
    {
	
	pdr_assemble_local_stiff_mat(Problem_id, Level_id, Comp_type,
				     nr_dof_ent, l_dof_ent_type,
				     l_dof_ent_id,l_dof_ent_nrdof, 
				     stiff_mat, rhs_vect, &rewrite);
    }	
	
	ielem_CPU++;

	if(nr_entries!=num_shap*num_shap*nreq*nreq){
	  printf("Error in indexing stiff_mat: nr_entries %d != %d\n",
		 nr_entries, num_shap*num_shap*nreq*nreq);
	}
	

       }
      } // end loop over elements within workgroup

      free(stiff_mat);
      free(rhs_vect);

    
    } // end loop over workgroups (parallel!)
    
#endif    
    
#ifdef NO_COPYING
#pragma omp parallel for default(none) firstprivate(nr_work_groups, nr_elems_per_work_group, nr_parts_of_stiff_mat) firstprivate(num_shap, work_group_size, kernel_index, nreq, num_dofs, el_data_out, first_elem_kercall, last_elem_kercall, Max_dofs_int_ent, Problem_id, L_int_ent_type, L_int_ent_id,L_dof_elem_to_struct, L_dof_face_to_struct, L_dof_edge_to_struct, L_dof_vert_to_struct, Comp_type, field_id, pdeg, Level_id  )
    for(group_id=0;group_id<nr_work_groups;group_id++){
      int ielem_CPU = first_elem_kercall+nr_elems_per_work_group*group_id;
      double *stiff_mat;
      double *rhs_vect;
      int ielem;
      for(ielem=0;ielem<nr_elems_per_work_group;ielem++){
       if(ielem_CPU<=last_elem_kercall){

	int intent = ielem_CPU;
	int l_dof_ent_id[PDC_MAX_DOF_PER_INT], l_dof_ent_nrdof[PDC_MAX_DOF_PER_INT];
	int l_dof_ent_posglob[PDC_MAX_DOF_PER_INT];
	int l_dof_ent_type[PDC_MAX_DOF_PER_INT];
	int nr_dof_ent = PDC_MAX_DOF_PER_INT;
	int nrdofs_int_ent = Max_dofs_int_ent;
	char rewrite;
	// FOR REWRITING STIFF_MAT AND LOAD_VEC THEY ARE NOT COPUTED
	// only data necessary for asembling are obtained
	pdr_comp_stiff_mat(Problem_id, L_int_ent_type[intent],
			   L_int_ent_id[intent], PDC_NO_COMP, NULL,
			   &nr_dof_ent,l_dof_ent_type,l_dof_ent_id,l_dof_ent_nrdof,
			   &nrdofs_int_ent, NULL, NULL, &rewrite);
	


	
	int aux_offset = (group_id*nr_elems_per_work_group+ielem)
	                              *(num_shap*num_shap+num_shap)*nreq*nreq;
	int index_GPU = aux_offset;
	int index_GPU_v = aux_offset + num_dofs*num_dofs;

	stiff_mat = &el_data_out[index_GPU];
	rhs_vect = &el_data_out[index_GPU_v];

	/* change the option compute SM and RHSV to rewrite SM and RHSV */
	int Comp_sm = Comp_type;
	if(Comp_sm!=PDC_NO_COMP) Comp_sm += 3;
	apr_get_stiff_mat_data(field_id, L_int_ent_id[intent], Comp_sm, 'N',
			       pdeg, 0, &nr_dof_ent, l_dof_ent_type,
			       l_dof_ent_id, l_dof_ent_nrdof,
			       &nrdofs_int_ent, stiff_mat, rhs_vect);

#pragma omp critical(assembling)
    {
	
      pdr_assemble_local_stiff_mat(Problem_id, Level_id, Comp_type,
				   nr_dof_ent, l_dof_ent_type,
				   l_dof_ent_id,l_dof_ent_nrdof, 
				   stiff_mat, rhs_vect, &rewrite);
	
    }
	
	ielem_CPU++;

	
       }
      } // end loop over elements within workgroup
    } // end loop over workgroups
    
#endif    

#ifdef TIME_TEST
      t_end = time_clock();
      printf("EXECUTION TIME: copying output buffer to global stiffness matrix: %lf\n",
	   t_end-t_begin);

    total_time += t_end-t_begin; 
#endif
      
    } // the end of loop over kernel calls

    
    // free mapped memory regions
#ifdef OPENCL_HSA
    clSVMFree(context,el_data_in);
    clSVMFree(context,el_data_out);
#else
	#ifndef PINNED
		free(el_data_in);
		free(el_data_out);
	#else
		clEnqueueUnmapMemObject(command_queue, cmPinnedBufIn, el_data_in, 0, NULL, NULL);
		clReleaseMemObject(cmPinnedBufIn);
		clEnqueueUnmapMemObject(command_queue, cmPinnedBufOut, el_data_out, 0, NULL, NULL);
		clReleaseMemObject(cmPinnedBufOut);
	#endif
		clReleaseMemObject(cmDevBufIn);
		clReleaseMemObject(cmDevBufOut);
#endif
    //clReleaseMemObject(execution_parameters_dev);
    if(gauss_dat_dev_bytes > 0) clReleaseMemObject(gauss_dat_dev);
    if(shape_fun_dev_bytes > 0) clReleaseMemObject(shape_fun_dev);
    
  } // the end of loop over colors
  
#ifdef TIME_TEST
  clFinish(command_queue);
  printf("\nTOTAL EXECUTION TIME: %lf\n", total_time);
#endif

  
  return(1);
}


///////////////////////////////////////

