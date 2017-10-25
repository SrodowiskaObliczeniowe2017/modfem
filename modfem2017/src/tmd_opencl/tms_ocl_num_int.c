/***********************************************************************
File tms_ocl_num_int - OpenCl routines supporting streaming processing 
                       on accelerators (common to all problem modules)

Contains definitions of routines:

local:
  utr_create_assemble_stiff_mat_elem_accel - to create element stiffness matrices
                           and assemble them to the global SM using ACCELERATOR
------------------------------
History:
	08.2016 - Krzysztof Banas, pobanas@cyf-kr.edu.pl, initial version
        10.2016 - Jan Bielanski, Kazimierz Chlon - assembling into crs on gpu
*************************************************************************/

#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<assert.h>
#include<signal.h>
#include<limits.h>

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

/* interface for linear algebra packages */
#include "lin_alg_intf.h"

/* interface for control parameters and some general purpose functions */
/* from problem dependent module */
#include "pdh_control_intf.h"

#include "pdh_intf.h"

#include "tmh_intf.h"

// interface of opencl implementation
#include "../tmd_opencl/tmh_ocl.h"


#ifdef TUNING
    FILE *optf,*resf;
    FILE *headf;  //header file only for result titles
    //#define COUNT_OPER
#endif

#define TIME_TEST
#ifdef TIME_TEST
double t_begin;
double t_end;
double total_time;
#endif

// Master switch: GPU versus PHI versus CPU - controlled by compilation options !!!!!
// MUST BE COMPATIBLE WITH: tmd_ocl/tms_ocl_intf.c and tmd_ocl/tmh_ocl.h 
//#define OPENCL_CPU
#define OPENCL_GPU
//#define OPENCL_PHI

//Opencl_HSA is for Heterogenous System Architecture with Shared Virtual Memory eg. APU

#ifdef OPENCL_HSA
	#define OPENCL_GPU
#endif

// Several important switches for different variants of the algorithm:
// SWITCH 1: float versus double (MUST BE COMPATIBLE WITH KERNEL SWITCH!!!!)
// data type for integration
//#define SCALAR float
#define SCALAR double

// SWITCH 2: one_el_one_thread strategy versus one_el_one_workgroup strategy
#define ONE_EL_ONE_THREAD
//#define ONE_EL_ONE_WORKGROUP
//#define ONE_EL_TWO_THREADS

// Less important switches - hacks for specific versions of kernels
// SWITCH 3: generic conv-diff (with plenty of coeffcients) versus Laplace
//#define GENERIC_CONV_DIFF
//#define LAPLACE
//#define HEAT
// artificial example - coefficients constant for all integration points
//#define TEST_SCALAR

// SWITCH 4: size for  work-group
#ifdef OPENCL_CPU
#define WORK_GROUP_SIZE 8
#elif defined OPENCL_GPU
#define WORK_GROUP_SIZE 64
#elif defined OPENCL_PHI
#define WORK_GROUP_SIZE 16
#endif

// this should always be defined for GPUs
#ifdef OPENCL_GPU
  #define COAL_WRITE
#endif


///////// Declarations of internal procedures

/*---------------------------------------------------------
  tmr_ocl_create_num_int_kernel - for a selected platform, device and kernel index
---------------------------------------------------------*/
int tmr_ocl_create_num_int_kernel(
  int Platform_index,
  int Device_index,
  int Kernel_index,
  char* Kernel_name,
  const char* Kernel_file,
  int Monitor
    );


int tmr_prepare_create_assemble_accel(
  int Problem_id,
  int Monitor
				       )
{
  
  // find out the module name
  char module_name[100];
  pdr_module_introduce(module_name);
  
  // get problem name
  char problem_name[100];
  pdr_problem_name(Problem_id, problem_name);
  
  // choose the platform 
  int platform_index = tmr_ocl_get_current_platform_index();
  
  int device_tmc_type = tmr_ocl_get_current_device_type();
  tmv_ocl_struct.current_device_type=device_tmc_type;
  
  // choose device_index
  int device_index = tmr_ocl_get_current_device(); //tmr_ocl_select_device(platform_index, device_tmc_type);
  
  // OpenCL device characteristics stored in data structure
  tmt_ocl_device_struct device_struct = 
    tmv_ocl_struct.list_of_platforms[platform_index].list_of_devices[device_index];
  
  //cl_platform_id platform = tmv_ocl_struct.list_of_platforms[platform_index].id;
  tmt_ocl_platform_struct platform_struct = 
    tmv_ocl_struct.list_of_platforms[platform_index];
  
  /*-------------- KERNEL CREATION PHASE--------------------*/
  
  int kernel_index;
  
#ifdef TIME_TEST
  t_begin = omp_get_wtime();
  //t_begin = time_clock();
#endif
  
  
  if(Monitor>TMC_PRINT_INFO){
    if(device_tmc_type==TMC_OCL_DEVICE_CPU){
      printf( 
	      "\nSelecting OpenCL device type %d (CPU) for platform %d\n", 
	      device_tmc_type, platform_index);
    }
    else if(device_tmc_type==TMC_OCL_DEVICE_GPU){
      printf( 
	      "\nSelecting OpenCL device type %d (GPU) for platform %d\n", 
	      device_tmc_type, platform_index);
    }
    else if(device_tmc_type==TMC_OCL_DEVICE_ACCELERATOR){
      printf( 
	      "\nSelecting OpenCL device type %d (ACCELERATOR) for platform %d\n", 
	      device_tmc_type, platform_index);
    }
    else {
      printf( 
	      "\nUnknown OpenCL device type %d for platform %d. Exiting!\n", 
	      device_tmc_type, platform_index);
      exit(-1);
    }
  }
  
  cl_kernel kernel;
  // if required context exists
  if(device_index >= 0){
    
    // select existing or create new kernel for numerical integration
    kernel_index = TMC_OCL_KERNEL_NUM_INT_INDEX;
    
    kernel = tmr_ocl_select_kernel(platform_index, device_index, kernel_index); 
    
    if(kernel==NULL){
      
      
      if(Monitor>TMC_PRINT_INFO){
	printf( 
		"\nCreating OpenCL kernel %d for device type %d and device index %d\n", 
		kernel_index, device_tmc_type, device_index);
      }
      
      char arg[255]={0};
      sprintf(arg,"%s/tmr_ocl_num_int_el.cl", tmv_ocl_struct.kernel_source_directory);
      tmr_ocl_create_num_int_kernel(platform_index, device_index, kernel_index,
				    "tmr_ocl_num_int_el", arg, Monitor);
      
      /////////////////////////
      // THE INFRASTRUCTURE FOR DIFFERENT KERNEL VERSIONS NOT USED (YET???)
      /////////////////////////


/*       // !!!!!!!!!!!!!!!!!!!!!!????????????????!!!!!!!!!!!!!!!!!!!!!!! */
/*       // APART FROM DIFFERENT KERNELS FOR DIFFERENT HARDWARE, THERE MAY BE */
/*       // DIFFERENT KERNEL VERSIONS FOR DIFFERENT ALGORITHMS; */
/* // numerical integration kernel versions depending on hardware */
/* #define TMC_OCL_KERNEL_NUM_INT_GENERIC   0 */
/* #define TMC_OCL_KERNEL_NUM_INT_CPU_OPT   1 */
/* #define TMC_OCL_KERNEL_NUM_INT_GPU_OPT   2 */
/* #define TMC_OCL_KERNEL_NUM_INT_CELL_OPT  3 */
/* #define TMC_OCL_KERNEL_NUM_INT_PHI_OPT   4 */

/* // numerical integration kernel versions depending on algorithm */
/* #define TMC_OCL_KERNEL_NUM_INT_DEFAULT     0 // which one is default? */
/* #define TMC_OCL_KERNEL_NUM_INT_ONE_EL_ONE_THREAD    1 */
/* #define TMC_OCL_KERNEL_NUM_INT_ONE_EL_ONE_WORKGROUP 2 */
      
      
/*       // HARDCODED DEFAULT */
/*       platform_struct.list_of_devices[device_index].num_int_kernel_version = 0; */
/*       int kernel_version_hw = 0; */
/*       int kernel_version_alg = 0;  */
      
/*       // BELOW IT IS READ FROM INPUT: */
/*       // IT SHOULD BE BASED ON INPUT PARAMETERS FOR OPENCL !!!!!!!!!!!!!!! */
/*       /\* if(Interactive_input == stdin){ *\/ */
      
/*       /\*   printf("Select kernel version (kernel_version_hw + 10*kernel_version_alg): \n"); *\/ */
/*       /\*   printf("see tmd_ocl/tmh_ocl.h for details\n"); *\/ */
/*       /\*   // HARDCODED DEFAULT *\/ */
/*       /\*   int i = 0; *\/ */
/*       /\*   fscanf(Interactive_input, "%d", &i); *\/ */
/*       /\*   platform_struct.list_of_devices[device_index].num_int_kernel_version = i; *\/ */
/*       /\*   int kernel_version_hw = i%10; *\/ */
/*       /\*   int kernel_version_alg = i/10;  *\/ */
      
/*       /\* } *\/ */
      
/*   // to simplify stupid extra-long names from tmd_ocl/tmh_ocl.h */
/* #define ONE_EL_ONE_THREAD_KERNEL TMC_OCL_KERNEL_NUM_INT_ONE_EL_ONE_THREAD */
/* #define ONE_EL_TWO_THREADS_KERNEL TMC_OCL_KERNEL_NUM_INT_ONE_EL_TWO_THREADS */
/* #define ONE_EL_ONE_WORKGROUP_KERNEL TMC_OCL_KERNEL_NUM_INT_ONE_EL_ONE_WORKGROUP */
  
/* #ifdef ONE_EL_ONE_THREAD */
/*   kernel_version_alg = TMC_OCL_KERNEL_NUM_INT_ONE_EL_ONE_THREAD; */
/* #elif defined(ONE_EL_TWO_THREADS) */
/*   kernel_version_alg = TMC_OCL_KERNEL_NUM_INT_ONE_EL_TWO_THREADS; */
/* #elif defined ONE_EL_ONE_WORKGROUP */
/*   kernel_version_alg = TMC_OCL_KERNEL_NUM_INT_ONE_EL_ONE_WORKGROUP; */
/* #else */
/*   printf("wrong kernel version specified!!!). Exiting.\n"); */
/*   exit(-1); */
/* #endif */
  
/*       // BASED ON INPUT PARAMETERS SELECT OPENCL KERNEL FOR NUMERICAL INTEGRATION !!!!!!! */
/*       if(kernel_index == TMC_OCL_KERNEL_NUM_INT_GENERIC){ */
/* 	tmr_ocl_create_num_int_kernel(platform_index,device_index,kernel_index, */
/* 				      "tmr_ocl_num_int_el", arg, Monitor); */
/*       } */
      /* else if(kernel_index == TMC_OCL_KERNEL_NUM_INT_CPU_OPT){ */
      
      /*   tmr_ocl_create_kernel(Interactive_output,platform_index,device_index,kernel_index, */
      /* 	   // kernel name:         , file: */
      /*     "tmr_ocl_num_int_el", "tmr_ocl_num_int_el_cpu.cl", Monitor); */
      /* } */
      /* else if(kernel_index == TMC_OCL_KERNEL_NUM_INT_GPU_OPT){ */
      
      /*   tmr_ocl_create_kernel(Interactive_output,platform_index,device_index,kernel_index, */
      /* 	   // kernel name:         , file: */
      /*     "tmr_ocl_num_int_el", "tmr_ocl_num_int_el_gpu.cl", Monitor); */
      /* } */
      /* else if(kernel_index == TMC_OCL_KERNEL_NUM_INT_CELL_OPT){ */
      
      /*   tmr_ocl_create_kernel(Interactive_output,platform_index,device_index,kernel_index, */
      /* 	   // kernel name:         , file: */
      /*     "tmr_ocl_num_int_el", "tmr_ocl_num_int_el_cell.cl", Monitor); */
      /* } */
      
#ifdef TIME_TEST
      t_end = omp_get_wtime();
      //t_end = time_clock();
      printf("EXECUTION TIME: creating OpenCL kernel for integration: %.12lf\n", t_end-t_begin);
#endif
      
    } // end if new kernel created
    else{

      if(Monitor>TMC_PRINT_INFO){
	printf( 
		"\nFound existing OpenCL kernel %d for device type %d and device index %d\n", 
		kernel_index, device_tmc_type, device_index);
      }

    }
      
  }
  else{
    printf("Requested device not supported for this platform.\n");
    exit(-1);
  }
      
  return(1);
}
      
/*---------------------------------------------------------
----------------------------------------------------------*/

int tmr_perform_create_assemble_accel(
  int Problem_id, 
  int Level_id, 
  int Comp_type,         /* in: indicator for the scope of computations: */
  //extern const int PDC_NO_COMP  ; /* do not compute stiff matrix and rhs vector */
  //extern const int PDC_COMP_SM  ; /* compute entries to stiff matrix only */
  //extern const int PDC_COMP_RHS ; /* compute entries to rhs vector only */
  //extern const int PDC_COMP_BOTH; /* compute entries for sm and rhsv */
  int Nr_elems,             //nr elems for one color
  int First_elem_in_color_index, //index of the first element in color
  int* L_elem_id,           //list elements in one color
  int* Asse_pos_first_dof_int_ent, //position of first entry
  int* Assembly_table,      //all Assembly_table
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
)
{
  int i,j,k;

  int field_id = pdr_ctrl_i_params(Problem_id, 3);
  int mesh_id = apr_get_mesh_id(field_id);
  int nreq = apr_get_nreq(field_id);

  // calculate parameters for execution and allocate buffers in accelerator memory
  
  // choose the platform 
  int platform_index = tmr_ocl_get_current_platform_index();
  
  int device_tmc_type = tmr_ocl_get_current_device_type();
  
  // choose device_index
  int device_index = tmr_ocl_get_current_device(); //tmr_ocl_select_device(platform_index, device_tmc_type);
  
  // OpenCL device characteristics stored in data structure
  tmt_ocl_device_struct device_struct = 
    tmv_ocl_struct.list_of_platforms[platform_index].list_of_devices[device_index];
  double global_mem_bytes = device_struct.global_mem_bytes;
  double global_max_alloc = device_struct.global_max_alloc;
  double shared_mem_bytes = device_struct.shared_mem_bytes;
  double constant_mem_bytes = device_struct.constant_mem_bytes;
  int max_num_comp_units = device_struct.max_num_comp_units;
  int max_work_group_size = device_struct.max_work_group_size;  

#ifdef GPU_ASSEMBLING

  // OpenCL peration counter
  #ifdef COUNT_OPER
  SCALAR count_oper[3];
  cl_mem ocl_count_oper;
  #endif 

  
  /*kc
  printf("Nr_elems=%d\n Assembly_table:\n",Nr_elems);
  for(i=0;i<Nr_elems;i++){
    printf("%d ",Assembly_table[i]);
  }
    printf("\nNum_dofs=%d \n\nAsse_pos_first_dof_int_ent:\n",Num_dofs);
  for(i=0;i<Num_dofs;i++){
    printf("%d ",Asse_pos_first_dof_int_ent[i]);
  }
  /*kc*/
  
  
    /*kc
  printf("Nr_elems=%d\n (L_elem_id,Asse_pos_first_dof_int_ent[L_elem_id[i]]):\n",Nr_elems);
  for(i=0;i<Nr_elems;i++){
    printf("(%d,%d) ",L_elem_id[i],Asse_pos_first_dof_int_ent[L_elem_id[i]]);
  }
  getchar();getchar();
  printf("\n\nAssembly_table:\n");
  for(i=0;i<Nr_elems;i++){
    for(j=Asse_pos_first_dof_int_ent[L_elem_id[i]];j<Asse_pos_first_dof_int_ent[L_elem_id[i]+1];j++){
      printf("(%d,%d) ",Asse_pos_first_dof_int_ent[L_elem_id[i]],Assembly_table[Asse_pos_first_dof_int_ent[L_elem_id[i]]+j]);
      
    }getchar();
  }
  /*kc*/
#endif
  
/*kbw*/
  printf("\n\nStarting OpenCL GPU computations: platform %d, device %d, nr_elems %d\n", 
	 platform_index, device_index, Nr_elems);
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
  
  // !!!!! FOR THE TIME BEING WE ALWAYS READ AND COMPILE tmr_ocl_num_int_el.cl or tmr_ocl_num_int_el_and_assembling.cl
  // when assembling is enabled
  // FROM THE CURRENT DIRECTORY AND ASSIGN DEFAULT VERSION = 0 
  
  cl_kernel kernel = tmr_ocl_select_kernel(platform_index, device_index, kernel_index); 
  
  int kernel_version_alg = 0 ;  
  
#define ONE_EL_ONE_THREAD_KERNEL 0
#define ONE_EL_TWO_THREADS_KERNEL 1
#define ONE_EL_ONE_WORKGROUP_KERNEL 2
  
#ifdef ONE_EL_ONE_THREAD
  kernel_version_alg = ONE_EL_ONE_THREAD_KERNEL;
#elif defined(ONE_EL_TWO_THREADS)
  kernel_version_alg = ONE_EL_TWO_THREADS_KERNEL;
#elif defined ONE_EL_ONE_WORKGROUP
  kernel_version_alg = ONE_EL_ONE_WORKGROUP_KERNEL;
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
  
  //  - coordinates and weights - for the reference element
  global_memory_req_in += 4*Ngauss;
  
  // 2. GEO_DOFS - 3 coordinates for each vertex 
  //               (geometrically multilinear elements only!)
  global_memory_req_one_el_in += 3*Num_geo_dofs;
  
  // 3. SHAPE FUNCTIONS
  // for the reference element
  //  - all functions and derivatives at all integration points
  global_memory_req_in += 4*Ngauss*Num_shap;
  
  // 4. PDE COEFFICIENTS
  int all_el_pde_coeff_size = Size_global_pde_data; // size for all elements
  int one_el_pde_coeff_size = Size_per_element_pde_data; 
  int one_int_p_pde_coeff_size = Size_per_int_point_pde_data;
  
  
  // global coefficients - the same for all elements
  //global_memory_req_in += all_el_pde_coeff_size;
  // COEFFICIENTS DIFFERENT FOR EACH ELEMENT BUT THE SAME FOR ALL INTEGRATION POINTS
  global_memory_req_one_el_in += one_el_pde_coeff_size;
  // COEFFICIENTS DIFFERENT FOR EACH ELEMENT AND EACH INTEGRATION POINT ?????????????
  global_memory_req_one_el_in += Ngauss*one_int_p_pde_coeff_size;

  
  // 5. COMPUTED STIFFNESS MATRIX AND LOAD VECTOR
  int one_el_stiff_mat_size = Num_dofs*Num_dofs;
  int one_el_load_vec_size = Num_dofs;
  global_memory_req_one_el_out += one_el_stiff_mat_size + one_el_load_vec_size;
  
  int global_memory_req = global_memory_req_in + global_memory_req_out;
  int global_memory_req_one_el = global_memory_req_one_el_in + global_memory_req_one_el_out;
  
/*kbw*/
  printf("\nGlobal memory requirements (before computing number of elements considered):\n");
  printf("\treference element %d, for each element %d\n",
	 global_memory_req, global_memory_req_one_el);
/*kew*/
  

  
// GLOBAL GPU MEMORY CONSIDERATIONS FOR SETTING THE NUMBER OF ELEMENTS PER KERNEL
  
  // nr_elems_per_kernel is related to the size of input and output data for kernel
  // therefore:
  // 1. nr_elems_per_kernel is optimized to make the best use of global GPU memory
  // 2. nr_elems_per_kernel is optimized to speed up transfer times for input and output
  // 3. nr_elems_per_kernel may be optimized for overlapping communication with computations
  
  
  // in the current version the procedure is called for a single color
  int nr_elems_per_color = Nr_elems;
  
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
  
  // ASSUMED WORK_GROUP_SIZE
  int work_group_size = WORK_GROUP_SIZE;
  
  
  // how many workgroups will be enough to make GPU cores (comp_units) busy?
#ifdef OPENCL_CPU
#define NR_WORK_GROUPS_PER_COMP_UNIT  1 // ????
#elif defined OPENCL_GPU
#define NR_WORK_GROUPS_PER_COMP_UNIT  4 // ????
#elif defined OPENCL_PHI
#define NR_WORK_GROUPS_PER_COMP_UNIT  1 // ????
#endif
  
  // usually for GPUs it is good to maximize the number of work_groups and threads
  int nr_work_groups = NR_WORK_GROUPS_PER_COMP_UNIT * max_num_comp_units;
  int nr_elems_per_work_group;
  
  // CHECK 2 nr_elems_per_kernel can be larger than nr_elems_per_color but not too much
  // i.e. we can allocate more memory for kernel execution but not much more
  // than necessary

  double arbitrary_multiplication_constant = 1.2;  
  if( kernel_version_alg==ONE_EL_ONE_THREAD_KERNEL) {
    if(nr_elems_per_kernel > arbitrary_multiplication_constant * nr_elems_per_color) {
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
  

  // currently we assume that the number of elements per color is small enough
  // but later when global assembly is performed on GPU THIS SHOULD BE CHANGED
  if(nr_elems_per_kernel<nr_elems_per_color){
    
    printf("Data do not fit into accelerator memory! Change algorithm...\n");
    exit(-1);
    
  }

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
  
#ifdef TIME_TEST
  t_begin = omp_get_wtime();
  //t_begin = time_clock();
#endif
  
  // WE ALLOCATE OPENCL BUFFERS FOR nr_elems_per_kernel ELEMENTS
  // this should be good for all kernel invocations and all colors !!!
  // (actual kernel calculations will be for nr_elems_this_kercall)
  
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
  
  const size_t gauss_dat_dev_bytes = 4 * Ngauss * sizeof(SCALAR);
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
  
  
  const size_t shape_fun_dev_bytes  = 4 * Num_shap * Ngauss * sizeof(SCALAR);
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
    (all_el_pde_coeff_size + global_memory_req_one_el_in * nr_elems_per_kernel) * sizeof(SCALAR);
  
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
  

#ifdef GPU_ASSEMBLING
  
  /* -------------------------------------------------------------------- GPU ASSEMBLING ------------------------------------------------------------------------- */

  // Bind assembly tables
  retval |= clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &tmv_ocl_crs_struct.ocl_asse_pos_first_dof_int_ent);
  retval |= clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &tmv_ocl_crs_struct.ocl_assembly_table);
  
  /* ------------------- BIND CRS STRUCTURE WITH KERNEL -------------------- */

  // Bind RHS vector
  retval |= clSetKernelArg(kernel, 6, sizeof(cl_mem), (void *) &tmv_ocl_crs_struct.ocl_rhs_val);

  // Bind CRS structure
  retval |= clSetKernelArg(kernel, 7, sizeof(cl_mem), (void *) &tmv_ocl_crs_struct.ocl_crs_val);

  /* tmv_ocl_crs_struct.ocl_crs_col_ind and tmv_ocl_crs_struct.ocl_crs_row_ptr not used in tmr_ocl_num_int_el kernel */
  // retval |= clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &tmv_ocl_crs_struct.ocl_crs_col_ind);
  // retval |= clSetKernelArg(kernel, 9, sizeof(cl_mem), (void *) &tmv_ocl_crs_struct.ocl_crs_row_ptr);


  #ifdef COUNT_OPER
  ocl_count_oper = clCreateBuffer(context, CL_MEM_WRITE_ONLY,3*sizeof(SCALAR), NULL, NULL);
  retval |= clSetKernelArg(kernel, 8, sizeof(cl_mem), (void *) &ocl_count_oper);
  #endif

  /* -------------------------------------------------------------------- GPU ASSEMBLING ------------------------------------------------------------------------- */

#else

  /* !BEGIN! CREATE BUFER FOR THE MATRIX WITH RESULTS !BEGIN! */
  
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

  /* !END! CREATE BUFER FOR THE MATRIX WITH RESULTS !END! */

#endif
  
  if (retval != CL_SUCCESS) {
    printf("Failed to Set the kernel argument %d.\n", 4);
    exit(0);
  }

#ifdef TIME_TEST
  clFinish(command_queue);
  t_end = omp_get_wtime();
  //t_end = time_clock();
  printf("\nEXECUTION TIME: initial settings on CPU and creating buffers on GPU %.12lf\n",
	 t_end-t_begin);
  total_time += t_end-t_begin; 
#endif

  
/*kbw*/
#ifdef GPU_ASSEMBLING
  printf("\nAllocated buffer on device for %d elements data:\n", nr_elems_per_kernel );
  printf("\tel_data_in bytes %.3lf[MB]\n",
	 (double)el_data_in_bytes*1.0e-6);
#else
  printf("\nAllocated buffer on device for %d elements data:\n", nr_elems_per_kernel );
  printf("\tel_data_in bytes %.3lf[MB], el_data_out bytes %.3lf [MB]\n",
	 (double)el_data_in_bytes*1.0e-6, (double)el_data_out_bytes*1.0e-6);
#endif
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

  
  // for each color we set the ID for the first and the last element
  // in this version the procedure is called for a single color
  int first_elem_color = 0;
  int last_elem_color = Nr_elems-1; // IDs from 0 to Nr_elems-1
  int nr_elems_this_color = Nr_elems;
  
  // AFTER ALLOCATING DATA BUFFERS ON GPU CARD IN A LOOP OVER KERNEL CALLS
  // BUFFERS ARE FILLED WITH DATA, DATA ARE SEND TO GPU, CALCULATIONS PERFORMED
  // AND OUTPUT DATA TRANSFERRED BACK TO HOST COMPUTER MEMORY
  
  // the number of kernel calls per single color depends on the number of elements
  int nr_kernel_calls=1;
  if(nr_elems_per_kernel < nr_elems_this_color){
    
    nr_kernel_calls = ceil(((double) nr_elems_this_color) / nr_elems_per_kernel);
    
  }
  
  int nr_elems_per_kercall = nr_elems_per_kernel;
  
  // DOES NOT WORK FOR SMALL NUMBER OF ELEMENTS FOR COLOR 
  /* if( kernel_version_alg==ONE_EL_ONE_THREAD_KERNEL) { */
  /*   // if the last kernel call gets too low number of elements */
  /*   if(nr_elems_this_color - (nr_kernel_calls-1)*nr_elems_per_kernel <  */
  /*      nr_work_groups*work_group_size){ */
      
  /*     // we decrease the number of elements per kernel call */
  /*     nr_elems_per_kercall = nr_elems_per_kernel-nr_work_groups*work_group_size; */
      
  /*   } */
  /* } */
  
/*kbw
  printf("\nnr_colors %d, nr_elems_per_color %d, nr_elems_this_color %d\n",
  nr_colors, nr_elems_per_color, nr_elems_this_color);
//kew*/

  

/*kbw
  printf("\nfirst_elem_color %d, last_elem_color %d\n",
  first_elem_color, last_elem_color);
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
      
      
      // if nr_elems_per_kernel > nr_elems_per_color + nr_work_groups
      // nr_elems_per_color < nr_elems_per_work_group*nr_work_groups 
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
    t_begin = omp_get_wtime();
    //t_begin = time_clock();
#endif
    
    // FILL MEMORY OBJECTS WITH DATA AND SEND INPUT DATA TO GPU MEMORY
    
    
    // !!!!!!!!!!!!!!!***************!!!!!!!!!!!!!!!!!
    // first  fill buffers for assumed and computed execution parameters
    // there are NR_EXEC_PARAMS execution parameters (ALLWAYS CHECK!!!)
    int execution_parameters_host[NR_EXEC_PARAMS];
    
    execution_parameters_host[0] = nr_elems_per_thread;
    execution_parameters_host[1] = nr_elems_this_kercall;

#ifdef GPU_ASSEMBLING
    execution_parameters_host[2] = First_elem_in_color_index;
    
    /*jbw
    printf("\nCPU: Number of elements for color: %d ; ID of the first element in color: %d\n\n",nr_elems_this_kercall, First_elem_in_color_index);
    getchar();
    /*jbw*/
#endif
 
    //execution_parameters_host[1] = nr_parts_of_stiff_mat;
    //execution_parameters_host[2] = nr_blocks_per_thread;
    //execution_parameters_host[3] = nreq;
    //execution_parameters_host[4] = pdeg_single;
    //execution_parameters_host[5] = num_shap;
    //execution_parameters_host[6] = base;
    //execution_parameters_host[7] = Ngauss;
    //execution_parameters_host[8] = nr_work_groups;
    //execution_parameters_host[9] = nr_elems_per_kernel; 
    //execution_parameters_host[9] = work_group_size; - taken from OpenCL in kernels
    
    // send data to device - non-blocking call (CL_FALSE)
    clEnqueueWriteBuffer(command_queue, execution_parameters_dev, CL_TRUE, 0,
			 execution_parameters_dev_bytes, execution_parameters_host,
			 0, NULL, NULL);
    
	
/*kbw
    for(i=0;i<5;i++){
      
      printf("Gauss_dat_host[%d] = %f\n", i, Gauss_dat_host[i]);
      //printf("Gauss_dat_host[%d] = %lf\n", i, Gauss_dat_host[i]);
      
    }
    for(i=Ngauss-5;i<Ngauss;i++){
      
      printf("Gauss_dat_host[%d] = %f\n", i, Gauss_dat_host[i]);
      //printf("Gauss_dat_host[%d] = %lf\n", i, Gauss_dat_host[i]);
      
    }
//kew*/

    // send data to device - blocking version to ensure Gauss_dat_host is still valid
    clEnqueueWriteBuffer(command_queue, gauss_dat_dev, CL_TRUE, 0,
			 gauss_dat_dev_bytes, Gauss_dat_host, 0, NULL, NULL);
    
    
    
      
      
/*kbw
    for(i=0;i<5;i++){
      
      printf("Shape_fun_dat_host[%d] = %f\n", i, Shape_fun_dat_host[i]);
      //printf("Shape_fun_dat_host[%d] = %lf\n", i, Shape_fun_dat_host[i]);
      
    }
    for(i=Ngauss*Num_shap-5;i<Ngauss*Num_shap;i++){
      
      printf("Shape_fun_dat_host[%d] = %f\n", i, Shape_fun_dat_host[i]);
      //printf("Shape_fun_dat_host[%d] = %lf\n", i, Shape_fun_dat_host[i]);
      
    }
//kew*/

    // send data to device - CL_TRUE i.e. blocking to free Shape_fun_dat_host
    clEnqueueWriteBuffer(command_queue, shape_fun_dev, CL_TRUE, 0,
			 shape_fun_dev_bytes, Shape_fun_dat_host, 0, NULL, NULL);
    
    
    
    
    
    
    // in create_assemble we have to send to GPU coefficients for all elements
    // processed by kernel and for all integration points
    // WE REWRITE COEFFICIENT MATRICESS TO THE FORM THAT IS SUITABLE FOR THREADS
    // EXECUTING KERNELS AND ACCESSING GLOBAL AND LOCAL MEMORY
    int packed_bytes = 0;
#ifdef ONE_EL_ONE_WORKGROUP
    // for one_el_one_workgroup strategy we keep geo_dofs and coeffs for each element together
    int position_all = 0; // actual position in el_data_in and also the number of packed items
#endif
    
#ifdef ONE_EL_ONE_THREAD
    // for one_el_one_thread strategy we keep geo_dofs and coeffs for each element separately
    // i.e. geo_dofs for all elements and next coeffs for all elements
    int position_geo_dofs = 0;
    int position_coeff = 0;
    int offset_coeff = nr_elems_this_kercall*Num_geo_dofs*3;
#endif

    // here: write coeff for all elements
    for(i=0; i<all_el_pde_coeff_size; i++) {
      el_data_in[position_coeff+offset_coeff] = El_pde_dat_host[i];
      position_coeff++;
    }
    packed_bytes += Size_global_pde_data*sizeof(SCALAR);
    
    // here: simple rewriting from provided GEO_DAT and PDE_DAT matrices
    int ielem;
   
    for(ielem=first_elem_kercall; ielem<=last_elem_kercall; ielem++){
      
#ifdef ONE_EL_ONE_THREAD
      for(i=0;i<Num_geo_dofs;i++){
	el_data_in[position_geo_dofs] = El_geo_dat_host[position_geo_dofs];
	el_data_in[position_geo_dofs+1] = El_geo_dat_host[position_geo_dofs+1];
	el_data_in[position_geo_dofs+2] = El_geo_dat_host[position_geo_dofs+2];
	position_geo_dofs += 3;

	//printf("\t--> GEO DATA %d: %lf ; %lf ; %lf\n",i,el_data_in[position_geo_dofs-3],el_data_in[position_geo_dofs-2],el_data_in[position_geo_dofs-1]);
      }
      //getchar(); getchar();
#endif
      
#ifdef ONE_EL_ONE_WORKGROUP
      for(i=0;i<Num_geo_dofs;i++){
	el_data_in[position_all] = El_geo_dat_host[position_geo_dofs];
	el_data_in[position_all+1] = El_geo_dat_host[position_geo_dofs+1];
	el_data_in[position_all+2] = El_geo_dat_host[position_geo_dofs+2];
	position_geo_dofs += 3;
	position_all += 3;
      }
#endif
      
      packed_bytes += Num_geo_dofs*3*sizeof(SCALAR);
      
      for(i=0; i< Size_per_element_pde_data; i++){
	
	el_data_in[position_coeff+offset_coeff] = El_pde_dat_host[position_coeff];
	position_coeff++;
      }

      packed_bytes += Size_per_element_pde_data*sizeof(SCALAR);
      
      int ki;
      for(ki=0;ki<Ngauss;ki++){
	
	for(i=0; i< Size_per_int_point_pde_data; i++){
	  
	  el_data_in[position_coeff+offset_coeff] = El_pde_dat_host[position_coeff];
	  position_coeff++;

	}

	packed_bytes += Size_per_int_point_pde_data*sizeof(SCALAR);
	
      }
      
    }
   
    // !!!OPT_PDT!!!
    
    // for many problems the coefficients can be sent not as full matrices but
    // in the form of problem specific parameters for which matrices are calculated
    // by kernels
    
    
/*kbw
      for(i=0;i<50;i++){
	
	printf("el_data_in[%d] = %f\n", i, el_data_in[i]);
	//printf("el_data_in[%d] = %lf\n", i, el_data_in[i]);
	
      }
      for(i=offset_coeff+position_coeff-5;i<offset_coeff+position_coeff;i++){
	printf("el_data_in[%d] = %f\n", i, el_data_in[i]);
	//printf("el_data_in[%d] = %lf\n", i, el_data_in[i]);
	
      }
//kew*/

#ifdef TIME_TEST
    clFinish(command_queue);
    t_end = omp_get_wtime();
    //t_end = time_clock();
    printf("\nEXECUTION TIME: filling buffers and sending to GPU memory %.12lf\n", 
	   t_end-t_begin);
    total_time += t_end-t_begin;
    t_begin = omp_get_wtime();
    //t_begin = time_clock();
#endif

      
    // writing input data to global GPU memory

#ifndef OPENCL_HSA
    clEnqueueWriteBuffer(command_queue, cmDevBufIn, CL_FALSE, 0,
			 el_data_in_bytes, el_data_in, 0, NULL, NULL);
#endif
    
#ifdef TIME_TEST
    clFinish(command_queue);
    t_end = omp_get_wtime();
    //t_end = time_clock();
    printf("\nEXECUTION TIME: sending el_data_in to GPU memory %.12lf (speed %lf [GB/s])\n", 
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
    t_begin = omp_get_wtime();
    //t_begin = time_clock();
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

    t_end = omp_get_wtime();
    //t_end = time_clock();
    
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
    printf("EXECUTION TIME: executing kernel: %.12lf (internal %lf)\n",
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
    t_begin = omp_get_wtime();
    //t_begin = time_clock();
#endif
    
    

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
      for(i=Ngauss-5;i<Ngauss;i++){
	
	//printf("gauss_dat_host[%d] = %f\n", i, gauss_dat_host[i]);
	printf("gauss_dat_host[%d] = %lf\n", i, gauss_dat_host[i]);
	
      }
    }     
    
    if(shape_fun_dev_bytes!=0){
      
      int shape_fun_host_bytes = shape_fun_dev_bytes;
      SCALAR *shape_fun_host = (SCALAR*)malloc(shape_fun_host_bytes);
      memset(shape_fun_host, 0, shape_fun_host_bytes);
      clEnqueueReadBuffer(command_queue, shape_fun_dev, CL_TRUE, 0,
			  shape_fun_dev_bytes, shape_fun_host, 0, NULL, NULL);
      
      for(i=0;i<5;i++){
	
	//printf("shape_fun_host[%d] = %f\n", i, shape_fun_host[i]);
	printf("shape_fun_host[%d] = %lf\n", i, shape_fun_host[i]);
	
      }
      for(i=Num_shap-5;i<Num_shap;i++){
	//printf("shape_fun_host[%d] = %f\n", i, shape_fun_host[i]);
	printf("shape_fun_host[%d] = %lf\n", i, shape_fun_host[i]);
	
      }
      free(  shape_fun_host);
      
    }
    
    memset(el_data_in, 0, el_data_in_bytes);
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
    
#ifndef GPU_ASSEMBLING
    
#ifdef TIME_TEST
    t_begin = omp_get_wtime();
    //t_begin = time_clock();
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
    /*kew*/
    
#ifdef TIME_TEST
    clFinish(command_queue);
    t_end = omp_get_wtime();
    //t_end = time_clock();
    printf("EXECUTION TIME: copying output buffer: %.12lf (total with kernel %.12lf)\n",
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
    t_begin = omp_get_wtime();
    //t_begin = time_clock();
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

#endif
    
//////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
// ALL THE CODE BELOW SHOULD BE MOVED TO ACCELERATOR
// GPU SHOULD ASSEMBLE ITS PART OF STIFFNESS MATRIX AND LOAD VECTOR
// AND RETURN THEM TO utr_create_assemble_stiff_mat_elem_accel
//
// THEN ASSEMBLED PARTS OF GLOBAL STIFFNESS MATRIX AND LOAD VECTOR
// OBTAINED ON CPU AND GPU SHOULD BE COMBINED INTO THE FINAL FORM
//////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


    // // perform streaming assembly (currently on CPU)
    // int int_ent;
    // for(int_ent = 0; int_ent<0; int_ent++){
      
      // double *stiff_mat = NULL;
      // double *rhs_vect = NULL;
      // char rewrite;
      // int l_dof_ent_id[PDC_MAX_DOF_PER_INT], l_dof_ent_nrdof[PDC_MAX_DOF_PER_INT];
      // int l_dof_ent_posglob[PDC_MAX_DOF_PER_INT];
      // int l_dof_ent_type[PDC_MAX_DOF_PER_INT];
      // int nr_dof_ent = PDC_MAX_DOF_PER_INT;
      
      // if(Assembly_table != NULL){     
	
	// int position = Asse_pos_first_dof_int_ent[int_ent];
	
// /*kbw
	    // printf("Assembling int_ent %d, using assembly_table at position %d\n",
		   // int_ent, position);
// /*kew*/

	// pdr_assemble_local_stiff_mat_with_table(Problem_id, Level_id, Comp_type,
						// nr_dof_ent,
						// &Assembly_table[position],
						// stiff_mat, rhs_vect, &rewrite);
	
      // }
      // else{
	
	// pdr_assemble_local_stiff_mat(Problem_id, Level_id, Comp_type,
				     // nr_dof_ent, l_dof_ent_type,
				     // l_dof_ent_id,l_dof_ent_nrdof, 
				     // stiff_mat, rhs_vect, &rewrite);
	
      // }
    
    // }


#ifndef GPU_ASSEMBLING

    #define NO_TESTING_CORRECTNESS
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
	int wg,thr,elem; int num_shap = Num_shap;
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
    int wg,thr,elem;int num_shap = Num_shap;
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
    
    
#pragma omp parallel for default(none) firstprivate(nr_work_groups, nr_elems_per_work_group, nr_parts_of_stiff_mat) firstprivate(work_group_size, kernel_index, el_data_out, first_elem_kercall, last_elem_kercall, Max_dofs_int_ent, Problem_id, L_elem_id, Comp_type, field_id, nreq, Level_id, nr_blocks_per_thread, kernel_version_alg)
    /* #pragma omp parallel for default(none) firstprivate(nr_work_groups, nr_elems_per_work_group, nr_parts_of_stiff_mat, nr_iter_within_part, num_shap, work_group_size, kernel_index, nreq, num_dofs, stiff_mat, el_data_out) shared(first_elem_kercall) */
    for(group_id=0;group_id<nr_work_groups;group_id++){
      int ielem_CPU = first_elem_kercall+nr_elems_per_work_group*group_id;
      int el_id = L_elem_id[0];
      int pdeg = apr_get_el_pdeg(field_id, el_id, NULL);  
      int num_shap = apr_get_el_pdeg_numshap(field_id, el_id, &pdeg); 
      int num_dofs = num_shap * nreq;
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
	  
	  
	  pdr_comp_stiff_mat(Problem_id, PDC_ELEMENT,
			     L_elem_id[intent], PDC_COMP_BOTH, NULL,
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
			   aux_offset, group_id, ielem, L_elem_id[intent], ipart, iiter, idofs, ieq2);
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
    if(L_elem_id[intent]!=-1) {
      printf("pdr_create_assemble: before assemble: Solver_id %d, level_id %d, sol_typ %d\n", 
	     Problem_id, level_id, Comp_type);
      int ibl,jbl,pli,plj,nri,nrj,nrdof,jaux;
      printf("ient %d, int_ent_id %d, nr_dof_ent %d\n", 
	     intent, L_elem_id[intent], nrdofbl);
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
	  apr_get_stiff_mat_data(field_id, L_elem_id[intent], Comp_sm, 'N',
				 pdeg, 0, &nr_dof_ent, l_dof_ent_type,
				 l_dof_ent_id, l_dof_ent_nrdof,
				 &nrdofs_int_ent, stiff_mat, rhs_vect);
	  
	  
	  //#pragma omp critical(assembling)
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
#pragma omp parallel for default(none) firstprivate(nr_work_groups, nr_elems_per_work_group, nr_parts_of_stiff_mat) firstprivate(num_shap, work_group_size, kernel_index, nreq, num_dofs, el_data_out, first_elem_kercall, last_elem_kercall, Max_dofs_int_ent, Problem_id, L_elem_id,L_dof_elem_to_struct, L_dof_face_to_struct, L_dof_edge_to_struct, L_dof_vert_to_struct, Comp_type, field_id, pdeg, Level_id  )
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
	  pdr_comp_stiff_mat(Problem_id, PDC_ELEMENT,
			     L_elem_id[intent], PDC_NO_COMP, NULL,
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
	  apr_get_stiff_mat_data(field_id, L_elem_id[intent], Comp_sm, 'N',
				 pdeg, 0, &nr_dof_ent, l_dof_ent_type,
				 l_dof_ent_id, l_dof_ent_nrdof,
				 &nrdofs_int_ent, stiff_mat, rhs_vect);
	  
	  //#pragma omp critical(assembling)
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
    t_end = omp_get_wtime();
    //t_end = time_clock();
    printf("EXECUTION TIME: copying output buffer to global stiffness matrix: %.12lf\n",
	   t_end-t_begin);
    
    total_time += t_end-t_begin; 
#endif
  #endif  //NOT DEFINED GPU_ASSEMBLING


    #ifdef GPU_ASSEMBLING
      #ifdef COUNT_OPER
        clEnqueueReadBuffer(command_queue, ocl_count_oper, CL_TRUE, 0, 3*sizeof(SCALAR), count_oper, 0, NULL, NULL);

	printf("CRS --> Arthmetic operations: %.0lf, Local mem access: %.0lf, Global mem access: %.0lf\n",count_oper[0],count_oper[1],count_oper[2]);
        clReleaseMemObject(ocl_count_oper);
      #endif
    #endif


    
  } // the end of loop over kernel calls
  
  
    // free mapped memory regions
#ifdef OPENCL_HSA
  clSVMFree(context,el_data_in);
#ifndef GPU_ASSEMBLING
  clSVMFree(context,el_data_out);
#endif
#else
#ifndef PINNED
  free(el_data_in);
#ifndef GPU_ASSEMBLING
  free(el_data_out);
#endif
#else
  clEnqueueUnmapMemObject(command_queue, cmPinnedBufIn, el_data_in, 0, NULL, NULL);
  clReleaseMemObject(cmPinnedBufIn);
#ifndef GPU_ASSEMBLING
  clEnqueueUnmapMemObject(command_queue, cmPinnedBufOut, el_data_out, 0, NULL, NULL);
  clReleaseMemObject(cmPinnedBufOut);
#endif
#endif
  clReleaseMemObject(cmDevBufIn);
#ifndef GPU_ASSEMBLING
  clReleaseMemObject(cmDevBufOut);
#endif
#endif
  //clReleaseMemObject(execution_parameters_dev);
  if(gauss_dat_dev_bytes > 0) clReleaseMemObject(gauss_dat_dev);
  if(shape_fun_dev_bytes > 0) clReleaseMemObject(shape_fun_dev);
  
  
#ifdef TIME_TEST
  clFinish(command_queue);
  printf("\nTOTAL EXECUTION TIME: %lf\n", total_time);
#endif
  

 
  
  return(0);
} 


/*---------------------------------------------------------
  tmr_ocl_create_num_int_kernel - for a selected platform, device and kernel index
---------------------------------------------------------*/
int tmr_ocl_create_num_int_kernel(
  int Platform_index,
  int Device_index,
  int Kernel_index,
  char* Kernel_name,
  const char* Kernel_file,
  int Monitor
)
{

  FILE* Interactive_output = stdout;
  cl_int retval;

  // choose the platform
  tmt_ocl_platform_struct platform_struct = tmv_ocl_struct.list_of_platforms[Platform_index];

  // check the device !!!!!!!!!!!!!!!!! (or at least its index)
  if(Device_index < 0){
    printf("Wrong device_index %d passed to tmr_ocl_create_kernel! Exiting.\n",
	   Device_index);
    exit(-1);
  } 

  cl_device_id device = platform_struct.list_of_devices[Device_index].id;
  cl_context context = platform_struct.list_of_contexts[ 
			 platform_struct.list_of_devices[Device_index].context_index
							 ];

  if(Monitor>TMC_PRINT_INFO){
    printf("Program file is: %s\n", Kernel_file);
  }

#ifdef TUNING
  //for assembler
  system("rm -rf ~/.nv");
#endif

  // read source file into data structure
  const char* source = tmr_ocl_readSource(Interactive_output,Kernel_file);

  cl_program program = clCreateProgramWithSource(context, 1,
				      &source,
				      NULL, NULL);
  if (program == NULL)
    {
      printf("Failed to create CL program from source.\n");
      exit(-1);
    }

#ifdef TUNING

	  optf = fopen("options.txt", "r");
	  if(!optf) {
		 printf("Could not open options file!\n");
		 exit(-1);
	  }
	  resf = fopen("result.csv", "a+");
	  if(!resf) {
		 printf("Could not open results file!\n");
		 exit(-1);
	  }

	  #define NOPT 9

	  char opt[NOPT][40];
	  char buffer[40*(NOPT+10)];
	  int options[NOPT];
	  int tmp;
	  int i;
	  for(i=0;i<40*(NOPT+10);i++)
		  buffer[i]=0;



	  sprintf(&opt[0]," -D COAL_READ");
	  sprintf(&opt[1]," -D COAL_WRITE");
	  sprintf(&opt[2]," -D COMPUTE_ALL_SHAPE_FUN_DER");
	  sprintf(&opt[3]," -D USE_WORKSPACE_FOR_PDE_COEFF");
	  sprintf(&opt[4]," -D USE_WORKSPACE_FOR_GEO_DATA");
	  sprintf(&opt[5]," -D USE_WORKSPACE_FOR_SHAPE_FUN");
	  sprintf(&opt[6]," -D USE_WORKSPACE_FOR_STIFF_MAT");
	  sprintf(&opt[7]," -D WORKSPACE_PADDING=0");
	  sprintf(&opt[8]," -D WORKSPACE_PADDING=1");

//	  sprintf(&opt[0]," -D USE_WORKSPACE");
//	  sprintf(&opt[1]," -D USE_REGISTERS_FOR_COEFF");
//	  sprintf(&opt[2]," -D USE_SHAPE_FUN_WORKSPACE");
//	  sprintf(&opt[3]," -D USE_REGISTERS_FOR_SHAPE_FUN");
//	  sprintf(&opt[4]," -D USE_SHAPE_FUN_REF_DIRECTLY");
//	  sprintf(&opt[5]," -D STIFF_MAT_IN_SHARED");
//	  sprintf(&opt[6]," -D PADDING=0");
//	  sprintf(&opt[7]," -D PADDING=1");

	  //sprintf(&opt[7]," -D COAL_READ");
	  //sprintf(&opt[8]," -D COAL_WRITE");
	  //sprintf(&opt[9]," -D CONSTANT_COEFF");
	  //sprintf(&opt[9]," -D COUNT_OPER");

	  unsigned long line_count = 0;
	  int c;

	  while ( (c=fgetc(resf)) != EOF ) {
		 if ( c == '\n' )
				line_count++;
	  }

	  printf("Result file has %u lines\n", line_count);
	  if(line_count==0)
	  {
		  headf = fopen("header.csv", "w");
			  if(!headf) {
				 printf("Could not open header file!\n");
				 exit(-1);
		  }
		  system("touch num");
	  }

	  for(i=0;i<NOPT*line_count;i++)
		  fscanf(optf,"%d",&tmp);

	  for(i=0;i<NOPT;i++)
		  fscanf(optf,"%d",&options[i]);

	  printf("Options indicators:\n");
	  for(i=0;i<NOPT;i++)
	  {
		  if(line_count==0)
		  {
			  fprintf(headf, "%s,",opt[i]);
		  }
		  printf("%d ,",options[i]);
		  fprintf(resf,"%d ,",options[i]);
	  }
	  //fprintf(resf,";");
	  printf("Options:\n");
	  strcpy(buffer," ");
	  for(i=0;i<NOPT;i++)
	  {
		  if(options[i]==1)
		  {
			  strcat(buffer,opt[i]);
			  if(i==1)
			  {
				  system("touch coal_write");
			  }
		  }
		  else
			  strcat(buffer," ");
	  }

	#ifdef OPENCL_CPU
	 strcat(buffer," -D WORK_GROUP_SIZE=8");
	#elif defined OPENCL_GPU
	 strcat(buffer," -D WORK_GROUP_SIZE=64");
	#elif defined OPENCL_PHI
	 strcat(buffer," -D WORK_GROUP_SIZE=16");
	#endif

	#ifdef LAPLACE
	 strcat(buffer," -D LAPLACE");
	#elif TEST_SCALAR
	 strcat(buffer," -D TEST_NUMINT");
	#elif HEAT
	 strcat(buffer," -D HEAT");
	#endif

	#ifdef COUNT_OPER
	   strcat(buffer," -D COUNT_OPER");
	#endif

#ifdef OPENCL_HSA
	strcat(buffer," -cl-std=CL2.0");
#endif

	if(tmr_ocl_check_vendor()==1)
		strcat(buffer," -cl-nv-verbose");

	if(tmr_ocl_check_vendor()==3)
	{
		//strcat(buffer," -fno-bin-source -fno-bin-llvmir -fno-bin-amdil -save-temps=");
		strcat(buffer," -save-temps=");
		system("mkdir kernele");
		char name[50];
		int word;
		FILE* fp;
		for(i=0;i<50;i++)
		  name[i]=0;

		strcat(name,"kernele/");
		for(i=0;i<NOPT;i++)
		{
		  if(options[i]==1)
			  strcat(name,"1");
		  else
			  strcat(name,"0");
		}


#ifdef OPENCL_GPU
	#ifdef OPENCL_HSA
		strcat(name,"_HSA");
	#else
		strcat(name,"_GPU");
	#endif
#endif
#ifdef OPENCL_CPU
	  strcat(name,"_CPU");
#endif

		word=wc("tmr_ocl_num_int_el.cl","//QSS");
		if(word==1)
		  strcat(name,"_QSS");
		word=wc("tmr_ocl_num_int_el.cl","//SQS");
		if(word==1)
		  strcat(name,"_SQS");
		word=wc("tmr_ocl_num_int_el.cl","//SSQ");
		if(word==1)
		  strcat(name,"_SSQ");

		strcat(buffer,name);

	}

	printf("%s\n",buffer);

	retval = clBuildProgram(program, 0, NULL, buffer, NULL, NULL);

#else

#ifdef OPENCL_HSA
  retval = clBuildProgram(program, 0, NULL, "-cl-std=CL2.0", NULL, NULL);
#else

  if(tmr_ocl_check_vendor()==1)
  {
	  // TO GET INFO FROM NVIDIA COMPILER
	  retval = clBuildProgram(program, 0, NULL, "-cl-nv-verbose", NULL, NULL);
	  // TO FORCE NVIDIA COMPILER TO LIMIT THE NUMBER OF USED REGISTERS
	  // retval = clBuildProgram(program, 0, NULL, "-cl-nv-maxrregcount=32", NULL, NULL);
  }
  else
  {
	  // generic call - no compiler options
	  retval = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
	  //retval = clBuildProgram(program, 0, NULL, "-cl-nv-verbose", NULL, NULL);  
  }

#endif

#endif //TUNING


  char* buildLog; size_t size_of_buildLog; 
  clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 
			0, NULL, &size_of_buildLog); 
  buildLog = (char *)malloc(size_of_buildLog+1); 
  clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG, 
			size_of_buildLog, buildLog, NULL); 
  buildLog[size_of_buildLog]= '\0'; 
  if(Monitor>TMC_PRINT_INFO){
    printf("Kernel buildLog: %s\n", buildLog);
  }
  if (retval != CL_SUCCESS)
	{
	  printf("Error building program in tmr_ocl_create_kernel\n");
	  #ifdef TUNING
		  //fprintf(resf,"Kernel buildLog: %s", buildLog);
	  	  fprintf(resf,"Kernel error - too much shared data\n");
		  fclose(resf);
	  #endif
	  exit(-1);
	}


#ifdef TUNING

  if(tmr_ocl_check_vendor()==1)
  {
	  char *s;
	  int t=0;
	  int stack;
	  int stores;
	  int loads;
	  int registers;
	  int smem;
	  int cmem;
	  int cmem2;
	  int nr;

	  s = strstr(buildLog,"stack");
	  if (s != NULL)
	  {
		   t = (int)(s - buildLog);
		   //printf("Found string at index = %d\n", t );
		   sscanf(&buildLog[t-10],"%d",&stack);
		   //printf("stack=%d\n",stack);
	  }
	  else
	  {
		   //printf("String not found\n");
		  stack = 0;
	  }
	  if(line_count==0)
	  {
		  fprintf(headf, "bytes stack frame,");
	  }
	  fprintf(resf,"%d,",stack);
	  s = strstr(buildLog,"bytes spill stores");
	  if (s != NULL)
	  {
		   t = (int)(s - buildLog);
		   //printf("Found string at index = %d\n", t );
		   sscanf(&buildLog[t-10],"%*s %d",&stores);
		   //printf("stores=%d\n",stores);
	  }
	  else
	  {
		   //printf("String not found\n");
		  stores = 0;
	  }
	  if(line_count==0)
	  {
		  fprintf(headf, "bytes spill stores,");
	  }
	  fprintf(resf,"%d,",stores);
	  s = strstr(buildLog,"bytes spill loads");
	  if (s != NULL)
	  {
		  t = (int)(s - buildLog);
		  //printf("Found string at index = %d\n", t );
		  sscanf(&buildLog[t-10],"%*s %d",&loads);
		  //printf("loads=%d\n",loads);
	  }
	  else
	  {
		  //printf("String not found\n");
		  loads = 0;
	  }
	  if(line_count==0)
	  {
		  fprintf(headf, "bytes spill loads,");
	  }
	  fprintf(resf,"%d,",loads);
	  s = strstr(buildLog,"registers");
	  if (s != NULL)
	  {
		  t = (int)(s - buildLog);
		  //printf("Found string at index = %d\n", t );
		  sscanf(&buildLog[t-7],"%*s %d",&registers);
		  //printf("registers=%d\n",registers);
	  }
	  else
	  {
		  //printf("String not found\n");
		  registers = 0;
	  }
	  if(line_count==0)
	  {
		  fprintf(headf, "registers,");
	  }
	  fprintf(resf,"%d,",registers);
	  s = strstr(buildLog,"smem");
	  if (s != NULL)
	  {
		  t = (int)(s - buildLog);
		  //printf("Found string at index = %d\n", t );
		  sscanf(&buildLog[t-20],"%*s %d",&smem);
		  //printf("smem=%d\n",smem);
	  }
	  else
	  {
		  //printf("String not found\n");
		  smem = 0;
	  }
	  if(line_count==0)
	  {
		  fprintf(headf, "bytes smem,");
	  }
	  fprintf(resf,"%d,",smem);
	  s = strstr(buildLog,"cmem[0]");
	  if (s != NULL)
	  {
		  t = (int)(s - buildLog);
		  //printf("Found string at index = %d\n", t );
		  sscanf(&buildLog[t-10],"%d",&cmem);
		  //printf("cmem[0]=%d\n",cmem);
	  }
	  else
	  {
		  //printf("String not found\n");
		  cmem = 0;
	  }
	  if(line_count==0)
	  {
		  fprintf(headf, "bytes cmem[0],");
	  }
	  fprintf(resf,"%d,",cmem);
	  s = strstr(&buildLog[t+9],"cmem[");
	  if (s != NULL)
	  {
		  t = (int)(s - buildLog);
		  //printf("Found string at index = %d\n", t );
		  sscanf(&buildLog[t-10],"%d",&cmem2);
		  sscanf(&buildLog[t+5],"%d",&nr);
		  //printf("cmem[%d]=%d\n",nr,cmem2);
	  }
	  else
	  {
		  //printf("String not found\n");
		  nr=-1;
		  cmem2 = 0;
	  }
	  if(line_count==0)
	  {
		  fprintf(headf, "bytes cmem[%d],",nr);
	  }
	  fprintf(resf,"%d,",cmem2);

  }//end nvidia

#endif

  free(buildLog);

  // Create OpenCL kernel
  cl_kernel kernel = clCreateKernel(program, Kernel_name, NULL);
  if (kernel == NULL)
    {
      printf("Failed to create kernel.\n");
      exit(-1);
    }
  
  if(Monitor>TMC_PRINT_INFO){
    printf("Created kernel for platform %d, device %d, kernel index %d\n",
	   Platform_index, Device_index, Kernel_index);
  }
  
#ifdef TUNING

  if(tmr_ocl_check_vendor()==1)
  {
	  system("mkdir kernele");
	  char name[50];
	  FILE* fp;
	  size_t bin_sz;
	  int word;
	  for(i=0;i<50;i++)
		  name[i]=0;

	  strcat(name,"kernele/");
	  for(i=0;i<NOPT;i++)
	  {
		  if(options[i]==1)
			  strcat(name,"1");
		  else
			  strcat(name,"0");
	  }

	  word=wc("tmr_ocl_num_int_el.cl","//QSS");
	  if(word==1)
		  strcat(name,"_QSS");
	  word=wc("tmr_ocl_num_int_el.cl","//SQS");
	  if(word==1)
		  strcat(name,"_SQS");
	  word=wc("tmr_ocl_num_int_el.cl","//SSQ");
	  if(word==1)
		  strcat(name,"_SSQ");

	  strcat(name,"_kernel.ptx");

	  printf("Kernel: %s\n",name);

	  clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES, sizeof(size_t), &bin_sz, NULL);

	  // Read binary (PTX file) to memory buffer
	  unsigned char *bin = (unsigned char *)malloc(bin_sz);
	  clGetProgramInfo(program, CL_PROGRAM_BINARIES, sizeof(unsigned char *), &bin, NULL);

	  printf("TU!!!\n");


	  // Save PTX to add_vectors_ocl.ptx
	  fp = fopen(name, "wb");
	  fwrite(bin, sizeof(char), bin_sz, fp);
	  fclose(fp);
	  free(bin);


	  word=wc(name,"add.f64");
	  if(line_count==0)
	  {
		  fprintf(headf, "add.f64,");
	  }
	  fprintf(resf,"%d,",word);
	  //printf("add.f64=%d\n",word);
	  word=wc(name,"sub.f64");
	  if(line_count==0)
	  {
		  fprintf(headf, "sub.f64,");
	  }
	  fprintf(resf,"%d,",word);
	  word=wc(name,"mul.f64");
	  if(line_count==0)
	  {
		  fprintf(headf, "mul.f64,");
	  }
	  fprintf(resf,"%d,",word);
	  word=wc(name,"fma.rn.f64");
	  if(line_count==0)
	  {
		  fprintf(headf, "fma.rn.f64,");
	  }
	  fprintf(resf,"%d,",word);
	  word=wc(name,"neg.f64");
	  if(line_count==0)
	  {
		  fprintf(headf, "neg.f64,");
	  }
	  fprintf(resf,"%d,",word);
	  word=wc(name,"rcp.rn.f64");
	  if(line_count==0)
	  {
		  fprintf(headf, "rcp.rn.f64,");
	  }
	  fprintf(resf,"%d,",word);

		word=wc(name,"global.f64");
		if(line_count==0)
		{
			fprintf(headf, "global.f64,");
		}
		fprintf(resf,"%d,",word);
		word=wc(name,"shared.f64");
		if(line_count==0)
		{
			fprintf(headf, "shared.f64,");
		}
		fprintf(resf,"%d,",word);
		word=wc(name,"const.f64");
		if(line_count==0)
		{
			fprintf(headf, "const.f64,");
		}
		fprintf(resf,"%d,",word);
		word=wc(name,"local.f64");
		if(line_count==0)
		{
		  fprintf(headf, "local.f64,");
		}
		fprintf(resf,"%d,",word);


  }
  else if(tmr_ocl_check_vendor()==2)
  {
	  system("mkdir kernele");
	  char name[50];
	  int word;
	  for(i=0;i<50;i++)
		  name[i]=0;

	  strcat(name,"kernele/");
	  for(i=0;i<NOPT;i++)
	  {
		  if(options[i]==1)
			  strcat(name,"1");
		  else
			  strcat(name,"0");
	  }

	  word=wc("tmr_ocl_num_int_el.cl","//QSS");
	  if(word==1)
		  strcat(name,"_QSS");
	  word=wc("tmr_ocl_num_int_el.cl","//SQS");
	  if(word==1)
		  strcat(name,"_SQS");
	  word=wc("tmr_ocl_num_int_el.cl","//SSQ");
	  if(word==1)
		  strcat(name,"_SSQ");

#ifdef OPENCL_CPU
	  strcat(name,"_kernel_CPU.ptx");
#endif
#ifdef OPENCL_PHI
	  strcat(name,"_kernel_PHI.ptx");
#endif


	  printf("Kernel: %s\n",name);

	  char compile[500];

	  for(i=0;i<500;i++)
		  compile[i]=0;
	  strcat(compile,"ioc64 -cmd=build -input=tmr_ocl_num_int_el.cl -device=");
#ifdef OPENCL_CPU
	  strcat(compile,"cpu -asm=");
#endif
#ifdef OPENCL_PHI
	  strcat(compile,"co -asm=");
#endif
	  strcat(compile,name);
	  strcat(compile," -bo=\"");
	  strcat(compile,buffer);
	  strcat(compile,"\"");

	  printf("Compile: %s\n",compile);

	  system(compile);
//PHI
	  word=wc(name,"vsubrpd");
	  if(line_count==0)
	  {
		  fprintf(headf, "vsubrpd,");
	  }
	  fprintf(resf,"%d,",word);
	  //printf("add.f64=%d\n",word);
	  word=wc(name,"vfmadd213pd");
	  if(line_count==0)
	  {
		  fprintf(headf, "vfmadd213pd,");
	  }
	  fprintf(resf,"%d,",word);
	  word=wc(name,"vfnmadd231pd");
	  if(line_count==0)
	  {
		  fprintf(headf, "vfnmadd231pd,");
	  }
	  fprintf(resf,"%d,",word);
	  word=wc(name,"vfmadd231pd");
	  if(line_count==0)
	  {
		  fprintf(headf, "vfmadd231pd,");
	  }
	  fprintf(resf,"%d,",word);
	  word=wc(name,"vfmadd132pd");
	  if(line_count==0)
	  {
		  fprintf(headf, "vfmadd132pd,");
	  }
	  fprintf(resf,"%d,",word);
//common
	  word=wc(name,"vmulpd");
	  if(line_count==0)
	  {
		  fprintf(headf, "vmulpd,");
	  }
	  fprintf(resf,"%d,",word);
	  word=wc(name,"vaddpd");
	  if(line_count==0)
	  {
		  fprintf(headf, "vaddpd,");
	  }
	  fprintf(resf,"%d,",word);

//cpu
	  word=wc(name,"vsubsd");
	  if(line_count==0)
	  {
		  fprintf(headf, "vsubsd,");
	  }
	  fprintf(resf,"%d,",word);
	  word=wc(name,"vmulsd");
	  if(line_count==0)
	  {
		  fprintf(headf, "vmulsd,");
	  }
	  fprintf(resf,"%d,",word);
	  word=wc(name,"vaddsd");
	  if(line_count==0)
	  {
		  fprintf(headf, "vaddsd,");
	  }
	  fprintf(resf,"%d,",word);
	  word=wc(name,"vsubpd");
	  if(line_count==0)
	  {
		  fprintf(headf, "vsubpd,");
	  }
	  fprintf(resf,"%d,",word);

  }
  else if(tmr_ocl_check_vendor()==3)
  {

	  //spectre
	  /*
	  v_add_f64
	  v_mul_f64
	  v_div_scale_f64
	  v_rcp_f64
	  v_fma_f64
	  v_div_fmas_f64
	  v_div_fixup_f64
*/

  }



#endif


  platform_struct.list_of_devices[Device_index].program[Kernel_index] = program;
  platform_struct.list_of_devices[Device_index].kernel[Kernel_index] = kernel;
  
  return(1);
}


/*---------------------------------------------------------
  tmr_perform_creation_crs - create crs structure on gpu
---------------------------------------------------------*/

int tmr_perform_creation_crs(int Problem_id,
			   int Nr_int_ent,
			   int* Asse_pos_first_dof_int_ent,
			   int* Assembly_table,
			   int Max_dofs_int_ent
			   //,tmt_ocl_crs_struct *tmv_ocl_crs_struct - currently structure is
			   //                                          global defined
			   ){

  int i;
  
  int field_id = pdr_ctrl_i_params(Problem_id, 3);
  int nreq = apr_get_nreq(field_id);

  // Initialize data structures
  tmv_ocl_crs_struct.ocl_crs_val = NULL;
  tmv_ocl_crs_struct.ocl_crs_col_ind = NULL;
  tmv_ocl_crs_struct.ocl_crs_row_ptr = NULL;
  tmv_ocl_crs_struct.ocl_rhs_val = NULL;

  // ----------- Get platform configuration ----------------
  
  // choose the platform 
  int platform_index = tmr_ocl_get_current_platform_index();
  
  int device_tmc_type = tmr_ocl_get_current_device_type();
  
  // choose device_index
  int device_index = tmr_ocl_get_current_device(); //tmr_ocl_select_device(platform_index, device_tmc_type);

  // choose the context
  cl_context context = tmr_ocl_select_context(platform_index, device_index);

  // choose the command queue
  cl_command_queue command_queue = 
    tmr_ocl_select_command_queue(platform_index, device_index);

  // Create OpenCL memory objects (only create and copy necessary data without kernel assigning

  // For CPU better result is when we don't use pinned buffer
  #ifndef OPENCL_CPU
  #define PINNED
  #endif

  #ifdef OPENCL_HSA

  // TODO

  #else

#ifdef TIME_TEST
  t_begin = omp_get_wtime();
  //t_begin = time_clock();
#endif
    
  
  /* ------------- ALLOCATE MEMEORY FOR ASSEMBLY TABLE -------------------- */
  
  // Calculate the block size
  int block_size=(Max_dofs_int_ent+1)*nreq*Max_dofs_int_ent*nreq;
  
  //Calculate size of the assembly table
  int size_assemble_table_col=block_size*Nr_int_ent;   // compute size of assembly table (size of block for int ent * number of int ent)
  tmv_ocl_crs_struct.ocl_assembly_table_bytes = size_assemble_table_col* sizeof(int);
  tmv_ocl_crs_struct.ocl_asse_pos_first_dof_int_ent_bytes = Nr_int_ent* sizeof(int);
  
  

  // Create assembly table buffer on GPU
  tmv_ocl_crs_struct.ocl_asse_pos_first_dof_int_ent = clCreateBuffer(context, CL_MEM_READ_ONLY, tmv_ocl_crs_struct.ocl_asse_pos_first_dof_int_ent_bytes, NULL, NULL);
  tmv_ocl_crs_struct.ocl_assembly_table = clCreateBuffer(context, CL_MEM_READ_ONLY, tmv_ocl_crs_struct.ocl_assembly_table_bytes, NULL, NULL);
  
  // Write data on GPU
  clEnqueueWriteBuffer(command_queue, tmv_ocl_crs_struct.ocl_asse_pos_first_dof_int_ent, CL_TRUE, 0, tmv_ocl_crs_struct.ocl_asse_pos_first_dof_int_ent_bytes, Asse_pos_first_dof_int_ent, 0, NULL, NULL);
  clEnqueueWriteBuffer(command_queue, tmv_ocl_crs_struct.ocl_assembly_table, CL_TRUE, 0, tmv_ocl_crs_struct.ocl_assembly_table_bytes, Assembly_table, 0, NULL, NULL);
  
  /* ------------- ALLOCATE MEMEORY CRS STRUCTURE -------------------- */
 
  tmv_ocl_crs_struct.ocl_crs_val = clCreateBuffer(context, CL_MEM_READ_WRITE,
      				   tmv_ocl_crs_struct.ocl_crs_val_bytes, NULL, NULL);

  tmv_ocl_crs_struct.ocl_crs_col_ind = clCreateBuffer(context, CL_MEM_READ_ONLY,
      				   tmv_ocl_crs_struct.ocl_crs_col_ind_bytes, NULL, NULL);

  
  tmv_ocl_crs_struct.ocl_crs_row_ptr = clCreateBuffer(context, CL_MEM_READ_ONLY,
      				   tmv_ocl_crs_struct.ocl_crs_row_ptr_bytes, NULL, NULL);

  tmv_ocl_crs_struct.ocl_rhs_val = clCreateBuffer(context, CL_MEM_READ_WRITE,
      				   tmv_ocl_crs_struct.ocl_rhs_bytes, NULL, NULL);


  // Write CRS data column indices and rows pointer

  /* structure must be created before write */
  
  // Write column index table
  clEnqueueWriteBuffer(command_queue, tmv_ocl_crs_struct.ocl_crs_col_ind, CL_TRUE, 0,
  		       tmv_ocl_crs_struct.ocl_crs_col_ind_bytes, tmv_ocl_crs_struct.crs_col_ind,
		       0, NULL, NULL);

  // Write row pointer
  clEnqueueWriteBuffer(command_queue, tmv_ocl_crs_struct.ocl_crs_row_ptr, CL_TRUE, 0,
  		       tmv_ocl_crs_struct.ocl_crs_row_ptr_bytes, tmv_ocl_crs_struct.crs_row_ptr,
		       0, NULL, NULL);

  // Write val buffer on gpu for opencl 1.1
  clEnqueueWriteBuffer(command_queue, tmv_ocl_crs_struct.ocl_crs_val, CL_TRUE, 0,
		       tmv_ocl_crs_struct.ocl_crs_val_bytes, tmv_ocl_crs_struct.crs_val_gpu,
		       0, NULL, NULL);
  // Write val buffer on gpu for opencl 1.1
  clEnqueueWriteBuffer(command_queue, tmv_ocl_crs_struct.ocl_rhs_val, CL_TRUE, 0,
		       tmv_ocl_crs_struct.ocl_rhs_bytes, tmv_ocl_crs_struct.rhs_val_gpu,
		       0, NULL, NULL);
	
  // Write val buffer on gpu for new opencl > 1.1
  // double fill_data = 0.0;
  // clEnqueueFillBuffer(command_queue, tmv_ocl_crs_struct.ocl_crs_val,&fill_data, sizeof(double), 0,
  		      // tmv_ocl_crs_struct.ocl_crs_val_bytes, 0, NULL, NULL); 
  // Write val buffer on gpu for new opencl > 1.1
  // double fill_data = 0.0;
  // clEnqueueFillBuffer(command_queue, tmv_ocl_crs_struct.ocl_rhs_val,&fill_data, sizeof(double), 0,
  		      // tmv_ocl_crs_struct.ocl_rhs_bytes, 0, NULL, NULL); 

      

#ifdef TIME_TEST
  clFinish(command_queue);
  t_end = omp_get_wtime();
  //t_end = time_clock();
  //total_time+=t_end-t_begin;
  printf("CRS --> EXECUTION TIME: create and initialize assembly tables and crs on GPU: %.12lf\n",
	 t_end-t_begin);
#endif

  // TO DO ...
    
  #endif //ifndef HSA


  return(0);
  }

/*---------------------------------------------------------
  tmr_send_data_to_gpu_crs - send 
---------------------------------------------------------*/
int tmr_prepare_final_crs() {

#ifdef TIME_TEST
  double kernel_execution_time=0.0;
#endif

#ifdef COUNT_OPER
  SCALAR count_oper[3];
  cl_mem ocl_count_oper;
#endif

  // ----------- Get platform configuration ----------------
  

  int i;

  // choose the platform 
  int platform_index = tmr_ocl_get_current_platform_index();
  
  int device_tmc_type = tmr_ocl_get_current_device_type();
  
  // choose device_index
  int device_index = tmr_ocl_get_current_device(); //tmr_ocl_select_device(platform_index, device_tmc_type);

  // OpenCL device characteristics stored in data structure
  tmt_ocl_device_struct device_struct = 
    tmv_ocl_struct.list_of_platforms[platform_index].list_of_devices[device_index];
  int max_num_comp_units = device_struct.max_num_comp_units;
  int max_work_group_size = device_struct.max_work_group_size;

  // ASSUMED WORK_GROUP_SIZE
  int work_group_size = WORK_GROUP_SIZE;

  // usually for GPUs it is good to maximize the number of work_groups and threads
  int nr_work_groups = NR_WORK_GROUPS_PER_COMP_UNIT * max_num_comp_units;
  int nr_threads = nr_work_groups * work_group_size;

  // choose the context
  cl_context context = tmr_ocl_select_context(platform_index, device_index);

  // choose the command queue
  cl_command_queue command_queue = 
    tmr_ocl_select_command_queue(platform_index, device_index);

  // THERE MAY BE SEVERAL KERNELS FOR THE CODE, THE INDEX IN OPENCL DATA 
  // STRUCTURE FOR THE NUMERICAL INTEGRATION KERNEL IS ASSIGNED IN (tmd_ocl/tmh_ocl.h)
  int kernel_index=TMC_OCL_KERNEL_CRS_FINALIZE;

  // !!!!! FOR THE TIME BEING WE ALWAYS READ AND COMPILE tmr_ocl_num_int_el_and_assembling.cl
  cl_kernel kernel = tmr_ocl_select_kernel(platform_index, device_index, kernel_index); 

  #ifndef NR_EXEC_PARAMS
  #define NR_EXEC_PARAMS 16  // size of array with execution parameters
  #endif

  // WE ALLOCATE OPENCL BUFFERS FOR nr_elems_per_kernel ELEMENTS
  // this should be good for all kernel invocations and all colors !!!
  // (actual kernel calculations will be for nr_elems_this_kercall)
  
  cl_int retval = CL_SUCCESS;
  
#ifdef TIME_TEST
  t_begin = omp_get_wtime();
  //t_begin = time_clock();
#endif


  // there are NR_EXEC_PARAMS execution parameters (ALLWAYS CHECK!!!)
  
  const size_t execution_parameters_dev_bytes = NR_EXEC_PARAMS * sizeof(int);
  
  cl_mem execution_parameters_dev;
  
  execution_parameters_dev = clCreateBuffer(context, CL_MEM_READ_ONLY,
					    execution_parameters_dev_bytes, NULL, NULL);
  retval |= clSetKernelArg(kernel, 0, sizeof(cl_mem), (void *) &execution_parameters_dev);
  
  // Create OpenCL memory objects (only create and copy necessary data without kernel assigning

  // For CPU better result is when we don't use pinned buffer
  #ifndef OPENCL_CPU
  #define PINNED
  #endif

  #ifdef OPENCL_HSA

  // TODO

  #else

  /* ------------- ALLOCATE MEMEORY CRS STRUCTURE -------------------- */
 
  tmv_ocl_crs_struct.ocl_crs_external_val = clCreateBuffer(context, CL_MEM_READ_ONLY,
      				   tmv_ocl_crs_struct.ocl_crs_val_bytes, NULL, NULL);

  tmv_ocl_crs_struct.ocl_rhs_external_val = clCreateBuffer(context, CL_MEM_READ_ONLY,
      				   tmv_ocl_crs_struct.ocl_rhs_bytes, NULL, NULL);

  /* structure must be created before write */

  clEnqueueWriteBuffer(command_queue, tmv_ocl_crs_struct.ocl_crs_external_val, CL_TRUE, 0,
		       tmv_ocl_crs_struct.ocl_crs_val_bytes, tmv_ocl_crs_struct.crs_val_cpu,
		       0, NULL, NULL);

  clEnqueueWriteBuffer(command_queue, tmv_ocl_crs_struct.ocl_rhs_external_val, CL_TRUE, 0,
		       tmv_ocl_crs_struct.ocl_rhs_bytes, tmv_ocl_crs_struct.rhs_val_cpu,
		       0, NULL, NULL);

  #ifdef COUNT_OPER
  ocl_count_oper = clCreateBuffer(context, CL_MEM_WRITE_ONLY,3*sizeof(SCALAR), NULL, NULL);
  #endif

  /* ------------- BIND DATA -------------------- */
  
  // Bind externel RHS vector 
  retval |= clSetKernelArg(kernel, 1, sizeof(cl_mem), (void *) &tmv_ocl_crs_struct.ocl_rhs_external_val);

  // Bind externel CRS structure
  retval |= clSetKernelArg(kernel, 2, sizeof(cl_mem), (void *) &tmv_ocl_crs_struct.ocl_crs_external_val);

  
  // Bind RHS vector
  retval |= clSetKernelArg(kernel, 3, sizeof(cl_mem), (void *) &tmv_ocl_crs_struct.ocl_rhs_val);

  // Bind CRS structure
  retval |= clSetKernelArg(kernel, 4, sizeof(cl_mem), (void *) &tmv_ocl_crs_struct.ocl_crs_val);

  #ifdef COUNT_OPER
  retval |= clSetKernelArg(kernel, 5, sizeof(cl_mem), (void *) &ocl_count_oper);
  #endif

  /* ------------- PERFORM KERNEL -------------------- */

  int execution_parameters_host[NR_EXEC_PARAMS];
  execution_parameters_host[0] = tmv_ocl_crs_struct.Nrdof;
  execution_parameters_host[1] = tmv_ocl_crs_struct.Nnz;
  execution_parameters_host[2] = tmv_ocl_crs_struct.Nrdof/(nr_work_groups*work_group_size);
  execution_parameters_host[3] = tmv_ocl_crs_struct.Nrdof%(nr_work_groups*work_group_size);
  execution_parameters_host[4] = tmv_ocl_crs_struct.Nnz/(nr_work_groups*work_group_size);
  execution_parameters_host[5] = tmv_ocl_crs_struct.Nnz%(nr_work_groups*work_group_size);

  // send data to device - non-blocking call (CL_FALSE)
  clEnqueueWriteBuffer(command_queue, execution_parameters_dev, CL_TRUE, 0,
		       execution_parameters_dev_bytes, execution_parameters_host,
		       0, NULL, NULL);


#ifdef TIME_TEST
  clFinish(command_queue);
  t_end = omp_get_wtime();
  //t_end = time_clock();
  
  total_time+=t_end-t_begin;
  printf("CRS --> EXECUTION TIME: copying input crs buffer: %.12lf, total time: %.12lf\n",
	 t_end-t_begin, total_time);
#endif

  // FINALLY EXECUTE THE KERNEL
  // set kernel invocation parameters 
  size_t globalWorkSize[1] = { 0 };
  size_t localWorkSize[1] = { 0 };
  globalWorkSize[0] = nr_threads ;
  localWorkSize[0] = work_group_size ;
    
  cl_event ndrEvt;
    
#ifdef TIME_TEST
  t_begin = omp_get_wtime();
  //t_begin = time_clock();
#endif

  // Queue the kernel up for execution
  retval = clEnqueueNDRangeKernel(command_queue, kernel, 1, NULL,
				  globalWorkSize, localWorkSize,
				  0, NULL,  &ndrEvt);

  #ifdef TIME_TEST
    clWaitForEvents(1, &ndrEvt);
    clFinish(command_queue);

    t_end = omp_get_wtime();
    //t_end = time_clock();
    
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
    
    // double kernel_execution_time=0.0;
    kernel_execution_time = t_end-t_begin;
    printf("CRS --> EXECUTION TIME: executing addition kernel: %.12lf (internal %.12lf)\n",
	   kernel_execution_time, time_internal);
    
    total_time += t_end-t_begin;
    t_begin = omp_get_wtime();
    //t_begin = time_clock();
#endif


#ifdef GPU_DEBUG
    memset(tmv_ocl_crs_struct.crs_val_cpu,0,tmv_ocl_crs_struct.ocl_crs_val_bytes);
    memset(tmv_ocl_crs_struct.rhs_val_cpu,0,tmv_ocl_crs_struct.ocl_rhs_bytes);
#endif

    clEnqueueReadBuffer(command_queue, tmv_ocl_crs_struct.ocl_crs_val, CL_TRUE, 0, tmv_ocl_crs_struct.ocl_crs_val_bytes, tmv_ocl_crs_struct.crs_val_cpu, 0, NULL, NULL);
    clEnqueueReadBuffer(command_queue, tmv_ocl_crs_struct.ocl_rhs_val, CL_TRUE, 0, tmv_ocl_crs_struct.ocl_rhs_bytes, tmv_ocl_crs_struct.rhs_val_cpu, 0, NULL, NULL);
  
/*
  for(i=0; i<tmv_ocl_crs_struct.Nnz; i++) {
    tmv_ocl_crs_struct.crs_val_cpu[i]= tmv_ocl_crs_struct.crs_val_gpu[i];  
  }

  for(i=0; i<(tmv_ocl_crs_struct.ocl_rhs_bytes/sizeof(double)); i++) {
    tmv_ocl_crs_struct.rhs_val_cpu[i]=tmv_ocl_crs_struct.rhs_val_gpu[i];
  }
*/

  #ifdef COUNT_OPER
    clEnqueueReadBuffer(command_queue, ocl_count_oper, CL_TRUE, 0, 3*sizeof(SCALAR), count_oper, 0, NULL, NULL);
  #endif  
  
#ifdef TIME_TEST
  clFinish(command_queue);
  t_end = omp_get_wtime();
  //t_end = time_clock();
  total_time+=t_end-t_begin;
  printf("CRS --> EXECUTION TIME: copying output crs buffer: %.12lf, total time: %.12lf\n",
	 t_end-t_begin, total_time);
#endif

  #ifdef COUNT_OPER
  printf("CRS --> Arthmetic operations: %.0lf, Local mem access: %.0lf, Global mem access: %.0lf\n",count_oper[0],count_oper[1],count_oper[2]);
  clReleaseMemObject(ocl_count_oper);
  #endif  
  
  #endif
  return(0);

}


/*---------------------------------------------------------
  tmr_cleanup_crs - free tmr crs structure
---------------------------------------------------------*/
int tmr_cleanup_crs() {

  // Free CRS data
  if(tmv_ocl_crs_struct.crs_val_gpu != NULL) {
    free(tmv_ocl_crs_struct.crs_val_gpu);
    tmv_ocl_crs_struct.crs_val_gpu = NULL;
  }
  if(tmv_ocl_crs_struct.rhs_val_gpu != NULL) {
    free(tmv_ocl_crs_struct.rhs_val_gpu);
    tmv_ocl_crs_struct.rhs_val_gpu = NULL;
  }
  tmv_ocl_crs_struct.crs_col_ind = NULL;
  tmv_ocl_crs_struct.crs_row_ptr = NULL;

  // Release GPU memory objects
  clReleaseMemObject(tmv_ocl_crs_struct.ocl_asse_pos_first_dof_int_ent);
  clReleaseMemObject(tmv_ocl_crs_struct.ocl_assembly_table);
  clReleaseMemObject(tmv_ocl_crs_struct.ocl_crs_val);
  clReleaseMemObject(tmv_ocl_crs_struct.ocl_crs_external_val);
  clReleaseMemObject(tmv_ocl_crs_struct.ocl_crs_col_ind);
  clReleaseMemObject(tmv_ocl_crs_struct.ocl_crs_row_ptr);
  clReleaseMemObject(tmv_ocl_crs_struct.ocl_rhs_val);
  clReleaseMemObject(tmv_ocl_crs_struct.ocl_rhs_external_val);

  // Cleaning
  tmv_ocl_crs_struct.Nnz = 0;
  tmv_ocl_crs_struct.Nrdof = 0;
  tmv_ocl_crs_struct.ocl_crs_val_bytes = 0;
  tmv_ocl_crs_struct.ocl_crs_col_ind_bytes = 0;
  tmv_ocl_crs_struct.ocl_crs_row_ptr_bytes = 0;

  return(0);
}


/*---------------------------------------------------------
  tmr_get_crs_structure - get crs structure from lsd
---------------------------------------------------------*/
int tmr_get_crs_structure(
  int Nrdof_glob,
  int nnz,
  int* crs_col_ind,
  int* crs_row_ptr,
  double* crs_val,
  double* rhs
 ){
   
  
  tmv_ocl_crs_struct.crs_col_ind = crs_col_ind;
  tmv_ocl_crs_struct.crs_row_ptr = crs_row_ptr;
  tmv_ocl_crs_struct.crs_val_cpu = crs_val;
  tmv_ocl_crs_struct.rhs_val_cpu = rhs;
  tmv_ocl_crs_struct.crs_val_gpu = (SCALAR*) malloc (nnz*sizeof(SCALAR));
  tmv_ocl_crs_struct.rhs_val_gpu = (SCALAR*) malloc (Nrdof_glob*sizeof(SCALAR));


  // For opencl 1.1
  memset(tmv_ocl_crs_struct.crs_val_gpu,0.0,nnz*sizeof(SCALAR));
  memset(tmv_ocl_crs_struct.rhs_val_gpu,0.0,Nrdof_glob*sizeof(SCALAR));
  
  tmv_ocl_crs_struct.Nnz = nnz;
  tmv_ocl_crs_struct.Nrdof = Nrdof_glob;
  tmv_ocl_crs_struct.ocl_crs_val_bytes = nnz*sizeof(SCALAR);
  tmv_ocl_crs_struct.ocl_crs_col_ind_bytes = (nnz+1)*sizeof(int);
  tmv_ocl_crs_struct.ocl_crs_row_ptr_bytes = (Nrdof_glob+1)*sizeof(int);
  tmv_ocl_crs_struct.ocl_rhs_bytes = Nrdof_glob*sizeof(SCALAR);
  
   return(0);
 }
  
  
