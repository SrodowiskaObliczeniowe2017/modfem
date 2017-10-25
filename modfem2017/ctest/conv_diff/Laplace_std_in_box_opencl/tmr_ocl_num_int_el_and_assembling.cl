/* QSS with ASSEMBLING kernel
 *
 * OpenCL Kernel for following probems:
 * - LAPLACE
 * - TEST_NUMINT
 * - HEAT
 * Other notes:
 * Loop model: QSS
 * Implemented features:
 * - CRS ASSEMBLING
 * - simple assembling to local matrix
 * To do:
 * - TEST HEAT PROBLEM
 * - TEST DIFFERENT CONFIGURATIONS OF SHARED / L1 MEMORY
 * - LOOP UNROLLING = REGISTER BLOCKING FOR INNERMOST LOOPS
 *
 * 02.2013 - Krzysztof Banas, initial version
 * 10.2016 - Jan Bielanski, Kazimierz Chlon - assembling into crs on gpu
 * 10.2016 - Jan Bielanski - support heat problem
 */

#if defined(cl_amd_fp64)
  #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#elif defined(cl_khr_fp64)
  #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#else
  #error "Double precision floating point not supported by OpenCL implementation."
#endif

//----------------------------------------------------
// TWO MASTER SWITCHES (float<->double, work_group_size)
//#define FLOAT
#ifdef FLOAT
  #define SCALAR float
  #define zero 0.0f
  #define one 1.0f
  #define two 2.0f
  #define half 0.5f
  #define one_fourth 0.25f
  #define one_sixth (0.16666666667f)
#else
  #define SCALAR double
  #define zero 0.0
  #define one 1.0
  #define two 2.0
  #define half 0.5
  #define one_fourth 0.25
  #define one_sixth (0.166666666666666667)
#endif

//#define WORK_GROUP_SIZE 16 // for XEON_PHI
#define WORK_GROUP_SIZE 64 // for GPUs
//#define WORK_GROUP_SIZE 8 // for CPUs
//----------------------------------------------------


//----------------------------------------------------
// SEVERAL PROBLEM, ELEMENT AND APPROXIMATION DEPENDENT SWITCHES 
// (nr_exec_params, nreq, num_shap, num_gauss, num_geo_dofs
#define NR_EXEC_PARAMS 16  // size of array with execution parameters
// here: the smallest work-group for reading data is selected
// exec_params are read from global to shared memory and used when needed
// if shared memory resources are scarce this can be reduced

// FOR SCALAR PROBLEMS !!!!!!!!!!!!!!!!!!!!!
#define nreq 1
// FOR NS_SUPG PROBLEM !!!!!!!!!!!!!!!!!!!!!
//#define nreq 4

// FOR LINEAR PRISMS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#define num_shap 6
#define num_gauss 6
#define num_geo_dofs 6
#define weight_linear_prism (one_sixth)
#define weight_gauss weight_linear_prism

// FOR LINEAR TETRAHEDRA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//#define num_shap 4
//#define num_gauss 4
//#define num_geo_dofs 4
//#define weight_linear_tetra (one_fourth*one_sixth)
//#define weight_gauss weight_linear_tetra

#define num_dofs (num_shap*nreq)
#define EL_GEO_DAT_SIZE (3*num_geo_dofs)

// J_AND_DETJ_SIZE=10 - for NOJAC variants
//#define J_AND_DETJ_SIZE 10

// the number of coefficients sent for elements
// either coefficients constant for the whole element
// or different for every integration point
#define LAPLACE
//#define TEST_NUMINT
//#define HEAT
#ifdef LAPLACE
  #define NR_COEFFS_SENT_PER_ELEMENT 0
  #define NR_COEFFS_SENT_PER_INT_POINT 1
  #define NR_COEFFS_IN_SM_CALCULATIONS num_gauss
#elif defined TEST_NUMINT
  #define NR_COEFFS_SENT_PER_ELEMENT 20
  #define NR_COEFFS_SENT_PER_INT_POINT 0
  #define NR_COEFFS_IN_SM_CALCULATIONS 20
#elif defined HEAT
  #define NR_COEFFS_SENT_PER_ELEMENT (num_dofs+num_dofs+1)
  #define NR_COEFFS_SENT_PER_INT_POINT 5
  #define NR_COEFFS_IN_SM_CALCULATIONS 20
  #define NR_PDE_COEFFS_READ (num_dofs+1 + NR_COEFFS_SENT_PER_INT_POINT*num_gauss)
#endif

#define NR_PDE_COEFFS_SENT (NR_COEFFS_SENT_PER_ELEMENT + NR_COEFFS_SENT_PER_INT_POINT*num_gauss)
//----------------------------------------------------


//----------------------------------------------------
// SWITCHES FOR DIFFERENT OPTIMIZATION OPTIONS !!!!!!!!!!!!!!!!!!!!!!!!

//gpu assembling - enable assembling on gpu
#define ASSEMBLING

//load vector computing - not defined only for some tests
#define LOAD_VEC_COMP

// COAL_READ - both PDE_COEFF and GEO_DATA are read in a coalesced way (using a single workspace)
// coalesced reading may be good for GPUs (requires large workspace and several barriers)
// (the workspace can be further used by GEO_DAT or PDE_COEFF or SHAPE_FUN or STIFF_MAT !!!)
#define COAL_READ

// coalesced writing requires host code to adapt to the order of data!!!
// coalesced writing should be switched on for GPUs
#define COAL_WRITE

// COMPUTE_ALL_SHAPE_FUN_DER - to compute all shape functions and their derivatives
//                         before entering the loops over shape functions
#define COMPUTE_ALL_SHAPE_FUN_DER

//#define COUNT_OPER

// ********** THE SET OF OPTIONS FOR USING THE WORKSPACE IN SHARED MEMORY
// ********** AT MOST 1 SHOULD BE DEFINED !!!!!!!!!!!!!!!
// USE_WORKSPACE_FOR_PDE_COEFF - to use shared memory for pde_coeff during SM calculations
// otherwise - registers are used
//#define USE_WORKSPACE_FOR_PDE_COEFF

// USE_WORKSPACE_FOR_GEO_DATA - to use shared memory for geo_dat during SM calculations
// otherwise - registers are used
//#define USE_WORKSPACE_FOR_GEO_DATA


// USE_WORKSPACE_FOR_SHAPE_FUN - to use shared memory for shape functions 
//                               and their derivatives during SM calculations
// otherwise - registers are used
//#define USE_WORKSPACE_FOR_SHAPE_FUN


// USE_WORKSPACE_FOR_STIFF_MAT - to use shared memory for SM during SM calculations
// otherwise - registers are used
//#define USE_WORKSPACE_FOR_STIFF_MAT

// FOR EACH ARCHITECTURE PADDING SHOULD BE TESTED TO DETECT SHARED MEMORY BANK CONFLICTS
#define WORKSPACE_PADDING 0
#define PADDING WORKSPACE_PADDING

// THE END OF: SWITCHES FOR DIFFERENT OPTIMIZATION OPTIONS
//----------------------------------------------------

/* AUXILLIARY FUNCTIONS FOR COMMON PART OF CODE */

// Initialize PDE_COEFF
void initialize_pde_coeff(__local SCALAR *workspace, __private SCALAR *pde_coeff, __private int *offset)
{
  #ifdef LAPLACE
    pde_coeff[0]=workspace[(*offset)+0];
    pde_coeff[1]=workspace[(*offset)+1];
    pde_coeff[2]=workspace[(*offset)+2];
    pde_coeff[3]=workspace[(*offset)+3];
    pde_coeff[4]=workspace[(*offset)+4];
    pde_coeff[5]=workspace[(*offset)+5];
  #elif defined(TEST_NUMINT)
    pde_coeff[0]=workspace[(*offset)+0]; //coeff00
    pde_coeff[1]=workspace[(*offset)+1]; //coeff01
    pde_coeff[2]=workspace[(*offset)+2]; //coeff02
    pde_coeff[3]=workspace[(*offset)+3]; //coeff10
    pde_coeff[4]=workspace[(*offset)+4]; //coeff11
    pde_coeff[5]=workspace[(*offset)+5]; //coeff12
    pde_coeff[6]=workspace[(*offset)+6]; //coeff20
    pde_coeff[7]=workspace[(*offset)+7]; //coeff21
    pde_coeff[8]=workspace[(*offset)+8]; //coeff22
    pde_coeff[9]=workspace[(*offset)+9]; //coeff30
    pde_coeff[10]=workspace[(*offset)+10]; //coeff31
    pde_coeff[11]=workspace[(*offset)+11]; //coeff32
    pde_coeff[12]=workspace[(*offset)+12]; //coeff03
    pde_coeff[13]=workspace[(*offset)+13]; //coeff13
    pde_coeff[14]=workspace[(*offset)+14]; //coeff23
    pde_coeff[15]=workspace[(*offset)+15]; //coeff33
    pde_coeff[16]=workspace[(*offset)+16]; //coeff04
    pde_coeff[17]=workspace[(*offset)+17]; //coeff14
    pde_coeff[18]=workspace[(*offset)+18]; //coeff24
    pde_coeff[19]=workspace[(*offset)+19]; //coeff34
  #elif defined(HEAT)
    pde_coeff[0]=zero; //coeff00
    pde_coeff[1]=zero; //coeff01
    pde_coeff[2]=zero; //coeff02
    pde_coeff[3]=zero; //coeff10
    pde_coeff[4]=zero; //coeff11
    pde_coeff[5]=zero; //coeff12
    pde_coeff[6]=zero; //coeff20
    pde_coeff[7]=zero; //coeff21
    pde_coeff[8]=zero; //coeff22
    pde_coeff[9]=zero; //coeff30
    pde_coeff[10]=zero; //coeff31
    pde_coeff[11]=zero; //coeff32
    pde_coeff[12]=zero; //coeff03
    pde_coeff[13]=zero; //coeff13
    pde_coeff[14]=zero; //coeff23
    pde_coeff[15]=zero; //coeff33
    pde_coeff[16]=zero; //coeff04
    pde_coeff[17]=zero; //coeff14
    pde_coeff[18]=zero; //coeff24
    pde_coeff[19]=zero; //coeff34
  #endif
}
/* THE END: AUXILLIARY FUNCTIONS FOR COMMON PART OF CODE */
//----------------------------------------------------


/* KERNEL WITH NUMERICAL INTEGRATION IN ELEMENTS AND ASSEBLING OF STIFFNESS MATRIX */
kernel void tmr_ocl_num_int_el(
  // execution_parameters can be read directly from constant memory, assuming it is cached and
  // further accesses are realized from cache
  __constant int* execution_parameters,
  //__global int* execution_parameters,
  // gauss data can be read directly from constant memory, assuming it is cached and
  // further accesses are realized from cache
  __constant SCALAR* gauss_dat, // integration points data of elements having given p
  //__global SCALAR* gauss_dat, // integration points data of elements having given p
  // shape function values can be read directly from constant memory, assuming it is cached and
  // further accesses are realized from cache
  __constant SCALAR* shpfun_ref, // shape functions on a reference element
  //__global SCALAR* shpfun_ref, // shape functions on a reference element

  __global SCALAR* el_data_in, // data for integration of NR_ELEMS_THIS_KERCALL elements

  //__global SCALAR* el_data_in_geo, // geo data for integration of NR_ELEMS_THIS_KERCALL elements
  //__global SCALAR* el_data_in_coeff, // coeff data for integration of elements

#ifdef ASSEMBLING  
  __global int* assembly_position_table,  //in: assembly table
  __global int* assembly_table,  //in: assembly table
  
  __global SCALAR* rhs,  //out: 
  __global SCALAR* crs  //out:

  // __global int *crs_col_ind, //in: not used in this kernel
  // __global int *crs_row_ptr //in: not used in this kernel

#else
  __global SCALAR* stiff_mat_out // out: result of integration of NR_ELEMS_THIS_KERCALL elements
#endif

#ifdef COUNT_OPER
  , // buffer give access to kernel profiling data such as number of arthmetic operations and shared/global memory access
  __global SCALAR* kernel_profile_data_out

#endif
  
){

#ifdef COUNT_OPER
  SCALAR nr_oper=0.0;
  SCALAR nr_access_shared=0.0;
  SCALAR nr_global_access=0.0;
#endif

  const int group_id = get_group_id(0);
  const int thread_id = get_local_id(0);
  //const int work_group_size = get_local_size(0);
  const int nr_work_groups = get_num_groups(0);


//----------------------------------------------------
// DEFINITIONS DEPENDENT ON OPTIMIZATION OPTIONS

#ifdef USE_WORKSPACE_FOR_PDE_COEFF
  #ifdef HEAT
    #define WORKSPACE_SIZE_FOR_PDE_COEFF (NR_COEFFS_IN_SM_CALCULATIONS+NR_PDE_COEFFS_READ)
  #else
    #define WORKSPACE_SIZE_FOR_PDE_COEFF (NR_COEFFS_IN_SM_CALCULATIONS)
  #endif
#else
  #define WORKSPACE_SIZE_FOR_PDE_COEFF 0
  #ifdef HEAT
    SCALAR pde_coeff[NR_COEFFS_IN_SM_CALCULATIONS+NR_PDE_COEFFS_READ];
  #elif defined LAPLACE
    SCALAR pde_coeff[num_shap];
  #else
    SCALAR pde_coeff[NR_COEFFS_IN_SM_CALCULATIONS];
  #endif
#endif

#ifdef USE_WORKSPACE_FOR_GEO_DATA
  #define WORKSPACE_SIZE_FOR_GEO_DATA (3*num_geo_dofs)
#else
  #define WORKSPACE_SIZE_FOR_GEO_DATA 0
  SCALAR geo_dat[3*num_geo_dofs];
#endif

#ifdef USE_WORKSPACE_FOR_SHAPE_FUN
  #define WORKSPACE_SIZE_FOR_SHAPE_FUN (3*num_shap)
#else
  #define WORKSPACE_SIZE_FOR_SHAPE_FUN 0
  #ifdef COMPUTE_ALL_SHAPE_FUN_DER
  SCALAR tab_fun_u_derx[num_shap];
  SCALAR tab_fun_u_dery[num_shap];
  SCALAR tab_fun_u_derz[num_shap];
  #endif
#endif

#ifdef USE_WORKSPACE_FOR_STIFF_MAT
  #define WORKSPACE_SIZE_FOR_STIFF_MAT (num_dofs*(num_dofs+1))
#else
  #define WORKSPACE_SIZE_FOR_STIFF_MAT 0
  SCALAR stiff_mat[num_dofs*num_dofs];
  SCALAR load_vec[num_dofs];
#endif

#ifdef COAL_READ
  #if (NR_PDE_COEFFS_SENT) > (3*num_geo_dofs)
    #define WORKSPACE_SIZE_FOR_READING NR_PDE_COEFFS_SENT
  #else
    #define WORKSPACE_SIZE_FOR_READING (3*num_geo_dofs)
  #endif
#else
  #define WORKSPACE_SIZE_FOR_READING 0
#endif

#define TMP_SIZE (WORKSPACE_SIZE_FOR_PDE_COEFF+WORKSPACE_SIZE_FOR_GEO_DATA+WORKSPACE_SIZE_FOR_SHAPE_FUN+WORKSPACE_SIZE_FOR_STIFF_MAT)

#if (TMP_SIZE) > (WORKSPACE_SIZE_FOR_READING)
  #define WORKSPACE_SIZE ((TMP_SIZE+WORKSPACE_PADDING)*WORK_GROUP_SIZE)
#else
  #define WORKSPACE_SIZE ((WORKSPACE_SIZE_FOR_READING+WORKSPACE_PADDING)*WORK_GROUP_SIZE)
#endif

#if (WORKSPACE_SIZE) > 0
  __local SCALAR workspace[WORKSPACE_SIZE]; //
#endif

// THE END OF: DEFINITIONS DEPENDENT ON OPTIMIZATION OPTIONS
//----------------------------------------------------


  // ASSUMPTION: one element = one thread

  int nr_elems_per_thread = execution_parameters[0];
  int nr_elems_this_kercall = execution_parameters[1];
#ifdef ASSEMBLING  
  int first_elem_index = execution_parameters[2]; // The index of the first element in assebly table
#endif

  int ielem;
  int offset;


//-------------------------------------------------------------
//******************* loop over elements processed by a thread *********************
  for(ielem = 0; ielem < nr_elems_per_thread; ielem++){

    int element_index = group_id * nr_elems_per_thread * WORK_GROUP_SIZE +
                                                 ielem * WORK_GROUP_SIZE +
                                                               thread_id ;
    int i;


//-------------------------------------------------------------
// ******************* READING INPUT DATA *********************

#ifdef COAL_READ

  #ifdef USE_WORKSPACE_FOR_GEO_DATA // workspace is used for GEO_DATA hence we read PDE_COEFF
                                    // and immediately rewrite them to registers

    offset = nr_elems_this_kercall * EL_GEO_DAT_SIZE + (element_index - thread_id) * NR_PDE_COEFFS_SENT;
    // !!! OLD BACKUP
    //offset = (element_index - thread_id) * NR_PDE_COEFFS_SENT;
    // !!! OLD BACKUP

    barrier(CLK_LOCAL_MEM_FENCE); // I don't know why but without barrier here, one per ten runs gives bad results

    // TRY TO UNROLL THIS LOOP TO INCREASE MEMORY PARALLELISM !!!
    #ifdef HEAT
    workspace[i*WORK_GROUP_SIZE+thread_id] = zero;
    #else
    for(i=0;i<NR_PDE_COEFFS_SENT;i++) {

      workspace[i*WORK_GROUP_SIZE+thread_id] =
	el_data_in[offset+i*WORK_GROUP_SIZE+thread_id];

      // !!! OLD BACKUP
      //workspace[i*WORK_GROUP_SIZE+thread_id] =
      //el_data_in_coeff[offset+i*WORK_GROUP_SIZE+thread_id];
      // !!! OLD BACKUP

    }

    #ifdef COUNT_OPER
    nr_global_access += NR_PDE_COEFFS_SENT; // we neglect shared memory accesses at this stage
    #endif
    #endif

    barrier(CLK_LOCAL_MEM_FENCE); // !!!!!!!!!!!!!!!!!!!!!!

    offset=thread_id*(NR_PDE_COEFFS_SENT);

    // Initialize PDE_COEFF values
    initialize_pde_coeff(workspace,pde_coeff,&offset);

    #ifdef COUNT_OPER
      #ifdef HEAT
         nr_access_shared += 0;
      #else
         nr_access_shared += NR_PDE_COEFFS_SENT; // we count accesses because of the barriers
      #endif
    #endif

    barrier(CLK_LOCAL_MEM_FENCE); // !!!!!!!!!!!!!!!!!!!!!!

    // after rewriting PDE_COEFF to registers we read GEO_DATA to workspace
    offset = (element_index-thread_id)*(EL_GEO_DAT_SIZE);

    // TRY TO UNROLL THIS LOOP TO INCREASE MEMORY PARALLELISM !!!
    for(i = 0; i < EL_GEO_DAT_SIZE; i++){

      workspace[i*WORK_GROUP_SIZE+thread_id] =
	el_data_in[offset+i*WORK_GROUP_SIZE+thread_id];
      
      // !!! OLD BACKUP
      //workspace[i*WORK_GROUP_SIZE+thread_id] = 
      //el_data_in_geo[offset+i*WORK_GROUP_SIZE+thread_id];
      // !!! OLD BACKUP

    }

    #ifdef COUNT_OPER
    nr_global_access += EL_GEO_DAT_SIZE; // we neglect shared memory accesses at this stage
    #endif


  #else // if workspace not used for geo_data
        // (hence used for pde_coeff or something else)

    // first we read geo data to workspace 
    offset = (element_index-thread_id)*(EL_GEO_DAT_SIZE);

    barrier(CLK_LOCAL_MEM_FENCE); // I don't know why but without barrier here, one per ten runs gives bad results

    // TRY TO UNROLL THIS LOOP TO INCREASE MEMORY PARALLELISM !!!
    for(i = 0; i < EL_GEO_DAT_SIZE; i++){

      workspace[i*WORK_GROUP_SIZE+thread_id] =
	el_data_in[offset+i*WORK_GROUP_SIZE+thread_id];  
      
      // !!! OLD BACKUP
      //workspace[i*WORK_GROUP_SIZE+thread_id] = 	
      //el_data_geo[offset+i*WORK_GROUP_SIZE+thread_id];
      // !!! OLD BACKUP
  
    }

    #ifdef COUNT_OPER
    nr_global_access += EL_GEO_DAT_SIZE; // we neglect shared memory accesses at this stage
    #endif

    barrier(CLK_LOCAL_MEM_FENCE); // !!!!!!!!!!!!!!!!!!!!!!

    // we rewrite geo_data to registers
    offset=thread_id*EL_GEO_DAT_SIZE;
    
    for(i=0;i<num_geo_dofs;i++){  
      
      geo_dat[3*i] = workspace[offset+3*i];  //node coor
      geo_dat[3*i+1] = workspace[offset+3*i+1];
      geo_dat[3*i+2] = workspace[offset+3*i+2];
      
    }
    
    #ifdef COUNT_OPER
    nr_access_shared += 3*num_geo_dofs; // we count accesses because of the barriers
    #endif

    barrier(CLK_LOCAL_MEM_FENCE); // !!!!!!!!!!!!!!!!!!!!!!


    // after rewriting GEO_DATA to registers we read  PDE_COEFF to workspace
    offset = nr_elems_this_kercall * EL_GEO_DAT_SIZE + (element_index - thread_id) * NR_PDE_COEFFS_SENT;

    // !!! OLD BACKUP
    //offset = (element_index - thread_id) * NR_PDE_COEFFS_SENT;
    // !!! OLD BACKUP

    // TRY TO UNROLL THIS LOOP TO INCREASE MEMORY PARALLELISM !!!
    #ifdef HEAT
    workspace[i*WORK_GROUP_SIZE+thread_id] = zero;
    #else
    for(i=0;i<NR_PDE_COEFFS_SENT;i++) {

      workspace[i*WORK_GROUP_SIZE+thread_id] =
	el_data_in[offset+i*WORK_GROUP_SIZE+thread_id];
      
      // !!! OLD BACKUP
      //workspace[i*WORK_GROUP_SIZE+thread_id] = 
      //el_data_in_coeff[offset+i*WORK_GROUP_SIZE+thread_id];
      // !!! OLD BACKUP

    }

    #ifdef COUNT_OPER
    nr_global_access += NR_PDE_COEFFS_SENT; // we neglect shared memory accesses at this stage
    #endif
    #endif

    // if we do not leave pde coeffs in workspace we rewrite them to registers
    #ifndef USE_WORKSPACE_FOR_PDE_COEFF

    barrier(CLK_LOCAL_MEM_FENCE); // !!!!!!!!!!!!!!!!!!!!!!

    offset=thread_id*(NR_PDE_COEFFS_SENT);

    // Initialize PDE_COEFF values
    initialize_pde_coeff(workspace,pde_coeff,&offset);

    #ifdef COUNT_OPER
      #ifdef HEAT
    nr_access_shared += 0;
      #else
        nr_access_shared += NR_PDE_COEFFS_SENT; // we count accesses because of the barriers
      #endif
    #endif

    #endif // end if we do not leave pde coeffs in workspace and rewrite them to registers

  #endif // end if workspace not used for geo_data (hence used for pde coeffs or something else)

#else // if not COAL_READ (i.e. only one or none of PDE_COEFF and GEO_DAT read in a coalsced way)


  #ifdef USE_WORKSPACE_FOR_PDE_COEFF

    offset = nr_elems_this_kercall * EL_GEO_DAT_SIZE + (element_index - thread_id) * NR_PDE_COEFFS_SENT;
	
    // !!! OLD BACKUP	
    //offset = (element_index - thread_id) * NR_PDE_COEFFS_SENT;
    // !!! OLD BACKUP
    
    // TRY TO UNROLL THIS LOOP TO INCREASE MEMORY PARALLELISM !!!
    #ifdef HEAT
    workspace[i*WORK_GROUP_SIZE+thread_id] = zero;
    #else
    for(i=0;i<NR_PDE_COEFFS_SENT;i++) {

      workspace[i*WORK_GROUP_SIZE+thread_id] = el_data_in[offset+i*WORK_GROUP_SIZE+thread_id];
      
      // !!! OLD BACKUP	
      //workspace[i*WORK_GROUP_SIZE+thread_id] = 
      //el_data_in_coeff[offset+i*WORK_GROUP_SIZE+thread_id];
      // !!! OLD BACKUP	

    }

    #ifdef COUNT_OPER
    nr_global_access += NR_PDE_COEFFS_SENT; // we neglect shared memory accesses at this stage
    #endif
    #endif

  #else // if not USE_WORKSPACE_FOR_PDE_COEFF

    offset = nr_elems_this_kercall * EL_GEO_DAT_SIZE + element_index * (NR_PDE_COEFFS_SENT);
    
    // !!! OLD BACKUP	
    //offset = element_index * (NR_PDE_COEFFS_SENT);
    // !!! OLD BACKUP	

    // Initialize PDE_COEFF values
    initialize_pde_coeff(workspace,pde_coeff,&offset);
    
    #ifdef COUNT_OPER
      #ifdef HEAT
    nr_global_access += 0;
      #else
        nr_global_access += NR_PDE_COEFFS_SENT;
      #endif
    #endif

  #endif // end if not USE_WORKSPACE_FOR_PDE_COEFF

  #ifdef USE_WORKSPACE_FOR_GEO_DATA

    offset = (element_index-thread_id)*(EL_GEO_DAT_SIZE);

    barrier(CLK_LOCAL_MEM_FENCE); // I don't know why but without barrier here, one per ten runs gives bad results

    // TRY TO UNROLL THIS LOOP TO INCREASE MEMORY PARALLELISM !!!
    for(i = 0; i < EL_GEO_DAT_SIZE; i++){

      workspace[i*WORK_GROUP_SIZE+thread_id] = el_data_in[offset+i*WORK_GROUP_SIZE+thread_id];
      
      // !!! OLD BACKUP
      //workspace[i*WORK_GROUP_SIZE+thread_id] = el_data_in_geo[offset+i*WORK_GROUP_SIZE+thread_id];
      // !!! OLD BACKUP

    }

  #else // if not USE_WORKSPACE_FOR_GEO_DATA

    // we read geo data to registers
    offset = (element_index)*(EL_GEO_DAT_SIZE);

    // TRY TO UNROLL THIS LOOP TO INCREASE MEMORY PARALLELISM !!!
    for(i = 0; i < EL_GEO_DAT_SIZE; i++){

      geo_dat[i] = el_data_in[offset+i];
      
      // !!! OLD BACKUP
      //geo_dat[i] = el_data_in_geo[offset+i];
      // !!! OLD BACKUP
      
    }

  #endif  // if not USE_WORKSPACE_FOR_GEO_DATA

  #ifdef COUNT_OPER
    nr_global_access += EL_GEO_DAT_SIZE; // we neglect shared memory accesses at this stage
  #endif

#endif // end if not COAL_READ

// ******* THE END OF: READING INPUT DATA *********************
//-------------------------------------------------------------


#if (WORKSPACE_SIZE) > 0
    barrier(CLK_LOCAL_MEM_FENCE); // !!!!!!!!!!!!!!!!!!!!!!
#endif

//-------------------------------------------------------------
//******************** INITIALIZING SM AND LV ******************//

#ifdef USE_WORKSPACE_FOR_STIFF_MAT

    // stiff_mat_workspace holds SM and LV
    for(i = 0; i < num_dofs*(num_dofs+1); i++) {
      workspace[thread_id*(num_dofs*(num_dofs+1)+PADDING)+i] = zero;
    }

#ifdef COUNT_OPER
    nr_access_shared += num_dofs*(num_dofs+1);
#endif

#else // if not  USE_WORKSPACE_FOR_STIFF_MAT


    for(i = 0; i < num_dofs*num_dofs; i++) stiff_mat[i] = zero;

  #ifdef LOAD_VEC_COMP
    for(i = 0; i < num_dofs; i++) load_vec[i] = zero;
  #endif

#endif // end if not  USE_WORKSPACE_FOR_STIFF_MAT

//******************** END OF: INITIALIZING SM AND LV ******************//
//-------------------------------------------------------------



//-------------------------------------------------------------
//************************* INITIALIZE HEAT DATA ************************//
        
    #ifdef USE_WORKSPACE_FOR_PDE_COEFF
      // read sol_dofs_n to local memory
      #ifdef HEAT

	/* EL_DATA_IN_COEFF:
	 * [e+0]                               - h_k
	 * [e+1] ... [e+1+num_shap]            - sol_dofs_k
	 * [e+1+num_shap] ... [e+1+2*num_shap] - sol_dofs_n
	 * [e+1+2*num_shap+i*num_gauss+0]      - vx
	 * [e+1+2*num_shap+i*num_gauss+1]      - vy
	 * [e+1+2*num_shap+i*num_gauss+2]      - vz
	 * [e+1+2*num_shap+i*num_gauss+3]      - thermal conductivity
	 * [e+1+2*num_shap+i*num_gauss+4]      - density times specific heat
	 */

	//Read pde coeffs from global memory
	offset = nr_elems_this_kercall * EL_GEO_DAT_SIZE + (element_index - thread_id) * NR_PDE_COEFFS_SENT;
	workspace[thread_id*(NR_PDE_COEFFS_READ+NR_COEFFS_IN_SM_CALCULATIONS)+NR_COEFFS_IN_SM_CALCULATIONS+0] = el_data_in[offset+0];
	
	for(i=1; i<NR_PDE_COEFFS_READ; i++) {	  
	  workspace[thread_id*(NR_PDE_COEFFS_READ+NR_COEFFS_IN_SM_CALCULATIONS)+NR_COEFFS_IN_SM_CALCULATIONS+i] = el_data_in[offset+num_shap+i];
	}

	#ifdef COUNT_OPER
	  nr_global_access += NR_PDE_COEFFS_READ;
	#endif
	
      #endif
    #else	
      // read sol_dofs_n to private memory
      #ifdef HEAT
        SCALAR h_k; //h_k
	// SCALAR sol_dofs[num_shap+num_shap]; // sol_dofs_k + sol_dofs_n [for future code]
	SCALAR sol_dofs[num_shap]; // sol_dofs_n
	SCALAR coeff_in_gauss[num_gauss*5]; // coeff in gauss points [v0,v1,v2,tc,dtsh]G0 ... [v0,v1,v2,tc,dtsh]GN
      
	offset = nr_elems_this_kercall * EL_GEO_DAT_SIZE + (element_index - thread_id) * NR_PDE_COEFFS_SENT;

	h_k=el_data_in[offset+0];
	
	// [for future code]
	//for(i=0; i<(num_shap+num_shap); i++) {
	//  sol_dofs[i] = el_data_in_coeff[offset+i+1];
	//}
	for(i=0; i<num_shap; i++) {
	  sol_dofs[i] = el_data_in[offset+i+num_shap+1];
	}
	for(i=0; i<(num_gauss*5); i++) {
	  coeff_in_gauss[i] = el_data_in[offset+i+num_shap+num_shap+1];
	}
        #ifdef COUNT_OPER
          nr_global_access += 5*num_gauss+num_shap+1;
	  //nr_global_access += 5*num_gauss+2*num_shap+1; [for future code]
        #endif
      #endif
    #endif

//************************* END OF: INITIALIZE HEAT DATA ************************//
//-------------------------------------------------------------


//-------------------------------------------------------------
//************************* LOOP OVER INTEGRATION POINTS ************************//

    // in a loop over gauss points
    int igauss;
    int idof, jdof;
    for(igauss = 0; igauss < num_gauss; igauss++){


      // integration data read from cached constant or shared  memory
      SCALAR daux = gauss_dat[4*igauss];
      SCALAR faux = gauss_dat[4*igauss+1];
      SCALAR eaux = gauss_dat[4*igauss+2];
      //SCALAR vol = gauss_dat[4*igauss+3]; // vol = weight
      SCALAR vol = weight_gauss; // vol = weight CONSTANT FOR LINEAR PRISM!!!

#ifdef COUNT_OPER
    nr_access_shared += 4;
#endif


//-------------------------------------------------------------
//************************* JACOBIAN TERMS CALCULATIONS *************************//

      // when geometrical shape functions are not necessary 
      // (only derivatives are used for Jacobian calculations)
      SCALAR temp1 = zero;
      SCALAR temp2 = zero;
      SCALAR temp3 = zero;
      SCALAR temp4 = zero;
      SCALAR temp5 = zero;
      SCALAR temp6 = zero;
      SCALAR temp7 = zero;
      SCALAR temp8 = zero;
      SCALAR temp9 = zero;

#ifdef COMPUTE_ALL_SHAPE_FUN_DER
      { // block to indicate the scope of jac_x registers
#endif

      SCALAR jac_0 = zero;
      SCALAR jac_1 = zero;
      SCALAR jac_2 = zero;
      SCALAR jac_3 = zero;
      SCALAR jac_4 = zero;
      SCALAR jac_5 = zero;
      SCALAR jac_6 = zero;
      SCALAR jac_7 = zero;
      SCALAR jac_8 = zero;

      // derivatives of geometrical shape functions
      { // block to indicate the scope of jac_data

        // derivatives of geometrical shape functions are stored in jac_data
	SCALAR jac_data[3*num_geo_dofs];
	jac_data[0] = -(one-eaux)*half;
	jac_data[1] =  (one-eaux)*half;
	jac_data[2] =  zero;
	jac_data[3] = -(one+eaux)*half;
	jac_data[4] =  (one+eaux)*half;
	jac_data[5] =  zero;
	jac_data[6] = -(one-eaux)*half;
	jac_data[7] =  zero;
	jac_data[8] =  (one-eaux)*half;
	jac_data[9] = -(one+eaux)*half;
	jac_data[10] =  zero;
	jac_data[11] =  (one+eaux)*half;
	jac_data[12] = -(one-daux-faux)*half;
	jac_data[13] = -daux*half;
	jac_data[14] = -faux*half;
	jac_data[15] =  (one-daux-faux)*half;
	jac_data[16] =  daux*half;
	jac_data[17] =  faux*half;


#ifdef COUNT_OPER
	nr_oper += 14; // after optimization
#endif


	/* Jacobian matrix J */
#ifdef USE_WORKSPACE_FOR_GEO_DATA
	offset=thread_id*EL_GEO_DAT_SIZE;
#endif

	for(i=0;i<num_geo_dofs;i++){

	  jac_1 = jac_data[i];
	  jac_2 = jac_data[num_geo_dofs+i];
	  jac_3 = jac_data[2*num_geo_dofs+i];

#ifdef USE_WORKSPACE_FOR_GEO_DATA

	  jac_4 = workspace[offset+3*i];  //node coor
	  jac_5 = workspace[offset+3*i+1];
	  jac_6 = workspace[offset+3*i+2];

#ifdef COUNT_OPER
	  nr_access_shared += 3;
#endif

#else // if not USE_WORKSPACE_FOR_GEO_DATA (geo data in registers)

	  jac_4 = geo_dat[3*i];  //node coor
	  jac_5 = geo_dat[3*i+1];
	  jac_6 = geo_dat[3*i+2];

#endif // end if not USE_GEO_DAT_WORKSPACE

	  temp1 += jac_4 * jac_1;
	  temp2 += jac_4 * jac_2;
	  temp3 += jac_4 * jac_3;
	  temp4 += jac_5 * jac_1;
	  temp5 += jac_5 * jac_2;
	  temp6 += jac_5 * jac_3;
	  temp7 += jac_6 * jac_1;
	  temp8 += jac_6 * jac_2;
	  temp9 += jac_6 * jac_3;

	}

      } // the end of scope for jac_data

#ifdef COUNT_OPER
      nr_oper += 18*num_geo_dofs; // after optimization
#endif


      jac_0 = (temp5*temp9 - temp8*temp6);
      jac_1 = (temp8*temp3 - temp2*temp9);
      jac_2 = (temp2*temp6 - temp3*temp5);

      daux = temp1*jac_0 + temp4*jac_1 + temp7*jac_2;

      /* Jacobian calculations - |J| and inverse of the Jacobian matrix*/
      vol *= daux; // vol = weight * det J

      faux = one/daux;

      jac_0 *= faux;
      jac_1 *= faux;
      jac_2 *= faux;

      jac_3 = (temp6*temp7 - temp4*temp9)*faux;
      jac_4 = (temp1*temp9 - temp7*temp3)*faux;
      jac_5 = (temp3*temp4 - temp1*temp6)*faux;

      jac_6 = (temp4*temp8 - temp5*temp7)*faux;
      jac_7 = (temp2*temp7 - temp1*temp8)*faux;
      jac_8 = (temp1*temp5 - temp2*temp4)*faux;

#ifdef COUNT_OPER
 nr_oper += 43; // after optimization, includes 1 inverse and 6 sign changes
 // total: 14+18*num_geo_dofs+43 = 165 (for prisms)
#endif

//************* THE END OF: JACOBIAN TERMS CALCULATIONS *************************//
//-------------------------------------------------------------


//-------------------------------------------------------------
//***** SEPARATE COMPUTING OF ALL GLOBAL DERIVATIVES OF ALL SHAPE FUNCTIONS *****//

#ifdef COMPUTE_ALL_SHAPE_FUN_DER

 //************ loop for computing ALL shape function values at integration point **********//
      for(idof = 0; idof < num_shap; idof++){

	// read proper values of shape functions and their derivatives
	temp1 = shpfun_ref[igauss*4*num_shap+4*idof+1];
	temp2 = shpfun_ref[igauss*4*num_shap+4*idof+2];
	temp3 = shpfun_ref[igauss*4*num_shap+4*idof+3];

  #ifdef COUNT_OPER
	nr_access_shared += 3; // 3 reads from constant cache
  #endif

	// compute derivatives wrt global coordinates
	// 15 operations

  #ifdef USE_WORKSPACE_FOR_SHAPE_FUN

	workspace[thread_id*(3*num_shap+PADDING)+3*idof]   = temp1*jac_0+temp2*jac_3+temp3*jac_6;
	workspace[thread_id*(3*num_shap+PADDING)+3*idof+1] = temp1*jac_1+temp2*jac_4+temp3*jac_7;
	workspace[thread_id*(3*num_shap+PADDING)+3*idof+2] = temp1*jac_2+temp2*jac_5+temp3*jac_8;

#ifdef COUNT_OPER
	nr_access_shared += 3; //  3 writes to shared memory
#endif

  #else  // if not USE_WORKSPACE_FOR_SHAPE_FUN

	tab_fun_u_derx[idof] = temp1*jac_0+temp2*jac_3+temp3*jac_6;
	tab_fun_u_dery[idof] = temp1*jac_1+temp2*jac_4+temp3*jac_7;
	tab_fun_u_derz[idof] = temp1*jac_2+temp2*jac_5+temp3*jac_8;

  #endif // end if not USE_WORKSPACE_FOR_SHAPE_FUN

      } // end loop over shape functions for which global derivatives were computed

#ifdef COUNT_OPER
 nr_oper += 15*num_shap; 
#endif

      } // the end of block to indicate the scope of jac_x registers

#endif // end if COMPUTE_ALL_SHAPE_FUN_DER

//*** THE END OF: SEPARATE COMPUTING OF ALL GLOBAL DERIVATIVES OF ALL SHAPE FUNCTIONS ***//
//-------------------------------------------------------------


//-------------------------------------------------------------
//***** SUBSTITUTING ACTUAL COEFFICIENTS FOR SM AND LV CALCULATIONS *****//

#ifdef USE_WORKSPACE_FOR_PDE_COEFF
      
#ifdef HEAT

	offset = thread_id*(NR_COEFFS_IN_SM_CALCULATIONS+NR_PDE_COEFFS_READ);

	// for non-constant, non-linear coefficients a place for call to problem dependent
	// function calculating actual PDE coefficients based on data in coeff 
	// workspace or registers and  storing data back in workspace or in registers
      
	/* EL_DATA_IN_COEFF:
	 * [e+0]                               - h_k
	 * [e+1] ... [e+1+num_shap]            - sol_dofs_k
	 * [e+1+num_shap] ... [e+1+2*num_shap] - sol_dofs_n
	 * [e+1+2*num_shap+i*num_gauss+0]      - vx
	 * [e+1+2*num_shap+i*num_gauss+1]      - vy
	 * [e+1+2*num_shap+i*num_gauss+2]      - vz
	 * [e+1+2*num_shap+i*num_gauss+3]      - thermal conductivity
	 * [e+1+2*num_shap+i*num_gauss+4]      - density times specific heat
	 *
	 *
	 * COEFF description:
	 * pde_coeff -> 0 - 2   = axx,axy,axz
	 * pde_coeff -> 3 - 5   = xyx,ayy,ayz
	 * pde_coeff -> 6 - 8   = azx,azy,azz
	 * pde_coeff -> 9 - 11  = bx,by,bz
	 * pde_coeff -> 12 - 14 = tx,ty,tz
	 * pde_coeff -> 15      = cval/mval
	 * pde_coeff -> 16 - 18 = qx,qy,qz
	 * pde_coeff -> 19      = sval/lval
	 *
	 */
        
	// copy data to temporary registers       
	// temp1 = (double) el_data_in[0]; //Implicitness coeff
	// temp2 = (double) el_data_in[1]; //Delta T current time step
	// temp3 -> m_k
	// temp4 -> norm_u
	// temp5 -> peclet_local
	// temp6 -> tau_therm      
       
	// calculations auxiliary variables
	//temp1 *= 0.001; //Implicitness parameter [this is being sent as integer]
	//temp2 *= 0.000001; //Delta T [this is being sent as integer]
	//temp3 = 1.0/3.0; //m_k
       
	//norm_u
	temp5 = workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+0]; temp6 = workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+1]; temp7 = workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+2];
	temp4 = sqrt(temp5*temp5+temp6*temp6+temp7*temp7); // sqrt(Vx^2+Vy^2+Vz^2)
	
        #ifdef COUNT_OPER
	  nr_oper += 7;
	  nr_access_shared += 3;
	#endif
	
	
	//peclet_local
	if(workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+3] > 1.e-6) {
	  temp5 = (temp3*temp4*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+0])/(2.0*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+3]);
	  #ifdef COUNT_OPER
	    nr_oper += 7;
	    nr_access_shared += 2;
	  #endif
	} else {
	  temp5 = 1.0e6;
	  #ifdef COUNT_OPER
	    nr_oper += 1;
	  #endif
	}
       
	//tau_therm
	if(temp5 < 1.0) {
	  temp6 = (workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+0]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+0]*temp3)/(4.0*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+3]);
	  #ifdef COUNT_OPER
	    nr_oper += 7;
	    nr_access_shared += 3;
	  #endif
	} else {
	  temp6 = (4.0*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+0])/(2.0*temp4);
	  #ifdef COUNT_OPER
	    nr_oper += 6;
	    nr_access_shared += 1;
	  #endif
	}

	/*! ------------------ CALCULATE VALUES FROM PREVIEW STEP -------------------! */
	SCALAR un_val; SCALAR un_xyz[3];
	un_val = zero; un_xyz[0] = zero; un_xyz[1] = zero; un_xyz[2] = zero;
	#ifdef COUNT_OPER
	  nr_oper += 4;
	#endif

	//Loop over shape functions
	for(idof = 0; idof < num_shap; idof++){
	  
          #ifdef COMPUTE_ALL_SHAPE_FUN_DER	  
	  
	    // read proper values of shape functions and their derivatives
	    SCALAR shp_fun_u = shpfun_ref[igauss*4*num_shap+4*idof];
	    SCALAR fun_u_derx = tab_fun_u_derx[idof];
	    SCALAR fun_u_dery = tab_fun_u_dery[idof];
	    SCALAR fun_u_derz = tab_fun_u_derz[idof];

            #ifdef COUNT_OPER
               nr_access_shared += 1;
            #endif	  
	  
          #else // if not COMPUTE_ALL_SHAPE_FUN_DER

            // read proper values of shape functions and their derivatives
            SCALAR shp_fun_u = shpfun_ref[igauss*4*num_shap+4*idof];
	    SCALAR fun_u_derx  = shpfun_ref[igauss*4*num_shap+4*idof+1];
	    SCALAR fun_u_dery = shpfun_ref[igauss*4*num_shap+4*idof+2];
	    SCALAR fun_u_derz = shpfun_ref[igauss*4*num_shap+4*idof+3];
	    
            #ifdef COUNT_OPER
              nr_access_shared += 4; // constant cache accesses
	      nr_oper += 15; // after optimization
	      // total: 13+5+18*num_geo_dofs+15+36+15*num_shap = 177+90 = 267 (for prisms)
            #endif
          #endif
	      
	  // Computed solution from previous time step
	  un_val += workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+num_shap+1+idof] * shp_fun_u; //sol_dofs[num_dofs + idof] * shp_fun_u; // if sol_dofs_k used
	  // Gradient of computed solution Un_val X
	  un_xyz[0] += workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+num_shap+1+idof] * fun_u_derz; //sol_dofs[num_dofs + idof] * fun_u_derx; // if sol_dofs_k used
	  // Gradient of computed solution Un_val Y
	  un_xyz[1] += workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+num_shap+1+idof] * fun_u_derz; //sol_dofs[num_dofs + idof] * fun_u_dery; // if sol_dofs_k used
	  // Gradient of computed solution Un_val Z
	  un_xyz[2] += workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+num_shap+1+idof] * fun_u_derz; //sol_dofs[num_dofs + idof] * fun_u_derz; // if sol_dofs_k used

	  #ifdef COUNT_OPER
	    nr_oper += 8;
	  #endif
	}

        // Copy Implicitness coeff to register
        temp1 = ((SCALAR) el_data_in[0]) * 0.001; //Implicitness coeff
	// Copy Delta T current time step to register
	temp2 = ((SCALAR) el_data_in[1]) * 0.000001; //Delta T current time step

	// Calculate m_k val
	temp3 = 1.0/3.0;

	#ifdef COUNT_OPER
	  nr_oper += 6;
	#endif
       
        /*! ------------------ CALCULATE STABILIZATION COEFFS -------------------! */
       
	// AXX, AXY, AXZ
	workspace[offset+0] = temp6*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+0]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+0]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4] + temp1*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+3]; //AXX tau_therm*Vx^2*density_time_specific_heat + implicitness_coeff*tc /* 6 nr_oper */
	workspace[offset+1] = temp6*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+0]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+1]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4]; //AXY tau_therm*Vx*Vy*density_time_specific_heat /* 4 nr_oper */
	workspace[offset+2] = temp6*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+0]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+2]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4]; //AXZ tau_therm*Vx*Vz*density_time_specific_heat /* 4 nr_oper */
	// AYX, AYY, AYZ
	workspace[offset+3] = temp6*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+1]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+0]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4]; //AYX tau_therm*Vy*Vx*density_time_specific_heat /* 4 nr_oper */
	workspace[offset+4] = temp6*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+1]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+1]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4] + temp1*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+3]; //AYY tau_therm*Vy*Vy*density_time_specific_heat + implicitness_coeff*tc /* 6 nr_oper */
	workspace[offset+5] = temp6*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+1]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+2]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4]; //AYZ tau_therm*Vy*Vz*density_time_specific_heat /* 4 nr_oper */
	//AZX, AZY, AZZ
	workspace[offset+6] = temp6*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+2]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+0]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4]; //AZX tau_therm*Vz*Vx*density_time_specific_heat /* 4 nr_oper */
	workspace[offset+7] = temp6*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+2]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+1]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4]; //AZY tau_therm*Vz*Vy*density_time_specific_heat /* 4 nr_oper */
	workspace[offset+8] = temp6*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+2]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+2]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4] + temp1*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+3]; //AZZ tau_therm*Vz*Vz*density_time_specific_heat + implicitness_coeff*tc /* 6 nr_oper */

	// Calculate BX, BY, BZ
	workspace[offset+9] = temp1*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+0]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4]; //BX implicitness_coeff*Vx*density_time_specific_heat /* 3 nr_oper */
	workspace[offset+10] = temp1*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+1]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4]; //BY implicitness_coeff*Vy*density_time_specific_heat /* 3 nr_oper */
	workspace[offset+11] = temp1*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+2]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4]; //BZ implicitness_coeff*Vz*density_time_specific_heat /* 3 nr_oper */
       
	//TX, TY, TZ
	workspace[offset+12] = temp6*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+0]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4]/temp2; //TX tau_therm*Vx*density_time_specific_heat/delta_t /* 4 nr_oper */
	workspace[offset+13] = temp6*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+1]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4]/temp2; //TY tau_therm*Vy*density_time_specific_heat/delta_t /* 4 nr_oper */
	workspace[offset+14] = temp6*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+2]*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4]/temp2; //TZ tau_therm*Vz*density_time_specific_heat/delta_t /* 4 nr_oper */
       
	// Calculate cval/mval
	workspace[offset+15] = workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4]/temp2; //CVAL,MVAL density_time_specific_heat/delta_t /* 2 nr_oper */
       
	// QX, QY, QZ
	if(temp1<1.0) { //Implicitness_coeff < 1.0
	  workspace[offset+16] = (temp1-1.0)*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+3]*un_xyz[0]; //Qx[0] += (Implicitness_coeff-1.0)*Thermal_diffusivity*Un_x[0] /* 3 nr_oper */
	  workspace[offset+17] = (temp1-1.0)*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+3]*un_xyz[1]; //Qy[0] += (Implicitness_coeff-1.0)*Thermal_diffusivity*Un_y[0] /* 3 nr_oper */
	  workspace[offset+18] = (temp1-1.0)*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+3]*un_xyz[2]; //Qz[0] += (Implicitness_coeff-1.0)*Thermal_diffusivity*Un_z[0] /* 3 nr_oper */
	  #ifdef COUNT_OPER
	    nr_oper += 9;
	    nr_access_shared += 3;
	  #endif
	}
	workspace[offset+16] = un_val*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4]/temp2*temp6*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+0]; //Qx[0] += Un_val[0] * Density_times_specific_heat / Delta_t * tau_therm * Velocity[0] /* 5 nr_oper */
	workspace[offset+17] = un_val*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4]/temp2*temp6*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+1]; //Qy[0] += Un_val[0] * Density_times_specific_heat / Delta_t * tau_therm * Velocity[1] /* 5 nr_oper */
	workspace[offset+18] = un_val*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4]/temp2*temp6*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+2]; //Qy[0] += Un_val[0] * Density_times_specific_heat / Delta_t * tau_therm * Velocity[2] /* 5 nr_oper */
	
	// SVAL
	if(temp1<1.0) { //Implicitness_coeff < 1.0
	  // Sval[0] = (Implicitness_coeff-1.0) *			\
	  // (Velocity[0]*Un_x[0] + Velocity[1]*Un_y[0] + Velocity[2]*Un_z[0]) *
	  // Density_times_specific_heat
	  workspace[offset+19] = (temp1-1.0)*(workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+0]*un_xyz[0]+workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+1]*un_xyz[1]+workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+2]*un_xyz[2])*workspace[offset+NR_COEFFS_IN_SM_CALCULATIONS+1+num_dofs+igauss*5+4]; /* 9 nr_oper */
	  #ifdef COUNT_OPER
	    nr_oper += 9;
	    nr_access_shared += 1;
	  #endif
	} else {
	  workspace[offset+19] = zero;
	}

	#ifdef COUNT_OPER
	  nr_oper += 78;
	  nr_access_shared += 20;
	#endif

	/* Print HEAT COEFF
	   if(get_global_id(0) == 0 && igauss == 0) {
            printf("ielem = %d\n",ielem); 
            printf("Un_val[0] = %lf ; Un_x[0] = %lf ; Un_y[0] = %lf ; Un_z[0] = %lf\n",un_val,un_xyz[0],un_xyz[1],un_xyz[2]);
            printf("AXX = %lf ; AXY = %lf ; AXZ = %lf\n",pde_coeff[0],pde_coeff[1],pde_coeff[2]);
            printf("AYX = %lf ; AYY = %lf ; AYZ = %lf\n",pde_coeff[3],pde_coeff[4],pde_coeff[5]);
            printf("AZX = %lf ; AZY = %lf ; AZZ = %lf\n",pde_coeff[6],pde_coeff[7],pde_coeff[8]);
            printf("BX = %lf ; BY = %lf ; BZ = %lf\n",pde_coeff[9],pde_coeff[10],pde_coeff[11]);
            printf("TX = %lf ; TY = %lf ; TZ = %lf\n",pde_coeff[12],pde_coeff[13],pde_coeff[14]);
            printf("QX = %lf ; QY = %lf ; QZ = %lf\n",pde_coeff[16],pde_coeff[17],pde_coeff[18]);
            printf("CVAL/MVAL = %lf\n",pde_coeff[15]);
            printf("SVAL = %lf\n",pde_coeff[19]);
	   }
        /* Print HEAT COEFF */
      
#endif


  #ifdef LAPLACE

      // offset for reading data
      offset=thread_id*(NR_PDE_COEFFS_SENT);

      pde_coeff[12] = workspace[offset+igauss]; //coeff03

    #ifdef COUNT_OPER
    nr_access_shared += 1;
    #endif

  #endif // end if LAPLACE

  // offset for computations
  #ifdef HEAT
    offset=thread_id*(NR_COEFFS_IN_SM_CALCULATIONS+NR_PDE_COEFFS_READ);
  #else
    offset=thread_id*(NR_COEFFS_IN_SM_CALCULATIONS);
  #endif

#else // if not USE_WORKSPACE_FOR_PDE_COEFF


  #ifdef HEAT

        offset = thread_id*NR_PDE_COEFFS_SENT;

	// for non-constant, non-linear coefficients a place for call to problem dependent
	// function calculating actual PDE coefficients based on data in coeff 
	// workspace or registers and  storing data back in workspace or in registers
      
	/* EL_DATA_IN_COEFF:
	 * [e+0]                               - h_k
	 * [e+1] ... [e+1+num_shap]            - sol_dofs_k
	 * [e+1+num_shap] ... [e+1+2*num_shap] - sol_dofs_n
	 * [e+1+2*num_shap+i*num_gauss+0]      - vx
	 * [e+1+2*num_shap+i*num_gauss+1]      - vy
	 * [e+1+2*num_shap+i*num_gauss+2]      - vz
	 * [e+1+2*num_shap+i*num_gauss+3]      - thermal conductivity
	 * [e+1+2*num_shap+i*num_gauss+4]      - density times specific heat
	 *
	 *
	 * COEFF description:
	 * pde_coeff -> 0 - 2   = axx,axy,axz
	 * pde_coeff -> 3 - 5   = xyx,ayy,ayz
	 * pde_coeff -> 6 - 8   = azx,azy,azz
	 * pde_coeff -> 9 - 11  = bx,by,bz
	 * pde_coeff -> 12 - 14 = tx,ty,tz
	 * pde_coeff -> 15      = cval/mval
	 * pde_coeff -> 16 - 18 = qx,qy,qz
	 * pde_coeff -> 19      = sval/lval
	 *
	 */
        
	// copy data to temporary registers       
	// temp1 = (double) el_data_in[0]; //Implicitness coeff
	// temp2 = (double) el_data_in[1]; //Delta T current time step
	// temp3 -> m_k
	// temp4 -> norm_u
	// temp5 -> peclet_local
	// temp6 -> tau_therm      
       
	//coeff_in_gauss[igauss*5+0] ; //Velocity X
	//coeff_in_gauss[igauss*5+1] ; //Velocity Y
	//coeff_in_gauss[igauss*5+2] ; //Velocity Z
       
	//coeff_in_gauss[igauss*5+3] ; //Thermal condactivity
	//coeff_in_gauss[igauss*5+4] ; //Density time specific heat
       
	// calculations auxiliary variables
	//temp1 *= 0.001; //Implicitness parameter [this is being sent as integer]
	//temp2 *= 0.000001; //Delta T [this is being sent as integer]
	//temp3 = 1.0/3.0; //m_k
       
	//norm_u
	temp4 = sqrt(coeff_in_gauss[igauss*5+0]*coeff_in_gauss[igauss*5+0] + coeff_in_gauss[igauss*5+1]*coeff_in_gauss[igauss*5+1] + coeff_in_gauss[igauss*5+2]*coeff_in_gauss[igauss*5+2]); // sqrt(Vx^2+Vy^2+Vz^2)
	
        #ifdef COUNT_OPER
	  nr_oper += 7;
	#endif
	
	
	//peclet_local
	if(coeff_in_gauss[igauss*5+3] > 1.e-6) {
	  temp5 = (temp3*temp4*h_k)/(2.0*coeff_in_gauss[igauss*5+3]);
	  #ifdef COUNT_OPER
	    nr_oper += 7;
	  #endif
	} else {
	  temp5 = 1.0e6;
	  #ifdef COUNT_OPER
	    nr_oper += 1;
	  #endif
	}
       
	//tau_therm
	if(temp5 < 1.0) {
	  temp6 = (h_k*h_k*temp3)/(4.0*coeff_in_gauss[igauss*5+3]);
	  #ifdef COUNT_OPER
	    nr_oper += 7;
	  #endif
	} else {
	  temp6 = (4.0*h_k)/(2.0*temp4);
	  #ifdef COUNT_OPER
	    nr_oper += 6;
	  #endif
	}

	/*! ------------------ CALCULATE VALUES FROM PREVIEW STEP -------------------! */
	SCALAR un_val; SCALAR un_xyz[3];
	un_val = zero; un_xyz[0] = zero; un_xyz[1] = zero; un_xyz[2] = zero;
	#ifdef COUNT_OPER
	  nr_oper += 4;
	#endif

	//Loop over shape functions
	for(idof = 0; idof < num_shap; idof++){
	  
          #ifdef COMPUTE_ALL_SHAPE_FUN_DER	  
             #ifdef USE_WORKSPACE_FOR_SHAPE_FUN
	     
	       SCALAR shp_fun_u = shpfun_ref[igauss*4*num_shap+4*idof];
	       SCALAR fun_u_derx = workspace[thread_id*(3*num_shap+PADDING)+3*idof];
	       SCALAR fun_u_dery = workspace[thread_id*(3*num_shap+PADDING)+3*idof+1];
	       SCALAR fun_u_derz = workspace[thread_id*(3*num_shap+PADDING)+3*idof+2];
	  
               #ifdef COUNT_OPER
                 nr_access_shared += 4; // including 1 for constant cache
               #endif
             #else // if not USE_WORKSPACE_FOR_SHAPE_FUN
	  
	       // read proper values of shape functions and their derivatives
	       SCALAR shp_fun_u = shpfun_ref[igauss*4*num_shap+4*idof];
	       SCALAR fun_u_derx = tab_fun_u_derx[idof];
	       SCALAR fun_u_dery = tab_fun_u_dery[idof];
	       SCALAR fun_u_derz = tab_fun_u_derz[idof];

               #ifdef COUNT_OPER
                 nr_access_shared += 1;
               #endif	  
            #endif // end if not USE_SHAPE_FUN_WORKSPACE	  
          #else // if not COMPUTE_ALL_SHAPE_FUN_DER

            // read proper values of shape functions and their derivatives
            SCALAR shp_fun_u = shpfun_ref[igauss*4*num_shap+4*idof];
	    SCALAR fun_u_derx  = shpfun_ref[igauss*4*num_shap+4*idof+1];
	    SCALAR fun_u_dery = shpfun_ref[igauss*4*num_shap+4*idof+2];
	    SCALAR fun_u_derz = shpfun_ref[igauss*4*num_shap+4*idof+3];
	    
            #ifdef COUNT_OPER
              nr_access_shared += 4; // constant cache accesses
	      nr_oper += 15; // after optimization
	      // total: 13+5+18*num_geo_dofs+15+36+15*num_shap = 177+90 = 267 (for prisms)
            #endif
          #endif
	      
	  // Computed solution from previous time step
	  un_val += sol_dofs[idof] * shp_fun_u; //sol_dofs[num_dofs + idof] * shp_fun_u; // if sol_dofs_k used
	  // Gradient of computed solution Un_val X
	  un_xyz[0] += sol_dofs[idof] * fun_u_derz; //sol_dofs[num_dofs + idof] * fun_u_derx; // if sol_dofs_k used
	  // Gradient of computed solution Un_val Y
	  un_xyz[1] += sol_dofs[idof] * fun_u_derz; //sol_dofs[num_dofs + idof] * fun_u_dery; // if sol_dofs_k used
	  // Gradient of computed solution Un_val Z
	  un_xyz[2] += sol_dofs[idof] * fun_u_derz; //sol_dofs[num_dofs + idof] * fun_u_derz; // if sol_dofs_k used

	  #ifdef COUNT_OPER
	    nr_oper += 8;
	  #endif
	}

	// Copy Implicitness coeff to register
        temp1 = ((SCALAR) el_data_in[0]) * 0.001; //Implicitness coeff
	// Copy Delta T current time step to register
	temp2 = ((SCALAR) el_data_in[1]) * 0.000001; //Delta T current time step
	
	// Calculate m_k val
	temp3 = 1.0/3.0;

	#ifdef COUNT_OPER
	  nr_oper += 6;
	#endif
       
        /*! ------------------ CALCULATE STABILIZATION COEFFS -------------------! */
       
	// AXX, AXY, AXZ
	pde_coeff[0] = temp6*coeff_in_gauss[igauss*5+0]*coeff_in_gauss[igauss*5+0]*coeff_in_gauss[igauss*5+4] + temp1*coeff_in_gauss[igauss*5+3]; //AXX tau_therm*Vx^2*density_time_specific_heat + implicitness_coeff*tc /* 6 nr_oper */
	pde_coeff[1] = temp6*coeff_in_gauss[igauss*5+0]*coeff_in_gauss[igauss*5+1]*coeff_in_gauss[igauss*5+4]; //AXY tau_therm*Vx*Vy*density_time_specific_heat /* 4 nr_oper */
	pde_coeff[2] = temp6*coeff_in_gauss[igauss*5+0]*coeff_in_gauss[igauss*5+2]*coeff_in_gauss[igauss*5+4]; //AXZ tau_therm*Vx*Vz*density_time_specific_heat /* 4 nr_oper */
	// AYX, AYY, AYZ
	pde_coeff[3] = temp6*coeff_in_gauss[igauss*5+1]*coeff_in_gauss[igauss*5+0]*coeff_in_gauss[igauss*5+4]; //AYX tau_therm*Vy*Vx*density_time_specific_heat /* 4 nr_oper */
	pde_coeff[4] = temp6*coeff_in_gauss[igauss*5+1]*coeff_in_gauss[igauss*5+1]*coeff_in_gauss[igauss*5+4] + temp1*coeff_in_gauss[igauss*5+3]; //AYY tau_therm*Vy*Vy*density_time_specific_heat + implicitness_coeff*tc /* 6 nr_oper */
	pde_coeff[5] = temp6*coeff_in_gauss[igauss*5+1]*coeff_in_gauss[igauss*5+2]*coeff_in_gauss[igauss*5+4]; //AYZ tau_therm*Vy*Vz*density_time_specific_heat /* 4 nr_oper */
	//AZX, AZY, AZZ
	pde_coeff[6] = temp6*coeff_in_gauss[igauss*5+2]*coeff_in_gauss[igauss*5+0]*coeff_in_gauss[igauss*5+4]; //AZX tau_therm*Vz*Vx*density_time_specific_heat /* 4 nr_oper */
	pde_coeff[7] = temp6*coeff_in_gauss[igauss*5+2]*coeff_in_gauss[igauss*5+1]*coeff_in_gauss[igauss*5+4]; //AZY tau_therm*Vz*Vy*density_time_specific_heat /* 4 nr_oper */
	pde_coeff[8] = temp6*coeff_in_gauss[igauss*5+2]*coeff_in_gauss[igauss*5+2]*coeff_in_gauss[igauss*5+4] + temp1*coeff_in_gauss[igauss*5+3]; //AZZ tau_therm*Vz*Vz*density_time_specific_heat + implicitness_coeff*tc /* 6 nr_oper */

	// Calculate BX, BY, BZ
	pde_coeff[9] = temp1*coeff_in_gauss[igauss*5+0]*coeff_in_gauss[igauss*5+4]; //BX implicitness_coeff*Vx*density_time_specific_heat /* 3 nr_oper */
	pde_coeff[10] = temp1*coeff_in_gauss[igauss*5+1]*coeff_in_gauss[igauss*5+4]; //BY implicitness_coeff*Vy*density_time_specific_heat /* 3 nr_oper */
	pde_coeff[11] = temp1*coeff_in_gauss[igauss*5+2]*coeff_in_gauss[igauss*5+4]; //BZ implicitness_coeff*Vz*density_time_specific_heat /* 3 nr_oper */
       
	//TX, TY, TZ
	pde_coeff[12] = temp6*coeff_in_gauss[igauss*5+0]*coeff_in_gauss[igauss*5+4]/temp2; //TX tau_therm*Vx*density_time_specific_heat/delta_t /* 4 nr_oper */
	pde_coeff[13] = temp6*coeff_in_gauss[igauss*5+1]*coeff_in_gauss[igauss*5+4]/temp2; //TY tau_therm*Vy*density_time_specific_heat/delta_t /* 4 nr_oper */
	pde_coeff[14] = temp6*coeff_in_gauss[igauss*5+2]*coeff_in_gauss[igauss*5+4]/temp2; //TZ tau_therm*Vz*density_time_specific_heat/delta_t /* 4 nr_oper */
       
	// Calculate cval/mval
	pde_coeff[15] = coeff_in_gauss[igauss*5+4]/temp2; //CVAL,MVAL density_time_specific_heat/delta_t /* 2 nr_oper */
       
	// QX, QY, QZ
	if(temp1<1.0) { //Implicitness_coeff < 1.0
	  pde_coeff[16] = (temp1-1.0)*coeff_in_gauss[igauss*5+3]*un_xyz[0]; //Qx[0] += (Implicitness_coeff-1.0)*Thermal_diffusivity*Un_x[0] /* 3 nr_oper */
	  pde_coeff[17] = (temp1-1.0)*coeff_in_gauss[igauss*5+3]*un_xyz[1]; //Qy[0] += (Implicitness_coeff-1.0)*Thermal_diffusivity*Un_y[0] /* 3 nr_oper */
	  pde_coeff[18] = (temp1-1.0)*coeff_in_gauss[igauss*5+3]*un_xyz[2]; //Qz[0] += (Implicitness_coeff-1.0)*Thermal_diffusivity*Un_z[0] /* 3 nr_oper */
	  #ifdef COUNT_OPER
	    nr_oper += 9;
	  #endif
	}
	pde_coeff[16] = un_val*coeff_in_gauss[igauss*5+4]/temp2*temp6*coeff_in_gauss[igauss*5+0]; //Qx[0] += Un_val[0] * Density_times_specific_heat / Delta_t * tau_therm * Velocity[0] /* 5 nr_oper */
	pde_coeff[17] = un_val*coeff_in_gauss[igauss*5+4]/temp2*temp6*coeff_in_gauss[igauss*5+1]; //Qy[0] += Un_val[0] * Density_times_specific_heat / Delta_t * tau_therm * Velocity[1] /* 5 nr_oper */
	pde_coeff[18] = un_val*coeff_in_gauss[igauss*5+4]/temp2*temp6*coeff_in_gauss[igauss*5+2]; //Qy[0] += Un_val[0] * Density_times_specific_heat / Delta_t * tau_therm * Velocity[2] /* 5 nr_oper */
	
	// SVAL
	if(temp1<1.0) { //Implicitness_coeff < 1.0
	  // Sval[0] = (Implicitness_coeff-1.0) *			\
	  // (Velocity[0]*Un_x[0] + Velocity[1]*Un_y[0] + Velocity[2]*Un_z[0]) *
	  // Density_times_specific_heat
	  pde_coeff[19] = (temp1-1.0)*(coeff_in_gauss[igauss*5+0]*un_xyz[0]+coeff_in_gauss[igauss*5+1]*un_xyz[1]+coeff_in_gauss[igauss*5+2]*un_xyz[2])*coeff_in_gauss[igauss*5+4]; /* 9 nr_oper */
	  #ifdef COUNT_OPER
	    nr_oper += 9;
	  #endif
	}

	#ifdef COUNT_OPER
	  nr_oper += 78;
	#endif

	/* Print HEAT COEFF
	   if(get_global_id(0) == 0 && igauss == 0) {
            printf("ielem = %d\n",ielem); 
            printf("Un_val[0] = %lf ; Un_x[0] = %lf ; Un_y[0] = %lf ; Un_z[0] = %lf\n",un_val,un_xyz[0],un_xyz[1],un_xyz[2]);
            printf("AXX = %lf ; AXY = %lf ; AXZ = %lf\n",pde_coeff[0],pde_coeff[1],pde_coeff[2]);
            printf("AYX = %lf ; AYY = %lf ; AYZ = %lf\n",pde_coeff[3],pde_coeff[4],pde_coeff[5]);
            printf("AZX = %lf ; AZY = %lf ; AZZ = %lf\n",pde_coeff[6],pde_coeff[7],pde_coeff[8]);
            printf("BX = %lf ; BY = %lf ; BZ = %lf\n",pde_coeff[9],pde_coeff[10],pde_coeff[11]);
            printf("TX = %lf ; TY = %lf ; TZ = %lf\n",pde_coeff[12],pde_coeff[13],pde_coeff[14]);
            printf("QX = %lf ; QY = %lf ; QZ = %lf\n",pde_coeff[16],pde_coeff[17],pde_coeff[18]);
            printf("CVAL/MVAL = %lf\n",pde_coeff[15]);
            printf("SVAL = %lf\n",pde_coeff[19]);
	   }
        /* Print HEAT COEFF */    

  #endif

  #ifdef LAPLACE

      pde_coeff[12] = pde_coeff[igauss];

      // !!!OLD BACKUP
      //SCALAR coeff03 = pde_coeff[igauss]; // DOES NOT WORK FOR AMD !!!!!!!!!
      //SCALAR coeff03 = 0.0;
      //switch(igauss){
      //case 0:
      //coeff03 = pde_coeff[0];
      //break;
      //case 1:
      //coeff03 = pde_coeff[1];
      //break;
      //case 2:
      //coeff03 = pde_coeff[2];
      //break;
      //case 3:
      //coeff03 = pde_coeff[3];
      //break;
      //case 4:
      //coeff03 = pde_coeff[4];
      //break;
      //case 5:
      //coeff03 = pde_coeff[5];
      //break;
      //}
      // !!!OLD BACKUP

  #endif // end if LAPLACE

#endif // end if not USE_WORKSPACE_FOR_PDE_COEFF

//*** THE END OF: SUBSTITUTING ACTUAL COEFFICIENTS FOR SM AND LV CALCULATIONS ***//
//-------------------------------------------------------------


//-------------------------------------------------------------
//********************* first loop over shape functions ***********************//
      for(idof = 0; idof < num_shap; idof++){
	
	{ // beginning of using registers for u  (shp_fun_u, fun_u_der.)
	  

//-------------------------------------------------------------
//****** SUBSTITUTING OR COMPUTING GLOBAL DERIVATIVES OF IDOF SHAPE FUNCTION ******//

#ifdef COMPUTE_ALL_SHAPE_FUN_DER
	  
  #ifdef USE_WORKSPACE_FOR_SHAPE_FUN
	  
	  SCALAR shp_fun_u = shpfun_ref[igauss*4*num_shap+4*idof];
	  SCALAR fun_u_derx = workspace[thread_id*(3*num_shap+PADDING)+3*idof];
	  SCALAR fun_u_dery = workspace[thread_id*(3*num_shap+PADDING)+3*idof+1];
	  SCALAR fun_u_derz = workspace[thread_id*(3*num_shap+PADDING)+3*idof+2];

    #ifdef COUNT_OPER
          nr_access_shared += 4; // including 1 for constant cache
    #endif

  #else // if not USE_WORKSPACE_FOR_SHAPE_FUN

	  // read proper values of shape functions and their derivatives
          SCALAR shp_fun_u = shpfun_ref[igauss*4*num_shap+4*idof];
          SCALAR fun_u_derx = tab_fun_u_derx[idof];
          SCALAR fun_u_dery = tab_fun_u_dery[idof];
          SCALAR fun_u_derz = tab_fun_u_derz[idof];

    #ifdef COUNT_OPER
          nr_access_shared += 1;
    #endif
	  
  #endif // end if not USE_SHAPE_FUN_WORKSPACE

#else // if not COMPUTE_ALL_SHAPE_FUN_DER

	  // read proper values of shape functions and their derivatives
	  SCALAR shp_fun_u = shpfun_ref[igauss*4*num_shap+4*idof];
	  temp1 = shpfun_ref[igauss*4*num_shap+4*idof+1];
	  temp2 = shpfun_ref[igauss*4*num_shap+4*idof+2];
	  temp3 = shpfun_ref[igauss*4*num_shap+4*idof+3];
	  
	  
	  // compute derivatives wrt global coordinates
	  // 15 operations
	  SCALAR fun_u_derx = temp1*jac_0 + temp2*jac_3 + temp3*jac_6;
	  SCALAR fun_u_dery = temp1*jac_1 + temp2*jac_4 + temp3*jac_7;
	  SCALAR fun_u_derz = temp1*jac_2 + temp2*jac_5 + temp3*jac_8;
	  
  #ifdef COUNT_OPER
	  nr_access_shared += 4; // constant cache accesses
	  nr_oper += 15; // after optimization
	  // total: 13+5+18*num_geo_dofs+15+36+15*num_shap = 177+90 = 267 (for prisms)
  #endif

#endif // end if not COMPUTE_ALL_SHAPE_FUN_DER

//*** THE END OF: SUBSTITUTING OR COMPUTING GLOBAL DERIVATIVES OF IDOF SHAPE FUNCTION ***//
//-------------------------------------------------------------

//-------------------------------------------------------------
//*** ACTUAL INTERMEDIATE CALCULATIONS FOR IDOF SHAPE FUNCTION ***//

#ifdef LAPLACE

	  temp4=fun_u_derx;
	  temp5=fun_u_dery;
	  temp6=fun_u_derz;

#elif defined(TEST_NUMINT)

  #ifdef USE_WORKSPACE_FOR_PDE_COEFF

	  temp4 = workspace[offset+0]  * fun_u_derx +
	          workspace[offset+1]  * fun_u_dery +
	          workspace[offset+2]  * fun_u_derz +
	          workspace[offset+12] * shp_fun_u ;
	  
	  temp5 = workspace[offset+3]  * fun_u_derx +
	          workspace[offset+4]  * fun_u_dery +
	          workspace[offset+5]  * fun_u_derz +
	          workspace[offset+13] * shp_fun_u;
	  
	  temp6 = workspace[offset+6]  * fun_u_derx +
	          workspace[offset+7]  * fun_u_dery +
	          workspace[offset+8]  * fun_u_derz +
	          workspace[offset+14] * shp_fun_u;
	  
	  temp7 = workspace[offset+9]  * fun_u_derx +
	          workspace[offset+10] * fun_u_dery +
	          workspace[offset+11] * fun_u_derz +
	          workspace[offset+15] * shp_fun_u;
	  
    #ifdef COUNT_OPER
	  nr_access_shared += 16;
    #endif

  #else // if not USE_WORKSPACE_FOR_PDE_COEFF

	  temp4 = pde_coeff[0]*fun_u_derx + pde_coeff[1]*fun_u_dery + pde_coeff[2]*fun_u_derz + pde_coeff[12]*shp_fun_u;
	  temp5 = pde_coeff[3]*fun_u_derx + pde_coeff[4]*fun_u_dery + pde_coeff[5]*fun_u_derz + pde_coeff[13]*shp_fun_u;
	  temp6 = pde_coeff[6]*fun_u_derx + pde_coeff[7]*fun_u_dery + pde_coeff[8]*fun_u_derz + pde_coeff[14]*shp_fun_u;
	  temp7 = pde_coeff[9]*fun_u_derx + pde_coeff[10]*fun_u_dery + pde_coeff[11]*fun_u_derz + pde_coeff[15]*shp_fun_u;

	  //!!! OLD BACKUP
	  //temp4 = coeff00*fun_u_derx + coeff01*fun_u_dery + coeff02*fun_u_derz + coeff03*shp_fun_u;
	  //temp5 = coeff10*fun_u_derx + coeff11*fun_u_dery + coeff12*fun_u_derz + coeff13*shp_fun_u;
	  //temp6 = coeff20*fun_u_derx + coeff21*fun_u_dery + coeff22*fun_u_derz + coeff23*shp_fun_u;
	  //temp7 = coeff30*fun_u_derx + coeff31*fun_u_dery + coeff32*fun_u_derz + coeff33*shp_fun_u;
	  //!!! OLD BACKUP

  #endif // if not USE_WORKSPACE_FOR_PDE_COEFF
	
  #ifdef COUNT_OPER
	  nr_oper += 7*4;
  #endif

#elif defined(HEAT)

          #ifdef USE_WORKSPACE_FOR_PDE_COEFF

	    temp4 = workspace[offset+0]  * fun_u_derx +
	            workspace[offset+1]  * fun_u_dery +
	            workspace[offset+2]  * fun_u_derz +
	            workspace[offset+12] * shp_fun_u ;
	  
	    temp5 = workspace[offset+3]  * fun_u_derx +
	            workspace[offset+4]  * fun_u_dery +
	           workspace[offset+5]  * fun_u_derz +
	           workspace[offset+13] * shp_fun_u;
	  
	    temp6 = workspace[offset+6]  * fun_u_derx +
	            workspace[offset+7]  * fun_u_dery +
	            workspace[offset+8]  * fun_u_derz +
	            workspace[offset+14] * shp_fun_u;
	  
	    temp7 = workspace[offset+9]  * fun_u_derx +
	            workspace[offset+10] * fun_u_dery +
	            workspace[offset+11] * fun_u_derz +
	            workspace[offset+15] * shp_fun_u;
	  
            #ifdef COUNT_OPER
	      nr_access_shared += 16;
            #endif

          #else // if not USE_WORKSPACE_FOR_PDE_COEFF

	    temp4 = pde_coeff[0]*fun_u_derx + pde_coeff[1]*fun_u_dery + pde_coeff[2]*fun_u_derz + pde_coeff[12]*shp_fun_u;
	    temp5 = pde_coeff[3]*fun_u_derx + pde_coeff[4]*fun_u_dery + pde_coeff[5]*fun_u_derz + pde_coeff[13]*shp_fun_u;
	    temp6 = pde_coeff[6]*fun_u_derx + pde_coeff[7]*fun_u_dery + pde_coeff[8]*fun_u_derz + pde_coeff[14]*shp_fun_u;
	    temp7 = pde_coeff[9]*fun_u_derx + pde_coeff[10]*fun_u_dery + pde_coeff[11]*fun_u_derz + pde_coeff[15]*shp_fun_u;

          #endif // if not USE_WORKSPACE_FOR_PDE_COEFF
	
          #ifdef COUNT_OPER
	    nr_oper += 7*4;
          #endif	  

#endif

//*** THE END OF: ACTUAL INTERMEDIATE CALCULATIONS FOR IDOF SHAPE FUNCTION ***//
//-------------------------------------------------------------

//-------------------------------------------------------------
//*** ACTUAL CALCULATIONS FOR LOAD VECTOR (AND IDOF SHAPE FUNCTION) ***//

#ifdef LOAD_VEC_COMP

  #ifdef USE_WORKSPACE_FOR_STIFF_MAT

    #ifdef COUNT_OPER
	  nr_access_shared += 2;
    #endif

	  workspace[thread_id*(num_dofs*(num_dofs+1)+PADDING)+num_dofs*num_dofs+idof] += (

  #else

	  load_vec[idof] += (

  #endif

  #ifdef LAPLACE

			     pde_coeff[12] * shp_fun_u

  #elif defined(TEST_NUMINT)

    #ifdef USE_WORKSPACE_FOR_PDE_COEFF
			     
			     workspace[offset+16] * fun_u_derx +
			     workspace[offset+17] * fun_u_dery +
			     workspace[offset+18] * fun_u_derz +
			     workspace[offset+19] * shp_fun_u

    #else // if not USE_WORKSPACE_FOR_PDE_COEFF

			     pde_coeff[16] * fun_u_derx +
			     pde_coeff[17] * fun_u_dery +
			     pde_coeff[18] * fun_u_derz +
			     pde_coeff[19] * shp_fun_u


    #endif // end if not USE_WORKSPACE_FOR_PDE_COEFF

  #elif defined(HEAT)

              #ifdef USE_WORKSPACE_FOR_PDE_COEFF
			     
	        workspace[offset+16] * fun_u_derx +
		workspace[offset+17] * fun_u_dery +
		workspace[offset+18] * fun_u_derz +
		workspace[offset+19] * shp_fun_u

              #else // if not USE_WORKSPACE_FOR_PDE_COEFF

	        pde_coeff[16] * fun_u_derx +
		pde_coeff[17] * fun_u_dery +
		pde_coeff[18] * fun_u_derz +
		pde_coeff[19] * shp_fun_u

              #endif // end if not USE_WORKSPACE_FOR_PDE_COEFF

  #endif
			     
			     ) * vol;
	  
  #ifdef COUNT_OPER
    #ifdef LAPLACE

	  nr_oper += 3;

    #elif defined(TEST_NUMINT)

      #ifdef USE_WORKSPACE_FOR_PDE_COEFF
	  nr_access_shared += 4; 
      #endif
	  nr_oper += 9;

    #elif defined(HEAT)

	  #ifdef USE_WORKSPACE_FOR_PDE_COEFF
	      nr_access_shared += 4; 
          #endif
	      nr_oper += 9;

    #endif

  #endif
    
#endif // end if computing RHS vector

//*** THE END OF: ACTUAL CALCULATIONS FOR LOAD VECTOR (AND IDOF SHAPE FUNCTION) ***//
//-------------------------------------------------------------

	  } // the end of using registers for u (shp_fun_u, fun_u_der.)

//-------------------------------------------------------------
// ************************* second loop over shape functions ****************************//
        for(jdof = 0; jdof < num_shap; jdof++){
	  
//-------------------------------------------------------------
//****** SUBSTITUTING OR COMPUTING GLOBAL DERIVATIVES OF JDOF SHAPE FUNCTION ******//

#ifdef COMPUTE_ALL_SHAPE_FUN_DER
	  
  #ifdef USE_WORKSPACE_FOR_SHAPE_FUN

	  SCALAR shp_fun_v = shpfun_ref[igauss*4*num_shap+4*jdof];
	  SCALAR fun_v_derx = workspace[thread_id*(3*num_shap+PADDING)+3*jdof];
	  SCALAR fun_v_dery = workspace[thread_id*(3*num_shap+PADDING)+3*jdof+1];
	  SCALAR fun_v_derz = workspace[thread_id*(3*num_shap+PADDING)+3*jdof+2];

    #ifdef COUNT_OPER
	  nr_access_shared += 4;  // including 1 for constant cache
    #endif

  #else // if not USE_WORKSPACE_FOR_SHAPE_FUN
 
	  SCALAR shp_fun_v = shpfun_ref[igauss*4*num_shap+4*jdof];
	  SCALAR fun_v_derx = tab_fun_u_derx[jdof];
	  SCALAR fun_v_dery = tab_fun_u_dery[jdof];
	  SCALAR fun_v_derz = tab_fun_u_derz[jdof];

    #ifdef COUNT_OPER
	  nr_access_shared += 1; // constant cache access
    #endif
	  
  #endif // end if not USE_WORKSPACE_FOR_SHAPE_FUN
	  
#else // if not COMPUTE_ALL_SHAPE_FUN_DER

	// read proper values of shape functions and their derivatives
	SCALAR shp_fun_v = shpfun_ref[igauss*4*num_shap+4*jdof];
	temp1 = shpfun_ref[igauss*4*num_shap+4*jdof+1];
	temp2 = shpfun_ref[igauss*4*num_shap+4*jdof+2];
	temp3 = shpfun_ref[igauss*4*num_shap+4*jdof+3];
	
	// compute derivatives wrt global coordinates
	// 15 operations
	SCALAR fun_v_derx = temp1*jac_0 + temp2*jac_3 + temp3*jac_6;
	SCALAR fun_v_dery = temp1*jac_1 + temp2*jac_4 + temp3*jac_7;
	SCALAR fun_v_derz = temp1*jac_2 + temp2*jac_5 + temp3*jac_8;
	
  #ifdef COUNT_OPER
	nr_access_shared += 4; // constant cache accesses
	nr_oper += 15; 
  #endif
	  
#endif // end if not COMPUTE_ALL_SHAPE_FUN_DER

//*** THE END OF: SUBSTITUTING OR COMPUTING GLOBAL DERIVATIVES OF IDOF SHAPE FUNCTION ***//
//-------------------------------------------------------------

//-------------------------------------------------------------
//********* ACTUAL FINAL CALCULATIONS FOR SM ENTRY  *********//

#ifdef USE_WORKSPACE_FOR_STIFF_MAT

  #ifdef COUNT_OPER
	nr_access_shared += 2;
  #endif

	workspace[thread_id*(num_dofs*(num_dofs+1)+PADDING)+idof*num_dofs+jdof] += (

#else

	stiff_mat[idof*num_dofs+jdof] += (

#endif

#ifdef LAPLACE

      	    temp4 * fun_v_derx +
       	    temp5 * fun_v_dery +
       	    temp6 * fun_v_derz

#elif defined(TEST_NUMINT)

      	    temp4 * fun_v_derx +
       	    temp5 * fun_v_dery +
       	    temp6 * fun_v_derz +
       	    temp7 * shp_fun_v

#elif defined(HEAT)

      	    temp4 * fun_v_derx +
       	    temp5 * fun_v_dery +
       	    temp6 * fun_v_derz +
       	    temp7 * shp_fun_v

#endif

					    ) * vol;

#ifdef COUNT_OPER
  #ifdef LAPLACE
	nr_oper += 7; 
  #elif defined(TEST_NUMINT)
	nr_oper += 9;
  #elif defined(HEAT)
	nr_oper += 9;
  #endif
#endif

//*** THE END OF: ACTUAL FINAL CALCULATIONS FOR SM ENTRY  ***//
//-------------------------------------------------------------

	}//jdof

//******* THE END OF: first loop over shape functions *******//
//-------------------------------------------------------------

      }//idof

//******* THE END OF: second loop over shape functions *******//
//-------------------------------------------------------------

    }//gauss

// ******** THE END OF: loop over integration points ********//
//-------------------------------------------------------------

//-------------------------------------------------------------
//******** WRITING OUTPUT SM AND LV TO GLOBAL MEMORY ********//


#ifdef ASSEMBLING

      // Write index
      int write_index = 0;

      // Write data into crs, perform assembling
      if(element_index < nr_elems_this_kercall) { // perform only for real elements
	i=0;
	for(idof=0; idof < num_shap; idof++)
	  {
	    for(jdof=0; jdof < num_shap; jdof++)
	      {
		// Calculate write index, read from the assembly table and position
		write_index = assembly_table[assembly_position_table[first_elem_index+element_index]+i];
			
		/*jbw
		printf("GPU tmr_ocl_num_int_el.cl %d element_index=%d ; assembly_position_table[%d]=%d ; assembly_table[%d] = %d\n",
		       __LINE__,element_index,
		       (first_elem_index+element_index),assembly_position_table[first_elem_index+element_index],
		       (assembly_position_table[first_elem_index+element_index]+i),write_index);
		/*jbw*/

               #ifdef USE_WORKSPACE_FOR_STIFF_MAT

		/*jbw
		int workspace_index = thread_id*(num_dofs*(num_dofs+1)+PADDING)+idof*num_dofs+jdof;
		double workspace_val = workspace[workspace_val];
		printf("GPU tmr_ocl_num_int_el.cl %d element_index=%d ; assembly_position_table[%d]=%d ; assembly_table[%d] = %d ; written value workspace[%d] = %lf\n",
		       __LINE__,element_index,
		       (first_elem_index+element_index),assembly_position_table[first_elem_index+element_index],
		       (assembly_position_table[first_elem_index+element_index]+i),write_index,workspace_index,workspace_val);
		/*jbw*/

		crs[write_index] += workspace[thread_id*(num_dofs*(num_dofs+1)+PADDING)+idof*num_dofs+jdof];

	       #else

		/*jbw
		if(element_index == 0) {
		int stiff_mat_index = idof*num_dofs+jdof;
		double stiff_mat_val = stiff_mat[stiff_mat_index];
		printf("GPU tmr_ocl_num_int_el.cl GPU element_index=%d ; assembly_position_table[%d]=%d ; i=%d ;assembly_table[%d] = %d ; written value stiff_mat[%d] = %lf\n",
		       element_index,
		       (first_elem_index+element_index),assembly_position_table[first_elem_index+element_index],i,
		       (assembly_position_table[first_elem_index+element_index]+i),write_index,stiff_mat_index,stiff_mat_val);
		}
		/*jbw*/

		crs[write_index] += stiff_mat[idof*num_dofs+jdof];

               #endif
			
		i++;
	      }
	  }

	  #ifdef COUNT_OPER
	    nr_global_access += 2*num_dofs*num_dofs; // we neglect shared memory accesses at this stage (read and write)
          #endif

	  #ifdef LOAD_VEC_COMP
	    /* currently disabled write load vector */

	    for(i=0; i < num_shap; i++) { // write load vector

	      // Calculate write index, read from the assembly table and position
	      write_index = assembly_table[assembly_position_table[first_elem_index+element_index]+num_shap*num_shap+i];

              #ifdef USE_WORKSPACE_FOR_STIFF_MAT

		rhs[write_index] += workspace[thread_id*(num_dofs*(num_dofs+1)+PADDING)+num_dofs*num_dofs+i];

              #else

	        rhs[write_index] += load_vec[i];

             #endif

	    }

            #ifdef COUNT_OPER
	      nr_global_access += 2*num_dofs; // we neglect shared memory accesses at this stage (read and write)
            #endif
	  
	   /* currently disabled write load vector */
        #endif // end if LOAD_VEC_COMP
      }

#else // if assembling disabled
	    
#ifdef COAL_WRITE

    // write stiffness matrix - in a coalesced way
    offset = (element_index-thread_id)*(num_shap*num_shap+num_shap);
    i=0;
    for(idof=0; idof < num_shap; idof++)
      {
    	for(jdof=0; jdof < num_shap; jdof++)
	  {

  #ifdef USE_WORKSPACE_FOR_STIFF_MAT

	    stiff_mat_out[offset+i*WORK_GROUP_SIZE+thread_id] =
	      workspace[thread_id*(num_dofs*(num_dofs+1)+PADDING)+idof*num_dofs+jdof];

  #else

	    stiff_mat_out[offset+i*WORK_GROUP_SIZE+thread_id] = stiff_mat[idof*num_dofs+jdof];

  #endif

	    i++;
	  }
      }

  #ifdef COUNT_OPER
    nr_global_access += num_dofs*num_dofs; // we neglect shared memory accesses at this stage
  #endif


  #ifdef LOAD_VEC_COMP

    for(i=0; i < num_shap; i++){
      // write load vector

    #ifdef USE_WORKSPACE_FOR_STIFF_MAT

      stiff_mat_out[offset+(num_shap*num_shap+i)*WORK_GROUP_SIZE+thread_id] =
	workspace[thread_id*(num_dofs*(num_dofs+1)+PADDING)+num_dofs*num_dofs+i];

    #else

      stiff_mat_out[offset+(num_shap*num_shap+i)*WORK_GROUP_SIZE+thread_id] = load_vec[i];

    #endif

    }

    #ifdef COUNT_OPER
    nr_global_access += num_dofs; // we neglect shared memory accesses at this stage
    #endif

  #endif // end if not LOAD_VEC_COMP

#else // if not coalesced write

    // write stiffness matrix - threads compute subsequent elements
    offset = element_index*(num_shap*num_shap+num_shap);
    i=0;
    for(idof=0; idof < num_shap; idof++)
      {
	for(jdof=0; jdof < num_shap; jdof++)
	  {

  #ifdef USE_WORKSPACE_FOR_STIFF_MAT

	    stiff_mat_out[offset+i] =
	      workspace[thread_id*(num_dofs*(num_dofs+1)+PADDING)+idof*num_dofs+jdof];

  #else

	    stiff_mat_out[offset+i] = stiff_mat[idof*num_dofs+jdof];

  #endif

	    i++;
	  }
      }

  #ifdef COUNT_OPER
    nr_global_access += num_dofs*num_dofs; // we neglect shared memory accesses at this stage
  #endif

  #ifdef LOAD_VEC_COMP

    for(i=0; i < num_shap; i++){
      // write load vector

    #ifdef USE_WORKSPACE_FOR_STIFF_MAT

      stiff_mat_out[offset+num_shap*num_shap+i] =
	workspace[thread_id*(num_dofs*(num_dofs+1)+PADDING)+num_dofs*num_dofs+i];

    #else

      stiff_mat_out[offset+num_shap*num_shap+i] = load_vec[i];

    #endif

    }

    #ifdef COUNT_OPER
    nr_global_access += num_dofs; // we neglect shared memory accesses at this stage
    #endif

  #endif // end if LOAD_VEC_COMP

#endif // end if not COAL_WRITE

barrier(CLK_GLOBAL_MEM_FENCE);

#endif // end if not ASSEMBLING
    
// *** THE END OF: WRITING OUTPUT SM AND LV TO GLOBAL MEMORY ***//
//-------------------------------------------------------------

  } // the end of loop over elements

// ************* THE END OF: LOOP OVER ELEMENTS *************//
//-------------------------------------------------------------

#ifdef COUNT_OPER

  if(group_id==0 && thread_id==0){
    kernel_profile_data_out[0] = nr_oper;
    kernel_profile_data_out[1] = nr_access_shared;
    kernel_profile_data_out[2] = nr_global_access;
  }

#endif


};


kernel void tmr_ocl_prepare_final_crs(
  // execution_parameters can be read directly from constant memory, assuming it is cached and
  // further accesses are realized from cache
  __constant int* execution_parameters,

  // External CRS data
  __global SCALAR* rhs_external,  //in: 
  __global SCALAR* crs_external,  //in:
  
  // Internal CRS data
  __global SCALAR* rhs,  //out: 
  __global SCALAR* crs  //out:  

  // __global int *crs_col_ind, //in: not used in this kernel
  // __global int *crs_row_ptr //in: not used in this kernel

  #ifdef COUNT_OPER
  , // buffer give access to kernel profiling data such as number of arthmetic operations and shared/global memory access
  __global SCALAR* kernel_profile_data_out

  #endif
  
){

#ifdef COUNT_OPER
  SCALAR nr_oper=0.0;
  SCALAR nr_access_shared=0.0;
  SCALAR nr_global_access=0.0;
#endif

  const int group_id = get_group_id(0);
  const int thread_id = get_local_id(0);
  const int thread_global_id = get_global_id(0);
  const int work_group_size = get_local_size(0);
  const int nr_work_groups = get_num_groups(0);

  // Get CRS and RHS size
  int rhs_size = execution_parameters[0];
  int crs_size = execution_parameters[1];

  int i;

#ifdef COAL_WRITE

  int index_start;
  int index_end;
  int stride;

  index_start = thread_global_id;
  stride = nr_work_groups*work_group_size;

  //Perform addition input and internal crs
  index_end = crs_size;

  for(i=index_start; i<index_end; i+=stride) {
    crs[i] += crs_external[i];

#ifdef COUNT_OPER
    nr_oper = nr_oper+2;
    nr_global_access = nr_global_access+3;
#endif   
  }

  //Perform addition input and internal rhs
  index_end = rhs_size;
  for(i=index_start; i<index_end; i+=stride) {
    rhs[i] += rhs_external[i];

#ifdef COUNT_OPER
    nr_oper = nr_oper+2;
    nr_global_access = nr_global_access+3;
#endif 
  }

#else
    
  int crs_nr_index_per_thread;
  int crs_offset;
  int rhs_nr_index_per_thread;
  int rhs_offset;
  int index;

  // Get number of elements per thread and offset size
  rhs_nr_index_per_thread = execution_parameters[2];
  rhs_offset = execution_parameters[3];
  crs_nr_index_per_thread = execution_parameters[4];
  crs_offset = execution_parameters[5];
  
  
  //Perform addition input and internal crs
  index = (group_id*work_group_size+thread_id)*crs_nr_index_per_thread;

  for(i=0; i<crs_nr_index_per_thread; i++, index++) {

	//printf("tmr_ocl_prepare_final_crs -> CRS -> group_id=%d ; thread_id=%d ; nr_index_per_thread=%d/%d (offset=%d); index=%d ; CRS_internal = %lf ; CRS_external=%lf\n",group_id,thread_id,crs_nr_index_per_thread,(nr_work_groups*work_group_size),crs_offset,index,crs[index],crs_external[index]);
	
	crs[index] += crs_external[index];

#ifdef COUNT_OPER
	nr_oper = nr_oper+2;
	nr_global_access = nr_global_access+3;
#endif 
  }

  index = (nr_work_groups*work_group_size*crs_nr_index_per_thread)+group_id*work_group_size+thread_id;
  if(index < crs_size) {

    //printf("tmr_ocl_prepare_final_crs -> CRS (OFFSET %d) -> group_id=%d ; thread_id=%d ; nr_index_per_thread=%d/%d (offset=%d); index=%d ; CRS_internal = %lf ; CRS_external=%lf\n",crs_offset,group_id,thread_id,crs_nr_index_per_thread,(nr_work_groups*work_group_size),crs_offset,index,crs[index],crs_external[index]);
    
	crs[index] += crs_external[index];

#ifdef COUNT_OPER
	nr_oper = nr_oper+2;
	nr_global_access = nr_global_access+3;
#endif 
  }


  //Perform addition input and internal rhs
  index = (group_id*work_group_size+thread_id)*rhs_nr_index_per_thread;

  for(i=0; i<rhs_nr_index_per_thread; i++, index++) {

    //printf("tmr_ocl_prepare_final_crs -> RHS -> group_id=%d ; thread_id=%d ; nr_index_per_thread=%d/%d (offset=%d); index=%d ; RHS_internal = %lf ; RHS_external=%lf\n",group_id,thread_id,rhs_nr_index_per_thread,(nr_work_groups*work_group_size),rhs_offset,index,rhs[index],rhs_external[index]);
	
	rhs[index] += rhs_external[index];
	
#ifdef COUNT_OPER
	nr_oper = nr_oper+2;
	nr_global_access = nr_global_access+3;
#endif 	
  }

  index = (nr_work_groups*work_group_size*rhs_nr_index_per_thread)+group_id*work_group_size+thread_id;
  if(index < rhs_size) {

    //printf("tmr_ocl_prepare_final_crs -> RHS (OFFSET %d) -> group_id=%d ; thread_id=%d ; nr_index_per_thread=%d/%d (offset=%d); index=%d ; RHS_internal = %lf ; RHS_external=%lf\n",rhs_offset,group_id,thread_id,rhs_nr_index_per_thread,(nr_work_groups*work_group_size),rhs_offset,index,rhs[index],rhs_external[index]);
    
	rhs[index] += rhs_external[index];
	
#ifdef COUNT_OPER
	nr_oper = nr_oper+2;
	nr_global_access = nr_global_access+3;
#endif 	
  }

#endif

  #ifdef COUNT_OPER

  barrier(CLK_GLOBAL_MEM_FENCE);

  if(group_id==0 && thread_id==0){
    kernel_profile_data_out[0] = nr_oper;
    kernel_profile_data_out[1] = nr_access_shared;
    kernel_profile_data_out[2] = nr_global_access;
  }

  #endif
  
  // Memory synchronization before solve
  //barrier(CLK_GLOBAL_MEM_FENCE);	

}


