//QSS
//!!!!!!!!!!! TODO: SEPARATE GEO_DATA AND PDE_COEFF !!!!!!!!!!!!!!!!

//!!!!!!!!!!! TEST DIFFERENT CONFIGURATIONS OF SHARED / L1 MEMORY

//!!!!!!!!!!! LOOP UNROLLING = REGISTER BLOCKING FOR INNERMOST LOOPS ??????????

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

// FOR LINEAR TETRAHEDRA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//#define num_shap 4
//#define num_gauss 4
//#define num_geo_dofs 4
//#define weight_linear_tetra (one_fourth*one_sixth)

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
  #define NR_COEFFS_SENT_PER_ELEMENT (num_dofs+1)
  #define NR_COEFFS_SENT_PER_INT_POINT 5
  #define NR_COEFFS_IN_SM_CALCULATIONS 20
#endif

#define NR_PDE_COEFFS_SENT (NR_COEFFS_SENT_PER_ELEMENT + NR_COEFFS_SENT_PER_INT_POINT*num_gauss)
//----------------------------------------------------


//----------------------------------------------------
// SWITCHES FOR DIFFERENT OPTIMIZATION OPTIONS !!!!!!!!!!!!!!!!!!!!!!!!

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
#define USE_WORKSPACE_FOR_SHAPE_FUN


// USE_WORKSPACE_FOR_STIFF_MAT - to use shared memory for SM during SM calculations
// otherwise - registers are used
//#define USE_WORKSPACE_FOR_STIFF_MAT

// FOR EACH ARCHITECTURE PADDING SHOULD BE TESTED TO DETECT SHARED MEMORY BANK CONFLICTS
#define WORKSPACE_PADDING 0
#define PADDING WORKSPACE_PADDING

// THE END OF: SWITCHES FOR DIFFERENT OPTIMIZATION OPTIONS
//----------------------------------------------------


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
  // TODO!!!!!!!!!!!!!!!!!!
  //__global SCALAR* el_data_in_geo, // geo data for integration of NR_ELEMS_THIS_KERCALL elements
  //__global SCALAR* el_data_in_coeff, // coeff data for integration of elements
  __global SCALAR* stiff_mat_out // result of integration of NR_ELEMS_THIS_KERCALL elements
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
  #define WORKSPACE_SIZE_FOR_PDE_COEFF (NR_COEFFS_IN_SM_CALCULATIONS)
#else
  #define WORKSPACE_SIZE_FOR_PDE_COEFF 0
  //SCALAR pde_coeff[NR_COEFFS_IN_SM_CALCULATIONS];
  #ifdef LAPLACE
  SCALAR pde_coeff[num_shap];
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

    offset= nr_elems_this_kercall * EL_GEO_DAT_SIZE +
      // TODO !!!!!!!!!!!!!!
      //offset= 
                       (element_index - thread_id) * NR_PDE_COEFFS_SENT;

    barrier(CLK_LOCAL_MEM_FENCE); // I don't know why but without barrier here, one per ten runs gives bad results

    // TRY TO UNROLL THIS LOOP TO INCREASE MEMORY PARALLELISM !!!
    for(i=0;i<NR_PDE_COEFFS_SENT;i++) {

      workspace[i*WORK_GROUP_SIZE+thread_id] = 
	                           el_data_in[offset+i*WORK_GROUP_SIZE+thread_id];
      // TODO!!!!!!!!!!!!!!
      //el_data_in_coeff[offset+i*WORK_GROUP_SIZE+thread_id];

    }

    #ifdef COUNT_OPER
    nr_global_access += NR_PDE_COEFFS_SENT; // we neglect shared memory accesses at this stage
    #endif

    barrier(CLK_LOCAL_MEM_FENCE); // !!!!!!!!!!!!!!!!!!!!!!

    offset=thread_id*(NR_PDE_COEFFS_SENT);

    #ifdef LAPLACE

    pde_coeff[0]=workspace[offset+0];
    pde_coeff[1]=workspace[offset+1];
    pde_coeff[2]=workspace[offset+2];
    pde_coeff[3]=workspace[offset+3];
    pde_coeff[4]=workspace[offset+4];
    pde_coeff[5]=workspace[offset+5];

    #elif defined(TEST_NUMINT)

    SCALAR coeff00=workspace[offset+0];
    SCALAR coeff01=workspace[offset+1];
    SCALAR coeff02=workspace[offset+2];
    SCALAR coeff10=workspace[offset+3];
    SCALAR coeff11=workspace[offset+4];
    SCALAR coeff12=workspace[offset+5];
    SCALAR coeff20=workspace[offset+6];
    SCALAR coeff21=workspace[offset+7];
    SCALAR coeff22=workspace[offset+8];
    SCALAR coeff30=workspace[offset+9];
    SCALAR coeff31=workspace[offset+10];
    SCALAR coeff32=workspace[offset+11];
    SCALAR coeff03=workspace[offset+12];
    SCALAR coeff13=workspace[offset+13];
    SCALAR coeff23=workspace[offset+14];
    SCALAR coeff33=workspace[offset+15];
    SCALAR coeff04=workspace[offset+16];
    SCALAR coeff14=workspace[offset+17];
    SCALAR coeff24=workspace[offset+18];
    SCALAR coeff34=workspace[offset+19];

    #elif defined(HEAT)

    #endif

    #ifdef COUNT_OPER
    nr_access_shared += NR_PDE_COEFFS_SENT; // we count accesses because of the barriers
    #endif

    barrier(CLK_LOCAL_MEM_FENCE); // !!!!!!!!!!!!!!!!!!!!!!

    // after rewriting PDE_COEFF to registers we read GEO_DATA to workspace
    offset = (element_index-thread_id)*(EL_GEO_DAT_SIZE);

    // TRY TO UNROLL THIS LOOP TO INCREASE MEMORY PARALLELISM !!!
    for(i = 0; i < EL_GEO_DAT_SIZE; i++){

      workspace[i*WORK_GROUP_SIZE+thread_id] = 
	                         el_data_in[offset+i*WORK_GROUP_SIZE+thread_id];
	// TODO !!!!
	// el_data_in_geo[offset+i*WORK_GROUP_SIZE+thread_id];

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
	// TODO !!!!
	// el_data_in_geo[offset+i*WORK_GROUP_SIZE+thread_id];
      
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
    offset= nr_elems_this_kercall * EL_GEO_DAT_SIZE +
      // TODO !!!!!!!!!!!!!!
      //offset= 
                    (element_index - thread_id) * NR_PDE_COEFFS_SENT;

    // TRY TO UNROLL THIS LOOP TO INCREASE MEMORY PARALLELISM !!!
    for(i=0;i<NR_PDE_COEFFS_SENT;i++) {

      workspace[i*WORK_GROUP_SIZE+thread_id] = 
	el_data_in[offset+i*WORK_GROUP_SIZE+thread_id];
	// TODO !!!!
	//el_data_in_coeff[offset+i*WORK_GROUP_SIZE+thread_id];

    }

    #ifdef COUNT_OPER
    nr_global_access += NR_PDE_COEFFS_SENT; // we neglect shared memory accesses at this stage
    #endif

    // if we do not leave pde coeffs in workspace we rewrite them to registers
    #ifndef USE_WORKSPACE_FOR_PDE_COEFF

    barrier(CLK_LOCAL_MEM_FENCE); // !!!!!!!!!!!!!!!!!!!!!!

    offset=thread_id*(NR_PDE_COEFFS_SENT);

      #ifdef LAPLACE

    pde_coeff[0]=workspace[offset+0];
    pde_coeff[1]=workspace[offset+1];
    pde_coeff[2]=workspace[offset+2];
    pde_coeff[3]=workspace[offset+3];
    pde_coeff[4]=workspace[offset+4];
    pde_coeff[5]=workspace[offset+5];

      #elif defined(TEST_NUMINT)

    SCALAR coeff00=workspace[offset+0];
    SCALAR coeff01=workspace[offset+1];
    SCALAR coeff02=workspace[offset+2];
    SCALAR coeff10=workspace[offset+3];
    SCALAR coeff11=workspace[offset+4];
    SCALAR coeff12=workspace[offset+5];
    SCALAR coeff20=workspace[offset+6];
    SCALAR coeff21=workspace[offset+7];
    SCALAR coeff22=workspace[offset+8];
    SCALAR coeff30=workspace[offset+9];
    SCALAR coeff31=workspace[offset+10];
    SCALAR coeff32=workspace[offset+11];
    SCALAR coeff03=workspace[offset+12];
    SCALAR coeff13=workspace[offset+13];
    SCALAR coeff23=workspace[offset+14];
    SCALAR coeff33=workspace[offset+15];
    SCALAR coeff04=workspace[offset+16];
    SCALAR coeff14=workspace[offset+17];
    SCALAR coeff24=workspace[offset+18];
    SCALAR coeff34=workspace[offset+19];

      #elif defined(HEAT)

      #endif

      #ifdef COUNT_OPER
    nr_access_shared += NR_PDE_COEFFS_SENT; // we count accesses because of the barriers
      #endif

    #endif // end if we do not leave pde coeffs in workspace and rewrite them to registers

  #endif // end if workspace not used for geo_data (hence used for pde coeffs or something else)

#else // if not COAL_READ (i.e. only one or none of PDE_COEFF and GEO_DAT read in a coalsced way)


  #ifdef USE_WORKSPACE_FOR_PDE_COEFF

    offset= nr_elems_this_kercall * EL_GEO_DAT_SIZE +
      // TODO !!!!!!!!!!!!!!
      //offset= 
                       (element_index - thread_id) * NR_PDE_COEFFS_SENT;
    // TRY TO UNROLL THIS LOOP TO INCREASE MEMORY PARALLELISM !!!
    for(i=0;i<NR_PDE_COEFFS_SENT;i++) {

      workspace[i*WORK_GROUP_SIZE+thread_id] = 
	                           el_data_in[offset+i*WORK_GROUP_SIZE+thread_id];
      // TODO!!!!!!!!!!!!!!
      //el_data_in_coeff[offset+i*WORK_GROUP_SIZE+thread_id];

    }

    #ifdef COUNT_OPER
    nr_global_access += NR_PDE_COEFFS_SENT; // we neglect shared memory accesses at this stage
    #endif

  #else // if not USE_WORKSPACE_FOR_PDE_COEFF

    offset = nr_elems_this_kercall * EL_GEO_DAT_SIZE +
      // TODO !!!!!!!!!!!!!!
      //offset= 
	  			element_index * (NR_PDE_COEFFS_SENT);

    #ifdef LAPLACE

    pde_coeff[0]=el_data_in[offset+0];
    pde_coeff[1]=el_data_in[offset+1];
    pde_coeff[2]=el_data_in[offset+2];
    pde_coeff[3]=el_data_in[offset+3];
    pde_coeff[4]=el_data_in[offset+4];
    pde_coeff[5]=el_data_in[offset+5];

    #elif defined(TEST_NUMINT)

    SCALAR coeff00=el_data_in[offset+0];
    SCALAR coeff01=el_data_in[offset+1];
    SCALAR coeff02=el_data_in[offset+2];
    SCALAR coeff10=el_data_in[offset+3];
    SCALAR coeff11=el_data_in[offset+4];
    SCALAR coeff12=el_data_in[offset+5];
    SCALAR coeff20=el_data_in[offset+6];
    SCALAR coeff21=el_data_in[offset+7];
    SCALAR coeff22=el_data_in[offset+8];
    SCALAR coeff30=el_data_in[offset+9];
    SCALAR coeff31=el_data_in[offset+10];
    SCALAR coeff32=el_data_in[offset+11];
    SCALAR coeff03=el_data_in[offset+12];
    SCALAR coeff13=el_data_in[offset+13];
    SCALAR coeff23=el_data_in[offset+14];
    SCALAR coeff33=el_data_in[offset+15];
    SCALAR coeff04=el_data_in[offset+16];
    SCALAR coeff14=el_data_in[offset+17];
    SCALAR coeff24=el_data_in[offset+18];
    SCALAR coeff34=el_data_in[offset+19];


    #elif defined(HEAT)

    #endif

    #ifdef COUNT_OPER
    nr_global_access += NR_PDE_COEFFS_SENT;
    #endif

  #endif // end if not USE_WORKSPACE_FOR_PDE_COEFF

  #ifdef USE_WORKSPACE_FOR_GEO_DATA

    offset = (element_index-thread_id)*(EL_GEO_DAT_SIZE);

    barrier(CLK_LOCAL_MEM_FENCE); // I don't know why but without barrier here, one per ten runs gives bad results

    // TRY TO UNROLL THIS LOOP TO INCREASE MEMORY PARALLELISM !!!
    for(i = 0; i < EL_GEO_DAT_SIZE; i++){

      workspace[i*WORK_GROUP_SIZE+thread_id] = 
                    el_data_in[offset+i*WORK_GROUP_SIZE+thread_id];
      // TODO !!!!!!!!!
      // el_data_in_geo[offset+i*WORK_GROUP_SIZE+thread_id];

    }

  #else // if not USE_WORKSPACE_FOR_GEO_DATA

    // we read geo data to registers
    offset = (element_index)*(EL_GEO_DAT_SIZE);

    // TRY TO UNROLL THIS LOOP TO INCREASE MEMORY PARALLELISM !!!
    for(i = 0; i < EL_GEO_DAT_SIZE; i++){
      
      geo_dat[i] = el_data_in[offset+i];
      // TODO !!!!!!!!!
      // geo_dat[i] = el_data_in_geo[offset+i];
      
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
      SCALAR vol = weight_linear_prism;	// vol = weight CONSTANT FOR LINEAR PRISM!!!

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

#ifdef HEAT

      // for non-constant, non-linear coefficients a place for call to problem dependent
      // function calculating actual PDE coefficients based on data in coeff 
      // workspace or registers and  storing data back in workspace or in registers

#endif


#ifdef USE_WORKSPACE_FOR_PDE_COEFF


  #ifdef LAPLACE

      // offset for reading data
      offset=thread_id*(NR_PDE_COEFFS_SENT);

      SCALAR coeff03 = workspace[offset+igauss];

    #ifdef COUNT_OPER
    nr_access_shared += 1;
    #endif

  #endif // end if LAPLACE

      // offset for computations
      offset=thread_id*(NR_COEFFS_IN_SM_CALCULATIONS);

#else // if not USE_WORKSPACE_FOR_PDE_COEFF

  #ifdef LAPLACE

      //SCALAR coeff03 = pde_coeff[igauss]; // DOES NOT WORK FOR AMD !!!!!!!!!
      SCALAR coeff03 = 0.0;
      switch(igauss){
      case 0:
	coeff03 = pde_coeff[0];
	break;
      case 1:
	coeff03 = pde_coeff[1];
	break;
      case 2:
	coeff03 = pde_coeff[2];
	break;
      case 3:
	coeff03 = pde_coeff[3];
	break;
      case 4:
	coeff03 = pde_coeff[4];
 	break;
      case 5:
	coeff03 = pde_coeff[5];
	break;
      }

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

	  temp4 = coeff00*fun_u_derx + coeff01*fun_u_dery + coeff02*fun_u_derz + coeff03*shp_fun_u;
	  temp5 = coeff10*fun_u_derx + coeff11*fun_u_dery + coeff12*fun_u_derz + coeff13*shp_fun_u;
	  temp6 = coeff20*fun_u_derx + coeff21*fun_u_dery + coeff22*fun_u_derz + coeff23*shp_fun_u;
	  temp7 = coeff30*fun_u_derx + coeff31*fun_u_dery + coeff32*fun_u_derz + coeff33*shp_fun_u;

  #endif // if not USE_WORKSPACE_FOR_PDE_COEFF
	
  #ifdef COUNT_OPER
	  nr_oper += 7*4;
  #endif

#elif defined(HEAT)

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

			     coeff03 * shp_fun_u

  #elif defined(TEST_NUMINT)

    #ifdef USE_WORKSPACE_FOR_PDE_COEFF
			     
			     workspace[offset+16] * fun_u_derx +
			     workspace[offset+17] * fun_u_dery +
			     workspace[offset+18] * fun_u_derz +
			     workspace[offset+19] * shp_fun_u

    #else // if not USE_WORKSPACE_FOR_PDE_COEFF

			     coeff04 * fun_u_derx +
			     coeff14 * fun_u_dery +
			     coeff24 * fun_u_derz +
			     coeff34 * shp_fun_u


    #endif // end if not USE_WORKSPACE_FOR_PDE_COEFF

  #elif defined(HEAT)

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

// *** THE END OF: WRITING OUTPUT SM AND LV TO GLOBAL MEMORY ***//
//-------------------------------------------------------------

  } // the end of loop over elements

// ************* THE END OF: LOOP OVER ELEMENTS *************//
//-------------------------------------------------------------

#ifdef COUNT_OPER

  if(group_id==0 && thread_id==0){
    stiff_mat_out[0] = nr_oper;
    stiff_mat_out[1] = nr_access_shared;
    stiff_mat_out[2] = nr_global_access;
  }

#endif


};
