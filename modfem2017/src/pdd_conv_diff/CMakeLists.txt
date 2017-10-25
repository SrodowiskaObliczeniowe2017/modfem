# Uses:
#
# PDD_CONV_DIFF_EXTERNAL_INCLUDE_DIRS
# PDD_CONV_DIFF_EXTERNAL_LIBRARIES
#
# CREATE_PDD_CONV_DIFF_PRISM_STD
# CREATE_PDD_CONV_DIFF_PRISM_STD_QUAD
# CREATE_PDD_CONV_DIFF_HYBRID_STD
# CREATE_PDD_CONV_DIFF_HYBRID_STD_QUAD
# CREATE_PDD_CONV_DIFF_PRISM_DG
# CREATE_PDD_CONV_DIFF_HYBRID_DG

# -- Add target to global list --
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_CONV_DIFF_PRISM_STD")
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_CONV_DIFF_PRISM_STD_QUAD")
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_CONV_DIFF_HYBRID_STD")
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_CONV_DIFF_HYBRID_STD_QUAD")
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_CONV_DIFF_PRISM2D_STD")
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_CONV_DIFF_PRISM2D_STD_QUAD")
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_CONV_DIFF_PRISM_DG")
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_CONV_DIFF_HYBRID_DG")
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_CONV_DIFF_PRISM2D_DG")
set(MOD_FEM_EXEC_LIST "${MOD_FEM_EXEC_LIST}" PARENT_SCOPE)

# -- Exec files --
set(pdd_conv_diff_srcdir "${CMAKE_SOURCE_DIR}/pdd_conv_diff")
set(PDD_CONV_DIFF_PRISM_STD_EXE "MOD_FEM_conv_diff_prism_std")
set(PDD_CONV_DIFF_PRISM_STD_QUAD_EXE "MOD_FEM_conv_diff_prism_std_quad")
set(PDD_CONV_DIFF_PRISM_DG_EXE "MOD_FEM_conv_diff_prism_dg")
set(PDD_CONV_DIFF_HYBRID_STD_EXE "MOD_FEM_conv_diff_hybrid_std")
set(PDD_CONV_DIFF_HYBRID_STD_QUAD_EXE "MOD_FEM_conv_diff_hybrid_std_quad")
set(PDD_CONV_DIFF_HYBRID_DG_EXE "MOD_FEM_conv_diff_hybrid_dg")
set(PDD_CONV_DIFF_PRISM2D_STD_EXE "MOD_FEM_conv_diff_prism2d_std")
set(PDD_CONV_DIFF_PRISM2D_STD_QUAD_EXE "MOD_FEM_conv_diff_prism2d_std_quad")
set(PDD_CONV_DIFF_PRISM2D_DG_EXE "MOD_FEM_conv_diff_prism2d_dg")

# --- OpenCL support ---
if(${accel} STREQUAL "opencl")
  list(APPEND PDD_CONV_DIFF_STD_INCLUDE_DIRS
    ${pdd_conv_diff_srcdir}/linear_solver_interface/opencl/std_lin
    )
  list(APPEND PDD_CONV_DIFF_STD_QUAD_INCLUDE_DIRS
    ${pdd_conv_diff_srcdir}/linear_solver_interface/opencl/std_quad
    )
endif()

# -- Include directories for std conv diff --
list(APPEND PDD_CONV_DIFF_STD_INCLUDE_DIRS
  ${pdd_conv_diff_srcdir}/main
  ${pdd_conv_diff_srcdir}/main/util_std
  ${pdd_conv_diff_srcdir}/time_integration
  ${pdd_conv_diff_srcdir}/linear_solver_interface
  ${pdd_conv_diff_srcdir}/adaptation
  ${pdd_conv_diff_srcdir}/input_output
  ${pdd_conv_diff_srcdir}/weak_formulation
  ${pdd_conv_diff_srcdir}/weak_formulation/approx_std
  ${pdd_conv_diff_srcdir}/include
)
# -- Include directories for std quad conv diff --
list(APPEND PDD_CONV_DIFF_STD_QUAD_INCLUDE_DIRS
  ${pdd_conv_diff_srcdir}/main
  ${pdd_conv_diff_srcdir}/main/util_std
  ${pdd_conv_diff_srcdir}/time_integration
  ${pdd_conv_diff_srcdir}/linear_solver_interface
  ${pdd_conv_diff_srcdir}/adaptation
  ${pdd_conv_diff_srcdir}/input_output
  ${pdd_conv_diff_srcdir}/weak_formulation
  ${pdd_conv_diff_srcdir}/weak_formulation/approx_std
  ${pdd_conv_diff_srcdir}/include
)
# -- Include directories for dg conv_diff --
list(APPEND PDD_CONV_DIFF_DG_INCLUDE_DIRS
  ${pdd_conv_diff_srcdir}/main
  ${pdd_conv_diff_srcdir}/main/util_dg
  ${pdd_conv_diff_srcdir}/time_integration
  ${pdd_conv_diff_srcdir}/linear_solver_interface
  ${pdd_conv_diff_srcdir}/adaptation
  ${pdd_conv_diff_srcdir}/input_output
  ${pdd_conv_diff_srcdir}/weak_formulation
  ${pdd_conv_diff_srcdir}/weak_formulation/approx_dg
  ${pdd_conv_diff_srcdir}/include
)

if(PDD_CONV_DIFF_STD_INCLUDE_DIRS)
  set(PDD_CONV_DIFF_STD_INCLUDE_DIRS ${PDD_CONV_DIFF_STD_INCLUDE_DIRS} PARENT_SCOPE)
endif()

if(PDD_CONV_DIFF_STD_QUAD_INCLUDE_DIRS)
  set(PDD_CONV_DIFF_STD_QUAD_INCLUDE_DIRS ${PDD_CONV_DIFF_STD_QUAD_INCLUDE_DIRS} PARENT_SCOPE)
endif()

if(PDD_CONV_DIFF_DG_INCLUDE_DIRS)
  set(PDD_CONV_DIFF_DG_INCLUDE_DIRS ${PDD_CONV_DIFF_DG_INCLUDE_DIRS} PARENT_SCOPE)
endif()

    enable_testing ()
    # message(${CMAKE_BUILD_TYPE})

    foreach(folder ${MKL_INCLUDE_DIRS})
        set(KKK "${KKK}${folder}:")
    endforeach()

    foreach(folder ${MKL_LIBRARY_DIRS})
        set(KKK "${KKK}${folder}:")
    endforeach()

    set(KKK "${KKK}.")
    #message(${KKK})
    set(CTEST_BUILD_DIR $ENV{PWD})
    set(CTEST_DIR ${CMAKE_SOURCE_DIR}/../ctest)

# Create MOD_FEM_CONV_DIFF_PRISM_STD with OpenMP/OpenCL support
if(CREATE_MOD_FEM_CONV_DIFF_PRISM_STD)
  # -- Source files list --
  collect_src_from_include_dirs(PDD_CONV_DIFF_STD_INCLUDE_DIRS)
  # -- Target --
  include_directories(${PDD_CONV_DIFF_STD_INCLUDE_DIRS} ${PDD_CONV_DIFF_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

  add_executable(${PDD_CONV_DIFF_PRISM_STD_EXE} ${SOURCE_FILES})
  if(${MODFEM_ACCEL} STREQUAL "opencl")
    target_link_libraries(${PDD_CONV_DIFF_PRISM_STD_EXE} ${PDD_CONV_DIFF_EXTERNAL_LIBRARIES} mml_prism apl_std_lin utl_system utl_util tml_ocl ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD} ${OPENCL_LIBRARIES})
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM_STD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
  else()
    target_link_libraries(${PDD_CONV_DIFF_PRISM_STD_EXE} ${PDD_CONV_DIFF_EXTERNAL_LIBRARIES} mml_prism apl_std_lin utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM_STD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
  endif()
  SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM_STD_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
  SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM_STD_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")

  if(${MODFEM_ACCEL} STREQUAL "opencl")
    #add_mf_tests_mpi(2 "conv_diff" "prism" "std" "Laplace_std_in_box_opencl" "Run;L2;H1")
    
    add_mf_tests_mpi(1 "conv_diff" "prism" "std" "Laplace_std_in_box_mpi_opencl" "Run;L2;H1" "conf.txt_s211")
    add_mf_tests_mpi(1 "conv_diff" "prism" "std" "Laplace_std_in_box_mpi_opencl" "Run;L2;H1" "conf.txt_s214")
    add_mf_tests_mpi(1 "conv_diff" "prism" "std" "Laplace_std_in_box_mpi_opencl" "Run;L2;H1" "conf.txt_s411")
    
    add_mf_tests_mpi(2 "conv_diff" "prism" "std" "Laplace_std_in_box_mpi_opencl" "Run;L2;H1" "conf.txt_s211")
    add_mf_tests_mpi(2 "conv_diff" "prism" "std" "Laplace_std_in_box_mpi_opencl" "Run;L2;H1" "conf.txt_s214")
    add_mf_tests_mpi(2 "conv_diff" "prism" "std" "Laplace_std_in_box_mpi_opencl" "Run;L2;H1" "conf.txt_s411")
    
    add_mf_tests_mpi(4 "conv_diff" "prism" "std" "Laplace_std_in_box_mpi_opencl" "Run;L2;H1" "conf.txt_s211")
    add_mf_tests_mpi(4 "conv_diff" "prism" "std" "Laplace_std_in_box_mpi_opencl" "Run;L2;H1" "conf.txt_s214")
    add_mf_tests_mpi(4 "conv_diff" "prism" "std" "Laplace_std_in_box_mpi_opencl" "Run;L2;H1" "conf.txt_s411")
    
  else()
    add_mf_tests("conv_diff" "prism" "std" "Laplace_std_in_cube" "Run;L2;H1;ZZ")
    add_mf_tests("conv_diff" "prism" "std" "Laplace_std" "Run;L2;H1;ZZ")

    add_mf_tests_mpi(2 "conv_diff" "prism" "std" "Laplace_std_in_cube" "Run;L2;H1")
    add_mf_tests_mpi(2 "conv_diff" "prism" "std" "Laplace_std" "Run;L2;H1")
  endif()
    
endif()

unset(SOURCE_FILES)

# Create MOD_FEM_CONV_DIFF_PRISM_STD_QUAD with OpenMP/OpenCL support
if(CREATE_MOD_FEM_CONV_DIFF_PRISM_STD_QUAD)
  # -- Source files list --
  collect_src_from_include_dirs(PDD_CONV_DIFF_STD_QUAD_INCLUDE_DIRS)
  # -- Target --
  include_directories(${PDD_CONV_DIFF_STD_QUAD_INCLUDE_DIRS} ${PDD_CONV_DIFF_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

  add_executable(${PDD_CONV_DIFF_PRISM_STD_QUAD_EXE} ${SOURCE_FILES})
  if(${MODFEM_ACCEL} STREQUAL "opencl")
    target_link_libraries(${PDD_CONV_DIFF_PRISM_STD_QUAD_EXE} ${PDD_CONV_DIFF_EXTERNAL_LIBRARIES} mml_prism apl_std_quad utl_system utl_util tml_ocl ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD} ${OPENCL_LIBRARIES})
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM_STD_QUAD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
  else()
    target_link_libraries(${PDD_CONV_DIFF_PRISM_STD_QUAD_EXE} ${PDD_CONV_DIFF_EXTERNAL_LIBRARIES} mml_prism apl_std_quad utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM_STD_QUAD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
  endif()
  SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM_STD_QUAD_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
  SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM_STD_QUAD_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
endif()

unset(SOURCE_FILES)

# Create MOD_FEM_CONV_DIFF_PRISM_DG with OpenMP support
if(CREATE_MOD_FEM_CONV_DIFF_PRISM_DG)
  if(NOT ${MODFEM_ACCEL} STREQUAL "opencl")
    # -- Source files list --
    collect_src_from_include_dirs(PDD_CONV_DIFF_DG_INCLUDE_DIRS)
    # -- Target --
    include_directories(${PDD_CONV_DIFF_DG_INCLUDE_DIRS} ${PDD_CONV_DIFF_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

    add_executable(${PDD_CONV_DIFF_PRISM_DG_EXE} ${SOURCE_FILES})
    target_link_libraries(${PDD_CONV_DIFF_PRISM_DG_EXE} ${PDD_CONV_DIFF_EXTERNAL_LIBRARIES} mml_prism apl_dg_prism utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_DG})
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM_DG_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE}  RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM_DG_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM_DG_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
  else()
    message(WARNING "OpenCL for PDD_CONV_DIFF_PRISM_DG not implemented yet!!!")
  endif()

  #solver tests
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s110")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s111")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s112")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s113")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s114")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s115")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s116")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s117")

  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s120")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s121")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s122")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s123")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s124")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s125")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s126")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s127")
  
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s220")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s221")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s222")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s223")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s224")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s225")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s226")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s227")
  
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s210")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s211")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s212")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s213")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s214")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s215")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s216")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s217")
  
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s411")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_s421")
  
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_r500")

  #PDEG tests:
  
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p101")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p102")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p103")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p104")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p105")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p201")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p202")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p203")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p204")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p205")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p301")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p302")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p303")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p304")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p305")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p401")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p402")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p403")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p404")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p405")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p501")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p502")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p503")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p504")
  add_mf_tests("conv_diff" "prism" "dg" "Laplace_dg_in_cube" "Run;L2;H1" "conf.txt_p505")
  
endif()

unset(SOURCE_FILES)

# Create MOD_FEM_CONV_DIFF_HYBRID_STD with OpenMP/OpenCL support
if(CREATE_MOD_FEM_CONV_DIFF_HYBRID_STD)
    # -- Source files list --
    collect_src_from_include_dirs(PDD_CONV_DIFF_STD_INCLUDE_DIRS)
    # -- Target --
    include_directories(${PDD_CONV_DIFF_STD_INCLUDE_DIRS} ${PDD_CONV_DIFF_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

    add_executable(${PDD_CONV_DIFF_HYBRID_STD_EXE} ${SOURCE_FILES})
    if(${MODFEM_ACCEL} STREQUAL "opencl")
      target_link_libraries(${PDD_CONV_DIFF_HYBRID_STD_EXE} ${PDD_CONV_DIFF_EXTERNAL_LIBRARIES} mml_t4_prism apl_std_lin utl_system utl_util tml_ocl ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD} ${OPENCL_LIBRARIES})
      SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_HYBRID_STD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
    else()
      target_link_libraries(${PDD_CONV_DIFF_HYBRID_STD_EXE} ${PDD_CONV_DIFF_EXTERNAL_LIBRARIES} mml_t4_prism apl_std_lin utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
      SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_HYBRID_STD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
    endif()
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_HYBRID_STD_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_HYBRID_STD_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")


    add_mf_tests("conv_diff" "hybrid" "std" "Laplace_std_in_cube_tetra" "Run;L2;H1;ZZ")
    add_mf_tests("conv_diff" "hybrid" "std" "Laplace_std_in_cube" "Run;L2;H1;ZZ")
    add_mf_tests("conv_diff" "hybrid" "std" "Laplace_std" "Run;L2;H1;ZZ")
    add_mf_tests("conv_diff" "hybrid" "std" "Laplace_std_tetra" "Run;L2;H1;ZZ")

    add_mf_tests_mpi(2 "conv_diff" "hybrid" "std" "Laplace_std_in_cube" "Run;L2;H1")
    add_mf_tests_mpi(2 "conv_diff" "hybrid" "std" "Laplace_std_in_cube_tetra" "Run;L2;H1")
    # add_mf_tests_mpi(6 "conv_diff" "hybrid" "std" "Laplace_std_in_cube_tetra" "Run;L2;H1")

    add_mf_tests_mpi(2 "conv_diff" "hybrid" "std" "Laplace_std" "Run;L2;H1")
    add_mf_tests_mpi(2 "conv_diff" "hybrid" "std" "Laplace_std_tetra" "Run;L2;H1")
    add_mf_tests_mpi(3 "conv_diff" "hybrid" "std" "Laplace_std" "Run;L2;H1")
    add_mf_tests_mpi(3 "conv_diff" "hybrid" "std" "Laplace_std_tetra" "Run;L2;H1")
    add_mf_tests_mpi(4 "conv_diff" "hybrid" "std" "Laplace_std" "Run;L2;H1")
    add_mf_tests_mpi(4 "conv_diff" "hybrid" "std" "Laplace_std_tetra" "Run;L2;H1")
    add_mf_tests_mpi(5 "conv_diff" "hybrid" "std" "Laplace_std" "Run;L2;H1")
    add_mf_tests_mpi(5 "conv_diff" "hybrid" "std" "Laplace_std_tetra" "Run;L2;H1")
    add_mf_tests_mpi(6 "conv_diff" "hybrid" "std" "Laplace_std" "Run;L2;H1")
    add_mf_tests_mpi(6 "conv_diff" "hybrid" "std" "Laplace_std_tetra" "Run;L2;H1")
    add_mf_tests_mpi(7 "conv_diff" "hybrid" "std" "Laplace_std" "Run;L2;H1")
    add_mf_tests_mpi(7 "conv_diff" "hybrid" "std" "Laplace_std_tetra" "Run;L2;H1")
    add_mf_tests_mpi(8 "conv_diff" "hybrid" "std" "Laplace_std" "Run;L2;H1")
    add_mf_tests_mpi(8 "conv_diff" "hybrid" "std" "Laplace_std_tetra" "Run;L2;H1")
endif()

unset(SOURCE_FILES)

# Create MOD_FEM_CONV_DIFF_HYBRID_STD_QUAD with OpenMP/OpenCL support
if(CREATE_MOD_FEM_CONV_DIFF_HYBRID_STD_QUAD)
  if(NOT ${MODFEM_ACCEL} STREQUAL "opencl")
    # -- Source files list --
    collect_src_from_include_dirs(PDD_CONV_DIFF_STD_QUAD_INCLUDE_DIRS)
    # -- Target --
    include_directories(${PDD_CONV_DIFF_STD_QUAD_INCLUDE_DIRS} ${PDD_CONV_DIFF_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

    add_executable(${PDD_CONV_DIFF_HYBRID_STD_QUAD_EXE} ${SOURCE_FILES})
    if(${MODFEM_ACCEL} STREQUAL "opencl")
      target_link_libraries(${PDD_CONV_DIFF_HYBRID_STD_QUAD_EXE} ${PDD_CONV_DIFF_EXTERNAL_LIBRARIES} mml_t4_prism apl_std_quad utl_system utl_util tml_ocl ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD} ${OPENCL_LIBRARIES})
      SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_HYBRID_STD_QUAD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
    else()
      target_link_libraries(${PDD_CONV_DIFF_HYBRID_STD_QUAD_EXE} ${PDD_CONV_DIFF_EXTERNAL_LIBRARIES} mml_t4_prism apl_std_quad utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
      SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_HYBRID_STD_QUAD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
    endif()
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_HYBRID_STD_QUAD_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_HYBRID_STD_QUAD_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
  else()
    message(WARNING "OpenCL for PDD_CONV_DIFF_HYBRID_QUAD_STD not implemented yet!!!")
  endif()
endif()

unset(SOURCE_FILES)

# Create MOD_FEM_CONV_DIFF_HYBRID_DG with OpenMP support
if(CREATE_MOD_FEM_CONV_DIFF_HYBRID_DG)
  if(NOT ${MODFEM_ACCEL} STREQUAL "opencl")
    # -- Source files list --
    collect_src_from_include_dirs(PDD_CONV_DIFF_DG_INCLUDE_DIRS)
    # -- Target --
    include_directories(${PDD_CONV_DIFF_DG_INCLUDE_DIRS} ${PDD_CONV_DIFF_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

    add_executable(${PDD_CONV_DIFF_HYBRID_DG_EXE} ${SOURCE_FILES})
    target_link_libraries(${PDD_CONV_DIFF_HYBRID_DG_EXE} ${PDD_CONV_DIFF_EXTERNAL_LIBRARIES} mml_t4_prism apl_dg_prism utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_DG})
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_HYBRID_DG_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_HYBRID_DG_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_HYBRID_DG_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
  else()
    message(WARNING "OpenCL for CREATE_MOD_FEM_CONV_DIFF_HYBRID_DG not implemented yet!!!")
  endif()
endif()

unset(SOURCE_FILES)

# Create MOD_FEM_CONV_DIFF_PRISM2D_STD with OpenMP/OpenCL support
if(CREATE_MOD_FEM_CONV_DIFF_PRISM2D_STD)
  # -- Source files list --
  collect_src_from_include_dirs(PDD_CONV_DIFF_STD_INCLUDE_DIRS)
  # -- Target --
  include_directories(${PDD_CONV_DIFF_STD_INCLUDE_DIRS} ${PDD_CONV_DIFF_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

  add_executable(${PDD_CONV_DIFF_PRISM2D_STD_EXE} ${SOURCE_FILES})
  if(${MODFEM_ACCEL} STREQUAL "opencl")
    target_link_libraries(${PDD_CONV_DIFF_PRISM2D_STD_EXE} ${PDD_CONV_DIFF_EXTERNAL_LIBRARIES} mml_prism_2D apl_std_lin utl_system utl_util tml_ocl ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD} ${OPENCL_LIBRARIES})
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM2D_STD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
  else()
    target_link_libraries(${PDD_CONV_DIFF_PRISM2D_STD_EXE} ${PDD_CONV_DIFF_EXTERNAL_LIBRARIES} mml_prism_2D apl_std_lin utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM2D_STD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
  endif()
  SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM2D_STD_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
  SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM2D_STD_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
  
    #solver tests
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s11-4")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s11-2")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s110")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s111")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s112")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s113")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s114")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s115")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s116")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s117")

  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s12-4")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s12-2")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s120")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s121")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s122")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s123")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s124")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s125")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s126")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s127")
  
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s22-4")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s22-2")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s220")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s221")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s222")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s223")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s224")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s225")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s226")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s227")
  
  
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s21-4")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s21-2")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s210")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s211")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s212")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s213")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s214")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s215")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s216")
  add_mf_tests("conv_diff" "prism" "std" "Laplace_std2d_in_cube" "Run;L2;H1" "conf.txt_s217")
  
  
endif()

unset(SOURCE_FILES)

# Create MOD_FEM_CONV_DIFF_PRISM2D_STD_QUAD with OpenMP/OpenCL support
if(CREATE_MOD_FEM_CONV_DIFF_PRISM2D_STD_QUAD)
  # -- Source files list --
  collect_src_from_include_dirs(PDD_CONV_DIFF_STD_QUAD_INCLUDE_DIRS)
  # -- Target --
  include_directories(${PDD_CONV_DIFF_STD_QUAD_INCLUDE_DIRS} ${PDD_CONV_DIFF_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

  add_executable(${PDD_CONV_DIFF_PRISM2D_STD_QUAD_EXE} ${SOURCE_FILES})
  if(${MODFEM_ACCEL} STREQUAL "opencl")
    target_link_libraries(${PDD_CONV_DIFF_PRISM2D_STD_QUAD_EXE} ${PDD_CONV_DIFF_EXTERNAL_LIBRARIES} mml_prism_2D apl_std_quad utl_system utl_util tml_ocl ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD} ${OPENCL_LIBRARIES})
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM2D_STD_QUAD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
  else()
    target_link_libraries(${PDD_CONV_DIFF_PRISM2D_STD_QUAD_EXE} ${PDD_CONV_DIFF_EXTERNAL_LIBRARIES} mml_prism_2D apl_std_quad utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM2D_STD_QUAD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
  endif()
  SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM2D_STD_QUAD_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
  SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM2D_STD_QUAD_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
endif()

unset(SOURCE_FILES)

# Create MOD_FEM_CONV_DIFF_PRISM2D_DG with OpenMP support
if(CREATE_MOD_FEM_CONV_DIFF_PRISM2D_DG)
  if(NOT ${MODFEM_ACCEL} STREQUAL "opencl")
    # -- Source files list --
    collect_src_from_include_dirs(PDD_CONV_DIFF_DG_INCLUDE_DIRS)
    # -- Target --
    include_directories(${PDD_CONV_DIFF_DG_INCLUDE_DIRS} ${PDD_CONV_DIFF_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

    add_executable(${PDD_CONV_DIFF_PRISM2D_DG_EXE} ${SOURCE_FILES})
    target_link_libraries(${PDD_CONV_DIFF_PRISM2D_DG_EXE} ${PDD_CONV_DIFF_EXTERNAL_LIBRARIES} mml_prism_2D apl_dg_prism utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_DG})
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM2D_DG_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM2D_DG_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${PDD_CONV_DIFF_PRISM2D_DG_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
  else()
    message(WARNING "OpenCL for MOD_FEM_CONV_DIFF_PRISM2D_DG not implemented yet!!!")
  endif()
endif()

