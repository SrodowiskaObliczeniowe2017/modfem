# -- Add target to global list --
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_NS_SUPG_PRISM_STD")
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_NS_SUPG_HYBRID_STD")
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_NS_SUPG_REMESH_STD")
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_NS_SUPG_PRISM2D_STD")
#list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_NS_SUPG_PRISM_STD_TURBULENT")
#list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_NS_SUPG_HYBRID_STD_TURBULENT")
set(MOD_FEM_EXEC_LIST "${MOD_FEM_EXEC_LIST}" PARENT_SCOPE)

# -- Exec files --
set(pdd_ns_supg_srcdir "${CMAKE_SOURCE_DIR}/pdd_ns_supg")
set(PDD_NS_SUPG_PRISM_STD_EXE "MOD_FEM_ns_supg_prism_std")
#set(PDD_NS_SUPG_PRISM_STD_TURBULENT_EXE "MOD_FEM_ns_supg_prism_std_turbulent")
set(PDD_NS_SUPG_HYBRID_STD_EXE "MOD_FEM_ns_supg_hybrid_std")
#set(PDD_NS_SUPG_HYBRID_STD_TURBULENT_EXE "MOD_FEM_ns_supg_hybrid_std_turbulent")
set(PDD_NS_SUPG_PRISM2D_STD_EXE "MOD_FEM_ns_supg_prism2d_std")
#set(PDD_NS_SUPG_REMESH_EXE "MOD_FEM_ns_supg_remesh")
set(PDD_NS_SUPG_REMESH_STD_EXE "MOD_FEM_ns_supg_remesh_std")
#set(PDD_NS_SUPG_FCPM_STD_EXE "MOD_FEM_ns_supg_fcp_std")

# -- Include directories --
list(APPEND PDD_NS_SUPG_INCLUDE_DIRS
  ${pdd_ns_supg_srcdir}/adaptation
  ${pdd_ns_supg_srcdir}/input_output
  ${pdd_ns_supg_srcdir}/materials
  ${pdd_ns_supg_srcdir}/linear_solver_interface
  ${pdd_ns_supg_srcdir}/main
  ${pdd_ns_supg_srcdir}/include
  ${pdd_ns_supg_srcdir}/time_integration
  ${pdd_ns_supg_srcdir}/weak_formulation
  ${pdd_ns_supg_srcdir}
)
if(PDD_NS_SUPG_INCLUDE_DIRS)
  set(PDD_NS_SUPG_INCLUDE_DIRS ${PDD_NS_SUPG_INCLUDE_DIRS} PARENT_SCOPE)
endif()
# -- Source files list --
collect_src_from_include_dirs(PDD_NS_SUPG_INCLUDE_DIRS)

# -- Target --
include_directories(${PDD_NS_SUPG_INCLUDE_DIRS} ${PDD_NS_SUPG_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

enable_testing ()

foreach(folder ${MKL_INCLUDE_DIRS})
    set(KKK "${KKK}${folder}:")
endforeach()

foreach(folder ${MKL_LIBRARY_DIRS})
    set(KKK "${KKK}${folder}:")
endforeach()

set(KKK "${KKK}.")

set(CTEST_BUILD_DIR $ENV{PWD})
set(CTEST_DIR ${CMAKE_SOURCE_DIR}/../ctest)

if(CREATE_MOD_FEM_NS_SUPG_PRISM_STD)
  if(NOT ${MODFEM_ACCEL} STREQUAL "opencl")
    add_executable(${PDD_NS_SUPG_PRISM_STD_EXE} ${SOURCE_FILES})
    #target_link_libraries(${PDD_NS_SUPG_PRISM_STD_EXE} ${PDD_NS_SUPG_EXTERNAL_LIBRARIES} mml_prism apl_std_lin utl_system utl_util tml_omp ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
    target_link_libraries(${PDD_NS_SUPG_PRISM_STD_EXE} ${PDD_NS_SUPG_EXTERNAL_LIBRARIES} mml_prism apl_std_lin utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_PRISM_STD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_PRISM_STD_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_PRISM_STD_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
  else()
    message(WARNING "OpenCL for PDD_NS_SUPG_PRISM_STD not implemented yet!!!")
  endif()

  add_mf_tests("ns_supg" "prism" "std" "LDC_jk_1x1" "Run;Diff;ZZ")
  add_mf_tests_mpi(2 "ns_supg" "prism" "std" "LDC_jk_1x1" "Run")
  add_mf_tests("ns_supg" "prism" "std" "bfs_jk_Re100" "Run;Diff")
endif()

if(CREATE_MOD_FEM_NS_SUPG_HYBRID_STD)
  if(NOT ${MODFEM_ACCEL} STREQUAL "opencl")
    add_executable(${PDD_NS_SUPG_HYBRID_STD_EXE} ${SOURCE_FILES})
    target_link_libraries(${PDD_NS_SUPG_HYBRID_STD_EXE} ${PDD_NS_SUPG_EXTERNAL_LIBRARIES} mml_t4_prism apl_std_lin utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD} )
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_HYBRID_STD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_HYBRID_STD_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_HYBRID_STD_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}" )
    #message("${PDD_NS_SUPG_HYBRID_STD_EXE} LINK_FLAGS ARE: ${LINK_FLAGS}")
    #message("${PDD_NS_SUPG_HYBRID_STD_EXE} COMPILE_DEFINITIONS ARE: ${MODFEM_MPI_CODE_DEFINE}")
    #message("${PDD_NS_SUPG_HYBRID_STD_EXE} COMPILE_FLAGS ARE: ${COMPILE_FLAGS}")
  else()
    message(WARNING "OpenCL for PDD_NS_SUPG_HYBRID_STD not implemented yet!!!")
  endif()

    add_mf_tests("ns_supg" "hybrid" "std" "LDC_nas_20x20" "Run;Diff")
    add_mf_tests("ns_supg" "hybrid" "std" "bfs_nas_Re100" "Run;Diff")
    add_mf_tests("ns_supg" "hybrid" "std" "VonKarman" "Run")
    #add_mf_tests("ns_supg" "hybrid" "std" "bfs_jk_Re100" "Run;Diff")

endif()

if(CREATE_MOD_FEM_NS_SUPG_REMESH_STD)
  if(NOT ${MODFEM_ACCEL} STREQUAL "opencl")
    add_executable(${PDD_NS_SUPG_REMESH_STD_EXE} ${SOURCE_FILES})
    target_link_libraries(${PDD_NS_SUPG_REMESH_STD_EXE} ${PDD_NS_SUPG_EXTERNAL_LIBRARIES} mml_remesh apl_std_lin utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_REMESH_STD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_REMESH_STD_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_REMESH_STD_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
  else()
    message(WARNING "OpenCL for PDD_NS_SUPG_REMESH_STD not implemented yet!!!")
  endif()
    add_mf_tests("ns_supg" "remesh" "std" "LDC_nas_20x20" "Run;Diff")
    add_mf_tests("ns_supg" "remesh" "std" "bfs_nas_Re100" "Run;Diff")
endif()


if(CREATE_MOD_FEM_NS_SUPG_PRISM2D_STD)
  if(NOT ${MODFEM_ACCEL} STREQUAL "opencl")
    add_executable(${PDD_NS_SUPG_PRISM2D_STD_EXE} ${SOURCE_FILES})
    target_link_libraries(${PDD_NS_SUPG_PRISM2D_STD_EXE} ${PDD_NS_SUPG_EXTERNAL_LIBRARIES} mml_prism_2D apl_std_lin utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_PRISM2D_STD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_PRISM2D_STD_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_PRISM2D_STD_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
  else()
    message(WARNING "OpenCL for PDD_NS_SUPG_PRISM2D_STD not implemented yet!!!")
  endif()
    add_mf_tests("ns_supg" "prism" "std" "LDC_jk_prism2d_std" "Run;ZZ" "conf.txt_s224")
    add_mf_tests("ns_supg" "prism" "std" "LDC_jk_prism2d_std" "Run;ZZ" "conf.txt_s411")
    add_mf_tests("ns_supg" "prism" "std" "LDC_jk_prism2d_std" "Run;ZZ" "conf.txt_s4011")
    add_mf_tests("ns_supg" "prism" "std" "LDC_jk_prism2d_std" "Run;ZZ" "conf.txt_s4111")
    add_mf_tests("ns_supg" "prism" "std" "LDC_jk_prism2d_std" "Run;ZZ" "conf.txt_s4211")
    add_mf_tests("ns_supg" "prism" "std" "LDC_jk_prism2d_std" "Run;ZZ" "conf.txt_s4311")
    add_mf_tests("ns_supg" "prism" "std" "LDC_jk_prism2d_std" "Run;ZZ" "conf.txt_s5011")
    add_mf_tests("ns_supg" "prism" "std" "LDC_jk_prism2d_std" "Run;ZZ" "conf.txt_s5111")
    add_mf_tests("ns_supg" "prism" "std" "LDC_jk_prism2d_std" "Run;ZZ" "conf.txt_s5211")
    add_mf_tests("ns_supg" "prism" "std" "LDC_jk_prism2d_std" "Run;ZZ" "conf.txt_s5311")
endif()

if(CREATE_MOD_FEM_NS_SUPG_PRISM_STD_TURBULENT)
  if(NOT ${MODFEM_ACCEL} STREQUAL "opencl")
    add_executable(${PDD_NS_SUPG_PRISM_STD_TURBULENT_EXE} ${SOURCE_FILES})
    target_link_libraries(${PDD_NS_SUPG_PRISM_STD_TURBULENT_EXE} ${PDD_NS_SUPG_PRISM_STD_TURBULENT_EXTERNAL_LIBRARIES} mml_prism apl_std_lin utl_system utl_util pdl_turbulent_lib pdl_cbs_lib ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_PRISM_STD_TURBULENT_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_PRISM_STD_TURBULENT_EXE} PROPERTIES COMPILE_DEFINITIONS TURBULENTFLOW)
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_PRISM_STD_TURBULENT_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_PRISM_STD_TURBULENT_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
  else()
    message(WARNING "OpenCL for PDD_NS_SUPG_PRISM_STD_TURBULENT not implemented yet!!!")
  endif()
endif()

if(CREATE_MOD_FEM_NS_SUPG_HYBRID_STD_TURBULENT)
  if(NOT ${MODFEM_ACCEL} STREQUAL "opencl")
    add_executable(${PDD_NS_SUPG_HYBRID_STD_TURBULENT_EXE} ${SOURCE_FILES})
    target_link_libraries(${PDD_NS_SUPG_HYBRID_STD_TURBULENT_EXE} ${PDD_NS_SUPG_HYBRID_STD_TURBULENT_EXTERNAL_LIBRARIES} mml_t4_prism apl_std_lin utl_system utl_util pdl_turbulent_lib pdl_cbs_lib ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_HYBRID_STD_TURBULENT_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_HYBRID_STD_TURBULENT_EXE} PROPERTIES COMPILE_DEFINITIONS TURBULENTFLOW)
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_HYBRID_STD_TURBULENT_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_HYBRID_STD_TURBULENT_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
  else()
    message(WARNING "OpenCL for PDD_NS_SUPG_HYBRID_STD_TURBULENT not implemented yet!!!")
  endif()
endif()

if(CREATE_MOD_FEM_NS_SUPG_FCPM_STD)
  if(NOT ${MODFEM_ACCEL} STREQUAL "opencl")
    add_executable(${PDD_NS_SUPG_FCPM_STD_EXE} ${SOURCE_FILES})
    target_link_libraries(${PDD_NS_SUPG_FCPM_STD_EXE} ${PDD_NS_SUPG_EXTERNAL_LIBRARIES} mml_fcpm apl_std_lin utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD} )
    target_link_libraries(${PDD_NS_SUPG_FCPM_STD_EXE} ${MODFEM_USER_LIBS})
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_FCPM_STD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_FCPM_STD_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${PDD_NS_SUPG_FCPM_STD_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}" )
    #message("${PDD_NS_SUPG_FCPM_STD_EXE} LINK_FLAGS ARE: ${LINK_FLAGS}")
    #message("${PDD_NS_SUPG_FCPM_STD_EXE} COMPILE_DEFINITIONS ARE: ${MODFEM_MPI_CODE_DEFINE}")
    #message("${PDD_NS_SUPG_FCPM_STD_EXE} COMPILE_FLAGS ARE: ${COMPILE_FLAGS}")
  else()
    message(WARNING "OpenCL for PDD_NS_SUPG_FCPM_STD not implemented yet!!!")
  endif()
endif()

##############################################################################################################

