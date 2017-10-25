# Uses:
#
# PDD_HEAT_INCLUDE_DIRS
# PDD_HEAT_EXTERNAL_LIBRARIES
#
# CREATE_PDD_HEAT_PRISM_STD
# CREATE_PDD_HEAT_PRISM_STD_QUAD
# CREATE_PDD_HEAT_HYBRID_STD
# CREATE_PDD_HEAT_HYBRID_STD_QUAD
#

# -- Add target to global list --
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_HEAT_PRISM_STD")
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_HEAT_PRISM_STD_QUAD")
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_HEAT_PRISM2D_STD")
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_HEAT_PRISM2D_STD_QUAD")
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_HEAT_HYBRID_STD")
list(APPEND MOD_FEM_EXEC_LIST "CREATE_MOD_FEM_HEAT_HYBRID_STD_QUAD")
set(MOD_FEM_EXEC_LIST "${MOD_FEM_EXEC_LIST}" PARENT_SCOPE)

# -- Exec files --
set(pdd_heat_srcdir "${CMAKE_SOURCE_DIR}/pdd_heat")
set(PDD_HEAT_PRISM_STD_EXE "MOD_FEM_heat_prism_std")
set(PDD_HEAT_PRISM_STD_QUAD_EXE "MOD_FEM_heat_prism_std_quad")
set(PDD_HEAT_PRISM2D_STD_EXE "MOD_FEM_heat_prism2d_std")
set(PDD_HEAT_PRISM2D_STD_QUAD_EXE "MOD_FEM_heat_prism2d_std_quad")
set(PDD_HEAT_HYBRID_STD_EXE "MOD_FEM_heat_hybrid_std")
set(PDD_HEAT_HYBRID_STD_QUAD_EXE "MOD_FEM_heat_hybrid_std_quad")

# --- OpenCL support ---
if(${accel} STREQUAL "opencl")
  list(APPEND PDD_HEAT_STD_INCLUDE_DIRS
    ${pdd_heat_srcdir}/linear_solver_interface/opencl/std_lin
    )
  list(APPEND PDD_HEAT_STD_QUAD_INCLUDE_DIRS
    ${pdd_heat_srcdir}/linear_solver_interface/opencl/std_quad
    )
endif()


# -- Include directories for heat library --
list(APPEND PDD_HEAT_INCLUDE_DIRS
  ${pdd_heat_srcdir}/input_output
  ${pdd_heat_srcdir}/materials
  ${pdd_heat_srcdir}/weak_formulation
  ${pdd_heat_srcdir}/include
  ${pdd_heat_srcdir}/
)

# -- Include directories for std heat --
list(APPEND PDD_HEAT_STD_INCLUDE_DIRS
  ${pdd_heat_srcdir}/main
  ${pdd_heat_srcdir}/time_integration
  ${pdd_heat_srcdir}/linear_solver_interface
  ${pdd_heat_srcdir}/adaptation
  ${pdd_heat_srcdir}/input_output
  ${pdd_heat_srcdir}/materials
  ${pdd_heat_srcdir}/weak_formulation
  ${pdd_heat_srcdir}/include
  )

# -- Include directories for std quad heat --
list(APPEND PDD_HEAT_STD_QUAD_INCLUDE_DIRS
  ${pdd_heat_srcdir}/main
  ${pdd_heat_srcdir}/time_integration
  ${pdd_heat_srcdir}/linear_solver_interface
  ${pdd_heat_srcdir}/adaptation
  ${pdd_heat_srcdir}/input_output
  ${pdd_heat_srcdir}/materials
  ${pdd_heat_srcdir}/weak_formulation
  ${pdd_heat_srcdir}/include
  )

if(PDD_HEAT_INCLUDE_DIRS)
  set(PDD_HEAT_INCLUDE_DIRS ${PDD_HEAT_INCLUDE_DIRS} PARENT_SCOPE)
endif()

if(PDD_HEAT_STD_INCLUDE_DIRS)
  set(PDD_HEAT_STD_INCLUDE_DIRS ${PDD_HEAT_STD_INCLUDE_DIRS} PARENT_SCOPE)
endif()

if(PDD_HEAT_STD_QUAD_INCLUDE_DIRS)
  set(PDD_HEAT_STD_QUAD_INCLUDE_DIRS ${PDD_HEAT_STD_QUAD_INCLUDE_DIRS} PARENT_SCOPE)
endif()

enable_testing ()

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


# ----------------------- Create HEAT library only ------------------------------- #

# -- Source files list --
collect_src_from_include_dirs(PDD_HEAT_INCLUDE_DIRS)

# -- Target --
include_directories(${PDD_HEAT_INCLUDE_DIRS} ${PDD_HEAT_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

#if(CREATE_PDD_HEAT_PRISM)
  add_library(pdl_heat ${SOURCE_FILES})
  target_link_libraries(pdl_heat ${PDD_HEAT_EXTERNAL_LIBRARIES})
  SET_TARGET_PROPERTIES(pdl_heat PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
  SET_TARGET_PROPERTIES(pdl_heat PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
  SET_TARGET_PROPERTIES(pdl_heat PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
#endif()
unset(SOURCE_FILES)

# -------- Create MOD_FEM_HEAT_PRISM_STD with OpenMP/OpenCL support --------- #
if(CREATE_MOD_FEM_HEAT_PRISM_STD)
  # -- Source files list --
  collect_src_from_include_dirs(PDD_HEAT_STD_INCLUDE_DIRS)
  # -- Target --
  include_directories(${PDD_HEAT_STD_INCLUDE_DIRS} ${PDD_HEAT_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

  add_executable(${PDD_HEAT_PRISM_STD_EXE} ${SOURCE_FILES})
  if(${MODFEM_ACCEL} STREQUAL "opencl")
    target_link_libraries(${PDD_HEAT_PRISM_STD_EXE} ${PDD_HEAT_EXTERNAL_LIBRARIES} mml_prism apl_std_lin utl_system utl_util tml_ocl ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD} ${OPENCL_LIBRARIES})
    SET_TARGET_PROPERTIES(${PDD_HEAT_PRISM_STD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
  else()
    target_link_libraries(${PDD_HEAT_PRISM_STD_EXE} ${PDD_HEAT_EXTERNAL_LIBRARIES} mml_prism apl_std_lin utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
    SET_TARGET_PROPERTIES(${PDD_HEAT_PRISM_STD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
  endif()
  SET_TARGET_PROPERTIES(${PDD_HEAT_PRISM_STD_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
  SET_TARGET_PROPERTIES(${PDD_HEAT_PRISM_STD_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")

  add_mf_tests("heat" "prism" "std" "cubeP_conv_prism" "Run;ZZ")

endif()
unset(SOURCE_FILES)

# ----- Create MOD_FEM_HEAT_PRISM_STD_QUAD with OpenMP/OpenCL support ------- #
if(CREATE_MOD_FEM_HEAT_PRISM_STD_QUAD)
  # -- Source files list --
  collect_src_from_include_dirs(PDD_HEAT_STD_QUAD_INCLUDE_DIRS)
  # -- Target --
  include_directories(${PDD_HEAT_STD_QUAD_INCLUDE_DIRS} ${PDD_HEAT_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

  add_executable(${PDD_HEAT_PRISM_STD_QUAD_EXE} ${SOURCE_FILES})
  if(${MODFEM_ACCEL} STREQUAL "opencl")
    target_link_libraries(${PDD_HEAT_PRISM_STD_QUAD_EXE} ${PDD_HEAT_EXTERNAL_LIBRARIES} mml_prism apl_std_quad utl_system utl_util tml_ocl ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD} ${OPENCL_LIBRARIES})
    SET_TARGET_PROPERTIES(${PDD_HEAT_PRISM_STD_QUAD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
  else()
    target_link_libraries(${PDD_HEAT_PRISM_STD_QUAD_EXE} ${PDD_HEAT_EXTERNAL_LIBRARIES} mml_prism apl_std_quad utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
    SET_TARGET_PROPERTIES(${PDD_HEAT_PRISM_STD_QUAD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
  endif()
  SET_TARGET_PROPERTIES(${PDD_HEAT_PRISM_STD_QUAD_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
  SET_TARGET_PROPERTIES(${PDD_HEAT_PRISM_STD_QUAD_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
endif()
unset(SOURCE_FILES)

# ------- Create MOD_FEM_HEAT_PRISM2D_STD with OpenMP/OpenCL support -------- #
if(CREATE_MOD_FEM_HEAT_PRISM2D_STD)
  if(NOT ${MODFEM_ACCEL} STREQUAL "opencl")
    # -- Source files list --
    collect_src_from_include_dirs(PDD_HEAT_STD_INCLUDE_DIRS)
    # -- Target --
    include_directories(${PDD_HEAT_STD_INCLUDE_DIRS} ${PDD_HEAT_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

    add_executable(${PDD_HEAT_PRISM2D_STD_EXE} ${SOURCE_FILES})
    if(${MODFEM_ACCEL} STREQUAL "opencl")
      #TODO
    else()
      target_link_libraries(${PDD_HEAT_PRISM2D_STD_EXE} ${PDD_HEAT_EXTERNAL_LIBRARIES} mml_prism_2D apl_std_lin utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
      SET_TARGET_PROPERTIES(${PDD_HEAT_PRISM2D_STD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
    endif()
    SET_TARGET_PROPERTIES(${PDD_HEAT_PRISM2D_STD_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${PDD_HEAT_PRISM2D_STD_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
  else()
    message(WARNING "OpenCL for MOD_FEM_HEAT_PRISM2D_STD not implemented yet!!!")
  endif()
endif()
unset(SOURCE_FILES)

# ---- Create MOD_FEM_HEAT_PRISM2D_STD_QUAD with OpenMP/OpenCL support ------ #
if(CREATE_MOD_FEM_HEAT_PRISM2D_STD_QUAD)
  if(NOT ${MODFEM_ACCEL} STREQUAL "opencl")
    # -- Source files list --
    collect_src_from_include_dirs(PDD_HEAT_STD_QUAD_INCLUDE_DIRS)
    # -- Target --
    include_directories(${PDD_HEAT_STD_QUAD_INCLUDE_DIRS} ${PDD_HEAT_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

    add_executable(${PDD_HEAT_PRISM2D_STD_QUAD_EXE} ${SOURCE_FILES})
    if(${MODFEM_ACCEL} STREQUAL "opencl")
      #TODO
    else()
      target_link_libraries(${PDD_HEAT_PRISM2D_STD_QUAD_EXE} ${PDD_HEAT_EXTERNAL_LIBRARIES} mml_prism_2D apl_std_quad utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
      SET_TARGET_PROPERTIES(${PDD_HEAT_PRISM2D_STD_QUAD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
    endif()
    SET_TARGET_PROPERTIES(${PDD_HEAT_PRISM2D_STD_QUAD_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${PDD_HEAT_PRISM2D_STD_QUAD_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
  else()
    message(WARNING "OpenCL for PDD_HEAT_PRISM_STD2D_QUAD not implemented yet!!!")
  endif()
endif()
unset(SOURCE_FILES)

# -------- Create MOD_FEM_HEAT_HYBRID_STD with OpenMP/OpenCL support --------- #
if(CREATE_MOD_FEM_HEAT_HYBRID_STD)
  # -- Source files list --
  collect_src_from_include_dirs(PDD_HEAT_STD_INCLUDE_DIRS)
  # -- Target --
  include_directories(${PDD_HEAT_STD_INCLUDE_DIRS} ${PDD_HEAT_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
    
  add_executable(${PDD_HEAT_HYBRID_STD_EXE} ${SOURCE_FILES})
  if(${MODFEM_ACCEL} STREQUAL "opencl")
    target_link_libraries(${PDD_HEAT_HYBRID_STD_EXE} ${PDD_HEAT_EXTERNAL_LIBRARIES} mml_t4_prism apl_std_lin utl_system utl_util tml_ocl ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD} ${OPENCL_LIBRARIES})
    SET_TARGET_PROPERTIES(${PDD_HEAT_HYBRID_STD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
  else()
    target_link_libraries(${PDD_HEAT_HYBRID_STD_EXE} ${PDD_HEAT_EXTERNAL_LIBRARIES} mml_t4_prism apl_std_lin utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
    SET_TARGET_PROPERTIES(${PDD_HEAT_HYBRID_STD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
  endif()
  SET_TARGET_PROPERTIES(${PDD_HEAT_HYBRID_STD_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
  SET_TARGET_PROPERTIES(${PDD_HEAT_HYBRID_STD_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")

  #add_mf_tests("heat" "hybrid" "std" "cubeP_conv_prism" "Run")

endif()
unset(SOURCE_FILES)


# ---- Create MOD_FEM_HEAT_HYBRID_STD_QUAD with OpenMP/OpenCL support --------- #
if(CREATE_MOD_FEM_HEAT_HYBRID_STD_QUAD)
  if(NOT ${MODFEM_ACCEL} STREQUAL "opencl")
    # -- Source files list --
    collect_src_from_include_dirs(PDD_HEAT_STD_QUAD_INCLUDE_DIRS)
    # -- Target --
    include_directories(${PDD_HEAT_STD_QUAD_INCLUDE_DIRS} ${PDD_HEAT_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

    add_executable(${PDD_HEAT_HYBRID_STD_QUAD_EXE} ${SOURCE_FILES})
    if(${MODFEM_ACCEL} STREQUAL "opencl")
      #TODO
    else()
      target_link_libraries(${PDD_HEAT_HYBRID_STD_QUAD_EXE} ${PDD_HEAT_EXTERNAL_LIBRARIES} mml_t4_prism apl_std_quad utl_system utl_util ${MODFEM_DIRECT_SOLVER_MODULE} ${MODFEM_ITER_SOLVER_MODULE} ${MODFEM_MPI_LIBS_STD})
      SET_TARGET_PROPERTIES(${PDD_HEAT_HYBRID_STD_QUAD_EXE} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE} RELEASE_POSTFIX "${RELEASE_POSTFIX}" DEBUG_POSTFIX "${DEBUG_POSTFIX}")
    endif()
    SET_TARGET_PROPERTIES(${PDD_HEAT_HYBRID_STD_QUAD_EXE} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
    SET_TARGET_PROPERTIES(${PDD_HEAT_HYBRID_STD_QUAD_EXE} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
  else()
    message(WARNING "OpenCL for PDD_HEAT_HYBRID_STD_QUAD not implemented yet!!!")
  endif()
endif()
unset(SOURCE_FILES)
