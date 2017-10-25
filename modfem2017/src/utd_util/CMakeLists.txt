set(utd_util_srcdir "${CMAKE_SOURCE_DIR}/utd_util")

# -- Include directories --
list(APPEND UTD_UTIL_INCLUDE_DIRS
  ${utd_util_srcdir}
)

if(UTD_UTIL_INCLUDE_DIRS)
  set(UTD_UTIL_INCLUDE_DIRS ${UTD_UTIL_INCLUDE_DIRS} PARENT_SCOPE)
endif()

# -- Source files list --
#collect_src_from_include_dirs(UTD_UTIL_INCLUDE_DIRS)



# -- Target --
include_directories( ${UTD_UTIL_INCLUDE_DIRS} ${UTD_UTIL_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

add_library(utl_log "uts_log")
SET_TARGET_PROPERTIES(utl_log PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(utl_log PROPERTIES LINK_FLAGS "${LINK_FLAGS}")

add_library(utl_io_results "uts_io_results")
target_link_libraries(utl_io_results utl_log ${UTD_UTIL_EXTERNAL_LIBRARIES})
SET_TARGET_PROPERTIES(utl_io_results PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(utl_io_results PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(utl_io_results PROPERTIES LINK_FLAGS "${LINK_FLAGS}")

add_library(utl_bc "uts_bc")
target_link_libraries(utl_bc utl_log ${UTD_UTIL_EXTERNAL_LIBRARIES})
SET_TARGET_PROPERTIES(utl_bc PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(utl_bc PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(utl_bc PROPERTIES LINK_FLAGS "${LINK_FLAGS}")

add_library(utl_mat "uts_mat")
target_link_libraries(utl_mat utl_log ${UTD_UTIL_EXTERNAL_LIBRARIES})
SET_TARGET_PROPERTIES(utl_mat PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(utl_mat PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(utl_mat PROPERTIES LINK_FLAGS "${LINK_FLAGS}")

add_library(utl_mesh "uts_mesh")
target_link_libraries(utl_mesh utl_mat utl_log utl_mesh_remesh ${UTD_UTIL_EXTERNAL_LIBRARIES})
SET_TARGET_PROPERTIES(utl_mesh PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(utl_mesh PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(utl_mesh PROPERTIES LINK_FLAGS "${LINK_FLAGS}")

add_library(utl_io "uts_io_intf" "uts_io_compression" "miniz")
target_link_libraries(utl_io utl_io_results utl_log utl_bc utl_system utl_mesh ${UTD_UTIL_EXTERNAL_LIBRARIES})
SET_TARGET_PROPERTIES(utl_io PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(utl_io PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(utl_io PROPERTIES LINK_FLAGS "${LINK_FLAGS}")


# --- OpenCL support ---
if(${accel} STREQUAL "opencl")
set(SOURCE_FILES "uts_accel_intf")
endif()
# rest of source files
set(SOURCE_FILES ${SOURCE_FILES}
"uts_adapt"
"uts_coloring"
"uts_io_intf"
"uts_ls_intf"
"uts_util"
"uts_write_paraview"
"base64")

add_library(utl_util ${SOURCE_FILES} )
target_link_libraries(utl_util utl_io_results utl_system utl_bc utl_mat utl_mesh_remesh utl_io apl_std_lin ${UTD_UTIL_EXTERNAL_LIBRARIES})
SET_TARGET_PROPERTIES(utl_util PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(utl_util PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(utl_util PROPERTIES LINK_FLAGS "${LINK_FLAGS}")




enable_testing()

if(MODFEM_MPI STREQUAL "mpi" OR "${CMAKE_CXX_COMPILER}" MATCHES ".*mpi.*" )
    if(MODFEM_NEW_MPI)
            add_executable(MOD_FEM_test_util "uts_test.cpp")
            target_link_libraries(MOD_FEM_test_util utl_log utl_mat utl_util mml_t4_prism utl_system utl_io apl_std_lin pcl_mpi_safe mmpl_adapter ${UTD_UTIL_EXTERNAL_LIBRARIES})
            add_test(NAME MOD_FEM_test_util
            #COMMAND mpirun -np 2 MOD_FEM_test_util # here is some problem with parallel file access
            COMMAND MOD_FEM_test_util
            WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/..")
        endif()
else()
        add_executable(MOD_FEM_test_util "uts_test.cpp")
        target_link_libraries(MOD_FEM_test_util utl_log utl_mat utl_util mml_t4_prism utl_system utl_io apl_std_lin ${UTD_UTIL_EXTERNAL_LIBRARIES})
        add_test(NAME MOD_FEM_test_util
        COMMAND MOD_FEM_test_util
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/..")
endif()


