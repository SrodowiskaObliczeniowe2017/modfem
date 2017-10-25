set(lsd_mkb_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb")
set(lsd_mkb_core_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb/lsd_mkb_core")

set(lad_block_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb/lad_block")
set(lad_amg_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb/amg_ext")
#set(lsd_ns_supg_ext_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb/lsd_ns_supg_ext")
set(lsd_ns_supg_ext_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb/amg_ext")
set(lad_crs_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb/lad_crs")
set(lad_crs_generic_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb/lad_crs_generic")
set(lad_bcrs_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb/lad_bcrs")
set(lsd_mkb_pardiso_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb/lsd_mkb_pardiso" "${CMAKE_SOURCE_DIR}/lsd_mkb/lad_crs")
set(lsd_mkb_superlu_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb/lsd_mkb_superlu" "${CMAKE_SOURCE_DIR}/lsd_mkb/lad_crs" "${CMAKE_SOURCE_DIR}/lsd_mkb/lsd_mkb_superlu/superlu_threads")
set(lsd_mkb_superlu_seq_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb/lsd_mkb_superlu" "${CMAKE_SOURCE_DIR}/lsd_mkb/lad_crs" "${CMAKE_SOURCE_DIR}/lsd_mkb/lsd_mkb_superlu/superlu_seq")
set(lsd_mkb_viennacl_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb/lsd_mkb_viennacl" "${CMAKE_SOURCE_DIR}/lsd_mkb/lad_crs")

#set(lad_block_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb/lad_block")
#set(lad_amg_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb/amg_ext")
#set(lsd_ns_supg_ext_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb/lad_amg_ns_supg")
#set(lsd_ns_supg_ext_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb/amg_ext")


# -- Include directories --
list(APPEND LSD_MKB_INCLUDE_DIRS
  ${lsd_mkb_srcdir}
  ${lsd_mkb_core_srcdir}
  ${lad_block_srcdir}
  ${lad_crs_srcdir}
  ${lad_crs_generic_srcdir}
  ${lad_bcrs_srcdir}
  ${lad_amg_srcdir}
  ${lsd_ns_supg_ext_srcdir}
)


#if(MODFEM_ITER_SOLVER_MODULE STREQUAL "sil_mkb")
#  set(lad_block_srcdir "${CMAKE_SOURCE_DIR}/lsd_mkb/lad_amg")
#endif()

if(MODFEM_MKB_DIRECT_SOLVER_MODULE STREQUAL "PARDISO")
  message(" --> MKB INTERFACE selected direct solver: PARDISO")
  list(APPEND LSD_MKB_INCLUDE_DIRS ${lsd_mkb_pardiso_srcdir})
  add_definitions(-DPARDISO_MKB_DIRECT)
elseif(MODFEM_MKB_DIRECT_SOLVER_MODULE STREQUAL "VIENNACL")
  message(" --> MKB INTERFACE selected direct solver : VIENNACL")
  list(APPEND LSD_MKB_INCLUDE_DIRS ${lsd_mkb_viennacl_srcdir} ${VIENNACL_INCLUDE_DIRS})
  include_directories( ${SID_VIENNACL_CRS_INCLUDE_DIRS} ${SID_VIENNACL_CRS_EXTERNAL_INCLUDE_DIRS} ${VIENNACL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
  add_definitions(-DVIENNACL_MKB_DIRECT)
else()
  message(" --> MKB INTERFACE selected direct solver : SUPERLU")
  list(APPEND LSD_MKB_INCLUDE_DIRS ${lsd_mkb_superlu_srcdir})
#  list(APPEND LSD_MKB_INCLUDE_DIRS ${lsd_mkb_superlu_seq_srcdir})
  add_definitions(-DSUPERLU_MKB_DIRECT)
endif()

# SuperLU added as default direct solver for MKB (above)
# if(MODFEM_MKB_DIRECT_SOLVER_MODULE STREQUAL "SuperLU")
#   list(APPEND LSD_MKB_INCLUDE_DIRS ${lsd_mkb_superlu_srcdir})
#   add_definitions(-DSUPERLU_MKB_DIRECT)
# endif()

if(LSD_MKB_INCLUDE_DIRS)
  set(LSD_MKB_INCLUDE_DIRS ${LSD_MKB_INCLUDE_DIRS} PARENT_SCOPE)
endif()

# -- Source files list --
collect_src_from_include_dirs(LSD_MKB_INCLUDE_DIRS)
if(MSVC)  
  SET_SOURCE_FILES_PROPERTIES(${SOURCE_FILES} PROPERTIES LANGUAGE C)
  #SET_SOURCE_FILES_PROPERTIES("las_bcrs_intf.c" PROPERTIES LANGUAGE CXX)
  #SET_SOURCE_FILES_PROPERTIES("las_block_intf.c" PROPERTIES LANGUAGE CXX)
endif()

# -- Target --
include_directories( ${LSD_MKB_INCLUDE_DIRS} ${LSD_MKB_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
add_library(lsl_mkb ${SOURCE_FILES} )


if(MODFEM_ITER_SOLVER_MODULE STREQUAL "sil_mkb")
  target_link_libraries(lsl_mkb ${LSD_MKB_EXTERNAL_LIBRARIES} sil_mkb )
  if(MODFEM_USE_PETSC)
  	add_definitions(-DAMG_LAR_EXTENSIONS)
  	add_definitions(-DAMG_MKB_EXTENSIONS)
  	add_definitions(-DAMG_MKB_CORE_EXTENSIONS)
  	add_definitions(-DNS_SUPG_MKB_EXTENSIONS)
  	add_definitions(-DNS_SUPG_MKB_CORE_EXTENSIONS)
  	target_link_libraries(lsl_mkb ${LSD_MKB_EXTERNAL_LIBRARIES} amg_mkb )
  endif()
endif()


SET_TARGET_PROPERTIES(lsl_mkb PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(lsl_mkb PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(lsl_mkb PROPERTIES LINK_FLAGS "${LINK_FLAGS}")

#IF (WIN32)
#target_link_libraries(lsd_mkb sid_mkb blas lapack pthreadVC2)
#ELSE (WIN32)
#target_link_libraries(lsd_mkb sid_mkb )
#ENDIF (WIN32)
