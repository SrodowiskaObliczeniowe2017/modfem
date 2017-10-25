set(tmd_ocl_srcdir "${CMAKE_SOURCE_DIR}/tmd_opencl")

# -- Include directories --
list(APPEND TMD_OCL_INCLUDE_DIRS
  ${tmd_ocl_srcdir}
)

if(TMD_OCL_INCLUDE_DIRS)
  set(TMD_OCL_INCLUDE_DIRS ${TMD_OCL_INCLUDE_DIRS} PARENT_SCOPE)
endif()

# -- Source files list --
collect_src_from_include_dirs(TMD_OCL_INCLUDE_DIRS)

# -- Target --
include_directories( ${TMD_OCL_INCLUDE_DIRS} ${TMD_OCL_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
add_library(tml_ocl ${SOURCE_FILES} )

target_link_libraries(tml_ocl ${TMD_OCL_EXTERNAL_LIBRARIES})

SET_TARGET_PROPERTIES(tml_ocl PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(tml_ocl PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(tml_ocl PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
