set(apd_dg_prism_srcdir "${CMAKE_SOURCE_DIR}/apd_dg_prism")

# -- Include directories --
list(APPEND APD_DG_PRISM_INCLUDE_DIRS
  ${apd_dg_prism_srcdir}
)

if(APD_DG_PRISM_INCLUDE_DIRS)
  set(APD_DG_PRISM_INCLUDE_DIRS ${APD_DG_PRISM_INCLUDE_DIRS} PARENT_SCOPE)
endif()

# -- Source files list --
collect_src_from_include_dirs(APD_DG_PRISM_INCLUDE_DIRS)

# -- Target --
#include_directories( ${APD_DG_PRISM_INCLUDE_DIRS} ${APD_DG_PRISM_EXTERNAL_INCLUDE_DIRS})
include_directories( ${APD_DG_PRISM_INCLUDE_DIRS} ${APD_DG_PRISM_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
add_library( apl_dg_prism ${SOURCE_FILES} )

target_link_libraries(apl_dg_prism ${APD_DG_PRISM_EXTERNAL_LIBRARIES})
#target_link_libraries(apl_dg_prism utd_util)

SET_TARGET_PROPERTIES(apl_dg_prism PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(apl_dg_prism PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(apl_dg_prism PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
