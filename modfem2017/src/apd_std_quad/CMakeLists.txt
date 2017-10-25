set(apd_std_quad_srcdir "${CMAKE_SOURCE_DIR}/apd_std_quad")

# -- Include directories --
list(APPEND APD_STD_QUAD_INCLUDE_DIRS
  ${apd_std_quad_srcdir}
)

if(APD_STD_QUAD_INCLUDE_DIRS)
  set(APD_STD_QUAD_INCLUDE_DIRS ${APD_STD_QUAD_INCLUDE_DIRS} PARENT_SCOPE)
endif()

# -- Source files list --
collect_src_from_include_dirs(APD_STD_QUAD_INCLUDE_DIRS)

# -- Target --
#include_directories( ${APD_STD_QUAD_INCLUDE_DIRS} ${APD_STD_QUAD_EXTERNAL_INCLUDE_DIRS})
include_directories( ${APD_STD_QUAD_INCLUDE_DIRS} ${APD_STD_QUAD_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
add_library( apl_std_quad ${SOURCE_FILES} )

target_link_libraries(apl_std_quad ${APD_STD_QUAD_EXTERNAL_LIBRARIES})
target_link_libraries(apl_std_quad utl_util)

SET_TARGET_PROPERTIES(apl_std_quad PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(apl_std_quad PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(apl_std_quad PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
