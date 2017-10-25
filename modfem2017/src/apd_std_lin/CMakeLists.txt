set(apd_std_lin_srcdir "${CMAKE_SOURCE_DIR}/apd_std_lin")

# -- Include directories --
list(APPEND APD_STD_LIN_INCLUDE_DIRS
  ${apd_std_lin_srcdir}
)

if(APD_STD_LIN_INCLUDE_DIRS)
  set(APD_STD_LIN_INCLUDE_DIRS ${APD_STD_LIN_INCLUDE_DIRS} PARENT_SCOPE)
endif()

# -- Source files list --
collect_src_from_include_dirs(APD_STD_LIN_INCLUDE_DIRS)

# -- Target --
#include_directories( ${APD_STD_LIN_INCLUDE_DIRS} ${APD_STD_LIN_EXTERNAL_INCLUDE_DIRS})
include_directories( ${APD_STD_LIN_INCLUDE_DIRS} ${APD_STD_LIN_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
add_library( apl_std_lin ${SOURCE_FILES} )

target_link_libraries(apl_std_lin ${APD_STD_LIN_EXTERNAL_LIBRARIES})
target_link_libraries(apl_std_lin utl_util)

SET_TARGET_PROPERTIES(apl_std_lin PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(apl_std_lin PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(apl_std_lin PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
