set(appd_std_lin_srcdir "${CMAKE_SOURCE_DIR}/appd_std_lin")

# -- Include directories --
list(APPEND APPD_STD_LIN_INCLUDE_DIRS
  ${appd_std_lin_srcdir}
)

if(APPD_STD_LIN_INCLUDE_DIRS)
  set(APPD_STD_LIN_INCLUDE_DIRS ${APPD_STD_LIN_INCLUDE_DIRS} PARENT_SCOPE)
endif()

# -- Source files list --
collect_src_from_include_dirs(APPD_STD_LIN_INCLUDE_DIRS)

# -- Target --
include_directories( ${APPD_STD_LIN_INCLUDE_DIRS} ${APPD_STD_LIN_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
add_library( appl_std_lin ${SOURCE_FILES} )

target_link_libraries(appl_std_lin ${APPD_STD_LIN_EXTERNAL_LIBRARIES})
target_link_libraries(appl_std_lin utl_util)

SET_TARGET_PROPERTIES(appl_std_lin PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(appl_std_lin PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(appl_std_lin PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
