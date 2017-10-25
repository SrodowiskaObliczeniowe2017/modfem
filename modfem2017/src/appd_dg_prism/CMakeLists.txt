set(appd_dg_prism_srcdir "${CMAKE_SOURCE_DIR}/appd_dg_prism")

# -- Include directories --
list(APPEND APPD_DG_PRISM_INCLUDE_DIRS
  ${appd_dg_prism_srcdir}
)

if(APPD_DG_PRISM_INCLUDE_DIRS)
  set(APPD_DG_PRISM_INCLUDE_DIRS ${APPD_DG_PRISM_INCLUDE_DIRS} PARENT_SCOPE)
endif()

# -- Source files list --
collect_src_from_include_dirs(APPD_DG_PRISM_INCLUDE_DIRS)

# -- Target --
include_directories( ${APPD_DG_PRISM_INCLUDE_DIRS} ${APPD_DG_PRISM_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
add_library( appl_dg_prism ${SOURCE_FILES} )

target_link_libraries(appl_dg_prism ${APPD_DG_PRISM_EXTERNAL_LIBRARIES})

SET_TARGET_PROPERTIES(appl_dg_prism PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(appl_dg_prism PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(appl_dg_prism PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
