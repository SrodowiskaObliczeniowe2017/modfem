# -- Exec files --
set(pdd_vof_srcdir "${CMAKE_SOURCE_DIR}/pdd_vof")
# -- Include directories --
list(APPEND PDD_VOF_INCLUDE_DIRS
  ${pdd_vof_srcdir}/input_output
  ${pdd_vof_srcdir}/materials
  ${pdd_vof_srcdir}/weak_formulation
  ${pdd_vof_srcdir}/include
  ${pdd_vof_srcdir}/
)
if(PDD_VOF_INCLUDE_DIRS)
  set(PDD_VOF_INCLUDE_DIRS ${PDD_VOF_INCLUDE_DIRS} PARENT_SCOPE)
endif()

# -- Source files list --
collect_src_from_include_dirs(PDD_VOF_INCLUDE_DIRS)
# -- Target --
include_directories(${PDD_VOF_INCLUDE_DIRS} ${PDD_VOF_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

add_library(pdl_vof ${SOURCE_FILES})
target_link_libraries(pdl_vof ${PDD_VOF_EXTERNAL_LIBRARIES})
SET_TARGET_PROPERTIES(pdl_vof PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(pdl_vof PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(pdl_vof PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
