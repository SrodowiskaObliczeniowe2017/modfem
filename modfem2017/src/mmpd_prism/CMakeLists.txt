set(mmpd_prism_srcdir "${CMAKE_SOURCE_DIR}/mmpd_prism")

# -- Include directories --
list(APPEND mmpd_PRISM_INCLUDE_DIRS
  ${mmpd_prism_srcdir}
  ${mmpd_prism_srcdir}/ddd_manager
)

if(mmpd_PRISM_INCLUDE_DIRS)
  set(mmpd_PRISM_INCLUDE_DIRS ${mmpd_PRISM_INCLUDE_DIRS} PARENT_SCOPE)
endif()

# -- Source files list --
collect_src_from_include_dirs(mmpd_PRISM_INCLUDE_DIRS)

# -- Target --
include_directories(${mmpd_PRISM_INCLUDE_DIRS} ${mmpd_PRISM_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

add_library(mmpl_prism ${SOURCE_FILES})
target_link_libraries(mmpl_prism ${mmpd_PRISM_EXTERNAL_LIBRARIES})
SET_TARGET_PROPERTIES(mmpl_prism PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(mmpl_prism PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(mmpl_prism PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
