set(sid_lapack_srcdir "${CMAKE_SOURCE_DIR}/sid_lapack")

# -- Include directories --
list(APPEND SID_LAPACK_INCLUDE_DIRS
  ${sid_lapack_srcdir}
)

if(SID_LAPACK_INCLUDE_DIRS)
  set(SID_LAPACK_INCLUDE_DIRS ${SID_LAPACK_INCLUDE_DIRS} PARENT_SCOPE)
endif()

# -- Source files list --
collect_src_from_include_dirs(SID_LAPACK_INCLUDE_DIRS)

# -- Target --
include_directories( ${SID_LAPACK_INCLUDE_DIRS} ${SID_LAPACK_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
add_library( sil_lapack ${SOURCE_FILES} )

target_link_libraries(sil_lapack ${SID_LAPACK_EXTERNAL_LIBRARIES})

SET_TARGET_PROPERTIES(sil_lapack PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(sil_lapack PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(sil_lapack PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
