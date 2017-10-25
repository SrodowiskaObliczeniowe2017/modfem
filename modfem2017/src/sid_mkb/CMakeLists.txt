set(sid_mkb_srcdir "${CMAKE_SOURCE_DIR}/sid_mkb")

# -- Include directories --
list(APPEND SID_MKB_INCLUDE_DIRS
  ${sid_mkb_srcdir}
)

if(SID_MKB_INCLUDE_DIRS)
  set(SID_MKB_INCLUDE_DIRS ${SID_MKB_INCLUDE_DIRS} PARENT_SCOPE)
endif()

# -- Source files list --
collect_src_from_include_dirs(SID_MKB_INCLUDE_DIRS)

# -- Target --
include_directories( ${SID_MKB_INCLUDE_DIRS} ${SID_MKB_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
add_library(sil_mkb ${SOURCE_FILES} )

#target_link_libraries(sil_mkb ${SID_MKB_EXTERNAL_LIBRARIES}  lsl_mkb )
target_link_libraries(sil_mkb ${SID_MKB_EXTERNAL_LIBRARIES}  lsl_mkb )

SET_TARGET_PROPERTIES(sil_mkb PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(sil_mkb PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(sil_mkb PROPERTIES LINK_FLAGS "${LINK_FLAGS}")


