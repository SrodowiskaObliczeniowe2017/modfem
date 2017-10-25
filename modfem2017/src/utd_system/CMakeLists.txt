if(WIN32)
  set(utd_system_srcdir "${CMAKE_SOURCE_DIR}/utd_system/WIN")
else()
  set(utd_system_srcdir "${CMAKE_SOURCE_DIR}/utd_system/UNIX")
endif()

# -- Include directories --
list(APPEND UTD_SYSTEM_INCLUDE_DIRS
  ${utd_system_srcdir}
)

if(UTD_SYSTEM_INCLUDE_DIRS)
  set(UTD_SYSTEM_INCLUDE_DIRS ${UTD_SYSTEM_INCLUDE_DIRS} PARENT_SCOPE)
endif()

# -- Source files list --
collect_src_from_include_dirs(UTD_SYSTEM_INCLUDE_DIRS)

# -- Target --
include_directories( ${UTD_SYSTEM_INCLUDE_DIRS} ${UTD_SYSTEM_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

add_library( utl_system ${SOURCE_FILES} )

target_link_libraries(utl_system ${UTD_SYTEM_EXTERNAL_LIBRARIES})

SET_TARGET_PROPERTIES(utl_system PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(utl_system PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(utl_system PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
