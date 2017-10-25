set(ddd_parmetis_srcdir "${CMAKE_SOURCE_DIR}/ddd_parmetis")

# -- Include directories --
list(APPEND ddd_parmetis_INCLUDE_DIRS
  ${ddd_parmetis_srcdir}
  ${ddd_parmetis_srcdir}/ddd_metis_adapter
)

if(ddd_parmetis_INCLUDE_DIRS)
  set(ddd_parmetis_INCLUDE_DIRS ${mmpd_ADAPTER_INCLUDE_DIRS} PARENT_SCOPE)
endif()

list(APPEND ddd_metis_EXTERNAL_LIBRARIES "-Wl,--start-group ${METIS_LIBRARY} ${PARMETIS_LIBRARY} -Wl,--end-group")

# -- Source files list --
#collect_src_from_include_dirs(ddd_parmetis_INCLUDE_DIRS)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# -- Target --
include_directories(${ddd_parmetis_INCLUDE_DIRS} ${ddd_parmetis_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

add_library(ddl_parmetis "dds_metis_parmetis")
target_link_libraries(ddl_parmetis ${ddd_parmetis_EXTERNAL_LIBRARIES} ${METIS_LIBRARY} ${PARMETIS_LIBRARY})
SET_TARGET_PROPERTIES(ddl_parmetis PROPERTIES COMPILE_DEFINITIONS "${MODFEM_MPI_CODE_DEFINE}")
SET_TARGET_PROPERTIES(ddl_parmetis PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(ddl_parmetis PROPERTIES LINK_FLAGS "${LINK_FLAGS}")

if("${CMAKE_CXX_COMPILER}" MATCHES ".*mpi.*")
target_link_libraries(ddl_parmetis mmpl_adapter)
endif()

enable_testing()

#add_executable(MOD_FEM_test_ddl_parmetis "dds_parmetis_test")
#target_link_libraries(MOD_FEM_test_ddl_parmetis
#ddl_parmetis utl_io mml_t4_prism pcl_mpi_safe ${Boost_LIBRARIES})

#add_test(NAME MOD_FEM_test_ddl_parmetis
#COMMAND MOD_FEM_test_ddl_parmetis
#WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/..")

