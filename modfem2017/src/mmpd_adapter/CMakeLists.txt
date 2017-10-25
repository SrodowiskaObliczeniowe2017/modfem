set(mmpd_adapter_srcdir "${CMAKE_SOURCE_DIR}/mmpd_adapter")

# -- Include directories --
list(APPEND mmpd_ADAPTER_INCLUDE_DIRS
  ${mmpd_adapter_srcdir}
  ${mmpd_adapter_srcdir}/fpcm
)

if(mmpd_ADAPTER_INCLUDE_DIRS)
  set(mmpd__INCLUDE_DIRS ${mmpd_ADAPTER_INCLUDE_DIRS} PARENT_SCOPE)
endif()

#list(APPEND mmpd_ADAPTER_EXTERNAL_LIBRARIES "-Wl,--start-group ${METIS_LIBRARY} ${PARMETIS_LIBRARY} -Wl,--end-group")

# -- Source files list --
#collect_src_from_include_dirs(mmpd_ADAPTER_INCLUDE_DIRS)
#list(REMOVE_ITEM SOURCE_FILES mmps_adapter_test.cpp)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

# -- Target --
include_directories(${mmpd_ADAPTER_INCLUDE_DIRS} ${mmpd_ADAPTER_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

add_library(mmpl_adapter "mmph_adapter" "mmps_adapter_intf" "Transferer" "TransfererWithOwnership" "CompositeTransfererWithOwnership" "EntitiesOwnership"
"mmps_ipid" "fpcm/compressed_mesh")
target_link_libraries(mmpl_adapter ddl_parmetis ${mmpd_ADAPTER_EXTERNAL_LIBRARIES})
SET_TARGET_PROPERTIES(mmpl_adapter PROPERTIES COMPILE_DEFINITIONS "${MODFEM_MPI_CODE_DEFINE}")
SET_TARGET_PROPERTIES(mmpl_adapter PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(mmpl_adapter PROPERTIES LINK_FLAGS "${LINK_FLAGS}")

enable_testing()

add_executable(MOD_FEM_test_mmpl_adapter "mmps_adapter_test.cpp")
target_link_libraries(MOD_FEM_test_mmpl_adapter utl_log utl_io mml_t4_prism pcl_mpi_safe mmpl_adapter  ${Boost_LIBRARIES})

add_test(NAME MOD_FEM_test_mmpl_adapter
COMMAND mpirun -np 2 ${CMAKE_BINARY_DIR}/MOD_FEM_test_mmpl_adapter
WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/..")

