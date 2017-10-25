add_library(pcl_mpi_safe pcs_mpi_safe_intf.cpp ../include/pch_intf.h)

enable_testing()

add_executable(MOD_FEM_test_pcl "pcs_test.cpp")
target_link_libraries(MOD_FEM_test_pcl utl_util pcl_mpi_safe)


add_test(NAME MOD_FEM_test_pcl
COMMAND mpirun -np 2 ${CMAKE_BINARY_DIR}/MOD_FEM_test_pcl
WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/..")
