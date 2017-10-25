set(amg_mkb_testdir "${CMAKE_SOURCE_DIR}/amg_mkb/test")

#test
list(APPEND AMG_MKB_INCLUDE_DIRS
  ${amg_mkb_testdir}/
)

include_directories(${PETSC_HOME}/include ${PETSC_HOME}/${PETSC_ARCH}/include)
find_package(MPI REQUIRED)
include_directories( ${AMG_MKB_INCLUDE_DIRS} ${AMG_MKB_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include ${MPI_INCLUDE_PATH})

add_executable(spai_test ${amg_mkb_testdir}/util/ApproximateInverseTest.cpp)
add_test(amg_test spai_test)
target_link_libraries(spai_test amg_mkb boost_unit_test_framework) 

add_executable(rscf_splitter_test ${amg_mkb_testdir}/amg/RSCFSplitterTest.cpp)
add_test(amg_test rscf_splitter_test)
target_link_libraries(rscf_splitter_test amg_mkb boost_unit_test_framework)

add_executable(amg_solver_structure_test ${amg_mkb_testdir}/amg/AMGSolverStructureTest.cpp)
add_test(amg_test amg_solver_structure_test)
target_link_libraries(amg_solver_structure_test amg_mkb boost_unit_test_framework)