set(amg_mkb_srcdir "${CMAKE_SOURCE_DIR}/amg_mkb")

# -- Include directories --
list(APPEND AMG_MKB_INCLUDE_DIRS
  ${amg_mkb_srcdir}
  ${amg_mkb_srcdir}/amg
  ${amg_mkb_srcdir}/amg/Interpolation
  ${amg_mkb_srcdir}/lad_petsc
  ${amg_mkb_srcdir}/lad_petsc/ns_supg
  ${amg_mkb_srcdir}/lsd_amg
)

add_subdirectory(test)
 
if(AMG_MKB_INCLUDE_DIRS)
  set(AMG_MKB_INCLUDE_DIRS ${AMG_MKB_INCLUDE_DIRS} PARENT_SCOPE)
endif()

# -- Source files list --
collect_src_from_include_dirs(AMG_MKB_INCLUDE_DIRS)

# -- Target --
find_package(MPI REQUIRED)
include_directories( ${AMG_MKB_INCLUDE_DIRS} ${AMG_MKB_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include ${MPI_INCLUDE_PATH})
include_directories(${PETSC_HOME}/include ${PETSC_HOME}/${PETSC_ARCH}/include)
add_library(amg_mkb ${SOURCE_FILES} )

link_directories(${PETSC_HOME}/${PETSC_ARCH}/lib)
target_link_libraries(amg_mkb ${AMG_MKB_EXTERNAL_LIBRARIES} sil_mkb ${PETSC_HOME}/${PETSC_ARCH}/lib/libpetsc.so ${BLAS_LIB_NAME} lapack)

SET_TARGET_PROPERTIES(amg_mkb PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(amg_mkb PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(amg_mkb PROPERTIES LINK_FLAGS "${AMG_MKB_LINK_FLAGS}")



