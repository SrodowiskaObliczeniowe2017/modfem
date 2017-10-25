
set(mmd_remesh_srcdir "${CMAKE_SOURCE_DIR}/mmd_remesh")

list(APPEND MMD_REMESH_INCLUDE_DIRS
  ${mmd_remesh_srcdir}
)

if(MMD_REMESH_INCLUDE_DIRS)
  set(MMD_REMESH_INCLUDE_DIRS ${MMD_REMESH_INCLUDE_DIRS} PARENT_SCOPE)
endif()
set(SOURCE_FILES
${mmd_remesh_srcdir}/teksty
${mmd_remesh_srcdir}/wektor
${mmd_remesh_srcdir}/plansza
${mmd_remesh_srcdir}/rozrost
${mmd_remesh_srcdir}/losowanie
${mmd_remesh_srcdir}/delanouy3d
${mmd_remesh_srcdir}/figury
${mmd_remesh_srcdir}/ticooMesh3D
${mmd_remesh_srcdir}/ziarno
)


add_library(utl_mesh_remesh ${SOURCE_FILES})
target_link_libraries(utl_mesh_remesh ${MMD_REMESH_EXTERNAL_LIBRARIES})

set(SOURCE_FILES
${mmd_remesh_srcdir}/mms_prism_t4_remesh_intf
)

# -- Source files list --
#collect_src_from_include_dirs(MMD_REMESH_INCLUDE_DIRS)

# -- Target --
include_directories(${MMD_REMESH_INCLUDE_DIRS} ${MMD_REMESH_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)


# Library
add_library(mml_remesh ${SOURCE_FILES})
target_link_libraries(mml_remesh utl_mesh_remesh ${MMD_REMESH_EXTERNAL_LIBRARIES}) #which external libraries target requires is decided in main CMakeLists.txt

SET_TARGET_PROPERTIES(mml_remesh PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(mml_remesh PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(mml_remesh PROPERTIES LINK_FLAGS "${LINK_FLAGS}")







