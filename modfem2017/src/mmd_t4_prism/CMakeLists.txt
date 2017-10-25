set(mmd_t4_prism_srcdir "${CMAKE_SOURCE_DIR}/mmd_t4_prism")

# -- Include directories --
list(APPEND MMD_T4_PRISM_INCLUDE_DIRS
  ${mmd_t4_prism_srcdir}
  ${mmd_t4_prism_srcdir}/GeometryModule
  ${mmd_t4_prism_srcdir}/MeshModule
  ${mmd_t4_prism_srcdir}/MeshWrite
  ${mmd_t4_prism_srcdir}/MeshRead
)

if(MMD_T4_PRISM_INCLUDE_DIRS)
  set(MMD_T4_PRISM_INCLUDE_DIRS ${MMD_T4_PRISM_INCLUDE_DIRS} PARENT_SCOPE)
endif()

# -- Source files list --
#collect_src_from_include_dirs(MMD_T4_PRISM_INCLUDE_DIRS)
#No automatic collecting here
set(SOURCE_FILES
${mmd_t4_prism_srcdir}/mms_prism_t4_intf.cpp
${mmd_t4_prism_srcdir}/MeshModule/hHybridMesh
${mmd_t4_prism_srcdir}/MeshModule/hHybridMeshWithContacts
${mmd_t4_prism_srcdir}/MeshModule/mmr_vec3
${mmd_t4_prism_srcdir}/MeshModule/hObj
${mmd_t4_prism_srcdir}/MeshModule/hParent
${mmd_t4_prism_srcdir}/MeshModule/Vertex
${mmd_t4_prism_srcdir}/MeshModule/Edge
${mmd_t4_prism_srcdir}/MeshModule/Face3
${mmd_t4_prism_srcdir}/MeshModule/Face4
${mmd_t4_prism_srcdir}/MeshModule/ElemT4
${mmd_t4_prism_srcdir}/MeshModule/ElemPrism
${mmd_t4_prism_srcdir}/MeshModule/mesh_configuration
${mmd_t4_prism_srcdir}/MeshModule/SuperFastHash
${mmd_t4_prism_srcdir}/MeshWrite/MeshFileWriter
${mmd_t4_prism_srcdir}/MeshWrite/DmpFileExporter
${mmd_t4_prism_srcdir}/MeshWrite/MileninExporter
${mmd_t4_prism_srcdir}/MeshWrite/KazFileExporter
${mmd_t4_prism_srcdir}/MeshWrite/BinaryFileWriter
${mmd_t4_prism_srcdir}/MeshRead/IMeshReader
${mmd_t4_prism_srcdir}/MeshRead/MeshFileImporter
${mmd_t4_prism_srcdir}/MeshRead/NasFileImporter
${mmd_t4_prism_srcdir}/MeshRead/AdinaDatImporter
${mmd_t4_prism_srcdir}/MeshRead/DmpFileImporter
${mmd_t4_prism_srcdir}/MeshRead/KazFileImporter
${mmd_t4_prism_srcdir}/MeshRead/BinaryFileReader
${mmd_t4_prism_srcdir}/MeshRead/SimplePointReader
${mmd_t4_prism_srcdir}/MeshRead/InFileImporter
)


# -- Target --
include_directories(${MMD_T4_PRISM_INCLUDE_DIRS} ${MMD_T4_PRISM_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)

# Library
add_library(mml_t4_prism ${SOURCE_FILES})
target_link_libraries(mml_t4_prism utl_log ${MMD_T4_PRISM_EXTERNAL_LIBRARIES}) #which external libraries target requires is decided in main CMakeLists.txt
SET_TARGET_PROPERTIES(mml_t4_prism PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(mml_t4_prism PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(mml_t4_prism PROPERTIES LINK_FLAGS "${LINK_FLAGS}")

# Test
#add_executable(test_HybridMesh HybridMeshUnitTest)
#target_link_libraries(test_HybridMesh mmd_t4_prism)
#SET_TARGET_PROPERTIES(test_HybridMesh PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
#SET_TARGET_PROPERTIES(test_HybridMesh PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
#add_test(mmd_t4_prism_test test_HybridMesh mms_intf)
