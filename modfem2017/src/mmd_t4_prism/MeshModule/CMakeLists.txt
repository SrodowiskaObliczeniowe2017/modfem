add_library(hHybridMesh 
mmr_vec3
SuperFastHash
StaticPool.hpp ArrayPool.hpp
hHybridMesh.h hHybridMesh.cpp 
EntityAttributes.hpp
hObj.h hObj.cpp
../Field.hpp
../Common.h
hParent.hpp hParent.cpp
Vertex.h Vertex.cpp
Edge.h Edge.cpp
Face3.h Face3.cpp
Face4.h Face4.cpp
ElemT4.h ElemT4.cpp
ElemPrism.h ElemPrism.cpp
MeshEntity.hpp
mesh_configuration.h mesh_configuration.cpp
)

target_link_libraries(hHybridMesh MeshRead MeshWrite)
