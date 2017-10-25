set( STD_SRC 
      ../include/ApproxModule.h
      ../ApproxModule/ApproxManager.h
      ./ApproxSTD.h
      ./ApproxSTD.cpp )

if(WIN32)
  add_definitions(-D_CONSOLE -D_CRT_SECURE_NO_WARNINGS)
else(WIN32)
  #add_definitions(-g -O2 -fPIC -Wall -W)
endif(WIN32)

add_definitions(-DNODEBUG -DNDEBUG)
add_library( ApproxSTD MODULE ${STD_SRC} )
add_dependencies(ApproxSTD apd_std_lin mmd_prism)
target_link_libraries( ApproxSTD apd_std_lin mmd_prism )
set_target_properties( ApproxSTD PROPERTIES VERSION 1.0.0 SOVERSION 1 )
#install( TARGETS ApproxSTD
#         LIBRARY DESTINATION lib
#         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
