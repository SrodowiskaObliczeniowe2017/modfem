set( DG_SRC
	../include/ApproxModule.h
    ../ApproxModule/ApproxManager.h 
    ./ApproxDG.h
    ./ApproxDG.cpp )

if(WIN32)
	string(REPLACE "/D_WINDOWS" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
	add_definitions(-D_CONSOLE -D_CRT_SECURE_NO_WARNINGS)
endif(WIN32)


add_library( ApproxDG MODULE ${DG_SRC} )
add_dependencies(ApproxDG apd_dg_prism mmd_prism)
target_link_libraries( ApproxDG apd_dg_prism mmd_prism)
set_target_properties(ApproxDG PROPERTIES VERSION 1.0.0 SOVERSION 1)
#install(TARGETS ApproxDG
#        LIBRARY DESTINATION lib
#        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
        

