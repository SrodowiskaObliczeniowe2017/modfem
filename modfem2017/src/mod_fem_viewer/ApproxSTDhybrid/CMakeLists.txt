set( STDh_SRC 
      ../include/ApproxModule.h
      ../ApproxModule/ApproxManager.h
      ApproxSTDhybrid.h
      ApproxSTDhybrid.cpp )



if(WIN32)
  if(MSVC)
	string(REPLACE "/D_WINDOWS" "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
	add_definitions(-D_CONSOLE -D_CRT_SECURE_NO_WARNINGS)
  endif(MSCV)
endif(WIN32)


add_library( ApproxSTDhybrid MODULE ${STDh_SRC} )
add_dependencies(ApproxSTDhybrid apd_std_lin mmd_t4_prism utd_util)
target_link_libraries( ApproxSTDhybrid apd_std_lin mmd_t4_prism utd_util)
set_target_properties(ApproxSTDhybrid PROPERTIES VERSION 1.0.0 SOVERSION 1)
#install(TARGETS ApproxSTDhybrid
#        LIBRARY DESTINATION lib
#        PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
