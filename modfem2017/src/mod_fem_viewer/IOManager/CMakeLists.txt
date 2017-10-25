include_directories( inc )

AUX_SOURCE_DIRECTORY( inc IO_INC )
AUX_SOURCE_DIRECTORY( src IO_SRC )


if(WIN32)
  add_definitions(-D_CONSOLE -D_CRT_SECURE_NO_WARNINGS)
else(WIN32)
  #add_definitions(-g -O2 -fPIC -Wall -W)
endif(WIN32)

add_library( IOManager SHARED ${IO_INC} ${IO_SRC} )

#target_link_libraries( IOManager )
set_target_properties( IOManager PROPERTIES VERSION 1.0.0 SOVERSION 1 )
#install( TARGETS IOManager
#         LIBRARY DESTINATION lib
#         PERMISSIONS OWNER_READ GROUP_READ WORLD_READ)
