include_directories(inc)

AUX_SOURCE_DIRECTORY(inc MF_INC)
AUX_SOURCE_DIRECTORY(src MF_SRC)

remove_definitions(-D_WINDOWS)
add_definitions(-D_CONSOLE)
add_library(mod_fem ${MF_INC} ${MF_SRC})