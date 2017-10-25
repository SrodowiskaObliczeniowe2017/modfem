# Will define:
# - LIBCONFIG_INCLUDE_DIRS - where to find Libconfig headers
# - LIBCONFIG_LIBRARIES - List of libraries when using Libconfig
# - LIBCONFIG_FOUND - True if Libconfig found

include(LibFindMacros)

set(LIBCONFIG_LIBRARY)
set(LIBCONFIG_PROCESS_LIBS)
set(LIBCONFIG_PROCESS_INCLUDES)

link_directories(${LIBCONFIG_LIBRARY_DIRS})

find_path(LIBCONFIG_INCLUDE_DIR
          NAMES libconfig.h
          PATHS ${LIBCONFIG_INCLUDE_DIRS})

# mkl_intel
find_library(LIBCONFIG_LIBRARY ${LIBCONFIG_LIB_NAME}
               PATHS ${LIBCONFIG_LIBRARY_DIRS})

if(LIBCONFIG_LIBRARY)
  mark_as_advanced(LIBCONFIG_LIBRARY)
endif()

set(LIBCONFIG_PROCESS_LIBS LIBCONFIG_LIBRARY)
set(LIBCONFIG_PROCESS_INCLUDES LIBCONFIG_INCLUDE_DIR)
libfind_process(LIBCONFIG)

message("Libconfig library: ${LIBCONFIG_LIBRARIES}")


