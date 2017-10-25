# Will define:
# - VOROPP_INCLUDE_DIRS - where to find Voro++ headers
# - VOROPP_LIBRARIES - List of libraries when using Voro++
# - VOROPP_FOUND - True if Voro++ found

include(LibFindMacros)

set(VOROPP_LIBRARY)
set(VOROPP_PROCESS_LIBS)
set(VOROPP_PROCESS_INCLUDES)

link_directories(${VOROPP_LIBRARY_DIRS})

find_path(VOROPP_INCLUDE_DIR
          NAMES voro++.hh
          PATHS ${VOROPP_INCLUDE_DIRS}
		/usr/local/include
		/usr/local/include/voro++)

# mkl_intel
find_library(VOROPP_LIBRARY ${VOROPP_LIB_NAME}
               PATHS ${VOROPP_LIBRARY_DIRS})

if(VOROPP_LIBRARY)
  mark_as_advanced(VOROPP_LIBRARY)
endif()

set(VOROPP_PROCESS_LIBS VOROPP_LIBRARY)
set(VOROPP_PROCESS_INCLUDES VOROPP_INCLUDE_DIR)
libfind_process(VOROPP)

message("Voro++ library: ${VOROPP_LIBRARIES}")


