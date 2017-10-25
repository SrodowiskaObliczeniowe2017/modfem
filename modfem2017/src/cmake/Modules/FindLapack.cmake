# Will define:
# - LAPACK_INCLUDE_DIRS - where to find LAPACK headers
# - LAPACK_LIBRARIES - List of libraries when using LAPACK
# - LAPACK_FOUND - True if LAPACK


include(LibFindMacros)

set(MKL_CORE_LIBRARY)
set(MKL_INTEL_LIBRARY)
set(MKL_LAPACK_LIBRARY)
set(MKL_SOLVER_LIBRARY)
set(MKL_THREAD_LIBRARY)
set(IOMP5_LIBRARY)
set(DYN_LIBS_BEFORE_GROUP)
set(DYN_LIBS_AFTER_GROUP)
set(LAPACK_PROCESS_LIBS)
set(LAPACK_PROCESS_INCLUDES)
set(M_LIBRARY)
set(DL_LIBRARY)


# find system library
find_library(M_LIBRARY m)
find_library(DL_LIBRARY dl)

if(WIN32)
  # supressing errors
  set(M_LIBRARY "")
  set(DL_LIBRARY "")
endif()
  
if(MODFEM_BLASLAPACK STREQUAL "MKL")
  link_directories(${MKL_LIBRARY_DIRS})

  find_path(LAPACK_INCLUDE_DIR
    NAMES mkl.h
    PATHS ${MKL_INCLUDE_DIRS})

  # mkl_intel
  find_library(MKL_INTEL_LIBRARY ${MKL_INTEL_LIB_NAME}
    PATHS ${MKL_LIBRARY_DIRS})

  # mkl_lapack
  find_library(MKL_LAPACK_LIBRARY ${MKL_LAPACK_LIB_NAME}
    PATHS ${MKL_LIBRARY_DIRS})

  if(MKL_LAPACK_LIBRARY)
    mark_as_advanced(MKL_LAPACK_LIBRARY)
  endif()

  # mkl_solver
  find_library(MKL_SOLVER_LIBRARY ${MKL_SOLVER_LIB_NAME}
    PATHS ${MKL_LIBRARY_DIRS}) #dla sequential :mkl_solver_lp64_sequential

  if(MKL_SOLVER_LIBRARY)
    mark_as_advanced(MKL_SOLVER_LIBRARY)
  endif()

  ################
  set(CMAKE_FIND_LIBRARY_SUFFIXES_ORIGINAL ${CMAKE_FIND_LIBRARY_SUFFIXES})
  
  if(MODFEM_USE_STATIC)
    if(UNIX)
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a")
      #message("STATIC 1")
    endif()
  endif()
 
  # iomp5
  if(MKL_IOMP5_LIB_NAME)
    find_library(IOMP5_LIBRARY ${MKL_IOMP5_LIB_NAME}
      PATHS ${MKL_LIBRARY_DIRS} NO_SYSTEM_ENVIRONMENT_PATH)
    if(WIN32)
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".dll")
      find_library(IOMP5_LIBRARY_DLL ${MKL_IOMP5_LIB_NAME}
	PATHS ${MKL_LIBRARY_DIRS} NO_SYSTEM_ENVIRONMENT_PATH)
    endif()
  endif()

  set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_ORIGINAL})
  #############

  # mkl_intel_thread / mkl_sequential
  find_library(MKL_THREAD_LIBRARY ${MKL_THREAD_LIB_NAME}
    PATHS ${MKL_LIBRARY_DIRS}) #albo samo mkl_sequential (dla WIN32 _dll i samo mkl_intel_thread_dll dla mt)

  # mkl_core
  find_library(MKL_CORE_LIBRARY ${MKL_CORE_LIB_NAME}
    PATHS ${MKL_LIBRARY_DIRS}) #pamietac wszedzie o _dll dla WIN32

  ###############	     
  # add libraries
  if(MODFEM_USE_STATIC)
    set(DYN_LIBS_BEFORE_GROUP MKL_SOLVER_LIBRARY)
  endif()

  if(MODFEM_USE_STATIC)
    list(APPEND INGROUP MKL_INTEL_LIBRARY MKL_THREAD_LIBRARY MKL_CORE_LIBRARY M_LIBRARY DL_LIBRARY)
  else()
    list(APPEND INGROUP MKL_SOLVER_LIBRARY MKL_INTEL_LIBRARY MKL_THREAD_LIBRARY MKL_CORE_LIBRARY)
  endif()

  if(IOMP5_LIBRARY) #lub MKL_IOMP5_LIB_NAME
    list(APPEND DYN_LIBS_AFTER_GROUP IOMP5_LIBRARY)
  endif()

  if(MSVC)
    find_library( INTEL_IFCORE libifcoremd PATHS ${MKL_LIBRARY_DIRS})
    list(APPEND DYN_LIBS_AFTER_GROUP INTEL_IFCORE)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".dll")
    find_library(INTEL_IFCORE_DLL libifcoremd PATHS ${MKL_LIBRARY_DIRS})
    find_library( INTEL_M_DLL libmmd PATHS ${MKL_LIBRARY_DIRS})
    set(CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_ORIGINAL})
  endif()

  if(MODFEM_USE_STATIC)
    set(LAPACK_PROCESS_LIBS ${DYN_LIBS_BEFORE_GROUP} STARTGROUP ${INGROUP} ENDGROUP ${DYN_LIBS_AFTER_GROUP})
  else()
    set(LAPACK_PROCESS_LIBS ${DYN_LIBS_BEFORE_GROUP} ${INGROUP} ${DYN_LIBS_AFTER_GROUP})
  endif()
 
elseif(MODFEM_BLASLAPACK STREQUAL "GENERIC")

  # Find lapack
  if(EXISTS ${LAPACK_DIRS})
    find_library(LAPACK_LIBRARY NAMES ${LAPACK_LIB_NAME} PATHS ${LAPACK_DIRS} NO_DEFAULT_PATH)
    if(${LAPACK_LIBRARY} STREQUAL "LAPACK_LIBRARY-NOTFOUND")
      message(WARNING "LAPACK library not found, try to find system library")
      find_library(LAPACK_LIBRARY NAMES ${LAPACK_LIB_NAME})
    endif()
  else()
    find_library(LAPACK_LIBRARY NAMES ${LAPACK_LIB_NAME}) 
  endif()
  if(${LAPACK_LIBRARY} STREQUAL "LAPACK_LIBRARY-NOTFOUND")
    message(FATAL_ERROR "Can not find LAPACK library!!!")
  endif()

  # Find blas
  if(EXISTS ${BLAS_DIRS})
    find_library(BLAS_LIBRARY NAMES ${BLAS_LIB_NAME} PATHS ${BLAS_DIRS} NO_DEFAULT_PATH)
    if(${BLAS_LIBRARY} STREQUAL "BLAS_LIBRARY-NOTFOUND")
      message(WARNING "BLAS library not found, try to find system library")
      find_library(BLAS_LIBRARY NAMES ${BLAS_LIB_NAME})
    endif()
  else()
    find_library(BLAS_LIBRARY NAMES ${BLAS_LIB_NAME})
  endif()
  if(${BLAS_LIBRARY} STREQUAL "BLAS_LIBRARY-NOTFOUND")
    message(FATAL_ERROR "Can not find BLAS library!!!")
  endif()

  # Find LIBF2C
  if(EXISTS ${LIBF2C_DIRS})
	find_library(LIBF2C_LIBRARY NAMES ${LIBF2C_LIB_NAME} PATHS ${LIBF2C_DIRS} NO_DEFAULT_PATH)
	message("LIBF2C_LIB_NAME: ${LIBF2C_LIB_NAME},  LIBF2C_DIRS:${LIBF2C_DIRS}" )
  else()
	find_library(LIBF2C_LIBRARY NAMES ${LIBF2C_LIB_NAME})
  endif()
   
  # End configuration
  set(LAPACK_PROCESS_LIBS LAPACK_LIBRARY BLAS_LIBRARY LIBF2C_LIBRARY)

elseif(MODFEM_BLASLAPACK STREQUAL "ACML")

  link_directories(${ACML_LIBRARY_DIRS})
  find_path(LAPACK_INCLUDE_DIR
    NAMES acml.h
    PATHS ${ACML_INCLUDE_DIRS})

  set(LAPACK_LIBRARY)
  if(MODFEM_USE_STATIC)
    find_library(LAPACK_LIBRARY acml PATHS ${ACML_LIBRARY_DIRS})
  else()
    find_library(LAPACK_LIBRARY acml PATHS ${ACML_LIBRARY_DIRS})
  endif()
  
  set(LAPACK_PROCESS_LIBS LAPACK_LIBRARY)
  
endif()

set(LAPACK_PROCESS_INCLUDES LAPACK_INCLUDE_DIR)
libfind_process(LAPACK)

if( MODFEM_BLASLAPACK STREQUAL "MKL" )
  find_package(Threads)
  list(APPEND LAPACK_LIBRARIES  ${CMAKE_THREAD_LIBS_INIT} ${M_LIBRARY} ${DL_LIBRARY})
elseif( MODFEM_BLASLAPACK STREQUAL "GENERIC" )
  list(APPEND LAPACK_LIBRARIES ${M_LIBRARY} ${DL_LIBRARY})
elseif( MODFEM_BLASLAPACK STREQUAL "ACML")
  list(APPEND LAPACK_LIBRARIES ${M_LIBRARY} ${DL_LIBRARY})
endif()

message("Lapack/MKL/ACML libraries: ${LAPACK_LIBRARIES}")
