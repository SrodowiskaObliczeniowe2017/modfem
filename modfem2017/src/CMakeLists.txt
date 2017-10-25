# THIS IS NEW CMAKELISTS FILE !!!
cmake_minimum_required(VERSION 2.6)

message("")
message("*** Entering main CMakeLists ***")
message("")

set(UTS_TEST FALSE)



# -----------------------------------------------------------------------------
# ----------------------------- Project configuration -------------------------
# -----------------------------------------------------------------------------

# -------------------------
# -- exec list and tools --
# -------------------------
set(MOD_FEM_EXEC_LIST)
set(MOD_FEM_TOOLS_LIST)

# Environment
set(ENV_PLATFORM $ENV{MOD_FEM_ARCH_CMAKE})
# perhaps MOD_FEM_ARCH was passed as argument?
if(NOT DEFINED ENV_PLATFORM)
  set(ENV_PLATFORM ${MOD_FEM_ARCH_CMAKE})
endif()

if(NOT DEFINED ENV_PLATFORM)

  # 1.2 PRO USER MODE (correct folder name)
  # Correct folder name has form of:
  # "MF_ARCH"_"MF_MPI"_"MF_ACCEL"_"MF_CC"_"MF_CXX"
  # ^.*_.*_.*_.*_.*$
  if( ${CMAKE_BINARY_DIR} MATCHES "^.*_.*_.*_.*_.*$")
    set(MF_USER_MODE "pro")
    SET(MF_DIR_NAME ${CMAKE_BINARY_DIR})
    get_filename_component(MF_DIR_NAME  ${MF_DIR_NAME} NAME)
    message(STATUS "Dir name: ${MF_DIR_NAME}")
    STRING(LENGTH ${MF_DIR_NAME} MF_DIR_NAME_LENGTH)

    STRING(REGEX MATCH "[^_]+$" MF_CXX ${MF_DIR_NAME})
    message(STATUS "MF_CXX=${MF_CXX}")
    STRING(LENGTH ${MF_CXX} TMP_LEN)
    MATH(EXPR LEN "${MF_DIR_NAME_LENGTH} - ${TMP_LEN} - 1")
    STRING(SUBSTRING ${MF_DIR_NAME} 0 ${LEN} MF_DIR_NAME )

    STRING(REGEX MATCH "[^_]+$" MF_CC ${MF_DIR_NAME})
    message(STATUS "MF_CC=${MF_CC}")
    STRING(LENGTH ${MF_CC} TMP_LEN)
    STRING(LENGTH ${MF_DIR_NAME} MF_DIR_NAME_LENGTH)
    MATH(EXPR LEN "${MF_DIR_NAME_LENGTH} - ${TMP_LEN} - 1")
    STRING(SUBSTRING ${MF_DIR_NAME} 0 ${LEN} MF_DIR_NAME )

    STRING(REGEX MATCH "[^_]+$" MF_ACCEL ${MF_DIR_NAME})
    message(STATUS "MF_ACCEL=${MF_ACCEL}")
    STRING(LENGTH ${MF_ACCEL} TMP_LEN)
    STRING(LENGTH ${MF_DIR_NAME} MF_DIR_NAME_LENGTH)
    MATH(EXPR LEN "${MF_DIR_NAME_LENGTH} - ${TMP_LEN} - 1")
    STRING(SUBSTRING ${MF_DIR_NAME} 0 ${LEN} MF_DIR_NAME )

    STRING(REGEX MATCH "[^_]+$" MF_MPI ${MF_DIR_NAME})
    message(STATUS "MF_MPI=${MF_MPI}")
    STRING(LENGTH ${MF_MPI} TMP_LEN)
    STRING(LENGTH ${MF_DIR_NAME} MF_DIR_NAME_LENGTH)
    MATH(EXPR LEN "${MF_DIR_NAME_LENGTH} - ${TMP_LEN} - 1")
    STRING(SUBSTRING ${MF_DIR_NAME} 0 ${LEN} MF_DIR_NAME )
    SET(MF_ARCH ${MF_DIR_NAME})

    set(MODFEM_MPI ${MF_MPI})
    set(MODFEM_ACCEL ${MF_ACCEL})
    set(ENV_PLATFORM ${MF_ARCH})

    set(CMAKE_CXX_COMPILER ${MF_CXX})
    set(CMAKE_C_COMPILER ${MF_CC})
    
  endif()
else()
  if(DEFINED CC)
    set(CMAKE_C_COMPILER ${CC})
  endif()
  if(DEFINED CXX)
    set(CMAKE_CXX_COMPILER ${CXX})
  endif()
endif()

# Additional folders to search for CMake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Platforms/")

# -- PLATFORM --
# --------------
if(DEFINED ENV_PLATFORM)
  if(ENV_PLATFORM)
    if(EXISTS "${CMAKE_SOURCE_DIR}/cmake/Platforms/${ENV_PLATFORM}.cmake")
      include(${ENV_PLATFORM} CUSTOMPLATFORM_EXISTS)
      message(STATUS "MOD_FEM_ARCH_CMAKE is '${ENV_PLATFORM}'") # STATUS ????
      message(STATUS "File is ${CMAKE_SOURCE_DIR}/cmake/Platforms/${ENV_PLATFORM}.cmake") # STATUS ????
    else()
      message("MOD_FEM_ARCH_CMAKE is '${ENV_PLATFORM}'")
      message(FATAL_ERROR "File ${CMAKE_SOURCE_DIR}/cmake/Platforms/${ENV_PLATFORM}.cmake doesn't exist!")
    endif()  
  else()
    message(FATAL_ERROR "MOD_FEM_ARCH_CMAKE is not defined!")
  endif()
else()
  message(FATAL_ERROR "MOD_FEM_ARCH_CMAKE is not defined!")
endif()


# Set default CMAKE_BUILD_TYPE
if(DEFINED CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose build type")
else()
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose build type")
endif()
mark_as_advanced(CLEAR CMAKE_VERBOSE_MAKEFILE)

# Define project with languages
project(mod_fem_code)
enable_language(CXX)
enable_language(C)

# -- OPERATIG SYSTEM DEPENDENT --
# --------------
include("ModFEM_OS.cmake")

#Create postfix
set(MF_POSTFIX "")
if(MODFEM_MPI STREQUAL "mpi")
  set(MF_POSTFIX "${MF_POSTFIX}_mpi")
endif()
if(MODFEM_ACCEL STREQUAL "opencl")
  set(MF_POSTFIX "${MF_POSTFIX}_ocl")
endif()

#Binaries appendixes
set(RELEASE_POSTFIX "${MF_POSTFIX}" )
set(CMAKE_RELEASE_POSTFIX "${MF_POSTFIX}" )

set(DEBUG_POSTFIX "${MF_POSTFIX}_d" )
set(CMAKE_DEBUG_POSTFIX "${MF_POSTFIX}_d" )

set(RELWITHDEBINFO_POSTFIX "${MF_POSTFIX}_i" )
set(CMAKE_RELWITHDEBINFO_POSTFIX "${MF_POSTFIX}_i" )

# Builds directories
set(ARCHIVE_OUTPUT_PATH ${CMAKE_BUILD_TYPE} CACHE PATH "Subdirectory for archives" FORCE)
set(LIBRARY_OUTPUT_PATH ${CMAKE_BUILD_TYPE} CACHE PATH "Subdirectory for libraries" FORCE)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR} CACHE PATH "Subdirectory for executables" FORCE)
mark_as_advanced(ARCHIVE_OUTPUT_PATH LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

# Default compilation flags
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
set(CMAKE_C_FLAGS_RELEASE "-O2 -DNDEBUG")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}") 
if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL GNU )
  set(CMAKE_CXX_FLAGS     "-Wno-unused-variable -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-comment ${CMAKE_CXX_FLAGS} -pipe -ldl")
  set(CMAKE_C_FLAGS       "-Wno-unused-variable -Wno-unused-parameter -Wno-unused-but-set-variable -Wno-comment ${CMAKE_C_FLAGS}  -pipe -ldl")
elseif( ${CMAKE_CXX_COMPILER_ID} STREQUAL Intel )
  set(CMAKE_CXX_FLAGS     "${CMAKE_CXX_FLAGS}")
elseif(MSVC)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /Qlowercase")  # For correct call to ilut,...
endif()
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Dfloat='double'")

#Time measure flags
if(DEFINED TIME_TEST)
  if(TIME_TEST)
    add_definitions(-DTIME_TEST)
  endif()
endif()
if(DEFINED TIME_TEST_2)
  if(TIME_TEST_2)
    add_definitions(-DTIME_TEST_2)
  endif()
endif()

# OPENCL_TASK - hold LAPLACE or TEST_SCALAR flag

message("MOD_FEM_ARCH_CMAKE is '${MODFEM_ACCEL}'")
if(${MODFEM_ACCEL} STREQUAL "opencl")
  if(DEFINED LAPLACE AND NOT DEFINED TEST_SCALAR)
    if(LAPLACE)
      add_definitions(-DLAPLACE)
    endif()
  elseif(DEFINED TEST_SCALAR AND NOT DEFINED LAPLACE)
    if(TEST_SCALAR)
      add_definitions(-DTEST_SCALAR)
    endif()
  elseif(DEFINED LAPLACE AND DEFINED TEST_SCALAR)
    if(LAPLACE AND NOT TEST_SCALAR)
      add_definitions(-DLAPLACE)
    elseif(TEST_SCALAR AND NOT LAPLACE)
      add_definitions(-DTEST_SCALAR)
    endif()
  endif()
endif()


if(DEFINED LAPLACE AND NOT DEFINED TEST_SCALAR)
  if(LAPLACE)
    set(MF_OPENCL_TASK "-DLAPLACE=\"ON\"")
  endif()
  unset(LAPLACE CACHE)
elseif(DEFINED TEST_SCALAR AND NOT DEFINED LAPLACE)
  if(TEST_SCALAR)
    set(MF_OPENCL_TASK "-DTEST_SCALAR=\"ON\"")
  endif()
  unset(TEST_SCALAR CACHE)
elseif(DEFINED LAPLACE AND DEFINED TEST_SCALAR)
  if(LAPLACE AND NOT TEST_SCALAR)
    set(MF_OPENCL_TASK "-DLAPLACE=\"ON\"")
  elseif(TEST_SCALAR AND NOT LAPLACE)
    set(MF_OPENCL_TASK "-DTEST_SCALAR=\"ON\"")
  endif()
  unset(LAPLACE CACHE)
  unset(TEST_SCALAR CACHE)
endif()


# Some flags in release/debug mode
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
#    set(CMAKE_CXX_FLAGS "-O2")
  else()
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  endif()

  if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  else()
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  endif()

elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "-O0 -g ${CMAKE_CXX_FLAGS}")
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel" )
    set(CMAKE_CXX_FLAGS "-O0 -g ${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  endif()

  if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_C_FLAGS "-O0 -g ${CMAKE_C_FLAGS}")
  elseif(CMAKE_C_COMPILER_ID STREQUAL "Intel" )
    set(CMAKE_C_FLAGS "-O0 -g ${CMAKE_C_FLAGS}")
  else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  endif()

  add_definitions(-DDEBUG)
endif()

# Extra debug flag
  if(DEFINED DEBUG_APM)
    add_definitions(-DDEBUG_APM)
  endif()
  if(DEFINED DEBUG_LSM)
    add_definitions(-DDEBUG_LSM)
  endif()
  if(DEFINED DEBUG_MMM)
    add_definitions(-DDEBUG_MMM)
  endif()
  if(DEFINED DEBUG_PCM)
    add_definitions(-DDEBUG_PCM)
  endif()
  if(DEFINED DEBUG_SIM)
    add_definitions(-DDEBUG_SIM)
  endif()
  if(DEFINED DEBUG_TMM)
    add_definitions(-DDEBUG_TMM)
  endif()

# SVN revision 
# -------------------------
include(FindSubversion OPTIONAL)
IF(Subversion_FOUND)
  get_filename_component(CMAKE_SOURCE_DIR_REALPATH ${CMAKE_SOURCE_DIR} REALPATH)
  set(ERR_VAR "")
  execute_process(COMMAND svn info 
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR_REALPATH}"
    ERROR_VARIABLE ERR_VAR
    )
  if(NOT ERR_VAR)
    Subversion_WC_INFO(${CMAKE_SOURCE_DIR_REALPATH} ModFEM)
    SET(SVN_REVISION "${ModFEM_WC_REVISION}")
    file(WRITE "${CMAKE_SOURCE_DIR_REALPATH}/include/svnversion.h" "const int SVNVERSION = ${SVN_REVISION} ;")
  endif()
ENDIF(Subversion_FOUND)

if(DEFINED SVN_REVISION)
    MESSAGE("Current revision is ${SVN_REVISION}")
else()
    MESSAGE("Current revision is not sync with local svn.")
endif()

# GIT COMMIT
# -------------------------
find_package(Git)
get_filename_component(CMAKE_SOURCE_DIR_REALPATH ${CMAKE_SOURCE_DIR} REALPATH)
IF(GIT_FOUND)
  EXECUTE_PROCESS(
    COMMAND ${GIT_EXECUTABLE} rev-parse --short HEAD
    WORKING_DIRECTORY "${CMAKE_SOURCE_DIR_REALPATH}"
    OUTPUT_VARIABLE "GIT_COMMIT"
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
    file(WRITE "${CMAKE_SOURCE_DIR_REALPATH}/include/gitversion.h" "const char GITVERSION[] = \"${GIT_COMMIT}\" ;")
  MESSAGE( STATUS "Git version: ${GIT_COMMIT}" )
ELSE(GIT_FOUND)
  MESSAGE("Git not found.")
  file(WRITE "${CMAKE_SOURCE_DIR_REALPATH}/include/gitversion.h" "const char GITVERSION[] = \"unversioned\" ;")
ENDIF(GIT_FOUND)



# -- MODFEM_TEST_COVERAGE --
# -------------------------
set(MODFEM_TEST_COVERAGE_DEFAULT FALSE)
set(MODFEM_TEST_COVERAGE ${MODFEM_TEST_COVERAGE_DEFAULT} CACHE BOOL "Compilation with coverage tools (-fprofile-arcs -ftest-coverage) - default: FALSE")

if(MODFEM_TEST_COVERAGE)
  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    if(NOT MODFEM_USE_STATIC)
      message(WARNING " -- Coverage testing does not work properly with shared libraries; Using static compilation")
      set(MODFEM_USE_STATIC ON)
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
    set(CMAKE_FORTRAN_FLAGS "${CMAKE_FORTRAN_FLAGS} -fprofile-arcs -ftest-coverage")
    set(LINK_FLAGS "${LINK_FLAGS}  -fprofile-arcs -ftest-coverage -static-libstdc++")
  else()
    message(WARNING " -- 'MODFEM_TEST_COVERAGE' is switched on while 'CMAKE_BUILD_TYPE' is not in 'Debug' mode; 'MODFEM_TEST_COVERAGE' is switched to false")
    set(MODFEM_TEST_COVERAGE FALSE)
  endif()
endif()

# -- MODFEM_USE_PROFILING --
# --------------------------
set(MODFEM_USE_PROFILING_DEFAULT FALSE)
set(MODFEM_USE_PROFILING ${MODFEM_USE_PROFILING_DEFAULT} CACHE BOOL "Enable profiling info (-pg) - default: ${MODFEM_USE_PROFILING_DEFAULT}")

if(MODFEM_USE_PROFILING)
  set(COMPILE_FLAGS "${COMPILE_FLAGS} -p")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -p")
else()
  if(DEFINED COMPILE_FLAGS)
    STRING(REPLACE "-p" "" COMPILE_FLAGS ${COMPILE_FLAGS})
  endif()
  if(DEFINED CMAKE_EXE_LINER_FLAGS)
    STRING(REPLACE "-p" "" CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS})
  endif()
endif()

# -- MODFEM_USE_STATIC --
# -----------------------
if(MODFEM_USE_STATIC)
  if(UNIX)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so")
    message("-- Searching libraries will be done for static linking.")
  endif()
  set(BUILD_SHARED_LIBS FALSE)
else()
  if(UNIX)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a")
    message("-- Searching libraries will be done for dynamic linking.")
  endif()
  set(BUILD_SHARED_LIBS TRUE)
endif()

# -- OpenMP SUPPORT --
# --------------------

if(CMAKE_C_COMPILER_ID STREQUAL "Clang")
    #set(CMAKE_C_FLAGS "-fopenmp=libomp ${CMAKE_C_FLAGS}")
    INCLUDE_DIRECTORIES("/usr/lib/gcc/x86_64-linux-gnu/5/include")
    LINK_DIRECTORIES("/usr/lib/gcc/x86_64-linux-gnu/5/")
    SET(CMAKE_EXE_LINKER_FLAGS "/usr/lib/gcc/x86_64-linux-gnu/5/libgomp.so ${CMAKE_EXE_LINKER_FLAGS}")
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    #set(CMAKE_CXX_FLAGS "-fopenmp=libomp ${CMAKE_CXX_FLAGS}")
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
  
  # Vienna OpenMP module
  if(MODFEM_ITER_SOLVER_MODULE STREQUAL "sil_viennacl_crs" OR MODFEM_DIRECT_SOLVER_MODULE STREQUAL "sil_viennacl_crs" OR MODFEM_MKB_DIRECT_SOLVER_MODULE STREQUAL "VIENNACL")
    add_definitions(-DVIENNACL_WITH_OPENMP)
  endif()


  
endif()

# -- OpenCL SUPPORT --
# --------------------
if(MODFEM_ACCEL)
  #Add OpenCL support
  if(${MODFEM_ACCEL} STREQUAL "opencl" OR (CREATE_MOD_FEM_VIEWER) OR (OPENCL_USE_LIBRARY_ONLY))
    #OpenCL found flag
    set(OpenCL_Inc_Found FALSE)
    set(OpenCL_Dir_Found FALSE)
    add_definitions(-DVIENNACL_WITH_OPENCL)
    
    #Find include dirs
    foreach(opencl_inc ${OPENCL_INCLUDE_DIRS})
      #C and C++ language headers
      if(EXISTS ${opencl_inc}/CL/cl.h)
	#ModFEM viewer c++ header
	if(CREATE_MOD_FEM_VIEWER)
	  if(NOT EXISTS ${opencl_inc}/CL/cl.hpp)
	    message(FATAL_ERROR "Can not build ModFEMViewer without 'cl.hpp' header") 
	  endif()
	endif()

	#OpenCL C header check
	set(OPENCL_INCLUDE_DIRS "${opencl_inc}")
	SET(OPENCL_HAS_CPP_BINDINGS TRUE)
	set(OpenCL_Inc_Found TRUE)
	break()
      elseif(EXISTS ${opencl_inc}/OpenCL/cl.h)
	#ModFEM viewer c++ header
	if(CREATE_MOD_FEM_VIEWER)
	  if(NOT EXISTS ${opencl_inc}/OpenCL/cl.hpp)
	    message(FATAL_ERROR "Can not build ModFEMViewer without 'cl.hpp' header") 
	  endif()
	endif()

	#OpenCL C header check
	set(OPENCL_INCLUDE_DIRS "${opencl_inc}")
	SET(OPENCL_HAS_CPP_BINDINGS TRUE)
	set(OpenCL_Inc_Found TRUE)
	break()
      endif()
    endforeach(opencl_inc)

    #Find library
    foreach(opencl_lib ${OPENCL_LIBRARY_DIRS})
      if(WIN32)
	if(EXISTS ${opencl_lib}/OpenCL.lib)
	  set(OPENCL_LIBRARIES ${opencl_lib}/OpenCL.lib)
	  set(OpenCL_Dir_Found TRUE)
	  break()
	endif()
      else()
	if(EXISTS ${opencl_lib}/libOpenCL.so)
	  #set(OPENCL_LIBRARIES ${opencl_lib}/libOpenCL.so)
    set(OPENCL_LIBRARIES "-L${opencl_lib}/ -lOpenCL")
	  set(OpenCL_Dir_Found TRUE)
	  break()
	endif()
      endif()
    endforeach(opencl_lib)

    #Error message
    if(NOT (OpenCL_Inc_Found AND OpenCL_Dir_Found))   
      message(FATAL_ERROR "OpenCL not found in this system (check your cmake file)!!!")
    endif()
    
    include_directories("${OPENCL_INCLUDE_DIRS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OPENCL_LIBRARIES}")
  endif()
else()
  message(FATAL_ERROR "MODFEM_ACCEL not defined. Allowed values: none, opencl, cuda.")
endif()
 
# -----------------------------------------------------------------------------
# ----------------------------- Third-party libraries -------------------------
# -----------------------------------------------------------------------------

include(FindPackageHandleStandardArgs)  # Used in local find_package

# -----------------------------------------------------------------------------
# ------------------------- Functions -----------------------------------------
# -----------------------------------------------------------------------------
function(collect_src_from_include_dirs INCLUDE_DIRS)
  foreach(folder ${${INCLUDE_DIRS}})
    file(
      GLOB
      source_files_one_dir
      ${folder}/*.f
      ${folder}/*.F
#      ${folder}/*.h
#      ${folder}/*.hh
#      ${folder}/*.hxx
#      ${folder}/*.hpp
      ${folder}/*.cc
      ${folder}/*.cxx
      ${folder}/*.cpp
      ${folder}/*.c
    )
    list(APPEND SOURCE_FILES_LIST ${source_files_one_dir})
  endforeach()
  set(SOURCE_FILES ${SOURCE_FILES_LIST} PARENT_SCOPE)
if(MSVC)  
  SET_SOURCE_FILES_PROPERTIES(${SOURCE_FILES_LIST} PROPERTIES LANGUAGE CXX)
endif()
endfunction()

# -----------------------------------------------------------------------------
# TESTS functions
# -----------------------------------------------------------------------------
# problem_name = name of the problem (without folder prefix: e.g. conv_diff etc.)
# ap = std | dg
# mm = prism|hybrid|remesh
# test_ver = Run|L2|H1|ZZ|Diff|Ref etc..
function(add_single_mf_test n_proc problem_name mm ap test_name test_ver)

#message("Kazik_DEBUG: Setting make test: ${ARGV6}")
    SET(config_file "")
    if(${ARGC} EQUAL 7) 
       SET(config_file "${ARGV6}")
    endif()

    if( (n_proc GREATER 0) AND (MODFEM_MPI STREQUAL "nompi") )
        RETURN() # parallel tests are ignored in not parallel env.
    elseif(n_proc EQUAL 0 AND MODFEM_MPI STREQUAL "mpi")
        RETURN() # sequential tests are off in parallel mode
    endif()

    SET(test_exacutable MOD_FEM_${problem_name}_${mm}_${ap}${MF_POSTFIX})
    
    if(${ARGC} EQUAL 7)
      SET(full_test_name ${problem_name}_${mm}_${ap}_${test_name}_${test_ver}_${config_file} )
    else()
      SET(full_test_name ${problem_name}_${mm}_${ap}_${test_name}_${test_ver} )
    endif()
    
    if(NOT n_proc EQUAL 0) # 0 = sequential code = no parallelism
        SET(full_test_name "${full_test_name}_mpi_${n_proc}")
        message("CMAKE_DEBUG: Setting mpi test: ${full_test_name}")
    else()
        message("CMAKE_DEBUG: Setting test: ${full_test_name}")
    endif()

if(${test_ver} STREQUAL "Run")
    add_test(${full_test_name} ${CTEST_DIR}/${problem_name}/${problem_name}_${test_ver}.py
             ${test_exacutable}
             ${CTEST_DIR}/${problem_name}/${test_name}/
             ${CTEST_BUILD_DIR}
             ${KKK}
             ${n_proc}
             ${config_file})
else()
    add_test(${full_test_name} ${CTEST_DIR}/${problem_name}/${problem_name}_${test_ver}.py
             ${CTEST_DIR}/${problem_name}/${test_name}/)
endif()
    set (passRegex "^Test passed")
    set (failRegex "Fail")
    set_property(TEST ${full_test_name} PROPERTY PASS_REGULAR_EXPRESSION "${passRegex}")
    set_property(TEST ${full_test_name} PROPERTY FAIL_REGULAR_EXPRESSION "${failRegex}")
if(CMAKE_BUILD_TYPE STREQUAL "Debug") # need more time in Debug mode
    set_tests_properties ( ${full_test_name} PROPERTIES TIMEOUT 600)
else()  # standard Timeout - if tests exceed it, that means something changed badly :/
    set_tests_properties ( ${full_test_name} PROPERTIES TIMEOUT 300)
endif()

    if(${test_ver} STREQUAL "Run")
        set_tests_properties(${full_test_name} PROPERTIES DEPENDS ${problem_name}_Make)
    elseif(${test_ver} STREQUAL "L2"
        OR ${test_ver} STREQUAL "H1"
        OR ${test_ver} STREQUAL "ZZ"
        OR ${test_ver} MATCHES "Diff.*"
        OR ${test_ver} STREQUAL "Run")
        set_tests_properties(${full_test_name} PROPERTIES DEPENDS ${problem_name}_${mm}_${ap}_Run)
    elseif(${test_ver} STREQUAL "Ref")
        set_tests_properties(${problem_name}_Make PROPERTIES DEPENDS ${full_test_name})
    endif()

endfunction()
# -----------------------------------------------------------------------------
# problem_name = name of the problem (without folder prefix: e.g. conv_diff etc.)
# ap = std | dg
# mm = prism|hybrid|remesh
# test_versions = list of all required tests, eq.: "Run;H1;L2"
# ARGV6 = (optional) config file name
function(add_mf_tests_mpi n_proc problem_name mm ap test_name test_versions)

    #message("problem name = ${problem_name}")
    if(NOT ${problem_name}_Make_test)
        message("CMAKE_DEBUG: Setting make test: ${problem_name}_Make")

        SET(${problem_name}_Make_test TRUE CACHE BOOL "${problem_name}_Make_test")

        add_test(${problem_name}_Make ${CTEST_DIR}/build.sh ${CTEST_BUILD_DIR}/pdd_${problem_name} ${ENV_PLATFORM} ${KKK} )

    endif()

    foreach(it_test_ver IN LISTS test_versions)
      if(${ARGC} EQUAL 6)
        add_single_mf_test(${n_proc} ${problem_name} ${mm} ${ap} ${test_name} ${it_test_ver})
      else()
        add_single_mf_test(${n_proc} ${problem_name} ${mm} ${ap} ${test_name} ${it_test_ver} ${ARGV6})
      endif()
    endforeach()

endfunction()
# ---------------
# problem_name = name of the problem (without folder prefix: e.g. conv_diff etc.)
# ap = std | dg
# mm = prism|hybrid|remesh
# test_versions = list of all required tests, eq.: "Run;H1;L2"
# ARGV5 = (optional) config file name
function(add_mf_tests problem_name mm ap test_name test_versions)

    # for sequential (no parallelism) test n_proc is forced to 0
  if(${ARGC} EQUAL 5)
    add_mf_tests_mpi("0" "${problem_name}" "${mm}" "${ap}" "${test_name}" "${test_versions}")
  else()
    add_mf_tests_mpi("0" "${problem_name}" "${mm}" "${ap}" "${test_name}" "${test_versions}" "${ARGV5}")
  endif()

endfunction()
# -----------------------------------------------------------------------------

# -- libconfig --
# ---------------

find_package(Libconfig REQUIRED)

# utd_util
set(UTD_UTIL_EXTERNAL_INCLUDE_DIRS ${UTD_UTIL_EXTERNAL_INCLUDE_DIRS} ${LIBCONFIG_INCLUDE_DIRS}) 
list(APPEND UTD_UTIL_EXTERNAL_LIBRARIES ${LIBCONFIG_LIBRARIES}) # utd_util

# pdd_conv_diff
if( CREATE_MOD_FEM_HEAT_PRISM_STD OR
    CREATE_MOD_FEM_HEAT_PRISM_STD_QUAD OR
    CREATE_MOD_FEM_HEAT_PRISM2D_STD OR
    CREATE_MOD_FEM_HEAT_PRISM2D_STD_QUAD OR
    CREATE_MOD_FEM_HEAT_HYBRID_STD OR
    CREATE_MOD_FEM_HEAT_HYBRID_STD_QUAD
    )
  set(PDD_HEAT_EXTERNAL_INCLUDE_DIRS ${PDD_HEAT_EXTERNAL_INCLUDE_DIRS} ${LIBCONFIG_INCLUDE_DIRS})
  list(APPEND PDD_HEAT_EXTERNAL_LIBRARIES ${LIBCONFIG_LIBRARIES})
endif()

# pdd_plast_flow
if( CREATE_MOD_FEM_PLAST_FLOW_PRISM_STD OR
    CREATE_MOD_FEM_PLAST_FLOW_HYBRID_STD OR
    CREATE_MOD_FEM_PLAST_FLOW_PRISM2D_STD OR
    CREATE_MOD_FEM_PLAST_FLOW_REMESH_STD
    )
  set(PDD_PLAST_FLOW_EXTERNAL_INCLUDE_DIRS ${PDD_PLAST_FLOW_EXTERNAL_INCLUDE_DIRS} ${LIBCONFIG_INCLUDE_DIRS})
  list(APPEND PDD_PLAST_FLOW_EXTERNAL_LIBRARIES ${LIBCONFIG_LIBRARIES})
endif()

# pdd_plast_flow_heat
if( CREATE_MOD_FEM_PLAST_FLOW_HEAT_PRISM_STD OR
    CREATE_MOD_FEM_PLAST_FLOW_HEAT_HYBRID_STD OR
    CREATE_MOD_FEM_PLAST_FLOW_HEAT_PRISM2D_STD OR
    CREATE_MOD_FEM_PLAST_FLOW_HEAT_REMESH_STD
    )
  set(PDD_PLAST_FLOW_HEAT_EXTERNAL_INCLUDE_DIRS ${PDD_PLAST_FLOW_HEAT_EXTERNAL_INCLUDE_DIRS} ${LIBCONFIG_INCLUDE_DIRS})
  list(APPEND PDD_PLAST_FLOW_HEAT_EXTERNAL_LIBRARIES ${LIBCONFIG_LIBRARIES})
endif()

# pdd_forming
if( CREATE_MOD_FEM_FORMING_PRISM_STD OR
    CREATE_MOD_FEM_FORMING_PRISM_STD_QUAD OR
    CREATE_MOD_FEM_FORMING_HYBRID_STD OR
    CREATE_MOD_FEM_FORMING_HYBRID_STD_QUAD OR
    CREATE_MOD_FEM_FORMING_PRISM2D_STD OR
    CREATE_MOD_FEM_FORMING_PRISM2D_STD_QUAD OR
    CREATE_MOD_FEM_FORMING_REMESH_STD OR
    CREATE_MOD_FEM_FORMING_REMESH_STD_QUAD
    )
  set(PDD_FORMING_EXTERNAL_INCLUDE_DIRS ${PDD_FORMING_EXTERNAL_INCLUDE_DIRS} ${LIBCONFIG_INCLUDE_DIRS})
  list(APPEND PDD_FORMING_EXTERNAL_LIBRARIES ${LIBCONFIG_LIBRARIES})
endif()

# pdd_forming_heat
if( CREATE_MOD_FEM_FORMING_HEAT_PRISM_STD OR
    CREATE_MOD_FEM_FORMING_HEAT_HYBRID_STD OR
    CREATE_MOD_FEM_FORMING_HEAT_PRISM2D_STD OR
    CREATE_MOD_FEM_FORMING_HEAT_REMESH_STD
    )
  set(PDD_FORMING_HEAT_EXTERNAL_INCLUDE_DIRS ${PDD_FORMIG_HEAT_EXTERNAL_INCLUDE_DIRS} ${LIBCONFIG_INCLUDE_DIRS})
  list(APPEND PDD_FORMING_HEAT_EXTERNAL_LIBRARIES ${LIBCONFIG_LIBRARIES})
endif()

# pdd_forming_mixed
if( CREATE_MOD_FEM_FORMING_PRISM_STD_MIXED OR
    CREATE_MOD_FEM_FORMING_HYBRID_STD_MIXED OR
    CREATE_MOD_FEM_FORMING_REMESH_STD_MIXED OR
    CREATE_MOD_FEM_FORMING_PRISM2D_STD_MIXED
    )
  set(PDD_FORMING_MIXED_EXTERNAL_INCLUDE_DIRS ${PDD_FORMING_MIXED_EXTERNAL_INCLUDE_DIRS} ${LIBCONFIG_INCLUDE_DIRS})
  list(APPEND PDD_FORMING_MIXED_EXTERNAL_LIBRARIES ${LIBCONFIG_LIBRARIES})
endif()

# pdd_ns_supg
if( CREATE_MOD_FEM_NS_SUPG_PRISM_STD OR
    CREATE_MOD_FEM_NS_SUPG_HYBRID_STD OR
    CREATE_MOD_FEM_NS_SUPG_PRISM2D_STD OR
    CREATE_MOD_FEM_NS_SUPG_REMESH_STD
    )
  set(PDD_NS_SUPG_EXTERNAL_INCLUDE_DIRS ${PDD_NS_SUPG_EXTERNAL_INCLUDE_DIRS} ${LIBCONFIG_INCLUDE_DIRS})
  list(APPEND PDD_NS_SUPG_EXTERNAL_LIBRARIES ${LIBCONFIG_LIBRARIES})
endif()

# pdd_ns_supg_heat
if( CREATE_MOD_FEM_NS_SUPG_HEAT_PRISM_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_HYBRID_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_PRISM2D_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_REMESH_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_ALE_REMESH_STD
    )
  set(PDD_NS_SUPG_HEAT_EXTERNAL_INCLUDE_DIRS ${PDD_NS_SUPG_HEAT_EXTERNAL_INCLUDE_DIRS} ${LIBCONFIG_INCLUDE_DIRS})
  list(APPEND PDD_NS_SUPG_HEAT_EXTERNAL_LIBRARIES ${LIBCONFIG_LIBRARIES})
endif()
    
# pdd_ns_supg_ale
if( CREATE_MOD_FEM_NS_SUPG_ALE_REMESH_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_ALE_REMESH_STD
    )
  set(PDD_NS_SUPG_ALE_EXTERNAL_INCLUDE_DIRS ${PDD_NS_SUPG_ALE_EXTERNAL_INCLUDE_DIRS} ${LIBCONFIG_INCLUDE_DIRS})
  list(APPEND PDD_NS_SUPG_ALE_EXTERNAL_LIBRARIES ${LIBCONFIG_LIBRARIES})
endif()

# pdd_ns_supg_heat_vof
if( CREATE_MOD_FEM_NS_SUPG_HEAT_VOF_PRISM_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_VOF_HYBRID_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_VOF_PRISM2D_STD
    )
    set(PDD_NS_SUPG_HEAT_VOF_EXTERNAL_INCLUDE_DIRS ${PDD_NS_SUPG_HEAT_VOF_EXTERNAL_INCLUDE_DIRS} ${LIBCONFIG_INCLUDE_DIRS})
    set(PDD_VOF_EXTERNAL_INCLUDE_DIRS ${PDD_VOF_EXTERNAL_INCLUDE_DIRS} ${LIBCONFIG_INCLUDE_DIRS})
  list(APPEND PDD_NS_SUPG_HEAT_VOF_EXTERNAL_LIBRARIES ${LIBCONFIG_LIBRARIES})
  list(APPEND PDD_VOF_EXTERNAL_LIBRARIES ${LIBCONFIG_LIBRARIES})

endif()

# ---------------
# -- voro++ --
# ---------------

  if( CREATE_MOD_FEM_NS_SUPG_HEAT_VOF_PRISM_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_VOF_HYBRID_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_VOF_PRISM2D_STD
    )
    find_package(VoroPP REQUIRED)
    #set EXTERNAL_INCLUDE_DIRS of all supertargets that use library
    set(PDD_NS_SUPG_EXTERNAL_INCLUDE_DIRS ${PDD_NS_SUPG_EXTERNAL_INCLUDE_DIRS} ${VOROPP_INCLUDE_DIRS})
    set(PDD_NS_SUPG_HEAT_EXTERNAL_INCLUDE_DIRS ${PDD_NS_SUPG_HEAT_EXTERNAL_INCLUDE_DIRS} ${VOROPP_INCLUDE_DIRS})
    set(PDD_HEAT_EXTERNAL_INCLUDE_DIRS ${PDD_HEAT_EXTERNAL_INCLUDE_DIRS} ${VOROPP_INCLUDE_DIRS})
    set(PDD_NS_SUPG_HEAT_VOF_EXTERNAL_INCLUDE_DIRS ${PDD_NS_SUPG_HEAT_VOF_EXTERNAL_INCLUDE_DIRS} ${VOROPP_INCLUDE_DIRS})
    set(PDD_VOF_EXTERNAL_INCLUDE_DIRS ${PDD_VOF_EXTERNAL_INCLUDE_DIRS} ${VOROPP_INCLUDE_DIRS})
    set(PDD_NS_SUPG_ALE_EXTERNAL_INCLUDE_DIRS ${PDD_NS_SUPG_ALE_EXTERNAL_INCLUDE_DIRS} ${VOROPP_INCLUDE_DIRS})
    #set EXTERNAL_LIBRARIES of all supertargets that use library
    list(APPEND PDD_NS_SUPG_EXTERNAL_LIBRARIES ${VOROPP_LIBRARIES})                 #pdd_ns_supg
    list(APPEND PDD_NS_SUPG_HEAT_EXTERNAL_LIBRARIES ${VOROPP_LIBRARIES})            #pdd_ns_supg_heat
    list(APPEND PDD_NS_SUPG_ALE_EXTERNAL_LIBRARIES ${VOROPP_LIBRARIES})             #pdd_ns_supg_ale
    list(APPEND PDD_HEAT_EXTERNAL_LIBRARIES ${VOROPP_LIBRARIES})                    #pdd_heat
    list(APPEND PDD_NS_SUPG_PRISM_STD_TURBULENT_EXTERNAL_LIBRARIES ${VOROPP_LIBRARIES}) #pdd_ns_supg_prism_std_turbulent
    list(APPEND PDD_NS_SUPG_HYBRID_STD_TURBULENT_EXTERNAL_LIBRARIES ${VOROPP_LIBRARIES}) #pdd_ns_supg_hybrid_std_turbulent
    list(APPEND PDD_NS_SUPG_HEAT_VOF_EXTERNAL_LIBRARIES ${VOROPP_LIBRARIES})        #pdd_ns_supg_heat_vof
    list(APPEND PDD_VOF_EXTERNAL_LIBRARIES ${VOROPP_LIBRARIES})		       #pdd_vof
    # -- pdd_vof ----
  endif()

# ------------------------------
# -- PARDISO && VIENNA SOLVER --
# ------------------------------

# PARDISO
set(SID_PARDISO_EXTERNAL_INCLUDE_DIRS ${SID_PARDISO_EXTERNAL_INCLUDE_DIRS})
list(APPEND SID_PARDISO_EXTERNAL_LIBRARIES)  #sid_pardiso

# VIENNA
set(SID_VIENNACL_CRS_EXTERNAL_INCLUDE_DIRS ${SID_VIENNACL_CRS_EXTERNAL_INCLUDE_DIRS} ${VIENNACL_INCLUDE_DIRS})
list(APPEND SID_VIENNACL_CRS_EXTERNAL_LIBRARIES ${VIENNACL_LIBRARIES})         #sid_viennacl_crs

  
  
# --------------------------------------------
# -- BLAS and Lapack (from MKL, ACML or standard) --
# --------------------------------------------
if((${MODFEM_BLASLAPACK} STREQUAL "GENERIC") AND (NOT WIN32))

  set(MOD_FEM_FORTRAN gfortran)
  set(CMAKE_EXE_LINKER_FLAGS "-pthread ${CMAKE_EXE_LINKER_FLAGS} -lgfortran")

else()
  set(MOD_FEM_FORTRAN)
endif()

find_package(Lapack REQUIRED) #lapack/blas always required!
#set EXTERNAL_INCLUDE_DIRS of all supertargets that use library
set(SID_VIENNACL_CRS_EXTERNAL_INCLUDE_DIRS ${SID_VIENNACL_CRS_EXTERNAL_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS})
set(SID_PARDISO_EXTERNAL_INCLUDE_DIRS ${SID_PARDISO_EXTERNAL_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS})
set(SID_MKB_EXTERNAL_INCLUDE_DIRS ${SID_MKB_EXTERNAL_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS})
set(LSD_MKB_EXTERNAL_INCLUDE_DIRS ${SID_MKB_EXTERNAL_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS})
set(SID_LAPACK_EXTERNAL_INCLUDE_DIRS ${SID_LAPACK_EXTERNAL_INCLUDE_DIRS} ${LAPACK_INCLUDE_DIRS})


#set EXTERNAL_LIBRARIES of all targets that use library
list(APPEND SID_VIENNACL_CRS_EXTERNAL_LIBRARIES ${LAPACK_LIBRARIES})         #sid_viennacl_crs
list(APPEND SID_PARDISO_EXTERNAL_LIBRARIES ${LAPACK_LIBRARIES})              #sid_pardiso
list(APPEND SID_MKB_EXTERNAL_LIBRARIES ${LAPACK_LIBRARIES})                  #sid_mkb
list(APPEND LSD_MKB_EXTERNAL_LIBRARIES ${LAPACK_LIBRARIES})                  #lsd_mkb
list(APPEND SID_LAPACK_EXTERNAL_LIBRARIES ${LAPACK_LIBRARIES})               #sid_lapack

# -----------
# -- Boost --
# -----------
  
set(Boost_USE_MULTITHREADED ON)
find_package(Boost REQUIRED COMPONENTS filesystem system regex )
find_library(Zlib_library z)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIR} ${ZLIB_INCLUDE_DIRS})
  link_directories(${Boost_LIBRARY_DIRS} ${ZLIB_LIBRARY_DIRS})

  set(UTD_UTIL_EXTERNAL_INCLUDE_DIRS ${UTD_UTIL_EXTERNAL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS})
  list(APPEND UTD_UTIL_EXTERNAL_LIBRARIES ${Boost_LIBRARIES} ${Zlib_library})
  add_definitions("-DHAS_BOOST")
endif()  

# -----------
# --  MPI  --
# -----------
if(MODFEM_NEW_MPI)
  add_definitions(-DMODFEM_NEW_MPI=true)
  #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMODFEM_NEW_MPI=true")
  #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMODFEM_NEW_MPI=true")
endif()
if(DEFINED MODFEM_MPI)
  if(${MODFEM_MPI} STREQUAL "mpi")
	find_package(MPI REQUIRED)
	include_directories(${MPI_INCLUDE_PATH})
	#include_directories(${METIS_INCLUDE_DIR})
	#include_directories(${PARMETIS_INCLUDE_DIR})
    set(MODFEM_MPI_CODE_DEFINE PARALLEL)
    add_definitions(-DPARALLEL)
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DPARALLEL")
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPARALLEL")
    if(MODFEM_NEW_MPI)
      message("MODFEM_NEW_MPI set to TRUE")
      find_package(METIS REQUIRED)
      find_package(PARMETIS REQUIRED)
      include_directories(${METIS_INCLUDE_DIR} ${PARMETIS_INCLUDE_DIR})
      link_directories(${METIS_LIBRARY_DIR} ${PARMETIS_LIBRARY_DIR})
      list(APPEND MODFEM_MPI_LIBS_STD mmpl_adapter appl_std_lin pcl_mpi_safe ${PARMETIS_LIBRARY} ${METIS_LIBRARY} ${MOD_FEM_FORTRAN})
      list(APPEND MODFEM_MPI_LIBS_DG mmpl_adapter appl_dg_prism pcl_mpi_safe ${PARMETIS_LIBRARY} ${METIS_LIBRARY} ${MOD_FEM_FORTRAN})
      add_subdirectory(ddd_parmetis)
      add_subdirectory(mmpd_adapter)
      add_subdirectory(pcd_mpi_safe)
      #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    else()
      message("MODFEM_NEW_MPI set to FALSE")
      list(APPEND MODFEM_MPI_LIBS_STD mmpl_prism appl_std_lin pcl_mpi ${MOD_FEM_FORTRAN})
      list(APPEND MODFEM_MPI_LIBS_DG mmpl_prism appl_dg_prism pcl_mpi ${MOD_FEM_FORTRAN})
      add_subdirectory(pcd_mpi)
    endif() 
          
  elseif(${MODFEM_MPI} STREQUAL "nompi")
    set(MODFEM_MPI_CODE_DEFINE NPARALLEL)
  else()
    message("MODFEM_MPI: ${MODFEM_MPI} ")
    message(FATAL_ERROR "MODFEM_MPI incorecctly defined. Allowed values: mpi, nompi.")
  endif()
else()
 message(FATAL_ERROR "MODFEM_MPI not defined. Allowed values: mpi, nompi.")
endif()

# -- check for wrong config file settings
if(MODFEM_ITER_SOLVER_MODULE STREQUAL "sil_mkb")
  if(NOT MODFEM_USE_STATIC)
    message(FATAL_ERROR "For sil_mkb (due to cyclic dependency) you need to set MODFEM_USE_STATIC to TRUE. Please edit your config file.")
  endif()
endif()

if(MODFEM_ITER_SOLVER_MODULE STREQUAL "sil_viennacl_crs")
    if(NOT MODFEM_USE_STATIC)
      message(FATAL_ERROR "For sil_viennacl_crs (due to cyclic dependency) you need to set MODFEM_USE_STATIC to TRUE. Please edit your config file.")
    endif()
endif()

if(MODFEM_DIRECT_SOLVER_MODULE STREQUAL "sil_pardiso")
  if(NOT MODFEM_BLASLAPACK STREQUAL "MKL")
    message(FATAL_ERROR "For sil_pardiso MODFEM_BLASLAPACK needs to be set to MKL. Please edit your config file.")
  endif()
endif()


# -- configure opencl
if(MODFEM_ACCEL)
  if(NOT ${MODFEM_ACCEL} STREQUAL "none")
    #Link OpenMP
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    
    #Set multithread flag    
    add_definitions(-DMULTITHREADED)

    #Set appropriate opencl flag
    if(${MODFEM_ACCEL} STREQUAL "opencl")
      set(accel "opencl")

      #enable gpu assembling
      if(OPENCL_GPU_ASSEMBLING)
	add_definitions(-DGPU_ASSEMBLING)
      endif()

      #gpu debug flag
      if(GPU_DEBUG)
	add_definitions(-DGPU_DEBUG)
      endif()
      
      #Check machine
      if(${OPENCL_MACHINE} STREQUAL "cpu")
	add_definitions(-DOPENCL_CPU)
      elseif(${OPENCL_MACHINE} STREQUAL "gpu")
	add_definitions(-DOPENCL_GPU)
      elseif(${OPENCL_MACHINE} STREQUAL "phi")
	add_definitions(-DOPENCL_PHI)
      elseif(${OPENCL_MACHINE} STREQUAL "hsa")
	add_definitions(-DOPENCL_HSA)
      else()
	message(FATAL_ERROR "OPENCL_MACHINE not defined. Allowed values is: cpu, gpu, phi, hsa.")
      endif()

      #Auto-Tuning
      if(OPENCL_AUTO_TUNING)
	add_definitions(-DTUNING=true)
      endif()
      
    #Set appropriate opencl flag
    elseif(${MODFEM_ACCEL} STREQUAL "cuda")
      #add_subdirectory(tmd_cuda)
      message(FATAL_ERROR "MODFEM_ACCEL: cuda not yet implemented.")
    else()
      set(accel "openmp")
    endif()
  else()
    set(accel "none")

    #Link OpenMP for modfem viewer
    if(CREATE_MOD_FEM_VIEWER)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    endif()
    
  endif()
else()
  message(FATAL_ERROR "MODFEM_ACCEL not defined. Allowed values: none, opencl, cuda.")
endif()


########################################################
# ----------------------- TARGETS -------------------- #
########################################################

#######################################
# -- APPROX, MESH AND UTIL TARGETS -- #
#######################################

# -- utd_util ----
add_subdirectory(utd_util) 
# -- utd_system ----
add_subdirectory(utd_system)
# -- mmd_prism ------
add_subdirectory(mmd_prism)
# -- mmd_t4_prism ------
add_subdirectory(mmd_t4_prism)
# -- mmd_prism_2D ----
add_subdirectory(mmd_prism_2D)
# -- mmd_remesh ------
add_subdirectory(mmd_remesh)
# -- mmd_fcpm
# add_subdirectory(mmd_fcpm)

# -- apd_std_quad ----
add_subdirectory(apd_std_quad)
# -- apd_std_lin ----
add_subdirectory(apd_std_lin)

# -- apd_dg_prism ----
add_subdirectory(apd_dg_prism)

if( CREATE_MOD_FEM_NS_SUPG_HEAT_VOF_PRISM_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_VOF_HYBRID_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_VOF_PRISM2D_STD
    )
    add_subdirectory(pdd_vof)
endif()

# -- mod_fem_viewer ---
list(APPEND MOD_FEM_TOOLS_LIST "CREATE_MOD_FEM_VIEWER")

if(CREATE_MOD_FEM_VIEWER)
  if(MOD_FEM_VIEWER_USE_AS_LIB)
    set(MOD_FEM_VIEWER_LIB "mod_fem_viewer")
    add_definitions(-DMOD_FEM_VIEWER)
  else()
    set(MOD_FEM_VIEWER_APP "FemViewer")
  endif()

  # Support for C++ 11 and link required libraries
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  if ( ${CMAKE_CXX_COMPILER_ID} STREQUAL GNU )
    set(LINK_FLAGS "${LINK_FLAGS} -lX11 -lXext -lXxf86vm -lXmu -lXi -lXrandr -lGLU -lGL")
  endif()
  
  add_subdirectory(mod_fem_viewer)
  # add mod_fem_viewer to targets
  list(APPEND PDD_CONV_DIFF_EXTERNAL_LIBRARIES ${MOD_FEM_VIEWER_LIB})            #pdd_conv_diff
  list(APPEND PDD_PLAST_FLOW_EXTERNAL_LIBRARIES  ${MOD_FEM_VIEWER_LIB})	 	 #pdd_plast_flow
  list(APPEND PDD_PLAST_FLOW_HEAT_EXTERNAL_LIBRARIES  ${MOD_FEM_VIEWER_LIB})	 #pdd_plast_flow_heat
  list(APPEND PDD_FORMING_EXTERNAL_LIBRARIES  ${MOD_FEM_VIEWER_LIB})	 	 #pdd_forming
  list(APPEND PDD_FORMING_HEAT_EXTERNAL_LIBRARIES  ${MOD_FEM_VIEWER_LIB})	 #pdd_forming_heat
  list(APPEND PDD_FORMING_MIXED_EXTERNAL_LIBRARIES ${MOD_FEM_VIEWER_LIB})       #pdd_forming_mixed
  list(APPEND PDD_NS_SUPG_EXTERNAL_LIBRARIES  ${MOD_FEM_VIEWER_LIB})	 	 #pdd_ns_supg
  list(APPEND PDD_HEAT_EXTERNAL_LIBRARIES  ${MOD_FEM_VIEWER_LIB})	 	 #pdd_heat
  list(APPEND PDD_NS_SUPG_HEAT_EXTERNAL_LIBRARIES  ${MOD_FEM_VIEWER_LIB})	 #pdd_ns_supg_heat
  
endif()

# -- modfem mpi ---
if(${MODFEM_MPI} STREQUAL "mpi")
  #add_subdirectory(mmpd_adapter)
  add_subdirectory(mmpd_prism)
  add_subdirectory(appd_std_lin)
  add_subdirectory(appd_dg_prism)
 # add_subdirectory(pcd_mpi) #this was done above
endif()

# -- modfem openmp && opencl --- 
if(MODFEM_ACCEL)
  if(${MODFEM_ACCEL} STREQUAL "opencl")
    # -- tmd_ocl ----
    add_subdirectory(tmd_opencl)
  elseif(${MODFEM_ACCEL} STREQUAL "cuda")
    #add_subdirectory(tmd_cuda)
    message(FATAL_ERROR "MODFEM_ACCEL: cuda not yet implemented.")
  endif()
else()
  message(FATAL_ERROR "MODFEM_ACCEL not defined. Allowed values: none, opencl, cuda.")
endif()


########################
# -- SOLVER TARGETS -- #
########################

if(MODFEM_DIRECT_SOLVER_MODULE STREQUAL "sil_pardiso")
  # -- sid_pardiso ----
  add_subdirectory(sid_pardiso)
elseif(MODFEM_DIRECT_SOLVER_MODULE STREQUAL "sil_lapack")
  # -- sid_lapack ----
  add_subdirectory(sid_lapack)
else()
  message(FATAL_ERROR "MODFEM_DIRECT_SOLVER_MODULE is wrong. Please check your configuration file.")
endif()

if(MODFEM_ITER_SOLVER_MODULE STREQUAL "sil_mkb")
  if(MODFEM_USE_PETSC)
  	add_subdirectory(amg_mkb)
  endif()
  add_subdirectory(lsd_mkb)
  add_subdirectory(sid_mkb)

elseif(MODFEM_ITER_SOLVER_MODULE STREQUAL "sil_viennacl_crs")
  add_subdirectory(sid_viennacl_crs)
  
else()
  message(FATAL_ERROR "MODFEM_ITER_SOLVER_MODULE is wrong. Please check your configuration file.")
endif()

#####################
# -- PDD TARGETS -- #
#####################
# -- (which of them are created is decided in their CMake files)
if(     CREATE_MOD_FEM_CONV_DIFF_PRISM_STD OR
    CREATE_MOD_FEM_CONV_DIFF_PRISM_STD_QUAD OR
    CREATE_MOD_FEM_CONV_DIFF_HYBRID_STD OR
    CREATE_MOD_FEM_CONV_DIFF_HYBRID_STD_QUAD OR
    CREATE_MOD_FEM_CONV_DIFF_PRISM2D_STD OR
    CREATE_MOD_FEM_CONV_DIFF_PRISM2D_STD_QUAD OR
    CREATE_MOD_FEM_CONV_DIFF_PRISM_DG OR
    CREATE_MOD_FEM_CONV_DIFF_HYBRID_DG OR
    CREATE_MOD_FEM_CONV_DIFF_PRISM2D_DG)
    add_subdirectory(pdd_conv_diff)
endif()

if( CREATE_MOD_FEM_PLAST_FLOW_PRISM_STD OR
    CREATE_MOD_FEM_PLAST_FLOW_HYBRID_STD OR
    CREATE_MOD_FEM_PLAST_FLOW_PRISM2D_STD OR
    CREATE_MOD_FEM_PLAST_FLOW_REMESH_STD)
# currently pdd_plast_flow_supg assumed as default implementation for plast_flow
    add_subdirectory(pdd_plast_flow_supg)
endif()

if( CREATE_MOD_FEM_FORMING_PRISM_STD OR
    CREATE_MOD_FEM_FORMING_PRISM_STD_QUAD OR
    CREATE_MOD_FEM_FORMING_HYBRID_STD OR
    CREATE_MOD_FEM_FORMING_HYBRID_STD_QUAD OR
    CREATE_MOD_FEM_FORMING_PRISM2D_STD OR
    CREATE_MOD_FEM_FORMING_PRISM2D_STD_QUAD OR
    CREATE_MOD_FEM_FORMING_REMESH_STD OR
    CREATE_MOD_FEM_FORMING_REMESH_STD_QUAD)
    add_subdirectory(pdd_forming)
endif()

if( CREATE_MOD_FEM_FORMING_PRISM_STD_MIXED OR
    CREATE_MOD_FEM_FORMING_HYBRID_STD_MIXED OR
    CREATE_MOD_FEM_FORMING_PRISM2D_STD_MIXED OR
    CREATE_MOD_FEM_FORMING_REMESH_STD_MIXED)
    add_subdirectory(pdd_forming_mixed)
endif()

if( CREATE_MOD_FEM_NS_SUPG_PRISM_STD OR
    CREATE_MOD_FEM_NS_SUPG_HYBRID_STD OR
    CREATE_MOD_FEM_NS_SUPG_PRISM2D_STD OR
    CREATE_MOD_FEM_NS_SUPG_REMESH_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_PRISM_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_HYBRID_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_PRISM2D_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_REMESH_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_ALE_REMESH_STD)
    add_subdirectory(pdd_ns_supg)
endif()

if( CREATE_MOD_FEM_HEAT_PRISM_STD OR
    CREATE_MOD_FEM_HEAT_PRISM_STD_QUAD OR
    CREATE_MOD_FEM_HEAT_PRISM2D_STD OR
    CREATE_MOD_FEM_HEAT_PRISM2D_STD_QUAD OR
    CREATE_MOD_FEM_HEAT_HYBRID_STD OR
    CREATE_MOD_FEM_HEAT_HYBRID_STD_QUAD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_PRISM_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_HYBRID_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_PRISM2D_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_REMESH_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_ALE_REMESH_STD)
    add_subdirectory(pdd_heat)
endif()


# OBSOLETE
# add_subdirectory(pdd_elast)

if( CREATE_MOD_FEM_NS_SUPG_HEAT_ALE_REMESH_STD)
add_subdirectory(pdd_ns_supg_ale)
endif()

if(     CREATE_MOD_FEM_NS_SUPG_HEAT_PRISM_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_HYBRID_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_PRISM2D_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_REMESH_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_ALE_REMESH_STD)
    add_subdirectory(pdd_ns_supg_heat)
endif()

if(    CREATE_MOD_FEM_NS_SUPG_HEAT_VOF_PRISM_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_VOF_HYBRID_STD OR
    CREATE_MOD_FEM_NS_SUPG_HEAT_VOF_PRISM2D_STD)
    add_subdirectory(pdd_ns_supg_heat_vof)
endif()

# -- summary of settings
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_Fortran_FLAGS_ALL "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}")
  set(CMAKE_CXX_FLAGS_ALL "${CMAKE_CXX_FLAGS_ALL} ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
  set(CMAKE_C_FLAGS_ALL "${CMAKE_C_FLAGS_ALL} ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
  set(CMAKE_EXE_LINKER_FLAGS_ALL "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELEASE}")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  set(CMAKE_Fortran_FLAGS_ALL "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELWITHDEBINFO}")
  set(CMAKE_CXX_FLAGS_ALL "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELWITHDEBINFO}")
  set(CMAKE_C_FLAGS_ALL "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELWITHDEBINFO}")
  set(CMAKE_EXE_LINKER_FLAGS_ALL "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_RELWITHDEBINFO}")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_Fortran_FLAGS_ALL "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
  set(CMAKE_CXX_FLAGS_ALL "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
  set(CMAKE_C_FLAGS_ALL "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
  set(CMAKE_EXE_LINKER_FLAGS_ALL "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
else()
  set(CMAKE_Fortran_FLAGS_ALL "Unknown")
  set(CMAKE_CXX_FLAGS_ALL "Unknown")
  set(CMAKE_C_FLAGS_ALL "Unknown")
  set(CMAKE_EXE_LINKER_FLAGS_ALL "Unknown")
endif()

#get_property(SID_PARDISO_DEFINITIONS DIRECTORY sid_pardiso PROPERTY COMPILE_DEFINITIONS)

# -- list of definitions
get_directory_property( MODFEM_DEFINITIONS DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )

if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(ColourBold  "${Esc}[1m")
  set(Red         "${Esc}[31m")
  set(Green       "${Esc}[32m")
  set(Yellow      "${Esc}[33m")
  set(Blue        "${Esc}[34m")
  set(Magenta     "${Esc}[35m")
  set(Cyan        "${Esc}[36m")
  set(White       "${Esc}[37m")
  set(BoldRed     "${Esc}[1;31m")
  set(BoldGreen   "${Esc}[1;32m")
  set(BoldYellow  "${Esc}[1;33m")
  set(BoldBlue    "${Esc}[1;34m")
  set(BoldMagenta "${Esc}[1;35m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
endif()

message("")
message("=> Config file in use: ${ENV_PLATFORM}.cmake")
message("=> Iterative solver module in use: ${MODFEM_ITER_SOLVER_MODULE}")
message("=> Direct solver module in use: ${MODFEM_DIRECT_SOLVER_MODULE}")
message("=> BLAS/LAPACK from: ${MODFEM_BLASLAPACK}")
message("=> Compilation:")
message("   - CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")
message("   - MODFEM_USE_STATIC: ${MODFEM_USE_STATIC}")
message("   - Fortran:")
message("     - Compiler: ${CMAKE_Fortran_COMPILER}")
message("     - Flags: ${CMAKE_Fortran_FLAGS_ALL}")
message("   - C:")
message("     - Compiler: ${CMAKE_C_COMPILER}")
message("     - Flags: ${CMAKE_C_FLAGS_ALL}")
message("   - C++:")
message("     - Compiler: ${CMAKE_CXX_COMPILER}")
message("     - Flags: ${CMAKE_CXX_FLAGS_ALL}")
message("   - Link:")
message("     - Flags for exe: ${CMAKE_EXE_LINKER_FLAGS_ALL}")
message("     - Flags added manually: ${LINK_FLAGS}")
message("   - List of definitions:")
message("     - Definitions: ${MODFEM_DEFINITIONS}")  
message("   - MODFEM_MPI: ${MODFEM_MPI}")
message("   - MODFEM_ACCEL: ${MODFEM_ACCEL}")
if(${MODFEM_ACCEL} STREQUAL "opencl")
  message("     - OpenCL headers directory: ${OPENCL_INCLUDE_DIRS}")
  message("     - OpenCL libraries directories: ${OPENCL_LIBRARY_DIRS}")
  message("     - OpenCL library file: ${OPENCL_LIBRARIES}")
  if(${OPENCL_MACHINE} STREQUAL "cpu")
    message("     - OpenCL acceleration enabled for CPU (intel/amd)")
  elseif(${OPENCL_MACHINE} STREQUAL "gpu")
    message("     - OpenCL acceleration enabled for GPU (intel/nvidia/amd)")
  elseif(${OPENCL_MACHINE} STREQUAL "phi")
    message("     - OpenCL acceleration enabled for ACCELERATOR (intel xeon phi)")
  elseif(${OPENCL_MACHINE} STREQUAL "hsa")
    message("     - OpenCL acceleration enabled for AMD Heterogeneous System Architecture (amd)")
  endif()
  if(OPENCL_AUTO_TUNING STREQUAL "TRUE")
    message("     - OpenCL Auto-Tuning is ${Green}enabled${ColourReset}")
  else()
    message("     - OpenCL Auto-Tuning is ${Red}disabled${ColourReset}")
  endif()
endif()

#message("   - Definitions:")
#message("     - sid_pardiso: ${SID_PARDISO_DEFINITIONS}")

message("=> List of created auxiliary tools:")

list(LENGTH MOD_FEM_TOOLS_LIST number_of_tools)
math(EXPR number_of_tools "${number_of_tools} - 1")
foreach(iter RANGE ${number_of_tools})
  list(GET MOD_FEM_TOOLS_LIST ${iter} loc_name)
  if(${loc_name} STREQUAL "TRUE")
    message(STATUS "${Green}${loc_name}${ColourReset} flag is ${Green}${${loc_name}}${ColourReset}")
  elseif(${loc_name} STREQUAL "FALSE")
    message(STATUS "${Red}${loc_name}${ColourReset} flag is ${Red}${${loc_name}}${ColourReset}")
  else()
    message(STATUS "${Yellow}${loc_name}${ColourReset} flag is ${Yellow}UNDEFINED${ColourReset}")
  endif()
  unset(loc_name)
endforeach(iter)
unset(iter)
unset(number_of_tools)

message("=> List of created executable files:")

list(LENGTH MOD_FEM_EXEC_LIST number_of_exec)
math(EXPR number_of_exec "${number_of_exec} - 1")
foreach(iter RANGE ${number_of_exec})
  list(GET MOD_FEM_EXEC_LIST ${iter} loc_name)
  if(${loc_name} STREQUAL "TRUE")
    message(STATUS "${Green}${loc_name}${ColourReset} flag is ${Green}${${loc_name}}${ColourReset}")
  elseif(${loc_name} STREQUAL "FALSE")
    message(STATUS "${Red}${loc_name}${ColourReset} flag is ${Red}${${loc_name}}${ColourReset}")
  else()
    message(STATUS "${Yellow}${loc_name}${ColourReset} flag is ${Yellow}UNDEFINED${ColourReset}")
  endif()
  unset(loc_name)
endforeach(iter)
unset(iter)
unset(number_of_exec)

message("--------------------------------------------------------------------")
message("")

INCLUDE(CTest)

