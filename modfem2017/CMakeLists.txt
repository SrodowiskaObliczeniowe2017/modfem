# New CMAKE for ModFEM compilation it replace (config.sh) and add some new features
# Before use, please read short description:
# bin_cmake/README.txt
# 
# ----------------------------------------------
# History:
# 11.06.2015 Jan Bielanski, jbielan@agh.edu.pl, version 1.01 - support compilers with version postfix 
# 11.06.2015 Jan Bielanski, jbielan@agh.edu.pl, version 1.00 - added default MSVC compiler for windows  
# 09.06.2015 Jan Bielanski, jbielan@agh.edu.pl, version 0.90 - mode name fixes and many small changes
# 01.05.2015 Jan Bielanski, jbielan@agh.edu.pl, version 0.86/0.85 - bug fixes
# 30.04.2015 Jan Bielanski, jbielan@agh.edu.pl, version 0.80 - first version
#

# Launcher cmake script
cmake_minimum_required(VERSION 2.6)

#
# Initializing variables etc.
#
PROJECT(ModFEM_launcher)

#
# Create clean script and add clean command
#
file(WRITE ${CMAKE_BINARY_DIR}/clean-upper.cmake
"
set(cmake_generated \${CMAKE_BINARY_DIR}/CMakeCache.txt
                    \${CMAKE_BINARY_DIR}/bin_cmake/NOTFOUND
                    \${CMAKE_BINARY_DIR}/cmake_install.cmake  
                    \${CMAKE_BINARY_DIR}/Makefile
                    \${CMAKE_BINARY_DIR}/CMakeFiles
                    \${CMAKE_BINARY_DIR}/clean-upper.cmake
)
foreach(file \${cmake_generated})
  if(EXISTS \${file})
     file(REMOVE_RECURSE \${file})
  endif()
endforeach(file)
")

add_custom_target(clean-upper
  COMMAND ${CMAKE_BUILD_TOOL} clean
  COMMAND ${CMAKE_COMMAND} -P clean-upper.cmake
)

#
# Default configuration fo compiler and calculation mode
#
if(WIN32)
  set(DEFAULT_MF_CC "msvc")
  set(DEFAULT_MF_CXX "msvc")
else()
  set(DEFAULT_MF_CC "gcc")
  set(DEFAULT_MF_CXX "g++")
endif()
set(DEFAULT_MF_ACCEL "openmp") # Available flag [none,openmp,opencl]
set(DEFAULT_MF_MPI "nompi") # Available flag [nompi,mpi]
set(DEFAULT_MF_BUILD "Release") # Available flag [Release,RelWithDebInfo,Debug]

#
# Get list of platforms
#
SET(MF_ARCH_LIST)
FILE(GLOB  TMP_MF_ARCH_LIST  "src/cmake/Platforms/*.cmake")
foreach(a ${TMP_MF_ARCH_LIST})
  get_filename_component(TMP_MF_ARCH_FILENAME "${a}" NAME_WE)
  LIST(APPEND MF_ARCH_LIST ${TMP_MF_ARCH_FILENAME})
endforeach(a)

#
# ModFEM arch test function
#
SET(MF_ARCH_CONFIRMED -1)
function(CHECK_MF_ARCH ARCH_NAME)
  LIST(FIND MF_ARCH_LIST ${ARCH_NAME} FOUND)
  SET(MF_ARCH_CONFIRMED ${FOUND} PARENT_SCOPE)
  # Show available platforms if it was wrong
  if(${FOUND} EQUAL -1 OR ${ARCH_NAME} STREQUAL "REFERENCE_PLATFORM")
    SET(MF_ARCH_CONFIRMED "-1" PARENT_SCOPE)
    SET(MF_BAD_ARCH_MSG "MOD_FEM_ARCH is set to  \"${ARCH_NAME}\" but corresponding file ${ARCH_NAME}.cmake doesn't exist.")
    SET(MF_BAD_ARCH_MSG "${MF_BAD_ARCH_MSG} Create architecture file \"src/cmake/Platforms/${ARCH_NAME}.cmake\" exclude name \"REFERENCE_PLATFORM.cmake\" or use one of the following architectures:\n")
    foreach(a ${MF_ARCH_LIST})
      if(NOT ${a} STREQUAL "REFERENCE_PLATFORM")
	SET(MF_BAD_ARCH_MSG "${MF_BAD_ARCH_MSG}  -- ${a}\n")
      endif()
    endforeach(a)
    message(WARNING "${MF_BAD_ARCH_MSG}")
  endif()
endfunction(CHECK_MF_ARCH)

function(CHECK_WEAK_MF_ARCH ARCH_NAME)
  LIST(FIND MF_ARCH_LIST ${ARCH_NAME} FOUND)
  SET(MF_ARCH_CONFIRMED ${FOUND} PARENT_SCOPE)
endfunction(CHECK_WEAK_MF_ARCH)

#
# ModFEM configuration variables
#
set(MF_ARCH "" CACHE PATH "Architecture file name [see src/cmake/Platforms/]" FORCE)
if(DEFINED MF_CC)
  set(MF_CC "${MF_CC}")
else()
  set(MF_CC "" CACHE PATH "C Compiler name" FORCE)
endif()
if(DEFINED MF_CXX)
  set(MF_CXX "${MF_CXX}")
else()
  set(MF_CXX "" CACHE PATH "C++ Compiler name" FORCE)
endif()
if(DEFINED MF_ACCEL)
  set(MF_ACCEL "${MF_ACCEL}")
else()
  set(MF_ACCEL "" CACHE PATH "Hardware acceleration (w/multithreading) flag [none,openmp,opencl]" FORCE)
endif()
if(DEFINED MF_MPI)
  set(MF_MPI "${MF_MPI}")
else()
  set(MF_MPI "" CACHE PATH "Distributed memory support flag [nompi,mpi]" FORCE)
endif()
set(MF_BUILD "" CACHE PATH "Project build type [Release,RelWithDebInfo,Debug]" FORCE)
# Options for old make based compilation
# MF_DIR - hold ModFEM directory (path where is src subdirectory)
# MF_PROBLEM - hold information about problem which we will build


# MF_OPENCL_TASK - hold LAPLACE or TEST_SCALAR flag
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

#
# PART 1: Gathering the most important setup.
#

# 1.1 ANCIENT_MAKE MODE - use make and Makefile_explicit
# 1.2 SINGLE_PLATFORM MODE - use cmake
# 1.3 MULTI_PLATFORM MODE - use cmake

# ANCIENT_MAKE MODE
if(DEFINED ENV{MOD_FEM_ARCH} AND NOT DEFINED ENV{MOD_FEM_ARCH_CMAKE})
  message(STATUS "Selected ANCIENT_MAKE MODE use make & Mafile_explicit (found MOD_FEM_ARCH variable only).")
  
  #Set ModFEM arch
  set(MF_ARCH $ENV{MOD_FEM_ARCH})

  #Check directory
  if(NOT EXISTS "${MF_DIR}/src/Makefile_explicit")
    set(MF_DIR $ENV{MOD_FEM_DIR})
    if(NOT EXISTS "${MF_DIR}/src/Makefile_explicit")
      set(MF_DIR ${CMAKE_CURRENT_SOURCE_DIR})
      if(NOT EXISTS "${MF_DIR}/src/Makefile_explicit")
	message(FATAL_ERROR "Can not find Makefile_explicit")
      endif()
    endif()
  endif()
  set(MF_USER_MODE "ancient_make")

# SINGLE_PLATFORM MODE
elseif(DEFINED ENV{MOD_FEM_ARCH_CMAKE})
  message(STATUS "Selected SINGLE_PLATFORM MODE use cmake (found MOD_FEM_ARCH_CMAKE variable only).")
  
  #Set ModFEM arch
  set(MF_ARCH $ENV{MOD_FEM_ARCH_CMAKE})
  set(MF_USER_MODE "single_platform")

  #Check ModFEM arch
  CHECK_MF_ARCH(${MF_ARCH})

  if(NOT ${MF_ARCH_CONFIRMED} EQUAL -1)
    #
    # Attempts to detect variable values from cmake file name
    #

    # Temporary variables
    set(MF_TMP_ARCH "")
    set(MF_TMP1 "")
    set(MF_TMP2 "")
    set(MF_TMP3 "")

    #Set letters to low in platform name
    string(TOLOWER ${MF_ARCH} MF_TMP_ARCH)
    
    # MPI flag
    if(NOT MF_MPI)
      set(MF_MPI "false")
    else()
      if(NOT (${MF_MPI} STREQUAL "mpi" OR NOT ${MF_MPI} STREQUAL "nompi"))
	set(MF_MPI "false")
      endif()
    endif()

    #Try to read MPI flag
    if(${MF_MPI} STREQUAL "false")
      string(REGEX MATCH "[-_]mpi" MF_TMP1 ${MF_TMP_ARCH})
      string(REGEX MATCH "[-_]nompi" MF_TMP2 ${MF_TMP_ARCH})
      if(MF_TMP1)
	if(${MF_TMP1} STREQUAL "-mpi" OR ${MF_TMP1} STREQUAL "_mpi")
	  set(MF_MPI "mpi")
	endif()
      elseif(MF_TMP2)
	if(${MF_TMP2} STREQUAL "-nompi" OR ${MF_TMP2} STREQUAL "_nompi")
	  set(MF_MPI "nompi")
	endif()
      endif()
    endif()

    #Set default MPI if it not set
    if(${MF_MPI} STREQUAL "false")
      set(MF_MPI "${DEFAULT_MF_MPI}")
    endif()

    # Acceleration flag
    if(NOT MF_ACCEL)
      set(MF_ACCEL "false")
    else()
      if(NOT (${MF_ACCEL} STREQUAL "openmp" OR ${MF_ACCEL} STREQUAL "opencl"))
	if(NOT ${MF_ACCEL} STREQUAL "none")
	  set(MF_ACCEL "false")
	endif()
      endif()
    endif()

    #Try to read ACCELERATION flag
    string(REGEX MATCH "[-_]opencl" MF_TMP1 ${MF_TMP_ARCH})
    string(REGEX MATCH "[-_]openmp" MF_TMP2 ${MF_TMP_ARCH})
    string(REGEX MATCH "[-_]none" MF_TMP3 ${MF_TMP_ARCH})
    if(MF_TMP1)
      if(${MF_TMP1} STREQUAL "-opencl" OR ${MF_TMP1} STREQUAL "_opencl")
	set(MF_ACCEL "opencl")
      endif()
    elseif(MF_TMP2)
      if(${MF_TMP2} STREQUAL "-openmp" OR ${MF_TMP2} STREQUAL "_openmp")
	set(MF_ACCEL "openmp")
      endif()
    elseif(MF_TMP3)
      if(${MF_TMP3} STREQUAL "-none" OR ${MF_TMP3} STREQUAL "_none")
	set(MF_ACCEL "none")
      endif()
    endif()

    #Set default ACCELERATION if it not set
    if(${MF_ACCEL} STREQUAL "false")
      set(MF_ACCEL "${DEFAULT_MF_ACCEL}")
    endif()
  
    # Detect compiler
    if(NOT (MF_CC AND MF_CXX))
      if(WIN32)
	# - msvc
	unset(MF_TMP1)
	unset(MF_TMP2)
	string(REGEX MATCH "[-_]msvc" MF_TMP1 ${MF_TMP_ARCH})
	string(REGEX MATCH "[-_]msvc" MF_TMP2 ${MF_TMP_ARCH})
	if(MF_TMP1 AND MF_TMP2)
	  if((${MF_TMP1} STREQUAL "-msvc" OR ${MF_TMP1} STREQUAL "_msvc") AND (${MF_TMP2} STREQUAL "-msvc" OR ${MF_TMP2} STREQUAL "_msvc"))
	    #message("COMPILER: msvc/msvc")
	    set(MF_CC "msvc")
	    set(MF_CXX "msvc")
	  endif()
	endif()
	
      else()
	# - gcc and g++
	unset(MF_TMP1)
	unset(MF_TMP2)
	string(REGEX MATCH "[-_]gcc" MF_TMP1 ${MF_TMP_ARCH})
	string(REGEX MATCH "[-_]g[+][+]" MF_TMP2 ${MF_TMP_ARCH})
	if(MF_TMP1 AND MF_TMP2)
	  if((${MF_TMP1} STREQUAL "-gcc" OR ${MF_TMP1} STREQUAL "_gcc") AND (${MF_TMP2} STREQUAL "-g++" OR ${MF_TMP2} STREQUAL "_g++"))
	    #message("COMPILER: gcc/g++")
	    set(MF_CC "gcc")
	    set(MF_CXX "g++")
	  endif()
	endif()

	# - clang and clang++
	unset(MF_TMP1)
	unset(MF_TMP2)
	string(REGEX MATCH "[-_]clang" MF_TMP1 ${MF_TMP_ARCH})
	string(REGEX MATCH "[-_]clang[+][+]" MF_TMP2 ${MF_TMP_ARCH})
	if(MF_TMP1 AND MF_TMP2)
	  if((${MF_TMP1} STREQUAL "-clang" OR ${MF_TMP1} STREQUAL "_clang") AND (${MF_TMP2} STREQUAL "-clang++" OR ${MF_TMP2} STREQUAL "_clang++"))
	    #message("COMPILER: clang/clang++")
	    set(MF_CC "clang")
	    set(MF_CXX "clang++")
	  endif()
	endif()


	# - mpicc and mpic++
	if(${MF_MPI} STREQUAL "mpi")
	  unset(MF_TMP1)
	  unset(MF_TMP2)
	  string(REGEX MATCH "[-_]mpicc" MF_TMP1 ${MF_TMP_ARCH})
	  string(REGEX MATCH "[-_]mpic[+][+]" MF_TMP2 ${MF_TMP_ARCH})
	  if(NOT MF_TMP2)
	    string(REGEX MATCH "[-_]mpiCC" MF_TMP2 ${MF_TMP_ARCH})
	  endif()
	  if(NOT MF_TMP2)
	    string(REGEX MATCH "[-_]mpicxx" MF_TMP2 ${MF_TMP_ARCH})
	  endif()
	  
	  if(MF_TMP1 AND MF_TMP2)
	    if((${MF_TMP1} STREQUAL "-mpicc" OR ${MF_TMP1} STREQUAL "_mpicc") AND
		(${MF_TMP2} STREQUAL "-mpic++" OR ${MF_TMP2} STREQUAL "_mpic++" OR
		  ${MF_TMP2} STREQUAL "-mpiCC" OR ${MF_TMP2} STREQUAL "_mpiCC" OR
		  ${MF_TMP2} STREQUAL "-mpicxx" OR ${MF_TMP2} STREQUAL "_mpicxx"))
	      #message("COMPILER: gcc/g++")
	      set(MF_CC "mpicc")
	      set(MF_CXX "mpic++")
	    endif()
	  endif()
	endif()
      endif()
      
      # - icc and icpc
      unset(MF_TMP1)
      unset(MF_TMP2)
      string(REGEX MATCH "[-_]icc" MF_TMP1 ${MF_TMP_ARCH})
      string(REGEX MATCH "[-_]icpc" MF_TMP2 ${MF_TMP_ARCH})
      if(MF_TMP1 AND MF_TMP2)
	if((${MF_TMP1} STREQUAL "-icc" OR ${MF_TMP1} STREQUAL "_icc") AND (${MF_TMP2} STREQUAL "-icpc" OR ${MF_TMP2} STREQUAL "_icpc"))
	  #message("COMPILER: icc/icpc")
	  set(MF_CC "icc")
	  set(MF_CXX "icpc")
	endif()
      endif()

      # - mpiicc and mpiicpc
      if(${MF_MPI} STREQUAL "mpi")
	unset(MF_TMP1)
	unset(MF_TMP2)
	string(REGEX MATCH "[-_]mpiicc" MF_TMP1 ${MF_TMP_ARCH})
	string(REGEX MATCH "[-_]mpiicpc" MF_TMP2 ${MF_TMP_ARCH})
	if(MF_TMP1 AND MF_TMP2)
	  if((${MF_TMP1} STREQUAL "-mpiicc" OR ${MF_TMP1} STREQUAL "_mpiicc") AND (${MF_TMP2} STREQUAL "-mpiicpc" OR ${MF_TMP2} STREQUAL "_mpiicpc"))
	    #message("COMPILER: mpiicc/mpiicpc")
	    set(MF_CC "mpiicc")
	    set(MF_CXX "mpiicpc")
	  endif()
	endif()
      endif()  
    endif()

    # MPI Warning for unix systems
    if(NOT WIN32)
      if(${MF_MPI} STREQUAL "mpi")

	# Test C compiler
	set(C_MPI_COMPILER "ok")
	unset(MF_TMP1)

	if(NOT ${MF_CC} STREQUAL "mpicc" AND NOT ${MF_CC} STREQUAL "mpiicc")
	  string(REGEX MATCH "^.*-" MF_TMP1 ${MF_CC})
	  if(MF_TMP1)
	    if(NOT ${MF_TMP1} STREQUAL "mpicc-" AND NOT ${MF_TMP1} STREQUAL "mpiicc-")
	      set(C_MPI_COMPILER "fail")
	    endif()
	  else()
	    set(C_MPI_COMPILER "fail")
	  endif()
	endif()

	# Test C++ compiler
	set(CPP_MPI_COMPILER "ok")
	unset(MF_TMP1)

	if(NOT ${MF_CXX} STREQUAL "mpiCC" AND NOT ${MF_CXX} STREQUAL "mpicxx" AND NOT ${MF_CXX} STREQUAL "mpic++" AND NOT ${MF_CXX} STREQUAL "mpiicc")
	  string(REGEX MATCH "^.*-" MF_TMP1 ${MF_CXX})
	  if(MF_TMP1)
	    if(NOT ${MF_TMP1} STREQUAL "mpiCC-" AND NOT ${MF_TMP1} STREQUAL "mpicxx-" AND NOT ${MF_TMP1} STREQUAL "mpic++-" AND NOT ${MF_TMP1} STREQUAL "mpiicpc-")
	      set(CPP_MPI_COMPILER "fail")
	    endif()
	  else()
	    set(CPP_MPI_COMPILER "fail")
	  endif()
	endif()
	
	# Compiler info
	if(${C_MPI_COMPILER} STREQUAL "fail" OR ${CPP_MPI_COMPILER} STREQUAL "fail")
	  set(MSG_WARN "Compiler not support Message Passing Interface (MPI)\n")
	  set(MSG_WARN "${MSG_WARN}You should one of the following compilers:\n")
	  set(MSG_WARN "${MSG_WARN}-- GNU COMPILERS:   mpicc and mpic++ / mpicxx / mpiCC\n")
	  set(MSG_WARN "${MSG_WARN}-- INTEL COMPILERS: mpiicc and mpiicpc")
	  message(WARNING "${MSG_WARN}")
	  set(MF_CC "fail")
	  set(MF_CXX "fail")
	endif()
      endif()
    endif()

    #Read system variable if compiler not set
    if(NOT (MF_CC AND MF_CXX))
      #C compiler
      if(DEFINED $ENV{CC})
	set(MF_CC $ENV{CC})
      else()
	set(MF_CC "${DEFAULT_MF_CC}")
      endif()

      #C++ compiler
      if(DEFINED $ENV{CXX})
	set(MF_CXX $ENV{CXX})
      else()
	set(MF_CXX "${DEFAULT_MF_CXX}")
      endif()
    endif()

    #Check compiler
    if(NOT MF_CC)
      message(SEND_ERROR "ENVIRONMENT VARIABLE 'CC' NOT FOUND")
    endif()
    if(NOT MF_CXX)
      message(SEND_ERROR "ENVIRONMENT VARIABLE 'CXX' NOT FOUND")
    endif()
  
    #Create folder name
    set(MF_TMP_ARCH "${MF_ARCH}")
    string(REGEX REPLACE "[+]" "[+]" MF_TMP1 ${MF_CXX})
    string(REGEX REPLACE "_${MF_TMP1}+$" "" MF_TMP_ARCH ${MF_TMP_ARCH})
    string(REGEX REPLACE "_${MF_CC}+$" "" MF_TMP_ARCH ${MF_TMP_ARCH})
    string(REGEX REPLACE "_${MF_ACCEL}+$" "" MF_TMP_ARCH ${MF_TMP_ARCH})
    string(REGEX REPLACE "_${MF_MPI}+$" "" MF_TMP_ARCH ${MF_TMP_ARCH})
  
    SET(MF_BINARY_DIR "bin_cmake/${MF_TMP_ARCH}_${MF_MPI}_${MF_ACCEL}_${MF_CC}_${MF_CXX}")

    #Remove useless variables
    unset(MF_TMP_ARCH)
    unset(MF_TMP1)
    unset(MF_TMP2)
    unset(MF_TMP3)
  endif()

# MULTI_PLATFORM MODE
else()
  message(STATUS "Selected MULTI_PLATFORM MODE.")
  set(MF_USER_MODE "multi_platform")

  #List of bin_cmake directories
  SET(MF_CORRECT_DIR_LIST)
  SET(MF_CORRECT_ARCH_LIST)
  SET(MF_WRONG_DIR_LIST)
  
  #
  # Create list of correct architectures
  #
  FILE(GLOB TMP_MF_BIN_DIR_LIST RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/bin_cmake" "bin_cmake/*")
  foreach(d ${TMP_MF_BIN_DIR_LIST})
    if(IS_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin_cmake/${d}")

      #Check correction
      if( ${d} MATCHES "^.*_.*_.*_.*_.*$")
	set(MF_TMP1 "${d}")
	set(MF_TMP2 "")

	#Cut postfix
	string(LENGTH ${d} TMP_LEN1)
	string(REGEX MATCH "_[^_]+_[^_]+_[^_]+_[^_]+$" MF_TMP2 ${d})
	string(LENGTH ${MF_TMP2} TMP_LEN2)
	math(EXPR LEN "${TMP_LEN1} - ${TMP_LEN2}")

	#Create name without postfix
	string(SUBSTRING ${MF_TMP1} 0 ${LEN} MF_TMP2 )

	#Check name with and without postfix
	CHECK_WEAK_MF_ARCH(${MF_TMP1})
	if(${MF_ARCH_CONFIRMED} EQUAL -1)
	  CHECK_WEAK_MF_ARCH(${MF_TMP2})
	  if(${MF_ARCH_CONFIRMED} EQUAL -1)
	    LIST(APPEND MF_WRONG_DIR_LIST ${d})
	  else()
	    if(${MF_TMP2} STREQUAL "REFERENCE_PLATFORM")
	      LIST(APPEND MF_WRONG_DIR_LIST ${d})
	    else()
	      LIST(APPEND MF_CORRECT_DIR_LIST ${d})
	      LIST(APPEND MF_CORRECT_ARCH_LIST ${MF_TMP2})
	    endif()
	  endif()
	else()
	  if(${MF_TMP1} STREQUAL "REFERENCE_PLATFORM")
	    LIST(APPEND MF_WRONG_DIR_LIST ${d})
	  else()
	    LIST(APPEND MF_CORRECT_DIR_LIST ${d})
	    LIST(APPEND MF_CORRECT_ARCH_LIST ${MF_TMP1})
	  endif()
	endif()

	unset(MF_TMP1)
	unset(MF_TMP2)
	unset(MF_LEN1)
	unset(MF_LEN2)
      else()
	LIST(APPEND MF_WRONG_DIR_LIST ${d})
      endif()
    endif()
  endforeach(d)

  #Show wrong directories
  if(MF_WRONG_DIR_LIST)
    set(MF_BAD_ARCH_MSG "This directories will not be used:\n")
    foreach( d ${MF_WRONG_DIR_LIST})
      set(MF_BAD_ARCH_MSG "${MF_BAD_ARCH_MSG} --> ${d}\n")
    endforeach(d)
    message(WARNING "${MF_BAD_ARCH_MSG}")
  endif()
endif()

#
# PART 2: Check for missing variables.
#

# Set build type
if(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
  set(MF_BUILD ${CMAKE_BUILD_TYPE})
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(MF_BUILD ${CMAKE_BUILD_TYPE})
else()
  set(MF_BUILD ${DEFAULT_MF_BUILD})
endif()

# List of additional debug info and time tests
set(MF_DEBUG_FLAG_LIST)
set(MF_TIMETEST_FLAG_LIST)

#Debug info
if(DEBUG_APM)
  LIST(APPEND MF_DEBUG_FLAG_LIST "-DDEBUG_APM=ON")
  message(STATUS "Added DEBUG_APM flag")
  unset(DEBUG_APM)
endif()
if(DEBUG_LSM)
  LIST(APPEND MF_DEBUG_FLAG_LIST "-DDEBUG_LSM=ON")
  message(STATUS "Added DEBUG_LSM flag")
  unset(DEBUG_LSM)
endif()
if(DEBUG_MMM)
  LIST(APPEND MF_DEBUG_FLAG_LIST "-DDEBUG_MMM=ON")
  message(STATUS "Added DEBUG_MMM flag")
  unset(DEBUG_MMM)
endif()
if(DEBUG_PCM)
  LIST(APPEND MF_DEBUG_FLAG_LIST "-DDEBUG_PCM=ON")
  message(STATUS "Added DEBUG_PCM flag")
  unset(DEBUG_PCM)
endif()
if(DEBUG_SIM)
  LIST(APPEND MF_DEBUG_FLAG_LIST "-DDEBUG_SIM=ON")
  message(STATUS "Added DEBUG_SIM flag")
  unset(DEBUG_SIM)
endif()
if(DEBUG_TMM)
  LIST(APPEND MF_DEBUG_FLAG_LIST "-DDEBUG_TMM=ON")
  message(STATUS "Added DEBUG_TMM flag")
  unset(DEBUG_TMM)
endif()

#Time tests
if(TIME_TEST)
  LIST(APPEND MF_TIMETEST_FLAG_LIST "-DTIME_TEST")
  message(STATUS "Added TIME_TEST flag")
  unset(TIME_TEST)
endif()
if(TIME_TEST_2)
  LIST(APPEND MF_TIMETEST_FLAG_LIST "-DTIME_TEST_2")
  message(STATUS "Added TIME_TEST_2 flag")
  unset(TIME_TEST_2)
endif()

#
# PART 3: Lunch actual cmake-build OR lunch ccmake to provide missing variables
#

# Call make or cmake
if(${MF_USER_MODE} STREQUAL "ancient_make")

  #Set MOD_FEM_DIR variable
  if(NOT DEFINED ENV{MOD_FEM_DIR})
    set(ENV{MOD_FEM_DIR} ${MF_DIR})
  endif()

  message(STATUS "PLATFORM FILE: \"${MF_ARCH}\" MAKE OPTIONS: \"${MF_PROBLEM}\" WORKING_DIRECTORY: \"${MF_DIR}\" ")
  # Create directories for compilation objects
  execute_process(COMMAND "make" "--makefile=Makefile_explicit" WORKING_DIRECTORY "${MF_DIR}/src")
  if(MF_PROBLEM)
    # Run compilation with make
    execute_process(COMMAND "make" "--makefile=Makefile_explicit" "${MF_PROBLEM}" WORKING_DIRECTORY "${MF_DIR}/src")
  endif()
  
elseif(${MF_USER_MODE} STREQUAL "single_platform" AND NOT ${MF_ARCH_CONFIRMED} EQUAL -1)
  if(NOT ${MF_CC} STREQUAL "fail")
  
    #Info
    message("Configuration:")
    message(STATUS "MF_ARCH  --> ${MF_ARCH}")
    message(STATUS "MF_CC    --> ${MF_CC}")
    message(STATUS "MF_CXX   --> ${MF_CXX}")
    message(STATUS "MF_ACCEL --> ${MF_ACCEL}")
    message(STATUS "MF_MPI   --> ${MF_MPI}")
    message(STATUS "MF_BUILD --> ${MF_BUILD}")
    message(STATUS "COMPILATION FLAGS: --> ${MF_DEBUG_FLAG_LIST} ${MF_TIMETEST_FLAG_LIST}")

    #Remove binary dir if it exist
    if(EXISTS ${MF_BINARY_DIR})
      message(STATUS "Remove directory: ${MF_BINARY_DIR}")
      FILE(REMOVE_RECURSE "${MF_BINARY_DIR}")
    endif()
  
    #Create new directory
    FILE(MAKE_DIRECTORY "${MF_BINARY_DIR}")

    #Run cmake command
    message(STATUS "Run CMAKE")
    if(${MF_BUILD} STREQUAL "Release")
      execute_process(COMMAND "cmake" "--no-warn-unused-cli" "." "-DCMAKE_BUILD_TYPE:STRING=${MF_BUILD}" "-DMODFEM_ARCH:STRING=${MF_ARCH}" "-DMODFEM_MPI:STRING=${MF_MPI}" "-DMODFEM_ACCEL:STRING=${MF_ACCEL}" "-DCC:STRING=${MF_CC}" "-DCXX:STRING=${MF_CXX}" "${MF_TIMETEST_FLAG_LIST}" "${MF_OPENCL_TASK}" "${CMAKE_SOURCE_DIR}/src" WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/${MF_BINARY_DIR}")
    else()
      execute_process(COMMAND "cmake" "--no-warn-unused-cli" "." "-DCMAKE_BUILD_TYPE:STRING=${MF_BUILD}" "-DMODFEM_ARCH:STRING=${MF_ARCH}" "-DMODFEM_MPI:STRING=${MF_MPI}" "-DMODFEM_ACCEL:STRING=${MF_ACCEL}" "-DCC:STRING=${MF_CC}" "-DCXX:STRING=${MF_CXX}" "${MF_DEBUG_FLAG_LIST}" "${MF_OPENCL_TASK}" "${MF_TIMETEST_FLAG_LIST}" "${CMAKE_SOURCE_DIR}/src" WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/${MF_BINARY_DIR}")
    endif()

    #Build ModFEM
    message(STATUS "Run ModFEM Compilation")
    execute_process(COMMAND "make" WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/${MF_BINARY_DIR}")
  endif()
    
elseif(${MF_USER_MODE} STREQUAL "multi_platform")

  # Parent directory
  SET(MF_BINARY_PARENT_DIR "${CMAKE_BINARY_DIR}/bin_cmake")
  
  foreach(d ${MF_CORRECT_DIR_LIST})
    message(STATUS "Correct directory: ${MF_BINARY_PARENT_DIR}/${d}")
  endforeach(d)

  #
  # Run cmake and compilation for all correct directories
  #

  # Calculate list size
  list(LENGTH MF_CORRECT_DIR_LIST TMP_LEN)
  math(EXPR LEN "${TMP_LEN} - 1")

  # Loop over list elements
  if(NOT ${LEN} EQUAL -1) 
    foreach(val RANGE ${LEN})

      #Get directory and arch
      list(GET MF_CORRECT_DIR_LIST ${val} MF_TMP_DIR)
      list(GET MF_CORRECT_ARCH_LIST ${val} MF_TMP_ARCH)

      #Remove binary dir if it exist
      if(EXISTS "${MF_BINARY_PARENT_DIR}/${MF_TMP_DIR}")
	message(STATUS "Remove directory: ${MF_BINARY_PARENT_DIR}/${MF_TMP_DIR}")
	FILE(REMOVE_RECURSE "${MF_BINARY_PARENT_DIR}/${MF_TMP_DIR}")
      endif()
  
      #Create new directory
      FILE(MAKE_DIRECTORY "${MF_BINARY_PARENT_DIR}/${MF_TMP_DIR}")
    
      #Run cmake command
      message(STATUS "Run CMAKE for ${MF_TMP_DIR}")
      if(${MF_BUILD} STREQUAL "Release")
	execute_process(COMMAND "cmake" "--no-warn-unused-cli" "." "-DCMAKE_BUILD_TYPE:STRING=${MF_BUILD}" "${MF_TIMETEST_FLAG_LIST}" "${MF_OPENCL_TASK}" "${CMAKE_SOURCE_DIR}/src" WORKING_DIRECTORY "${MF_BINARY_PARENT_DIR}/${MF_TMP_DIR}")
      else()
	execute_process(COMMAND "cmake" "--no-warn-unused-cli" "." "-DCMAKE_BUILD_TYPE:STRING=${MF_BUILD}" "${MF_DEBUG_FLAG_LIST}" "${MF_TIMETEST_FLAG_LIST}" "${MF_OPENCL_TASK}" "${CMAKE_SOURCE_DIR}/src" WORKING_DIRECTORY "${MF_BINARY_PARENT_DIR}/${MF_TMP_DIR}")
      endif()

      # Run compilation
      message(STATUS "Run ModFEM Compilation for ${MF_TMP_DIR}")
      execute_process(COMMAND "make" WORKING_DIRECTORY "${MF_BINARY_PARENT_DIR}/${MF_TMP_DIR}")

    endforeach()
  endif()

  #List of platform if not found any
  if(NOT MF_CORRECT_ARCH_LIST)
    message("")
    message("Not found user platform directory, please create one using following rule:")
    message("  [NameOfYourCmakeConfigFile]_[nompi|mpi]_[none|opencl|cuda]_[c_compiler]_[cxx_compiler]")
    message("The list of correct CMake configuration files is below:")
    foreach(a ${MF_ARCH_LIST})
      if(NOT ${a} STREQUAL "REFERENCE_PLATFORM")
	message("-- ${a}")
      endif()
    endforeach(a)
  endif()
  
endif()

message("To finish proccess run:")
message(" -> make clean-upper")
message("or delete following files:")
message(" -> CMakeFiles")
message(" -> CMakeCache.txt")
message(" -> clean-upper.cmake")
message(" -> cmake_install.cmake")
message(" -> Makefile")
message("")

MESSAGE(STATUS "ModFEM launcher end")

