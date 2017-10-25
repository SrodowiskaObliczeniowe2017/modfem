	Use new CMAKE

New CMAKE support three ways of compilation ModFEM, cmake should run
in main ModFEM directory: modular_fem_code_2009

1. ANCIENT_MAKE MODE
 In this mode cmake support make command based on Makefile_explicit.
 To run this mode in system must set variable MOD_FEM_ARCH only!!!
 Cmake take only one parameter -DMF_PROBLEM="problem"
	Example:
 cmake . -DMF_PROBLEM="laplace_prism_std_krb"

2. SINGLE_PLATFORM MODE
 In this mode cmake run for only one architecture defined in
 variable MOD_FEM_ARCH_CMAKE (if MOD_FEM_ARCH has been set too
 this mode will prefered). During the work cmake will create
 correct folder in cmake_bin based on architecture name and passed
 values.
 Available options:
  -DCMAKE_BUILD_TYPE="bulid type" --> Build type [Release/RelWithDebInfo/Debug] default [Release]
  -DMF_CC="compiler"  --> C compiler
  -DMF_CXX="compiler" --> C++ compiler
  -DMF_MPI="value"    --> [mpi/nompi] - MPI settings default [nompi]
  -DMF_ACCEL="value"  --> [none/openmp/opencl] - acceleration settings default [openmp]
	Time mesare flags:
  -DTIME_TEST="switch"   --> [ON/OFF]  
  -DTIME_TEST_2="switch" --> [ON/OFF]
  	OpenCL test flag (only one flag could be enabled):
  -DTEST_SCALAR="switch" --> [ON/OFF]
  -DLAPLACE="switch"     --> [ON/OFF]
	Debug print (working in debug mode):
  -DDEBUG_APM="switch"   --> [ON/OFF]
  -DDEBUG_LSM="switch"   --> [ON/OFF]
  -DDEBUG_MMM="switch"   --> [ON/OFF]  
  -DDEBUG_PCM="switch"   --> [ON/OFF] 
  -DDEBUG_SIM="switch"   --> [ON/OFF]
  -DDEBUG_TMM="switch"   --> [ON/OFF]
 If primary variables is not set, cmake try to detect all values from
 the name of MOD_FEM_ARCH_CMAKE, in other way it use default settings.
	Example:
 export MOD_FEM_ARCH_CMAKE=computerJB-NoMPI-Lapack
 cmake . -DMF_CC="icc" -DMF_CXX="icpc" -DTIME_TEST=ON

3. MULTI_PLATFORM MODE
 In this mode user must create correct folder in bin_cmake directory,
 creation rule:

	[NameOfYourCmakeConfigFile]_[nompi|mpi]_[none|opencl|cuda]_[c_compiler]_[cxx_compiler]

 After tha user can run cmake, cmake will run for all created directories.
 During the work directories has been checked for correction. If directory is
 invalid then it will be deleted.
 In this mode user can only set cmake values:
  -DCMAKE_BUILD_TYPE="bulid type" --> Build type [Release/RelWithDebInfo/Debug] default [Release] 
	Time mesare flags:
  -DTIME_TEST="switch"   --> [ON/OFF]  
  -DTIME_TEST_2="switch" --> [ON/OFF]
    	OpenCL test flag (only one flag could be enabled):
  -DTEST_SCALAR="switch" --> [ON/OFF]
  -DLAPLACE="switch"     --> [ON/OFF]
	Debug print (working in debug mode):
  -DDEBUG_APM="switch"   --> [ON/OFF]
  -DDEBUG_LSM="switch"   --> [ON/OFF]
  -DDEBUG_MMM="switch"   --> [ON/OFF]  
  -DDEBUG_PCM="switch"   --> [ON/OFF] 
  -DDEBUG_SIM="switch"   --> [ON/OFF]
  -DDEBUG_TMM="switch"   --> [ON/OFF]
 Other settings will read from directory name.
	Example:
 cmake . -DCMAKE_BUILD_TYPE="Debug" -DTIME_TEST=ON

In mode 2 and 3 cmake will recreate directory.
COMPILATION with make will run directly after
CMake execution.

4. The older method of ModFEM compilation with using config.sh and other will still working.
5. Afrer run cmake followed files and directories should be delete:
  - modular_fem_code_2009/CMakeFiles
  - modular_fem_code_2009/CMakeCache.txt
  - modular_fem_code_2009/clean-upper.cmake
  - modular_fem_code_2009/cmake_install.cmake
  - modular_fem_code_2009/Makefile

	OR

  - run command:
  make clean-upper

# Changes:
01.05.2015 - bug fixes
30.04.2015 - write new cmake

# ------------------ OLD DOCUMENTATION BEFORE 30.04.2015 --------------------------------------------------------- #	
Create folder for your build here.
Folder name:

[NameOfYourCmakeConfigFile]_[nompi|mpi]_[none|opencl|cuda]_[c_compiler]_[cxx_compiler]

for ex.:

notebookPP_nompi_none_gcc_g++

(assuming notebookPP.cmake file is present in src/cmake/Platforms)

Enable debug/release mode:
cmake -DCMAKE_BUILD_TYPE=Debug .
cmake -DCMAKE_BUILD_TYPE=Release .

IMPORTANT:
In case of linking errors, try: "-ldl -liomp5 -lpthread" under CMAKE_EXE_LINKER_FLAGS (ccmake ../../src then 't')
