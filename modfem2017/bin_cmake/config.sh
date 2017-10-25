#!/bin/bash


# Here is fix to avoid never ending cycle of:
# 1. clear dir of build of all files
# 2. copy 'config.sh' from 'bin_cmake' into build dir
# 3. run './config.sh' in build dir
# just in one step.
# In order to make in work just call 'config.sh' from 'bin_cmake' directory
# with additional argument naming existing build folder.
# E.g.: ./config.sh pcKM_MKL_mpi_none_mpicc_mpic++
# Date: 04.11.2015
# Author: kamich@agh.edu.pl
if [ ! "$1" = "" ]; then
	echo " Clearing $1"
	echo "$(rm -rf $1/*)"
	echo " Copying config.sh into $1"
	echo "$(cp config.sh $1/config.sh)"
	echo "Entering $1"
	cd $1
	echo "Inside $(pwd)"
	echo "Running configuration script..."
	./config.sh
	exit
fi
#end of fix.

#BINDIR=BAD
BINDIR="$(basename `pwd`)"
#echo $BINDIR
varCXX=$(echo ${BINDIR//*_})
#echo "varCXX" $varCXX
woCXX=$(echo ${BINDIR%_$varCXX})
#echo "woCXX" $woCXX
varCC=$(echo ${woCXX//*_})
#echo "varCC" $varCC
woCXX_CC=$(echo ${woCXX%_$varCC})
#echo "woCXX_CC" $woCXX_CC
accel=$(echo ${woCXX_CC//*_})
#echo "accel" $accel
woCXX_CC_accel=$(echo ${woCXX_CC%_$accel})
#echo "woCXX_CC_accel" $woCXX_CC_accel
mpi=$(echo ${woCXX_CC_accel//*_})
#echo "mpi" $mpi
cmakefile=$(echo ${woCXX_CC_accel%_$mpi})
#echo "cmakefile" $cmakefile

#echo $compiler
#echo "Cmake config file is:" $cmakefile".cmake"
export MOD_FEM_ARCH_CMAKE="$cmakefile"
echo "MOD_FEM_ARCH_CMAKE is:" $MOD_FEM_ARCH_CMAKE
echo "Cmake config file is:" $cmakefile".cmake"

if [ ! -f ../../src/cmake/Platforms/"$cmakefile".cmake ]
then
    echo "File" $cmakefile".cmake not found in ../../src/cmake/Platforms."
    echo "Name your directory: [CMAKE_CONFIG_FILENAME]_[mpi|nompi]_[none|opencl|cuda]_[CC]_[CXX]"
    echo "Exiting."
    exit
fi

export CC="$varCC"
export CXX="$varCXX"

#if [ "$compiler" == "gcc" ]
#then
#    export CC=gcc
#    export CXX=g++
#elif [ "$compiler" == "icc" ]
#then
#    export CC=icc
#    export CXX=icpc
#else
#    export CC=gcc
#    export CXX=g++
#    echo "Compiler unknown/not present in folder name. Setting default."
#    echo "Name your directory: [CMAKE_CONFIG_FILENAME]_[gcc|icc]"
#fi

echo "C compiler is:" $CC
echo "C++ compiler is:" $CXX
echo "MODFEM_MPI:" $mpi
echo "MODFEM_ACCEL:" $accel

echo "Starting cmake..."
cmake -G "CodeBlocks - Unix Makefiles" -DMODFEM_MPI:STRING=$mpi -DMODFEM_ACCEL:STRING=$accel ../../src

## Below fix is moved from here to src/CMakeLists.txt!
## DEPRECATED:
## Here is fix for debian based systems and MKL (without this, no binary will ever built because of linkage failure).
## Date: 04.11.2015
## Author: kamich@agh.edu.pl
#if [ -f /etc/debian_version ]; then
#	echo " Altering linker input"
#	linker_files="$(find ./ -name 'link.txt')"
#	echo " Found linker files: $linker_files. Now updateing "
## For Intel MKL
#	old="libmkl_core.a -Wl,--end-group"
#	new="libmkl_core.a  -lm -ldl -lgfortran -Wl,--end-group"
#	replaced="$(sed -i "s/$old/$new/g" $linker_files)"
## For Generic BLAS/LAPCK
#	old="-llapack -lblas -lm -ldl"
#	new="-llapack -lblas -lgfortran -lquadmath -lm -ldl"
#	replaced="$(sed -i "s/$old/$new/g" $linker_files)"
#fi
##end of fix.

# Here is fix for auto-creation of scripts to run mpi, without typing a lot.
# Date: 04.11.2015
# Author: kamich@agh.edu.pl
if [ "$mpi" = "mpi" ]; then
	echo "Creating lunch script for mpi"
        lsmpi="$(echo mpirun -n 2 xterm -hold -e gdb -ex run --args MOD_FEM_ns_supg_hybrid_std_d ../../examples/pdd_ns_supg/LDC/nas_20x20 >> run_mpi_gdb.sh)"
        lsmpi="$(echo mpirun -n 2 xterm -hold -e valgrind MOD_FEM_ns_supg_hybrid_std_d ../../examples/pdd_ns_supg/LDC/nas_20x20 >> run_mpi_val.sh)"
        lsmpi="$(echo mpirun -n 2 xterm -hold -e valgrind --vgdb-error=0 MOD_FEM_ns_supg_hybrid_std_d ../../examples/pdd_ns_supg/LDC/nas_20x20 >> run_mpi_val_gdb.sh)"
fi
# end of fix.

#for x in $arr
#do
#    echo "> [$x]"
#done
