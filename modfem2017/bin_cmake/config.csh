#!/bin/csh

set BINFDIR=`pwd`
set BINDIR=`basename $BINFDIR`

set varCXX=`echo "$BINDIR" | sed "s/.*_//"`
#echo "varCXX" $varCXX

set woCXX = `echo "$BINDIR" | sed "s/_$varCXX//"`
#echo "woCXX" $woCXX

set varCC = `echo "$woCXX" | sed "s/.*_//"`
#echo "varCC" $varCC

set woCXX_CC=`echo "$woCXX" | sed "s/_$varCC//"`
#echo "woCXX_CC" $woCXX_CC

set accel=`echo "$woCXX_CC" | sed "s/.*_//"`
#echo "accel" $accel

set woCXX_CC_accel=`echo "$woCXX_CC" | sed "s/_$accel//"`
#echo "woCXX_CC_accel" $woCXX_CC_accel

set mpi=`echo "$woCXX_CC_accel" | sed "s/.*_//"`
#echo "mpi" $mpi

set cmakefile=`echo "$woCXX_CC_accel" | sed "s/_$mpi//"`
#echo "cmakefile" $cmakefile


setenv MOD_FEM_ARCH_CMAKE "$cmakefile"
echo "MOD_FEM_ARCH_CMAKE is:" $MOD_FEM_ARCH_CMAKE
echo "Cmake config file is:" $cmakefile".cmake"

if ( ! -f ../../src/cmake/Platforms/"$cmakefile".cmake ) then
    echo "File" $cmakefile".cmake not found in ../../src/cmake/Platforms."
    echo "Name your directory: [CMAKE_CONFIG_FILENAME]_[MPI]_[ACCEL]_[CC]_[CXX]"
    echo "Exiting."
    exit
endif

setenv CC "$varCC"
setenv CXX "$varCXX"

echo "C compiler is:" $CC
echo "C++ compiler is:" $CXX
echo "MODFEM_MPI:" $mpi
echo "MODFEM_ACCEL:" $accel

echo "Starting cmake..."
cmake -DMODFEM_MPI:STRING=$mpi -DMODFEM_ACCEL:STRING=$accel ../../src
