#!/bin/sh
export MOD_FEM_ARCH_CMAKE=$2

#LD_LIBRARY_PATH=/opt/intel/composer_xe_2011_sp1.9.293/compiler/lib/intel64:/opt/intel/composer_xe_2011_sp1.9.293/ipp/../compiler/lib/intel64:/opt/intel/composer_xe_2011_sp1.9.293/ipp/lib/intel64:/opt/intel/composer_xe_2011_sp1.9.293/compiler/lib/intel64:/opt/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64:/opt/intel/composer_xe_2011_sp1.9.293/tbb/lib/intel64//cc4.1.0_libc2.4_kernel2.6.16.21:/opt/intel/composer_xe_2011_sp1.9.293/debugger/lib/intel64:/opt/intel/composer_xe_2011_sp1.9.293/mpirt/lib/intel64:/opt/intel/composer_xe_2011_sp1.9.293/mkl/include:/opt/intel/composer_xe_2011_sp1.9.293/mkl:/opt/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64:/opt/intel/composer_xe_2011_sp1.9.293/compiler/lib/intel64


export LIBRARY_PATH=$LIBRARY_PATH:$LD_LIBRARY_PATH:$3
export LIBRARY_PATH=/opt/gnu/gcc/lib64/:$LIBRARY_PATH
export LD_LIBRARY_PATH=$LIBRARY_PATH

#LIBRARY_PATH=/opt/intel/composer_xe_2011_sp1.9.293/compiler/lib/intel64:/opt/intel/composer_xe_2011_sp1.9.293/ipp/../compiler/lib/intel64:/opt/intel/composer_xe_2011_sp1.9.293/ipp/lib/intel64:/opt/intel/composer_xe_2011_sp1.9.293/compiler/lib/intel64:/opt/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64:/opt/intel/composer_xe_2011_sp1.9.293/tbb/lib/intel64//cc4.1.0_libc2.4_kernel2.6.16.21:/opt/intel/composer_xe_2011_sp1.9.293/debugger/lib/intel64:/opt/intel/composer_xe_2011_sp1.9.293/mpirt/lib/intel64:/opt/intel/composer_xe_2011_sp1.9.293/mkl/include:/opt/intel/composer_xe_2011_sp1.9.293/mkl:/opt/intel/composer_xe_2011_sp1.9.293/mkl/lib/intel64:/opt/intel/composer_xe_2011_sp1.9.293/compiler/lib/intel64
#export LIBRARY_PATH

#echo $LD_LIBRARY_PATH
#echo $LIBRARY_PATH

#katalog do skompilowania
cd $1
if [ -f /usr/bin/gmake ]; then
    /usr/bin/gmake
else
    /usr/bin/make
fi
