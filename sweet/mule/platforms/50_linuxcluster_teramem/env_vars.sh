#
# Default configuration file
#


#
# Don't touch intel stuff, otherwise the mpicxx is not available
#
#module unload intel
#module load intel/19.1.0


module load salloc_conf/teramem
module load intel-oneapi-compilers intel-mkl intel-mpi


export MULE_JOB_SCHEDULER_NUM_JOB_LIMITATION=2
export MULE_JOB_SCHEDULER_SLEEP_SECS=10


#
# Compiler environment
#
export F90=ifort
export CC=icc
export CXX=icpc
export FC=$F90
export LD=ld
export MULE_LINK=$MULE_CXX

export MULE_MPICC=mpicc
export MULE_MPICXX=mpicxx
export MULE_MPIF90=mpif90


export MULE_MPILINK=mpicxx
# If we link with mpif90, we have to add stdc++ for C++
#export MULE_MPILIBS=stdc++
export MULE_MPILIBS=gfortran


export MULE_CC_COMPILER=intel
export MULE_CXX_COMPILER=intel
export MULE_F90_COMPILER=intel
