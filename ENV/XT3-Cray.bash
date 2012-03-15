# IPPL ----------------------------------------------------------------------
export IPPL_ARCH=XT3
export IPPL_ROOT=~/svnwork/ippl
# EXPDE end ----------------------------------------------------------------


module load craypat
module load gcc-catamount
module load pgi/6.0.5
module load xt-lustre-ss
export MPICH_CXX=~roberto/CSCS/MAD9P_GNU/qk-g++

export CXX=CC
export CPP=cc
export F77=ftn
export MYMPIRUN=yod

#export IOBUF_PARAMS='%stdout,dataH5.dat'
#export IOBUF_PARAMS='*:verbose'
#export IOBUF_VERBOSE=1


