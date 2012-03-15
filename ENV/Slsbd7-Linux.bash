# IPPL ----------------------------------------------------------------------
export IPPL_ARCH=LINUX-SER
export IPPL_ROOT=~/svnwork/ippl
# EXPDE end ----------------------------------------------------------------

export CXX=g++ #$HOME/install/mpich/bin/mpicxx
export CC=gcc #$HOME/install/mpich/bin/mpicc
export F77=g77 #$HOME/install/mpich/bin/mpif77

#export CXX=$HOME/install/mpich/bin/mpicxx
#export CC=$HOME/install/mpich/bin/mpicc
#export F77=$HOME/install/mpich/bin/mpif77





export MYMPIRUN=mpirun
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$H5HOME/install/lib
