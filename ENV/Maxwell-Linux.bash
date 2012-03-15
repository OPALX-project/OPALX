# IPPL ----------------------------------------------------------------------
export IPPL_ARCH=LINUX
export IPPL_ROOT=/home/adelmann/svnwork/ippl
# EXPDE end ----------------------------------------------------------------

export CXX=mpic++
export CC=mpicc
export F77=mpif77
export MYMPIRUN=mpirun
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$H5HOME/install/lib

source /home/adelmann/svnwork/femaXX/install_3pl/config_maxwell
export FEMAXX_PREFIX=/home/adelmann/svnwork/femaXX
