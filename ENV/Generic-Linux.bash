# HDF5 ---------------------------------------------------------------------
export H5HOME=~/svnwork/hdf5inst
# end HDF5 -----------------------------------------------------------------

# IPL ----------------------------------------------------------------------
export IPL_ARCH=LINUX
export IPL_ROOT=~/svnwork/ipl
# EXPDE end ----------------------------------------------------------------

# GenTrackE start ----------------------------------------------------------
export GTHOME=~/svnwork/gte
export GTCLASSIC=~/svnwork/classic/3.3/src
# GenTrackE end ------------------------------------------------------------

# CMEE start ---------------------------------------------------------------
export CMEEHOME=~/svnwork/cmee-1.1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CMEEHOME/src/.libs
# CMEE end -----------------------------------------------------------------

export CXX=$HOME/mpich/bin/mpicxx
export CC=$HOME/mpich/bin/mpicc
export F77=$HOME/mpich/bin/mpif77
export MYMPIRUN=$HOME/mpich/bin/mpirun

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/dshome/adelmann/svnwork/hdf5inst/lib
