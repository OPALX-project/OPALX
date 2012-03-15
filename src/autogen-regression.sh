#!/bin/bash
#
mkdir -p config 
aclocal --force
libtoolize --force --copy
automake --force --add-missing --copy
autoheader --force
autoconf --force
autoreconf
CXX=mpicxx ./configure --with-classic-includedir=$CLASSIC_ROOT/src --with-classic-libdir=$CLASSIC_ROOT/src \
    --with-doom-includedir=$DOOM_ROOT --with-doom-libdir=$DOOM_ROOT \
    --with-ippl-includedir=$IPPL_ROOT/src --with-ippl-libdir=$IPPL_ROOT/lib/$IPPL_ARCH \
    --with-h5part-includedir=$H5Part/src --with-h5part-libdir=$H5Part/src \
    --with-hdf5-includedir=/opt/hdf5/hdf5-1.6.10-openmpi-1.2.6-intel-11.1/include --with-hdf5-libdir=/opt/hdf5/hdf5-1.6.10-openmpi-1.2.6-intel-11.1/lib \
    --with-gsl-includedir=/opt/gsl/gsl-1.12/include --with-gsl-libdir=/opt/gsl/gsl-1.12/lib \
    --with-libdir="-L/opt/parmetis/parmetis-3.1 -L/opt/intel-mkl/mkl-10.0/lib/em64t"  \
    --with-libs="-lz -lm -lparmetis -lmetis -lirc" \
    --with-blas=mkl --with-lapack=mkl \
    --with-trilinos-includedir=/opt/trilinos/trilinos-10.2.0/include --with-trilinos-libdir=/opt/trilinos/trilinos-10.2.0/lib \
    --disable-ml-solver
#    --enable-ml-solver
#make -j 4


