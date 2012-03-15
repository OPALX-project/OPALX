#!/bin/bash
#
mkdir -p config 
aclocal --force
libtoolize --force --copy
automake --force --add-missing --copy
autoheader --force
autoconf --force
autoreconf
CXX=mpicxx ./configure  --without-wake-test --without-envelope-solver --with-classic-includedir=$HOME/svnwork/OPAL/classic/5.0/src --with-classic-libdir=$HOME/svnwork/OPAL/classic/5.0/src --with-doom-includedir=$HOME/svnwork/OPAL/doom/ \
--with-doom-libdir=$HOME/svnwork/OPAL/doom/ --with-ippl-includedir=$HOME/svnwork/ippl/src --with-ippl-libdir=$HOME/svnwork/ippl/lib/LINUX --with-h5part-includedir=$HOME/svnwork/H5Part/src \
--with-h5part-libdir=$HOME/svnwork/H5Part/src --with-hdf5-includedir=$HDF5HOME/include --with-hdf5-ibdir=$HDF5HOME/lib \
--with-libs="-lz -lm -lfftw3" 
make -j 20


