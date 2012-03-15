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
    --with-h5part-includedir=$H5Part/src --with-h5part-libdir=$H5Part/src
make -j 10


