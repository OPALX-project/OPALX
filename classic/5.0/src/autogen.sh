#!/bin/bash
#
mkdir -p config
aclocal --force
libtoolize --force --copy
automake --force --add-missing --copy
autoheader --force
autoconf --force
autoreconf
CXX=mpicxx ./configure --with-ippl-includedir=$IPPL_ROOT/src --with-h5part-includedir=$H5Part/src --with-h5part-libdir=$H5Part/src
make -j 10


