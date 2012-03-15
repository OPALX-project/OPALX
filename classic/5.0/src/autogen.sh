#!/bin/bash
#
mkdir -p config
aclocal --force
libtoolize --force --copy
automake --force --add-missing --copy
autoheader --force
autoconf --force
autoreconf
CXX=mpicxx ./configure --with-ippl-includedir=$IPPL_ROOT/src
make -j 10


