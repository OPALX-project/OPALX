#!/bin/bash
#
mkdir -p config 
aclocal --force
glibtoolize --force --copy
automake --force --add-missing --copy
autoheader --force
autoconf --force
autoreconf
CXX=mpicxx ./configure LDFLAGS=-static --with-ippl-includedir=$HOME/extlib/ippl/include --with-ippl-libdir=$HOME/extlib/ippl/lib
#
make -j 10


