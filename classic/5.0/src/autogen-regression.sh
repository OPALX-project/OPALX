#!/bin/bash
#
mkdir -p config
aclocal --force
libtoolize --force --copy
automake --force --add-missing --copy
autoheader --force
autoconf --force
autoreconf

./configure --with-ippl-includedir=/home2/ineichen/svnwork/ippl//src --with-h5part-includedir=/home2/ineichen/felsim/H5Part/src --with-h5part-libdir=/home2/ineichen/felsim/H5Part/src --with-hdf5-includedir=/opt/hdf5/hdf5-1.6.10-openmpi-1.2.6-intel-11.1/include --with-hdf5-libdir=/opt/hdf5/hdf5-1.6.10-openmpi-1.2.6-intel-11.1/lib --with-gsl-includedir=/opt/gsl/gsl-1.12/include --with-gsl-libdir=/opt/gsl/gsl-1.12/lib CC=mpicc CXX=mpicxx F77=mpif77
