#!/bin/bash

g++ -O2 -o AnimateSimulation makeFrames.cpp \
    -I/home2/kraus/felsim/include \
    -I/home2/kraus/felsim/extlib/root/include \
    -I/home2/kraus/felsim/svnwork/H5Part-S/build/include \
    -I/home2/kraus/felsim/svnwork/hdf5-1.6.5-S/build/include \
    -L/home2/kraus/felsim/lib \
    -L/home2/kraus/felsim/extlib/root/lib \
    -L/home2/kraus/felsim/svnwork/H5Part-S/build/lib \
    -L/home2/kraus/felsim/svnwork/hdf5-1.6.5-S/build/lib \
    -lz \
    -lh5root \
    -lH5Part \
    -lhdf5 \
    -lSpectrum \
    -lCore \
    -lCint \
    -lRIO \
    -lNet \
    -lHist \
    -lGraf \
    -lGraf3d \
    -lGpad \
    -lTree \
    -lRint \
    -lPostscript \
    -lMatrix \
    -lPhysics \
    -lMathCore \
    -lThread \
    -lGui \
    -pthread \
    -lm \
    -ldl \
    -rdynamic
    