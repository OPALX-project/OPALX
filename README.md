# This is OPAL-X


# BUILDING OPAL-X on Merlin


## Modules needed OPENMP build

cmake/3.25.2

openmpi/4.1.5_slurm

fftw/3.3.10_merlin6    

gsl/2.7                

H5hut/2.0.0rc6_slurm

gcc/12.3.0             

boost/1.82.0_slurm     

gtest/1.13.0-1         

hdf5/1.10.8_slurm     

gnutls/3.5.19


## Clone repo and build opal-x with OPENMP 

% git clone git@gitlab.psi.ch:OPAL/opal-x/src.git opal-x

% cd opal-x

% ./gen_OPALrevision

% mkdir build_openmp && cd build_openmp

% cmake .. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_STANDARD=20  -DENABLE_SOLVERS=ON  -DENABLE_FFT=ON -DIPPL_PLATFORMS=openmp



The documentation has been moved to the [Wiki](https://gitlab.psi.ch/OPAL/src/wikis/home).
