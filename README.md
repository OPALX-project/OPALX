# This is OPAL-X


# BUILDING OPAL-X on Merlin


## Modules needed

load_ipplomp() { 
module purge
module use unstable
module load cmake/3.25.2
module load gcc/11.4.0      # does not work for openmp and mixed module load gcc/12.3.0
module load openmpi
module load boost
module load fftw
module load gtest
module load gsl
module load hdf5
module load H5hut
module load boost
module load cuda/12.1.1
export BOOST_ROOT=$BOOST_DIR
export OMP_PROC_BIND=spread
}

## build for Kokkos, FEFFTe and IPPL

% git clone git@github.com:IPPL-framework/ippl-build-scripts.git

set the correct environment, adjust the path in ITB_DOWNLOAD_DIR

% export ITB_DOWNLOAD_DIR=/data/project/general/isodarUQ/adelmann
% export ITB_SRC_DIR=${ITB_DOWNLOAD_DIR}/downloads
% export ITB_PREFIX=${ITB_DOWNLOAD_DIR}/install
% mkdir -p ${ITB_SRC_DIR}
% mkdir -p ${ITB_PREFIX}

Now you can build serial and with openmp

% ./999-build-everything -t serial  --kokkos --heffte --ippl --export -u
% ./999-build-everything -t openmp --enable-openmp --kokkos --heffte --ippl --export -u


here is still a problem with cuda (./999-build-everything -t mixed --enable-cuda --enable-openmp --kokkos --heffte --ippl --export --arch=PASCAL6 -u)

## Checkout OPAL-X

% git clone git@gitlab.psi.ch:OPAL/opal-x/src.git 

Change to build-script

401-build-opal -r openmp --export -u 




The documentation has been moved to the [Wiki](https://gitlab.psi.ch/OPAL/src/wikis/home).
