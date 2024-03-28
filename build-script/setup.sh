#!/bin/bash -l
export OTB_SRC_DIR=${HOME}/opal-x/downloads
#
# the build type will be appended by 401-build-opal 
#
export OTB_PREFIX=${HOME}/opal-x/install
export HIK_INSTALL_DIR=${HOME}/opal-x/ippl-build-scripts
export PATH=${OTB_SRC_DIR}/opal-x/build_openmp/src:${PATH}