# Third-Party Licenses

This directory contains the license texts for dependencies fetched by OPALX at build time.

## Covered Dependencies

| Dependency | Typical Linking | License Files in This Directory | Repository | Upstream License |
|-----------|------------------|----------------------------------|------------|------------------|
| OPAL (legacy codebase) | Incorporated source code | `OPAL/LICENSE` | https://github.com/OPALX-project/OPAL | https://github.com/OPALX-project/OPAL/blob/master/LICENSES/OPAL/LICENSE |
| OpenMPI (optional MPI implementation) | Usually dynamic system library | `OpenMPI/LICENSE` | https://github.com/open-mpi/ompi | https://github.com/open-mpi/ompi/blob/main/LICENSE |
| MPICH (optional MPI implementation) | Usually dynamic system library | `MPICH/LICENSE` | https://github.com/pmodels/mpich | https://github.com/pmodels/mpich/blob/main/COPYRIGHT |
| Heffte | Static (via IPPL) | `Heffte/LICENSE` | https://github.com/icl-utk-edu/heffte | https://github.com/icl-utk-edu/heffte/blob/master/LICENSE |
| Kokkos | Static (via IPPL) | `Kokkos/LICENSE` | https://github.com/kokkos/kokkos | https://github.com/kokkos/kokkos/blob/master/LICENSE |
| IPPL | Static | `IPPL/LICENSE` | https://github.com/IPPL-framework/ippl | https://github.com/IPPL-framework/ippl/blob/master/LICENSE |
| HDF5 | Static by default in this project | `HDF5/LICENSE` | https://github.com/HDFGroup/hdf5 | https://github.com/HDFGroup/hdf5/blob/develop/COPYING |
| H5hut | Static by default in this project | `H5hut/LICENSE` | https://github.com/eth-cscs/h5hut | https://github.com/eth-cscs/h5hut/blob/master/COPYING |
| ZLIB | Usually dynamic system library | `ZLIB/LICENSE` | https://github.com/madler/zlib | https://github.com/madler/zlib/blob/develop/LICENSE |
| GoogleTest (tests only) | Static in test builds | `GTest/LICENSE` | https://github.com/google/googletest | https://github.com/google/googletest/blob/main/LICENSE |

## Notes on Distribution

`BUILD_SHARED_LIBS=OFF` is the default in this project, so third-party code is typically statically linked into OPALX artifacts.

For MPI specifically: include the license for the MPI implementation used to build or package your binaries.
- If you distribute binaries together with OpenMPI runtime files, include `OpenMPI/LICENSE`.
- If you distribute binaries together with MPICH runtime files, include `MPICH/LICENSE`.
- If MPI comes entirely from the target system (not redistributed by you), these files are still useful as references but may not be strictly required in your release bundle.

If you publish binaries, include this `LICENSES/` directory with the release and keep the repository root `LICENSE` file.

## Traceability

These license files were copied from the downloaded source trees.
The OPAL license (`OPAL/LICENSE`) is included because OPALX contains substantial legacy OPAL code.
