# Third-Party Licenses

This directory contains the licenses for all third-party dependencies statically linked into OPALX binaries. 

## License Summary

| Dependency | License | Type | Notes |
|-----------|---------|------|-------|
| **Heffte** | BSD 3-Clause | Static | FFT library for distributed arrays |
| **Kokkos** | Apache 2.0 + LLVM Exceptions | Static | Performance library for parallel computing |
| **IPPL** | GPLv3 | Static | The IPPL framework (required for OPALX) |
| **HDF5** | BSD + proprietary (LBNL) | Static | Scientific data format library |
| **H5Hut** | BSD 3-Clause | Static | I/O and parallel file format library |
| **ZLIB** | zlib License | Static | Compression library |
| **MPI** | Varies | Dynamic | Message Passing Interface (system library) |
| **GCC Runtime** | GCC Runtime Library Exception | Dynamic | C++ standard library |

## Important Notes

### Static Linking Implications

Since OPALX is built with `BUILD_SHARED_LIBS=OFF` by default, all third-party libraries listed above are **statically linked** into the binary. This means:

1. **License Compliance**: When distributing OPALX binaries, you must include copies of these licenses.
2. **No Source Distribution Required**: Static linking does not require distributing source code for permissively-licensed libraries (BSD, Apache 2.0).
3. **GPLv3 Compliance**: IPPL is licensed under GPLv3. You should be aware of the copyleft implications when distributing OPALX binaries.

### Dynamic Linking Alternative

If you build OPALX with `BUILD_SHARED_LIBS=ON`, most dependencies are dynamically linked:

- You do **not** need to include third-party library licenses in your binary distribution.
- End users must install required libraries on their systems.
- System libraries (MPI, compiler runtime) remain dynamically linked in either case.

## Build Configuration

Check the current linking configuration:

```bash
cd build
cmake -L | grep BUILD_SHARED_LIBS
```

Or examine `CMakeCache.txt`:

```bash
grep "BUILD_SHARED_LIBS" build/CMakeCache.txt
```

## Distribution Guidelines

When uploading OPALX binaries to GitHub:

1. **Always include this LICENSES directory** with the OPALX binary release.
2. Add a **LICENSE** file at the root of the repository (OPALX's own license).
3. Include a **THIRD_PARTY_LICENSES.md** or similar file referencing this directory.
4. Document the build configuration used (static vs. dynamic linking).

### Example Release Description

```
OPALX Binary v1.0.0

This binary release includes statically-linked third-party libraries. 
Please see the LICENSES/ directory for license information of all 
dependencies bundled with this binary.
```

## License Details

- **Heffte/**: BSD 3-Clause license
- **Kokkos/**: Apache 2.0 license with LLVM exceptions  
- **IPPL/**: GNU General Public License v3
- **HDF5/**: Available in existing LICENSES/hdf5/ directory
- **H5Hut/**: Available in existing LICENSES/H5hut/ directory
- **ZLIB/**: zlib license
- **Other libraries/**: See existing subdirectories in LICENSES/

For more information on each library:

- Heffte: https://github.com/icl-utk-edu/heffte
- Kokkos: https://github.com/kokkos/kokkos
- IPPL: https://github.com/IPPL-framework/ippl
- HDF5: https://www.hdfgroup.org/
- H5Hut: https://github.com/eth-cscs/h5hut
- ZLIB: https://zlib.net/
