# AdaptBins

This repository contains my progress on the AdaptBins class created for OPALX. The testing is done using an adapted version of the "PICnd.cpp" from IPPL. The AdaptBins files are later supposed to be included in the OPALX repository. 

Also add to `CMakeLists.txt` inside the `test/` folder:
```
add_subdirectory (binning)
```

In the `CMakeLists.txt` on above can contain the following to disable all the other tests:
```
# Define an option for ONLY_BINNING
option(ONLY_BINNING "Enable ONLY_BINNING mode" OFF)

if(NOT ONLY_BINNING)
  add_subdirectory (kokkos)
  add_subdirectory (types)
  add_subdirectory (field)
  if (ENABLE_FFT)
      add_subdirectory (FFT)
  endif()
  if (ENABLE_SOLVERS)
      add_subdirectory (solver)
  endif()
  add_subdirectory (particle)
  add_subdirectory (region)
  add_subdirectory (random)
  add_subdirectory (serialization)
endif()
add_subdirectory (binning)
```
