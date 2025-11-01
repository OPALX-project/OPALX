# cmake/InstallOPALX.cmake
# -------------------------------------------------------
# Installation logic for the OPALX library
# -------------------------------------------------------
# Install headers (optional — keep if your build needs them)
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/opalx
  FILES_MATCHING
  PATTERN "*.h" PATTERN "*.hpp" PATTERN "*.hh" PATTERN "*.H" PATTERN "*.cuh" PATTERN "*.tpp"
  EXCLUDE
  PATTERN "CMakeFiles" EXCLUDE
  PATTERN "CMakeLists.txt" EXCLUDE
  PATTERN "*.cpp" EXCLUDE
  PATTERN "*.cc" EXCLUDE
  PATTERN "*.cu" EXCLUDE)

# Install generated version header
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/OPALXVersions.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/opalx)

# -------------------------------------------------------
# Install the actual executable
install(TARGETS opalx
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})