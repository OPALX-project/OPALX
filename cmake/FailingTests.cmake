# -----------------------------------------------------------------------------
# FailingTests.cmake
#
# Maintains a list of tests that are excluded from testing due to the fact that they are known to
# fail and spoil an otherwise green dashboard.
#
# Fixing these tests should be considered a priority
#
# -----------------------------------------------------------------------------

if(OPALX_SKIP_FAILING_TESTS)
  message(FATAL_ERROR
    "OPALX_SKIP_FAILING_TESTS is currently not implemented. "
    "It is reserved for future work to exclude known failing tests from the build/test set. "
    "Please configure with -DOPALX_SKIP_FAILING_TESTS=OFF for now.")
endif()