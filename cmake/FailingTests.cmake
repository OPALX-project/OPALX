# -----------------------------------------------------------------------------
# FailingTests.cmake
#
# Maintains a list of tests that are excluded from testing due to the fact that they are known to
# fail and spoil an otherwise green dashboard.
#
# Fixing these tests should be considered a priority
#
# -----------------------------------------------------------------------------

if(BUILD_TESTING AND OPALX_SKIP_FAILING_TESTS)
  #set(OPALX_DISABLED_TEST_LIST
  #    AssembleRHS)
endif()