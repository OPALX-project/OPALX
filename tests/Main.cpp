#include "gtest/gtest.h"

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  int test_out = RUN_ALL_TESTS();
  return test_out;
}

