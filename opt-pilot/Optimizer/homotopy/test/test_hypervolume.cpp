#include "gtest/gtest.h"

#include <string>
#include "wfgHypervolume/hypervolume.h"

// processes each front from the file 
TEST(HYPERVOLUMETest, Alive) {
	EXPECT_NEAR(0.444181, Hypervolume::FromFile("./data/hypervolume.dat"), 1e-6);
}