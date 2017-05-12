/* 
 *  Copyright (c) 2017, Chris Rogers
 *  All rights reserved.
 *  Redistribution and use in source and binary forms, with or without 
 *  modification, are permitted provided that the following conditions are met: 
 *  1. Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer. 
 *  2. Redistributions in binary form must reproduce the above copyright notice, 
 *     this list of conditions and the following disclaimer in the documentation 
 *     and/or other materials provided with the distribution.
 *  3. Neither the name of STFC nor the names of its contributors may be used to 
 *     endorse or promote products derived from this software without specific 
 *     prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
 *  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
 *  POSSIBILITY OF SUCH DAMAGE.
 */

#include "gtest/gtest.h"
#include "Fields/Interpolation/NDGrid.h"

// namespace ndgridtest {
class NDGridTest : public ::testing::Test { 
public: 
    NDGridTest() : grid_m(NULL) { 
    }
 
    void SetUp( ) { 
        std::vector< std::vector<double> > gridCoordinates(2);
        gridCoordinates[0] = std::vector<double>(2, 0.);
        gridCoordinates[0][1] = 3.;	
        gridCoordinates[1] = std::vector<double>(3, 1.);
        gridCoordinates[1][1] = 5.;	
        gridCoordinates[1][2] = 9.;
        grid_m = new interpolation::NDGrid(gridCoordinates);
    }
 
    void TearDown( ) { 
        delete grid_m;
        grid_m = NULL;
    }

    ~NDGridTest() {
    }

private:
    interpolation::NDGrid* grid_m;

};

TEST_F(NDGridTest, DefaultConstructorTest) {
    interpolation::NDGrid grid;
    EXPECT_EQ(grid.begin(), grid.end());
}

// also tests size
TEST_F(NDGridTest, Constructor1Test) {
    int size[] = {5, 6, 7, 8};
    double spacing[] = {1., 2., 3., 4.};
    double min[] = {-1., -2., -3., -4.};
    interpolation::NDGrid grid(4, size, spacing, min);

    ASSERT_EQ(grid.getPositionDimension(), 4);
    for (int i = 0; i < 4; ++i) {
        ASSERT_EQ(grid.size(i), size[i]);
        EXPECT_NEAR(grid.coord(1, i), min[i], 1e-12) << "Failed for i " << i;
        EXPECT_NEAR(grid.coord(size[i], i), min[i]+spacing[i]*(size[i]-1), 1e-12) << "Failed for i " << i;
    }
}

TEST_F(NDGridTest, Constructor2Test) {
    std::vector<int> size(2);
    size[0] = 2;
    size[1] = 3;
    std::vector<const double*> gridCoordinates(2, NULL);
    const double vec1[] = {0., 2.};
    const double vec2[] = {1., 5., 9.};
    gridCoordinates[0] = vec1;
    gridCoordinates[1] = vec2;
    interpolation::NDGrid grid1(size, gridCoordinates);
    ASSERT_EQ(grid1.getPositionDimension(), 2);
    for (int i = 0; i < 2; ++i) {
        ASSERT_EQ(grid1.size(i), size[i]);
        EXPECT_NEAR(grid1.coord(1, i), gridCoordinates[i][0], 1e-12) << "Failed for i " << i;
        EXPECT_NEAR(grid1.coord(i+2, i), gridCoordinates[i][size[i]-1], 1e-12) << "Failed for i " << i;
    }   
}

TEST_F(NDGridTest, Constructor3Test) {
    std::vector< std::vector<double> > gridCoordinates(2);
    gridCoordinates[0] = std::vector<double>(2, 0.);
    gridCoordinates[1] = std::vector<double>(3, 1.);
    gridCoordinates[1][2] = 9.;
    interpolation::NDGrid grid1(gridCoordinates);
    ASSERT_EQ(grid1.getPositionDimension(), 2);
    for (int i = 0; i < 2; ++i) {
        size_t size = gridCoordinates[i].size();
        ASSERT_EQ(grid1.size(i), size);
        EXPECT_NEAR(grid1.coord(1, i), gridCoordinates[i][0], 1e-12) << "Failed for i " << i;
        EXPECT_NEAR(grid1.coord(i+2, i), gridCoordinates[i][size-1], 1e-12) << "Failed for i " << i;
    }   
}

TEST_F(NDGridTest, CoordTest) {
    std::vector< std::vector<double> > gridCoordinates(2);
    gridCoordinates[0] = std::vector<double>(2, 0.);
    gridCoordinates[1] = std::vector<double>(3, 1.);
    gridCoordinates[1][2] = 9.;
    interpolation::NDGrid grid_var(gridCoordinates);
    const interpolation::NDGrid grid_const(gridCoordinates);
    for (int i = 0; i < 2; ++i) {
            EXPECT_NEAR(grid_var.coord(1, i), gridCoordinates[i][0], 1e-12) << "Failed for i " << i;
            EXPECT_NEAR(grid_const.coord(1, i), gridCoordinates[i][0], 1e-12) << "Failed for i " << i;
    }   
}

TEST_F(NDGridTest, CoordVectorTest) {  // and newCoordArray
    std::vector< std::vector<double> > gridCoordinates(2);
    gridCoordinates[0] = std::vector<double>(2, 0.);
    gridCoordinates[1] = std::vector<double>(3, 1.);
    gridCoordinates[1][2] = 9.;	
    interpolation::NDGrid grid(gridCoordinates);
    for (int i = 0; i < 2; ++i) {
        std::vector<double> coords_v = grid.coordVector(i);
        double* coords_a = grid.newCoordArray(i);
        ASSERT_EQ(coords_v.size(), gridCoordinates[i].size());
        for (size_t j = 0; j < gridCoordinates[i].size(); j++) {
            EXPECT_NEAR(coords_v[j], gridCoordinates[i][j], 1e-12);
            EXPECT_NEAR(coords_a[j], gridCoordinates[i][j], 1e-12);
        }
        delete coords_a;
    }
}
//} // namespace ndgridtest

