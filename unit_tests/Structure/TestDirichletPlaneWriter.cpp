//
// Unit tests for DirichletPlaneWriter.
//

#include <gtest/gtest.h>

#include "Ippl.h"
#include "Structure/DirichletPlaneWriter.h"
#include "Utilities/OpalException.h"

#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

class DirichletPlaneWriterTest : public ::testing::Test {
protected:
    static void SetUpTestSuite() {
        int argc = 0;
        char** argv = nullptr;
        ippl::initialize(argc, argv);
    }

    static void TearDownTestSuite() {
        ippl::finalize();
    }

    static std::string readFile(const std::string& path) {
        std::ifstream f(path);
        std::stringstream ss;
        ss << f.rdbuf();
        return ss.str();
    }
};

TEST_F(DirichletPlaneWriterTest, ConstructorCreatesOutputDirectory) {
    const std::string dir = "dirichlet_writer_test_dir";
    std::filesystem::remove_all(dir);

    DirichletPlaneWriter writer(dir);

    EXPECT_TRUE(std::filesystem::exists(dir));
    EXPECT_TRUE(std::filesystem::is_directory(dir));
}

TEST_F(DirichletPlaneWriterTest, WritePlaneCreatesAsciiDumpWithMetadata) {
    const std::string dir = "dirichlet_writer_dump_dir";
    std::filesystem::remove_all(dir);

    DirichletPlaneWriter writer(dir);

    const std::size_t nx = 2;
    const std::size_t ny = 3;
    const std::vector<double> x = {0.0, 0.1};
    const std::vector<double> y = {-0.2, 0.0, 0.2};
    const std::vector<double> phi = {
            0.0, 1.0, 2.0,
            3.0, 4.0, 5.0};

    writer.writePlane(
            /*step=*/7,
            /*time=*/1.5e-9,
            /*zPlane=*/0.25,
            x,
            y,
            phi,
            nx,
            ny,
            "legacy");

    const std::string expectedFile = dir + "/phi_dirichlet_legacy_step-00000007.dat";
    ASSERT_TRUE(std::filesystem::exists(expectedFile));

    const std::string content = readFile(expectedFile);
    EXPECT_NE(content.find("# Dirichlet-plane potential dump"), std::string::npos);
    EXPECT_NE(content.find("step=7"), std::string::npos);
    EXPECT_NE(content.find("zPlane=0.25"), std::string::npos);
    EXPECT_NE(content.find("columns: i j x[m] y[m] phi[V]"), std::string::npos);
}

TEST_F(DirichletPlaneWriterTest, InvalidDimensionsThrow) {
    const std::string dir = "dirichlet_writer_invalid_dims";
    std::filesystem::remove_all(dir);

    DirichletPlaneWriter writer(dir);

    const std::vector<double> x = {0.0};
    const std::vector<double> y = {0.0};
    const std::vector<double> phi = {1.0};

    EXPECT_THROW(
            writer.writePlane(
                    /*step=*/0,
                    /*time=*/0.0,
                    /*zPlane=*/0.0,
                    x,
                    y,
                    phi,
                    /*nx=*/2,
                    /*ny=*/1,
                    "legacy"),
            OpalException);
}
