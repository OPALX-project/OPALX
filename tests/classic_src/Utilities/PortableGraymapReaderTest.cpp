#include "gtest/gtest.h"
#include "OpalSourcePath.h"

#include "Utilities/PortableGraymapReader.h"
#include "opal_test_utilities/SilenceTest.h"

#include <iostream>

TEST(PGMReaderTest, SimpleAscii8Test) {
    OpalTestUtilities::SilenceTest silencer;

    std::string opalSourcePath = OPAL_SOURCE_DIR;
    std::string pathToBitmapFile = opalSourcePath + "/tests/classic_src/Utilities/Untitled8.pgm";
    PortableGraymapReader reader(pathToBitmapFile);

    float pixel = reader.getPixel(199, 299);
    float expected = 103 / (float)255;
    ASSERT_NEAR(expected, pixel, 1e-7);

    pixel = reader.getPixel(200, 299);
    expected = 68 / (float)255;
    ASSERT_NEAR(expected, pixel, 1e-7);

    pixel = reader.getPixel(200, 300);
    expected = 219 / (float)255;
    ASSERT_NEAR(expected, pixel, 1e-7);

    pixel = reader.getPixel(199, 300);
    expected = 155 / (float)255;
    ASSERT_NEAR(expected, pixel, 1e-7);
}

TEST(PGMReaderTest, SimpleAscii16Test) {
    OpalTestUtilities::SilenceTest silencer;

    std::string opalSourcePath = OPAL_SOURCE_DIR;
    std::string pathToBitmapFile = opalSourcePath + "/tests/classic_src/Utilities/Untitled16.pgm";
    PortableGraymapReader reader(pathToBitmapFile);

    float pixel = reader.getPixel(199, 299);
    float expected = 26610 / (float)65535;
    ASSERT_NEAR(expected, pixel, 1e-7);

    pixel = reader.getPixel(200, 299);
    expected = 17689 / (float)65535;
    ASSERT_NEAR(expected, pixel, 1e-7);

    pixel = reader.getPixel(200, 300);
    expected = 56453 / (float)65535;
    ASSERT_NEAR(expected, pixel, 1e-7);

    pixel = reader.getPixel(199, 300);
    expected = 39879 / (float)65535;
    ASSERT_NEAR(expected, pixel, 1e-7);
}

TEST(PGMReaderTest, SimpleBinary8Test) {
    OpalTestUtilities::SilenceTest silencer;

    std::string opalSourcePath = OPAL_SOURCE_DIR;
    std::string pathToBitmapFile = opalSourcePath + "/tests/classic_src/Utilities/Untitled_binary8.pgm";
    PortableGraymapReader reader(pathToBitmapFile);

    float pixel = reader.getPixel(199, 299);
    float expected = 103 / (float)255;
    // std::cout << pixel << "\t" << expected << std::endl;
    ASSERT_NEAR(expected, pixel, 1e-7);

    pixel = reader.getPixel(200, 299);
    expected = 68 / (float)255;
    // std::cout << pixel << "\t" << expected << std::endl;
    ASSERT_NEAR(expected, pixel, 1e-7);

    pixel = reader.getPixel(200, 300);
    expected = 219 / (float)255;
    // std::cout << pixel << "\t" << expected << std::endl;
    ASSERT_NEAR(expected, pixel, 1e-7);

    pixel = reader.getPixel(199, 300);
    expected = 155 / (float)255;
    // std::cout << pixel << "\t" << expected << std::endl;
    ASSERT_NEAR(expected, pixel, 1e-7);
}

TEST(PGMReaderTest, SimpleBinary16Test) {
    OpalTestUtilities::SilenceTest silencer;

    std::string opalSourcePath = OPAL_SOURCE_DIR;
    std::string pathToBitmapFile = opalSourcePath + "/tests/classic_src/Utilities/Untitled_binary16.pgm";
    PortableGraymapReader reader(pathToBitmapFile);

    float pixel = reader.getPixel(199, 299);
    float expected = 26610 / (float)65535;
    ASSERT_NEAR(expected, pixel, 1e-7);

    pixel = reader.getPixel(200, 299);
    expected = 17689 / (float)65535;
    ASSERT_NEAR(expected, pixel, 1e-7);

    pixel = reader.getPixel(200, 300);
    expected = 56453 / (float)65535;
    ASSERT_NEAR(expected, pixel, 1e-7);

    pixel = reader.getPixel(199, 300);
    expected = 39879 / (float)65535;
    ASSERT_NEAR(expected, pixel, 1e-7);
}