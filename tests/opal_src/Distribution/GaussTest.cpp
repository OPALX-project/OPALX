#include "gtest/gtest.h"

#include "opal.h"

Ippl *ippl;
Inform *gmsg;

#include "AbstractObjects/OpalData.h"
#include "OpalConfigure/Configure.h"
#include "Utilities/OpalException.h"
#include "Distribution/Distribution.h"
#include "OpalParser/OpalParser.h"

#include "Parser/FileStream.h"
#include "Physics/Physics.h"

#include "gsl/gsl_statistics_double.h"

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sstream>

namespace {
    std::string burnAfterReading(std::ostringstream &ostr) {
        std::string returnValue = ostr.str();
        ostr.str("");

        return returnValue;
    }
}

TEST(GaussTest, FullSigmaTest1) {

    std::streambuf *defaultCout;
    std::streambuf *defaultCerr;
    std::ostringstream debugOutput;
    defaultCout = std::cout.rdbuf(debugOutput.rdbuf());
    defaultCerr = std::cerr.rdbuf(debugOutput.rdbuf());

    char inputFileName[] = "GaussDistributionTest.in";
    std::string input = "OPTION, ECHO=FALSE;\n"
        "OPTION, CZERO=FALSE;\n"
        "TITLE, STRING=\"gauss distribution unit test\";\n"
        "DIST1: DISTRIBUTION, DISTRIBUTION = \"GAUSS\", \n"
        "SIGMAX = 1.978e-3, SIGMAY = 2.498e-3, SIGMAZ = 1.537e-3, \n"
        "SIGMAPX = 0.7998, SIGMAPY = 0.6212, SIGMAPZ = 0.9457, \n"
        "R = {-0.40993, 0, 0, 0.14935, 0.72795,   0, 0, 0.59095, -0.3550,   0.77208, 0, 0,   0, 0,  0.12051}, \n"
        "EKIN = 0.63, \n"
        "EMITTED = FALSE;\n";

    int narg = 7;
    char exe_name[] = "opal_unit_tests";
    char commlib[] = "--nocomminit";
    char info[] = "--info";
    char info0[] = "0";
    char warn[] = "--warn";
    char warn0[] = "0";

    char **arg = new char*[7];
    arg[0] = exe_name;
    arg[1] = inputFileName;
    arg[2] = commlib;
    arg[3] = info;
    arg[4] = info0;
    arg[5] = warn;
    arg[6] = warn0;

    if (!ippl)
        ippl = new Ippl(narg, arg, Ippl::KEEP, MPI_COMM_WORLD);
    gmsg = new Inform("OPAL ");

    std::ofstream inputFile(inputFileName);
    inputFile << input << std::endl;
    inputFile.close();

    OpalData *OPAL = OpalData::getInstance();
    Configure::configure();
    OPAL->storeInputFn(inputFileName);

    FileStream *is = 0;
    try {
        is = new FileStream(inputFileName);
    } catch(...) {
        is = 0;
        std::cout.rdbuf(defaultCout);
        std::cerr.rdbuf(defaultCerr);
        throw new OpalException("FullSigmaTest", "Could not read string");
    }

    OpalParser *parser = new OpalParser();
    if (is) {
        try {
            parser->run(is);
        } catch (...) {
            std::cout.rdbuf(defaultCout);;
            std::cerr.rdbuf(defaultCerr);;
            throw new OpalException("FullSigmaTest", "Could not parse input");
        }
    }

    Object *distObj;
    try {
        distObj = OPAL->find("DIST1");
    } catch(...) {
        distObj = 0;
        std::cout.rdbuf(defaultCout);;
        std::cerr.rdbuf(defaultCerr);;
        throw new OpalException("FullSigmaTest", "Could not find distribution");
    }

    if (distObj) {
        Distribution *dist = dynamic_cast<Distribution*>(distObj);

        dist->SetDistType();
        dist->CheckIfEmitted();

        size_t numParticles = 1000000;
        dist->Create(numParticles, Physics::m_p);


        double R11 = gsl_stats_variance(&(dist->xDist_m[0]), 1, dist->xDist_m.size()) * 1e6;
        double R21 = gsl_stats_covariance(&(dist->xDist_m[0]), 1, &(dist->pxDist_m[0]), 1, dist->xDist_m.size()) * 1e3;
        double R22 = gsl_stats_variance(&(dist->pxDist_m[0]), 1, dist->pxDist_m.size());

        double R51 = gsl_stats_covariance(&(dist->xDist_m[0]), 1, &(dist->tOrZDist_m[0]), 1, dist->xDist_m.size()) * 1e6;
        double R52 = gsl_stats_covariance(&(dist->pxDist_m[0]), 1, &(dist->tOrZDist_m[0]), 1, dist->pxDist_m.size()) * 1e3;
        double R61 = gsl_stats_covariance(&(dist->xDist_m[0]), 1, &(dist->pzDist_m[0]), 1, dist->xDist_m.size()) * 1e3;
        double R62 = gsl_stats_covariance(&(dist->pxDist_m[0]), 1, &(dist->pzDist_m[0]), 1, dist->pxDist_m.size());

        const double expectedR11 = 3.914;
        const double expectedR21 = -0.6486;
        const double expectedR22 = 0.6396;
        const double expectedR51 = 0.4542;
        const double expectedR52 = 0.7265;
        const double expectedR61 = 1.362;
        const double expectedR62 = -0.2685;

        EXPECT_LT(std::abs(expectedR11 - R11),  0.05 * expectedR11) << ::burnAfterReading(debugOutput);
        EXPECT_LT(std::abs(expectedR21 - R21), -0.05 * expectedR21) << ::burnAfterReading(debugOutput);
        EXPECT_LT(std::abs(expectedR22 - R22),  0.05 * expectedR22) << ::burnAfterReading(debugOutput);
        EXPECT_LT(std::abs(expectedR51 - R51),  0.05 * expectedR51) << ::burnAfterReading(debugOutput);
        EXPECT_LT(std::abs(expectedR52 - R52),  0.05 * expectedR52) << ::burnAfterReading(debugOutput);
        EXPECT_LT(std::abs(expectedR61 - R61),  0.05 * expectedR61) << ::burnAfterReading(debugOutput);
        EXPECT_LT(std::abs(expectedR62 - R62), -0.05 * expectedR62) << ::burnAfterReading(debugOutput);
    }

    std::cout.rdbuf(defaultCout);;
    std::cerr.rdbuf(defaultCerr);;

    OpalData::deleteInstance();
    delete parser;
    delete gmsg;
    //    delete ippl;
    delete[] arg;

    std::remove(inputFileName);
}

TEST(GaussTest, FullSigmaTest2) {

    std::streambuf *defaultCout;
    std::streambuf *defaultCerr;
    std::ostringstream debugOutput;
    defaultCout = std::cout.rdbuf(debugOutput.rdbuf());
    defaultCerr = std::cerr.rdbuf(debugOutput.rdbuf());

    char inputFileName[] = "GaussDistributionTest.in";
    std::string input = "OPTION, ECHO=FALSE;\n"
        "OPTION, CZERO=FALSE;\n"
        "TITLE, STRING=\"gauss distribution unit test\";\n"
        "DIST1: DISTRIBUTION, DISTRIBUTION = \"GAUSS\", \n"
        "SIGMAX = 1.978e-3, SIGMAY = 2.498e-3, SIGMAZ = 1.537e-3, \n"
        "SIGMAPX = 0.7998, SIGMAPY = 0.6212, SIGMAPZ = 0.9457, \n"
        "CORRX= -0.40993, CORRY=0.77208, CORRZ=0.12051, \n"
        "R51=0.14935, R52=0.59095, R61=0.72795, R62=-0.3550, \n"
        "EKIN = 0.63, \n"
        "EMITTED = FALSE;\n";


    int narg = 7;
    char exe_name[] = "opal_unit_tests";
    char commlib[] = "--nocomminit";
    char info[] = "--info";
    char info0[] = "0";
    char warn[] = "--warn";
    char warn0[] = "0";

    char **arg = new char*[7];
    arg[0] = exe_name;
    arg[1] = inputFileName;
    arg[2] = commlib;
    arg[3] = info;
    arg[4] = info0;
    arg[5] = warn;
    arg[6] = warn0;

    if (!ippl)
        ippl = new Ippl(narg, arg, Ippl::KEEP, MPI_COMM_WORLD);
    gmsg = new Inform("OPAL ");

    std::ofstream inputFile(inputFileName);
    inputFile << input << std::endl;
    inputFile.close();

    OpalData *OPAL = OpalData::getInstance();
    Configure::configure();
    OPAL->storeInputFn(inputFileName);

    FileStream *is = 0;
    try {
        is = new FileStream(inputFileName);
    } catch(...) {
        is = 0;
        std::cout.rdbuf(defaultCout);;
        std::cerr.rdbuf(defaultCerr);;
        throw new OpalException("FullSigmaTest", "Could not read string");
    }

    OpalParser *parser = new OpalParser();
    if (is) {
        try {
            parser->run(is);
        } catch (...) {
            std::cout.rdbuf(defaultCout);;
            std::cerr.rdbuf(defaultCerr);;
            throw new OpalException("FullSigmaTest", "Could not parse input");
        }
    }

    Object *distObj;
    try {
        distObj = OPAL->find("DIST1");
    } catch(...) {
        distObj = 0;
        std::cout.rdbuf(defaultCout);;
        std::cerr.rdbuf(defaultCerr);;
        throw new OpalException("FullSigmaTest", "Could not find distribution");
    }

    if (distObj) {
        Distribution *dist = dynamic_cast<Distribution*>(distObj);

        dist->SetDistType();
        dist->CheckIfEmitted();

        size_t numParticles = 1000000;
        dist->Create(numParticles, Physics::m_p);


        double R11 = gsl_stats_variance(&(dist->xDist_m[0]), 1, dist->xDist_m.size()) * 1e6;
        double R21 = gsl_stats_covariance(&(dist->xDist_m[0]), 1, &(dist->pxDist_m[0]), 1, dist->xDist_m.size()) * 1e3;
        double R22 = gsl_stats_variance(&(dist->pxDist_m[0]), 1, dist->pxDist_m.size());

        double R51 = gsl_stats_covariance(&(dist->xDist_m[0]), 1, &(dist->tOrZDist_m[0]), 1, dist->xDist_m.size()) * 1e6;
        double R52 = gsl_stats_covariance(&(dist->pxDist_m[0]), 1, &(dist->tOrZDist_m[0]), 1, dist->pxDist_m.size()) * 1e3;
        double R61 = gsl_stats_covariance(&(dist->xDist_m[0]), 1, &(dist->pzDist_m[0]), 1, dist->xDist_m.size()) * 1e3;
        double R62 = gsl_stats_covariance(&(dist->pxDist_m[0]), 1, &(dist->pzDist_m[0]), 1, dist->pxDist_m.size());

        const double expectedR11 = 3.914;
        const double expectedR21 = -0.6486;
        const double expectedR22 = 0.6396;
        const double expectedR51 = 0.4542;
        const double expectedR52 = 0.7265;
        const double expectedR61 = 1.362;
        const double expectedR62 = -0.2685;

        EXPECT_LT(std::abs(expectedR11 - R11),  0.05 * expectedR11) << ::burnAfterReading(debugOutput);
        EXPECT_LT(std::abs(expectedR21 - R21), -0.05 * expectedR21) << ::burnAfterReading(debugOutput);
        EXPECT_LT(std::abs(expectedR22 - R22),  0.05 * expectedR22) << ::burnAfterReading(debugOutput);
        EXPECT_LT(std::abs(expectedR51 - R51),  0.05 * expectedR51) << ::burnAfterReading(debugOutput);
        EXPECT_LT(std::abs(expectedR52 - R52),  0.05 * expectedR52) << ::burnAfterReading(debugOutput);
        EXPECT_LT(std::abs(expectedR61 - R61),  0.05 * expectedR61) << ::burnAfterReading(debugOutput);
        EXPECT_LT(std::abs(expectedR62 - R62), -0.05 * expectedR62) << ::burnAfterReading(debugOutput);
    }

    std::cout.rdbuf(defaultCout);;
    std::cerr.rdbuf(defaultCerr);;

    OpalData::deleteInstance();
    delete parser;
    delete gmsg;
    //    delete ippl;
    delete[] arg;

    std::remove(inputFileName);
}