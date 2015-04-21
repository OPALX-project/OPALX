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

#include <iostream>
#include <fstream>
#include <string>
#include <cstring>

std::string input = "OPTION, ECHO=TRUE;\n"
    "OPTION, CZERO=FALSE;\n"
    "TITLE, STRING=\"gauss distribution unit test\";\n"
    "DIST1: DISTRIBUTION, DISTRIBUTION = \"GAUSS\", \n"
    "SIGMAX = 3.914e-6, SIGMAY = 6.239e-6, SIGMAZ = 2.363e-6, \n"
    "SIGMAPX = 0.6396, SIGMAPY = 0.3859, SIGMAPZ = 0.8944, \n"
    "R = {-0.6486e-3, 0, 0, 0.4542e-6, 1.362e-3, 0, 0, 0.7265e-3, -0.2685, 1.198e-3, 0, 0, 0, 0, 0.1752e-3}, \n"
    "EKIN = 0.63, \n"
    "EMITTED = FALSE;\n";

TEST(GaussTest, FullSigmaTest) {

    char inputFileName[] = "GaussDistributionTest.in";

    int narg = 8;
    char exe_name[] = "opal_unit_tests";
    char commlib[] = "--commlib";
    char mpicomm[] = "mpi";
    char info[] = "--info";
    char info0[] = "0";
    char warn[] = "--warn";
    char warn0[] = "0";

    char **arg = new char*[8];
    arg[0] = exe_name;
    arg[1] = inputFileName;
    arg[2] = commlib;
    arg[3] = mpicomm;
    arg[4] = info;
    arg[5] = info0;
    arg[6] = warn;
    arg[7] = warn0;

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
        throw new OpalException("FullSigmaTest", "Could not read string");
    }

    OpalParser *parser = new OpalParser();
    if (is) {
        try {
            parser->run(is);
        } catch (...) {
            throw new OpalException("FullSigmaTest", "Could not parse input");
        }
    }

    Object *distObj;
    try {
        distObj = OPAL->find("DIST1");
    } catch(...) {
        distObj = 0;
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

        std::cout << std::setprecision(4)
                  << std::setw(11) << R11
                  << std::setw(11) << R21
                  << std::setw(11) << R22
                  << std::endl
                  << std::setw(11) << R51
                  << std::setw(11) << R52
                  << std::setw(11) << R61
                  << std::setw(11) << R62
                  << std::endl;

        //EXPECT_NEAR();
    }

    OpalData::deleteInstance();
    delete parser;
    delete gmsg;
    delete ippl;
    delete[] arg;
}