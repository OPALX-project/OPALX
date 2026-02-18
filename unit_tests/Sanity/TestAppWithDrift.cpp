#include "gtest/gtest.h"

extern "C" {
#include "H5hut.h"
}
#include "AbstractObjects/OpalData.h"
#include "Fields/Fieldmap.h"
#include "OpalConfigure/Configure.h"
#include "OpalParser/FileStream.h"
#include "OpalParser/OpalParser.h"
#include "Utilities/Util.h"
#include "GSLErrorHandling.h"
#include "Utility/Inform.h"
#include "Utility/IpplTimings.h"
#include "Utilities/GSLCompat.h"
#include <filesystem>
#include <system_error>
#include <iostream>
#include <string>

class TemporaryDirectory {
public:
    TemporaryDirectory() {
        const auto base = std::filesystem::temp_directory_path();
        path_ = base / randomName();
        std::filesystem::create_directory(path_);
    }

    ~TemporaryDirectory() {
        std::error_code ec;
        std::filesystem::remove_all(path_, ec);
    }

    const std::filesystem::path& path() const { return path_; }

private:
    static std::string randomName() {
        static std::mt19937_64 rng{std::random_device{}()};
        static std::uniform_int_distribution<unsigned long long> dist;
        return "testdir_" + std::to_string(dist(rng));
    }

    std::filesystem::path path_;
};

TEST(TestAppWithDrift, Drift1) {
    TemporaryDirectory tmpDir;
    // Create the input file
    auto fileName = tmpDir.path() / "Drift-3.in";
    {
        std::ofstream ofs(fileName);
        ofs << R"====(
/*
    Drinft-3.in produces a meaningful space charge drift that can be compared to
    OPAL.
*/

OPTION, PSDUMPFREQ      = 10;     // 6d data written every 10th time step (h5).
OPTION, STATDUMPFREQ    = 1;     // Beam Stats written every time step (stat).
OPTION, BOUNDPDESTROYFQ = 10;    // Delete lost particles, if any out of 10 \sigma
OPTION, AUTOPHASE       = 4;     // Autophase is on, and phase of max energy
                                 // gain will be found automatically for cavities.
OPTION, VERSION=10900;

Title, string="Test simple Gaussian distribution with space charge";

Value,{OPALVERSION};

// ----------------------------------------------------------------------------
// Global Parameters

REAL N                   = 32;
REAL rf_freq             = 1.3e9;    // RF frequency. (Hz)
REAL n_particles         = N*N*N;      // Number of particles in simulation.
REAL beam_bunch_charge   = 1e-11;     // Charge of bunch. (C)

// ----------------------------------------------------------------------------
// Initial Momentum Calculation

REAL Edes    = 1e-9; //initial energy in GeV
REAL gamma   = (Edes+EMASS)/EMASS;
REAL beta    = sqrt(1-(1/gamma^2));
REAL P0      = gamma*beta*EMASS;    //inital z momentum

//Printing initial energy and momentum to terminal output.
value, {Edes, P0};

D1: DRIFT, L = 1.25, ELEMEDGE = 0.0;
myLine: Line = (D1);

BEAM1:  BEAM, PARTICLE = ELECTRON, pc = P0, NPART = n_particles,
        BFREQ = rf_freq, BCURRENT = beam_bunch_charge * rf_freq * 1e6, CHARGE = -1;

FS1: Fieldsolver, TYPE=OPEN,
    NX=N,
    NY=N,
    NZ=N,
    PARFFTX = true,
    PARFFTY = true,
    PARFFTZ = false,
    BCFFTX=OPEN,
    BCFFTY=OPEN,
    BCFFTZ=OPEN, BBOXINCR = 5, GREENSF = STANDARD;

Dist1: DISTRIBUTION, TYPE=GAUSS,
    SIGMAX=1e-3,
    SIGMAY=1e-3,
    SIGMAZ=1e-3,
    SIGMAPX=1e-9,
    SIGMAPY=1e-9,
    SIGMAPZ=1e-9;

TRACK, LINE = myLine, BEAM = BEAM1, MAXSTEPS = 20, DT = {5e-11}, ZSTOP=1.25;
 RUN, METHOD = "PARALLEL", BEAM = BEAM1, FIELDSOLVER = FS1, DISTRIBUTION = Dist1;ENDTRACK;
QUIT;
)====";
    }
    // Run Opal
    // This is a cut down version of what is in Main.cpp.  It would be best
    // if Main.cpp could be refactored such that what is in the main()
    // function could be linked to from here rather than code being repeated.
    std::string appName = "opalx";
    const char* args[] = {appName.c_str(), fileName.c_str(), nullptr};
    int argCount = 2;
    ippl::initialize(argCount, const_cast<char**>(args));
    gmsg = new Inform("OPAL-X");
    H5SetVerbosityLevel(1); // 65535);
    gsl_set_error_handler(&handleGSLErrors);
    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("mainTimer");
    IpplTimings::startTimer(mainTimer);
    const OpalParser parser;
    std::cout.precision(16);
    std::cout.setf(std::ios::scientific, std::ios::floatfield);
    std::cerr.precision(16);
    std::cerr.setf(std::ios::scientific, std::ios::floatfield);
    OpalData* opal = OpalData::getInstance();
    if(ippl::Comm->rank() == 0) {
        if(!std::filesystem::exists(opal->getAuxiliaryOutputDirectory())) {
            std::error_code error_code;
            std::filesystem::create_directory(opal->getAuxiliaryOutputDirectory(), error_code);
        }
    }
    ippl::Comm->barrier();
    Configure::configure();
    opal->storeInputFn(fileName);
    FileStream* is;
    is = new FileStream(fileName);
    parser.run(is);
    IpplTimings::stopTimer(mainTimer);
    IpplTimings::print();
    IpplTimings::print(
            std::string("timing.dat"), OpalData::getInstance()->getProblemCharacteristicValues());
    ippl::Comm->barrier();
    Fieldmap::clearDictionary();
    delete gmsg;
    ippl::finalize();
    // Do some crude validation of the resulting files
    EXPECT_TRUE(std::filesystem::exists(tmpDir.path() / "Drift-3.h5"));
    EXPECT_TRUE(std::filesystem::exists(tmpDir.path() / "Drift-3.stat"));
}
