//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Structure/DataSink.h"

#include "OPALconfig.h"
#include "OPALrevision.h"
#include "Algorithms/bet/EnvelopeBunch.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"
#include "Utilities/Options.h"
#include "Fields/Fieldmap.h"
#include "Structure/BoundaryGeometry.h"
#include "Structure/H5PartWrapper.h"
#include "Structure/H5PartWrapperForPS.h"
#include "Utilities/Timer.h"

#include "H5hut.h"

#include <boost/filesystem.hpp>
#include <boost/regex.hpp>

#include <queue>
#include <sstream>

extern Inform *gmsg;

using namespace std;

//using Physics::m_e;

// backward compatibility with 1.99.5
#if defined(H5_O_FLUSHSTEP)
#define H5_FLUSH_STEP H5_O_FLUSHSTEP
#endif

DataSink::DataSink() :
    lossWrCounter_m(0),
    doHDF5_m(true),
    h5wrapper_m(NULL),
    H5call_m(0)
{ }

DataSink::DataSink(H5PartWrapper *h5wrapper, int restartStep):
    lossWrCounter_m(0),
    h5wrapper_m(h5wrapper),
    H5call_m(0)
{
    namespace fs = boost::filesystem;

    doHDF5_m = Options::enableHDF5;
    if (!doHDF5_m) {
        throw OpalException("DataSink::DataSink(int)",
                            "Can not restart when HDF5 is disabled");
    }

    H5PartTimer_m = IpplTimings::getTimer("H5PartTimer");
    StatMarkerTimer_m = IpplTimings::getTimer("StatMarkerTimer");

    /// Set flags whether SDDS file exists or have to write header first.
    firstWriteToStat_m = true;

    string fn = OpalData::getInstance()->getInputBasename();

    statFileName_m = fn + string(".stat");
    lBalFileName_m = fn + string(".lbal");

    unsigned int linesToRewind = 0;
    if (fs::exists(statFileName_m)) {
        firstWriteToStat_m = false;
        INFOMSG("Appending statistical data to existing data file: " << statFileName_m << endl);
        double spos = h5wrapper->getLastPosition();
        linesToRewind = rewindSDDStoSPos(spos);
    } else {
        INFOMSG("Creating new file for statistical data: " << statFileName_m << endl);
    }

    if (fs::exists(lBalFileName_m)) {
        INFOMSG("Appending load balance data to existing data file: " << lBalFileName_m << endl);
        rewindLinesLBal(linesToRewind);
    } else {
        INFOMSG("Creating new file for load balance data: " << lBalFileName_m << endl);
    }
}

DataSink::DataSink(H5PartWrapper *h5wrapper):
    lossWrCounter_m(0),
    h5wrapper_m(h5wrapper),
    H5call_m(0)
{
    /// Constructor steps:
    /// Get timers from IPPL.
    H5PartTimer_m = IpplTimings::getTimer("H5PartTimer");
    StatMarkerTimer_m = IpplTimings::getTimer("StatMarkerTimer");

    /// Set file write flags to true. These will be set to false after first
    /// write operation.
    firstWriteToStat_m = true;
    firstWriteH5Surface_m = true;
    /// Define file names.
    string fn = OpalData::getInstance()->getInputBasename();
    surfaceLossFileName_m = fn + string(".SurfaceLoss.h5");
    statFileName_m = fn + string(".stat");
    lBalFileName_m = fn + string(".lbal");

    doHDF5_m = Options::enableHDF5;

    h5wrapper_m->writeHeader();
}

DataSink::~DataSink() {
    h5wrapper_m = NULL;
}

void DataSink::storeCavityInformation() {
    if (!doHDF5_m) return;

    h5wrapper_m->storeCavityInformation();
}

void DataSink::writePhaseSpace(PartBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail) {

    if (!doHDF5_m) return;

    IpplTimings::startTimer(H5PartTimer_m);
    std::map<std::string, double> additionalAttributes = {
        std::make_pair("sposHead", sposHead),
        std::make_pair("sposRef", sposRef),
        std::make_pair("sposTail", sposTail),
        std::make_pair("B-head_x", FDext[0](0)),
        std::make_pair("B-head_z", FDext[0](1)),
        std::make_pair("B-head_y", FDext[0](2)),
        std::make_pair("E-head_x", FDext[1](0)),
        std::make_pair("E-head_z", FDext[1](1)),
        std::make_pair("E-head_y", FDext[1](2)),
        std::make_pair("B-ref_x", FDext[2](0)),
        std::make_pair("B-ref_z", FDext[2](1)),
        std::make_pair("B-ref_y", FDext[2](2)),
        std::make_pair("E-ref_x", FDext[3](0)),
        std::make_pair("E-ref_z", FDext[3](1)),
        std::make_pair("E-ref_y", FDext[3](2)),
        std::make_pair("B-tail_x", FDext[4](0)),
        std::make_pair("B-tail_z", FDext[4](1)),
        std::make_pair("B-tail_y", FDext[4](2)),
        std::make_pair("E-tail_x", FDext[5](0)),
        std::make_pair("E-tail_z", FDext[5](1)),
        std::make_pair("E-tail_y", FDext[5](2))};

    h5wrapper_m->writeStep(beam, additionalAttributes);
    IpplTimings::stopTimer(H5PartTimer_m);

    return;
}



int DataSink::writePhaseSpace_cycl(PartBunch &beam, Vector_t FDext[], double meanEnergy,
                                   double refPr, double refPt, double refPz,
                                   double refR, double refTheta, double refZ,
                                   double azimuth, double elevation, bool local) {

    if (!doHDF5_m) return -1;
    if (beam.getLocalNum() < 3) return -1; // in single particle mode and tune calculation (2 particles) we do not need h5 data

    IpplTimings::startTimer(H5PartTimer_m);
    std::map<std::string, double> additionalAttributes = {
        std::make_pair("REFPR", refPr),
        std::make_pair("REFPT", refPt),
        std::make_pair("REFPZ", refPz),
        std::make_pair("REFR", refR),
        std::make_pair("REFTHETA", refTheta),
        std::make_pair("REFZ", refZ),
        std::make_pair("AZIMUTH", azimuth),
        std::make_pair("ELEVATION", elevation),
        std::make_pair("B-head_x", FDext[0](0)),
        std::make_pair("B-head_z", FDext[0](1)),
        std::make_pair("B-head_y", FDext[0](2)),
        std::make_pair("E-head_x", FDext[1](0)),
        std::make_pair("E-head_z", FDext[1](1)),
        std::make_pair("E-head_y", FDext[1](2)),
        std::make_pair("B-ref_x", FDext[2](0)),
        std::make_pair("B-ref_z", FDext[2](1)),
        std::make_pair("B-ref_y", FDext[2](2)),
        std::make_pair("E-ref_x", FDext[3](0)),
        std::make_pair("E-ref_z", FDext[3](1)),
        std::make_pair("E-ref_y", FDext[3](2)),
        std::make_pair("B-tail_x", FDext[4](0)),
        std::make_pair("B-tail_z", FDext[4](1)),
        std::make_pair("B-tail_y", FDext[4](2)),
        std::make_pair("E-tail_x", FDext[5](0)),
        std::make_pair("E-tail_z", FDext[5](1)),
        std::make_pair("E-tail_y", FDext[5](2))};

    h5wrapper_m->writeStep(beam, additionalAttributes);
    IpplTimings::stopTimer(H5PartTimer_m);

    ++ H5call_m;
    return H5call_m - 1;
}

void DataSink::writePhaseSpaceEnvelope(EnvelopeBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail) {

    if (!doHDF5_m) return;

    IpplTimings::startTimer(H5PartTimer_m);
    std::map<std::string, double> additionalAttributes = {
        std::make_pair("sposHead", sposHead),
        std::make_pair("sposRef", sposRef),
        std::make_pair("sposTail", sposTail),
        std::make_pair("B-head_x", FDext[0](0)),
        std::make_pair("B-head_z", FDext[0](1)),
        std::make_pair("B-head_y", FDext[0](2)),
        std::make_pair("E-head_x", FDext[1](0)),
        std::make_pair("E-head_z", FDext[1](1)),
        std::make_pair("E-head_y", FDext[1](2)),
        std::make_pair("B-ref_x", FDext[2](0)),
        std::make_pair("B-ref_z", FDext[2](1)),
        std::make_pair("B-ref_y", FDext[2](2)),
        std::make_pair("E-ref_x", FDext[3](0)),
        std::make_pair("E-ref_z", FDext[3](1)),
        std::make_pair("E-ref_y", FDext[3](2)),
        std::make_pair("B-tail_x", FDext[4](0)),
        std::make_pair("B-tail_z", FDext[4](1)),
        std::make_pair("B-tail_y", FDext[4](2)),
        std::make_pair("E-tail_x", FDext[5](0)),
        std::make_pair("E-tail_z", FDext[5](1)),
        std::make_pair("E-tail_y", FDext[5](2))};

    h5wrapper_m->writeStep(beam, additionalAttributes);
    IpplTimings::stopTimer(H5PartTimer_m);
}

void DataSink::stashPhaseSpaceEnvelope(EnvelopeBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail) {

    if (!doHDF5_m) return;

    /// Start timer.
    IpplTimings::startTimer(H5PartTimer_m);

    static_cast<H5PartWrapperForPS*>(h5wrapper_m)->stashPhaseSpaceEnvelope(beam,
                                                                           FDext,
                                                                           sposHead,
                                                                           sposRef,
                                                                           sposTail);
    H5call_m++;

    /// %Stop timer.
    IpplTimings::stopTimer(H5PartTimer_m);
}

void DataSink::dumpStashedPhaseSpaceEnvelope() {

    if (!doHDF5_m) return;

    static_cast<H5PartWrapperForPS*>(h5wrapper_m)->dumpStashedPhaseSpaceEnvelope();

    /// %Stop timer.
    IpplTimings::stopTimer(H5PartTimer_m);
}

void DataSink::writeStatData(PartBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail, double E) {
    doWriteStatData(beam, FDext, sposHead, sposRef, sposTail, E, std::vector<std::pair<std::string, unsigned int> >());
}

void DataSink::writeStatData(PartBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail, const std::vector<std::pair<std::string, unsigned int> >& losses) {
    doWriteStatData(beam, FDext, sposHead, sposRef, sposTail, beam.get_meanEnergy(), losses);
}


void DataSink::doWriteStatData(PartBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail, double E, const std::vector<std::pair<std::string, unsigned int> > &losses) {

    /// Function steps:

    /// Start timer.
    IpplTimings::startTimer(StatMarkerTimer_m);

    /// Set width of write fields in output files.
    unsigned int pwi = 10;

    /// Calculate beam statistics and gather load balance statistics.
//    beam.calcBeamParameters_cycl();
    beam.gatherLoadBalanceStatistics();

    size_t npOutside = 0;
    if (Options::beamHaloBoundary>0)
        npOutside = beam.calcNumPartsOutside(Options::beamHaloBoundary*beam.get_rrms());
    // *gmsg << "npOutside 1 = " << npOutside << " beamHaloBoundary= " << Options::beamHaloBoundary << " rrms= " << beam.get_rrms() << endl;

    double  pathLength = 0.0;
    if (OpalData::getInstance()->isInOPALCyclMode())
        pathLength = beam.getLPath();
    else
        pathLength = sposRef;

    /// Write data to files. If this is the first write to the beam statistics file, write SDDS
    /// header information.
    ofstream os_statData;
    ofstream os_lBalData;
    double Q = beam.getCharge();

    if(Ippl::myNode() == 0) {
        if(firstWriteToStat_m) {
            os_statData.open(statFileName_m.c_str(), ios::out);
            os_statData.precision(15);
            os_statData.setf(ios::scientific, ios::floatfield);
            writeSDDSHeader(os_statData, losses);

            os_lBalData.open(lBalFileName_m.c_str(), ios::out);
            os_lBalData.precision(15);
            os_lBalData.setf(ios::scientific, ios::floatfield);
            os_lBalData << "# " << Ippl::getNodes() << endl;

            firstWriteToStat_m = false;
        } else {
            os_statData.open(statFileName_m.c_str(), ios::app);
            os_statData.precision(15);
            os_statData.setf(ios::scientific, ios::floatfield);

            os_lBalData.open(lBalFileName_m.c_str(), ios::app);
            os_lBalData.precision(15);
            os_lBalData.setf(ios::scientific, ios::floatfield);
        }

        os_statData << beam.getT() << setw(pwi) << "\t"                                       // 1
                    << pathLength << setw(pwi) << "\t"                                        // 2

                    << beam.getTotalNum() << setw(pwi) << "\t"                                // 3
                    << Q << setw(pwi) << "\t"                                                 // 4

                    << E << setw(pwi) << "\t"                                                 // 5

                    << beam.get_rrms()(0) << setw(pwi) << "\t"                                // 6
                    << beam.get_rrms()(1) << setw(pwi) << "\t"                                // 7
                    << beam.get_rrms()(2) << setw(pwi) << "\t"                                // 8

                    << beam.get_prms()(0) << setw(pwi) << "\t"                                // 9
                    << beam.get_prms()(1) << setw(pwi) << "\t"                                // 10
                    << beam.get_prms()(2) << setw(pwi) << "\t"                                // 11

                    << beam.get_norm_emit()(0) << setw(pwi) << "\t"                           // 12
                    << beam.get_norm_emit()(1) << setw(pwi) << "\t"                           // 13
                    << beam.get_norm_emit()(2) << setw(pwi) << "\t"                           // 14

                    << beam.get_rmean()(0)  << setw(pwi) << "\t"                              // 15
                    << beam.get_rmean()(1)  << setw(pwi) << "\t"                              // 16
                    << beam.get_rmean()(2)  << setw(pwi) << "\t"                              // 17

                    << beam.get_maxExtend()(0) << setw(pwi) << "\t"                           // 18
                    << beam.get_maxExtend()(1) << setw(pwi) << "\t"                           // 19
                    << beam.get_maxExtend()(2) << setw(pwi) << "\t"                           // 20

            // Write out Courant Snyder parameters.
                    << beam.get_rprrms()(0) << setw(pwi) << "\t"                              // 21
                    << beam.get_rprrms()(1) << setw(pwi) << "\t"                              // 22
                    << beam.get_rprrms()(2) << setw(pwi) << "\t"                              // 23

                    << 0.0 << setw(pwi) << "\t"                                                // 24
                    << 0.0 << setw(pwi) << "\t"                                                // 25

            // Write out dispersion.
                    << beam.get_Dx() << setw(pwi) << "\t"                                      // 26
                    << beam.get_DDx() << setw(pwi) << "\t"                                     // 27
                    << beam.get_Dy() << setw(pwi) << "\t"                                      // 28
                    << beam.get_DDy() << setw(pwi) << "\t"                                     // 29


            // Write head/reference particle/tail field information.
                    << FDext[0](0) << setw(pwi) << "\t"                                         // 30 B-head x
                    << FDext[0](1) << setw(pwi) << "\t"                                         // 31 B-head y
                    << FDext[0](2) << setw(pwi) << "\t"                                         // 32 B-head z

                    << FDext[1](0) << setw(pwi) << "\t"                                         // 33 E-head x
                    << FDext[1](1) << setw(pwi) << "\t"                                         // 34 E-head y
                    << FDext[1](2) << setw(pwi) << "\t"                                         // 35 E-head z

                    << FDext[2](0) << setw(pwi) << "\t"                                         // 36 B-ref x
                    << FDext[2](1) << setw(pwi) << "\t"                                         // 37 B-ref y
                    << FDext[2](2) << setw(pwi) << "\t"                                         // 38 B-ref z

                    << FDext[3](0) << setw(pwi) << "\t"                                         // 39 E-ref x
                    << FDext[3](1) << setw(pwi) << "\t"                                         // 40 E-ref y
                    << FDext[3](2) << setw(pwi) << "\t"                                         // 41 E-ref z

                    << FDext[4](0) << setw(pwi) << "\t"                                         // 42 B-tail x
                    << FDext[4](1) << setw(pwi) << "\t"                                         // 43 B-tail y
                    << FDext[4](2) << setw(pwi) << "\t"                                         // 44 B-tail z

                    << FDext[5](0) << setw(pwi) << "\t"                                         // 45 E-tail x
                    << FDext[5](1) << setw(pwi) << "\t"                                         // 46 E-tail y
                    << FDext[5](2) << setw(pwi) << "\t"                                         // 47 E-tail z

                    << beam.getdE() << setw(pwi) << "\t"                                        // 48 dE energy spread

                    << npOutside << setw(pwi) << "\t";                                          // 49 number of particles outside n*sigma

        if(Ippl::getNodes() == 1 && beam.getLocalNum() > 0) {
            os_statData << beam.R[0](0) << setw(pwi) << "\t";                                   // 50 R0_x
            os_statData << beam.R[0](1) << setw(pwi) << "\t";                                   // 51 R0_y
            os_statData << beam.R[0](2) << setw(pwi) << "\t";                                   // 52 R0_z
            os_statData << beam.P[0](0) << setw(pwi) << "\t";                                   // 53 P0_x
            os_statData << beam.P[0](1) << setw(pwi) << "\t";                                   // 54 P0_y
            os_statData << beam.P[0](2) << setw(pwi) << "\t";                                   // 55 P0_z
            os_statData << beam.R[0](0) << setw(pwi) << "\t";                                   // 50 R0_x
            os_statData << beam.R[1](1) << setw(pwi) << "\t";                                   // 51 R0_y
            os_statData << beam.R[1](2) << setw(pwi) << "\t";                                   // 52 R0_z
            os_statData << beam.P[1](0) << setw(pwi) << "\t";                                   // 53 P0_x
            os_statData << beam.P[1](1) << setw(pwi) << "\t";                                   // 54 P0_y
            os_statData << beam.P[1](2) << setw(pwi) << "\t";                                   // 55 P0_z
            os_statData << beam.R[2](0) << setw(pwi) << "\t";                                   // 50 R0_x
            os_statData << beam.R[2](1) << setw(pwi) << "\t";                                   // 51 R0_y
            os_statData << beam.R[2](2) << setw(pwi) << "\t";                                   // 52 R0_z
            os_statData << beam.P[2](0) << setw(pwi) << "\t";                                   // 53 P0_x
            os_statData << beam.P[2](1) << setw(pwi) << "\t";                                   // 54 P0_y
            os_statData << beam.P[2](2) << setw(pwi) << "\t";                                   // 55 P0_z
            os_statData << beam.R[3](0) << setw(pwi) << "\t";                                   // 50 R0_x
            os_statData << beam.R[3](1) << setw(pwi) << "\t";                                   // 51 R0_y
            os_statData << beam.R[3](2) << setw(pwi) << "\t";                                   // 52 R0_z
            os_statData << beam.P[3](0) << setw(pwi) << "\t";                                   // 53 P0_x
            os_statData << beam.P[3](1) << setw(pwi) << "\t";                                   // 54 P0_y
            os_statData << beam.P[3](2) << setw(pwi) << "\t";                                   // 55 P0_z
        }


        for(size_t i = 0; i < losses.size(); ++ i) {
            os_statData << losses[i].second << setw(pwi) << "\t";
        }
        os_statData   << endl;

        for(int p = 0; p < Ippl::getNodes(); p++)
            os_lBalData << beam.getLoadBalance(p)  << setw(pwi) << "\t";
        os_lBalData << endl;

        os_statData.close();
        os_lBalData.close();
    }

    /// %Stop timer.
    IpplTimings::stopTimer(StatMarkerTimer_m);
}

void DataSink::writeStatData(EnvelopeBunch &beam, Vector_t FDext[], double sposHead, double sposRef, double sposTail) {
    /// Function steps:
    /// Start timer.
    IpplTimings::startTimer(StatMarkerTimer_m);

    /// Set width of write fields in output files.
    unsigned int pwi = 10;

    /// Calculate beam statistics and gather load balance statistics.
    beam.calcBeamParameters();
    beam.gatherLoadBalanceStatistics();

    /// Write data to files. If this is the first write to the beam statistics file, write SDDS
    /// header information.
    ofstream os_statData;
    ofstream os_lBalData;
    double en = beam.get_meanEnergy() * 1e-6;

    if(Ippl::myNode() == 0) {
        if(firstWriteToStat_m) {
            os_statData.open(statFileName_m.c_str(), ios::out);
            os_statData.precision(15);
            os_statData.setf(ios::scientific, ios::floatfield);
            writeSDDSHeader(os_statData);

            os_lBalData.open(lBalFileName_m.c_str(), ios::out);
            os_lBalData.precision(15);
            os_lBalData.setf(ios::scientific, ios::floatfield);
            os_lBalData << "# " << Ippl::getNodes() << endl;

            firstWriteToStat_m = false;
        } else {
            os_statData.open(statFileName_m.c_str(), ios::app);
            os_statData.precision(15);
            os_statData.setf(ios::scientific, ios::floatfield);

            os_lBalData.open(lBalFileName_m.c_str(), ios::app);
            os_lBalData.precision(15);
            os_lBalData.setf(ios::scientific, ios::floatfield);
        }


        os_statData << beam.getT() << setw(pwi) << "\t"                                       // 1
                    << sposRef << setw(pwi) << "\t"                                           // 2

                    << beam.getTotalNum() << setw(pwi) << "\t"                                // 3
                    << beam.getTotalNum() * beam.getChargePerParticle() << setw(pwi) << "\t"  // 4
                    << en << setw(pwi) << "\t"                                                // 5

                    << beam.get_rrms()(0) << setw(pwi) << "\t"                                // 6
                    << beam.get_rrms()(1) << setw(pwi) << "\t"                                // 7
                    << beam.get_rrms()(2) << setw(pwi) << "\t"                                // 8

                    << beam.get_prms()(0) << setw(pwi) << "\t"                                // 9
                    << beam.get_prms()(1) << setw(pwi) << "\t"                                // 10
                    << beam.get_prms()(2) << setw(pwi) << "\t"                                // 11

                    << beam.get_norm_emit()(0) << setw(pwi) << "\t"                           // 12
                    << beam.get_norm_emit()(1) << setw(pwi) << "\t"                           // 13
                    << beam.get_norm_emit()(2) << setw(pwi) << "\t"                           // 14

                    << beam.get_rmean()(0)  << setw(pwi) << "\t"                              // 15
                    << beam.get_rmean()(1)  << setw(pwi) << "\t"                              // 16
                    << beam.get_rmean()(2)  << setw(pwi) << "\t"                              // 17

                    << beam.get_maxExtend()(0) << setw(pwi) << "\t"                           // 18
                    << beam.get_maxExtend()(1) << setw(pwi) << "\t"                           // 19
                    << beam.get_maxExtend()(2) << setw(pwi) << "\t"                           // 20

            // Write out Courant Snyder parameters.
                    << 0.0  << setw(pwi) << "\t"                                              // 21
                    << 0.0  << setw(pwi) << "\t"                                              // 22

                    << 0.0 << setw(pwi) << "\t"                                               // 23
                    << 0.0 << setw(pwi) << "\t"                                               // 24

            // Write out dispersion.
                    << beam.get_Dx() << setw(pwi) << "\t"                                     // 25
                    << beam.get_DDx() << setw(pwi) << "\t"                                    // 26
                    << beam.get_Dy() << setw(pwi) << "\t"                                     // 27
                    << beam.get_DDy() << setw(pwi) << "\t"                                    // 28


            // Write head/reference particle/tail field information.
                    << FDext[0](0) << setw(pwi) << "\t"                                       // 29 B-head x
                    << FDext[0](1) << setw(pwi) << "\t"                                       // 30 B-head y
                    << FDext[0](2) << setw(pwi) << "\t"                                       // 31 B-head z

                    << FDext[1](0) << setw(pwi) << "\t"                                       // 32 E-head x
                    << FDext[1](1) << setw(pwi) << "\t"                                       // 33 E-head y
                    << FDext[1](2) << setw(pwi) << "\t"                                       // 34 E-head z

                    << FDext[2](0) << setw(pwi) << "\t"                                       // 35 B-ref x
                    << FDext[2](1) << setw(pwi) << "\t"                                       // 36 B-ref y
                    << FDext[2](2) << setw(pwi) << "\t"                                       // 37 B-ref z

                    << FDext[3](0) << setw(pwi) << "\t"                                       // 38 E-ref x
                    << FDext[3](1) << setw(pwi) << "\t"                                       // 39 E-ref y
                    << FDext[3](2) << setw(pwi) << "\t"                                       // 40 E-ref z

                    << FDext[4](0) << setw(pwi) << "\t"                                       // 41 B-tail x
                    << FDext[4](1) << setw(pwi) << "\t"                                       // 42 B-tail y
                    << FDext[4](2) << setw(pwi) << "\t"                                       // 43 B-tail z

                    << FDext[5](0) << setw(pwi) << "\t"                                       // 44 E-tail x
                    << FDext[5](1) << setw(pwi) << "\t"                                       // 45 E-tail y
                    << FDext[5](2) << setw(pwi) << "\t"                                       // 46 E-tail z

                    << beam.get_dEdt() << setw(pwi) << "\t"                                   // 47 dE energy spread

                    << endl;


        for(int p = 0; p < Ippl::getNodes(); p++)
            os_lBalData << beam.getLoadBalance(p)  << setw(pwi) << "\t";
        os_lBalData << endl;

        os_statData.close();
        os_lBalData.close();
    }

    /// %Stop timer.
    IpplTimings::stopTimer(StatMarkerTimer_m);
}

void DataSink::writeSDDSHeader(ofstream &outputFile) {
    writeSDDSHeader(outputFile,
                    std::vector<std::pair<std::string, unsigned int> >());
}
void DataSink::writeSDDSHeader(ofstream &outputFile,
                               const std::vector<std::pair<std::string, unsigned int> > &losses) {
    OPALTimer::Timer simtimer;

    string dateStr(simtimer.date());
    string timeStr(simtimer.time());

    outputFile << "SDDS1" << endl;
    outputFile << "&description text=\"Statistics data " << OpalData::getInstance()->getInputFn() << " " << dateStr << " " << timeStr << "\" " << endl;
    outputFile << ", contents=\"stat parameters\" &end" << endl;

    outputFile << "&parameter name=processors, type=long, ";
    outputFile << "description=\"Number of Cores used\" &end" << endl;

    outputFile << "&parameter name=revision, type=string, "
               << "description=\"git revision of opal\" &end\n";

    outputFile << "&column name=t, type=double, units=s, ";
    outputFile << "description=\"1 Time\" &end" << endl;
    outputFile << "&column name=s, type=double, units=m, ";
    outputFile << "description=\"2 Average Longitudinal Position\" &end" << endl;

    outputFile << "&column name=numParticles, type=long, units=1, ";
    outputFile << "description=\"3 Number of Macro Particles\" &end" << endl;
    outputFile << "&column name=charge, type=double, units=1, ";
    outputFile << "description=\"4 Bunch Charge\" &end" << endl;

    outputFile << "&column name=energy, type=double, units=MeV, ";
    outputFile << "description=\"5 Mean Energy\" &end" << endl;

    outputFile << "&column name=rms_x, type=double, units=m , ";
    outputFile << "description=\"6 RMS Beamsize in x  \" &end" << endl;
    outputFile << "&column name=rms_y, type=double, units=m , ";
    outputFile << "description=\"7 RMS Beamsize in y  \" &end" << endl;
    outputFile << "&column name=rms_s, type=double, units=m , ";
    outputFile << "description=\"8 RMS Beamsize in s  \" &end" << endl;

    outputFile << "&column name=rms_px, type=double, units=1 , ";
    outputFile << "description=\"9 RMS Momenta in x  \" &end" << endl;
    outputFile << "&column name=rms_py, type=double, units=1 , ";
    outputFile << "description=\"10 RMS Momenta in y  \" &end" << endl;
    outputFile << "&column name=rms_ps, type=double, units=1 , ";
    outputFile << "description=\"11 RMS Momenta in s  \" &end" << endl;

    outputFile << "&column name=emit_x, type=double, units=m , ";
    outputFile << "description=\"12 Normalized Emittance x  \" &end" << endl;
    outputFile << "&column name=emit_y, type=double, units=m , ";
    outputFile << "description=\"13 Normalized Emittance y  \" &end" << endl;
    outputFile << "&column name=emit_s, type=double, units=m , ";
    outputFile << "description=\"14 Normalized Emittance s  \" &end" << endl;

    outputFile << "&column name=mean_x, type=double, units=m , ";
    outputFile << "description=\"15 Mean Beam Position in x  \" &end" << endl;
    outputFile << "&column name=mean_y, type=double, units=m , ";
    outputFile << "description=\"16 Mean Beam Position in y  \" &end" << endl;
    outputFile << "&column name=mean_s, type=double, units=m , ";
    outputFile << "description=\"17 Mean Beam Position in s  \" &end" << endl;

    outputFile << "&column name=max_x, type=double, units=m , ";
    outputFile << "description=\"18 Max Beamsize in x  \" &end" << endl;
    outputFile << "&column name=max_y, type=double, units=m , ";
    outputFile << "description=\"19 Max Beamsize in y  \" &end" << endl;
    outputFile << "&column name=max_s, type=double, units=m , ";
    outputFile << "description=\"20 Max Beamsize in s  \" &end" << endl;

    outputFile << "&column name=xpx, type=double, units=1 , ";
    outputFile << "description=\"21 Correlation xpx  \" &end" << endl;
    outputFile << "&column name=ypy, type=double, units=1 , ";
    outputFile << "description=\"22 Correlation ypyy  \" &end" << endl;

    outputFile << "&column name=zpz, type=double, units=1 , ";
    outputFile << "description=\"23 Correlation zpz  \" &end" << endl;

    outputFile << "&column name=notused1, type=double, units=1 , ";
    outputFile << "description=\"24 notused1 in y  \" &end" << endl;

    outputFile << "&column name=notused2, type=double, units=1 , ";
    outputFile << "description=\"25 notused2 in y  \" &end" << endl;

    if (Ippl::getNodes()>1) {
      outputFile << "&column name=Dx, type=double, units=m , ";
      outputFile << "description=\"26 Dispersion in x  \" &end" << endl;
      outputFile << "&column name=DDx, type=double, units=1 , ";
      outputFile << "description=\"27 Derivative of dispersion in x  \" &end" << endl;

      outputFile << "&column name=Dy, type=double, units=m , ";
      outputFile << "description=\"28 Dispersion in y  \" &end" << endl;
      outputFile << "&column name=DDy, type=double, units=1 , ";
      outputFile << "description=\"29 Derivative of dispersion in y  \" &end" << endl;
    }
    else {
      outputFile << "&column name=Dx1, type=double, units=m , ";
      outputFile << "description=\"26 Dispersion in x particle 1  \" &end" << endl;
      outputFile << "&column name=Dx2, type=double, units=m , ";
      outputFile << "description=\"27 Dispersion in x particle 2  \" &end" << endl;

      outputFile << "&column name=Dx3, type=double, units=m , ";
      outputFile << "description=\"28 Dispersion in x particle 3  \" &end" << endl;
      outputFile << "&column name=Dx4, type=double, units=m , ";
      outputFile << "description=\"29 Dispersion in x particle 4  \" &end" << endl;
    }
    outputFile << "&column name=Bx_head, type=double, units=T , ";
    outputFile << "description=\"30 Bx-Field component of head particle  \" &end" << endl;
    outputFile << "&column name=By_head, type=double, units=T , ";
    outputFile << "description=\"31 By-Field component of head particle  \" &end" << endl;
    outputFile << "&column name=Bz_head, type=double, units=T , ";
    outputFile << "description=\"32 Bz-Field component of head particle  \" &end" << endl;

    outputFile << "&column name=Ex_head, type=double, units=MV/m , ";
    outputFile << "description=\"33 Ex-Field component of head particle  \" &end" << endl;
    outputFile << "&column name=Ey_head, type=double, units=MV/m , ";
    outputFile << "description=\"34 Ey-Field component of head particle  \" &end" << endl;
    outputFile << "&column name=Ez_head, type=double, units=MV/m , ";
    outputFile << "description=\"35 Ez-Field component of head particle  \" &end" << endl;

    outputFile << "&column name=Bx_ref, type=double, units=T , ";
    outputFile << "description=\"36 Bx-Field component of ref particle  \" &end" << endl;
    outputFile << "&column name=By_ref, type=double, units=T , ";
    outputFile << "description=\"37 By-Field component of ref particle  \" &end" << endl;
    outputFile << "&column name=Bz_ref, type=double, units=T , ";
    outputFile << "description=\"38 Bz-Field component of ref particle  \" &end" << endl;

    outputFile << "&column name=Ex_ref, type=double, units=MV/m , ";
    outputFile << "description=\"39 Ex-Field component of ref particle  \" &end" << endl;
    outputFile << "&column name=Ey_ref, type=double, units=MV/m , ";
    outputFile << "description=\"40 Ey-Field component of ref particle  \" &end" << endl;
    outputFile << "&column name=Ez_ref, type=double, units=MV/m , ";
    outputFile << "description=\"41 Ez-Field component of ref particle  \" &end" << endl;

    outputFile << "&column name=Bx_tail, type=double, units=T , ";
    outputFile << "description=\"42 Bx-Field component of tail particle  \" &end" << endl;
    outputFile << "&column name=By_tail, type=double, units=T , ";
    outputFile << "description=\"43 By-Field component of tail particle  \" &end" << endl;
    outputFile << "&column name=Bz_tail, type=double, units=T , ";
    outputFile << "description=\"44 Bz-Field component of tail particle  \" &end" << endl;

    outputFile << "&column name=Ex_tail, type=double, units=MV/m , ";
    outputFile << "description=\"45 Ex-Field component of tail particle  \" &end" << endl;
    outputFile << "&column name=Ey_tail, type=double, units=MV/m , ";
    outputFile << "description=\"46 Ey-Field component of tail particle  \" &end" << endl;
    outputFile << "&column name=Ez_tail, type=double, units=MV/m , ";
    outputFile << "description=\"47 Ez-Field component of tail particle  \" &end" << endl;

    outputFile << "&column name=dE, type=double, units=MeV , ";
    outputFile << "description=\"48 energy spread of the beam  \" &end" << endl;

    outputFile << "&column name=partsOutside, type=double, units=1 , ";
    outputFile << "description=\"49 outside n*sigma of the beam  \" &end" << endl;

    unsigned int columnStart = 49;
    if(Ippl::getNodes() == 1) {
        outputFile << "&column name=R0_x, type=double, units=m , ";
        outputFile << "description=\"50 R0 Particle Position in x  \" &end" << endl;
        outputFile << "&column name=R0_y, type=double, units=m , ";
        outputFile << "description=\"51 R0 Particle Position in y  \" &end" << endl;
        outputFile << "&column name=R0_s, type=double, units=m , ";
        outputFile << "description=\"52 R0 Particle Position in s  \" &end" << endl;
        outputFile << "&column name=P0_x, type=double, units=1 , ";
        outputFile << "description=\"53 P0 Particle Position in px  \" &end" << endl;
        outputFile << "&column name=P0_y, type=double, units=1 , ";
        outputFile << "description=\"54 P0 Particle Position in py  \" &end" << endl;
        outputFile << "&column name=P0_s, type=double, units=1 , ";
        outputFile << "description=\"55 P0 Particle Position in ps  \" &end" << endl;

        outputFile << "&column name=R1_x, type=double, units=m , ";
        outputFile << "description=\"56 R1 Particle Position in x  \" &end" << endl;
        outputFile << "&column name=R1_y, type=double, units=m , ";
        outputFile << "description=\"57 R1 Particle Position in y  \" &end" << endl;
        outputFile << "&column name=R1_s, type=double, units=m , ";
        outputFile << "description=\"58 R1 Particle Position in s  \" &end" << endl;
        outputFile << "&column name=P1_x, type=double, units=1 , ";
        outputFile << "description=\"59 P1 Particle Position in px  \" &end" << endl;
        outputFile << "&column name=P1_y, type=double, units=1 , ";
        outputFile << "description=\"60 P1 Particle Position in py  \" &end" << endl;
        outputFile << "&column name=P1_s, type=double, units=1 , ";
        outputFile << "description=\"61 P1 Particle Position in ps  \" &end" << endl;

        outputFile << "&column name=R2_x, type=double, units=m , ";
        outputFile << "description=\"62 R2 Particle Position in x  \" &end" << endl;
        outputFile << "&column name=R2_y, type=double, units=m , ";
        outputFile << "description=\"63 R2 Particle Position in y  \" &end" << endl;
        outputFile << "&column name=R2_s, type=double, units=m , ";
        outputFile << "description=\"64 R2 Particle Position in s  \" &end" << endl;
        outputFile << "&column name=P2_x, type=double, units=1 , ";
        outputFile << "description=\"65 P2 Particle Position in px  \" &end" << endl;
        outputFile << "&column name=P2_y, type=double, units=1 , ";
        outputFile << "description=\"66 P2 Particle Position in py  \" &end" << endl;
        outputFile << "&column name=P2_s, type=double, units=1 , ";
        outputFile << "description=\"67 P2 Particle Position in ps  \" &end" << endl;

        outputFile << "&column name=R3_x, type=double, units=m , ";
        outputFile << "description=\"68 R3 Particle Position in x  \" &end" << endl;
        outputFile << "&column name=R3_y, type=double, units=m , ";
        outputFile << "description=\"69 R3 Particle Position in y  \" &end" << endl;
        outputFile << "&column name=R3_s, type=double, units=m , ";
        outputFile << "description=\"70 R3 Particle Position in s  \" &end" << endl;
        outputFile << "&column name=P3_x, type=double, units=1 , ";
        outputFile << "description=\"71 P3 Particle Position in px  \" &end" << endl;
        outputFile << "&column name=P3_y, type=double, units=1 , ";
        outputFile << "description=\"72 P3 Particle Position in py  \" &end" << endl;
        outputFile << "&column name=P3_s, type=double, units=1 , ";
        outputFile << "description=\"73 P3 Particle Position in ps  \" &end" << endl;

        columnStart = 73;
    }

    for (size_t i = 0; i < losses.size(); ++ i) {
        outputFile << "&column name=" << losses[i].first << ", type=long, units=1, ";
        outputFile << "description=\"" << columnStart++ << " " << losses[i].second << "\" &end" << endl;
    }
    outputFile << "&data mode=ascii, no_row_counts=1 &end" << endl;

    outputFile << Ippl::getNodes() << endl;
    outputFile << PACKAGE_NAME << " " << PACKAGE_VERSION << " git rev. " << GIT_VERSION << endl;
}


void DataSink::writePartlossZASCII(PartBunch &beam, BoundaryGeometry &bg, string fn) {

    size_t temp = lossWrCounter_m ;

    string ffn = fn + convertToString(temp) + string("Z.dat");
    std::unique_ptr<Inform> ofp(new Inform(NULL, ffn.c_str(), Inform::OVERWRITE, 0));
    Inform &fid = *ofp;
    setInform(fid);
    fid.precision(6);

    string ftrn =  fn + string("triangle") + convertToString(temp) + string(".dat");
    std::unique_ptr<Inform> oftr(new Inform(NULL, ftrn.c_str(), Inform::OVERWRITE, 0));
    Inform &fidtr = *oftr;
    setInform(fidtr);
    fidtr.precision(6);

    Vector_t Geo_nr = bg.getnr();
    Vector_t Geo_hr = bg.gethr();
    Vector_t Geo_mincoords = bg.getmincoords();
    double t = beam.getT();
    double t_step = t * 1.0e9;
    double* prPartLossZ = new double[bg.getnr() (2)];
    double* sePartLossZ = new double[bg.getnr() (2)];
    double* fePartLossZ = new double[bg.getnr() (2)];
    fidtr << "# Time/ns" << std::setw(18) << "Triangle_ID" << std::setw(18)
          << "Xcoordinates (m)" << std::setw(18)
          << "Ycoordinates (m)" << std::setw(18)
          << "Zcoordinates (m)" << std::setw(18)
          << "Primary part. charge (C)" << std::setw(40)
          << "Field emit. part. charge (C)" << std::setw(40)
          << "Secondary emit. part. charge (C)" << std::setw(40) << endl ;
    for(int i = 0; i < Geo_nr(2) ; i++) {
        prPartLossZ[i] = 0;
        sePartLossZ[i] = 0;
        fePartLossZ[i] = 0;
        for(int j = 0; j < bg.getNumBFaces(); j++) {
            if(((Geo_mincoords[2] + Geo_hr(2)*i) < bg.TriBarycenters_m[j](2))
               && (bg.TriBarycenters_m[j](2) < (Geo_hr(2)*i + Geo_hr(2) + Geo_mincoords[2]))) {
                prPartLossZ[i] += bg.TriPrPartloss_m[j];
                sePartLossZ[i] += bg.TriSePartloss_m[j];
                fePartLossZ[i] += bg.TriFEPartloss_m[j];
            }

        }
    }
    for(int j = 0; j < bg.getNumBFaces(); j++) {
        fidtr << t_step << std::setw(18) << j << std::setw(18)// fixme: maybe gether particle loss data, i.e., do a reduce() for each triangle in each node befor write to file.
              << bg.TriBarycenters_m[j](0) << std::setw(18)
              << bg.TriBarycenters_m[j](1) << std::setw(18)
              << bg.TriBarycenters_m[j](2) <<  std::setw(40)
              << -bg.TriPrPartloss_m[j] << std::setw(40)
              << -bg.TriFEPartloss_m[j] <<  std::setw(40)
              << -bg.TriSePartloss_m[j] << endl;
    }
    fid << "# Delta_Z/m" << std::setw(18)
        << "Zcoordinates (m)" << std::setw(18)
        << "Primary part. charge (C)" << std::setw(40)
        << "Field emit. part. charge (C)" << std::setw(40)
        << "Secondary emit. part. charge (C)" << std::setw(40) << "t" << endl ;


    for(int i = 0; i < Geo_nr(2) ; i++) {
        double primaryPLoss = -prPartLossZ[i];
        double secondaryPLoss = -sePartLossZ[i];
        double fieldemissionPLoss = -fePartLossZ[i];
        reduce(primaryPLoss, primaryPLoss, OpAddAssign());
        reduce(secondaryPLoss, secondaryPLoss, OpAddAssign());
        reduce(fieldemissionPLoss, fieldemissionPLoss, OpAddAssign());
        fid << Geo_hr(2) << std::setw(18)
            << Geo_mincoords[2] + Geo_hr(2)*i << std::setw(18)
            << primaryPLoss << std::setw(40)
            << fieldemissionPLoss << std::setw(40)
            << secondaryPLoss << std::setw(40) << t << endl;
    }
    lossWrCounter_m++;
}

void DataSink::writeSurfaceInteraction(PartBunch &beam, long long &step, BoundaryGeometry &bg, string fn) {

    if (!doHDF5_m) return;

    h5_int64_t rc;
    /// Start timer.
    IpplTimings::startTimer(H5PartTimer_m);
    if(firstWriteH5Surface_m) {
        firstWriteH5Surface_m = false;

#if defined (USE_H5HUT2)
	h5_prop_t props = H5CreateFileProp ();
	MPI_Comm comm = Ippl::getComm();
	H5SetPropFileMPIOCollective (props, &comm);
	H5fileS_m = H5OpenFile (surfaceLossFileName_m.c_str(), H5_O_WRONLY, props);
        if(H5fileS_m == (h5_file_t)H5_ERR) {
            throw OpalException("DataSink::writeSurfaceInteraction",
                                "failed to open h5 file '" + surfaceLossFileName_m + "' for surface loss");
        }
#else
        H5fileS_m = H5OpenFile(surfaceLossFileName_m.c_str(), H5_FLUSH_STEP | H5_O_WRONLY, Ippl::getComm());
        if(H5fileS_m == (void*)H5_ERR) {
            throw OpalException("DataSink::writeSurfaceInteraction",
                                "failed to open h5 file '" + surfaceLossFileName_m + "' for surface loss");
        }
#endif

    }
    int nTot = bg.getNumBFaces();

    int N_mean = static_cast<int>(floor(nTot / Ippl::getNodes()));
    int N_extra = static_cast<int>(nTot - N_mean * Ippl::getNodes());
    int pc = 0;
    int count = 0;
    if(Ippl::myNode() == 0) {
        N_mean += N_extra;
    }
    std::unique_ptr<char[]> varray(new char[(N_mean)*sizeof(double)]);
    double *farray = reinterpret_cast<double *>(varray.get());
    h5_int64_t *larray = reinterpret_cast<h5_int64_t *>(varray.get());


    rc = H5SetStep(H5fileS_m, step);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    rc = H5PartSetNumParticles(H5fileS_m, N_mean);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    double    qi = beam.getChargePerParticle();
    rc = H5WriteStepAttribFloat64(H5fileS_m, "qi", &qi, 1);

    std::unique_ptr<double[]> tmploss(new double[nTot]);
    for(int i = 0; i < nTot; i++)
        tmploss[i] = 0.0;
    //memset( tmploss, 0.0, nTot * sizeof(double));
    reduce(bg.TriPrPartloss_m, bg.TriPrPartloss_m + nTot, tmploss.get(), OpAddAssign()); // may be removed if we have parallelized the geometry .

    for(int i = 0; i < nTot; i++) {
        if(pc == Ippl::myNode()) {
            if(count < N_mean) {
                if(pc != 0) {
                    //farray[count] =  bg.TriPrPartloss_m[Ippl::myNode()*N_mean+count];
                    if(((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::Absorption)) == (BGphysics::Absorption)) && (((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) && ((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission)))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[pc * N_mean + count + N_extra];
                    }
                    count ++;
                } else {
                    if(((bg.TriBGphysicstag_m[count] & (BGphysics::Absorption)) == (BGphysics::Absorption)) && (((bg.TriBGphysicstag_m[count] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) && ((bg.TriBGphysicstag_m[count] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission)))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[count];
                    }
                    count ++;

                }
            }
        }
        pc++;
        if(pc == Ippl::getNodes())
            pc = 0;
    }
    rc = H5PartWriteDataFloat64(H5fileS_m, "PrimaryLoss", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(int i = 0; i < nTot; i++)
        tmploss[i] = 0.0;
    reduce(bg.TriSePartloss_m, bg.TriSePartloss_m + nTot, tmploss.get(), OpAddAssign()); // may be removed if we parallelize the geometry as well.
    count = 0;
    pc = 0;
    for(int i = 0; i < nTot; i++) {
        if(pc == Ippl::myNode()) {
            if(count < N_mean) {

                if(pc != 0) {

                    if(((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::Absorption)) == (BGphysics::Absorption)) && (((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) && ((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission)))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[pc * N_mean + count + N_extra];
                    }
                    count ++;
                } else {
                    if(((bg.TriBGphysicstag_m[count] & (BGphysics::Absorption)) == (BGphysics::Absorption)) && (((bg.TriBGphysicstag_m[count] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) && ((bg.TriBGphysicstag_m[count] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission)))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[count];
                    }
                    count ++;

                }
            }
        }
        pc++;
        if(pc == Ippl::getNodes())
            pc = 0;
    }
    rc = H5PartWriteDataFloat64(H5fileS_m, "SecondaryLoss", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(int i = 0; i < nTot; i++)
        tmploss[i] = 0.0;

    reduce(bg.TriFEPartloss_m, bg.TriFEPartloss_m + nTot, tmploss.get(), OpAddAssign()); // may be removed if we parallelize the geometry as well.
    count = 0;
    pc = 0;
    for(int i = 0; i < nTot; i++) {
        if(pc == Ippl::myNode()) {
            if(count < N_mean) {

                if(pc != 0) {

                    if(((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::Absorption)) == (BGphysics::Absorption)) && (((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) && ((bg.TriBGphysicstag_m[pc * N_mean + count + N_extra] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission)))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[pc * N_mean + count + N_extra];
                    }
                    count ++;
                } else {
                    if(((bg.TriBGphysicstag_m[count] & (BGphysics::Absorption)) == (BGphysics::Absorption)) && (((bg.TriBGphysicstag_m[count] & (BGphysics::FNEmission)) != (BGphysics::FNEmission)) && ((bg.TriBGphysicstag_m[count] & (BGphysics::SecondaryEmission)) != (BGphysics::SecondaryEmission)))) {
                        farray[count] = 0.0;
                    } else {
                        farray[count] =  tmploss[count];
                    }
                    count ++;

                }

            }
        }
        pc++;
        if(pc == Ippl::getNodes())
            pc = 0;
    }
    rc = H5PartWriteDataFloat64(H5fileS_m, "FNEmissionLoss", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    count = 0;
    pc = 0;
    for(int i = 0; i < nTot; i++) {
        if(pc == Ippl::myNode()) {
            if(count < N_mean) {
                if(pc != 0)
                    larray[count] =  Ippl::myNode() * N_mean + count + N_extra; // node 0 will be 0*N_mean+count+N_extra also correct.
                else
                    larray[count] = count;
                count ++;
            }
        }
        pc++;
        if(pc == Ippl::getNodes())
            pc = 0;
    }
    rc = H5PartWriteDataInt64(H5fileS_m, "TriangleID", larray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    /// %Stop timer.
    IpplTimings::stopTimer(H5PartTimer_m);

}

void DataSink::writeImpactStatistics(PartBunch &beam, long long &step, size_t &impact, double &sey_num,
                                     size_t numberOfFieldEmittedParticles, bool nEmissionMode, string fn) {

    double charge = 0.0;
    size_t Npart = 0;
    double Npart_d = 0.0;
    if(!nEmissionMode) {
        charge = -1.0 * beam.getCharge();
        //reduce(charge, charge, OpAddAssign());
        Npart_d = -1.0 * charge / beam.getChargePerParticle();
    } else {
        Npart = beam.getTotalNum();
    }
    if(Ippl::myNode() == 0) {
        string ffn = fn + string(".dat");

        std::unique_ptr<Inform> ofp(new Inform(NULL, ffn.c_str(), Inform::APPEND, 0));
        Inform &fid = *ofp;
        setInform(fid);

        fid.precision(6);
        fid << setiosflags(ios::scientific);
        double t = beam.getT() * 1.0e9;
        if(!nEmissionMode) {

            if(step == 0) {
                fid << "#Time/ns"  << std::setw(18) << "#Geometry impacts" << std::setw(18) << "tot_sey" << std::setw(18)
                    << "TotalCharge" << std::setw(18) << "PartNum" << " numberOfFieldEmittedParticles " << endl ;
            }
            fid << t << std::setw(18) << impact << std::setw(18) << sey_num << std::setw(18) << charge
                << std::setw(18) << Npart_d << std::setw(18) << numberOfFieldEmittedParticles << endl ;
        } else {

            if(step == 0) {
                fid << "#Time/ns"  << std::setw(18) << "#Geometry impacts" << std::setw(18) << "tot_sey" << std::setw(18)
                    << "ParticleNumber" << " numberOfFieldEmittedParticles " << endl;
            }
            fid << t << std::setw(18) << impact << std::setw(18) << sey_num
                << std::setw(18) << double(Npart) << std::setw(18) << numberOfFieldEmittedParticles << endl ;
        }
    }
}

void DataSink::writeGeomToVtk(BoundaryGeometry &bg, string fn) {
    if(Ippl::myNode() == 0) {
        bg.writeGeomToVtk (fn);
    }
}

/** \brief Find out which if we write HDF5 or not
 *
 *
 *
 */
bool DataSink::doHDF5() {
    return doHDF5_m;
}


/** \brief
 *  delete the last 'numberOfLines' lines of the file 'fileName'
 */
void DataSink::rewindLines(const std::string &fileName, size_t numberOfLines) const {
    if (numberOfLines == 0) return;

    std::string line;
    std::queue<std::string> allLines;
    std::fstream fs;

    fs.open (fileName.c_str(), std::fstream::in);

    if (!fs.is_open()) return;

    while (getline(fs, line)) {
        allLines.push(line);
    }
    fs.close();


    fs.open (fileName.c_str(), std::fstream::out);

    if (!fs.is_open()) return;

    while (allLines.size() > numberOfLines) {
        fs << allLines.front() << "\n";
        allLines.pop();
    }
    fs.close();
}

/** \brief
 *  rewind the SDDS file such that the spos of the last step is less or equal to maxSPos
 */
unsigned int DataSink::rewindLinesSDDS(const std::string &fileName, double maxSPos) const {
    std::string line;
    std::queue<std::string> allLines;
    std::fstream fs;
    unsigned int numParameters = 0;
    unsigned int numColumns = 0;
    unsigned int sposColumnNr = 0;
    double spos, lastSPos = -1.0;

    boost::regex parameters("^&parameter name=");
    boost::regex column("^&column name=([a-zA-Z0-9\\$_]+),");
    boost::regex data("^&data mode=ascii");
    boost::smatch match;

    std::istringstream linestream;
    fs.open (fileName.c_str(), std::fstream::in);

    if (!fs.is_open()) return 0;

    while (getline(fs, line)) {
        allLines.push(line);
    }
    fs.close();


    fs.open (fileName.c_str(), std::fstream::out);

    if (!fs.is_open()) return 0;

    do {
        line = allLines.front();
        allLines.pop();
        fs << line << "\n";
        if (boost::regex_search(line, match, parameters)) {
            ++numParameters;
        } else if (boost::regex_search(line, match, column)) {
            ++numColumns;
            if (match[1] == "s") {
                sposColumnNr = numColumns;
            }
        }
    } while (!boost::regex_search(line, match, data));

    for (unsigned int i = 0; i < numParameters; ++ i) {
        fs << allLines.front() << "\n";
        allLines.pop();
    }

    while (allLines.size() > 0) {
        line = allLines.front();
        linestream.str(line);

        for (unsigned int i = 0; i < sposColumnNr; ++ i) {
            linestream >> spos;
        }

        if ((spos - maxSPos) > 1e-20 * Physics::c || (spos - lastSPos) < 1e-20 * Physics::c) break;

        allLines.pop();

        fs << line << "\n";

        lastSPos = spos;
    }

    fs.close();

    return allLines.size();
}

/***************************************************************************
 * $RCSfile: DataSink.cpp,v $   $Author: adelmann $
 * $Revision: 1.3 $   $Date: 2004/06/02 19:38:54 $
 ***************************************************************************/
