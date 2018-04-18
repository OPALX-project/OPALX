#ifndef OPAL_STATISTICALERRORS_H
#define OPAL_STATISTICALERRORS_H

//
//  Copyright & License: See Copyright.readme in src directory
//

/*!
  Class documentation
*/

#define SE_VISITELEMENT(elem) virtual void visit##elem(const elem &) { }

#include "Algorithms/Tracker.h"

#ifdef WITH_UNIT_TESTS
#include <gtest/gtest_prod.h>
#endif

#include <iostream>
#include <list>

class BMultipoleField;
class PartBunch;
class AlignWrapper;
class BeamBeam;
class Collimator;
class Corrector;
class CyclotronValley;
class Degrader;
class Diagnostic;
class Drift;
class ElementBase;
class FlexibleCollimator;
class Lambertson;
class Marker;
class Monitor;
class Multipole;
class ParallelPlate;
class Probe;
class RBend;
class RFCavity;
class RFQuadrupole;
class SBend;
class Separator;
class Septum;
class Solenoid;
class TravelingWave;
struct OpalInputInterpreter;


class StatisticalErrors: public Tracker {

public:
    /// Constructor.
    explicit StatisticalErrors(const Beamline &beamline,
                               const PartData &reference,
                               bool revBeam,
                               bool revTrack,
                               const std::string &method,
                               unsigned int nodesPerInstance,
                               unsigned int numInstances,
                               const std::vector<std::string> &objectives);

    virtual ~StatisticalErrors();

    SE_VISITELEMENT(AlignWrapper)
    SE_VISITELEMENT(Beamline)
    SE_VISITELEMENT(BeamBeam)
    SE_VISITELEMENT(Collimator)
    SE_VISITELEMENT(Corrector)
    SE_VISITELEMENT(CyclotronValley)
    SE_VISITELEMENT(Degrader)
    SE_VISITELEMENT(Diagnostic)
    SE_VISITELEMENT(Drift)
    SE_VISITELEMENT(FlexibleCollimator)
    SE_VISITELEMENT(Lambertson)
    SE_VISITELEMENT(Marker)
    SE_VISITELEMENT(Monitor)
    SE_VISITELEMENT(Multipole)
    SE_VISITELEMENT(ParallelPlate)
    SE_VISITELEMENT(Probe)
    SE_VISITELEMENT(RBend)
    SE_VISITELEMENT(RFCavity)
    SE_VISITELEMENT(RFQuadrupole)
    SE_VISITELEMENT(SBend)
    SE_VISITELEMENT(Separator)
    SE_VISITELEMENT(Septum)
    SE_VISITELEMENT(Solenoid)
    SE_VISITELEMENT(TravelingWave)

    virtual void execute();

private:

    StatisticalErrors();
    StatisticalErrors(const StatisticalErrors &);

    void operator=(const StatisticalErrors &);

    void runSimulation(const std::string &inputFileName, MPI_Comm comm);
    void stashEnvironment();
    void popEnvironment();

    static std::string getNextDirectoryName(bool reference);
    void createDirectory(const std::string &name);
    void linkFiles(const std::string &path);
    void removeLinks(const std::string &path);

    void calcIdealDivision(unsigned int numNodes);
    void formGroups();
    void shutReplicasDown(int tag,
                          const std::vector<std::string> &currentJobs);
    std::pair<int, int> identifyFreeReplica(int tag);
    std::string assignReplicaNewJob(const OpalInputInterpreter &interpreter,
                                    const std::pair<int, int> &source,
                                    bool referenceRun,
                                    std::string nextDir,
                                    int tag);

    struct PRMessage {
        int status;
        std::string nextDirectory;
    };

    PRMessage getCrunchJobFromPrimary(int tag);
    void collectOutputData(const std::string &directory);
    void processOutputData(const std::string &directory);
    // std::string readSDDSFile(const std::string &directory);
    void writeSDDSFile();
    void writeSDDSHeader(std::ofstream &out,
                         std::vector<std::string> &order);
    void writeSDDSData(std::ofstream &out,
                       const std::map<std::string, std::vector<double> > &data,
                       const std::vector<std::string> &order);
    std::vector<std::pair<double, double> > computeStatistics(const std::string &objective,
                                                              const std::vector<double> &referencePositions);
    std::vector<double> getNextDataset(std::ifstream &dataStream,
                                       std::ifstream &positionStream,
                                       int dataType,
                                       const std::vector<double> &referencePositions);

    std::vector<double > interpolateSDDSData(const std::vector<double> &data,
                                             const std::vector<double> &oldpositions,
                                             const std::vector<double> &newpositions);
    std::vector<double> readFloatData(std::ifstream &in);
    std::vector<double> readDoubleData(std::ifstream &in);
    std::vector<double> readShortData(std::ifstream &in);
    std::vector<double> readLongData(std::ifstream &in);

    // SDDS::file parseSDDSFile(std::string & contents);
    void removeOldDataFiles();
    void removeBinaryDataFiles();
    std::map<std::string, std::vector<double> > distributeProcessingJobs();
    void getProcessingJobFromPrimary();

    std::vector<double> getReferenceSamplingPositions();
    void sendReferenceSamplingPositions(std::vector<double> positionSampling);
    void sendObjectiveForProcessing(int tag,
                                    int source,
                                    unsigned int i,
                                    char* messageBuffer);
    std::pair<int, int> getNextObjectiveToProcess(int tag, char *messageBuffer);
    void sendProcessedDataToPrimary(int tag,
                                    int objectiveID,
                                    unsigned int length,
                                    char *messageBuffer);
    void copyDataToBuffer(char *messageBuffer,
                          const std::vector<std::pair<double, double> > &statistics);

    struct DataSource {
        int source;
        int ID;
        std::vector<double> data;
    };

    DataSource getProcessedDataFromReplica(int tag,
                                           unsigned int i,
                                           unsigned int length,
                                           char* messageBuffer);

    std::string method_m;
    unsigned int nodesPerInstance_m;
    unsigned int numInstances_m;

    MPI_Comm subComm_m;
    MPI_Group subGroup_m;
    std::vector<int> groupIDs_m;
    std::vector<int> groupRanks_m;
    std::vector<int> localMasters_m;
    std::vector<std::string> objectives_m;

    std::string directoryBaseName_m;
};

#endif // OPAL_STATISTICALERRORS_H