//
//  Copyright & License: See Copyright.readme in src directory
//

#include "Algorithms/StatisticalErrors.h"

#include "Algorithms/StatisticalErrorsUtilities.h"
#include "AbstractObjects/OpalData.h"
#include "OpalConfigure/Configure.h"
#include "OpalParser/OpalParser.h"
#include "Parser/FileStream.h"
#include "Structure/OpalInputInterpreter.h"
#include "Structure/IpplInfoWrapper.h"
#include "Structure/SDDSParser.h"
#include "Track/Track.h"
#include "Utilities/OpalException.h"
#include "Utilities/Options.h"
#include "Utilities/Timer.h"
#include "Utilities/Util.h"

#include "OPALconfig.h"
#include "Ippl.h"

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/filesystem.hpp>
#include <boost/regex.hpp>
#include <boost/algorithm/string/trim.hpp>

#include <sstream>
#include <string>
#include <iostream>
#include <fstream>

#define LENGTH_DIRECTORY_NAME 29

extern Ippl *ippl;
extern Inform *gmsg;

namespace{

    enum { NORMAL,
           DONE };
}

StatisticalErrors::StatisticalErrors(const Beamline &beamline,
                                     const PartData &reference,
                                     bool revBeam,
                                     bool revTrack,
                                     const std::string &method,
                                     unsigned int nodesPerInstance,
                                     unsigned int numInstances,
                                     const std::vector<std::string> &objectives):
    Tracker(beamline, reference, revBeam, revTrack),
    method_m(method),
    nodesPerInstance_m(nodesPerInstance),
    numInstances_m(numInstances),
    objectives_m(objectives)
{ }


StatisticalErrors::~StatisticalErrors() {

}

void StatisticalErrors::execute() {
    const bool infoState = Options::info;
    Ippl::Comm->barrier();
    formGroups();

    std::string inputFileName;

    int tag = Ippl::Comm->next_tag(IPPL_APP_TAG5, IPPL_APP_CYCLE);
    if (Ippl::myNode() == 0) {
        bool referenceDone = false;
        std::vector<std::string> runningJobs(localMasters_m.size(),"");
        inputFileName = OpalData::getInstance()->getInputFn();

        OpalInputInterpreter interpreter(inputFileName);
        interpreter.replaceString("STATISTICAL-ERRORS\\((.*?),\\s*\\d*\\s*,\\s*\\d*\\s*\\)",
                                  "${1}");

        removeOldDataFiles();

        std::string directoryBaseName = getNextDirectoryName(false);
        try {
            for (unsigned int instance = 0; instance < numInstances_m; ++ instance) {
                std::pair<int, int> source = identifyFreeReplica(tag);
                std::string replicasFormerDirectory = runningJobs[source.second - 1];

                char runNumber[6];
                sprintf(runNumber, "%05u", instance);
                std::string nextDir = directoryBaseName + "_run_" + std::string(runNumber);
                if (!referenceDone) {
                    nextDir = getNextDirectoryName(true);
                }
                runningJobs[source.second - 1] = assignReplicaNewJob(interpreter, source, !referenceDone, nextDir, tag);

                collectOutputData(replicasFormerDirectory);
                removeLinks(replicasFormerDirectory);

                referenceDone = true;
            }
        } catch (OpalException &ex) {
            shutReplicasDown(tag, runningJobs);
            throw OpalException(ex.where(), ex.what());
        } catch (...) {
            shutReplicasDown(tag, runningJobs);
            throw OpalException("StatisticalErrors::execute",
                                "something went wrong");
        }
        shutReplicasDown(tag, runningJobs);

    } else {
        while (true) {
            PRMessage message = getCrunchJobFromPrimary(tag);

            if (message.status == DONE) break;

            std::string nextDir = message.nextDirectory;
            std::string runNumber = "";
            if (nextDir != "reference") {
                runNumber = "_run_" + nextDir.substr(LENGTH_DIRECTORY_NAME - 6, 5);
            }
            std::string inputFn = OpalData::getInstance()->getInputFn();
            size_t startExtension = inputFn.find_last_of(".");
            std::string baseName = inputFn.substr(0, startExtension);
            std::string extension = inputFn.substr(startExtension);
            inputFileName = nextDir + "/" + baseName + runNumber + extension;

            Inform *origGmsg = gmsg;
            gmsg = 0;

            stashEnvironment();
            runSimulation(inputFileName, subComm_m);
            popEnvironment();

            gmsg = origGmsg;
            Options::info = infoState;
        }
    }
    Ippl::Comm->barrier();

    writeSDDSFile();

    removeBinaryDataFiles();
    Ippl::Comm->barrier();
}

void StatisticalErrors::runSimulation(const std::string &inputFileName, MPI_Comm comm) {
    IpplInfoWrapper *newippl = new IpplInfoWrapper(inputFileName, comm);
    std::string::size_type startExtension    = inputFileName.find_last_of('.');
    std::string::size_type startRelativePath = inputFileName.find_last_of('/');
    std::string relativePath("");
    if (startRelativePath != std::string::npos) {
        relativePath = inputFileName.substr(0, startRelativePath + 1);
    }
    std::string outputFileName = inputFileName.substr(0,startExtension) + ".out";
    std::ofstream output(outputFileName.c_str());

    IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("mainTimer");

    gmsg = new Inform("OPAL ", output);
    IpplInfo::Info->setDestination(output);
    IpplInfo::Error->setDestination(output);
    IpplInfo::Warn->setDestination(output);

    OpalData *opal = OpalData::getInstance();

    Configure::configure();
    opal->storeInputFn(inputFileName);

    IpplTimings::startTimer(mainTimer);

    FileStream *is = 0;
    try {
        is = new FileStream(inputFileName);
    } catch(...) {
        is = 0;
        throw new OpalException("StatisticalErrors::runSimulation()", "Could not open inputfile: " + inputFileName);
    }

    OpalParser *parser = new OpalParser();
    if(is) parser->run(is);
    IpplTimings::stopTimer(mainTimer);
    IpplTimings::print(relativePath + "timing.dat");

    Ippl::Comm->barrier();
    processOutputData(inputFileName.substr(0, startExtension) + ".stat");

    OpalData::deleteInstance();
    IpplInfo::Info->setDestination(std::cout);
    IpplInfo::Error->setDestination(std::cout);
    IpplInfo::Warn->setDestination(std::cout);
    // Fieldmap::clearDictionary();
    delete parser;
    // delete is; this is apparently already done in the parser.
    delete newippl;
    delete gmsg;
    gmsg = 0;

    output.close();
}

void StatisticalErrors::stashEnvironment() {
    Ippl::stash();
    IpplTimings::stash();
    Track::stash();
    OpalData::stashInstance();
}

void StatisticalErrors::popEnvironment() {
    Ippl::pop();
    IpplTimings::pop();
    OpalData::popInstance();
    Track::pop();
}

std::string StatisticalErrors::getNextDirectoryName(bool reference) {
    namespace fs = boost::filesystem;

    std::string nameFormat("%%%%%%%%-%%%%-%%%%");
    if (nameFormat.length() + 11 != LENGTH_DIRECTORY_NAME) {
        throw OpalException("StatisticalErrors::getNextDirectoryName",
                            "length of directory name not equal to " + std::to_string(LENGTH_DIRECTORY_NAME));
    }

    fs::path model(nameFormat);
    fs::path path;

    if (reference) {
        path = fs::path("reference");

        if (fs::exists(path)) {
            fs::remove_all(path);
        }
    } else {
        path = fs::unique_path(model);
    }

    return path.native();
}

void StatisticalErrors::createDirectory(const std::string &name) {
    namespace fs = boost::filesystem;

    if (!fs::exists(name))
        fs::create_directory(name);
}

void StatisticalErrors::linkFiles(const std::string &directory) {
    namespace fs = boost::filesystem;

    if (!fs::exists(directory) || !fs::is_directory(fs::path(directory))) return;

    fs::path currentPath = fs::current_path();
    fs::path destination = currentPath;
    destination /= fs::path(directory);

    fs::directory_iterator it(currentPath);
    fs::directory_iterator end;

    std::string baseName = OpalData::getInstance()->getInputBasename();
    std::string inputFile = OpalData::getInstance()->getInputFn();
    std::string statFile = baseName + ".stat";
    std::string lbalFile = baseName + ".lbal";
    std::string H5File = baseName + ".h5";
    std::string outputFile = baseName + ".out";
    std::string dataDirectory = "data";
    std::string referenceDirectory = "reference";
    std::string errorMsg = "errormsg.txt";

    boost::smatch matchResult;
    boost::regex baseNameRe("^" + baseName + "\\..*$");
    boost::regex emacsBackupRe(".*~$");
    boost::regex crunchDirectoryRe("^[0-9a-fA-F]{8}-[0-9a-fA-F]{4}-[0-9a-fA-F]{4}_run_[0-9]{5}$");

    for (; it != end; ++ it) {
        std::string entry = it->path().filename().native();
        if (entry == errorMsg ||
            entry == dataDirectory ||
            entry == referenceDirectory ||
            boost::regex_match(entry, matchResult, crunchDirectoryRe) ||
            boost::regex_match(entry, matchResult, emacsBackupRe)) continue;

        if (boost::regex_match(entry, matchResult, baseNameRe) &&
            (entry == inputFile ||
             entry == statFile ||
             entry == lbalFile ||
             entry == H5File ||
             entry == outputFile)) continue;

        fs::path linksTo = destination;
        linksTo /= fs::path(entry);

        fs::create_symlink(it->path(), linksTo);
    }
}

void StatisticalErrors::removeLinks(const std::string &directory) {
    namespace fs = boost::filesystem;

    if (!fs::exists(directory) || !fs::is_directory(fs::path(directory))) return;

    fs::path currentPath = fs::current_path();
    fs::path source = currentPath;
    source /= fs::path(directory);

    fs::directory_iterator it(source);
    fs::directory_iterator end;

    for (; it != end; ++ it) {
        if (fs::is_symlink(it->path())) {
            fs::remove(it->path());
        }
    }
}

void StatisticalErrors::formGroups() {
    MPI_Group fullGroup;
    unsigned int numNodes = Ippl::getNodes();

    Ippl::Comm->barrier();
    if (Ippl::myNode() == 0) {
        groupRanks_m.push_back(0);

        calcIdealDivision(numNodes);
        std::vector<int> localRanks(numNodes, 0);

        MPI_Bcast(&groupIDs_m[0], numNodes, MPI_INT, 0, Ippl::getComm());

        MPI_Comm_group(Ippl::getComm(), &fullGroup);
        MPI_Group_incl(fullGroup, 1, &groupRanks_m[0], &subGroup_m);
        MPI_Comm_create(Ippl::getComm(), subGroup_m, &subComm_m);

        MPI_Gather(MPI_IN_PLACE,
                   1,
                   MPI_INT,
                   &localRanks[0],
                   1,
                   MPI_INT,
                   0,
                   Ippl::getComm());

        for (unsigned int i = 1; i < numNodes; ++ i) {
            if (localRanks[i] == 0) {
                localMasters_m.push_back(i);
            }
        }
    } else {
        std::vector<int> groupDivision(numNodes, 0);
        MPI_Bcast(&groupDivision[0], numNodes, MPI_INT, 0, Ippl::getComm());
        int myGroupID = groupDivision[Ippl::myNode()];

        for (unsigned int i = 1; i < numNodes; ++ i) {
            if (groupDivision[i] == myGroupID) {
                groupRanks_m.push_back(i);
            }
        }

        MPI_Comm_group(Ippl::getComm(), &fullGroup);
        MPI_Group_incl(fullGroup, groupRanks_m.size(), &groupRanks_m[0], &subGroup_m);
        MPI_Comm_create(Ippl::getComm(), subGroup_m, &subComm_m);

        MPI_Comm_rank(subComm_m, &groupDivision[0]);

        MPI_Gather(&groupDivision[0],
                   1,
                   MPI_INT,
                   &groupDivision[0],
                   1,
                   MPI_INT,
                   0,
                   Ippl::getComm());
    }
}

void StatisticalErrors::calcIdealDivision(unsigned int numNodes) {
    groupIDs_m.resize(numNodes, 0);
    int groupID = 1;
    for (unsigned int i = 1; i < numNodes; i += nodesPerInstance_m, ++ groupID) {
        unsigned int end = std::min(i + nodesPerInstance_m, numNodes);
        for (unsigned int j = i; j < end; ++ j) {
            groupIDs_m[j] = groupID;
        }
    }
}

void StatisticalErrors::shutReplicasDown(int tag,
                                         const std::vector<std::string> &currentJobs) {
    MPI_Status mpiStatus;
    MPI_Request mpiRequest;
    std::vector<char> sendBuffer(sizeof(int) + LENGTH_DIRECTORY_NAME);
    int runStatus = DONE;

    memcpy(&sendBuffer[0], reinterpret_cast<const char*>(&runStatus), sizeof(int));
    for (unsigned int i = 0; i < localMasters_m.size(); ++ i) {
        MPI_Recv(&runStatus, 1, MPI_INT, MPI_ANY_SOURCE, tag, Ippl::getComm(), &mpiStatus);
        MPI_Isend(&sendBuffer[0],
                  sizeof(int) + LENGTH_DIRECTORY_NAME,
                  MPI_CHAR,
                  mpiStatus.MPI_SOURCE,
                  tag,
                  Ippl::getComm(),
                  &mpiRequest);

        unsigned int groupID = groupIDs_m[mpiStatus.MPI_SOURCE];
        std::string replicasFormerDirectory = currentJobs[groupID - 1];

        collectOutputData(replicasFormerDirectory);
        removeLinks(replicasFormerDirectory);
    }
}

std::pair<int, int> StatisticalErrors::identifyFreeReplica(int tag) {
    int runStatus;
    MPI_Status mpiStatus;

    MPI_Recv(&runStatus, 1, MPI_INT, MPI_ANY_SOURCE, tag, Ippl::getComm(), &mpiStatus);
    std::pair<int, int> source;
    source.first = mpiStatus.MPI_SOURCE;
    source.second = groupIDs_m[source.first];

    return source;
}

std::string StatisticalErrors::assignReplicaNewJob(const OpalInputInterpreter &interpreter,
                                                   const std::pair<int, int> &source,
                                                   bool referenceRun,
                                                   std::string nextDir,
                                                   int tag) {

    int runStatus = NORMAL;
    MPI_Request mpiRequest;
    std::vector<char> sendBuffer(sizeof(int) + LENGTH_DIRECTORY_NAME);
    // std::string nextDir = getNextDirectoryName(referenceRun);

    createDirectory(nextDir);
    linkFiles(nextDir);
    std::string runNumber = "";
    if (!referenceRun) {
        runNumber = "_run_" + nextDir.substr(LENGTH_DIRECTORY_NAME - 6, 5);
    }
    std::string inputFn = OpalData::getInstance()->getInputFn();
    size_t startExtension = inputFn.find_last_of(".");
    std::string baseName = inputFn.substr(0, startExtension);
    std::string extension = inputFn.substr(startExtension);
    std::string outputFileName = nextDir + "/" + baseName + runNumber + extension;

    std::ofstream ofh(outputFileName.c_str());

    if (referenceRun) {
        ofh << interpreter.processASTReference();
    } else {
        ofh << interpreter.processAST();
    }

    ofh.close();

    memcpy(&sendBuffer[0], reinterpret_cast<const char*>(&runStatus), sizeof(int));
    strcpy(&sendBuffer[sizeof(int)], nextDir.c_str());

    MPI_Isend(&sendBuffer[0],
              sizeof(int) + LENGTH_DIRECTORY_NAME,
              MPI_CHAR,
              source.first,
              tag,
              Ippl::getComm(),
              &mpiRequest);

    return nextDir;
}

StatisticalErrors::PRMessage StatisticalErrors::getCrunchJobFromPrimary(int tag) {
    MPI_Status mpiStatus;
    MPI_Request mpiRequest;
    std::vector<char> buffer(sizeof(int) + LENGTH_DIRECTORY_NAME);
    PRMessage message;
    int myLocalRank;
    MPI_Comm_rank(subComm_m, &myLocalRank);

    if (myLocalRank == 0) {
        message.status = DONE;
        MPI_Isend(&message.status, 1, MPI_INT, 0, tag, Ippl::getComm(), &mpiRequest);
        MPI_Recv(&buffer[0],
                 sizeof(int) + LENGTH_DIRECTORY_NAME,
                 MPI_CHAR,
                 0,
                 tag,
                 Ippl::getComm(),
                 &mpiStatus);
    }
    MPI_Bcast(&buffer[0],
              sizeof(int) + LENGTH_DIRECTORY_NAME,
              MPI_CHAR,
              0,
              subComm_m);

    message.status = *reinterpret_cast<int*>(&buffer[0]);
    message.nextDirectory = std::string(&buffer[sizeof(int)]);

    return message;
}

void StatisticalErrors::collectOutputData(const std::string &directory) {
    namespace fs = boost::filesystem;

    if (!fs::exists(directory)) return;

    try {
        std::string outdataDirName = "data/";
        std::string indataDirName = directory + "/data/";

        for (std::string col: objectives_m) {
            std::string contents;
            std::string inFname = indataDirName + col + ".dat";
            std::ifstream in(inFname, std::ios::binary);
            std::string outFname = outdataDirName + col + ".dat";
            std::ofstream out;


            in.seekg(0, std::ios::end);
            int sizeFile = in.tellg();
            size_t dataSize = sizeFile - sizeof(int);
            in.seekg(0, std::ios::beg);

            contents.resize(dataSize);

            if (!fs::exists(outFname)) {
                char typeSize[sizeof(int)];
                in.read(typeSize, sizeof(int));

                out.open(outFname, std::ios::binary);
                out.write(typeSize, sizeof(int));
            } else {
                in.seekg(sizeof(int));

                out.open(outFname, std::ios::app | std::ios::binary);
            }

            in.read(&contents[0], contents.size());
            out.write(&contents[0], contents.size());

            out.close();
            in.close();

            fs::remove(inFname);
        }
    } catch (OpalException &ex) {
        ERRORMSG(__FILE__ << ": " << __LINE__ << ": " << ex.where());
    }
}

void StatisticalErrors::processOutputData(const std::string &sddsFileName) {
    namespace fs = boost::filesystem;

    if (Ippl::myNode() != 0 || !fs::exists(sddsFileName)) return;

    // std::string contents = readSDDSFile(sddsFileName);
    SDDS::SDDSParser parser;
    try {
        parser = SDDS::SDDSParser(sddsFileName);
        parser.run();
    } catch (OpalException &ex) {
        ERRORMSG(__FILE__ << ": " << __LINE__ << ": could not parse sdds file '" << sddsFileName << "'\n" << endl);
    }

    try {
        size_t startExtension = sddsFileName.find_last_of('/');
        std::string dataDirName = sddsFileName.substr(0, startExtension) + "/data/";
        fs::create_directory(dataDirName);

        for (std::string col: objectives_m) {

            SDDS::ast::columnData_t values = parser.getColumnData(col);
            int length = values.size();
            std::string fname = dataDirName + col + ".dat";
            std::ofstream out(fname, std::ios::binary);

            boost::apply_visitor(DataTypeWriter(out), values.front());
            out.write(reinterpret_cast<const char*>(&length), sizeof(int));

            BinaryWriter converter(out);
            std::for_each(values.begin(), values.end(),
                          boost::apply_visitor(converter));

            out.close();
        }
    } catch (OpalException &ex) {
        ERRORMSG(__FILE__ << ": " << __LINE__ << ": " << ex.where());
    }
}

void StatisticalErrors::writeSDDSFile() {

    if (Ippl::myNode() == 0) {
        std::map<std::string, std::vector<double> > data = distributeProcessingJobs();

        std::string statName = OpalData::getInstance()->getInputBasename() + ".stat";
        std::ofstream statOut(statName);
        std::vector<std::string> order;

        writeSDDSHeader(statOut, order);
        writeSDDSData(statOut, data, order);
    } else {
        getProcessingJobFromPrimary();
    }
}

void StatisticalErrors::writeSDDSHeader(std::ofstream &out,
                                        std::vector<std::string> &order) {
    OPALTimer::Timer simtimer;
    std::string referenceStatFileName = "reference/" + OpalData::getInstance()->getInputBasename() + ".stat";
    SDDS::SDDSParser parser(referenceStatFileName);
    SDDS::file parsedRefData = parser.run();
    SDDS::columnList &columns = parsedRefData.sddsColumns_m;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());
    out << "SDDS1\n"
        << "&description text=\"Statistics data " << OpalData::getInstance()->getInputFn()
        << " " << dateStr << " " << timeStr << "\",\n contents=\"stat parameters\" &end\n"
        << "&parameter name=processors, type=long, description=\"Number of Cores used\" &end\n"
        << "&parameter name=revision, type=string, description=\"git revision of opal\" &end\n";

    out << "&column name=s, type=double, units=m, description=\"1 longitudinal position\" &end\n";
    order.push_back("s");

    for (unsigned int i = 0; i < objectives_m.size(); ++ i) {
        if (objectives_m[i] == "s") continue;
        std::string curObjective = objectives_m[i];

        SDDS::columnList::iterator it = columns.begin();
        for (; it != columns.end(); ++ it) {
            if (*(*it).name_m == curObjective) break;
        }
        if (it == columns.end()) continue;

        out << "&column name=mean_" << curObjective << ", type=double, units="
            << *(*it).units_m << ", description=\"" << 2 * order.size() << " mean of " << curObjective << " over all runs\" &end\n"
            << "&column name=rms_" << curObjective << ", type=double, units="
            << *(*it).units_m << ", description=\"" << 2 * order.size() + 1 << " rms of " << curObjective << " over all runs\" &end\n";

        order.push_back(curObjective);
    }

    out << "&data mode=ascii, no_row_counts=1 &end\n";
}

void StatisticalErrors::writeSDDSData(std::ofstream &out,
                                      const std::map<std::string, std::vector<double> > &data,
                                      const std::vector<std::string> &order) {
    out << Ippl::getNodes() << "\n"
        << PACKAGE_NAME << " " << OPAL_VERSION_STR << " # git rev. " << Util::getGitRevision() << std::endl;

    size_t sizeData = data.find("s")->second.size();
    size_t numColumns = order.size();
    for (unsigned int i = 0; i < sizeData; ++ i) {
        std::string colName = "s";
        auto it = data.find(colName);
        if (it == data.end()) break;
        out << std::setw(13) << std::setprecision(6) << (*it).second[i];

        for (unsigned int j = 0; j < numColumns; ++ j) {
            colName = "mean_" + order[j];
            it = data.find(colName);
            if (it == data.end()) continue;
            out << std::setw(13) << std::setprecision(6) << (*it).second[i];

            colName = "rms_" + order[j];
            it = data.find(colName);
            out << std::setw(13) << std::setprecision(6) << (*it).second[i];
        }
        out << "\n";
    }
}

std::map<std::string, std::vector<double> > StatisticalErrors::distributeProcessingJobs() {
    std::map<std::string, std::vector<double> > statistics;

    int tag = Ippl::Comm->next_tag(IPPL_APP_TAG5, IPPL_APP_CYCLE);
    if (Ippl::myNode() == 0) {
        DataSource objectiveData;
        std::vector<double> positionSampling = getReferenceSamplingPositions();
        sendReferenceSamplingPositions(positionSampling);

        statistics.insert(std::make_pair("s", positionSampling));

        size_t bufferSize = 2 * (sizeof(int) + positionSampling.size() * sizeof(double));
        std::vector<char> messageBuffer(bufferSize);
        unsigned int numMessages = objectives_m.size() + Ippl::getNodes() - 1;

        for (unsigned int i = 0; i < numMessages; ++ i) {
            if (i < objectives_m.size() && objectives_m[i] == "s") continue;
            objectiveData = getProcessedDataFromReplica(tag, i, messageBuffer.size(), &messageBuffer[0]);

            if (objectiveData.data.size() > 0) {
                std::string statName = "mean_" + objectives_m[objectiveData.ID];
                std::vector<double> meanData(objectiveData.data.begin(),
                                             objectiveData.data.begin() + positionSampling.size());
                statistics.insert(std::make_pair(statName, meanData));

                statName = "rms_" + objectives_m[objectiveData.ID];
                std::vector<double> statData(objectiveData.data.begin() + positionSampling.size(),
                                             objectiveData.data.end());
                statistics.insert(std::make_pair(statName, statData));
            }

            sendObjectiveForProcessing(tag, objectiveData.source, i, &messageBuffer[0]);
        }
    }

    return statistics;
}

void StatisticalErrors::getProcessingJobFromPrimary() {

    if (Ippl::myNode() == 0) throw OpalException("StatisticalErrors::receiveProcessingJob",
                                                 "Primary should not be here");
    int tag = Ippl::Comm->next_tag(IPPL_APP_TAG5, IPPL_APP_CYCLE);
    int objID = -1;
    std::vector<char> messageBuffer;
    std::vector<std::pair<double, double> > statistics;
    std::vector<double> pdata = getReferenceSamplingPositions();
    int length = pdata.size();

    messageBuffer.resize(2 * (sizeof(int) + length * sizeof(double)));
    while (true) {
        if (objID > 0) copyDataToBuffer(&messageBuffer[2 * sizeof(int)], statistics);
        sendProcessedDataToPrimary(tag, objID, length, &messageBuffer[0]);

        std::pair<int, int> nextJob = getNextObjectiveToProcess(tag, &messageBuffer[0]);
        if (nextJob.first == DONE) break;
        objID = nextJob.second;

        statistics = computeStatistics(objectives_m[objID], pdata);
    }
}

std::vector<std::pair<double, double> >
StatisticalErrors::computeStatistics(const std::string &objective,
                                     const std::vector<double> &referencePositions) {
    namespace fs = boost::filesystem;

    std::string posName = "data/s.dat";
    std::string inName = "data/" + objective + ".dat";
    if (!fs::exists(posName)) {
        throw OpalException("StatisticalErrors::computeStatistics",
                            "can't find file with sampling positions");
    }
    if (!fs::exists(inName)) {
        throw OpalException("StatisticalErrors::computeStatistics",
                            "can't find file with data of objective '" + objective + "'");
    }

    std::ifstream in(inName, std::ios::binary);
    std::ifstream pos(posName, std::ios::binary);

    char testBitData, testBitPos;
    int dataType, numDataSets = 0;
    std::vector<double> dataSum, squaredDataSum;
    std::vector<std::pair<double, double>> momenta;
    dataSum.resize(referencePositions.size(), 0.0);
    squaredDataSum.resize(referencePositions.size(), 0.0);
    momenta.resize(referencePositions.size(), std::make_pair(0.0, 0.0));

    pos.seekg(sizeof(int));
    in.read(reinterpret_cast<char*>(&dataType), sizeof(int));

    pos.get(testBitPos);
    in.get(testBitData);

    while (!in.eof() && !pos.eof()) {
        in.putback(testBitData);
        pos.putback(testBitPos);

        std::vector<double> currentData = getNextDataset(in, pos, dataType, referencePositions);
        for (unsigned int i = 0; i < dataSum.size(); ++ i) {
            dataSum[i] += currentData[i];
            squaredDataSum[i] += std::pow(currentData[i], 2.0);
        }

        ++ numDataSets;

        pos.get(testBitPos);
        in.get(testBitData);
    }

    for (unsigned int i = 0; i < dataSum.size(); ++ i) {
        momenta[i].first = dataSum[i] / numDataSets;
        momenta[i].second = sqrt(std::max(0.0, squaredDataSum[i] / numDataSets - std::pow(momenta[i].first, 2.0)));
    }

    return momenta;
}

std::vector<double> StatisticalErrors::interpolateSDDSData(const std::vector<double> &data,
                                                           const std::vector<double> &oldPositions,
                                                           const std::vector<double> &newPositions) {
    unsigned int j = 0;
    size_t newSize = newPositions.size(), oldSize = oldPositions.size();
    std::vector<double> interpolated(newSize);

    for (unsigned int i = 0; i < newSize; ++ i) {
        while (j < oldSize && oldPositions[j] <= newPositions[i]) {
            ++ j;
        }
        if (j == 0 || j == oldSize) {
            interpolated[i] = 0.0;
        } else {
            double lower = data[j - 1];
            double upper = data[j];
            double tau = (newPositions[i] - oldPositions[j - 1]) / (oldPositions[j] - oldPositions[j - 1]);

            interpolated[i] = lower + tau * (upper - lower);
        }
    }

    return interpolated;
}

std::vector<double> StatisticalErrors::readFloatData(std::ifstream &in) {
    int dataLength;
    std::vector<double> data;
    float singleSampling;

    in.read(reinterpret_cast<char*>(&dataLength), sizeof(int));
    data.resize(dataLength);
    for (int i = 0; i < dataLength; ++ i) {
        in.read(reinterpret_cast<char*>(&singleSampling), sizeof(float));
        data[i] = singleSampling;
    }

    return data;
}

std::vector<double> StatisticalErrors::readDoubleData(std::ifstream &in) {
    int dataLength;
    std::vector<double> data;
    double singleSampling;

    in.read(reinterpret_cast<char*>(&dataLength), sizeof(int));
    data.resize(dataLength);
    for (int i = 0; i < dataLength; ++ i) {
        in.read(reinterpret_cast<char*>(&singleSampling), sizeof(double));
        data[i] = singleSampling;
    }

    return data;
}

std::vector<double> StatisticalErrors::readShortData(std::ifstream &in) {
    int dataLength;
    std::vector<double> data;
    short singleSampling;

    in.read(reinterpret_cast<char*>(&dataLength), sizeof(int));
    data.resize(dataLength);
    for (int i = 0; i < dataLength; ++ i) {
        in.read(reinterpret_cast<char*>(&singleSampling), sizeof(short));
        data[i] = singleSampling;
    }

    return data;
}

std::vector<double> StatisticalErrors::readLongData(std::ifstream &in) {
    int dataLength;
    std::vector<double> data;
    long singleSampling;

    in.read(reinterpret_cast<char*>(&dataLength), sizeof(int));
    data.resize(dataLength);
    for (int i = 0; i < dataLength; ++ i) {
        in.read(reinterpret_cast<char*>(&singleSampling), sizeof(long));
        data[i] = singleSampling;
    }

    return data;
}

void StatisticalErrors::removeOldDataFiles() {
    removeBinaryDataFiles();
}

void StatisticalErrors::removeBinaryDataFiles() {
    namespace fs = boost::filesystem;

    if (Ippl::myNode() > 0) return;

    for (std::string col: objectives_m) {
        std::string fname = "data/" + col + ".dat";

        if (fs::exists(fname)) {
            fs::remove(fname);
        }
    }
}

std::vector<double> StatisticalErrors::getNextDataset(std::ifstream &dataStream,
                                                      std::ifstream &positionStream,
                                                      int dataType,
                                                      const std::vector<double> &referencePositions) {
    std::vector<double> currentData;
    std::vector<double> currentPositions = readDoubleData(positionStream);

    switch(dataType) {
    case SDDS::ast::FLOAT:
        currentData = readFloatData(dataStream);
        break;
    case SDDS::ast::DOUBLE:
        currentData = readDoubleData(dataStream);
        break;
    case SDDS::ast::SHORT:
        currentData = readShortData(dataStream);
        break;
    case SDDS::ast::LONG:
        currentData = readLongData(dataStream);
        break;
    default:
        throw OpalException("StatisticalErrors::getNextDataset",
                            "can't compute statistics of char or string data type");
    }

    return interpolateSDDSData(currentData,
                               currentPositions,
                               referencePositions);
}

std::vector<double> StatisticalErrors::getReferenceSamplingPositions() {
    std::vector<double> pdata;
    if (Ippl::myNode() == 0) {
        std::string referenceStatFileName = "reference/" + OpalData::getInstance()->getInputBasename() + ".stat";
        SDDS::SDDSParser parser(referenceStatFileName);
        SDDS::file parsedRefData = parser.run();
        SDDS::ast::columnData_t positions = parser.getColumnData("s");
        int length = positions.size();
        pdata.resize(length);

        auto it = positions.begin();
        for (int i = 0; i < length; ++ i, ++ it) {
            pdata[i] = boost::get<double>(*it);
        }
    } else {
        int length;

        MPI_Bcast(&length, 1, MPI_INT, 0, Ippl::getComm());
        pdata.resize(length, 0.0);
        MPI_Bcast(&pdata[0], length, MPI_DOUBLE, 0, Ippl::getComm());
    }

    return pdata;
}

void StatisticalErrors::sendReferenceSamplingPositions(std::vector<double> positionSampling) {
    int length = positionSampling.size();

    MPI_Bcast(&length, 1, MPI_INT, 0, Ippl::getComm());
    MPI_Bcast(&positionSampling[0], length, MPI_DOUBLE, 0, Ippl::getComm());
}

StatisticalErrors::DataSource StatisticalErrors::getProcessedDataFromReplica(int tag,
                                                                             unsigned int i,
                                                                             unsigned int length,
                                                                             char* messageBuffer) {
    MPI_Status mpiStatus;
    DataSource objectiveData;
    objectiveData.ID = -1;

    MPI_Recv(messageBuffer,
             2 * sizeof(int),
             MPI_CHAR,
             MPI_ANY_SOURCE,
             tag,
             Ippl::getComm(),
             &mpiStatus);
    objectiveData.source = mpiStatus.MPI_SOURCE;
    memcpy(reinterpret_cast<char*>(&objectiveData.ID), messageBuffer + sizeof(int), sizeof(int));

    if (objectiveData.ID >= 0) {
        unsigned int k = 2 * sizeof(int);
        MPI_Recv(messageBuffer + 2 * sizeof(int),
                 2 * length * sizeof(double),
                 MPI_CHAR,
                 objectiveData.source,
                 tag,
                 Ippl::getComm(),
                 &mpiStatus);

        objectiveData.data.resize(2 * length);
        for (unsigned int j = 0; j < 2 * length; ++ j, k += sizeof(double)) {
            memcpy(reinterpret_cast<char*>(&objectiveData.data[j]), messageBuffer + k, sizeof(double));
        }
    }
    return objectiveData;
}

void StatisticalErrors::sendProcessedDataToPrimary(int tag,
                                                   int objectiveID,
                                                   unsigned int length,
                                                   char *messageBuffer) {
    MPI_Request mpiRequest;
    int runStatus = DONE;
    memcpy(messageBuffer, reinterpret_cast<char*>(&runStatus), sizeof(int));
    memcpy(messageBuffer + sizeof(int), reinterpret_cast<char*>(&objectiveID), sizeof(int));

    MPI_Isend(messageBuffer,
              2 * sizeof(int),
              MPI_CHAR,
              0,
              tag,
              Ippl::getComm(),
              &mpiRequest);

    if (objectiveID >= 0) {
        MPI_Isend(messageBuffer + 2 * sizeof(int),
                  2 * length * sizeof(double),
                  MPI_CHAR,
                  0,
                  tag,
                  Ippl::getComm(),
                  &mpiRequest);
    }
}

void StatisticalErrors::copyDataToBuffer(char *messageBuffer,
                                         const std::vector<std::pair<double, double> > &statistics) {
    char *mptr = messageBuffer;

    unsigned int length = statistics.size();
    for (unsigned int i = 0; i < length; ++ i, mptr += sizeof(double)) {
        memcpy(mptr, reinterpret_cast<const char*>(&statistics[i].first), sizeof(double));
    }
    for (unsigned int i = 0; i < length; ++ i, mptr += sizeof(double)) {
        memcpy(mptr, reinterpret_cast<const char*>(&statistics[i].second), sizeof(double));
    }

}

void StatisticalErrors::sendObjectiveForProcessing(int tag,
                                                   int source,
                                                   unsigned int i,
                                                   char* messageBuffer) {
    MPI_Request mpiRequest;
    int runStatus = (i < objectives_m.size()) ? NORMAL: DONE;

    memcpy(messageBuffer, reinterpret_cast<char*>(&runStatus), sizeof(int));
    memcpy(messageBuffer + sizeof(int), reinterpret_cast<char*>(&i), sizeof(int));

    MPI_Isend(messageBuffer,
              2 * sizeof(int),
              MPI_CHAR,
              source,
              tag,
              Ippl::getComm(),
              &mpiRequest);

}

std::pair<int, int> StatisticalErrors::getNextObjectiveToProcess(int tag, char *messageBuffer) {
    MPI_Status mpiStatus;
    std::pair<int, int> nextJob;

    MPI_Recv(messageBuffer,
             2 * sizeof(int),
             MPI_CHAR,
             0,
             tag,
             Ippl::getComm(),
             &mpiStatus);


    memcpy(reinterpret_cast<char*>(&nextJob.first), messageBuffer, sizeof(int));
    memcpy(reinterpret_cast<char*>(&nextJob.second), messageBuffer + sizeof(int), sizeof(int));

    return nextJob;
}