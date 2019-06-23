#include "GridLBalWriter.h"

void GridLBalWriter::writeHeader() {
    AmrPartBunch* amrbeam = dynamic_cast<AmrPartBunch*>(beam);

    if ( !amrbeam )
        throw OpalException("DataSink::writeGridLBalHeader()",
                            "Can not write grid load balancing for non-AMR runs.");

    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());

    os_m << "SDDS1" << std::endl;
    
    std::stringstream ss;
    ss << "Grid load balancing statistics '"
       << OpalData::getInstance()->getInputFn() << "' "
       << dateStr << "" << timeStr;
    
    this->addDescription(ss.str(), "grid lbal parameters");
    
    this->addParameter("processors", "long", "Number of Cores used");
    
    this->addParameter("revision", "string", "git revision of opal");
    
    this->addParameter("flavor", "string", "OPAL flavor that wrote file");

    this->addColumn("t", "double", "ns", "Time");

    int nLevel = (amrbeam->getAmrObject())->maxLevel() + 1;

    for (int lev = 0; lev < nLevel; ++lev) {
        std::stringstream tmp1;
        tmp1 << "level-" << lev;
        
        std::stringstream tmp2;
        tmp2 << "Number of boxes at level " << lev;
        
        this->addColumn(tmp1.str(), "long", "1", tmp2.str());
    }

    for (int p = 0; p < Ippl::getNodes(); ++p) {
        std::stringstream tmp1;
        tmp1 << "processor-" << p;
        
        std::stringstream tmp2;
        tmp2 << "Number of grid points per processor " << p;
        
        this->addColumn(tmp1.str(), "long", "1", tmp2.str());
    }

    this->addInfo("ascii", 1);

    os_m << Ippl::getNodes() << std::endl
         << OPAL_PROJECT_NAME << " " << OPAL_PROJECT_VERSION << " git rev. #" << Util::getGitRevision() << std::endl
         << (OpalData::getInstance()->isInOPALTMode()? "opal-t":
                (OpalData::getInstance()->isInOPALCyclMode()? "opal-cycl": "opal-env")) << std::endl;
}


void GridLBalWriter::writeData(PartBunchBase<double, 3> *beam) {
    AmrPartBunch* amrbeam = dynamic_cast<AmrPartBunch*>(beam);

    if ( !amrbeam )
        throw OpalException("DataSink::writeGridLBalData()",
                            "Can not write grid load balancing for non-AMR runs.");

    this->writeValue(amrbeam->getT() * 1e9);    // 1

    std::map<int, long> gridPtsPerCore;

    int nLevel = (amrbeam->getAmrObject())->maxLevel() + 1;
    std::vector<int> gridsPerLevel;

    amrbeam->getAmrObject()->getGridStatistics(gridPtsPerCore, gridsPerLevel);

    for (int lev = 0; lev < nLevel; ++lev) {
        this->writeValue(gridsPerLevel[lev]);
    }

    int nProcs = Ippl::getNodes();
    for (int p = 0; p < nProcs; ++p) {
        this->writeValue(gridPtsPerCore[p]);
    }
}
