#include "GridLBalWriter.h"

void GridLBalWriter::fillHeader(PartBunchBase<double, 3> *beam) {

    static bool isNotFirst = false;
    if ( isNotFirst ) {
        return;
    }
    isNotFirst = true;

    columns_m.addColumn("t", "double", "ns", "Time");

    int nLevel = (amrbeam->getAmrObject())->maxLevel() + 1;

    for (int lev = 0; lev < nLevel; ++lev) {
        std::stringstream tmp1;
        tmp1 << "level-" << lev;

        std::stringstream tmp2;
        tmp2 << "Number of boxes at level " << lev;

        columns_m.addColumn(tmp1.str(), "long", "1", tmp2.str());
    }

    for (int p = 0; p < Ippl::getNodes(); ++p) {
        std::stringstream tmp1;
        tmp1 << "processor-" << p;

        std::stringstream tmp2;
        tmp2 << "Number of grid points per processor " << p;

        columns_m.addColumn(tmp1.str(), "long", "1", tmp2.str());
    }

    if ( mode_m == std::ios::app )
        return;

    AmrPartBunch* amrbeam = dynamic_cast<AmrPartBunch*>(beam);

    if ( !amrbeam )
        throw OpalException("DataSink::writeGridLBalHeader()",
                            "Can not write grid load balancing for non-AMR runs.");

    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());

    std::stringstream ss;
    ss << "Grid load balancing statistics '"
       << OpalData::getInstance()->getInputFn() << "' "
       << dateStr << "" << timeStr;

    this->addDescription(ss.str(), "grid lbal parameters");

    this->addDefaultParameters();

    this->addInfo("ascii", 1);
}


void GridLBalWriter::writeData(PartBunchBase<double, 3> *beam) {
    AmrPartBunch* amrbeam = dynamic_cast<AmrPartBunch*>(beam);

    if ( !amrbeam )
        throw OpalException("DataSink::writeGridLBalData()",
                            "Can not write grid load balancing for non-AMR runs.");

    amrbeam->getAmrObject()->getGridStatistics(gridPtsPerCore, gridsPerLevel);

    if ( Ippl::myNode() != 0 )
        return;

    this->fillHeader(beam);

    this->open();

    this->writeHeader();

    columns_m.addColumnValue("t", beam->getT() * 1e9); // 1

    std::map<int, long> gridPtsPerCore;

    int nLevel = (amrbeam->getAmrObject())->maxLevel() + 1;
    std::vector<int> gridsPerLevel;

    for (int lev = 0; lev < nLevel; ++lev) {
        std::stringstream ss;
        ss << "level-" << lev;
        columns_m.addColumnValue(ss.str(), gridsPerLevel[lev]);
    }

    int nProcs = Ippl::getNodes();
    for (int p = 0; p < nProcs; ++p) {
        std::stringstream ss;
        ss << "processor-" << p;
        columns_m.addColumnValue(ss.str(), gridPtsPerCore[p]);
    }

    this->writeRow();

    this->newline();

    this->close();
}
