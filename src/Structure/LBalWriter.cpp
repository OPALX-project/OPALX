#include "LBalWriter.h"

#include "OPALconfig.h"
#include "Utilities/Util.h"
#include "Utilities/Timer.h"


LBalWriter::LBalWriter(const std::string& fname, bool restart)
    : SDDSWriter(fname, restart)
{ }


void LBalWriter::fillHeader_m(PartBunchBase<double, 3> *beam) {
    if ( mode_m == std::ios::app )
        return;

    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());

    std::stringstream ss;
    ss << "Processor statistics '"
       << OpalData::getInstance()->getInputFn() << "' "
       << dateStr << "" << timeStr;

    this->addDescription(ss.str(), "lbal parameters");

    this->addDefaultParameters();

    this->addColumn("t", "double", "ns", "Time");

    for (int p = 0; p < Ippl::getNodes(); ++p) {
        std::stringstream tmp1;
        tmp1 << "processor-" << p;

        std::stringstream tmp2;
        tmp2 << "Number of particles of processor " << p;

        this->addColumn(ss.str(), "long", "1", tmp2.str());
    }

#ifdef ENABLE_AMR
    if ( AmrPartBunch* amrbeam = dynamic_cast<AmrPartBunch*>(beam) ) {

        int nLevel = (amrbeam->getAmrObject())->maxLevel() + 1;

        for (int lev = 0; lev < nLevel; ++lev) {
            std::stringstream tmp1;
            tmp1 << "level-" << lev;

            std::stringstream tmp2;
            tmps2 << "Number of particles at level " << lev;
            this->addColumn(tmp1.str(), "long", "1", tmp2.str());
        }
    }
#endif

    this->addInfo("ascii", 1);
}


void LBalWriter::write(PartBunchBase<double, 3> *beam) {

    beam->gatherLoadBalanceStatistics();

#ifdef ENABLE_AMR
    if ( AmrPartBunch* amrbeam = dynamic_cast<AmrPartBunch*>(beam) ) {
        amrbeam->gatherLevelStatistics();
    }
#endif

    if ( Ippl::myNode() != 0 )
        return;

    this->writeValue(beam->getT() * 1e9);   // 1

    size_t nProcs = Ippl::getNodes();
    for (size_t p = 0; p < nProcs; ++ p) {
        this->writeValue(beam->getLoadBalance(p));
    }

#ifdef ENABLE_AMR
    if ( AmrPartBunch* amrbeam = dynamic_cast<AmrPartBunch*>(beam) ) {
        int nLevel = (amrbeam->getAmrObject())->maxLevel() + 1;
        for (int lev = 0; lev < nLevel; ++lev) {
            this->writeValue(amrbeam->getLevelStatistics(lev));
        }
    }
#endif
    this->newline();
}