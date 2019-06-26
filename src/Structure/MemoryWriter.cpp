#include "MemoryWriter.h"

#include "OPALconfig.h"
#include "Utilities/Util.h"
#include "Utilities/Timer.h"
#include "AbstractObjects/OpalData.h"
#include "Ippl.h"

MemoryWriter::MemoryWriter(const std::string& fname, bool restart)
    : SDDSWriter(fname, restart)
{ }


void MemoryWriter::fillHeader_m() {

    static bool isNotFirst = false;
    if ( isNotFirst ) {
        return;
    }
    isNotFirst = true;

    this->addColumn("t", "double", "ns", "Time");

    this->addColumn("s", "double", "m", "Path length");

    this->addColumn("memory", "double", memory->getUnit(), "Total Memory");

    for (int p = 0; p < Ippl::getNodes(); ++p) {
        std::stringstream tmp1;
        tmp1 << "processor-" << p;

        std::stringstream tmp2;
        tmp2 << "Memory per processor " << p;
        this->addColumn(tmp1.str(), "double", memory->getUnit(), tmp2.str());
    }

    if ( mode_m == std::ios::app )
        return;

    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());

    IpplMemoryUsage::IpplMemory_p memory = IpplMemoryUsage::getInstance();

    std::stringstream ss;

    ss << "Memory statistics '"
       << OpalData::getInstance()->getInputFn() << "' "
       << dateStr << "" << timeStr;

    this->addDescription(ss.str(), "memory parameters");

    this->addDefaultParameters();


    this->addInfo("ascii", 1);
}


void MemoryWriter::write(PartBunchBase<double, 3> *beam)
{
    IpplMemoryUsage::IpplMemory_p memory = IpplMemoryUsage::getInstance();
    memory->sample();

    if (Ippl::myNode() != 0) {
        return;
    }

    double  pathLength = 0.0;
    if (OpalData::getInstance()->isInOPALCyclMode())
        pathLength = beam->getLPath();
    else
        pathLength = beam->get_sPos();

    fillHeader_m();

    this->open();

    this->writeHeader();

    SDDSDataRow row(this->getColumnNames());
    row.addColumn("t", beam->getT() * 1e9);             // 1
    row.addColumn("s", pathLength);                     // 2

    int nProcs = Ippl::getNodes();
    double total = 0.0;
    for (int p = 0; p < nProcs; ++p) {
        total += memory->getMemoryUsage(p);
    }

    row.addColumn("memory", total);

    for (int p = 0; p < nProcs; p++) {
        std::stringstream ss;
        ss << "processor-" << p;
        row.addColumn(ss.str(),  memory->getMemoryUsage(p));
    }

    this->writeRow(row);

    this->newline();

    this->close();
}
