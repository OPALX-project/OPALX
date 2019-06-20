#include "MemoryWriter.h"

void MemoryWriter::writeHeader() {
    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());

    IpplMemoryUsage::IpplMemory_p memory = IpplMemoryUsage::getInstance();

    os_m << "SDDS1" << std::endl;
    
    std::stringstream ss;
    
    ss << "Memory statistics '"
       << OpalData::getInstance()->getInputFn() << "' "
       << dateStr << "" << timeStr;
    
    this->addDescription(ss.str(), "memory parameters");
    
    this->addParameter("processors", "long", "Number of Cores used");
    
    this->addParameter("revision", "string", "git revision of opal");
    
    this->addParameter("flavor", "string", "OPAL flavor that wrote file");
    
    this->addColumn("t", "double", "ns", "Time");
    
    this->addColumn("memory", "double", memory->getUnit(), "Total Memory");

    for (int p = 0; p < Ippl::getNodes(); ++p) {
        std::stringstream tmp1;
        tmp1 << "processor-" << p;
        
        std::stringstream tmp2;
        tmp2 << "Memory per processor " << p;
        this->addColumn(tmp1.str(), "double", memory->getUnit(), tmp2.str());
    }
    
    this->addData("ascii", 1);
    
    os_m << Ippl::getNodes() << std::endl
         << OPAL_PROJECT_NAME << " " << OPAL_PROJECT_VERSION << " git rev. #" << Util::getGitRevision() << std::endl
         << (OpalData::getInstance()->isInOPALTMode()? "opal-t":
                (OpalData::getInstance()->isInOPALCyclMode()? "opal-cycl": "opal-env")) << std::endl;
}


void MemoryWriter::writeData(PartBunchBase<double, 3> *beam) {
    
    this->writeValue(beam->getT() * 1e9);   // 1

    IpplMemoryUsage::IpplMemory_p memory = IpplMemoryUsage::getInstance();

    int nProcs = Ippl::getNodes();
    double total = 0.0;
    for (int p = 0; p < nProcs; ++p) {
        total += memory->getMemoryUsage(p);
    }

    this->writeValue(total);
    
    for (int p = 0; p < nProcs; p++) {
        this->writeValue(memory->getMemoryUsage(p));
}
