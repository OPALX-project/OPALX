#include "MemoryProfiler.h"

#ifdef __linux__
    #include <sys/types.h>
    #include <unistd.h>
#endif

#include "OPALconfig.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Timer.h"
#include "Utilities/Util.h"
#include "Utilities/OpalException.h"

#include <boost/filesystem.hpp>

#include "Ippl.h"

#include <iomanip>


MemoryProfiler::MemoryProfiler(const std::string& fname, bool restart)
    : SDDSWriter(fname, restart)
{
    procinfo_m = {
        {"VmPeak:", VirtualMemory::VMPEAK},
        {"VmSize:", VirtualMemory::VMSIZE},
        {"VmHWM:",  VirtualMemory::VMHWM},
        {"VmRSS:",  VirtualMemory::VMRSS},
        {"VmStk:",  VirtualMemory::VMSTK},
        {"VmData:", VirtualMemory::VMDATA},
        {"VmExe:",  VirtualMemory::VMEXE},
        {"VmLck:",  VirtualMemory::VMLCK},
        {"VmPin:",  VirtualMemory::VMPIN},
        {"VmLib:",  VirtualMemory::VMLIB},
        {"VmPTE:",  VirtualMemory::VMPTE},
        {"VmPMD:",  VirtualMemory::VMPMD},
        {"VmSwap:", VirtualMemory::VMSWAP}
    };

    vmem_m.resize(procinfo_m.size());
    unit_m.resize(procinfo_m.size());
}


MemoryProfiler::~MemoryProfiler() {
    if ( ofstream_m.is_open() )
        ofstream_m.close();
}


void MemoryProfiler::header_m() {
    if ( mode_m == std::ios::app )
        return;

    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());

    std::stringstream ss;
    ss << "Memory statistics '"
       << OpalData::getInstance()->getInputFn() << "' "
       << dateStr << "" << timeStr;

    this->addDescription(ss.str(), "memory info");

    this->addDefaultParameters();

    this->addColumn("t", "double", "ns", "Time");

    this->addColumn("s", "double", "m", "Path length");

    // peak virtual memory size
    this->addColumn("VmPeak-Min", "double", unit_m[VirtualMemory::VMPEAK],
                    "Minimum peak virtual memory size");

    this->addColumn("VmPeak-Max", "double", unit_m[VirtualMemory::VMPEAK],
                    "Maximum peak virtual memory size");

    this->addColumn("VmPeak-Avg", "double", unit_m[VirtualMemory::VMPEAK],
                    "Average peak virtual memory size");

    // virtual memory size
    this->addColumn("VmSize-Min", "double", unit_m[VirtualMemory::VMSIZE],
                    "Minimum virtual memory size");

    this->addColumn("VmSize-Max", "double", unit_m[VirtualMemory::VMSIZE],
                    "Maximum virtual memory size");

    this->addColumn("VmSize-Avg", "double", unit_m[VirtualMemory::VMSIZE],
                    "Average virtual memory size");

    // peak resident set size ("high water mark")
    this->addColumn("VmHWM-Min", "double", unit_m[VirtualMemory::VMHWM],
                    "Minimum peak resident set size");

    this->addColumn("VmHWM-Max", "double", unit_m[VirtualMemory::VMHWM],
                    "Maximum peak resident set size");

    this->addColumn("VmHWM-Avg", "double", unit_m[VirtualMemory::VMHWM],
                    "Average peak resident set size");

    // resident set size
    this->addColumn("VmRSS-Min", "double", unit_m[VirtualMemory::VMRSS],
                    "Minimum resident set size");

    this->addColumn("VmRSS-Max", "double", unit_m[VirtualMemory::VMRSS],
                    "Maximum resident set size");

    this->addColumn("VmRSS-Avg", "double", unit_m[VirtualMemory::VMRSS],
                    "Average resident set size");

    // stack size
    this->addColumn("VmStk-Min", "double", unit_m[VirtualMemory::VMSTK],
                    "Minimum stack size");

    this->addColumn("VmStk-Max", "double", unit_m[VirtualMemory::VMSTK],
                    "Maximum stack size");

    this->addColumn("VmStk-Avg", "double", unit_m[VirtualMemory::VMSTK],
                    "Average stack size");

    this->addInfo("ascii", 1);
}


void MemoryProfiler::update_m() {
#ifdef __linux__
    static pid_t pid = getpid();
    std::string fname = "/proc/" + std::to_string(pid) + "/status";

    if ( !boost::filesystem::exists(fname) ) {
        throw OpalException("MemoryProfiler::update_m()",
                            "File '" + fname + "' doesn't exist.");
    }

    std::ifstream ifs(fname.c_str());

    if ( !ifs.is_open() ) {
        throw OpalException("MemoryProfiler::update_m()",
                            "Failed to open '" + fname + "'.");
    }

    std::string token = "";
    while (ifs >> token) {
        if ( !procinfo_m.count(token) ) {
            continue;
        }
        int idx = procinfo_m[token];
        ifs >> vmem_m[idx]
            >> unit_m[idx];
    }

    ifs.close();
#endif
}


void MemoryProfiler::compute_m(vm_t& vmMin,
                               vm_t& vmMax,
                               vm_t& vmAvg)
{
    if ( Ippl::getNodes() == 1 ) {
        for (unsigned int i = 0; i < vmem_m.size(); ++i) {
            vmMin[i] = vmMax[i] = vmAvg[i] = vmem_m[i];
        }
        return;
    }

    new_reduce(vmem_m.data(), vmAvg.data(), vmem_m.size(), std::plus<double>());

    double inodes = 1.0 / double(Ippl::getNodes());
    for (auto& vm : vmAvg) {
        vm *= inodes;
    }

    new_reduce(vmem_m.data(), vmMin.data(), vmem_m.size(), std::less<double>());
    new_reduce(vmem_m.data(), vmMax.data(), vmem_m.size(), std::greater<double>());
}


void MemoryProfiler::write(PartBunchBase<double, 3> *beam) {

    this->update_m();

    vm_t vmMin(vmem_m.size());
    vm_t vmMax(vmem_m.size());
    vm_t vmAvg(vmem_m.size());

    this->compute_m(vmMin, vmMax, vmAvg);

    if ( Ippl::myNode() != 0 ) {
        return;
    }

    double  pathLength = 0.0;
    if (OpalData::getInstance()->isInOPALCyclMode())
        pathLength = beam->getLPath();
    else
        pathLength = beam->get_sPos();

    header_m();

    this->open();

    this->writeHeader();

    SDDSDataRow row(this->getColumnNames());
    row.addColumn("t", beam->getT() * 1e9);             // 1
    row.addColumn("s", pathLength);                     // 2
    
    row.addColumn("VmPeak-Min", vmMin[VMPEAK]);
    row.addColumn("VmPeak-Max", vmMax[VMPEAK]);
    row.addColumn("VmPeak-Avg", vmAvg[VMPEAK]);

    row.addColumn("VmSize-Min", vmMin[VMSIZE]);
    row.addColumn("VmSize-Max", vmMax[VMSIZE]);
    row.addColumn("VmSize-Avg", vmAvg[VMSIZE]);

    row.addColumn("VmHWM-Min", vmMin[VMHWM]);
    row.addColumn("VmHWM-Max", vmMax[VMHWM]);
    row.addColumn("VmHWM-Avg", vmAvg[VMHWM]);

    row.addColumn("VmRSS-Min", vmMin[VMRSS]);
    row.addColumn("VmRSS-Max", vmMax[VMRSS]);
    row.addColumn("VmRSS-Avg", vmAvg[VMRSS]);

    row.addColumn("VmStk-Min", vmMin[VMSTK]);
    row.addColumn("VmStk-Max", vmMax[VMSTK]);
    row.addColumn("VmStk-Avg", vmAvg[VMSTK]);

    this->writeRow(row);

    this->newline();

    this->close();
}
