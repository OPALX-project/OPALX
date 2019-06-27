#include "MemoryProfiler.h"

#ifdef __linux__
#include <sys/types.h>
#include <unistd.h>
#endif

#include "Utilities/Timer.h"
#include "Utilities/OpalException.h"

#include <boost/filesystem.hpp>

#include "Ippl.h"

#include <sstream>

namespace {
    std::string to_string(long double d) {
        std::ostringstream ss;
        ss.precision(15);
        ss << d;
        return ss.str();
    }
}

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


void MemoryProfiler::header_m() {

    static bool isNotFirst = false;
    if ( isNotFirst ) {
        return;
    }
    isNotFirst = true;

    columns_m.addColumn("t", "double", "ns", "Time");

    columns_m.addColumn("s", "double", "m", "Path length");

    // peak virtual memory size
    columns_m.addColumn("VmPeak-Min", "double", unit_m[VirtualMemory::VMPEAK],
                        "Minimum peak virtual memory size");

    columns_m.addColumn("VmPeak-Max", "double", unit_m[VirtualMemory::VMPEAK],
                        "Maximum peak virtual memory size");

    columns_m.addColumn("VmPeak-Avg", "double", unit_m[VirtualMemory::VMPEAK],
                        "Average peak virtual memory size");

    // virtual memory size
    columns_m.addColumn("VmSize-Min", "double", unit_m[VirtualMemory::VMSIZE],
                        "Minimum virtual memory size");

    columns_m.addColumn("VmSize-Max", "double", unit_m[VirtualMemory::VMSIZE],
                        "Maximum virtual memory size");

    columns_m.addColumn("VmSize-Avg", "double", unit_m[VirtualMemory::VMSIZE],
                        "Average virtual memory size");

    // peak resident set size ("high water mark")
    columns_m.addColumn("VmHWM-Min", "double", unit_m[VirtualMemory::VMHWM],
                        "Minimum peak resident set size");

    columns_m.addColumn("VmHWM-Max", "double", unit_m[VirtualMemory::VMHWM],
                        "Maximum peak resident set size");

    columns_m.addColumn("VmHWM-Avg", "double", unit_m[VirtualMemory::VMHWM],
                        "Average peak resident set size");

    // resident set size
    columns_m.addColumn("VmRSS-Min", "double", unit_m[VirtualMemory::VMRSS],
                        "Minimum resident set size");

    columns_m.addColumn("VmRSS-Max", "double", unit_m[VirtualMemory::VMRSS],
                        "Maximum resident set size");

    columns_m.addColumn("VmRSS-Avg", "double", unit_m[VirtualMemory::VMRSS],
                        "Average resident set size");

    // stack size
    columns_m.addColumn("VmStk-Min", "double", unit_m[VirtualMemory::VMSTK],
                        "Minimum stack size");

    columns_m.addColumn("VmStk-Max", "double", unit_m[VirtualMemory::VMSTK],
                        "Maximum stack size");

    columns_m.addColumn("VmStk-Avg", "double", unit_m[VirtualMemory::VMSTK],
                        "Average stack size");

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

    columns_m.addColumnValue("t", beam->getT() * 1e9);             // 1
    columns_m.addColumnValue("s", pathLength);                     // 2

    // boost::variant can't overload double and long double. By using a
    // string this shortcoming can be bypassed.
    columns_m.addColumnValue("VmPeak-Min", ::to_string(vmMin[VMPEAK]));
    columns_m.addColumnValue("VmPeak-Max", ::to_string(vmMax[VMPEAK]));
    columns_m.addColumnValue("VmPeak-Avg", ::to_string(vmAvg[VMPEAK]));

    columns_m.addColumnValue("VmSize-Min", ::to_string(vmMin[VMSIZE]));
    columns_m.addColumnValue("VmSize-Max", ::to_string(vmMax[VMSIZE]));
    columns_m.addColumnValue("VmSize-Avg", ::to_string(vmAvg[VMSIZE]));

    columns_m.addColumnValue("VmHWM-Min", ::to_string(vmMin[VMHWM]));
    columns_m.addColumnValue("VmHWM-Max", ::to_string(vmMax[VMHWM]));
    columns_m.addColumnValue("VmHWM-Avg", ::to_string(vmAvg[VMHWM]));

    columns_m.addColumnValue("VmRSS-Min", ::to_string(vmMin[VMRSS]));
    columns_m.addColumnValue("VmRSS-Max", ::to_string(vmMax[VMRSS]));
    columns_m.addColumnValue("VmRSS-Avg", ::to_string(vmAvg[VMRSS]));

    columns_m.addColumnValue("VmStk-Min", ::to_string(vmMin[VMSTK]));
    columns_m.addColumnValue("VmStk-Max", ::to_string(vmMax[VMSTK]));
    columns_m.addColumnValue("VmStk-Avg", ::to_string(vmAvg[VMSTK]));

    this->writeRow();

    this->newline();

    this->close();
}