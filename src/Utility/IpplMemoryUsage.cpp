// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 * This program was prepared by PSI. 
 * All rights in the program are reserved by PSI.
 * Neither PSI nor the author(s)
 * makes any warranty, express or implied, or assumes any liability or
 * responsibility for the use of this software
 *
 * Visit www.amas.web.psi for more details
 *
 ***************************************************************************/

// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

// include files
#include "Utility/IpplMemoryUsage.h"
#include "Utility/Inform.h"
#include "Utility/IpplInfo.h"
#include "Message/GlobalComm.h"
#include "PETE/IpplExpressions.h"

#include <fstream>
#include <iostream>
#include <string>
#include <unistd.h>
#include <ios>
#include <iostream>
#include <fstream>

// static data members of IpplMemoryUsage class
IpplMemoryUsage::MemoryList_t IpplMemoryUsage::MemoryList;
IpplMemoryUsage::MemoryMap_t  IpplMemoryUsage::MemoryMap;

void IpplMemoryInfo::doSample() { 
	
    //////////////////////////////////////////////////////////////////////////////
    //
    // process_mem_usage(double &, double &) - takes two doubles by reference,
    // attempts to read the system-dependent data for a process' virtual memory
    // size and resident set size, and return the results in KB.
    //
    // On failure, returns 0.0, 0.0
    /* from the man pages
      vsize %lu
                     Virtual memory size in bytes.

      rss %ld
                     Resident  Set Size: number of pages the process has in real memory, minus
                     3 for administrative purposes. This is just the pages which count towards
                     text,  data,  or stack space.  This does not include pages which have not
                     been demand-loaded in, or which are swapped out.




    */
    
    // 'file' stat seems to give the most reliable results
    std::ifstream stat_stream("/proc/self/stat", std::ios_base::in);
    
    stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
		>> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
		>> utime >> stime >> cutime >> cstime >> priority >> nice
		>> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest
    
    page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    
    vm_usage     = vsize / 1024.0 / 1024.0;
    resident_set = rss * page_size_kb / 1024.0;
    
    reduce(vm_usage,vm_usage,OpAddAssign());
    reduce(resident_set,resident_set,OpAddAssign());
    
    maxM = 0;
    minM = 0;
    reduce(vm_usage,maxM,OpMaxAssign());
    reduce(vm_usage,minM,OpMinAssign());
}


//////////////////////////////////////////////////////////////////////
// default constructor
IpplMemoryUsage::IpplMemoryUsage() { }


//////////////////////////////////////////////////////////////////////
// destructor
IpplMemoryUsage::~IpplMemoryUsage() { }


//////////////////////////////////////////////////////////////////////
// create a timer, or get one that already exists
IpplMemoryUsage::MemRef IpplMemoryUsage::getMemObserver(const char *nm, double memPerCore) {
    std::string s(nm);
    MemoryInfo *tptr = 0;
    MemoryMap_t::iterator loc = MemoryMap.find(s);
    if (loc == MemoryMap.end()) {
	tptr = new MemoryInfo;
	tptr->indx = MemoryList.size();
	tptr->name = s;
	tptr->memPerCore = memPerCore;
	MemoryMap.insert(MemoryMap_t::value_type(s,tptr));
	MemoryList.push_back(my_auto_ptr<MemoryInfo>(tptr));
    } else {
	tptr = (*loc).second;
    }
    return tptr->indx;
}


//////////////////////////////////////////////////////////////////////
// sample memory
void IpplMemoryUsage::sample(MemRef t, const char *nm) {
    if (t < 0 || (unsigned int) t >= MemoryList.size())
	return;
    MemoryList[t]->doSample();
}

/////////////////////////////////////////////////////////////////////
// print out the timing results

void IpplMemoryUsage::doPrint(Inform &msg) {
    if (MemoryList.size() < 1)
	return;
    
    // report the average time for each timer

    {
        MemoryInfo *tptr = MemoryList[0].get();

        msg << "------------------------------------------------------------------------";
        msg << endl;
        msg << "     Memory Statistics results for " << Ippl::getNodes() << " node(s):";
        msg << " page size (kB) = " << tptr->page_size_kb << endl;
        msg << "     Total memory (GB)= " << tptr->memPerCore*Ippl::getNodes() << endl;
        msg << "------------------------------------------------------------------------";
        msg << endl;

        msg << tptr->name.c_str() << " ";
        for (int j=strlen(tptr->name.c_str()); j < 13; ++j)
            msg << ".";
        msg << " vm = ";
        msg.width(10);
        msg << tptr->vm_usage*0.001 << " (GB)\tdeltavm (GB)";
        msg << endl << endl;
    }

    for (unsigned int i=1; i < MemoryList.size(); ++i) {
	MemoryInfo *tptrprev = MemoryList[i-1].get();
	MemoryInfo *tptr = MemoryList[i].get();
	
	msg << tptr->name.c_str() << " ";
	for (int j=strlen(tptr->name.c_str()); j < 13; ++j)
	    msg << ".";
	msg << " vm = ";
	msg.width(10);
	msg << tptr->vm_usage*0.001 << " (GB)\t" << (tptr->vm_usage - tptrprev->vm_usage)*0.001 << endl;
	msg << endl;
    }
    msg << "------------------------------------------------------------------------";
    msg << endl;
}

void IpplMemoryUsage::print() {

    Inform *msg = new Inform("Memstat");
    doPrint(*msg);
    delete msg;
}

//////////////////////////////////////////////////////////////////////
// save the timing results into a file
void IpplMemoryUsage::print(std::string fn) {
  std::ofstream *memstat_stream;
  Inform *msg;
  memstat_stream = new std::ofstream;
  memstat_stream->open( fn.c_str(), std::ios::out );
  msg = new Inform( 0, *memstat_stream, 0 );
  doPrint(*msg);  
  memstat_stream->close();
  delete msg;
  delete memstat_stream;
}


/***************************************************************************
 * $RCSfile: IpplMemoryUsage.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 ***************************************************************************/


/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

