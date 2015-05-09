// -*- C++ -*-
/***************************************************************************
 *
 * The IPPL Framework
 * 
 *
 * Visit http://people.web.psi.ch/adelmann/ for more details
 *
 ***************************************************************************/

#ifndef IPPL_MEMPRYUSAGE_H
#define IPPL_MEMPRYUSAGE_H

/*************************************************************************
 * IpplMemoryUsage - a simple singleton class which lets the user watch the
 *   memory consumption of a program and print statistics at the end of the program.
 *
 * General usage
 *  1) create a (named) memory observer object and record the current memory status
 *     IpplMemoryUsage::MemoryRef val = IpplMemoryUsage::getMemObserver("name");
 *     This will either create a new one, or return a ref to an existing one
 *
 *  2) sample the memory consumtion:
 *     IpplMemoryUsage::sample(val,"name");
 *
 *  3) print out the results:
 *     IpplMemoryUsage::print();
 *
 *************************************************************************/

// include files
#include "Utility/my_auto_ptr.h"
#include "Ippl.h"
#ifdef IPPL_STDSTL
#include <vector>
#include <map>
#else
#include <vector.h>
#include <map.h>
#endif // IPPL_STDSTL

// a simple class used to store memory values
class IpplMemoryInfo
{
public:
    // typedef for reference to memory information
    typedef int MemRef;

    // constructor
    IpplMemoryInfo() : name(""), indx(-1) {
   
    }

    // destructor
    ~IpplMemoryInfo() { }

    void doSample(); 

    // the name of this memory observer
    std::string name;

    // memory per core in GB
    double memPerCore;

    // dummy vars for leading entries in stat that we don't care about
    //
    std::string pid, comm, state, ppid, pgrp, session, tty_nr;
    std::string tpgid, flags, minflt, cminflt, majflt, cmajflt;
    std::string utime, stime, cutime, cstime, priority, nice;
    std::string O, itrealvalue, starttime;
    
    // the two fields we want
    unsigned long vsize;
    long rss;

    long page_size_kb;

    double vm_usage;
    double resident_set;

    double maxM;
    double minM;

    // an index value for this memory observer
    MemRef indx;
};


class IpplMemoryUsage
{
public:
    // typedef for reference to memory information
    typedef int MemRef;
    
    // a typedef for the memory information object
    typedef IpplMemoryInfo MemoryInfo;

public:
    // Default constructor
    IpplMemoryUsage();
    
    // Destructor - clear out the existing memory info
    ~IpplMemoryUsage();
 

    // create a initial memory footprint, or get one that already exists
    static MemRef getMemObserver(const char *, double memPerCore);

    // sample the current memory consumption
    static void sample(MemRef, const char *);

    // return a MemoryInfo struct by asking for the name
    static MemoryInfo *infoMemory(const char *nm) {
	return MemoryMap[std::string(nm)];
    }

    //
    // I/O methods
    //

    // print the results to standard out
    static void print();

    // print the results to a file
    static void print(std::string fn);


private:

    static void doPrint(Inform &msg);

    // type of storage for list of MemoryInfo
    typedef std::vector<my_auto_ptr<MemoryInfo> > MemoryList_t;
    typedef std::map<std::string, MemoryInfo *> MemoryMap_t;

    // a list of memory info structs
    static MemoryList_t MemoryList;

    // a map of memory infos, keyed by string
    static MemoryMap_t MemoryMap;
};

#endif

/***************************************************************************
 * $RCSfile: IpplMemoryUsage.h,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 ***************************************************************************/

/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $ 
 ***************************************************************************/

