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
#include "Utility/IpplTimings.h"
#include "Utility/Inform.h"
#include "Utility/IpplInfo.h"
#include "Message/GlobalComm.h"
#include "PETE/IpplExpressions.h"

#include <boost/algorithm/string/predicate.hpp>

#include <fstream>
#include <iostream>
#include <algorithm>
// static data members of IpplTimings class

Timing* IpplTimings::instance = new Timing();
std::stack<Timing*> IpplTimings::stashedInstance;

//////////////////////////////////////////////////////////////////////
// default constructor
Timing::Timing():
    TimerList(),
    TimerMap()
{ }


//////////////////////////////////////////////////////////////////////
// destructor
Timing::~Timing() {
    for (TimerMap_t::iterator it = TimerMap.begin(); it != TimerMap.end(); ++ it) {
        it->second = 0;
    }
    TimerMap.clear();

    TimerList.clear();
}


//////////////////////////////////////////////////////////////////////
// create a timer, or get one that already exists
Timing::TimerRef Timing::getTimer(const char *nm) {
    std::string s(nm);
    TimerInfo *tptr = 0;
    TimerMap_t::iterator loc = TimerMap.find(s);
    if (loc == TimerMap.end()) {
        tptr = new TimerInfo;
        tptr->indx = TimerList.size();
        tptr->name = s;
        TimerMap.insert(TimerMap_t::value_type(s,tptr));
        TimerList.push_back(my_auto_ptr<TimerInfo>(tptr));
    } else {
        tptr = (*loc).second;
    }
    return tptr->indx;
}


//////////////////////////////////////////////////////////////////////
// start a timer
void Timing::startTimer(TimerRef t) {
    if (t >= TimerList.size())
        return;
    TimerList[t]->start();
}


//////////////////////////////////////////////////////////////////////
// stop a timer, and accumulate it's values
void Timing::stopTimer(TimerRef t) {
    if (t >= TimerList.size())
        return;
    TimerList[t]->stop();
}


//////////////////////////////////////////////////////////////////////
// clear a timer, by turning it off and throwing away its time
void Timing::clearTimer(TimerRef t) {
    if (t >= TimerList.size())
        return;
    TimerList[t]->clear();
}

#ifdef IPPL_XT3
//////////////////////////////////////////////////////////////////////
// print out the timing results
void Timing::print() {
    int i,j;
    if (TimerList.size() < 1)
	return;

    // report the average time for each timer
    Inform msg("Timings");
    msg << level1 << "-----------------------------------------------------------------\n";
    msg << "     Timing (dclock) results for " << Ippl::getNodes() << " nodes:" << "\n";
    msg << "-----------------------------------------------------------------";
    msg << "\n";
    for (i=0; i<1; ++i){
	TimerInfo *tptr = TimerList[i].get();
	double walltotal = 0.0;
	reduce(tptr->wallTime, walltotal, OpMaxAssign());
	msg << tptr->name.c_str() << " ";
	for (j=strlen(tptr->name.c_str()); j < 10; ++j)
	    msg << ".";
	msg << " Wall tot = ";
	msg.width(10);
	msg << walltotal;
	msg << "\n\n";
    }

    for (i=1; i < TimerList.size(); ++i) {
	TimerInfo *tptr = TimerList[i].get();
	double wallmax = 0.0, wallmin = 0.0;
	double wallavg = 0.0 ;
	reduce(tptr->wallTime, wallmax, OpMaxAssign());
	reduce(tptr->wallTime, wallmin, OpMinAssign());
	reduce(tptr->wallTime, wallavg, OpAddAssign());
	msg << tptr->name.c_str() << " ";
	for (j=strlen(tptr->name.c_str()); j < 10; ++j)
	    msg << ".";
	msg << " Wall max = ";
	msg.width(10);
	msg << wallmax << "\n";
	for (j = 0; j < 21; ++j)
	    msg << " ";
	msg << " Wall avg = ";
	msg.width(10);
	msg << wallavg / Ippl::getNodes() << "\n";
	for (j = 0; j < 21; ++j)
	    msg << " ";
	msg << " Wall min = ";
	msg.width(10);
	msg << wallmin << "\n\n";
    }
    msg << "-----------------------------------------------------------------";
    msg << endl;
}

#else

//////////////////////////////////////////////////////////////////////
// print out the timing results
void Timing::print() {
    if (TimerList.size() < 1)
        return;

    // report the average time for each timer
    Inform msg("Timings");
    msg << level1
        << "-----------------------------------------------------------------";
    msg << "\n";
    msg << "     Timing results for " << Ippl::getNodes() << " nodes:" << "\n";
    msg << "-----------------------------------------------------------------";
    msg << "\n";

    {
        TimerInfo *tptr = TimerList[0].get();
        double walltotal = 0.0, cputotal = 0.0;
        reduce(tptr->wallTime, walltotal, OpMaxAssign());
        reduce(tptr->cpuTime, cputotal, OpMaxAssign());
        size_t lengthName = std::min(tptr->name.length(), 19lu);
        msg << tptr->name.substr(0,lengthName)
            << std::string().assign(20 - lengthName,'.')
            << " Wall tot = " << std::setw(10) << walltotal << ","
            << " CPU tot = " << std::setw(10) << cputotal << "\n"
            << "\n";
    }

    auto begin = ++ TimerList.begin();
    auto end = TimerList.end();
    std::sort(begin, end, [](const my_auto_ptr<TimerInfo>& a, const my_auto_ptr<TimerInfo>& b)
              {
                  return boost::ilexicographical_compare(a->name, b->name);
              });

    for (unsigned int i=1; i < TimerList.size(); ++i) {
        TimerInfo *tptr = TimerList[i].get();
        double wallmax = 0.0, cpumax = 0.0, wallmin = 0.0, cpumin = 0.0;
        double wallavg = 0.0, cpuavg = 0.0;
        reduce(tptr->wallTime, wallmax, OpMaxAssign());
        reduce(tptr->cpuTime,  cpumax,  OpMaxAssign());
        reduce(tptr->wallTime, wallmin, OpMinAssign());
        reduce(tptr->cpuTime,  cpumin,  OpMinAssign());
        reduce(tptr->wallTime, wallavg, OpAddAssign());
        reduce(tptr->cpuTime,  cpuavg,  OpAddAssign());
        size_t lengthName = std::min(tptr->name.length(), 19lu);

        msg << tptr->name.substr(0,lengthName)
            << std::string().assign(20 - lengthName, '.')
            << " Wall max = " << std::setw(10) << wallmax << ","
            << " CPU max = " << std::setw(10) << cpumax << "\n"
            << std::string().assign(20,' ')
            << " Wall avg = " << std::setw(10) << wallavg / Ippl::getNodes() << ","
            << " CPU avg = " << std::setw(10) << cpuavg / Ippl::getNodes() << "\n"
            << std::string().assign(20,' ')
            << " Wall min = " << std::setw(10) << wallmin << ","
            << " CPU min = " << std::setw(10) << cpumin << "\n"
            << "\n";
    }
    msg << "-----------------------------------------------------------------";
    msg << endl;
}
#endif

//////////////////////////////////////////////////////////////////////
// save the timing results into a file
void Timing::print(const std::string &fn) {

    std::ofstream *timer_stream;
    Inform *msg;

    if (TimerList.size() < 1)
        return;

    timer_stream = new std::ofstream;
    timer_stream->open( fn.c_str(), std::ios::out );
    msg = new Inform( 0, *timer_stream, 0 );
    Inform gmsg("bla ");
    // report the average time for each timer
    // Inform msg("Timings");
    /*
     *msg << "---------------------------------------------------------------------------";
     *msg << "\n";
     *msg << "     Timing results for " << Ippl::getNodes() << " nodes:" << "\n";
     *msg << "---------------------------------------------------------------------------";
     *msg << " name nodes (cputot cpumax) (walltot wallmax) cpumin wallmin cpuav wallav  ";
     *msg << "\n";
     */

    *msg << std::setw(27) << "num Nodes"
         << std::setw(10) << "CPU tot"
         << std::setw(11) << "Wall tot\n"
         << std::string().assign(47,'=')
         << "\n";
    {
        TimerInfo *tptr = TimerList[0].get();
        double walltotal = 0.0, cputotal = 0.0;
        reduce(tptr->wallTime, walltotal, OpMaxAssign());
        reduce(tptr->cpuTime, cputotal, OpMaxAssign());
        size_t lengthName = std::min(tptr->name.length(), 19lu);
        *msg << tptr->name.substr(0,lengthName);
        for (int j=lengthName; j < 20; ++j) {
            *msg << ".";
        }
        *msg  << " " << std::setw(6)  << Ippl::getNodes()
              << " " << std::setw(9) << std::setprecision(4) << cputotal
              << " " << std::setw(9) << std::setprecision(4) << walltotal
              << "\n";
    }

    auto begin = ++ TimerList.begin();
    auto end = TimerList.end();
    std::sort(begin, end, [](const my_auto_ptr<TimerInfo>& a, const my_auto_ptr<TimerInfo>& b)
              {
                  return boost::ilexicographical_compare(a->name, b->name);
              });

    *msg << "\n"
         << std::setw(27) << "num Nodes"
         << std::setw(10) << "CPU max"
         << std::setw(10) << "Wall max"
         << std::setw(10) << "CPU min"
         << std::setw(10) << "Wall min"
         << std::setw(10) << "CPU avg"
         << std::setw(11) << "Wall avg\n"
         << std::string().assign(87,'=')
         << "\n";
    for (unsigned int i=0; i < TimerList.size(); ++i) {
        TimerInfo *tptr = TimerList[i].get();
        double wallmax = 0.0, cpumax = 0.0, wallmin = 0.0, cpumin = 0.0;
        double wallavg = 0.0, cpuavg = 0.0;
        reduce(tptr->wallTime, wallmax, OpMaxAssign());
        reduce(tptr->cpuTime,  cpumax,  OpMaxAssign());
        reduce(tptr->wallTime, wallmin, OpMinAssign());
        reduce(tptr->cpuTime,  cpumin,  OpMinAssign());
        reduce(tptr->wallTime, wallavg, OpAddAssign());
        reduce(tptr->cpuTime,  cpuavg,  OpAddAssign());
        size_t lengthName = std::min(tptr->name.length(), 19lu);
        *msg << tptr->name.substr(0,lengthName);
        for (int j=lengthName; j < 20; ++j) {
            *msg << ".";
        }
        *msg << " " << std::setw(6) << Ippl::getNodes()
             << " " << std::setw(9) << std::setprecision(4) << cpumax
             << " " << std::setw(9) << std::setprecision(4) << wallmax
             << " " << std::setw(9) << std::setprecision(4) << cpumin
             << " " << std::setw(9) << std::setprecision(4) << wallmin
             << " " << std::setw(9) << std::setprecision(4) << cpuavg / Ippl::getNodes()
             << " " << std::setw(9) << std::setprecision(4) << wallavg / Ippl::getNodes()
             << endl;
    }
    timer_stream->close();
    delete msg;
    delete timer_stream;
}

IpplTimings::IpplTimings() { }
IpplTimings::~IpplTimings() { }

void IpplTimings::stash() {
    PAssert_EQ(stashedInstance.size(), 0);

    stashedInstance.push(instance);
    instance = new Timing();
}

void IpplTimings::pop() {
    PAssert_GT(stashedInstance.size(), 0);

    delete instance;
    instance = stashedInstance.top();
    stashedInstance.pop();
}

/***************************************************************************
 * $RCSfile: IpplTimings.cpp,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:33 $
 ***************************************************************************/


/***************************************************************************
 * $RCSfile: addheaderfooter,v $   $Author: adelmann $
 * $Revision: 1.1.1.1 $   $Date: 2003/01/23 07:40:17 $
 * IPPL_VERSION_ID: $Id: addheaderfooter,v 1.1.1.1 2003/01/23 07:40:17 adelmann Exp $
 ***************************************************************************/
