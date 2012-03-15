// ------------------------------------------------------------------------
// $RCSfile: Main.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.9.2.2 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Main program for OPAL
//
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 20:10:10 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#ifdef DEBUG_ALLOCATOR
#include "Allocator.h"
#endif

#include "Ippl.h"

#include "AbstractObjects/DoomDB.h"
#include "AbstractObjects/OpalData.h"
#include "OpalConfigure/Configure.h"
#include "FixedAlgebra/FTps.h"
#include "OpalParser/OpalParser.h"
#include "Parser/FileStream.h"
#include "Parser/TerminalStream.h"
#include "Utilities/ClassicException.h"
#include "Utilities/ParseError.h"
#include "Utilities/Timer.h"

#include <iostream>
#include <vector>

#include <new>
#include <exception>

#include "config.h"

using namespace std;
using namespace OPALTimer;
//  DTA
#define NC 5

void printStringVector(const vector<string> &strings) {
    unsigned int iend = strings.size(), nc = 0;
    for(unsigned int i = 0; i < iend; ++i) {
        std::cout << "  " << strings[i];
        if(++nc == NC) { nc = 0; std::cout << std::endl; }
    }
    if(nc != 0) std::cout << std::endl;
    std::cout << std::endl;
}
// /DTA

// Global data.
// ------------------------------------------------------------------------

//: The global OPAL data structure.
//  Contains mainly the directory of known objects,
//  but also the directories used to maintain tables and expressions
//  up-to-date, as well as the run title.
OpalData OPAL;

//: The Doom data base.
DoomDB DOOM_DB;

//: A global Inform object

Inform *gmsg;
Inform *gmsg2all;

//: The OPAL main program.
int main(int argc, char *argv[]) {
    Ippl ippl(argc, argv);

    OPALTimer::Timer simtimer;

    string dateStr(simtimer.date());
    string timeStr(simtimer.time());

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("mainTimer");
    IpplTimings::startTimer(mainTimer);

    gmsg = new  Inform("OPAL ");
    Inform hmsg("");
    gmsg2all = new  Inform("OPAL", INFORM_ALL_NODES);
    string mySpace("            ");

    if(Ippl::myNode() == 0) remove("errormsg.txt");

    /*

    hmsg << mySpace << string("   _______  _______  _______  _  ")<< endl;
    hmsg << mySpace << string("  (  ___  )(  ____ )(  ___  )( ")<< endl;
    hmsg << mySpace << string("  | (   ) || (    )|| (   ) || (")<< endl;
    hmsg << mySpace << string("  | |   | || (____)|| (___) || |")<< endl;
    hmsg << mySpace << string("  | |   | ||  _____)|  ___  || |")<< endl;
    hmsg << mySpace << string("  | |   | || (      | (   ) || |")<< endl;
    //  hmsg << mySpace << string("  | (___) || )      | )   ( || (____/\")<< endl;
    //
    hmsg << mySpace << string("  (_______)|//       |//     /\|(_______// ")<< endl;
     */

    hmsg << mySpace <<  "   ____  _____       ___ " << endl;
    hmsg << mySpace <<  "  / __ \\|  __ \\ /\\   | | " << endl;
    hmsg << mySpace <<  " | |  | | |__) /  \\  | |" << endl;
    hmsg << mySpace <<  " | |  | |  ___/ /\\ \\ | |" << endl ;
    hmsg << mySpace <<  " | |__| | |  / ____ \\| |____" << endl;
    hmsg << mySpace <<  "  \\____/|_| /_/    \\_\\______|" << endl;


    *gmsg << endl << "This is OPAL (Object Oriented Parallel Accelerator Library) Version " << PACKAGE_VERSION << " SVN version " << SVN_VERSION << "  (c) PSI, http://amas.web.psi.ch" << endl
          << endl;
    *gmsg << "Please send cookies, goodies or other motivations (wine and beer ... ) to " << PACKAGE_BUGREPORT << endl;
    *gmsg << "Time: " << timeStr << " date: " << dateStr << endl;

    const OpalParser parser;

    //  DTA
    std::cout.precision(16);
    std::cout.setf(std::ios::scientific, std::ios::floatfield);
    std::cout.setf(std::ios::showpos);
    std::cerr.precision(16);
    std::cerr.setf(std::ios::scientific, std::ios::floatfield);
    std::cerr.setf(std::ios::showpos);
    // /DTA

    // Set global truncation orders.
    FTps<double, 2>::setGlobalTruncOrder(20);
    FTps<double, 4>::setGlobalTruncOrder(15);
    FTps<double, 6>::setGlobalTruncOrder(10);

    try {
        Configure::configure();
        int arg = 1;

        if(arg < argc  &&  strcmp(argv[arg], "-db") == 0) {
            if(argc > 2) {
                DOOM_DB.firstOpen(argv[arg+1]);
                arg += 2;
            } else {
                *gmsg << "Data base name missing." << endl;
                exit(1);
            }
        }

        // Read startup file.
        FileStream::setEcho(true);

#ifndef __LIBCATAMOUNT__
        string startup = getenv("HOME");
        startup += "/init.opal";
        FileStream::setEcho(false);
        FileStream *is;

        try {
            is = new FileStream(startup);
        } catch(...) {
            is = 0;
            *gmsg << "No startup file \"" << startup << "\" found." << endl;
        }

        if(is) {
            *gmsg << "Reading startup file \"" << startup << "\"." << endl;
            parser.run(is);
            *gmsg << "Finished reading startup file." << endl;
        }
#else
        *gmsg << "On Cray can not read startup file" << endl;
#endif

        if(argc > 1)
            OPAL.storeInputFn(string(argv[1]));

        if(argc > 3) {
            if(argc > 5) {
                // will write dumping date into a new h5 file
                for(int ii = 2; ii < 6; ii = ii + 2)
                    // The sequence of the two arguments is free
                    if(string(argv[ii]) == string("-restart")) {
                        OPAL.setRestartRun();
                        OPAL.setRestartStep(atoi(argv[ii+1]));
                    } else if((string(argv[ii]) == string("-restartfn"))) {
                        OPAL.setRestartFileName(argv[ii+1]);
                    }
            } else {
                // will append dumping date into old h5 file
                if(string(argv[2]) == string("-restart")) {
                    OPAL.setRestartRun();
                    OPAL.setRestartStep(atoi(argv[3]));
                }
            }
        }

        if(arg >= argc) {
            // Run commands from standard input
            parser.run(new TerminalStream("OPAL"));
        } else {
            FileStream *is;

            try {
                is = new FileStream(argv[arg]);
            } catch(...) {
                is = 0;
                *gmsg << "Input file \"" << argv[arg] << "\" not found." << endl;
            }

            if(is) {
                *gmsg << "Reading input stream \"" << argv[arg] << "\"." << endl;
                parser.run(is);
                *gmsg << "End of input stream \"" << argv[arg] << "\"." << endl;
            }
        }

        ////  DTA
        //  std::cout << std::endl << std::endl;
        //  vector<string> names=OPAL.getAllNames();
        //  printStringVector(names);
        //// /DTA

        DOOM_DB.shut();
        IpplTimings::stopTimer(mainTimer);
        IpplTimings::print();
        IpplTimings::print(string("timing.dat"));

        if(Ippl::myNode() == 0) {
            ifstream errormsg("errormsg.txt");
            if(errormsg.good()) {
                char buffer[256];
                string closure("*                                                                                  *\n");
                *gmsg << "\n"
                      << "* **********************************************************************************\n"
                      << "* ************** E R R O R * * M E S S A G E S *************************************\n"
                      << "* **********************************************************************************"
                      << endl;
                errormsg.getline(buffer, 256);
                while(errormsg.good()) {
                    if(errormsg.gcount() == 1) {
                        *gmsg << closure;
                    } else {
                        *gmsg << buffer << closure.substr(errormsg.gcount() - 1);
                    }
                    errormsg.getline(buffer, 256);
                }
                *gmsg << closure
                      << "* **********************************************************************************\n"
                      << "* **********************************************************************************"
                      << endl;
            }
            errormsg.close();
        }

        delete gmsg;
        delete gmsg2all;
        exit(0);
        //   return 0;
    } catch(ClassicException &ex) {
        *gmsg << endl << "*** User error detected by function \""
              << ex.where() << "\"" << endl;
        abort();
    } catch(bad_alloc) {
        *gmsg << "Sorry, virtual memory exhausted." << endl;
        abort();
    } catch(...) {
        *gmsg << "Unexpected error." << endl;
        abort();
    }
}
