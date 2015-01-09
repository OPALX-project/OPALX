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

//~ #ifdef DEBUG_ALLOCATOR
//~ #include "Allocator.h"
//~ #endif

#include "opal.h"

Ippl *ippl;
Inform *gmsg;

#include "H5hut.h"

#include "AbstractObjects/OpalData.h"
#include "OpalConfigure/Configure.h"
#include "OpalParser/OpalParser.h"
#include "Parser/FileStream.h"
#include "Parser/TerminalStream.h"
#include "Utilities/Timer.h"
#include "Fields/Fieldmap.hh"
#include "FixedAlgebra/FTps.h"


#include "config.h"

#ifdef HAVE_AMR_SOLVER
#include <ParallelDescriptor.H>
#endif

//  DTA
#define NC 5
#define MY_MASK 0755

void printStringVector(const std::vector<std::string> &strings) {
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
//OpalData OPAL;

//: A global Inform object

//: The OPAL main program.
int main(int argc, char *argv[]) {
    ippl = new Ippl(argc, argv);
    gmsg = new  Inform("OPAL");

#ifdef HAVE_AMR_SOLVER
    *gmsg << "Initializing BoxLib with inputs file " << argv[1] << endl;
    BoxLib::Initialize(argc, argv, true, Ippl::getComm());
    *gmsg << "Done initializing BoxLib with inputs file " << argv[1] << endl;
#endif

    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());

    H5SetVerbosityLevel(0); //65535);

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("mainTimer");
    IpplTimings::startTimer(mainTimer);

    Inform hmsg("");
    std::string mySpace("            ");

    if(Ippl::myNode() == 0) remove("errormsg.txt");

    /*

       hmsg << mySpace << std::string("   _______  _______  _______  _  ")<< endl;
       hmsg << mySpace << std::string("  (  ___  )(  ____ )(  ___  )( ")<< endl;
       hmsg << mySpace << std::string("  | (   ) || (    )|| (   ) || (")<< endl;
       hmsg << mySpace << std::string("  | |   | || (____)|| (___) || |")<< endl;
       hmsg << mySpace << std::string("  | |   | ||  _____)|  ___  || |")<< endl;
       hmsg << mySpace << std::string("  | |   | || (      | (   ) || |")<< endl;
    //  hmsg << mySpace << std::string("  | (___) || )      | )   ( || (____/\")<< endl;
    //
    hmsg << mySpace << std::string("  (_______)|//       |//     /\|(_______// ")<< endl;
    */

    hmsg << mySpace <<  "   ____  _____       ___ " << endl;
    hmsg << mySpace <<  "  / __ \\|  __ \\ /\\   | | " << endl;
    hmsg << mySpace <<  " | |  | | |__) /  \\  | |" << endl;
    hmsg << mySpace <<  " | |  | |  ___/ /\\ \\ | |" << endl ;
    hmsg << mySpace <<  " | |__| | |  / ____ \\| |____" << endl;
    hmsg << mySpace <<  "  \\____/|_| /_/    \\_\\______|" << endl;


    *gmsg << endl << "This is OPAL (Object Oriented Parallel Accelerator Library) Version " << PACKAGE_VERSION << "  (c) PSI, http://amas.web.psi.ch" << endl
          << endl;
    *gmsg << "Please send cookies, goodies or other motivations (wine and beer ... ) to the OPAL developers " << PACKAGE_BUGREPORT << endl;
    *gmsg << "Time: " << timeStr << " date: " << dateStr << endl;


    /*
      Make a directory data for some of the output
    */
    if(Ippl::myNode() == 0) {
        int temp = umask(0);
        if((temp = mkdir("data", MY_MASK)) != 0) {
            *gmsg << "NOTE " << errno << ": unable to mkdir; " << strerror(errno) << endl;
        }
    }

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

    OpalData *OPAL = OpalData::getInstance();

    try {
        Configure::configure();
        int arg = 1;

        // Read startup file.
        FileStream::setEcho(true);

#ifndef __LIBCATAMOUNT__
        std::string startup = getenv("HOME");
        startup += "/init.opal";
        FileStream::setEcho(false);
        FileStream *is;

        try {
            is = new FileStream(startup);
        } catch(...) {
            is = 0;
            *gmsg << "No startup file \"" << startup << "\" found. Note: this is not mandatory for an OPAL simulation!" << endl;
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
            OPAL->storeInputFn(std::string(argv[1]));

        if(argc > 3) {
            if(argc > 5) {
                // will write dumping date into a new h5 file
                for(int ii = 2; ii < 6; ii = ii + 2)
                    // The sequence of the two arguments is free
                    if(std::string(argv[ii]) == std::string("-restart")) {
                        OPAL->setRestartRun();
                        OPAL->setRestartStep(atoi(argv[ii+1]));
                        OPAL->setRestartFileName(argv[1]);
                    } else if((std::string(argv[ii]) == std::string("-restartfn"))) {
                        OPAL->setRestartFileName(argv[ii+1]);
                    }
            } else {
                // will append dumping date into old h5 file
                if(std::string(argv[2]) == std::string("-restart")) {
                    OPAL->setRestartRun();
                    OPAL->setRestartStep(atoi(argv[3]));
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
                *gmsg << "* Reading input stream \"" << argv[arg] << "\"." << endl;
                parser.run(is);
                *gmsg << "* End of input stream \"" << argv[arg] << "\"." << endl;
            }
        }

        ////  DTA
        //  std::cout << std::endl << std::endl;
        //  vector<std::string> names=OPAL->getAllNames();
        //  printStringVector(names);
        //// /DTA


        IpplTimings::stopTimer(mainTimer);
        IpplTimings::print();
        IpplTimings::print(std::string("timing.dat"));

        if(Ippl::myNode() == 0) {
            std::ifstream errormsg("errormsg.txt");
            if(errormsg.good()) {
                char buffer[256];
                std::string closure("*                                                                                  *\n");
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

        Ippl::Comm->barrier();
        // cleanup/free global data
        Fieldmap::clearDictionary();
        OpalData::deleteInstance();
        delete gmsg;
        delete ippl;
        delete Ippl::Info;
        delete Ippl::Warn;
        delete Ippl::Error;
        delete Ippl::Debug;
        return 0;

    } catch(ClassicException &ex) {
        *gmsg << endl << "*** User error detected by function \""
              << ex.where() << "\"" << endl;
        abort();
    } catch(std::bad_alloc) {
        *gmsg << "Sorry, virtual memory exhausted." << endl;
        abort();
    } catch(...) {
        *gmsg << "Unexpected error." << endl;
        abort();
    }
}
