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
#include "Fields/Fieldmap.h"
#include "FixedAlgebra/FTps.h"

#include "BasicActions/Option.h"
#include "Utilities/Options.h"
#include "Utilities/Options.h"
#include "Utilities/OpalException.h"

#include "OPALconfig.h"

#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#ifdef HAVE_AMR_SOLVER
#include <ParallelDescriptor.H>
#endif

#include <gsl/gsl_errno.h>

#include <boost/filesystem.hpp>

#include <cstring>

//  DTA
#define NC 5

void printStringVector(const std::vector<std::string> &strings) {
    unsigned int iend = strings.size(), nc = 0;
    for(unsigned int i = 0; i < iend; ++i) {
        std::cout << "  " << strings[i];
        if(++nc == NC) { nc = 0; std::cout << std::endl; }
    }
    if(nc != 0) std::cout << std::endl;
    std::cout << std::endl;
}

void errorHandlerGSL(const char *reason,
                     const char *file,
                     int line,
                     int gsl_errno);
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

    namespace fs = boost::filesystem;

#ifdef HAVE_AMR_SOLVER
    // false: build no parmparse, we use the OPAL parser instead.
    BoxLib::Initialize(argc, argv, false, Ippl::getComm());
#endif

    OPALTimer::Timer simtimer;

    std::string dateStr(simtimer.date());
    std::string timeStr(simtimer.time());

    H5SetVerbosityLevel(0); //65535);

    gsl_set_error_handler(&errorHandlerGSL);

    static IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("mainTimer");
    IpplTimings::startTimer(mainTimer);

    Inform hmsg("");
    std::string mySpace("            ");

    if(Ippl::myNode() == 0) remove("errormsg.txt");

    hmsg << mySpace <<  "   ____  _____       ___ " << endl;
    hmsg << mySpace <<  "  / __ \\|  __ \\ /\\   | | " << endl;
    hmsg << mySpace <<  " | |  | | |__) /  \\  | |" << endl;
    hmsg << mySpace <<  " | |  | |  ___/ /\\ \\ | |" << endl ;
    hmsg << mySpace <<  " | |__| | |  / ____ \\| |____" << endl;
    hmsg << mySpace <<  "  \\____/|_| /_/    \\_\\______|" << endl;


    *gmsg << endl
          << "This is OPAL (Object Oriented Parallel Accelerator Library) Version "
	  << OPAL_PACKAGE_VERSION << "\n\n"
          << "                (c) PSI, http://amas.web.psi.ch" << endl
          << endl;

#ifdef OPAL_DKS
    *gmsg << "OPAL compiled with DKS (Dynamic Kernel Scheduler) Version " 
	  << DKS_VERSION;
    if (IpplInfo::DKSEnabled)
      *gmsg << " GPU present" << endl << endl;
    else
      *gmsg << " GPU not present" << endl << endl;
#endif

    *gmsg << "Please send cookies, goodies or other motivations (wine and beer ... ) \nto the OPAL developers "
	  << OPAL_PACKAGE_BUGREPORT << "\n" << endl;
    *gmsg << "Time: " << timeStr << " date: " << dateStr << "\n" << endl;

    /*
      Make a directory data for some of the output
    */
    if(Ippl::myNode() == 0) {
        if (!fs::exists("data")) {
            boost::system::error_code error_code;
            if (!fs::create_directory("data", error_code)) {
                std::cerr << error_code.message() << std::endl;
                // use error code to prevent create_directory from throwing an exception
            }
        }
    }
    Ippl::Comm->barrier();
    if (!fs::is_directory("data")) {
        std::cerr << "unable to create directory; aborting" << std::endl;
        abort();
    }

    const OpalParser parser;

    //  DTA
    std::cout.precision(16);
    std::cout.setf(std::ios::scientific, std::ios::floatfield);
    std::cerr.precision(16);
    std::cerr.setf(std::ios::scientific, std::ios::floatfield);
    // /DTA

    // Set global truncation orders.
    FTps<double, 2>::setGlobalTruncOrder(20);
    FTps<double, 4>::setGlobalTruncOrder(15);
    FTps<double, 6>::setGlobalTruncOrder(10);

    OpalData *opal = OpalData::getInstance();

    try {
        Configure::configure();
        int arg = 1;

        // Read startup file.
        FileStream::setEcho(Options::echo);

        char *startup = getenv("HOME");
	boost::filesystem::path p = strncat(startup, "/init.opal", 20);
	if (startup != NULL && is_regular_file(p)) {
	  
	    FileStream::setEcho(false);
	    FileStream *is;

            try {
                is = new FileStream(startup);
            } catch(...) {
                is = 0;
                ERRORMSG("Could not open startup file \"" << startup << "\".\n"
                         << "Note: this is not mandatory for an OPAL simulation!\n");
            }

            if(is) {
                *gmsg << "Reading startup file \"" << startup << "\"." << endl;
                parser.run(is);
                *gmsg << "Finished reading startup file." << endl;
            }
            FileStream::setEcho(Options::echo);
        } else {
            *gmsg << "Couldn't find startup file \"" << startup << "\".\n"
                  << "Note: this is not mandatory for an OPAL simulation!\n" << endl;
        }

        if(argc > 1) {
            opal->storeInputFn(std::string(argv[1]));
            opal->setRestartFileName(opal->getInputBasename() + std::string(".h5"));
        }

        if(argc > 3) {
            if(argc > 5) {
                // will write dumping date into a new h5 file
                for(int ii = 2; ii < 6; ii = ii + 2)
                    // The sequence of the two arguments is free
                    if(std::string(argv[ii]) == std::string("-restart")) {
                        opal->setRestartRun();
                        opal->setRestartStep(atoi(argv[ii+1]));
                        opal->setRestartFileName(argv[1]);
                    } else if((std::string(argv[ii]) == std::string("-restartfn"))) {
                        opal->setRestartFileName(argv[ii+1]);
                    }
            } else {
                // will append dumping date into old h5 file
                if(std::string(argv[2]) == std::string("-restart")) {
                    opal->setRestartRun();
                    opal->setRestartStep(atoi(argv[3]));
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


        IpplTimings::stopTimer(mainTimer);

        IpplTimings::print();

	IpplTimings::print(std::string("timing.dat"));

        if(Ippl::myNode() == 0) {
            std::ifstream errormsg("errormsg.txt");
            if(errormsg.good()) {
                char buffer[256];
                std::string closure("                                                                                 *\n");
                ERRORMSG("\n"
                         << "* **********************************************************************************\n"
                         << "* ************** W A R N I N G / E R R O R * * M E S S A G E S *********************\n"
                         << "* **********************************************************************************"
                         << endl);
                errormsg.getline(buffer, 256);
                while(errormsg.good()) {
                    ERRORMSG("* ");
                    if(errormsg.gcount() == 1) {
                        ERRORMSG(closure);
                    } else if ((size_t)errormsg.gcount() <= closure.size()) {
                        ERRORMSG(buffer << closure.substr(errormsg.gcount() - 1));
                    } else {
                        ERRORMSG(buffer << endl);
                    }
                    errormsg.getline(buffer, 256);
                }
                ERRORMSG("* " << closure
                         << "* **********************************************************************************\n"
                         << "* **********************************************************************************"
                         << endl);
            }
            errormsg.close();
        }

        Ippl::Comm->barrier();
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

    } catch(std::bad_alloc &) {
      *gmsg << "Sorry, virtual memory exhausted." << endl;
      abort();

    } catch(std::exception const& e) {
      *gmsg << "Exception: " << e.what() << "\n";
       abort();

    } catch(...) {
      *gmsg << "Unexpected exception." << endl;
      abort();
    }
}

void errorHandlerGSL(const char *reason,
                     const char *file,
                     int,
                     int) {
    throw OpalException(file, reason);
}
