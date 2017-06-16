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
#include "Utilities/Util.h"

#include "OPALconfig.h"

#ifdef HAVE_AMR_SOLVER
#include <ParallelDescriptor.H>
#endif

#include <gsl/gsl_errno.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <cstring>
#include <set>
#include <algorithm>

/*
  
  Includes related the to optimizer
  
 */

#include "boost/smart_ptr.hpp"

#include "Pilot/Pilot.h"
#include "Util/CmdArguments.h"
#include "Util/OptPilotException.h"

#include "Optimizer/EA/FixedPisaNsga2.h"
#include "Optimizer/EA/BlendCrossover.h"
#include "Optimizer/EA/IndependentBitMutation.h"

#include "Util/OpalInputFileParser.h"
#include "Simulation/OpalSimulation.h"

#include "Comm/CommSplitter.h"
#include "Comm/Topology/NoCommTopology.h"
#include "Comm/Splitter/ManyMasterSplit.h"
#include "Comm/MasterGraph/SocialNetworkGraph.h"

#include "Expression/Parser/function.hpp"
#include "Expression/FromFile.h"
#include "Expression/SumErrSq.h"
#include "Expression/SDDSVariable.h"
#include "Expression/RadialPeak.h"
#include "Expression/SumErrSqRadialPeak.h"
#include "Expression/ProbeVariable.h"

struct sameSDDSVariable {
    sameSDDSVariable(const std::string & base_filename) {
        size_t pos = base_filename.find_last_of("/");
        std::string tmplfile = base_filename;
        if(pos != std::string::npos)
            tmplfile = base_filename.substr(pos+1);
        pos = tmplfile.find(".");
        // std::string simName =
        stat_filename_ = tmplfile.substr(0,pos) + ".stat";
    }

    Expressions::Result_t operator()(client::function::arguments_t args) {
        args.push_back(stat_filename_);
        return var_(args);
    }

private:
    client::function::argument_t stat_filename_;
    SDDSVariable var_;
};


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


bool haveOpimiseRun(int argc, char *argv[]) {
  
  namespace fs = boost::filesystem;

  int arg = -1;
  std::string fname;

  for(int ii = 1; ii < argc; ++ ii) {
    std::string argStr = std::string(argv[ii]);
    if (argStr == std::string("--input")) {
      ++ ii;
      arg = ii;
      INFOMSG(argv[ii] << endl);
      continue;
    } else if (argStr == std::string("-restart") ||
	       argStr == std::string("--restart")) {
      return false;
    } else if (argStr == std::string("-restartfn") ||
	       argStr == std::string("--restartfn")) {
      return false;
    } else if (argStr == std::string("-version") ||
	       argStr == std::string("--version")) {
      return false;
    } else if (argStr == std::string("-help") ||
	       argStr == std::string("--help")) {
      return false;
    } else {
      if (arg == -1 &&
	  (ii == 1 || ii + 1 == argc) &&
	  argv[ii][0] != '-') {
	arg = ii;
	continue;
      } else {
	continue;
      }
    }
  }

  if (arg == -1) {
    INFOMSG("No input file provided!" << endl);
    exit(1);
  }

  fname = std::string(argv[arg]);
  if (!fs::exists(fname)) {
    INFOMSG("Input file \"" << fname << "\" doesn't exist!" << endl);
    exit(1);
  }

  std::ifstream inFile;
  inFile.open(fname);

  std::stringstream strStream;
  strStream << inFile.rdbuf();
  std::string str = strStream.str();

  std::transform(str.begin(), str.end(),str.begin(), ::toupper);

  const std::string s1("OPTIMIZE");
  const std::string s2("OBJECTIVE");
  const std::string s3("DVAR");
  
  bool res = (boost::algorithm::contains(str, s1) && 
	      boost::algorithm::contains(str, s2) &&
	      boost::algorithm::contains(str, s3));

  inFile.close();

  return res;
}

int mainOPALOptimiser(int argc, char *argv[]) {

  gmsg = new  Inform("OPAL");

  *gmsg << "We would start the optimiser here .... " << endl;

    MPI_Init(&argc, &argv);

    // Setup/Configuration
    //////////////////////////////////////////////////////////////////////////
    typedef OpalInputFileParser Input_t;
    typedef OpalSimulation Sim_t;

    typedef FixedPisaNsga2< BlendCrossover, IndependentBitMutation > Opt_t;

    //    typedef PisaVariator< BlendCrossover, IndependentBitMutation > Opt_t;

    typedef CommSplitter< ManyMasterSplit< NoCommTopology > > Comm_t;
    typedef SocialNetworkGraph< NoCommTopology > SolPropagationGraph_t;

    typedef Pilot<Input_t, Opt_t, Sim_t, SolPropagationGraph_t, Comm_t>
        pilot_t;

    // prepare function dictionary and add all available functions in
    // expressions
    functionDictionary_t funcs;
    client::function::type ff;
    ff = FromFile();
    funcs.insert(std::pair<std::string, client::function::type>
            ("fromFile", ff));
    ff = SumErrSq();
    funcs.insert(std::pair<std::string, client::function::type>
            ("sumErrSq", ff));
    ff = SDDSVariable();
    funcs.insert(std::pair<std::string, client::function::type>
            ("sddsVariableAt", ff));

    ff = RadialPeak();
    funcs.insert(std::pair<std::string, client::function::type>
            ("radialPeak", ff));
    ff = SumErrSqRadialPeak();
    funcs.insert(std::pair<std::string, client::function::type>
            ("sumErrSqRadialPeak", ff));

    ff = ProbeVariable();
    funcs.insert(std::pair<std::string, client::function::type>
            ("probVariableWithID", ff));

    //////////////////////////////////////////////////////////////////////////






  return 0;
}


int mainOPAL(int argc, char *argv[]) {
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
          << "This is OPAL (Object Oriented Parallel Accelerator Library) Version " << OPAL_VERSION_STR << "\n\n"
          << "                (c) PSI, http://amas.web.psi.ch" << endl
          << endl;

#ifdef OPAL_DKS
    *gmsg << "OPAL compiled with DKS (Dynamic Kernel Scheduler) Version "
	  << DKS_VERSION << endl << endl;
#endif

    *gmsg << "Please send cookies, goodies or other motivations (wine and beer ... ) \nto the OPAL developers " << PACKAGE_BUGREPORT << "\n" << endl;
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

        if(argc <= 1) {
            // Run commands from standard input
            parser.run(new TerminalStream("OPAL"));
        } else {
            int arg = -1;
            std::string fname;
            // if(argc > 3) {
            //     if(argc > 5) {
            //         // will write dumping date into a new h5 file
            for(int ii = 1; ii < argc; ++ ii) {
                std::string argStr = std::string(argv[ii]);
                // The sequence of the two arguments is free
                if (argStr == std::string("--input")) {
                    ++ ii;
                    arg = ii;
                    INFOMSG(argv[ii] << endl);
                    continue;
                } else if (argStr == std::string("-restart") ||
                           argStr == std::string("--restart")) {
                    opal->setRestartRun();
                    opal->setRestartStep(atoi(argv[++ ii]));
                    opal->setRestartFileName(argv[1]);
                    continue;
                } else if (argStr == std::string("-restartfn") ||
                           argStr == std::string("--restartfn")) {
                    opal->setRestartFileName(argv[++ ii]);
                    continue;
                } else if (argStr == std::string("-version") ||
                           argStr == std::string("--version")) {
                    INFOMSG("OPAL Version " << OPAL_VERSION_STR << ", git rev. " << Util::getGitRevision() << endl);
                    IpplInfo::printVersion(true);
                    std::string options = (IpplInfo::compileOptions() +
                                           std::string(" ") +
                                           std::string(PACKAGE_COMPILE_OPTIONS) +
                                           std::string(" "));
                    std::set<std::string> uniqOptions;
                    while (options.length() > 0) {
                        size_t n = options.find_first_of(' ');
                        while (n == 0) {
                            options = options.substr(n + 1);
                            n = options.find_first_of(' ');
                        }

                        uniqOptions.insert(options.substr(0, n));
                        options = options.substr(n + 1);
                    }
                    for (auto it: uniqOptions) {
                        options += it + " ";
                    }

                    std::string header("Compile-time options: ");
                    while (options.length() > 58) {
                        std::string line = options.substr(0, 58);
                        size_t n = line.find_last_of(' ');
                        INFOMSG(header << line.substr(0, n) << "\n");

                        header = std::string(22, ' ');
                        options = options.substr(n + 1);
                    }
                    INFOMSG(header << options << endl);
                    exit(0);
                } else if (argStr == std::string("-help") ||
                           argStr == std::string("--help")) {
                    IpplInfo::printHelp(argv);
                    INFOMSG("   --version            : Print a brief version summary.\n");
                    INFOMSG("   --input <fname>      : Specifies the input file <fname>.\n");
                    INFOMSG("   --restart <n>        : Performes a restart from step <n>.\n");
                    INFOMSG("   --restartfn <fname>  : Uses the file <fname> to restart from.\n");
                    INFOMSG("   --help               : Display this command-line summary.\n");
                    INFOMSG(endl);
                    exit(0);
                } else {
                    if (arg == -1 &&
                        (ii == 1 || ii + 1 == argc) &&
                        argv[ii][0] != '-') {
                        arg = ii;
                        continue;
                    } else {
                        INFOMSG("Unknown argument \"" << argStr << "\"" << endl);
                        IpplInfo::printHelp(argv);
                        INFOMSG("   --version            : Print a brief version summary.\n");
                        INFOMSG("   --input <fname>      : Specifies the input file <fname>.\n");
                        INFOMSG("   --restart <n>        : Performes a restart from step <n>.\n");
                        INFOMSG("   --restartfn <fname>  : Uses the file <fname> to restart from.\n");
                        INFOMSG("   --help               : Display this command-line summary.\n");
                        INFOMSG(endl);
                        exit(1);
                    }
                }
            }

            if (arg == -1) {
                INFOMSG("No input file provided!" << endl);
                exit(1);
            }

            fname = std::string(argv[arg]);
            if (!fs::exists(fname)) {
                INFOMSG("Input file \"" << fname << "\" doesn't exist!" << endl);
                exit(1);
            }

            opal->storeInputFn(fname);
            opal->setRestartFileName(opal->getInputBasename() + std::string(".h5"));


            FileStream *is;

            try {
                is = new FileStream(fname);
            } catch(...) {
                is = 0;
                *gmsg << "Input file \"" << fname << "\" not found." << endl;
            }

            if(is) {
                *gmsg << "* Reading input stream \"" << fname << "\"." << endl;
                parser.run(is);
                *gmsg << "* End of input stream \"" << fname << "\"." << endl;
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
      *gmsg << endl << "*** User error detected by function \"" << ex.where() << "\":\n"
            << ex.what() << endl;
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



// The OPAL main program.

int main(int argc, char *argv[]) {

  int res; 

  if ((argc <= 1) || !haveOpimiseRun(argc, argv))
    res = mainOPAL(argc, argv);
  else
    res = mainOPALOptimiser(argc, argv);
  return res;

}


