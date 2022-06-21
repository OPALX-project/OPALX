#include <gsl/gsl_errno.h>

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

#ifdef ENABLE_AMR
#include <AMReX_ParallelDescriptor.H>
#endif
/*
  Includes related to the optimizer
*/
#include "boost/smart_ptr.hpp"

#include "Pilot/Pilot.h"
#include "Util/CmdArguments.h"
#include "Util/OptPilotException.h"

#include "Optimizer/EA/FixedPisaNsga2.h"
#include "Optimizer/EA/BlendCrossover.h"
#include "Optimizer/EA/IndependentBitMutation.h"

#include "Optimize/OpalSimulation.h"

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

#include <gsl/gsl_errno.h>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <mpi.h>

#include <cstring>
#include <set>
#include <algorithm>
#include "Utilities/OpalException.h"
#define PYOPAL_GLOBALS_C
#include "PyOpal/Globals.h"

namespace {
    void errorHandlerGSL(const char *reason,
                         const char *file,
                         int line,
                         int gsl_errno) {
        throw OpalException(file, reason);
    }
}

namespace PyOpal {
namespace Globals {
void Initialise() {
    if (gmsg == NULL) {
        std::cerr << "GMSG init " << gmsg << std::endl;
        gmsg = new Inform("OPAL");
        std::cerr << "GMSG inited " << gmsg << std::endl;
    }
    if (ippl == NULL) {
        int argc = 3;

        char* argvr[] = {
            (char*)("pyopal"), 
            (char*)("--processes"),
            (char*)("3"),
            //(char*)("--commlib"),
            //(char*)("serial"),
            nullptr
        };
        char** argv = argvr;
        std::cerr << "INITIALISING IPPL" << std::endl;
        // Ippl is a typedef of IpplInfo in ippl/Utilities
        //MPI_Init(&argc, &argv);
        ippl = new Ippl(argc, argv);
        std::cerr << "INITIALISED  IPPL" << std::endl;
    }
    gsl_set_error_handler(&errorHandlerGSL);
}
}
}
