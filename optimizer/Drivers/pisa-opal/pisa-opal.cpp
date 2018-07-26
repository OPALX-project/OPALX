/** A very simple driver using OPAL as forward solver and PISA (NSGA2) as
 *  optimizer.
 *  Only runs one optimization problem using MPI_COMM_WORLD for pilot,
 *  optimizer and workers.
 */

#include <mpi.h>
#include "boost/smart_ptr.hpp"

#include "Pilot/Pilot.h"
#include "Util/CmdArguments.h"
#include "Util/OptPilotException.h"


// Setup/Configuration
//////////////////////////////////////////////////////////////////////////////
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
#include "Expression/MaxNormRadialPeak.h"
#include "Expression/NumberOfPeaks.h"
#include "Expression/ProbeVariable.h"
//////////////////////////////////////////////////////////////////////////////

// Define statics from OPAL main
#include "opal.h"
Ippl *ippl;
Inform *gmsg;

int main(int argc, char** argv) {
    ippl = new Ippl(argc, argv);
    gmsg = new  Inform("OPAL");

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
    ff = MaxNormRadialPeak();
    funcs.insert(std::pair<std::string, client::function::type>
            ("maxNormRadialPeak", ff));
    ff = NumberOfPeaks();
    funcs.insert(std::pair<std::string, client::function::type>
            ("numberOfPeaks", ff));
    ff = SumErrSqRadialPeak();
    funcs.insert(std::pair<std::string, client::function::type>
            ("sumErrSqRadialPeak", ff));

    ff = ProbeVariable();
    funcs.insert(std::pair<std::string, client::function::type>
            ("probVariableWithID", ff));

    //////////////////////////////////////////////////////////////////////////

    try {
        CmdArguments_t args(new CmdArguments(argc, argv));

        std::string fname = args->getArg<std::string>("inputfile", true);
        ff = sameSDDSVariable(fname);
        funcs.insert(std::pair<std::string, client::function::type>
                     ("sameSDDSVariableAt", ff));

        boost::shared_ptr<Comm_t>  comm(new Comm_t(args, MPI_COMM_WORLD));
        boost::scoped_ptr<pilot_t> pi(new pilot_t(args, comm, funcs));

    } catch (OptPilotException &e) {
        std::cout << "Exception caught: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -100);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
