/** Stand alone PISA optimizer used for benchmarks.
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


#include "Util/PlainInputFileParser.h"
#include "Simulation/FunctionEvaluator.h"

#include "Comm/CommSplitter.h"
#include "Comm/Topology/NoCommTopology.h"
#include "Comm/Splitter/ManyMasterSplit.h"
#include "Comm/MasterGraph/SocialNetworkGraph.h"

#include "Expression/Parser/function.hpp"
#include "Expression/FromFile.h"
#include "Expression/SumErrSq.h"

//Test Problems
#include "Problems/FON.h"
//////////////////////////////////////////////////////////////////////////////


int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);

    // Setup/Configuration
    //////////////////////////////////////////////////////////////////////////
    typedef PlainInputFileParser Input_t;
    typedef FunctionEvaluator Sim_t;
    typedef FixedPisaNsga2< BlendCrossover, IndependentBitMutation >
        Opt_t;

    typedef SocialNetworkGraph< NoCommTopology > SolPropagationGraph_t;
    typedef CommSplitter< ManyMasterSplit<NoCommTopology> > Comm_t;

    typedef Pilot<Input_t, Opt_t, Sim_t, SolPropagationGraph_t, Comm_t>
        pilot_t;

    functionDictionary_t funcs;
    client::function::type ff;
    ff = FON();
    funcs.insert(std::pair<std::string, client::function::type> ("FON", ff));

    //////////////////////////////////////////////////////////////////////////


    try {
        CmdArguments_t args(new CmdArguments(argc, argv));

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
