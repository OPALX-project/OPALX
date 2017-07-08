/** A very simple driver using OPAL as forward solver and PISA (NSGA2) as
 *  optimizer.
 *  Only runs one optimization problem using MPI_COMM_WORLD for pilot,
 *  optimizer and workers.
 */

#include <mpi.h>
#include "boost/smart_ptr.hpp"

#include "Util/OpalInputFileParser.h"
#include "Optimizer/PisaVariator.h"
#include "Simulation/OpalSimulation.h"
#include "Comm/Splitter/OneMasterSplit.h"
#include "Comm/Topology/HwLocCommTopology.h"


#include "Util/CmdArguments.h"
#include "Pilot/Pilot.h"
#include "Util/OptPilotException.h"


int main(int argc, char** argv) {

    typedef OpalInputFileParser Input_t;
    typedef PisaVariator Opt_t;
    typedef OpalSimulation Sim_t;
    //typedef OneMasterSplit SplitterPolicy;
    typedef HwLocCommTopology TopoDiscoveryPolicy;

    //FIXME: use template <class> class SplitterPolicy to be able to specify
    //       just OneMasterSplit as policy (no template stacking).
    typedef Pilot< Input_t, Opt_t, Sim_t,
                   OneMasterSplit<HwLocCommTopology>,
                   TopoDiscoveryPolicy > pilot_t;

    MPI_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;

    try {
        CmdArguments_t args(new CmdArguments(argc, argv));

        boost::scoped_ptr< pilot_t > pi(new pilot_t(args, comm));

    } catch (OptPilotException &e) {
        std::cout << "Exception caught: " << e.what() << std::endl;
        MPI_Abort(comm, -1);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
