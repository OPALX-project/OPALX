/** A very simple driver using OPAL as forward solver and PISA (NSGA2) as
 *  optimizer.
 *  Only runs one optimization problem using MPI_COMM_WORLD for pilot,
 *  optimizer and workers.
 */

#include <mpi.h>
#include "boost/smart_ptr.hpp"

#include "Optimizer/PisaVariator.h"
#include "Util/OpalInputFileParser.h"
#include "Simulation/OpalSimulation.h"

#include "Comm/Splitter/OneMasterSplit.h"
#include "Comm/Topology/NoCommTopology.h"

#include "Pilot/Pilot.h"
#include "Util/CmdArguments.h"
#include "Util/OptPilotException.h"


int main(int argc, char** argv) {

    typedef OpalInputFileParser Input_t;
    typedef PisaVariator Opt_t;
    typedef OpalSimulation Sim_t;

    typedef CommSplitter<OneMasterSplit<NoCommTopology>, NoCommTopology> Comm_t;

    typedef Pilot< Input_t, Opt_t, Sim_t, Comm_t> pilot_t;

    MPI_Init(&argc, &argv);

    try {
        CmdArguments_t args(new CmdArguments(argc, argv));

        boost::shared_ptr<Comm_t>  comm(new Comm_t(args, MPI_COMM_WORLD));
        boost::scoped_ptr< pilot_t > pi(new pilot_t(args, comm));

    } catch (OptPilotException &e) {
        std::cout << "Exception caught: " << e.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, -100);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();

    return 0;
}
