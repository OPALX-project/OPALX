#ifndef __FUNCTION_EVALUATOR_H__
#define __FUNCTION_EVALUATOR_H__

#include <string>
#include <map>
#include <vector>

#include "mpi.h"

#include "Util/Types.h"
#include "Util/CmdArguments.h"
#include "Simulation/Simulation.h"

/**
 *  \class FunctionEvaluator
 *  \brief Simply evaluates the set of objectives for a given point in design
 *         space.
 */
class FunctionEvaluator : public Simulation {

public:

    /**
     *  Setup.
     *
     *  @param[in] objectives of the optimization problem
     *  @param[in] constraints of the optimization problem
     *  @param[in] params
     *  @param[in] name of the simulation
     *  @param[in] comm MPI communicator used for running the simulation
     *  @param[in] args command line arguments passed to the framework
     */
    FunctionEvaluator(Expressions::Named_t objectives,
                      Expressions::Named_t constraints,
                      Param_t params, std::string name, MPI_Comm comm,
                      CmdArguments_t args);

    virtual ~FunctionEvaluator();

    /// run simulation: here we can relax and watch the show...
    void run() {}

    /// evaluate the expression using the params passed to the class
    void collectResults();

    /// returns container containing all requested variables with results
    reqVarContainer_t getResults() { return requestedVars_; }

private:

    /// holds solutions returned to the optimizer
    reqVarContainer_t requestedVars_;

    Expressions::Named_t objectives_;
    Expressions::Named_t constraints_;
    Param_t params_;

    MPI_Comm comm_;
};

#endif