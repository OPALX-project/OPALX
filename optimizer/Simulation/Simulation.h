#ifndef __SIMULATION_H__
#define __SIMULATION_H__

#include "Util/Types.h"
#include "Util/CmdArguments.h"

/**
 *  \class Simulation
 *  \brief Abstract class every simulation has to implement to be able to
 *  work with the optimization pilot.
 *
 *  To be supported the simulation code has to be called as a library taking
 *  as an argument an MPI communicator specifying on which processors the
 *  simulation runs.
 *
 *  TODO:
 *    - formulate the 'simulation' part as a pipeline.. workers can simply put
 *      together their own pipelines situation-dependent.
 */
class Simulation {

public:

    Simulation(CmdArguments_t args)
        : args_(args)
    {}

    ~Simulation()
    {}

    CmdArguments_t getArgs() { return args_; }

    /**
     *  Run the simulation (different for every simulation we specify).
     *  This method can block or return immediately, the worker is calling
     *  collectResults() to get the results. This means that the
     *  implementation of the Simulation class has to make sure that
     *  collectResults() waits until the data is available.
     */
    virtual void run() = 0;

    /**
     *  Collect all results from after the simulation has been executed. Make
     *  sure that data is available before returning from function.
     *
     *  @see run()
     */
    virtual void collectResults() = 0;

    /**
     *  Get all requested information.
     */
    virtual reqVarContainer_t getResults() = 0;

private:
    CmdArguments_t args_;

};

#endif
