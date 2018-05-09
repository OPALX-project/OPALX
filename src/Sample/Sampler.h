#ifndef __OPAL_SAMPLER_H__
#define __OPAL_SAMPLER_H__

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <fstream>


#include "Comm/types.h"
#include "Util/Types.h"
#include "Util/CmdArguments.h"

#include "Optimizer/Optimizer.h"
#include "Sample/SIndividual.h"
#include "Sample/SamplingOperator.h"

#include <boost/smart_ptr.hpp>
#include <boost/chrono.hpp>

#include <memory>

#include <queue>


/**
 *  \class Sampler
 *  \brief Implementing sampling
 */
// template<
//       template <class> class SamplingOperator
// >
class Sampler : public Optimizer {

public:

    /** This constructor should never be called.
     *  It's provided due to the inheritance of SamplePilot from
     *  Pilot
     * 
     */
    Sampler(Expressions::Named_t objectives,
            Expressions::Named_t constraints,
            DVarContainer_t dvars,
            size_t dim, Comm::Bundle_t comms,
            CmdArguments_t args);
    
    
    /**
     *  Retrieves all (for the sampler) relevant arguments specified on the
     *  command line, initializes the variator and sets up statistics and
     *  debug traces.
     *
     *  @param[in] svars of sampling
     *  @param[in] comms available to the sampler
     *  @param[in] args the user passed on the command line
     */
    Sampler(
                   Expressions::Named_t type,
                   DVarContainer_t dvars, Comm::Bundle_t comms,
                   CmdArguments_t args);

    ~Sampler();

    /// Starting selection algorithm and variator PISA state machine
    void initialize();
    
    /// type used in solution state exchange with other optimizers
    typedef std::vector< SIndividual > SolutionState_t;

protected:

    // implementing poller hooks
    bool onMessage(MPI_Status status, size_t length);
    void postPoll();

    void setupPoll() {}
    void prePoll() {
        // std::ostringstream debug;
        // debug << "IN PRE POLL: ";
        // debug << getStateString(curState_m) << std::endl;
        // progress_->log(debug);
    }
    void onStop() {}

    // helper sending evaluation requests to the pilot
    void dispatch_forward_solves();


private:
    
    int gid;

    int my_local_pid_;

    typedef SIndividual  Individual_t;

    /// communicator bundle for the optimizer
    Comm::Bundle_t comms_;

//     /// mapping from unique job ID to individual
    std::map<size_t, boost::shared_ptr<SIndividual > > jobmapping_m;
    
    std::vector< std::unique_ptr<SamplingOperator> > samplingOp_m;
    
    std::queue<boost::shared_ptr<Individual_t> > individuals_m;

    /// bounds on each specified gene
    bounds_t dVarBounds_m;
    
    Expressions::Named_t type_m;
    
    /// design variables
    DVarContainer_t dvars_m;
    
    
    int nsamples_m;
    

    /// command line arguments specified by the user
    CmdArguments_t args_;

    /// current generation
    int act_sample_m;
    
    int done_sample_m;
    
    enum State {
        SUBMIT,
        STOP,
        TERMINATE
    };
    
    State curState_m;
    

    /// result file name
    std::string resultFile_m;
    std::string resultDir_m;


    // dump frequency
    int dump_freq_;

    /// Dumps index, objective values and bit string of all individuals in
    /// global_population.
    void dumpPopulationToJSON();
    
    
    void runStateMachine();
    
    void initSamplingMethods_m();
    
    void createNewIndividual_m();
};

#include "Sample/Sampler.tcc"

#endif