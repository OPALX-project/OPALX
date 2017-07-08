#ifndef __HOMOTOPY_H__
#define __HOMOTOPY_H__

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <utility>
#include <fstream>
#include <queue>

#include "Comm/types.h"
#include "Util/Types.h"
#include "Util/CmdArguments.h"
#include "Util/Statistics.h"

#include "Optimizer/Optimizer.h"

#include <boost/smart_ptr.hpp>
#include <boost/chrono.hpp>

#include "Util/Trace/Trace.h"

#include "HomotopyTypes.h"
using namespace Homotopy;


/**
 *  \class Homotopy
 *  \brief ...
 */
class HomotopyImpl : public Optimizer {

public:

    /**
     *  Setup PISA variator.
     *
     *  @param[in] objectives of optimization problem
     *  @param[in] constraints of optimization problem
     *  @param[in] dvars of optimization problem
     *  @param[in] dim number of objectives
     *  @param[in] comms access to various communicators
     *  @param[in] args access to command line arguments
     */
    HomotopyImpl(Expressions::Named_t objectives, Expressions::Named_t constraints,
             DVarContainer_t dvars, size_t dim, Comm::Bundle_t comms,
             CmdArguments_t args);

    /// Destructor.
    ~HomotopyImpl();

    void initialize();

    //FIXME: types of solution states
    typedef typename Homotopy::nodeEnvelope_t solution_t;
    typedef std::vector< typename Homotopy::nodeEnvelope_t > SolutionState_t;

protected:

    bool onMessage(MPI_Status status, size_t length);
    void postPoll();

    void setupPoll() {}
    void prePoll()   {}
    void onStop()    {}


private:

    //Functions required for homotopy
    bool onMessage(MPI_Status status, size_t length, std::queue< std::pair< size_t, objVars_t> >& q, std::queue< typename Homotopy::nodeEnvelope_t > &g);
    void dispatch_forward_solves( size_t jid, const Homotopy::designVars_t & vars);
    void exchangeSolutionStates( std::queue< typename Homotopy::NeighborMessage_t > & toSend );
    void SendStop();
    //TODO: encapsulate these details in a request manager
    std::map<size_t , std::pair< MPI_Request*, std::pair<size_t*, char*> > > reqs;
    std::map<size_t , std::pair< MPI_Request*, size_t*> > reqs_ids;
    std::map<size_t , Homotopy::designVars_t > debugVals;

    class ExchangeRequest {
    public:
        MPI_Request * buf_size_req;
        MPI_Request * buffer_req;
        size_t buf_size;
        char * buffer;

        ExchangeRequest(std::ostringstream &os);
        bool isComplete();
        ~ExchangeRequest();
    };
    std::vector< ExchangeRequest* > ExchangeState;


    int my_local_pid_;

    boost::scoped_ptr<Statistics<size_t> > statistics_;

    Comm::Bundle_t comms_;

    /// buffer holding all finished job id's
    std::queue<unsigned int> finishedBuffer_;

    /// mapping from unique job ID to solutions
    std::map<size_t, boost::shared_ptr<solution_t> > jobmapping_;

    // Optimization problem
    /// bounds on each specified gene
    bounds_t dVarBounds_;
    /// objectives
    Expressions::Named_t objectives_;
    /// constraints
    Expressions::Named_t constraints_;
    /// design variables
    DVarContainer_t dvars_;

    CmdArguments_t args_;

    /// number of objectives
    size_t dim_;
    /// number of iterations
    size_t iter_;

    /// result file name
    std::string resultFile_;
    std::string resultDir_;
    std::string file_param_descr_;

    double epsilon_;

    boost::chrono::system_clock::time_point run_clock_start_;
    boost::chrono::system_clock::time_point last_clock_;

    // DEBUG output helpers
    boost::scoped_ptr<Trace> job_trace_;
    boost::scoped_ptr<Trace> progress_;

    /// how often do we exchange solutions with other optimizers
    size_t exchangeSolStateFreq_m;


    void run();

    /// if necessary exchange solution state with other optimizers
    void exchangeSolutionStates();

    void dispatch_forward_solves();

    void dumpStateToFile();
};

#endif
