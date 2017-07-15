#ifndef __PISA_VARIATOR_H__
#define __PISA_VARIATOR_H__

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
#include "Util/Statistics.h"

#include "Optimizer/Optimizer.h"
#include "Optimizer/EA/Individual.h"
#include "Optimizer/EA/Variator.h"

#include <boost/smart_ptr.hpp>
#include <boost/chrono.hpp>

#include "Util/Trace/Trace.h"


/**
 *  \class PisaVariator
 *  \brief Implementing the Variator for the PISA state machine.
 *  \deprecated Please use FixedPisaNsga2
 *
 *  @see http://www.tik.ee.ethz.ch/pisa/
 *
 *  The convergence behavior of the optimizer can be steered in 3 ways,
 *  corresponding command line arguments are given in brackets:
 *    - limit the number of generations (maxGenerations),
 *    - specify a target hypervolume (expected-hypervol) and tolerance
 *      (epsilon)
 *    - specify a minimal hypervolume progress (conv-hvol-prog), relative to
 *      the last generation, ((prev - new)/prev) that has to be attained to
 *      continue optimizing.
 */
template<
      template <class> class CrossoverOperator
    , template <class> class MutationOperator
>
class PisaVariator : public Optimizer {

public:

    /**
     *  Retrieves all (for the optimizer) relevant arguments specified on the
     *  command line, initializes the variator and sets up statistics and
     *  debug traces.
     *
     *  @param[in] objectives of optimization problem
     *  @param[in] constraints of optimization problem
     *  @param[in] dvars of optimization problem
     *  @param[in] dim number of objectives
     *  @param[in] comms available to the optimizer
     *  @param[in] args the user passed on the command line
     */
    PisaVariator(Expressions::Named_t objectives,
                 Expressions::Named_t constraints,
                 DVarContainer_t dvars, size_t dim, Comm::Bundle_t comms,
                 CmdArguments_t args);

    ~PisaVariator();

    /// Starting selection algorithm and variator PISA state machine
    void initialize();

    /// type used in solution state exchange with other optimizers
    typedef std::vector< Individual > SolutionState_t;


protected:

    /// Write the variator config file
    void writeVariatorCfg();

    // implementing poller hooks
    bool onMessage(MPI_Status status, size_t length);
    void postPoll();

    void setupPoll() {}
    void prePoll() {}
    void onStop() {}

    // helper sending evaluation requests to the pilot
    void dispatch_forward_solves();


private:

    int my_local_pid_;

    /// all PISA states
    enum PisaState_t {
          Initialize         = 0
        , InitializeSelector = 1
        , Variate            = 2
        , Select             = 3
        , Stop               = 4
        , VariatorStopped    = 5
        , SelectorStopped    = 7
        , Reset              = 8
        , ReadyForReset      = 9
        , Restart            = 11
    };

    /// the current state of the state machine
    PisaState_t curState_m;

    /// collect some statistics of rejected and accepted individuals
    boost::scoped_ptr<Statistics<size_t> > statistics_;

    /// type of our variator
    typedef Individual Individual_t;
    typedef Variator< Individual_t, CrossoverOperator, MutationOperator >
        Variator_t;
    boost::scoped_ptr<Variator_t> variator_m;

    /// communicator bundle for the optimizer
    Comm::Bundle_t comms_;

    /// buffer holding all finished job id's
    std::queue<unsigned int> finishedBuffer_m;

    /// mapping from unique job ID to individual
    std::map<size_t, boost::shared_ptr<Individual> > jobmapping_m;

    /// string of executable of the selection algorithm
    std::string execAlgo_m;

    /// indicating if initial population has been created
    bool notInitialized_m;


    /// bounds on each specified gene
    bounds_t dVarBounds_m;
    /// objectives
    Expressions::Named_t objectives_m;
    /// constraints
    Expressions::Named_t constraints_m;
    /// design variables
    DVarContainer_t dvars_m;

    /// command line arguments specified by the user
    CmdArguments_t args_;

    /// size of initial population
    size_t alpha_m;
    /// number of parents the selector chooses
    size_t mu_m;
    /// number of children the variator produces
    size_t lambda_m;
    /// number of objectives
    size_t dim_m;
    /// current generation
    size_t act_gen;
    /// maximal generation (stopping criterion)
    size_t maxGenerations_m;

    /// parameter file name
    std::string paramfile;
    /// base of FSM file names
    std::string filenamebase;
    /// var file name
    std::string var_file;
    /// sel file name
    std::string sel_file;
    /// config file name
    std::string cfg_file;
    /// ini file name
    std::string ini_file;
    /// archive file name
    std::string arc_file;
    /// sta file name
    std::string sta_file;
    /// result file name
    std::string resultFile_m;
    std::string resultDir_m;

    // dump frequency
    int dump_freq_;

    /// convergence accuracy if maxGenerations not set
    double hvol_eps_;
    double expected_hvol_;
    double current_hvol_;
    double conv_hvol_progress_;
    double hvol_progress_;

    /// file header for result files contains this parameter description
    std::string file_param_descr_;

    boost::chrono::system_clock::time_point run_clock_start_;
    boost::chrono::system_clock::time_point last_clock_;

    // DEBUG output helpers
    boost::scoped_ptr<Trace> job_trace_;
    boost::scoped_ptr<Trace> progress_;


    // entry point for starting the selector side of the PISA state machine
    void startSelector(std::string filename_base);

    /// executes one loop of the PISA state machine
    void runStateMachine();

    /// passes num_individuals to the selector
    void toSelectorAndCommit(int num_individuals, std::string filename);

    /// how often do we exchange solutions with other optimizers
    size_t exchangeSolStateFreq_m;

    /// if necessary exchange solution state with other optimizers
    void exchangeSolutionStates();


    /// Read archive written by the selector and remove individuals that did
    /// not survive
    void readArchive(std::set<unsigned int> &survivors);

    /// returns the current state
    int readState();

    /**
     *  Retrieve selected indivduals
     *  @param[inout] set of selected individual ids
     */
    void readSelectedIndividuals(std::vector<unsigned int> &parents);

    /// Reset a specific state machine file
    void clearFile(std::string filename);

    /**
     *  Write evaluated individuals for selection to file
     *  @param[in] file name to dump individuals
     *  @param[in] set of individual ids passed to selector
     *  @return 0 if successful, 1 otherwise
     */
    int writeIndividualsForSelector(std::string ofile,
                                    std::set<unsigned int> ids);

    /**
     *  Write state of FSM.
     *  @param[in] state to be written
     *  @return 0 if successful, 1 otherwise
     */
    int writeState(int state);

    /// Check if file is ready for writing
    bool isFileReady(std::string filename);

    /// Check if selector selected new parents for recombination
    bool areNewParentsAvailable();

    /// Dumps index, objective values and bit string of all individuals in
    /// global_population.
    void dumpPopulationToFile();
    void dumpPopulationToJSON();
};

#include "Optimizer/EA/PisaVariator.tcc"

#endif
