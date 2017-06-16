#include <iostream>
#include <sstream>
#include <cassert>
#include <cmath>
#include <cstring>
#include <set>
#include <cstdlib>
#include <string>
#include <limits>

#include <sys/stat.h>

#include "boost/algorithm/string.hpp"
#include "boost/foreach.hpp"
#define foreach BOOST_FOREACH

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/map.hpp>

#include "Util/OptPilotException.h"
#include "Util/MPIHelper.h"

#include "Util/Trace/TraceComponent.h"
#include "Util/Trace/Timestamp.h"
#include "Util/Trace/FileSink.h"


template< template <class> class CO, template <class> class MO >
PisaVariator<CO, MO>::PisaVariator(Expressions::Named_t objectives,
                           Expressions::Named_t constraints,
                           DVarContainer_t dvars,
                           size_t dim, Comm::Bundle_t comms,
                           CmdArguments_t args)
             : Optimizer(comms.opt)
             , statistics_(new Statistics<size_t>("individuals"))
             , comms_(comms)
             , objectives_m(objectives)
             , constraints_m(constraints)
             , dvars_m(dvars)
             , args_(args)
             , dim_m(dim)
             , act_gen(0)
{
    my_local_pid_ = 0;
    MPI_Comm_rank(comms_.opt, &my_local_pid_);

    //FIXME: proper rand gen initialization (use boost?!)
    srand(time(NULL) + comms_.island_id);

    dump_freq_ = args->getArg<int>("dump-freq", 1, false);
    maxGenerations_m = args->getArg<int>("maxGenerations", true);
    resultFile_m = args->getArg<std::string>("outfile", "-th_generation.dat", false);
    resultDir_m = args->getArg<std::string>("outdir", "generations", false);

    // create output directory if it does not exists
    struct stat dirInfo;
    if(stat(resultDir_m.c_str(),&dirInfo) != 0)
        mkdir((const char*)(resultDir_m.c_str()), 0777);

    // solution exchange frequency
    exchangeSolStateFreq_m = args->getArg<size_t>("sol-synch", 0, false);

    // convergence arguments
    hvol_eps_           = args->getArg<double>("epsilon", 1e-3, false);
    expected_hvol_      = args->getArg<double>("expected-hypervol", 0.0,
                                               false);
    conv_hvol_progress_ = args->getArg<double>("conv-hvol-prog", 0.0, false);
    current_hvol_       = std::numeric_limits<double>::max();
    hvol_progress_      = std::numeric_limits<double>::max();

    //XXX: we can also set alpha_m to number of workers
    alpha_m = args->getArg<int>("initialPopulation", true);
    lambda_m = 2;
    mu_m = 2;

    file_param_descr_ = "%ID,";

    Expressions::Named_t::iterator it;
    for(it = objectives_m.begin(); it != objectives_m.end(); it++)
        file_param_descr_ += '%' + it->first + ',';

    file_param_descr_ += " DVAR: ";

    DVarContainer_t::iterator itr;
    for(itr = dvars_m.begin(); itr != dvars_m.end(); itr++) {
        file_param_descr_ += '%' + boost::get<VAR_NAME>(itr->second) + ',';
        dVarBounds_m.push_back(
                std::pair<double, double>
                    (boost::get<LOWER_BOUND>(itr->second),
                     boost::get<UPPER_BOUND>(itr->second)));
    }
    file_param_descr_ = file_param_descr_.substr(0,
                                                 file_param_descr_.size()-1);


    // setup variator
    variator_m.reset(new Variator_t(alpha_m, mu_m, lambda_m, dim, dVarBounds_m,
                                    args));


    // Traces and statistics
    std::ostringstream trace_filename;
    trace_filename << "opt.trace." << comms_.island_id;
    job_trace_.reset(new Trace("Optimizer Job Trace"));
    job_trace_->registerComponent( "sink",
            boost::shared_ptr<TraceComponent>(
                new FileSink(trace_filename.str())));

    std::ostringstream prog_filename;
    prog_filename << "opt.progress." << comms_.island_id;
    progress_.reset(new Trace("Optimizer Progress"));
    progress_->registerComponent( "timestamp",
            boost::shared_ptr<TraceComponent>(new Timestamp()));
    progress_->registerComponent( "sink",
            boost::shared_ptr<TraceComponent>(
                new FileSink(prog_filename.str())));

    statistics_->registerStatistic("accepted", 0);
    statistics_->registerStatistic("infeasible", 0);
}

template<
      template <class> class CO
    , template <class> class MO
>
PisaVariator<CO, MO>::~PisaVariator()
{}

template<
      template <class> class CO
    , template <class> class MO
>
void PisaVariator<CO, MO>::initialize() {

    std::ostringstream filename;
    filename << "PISA_" << comms_.island_id << "_";
    std::string filename_base = filename.str();
    filename << "cfg";

    std::ofstream pisa_config_file;
    pisa_config_file.open(filename.str().c_str(), std::ios::out);
    pisa_config_file << "alpha "  << alpha_m             << std::endl;
    pisa_config_file << "mu "     << mu_m                << std::endl;
    pisa_config_file << "lambda " << lambda_m            << std::endl;
    pisa_config_file << "dim "    << objectives_m.size() << std::endl;
    pisa_config_file.close();

    startSelector(filename_base);

    var_file = filename_base + "var";
    sel_file = filename_base + "sel";
    cfg_file = filename_base + "cfg";
    ini_file = filename_base + "ini";
    arc_file = filename_base + "arc";
    sta_file = filename_base + "sta";

    /* creating files and writing 0 in there */
    std::ofstream file;
    file.open(var_file.c_str(), std::ios::out);
    file << 0;
    file.close();

    file.open(sel_file.c_str(), std::ios::out);
    file << 0;
    file.close();

    file.open(ini_file.c_str(), std::ios::out);
    file << 0;
    file.close();

    file.open(arc_file.c_str(), std::ios::out);
    file << 0;
    file.close();


    curState_m = Initialize;
    writeState(curState_m);

    notInitialized_m = true;

    // start poll loop
    run_clock_start_ = boost::chrono::system_clock::now();
    last_clock_ = boost::chrono::system_clock::now();
    run();
    system("killall nsga2");

    current_hvol_ =
        variator_m->population()->computeHypervolume(comms_.island_id);

    boost::chrono::duration<double> total =
        boost::chrono::system_clock::now() - run_clock_start_;
    std::ostringstream stats;
    stats << "__________________________________________" << std::endl;
    stats << "GENERATION " <<  act_gen << std::endl;
    stats << "TOTAL = " << total.count() << "s" << std::endl;
    stats << "HYPERVOLUME = " << current_hvol_ << std::endl;
    stats << "time per accepted ind   = "
          << total.count()/statistics_->getStatisticValue("accepted")
          << std::endl;
    stats << "time per infeasible ind = "
          << total.count()/statistics_->getStatisticValue("infeasible")
          << std::endl;
    statistics_->dumpStatistics(stats);
    stats << "__________________________________________" << std::endl;
    progress_->log(stats);

}


template<
      template <class> class CO
    , template <class> class MO
>
void PisaVariator<CO, MO>::startSelector(std::string filename_base) {

    // kill old nsga2 executables still running
    std::ostringstream sp;
    sp << "./nsga2/nsga2_param.txt " << filename_base << " 0.1";

    std::string selector_exe =
        args_->getArg<std::string>("selector", "./nsga2/nsga2", false);
    std::string selector_params =
        args_->getArg<std::string>("selector-params", sp.str(), false);

    struct stat fileInfo;
    if(stat(selector_exe.c_str(), &fileInfo) != 0) {
        std::cout << "Could not find selector executable (" << selector_exe
                << ").." << std::endl;
        throw OptPilotException("PisaVariator::initialize",
                                "Invalid selector executable");
    }

    std::vector<std::string> exe_name;
    boost::split(exe_name, selector_exe, boost::is_any_of("/"),
                 boost::token_compress_on);

    std::ostringstream kill;
    kill << "killall " << exe_name.back();
    system(kill.str().c_str());

    std::ostringstream run;
    run << selector_exe << " " << selector_params << "&";
    system(run.str().c_str());
}


template< template <class> class CO, template <class> class MO >
bool PisaVariator<CO, MO>::onMessage(MPI_Status status, size_t length) {

    typedef typename PisaVariator::Individual_t individual;

    MPITag_t tag = MPITag_t(status.MPI_TAG);
    switch(tag) {

    case EXCHANGE_SOL_STATE_RES_SIZE_TAG: {

        size_t buf_size = length;
        size_t pilot_rank = status.MPI_SOURCE;

        std::ostringstream dump;
        dump << "new results from other cores " << buf_size << std::endl;
        job_trace_->log(dump);

        char *buffer = new char[buf_size];
        MPI_Recv(buffer, buf_size, MPI_CHAR, pilot_rank,
                 MPI_EXCHANGE_SOL_STATE_RES_TAG, comms_.opt, &status);

        dump.clear();
        dump.str(std::string());
        dump << "got results from other cores " << buf_size << std::endl;
        job_trace_->log(dump);

        SolutionState_t new_states;
        std::istringstream is(buffer);
        boost::archive::text_iarchive ia(is);
        ia >> new_states;
        delete[] buffer;

        std::set<unsigned int> new_state_ids;
        foreach(individual ind, new_states) {

            // only insert individual if not already in population
            if(variator_m->population()->isRepresentedInPopulation(ind.genes))
                continue;

            boost::shared_ptr<individual> new_ind(new individual);
            new_ind->genes = std::vector<double>(ind.genes);
            new_ind->objectives = std::vector<double>(ind.objectives);

            //XXX:   can we pass more than lambda_m files to selector?
            unsigned int id =
                variator_m->population()->add_individual(new_ind);
            finishedBuffer_m.push(id);

            dump.clear();
            dump.str(std::string());
            dump << "individual (ID: " << id
                 << ") successfully migrated from another island" << std::endl;
            job_trace_->log(dump);
        }

        return true;
    }

    case REQUEST_FINISHED: {

        unsigned int jid = static_cast<unsigned int>(length);
        typename std::map<size_t, boost::shared_ptr<individual> >::iterator it;
        it = jobmapping_m.find(jid);

        std::ostringstream dump;
        dump << "job with ID " << jid << " delivered results" << std::endl;
        job_trace_->log(dump);

        if(it == jobmapping_m.end()) {
            dump << "\t |-> NOT FOUND!" << std::endl;
            job_trace_->log(dump);
            std::cout << "NON-EXISTING JOB with ID = " << jid << std::endl;
            throw OptPilotException("PisaVariator::onMessage",
                    "non-existing job");
        }

        boost::shared_ptr<individual> ind = it->second;
        jobmapping_m.erase(it);

        //size_t dummy = 1;
        //MPI_Send(&dummy, 1, MPI_UNSIGNED_LONG, status.MPI_SOURCE,
        //         MPI_WORKER_FINISHED_ACK_TAG, comms_.listen);

        reqVarContainer_t res;
        MPI_Recv_reqvars(res, status.MPI_SOURCE, comms_.opt);

        ind->objectives.clear();

        //XXX: check order of genes
        reqVarContainer_t::iterator itr;
        std::map<std::string, double> vars;
        for(itr = res.begin(); itr != res.end(); itr++) {
            if(itr->second.is_valid) {
                if(itr->second.value.size() == 1)
                    ind->objectives.push_back(itr->second.value[0]);
            } else {
                std::ostringstream dump;
                dump << "invalid individual: objective/constraint \""
                     << itr->first << "\" failed to evaluate (to true)!"
                     << std::endl;
                job_trace_->log(dump);

                variator_m->infeasible(ind);
                statistics_->changeStatisticBy("infeasible", 1);
                dispatch_forward_solves();
                return true;
            }
        }

        finishedBuffer_m.push(jid);
        statistics_->changeStatisticBy("accepted", 1);

        return true;
    }

    default: {
        std::cout << "(PisaVariator) Error: unexpected MPI_TAG: "
                  << status.MPI_TAG << std::endl;
        return false;
    }
    }
}


template< template <class> class CO, template <class> class MO >
void PisaVariator<CO, MO>::postPoll() {

    // whenever lambda children are ready and we are in the variator phase
    // run the selector
    if(finishedBuffer_m.size() >= lambda_m && curState_m == Variate) {
        //std::cout << "▉";
        std::cout << "░" << std::flush;

        double hvol =
            variator_m->population()->computeHypervolume(comms_.island_id);
        hvol_progress_ = fabs(current_hvol_ - hvol) / current_hvol_;
        current_hvol_ = hvol;

        boost::chrono::duration<double> total =
            boost::chrono::system_clock::now() - run_clock_start_;
        boost::chrono::duration<double> dt    =
            boost::chrono::system_clock::now() - last_clock_;
        last_clock_ = boost::chrono::system_clock::now();
        std::ostringstream stats;
        stats << "__________________________________________" << std::endl;
        stats << "Arriving at generation " << act_gen + 1     << std::endl;
        stats << "dt = " << dt.count() << "s, total = " << total.count()
              << "s" << std::endl;
        stats << "Hypervolume = " << current_hvol_            << std::endl;
        stats << "__________________________________________" << std::endl;
        progress_->log(stats);

        // dump parents of current generation for visualization purposes
        if((act_gen + 1) % dump_freq_ == 0) {
            dumpPopulationToFile();
            dumpPopulationToJSON();
        }

        //XXX we can only use lambda_m here (selector does not support
        //    variable size). Change selector?
        toSelectorAndCommit(lambda_m, var_file);

        exchangeSolutionStates();

        //XXX: at the end of the Variate state (previous runStateMachine()
        //     call) we can safely change the state.
        curState_m = Select;
        writeState(curState_m);
        act_gen++;
    }

    // state will be reset to whatever is in the state file
    runStateMachine();
}


template<
      template <class> class CO
    , template <class> class MO
>
void PisaVariator<CO, MO>::exchangeSolutionStates() {

    typedef typename PisaVariator::Individual_t individual;

    size_t num_masters = args_->getArg<size_t>("num-masters", 1, false);

    if(num_masters <= 1 ||
       exchangeSolStateFreq_m == 0 ||
       act_gen % exchangeSolStateFreq_m != 0)
        return;

    int pilot_rank = comms_.master_local_pid;

    std::ostringstream os;
    boost::archive::text_oarchive oa(os);

    SolutionState_t population;
    typename std::map<unsigned int, boost::shared_ptr<individual> >::iterator itr;
    for(itr  = variator_m->population()->begin();
        itr != variator_m->population()->end(); itr++) {

        individual ind;
        ind.genes = std::vector<double>(itr->second->genes);
        ind.objectives = std::vector<double>(itr->second->objectives);
        population.push_back(ind);
    }

    oa << population;

    size_t buf_size = os.str().length();

    std::ostringstream dump;
    dump << "sending my buffer size " << buf_size << " bytes to PILOT"
         << std::endl;
    job_trace_->log(dump);

    MPI_Send(&buf_size, 1, MPI_UNSIGNED_LONG, pilot_rank,
             EXCHANGE_SOL_STATE_TAG, comms_.opt);

    char *buffer = new char[buf_size];
    memcpy(buffer, os.str().c_str(), buf_size);
    MPI_Send(buffer, buf_size, MPI_CHAR, pilot_rank,
             MPI_EXCHANGE_SOL_STATE_DATA_TAG, comms_.opt);
    delete[] buffer;

    dump.clear();
    dump.str(std::string());
    dump << "Sent " << buf_size << " bytes to PILOT" << std::endl;
    job_trace_->log(dump);
}


template< template <class> class CO, template <class> class MO >
void PisaVariator<CO, MO>::toSelectorAndCommit(int num_individuals,
                                       std::string filename) {

    std::set<unsigned int> ids;
    for(int i=0; i<num_individuals; i++) {
        unsigned int id = finishedBuffer_m.front();
        ids.insert(id);
        finishedBuffer_m.pop();
    }

    writeIndividualsForSelector(filename, ids);
    variator_m->population()->commit_individuals(ids);
}


template< template <class> class CO, template <class> class MO >
void PisaVariator<CO, MO>::dispatch_forward_solves() {

    typedef typename PisaVariator::Individual_t individual;

    while(variator_m->hasMoreIndividualsToEvaluate()) {

        //reqVarContainer_t reqs;
        //Expressions::Named_t::iterator it;
        //for(it = objectives_m.begin(); it != objectives_m.end(); it++) {
            //std::set<std::string> vars = it->second.getReqVars();
            //std::set<std::string>::iterator setitr;
            //for(setitr = vars.begin(); setitr != vars.end(); setitr++) {
                //if(reqs.count(*setitr) == 0) {
                    //reqVarInfo_t tmp;
                    //tmp.type = EVALUATE;
                    //tmp.value = 0.0;
                    //reqs.insert(std::pair<std::string, reqVarInfo_t>
                        //(*setitr, tmp));
                //}
            //}
        //}

        boost::shared_ptr<individual> ind =
            variator_m->popIndividualToEvaluate();
        Param_t params;
        DVarContainer_t::iterator itr;
        size_t i = 0;
        for(itr = dvars_m.begin(); itr != dvars_m.end(); itr++, i++) {
            params.insert(
                std::pair<std::string, double>
                    (boost::get<VAR_NAME>(itr->second),
                     ind->genes[i]));
        }

        size_t jid = static_cast<size_t>(ind->id);
        int pilot_rank = comms_.master_local_pid;

        // now send the request to the pilot
        MPI_Send(&jid, 1, MPI_UNSIGNED_LONG, pilot_rank, OPT_NEW_JOB_TAG, comms_.opt);

        MPI_Send_params(params, pilot_rank, comms_.opt);

        //MPI_Send_reqvars(reqs, pilot_rank, comms_.opt);

        jobmapping_m.insert(
                std::pair<size_t, boost::shared_ptr<individual> >(jid, ind));

        std::ostringstream dump;
        dump << "dispatched simulation with ID " << jid << std::endl;
        job_trace_->log(dump);
    }
}


template<
      template <class> class CO
    , template <class> class MO
>
void PisaVariator<CO, MO>::runStateMachine() {

    curState_m = (PisaState_t)readState();
    switch(curState_m) {

    // State 2 of the FSM: we read mu parents using readSelectedIndividuals()
    // and create lambda children by variation of the individuals specified
    // in sel file.
    case Variate: {

        if( (maxGenerations_m > 0 && act_gen >= maxGenerations_m) ||
            (hvol_progress_ < conv_hvol_progress_) ||
            (expected_hvol_ > 0.0 && fabs(current_hvol_ - expected_hvol_) < hvol_eps_)
          ) {
            curState_m = Stop;
            writeState(curState_m);
        } else {
            if(areNewParentsAvailable() && isFileReady(var_file)) {
                std::vector<unsigned int> parents;
                readSelectedIndividuals(parents);

                // remove non-surviving individuals
                std::set<unsigned int> survivors;
                try {
                    readArchive(survivors);
                } catch (OptPilotException &e) {
                    std::cout << "Exception while reading archive: "
                              << e.what() << std::endl;
                    return;
                }
                variator_m->population()->keepSurvivors(survivors);

                // only enqueue new individuals to pilot, the poll loop will
                // feed results back to the selector.
                //FIXME: variate works on staging set!!!
                variator_m->variate(parents);
                dispatch_forward_solves();

                clearFile(sel_file);
                clearFile(arc_file);
                curState_m = Variate;
                writeState(curState_m);
            }
        }
        break;
    }

    // State 0 of the FSM: generate initial population and write
    // information about initial population to ini file.
    case Initialize: {
        if(notInitialized_m) {
            variator_m->initial_population();
            dispatch_forward_solves();
            notInitialized_m = false;
        }

        //XXX: wait till the first alpha_m individuals are ready then start
        //     the selector
        if(finishedBuffer_m.size() >= alpha_m) {

            act_gen = 1;
            toSelectorAndCommit(alpha_m, ini_file);
            curState_m = InitializeSelector;
            writeState(curState_m);
            //FIXME: and double the population again to have all workers busy
            variator_m->initial_population();
            dispatch_forward_solves();
        }
        break;
    }

    // State 4 of the FSM: stopping state for the variator.
    case Stop: {
        std::set<unsigned int> survivors;
        try {
            readArchive(survivors);
            variator_m->population()->keepSurvivors(survivors);
            dumpPopulationToFile();
        } catch(OptPilotException &e) {
            std::cout << "Exception while reading archive: "
                      << e.what() << std::endl;
        }

        variator_m->population()->clean_population();
        curState_m = VariatorStopped;
        writeState(curState_m);

        // notify pilot that we have converged
        int dummy = 0;
        MPI_Request req;
        MPI_Isend(&dummy, 1, MPI_INT, comms_.master_local_pid,
                  MPI_OPT_CONVERGED_TAG, comms_.opt, &req);

        break;
    }

    // State 7 of the FSM: stopping state for the selector.
    case SelectorStopped: {
        curState_m = Stop;
        writeState(curState_m);
        break;
    }

    // State 8 of the FSM: reset the variator, restart in state 0.
    case Reset: {
        act_gen = 1;
        std::set<unsigned int> survivors;
        try {
            readArchive(survivors);
            variator_m->population()->keepSurvivors(survivors);
            dumpPopulationToFile();
        } catch(OptPilotException &e) {
            std::cout << "Exception while reading archive: "
                      << e.what() << std::endl;
        }

        variator_m->population()->clean_population();
        curState_m = ReadyForReset;
        writeState(curState_m);
        break;
    }

    // State 11 of the FSM: selector has just reset and is ready to
    // start again in state 1.
    case Restart: {
        curState_m = Initialize;
        writeState(curState_m);
        break;
    }

    default:
        break;

    }
}


template<
      template <class> class CO
    , template <class> class MO
>
void PisaVariator<CO, MO>::readArchive(std::set<unsigned int> &survivors) {

    std::ifstream file;
    file.open(arc_file.c_str(), std::ios::in);

    int len = 0;
    file >> len;
    if(len <= 0) {
        file.close();
        throw OptPilotException("PisaVariator::readArchive()",
                    "size<=0 in arc file! (We need to keep one individual!)");
    }

    unsigned int tmp = 0;
    for(int i = 0; i < len; i++) {
        file >> tmp;
        survivors.insert(tmp);
    }

    char tag[4];
    file >> tag;
    file.close();

    if (strcmp(tag, "END") != 0)
        throw OptPilotException("PisaVariator::readArchive()",
                                "No END tag in arc file!");
}


template< template <class> class CO, template <class> class MO >
void PisaVariator<CO, MO>::readSelectedIndividuals(
        std::vector<unsigned int> &parents) {

    int size = 0;
    std::ifstream file;
    file.open(sel_file.c_str(), std::ios::in);
    file >> size;

    for(int i = 0; i < size; i++) {
        unsigned int parentID = 0;
        file >> parentID;
        parents.push_back(parentID);
    }

    char tag[4];
    file >> tag;
    file.close();

    if (strcmp(tag, "END") != 0)
        throw OptPilotException("PisaVariator::readSelectedIndividuals()",
                                "No END tag in sel file!");
}


template< template <class> class CO, template <class> class MO >
int PisaVariator<CO, MO>::writeIndividualsForSelector(std::string ofile,
                                              std::set<unsigned int> ids) {

    typedef typename PisaVariator::Individual_t individual;

    std::ostringstream outdata;
    outdata.precision(15);
    outdata << (ids.size() * (dim_m + 1)) << std::endl;

    std::set<unsigned int>::iterator it;
    for(it = ids.begin(); it != ids.end(); it++) {
        boost::shared_ptr<individual> ind =
            variator_m->population()->get_staging(*it);
        outdata << *it << " ";
        for(size_t j = 0; j < ind->objectives.size(); j++)
            outdata << ind->objectives[j] << " ";
        outdata << std::endl;
    }

    outdata << "END";

    std::ofstream file;
    file.open(ofile.c_str(), std::ios::out);
    file.precision(15);
    file << outdata.str();
    file.flush();
    file.rdbuf()->pubsync();
    file.close();

    return 0;
}


template< template <class> class CO, template <class> class MO >
int PisaVariator<CO, MO>::writeState(int state) {

    assert((0 <= state) && (state <= 11));

    std::ofstream file;
    file.open(sta_file.c_str(), std::ios::out);
    file << state;
    file.close();
    return 0;
}


template< template <class> class CO, template <class> class MO >
int PisaVariator<CO, MO>::readState() {

    int state = -1;

    std::ifstream file;
    file.open(sta_file.c_str(), std::ios::in);
    if (file.good()) {
        file >> state;
        file.close();
        if(state == -1) return readState();
        if(state < 0 || state > 11) {
            std::cout << "invalid state read: " << state << std::endl;
            throw OptPilotException("PisaVariator::readState()",
                                    "invalid state read");
        }
    }

    return state;
}


template< template <class> class CO, template <class> class MO >
bool PisaVariator<CO, MO>::areNewParentsAvailable() {

    int control_element = 1;

    std::ifstream file;
    file.open(sel_file.c_str(), std::ios::in);
    file >> control_element;
    file.close();

    if(control_element == 0)
        return false;
    else
        return true;
}


template< template <class> class CO, template <class> class MO >
bool PisaVariator<CO, MO>::isFileReady(std::string filename) {

    int control_element = 1;

    std::ifstream file;
    file.open(filename.c_str(), std::ios::in);
    file >> control_element;
    file.close();

    if(control_element == 0)
        return true;
    else
        return false;
}


template< template <class> class CO, template <class> class MO >
void PisaVariator<CO, MO>::clearFile(std::string filename) {
    std::ofstream file;
    file.open(filename.c_str(), std::ios::out);
    file << 0;
    file.close();
}


template< template <class> class CO, template <class> class MO >
void PisaVariator<CO, MO>::dumpPopulationToFile() {

    // only dump old data format if the user requests it
    if(! args_->getArg<bool>("dump-dat", false, false)) return;

    typedef typename PisaVariator::Individual_t individual;
    boost::shared_ptr<individual> temp;

    std::ofstream file;
    std::ostringstream filename;
    filename << resultDir_m << "/" << act_gen << "_" << resultFile_m
             << "_" << comms_.island_id;
    file.open(filename.str().c_str(), std::ios::out);
    file << file_param_descr_ << std::endl;

    typename std::map<unsigned int, boost::shared_ptr<individual> >::iterator it;
    for(it = variator_m->population()->begin();
        it != variator_m->population()->end(); it++) {

        file << it->first << " ";

        temp = it->second;
        for(size_t i=0; i<temp->objectives.size(); i++)
            file << temp->objectives[i] << " ";

        for(size_t i=0; i<temp->genes.size(); i++)
            file << temp->genes[i] << " ";

        file << std::endl;
    }

    file.close();

}


template< template <class> class CO, template <class> class MO >
void PisaVariator<CO, MO>::dumpPopulationToJSON() {

    typedef typename PisaVariator::Individual_t individual;
    boost::shared_ptr<individual> temp;

    std::ofstream file;
    std::ostringstream filename;
    filename << resultDir_m << "/" << act_gen << "_" << resultFile_m
             << "_" << comms_.island_id << ".json";
    file.open(filename.str().c_str(), std::ios::out);

    file << "{" << std::endl;
    file << "\t" << "\"name\": " << "\"opt-pilot\"," << std::endl;
    
    file << "\t" << "\"dvar-bounds\": {" << std::endl;
    DVarContainer_t::iterator itr = dvars_m.begin();
    for ( bounds_t::iterator it = dVarBounds_m.begin();
          it != dVarBounds_m.end(); ++it, ++itr )
    {
         file << "\t\t\"" << boost::get<VAR_NAME>(itr->second) << "\": "
              << "[ " << it->first << ", " << it->second << " ]";
         
         if ( it != dVarBounds_m.end() - 1 )
             file << ",";
         file << "\n";
    }
    file << "\t}\n\t," << std::endl;
    
    file << "\t" << "\"constraints\": [" << std::endl;
    for ( Expressions::Named_t::iterator it = constraints_m.begin();
          it != constraints_m.end(); ++it )
    {
        file << "\t\t\"" << it->second->toString() << "\"";
        
        if ( it != std::prev(constraints_m.end(), 1) )
            file << ",";
        file << "\n";
    }
    file << "\t]\n\t," << std::endl;
    
    file << "\t" << "\"solutions\": " << "[" << std::endl;

    typename std::map<unsigned int, boost::shared_ptr<individual> >::iterator it;
    DVarContainer_t::iterator itr;
    for(it = variator_m->population()->begin();
        it != variator_m->population()->end(); it++) {

        if(it != variator_m->population()->begin())
            file << "\t" << "," << std::endl;

        file << "\t" << "{" << std::endl;
        file << "\t" << "\t" << "\"ID\": " << it->first << "," << std::endl;
    
        Expressions::Named_t::iterator expr_it;
        expr_it = objectives_m.begin();
        temp = it->second;
        
        file << "\t\t\"obj\":\n" << "\t\t{\n";
        for(size_t i=0; i < temp->objectives.size(); i++, expr_it++) {
            file << "\t" << "\t" << "\t" << "\"" << expr_it->first << "\": "
                 << temp->objectives[i];
            if( i + 1 != temp->objectives.size())
                file << ",";
            file << std::endl;
        }
        file << "\t\t}\n" << "\t\t,\n" << "\t\t\"dvar\":\n" << "\t\t{\n";
        size_t i = 0;
        for(itr = dvars_m.begin(); itr != dvars_m.end(); ++i, ++itr) {
            file << "\t\t\t\"" << boost::get<VAR_NAME>(itr->second) << "\": "
                 << temp->genes[i];
            if ( i + 1 != temp->genes.size())
                file << ",";
            file << "\n";
        }
        file << "\t\t}\n";

        file << "\t" << "}" << std::endl;
    }

    file << "\t" << "]" << std::endl;
    file << "}" << std::endl;
    file.close();
}
