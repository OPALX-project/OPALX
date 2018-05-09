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

#include "Sample/Uniform.h"

// #include <memory>

// #include "boost/algorithm/string.hpp"
// #include "boost/foreach.hpp"
// #define foreach BOOST_FOREACH

// #include <boost/archive/text_oarchive.hpp>
// #include <boost/archive/text_iarchive.hpp>
// #include <boost/serialization/map.hpp>

#include "Util/OptPilotException.h"
#include "Util/MPIHelper.h"

// #include "Util/Trace/TraceComponent.h"
// #include "Util/Trace/Timestamp.h"
// #include "Util/Trace/FileSink.h"

Sampler::Sampler(Expressions::Named_t objectives,
                 Expressions::Named_t constraints,
                 DVarContainer_t dvars,
                 size_t dim, Comm::Bundle_t comms,
                 CmdArguments_t args)
    : Optimizer(comms.opt)
{
    throw OptPilotException("Sampler::Sampler",
                            "We shouldn't get here!");
}


// template< template <class> class SO >
Sampler/*<SO>*/::Sampler(Expressions::Named_t type,
                   DVarContainer_t dvars, Comm::Bundle_t comms,
                   CmdArguments_t args)
    : Optimizer(comms.opt)
    , comms_(comms)
    , type_m(type)
    , dvars_m(dvars)
    , args_(args)
{
    my_local_pid_ = 0;
    MPI_Comm_rank(comms_.opt, &my_local_pid_);
    
    //FIXME: proper rand gen initialization (use boost?!)
    srand(time(NULL) + comms_.island_id);

    // create output directory if it does not exists
    struct stat dirInfo;
    if(stat(resultDir_m.c_str(),&dirInfo) != 0)
        mkdir((const char*)(resultDir_m.c_str()), 0777);
}


/*template< template <class> class SO >*/
Sampler/*<SO>*/::~Sampler()
{}


/*template< template <class> class SO >*/
void Sampler/*<SO>*/::initialize() {
    
    nsamples_m = args_->getArg<int>("nsamples", true);
    act_sample_m = 0;
    done_sample_m = 0;
    curState_m = SUBMIT;
    gid = 0;
    
    this->initSamplingMethods_m();
    
    
    // start poll loop
    run();
}


void Sampler::initSamplingMethods_m() {
    for ( Expressions::Named_t::iterator it = type_m.begin();
          it != type_m.end(); ++it )
    {
        std::string s = it->second->toString();
        
        if ( s == "UNIFORM_INT" ) {
            samplingOp_m.push_back( std::unique_ptr< Uniform<int> >(new Uniform<int>(42, 0, 1) ) );  //FIXME        // C++14 --> std::make_unique
        }
        else if ( s == "UNIFORM_DOUBLE" ) {
            samplingOp_m.push_back( std::unique_ptr<Uniform<double> >(new Uniform<double>(42, 0.0, 1.0) ) );
        } else {
            throw OptPilotException("Sampler::initSamplingMethods_m",
                                    "Sampling method '" + s + "' not supported.");
        }
    }
}


/*template< template <class> class SO >*/
bool Sampler/*<SO>*/::onMessage(MPI_Status status, size_t length) {
    
    typedef typename Sampler::Individual_t individual;

    MPITag_t tag = MPITag_t(status.MPI_TAG);
    switch(tag) {
        case REQUEST_FINISHED: {
            unsigned int jid = static_cast<unsigned int>(length);
            typename std::map<size_t, boost::shared_ptr<individual> >::iterator it;
            it = jobmapping_m.find(jid);
    
            if(it == jobmapping_m.end()) {
                std::cout << "NON-EXISTING JOB with ID = " << jid << std::endl;
                throw OptPilotException("Sampler::onMessage",
                        "non-existing job");
            }
    
    
            boost::shared_ptr<individual> ind = it->second;
            
            jobmapping_m.erase(it);
            
            done_sample_m++;
            
            return true;
        }
        default: {
            std::cout << "(Sampler) Error: unexpected MPI_TAG: "
                      << status.MPI_TAG << std::endl;
            return false;
        }
    }
}


/*template< template <class> class SO >*/
void Sampler/*<SO>*/::postPoll() {
    
    if ( act_sample_m < nsamples_m ) {
        this->createNewIndividual_m();
    }
    
    runStateMachine();
}


void Sampler::createNewIndividual_m() {
    
    std::cout << "createNewIndividual_m" << std::endl;
    
    std::vector<std::string> dNames;
    
    DVarContainer_t::iterator itr;
    for(itr = dvars_m.begin(); itr != dvars_m.end(); itr++) {
        std::string dName = boost::get<VAR_NAME>(itr->second);
        dNames.push_back(dName);
    }
    
    boost::shared_ptr<Individual_t> ind = boost::shared_ptr<Individual_t>( new Individual_t(dNames));
    
    for (uint i = 0; i < samplingOp_m.size(); ++i) {
        
        samplingOp_m[i]->create(ind, i);
        
    }
    
    ind->id = gid++;
    
    individuals_m.push(ind);
}


/*template< template <class> class SO >*/
void Sampler/*<SO>*/::dumpPopulationToJSON() {
    
    std::cout << "dumpPopulationToJSON" << std::endl;
}

void Sampler::runStateMachine() {
    
    std::cout << "runStateMachine" << std::endl;
    
    switch(curState_m) {
        
        case SUBMIT: {
            if ( done_sample_m == nsamples_m) {
                curState_m = STOP;
            } else {
                if ( act_sample_m != nsamples_m ) {
                    std::cout << "SUBMIT" << std::endl;
                    dispatch_forward_solves();
                }
            }
            break;
        }
        case STOP: {
            
            curState_m = TERMINATE;
            
            // notify pilot that we have converged
            int dummy = 0;
            MPI_Request req;
            MPI_Isend(&dummy, 1, MPI_INT, comms_.master_local_pid,
                      MPI_OPT_CONVERGED_TAG, comms_.opt, &req);
    
            break;
        }
        
        case TERMINATE: {
            break;
        }
    }
}




void Sampler::dispatch_forward_solves() {
    
    std::cout << "dispatch_forward_solves" << std::endl;

    typedef typename Sampler::Individual_t individual;

    while ( !individuals_m.empty() ) {
        boost::shared_ptr<Individual_t> ind = individuals_m.front();
        
        individuals_m.pop();
        
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
        
        
        act_sample_m++;
        
        // now send the request to the pilot
        MPI_Send(&jid, 1, MPI_UNSIGNED_LONG, pilot_rank, OPT_NEW_JOB_TAG, comms_.opt);

        MPI_Send_params(params, pilot_rank, comms_.opt);
        
        jobmapping_m.insert(
                std::pair<size_t, boost::shared_ptr<individual> >(jid, ind));
    }
}