#include <iostream>
#include <sstream>
#include <cassert>
#include <cmath>
#include <cstring>
#include <set>
#include <cstdlib>
#include <string>
#include <limits>

#include "Sample/Uniform.h"


#include "Util/OptPilotException.h"
#include "Util/MPIHelper.h"

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


Sampler::Sampler(const std::vector< std::shared_ptr<SamplingMethod> >& sampleMethods,
                 DVarContainer_t dvars, Comm::Bundle_t comms,
                 CmdArguments_t args)
    : Optimizer(comms.opt)
    , sampleMethods_m(sampleMethods)
    , comms_(comms)
    , dvars_m(dvars)
    , args_(args)
{
    my_local_pid_ = 0;
    MPI_Comm_rank(comms_.opt, &my_local_pid_);
}


Sampler::~Sampler()
{}


void Sampler::initialize() {
    
    nSamples_m = args_->getArg<int>("nsamples", true);
    act_sample_m = 0;
    done_sample_m = 0;
    curState_m = SUBMIT;
    
    int nMasters = args_->getArg<int>("num-masters", true);
    
    if ( nMasters > 1 )
        throw OptPilotException("Sampler::initialize",
                                "Only single master execution currently supported.");
        
    
    // unique job id, FIXME does not work with more than 1 sampler
    gid = 0;
    
    // start poll loop
    run();
}


bool Sampler::onMessage(MPI_Status status, size_t length) {
    
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


void Sampler::postPoll() {
    
    if ( act_sample_m < nSamples_m ) {
        this->createNewIndividual_m();
    }
    
    runStateMachine();
}


void Sampler::createNewIndividual_m() {
    
    std::vector<std::string> dNames;
    
    DVarContainer_t::iterator itr;
    for(itr = dvars_m.begin(); itr != dvars_m.end(); itr++) {
        std::string dName = boost::get<VAR_NAME>(itr->second);
        dNames.push_back(dName);
    }
    
    boost::shared_ptr<Individual_t> ind = boost::shared_ptr<Individual_t>( new Individual_t(dNames));
    
    for (uint i = 0; i < sampleMethods_m.size(); ++i) {
        sampleMethods_m[i]->create(ind, i);
        
        std::cout << dNames[i] << " " << ind->genes[i] << std::endl;
    }
    
    // FIXME does not work with more than 1 master
    ind->id = gid++;
    
    individuals_m.push(ind);
}


void Sampler::dumpPopulationToJSON() {
    
    std::cout << "dumpPopulationToJSON" << std::endl;
}


void Sampler::runStateMachine() {
    
    switch(curState_m) {
        
        case SUBMIT: {
            if ( done_sample_m == nSamples_m) {
                curState_m = STOP;
            } else {
                if ( act_sample_m != nSamples_m ) {
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