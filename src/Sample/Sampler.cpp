#include <iostream>
#include <string>

#include "Sample/Sampler.h"

#include "Util/OptPilotException.h"
#include "Util/MPIHelper.h"

#include <boost/property_tree/json_parser.hpp>

#include <boost/filesystem.hpp>

Sampler::Sampler(Expressions::Named_t objectives,
                 Expressions::Named_t constraints,
                 DVarContainer_t dvars,
                 size_t dim, Comm::Bundle_t comms,
                 CmdArguments_t args,
                 std::vector<double> hypervolRef)
    : Optimizer(comms.opt)
{
    throw OptPilotException("Sampler::Sampler",
                            "We shouldn't get here!");
}


Sampler::Sampler(const std::map<std::string,
                                std::shared_ptr<SamplingMethod>
                    >& sampleMethods,
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

    resultFile_m = args->getArg<std::string>("outfile", "samples.json", false);
    resultDir_m = args->getArg<std::string>("outdir", "samples", false);

    if ( !boost::filesystem::exists(resultDir_m) ) {
        boost::filesystem::create_directory(resultDir_m);
    }

    DVarContainer_t::iterator itr;
    for(itr = dvars_m.begin(); itr != dvars_m.end(); itr++) {
        dVarBounds_m.push_back(
                std::pair<double, double>
                    (boost::get<LOWER_BOUND>(itr->second),
                     boost::get<UPPER_BOUND>(itr->second)));
    }
}


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
    MPITag_t tag = MPITag_t(status.MPI_TAG);
    switch(tag) {
        case REQUEST_FINISHED: {
            unsigned int jid = static_cast<unsigned int>(length);
            typename std::map<size_t, boost::shared_ptr<Individual_t> >::iterator it;
            it = jobmapping_m.find(jid);

            if(it == jobmapping_m.end()) {
                std::cout << "NON-EXISTING JOB with ID = " << jid << std::endl;
                throw OptPilotException("Sampler::onMessage",
                        "non-existing job");
            }


            boost::shared_ptr<Individual_t> ind = it->second;

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

    for(itr = dvars_m.begin(); itr != dvars_m.end(); itr++) {
        std::string dName = boost::get<VAR_NAME>(itr->second);
        int i = ind->getIndex(dName);
        sampleMethods_m[dName]->create(ind, i);
    }

    // FIXME does not work with more than 1 master
    ind->id = gid++;

    addIndividualToJSON_m(ind);

    individuals_m.push(ind);
}

void Sampler::dumpIndividualsToJSON_m() {

    boost::property_tree::ptree tree;

    tree.put("name", "sampler");

    std::stringstream bounds;
    DVarContainer_t::iterator itr = dvars_m.begin();
    for (bounds_t::iterator it = dVarBounds_m.begin();
         it != dVarBounds_m.end(); ++it, ++itr)
    {
        std::string dvar = boost::get<VAR_NAME>(itr->second);
        bounds << "[ " << it->first << ", " << it->second << " ]";
        tree.put("dvar-bounds." + dvar, bounds.str());
        bounds.str("");
    }

    tree.add_child("samples", samples_m);


    std::ostringstream filename;
    filename << resultDir_m << "/" << resultFile_m
             << "_samples.json";

    boost::property_tree::write_json(filename.str(), tree);
}

void Sampler::addIndividualToJSON_m(const boost::shared_ptr<Individual_t>& ind) {
    boost::property_tree::ptree sample;

    sample.put("ID", ind->id);

    DVarContainer_t::iterator itr;
    for(itr = dvars_m.begin(); itr != dvars_m.end(); itr++) {
        std::string name = boost::get<VAR_NAME>(itr->second);
        int i = ind->getIndex(name);
        sample.put("dvar." + name, ind->genes[i]);
    }

    samples_m.push_back(std::make_pair("", sample));
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

            dumpIndividualsToJSON_m();

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

    while ( !individuals_m.empty() ) {
        boost::shared_ptr<Individual_t> ind = individuals_m.front();

        individuals_m.pop();

        Param_t params;
        DVarContainer_t::iterator itr;

        for(itr = dvars_m.begin(); itr != dvars_m.end(); itr++) {
            std::string dName = boost::get<VAR_NAME>(itr->second);
            int i = ind->getIndex(dName);
            params.insert(
                std::pair<std::string, double>
                    (dName, ind->genes[i]));
        }

        size_t jid = static_cast<size_t>(ind->id);
        int pilot_rank = comms_.master_local_pid;


        act_sample_m++;

        // now send the request to the pilot
        MPI_Send(&jid, 1, MPI_UNSIGNED_LONG, pilot_rank, OPT_NEW_JOB_TAG, comms_.opt);

        MPI_Send_params(params, pilot_rank, comms_.opt);

        jobmapping_m.insert(
                std::pair<size_t, boost::shared_ptr<Individual_t> >(jid, ind));
    }
}