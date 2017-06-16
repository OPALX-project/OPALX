#include <string>

#include <sys/stat.h>

#include "boost/algorithm/string.hpp"
#include "boost/foreach.hpp"
#define foreach BOOST_FOREACH

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/vector.hpp>

#include "Optimizer/HomotopyImpl.h"
#include "Util/OptPilotException.h"
#include "Util/MPIHelper.h"

#include "Util/Trace/TraceComponent.h"
#include "Util/Trace/Timestamp.h"
#include "Util/Trace/FileSink.h"

HomotopyImpl::ExchangeRequest::ExchangeRequest(std::ostringstream &os) {
    this->buf_size_req = new MPI_Request();
    this->buffer_req = new MPI_Request();
    this->buf_size = os.str().length();
    buffer = new char[buf_size];
    memcpy(buffer, os.str().c_str(), buf_size);
}

bool HomotopyImpl::ExchangeRequest::isComplete() {
    int bufferflag = 0;
    MPI_Test(this->buffer_req, &bufferflag, MPI_STATUS_IGNORE);
    if( bufferflag ) {
        int sizeflag = 0;
        MPI_Test(this->buf_size_req, &sizeflag, MPI_STATUS_IGNORE);
        if( sizeflag ) {
            return true;
        }
    }
    return false;
}

HomotopyImpl::ExchangeRequest::~ExchangeRequest() {
    delete buf_size_req;
    delete buffer_req;
    delete[] buffer;
}

HomotopyImpl::HomotopyImpl(Expressions::Named_t objectives,
                           Expressions::Named_t constraints,
                           DVarContainer_t dvars,
                           size_t dim, Comm::Bundle_t comms,
                           CmdArguments_t args)
             : Optimizer(comms.opt)
             , statistics_(new Statistics<size_t>("homotopy"))
             , comms_(comms)
             , objectives_(objectives)
             , constraints_(constraints)
             , dvars_(dvars)
             , args_(args)
             , dim_(dim)
{
    iter_ = 0;
    my_local_pid_ = 0;
    MPI_Comm_rank(comms_.opt, &my_local_pid_);

    // use some sane defaults if not specified
    resultFile_ = args->getArg<std::string>("outfile", "res", false);
    resultDir_ = args->getArg<std::string>("outdir", "generations", false);

    exchangeSolStateFreq_m = args->getArg<size_t>("sol-synch", 0, false);

    Expressions::Named_t::iterator it;
    for(it = objectives_.begin(); it != objectives_.end(); it++)
        file_param_descr_ += '%' + it->first + ',';

    DVarContainer_t::iterator itr;
    for(itr = dvars_.begin(); itr != dvars_.end(); itr++) {
        file_param_descr_ += '%' + boost::get<VAR_NAME>(itr->second) + ',';
        dVarBounds_.push_back(
                std::pair<double, double>
                    (boost::get<LOWER_BOUND>(itr->second),
                     boost::get<UPPER_BOUND>(itr->second)));
    }
    file_param_descr_ = file_param_descr_.substr(0,
                                                 file_param_descr_.size()-1);


    // Traces and statistics
    std::ostringstream trace_filename;
    trace_filename << "homotopy.trace." << comms_.island_id;
    job_trace_.reset(new Trace("Homotopy Job Trace"));
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
}

HomotopyImpl::~HomotopyImpl()
{}

void HomotopyImpl::initialize() {

    // start poll loop
    run_clock_start_ = boost::chrono::system_clock::now();
    last_clock_ = boost::chrono::system_clock::now();

    run();

    boost::chrono::duration<double> total =
        boost::chrono::system_clock::now() - run_clock_start_;
    std::ostringstream stats;
    stats << "__________________________________________" << std::endl;
    stats << "Iterations " <<  iter_ << std::endl;
    stats << "TOTAL = " << total.count() << "s" << std::endl;
    statistics_->dumpStatistics(stats);
    stats << "__________________________________________" << std::endl;
    progress_->log(stats);

}

#include "Problems.h"
#include "homotopy.hpp"
#include "mpiComms.hpp"

void HomotopyImpl::run() {

    //Get the number of points/iterations to sample and tolerance/fd_step
    size_t num_points = args_->getArg<size_t>("num-points", 3, false);
    size_t iterations = args_->getArg<size_t>("iterations", 1, false);
    double tolerance = args_->getArg<double>("tolerance", 1e-3, false);
    double fd_step = args_->getArg<double>("fd-step", 1e-6, false);
    int preproject = args_->getArg<int>("preprojection", 0, false);

    //Declare an abstract problem
    //TODO: Incorporate propper constraint handling
    Problem::Interface P;
    P.dimDesign = this->dvars_.size();
    P.dimObj = this->objectives_.size();
    DVarContainer_t::iterator i;
    for(i=this->dvars_.begin(); i!=this->dvars_.end(); i++){
        //Unpack the lower/upper bounds (stored as boost::tuples)
        P.lowerBounds.push_back( i->second.get<1>() );
        P.upperBounds.push_back( i->second.get<2>() );
    }    

    //Prepare a comms handler
    Homotopy::Communication::AdHoc< MPIimpl > Comm;
    Comm.comm_.AddComm( &(comms_.opt) );
    Comm.comm_.AddComm( &(comms_.world) );
    Comm.Dispatcher = boost::bind( &HomotopyImpl::dispatch_forward_solves, this, _1, _2 );
    Comm.Handler = boost::bind( &HomotopyImpl::onMessage, this, _1, _2, _3, _4 );
    Comm.Exchanger = boost::bind( &HomotopyImpl::exchangeSolutionStates, this, _1 );
    Comm.initialize();

    //Instantiate the homotopy object
    Pareto::homotopy h( &P, tolerance, fd_step, Comm );
    h.UsePreProjection = preproject?true:false;

    //HACK: to figure out where the optimizers are
    int id; int size;
    MPI_Comm_rank(comms_.world, &id);
    MPI_Comm_size (comms_.world, &size);

    //Sample the front
    h.GetFront(num_points, iterations, id, size);

    //Shutdown the processors
    this->SendStop();
    Comm.shutdown();
}

// notify pilot that we have converged
void HomotopyImpl::SendStop() {

    int dummy = 0;
    //MPI_Request req;
    //MPI_Isend(&dummy, 1, MPI_INT, comms_.master_local_pid, OPT_CONVERGED_TAG, comms_.opt, &req);
    MPI_Send(&dummy, 1, MPI_INT, comms_.master_local_pid, OPT_CONVERGED_TAG, comms_.opt);

}

bool HomotopyImpl::onMessage(MPI_Status status, size_t length) {
  return true;
}

bool HomotopyImpl::onMessage(MPI_Status status, size_t length, std::queue< std::pair< size_t, objVars_t> >& q, std::queue< typename Homotopy::nodeEnvelope_t > &g) {

    MPITag_t tag = MPITag_t(status.MPI_TAG);
    switch(tag) {

    //case EXCHANGE_SOL_STATE_RES_SIZE_TAG: {
    case EXCHANGE_SOL_STATE_TAG: {

        size_t buf_size = length;
        size_t pilot_rank = status.MPI_SOURCE;

        std::ostringstream dump;
        dump << "new results from other cores " << buf_size << std::endl;
        job_trace_->log(dump);

        char *buffer = new char[buf_size];
        MPI_Status stat;
        MPI_Recv(buffer, buf_size, MPI_CHAR, pilot_rank,
                 MPI_EXCHANGE_SOL_STATE_DATA_TAG, comms_.world, &stat);

        //Diagnostics
        //for(int i=0; i<buf_size; i++) std::cout << std::hex << buffer[i];
        //std::cout << std::endl;

        dump.clear();
        dump.str(std::string());
        dump << "got results from other cores " << buf_size << std::endl;
        job_trace_->log(dump);

        SolutionState_t new_states;
        try {
            std::istringstream is(buffer);
            boost::archive::text_iarchive ia(is);
            ia >> new_states;
            
        } catch ( boost::archive::archive_exception &e ) {
            dump << "boost::archive::archive_exception " << e.what() << std::endl;
            unsigned i;
            for(i=0; i<buf_size; i++) std::cout << std::hex << buffer[i] << " " ; 
                std::cout << std::endl;
        }
        delete[] buffer;

        foreach(solution_t ind, new_states) {

            //handle the received solution
            g.push( ind );

            dump.clear();
            dump.str(std::string());
            dump << "solution (ID: " << ind.first
                 << ") successfully migrated from another island" << std::endl;
            job_trace_->log(dump);
        }

        return true;
    }

    case REQUEST_FINISHED: {

        size_t jid = length;
        std::map<size_t, boost::shared_ptr<solution_t> >::iterator it;

        std::ostringstream dump;
        dump << "job with ID " << jid << " delivered results" << std::endl;

        reqVarContainer_t res;
        MPI_Recv_reqvars(res, status.MPI_SOURCE, comms_.opt);

        objVars_t recieved_solution;

        reqVarContainer_t::iterator itr;
        std::map<std::string, double> vars;
        for(itr = res.begin(); itr != res.end(); itr++) {
            //FIXME: proper constraint handling
            if(itr->second.is_valid) {
                if(itr->second.value.size() == 1){
                    recieved_solution.push_back(itr->second.value[0]);
			dump << "Received value: " << itr->second.value[0] << std::endl;
		}
		else{
			dump << "Received value: " << itr->second.value[0] << " with size != 1 " << std::endl;

		}

            } else {
	        dump << "job with ID " << jid << " has invalid results" << std::endl;
		dump << "Design Vars: " ;
		dump.precision( 15 );
		for(unsigned vals=0; vals<this->debugVals[jid].size(); vals++) dump << this->debugVals[jid][vals] << " ";
		dump << std::endl;
            }
        }

        //add results to outgoing queue
        if( recieved_solution.size() > 0 ) q.push( std::pair< size_t, objVars_t > ( jid, recieved_solution ) );

        job_trace_->log(dump);

        //free resources associated with request
        if ( reqs.find(jid) != reqs.end() ){
            int flag = 0;
            MPI_Test(reqs[jid].first, &flag, MPI_STATUS_IGNORE);
            if( flag ) {
                //erase the jobID information
                delete  reqs_ids[jid].second;
                delete  reqs_ids[jid].first;
                reqs_ids.erase( jid );

                //erase the message data buffer
                delete [] reqs[jid].second.second;
                delete reqs[jid].second.first;
                delete reqs[jid].first;
                reqs.erase( jid );
            }
        }

        return true;
    }

    default: {
        std::cout << "(Homotopy) Error: unexpected MPI_TAG: "
                  << status.MPI_TAG << std::endl;
        return false;
    }
    }
}


void HomotopyImpl::postPoll() {

    iter_++;
    boost::chrono::duration<double> total =
        boost::chrono::system_clock::now() - run_clock_start_;
    boost::chrono::duration<double> dt    =
        boost::chrono::system_clock::now() - last_clock_;
    last_clock_ = boost::chrono::system_clock::now();
    std::ostringstream stats;
    stats << "__________________________________________" << std::endl;
    stats << "Arriving at iteration " << iter_            << std::endl;
    stats << "dt = " << dt.count() << "s, total = " << total.count()
        << "s" << std::endl;
    stats << "__________________________________________" << std::endl;
    progress_->log(stats);

    // dump parents of current generation for visualization purposes
    dumpStateToFile();

    exchangeSolutionStates();
}

void HomotopyImpl::exchangeSolutionStates() {}
void HomotopyImpl::exchangeSolutionStates( std::queue< typename Homotopy::NeighborMessage_t > & toSend ){

    size_t num_masters = args_->getArg<size_t>("num-masters", 1, false);

    if(num_masters <= 1 ) return;
    
    //Pop completed sends
    std::vector< ExchangeRequest * >::iterator er = ExchangeState.begin();
    while( er!=ExchangeState.end() ){
        if( (*er)->isComplete() ) {
            delete *er;
            er = ExchangeState.erase( er );
        }
        else er++;
    }

    std::vector< typename Homotopy::nodeEnvelope_t > sent = toSend.front().second;

    while( ! toSend.empty() ){
        std::ostringstream os;
        boost::archive::text_oarchive oa(os);

        //ATTN: this address is the Simplex Subset,
        //which may be mapped non-trivially
        size_t addr = (size_t) toSend.front().first;
        SolutionState_t solutions(toSend.front().second);
        oa << solutions;
        toSend.pop();

        //Get ExchangeRequest object
        ExchangeRequest *req = new ExchangeRequest(os);

        std::ostringstream dump;
        dump << "sending my buffer size " << req->buf_size << " bytes to PILOT"
             << std::endl;
        job_trace_->log(dump);

        //MPI_Send(req.buf_size, 1, MPI_LONG, addr,
        //         EXCHANGE_SOL_STATE_TAG, comms_.opt);
        MPI_Isend(&(req->buf_size), 1, MPI_LONG, addr,
                  EXCHANGE_SOL_STATE_TAG, comms_.world, req->buf_size_req);

        MPI_Isend(req->buffer, req->buf_size, MPI_CHAR, addr,
                  MPI_EXCHANGE_SOL_STATE_DATA_TAG, comms_.world, req->buffer_req);

        //printf("Send (size: %d) : %s\n", *req.buf_size, req.buffer);
        //printf("Send (size: %d) : %s\n", buf_size, buffer);
        this->ExchangeState.push_back( req );

        dump.clear();
        dump.str(std::string());
        dump << "Sent " << req->buf_size << " bytes to PILOT" << std::endl;
        job_trace_->log(dump);
    }
}

void HomotopyImpl::dispatch_forward_solves() {}
void HomotopyImpl::dispatch_forward_solves( size_t jid, const Homotopy::designVars_t &vars) {

    Param_t params;
    DVarContainer_t::iterator itr;

	//DEBUG
	this->debugVals[jid] = vars;
	//DEBUG

    //FIXME: fill data in Param_t structure
    //cout << "Filling param_t struct for ID " << jid << endl;
    int i=0;
    for(itr = dvars_.begin(); itr != dvars_.end(); itr++, i++) {
        params.insert(
                std::pair<std::string, double>
                (boost::get<VAR_NAME>(itr->second),
                 vars[i]));
    }

    //size_t jid = 0;
    int pilot_rank = comms_.master_local_pid;

    // now send the request to the pilot
    //cout << "Sending pilot ID " << jid << endl;
    //MPI_Send(&jid, 1, MPI_LONG, pilot_rank, OPT_NEW_JOB_TAG, comms_.opt);
    size_t *JID = new size_t();
    *JID = jid;
    MPI_Request * req_id = new MPI_Request();
    MPI_Isend(JID, 1, MPI_LONG, pilot_rank, OPT_NEW_JOB_TAG, comms_.opt, req_id);
    std::pair< MPI_Request*, size_t*> pid(req_id, JID);
    reqs_ids[jid] = pid;
    //cout << "Sending pilot params for ID " << jid << endl;
    //MPI_Send_params(params, pilot_rank, comms_.opt);

    //generate new request object
    MPI_Request * req = new MPI_Request();
    std::pair<size_t*, char*> data = MPI_ISend_params(params, pilot_rank, comms_.opt, req);
    std::pair< MPI_Request*, std::pair<size_t*, char*> > p(req, data);
    reqs[jid] = p;

    //cout << "Sent pilot params for ID " << jid << endl;

    //jobmapping_.insert(
            //std::pair<size_t, boost::shared_ptr<solution_t> >(jid, ind));

    std::ostringstream dump;
    dump << "dispatched simulation with ID " << jid << endl;
    job_trace_->log(dump);
}


void HomotopyImpl::dumpStateToFile() {

    std::ofstream file;
    std::ostringstream filename;
    filename << resultDir_ << "/" << iter_ << "_" << resultFile_
             << "_" << comms_.island_id;
    file.open(filename.str().c_str(), std::ios::out);
    file << file_param_descr_ << std::endl;

    //FIXME: dump current solution state (all sample points) to file

    file.close();
}

