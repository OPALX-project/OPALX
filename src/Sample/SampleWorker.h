#ifndef __SAMPLE_WORKER_H__
#define __SAMPLE_WORKER_H__

#include "Pilot/Worker.h"


/**
 *  \class SampleWorker
 *  \brief A worker MPI entity consists of a processor group that runs a
 *  simulation of type Sim_t. The main loop in run() accepts new jobs from the
 *  master process runs the simulation and reports back the results.
 *
 *  @see SamplePilot
 *  @see Worker
 *  @see MPIHelper.h
 *
 *  @tparam Sim_T type of simulation to run
 */
template <class Sim_t>
class SampleWorker : protected Worker<Sim_t> {

public:

    SampleWorker(Expressions::Named_t constraints,
           std::string simName, Comm::Bundle_t comms, CmdArguments_t args)
        : Worker<Sim_t>(constraints, simName, comms, args)
    {
        // simulation pointer requires at least 1 objective --> provide dummy
        this->objectives_ = {
            {"dummy", new Expressions::Expr_t("dummy") }
        };
        
        int my_local_pid = 0;
        MPI_Comm_rank(this->coworker_comm_, &my_local_pid);
        
        // distinction between leader and coworkers
        if(my_local_pid == this->leader_pid_)
            this->run();
        else
            this->runCoWorker();
    }

    ~SampleWorker()
    {}

protected:
    
    bool onMessage(MPI_Status status, size_t recv_value) override {
        
        if(status.MPI_TAG == MPI_WORK_JOBID_TAG) {

            this->is_idle_ = false;
            size_t job_id = recv_value;

            // get new job
            Param_t params;
            MPI_Recv_params(params, (size_t)this->pilot_rank_, this->comm_m);

            // and forward to coworkers (if any)
            if(this->num_coworkers_ > 1) {
                this->notifyCoWorkers(MPI_COWORKER_NEW_JOB_TAG);
                MPI_Bcast_params(params, this->leader_pid_, this->coworker_comm_);
            }

            try {
                typename Worker<Sim_t>::SimPtr_t sim(new Sim_t(this->objectives_,
                                       this->constraints_,
                                       params,
                                       this->simulation_name_,
                                       this->coworker_comm_,
                                       this->cmd_args_));
                
                sim->setFilename(job_id);
                
                // run simulation in a "blocking" fashion
                sim->run();
                
            } catch(OptPilotException &ex) {
                std::cout << "Exception while running simulation: "
                          << ex.what() << std::endl;
            }

            MPI_Send(&job_id, 1, MPI_UNSIGNED_LONG, this->pilot_rank_,
                     MPI_WORKER_FINISHED_TAG, this->comm_m);

            size_t dummy = 0;
            MPI_Recv(&dummy, 1, MPI_UNSIGNED_LONG, this->pilot_rank_,
                     MPI_WORKER_FINISHED_ACK_TAG, this->comm_m, &status);
            
            this->is_idle_ = true;
            return true;

        } else {
            std::stringstream os;
            os << "Unexpected MPI_TAG: " << status.MPI_TAG;
            std::cout << "(Worker) Error: " << os.str() << std::endl;
            throw OptPilotException("Worker::onMessage", os.str());
        }
    }
};

#endif
