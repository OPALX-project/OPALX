#ifndef __SPLIT_STRATEGY__
#define __SPLIT_STRATEGY__

#include "mpi.h"

#include "Util/CmdArguments.h"
#include "Util/OptPilotException.h"


enum commGroupColorings_t {Internal, ExternalToOpt, ExternalToWorker};

/// Defines an interface for splitter strategy implementations
class SplitStrategy {

public:

    SplitStrategy(CmdArguments_t args,
                  MPI_Comm comm = MPI_COMM_WORLD)
        : comm_(comm)
        , cmd_args_(args)
        , role_(UNASSIGNED)
    {

        MPI_Comm_rank(comm, &rank_);
        MPI_Comm_rank(MPI_COMM_WORLD, &global_rank_);
        MPI_Comm_size(comm, &num_procs_);

        group_id_ = 0;

        if(num_procs_ < 3)
            throw OptPilotException("SplitStrategy::SplitStartegy",
                                    "We need 3 or more cores to split!");
    }


    virtual ~SplitStrategy()
    {}


    /**
     *  Forcing concrete implementation to split and assign poller, optimizer
     *  and worker nodes.
     */
    virtual void split() = 0;

    MPI_Comm getComm()  const { return comm_; }

    int getRank()       const { return rank_; }
    int getGlobalRank() const { return global_rank_; }
    int getNP()         const { return num_procs_; }

    Role_t getRole()    const { return role_; }
    int getLeader()     const { return leader_; }
    int getPoller()     const { return poller_; }

    std::vector<int> getWorkers()     const { return workers_; }
    std::vector<int> getOptimizers()  const { return optimizers_; }
    std::vector<int> getCoworkers()   const { return coworkers_; }


private:

    /// communicator we are splitting
    MPI_Comm comm_;


protected:

    int rank_;
    int global_rank_;
    int num_procs_;
    int group_id_;

    CmdArguments_t cmd_args_;

    Role_t role_;

    /// defines comm splitting
    std::vector<unsigned int> colorings_;

    /// every core specifies a leader (master is its own leader)
    int leader_;

    /// every core can specifies a master
    int poller_;

    /// used in master <-> workers communicator
    std::vector<int> workers_;

    /// used in aster <-> optimizers communicator
    std::vector<int> optimizers_;

    /// every role has one or more pids to solve the task at hand
    // a worker leader has more worker-cores to solve the forward problem
    // a optimizer leader has more opt-cores to solve the opt problem
    // a master has other masters working on the same opt problem
    std::vector<int> coworkers_;
};

#endif
