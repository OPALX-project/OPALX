#ifndef __RANDOM_STRATEGY_H__
#define __RANDOM_STRATEGY_H__

#include <set>

#include "Comm/Splitter/SplitStrategy.h"

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

/**
 * Randomly select num_masters master cores.
 */
class RandomStrategy : protected SplitStrategy {

public:

    RandomStrategy(size_t num_masters, boost::shared_ptr<CommTopology> topology,
                   MPI_Comm comm = MPI_COMM_WORLD)
        : SplitStrategy(num_masters, topology, comm) {

        fillPids(num_procs_);
    }


    virtual ~RandomStrategy()
    {}


    void split() {
        // we start by randomly picking processor IDs for the specified
        // number of masters
        std::vector<int> master_pids;
        for(int i=0; i < num_masters_; i++)
            master_pids.push_back(nextRandomPid());

        int num_workers = pids_.size();
        int num_assigned_masters = 0;

        // assign workers uniform at random to each master
        for(int master_pid : master_pids) {

            if(master_pid == rank_) role_ = Role_t.POLLER;

            // compute number of workers for this master
            int fraction_workers = num_workers / num_masters_ ;
            if(num_assigned_masters < num_workers % num_masters_)
                fraction_workers++;

            std::vector<int> worker_pids;

            for(int j = 0; j < fraction_workers; j++) {
                int worker_pid = nextRandomPid();
                if(worker_pid == rank_) role_ = Role_t.WORKER;
                worker_pids.push_back(worker_pid);
            }

            if(role_ == Role_t.POLLER || role_ == Role_t.WORKER)
                addCommGroup(master_pid, worker_pids);

            num_assigned_masters++;
        }
    }


private:

    std::set<int> pids_;
    boost::random::mt19937 rng_;


    void fillPids(int num) {
        for(int i=0; i < num; i++)
            pids_.insert(i);
    }


    int nextRandomPid() {
        boost::random::uniform_int_distribution<> random_pid(0,pids_.size()-1);
        int pid_idx= random_pid(rng_);

        int pid = 0;
        std::set<int>::iterator itr;

        for(itr = pids_.begin(); itr != pids_.end(); itr++, pid_idx--) {
            if(pid_idx == 0) {
                pid = *itr;
                pids_.erase(itr);
                break;
            }
        }

        return pid;
    }

};

#endif
