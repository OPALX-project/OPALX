#ifndef __COMM_TYPES__
#define __COMM_TYPES__

#include <vector>
#include <map>
#include <set>

#include "mpi.h"

#include "boost/tuple/tuple.hpp"
#include "boost/variant/variant.hpp"

namespace Comm {

    typedef size_t id_t;
    typedef size_t localId_t;

    // defining different groups of processors (consecutive vs. sets).. This
    // is currently not used.
    typedef boost::tuple<size_t, size_t> blockProcessorGroup_t;
    typedef std::set<size_t>             setProcessorGroup_t;

    typedef boost::variant <
            std::vector<blockProcessorGroup_t>
            , std::vector<setProcessorGroup_t>
            >
        processorGroups_t;

    /// bundles all communicators for a specific role/pid
    struct Bundle_t {
        int island_id;
        int leader_pid;
        int master_pid;
        int master_local_pid;
        MPI_Comm worker;
        MPI_Comm opt;
        MPI_Comm coworkers;
        MPI_Comm world;
    };
}

#endif
