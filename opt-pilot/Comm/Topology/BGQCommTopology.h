#ifdef ARCH_BG

#ifndef __BGQ_COMM_TOPOLOGY__
#define __BGQ_COMM_TOPOLOGY__

#include "mpi.h"
#include <vector>

#include "Comm/Topology/CommTopology.h"

/**
 * Using personality structure on IBM system to extract topology information
 */
class BGQCommTopology : protected CommTopology {

public:

    BGQCommTopology(MPI_Comm comm = MPI_COMM_WORLD);
    virtual ~BGQCommTopology()
    {}

    /// Discover underlaying BG/Q topology using Personality_t
    void discover();

    std::vector<unsigned int> getTorusLinks() const { return torus_links_; }
    unsigned int getProcID() const { return proc_id_; }

private:

    unsigned int proc_id_;
    std::vector<unsigned int> torus_links_;
};

#endif

#else

class BGQCommTopology;

#endif
