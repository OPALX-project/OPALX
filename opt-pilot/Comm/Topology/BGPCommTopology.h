#ifdef ARCH_BG

#ifndef __BGP_COMM_TOPOLOGY__
#define __BGP_COMM_TOPOLOGY__

#include "mpi.h"

#include "Comm/Topology/CommTopology.h"

/**
 * Using personality structure on IBM system to extract topology information
 */
class BGPCommTopology : protected CommTopology {

public:

    BGPCommTopology(MPI_Comm comm = MPI_COMM_WORLD);
    virtual ~BGPCommTopology()
    {}

    /**
     *  Discover underlaying BG/P topology using BGP_Personality_t
     *
     *  Requires an MPI_Allreduce (inplace) to compute the total extension
     *  in each direction.
     */
    void discover();
};

#endif

#else

class BGPCommTopology;

#endif
