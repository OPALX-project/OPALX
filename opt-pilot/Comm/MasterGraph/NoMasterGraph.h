#ifndef __NO_MASTER_GRAPH__
#define __NO_MASTER_GRAPH__

#include <set>

/// A simple empty master graph (no neighbors, every island works isolated).
template < class TopoDiscoveryStrategy_t >
class NoMasterGraph : public TopoDiscoveryStrategy_t {

public:

    std::set<size_t> execute(size_t numMasters, size_t dimensions, size_t id,
                             int island_id) {
        std::set<size_t> empty;
        return empty;
    }
};

#endif
