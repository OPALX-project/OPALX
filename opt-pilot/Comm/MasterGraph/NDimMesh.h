#ifndef __NDIM_MESH__
#define __NDIM_MESH__

#include <set>

//FIXME:
//#include "Mesh.hpp"

template < class TopoDiscoveryStrategy_t >
class NDimMesh : public TopoDiscoveryStrategy_t {

public:

    std::set<size_t> execute(size_t numMasters, size_t dimensions, size_t id,
                             int island_id) {
        return Mesh::Simplex::getNeighborIDs(island_id);
    }
};

#endif
