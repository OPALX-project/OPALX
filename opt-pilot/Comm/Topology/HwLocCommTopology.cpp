#include "Comm/Topology/HwLocCommTopology.h"

#include <iostream>


void HwLocCommTopology::discover() {

    hwloc_topology_t topology;

    hwloc_topology_init(&topology);
    hwloc_topology_set_flags(topology, HWLOC_TOPOLOGY_FLAG_WHOLE_IO);
    hwloc_topology_load(topology);

    num_nodes_ = hwloc_get_nbobjs_by_depth(topology, 1);
    num_cores_ = hwloc_get_nbobjs_by_depth(topology, 5);
    num_pus_   = hwloc_get_nbobjs_by_depth(topology, 6);

    //TODO: walk parent links until locating non-I/O object to determine
    //      locality
    size_t i = 0;
    hwloc_obj_t res = hwloc_get_obj_by_type(topology, HWLOC_OBJ_OS_DEVICE, i);
    while(res) {
        printChildren(topology, res, 0);
        res = hwloc_get_obj_by_type(topology, HWLOC_OBJ_OS_DEVICE, ++i);
    }

    hwloc_topology_destroy(topology);
}


void HwLocCommTopology::printChildren(hwloc_topology_t topology,
                                      hwloc_obj_t obj, int depth) {

    char info[128];

    hwloc_obj_snprintf(info, sizeof(info), topology, obj, "#", 0);
    std::cout << 2*depth << " " << info << ": " << obj->name << std::endl;
    for(unsigned i = 0; i < obj->arity; i++)
        printChildren(topology, obj->children[i], depth + 1);

}


void HwLocCommTopology::printTopologyInformation(hwloc_topology_t topology) {

    int topodepth = hwloc_topology_get_depth(topology);
    for(int depth = 0; depth < topodepth; depth++) {
        std::cout << "Level " << depth << std::endl;
        for(size_t i = 0; i < hwloc_get_nbobjs_by_depth(topology, depth); i++) {
            char info[512];
            hwloc_obj_snprintf(info, sizeof(info), topology,
                    hwloc_get_obj_by_depth(topology, depth, i), "#", 0);

            std::cout << "\t" << info << std::endl;
        }
    }
}
