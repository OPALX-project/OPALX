#ifdef ARCH_BG

#include "Comm/Topology/BGQCommTopology.h"

//#include <spi/kernel_interface.h>
//#include <common/bgp_personality.h>
//#include <common/bgp_personality_inlines.h>

#include <spi/include/kernel/location.h>


BGQCommTopology::BGQCommTopology(MPI_Comm comm = MPI_COMM_WORLD)
    : CommTopology(comm) {
        num_dims_ = 5;
        coords_.resize(num_dims_);
        dims_.resize(num_dims_);
        torus_links_.resize(num_dims_);
    }


void BGQCommTopology::discover() {

    proc_id_     = Kernel_ProcessorID();       // 0-63
    my_core_id_  = Kernel_ProcessorCoreID();   // 0-15
    hwID_        = Kernel_ProcessorThreadID(); // 0-3

    Personality_t pers;
    Kernel_GetPersonality(&pers, sizeof(pers));

    dims_[0] = pers.Network_Config.Anodes;
    dims_[1] = pers.Network_Config.Bnodes;
    dims_[2] = pers.Network_Config.Cnodes;
    dims_[3] = pers.Network_Config.Dnodes;
    dims_[4] = pers.Network_Config.Enodes;

    coords_[0] = pers.Network_Config.Acoord;
    coords_[1] = pers.Network_Config.Bcoord;
    coords_[2] = pers.Network_Config.Ccoord;
    coords_[3] = pers.Network_Config.Dcoord;
    coords_[4] = pers.Network_Config.Ecoord;

    uint64_t Nflags = pers.Network_Config.NetFlags;
    for(int i=0; i < num_dims_; i++)
        torus_links_[i] = 0;

    if (Nflags & ND_ENABLE_TORUS_DIM_A) torus_links_[0] = 1;
    if (Nflags & ND_ENABLE_TORUS_DIM_B) torus_links_[1] = 1;
    if (Nflags & ND_ENABLE_TORUS_DIM_C) torus_links_[2] = 1;
    if (Nflags & ND_ENABLE_TORUS_DIM_D) torus_links_[3] = 1;
    if (Nflags & ND_ENABLE_TORUS_DIM_E) torus_links_[4] = 1;

    if (getRank() == 0) {
        std::cout << "block shape : < ";
        for(int i=0; i < num_dims_; i++)
            std::cout << dims_[i] << " ";
        std::cout << std::endl;

        std::cout << "torus links enabled : < ";
        for(int i=0; i < num_dims_; i++)
            std::cout << torus_links_[i] << " ";
        std::cout << std::endl;
    }

    std::cout << "rank " << core << " location < ";
    for(int i=0; i < num_dims_; i++)
        std::cout << coords_[i] << " ";

    std::cout << " > core " << my_core_id_ << " hwthread " << hwID_
              << " procid " << proc_id_ << std::endl;
}

#endif
