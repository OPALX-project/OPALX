#ifdef ARCH_BG

#include "Comm/Topology/BGPCommTopology.h"

#include <spi/kernel_interface.h>
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>


BGPCommTopology::BGPCommTopology(MPI_Comm comm = MPI_COMM_WORLD)
    : CommTopology(comm) {
        num_dims_ = 3;
        coords_.resize(num_dims_);
        dims_.resize(num_dims_);
    }


void BGPCommTopology::discover() {

    _BGP_Personality_t personality;
    Kernel_GetPersonality(&personality, sizeof(personality));

    coords_[0] = personality.Network_Config.Xcoord;
    coords_[1] = personality.Network_Config.Ycoord;
    coords_[2] = personality.Network_Config.Zcoord;
    hwID_ = Kernel_PhysicalProcessorID();

    //node_config = personality.Kernel_Config.ProcessConfig;

    char location[512];
    BGP_Personality_getLocationString(&personality, location);

    //std::cout << "MPI rank " << rank_ << " has torus coords " << coords_[0];
    //std::cout << ", " << coords_[1] << ", " << coords_[2] << " cpu = " << hwID_;
    //std::cout << " location = " << location << std::endl;

    my_core_ = hwID_;

    for(size_t i=0; i < num_dims_; i++)
        dims_[i] = coords_[i] + 1;

    MPI_Allreduce(MPI_IN_PLACE, &dims_.front(), num_dims_,
                  MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    if(rank_ == 0)
        std::cout << "dim = " << dims_[0] << ", " << dims_[1] << ", "
                  << dims_[2] << std::endl;
}

#endif
