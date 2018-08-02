#ifndef __GEOMETRIC_STRATEGY_H__
#define __GEOMETRIC_STRATEGY_H__

#include <set>
#include <vector>

#include "Comm/Splitter/SplitStrategy.h"

typedef std::vector<int> coordinates_t;

typedef struct {
    coordinates_t extensions;
    coordinates_t origin;
    MPI_Comm      comm;
    int           master_pid;
    std::vector<int> worker_pids;
} region_t;

/**
 * The geometric strategy partitions the network graph in num_masters equally
 * sized regions and then places a master at the centroid of each region.
 * Workers inside a region are assigned to the master at its centroid.
 */
class GeometricStrategy : protected SplitStrategy {

public:

    GeometricStrategy(size_t num_masters,
                      boost::shared_ptr<CommTopology> topology, MPI_Comm comm)
        : Strategy(num_masters, topology, comm)
    {}

    virtual ~GeometricStrategy()
    {}

    void split() {

        // 1. Initialize the region list with the initial topology
        region_t initial_region;
        for(int dim=0; dim < topology_.getNumDimensions(); dim++) {
            initial_region.origin.push_back(0);
            initial_region.extensions.push_back(dims_[dim]);
        }
        regions_.push_back(initial_region);

        // 2. Build num_masters_ partitions of the topology
        while(regions_.size() < num_masters_) {
            splitInLargestDim();
        }

        // 3. Compute centroids of partition and assign master
        for(region_t & region : regions_) {
            computeCentroid(region);
        }

        // 4. Create communicator groups
        for(region_t & region : regions_) {
            bool is_worker = false;

            coordinates_t iter = region.origin;
            for(int x = 0; x < region.extensions[0]; x++) {
                iter[0] = region.origin[0] + x;
                for(int y = 0; y < region.extensions[1]; y++) {
                    iter[1] = region.origin[1] + y;
                    for(int z = 0; z < region.extensions[2]; z++) {
                        iter[2] = region.origin[2] + z;

                        std::vector<int> worker_pids;
                        coordinatesToPID(iter, worker_pids);

                        if(worker_pids[0] != region.master_pid)
                            region.worker_pids.push_back(worker_pids[0]);

                        region.worker_pids.push_back(worker_pids[1]);
                        region.worker_pids.push_back(worker_pids[2]);
                        region.worker_pids.push_back(worker_pids[3]);

                        if(rank_ == worker_pids[0] ||
                           rank_ == worker_pids[1] ||
                           rank_ == worker_pids[2] ||
                           rank_ == worker_pids[3] ) {
                            is_worker = true;
                        }

                    }
                }
            }

            if(rank_ == region.master_pid || is_worker)
                addCommGroup(region.master_pid, region.worker_pids);

            MPI_Barrier(comm_);
        }
    }

protected:

    std::vector<region_t> regions_;

    void splitInLargestDim() {

        // first we determine the largest extension of all remaining regions
        int origin = 0;
        int split_in_direction = 0;
        int max_ext_dir = 0;

        for(region_t region : regions_) {
            for(int dim = 0; dim < topology_.getNumDimensions(); dim++) {
                if(region.extensions[dim] > max_ext_dir) {
                    split_in_direction = dim;
                    origin = region.origin[dim];
                    max_ext_dir = region.extensions[dim];
                }
            }
        }

        // now we can perform splits as long as the number of requested
        // regions is not reached
        int split_at_coordinate = origin;
        split_at_coordinate += (int)(max_ext_dir / 2.0);

        splitRegionAt(split_at_coordinate, split_in_direction);
    }

    void splitRegionAt(int split_at_coordinate, int split_in_direction) {

        int removed = 0;
        bool continue_splitting = true;
        std::vector<region_t> new_regions;
        std::vector<region_t>::iterator itr;
        for(itr = regions_.begin(); itr != regions_.end(); itr++) {

            if(regions_.size() + new_regions.size() - removed == num_masters_)
                continue_splitting = false;

            double origin = itr->origin[split_in_direction];
            double extension = itr->extensions[split_in_direction];
            if( continue_splitting &&
                isSplittable(origin, split_at_coordinate, extension) ) {

                // clone region
                region_t second = *itr;

                // and split
                second.origin[split_in_direction]     = split_at_coordinate;
                second.extensions[split_in_direction] = extension - (split_at_coordinate - origin);

                itr->extensions[split_in_direction]   = split_at_coordinate - origin;

                new_regions.push_back(*itr);
                new_regions.push_back(second);
                removed++;

            } else {

                new_regions.push_back(*itr);
                removed++;
            }
        }

        regions_ = new_regions;
    }

    bool isSplittable(double origin, double split_at_coordinate, double extension) {
        return (origin < split_at_coordinate &&
                split_at_coordinate < origin + extension);
    }

    void computeCentroid(region_t &region) {

        coordinates_t centroid;
        centroid.resize(3, 0);

        for(unsigned int pos=0; pos < region.extensions.size(); pos++) {
            int centroid_coordinate = (int)(region.origin[pos]*1.0 + region.extensions[pos]/2.0);
            centroid[pos] = centroid_coordinate;
        }

        std::vector<int> pids;
        coordinatesToPID(centroid, pids);
        region.master_pid = pids[0];
    }

    void coordinatesToPID(coordinates_t coordinate, std::vector<int> &pid) {

        pid.resize(4, 0);

        if(coords_[0] == coordinate[0] &&
           coords_[1] == coordinate[1] &&
           coords_[2] == coordinate[2] ) {
            pid[my_core_] = rank_;
        }

        MPI_Allreduce(MPI_IN_PLACE, &pid[0], 4, MPI_INT, MPI_SUM, comm_);
    }

    void printRegion(region_t region) {
        if(rank_ == 0) {
            std::cout << "region: ";
            std::cout << "[" << region.origin[0] << ", " << region.origin[1] << ", " << region.origin[2] << "], ";
            std::cout << "[" << region.extensions[0] << ", " << region.extensions[1] << ", " << region.extensions[2] << "]";
            std::cout << std::endl;
        }
    }

};

#endif
