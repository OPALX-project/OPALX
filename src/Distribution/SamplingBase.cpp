#include "SamplingBase.hpp"
#include <mpi.h>
#include <algorithm>
#include <vector>

size_t SamplingBase::computeLocalEmitCount(size_t totalToSample) const {
    if (!pc_m) {
        return 0;
    }
    const int nranks = ippl::Comm->size();
    const int rank   = ippl::Comm->rank();
    if (nranks <= 0) {
        return totalToSample;
    }
    if (totalToSample == 0) {
        return 0;
    }

    const size_t capacity = pc_m->R.getView().extent(0);
    const size_t localNum = pc_m->getLocalNum();
    const size_t spaceLeft =
        (localNum <= capacity) ? (capacity - localNum) : size_t(0);

    std::vector<unsigned long> spaceLeftAll(static_cast<size_t>(nranks), 0);
    unsigned long mySpace = static_cast<unsigned long>(spaceLeft);
    MPI_Allgather(&mySpace, 1, MPI_UNSIGNED_LONG, spaceLeftAll.data(), 1,
                  MPI_UNSIGNED_LONG, ippl::Comm->getCommunicator());

    const size_t nranks_u = static_cast<size_t>(nranks);
    const size_t base     = totalToSample / nranks_u;
    const size_t rem      = totalToSample % nranks_u;

    std::vector<size_t> nlocal(static_cast<size_t>(nranks), 0);
    size_t sum = 0;
    for (int r = 0; r < nranks; ++r) {
        size_t ideal = base + (static_cast<size_t>(r) < rem ? 1 : 0);
        size_t cap   = static_cast<size_t>(spaceLeftAll[static_cast<size_t>(r)]);
        nlocal[static_cast<size_t>(r)] = std::min(ideal, cap);
        sum += nlocal[static_cast<size_t>(r)];
    }

    size_t deficit = totalToSample - sum;
    while (deficit > 0) {
        int chosen = -1;
        for (int r = 0; r < nranks && deficit > 0; ++r) {
            size_t cap = static_cast<size_t>(spaceLeftAll[static_cast<size_t>(r)]);
            if (nlocal[static_cast<size_t>(r)] < cap) {
                chosen = r;
                break;
            }
        }
        if (chosen < 0) {
            break;
        }
        nlocal[static_cast<size_t>(chosen)] += 1;
        deficit -= 1;
    }

    return nlocal[static_cast<size_t>(rank)];
}
