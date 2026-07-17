#ifndef OPAL_LOAD_BALANCER_H
#define OPAL_LOAD_BALANCER_H

#include <cstddef>

/**
 * @brief Compatibility state retained while PartBunch still derives from IPPL PicManager.
 *
 * PIC domain updates and ORB repartitioning are owned by PicDomainManager. This shell remains
 * only because PicManager currently carries a load-balancer template parameter and pointer.
 */
template <typename T, unsigned Dim>
class LoadBalancer {
public:
    explicit LoadBalancer(double threshold) : threshold_m(threshold) {}

    double getLoadBalanceThreshold() const { return threshold_m; }
    void setLoadBalanceThreshold(double threshold) { threshold_m = threshold; }

    std::size_t getNumBalances() const { return 0; }

private:
    double threshold_m;
};

#endif
