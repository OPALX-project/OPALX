#ifndef OPAL_LOAD_BALANCER_H
#define OPAL_LOAD_BALANCER_H

#include <cmath>
#include <memory>
#include <vector>

#include "PartBunch/FieldContainer.hpp"
#include "PartBunch/ParticleContainer.hpp"

template <typename T, unsigned Dim>
using ORB = ippl::OrthogonalRecursiveBisection<Field<double, Dim>, T>;

template <typename T, unsigned Dim>
class LoadBalancer {
    using Base = ippl::ParticleBase<
            ippl::ParticleSpatialLayout<T, Dim,ippl::UniformCartesian<T,Dim>>, Kokkos::DefaultExecutionSpace::memory_space>;
    using FieldSolver_t = ippl::FieldSolverBase<T, Dim>;

private:
    double loadbalancethreshold_m;
    std::shared_ptr<FieldContainer<T, Dim>> fc_m;
    std::shared_ptr<ParticleContainer<T, Dim>> pc_m;
    std::shared_ptr<FieldSolver_t> fs_m;
    ORB<T, Dim> orb;
    size_type numBalances;

public:
    LoadBalancer(
            double lbs, std::shared_ptr<FieldContainer<T, Dim>>& fc,
            std::shared_ptr<ParticleContainer<T, Dim>>& pc, std::shared_ptr<FieldSolver_t>& fs)
        : loadbalancethreshold_m(lbs), fc_m(fc), pc_m(pc), fs_m(fs), numBalances(0) {}

    ~LoadBalancer() {}

    double getLoadBalanceThreshold() const { return loadbalancethreshold_m; }
    void setLoadBalanceThreshold(double threshold) { loadbalancethreshold_m = threshold; }

    size_type getNumBalances() const { return numBalances; }

    Field_t<Dim>* getRho() const { return fc_m ? &fc_m->getRho() : nullptr; }

    VField_t<T, Dim>* getE() const { return fc_m ? &fc_m->getE() : nullptr; }

    Field<T, Dim>* getPhi() { return fc_m ? &fc_m->getPhi() : nullptr; }

    std::shared_ptr<ParticleContainer<T, Dim>> getParticleContainer() const { return pc_m; }
    void setParticleContainer(std::shared_ptr<ParticleContainer<T, Dim>> pc) { pc_m = pc; }

    std::shared_ptr<FieldSolver_t> getFieldSolver() const { return fs_m; }
    void setFieldSolver(std::shared_ptr<FieldSolver_t> fs) { fs_m = fs; }

    void updateLayout(ippl::FieldLayout<Dim>* fl, ippl::UniformCartesian<T, Dim>* mesh) {
        // Update local fields

        static IpplTimings::TimerRef tupdateLayout = IpplTimings::getTimer("updateLayout");
        IpplTimings::startTimer(tupdateLayout);
        fc_m->updateFieldLayoutsAfterLayoutChange(fs_m ? fs_m->getStype() : "");

        // Update layout with new FieldLayout
        PLayout_t<T, Dim>* layout = &pc_m->getLayout();
        (*layout).updateLayout(*fl, *mesh);
        IpplTimings::stopTimer(tupdateLayout);
        static IpplTimings::TimerRef tupdatePLayout = IpplTimings::getTimer("updatePB");
        IpplTimings::startTimer(tupdatePLayout);
        pc_m->update();
        pc_m->markMomentsDirty();
        IpplTimings::stopTimer(tupdatePLayout);
    }

    void initializeORB(ippl::FieldLayout<Dim>* fl, ippl::UniformCartesian<T, Dim>* mesh) {
        orb = ORB<T, Dim>();
        orb.initialize(*fl, *mesh, fc_m->getRho());
    }

    bool repartition(ippl::FieldLayout<Dim>* fl, ippl::UniformCartesian<T, Dim>* mesh) {
        // Repartition the domains

        using Base = ippl::ParticleBase<
                ippl::ParticleSpatialLayout<T, Dim,ippl::UniformCartesian<T,Dim>>, Kokkos::DefaultExecutionSpace::memory_space>;
        typename Base::particle_position_type* R;
        R = &pc_m->R;

        /*
        Currently, only a all parallel decomposition is supported by IPPL. This might change with
        https://github.com/IPPL-framework/ippl/pull/560 - but until then, we enforce that all
        dimensions are parallel. A similar check is implemented in
        FieldSolverCmd::setDomainDecomposition() to ensure consistency with the user input.

        Another note: "isFirstRepartition" could be set to true (skips the scatter step inside ORB)
        if the rho field is reused from the previous scatter. However, should we use binning, then
        this rho field only contains partial particle information and cannot be used. That's why -
        for now - it's way easier to just keep it at false. This might be a possible optimization in
        the future.
        */
        // bool res = orb.binaryRepartition(*R, *fl, false, fl->isParallel());
        bool res = orb.binaryRepartition(*R, *fl, false);
        if (res != true) {
            if (ippl::Comm->rank() == 0) {
                std::cout << "Could not repartition!" << std::endl;
            }
            return false;
        }
        // Update
        this->updateLayout(fl, mesh);
        numBalances++;  // increment only when repartition was actually performed
        return true;
    }

    bool repartitionFromCurrentParticles(
            ippl::FieldLayout<Dim>* fl, ippl::UniformCartesian<T, Dim>* mesh) {
        initializeORB(fl, mesh);
        return repartition(fl, mesh);
    }

    bool balance(size_type totalP) {
        if (totalP == 0) {
            return false;
        }
        if (ippl::Comm->size() < 2) {
            return false;
        } else {
            int local = 0;
            std::vector<int> res(ippl::Comm->size());
            double equalPart = (double)totalP / ippl::Comm->size();
            double dev       = std::abs((double)pc_m->getLocalNum() - equalPart) / totalP;
            if (dev > loadbalancethreshold_m) {
                local = 1;
            }
            MPI_Allgather(
                    &local, 1, MPI_INT, res.data(), 1, MPI_INT, ippl::Comm->getCommunicator());

            for (unsigned int i = 0; i < res.size(); i++) {
                if (res[i] == 1) {
                    return true;
                }
            }
            return false;
        }
    }
};

#endif
