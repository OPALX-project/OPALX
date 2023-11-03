//
// Class PartBunch.hpp
//   Particle Bunch.
//   A representation of a particle bunch as a vector of particles.
//
// Copyright (c) 2008 - 2020, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPAL_PartBunch_HPP
#define OPAL_PartBunch_HPP

class Distribution;

#include "Algorithms/CoordinateSystemTrafo.h"
#include "Algorithms/PartBins.h"
#include "Algorithms/PartData.h"
#include "Ippl.h"
#include "OPALtypes.h"
#include "Particle/ParticleAttrib.h"
#include "Particle/ParticleLayout.h"
#include "Utilities/Util.h"

/*

  Will change

*/
#include "PoissonSolvers/FFTOpenPoissonSolver.h"
#include "PoissonSolvers/FFTPeriodicPoissonSolver.h"
#include "PoissonSolvers/P3MSolver.h"

#include "PoissonSolvers/PoissonCG.h"

using ippl::detail::ConditionalType, ippl::detail::VariantFromConditionalTypes;

template <typename T = double, unsigned Dim = 3>
using FFTSolver_t = ConditionalType<
    Dim == 2 || Dim == 3, ippl::FFTPeriodicPoissonSolver<VField_t<T, Dim>, Field_t<Dim>>>;

template <typename T = double, unsigned Dim = 3>
using P3MSolver_t = ConditionalType<Dim == 3, ippl::P3MSolver<VField_t<T, Dim>, Field_t<Dim>>>;

template <typename T = double, unsigned Dim = 3>
using OpenSolver_t =
    ConditionalType<Dim == 3, ippl::FFTOpenPoissonSolver<VField_t<T, Dim>, Field_t<Dim>>>;

template <typename T = double, unsigned Dim = 3>
using Solver_t =
    VariantFromConditionalTypes<FFTSolver_t<T, Dim>, P3MSolver_t<T, Dim>, OpenSolver_t<T, Dim>>;

template <class PLayout, typename T, unsigned Dim = 3>
class PartBunch : public ippl::ParticleBase<PLayout> {
    using Base = ippl::ParticleBase<PLayout>;

public:
    Solver_t<T, Dim> solver_m;

    ParticleAttrib<double> Q;             // holds the charge for particle
    ParticleAttrib<double> M;             // holds the mass for particle
    ParticleAttrib<int> bunchNum;         // holds the mass for particle
    ParticleAttrib<Vector_t<T, Dim>> P;   // momenta
    ParticleAttrib<Vector_t<T, Dim>> E;   // E field vector f
    ParticleAttrib<Vector_t<T, Dim>> Ef;  // e field vector for gun simulations
    ParticleAttrib<Vector_t<T, Dim>> Bf;  // b field vector for gun simulations
    ParticleAttrib<int> Bin;     // energy bin in which the particle is in, if zero -> rebin
    ParticleAttrib<double> dt;   // holds the dt timestep for particle
    ParticleAttrib<double> Phi;  // the electric potential
    ParticleAttrib<Vector_t<T, Dim>> Eftmp;  // e field vector for gun simulations
    // ParticleAttrib<ParticleType> PType;  // particle names

    /**
      Reference particle structures
     */

    Vector_t<T, Dim> RefPartR_m;
    Vector_t<T, Dim> RefPartP_m;

    CoordinateSystemTrafo toLabTrafo_m;

    // ParticleOrigin refPOrigin_m;
    // ParticleType refPType_m;

    /// Initialize the translation vector and rotation quaternion
    /// here. Cyclotron tracker will reset these values each timestep
    /// TTracker can just use 0 translation and 0 rotation (quat[1 0 0 0]).
    // Vector_t globalMeanR_m = Vector_t(0.0, 0.0, 0.0);
    // Quaternion_t globalToLocalQuaternion_m = Quaternion_t(1.0, 0.0, 0.0, 0.0);
    Vector_t<T, Dim> globalMeanR_m;
    // Quaternion_t globalToLocalQuaternion_m;

    /**
       The structure for particle binning
    */

    PartBins* pbin_m;

    /// if larger than 0, emitt particles for tEmission_m [s]
    double tEmission_m;

    /// holds the gamma of the bin
    std::unique_ptr<double[]> bingamma_m;

    // FIXME: this should go into the Bin class!
    //  holds number of emitted particles of the bin
    //  jjyang: opal-cycl use *nBin_m of pbin_m
    std::unique_ptr<size_t[]> binemitted_m;

    /// steps per turn for OPAL-cycl
    int stepsPerTurn_m;

    /// step in a TRACK command
    long long localTrackStep_m;

    /// if multiple TRACK commands
    long long globalTrackStep_m;

    /// current bunch number
    short numBunch_m;

    /// number of particles per bunch
    std::vector<size_t> bunchTotalNum_m;
    std::vector<size_t> bunchLocalNum_m;

    /// this parameter records the current steps since last bunch injection
    /// it helps to inject new bunches correctly in the restart run of OPAL-cycl
    /// it is stored during phase space dump.
    int SteptoLastInj_m;

    bool fixed_grid;
    /// Mesh enlargement
    double dh_m;  /// relative enlargement of the mesh

    const PartData* reference;

    double couplingConstant_m;
    double qi_m;
    double massPerParticle_m;

    // unit state of PartBunch
    UnitState_t unit_state_m;
    UnitState_t stateOfLastBoundP_m;

    /// holds the centroid of the beam
    double centroid_m[2 * Dim];

    /// 6x6 matrix of the moments of the beam
    // matrix_t moments_m;

    /// holds the timestep in seconds
    double dt_m;
    /// holds the actual time of the integration
    double t_m;
    /// the position along design trajectory
    double spos_m;

    VField_t<T, Dim> EFD_m;
    VField_t<T, Dim> eg_m;

    Field_t<Dim> EFDMag_m;
    Field_t<Dim> rho_m;

    Field<T, Dim> phi_m;

    typedef ippl::BConds<Field<T, Dim>, Dim> bc_type;
    bc_type allPeriodic;

    /*
      Data structure for particle load balance information
    */

    ORB<T, Dim> orb;

    std::unique_ptr<size_t[]> globalPartPerNode_m;

    Vector_t<T, Dim> nr_m;

    std::array<bool, Dim> isParallel_m;

    Vector_t<double, Dim> hr_m;
    Vector_t<double, Dim> rmin_m;
    Vector_t<double, Dim> rmax_m;

    // flag to tell if we are a DC-beam
    bool dcBeam_m;
    double periodLength_m;

    double Q_m;

public:
    PartBunch(
        PLayout& pl, Vector_t<double, Dim> hr, Vector_t<double, Dim> rmin,
        Vector_t<double, Dim> rmax, std::array<bool, Dim> decomp, double Qtot)
        : ippl::ParticleBase<PLayout>(pl), hr_m(hr), rmin_m(rmin), rmax_m(rmax), Q_m(Qtot) {
        // register the particle attributes
        this->addAttribute(Q);
        this->addAttribute(M);
        this->addAttribute(P);
        this->addAttribute(E);
        this->addAttribute(Phi);
        this->addAttribute(Ef);
        this->addAttribute(Eftmp);
        this->addAttribute(Bf);
        this->addAttribute(Bin);
        this->addAttribute(dt);
        // this->addAttribute(PType);
        this->addAttribute(bunchNum);

        /// \todo setupBCs(a);
        for (unsigned int i = 0; i < Dim; i++) {
            isParallel_m[i] = decomp[i];
        }
    }

    // PartBunch() = delete;

    PartBunch& operator=(const PartBunch&) = delete;

    ~PartBunch() {
    }

    void initialize(FieldLayout_t<Dim>& fl, Mesh_t<Dim>& mesh) {
        updateLayout(fl, mesh);
    }

    /*

      Physics

    */

    double getCouplingConstant() const {
        return couplingConstant_m;
    }

    void setCouplingConstant(double c) {
        couplingConstant_m = c;
    }

    void calcBeamParameters() {
        get_bounds(rmin_m, rmax_m);
        // momentsComputer_m.compute(*this);
    }

    void calcBeamParametersInitial();  // Calculate initial beam parameters before emission.

    // set the charge per simulation particle
    void setCharge(double q) {
        qi_m = q;
        if (this->getTotalNum() != 0)
            Q = qi_m;
        else
            *ippl::Warn
                << "Could not set total charge in PartBunch::setCharge based on this->getTotalNum"
                << endl;
    }

    // set the charge per simulation particle when total particle number equals 0
    void setChargeZeroPart(double q) {
        qi_m = q;
    }

    // set the mass per simulation particle
    void setMass(double mass) {
        massPerParticle_m = mass;
        if (this->getTotalNum() != 0)
            M = mass;
    }

    void setMassZeroPart(double mass) {
        massPerParticle_m = mass;
    }

    /// get the total charge per simulation particle
    double getCharge() const {
        return 0.0;  // ADA sum(this->Q);
    }

    /// get the macro particle charge
    double getChargePerParticle() const {
        return qi_m;
    }

    double getMassPerParticle() const {
        return massPerParticle_m;
    }

    ///@{ Access to reference data
    double getQ() const {
        return reference->getQ();
    }
    double getM() const {
        return reference->getM();
    }
    double getP() const {
        return reference->getP();
    }

    double getE() const {
        return reference->getE();
    }

    // ParticleOrigin getPOrigin() const;
    // ParticleType getPType() const;

    double getInitialBeta() const {
        return reference->getBeta();
    }

    double getInitialGamma() const {
        return reference->getGamma();
    }

    ///@}
    ///@{ Set reference data

    void resetQ(double q) {
        const_cast<PartData*>(reference)->setQ(q);
    }
    void resetM(double m) {
        const_cast<PartData*>(reference)->setM(m);
    }

    // void setPOrigin(ParticleOrigin);
    void setPType(const std::string& type) {
        // refPType_m = ParticleProperties::getParticleType(type);
    }

    ///@}
    double getdE() const {
        return 0.0;
    }
    double getGamma(int i) const {
        return 0.0;
    }

    double getBeta(int i) const {
        return 0.0;
    }

    void actT() {
    }

    const PartData* getReference() const {
        return reference;
    }

    double getEmissionDeltaT() {
        return 0.0;
    }

    /*

        Domain Decomposition and Repartitioning

      */

    void do_binaryRepart() {
    }

    void updateDomainLength(Vector_t<int, 3>& grid) {
        ippl::NDIndex<3> domain = getFieldLayout().getDomain();
        for (unsigned int i = 0; i < Dim; i++)
            grid[i] = domain[i].length();
    }

    void updateFields(const Vector_t<T, Dim>& /*hr*/, const Vector_t<T, Dim>& origin) {
        /* \todo
        this->getMesh().set_meshSpacing(&(hr_m[0]));
        this->getMesh().set_origin(origin);
        rho_m.initialize(this->getMesh(), getFieldLayout(), 1);
        eg_m.initialize(this->getMesh(), getFieldLayout(), 1);
        */
    }

    void initializeORB(FieldLayout_t<Dim>& fl, Mesh_t<Dim>& mesh) {
        orb.initialize(fl, mesh, EFDMag_m);
    }

    void repartition(FieldLayout_t<Dim>& fl, Mesh_t<Dim>& mesh) {
        // Repartition the domains
        bool fromAnalyticDensity = false;
        bool res                 = orb.binaryRepartition(this->R, fl, fromAnalyticDensity);

        if (res != true) {
            std::cout << "Could not repartition!" << std::endl;
            return;
        }
        // Update
        this->updateLayout(fl, mesh);
    }

    bool balance(unsigned int totalP) {  //, int timestep = 1) {
        int local = 0;
        std::vector<int> res(ippl::Comm->size());
        double threshold = 1.0;
        double equalPart = (double)totalP / ippl::Comm->size();
        double dev       = std::abs((double)this->getLocalNum() - equalPart) / totalP;
        if (dev > threshold) {
            local = 1;
        }
        MPI_Allgather(&local, 1, MPI_INT, res.data(), 1, MPI_INT, ippl::Comm->getCommunicator());

        for (unsigned int i = 0; i < res.size(); i++) {
            if (res[i] == 1) {
                return true;
            }
        }
        return false;
    }

    void gatherLoadBalanceStatistics() {
        for (int i = 0; i < ippl::Comm->size(); i++)
            globalPartPerNode_m[i] = 0;

        std::size_t localnum = this->getLocalNum();
        // ADA gather(&localnum, &globalPartPerNode_m[0], 1);
    }

    size_t getLoadBalance(int p) const {
        return globalPartPerNode_m[p];
    }

    /*

      Mesh and Field Layout related functions

    */

    void resizeMesh();

    const Mesh_t<Dim>& getMesh() const {
        return &getFieldLayout().get_Mesh();
    }

    FieldLayout_t<Dim>& getFieldLayout() {
        // PLayout_t* layout = static_cast<Layout_t*>(&getLayout());
        //   return dynamic_cast<FieldLayout_t<Dim>&>(layout->getLayout().getFieldLayout());
        return this->getFieldLayout();
    }

    //     void setFieldLayout(FieldLayout_t* fLayout);

    void updateLayout(FieldLayout_t<Dim>& fl, Mesh_t<Dim>& mesh) {
        // Update local fields
        static IpplTimings::TimerRef tupdateLayout = IpplTimings::getTimer("updateLayout");
        IpplTimings::startTimer(tupdateLayout);
        this->EFD_m.updateLayout(fl);
        this->EFDMag_m.updateLayout(fl);

        // Update layout with new FieldLayout
        PLayout& layout = this->getLayout();
        layout.updateLayout(fl, mesh);
        IpplTimings::stopTimer(tupdateLayout);
        static IpplTimings::TimerRef tupdatePLayout = IpplTimings::getTimer("updatePB");
        IpplTimings::startTimer(tupdatePLayout);
        this->update();
        IpplTimings::stopTimer(tupdatePLayout);
    }

    bool isGridFixed() const {
        return fixed_grid;
    }

    void boundp() {
        /*
          Assume rmin_m < 0.0
         */

        // if (!R.isDirty() && stateOfLastBoundP_ == unit_state_) return;
        // if (!(this->R.isDirty() || this->ID.isDirty()) && stateOfLastBoundP_m == unit_state_m)
        //    return;  //-DW

        stateOfLastBoundP_m = unit_state_m;

        if (!isGridFixed()) {
            const int dimIdx = (dcBeam_m ? 2 : 3);

            /**
                In case of dcBeam_m && hr_m < 0
                this is the first call to boundp and we
                have to set hr completely i.e. x,y and z.
             */

            // \todo this->updateDomainLength(nr_m);
            get_bounds(rmin_m, rmax_m);
            Vector_t<T, Dim> len = rmax_m - rmin_m;

            double volume = 1.0;
            for (int i = 0; i < dimIdx; i++) {
                double length = std::abs(rmax_m[i] - rmin_m[i]);
                if (length < 1e-10) {
                    rmax_m[i] += 1e-10;
                    rmin_m[i] -= 1e-10;
                } else {
                    rmax_m[i] += dh_m * length;
                    rmin_m[i] -= dh_m * length;
                }
                hr_m[i] = (rmax_m[i] - rmin_m[i]) / (nr_m[i] - 1);
            }
            if (dcBeam_m) {
                rmax_m[2] = rmin_m[2] + periodLength_m;
                hr_m[2]   = periodLength_m / (nr_m[2] - 1);
            }
            for (int i = 0; i < dimIdx; ++i) {
                volume *= std::abs(rmax_m[i] - rmin_m[i]);
            }

            // ada if (volume < 1e-21 && this->this->getTotalNum() > 1 && std::abs(sum(Q)) >
            // 0.0) {
            //    WARNMSG(level1 << "!!! Extremely high particle density detected !!!" << endl);
            // }
            // *ippl::Info << "It is a full boundp hz= " << hr_m << " rmax= " << rmax_m << "
            // rmin= " << rmin_m
            // << endl;

            // ADA if (hr_m[0] * hr_m[1] * hr_m[2] <= 0) {
            // throw GeneralClassicException("boundp() ", "h<0, can not build a mesh");
            // }

            Vector_t<T, Dim> origin =
                rmin_m - Vector_t<T, Dim>(hr_m[0] / 2.0, hr_m[1] / 2.0, hr_m[2] / 2.0);
            this->updateFields(hr_m, origin);
        }
        this->update();
        // this->R.resetDirtyFlag();
    }

    /** This is only temporary in order to get the collimator and pepperpot working */
    size_t boundp_destroyT() {
        // \todo this->updateDomainLength(nr_m);

        std::vector<size_t> tmpbinemitted;

        boundp();

        size_t ne             = 0;
        const size_t localNum = this->getLocalNum();

        double rzmean             = 0.0;  // momentsComputer_m.getMeanPosition()(2);
        double rzrms              = 0.0;  // momentsComputer_m.getStandardDeviationPosition()(2);
        const bool haveEnergyBins = weHaveEnergyBins();
        if (haveEnergyBins) {
            tmpbinemitted.resize(getNumberOfEnergyBins(), 0.0);
        }
        for (unsigned int i = 0; i < localNum; i++) {
            if (this->Bin(i) < 0) {  // ADA || (Options::remotePartDel > 0
                // && std::abs(R(i)(2) - rzmean) < Options::remotePartDel * rzrms)) {
                ne++;
                // \todo this->destroy(1, i);
            } else if (haveEnergyBins) {
                tmpbinemitted[this->Bin(i)]++;
            }
        }

        boundp();

        calcBeamParameters();
        gatherLoadBalanceStatistics();

        //        ippl::Comm->reduce(ne, ne, 1, std::plus<size_t>());
        return ne;
    }

    size_t destroyT() {
        std::unique_ptr<size_t[]> tmpbinemitted;

        const size_t localNum = this->getLocalNum();
        const size_t totalNum = this->getTotalNum();
        size_t ne             = 0;

        if (weHaveEnergyBins()) {
            tmpbinemitted = std::unique_ptr<size_t[]>(new size_t[getNumberOfEnergyBins()]);
            for (int i = 0; i < getNumberOfEnergyBins(); i++) {
                tmpbinemitted[i] = 0;
            }
            for (unsigned int i = 0; i < localNum; i++) {
                if (this->Bin(i) < 0) {
                    // ADA destroy(1, i);
                    ++ne;
                } else
                    tmpbinemitted[this->Bin(i)]++;
            }
        } else {
            Inform dmsg("destroy: ", INFORM_ALL_NODES);
            for (size_t i = 0; i < localNum; i++) {
                if ((this->Bin(i) < 0)) {
                    ne++;
                    // ADA destroy(1, i);
                }
            }
        }
        /* ADA
           if (ne > 0) {
           performDestroy(true);
           }
        */

        calcBeamParameters();
        gatherLoadBalanceStatistics();

        size_t newTotalNum = this->getLocalNum();

        //
        ippl::Comm->reduce(newTotalNum, newTotalNum, 1, std::plus<size_t>());

        /// \todo ippl::Comm->reduce(newTotalNum, newTotalNum, OpAddAssign());
        /// \todo setTotalNum(newTotalNum);

        return totalNum - newTotalNum;
    }

    void set_meshEnlargement(double dh) {
        dh_m = dh;
    }

    /*

      Boundary Conditions

    */

    void setBCAllOpen() {
    }

    void setBCForDCBeam() {
    }

    void setupBCs() {
        setBCAllPeriodic();
    }

    void setBCAllPeriodic() {
        this->setParticleBC(ippl::BC::PERIODIC);
    }

    /*
      Misc Stuff
    */

    void runTests() {  // Field Solver
        Vector_t<T, Dim> ll(-0.005);
        Vector_t<T, Dim> ur(0.005);

        setBCAllPeriodic();

        ippl::NDIndex<Dim> domain = getFieldLayout().getDomain();
        for (unsigned int i = 0; i < Dim; i++)
            nr_m[i] = domain[i].length();

        for (int i = 0; i < 3; i++)
            hr_m[i] = (ur[i] - ll[i]) / nr_m[i];

        getFieldLayout().get_Mesh().set_meshSpacing(&(hr_m[0]));
        getFieldLayout().get_Mesh().set_origin(ll);

        rho_m.initialize(getFieldLayout().getMesh(), getFieldLayout(), 1);
        eg_m.initialize(getFieldLayout().getMesh(), getFieldLayout(), 1);

        // this->fs_m->solver_m->test(this);
    }

    void resetInterpolationCache(bool clearCache = false) {
        //    if (clearCache)
        //      interpolationCache_m.destroy(interpolationCache_m.size(), 0, true);
        //
    }

    void swap(unsigned int i, unsigned int j) {
        if (i >= this->getLocalNum() || j >= this->getLocalNum() || i == j)
            return;

        std::swap(this->R(i), this->R(j));
        std::swap(this->P(i), this->P(j));
        std::swap(this->Q(i), this->Q(j));
        std::swap(this->M(i), this->M(j));
        std::swap(this->Phi(i), this->Phi(j));
        /*
          std::swap(Ef(i), Ef[j]);
          std::swap(Eftmp(i), Eftmp[j]);
          std::swap(Bf(i), Bf[j]);
          std::swap(Bin(i), Bin[j]);
          std::swap(dt(i), dt[j]);
          std::swap(PType(i), PType[j]);
          std::swap(POrigin(i), POrigin[j]);
          std::swap(TriID(i), TriID[j]);
          std::swap(cavityGapCrossed(i), cavityGapCrossed[j]);
          std::swap(bunchNum(i), bunchNum[j]);
        */
    }

    double getRho(int x, int y, int z) {
        // return rho_m[x][y][z].get();
        return 0.0;
    }

    void gatherStatistics(unsigned int totalP) {
        ippl::Comm->barrier();
        std::cout << "Rank " << ippl::Comm->rank() << " has "
                  << (double)this->getLocalNum() / totalP * 100.0
                  << " percent of the total particles " << std::endl;
        ippl::Comm->barrier();
    }

    // FIXME: unify methods, use convention that all particles have own dt
    void switchToUnitlessPositions(bool use_dt_per_particle = false) {
        bool hasToReset = false;
        //   if (!this->R.isDirty())
        //    hasToReset = true;

        for (size_t i = 0; i < this->getLocalNum(); i++) {
            double dt = getdT();
            if (use_dt_per_particle)
                dt = this->dt(i);

            // this->R(i) /= Vector_t(Physics::c * dt);
            this->R(i) /= Vector_t<T, Dim>(dt);
        }

        unit_state_m = unitless;

        // if (hasToReset)
        //     this->R.resetDirtyFlag();
    }

    // FIXME: unify methods, use convention that all particles have own dt
    void switchOffUnitlessPositions(bool use_dt_per_particle = false) {
    }

    /** \brief returns the number of particles outside of a box defined by x */
    size_t calcNumPartsOutside(Vector_t<T, Dim> x) {
        std::size_t localnum         = 0;
        const Vector_t<T, Dim> meanR = get_rmean();

        for (unsigned long k = 0; k < this->getLocalNum(); ++k)
            if (std::abs(this->R(k)(0) - meanR(0)) > x(0)
                || std::abs(this->R(k)(1) - meanR(1)) > x(1)
                || std::abs(this->R(k)(2) - meanR(2)) > x(2)) {
                ++localnum;
            }

        // ADA gather(&localnum, &globalPartPerNode_m[0], 1);

        size_t npOutside = std::accumulate(
            globalPartPerNode_m.get(), globalPartPerNode_m.get() + ippl::Comm->size(), 0,
            std::plus<size_t>());

        return npOutside;
    }

    /**
     * \method calcLineDensity()
     * \brief calculates the 1d line density (not normalized) and append it to a file.
     * \see ParallelTTracker
     * \warning none yet
     *
     * DETAILED TODO
     *
     */
    void calcLineDensity(
        unsigned int nBins, std::vector<double>& lineDensity, std::pair<double, double>& meshInfo) {
        Vector_t<T, Dim> rmin, rmax;
        get_bounds(rmin, rmax);

        if (nBins < 2) {
            Vector<int, 3> /*NDIndex<3>*/ grid;
            // \todo this->updateDomainLength(grid);
            nBins = grid[2];
        }

        double length = rmax(2) - rmin(2);
        double zmin = rmin(2) - dh_m * length, zmax = rmax(2) + dh_m * length;
        double hz       = (zmax - zmin) / (nBins - 2);
        double perMeter = 1.0 / hz;  //(zmax - zmin);
        zmin -= hz;

        lineDensity.resize(nBins, 0.0);
        std::fill(lineDensity.begin(), lineDensity.end(), 0.0);

        const unsigned int lN = this->getLocalNum();
        for (unsigned int i = 0; i < lN; ++i) {
            const double z   = this->R(i)(2) - 0.5 * hz;
            unsigned int idx = (z - zmin) / hz;
            double tau       = (z - zmin) / hz - idx;

            // ADA lineDensity[idx] += Q(i) * (1.0 - tau) * perMeter;
            // lineDensity[idx + 1] += Q(i) * tau * perMeter;
        }

        // ADA reduce(&(lineDensity[0]), &(lineDensity[0]) + nBins, &(lineDensity[0]),
        // OpAddAssign());

        meshInfo.first  = zmin;
        meshInfo.second = hz;
    }

    void setBeamFrequency(double v) {
        periodLength_m = 1.;  // ADA Physics::c / f;
    }

    /*
        Compatibility function

    */

    VectorPair_t getEExtrema() {
        const Vector_t<T, Dim> maxE = max(eg_m);
        //      const double maxL = max(dot(eg_m,eg_m));
        const Vector_t<T, Dim> minE = min(eg_m);
        // *ippl::Info << "MaxE= " << maxE << " MinE= " << minE << endl;
        return VectorPair_t(maxE, minE);
    }

    /*

      Field Computations

    */

    void computeSelfFields() {
    }

    /** /brief used for self fields with binned distribution */
    void computeSelfFields(int b) {
    }

    // void setSolver(FieldSolver* fs);

    bool hasFieldSolver() {
        // if (this->fs_m)
        //   return this->fs_m->hasValidSolver();
        // else
        return false;
    }

    bool getFieldSolverType() const {
        // if (fs_m) {
        //     return true;  // fs_m->getFieldSolverType();
        // } else {
        return false;  // FieldSolverType::NONE;
                       //}
    }

    /*
      I / O



    Inform& print(Inform& os) {
        return print(os);
    }
    */
    /*
    Inform& operator<<(Inform& os, PartBunch<PLayout, T, Dim>& p) {
        return p.print(os);
    }

    */

    Inform& print(Inform& os) {
        // if (this->getLocalNum() != 0) {  // to suppress Nans
        Inform::FmtFlags_t ff = os.flags();

        double pathLength = get_sPos();

        os << std::scientific;
        os << level1 << "\n";
        os << "* ************** B U N C H "
              "********************************************************* \n";
        os << "* NP              = " << this->getLocalNum() << "\n";

        // os << "* Qtot            = " << std::setw(17) << Util::getChargeString(std::abs(sum(Q)))
        //    << "         "
        // os << "Qi    = " << std::setw(17) << Util::getChargeString(std::abs(qi_m)) << "\n";
        // os << "* Ekin            = " << std::setw(17)
        //   << Util::getEnergyString(get_meanKineticEnergy()) << "         "
        //   << "dEkin = " << std::setw(17) << Util::getEnergyString(getdE()) << "\n";
        // os << "* rmax            = " << Util::getLengthString(rmax_m, 5) << "\n";
        // os << "* rmin            = " << Util::getLengthString(rmin_m, 5) << "\n";
        /*
        if (this->getTotalNum() >= 2) {  // to suppress Nans
            os << "* rms beam size   = " << Util::getLengthString(get_rrms(), 5) << "\n";
            os << "* rms momenta     = " << std::setw(12) << std::setprecision(5) << get_prms()
               << " [beta gamma]\n";
            os << "* mean position   = " << Util::getLengthString(get_rmean(), 5) << "\n";
            os << "* mean momenta    = " << std::setw(12) << std::setprecision(5) << get_pmean()
               << " [beta gamma]\n";
            os << "* rms emittance   = " << std::setw(12) << std::setprecision(5) << get_emit()
               << " (not normalized)\n";
            os << "* rms correlation = " << std::setw(12) << std::setprecision(5) << get_rprms()
               << "\n";
        }

        os << "* hr              = " << Util::getLengthString(get_hr(), 5) << "\n";
        os << "* dh              = " << std::setw(13) << std::setprecision(5) << dh_m * 100
           << " [%]\n";
        os << "* t               = " << std::setw(17) << Util::getTimeString(getT()) << "         "
           << "dT    = " << std::setw(17) << Util::getTimeString(getdT()) << "\n";
        os << "* spos            = " << std::setw(17) << Util::getLengthString(pathLength) << "\n";
        */
        os << "* "
              "********************************************************************************** "
           << endl;
        os.flags(ff);
        // }
        return os;
    }

    // save particles in case of one core
    std::unique_ptr<Inform> pmsg_m;
    std::unique_ptr<std::ofstream> f_stream;

    void dumpData(int iteration) {
        auto view = P.getView();

        double Energy = 0.0;

        Kokkos::parallel_reduce(
            "Particle Energy", view.extent(0),
            KOKKOS_LAMBDA(const int i, double& valL) {
                double myVal = dot(view(i), view(i)).apply();
                valL += myVal;
            },
            Kokkos::Sum<double>(Energy));

        Energy *= 0.5;
        double gEnergy = 0.0;

        MPI_Reduce(&Energy, &gEnergy, 1, MPI_DOUBLE, MPI_SUM, 0, ippl::Comm->getCommunicator());

        Inform csvout(NULL, "data/energy.csv", Inform::APPEND);
        csvout.precision(10);
        csvout.setf(std::ios::scientific, std::ios::floatfield);

        csvout << iteration << " " << gEnergy << endl;

        ippl::Comm->barrier();
    }

    /*

    Field Particle Interpolation

    */

    void gatherCIC() {
        // static IpplTimings::TimerRef gatherTimer = IpplTimings::getTimer("gather");
        // IpplTimings::startTimer(gatherTimer);
        gather(this->E, EFD_m, this->R);
        // IpplTimings::stopTimer(gatherTimer);
    }

    void scatterCIC(unsigned int totalP, int iteration) {
        // static IpplTimings::TimerRef scatterTimer = IpplTimings::getTimer("scatter");
        // IpplTimings::startTimer(scatterTimer);
        Inform m("scatter ");
        EFDMag_m = 0.0;
        scatter(Q, EFDMag_m, this->R);
        // IpplTimings::stopTimer(scatterTimer);

        static IpplTimings::TimerRef sumTimer = IpplTimings::getTimer("CheckCharge");
        IpplTimings::startTimer(sumTimer);
        double Q_grid = EFDMag_m.sum();

        unsigned int Total_particles = 0;
        unsigned int local_particles = this->getLocalNum();

        MPI_Reduce(
            &local_particles, &Total_particles, 1, MPI_UNSIGNED, MPI_SUM, 0,
            ippl::Comm->getCommunicator());

        double rel_error = std::fabs((Q_m - Q_grid) / Q_m);
        m << "Rel. error in charge conservation = " << rel_error << endl;

        if (ippl::Comm->rank() == 0) {
            if (Total_particles != totalP || rel_error > 1e-10) {
                std::cout << "Total particles in the sim. " << totalP << " "
                          << "after update: " << Total_particles << std::endl;
                std::cout << "Total particles not matched in iteration: " << iteration << std::endl;
                std::cout << "Q grid: " << Q_grid << "Q particles: " << Q_m << std::endl;
                std::cout << "Rel. error in charge conservation: " << rel_error << std::endl;
                exit(1);
            }
        }

        IpplTimings::stopTimer(sumTimer);
    }

    /*

      Energy Bin Stuff

    */

    bool getIfBeamEmitting() {
        return true;
    }

    int getLastEmittedEnergyBin() {
        return 1;
    }

    size_t getNumberOfEmissionSteps() {
        return 1;
    }

    int getNumberOfEnergyBins() {
        return 1;
    }

    void Rebin() {
    }

    void setEnergyBins(int numberOfEnergyBins) {
    }

    bool weHaveEnergyBins() {
        return false;
    }

    void setTEmission(double t) {
    }

    double getTEmission() {
        return 1.0;
    }

    bool weHaveBins() const {
        return false;
    }

    void setPBins(PartBins* pbin) {
    }

    /** \brief Emit particles in the given bin
        i.e. copy the particles from the bin structure into the
        particle container
    */

    size_t emitParticles(double eZ) {
        return 0;
    }

    void updateNumTotal() {
        size_type total_particles = 0;
        size_type local_particles = this->getLocalNum();

        MPI_Reduce(
            &local_particles, &total_particles, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0,
            ippl::Comm->getCommunicator());
        // need to be stored
    }

    void rebin() {
        this->Bin = 0;
        // ADA pbin_m->resetBins();
        // delete pbin_m; we did not allocate it!
        pbin_m = nullptr;
    }

    int getLastemittedBin() {
        if (pbin_m != nullptr)
            return pbin_m->getLastemittedBin();
        else
            return 0;
    }

    void setLocalBinCount(size_t num, int bin) {
        if (pbin_m != nullptr) {
            // ADA pbin_m->setPartNum(bin, num);
        }
    }

    /** \brief Compute the gammas of all bins */
    void calcGammas() {
    }

    /** \brief Get gamma of one bin */
    double getBinGamma(int bin) {
        return 1.0;
    }

    bool hasBinning() const {
        return (pbin_m != nullptr);
    }

    /** \brief Set the charge of one bin to the value of q and all other to zero */
    void setBinCharge(int bin, double q) {
        // this->Q = where(eq(this->Bin, bin), q, 0.0);
    }

    /** \brief Set the charge of all other the ones in bin to zero */
    void setBinCharge(int bin) {
        //   this->Q = where(eq(this->Bin, bin), this->qi_m, 0.0);
    }

    /// calculate average angle of longitudinal direction of bins
    double calcMeanPhi() {
    }

    /// reset Bin[] for each particle according to the method given in paper PAST-AB(064402)
    /// by G. Fubiani et al.
    bool resetPartBinID2(const double eta) {
        return true;
    }

    bool resetPartBinBunch() {
        return true;
    }

    /*
      Part Bunch Access Functions

    */

    double getPx(int i) {
        return 0.0;
    }
    double getPy(int i) {
        return 0.0;
    }
    double getPz(int i) {
        return 0.0;
    }

    double getPx0(int i) {
        return 0.0;
    }
    double getPy0(int i) {
        return 0.0;
    }

    double getX(int i) {
        return 0.0;
    }
    double getY(int i) {
        return 0.0;
    }
    double getZ(int i) {
        return 0.0;
    }

    double getX0(int i) {
        return 0.0;
    }
    double getY0(int i) {
        return 0.0;
    }

    void setZ(int i, double zcoo) {
    }

    void get_bounds(Vector_t<T, Dim>& rmin, Vector_t<T, Dim>& rmax) const {
        this->getLocalBounds(rmin, rmax);

        double min[2 * Dim];

        for (unsigned int i = 0; i < Dim; ++i) {
            min[2 * i]     = rmin[i];
            min[2 * i + 1] = -rmax[i];
        }

        // ippl::Comm->allreduce(min, 2 * Dim, std::less<double>());

        for (unsigned int i = 0; i < Dim; ++i) {
            rmin[i] = min[2 * i];
            rmax[i] = -min[2 * i + 1];
        }
    }

    void getLocalBounds(Vector_t<T, Dim>& rmin, Vector_t<T, Dim>& rmax) const {
        const size_t localNum = this->getLocalNum();
        if (localNum == 0) {
            double maxValue = 1e8;
            rmin            = Vector_t<T, Dim>({maxValue, maxValue, maxValue});
            rmax            = Vector_t<T, Dim>({-maxValue, -maxValue, -maxValue});
            return;
        }

        rmin = this->R(0);
        rmax = this->R(0);
        for (size_t i = 1; i < localNum; ++i) {
            for (unsigned short d = 0; d < 3u; ++d) {
                if (rmin(d) > this->R(i)(d))
                    rmin(d) = this->R(i)(d);
                else if (rmax(d) < this->R(i)(d))
                    rmax(d) = this->R(i)(d);
            }
        }
    }

    std::pair<Vector_t<T, Dim>, double> getBoundingSphere() {
        Vector_t<T, Dim> rmin, rmax;
        get_bounds(rmin, rmax);

        std::pair<Vector_t<T, Dim>, double> sphere;
        sphere.first  = 0.5 * (rmin + rmax);
        sphere.second = std::sqrt(dot(rmax - sphere.first, rmax - sphere.first));

        return sphere;
    }

    std::pair<Vector_t<T, Dim>, double> getLocalBoundingSphere() {
        Vector_t<T, Dim> rmin, rmax;
        getLocalBounds(rmin, rmax);

        std::pair<Vector_t<T, Dim>, double> sphere;
        sphere.first  = 0.5 * (rmin + rmax);
        sphere.second = std::sqrt(dot(rmax - sphere.first, rmax - sphere.first));

        return sphere;
    }

    void get_PBounds(Vector_t<T, Dim>& min, Vector_t<T, Dim>& max) const {
        bounds(this->P, min, max);
    }

    /*
      Misc Bunch related quantities

    */

    // get 2nd order momentum matrix
    // matrix_t getSigmaMatrix() const;

    void setdT(double dt) {
        dt_m = dt;
    }

    double getdT() const {
        return dt_m;
    }

    void setT(double t) {
        t_m = t;
    }

    void incrementT() {
        t_m += dt_m;
    }

    double getT() const {
        return t_m;
    }

    /**
     * get the spos of the primary beam
     *
     * @param none
     *
     */

    void set_sPos(double s) {
        spos_m = s;
    }

    double get_gamma() const {
        return 1.00;
    }

    double get_meanKineticEnergy() {
        return 0.0;
    }

    Vector_t<T, Dim> get_origin() const {
        return rmin_m;
    }
    Vector_t<T, Dim> get_maxExtent() const {
        return rmax_m;
    }

    Vector_t<T, Dim> get_centroid() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_rrms() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_rprms() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_rmean() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_prms() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_pmean() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_pmean_Distribution() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_emit() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_norm_emit() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_halo() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_68Percentile() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_95Percentile() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_99Percentile() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_99_99Percentile() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_normalizedEps_68Percentile() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_normalizedEps_95Percentile() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_normalizedEps_99Percentile() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_normalizedEps_99_99Percentile() const {
        return Vector_t<T, Dim>(0.0);
    }
    Vector_t<T, Dim> get_hr() const {
        return Vector_t<T, Dim>(0.0);
    }

    double get_Dx() const {
        return 0.0;
    }
    double get_Dy() const {
        return 0.0;
    }
    double get_DDx() const {
        return 0.0;
    }
    double get_DDy() const {
        return 0.0;
    }

    /*
      Some quantities related to integrations/tracking
     */

    void setStepsPerTurn(int n) {
        stepsPerTurn_m = n;
    }

    int getStepsPerTurn() const {
        return stepsPerTurn_m;
    }

    /// step in multiple TRACK commands
    void setGlobalTrackStep(long long n) {
        globalTrackStep_m = n;
    }

    long long getGlobalTrackStep() const {
        return globalTrackStep_m;
    }

    /// step in a TRACK command
    void setLocalTrackStep(long long n) {
        localTrackStep_m = n;
    }

    void incTrackSteps() {
        globalTrackStep_m++;
        localTrackStep_m++;
    }

    long long getLocalTrackStep() const {
        return localTrackStep_m;
    }

    void setNumBunch(short n) {
        numBunch_m = n;
        bunchTotalNum_m.resize(n);
        bunchLocalNum_m.resize(n);
    }

    short getNumBunch() const {
        return numBunch_m;
    }

    void setGlobalMeanR(Vector_t<T, Dim> globalMeanR) {
        globalMeanR_m = globalMeanR;
    }

    Vector_t<T, Dim> getGlobalMeanR() {
        return globalMeanR_m;
    }

    // void setGlobalToLocalQuaternion(Quaternion_t globalToLocalQuaternion);
    // Quaternion_t getGlobalToLocalQuaternion();

    void setSteptoLastInj(int n) {
        SteptoLastInj_m = n;
    }

    int getSteptoLastInj() const {
        return SteptoLastInj_m;
    }

    double calculateAngle(double x, double y) {
        return 0.0;
    }

    double get_sPos() const {
        return spos_m;
    }
};

typedef PartBunch<PLayout_t<double, 3>, double, 3> PartBunch_t;

#endif
