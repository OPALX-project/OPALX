//
// Class PartBunch
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
#ifndef OPAL_PartBunch_HH
#define OPAL_PartBunch_HH

class Distribution;
class PartBins;
class PartBinsCyc;
class PartData;

#include "Ippl.h"
#include "OPALtypes.h"

#include "Particle/ParticleAttrib.h"
#include "Particle/ParticleLayout.h"

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

    // CoordinateSystemTrafo toLabTrafo_m;

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

    ippl::e_dim_tag decomp_m[Dim];

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
        Vector_t<double, Dim> rmax, ippl::e_dim_tag decomp[Dim], double Qtot)
        /*
            : ippl::ParticleBase<PLayout>(pl)
        */
        : hr_m(hr), rmin_m(rmin), rmax_m(rmax), Q_m(Qtot) {
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

        // ADA setupBCs(a);
        for (unsigned int i = 0; i < Dim; i++) {
            decomp_m[i] = decomp[i];
        }
    }

    // PartBunch() = delete;

    PartBunch& operator=(const PartBunch&) = delete;

    ~PartBunch() {
    }

    void initialize(FieldLayout_t<Dim>& fLayout, Mesh_t<Dim>& mesh);

    /*

      Physics

    */

    double getCouplingConstant() const;
    void setCouplingConstant(double c);

    void calcBeamParameters();
    void calcBeamParametersInitial();  // Calculate initial beam parameters before emission.

    // set the charge per simulation particle
    void setCharge(double q);
    // set the charge per simulation particle when total particle number equals 0
    void setChargeZeroPart(double q);

    // set the mass per simulation particle
    void setMass(double mass);
    void setMassZeroPart(double mass);

    /// get the total charge per simulation particle
    double getCharge() const;

    /// get the macro particle charge
    double getChargePerParticle() const;

    double getMassPerParticle() const;

    ///@{ Access to reference data
    double getQ() const;
    double getM() const;
    double getP() const;
    double getE() const;
    // ParticleOrigin getPOrigin() const;
    //  ParticleType getPType() const;
    double getInitialBeta() const;
    double getInitialGamma() const;
    ///@}
    ///@{ Set reference data

    void resetQ(double q);
    void resetM(double m);
    // void setPOrigin(ParticleOrigin);
    void setPType(const std::string& type);
    ///@}
    double getdE() const;
    double getGamma(int i);
    double getBeta(int i);
    void actT();

    const PartData* getReference() const;

    double getEmissionDeltaT();

    /*

        Domain Decomposition and Repartitioning

      */

    void do_binaryRepart();

    void updateDomainLength(Vector_t<int, 3>& grid);

    void updateFields(const Vector_t<T, Dim>& /*hr*/, const Vector_t<T, Dim>& origin);

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

    void gatherLoadBalanceStatistics();
    size_t getLoadBalance(int p) const;

    /*

      Mesh and Field Layout related functions

    */

    void resizeMesh();

    const Mesh_t<Dim>& getMesh() const {
        return &getFieldLayout().get_Mesh();
    }

    FieldLayout_t<Dim>& getFieldLayout();
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

    bool isGridFixed() const;

    void boundp();

    /** This is only temporary in order to get the collimator and pepperpot working */
    size_t boundp_destroyT();

    size_t destroyT();

    void set_meshEnlargement(double dh);

    /*

      Boundary Conditions

    */

    void setBCAllOpen();

    void setBCForDCBeam();

    void setupBCs() {
        setBCAllPeriodic();
    }

    void setBCAllPeriodic() {
        this->setParticleBC(ippl::BC::PERIODIC);
    }

    /*
      Misc Stuff
    */

    void runTests();  // Field Solver

    void resetInterpolationCache(bool clearCache = false);

    void swap(unsigned int i, unsigned int j);

    double getRho(int x, int y, int z);

    void gatherStatistics(unsigned int totalP) {
        ippl::Comm->barrier();
        std::cout << "Rank " << ippl::Comm->rank() << " has "
                  << (double)this->getLocalNum() / totalP * 100.0
                  << " percent of the total particles " << std::endl;
        ippl::Comm->barrier();
    }

    // FIXME: unify methods, use convention that all particles have own dt
    void switchToUnitlessPositions(bool use_dt_per_particle = false);

    // FIXME: unify methods, use convention that all particles have own dt
    void switchOffUnitlessPositions(bool use_dt_per_particle = false);

    /** \brief returns the number of particles outside of a box defined by x */
    size_t calcNumPartsOutside(Vector_t<T, Dim> x);

    void calcLineDensity(
        unsigned int nBins, std::vector<double>& lineDensity, std::pair<double, double>& meshInfo);

    void setBeamFrequency(double v);

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

    void computeSelfFields();

    /** /brief used for self fields with binned distribution */
    void computeSelfFields(int b);

    // void setSolver(FieldSolver* fs);

    bool hasFieldSolver();

    // FieldSolverType getFieldSolverType() const;

    /*
      I / O

    */

    Inform& print(Inform& os);

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

    bool getIfBeamEmitting();

    int getLastEmittedEnergyBin();

    size_t getNumberOfEmissionSteps();

    int getNumberOfEnergyBins();

    void Rebin();

    void setEnergyBins(int numberOfEnergyBins);

    bool weHaveEnergyBins();

    void setTEmission(double t);

    double getTEmission();

    bool weHaveBins() const;

    void setPBins(PartBins* pbin);

    void setPBins(PartBinsCyc* pbin);

    /** \brief Emit particles in the given bin
        i.e. copy the particles from the bin structure into the
        particle container
    */
    size_t emitParticles(double eZ);

    void updateNumTotal();

    void rebin();

    int getLastemittedBin();

    void setLocalBinCount(size_t num, int bin);

    /** \brief Compute the gammas of all bins */
    void calcGammas();

    /** \brief Get gamma of one bin */
    double getBinGamma(int bin);

    bool hasBinning() const;

    /** \brief Set the charge of one bin to the value of q and all other to zero */
    void setBinCharge(int bin, double q);

    /** \brief Set the charge of all other the ones in bin to zero */
    void setBinCharge(int bin);

    /// calculate average angle of longitudinal direction of bins
    double calcMeanPhi();

    /// reset Bin[] for each particle according to the method given in paper PAST-AB(064402) by  G.
    /// Fubiani et al.
    bool resetPartBinID2(const double eta);

    bool resetPartBinBunch();

    /*
      Part Bunch Access Functions

    */

    double getPx(int i);
    double getPy(int i);
    double getPz(int i);

    double getPx0(int i);
    double getPy0(int i);

    double getX(int i);
    double getY(int i);
    double getZ(int i);

    double getX0(int i);
    double getY0(int i);

    void setZ(int i, double zcoo);

    void get_bounds(Vector_t<T, Dim>& rmin, Vector_t<T, Dim>& rmax) const;
    void get_PBounds(Vector_t<T, Dim>& rmin, Vector_t<T, Dim>& rmax) const;

    void getLocalBounds(Vector_t<T, Dim>& rmin, Vector_t<T, Dim>& rmax) const;

    std::pair<Vector_t<T, Dim>, double> getBoundingSphere();

    std::pair<Vector_t<T, Dim>, double> getLocalBoundingSphere();

    /*
      Misc Bunch related quantities

    */

    // get 2nd order momentum matrix
    // matrix_t getSigmaMatrix() const;

    void setdT(double dt);
    double getdT() const;

    void setT(double t);
    void incrementT();
    double getT() const;

    /**
     * get the spos of the primary beam
     *
     * @param none
     *
     */
    double get_sPos() const;

    void set_sPos(double s);

    double get_gamma() const;
    double get_meanKineticEnergy() const;
    Vector_t<T, Dim> get_origin() const;
    Vector_t<T, Dim> get_maxExtent() const;
    Vector_t<T, Dim> get_centroid() const;
    Vector_t<T, Dim> get_rrms() const;
    Vector_t<T, Dim> get_rprms() const;
    Vector_t<T, Dim> get_rmean() const;
    Vector_t<T, Dim> get_prms() const;
    Vector_t<T, Dim> get_pmean() const;
    Vector_t<T, Dim> get_pmean_Distribution() const;
    Vector_t<T, Dim> get_emit() const;
    Vector_t<T, Dim> get_norm_emit() const;
    Vector_t<T, Dim> get_halo() const;
    Vector_t<T, Dim> get_68Percentile() const;
    Vector_t<T, Dim> get_95Percentile() const;
    Vector_t<T, Dim> get_99Percentile() const;
    Vector_t<T, Dim> get_99_99Percentile() const;
    Vector_t<T, Dim> get_normalizedEps_68Percentile() const;
    Vector_t<T, Dim> get_normalizedEps_95Percentile() const;
    Vector_t<T, Dim> get_normalizedEps_99Percentile() const;
    Vector_t<T, Dim> get_normalizedEps_99_99Percentile() const;
    Vector_t<T, Dim> get_hr() const;

    double get_Dx() const;
    double get_Dy() const;
    double get_DDx() const;
    double get_DDy() const;

    /*
      Some quantities related to integrations/tracking


     */

    void setStepsPerTurn(int n);
    int getStepsPerTurn() const;

    /// step in multiple TRACK commands
    void setGlobalTrackStep(long long n);
    long long getGlobalTrackStep() const;

    /// step in a TRACK command
    void setLocalTrackStep(long long n);
    void incTrackSteps();
    long long getLocalTrackStep() const;
    void setNumBunch(short n);
    short getNumBunch() const;

    void setGlobalMeanR(Vector_t<T, Dim> globalMeanR);
    Vector_t<T, Dim> getGlobalMeanR();

    // void setGlobalToLocalQuaternion(Quaternion_t globalToLocalQuaternion);
    // Quaternion_t getGlobalToLocalQuaternion();

    void setSteptoLastInj(int n);
    int getSteptoLastInj() const;

    double calculateAngle(double x, double y);

    /*
        void writePerRank() {
            double lq = 0.0, lqm = 0.0;
            auto viewRho = this->EFDMag_m.getView();
            auto viewqm  = this->qm.getView();
            int nghost   = this->EFDMag_m.getNghost();

            using index_array_type = typename ippl::RangePolicy<Dim>::index_array_type;
            ippl::parallel_reduce(
                "Particle Charge", ippl::getRangePolicy(viewRho, nghost),
                KOKKOS_LAMBDA(const index_array_type& args, double& val) {
                    val += ippl::apply(viewRho, args);
                },
                lq);
            Kokkos::parallel_reduce(
                "Particle QM", viewqm.extent(0),
                KOKKOS_LAMBDA(const int i, double& val) { val += viewqm(i); }, lqm);
        }

        void initFields() {
            static IpplTimings::TimerRef initFieldsTimer = IpplTimings::getTimer("initFields");
            IpplTimings::startTimer(initFieldsTimer);
            Inform m("initFields ");

            ippl::NDIndex<Dim> domain = EFD_m.getDomain();

            for (unsigned int i = 0; i < Dim; i++) {
                nr_m[i] = domain[i].length();
            }

            double phi0 = 0.1;
            double pi   = Kokkos::numbers::pi_v<double>;
            // scale_fact so that particles move more
            double scale_fact = 1e5;  // 1e6

            Vector_t<T, Dim> hr = hr_m;

            auto view                        = EFD_m.getView();
            const FieldLayout_t<Dim>& layout = EFD_m.getLayout();
            const ippl::NDIndex<Dim>& lDom   = layout.getLocalNDIndex();
            const int nghost                 = EFD_m.getNghost();

            using index_array_type = typename ippl::RangePolicy<Dim>::index_array_type;
            ippl::parallel_for(
                "Assign EFD_m", ippl::getRangePolicy(view, nghost),
                KOKKOS_LAMBDA(const index_array_type& args) {
                    // local to global index conversion
                    Vector_t<T, Dim> vec = (0.5 + args + lDom.first() - nghost) * hr;

                    ippl::apply(view, args)[0] = -scale_fact * 2.0 * pi * phi0;
                    for (unsigned d1 = 0; d1 < Dim; d1++) {
                        ippl::apply(view, args)[0] *= Kokkos::cos(2 * ((d1 + 1) % 3) * pi *
      vec[d1]);
                    }
                    for (unsigned d = 1; d < Dim; d++) {
                        ippl::apply(view, args)[d] = scale_fact * 4.0 * pi * phi0;
                        for (int d1 = 0; d1 < (int)Dim - 1; d1++) {
                            ippl::apply(view, args)[d] *=
                                Kokkos::sin(2 * ((d1 + 1) % 3) * pi * vec[d1]);
                        }
                    }
                });

            EFDMag_m = dot(EFD_m, EFD_m);
            EFDMag_m = sqrt(EFDMag_m);
            IpplTimings::stopTimer(initFieldsTimer);
        }


      // @param tag
      //        2 -> uniform(0,1)
      //        1 -> normal(0,1)
      //        0 -> gridpoints
      void initPositions(FieldLayout_t<Dim>& fl, Vector_t<T, Dim>& hr, unsigned int nloc,
                         int tag = 2) {
          Inform m("initPositions ");
          typename ippl::ParticleBase<PLayout>::particle_position_type::HostMirror R_host =
              this->R.getHostMirror();

          std::mt19937_64 eng[Dim];
          for (unsigned i = 0; i < Dim; ++i) {
              eng[i].seed(42 + i * Dim);
              eng[i].discard(nloc * ippl::Comm->rank());
          }

          std::mt19937_64 engN[4 * Dim];
          for (unsigned i = 0; i < 4 * Dim; ++i) {
              engN[i].seed(42 + i * Dim);
              engN[i].discard(nloc * ippl::Comm->rank());
          }

          auto dom                = fl.getDomain();
          unsigned int gridpoints = 1;
          for (unsigned d = 0; d < Dim; d++) {
              gridpoints *= dom[d].length();
          }
          if (tag == 0 && nloc * ippl::Comm->size() != gridpoints) {
              if (ippl::Comm->rank() == 0) {
                  std::cerr << "Particle count must match gridpoint count to use gridpoint "
                               "locations. Switching to uniform distribution."
                            << std::endl;
              }
              tag = 2;
          }

          if (tag == 0) {
              m << "Positions are set on grid points" << endl;
              int N = fl.getDomain()[0].length();  // this only works for boxes
              const ippl::NDIndex<Dim>& lDom = fl.getLocalNDIndex();
              int size                       = ippl::Comm->size();
              using index_type               = typename ippl::RangePolicy<Dim>::index_type;
              Kokkos::Array<index_type, Dim> begin, end;
              for (unsigned d = 0; d < Dim; d++) {
                  begin[d] = 0;
                  end[d]   = N;
              }
              end[0] /= size;
              // Loops over particles
              using index_array_type = typename ippl::RangePolicy<Dim>::index_array_type;
              ippl::parallel_for(
                  "initPositions", ippl::createRangePolicy(begin, end),
                  KOKKOS_LAMBDA(const index_array_type& args) {
                      int l = 0;
                      for (unsigned d1 = 0; d1 < Dim; d1++) {
                          int next = args[d1];
                          for (unsigned d2 = 0; d2 < d1; d2++) {
                              next *= N;
                          }
                          l += next / size;
                      }
                      R_host(l) = (0.5 + args + lDom.first()) * hr;
                  });

          } else if (tag == 1) {
              m << "Positions follow normal distribution" << endl;
              std::vector<double> mu = {0.5, 0.6, 0.2, 0.5, 0.6, 0.2};
              std::vector<double> sd = {0.75, 0.3, 0.2, 0.75, 0.3, 0.2};
              std::vector<double> states(Dim);

              Vector_t<T, Dim> length = 1;

              std::uniform_real_distribution<double> dist_uniform(0.0, 1.0);

              double sum_coord = 0.0;
              for (unsigned long long int i = 0; i < nloc; i++) {
                  for (unsigned d = 0; d < Dim; d++) {
                      double u1 = dist_uniform(engN[d * 2]);
                      double u2 = dist_uniform(engN[d * 2 + 1]);
                      states[d] = sd[d] * std::sqrt(-2.0 * std::log(u1))
                                      * std::cos(2.0 * Kokkos::numbers::pi * u2)
                                  + mu[d];
                      R_host(i)[d] = std::fabs(std::fmod(states[d], length[d]));
                      sum_coord += R_host(i)[d];
                  }
              }
          } else {
              double rmin = 0.0, rmax = 1.0;
              m << "Positions follow uniform distribution U(" << rmin << "," << rmax << ")" << endl;
              std::uniform_real_distribution<double> unif(rmin, rmax);
              for (unsigned long int i = 0; i < nloc; i++) {
                  for (unsigned d = 0; d < Dim; d++) {
                      R_host(i)[d] = unif(eng[d]);
                  }
              }
          }

          // Copy to device
          Kokkos::deep_copy(this->R.getView(), R_host);
      }
      */
};

template <class PLayout, typename T, unsigned dim>
inline double PartBunch<PLayout, T, dim>::getRho(int x, int y, int z) {
    // return rho_m[x][y][z].get();
    return x + y + z;
}

// inline
/*
template <class PLayout, typename T, unsigned dim>
inline Inform& operator<<(Inform& os, PartBunch<PLayout, T, dim>& p) {
    return p.print(os);
}
*/
/*

  FROM PartBunchBase.h

*/

#include "PartBunch.cpp"

#endif
