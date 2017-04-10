#ifndef OPAL_PartBunch_HH
#define OPAL_PartBunch_HH

// ------------------------------------------------------------------------
// $RCSfile: PartBunch.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class PartBunch
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: Andreas Adelmann  and Co. $
//
// ------------------------------------------------------------------------

#include "Ippl.h"
#include "Algorithms/PBunchDefs.h"
#include "Algorithms/Particle.h"
#include "Algorithms/CoordinateSystemTrafo.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include "Algorithms/PartBins.h"
#include "Algorithms/PartBinsCyc.h"
#include "Algorithms/PartData.h"
#include "Algorithms/Quaternion.h"
#include "Utilities/SwitcherError.h"
#include "Physics/Physics.h"

#include <iosfwd>
#include <vector>


class Distribution;
class LossDataSink;
class FieldSolver;
class ListElem;

template <class T, int, int> class FMatrix;
template <class T, int> class FVector;

// Class PartBunch.
// ------------------------------------------------------------------------
/// Particle Bunch.
//  A representation of a particle bunch as a vector of particles.

// class PartBunch: public std::vector<Particle>, public IpplParticleBase< ParticleSpatialLayout<double, 3> > {
class PartBunch: public IpplParticleBase< ParticleSpatialLayout<double, 3> > {

public:
    /// Default constructor.
    //  Construct empty bunch.
    PartBunch(const PartData *ref);

    /// Conversion.
    PartBunch(const std::vector<Particle> &, const PartData *ref);

    PartBunch(const PartBunch &);
    ~PartBunch();

    void runTests();

    bool getIfBeamEmitting();
    int getLastEmittedEnergyBin();
    size_t getNumberOfEmissionSteps();
    int getNumberOfEnergyBins();
    void Rebin();
    void setEnergyBins(int numberOfEnergyBins);
    bool weHaveEnergyBins();

    enum UnitState_t { units = 0, unitless = 1 };

    //FIXME: unify methods, use convention that all particles have own dt
    void switchToUnitlessPositions(bool use_dt_per_particle = false);

    //FIXME: unify methods, use convention that all particles have own dt
    void switchOffUnitlessPositions(bool use_dt_per_particle = false);

    /** \brief After each Schottky scan we delete
        all the particles.

    */
    void cleanUpParticles();

    void resetIfScan();

    double getRho(int x, int y, int z);

    void do_binaryRepart();

    void calcLineDensity(unsigned int nBins, std::vector<double> &lineDensity, std::pair<double, double> &meshInfo);

    void setDistribution(Distribution *d,
                         std::vector<Distribution *> addedDistributions,
                         size_t &np,
                         bool scan);

    bool isGridFixed();

    /*

      Energy bins related functions

    */

    void setTEmission(double t);
    double getTEmission();
    bool doEmission();

    bool weHaveBins() const;

    void setPBins(PartBins *pbin);
    void setPBins(PartBinsCyc *pbin);

    /** \brief Emit particles in the given bin
        i.e. copy the particles from the bin structure into the
        particle container
    */
    size_t emitParticles(double eZ);

    void updateNumTotal();

    void rebin();

    int getNumBins();

    int getLastemittedBin();

    /** \brief Compute the gammas of all bins */
    void calcGammas();

    void calcGammas_cycl();

    /** \brief Get gamma of one bin */
    double getBinGamma(int bin);

    /** \brief Set the charge of one bin to the value of q and all other to zero */
    void setBinCharge(int bin, double q);

    /** \brief Set the charge of all other the ones in bin to zero */
    void setBinCharge(int bin);

    /** \brief calculates back the max/min of the efield on the grid */
    std::pair<Vector_t, Vector_t> getEExtrema();

    /** \brief returns the number of particles outside of a box defined by x */
    size_t calcNumPartsOutside(Vector_t x);

    /*

      Mesh and Field Layout related functions

    */

    const Mesh_t &getMesh() const;

    Mesh_t &getMesh();

    FieldLayout_t &getFieldLayout();

    void setBCAllPeriodic();
    void setBCAllOpen();

    void setBCForDCBeam();

    void boundp();

    /** delete particles which are too far away from the center of beam*/
    void boundp_destroy();

    /** This is only temporary in order to get the collimator and pepperpot workinh */
    size_t boundp_destroyT();
    size_t destroyT();

    /* Definiere virtuelle Funktionen, um die Koordinaten auszulesen
     *     */

    virtual double getPx(int i);
    virtual double getPy(int i);
    virtual double getPz(int i);

    virtual double getPx0(int i);
    virtual double getPy0(int i);

    virtual double getX(int i);
    virtual double getY(int i);
    virtual double getZ(int i);

    virtual double getX0(int i);
    virtual double getY0(int i);

    virtual void setZ(int i, double zcoo);

    void get_bounds(Vector_t &rmin, Vector_t &rmax);
    void getLocalBounds(Vector_t &rmin, Vector_t &rmax);
    std::pair<Vector_t, double> getBoundingSphere();
    std::pair<Vector_t, double> getLocalBoundingSphere();

    /*
      Compatibility function push_back

    */

    void push_back(Particle p);

    void set_part(FVector<double, 6> z, int ii);

    void set_part(Particle p, int ii);

    Particle get_part(int ii);

    /// Return maximum amplitudes.
    //  The matrix [b]D[/b] is used to normalise the first two modes.
    //  The maximum normalised amplitudes for these modes are stored
    //  in [b]axmax[/b] and [b]aymax[/b].
    void maximumAmplitudes(const FMatrix<double, 6, 6> &D,
                           double &axmax, double &aymax);

    Inform &print(Inform &os);


    void   setdT(double dt);
    double getdT() const;

    void   setT(double t);
    double getT() const;
    void   computeSelfFields();

    /** /brief used for self fields with binned distribution */
    void   computeSelfFields(int b);

    void computeSelfFields_cycl(double gamma);
    void computeSelfFields_cycl(int b);

    void resetInterpolationCache(bool clearCache = false);

    /**
     * get the spos of the primary beam
     *
     * @param none
     *
     */
    double get_sPos();
    void set_sPos(double s);

    double   get_gamma() const;

    double get_meanKineticEnergy() const;
    Vector_t get_origin() const;
    Vector_t get_maxExtent() const;
    Vector_t get_centroid() const;
    Vector_t get_rrms() const;
    Vector_t get_rmean() const;
    Vector_t get_prms() const;
    Vector_t get_rprms() const;
    Vector_t get_pmean() const;
    Vector_t get_pmean_Distribution() const;
    Vector_t get_emit() const;
    Vector_t get_norm_emit() const;
    Vector_t get_hr() const;

    double get_Dx() const;
    double get_Dy() const;

    double get_DDx() const;
    double get_DDy() const;

    void set_meshEnlargement(double dh);

    void gatherLoadBalanceStatistics();
    size_t getLoadBalance(int p) const;


    void get_PBounds(Vector_t &min, Vector_t &max) const;

    void calcBeamParameters();
    void calcBeamParametersInitial(); // Calculate initial beam parameters before emission.

    double getCouplingConstant() const;
    void setCouplingConstant(double c);

    // set the charge per simulation particle
    void setCharge(double q);
    // set the charge per simulation particle when total particle number equals 0
    void setChargeZeroPart(double q);

    // set the mass per simulation particle
    void setMass(double mass);


    /// \brief Need Ek for the Schottky effect calculation (eV)
    double getEkin() const;

    /// Need the work function for the Schottky effect calculation (eV)
    double getWorkFunctionRf() const;

    /// Need the laser energy for the Schottky effect calculation (eV)
    double getLaserEnergy() const;

    /// get the total charge per simulation particle
    double getCharge() const;

    /// get the macro particle charge
    double getChargePerParticle() const;

    void setSolver(FieldSolver *fs);

    bool hasFieldSolver();

    std::string getFieldSolverType() const;

    void setLPath(double s);
    double getLPath() const;

    void setStepsPerTurn(int n);
    int getStepsPerTurn() const;

    /// step in multiple TRACK commands
    void setGlobalTrackStep(long long n);
    long long getGlobalTrackStep() const;

    /// step in a TRACK command
    void setLocalTrackStep(long long n);
    void incTrackSteps();
    long long getLocalTrackStep() const;

    void setNumBunch(int n);
    int getNumBunch() const;

    void setGlobalMeanR(Vector_t globalMeanR);
    Vector_t getGlobalMeanR();
    void setGlobalToLocalQuaternion(Quaternion_t globalToLocalQuaternion);
    Quaternion_t getGlobalToLocalQuaternion();

    void setSteptoLastInj(int n);
    int getSteptoLastInj();

    /// calculate average angle of longitudinal direction of bins
    double calcMeanPhi();

    /// reset Bin[] for each particle according to the method given in paper PAST-AB(064402) by  G. Fubiani et al.
    bool resetPartBinID2(const double eta);


    double getQ() const;
    double getM() const;
    double getP() const;
    double getE() const;

    void resetQ(double q);
    void resetM(double m);

    double getdE();
    double getInitialBeta() const;
    double getInitialGamma() const;
    virtual double getGamma(int i);
    virtual double getBeta(int i);
    virtual void actT();
    const PartData *getReference() const;

    double getEmissionDeltaT();

    Quaternion_t getQKs3D();
    void         setQKs3D(Quaternion_t q);
    Vector_t     getKs3DRefr();
    void         setKs3DRefr(Vector_t r);
    Vector_t     getKs3DRefp();
    void         setKs3DRefp(Vector_t p);

    void iterateEmittedBin(int binNumber);

    void calcEMean();
    void correctEnergy(double avrgp);

    void swap(unsigned int i, unsigned int j);

    // Particle container attributes
    ParticleAttrib< Vector_t > P;      // particle momentum //  ParticleSpatialLayout<double, 3>::ParticlePos_t P;
    ParticleAttrib< double >   Q;      // charge per simulation particle, unit: C.
    ParticleAttrib< double >   M;      // mass per simulation particle, for multi-species particle tracking, unit:GeV/c^2.
    ParticleAttrib< double >   Phi;    // the electric potential
    ParticleAttrib< Vector_t > Ef;     // e field vector
    ParticleAttrib< Vector_t > Eftmp;  // e field vector for gun simulations

    ParticleAttrib< Vector_t > Bf;   // b field vector
    ParticleAttrib< int >      Bin;   // holds the bin in which the particle is in, if zero particle is marked for deletion
    ParticleAttrib< double >   dt;   // holds the dt timestep for particle

    ParticleAttrib< short >    PType; // we can distinguish dark current particles from primary particle
    ParticleAttrib< int >      TriID; // holds the ID of triangle that the particle hit. Only for BoundaryGeometry case.

    Vector_t RefPartR_m;
    Vector_t RefPartP_m;
    CoordinateSystemTrafo toLabTrafo_m;

    /// scalar potential
    Field_t rho_m;

    /// vector field on the grid
    VField_t  eg_m;

    /// avoid calls to Ippl::myNode()
    int myNode_m;

    /// avoid calls to Ippl::getNodes()
    int nodes_m;

    /// if the grid does not have to adapt
    bool fixed_grid;

    // The structure for particle binning
    PartBins *pbin_m;

    std::unique_ptr<LossDataSink> lossDs_m;

    // save particles in case of one core
    std::unique_ptr<Inform> pmsg_m;
    std::unique_ptr<std::ofstream> f_stream;

    /// timer for IC, can not be in Distribution.h
    IpplTimings::TimerRef distrReload_m;
    IpplTimings::TimerRef distrCreate_m;

    // For AMTS integrator in OPAL-T
    double dtScInit_m, deltaTau_m;

    /// if a local node has less than 2 particles  lowParticleCount_m == true
    bool lowParticleCount_m;

protected:
    /// timer for selfField calculation
    IpplTimings::TimerRef selfFieldTimer_m;
    IpplTimings::TimerRef compPotenTimer_m;
    IpplTimings::TimerRef boundpTimer_m;
    IpplTimings::TimerRef statParamTimer_m;

    IpplTimings::TimerRef histoTimer_m;

private:

    const PartData *reference;
    void calcMoments();    // Calculates bunch moments using only emitted particles.
    void calcMomentsInitial(); // Calcualtes bunch moments by summing over bins (not accurate when any particles have been emitted).

    double calculateAngle(double x, double y);
    double calculateAngle2(double x, double y);

    /*
      Member variables starts here
    */

    // unit state of PartBunch
    UnitState_t unit_state_;
    UnitState_t stateOfLastBoundP_;

    /// holds the centroid of the beam
    double centroid_m[2 * Dim];

    /// resize mesh to geometry specified
    void resizeMesh();

    /// 6x6 matrix of the moments of the beam
    FMatrix<double, 2 * Dim, 2 * Dim> moments_m;

    /// holds the timestep in seconds
    double dt_m;
    /// holds the actual time of the integration
    double t_m;
    /// mean energy of the bunch (MeV)
    double eKin_m;
    /// energy spread of the beam in keV
    double dE_m;
    /// the position along design trajectory
    double spos_m;

    /// Initialize the translation vector and rotation quaternion
    /// here. Cyclotron tracker will reset these values each timestep
    /// TTracker can just use 0 translation and 0 rotation (quat[1 0 0 0]).
    //Vector_t globalMeanR_m = Vector_t(0.0, 0.0, 0.0);
    //Quaternion_t globalToLocalQuaternion_m = Quaternion_t(1.0, 0.0, 0.0, 0.0);
    Vector_t globalMeanR_m;
    Quaternion_t globalToLocalQuaternion_m;

    /// for coordinate transformation to Ks
    /// Ks is the coordinate system to calculate statistics
    Quaternion_t QKs3D_m;
    /// holds the referernce particle
    Vector_t     Ks3DRefr_m;
    Vector_t     Ks3DRefp_m;

    /// maximal extend of particles
    Vector_t rmax_m;
    /// minimal extend of particles
    Vector_t rmin_m;

    /// rms beam size (m)
    Vector_t rrms_m;
    /// rms momenta
    Vector_t prms_m;
    /// mean position (m)
    Vector_t rmean_m;
    /// mean momenta
    Vector_t pmean_m;

    /// rms emittance (not normalized)
    Vector_t eps_m;

    /// rms normalized emittance
    Vector_t eps_norm_m;
    /// rms correlation
    Vector_t rprms_m;

    /// dispersion x & y
    double Dx_m;
    double Dy_m;

    /// derivative of the dispersion
    double DDx_m;
    double DDy_m;

    /// meshspacing of cartesian mesh
    Vector_t hr_m;
    /// meshsize of cartesian mesh
    Vektor<int, 3> nr_m;

    /// for defining the boundary conditions
    BConds<double, 3, Mesh_t, Center_t> bc_m;
    BConds<Vector_t, 3, Mesh_t, Center_t> vbc_m;

    /// stores the used field solver
    FieldSolver *fs_m;

    double couplingConstant_m;

    double qi_m;

    bool interpolationCacheSet_m;

    ParticleAttrib<CacheDataCIC<double, 3U> > interpolationCache_m;

    /// counter to store the distributin dump
    int distDump_m;

    PartBunch &operator=(const PartBunch &) = delete;

    ///
    int fieldDBGStep_m;

    /// Mesh enlargement
    double dh_m; /// in % how much the mesh is enlarged

    /// if larger than 0, emitt particles for tEmission_m [s]
    double tEmission_m;

    /// holds the gamma of the bin
    std::unique_ptr<double[]> bingamma_m;

    //FIXME: this should go into the Bin class!
    // holds number of emitted particles of the bin
    // jjyang: opal-cycl use *nBin_m of pbin_m
    std::unique_ptr<size_t[]> binemitted_m;

    /// path length from the start point
    double lPath_m;

    /// steps per turn for OPAL-cycl
    int stepsPerTurn_m;

    /// step in a TRACK command
    long long localTrackStep_m;

    /// if multiple TRACK commands
    long long globalTrackStep_m;

    /// current bunch number
    int numBunch_m;

    /// this parameter records the current steps since last bunch injection
    /// it helps to inject new bunches correctly in the restart run of OPAL-cycl
    /// it is stored during phase space dump.
    int SteptoLastInj_m;

    /*
      Data structure for particle load balance information
    */

    std::unique_ptr<size_t[]> partPerNode_m;
    std::unique_ptr<size_t[]> globalPartPerNode_m;
    size_t minLocNum_m;


    Distribution *dist_m;

    // flag to tell if we are a DC-beam
    bool dcBeam_m;


};


inline Quaternion_t PartBunch::getQKs3D() {
    return QKs3D_m;
}

inline void PartBunch::setQKs3D(Quaternion_t q) {
    QKs3D_m=q;
}

inline Vector_t PartBunch::getKs3DRefr() {
    return Ks3DRefr_m;
}

inline void PartBunch::setKs3DRefr(Vector_t r) {
    Ks3DRefr_m=r;
}

inline Vector_t PartBunch::getKs3DRefp() {
    return Ks3DRefp_m;
}

inline void PartBunch::setKs3DRefp(Vector_t p) {
    Ks3DRefp_m=p;
}

inline
void PartBunch::switchToUnitlessPositions(bool use_dt_per_particle) {

    if(unit_state_ == unitless)
        throw SwitcherError("PartBunch::switchToUnitlessPositions",
                            "Cannot make a unitless PartBunch unitless");

    bool hasToReset = false;
    if(!R.isDirty()) hasToReset = true;

    for(size_t i = 0; i < getLocalNum(); i++) {
        double dt = getdT();
        if(use_dt_per_particle)
            dt = this->dt[i];

        R[i] /= Vector_t(Physics::c * dt);
    }

    unit_state_ = unitless;

    if(hasToReset) R.resetDirtyFlag();
}

//FIXME: unify methods, use convention that all particles have own dt
inline
void PartBunch::switchOffUnitlessPositions(bool use_dt_per_particle) {

    if(unit_state_ == units)
        throw SwitcherError("PartBunch::switchOffUnitlessPositions",
                            "Cannot apply units twice to PartBunch");

    bool hasToReset = false;
    if(!R.isDirty()) hasToReset = true;

    for(size_t i = 0; i < getLocalNum(); i++) {
        double dt = getdT();
        if(use_dt_per_particle)
            dt = this->dt[i];

        R[i] *= Vector_t(Physics::c * dt);
    }

    unit_state_ = units;

    if(hasToReset) R.resetDirtyFlag();
}

inline
double PartBunch::getRho(int x, int y, int z) {
    return rho_m[x][y][z].get();
}

inline
void PartBunch::do_binaryRepart() {
    get_bounds(rmin_m, rmax_m);
    BinaryRepartition(*this);
    update();
    get_bounds(rmin_m, rmax_m);
    boundp();
}

inline
bool PartBunch::isGridFixed() {
    return fixed_grid;
}

inline
void   PartBunch::setTEmission(double t) {
    tEmission_m = t;
}

inline
double PartBunch::getTEmission() {
    return tEmission_m;
}

inline
void PartBunch::rebin() {
    this->Bin = 0;
    pbin_m->resetBins();
    // delete pbin_m; we did not allocate it!
    pbin_m = NULL;
}

inline
int PartBunch::getNumBins() {
    if(pbin_m != NULL)
        return pbin_m->getNBins();
    else
        return 0;
}

inline
int PartBunch::getLastemittedBin() {
    if(pbin_m != NULL)
        return pbin_m->getLastemittedBin();
    else
        return 0;
}

inline
double PartBunch::getBinGamma(int bin) {
    return bingamma_m[bin];
}

inline
void PartBunch::setBinCharge(int bin, double q) {
    this->Q = where(eq(this->Bin, bin), q, 0.0);
}

inline
void PartBunch::setBinCharge(int bin) {
    this->Q = where(eq(this->Bin, bin), this->Q, 0.0);
}

inline
std::pair<Vector_t, Vector_t> PartBunch::getEExtrema() {
    const Vector_t maxE = max(eg_m);
    //      const double maxL = max(dot(eg_m,eg_m));
    const Vector_t minE = min(eg_m);
    // INFOMSG("MaxE= " << maxE << " MinE= " << minE << endl);
    return std::pair<Vector_t, Vector_t>(maxE, minE);
}

inline
const Mesh_t &PartBunch::getMesh() const {
    return getLayout().getLayout().getMesh();
}

inline
Mesh_t &PartBunch::getMesh() {
    return getLayout().getLayout().getMesh();
}

inline
FieldLayout_t &PartBunch::getFieldLayout() {
    return dynamic_cast<FieldLayout_t &>(getLayout().getLayout().getFieldLayout());
}

inline
double PartBunch::getPx0(int) {
    return 0;
}

inline
double PartBunch::getPy0(int) {
    return 0;
}

inline
void PartBunch::setZ(int /* i */, double /* zcoo */) {};

inline
void PartBunch::get_bounds(Vector_t &rmin, Vector_t &rmax) {
    bounds(this->R, rmin, rmax);
}

inline
std::pair<Vector_t, double> PartBunch::getBoundingSphere() {
    Vector_t rmin, rmax;
    get_bounds(rmin, rmax);

    std::pair<Vector_t, double> sphere;
    sphere.first = 0.5 * (rmin + rmax);
    sphere.second = sqrt(dot(rmax - sphere.first, rmax - sphere.first));

    return sphere;
}

inline
std::pair<Vector_t, double> PartBunch::getLocalBoundingSphere() {
    Vector_t rmin, rmax;
    getLocalBounds(rmin, rmax);

    std::pair<Vector_t, double> sphere;
    sphere.first = 0.5 * (rmin + rmax);
    sphere.second = sqrt(dot(rmax - sphere.first, rmax - sphere.first));

    return sphere;
}

inline
void PartBunch::push_back(Particle p) {
    Inform msg("PartBunch ");

    create(1);
    size_t i = getTotalNum();

    R[i](0) = p[0];
    R[i](1) = p[1];
    R[i](2) = p[2];

    P[i](0) = p[3];
    P[i](1) = p[4];
    P[i](2) = p[5];

    update();
    msg << "Created one particle i= " << i << endl;
}

inline
void PartBunch::set_part(FVector<double, 6> z, int ii) {
    R[ii](0) = z[0];
    P[ii](0) = z[1];
    R[ii](1) = z[2];
    P[ii](1) = z[3];
    R[ii](2) = z[4];
    P[ii](2) = z[5];
}

inline
void PartBunch::set_part(Particle p, int ii) {
    R[ii](0) = p[0];
    P[ii](0) = p[1];
    R[ii](1) = p[2];
    P[ii](1) = p[3];
    R[ii](2) = p[4];
    P[ii](2) = p[5];
}

inline
Particle PartBunch::get_part(int ii) {
    Particle part;
    part[0] = R[ii](0);
    part[1] = P[ii](0);
    part[2] = R[ii](1);
    part[3] = P[ii](1);
    part[4] = R[ii](2);
    part[5] = P[ii](2);
    return part;
}

inline
void   PartBunch::setdT(double dt) {
    dt_m = dt;
}

inline
double PartBunch::getdT() const {
    return dt_m;
}

inline
void   PartBunch::setT(double t) {
    t_m = t;
}

inline
double PartBunch::getT() const {
    return t_m;
}

inline
void PartBunch::resetInterpolationCache(bool clearCache) {
    interpolationCacheSet_m = false;
    if(clearCache) {
        interpolationCache_m.destroy(interpolationCache_m.size(), 0, true);
    }
}

inline
double PartBunch::get_sPos() {
    if(sum(PType != ParticleType::REGULAR)) {
        const size_t n = getLocalNum();
        size_t numPrimBeamParts = 0;
        double z = 0.0;
        if(n != 0) {
            for(size_t i = 0; i < n; i++) {
                if(PType[i] == ParticleType::REGULAR) {
                    z += R[i](2);
                    numPrimBeamParts++;
                }
            }
        }
        reduce(z, z, OpAddAssign());
	reduce(numPrimBeamParts, numPrimBeamParts, OpAddAssign());
	if(numPrimBeamParts != 0)
            z = z / numPrimBeamParts;
        return z;
    } else {
        return spos_m;
    }
}

inline
void PartBunch::set_sPos(double s) {
    spos_m = s;
}


inline
double   PartBunch::get_gamma() const {
    return eKin_m / (getM()*1e-6) + 1.0;
}


inline
double PartBunch::get_meanKineticEnergy() const {
    return eKin_m;
}

inline
Vector_t PartBunch::get_origin() const {
    return rmin_m;
}

inline
Vector_t PartBunch::get_maxExtent() const {
    return rmax_m;
}

inline
Vector_t PartBunch::get_centroid() const {
    return rmean_m;
}

inline
Vector_t PartBunch::get_rrms() const {
    return rrms_m;
}

inline
Vector_t PartBunch::get_rprms() const {
    return rprms_m;
}

inline
Vector_t PartBunch::get_rmean() const {
    return rmean_m;
}

inline
Vector_t PartBunch::get_prms() const {
    return prms_m;
}

inline
Vector_t PartBunch::get_pmean() const {
    return pmean_m;
}

inline
Vector_t PartBunch::get_emit() const {
    return eps_m;
}

inline
Vector_t PartBunch::get_norm_emit() const {
    return eps_norm_m;
}

inline
Vector_t PartBunch::get_hr() const {
    return hr_m;
}

inline
double PartBunch::get_Dx() const {
    return Dx_m;
}

inline
double PartBunch::get_Dy() const {
    return Dy_m;
}

inline
double PartBunch::get_DDx() const {
    return DDx_m;
}

inline
double PartBunch::get_DDy() const {
    return DDy_m;
}

inline
void PartBunch::set_meshEnlargement(double dh) {
    dh_m = dh;
}

inline
size_t PartBunch::getLoadBalance(int p) const {
    return globalPartPerNode_m[p];
}

inline
void PartBunch::get_PBounds(Vector_t &min, Vector_t &max) const {
    bounds(this->P, min, max);
}

inline
double PartBunch::getCouplingConstant() const {
    return couplingConstant_m;
}

inline
void PartBunch::setCouplingConstant(double c) {
    couplingConstant_m = c;
}

inline
void PartBunch::setCharge(double q) {
    qi_m = q;
    if(getTotalNum() != 0)
        Q = qi_m;
    else
        WARNMSG("Could not set total charge in PartBunch::setCharge based on getTotalNum" << endl);
}

inline
void PartBunch::setChargeZeroPart(double q) {
    qi_m = q;
}

inline
void PartBunch::setMass(double mass) {
    M = mass;
}

inline
double PartBunch::getCharge() const {
    return sum(Q);
}

inline
double PartBunch::getChargePerParticle() const {
    return qi_m;
}

inline
void PartBunch::setLPath(double s) {
    lPath_m = s;
}

inline
double PartBunch::getLPath() const {
    return lPath_m;
}

inline
void PartBunch::setStepsPerTurn(int n) {
    stepsPerTurn_m = n;
}

inline
int PartBunch::getStepsPerTurn() const {
    return stepsPerTurn_m;
}

inline
void PartBunch::setGlobalTrackStep(long long n) {
    globalTrackStep_m = n;
}

inline
long long PartBunch::getGlobalTrackStep() const {
    return globalTrackStep_m;
}

inline
void PartBunch::setLocalTrackStep(long long n) {
    localTrackStep_m = n;
}

inline
void PartBunch::incTrackSteps() {
    globalTrackStep_m++; localTrackStep_m++;
}

inline
long long PartBunch::getLocalTrackStep() const {
    return localTrackStep_m;
}

inline
void PartBunch::setNumBunch(int n) {
    numBunch_m = n;
}

inline
int PartBunch::getNumBunch() const {
    return numBunch_m;
}

inline
void PartBunch::setGlobalMeanR(Vector_t globalMeanR) {
    globalMeanR_m = globalMeanR;
}

inline
Vector_t PartBunch::getGlobalMeanR() {
    return globalMeanR_m;
}

inline
void PartBunch::setGlobalToLocalQuaternion(Quaternion_t globalToLocalQuaternion) {

    globalToLocalQuaternion_m = globalToLocalQuaternion;
}

inline
Quaternion_t PartBunch::getGlobalToLocalQuaternion() {
    return globalToLocalQuaternion_m;
}

inline
void PartBunch::setSteptoLastInj(int n) {
    SteptoLastInj_m = n;
}

inline
int PartBunch::getSteptoLastInj() {
    return SteptoLastInj_m;
}

inline
double PartBunch::getQ() const {
    return reference->getQ();
}

inline
double PartBunch::getM() const {
    return reference->getM();
}

inline
double PartBunch::getP() const {
    return reference->getP();
}

inline
double PartBunch::getE() const {
    return reference->getE();
}

inline
void PartBunch::resetQ(double q)  {
    const_cast<PartData *>(reference)->setQ(q);
}

inline
void PartBunch::resetM(double m)  {
    const_cast<PartData *>(reference)->setM(m);
}

inline
double PartBunch::getdE() {
    return dE_m;
}

inline
double PartBunch::getInitialBeta() const {
    return reference->getBeta();
}

inline
double PartBunch::getInitialGamma() const {
    return reference->getGamma();
}

inline
double PartBunch::getGamma(int) {
    return 0;
}

inline
double PartBunch::getBeta(int) {
    return 0;
}

inline
void PartBunch::actT() {};

inline
const PartData *PartBunch::getReference() const {
    return reference;
}

inline
Inform &operator<<(Inform &os, PartBunch &p) {
    return p.print(os);
}


#endif // OPAL_PartBunch_HH