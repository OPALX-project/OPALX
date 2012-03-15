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

#include "Algorithms/Particle.h"
#include "FixedAlgebra/FMatrix.h"
#include "FixedAlgebra/FVector.h"
#include "FixedAlgebra/FVps.h"
#include "FixedAlgebra/LinearMap.h"
#include "Algorithms/PartBins.h"
#include "Algorithms/PartBinsCyc.h"
#include "Algorithms/PartData.h"

#include <iosfwd>
#include <vector>
#include "Ippl.h"

class PartBunch;
class Distribution;
class LossDataSink;
#include "Structure/FieldSolver.h"

#include "ListElem.h"

template <class T, int, int> class FMatrix;
template <class T, int> class FVector;

typedef IntCIC  IntrplCIC_t;
typedef IntNGP  IntrplNGP_t;
typedef IntSUDS IntrplSUDS_t;

typedef ParticleSpatialLayout<double, 3>::ParticlePos_t Ppos_t;
typedef ParticleSpatialLayout<double, 3>::ParticleIndex_t PID_t;

typedef ParticleAttrib<double> Pscalar_t;

typedef InterpolatorTraits<double, 3, IntrplCIC_t>::Cache_t Pcache_t;

typedef UniformCartesian<3, double> Mesh_t;

typedef ParticleSpatialLayout<double, 3>::SingleParticlePos_t Vector_t;

typedef ParticleSpatialLayout< double, 3, Mesh_t  > Layout_t;

typedef Cell Center_t;

typedef CenteredFieldLayout<3, Mesh_t, Center_t> FieldLayout_t;
typedef Field<double, 3, Mesh_t, Center_t>       Field_t;
typedef Field<Vector_t, 3, Mesh_t, Center_t>     VField_t;


// Class PartBunch.
// ------------------------------------------------------------------------
/// Particle Bunch.
//  A representation of a particle bunch as a vector of particles.

// class PartBunch: public std::vector<Particle>, public ParticleBase< ParticleSpatialLayout<double, 3> > {
class PartBunch: public ParticleBase< ParticleSpatialLayout<double, 3> > {

public:

    // Particle container attributes
    ParticleAttrib< Vector_t > X;      // local 'lab frame' coordinates;
    ParticleAttrib< Vector_t > P;      // particle momentum //  ParticleSpatialLayout<double, 3>::ParticlePos_t P;
    ParticleAttrib< double >   Q;        // charge
    ParticleAttrib< Vector_t > Ef;     // e field vector
    ParticleAttrib< Vector_t > Eftmp;  // e field vector for gun simulations

    ParticleAttrib< Vector_t > Bf;   // b field vector
    ParticleAttrib< int >      Bin;   // holds the bin in which the particle is in, if zero particle is marked for deletion
    ParticleAttrib< double >   dt;   // holds the dt timestep for particle
    ParticleAttrib< long >     LastSection; // last em-field section


    ParticleAttrib< short >    PType; // we can distinguish dark current particles from primary particle
    ParticleAttrib< int >      TriID; // holds the ID of triangle that the particle hit. Only for BoundaryGeometry case.


    Vector_t RefPart_R;
    Vector_t RefPart_P;

    /// scalar potential
    Field_t rho_m;

    /// scalar fields for projecttion i.e. line densities
    Field_t tmpFieldZ_m;

    /// vector field on the grid
    VField_t  eg_m;

    /// avoid calls to Ippl::myNode()
    int myNode_m;

    /// avoid calls to Ippl::getNodes()
    int nodes_m;

    /// if the grid does not have to adapt
    bool fixed_grid;


    /// Default constructor.
    //  Construct empty bunch.
    PartBunch(const PartData *ref);

    /// Conversion.
    PartBunch(const std::vector<Particle> &, const PartData *ref);

    PartBunch(const PartBunch &);
    ~PartBunch();


    void resetIfScan();

    inline double getRho(NDIndex<3> e) { return rho_m.localElement(e); }

    inline double getRho(int x, int y, int z) { return rho_m[x][y][z].get(); }

    inline bool itIsMyTurn(int *n) {
        bool res = (*n == myNode_m);
        *n++;
        if(*n == nodes_m)
            *n = 0;
        return res;
    }

    bool hasZeroNLP();

    inline void do_binaryRepart() {
        get_bounds(rmin_m, rmax_m);
        BinaryRepartition(*this);
        update();
        get_bounds(rmin_m, rmax_m);
        boundp();
    }

    /// per default the MT value of the field solver is used
    void set_nBinsLineDensity(int n) {nBinsLineDensity_m = n;}

    void calcLineDensity();
    void fillArray(double *lineDensity, const list<ListElem> &l);
    void getLineDensity(vector<double> &lineDensity);

    void setDistribution(Distribution *d, size_t np, bool scan);


    inline void setGridIsFixed() { fixed_grid = true; }
    inline bool isGridFixed() { return fixed_grid; }


    /*

    Energy bins related functions

    */

    void   setTEmission(double t) {tEmission_m = t;}
    double getTEmission() {return tEmission_m;}
    bool doEmission() {return (tEmission_m > 0.0);}

    bool weHaveBins() const {
        if(pbin_m)
            return pbin_m->weHaveBins();
        else
            return false;
    }

    double getRebinEnergy() {
        return pbin_m->getRebinEnergy();
    }

    void weHaveNOBins() {
        if(pbin_m)
            delete pbin_m;
    }


    void setPBins(PartBins *pbin);
    void setPBins(PartBinsCyc *pbin);

    /** \brief Emit particles in the given bin
        i.e. copy the particles from the bin structure into the
        particle container
    */

    size_t emitParticles();
    size_t emitParticlesOLD(int bin);

    double calcTimeDelay(const double &jifactor);
    void moveBunchToCathode(double &t);
    void printBinHist();

    void rebin() {
        this->Bin = 0;
        pbin_m->resetBins();
    }


    int getNumBins() {
        if(pbin_m)
            return pbin_m->getNBins();
        else
            return 0;
    }
    int getLastemittedBin() {
        if(pbin_m)
            return pbin_m->getLastemittedBin();
        else
            return 0;
    }


    void updatePartInBin(size_t countLost[]) {
        if(pbin_m)
            pbin_m->updatePartInBin(countLost);
    }

    /** \brief we need this because only node 0 is emitting */
    void updateBinStructure();

    /** \brief report if any particle has been lost */
    void reportParticleLoss();

    /** \brief Compute the gammas of all bins */
    void calcGammas();

    void calcGammas_cycl();


    /** \brief Get gamma of one bin */
    inline double getBinGamma(int bin) { return bingamma_m[bin]; }

    /** \brief Set the charge of one bin to the value of q and all other to zero*/
    void setBinCharge(int bin, double q) { this->Q = where(eq(this->Bin, bin), q, 0.0); }

    /** \brief gets back the maximum dE of all the bins */
    double getMaxdEBins();

    /*

    Mesh and Field Layout related functions

    */


    inline const Mesh_t &getMesh() const { return getLayout().getLayout().getMesh(); }

    inline Mesh_t &getMesh() { return getLayout().getLayout().getMesh(); }

    inline FieldLayout_t &getFieldLayout() {
        return dynamic_cast<FieldLayout_t &>(getLayout().getLayout().getFieldLayout());
    }

    inline void setBCAllOpen() {
        for(int i = 0; i < 2 * 3; ++i) {
            bc_m[i] = new ZeroFace<double, 3, Mesh_t, Center_t>(i);
            vbc_m[i] = new ZeroFace<Vector_t, 3, Mesh_t, Center_t>(i);
            getBConds()[i] = ParticleNoBCond;
        }
    }

    void boundp();
    void boundpNoRep();

    /** delete particles which are too far away from the center of beam*/
    void boundp_destroy();
    /** delete particles when hit stripper in multi-turn extraction cyclotron*/
    bool boundp_extraction(const double R_stripper, const double theta_stripper, const double delta_theta, const double time);
    /** This is only temporary in order to get the collimator and pepperpot workinh */
    size_t boundp_destroyT();

    /* Definiere virtuelle Funktionen, um die Koordinaten auszulesen
     *     */

    virtual double getPx(int i);
    virtual double getPy(int i);
    virtual double getPz(int i);

    virtual double getPx0(int i) {};
    virtual double getPy0(int i) {};

    virtual double getX(int i);
    virtual double getY(int i);
    virtual double getZ(int i);

    virtual double getX0(int i);
    virtual double getY0(int i);

    virtual void setZ(int i, double zcoo) {};

    virtual void setKR(Vector_t value, int i) {};

    virtual void setKT(Vector_t value, int i) {};

    virtual void BetOut(FILE *dat, FILE *sli) {};

    virtual void plotR() {};


    void get_bounds(Vector_t &rmin, Vector_t &rmax) {
        bounds(this->R, rmin, rmax);
    }


    /*
      Compatibility function push_back

    */

    inline void push_back(Particle p) {
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

    inline void set_part(FVector<double, 6> z, int ii) {
        R[ii](0) = z[0];
        P[ii](0) = z[1];
        R[ii](1) = z[2];
        P[ii](1) = z[3];
        R[ii](2) = z[4];
        P[ii](2) = z[5];
    }

    inline void set_part(Particle p, int ii) {
        R[ii](0) = p[0];
        P[ii](0) = p[1];
        R[ii](1) = p[2];
        P[ii](1) = p[3];
        R[ii](2) = p[4];
        P[ii](2) = p[5];
    }

    inline Particle get_part(int ii) {
        Particle part;
        part[0] = R[ii](0);
        part[1] = P[ii](0);
        part[2] = R[ii](1);
        part[3] = P[ii](1);
        part[4] = R[ii](2);
        part[5] = P[ii](2);
        return part;
    }

    /// Return bunch distribution.
    //  Return the bunch centroid in [b]centroid[/b],
    //  and the second moments in [b]moment[/b].
    void beamEllipsoid(FVector<double, 6>   &centroid,
                       FMatrix<double, 6, 6> &moment);

    /// Return maximum amplitudes.
    //  The matrix [b]D[/b] is used to normalise the first two modes.
    //  The maximum normalised amplitudes for these modes are stored
    //  in [b]axmax[/b] and [b]aymax[/b].
    void maximumAmplitudes(const FMatrix<double, 6, 6> &D,
                           double &axmax, double &aymax);

    Inform &print(Inform &os);


    void   setdT(double dt) { dt_m = dt; }
    double getdT() const { return dt_m; }

    void   setT(double t) { t_m = t; }
    double getT() const { return t_m; }
    void   computeSelfFields();

    /** /brief used for self fields with binned distribution */
    void computeSelfFields(int b);

    void computeSelfFields_cycl(double gamma);

    void computeSelfFields_cycl(int b);

    void calcWeightedAverages(Vector_t &CentroidPosition, Vector_t &CentroidMomentum) const;

    /** EulerAngle[0] = rotation about the y-axis, EulerAngle[1] = rotation about x-axis
     *  EulerAngle[2] = rotation about the y'-axis */

    void rotateAbout(const Vector_t &Center, const Vector_t &EulerAngles);

    void moveBy(const Vector_t &Center);

    void ResetLocalCoordinateSystem(const int &i, const Vector_t &Orientation, const double &origin);

    // OLD stuff should go away
    inline bool isZPeriodic() const { return false; } // used in Distribution
    inline double getGaBeLa() const { return 1.0; }   // used in Distribution

    /**
     * get the spos of the primary beam
     *
     * @param none
     *
     */
    inline double get_sPos() {
        const size_t n = getLocalNum();
        size_t numPrimBeamParts = 0;
        double z = 0.0;
        for(size_t i = 0; i < n; i++) {
            if(PType[i] == 0) {
                z += R[i](2);
                numPrimBeamParts++;
            }
        }

        if(numPrimBeamParts != 0)
            z /= numPrimBeamParts;

        reduce(z, z, OpAddAssign());
        z = z / Ippl::getNodes();

        return z;

        /*
             const size_t n = getTotalNum();
            if (n>0)
                return sum(R(2))/getTotalNum();
            else
                return 0.0;
        */
    }

    inline double   get_phase() const { return 1.0; }
    inline double   get_gamma() const { return 1.0; }

    inline double get_meanEnergy() const { return eKin_m; }
    //  inline double* get_energy() {return energy_m; }
    inline Vector_t get_origin() const { return rmin_m; }
    inline Vector_t get_maxExtend() const { return rmax_m; }
    inline Vector_t get_centroid() const { return rmean_m; }
    inline Vector_t get_rrms() const { return rrms_m; }
    inline Vector_t get_rmean() const { return rmean_m; }
    inline Vector_t get_prms() const { return prms_m; }
    inline Vector_t get_pmean() const { return pmean_m; }
    inline Vector_t get_emit() const { return eps_m; }
    inline Vector_t get_norm_emit() const { return eps_norm_m; }
    inline Vector_t get_hr() const { return hr_m;}

    inline double get_Dx() const { return Dx_m; }
    inline double get_Dy() const { return Dy_m; }

    inline double get_DDx() const { return DDx_m; }
    inline double get_DDy() const { return DDy_m; }




    inline void set_meshEnlargement(double dh) { dh_m = dh; }
    inline double get_meshEnlargement() const { return dh_m; }

    void gatherLoadBalanceStatistics();
    inline size_t getLoadBalance(int p) const { return globalPartPerNode_m[p]; }

    inline void get_PBounds(Vector_t &min, Vector_t &max) const {
        bounds(this->P, min, max);
    }

    void calcBeamParameters();
    void calcBeamParametersInitial(); // Calculate initial beam parameters before emission.
    void calcBeamParameters_cycl();

    inline double getCouplingConstant() const { return couplingConstant_m; }
    inline void setCouplingConstant(double c) { couplingConstant_m = c;}

    // set the charge per simulation particle
    inline void setCharge(double q) {
        qi_m = q;
        if(getTotalNum() != 0)
            Q = qi_m;
        else
            WARNMSG("Could not set total charge in PartBunch::setCharge based on getTotalNum" << endl);
    }

    // set the charge per simulation particle when total particle number equals 0
    inline void setChargeZeroPart(double q) {
        qi_m = q;
        
    }


    /// get the total charge per simulation particle
    inline double getCharge() const { return sum(Q); }

    /// get the macro particle charge
    inline double getChargePerParticle() const { return qi_m; }

    void setSolver(FieldSolver *fs);

    bool hasFieldSolver();

    // The structure for particle binning
    PartBins *pbin_m;

    inline void setLPath(double s) {lPath_m = s;}
    inline double getLPath() const {return lPath_m;}

    inline void setStepsPerTurn(int n) {stepsPerTurn_m = n;}
    inline int getStepsPerTurn() const {return stepsPerTurn_m;}

    inline void setTrackStep(long long n) {trackStep_m = n;}
    inline long long getTrackStep() const {return trackStep_m;}

    inline void setNumBunch(int n) {numBunch_m = n;}
    inline int getNumBunch() const {return numBunch_m;}

    inline void setSteptoLastInj(int n) {SteptoLastInj_m = n;}
    inline int getSteptoLastInj() {return SteptoLastInj_m;}

    /// calculate average angle of longitudinal direction of bins
    double calcMeanPhi();

    size_t getNumPartInBin(int BinID) const;

    /// reset Bin[] for each particle
    bool resetPartBinID();
    /// reset Bin[] for each particle according to the method given in paper PAST-AB(064402) by  G. Fubiani et al.
    bool resetPartBinID2(double eta);


    double getQ() const { return reference->getQ();}
    double getM() const { return reference->getM();}
    double getP() const { return reference->getP();}
    double getE() const { return reference->getE();}

    inline double getdE() { return dE_m;}
    double getBeta() const { return reference->getBeta();}
    double getGamma() const { return reference->getGamma();}
    virtual double getGamma(int i) {return 0;}
    virtual double getBeta(int i) {return 0;}
    virtual void actT() {};
    const PartData *getReference() const { return reference; }

    double getTBin();

private:

    const PartData *reference;
    void calcMoments();    // Calculates bunch moments using only emitted particles.
    void calcMomentsInitial(); // Calcualtes bunch moments by summing over bins (not accurate when any particles have been emitted).

    double calculateAngle(double x, double y);
    double calculateAngle2(double x, double y);

    /*
      Member variables starts here
    */

    /// hold the line-density
    double *lineDensity_m;
    /// how many bins the line-density has
    unsigned int nBinsLineDensity_m;

    /// holds the centroid of the beam
    double centroid_m[2*Dim];

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
    /// energy of the bunch
    double *energy_m;
    /// energy spread of the beam in keV
    double dE_m;

    /// maximal extend of particles
    Vector_t rmax_m;
    /// minimal extend of particles
    Vector_t rmin_m;

    /// rms beam size (m)
    Vector_t rrms_m;
    /// rms momenta
    Vector_t prms_m;
    /// mean beam size (m)
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

protected:
    /// timer for selfField calculation
    IpplTimings::TimerRef selfFieldTimer_m;
    IpplTimings::TimerRef compPotenTimer_m;
    IpplTimings::TimerRef boundpTimer_m;
    IpplTimings::TimerRef statParamTimer_m;
public:
    /// timer for IC, can not be in Distribution.h
    IpplTimings::TimerRef distrReload_m;
    IpplTimings::TimerRef distrCreate_m;


    /*
      For saving particles on a collimator or
      dead particles (absobed)
    */

    LossDataSink *lossDs_m;

private:
    /// Mesh enlargement
    double dh_m; /// in % how much the mesh is enlarged

    /// if larger than 0, emitt particles for tEmission_m [s]
    double tEmission_m;

    /// holds the gamma of the bin
    double *bingamma_m;

    //FIXME: this should go into the Bin class!
    // holds number of emitted particles of the bin
    // jjyang: opal-cycl use *nBin_m of pbin_m
    size_t *binemitted_m;

    /// path length from the start point
    double lPath_m;

    /// steps per turn for OPAL-cycl
    int stepsPerTurn_m;

    /// current track step
    long long trackStep_m;

    /// current bunch number
    int numBunch_m;

    /// this parameter records the current steps since last bunch injection
    /// it helps to inject new bunches correctly in the restart run of OPAL-cycl
    /// it is stored during phase space dump.
    int SteptoLastInj_m;

    /*
      Data structure for particle load balance information
    */

    double *partPerNode_m;
    double *globalPartPerNode_m;

    Distribution *dist_m;
};

inline Inform &operator<<(Inform &os, PartBunch &p) {
    return p.print(os);
}


#endif // OPAL_PartBunch_HH
