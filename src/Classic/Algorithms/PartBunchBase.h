#ifndef PART_BUNCH_BASE_H
#define PART_BUNCH_BASE_H

#include "Algorithms/AbstractBunch.h"

template<class Bunch>
class PartBunchBase : public AbstractBunch {
    
public:
    
    PartBunchBase(const PartData *ref);
    
    PartBunchBase(const std::vector<OpalParticle> &rhs,
                  const PartData *ref);
    
    PartBunchBase(const PartBunchBase &rhs):
    
    virtual ~PartBunchBase();
    
    /*
     * bunch-specific functions
     */
    
    void runTests() {
        static_cast<Bunch*>(this)->runTests();
    }
    
    double getRho(int x, int y, int z) {
        return static_cast<Bunch*>(this)->getRho(x, y, z);
    }
    
    void setBCAllPeriodic() {
        static_cast<Bunch*>(this)->setBCAllPeriodic();
    }
    
    void boundp() {
        static_cast<Bunch*>(this)->boundp();
    }

    void boundp_destroy() {
        static_cast<Bunch*>(this)->boundp_destroy();
    }

    size_t boundp_destroyT() {
        return static_cast<Bunch*>(this)->boundp_destroyT();
    }
    
    void setBCAllOpen() {
        static_cast<Bunch*>(this)->setBCAllOpen();
    }

    void setBCForDCBeam() {
        static_cast<Bunch*>(this)->setBCForDCBeam();
    }
    
    /*
     * bunch-common functions
     */
    
    bool getIfBeamEmitting();
    
    int getLastEmittedEnergyBin();
    
    size_t getNumberOfEmissionSteps();
    
    int getNumberOfEnergyBins();
    
    void Rebin();
    
    void setEnergyBins(int numberOfEnergyBins);
    
    bool weHaveEnergyBins();
    
    void switchToUnitlessPositions(bool use_dt_per_particle = false);
    
    void switchOffUnitlessPositions(bool use_dt_per_particle = false);
    
    void resetIfScan();
    
    void do_binaryRepart();
    
    void calcLineDensity(unsigned int nBins,
                                 std::vector<double> &lineDensity,
                                 std::pair<double, double> &meshInfo) {
        static_cast<Bunch*>(this)->calcLineDensity(nBins, lineDensity, meshInfo);
    }
    
    void setDistribution(Distribution *d,
                         std::vector<Distribution *> addedDistributions,
                         size_t &np,
                         bool scan);
    
    void setTEmission(double t);
    
    double getTEmission();
    
    bool doEmission();
    
    bool weHaveBins() const;
    
    void setPBins(PartBins *pbin);
    
    void setPBins(PartBinsCyc *pbin);
    
    size_t emitParticles(double eZ);

    void updateNumTotal();

    void rebin();

    int getNumBins();

    int getLastemittedBin();
    
    void calcGammas();

    void calcGammas_cycl();
    
    void setBinCharge(int bin, double q);
    
    void setBinCharge(int bin);
    
    std::pair<Vector_t, Vector_t> getEExtrema();
    
    size_t calcNumPartsOutside(Vector_t x);
    
    size_t destroyT();
//     
//     const Mesh_t &getMesh() const {
//         static_cast<Bunch*>(this)->
//     }
// 
//     Mesh_t &getMesh() {
//         static_cast<Bunch*>(this)->
//     }
// 
//     FieldLayout_t &getFieldLayout() {
//         static_cast<Bunch*>(this)->
//     }
// 
//     
// 
//     
    virtual double getPy(int i);
    
    virtual double getPx(int i);
    
    virtual double getPz(int i);

    virtual double getPx0(int i);
    
    virtual double getPy0(int i):

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

    Inform &print(Inform &os) {
        return static_cast<Bunch*>(this)->print;
    }

    void setdT(double dt);
    
    double getdT() const;

    void   setT(double t);
    double getT() const;
    
    void   computeSelfFields() {
        static_cast<Bunch*>(this)->computeSelfFields();
    }

    /** /brief used for self fields with binned distribution */
    void   computeSelfFields(int b) {
        static_cast<Bunch*>(this)->computeSelfFields(b);
    }

    void computeSelfFields_cycl(double gamma) {
        static_cast<Bunch*>(this)->computeSelfFields_cycl(gamma);
    }
    
    void computeSelfFields_cycl(int b) {
        static_cast<Bunch*>(this)->computeSelfFields_cycl(b);
    }
    
// 
//     void resetInterpolationCache(bool clearCache = false);
// 
//     /**
//      * get the spos of the primary beam
//      *
//      * @param none
//      *
//      */
//     double get_sPos();
//     void set_sPos(double s);
// 
//     double   get_gamma() const;
// 
//     double get_meanKineticEnergy() const;
//     Vector_t get_origin() const;
//     Vector_t get_maxExtent() const;
//     Vector_t get_centroid() const;
//     Vector_t get_rrms() const;
//     Vector_t get_rmean() const;
//     Vector_t get_prms() const;
//     Vector_t get_rprms() const;
//     Vector_t get_pmean() const;
//     Vector_t get_pmean_Distribution() const;
//     Vector_t get_emit() const;
//     Vector_t get_norm_emit() const;
//     Vector_t get_hr() const;
// 
//     double get_Dx() const;
//     double get_Dy() const;
// 
//     double get_DDx() const;
//     double get_DDy() const;
// 
//     void set_meshEnlargement(double dh);
// 
//     void gatherLoadBalanceStatistics();
//     size_t getLoadBalance(int p) const;
// 
// 
//     void get_PBounds(Vector_t &min, Vector_t &max) const;
// 
//     void calcBeamParameters();
//     void calcBeamParametersInitial(); // Calculate initial beam parameters before emission.
// 
//     double getCouplingConstant() const;
//     void setCouplingConstant(double c);
// 
//     // set the charge per simulation particle
//     void setCharge(double q);
//     // set the charge per simulation particle when total particle number equals 0
//     void setChargeZeroPart(double q);
// 
//     // set the mass per simulation particle
//     void setMass(double mass);
// 
// 
//     /// \brief Need Ek for the Schottky effect calculation (eV)
//     double getEkin() const;
// 
//     /// Need the work function for the Schottky effect calculation (eV)
//     double getWorkFunctionRf() const;
// 
//     /// Need the laser energy for the Schottky effect calculation (eV)
//     double getLaserEnergy() const;
// 
//     /// get the total charge per simulation particle
//     double getCharge() const;
// 
//     /// get the macro particle charge
//     double getChargePerParticle() const;
// 
//     void setSolver(FieldSolver *fs);
// 
//     bool hasFieldSolver();
// 
//     std::string getFieldSolverType() const;
// 
//     void setLPath(double s);
//     double getLPath() const;
// 
//     void setStepsPerTurn(int n);
//     int getStepsPerTurn() const;
// 
//     /// step in multiple TRACK commands
//     void setGlobalTrackStep(long long n);
//     long long getGlobalTrackStep() const;
// 
//     /// step in a TRACK command
//     void setLocalTrackStep(long long n);
//     void incTrackSteps();
//     long long getLocalTrackStep() const;
// 
//     void setNumBunch(int n);
//     int getNumBunch() const;
// 
//     void setGlobalMeanR(Vector_t globalMeanR);
//     Vector_t getGlobalMeanR();
//     void setGlobalToLocalQuaternion(Quaternion_t globalToLocalQuaternion);
//     Quaternion_t getGlobalToLocalQuaternion();
// 
//     void setSteptoLastInj(int n);
//     int getSteptoLastInj();
// 
//     /// calculate average angle of longitudinal direction of bins
//     double calcMeanPhi();
// 
//     /// reset Bin[] for each particle according to the method given in paper PAST-AB(064402) by  G. Fubiani et al.
//     bool resetPartBinID2(const double eta);
// 
// 
//     double getQ() const;
//     double getM() const;
//     double getP() const;
//     double getE() const;
// 
//     void resetQ(double q);
//     void resetM(double m);
// 
//     double getdE();
//     double getInitialBeta() const;
//     double getInitialGamma() const;
//     double getGamma(int i);
//     double getBeta(int i);
//     void actT();
//     const PartData *getReference() const;
// 
//     double getEmissionDeltaT();
// 
//     Quaternion_t getQKs3D();
//     void         setQKs3D(Quaternion_t q);
//     Vector_t     getKs3DRefr();
//     void         setKs3DRefr(Vector_t r);
//     Vector_t     getKs3DRefp();
//     void         setKs3DRefp(Vector_t p);
// 
//     void iterateEmittedBin(int binNumber);
// 
//     void calcEMean();
//     void correctEnergy(double avrgp);
// 
//     void swap(unsigned int i, unsigned int j);
    
//     PartBunchBase() : bunch_mp(static_cast<Bunch*>(this)) { }

protected:
    long long localTrackStep_m;
    
private:
    Bunch* bunch_mp;
};

#include "PartBunchBase.hpp"

#endif
