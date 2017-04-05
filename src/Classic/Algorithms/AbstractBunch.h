#ifndef ABSTRACT_BUNCH_H
#define ABSTRACT_BUNCH_H

class AbstractBunch {
    
public:
    
    virtual void runTests() = 0;

    virtual bool getIfBeamEmitting() = 0;
    
    virtual int getLastEmittedEnergyBin() = 0;
    
    virtual size_t getNumberOfEmissionSteps() = 0;
    
    virtual int getNumberOfEnergyBins() = 0;
    
    virtual void Rebin() = 0;
    
    virtual void setEnergyBins(int numberOfEnergyBins) = 0;
    
    virtual bool weHaveEnergyBins() = 0;
    
    virtual void switchToUnitlessPositions(bool use_dt_per_particle = false) = 0;
    
    virtual void switchOffUnitlessPositions(bool use_dt_per_particle = false) = 0;
    
    virtual void resetIfScan() = 0;
    
    virtual double getRho(int x, int y, int z) = 0;
    
    virtual void do_binaryRepart() = 0;
    
    virtual void calcLineDensity(unsigned int nBins,
                                 std::vector<double> &lineDensity,
                                 std::pair<double, double> &meshInfo) = 0;
    
    virtual void setDistribution(Distribution *d,
                         std::vector<Distribution *> addedDistributions,
                         size_t &np,
                         bool scan) = 0;
    
    virtual void setTEmission(double t) = 0;
    
    virtual double getTEmission() = 0;
    
    virtual bool doEmission() = 0;
    
    virtual bool weHaveBins() const = 0;
    
    virtual void setPBins(PartBins *pbin) = 0;
    
    virtual void setPBins(PartBinsCyc *pbin) = 0;
    
    virtual size_t emitParticles(double eZ) = 0;

    virtual void updateNumTotal() = 0;

    virtual void rebin() = 0;

    virtual int getNumBins() = 0;

    virtual int getLastemittedBin() = 0;
    
    virtual void calcGammas() = 0;

    virtual void calcGammas_cycl() = 0;
    
    virtual void setBinCharge(int bin, double q) = 0;
    
    virtual void setBinCharge(int bin) = 0;
    
    virtual std::pair<Vector_t, Vector_t> getEExtrema() = 0;
    
    virtual size_t calcNumPartsOutside(Vector_t x) = 0;
    
    virtual const Mesh_t &getMesh() const = 0;

    virtual Mesh_t &getMesh() = 0;

    virtual FieldLayout_t &getFieldLayout() = 0;

    virtual void setBCAllPeriodic() = 0;
    
    virtual void setBCAllOpen() = 0;

    virtual void setBCForDCBeam() = 0;

    virtual void boundp() = 0;

    virtual void boundp_destroy() = 0;

    virtual size_t boundp_destroyT() = 0;
    
    virtual size_t destroyT() = 0;
    
    virtual double getPy(int i) = 0;
    virtual double getPx(int i) = 0;
    virtual double getPz(int i) = 0;

    virtual double getPx0(int i) = 0;
    virtual double getPy0(int i) = 0;

    virtual double getX(int i) = 0;
    virtual double getY(int i) = 0;
    virtual double getZ(int i) = 0;

    virtual double getX0(int i) = 0;
    virtual double getY0(int i) = 0;

    virtual void setZ(int i, double zcoo) = 0;

    virtual void get_bounds(Vector_t &rmin, Vector_t &rmax) = 0;
    
    virtual void getLocalBounds(Vector_t &rmin, Vector_t &rmax) = 0;
    
    virtual std::pair<Vector_t, double> getBoundingSphere() = 0;
    
    virtual std::pair<Vector_t, double> getLocalBoundingSphere() = 0;
    
    
    virtual void push_back(Particle p) = 0;

    virtual void set_part(FVector<double, 6> z, int ii);

    virtual void set_part(Particle p, int ii);

    virtual Particle get_part(int ii);

    /// Return maximum amplitudes.
    //  The matrix [b]D[/b] is used to normalise the first two modes.
    //  The maximum normalised amplitudes for these modes are stored
    //  in [b]axmax[/b] and [b]aymax[/b].
    virtual void maximumAmplitudes(const FMatrix<double, 6, 6> &D,
                                   double &axmax, double &aymax) = 0;

    virtual Inform &print(Inform &os) = 0;


    virtual void   setdT(double dt) = 0;
    virtual double getdT() const = 0;

    virtual void   setT(double t) = 0;
    virtual double getT() const = 0;
    virtual void   computeSelfFields() = 0;

    /** /brief used for self fields with binned distribution */
    virtual void   computeSelfFields(int b) = 0;

    virtual void computeSelfFields_cycl(double gamma) = 0;
    virtual void computeSelfFields_cycl(int b) = 0;

    virtual void resetInterpolationCache(bool clearCache = false) = 0;

    /**
     * get the spos of the primary beam
     *
     * @param none
     *
     */
    virtual double get_sPos() = 0;
    virtual void set_sPos(double s) = 0;

    virtual double   get_gamma() const = 0;

    virtual double get_meanKineticEnergy() const = 0;
    virtual Vector_t get_origin() const = 0;
    virtual Vector_t get_maxExtent() const = 0;
    virtual Vector_t get_centroid() const = 0;
    virtual Vector_t get_rrms() const = 0;
    virtual Vector_t get_rmean() const = 0;
    virtual Vector_t get_prms() const = 0;
    virtual Vector_t get_rprms() const = 0;
    virtual Vector_t get_pmean() const = 0;
    virtual Vector_t get_pmean_Distribution() const = 0;
    virtual Vector_t get_emit() const = 0;
    virtual Vector_t get_norm_emit() const = 0;
    virtual Vector_t get_hr() const = 0;

    virtual double get_Dx() const = 0;
    virtual double get_Dy() const = 0;

    virtual double get_DDx() const = 0;
    virtual double get_DDy() const = 0;

    virtual void set_meshEnlargement(double dh) = 0;

    virtual void gatherLoadBalanceStatistics() = 0;
    virtual size_t getLoadBalance(int p) const = 0;


    virtual void get_PBounds(Vector_t &min, Vector_t &max) const = 0;

    virtual void calcBeamParameters() = 0;
    virtual void calcBeamParametersInitial() = 0; // Calculate initial beam parameters before emission.

    virtual double getCouplingConstant() const = 0;
    virtual void setCouplingConstant(double c) = 0;

    // set the charge per simulation particle
    virtual void setCharge(double q) = 0;
    // set the charge per simulation particle when total particle number equals 0
    virtual void setChargeZeroPart(double q) = 0;

    // set the mass per simulation particle
    virtual void setMass(double mass) = 0;


    /// \brief Need Ek for the Schottky effect calculation (eV)
    virtual double getEkin() const = 0;

    /// Need the work function for the Schottky effect calculation (eV)
    virtual double getWorkFunctionRf() const = 0;

    /// Need the laser energy for the Schottky effect calculation (eV)
    virtual double getLaserEnergy() const = 0;

    /// get the total charge per simulation particle
    virtual double getCharge() const = 0;

    /// get the macro particle charge
    virtual double getChargePerParticle() const = 0;

    virtual void setSolver(FieldSolver *fs) = 0;

    virtual bool hasFieldSolver() = 0;

    virtual std::string getFieldSolverType() const = 0;

    virtual void setLPath(double s) = 0;
    virtual double getLPath() const = 0;

    virtual void setStepsPerTurn(int n) = 0;
    virtual int getStepsPerTurn() const = 0;

    /// step in multiple TRACK commands
    virtual void setGlobalTrackStep(long long n) = 0;
    virtual long long getGlobalTrackStep() const = 0;

    /// step in a TRACK command
    virtual void setLocalTrackStep(long long n) = 0;
    virtual void incTrackSteps() = 0;
    virtual long long getLocalTrackStep() const = 0;

    virtual void setNumBunch(int n) = 0;
    virtual int getNumBunch() const = 0;

    virtual void setGlobalMeanR(Vector_t globalMeanR) = 0;
    virtual Vector_t getGlobalMeanR() = 0;
    virtual void setGlobalToLocalQuaternion(Quaternion_t globalToLocalQuaternion) = 0;
    virtual Quaternion_t getGlobalToLocalQuaternion() = 0;

    virtual void setSteptoLastInj(int n) = 0;
    virtual int getSteptoLastInj() = 0;

    /// calculate average angle of longitudinal direction of bins
    virtual double calcMeanPhi() = 0;

    /// reset Bin[] for each particle according to the method given in paper PAST-AB(064402) by  G. Fubiani et al.
    virtual bool resetPartBinID2(const double eta) = 0;


    virtual double getQ() const = 0;
    virtual double getM() const = 0;
    virtual double getP() const = 0;
    virtual double getE() const = 0;

    virtual void resetQ(double q) = 0;
    virtual void resetM(double m) = 0;

    virtual double getdE() = 0;
    virtual double getInitialBeta() const = 0;
    virtual double getInitialGamma() const = 0;
    virtual double getGamma(int i) = 0;
    virtual double getBeta(int i) = 0;
    virtual void actT() = 0;
    virtual const PartData *getReference() const = 0;

    virtual double getEmissionDeltaT() = 0;

    virtual Quaternion_t getQKs3D() = 0;
    virtual void         setQKs3D(Quaternion_t q) = 0;
    virtual Vector_t     getKs3DRefr() = 0;
    virtual void         setKs3DRefr(Vector_t r) = 0;
    virtual Vector_t     getKs3DRefp() = 0;
    virtual void         setKs3DRefp(Vector_t p) = 0;

    virtual void iterateEmittedBin(int binNumber) = 0;

    virtual void calcEMean() = 0;
    virtual void correctEnergy(double avrgp) = 0;

    virtual void swap(unsigned int i, unsigned int j) = 0;
    
    
    

    Vector_t RefPartR_m;
    Vector_t RefPartP_m;
    CoordinateSystemTrafo toLabTrafo_m;
    
    PartBins *pbin_m;

    std::unique_ptr<LossDataSink> lossDs_m;

    /// timer for IC, can not be in Distribution.h
    IpplTimings::TimerRef distrReload_m;
    IpplTimings::TimerRef distrCreate_m;
    
    // For AMTS integrator in OPAL-T
    double dtScInit_m, deltaTau_m;
    
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
    
protected:
    virtual void cleanUpParticles() = 0;
    
    virtual bool isGridFixed() = 0;
    
    virtual double getBinGamma(int bin) = 0;
};

#endif
