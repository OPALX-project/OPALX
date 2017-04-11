#ifndef OPAL_Distribution_HH
#define OPAL_Distribution_HH

// ------------------------------------------------------------------------
// $RCSfile: Distribution.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Distribution
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:44 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------
#include <iosfwd>
#include <fstream>
#include <forward_list>
#include <string>
#include <map>

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"

#include "Algorithms/Vektor.h"
#include "Beamlines/Beamline.h"

#include "Ippl.h"

#include "H5hut.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_qrng.h>

#ifdef WITH_UNIT_TESTS
#include <gtest/gtest_prod.h>
#endif

class Beam;

template <class T, unsigned Dim>
class PartBunchBase;

class PartBins;
class EnvelopeBunch;
class BoundaryGeometry;
class LaserProfile;
class H5PartWrapper;

namespace DistrTypeT
{
    enum DistrTypeT {NODIST,
                    FROMFILE,
                    GAUSS,
                    BINOMIAL,
                    FLATTOP,
                    SURFACEEMISSION,
                    SURFACERANDCREATE,
                    GUNGAUSSFLATTOPTH,
	            ASTRAFLATTOPTH,
		    MATCHEDGAUSS
                    };
}

namespace EmissionModelT
{
    enum EmissionModelT {NONE,
                         ASTRA,
                         NONEQUIL
                        };
}

namespace InputMomentumUnitsT
{
    enum InputMomentumUnitsT {NONE,
                              EV
                              };
}

/*
 * Class Distribution
 *
 * Defines the initial beam that is injected or emitted into the simulation.
 */

class Distribution: public Definition {

public:

    Distribution();
    virtual ~Distribution();

    virtual bool canReplaceBy(Object *object);
    virtual Distribution *clone(const std::string &name);
    virtual void execute();
    virtual void update();

    void createBoundaryGeometry(PartBunchBase<double, 3> *p, BoundaryGeometry &bg);
    void createOpalCycl(PartBunchBase<double, 3> *beam,
                        size_t numberOfParticles,
			double current, const Beamline &bl,
                        bool scan);
    void createOpalE(Beam *beam,
                     std::vector<Distribution *> addedDistributions,
                     EnvelopeBunch *envelopeBunch,
                     double distCenter,
                     double Bz0);
    void createOpalT(PartBunchBase<double, 3> *beam,
                     std::vector<Distribution *> addedDistributions,
                     size_t &numberOfParticles,
                     bool scan);
    void createOpalT(PartBunchBase<double, 3> *beam, size_t &numberOfParticles, bool scan);
    void createPriPart(PartBunchBase<double, 3> *beam, BoundaryGeometry &bg);
    void doRestartOpalT(PartBunchBase<double, 3> *p, size_t Np, int restartStep, H5PartWrapper *h5wrapper);
    void doRestartOpalCycl(PartBunchBase<double, 3> *p, size_t Np, int restartStep,
                        const int specifiedNumBunch, H5PartWrapper *h5wrapper);
    void doRestartOpalE(EnvelopeBunch *p, size_t Np, int restartStep, H5PartWrapper *h5wrapper);
    size_t emitParticles(PartBunchBase<double, 3> *beam, double eZ);
    double getPercentageEmitted() const;
    static Distribution *find(const std::string &name);

    void eraseXDist();
    void eraseBGxDist();
    void eraseYDist();
    void eraseBGyDist();
    void eraseTOrZDist();
    void eraseBGzDist();
    bool getIfDistEmitting();
    int getLastEmittedEnergyBin();
    double getMaxTOrZ();
    double getMinTOrZ();
    size_t getNumberOfEmissionSteps();
    int getNumberOfEnergyBins();
    double getEmissionDeltaT();
    double getEnergyBinDeltaT();
    double getWeight();
    std::vector<double>& getXDist();
    std::vector<double>& getBGxDist();
    std::vector<double>& getYDist();
    std::vector<double>& getBGyDist();
    std::vector<double>& getTOrZDist();
    std::vector<double>& getBGzDist();

    /// Return the embedded CLASSIC PartData.
    const PartData &getReference() const;
    double getTEmission();

    Vector_t get_pmean() const;
    double getEkin() const;
    double getLaserEnergy() const;
    double getWorkFunctionRf() const;

    size_t getNumberOfDarkCurrentParticles();
    double getDarkCurrentParticlesInwardMargin();
    double getEInitThreshold();
    double getWorkFunction();
    double getFieldEnhancement();
    size_t getMaxFNemissionPartPerTri();
    double getFieldFNThreshold();
    double getFNParameterA();
    double getFNParameterB();
    double getFNParameterY();
    double getFNParameterVYZero();
    double getFNParameterVYSecond();
    int    getSecondaryEmissionFlag();
    bool   getEmissionMode() ;

    std::string getTypeofDistribution();

    double getvSeyZero();//return sey_0 in Vaughan's model
    double getvEZero();//return the energy related to sey_0 in Vaughan's model
    double getvSeyMax();//return sey max in Vaughan's model
    double getvEmax();//return Emax in Vaughan's model
    double getvKenergy();//return fitting parameter denotes the roughness of surface for impact energy in Vaughan's model
    double getvKtheta();//return fitting parameter denotes the roughness of surface for impact angle in Vaughan's model
    double getvVThermal();//return thermal velocity of Maxwellian distribution of secondaries in Vaughan's model
    double getVw();//return velocity scalar for parallel plate benchmark;
    int getSurfMaterial();//material type for Furman-Pivi's model 0 for copper, 1 for stainless steel

    Inform &printInfo(Inform &os) const;

    bool Rebin();
    void setDistToEmitted(bool emitted);
    void setDistType();
    void shiftBeam(double &maxTOrZ, double &minTOrZ);
    double getEmissionTimeShift() const;

    // in case if OPAL-cycl in restart mode
    double GetBeGa() {return bega_m;}
    double GetPr() {return referencePr_m;}
    double GetPt() {return referencePt_m;}
    double GetPz() {return referencePz_m;}
    double GetR() {return referenceR_m;}
    double GetTheta() {return referenceTheta_m;}
    double GetZ() {return referenceZ_m;}

    double GetPhi() {return phi_m;}
    double GetPsi() {return psi_m;}
    bool GetPreviousH5Local() {return previousH5Local_m;}

    void setNumberOfDistributions(unsigned int n) { numberOfDistributions_m = n; }

    DistrTypeT::DistrTypeT getType() const;
private:
#ifdef WITH_UNIT_TESTS
    FRIEND_TEST(GaussTest, FullSigmaTest1);
    FRIEND_TEST(GaussTest, FullSigmaTest2);
#endif

    Distribution(const std::string &name, Distribution *parent);

    // Not implemented.
    Distribution(const Distribution &);
    void operator=(const Distribution &);


    //    void printSigma(SigmaGenerator<double,unsigned int>::matrix_type& M, Inform& out);

    void addDistributions();
    void applyEmissionModel(double lowEnergyLimit, double &px, double &py, double &pz, std::vector<double> &additionalRNs);
    void applyEmissModelAstra(double &px, double &py, double &pz, std::vector<double> &additionalRNs);
    void applyEmissModelNone(double &pz);
    void applyEmissModelNonEquil(double eZ, double &px, double &py, double &pz, std::vector<double> &additionalRNs);
    void create(size_t &numberOfParticles, double massIneV);
    void calcPartPerDist(size_t numberOfParticles);
    void checkEmissionParameters();
    void checkIfEmitted();
    void checkParticleNumber(size_t &numberOfParticles);
    void chooseInputMomentumUnits(InputMomentumUnitsT::InputMomentumUnitsT inputMoUnits);
    double convertBetaGammaToeV(double valueInbega, double mass);
    double converteVToBetaGamma(double valueIneV, double massIneV);
    double convertMeVPerCToBetaGamma(double valueInMeVPerC, double massIneV);
    void createDistributionBinomial(size_t numberOfParticles, double massIneV);
    void createDistributionFlattop(size_t numberOfParticles, double massIneV);
    void createDistributionFromFile(size_t numberOfParticles, double massIneV);
    void createDistributionGauss(size_t numberOfParticles, double massIneV);
    void createMatchedGaussDistribution(size_t numberOfParticles, double massIneV);
    void destroyBeam(PartBunchBase<double, 3> *beam);
    void fillEBinHistogram();
    void fillParticleBins();
    size_t findEBin(double tOrZ);
    void generateAstraFlattopT(size_t numberOfParticles);
    void generateBinomial(size_t numberOfParticles);
    void generateFlattopLaserProfile(size_t numberOfParticles);
    void generateFlattopT(size_t numberOfParticles);
    void generateFlattopZ(size_t numberOfParticles);
    void generateGaussZ(size_t numberOfParticles);
    void generateLongFlattopT(size_t numberOfParticles);
    void generateTransverseGauss(size_t numberOfParticles);
    void initializeBeam(PartBunchBase<double, 3> *beam);
    void injectBeam(PartBunchBase<double, 3> *beam);
    void printDist(Inform &os, size_t numberOfParticles) const;
    void printDistBinomial(Inform &os) const;
    void printDistFlattop(Inform &os) const;
    void printDistFromFile(Inform &os) const;
    void printDistGauss(Inform &os) const;
    void printDistMatchedGauss(Inform &os) const;
    void printDistSurfEmission(Inform &os) const;
    void printDistSurfAndCreate(Inform &os) const;
    void printEmissionModel(Inform &os) const;
    void printEmissionModelAstra(Inform &os) const;
    void printEmissionModelNone(Inform &os) const;
    void printEmissionModelNonEquil(Inform &os) const;
    void printEnergyBins(Inform &os) const;
    void reflectDistribution(size_t &numberOfParticles);
    void scaleDistCoordinates();
    void setAttributes();
    void setDistParametersBinomial(double massIneV);
    void setDistParametersFlattop(double massIneV);
    void setDistParametersGauss(double massIneV);
    void setEmissionTime(double &maxT, double &minT);
    void setFieldEmissionParameters();
    void setupEmissionModel(PartBunchBase<double, 3> *beam);
    void setupEmissionModelAstra(PartBunchBase<double, 3> *beam);
    void setupEmissionModelNone(PartBunchBase<double, 3> *beam);
    void setupEmissionModelNonEquil();
    void setupEnergyBins(double maxTOrZ, double minTOrZ);
    void setupParticleBins(double massIneV, PartBunchBase<double, 3> *beam);
    void shiftDistCoordinates(double massIneV);
    void writeOutFileHeader();
    void writeOutFileEmission();
    void writeOutFileInjection();
    void writeToFile();

    std::string distT_m;                 /// Distribution type. Declared as string
    DistrTypeT::DistrTypeT distrTypeT_m; /// and list type for switch statements.
    std::ofstream os_m;                  /// Output file to write distribution.
    void writeToFileCycl(PartBunchBase<double, 3> *beam, size_t Np);

    unsigned int numberOfDistributions_m;

    bool emitting_m;                     /// Distribution is an emitted, and is currently
                                         /// emitting, rather than an injected, beam.

    bool scan_m;                         /// If we are doing a scan, we need to
                                         /// destroy existing particles before
                                         /// each run.

    PartData particleRefData_m;          /// Reference data for particle type (charge,
                                         /// mass etc.)

    /// Vector of distributions to be added to this distribution.
    std::vector<Distribution *> addedDistributions_m;
    std::vector<size_t> particlesPerDist_m;

    /// Emission Model.
    EmissionModelT::EmissionModelT emissionModel_m;

    /// Emission parameters.
    double tEmission_m;
    double tBin_m;
    double currentEmissionTime_m;
    int currentEnergyBin_m;
    int currentSampleBin_m;
    int numberOfEnergyBins_m;       /// Number of energy bins the distribution
                                    /// is broken into. Used for an emitted beam.
    int numberOfSampleBins_m;       /// Number of samples to use per energy bin
                                    /// when emitting beam.
    PartBins *energyBins_m;         /// Distribution energy bins.
    gsl_histogram *energyBinHist_m; /// GSL histogram used to define energy bin
                                    /// structure.

    gsl_rng *randGen_m;             /// Random number generator

    // ASTRA and NONE photo emission model.
    double pTotThermal_m;           /// Total thermal momentum.
    Vector_t pmean_m;

    // NONEQUIL photo emission model.
    double cathodeWorkFunc_m;       /// Cathode material work function (eV).
    double laserEnergy_m;           /// Laser photon energy (eV).
    double cathodeFermiEnergy_m;    /// Cathode material Fermi energy (eV).
    double cathodeTemp_m;           /// Cathode temperature (K).
    double emitEnergyUpperLimit_m;  /// Upper limit on emission energy distribution (eV).

    std::vector<std::vector<double> > additionalRNs_m;

    size_t totalNumberParticles_m;
    size_t totalNumberEmittedParticles_m;

    // Beam coordinate containers.
    std::vector<double> xDist_m;
    std::vector<double> pxDist_m;
    std::vector<double> yDist_m;
    std::vector<double> pyDist_m;
    std::vector<double> tOrZDist_m;
    std::vector<double> pzDist_m;

    // Initial coordinates for file write.
    std::vector<double> xWrite_m;
    std::vector<double> pxWrite_m;
    std::vector<double> yWrite_m;
    std::vector<double> pyWrite_m;
    std::vector<double> tOrZWrite_m;
    std::vector<double> pzWrite_m;
    std::vector<size_t> binWrite_m;

    // for compatibility reasons
    double avrgpz_m;



    //Distribution parameters.
    InputMomentumUnitsT::InputMomentumUnitsT inputMoUnits_m;
    double sigmaTRise_m;
    double sigmaTFall_m;
    double tPulseLengthFWHM_m;
    Vector_t sigmaR_m;
    Vector_t sigmaP_m;
    Vector_t cutoffR_m;
    Vector_t cutoffP_m;
    Vector_t mBinomial_m;
    SymTenzor<double, 6> correlationMatrix_m;

    // Laser profile.
    std::string laserProfileFileName_m;
    std::string laserImageName_m;
    double laserIntensityCut_m;
    LaserProfile *laserProfile_m;

    /*
     * Dark current calculation parameters.
     */
    size_t darkCurrentParts_m;      /// Number of dark current particles.
    double darkInwardMargin_m;      /// Dark current particle initialization position.
                                    /// Inward along the triangle normal, positive.
                                    /// Inside the geometry.
    double eInitThreshold_m;        /// Field threshold (MV/m) beyond which particles
                                    /// will be initialized.
    double workFunction_m;          /// Work function of surface material (eV).
    double fieldEnhancement_m;      /// Field enhancement factor beta for Fowler-
                                    /// Nordheim emission.
    double fieldThrFN_m;            /// Field threshold for Fowler-Nordheim
                                    /// emission (MV/m).
    size_t maxFN_m;                 /// Max. number of electrons emitted from a
                                    /// single triangle for Fowler-Nordheim emission.
    double paraFNA_m;               /// Empirical constant A in Fowler-Nordheim
                                    /// emission model.
    double paraFNB_m;               /// Empirical constant B in Fowler-Nordheim
                                    /// emission model.
    double paraFNY_m;               /// Constant for image charge effect parameter y(E)
                                    /// in Fowler-Nordheim emission model.
    double paraFNVYSe_m;            /// Second order constant for v(y) function in
                                    /// Fowler-Nordheim emission model.
    double paraFNVYZe_m;            /// Zero order constant for v(y) function in
                                    /// Fowler-Nordheim emission model.
    int    secondaryFlag_m;         /// Select the secondary model type:
                                    ///     0           ==> no secondary emission.
                                    ///     1           ==> Furman-Pivi
                                    ///     2 or larger ==> Vaughan's model
    double ppVw_m;                  /// Velocity scalar for parallel plate benchmark.
    double vVThermal_m;             /// Thermal velocity of Maxwellian distribution
                                    /// of secondaries in Vaughan's model.


    // AAA This is for the matched distribution
    double ex_m;
    double ey_m;
    double et_m;

    double I_m;
    double E_m;
    double bg_m;                      /// beta gamma
    double M_m;                       /// mass in terms of proton mass
    std::string bfieldfn_m;           /// only temporarly





    // Some legacy members that need to be cleaned up.

    /// time binned distribution with thermal energy
    double tRise_m;
    double tFall_m;
    double sigmaRise_m;
    double sigmaFall_m;
    double cutoff_m;

    // Cyclotron stuff
    double referencePr_m;
    double referencePt_m;
    double referencePz_m;
    double referenceR_m;
    double referenceTheta_m;
    double referenceZ_m;
    double bega_m;

    // Cyclotron for restart in local mode
    double phi_m;
    double psi_m;
    bool previousH5Local_m;
};

inline Inform &operator<<(Inform &os, const Distribution &d) {
    return d.printInfo(os);
}

inline
Vector_t Distribution::get_pmean() const {
    return pmean_m;
}

inline
DistrTypeT::DistrTypeT Distribution::getType() const {
    return distrTypeT_m;
}

inline double Distribution::getPercentageEmitted() const {
    return (double)totalNumberEmittedParticles_m / (double)totalNumberParticles_m;
}

#endif // OPAL_Distribution_HH