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

#include "Ippl.h"

#include "H5hut.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_qrng.h>

class Beam;
class PartBunch;
class PartBins;
class EnvelopeBunch;
class BoundaryGeometry;
class LaserProfile;

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
                    ASTRAFLATTOPTH
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

    void CreateBoundaryGeometry(PartBunch &p, BoundaryGeometry &bg);
    void CreateOpalCycl(PartBunch &beam,
                        size_t numberOfParticles,
                        bool scan);
    void CreateOpalE(Beam *beam,
                     std::vector<Distribution *> addedDistributions,
                     EnvelopeBunch *envelopeBunch,
                     double distCenter,
                     double Bz0);
    void CreateOpalT(PartBunch &beam,
                     std::vector<Distribution *> addedDistributions,
                     size_t &numberOfParticles,
                     bool scan);
    void CreateOpalT(PartBunch &beam, size_t &numberOfParticles, bool scan);
    void CreatePriPart(PartBunch *beam, BoundaryGeometry &bg);
    void DoRestartOpalT(PartBunch &p, size_t Np, int restartStep);
    void DoRestartOpalCycl(PartBunch &p, size_t Np, int restartStep,
                        const int specifiedNumBunch);
    void DoRestartOpalE(EnvelopeBunch &p, size_t Np, int restartStep);
    size_t EmitParticles(PartBunch &beam, double eZ);
    static Distribution *find(const std::string &name);

    void EraseXDist();
    void EraseBGxDist();
    void EraseYDist();
    void EraseBGyDist();
    void EraseTOrZDist();
    void EraseBGzDist();
    bool GetIfDistEmitting();
    int GetLastEmittedEnergyBin();
    double GetMaxTOrZ();
    double GetMinTOrZ();
    size_t GetNumberOfEmissionSteps();
    int GetNumberOfEnergyBins();
    double GetEmissionDeltaT();
    double GetEnergyBinDeltaT();
    double GetWeight();
    std::vector<double>& GetXDist();
    std::vector<double>& GetBGxDist();
    std::vector<double>& GetYDist();
    std::vector<double>& GetBGyDist();
    std::vector<double>& GetTOrZDist();
    std::vector<double>& GetBGzDist();

    /// Return the embedded CLASSIC PartData.
    const PartData &GetReference() const;
    double GetTEmission();

    double GetEkin() const;
    double GetLaserEnergy() const;
    double GetWorkFunctionRf() const;

    size_t GetNumberOfDarkCurrentParticles();
    double GetDarkCurrentParticlesInwardMargin();
    double GetEInitThreshold();
    double GetWorkFunction();
    double GetFieldEnhancement();
    size_t GetMaxFNemissionPartPerTri();
    double GetFieldFNThreshold();
    double GetFNParameterA();
    double GetFNParameterB();
    double GetFNParameterY();
    double GetFNParameterVYZero();
    double GetFNParameterVYSecond();
    int    GetSecondaryEmissionFlag();
    bool   GetEmissionMode() ;

    std::string GetTypeofDistribution();

    double GetvSeyZero();//return sey_0 in Vaughan's model
    double GetvEZero();//return the energy related to sey_0 in Vaughan's model
    double GetvSeyMax();//return sey max in Vaughan's model
    double GetvEmax();//return Emax in Vaughan's model
    double GetvKenergy();//return fitting parameter denotes the roughness of surface for impact energy in Vaughan's model
    double GetvKtheta();//return fitting parameter denotes the roughness of surface for impact angle in Vaughan's model
    double GetvVThermal();//return thermal velocity of Maxwellian distribution of secondaries in Vaughan's model
    double GetVw();//return velocity scalar for parallel plate benchmark;
    int GetSurfMaterial();//material type for Furman-Pivi's model 0 for copper, 1 for stainless steel

    Inform &printInfo(Inform &os) const;

    bool Rebin();
    void SetDistToEmitted(bool emitted);
    void SetDistType();
    void ShiftBeam(double &maxTOrZ, double &minTOrZ);

    // in case if OPAL-cycl in restart mode
    double GetBeGa() {return bega_m;}
    double GetPr() {return referencePr_m;}
    double GetPt() {return referencePt_m;}
    double GetPz() {return referencePz_m;}
    double GetR() {return referenceR_m;}
    double GetTheta() {return referenceTheta_m;}
    double GetZ() {return referenceZ_m;}

private:

    Distribution(const std::string &name, Distribution *parent);

    // Not implemented.
    Distribution(const Distribution &);
    void operator=(const Distribution &);

    void AddDistributions();
    void ApplyEmissionModel(double eZ, double &px, double &py, double &pz);
    void ApplyEmissModelAstra(double &px, double &py, double &pz);
    void ApplyEmissModelNone(double &pz);
    void ApplyEmissModelNonEquil(double eZ, double &px, double &py, double &pz);
    void Create(size_t &numberOfParticles, double massIneV);
    void CalcPartPerDist(size_t numberOfParticles);
    void CheckEmissionParameters();
    void CheckIfEmitted();
    void CheckParticleNumber(size_t &numberOfParticles);
    void ChooseInputMomentumUnits(InputMomentumUnitsT::InputMomentumUnitsT inputMoUnits);
    double ConvertBetaGammaToeV(double valueInbega, double mass);
    double ConverteVToBetaGamma(double valueIneV, double massIneV);
    double ConvertMeVPerCToBetaGamma(double valueInMeVPerC, double massIneV);
    void CreateDistributionBinomial(size_t numberOfParticles, double massIneV);
    void CreateDistributionFlattop(size_t numberOfParticles, double massIneV);
    void CreateDistributionFromFile(size_t numberOfParticles, double massIneV);
    void CreateDistributionGauss(size_t numberOfParticles, double massIneV);
    void DestroyBeam(PartBunch &beam);
    void FillEBinHistogram();
    void FillParticleBins();
    size_t FindEBin(double tOrZ);
    void GenerateAstraFlattopT(size_t numberOfParticles);
    void GenerateBinomial(size_t numberOfParticles);
    void GenerateFlattopLaserProfile(size_t numberOfParticles);
    void GenerateFlattopT(size_t numberOfParticles);
    void GenerateFlattopZ(size_t numberOfParticles);
    void GenerateGaussZ(size_t numberOfParticles);
    void GenerateLongFlattopT(size_t numberOfParticles);
    void GenerateTransverseGauss(size_t numberOfParticles);
    void InitializeBeam(PartBunch &beam);
    void InjectBeam(PartBunch &beam);
    void PrintDist(Inform &os, size_t numberOfParticles) const;
    void PrintDistBinomial(Inform &os) const;
    void PrintDistFlattop(Inform &os) const;
    void PrintDistFromFile(Inform &os) const;
    void PrintDistGauss(Inform &os) const;
    void PrintDistSurfEmission(Inform &os) const;
    void PrintDistSurfAndCreate(Inform &os) const;
    void PrintEmissionModel(Inform &os) const;
    void PrintEmissionModelAstra(Inform &os) const;
    void PrintEmissionModelNone(Inform &os) const;
    void PrintEmissionModelNonEquil(Inform &os) const;
    void PrintEnergyBins(Inform &os) const;
    void ReflectDistribution(size_t &numberOfParticles);
    void ScaleDistCoordinates();
    void SetAttributes();
    void SetDistParametersBinomial(double massIneV);
    void SetDistParametersFlattop(double massIneV);
    void SetDistParametersGauss(double massIneV);
    void SetEmissionTime(double &maxT, double &minT);
    void SetFieldEmissionParameters();
    void SetupEmissionModel(PartBunch &beam);
    void SetupEmissionModelAstra(PartBunch &beam);
    void SetupEmissionModelNone(PartBunch &beam);
    void SetupEmissionModelNonEquil();
    void SetupEnergyBins(double maxTOrZ, double minTOrZ);
    void SetupParticleBins(double massIneV, PartBunch &beam);
    void ShiftDistCoordinates(double massIneV);
    void WriteOutFileHeader();
    void WriteOutFileEmission();
    void WriteOutFileInjection();
    void WriteToFile();

    std::string distT_m;                 /// Distribution type. Declared as string
    DistrTypeT::DistrTypeT distrTypeT_m; /// and list type for switch statements.
    std::ofstream os_m;                  /// Output file to write distribution.
    void writeToFileCycl(PartBunch &beam, size_t Np);

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

    gsl_rng *randGenEmit_m;         /// Random number generator for thermal emission
                                    /// models.

    // ASTRA and NONE photo emission model.
    double pTotThermal_m;           /// Total thermal momentum.

    // NONEQUIL photo emission model.
    double cathodeWorkFunc_m;       /// Cathode material work function (eV).
    double laserEnergy_m;           /// Laser photon energy (eV).
    double cathodeFermiEnergy_m;    /// Cathode material Fermi energy (eV).
    double cathodeTemp_m;           /// Cathode temperature (K).
    double emitEnergyUpperLimit_m;  /// Upper limit on emission energy distribution (eV).


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
    std::vector<double> distCorr_m;

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


    // Some legacy members that need to be cleaned up.

    /// time binned distribution with thermal energy
    double tRise_m;
    double tFall_m;
    double sigmaRise_m;
    double sigmaFall_m;
    double cutoff_m;

    double referencePr_m;
    double referencePt_m;
    double referencePz_m = 0.0;
    double referenceR_m;
    double referenceTheta_m;
    double referenceZ_m = 0.0;
    double bega_m;

};

inline Inform &operator<<(Inform &os, const Distribution &d) {
    return d.printInfo(os);
}

#endif // OPAL_Distribution_HH
