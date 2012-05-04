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
#include <string>
#include <map>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_qrng.h>

#include "H5hut.h"

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"
#include "ranlib.h"
#include "Algorithms/Vektor.h"

#define RANLIBX
#define sqr(x) x*x

class PartBunch;
class PartBins;
class EnvelopeBunch;
class BoundaryGeometry;
class LaserProfile;

enum DistrTypeT {NODIST,
                 GAUSS,
                 FROMFILE,
                 GUNGAUSSFLATTOPTH,
                 BINOMIAL,
                 SURFACEEMISSION,
                 SURFACERANDCREATE,
                 ASTRAFLATTOPTH
                };

// Class Distribution
// ------------------------------------------------------------------------
/// The Distribution definition.
//
//  A Distribution definition is used by most physics commands to define the
//  beam distribution. This includes particle charge, reference momentum, beam size,
//  beam momentum spread etc.

class Distribution: public Definition {

    //===============
    // Class methods.
    //===============
public:

    /// Exemplar constructor.
    Distribution();

    virtual ~Distribution();

    /// Test if replacement is allowed.
    //  Can replace only by another Distribution.
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual Distribution *clone(const std::string &name);

    /// Check the Distribution data.
    virtual void execute();

    /// Find named Distribution.
    static Distribution *find(const std::string &name);

    /// Return the embedded CLASSIC PartData.
    const PartData &getReference() const;

    /// Update the Distribution data.
    virtual void update();

    /// configure rng generator such that calls to sample
    /// gives the random values
    void setup(PartBunch &beam, size_t Np, bool scan);

    /**
     * Add a list of distributions. This takes a list
     * of distributions and adds them together. Currently
     * this only works for the GUNGAUSSFLATTOPTH distribution.
     * Also, the transverse properties of the beam are defined only
     * by the first distribution in the list. For subsequent
     * distributions, only the longitudinal properties are used.
     */
    bool addDistributions(PartBunch &beam, std::vector<Distribution *> distributions, size_t numberOfParticles);

    /// gets back (x,p)
    std::pair<Vector_t, Vector_t> sample(double dt, int binNumber);

    /// gets back (x,p)
    std::pair<Vector_t, Vector_t> sampleNEW(double dt, int binNumber);

    void create(PartBunch &p, size_t Np);
    void create(PartBunch &p, size_t Np, bool scan);

    void create(PartBunch &p, BoundaryGeometry &bg);

    void createPriPart(PartBunch *beam, BoundaryGeometry &bg);

    void createSlicedBunch(int sl, double charge, double gamma, double mass, double current, double center, double Bz0, EnvelopeBunch *p);

    void doRestart(PartBunch &p, size_t Np, int restartStep);

    void doRestart_cycl(PartBunch &p, size_t Np, int restartStep, const int specifiedNumBunch);

    void doRestartEnvelope(EnvelopeBunch &p, size_t Np, int restartStep);

    Inform &printInfo(Inform &os) const;

    double getTBin() { return tBin_m; }
    double getTEmission();

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

private:
    void sampleGauss(PartBunch &beam, size_t Np);

    void createTimeBins(const int Np);

    void createBinom(Vector_t emit, Vector_t alpha, Vector_t beta, Vector_t gamma,
                     Vector_t bincoef, PartBunch &beam, size_t particles,
                     bool isBinned);

    void binnDistribution(PartBunch &beam, size_t Np, std::string distType);

    void binnDistributionFromFile(PartBunch &beam, const std::string fn);

    double eVtoBetaGamma(const double &valueIneV, const double &mass) {
        double tmp = 1. + valueIneV / mass;
        return sqrt(tmp * tmp - 1.);
    }

    double betaGammatoeV(const double &valueInbega, const double &mass) {
        return (sqrt(valueInbega * valueInbega + 1.) - 1.) * mass ;
    }

    void writeToFile();



    // Not implemented.
    Distribution(const Distribution &);
    void operator=(const Distribution &);

    // Clone constructor.
    Distribution(const std::string &name, Distribution *parent);

    //===============
    // Class members.
    //===============
private:

    // The particle reference data.
    PartData reference;

    // The structure for particle binning
    PartBins *pbin_m;

    // If we scan we already have a particles and to not have to create them anymore
    bool scan_m;


    //`setup() fills the following variables

    std::string distT_m;
    DistrTypeT distrTypeT_m;   // will replace distT_m

    double corr_m[7];

    Vector_t sigx_m;
    Vector_t sigp_m;

    Vector_t binc_m;

    Vector_t emit_m;
    Vector_t alpha_m;
    Vector_t beta_m;
    Vector_t gamma_m;

    double gauss_offset_m[2];

    double avrgpt_m;
    double avrgt_m;

    int nBins_m; // Number of energy bins the distribution is broken up into.
    int sBins_m; // Number of samples used to create the time histogram of each energy bin.
    gsl_rng *rn_m;
    gsl_qrng *R_m;

    /// this is for the transverse dimension
    gsl_qrng *qrng_m;

    RANLIB_class *rGen_m;

    /// time binned distribution with thermal energy
    double transvCutOff_m;
    double tEmission_m;
    double tPulseLengthFWHM_m;
    double tRise_m;
    double tFall_m;
    double sigmaRise_m;
    double sigmaFall_m;
    double cutoff_m;
    double tBin_m;

    gsl_histogram *h_m;
    double *distributionTable_m;

    bool   astraMode_m;
    double workf_m;          // eV
    double siglaser_m;      // m
    double elaser_m;        // eV
    double fe_m;            // Fermi energy eV
    double ag_m;            // Acceleration gradient eV/m
    double ekin_m;          // eV
    double phimax_m;        // rad
    double schottky_m;;      // eV
    double ptot_m;          // beta gamma

    std::ofstream os_m;

    /**
       Data for Laser Profile read in
    */
    std::string laserProfileFn_m;
    std::string laserImage_m;
    double intensityCut_m;

    LaserProfile *lp_m;


    /**
       Data for Dark Current calculations
    */
    size_t darkCurrentParts_m;
    /**
       Data for Dark Current initialized positions. A little bit inward along the triangle normal, positive: inside the geometry.
    */
    double darkInwardMargin_m;

    double eInitThreshold_m;

    double workFunction_m;

    double fieldEnhancement_m;

    double fieldThrFN_m;

    size_t maxFN_m;

    double paraFNA_m;

    double paraFNB_m;

    double paraFNY_m;

    double paraFNVYSe_m;

    double paraFNVYZe_m;

    int    secondaryFlag_m;

    double ppVw_m;

    double vVThermal_m;



};

inline Inform &operator<<(Inform &os, const Distribution &d) {
    return d.printInfo(os);
}

#endif // OPAL_Distribution_HH
