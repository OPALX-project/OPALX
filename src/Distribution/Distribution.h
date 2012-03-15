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

#include "AbstractObjects/Definition.h"
#include "Algorithms/PartData.h"
#include "Algorithms/PartBunch.h"
#include "Algorithms/bet/EnvelopeBunch.h"
#include "AbstractObjects/Expressions.h"
#include "AbstractObjects/OpalData.h"
#include "Attributes/Attributes.h"
#include "Attributes/String.h"
#include "Attributes/Real.h"
#include "Expressions/SAutomatic.h"
#include "Expressions/SRefExpr.h"
#include "Physics/Physics.h"
#include "Utilities/OpalException.h"

#include "AbstractObjects/OpalData.h"
#include "AbstractObjects/ObjectFunction.h"
#include "Attributes/Real.h"
#include "Attributes/String.h"
#include "Attributes/RealArray.h"

#include "ValueDefinitions/StringConstant.h"
#include "ValueDefinitions/RealVariable.h"

#include "Structure/BoundaryGeometry.h"

#include "Ippl.h"

#include <hdf5.h>
#include "H5Part.h"

#include <string>
#include <map>

#include "Algorithms/PartBins.h"


#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_qrng.h>

#include "halton1d_sequence.hh"

#define RANLIBX
#include "ranlib.h"

#define sqr(x) pow(x,2.)

enum DistrTypeT {GAUSS, BINOMIAL, UNIFORMXYZ, FROMFILE, GUNGAUSSFLATTOPTH, GUNUNIFORM, GUNGAUSS3D, NODIST, SURFACEEMISSION, SURFACERANDCREATE, ASTRAFLATTOPTH};


class LaserProfile;

// Class Distribution
// ------------------------------------------------------------------------
/// The Distribution definition.
//  A Distribution definition is used by most physics commands to define the
//  particle charge and the reference momentum, together with some other
//  data.

class Distribution: public Definition {

public:

    /// Exemplar constructor.
    Distribution();

    virtual ~Distribution();

    /// Test if replacement is allowed.
    //  Can replace only by another Distribution.
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual Distribution *clone(const string &name);

    /// Check the Distribution data.
    virtual void execute();

    /// Find named Distribution.
    static Distribution *find(const string &name);

    /// Return the embedded CLASSIC PartData.
    const PartData &getReference() const;

    /// Update the Distribution data.
    virtual void update();

    /// configure rng generator such that calls to sample
    /// gives the random values
    void setup(PartBunch &beam, size_t Np, bool scan);

    /// gets back (x,p)
    pair<Vector_t, Vector_t> sample(double dt, int binNumber);

    void create(PartBunch &p, size_t Np);
    void create(PartBunch &p, size_t Np, bool scan);

    void create(PartBunch &p, BoundaryGeometry &bg);

    void createPriPart(PartBunch *beam, BoundaryGeometry &bg);

    void createSlicedBunch(int sl, double charge, double gamma, double mass, double current, double center, double Bz0, EnvelopeBunch *p);

    void doRestart(PartBunch &p, size_t Np, size_t restartStep);

    void doRestart_cycl(PartBunch &p, size_t Np, size_t restartStep, const int specifiedNumBunch);

    void doRestartEnvelope(EnvelopeBunch &p, size_t Np, size_t restartStep);

    /// Print the TFS descriptors for the beam.
    void tfsDescriptors(std::ostream &os) const;

    Inform &print(Inform &os) const;

    double getTBin() { return tBin_m; }
    double getTEmission();

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
    string getTypeofDistribution();
    
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

    void createUniformTUniformL(Vector_t emit, Vector_t alpha, Vector_t beta, Vector_t gamma,
                                PartBunch &beam, size_t particles, bool isBinned);

    void binnDistributionZ(PartBunch &beam, size_t Np, string distType);
    void binnDistributionT(PartBunch &beam, size_t Np, string distType);

    void binnDistributionFromFile(PartBunch &beam, const string fn);

    double eVtoBetaGamma(const double &valueIneV, const double &mass) {
        double tmp = 1. + valueIneV / mass;
        return sqrt(tmp * tmp - 1.);
    }

    double betaGammatoeV(const double &valueInbega, const double &mass) {
        return (sqrt(valueInbega * valueInbega + 1.) - 1.) * mass ;
    }

    // Not implemented.
    Distribution(const Distribution &);
    void operator=(const Distribution &);

    // Clone constructor.
    Distribution(const string &name, Distribution *parent);

    // The particle reference data.
    PartData reference;

    // The structure for particle binning
    PartBins *pbin_m;

    // If we scan we already have a particles and to not have to create them anymore
    bool scan_m;


    //`setup() fills the following variables

    string distT_m;
    DistrTypeT distrTypeT_m;   // will replace distT_m

    double corr_m[7];

    Vector_t sigx_m;
    Vector_t sigp_m;

    Vector_t binc_m;

    Vector_t emit_m;
    Vector_t alpha_m;
    Vector_t beta_m;
    Vector_t gamma_m;

    double Hs2a_m;
    double Hs2b_m;
    double gauss_offset_m[2];

    double Vs2a_m;
    double Vs2b_m;

    double Ls2a_m;
    double Ls2b_m;

    double avrgpt_m;
    double avrgt_m;

    double nBins_m;
    double sBins_m;
    gsl_rng* rn_m;
    gsl_qrng* R_m;

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
    string laserProfileFn_m;
    string laserImage_m;
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
    return d.print(os);
}



#endif // OPAL_Distribution_HH
