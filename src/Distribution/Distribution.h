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

#define RANLIBX
#include "ranlib.h"

#define sqr(x) pow(x,2.)

enum DistrTypeT {GAUSS, BINOMIAL, UNIFORMXYZ, FROMFILE, GUNGAUSSFLATTOPTH, GUNUNIFORM, GUNGAUSS3D, NODIST, SURFACEEMISSION};

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
    pair<Vector_t, Vector_t> sample(double dt);

    void create(PartBunch &p, size_t Np);
    void create(PartBunch &p, size_t Np, bool scan);

    void create(PartBunch &p, BoundaryGeometry &bg);

    void createSlicedBunch(double charge, double gamma, double mass, double current, double center, double Bz0, EnvelopeBunch *p);

    void doRestart(PartBunch &p, size_t Np, size_t restartStep);

    void doRestart_cycl(PartBunch &p, size_t Np, size_t restartStep, const int specifiedNumBunch);

    /// Print the TFS descriptors for the beam.
    void tfsDescriptors(std::ostream &os) const;

    Inform &print(Inform &os) const;

    double getTBin() { return tBin_m; }

    size_t getNumberOfDarkCurrentParticles();
    double getDarkCurrentParticlesInwardMargin();
    double getWorkFunction();
    double getFieldEnhancement();
    size_t getMaxFNemissionPartPerTri();

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

    Vector_t corr_m;
    Vector_t sigx_m;
    Vector_t sigp_m;

    Vector_t binc_m;

    Vector_t emit_m;
    Vector_t alpha_m;
    Vector_t beta_m;
    Vector_t gamma_m;

    double gauss_corr_m[7];
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

    double workFunction_m;

    double fieldEnhancement_m;

    size_t maxFN_m;

};

inline Inform &operator<<(Inform &os, const Distribution &d) {
    return d.print(os);
}



#endif // OPAL_Distribution_HH
