#ifndef distr_HH
#define distr_HH

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_erf.h>
#include <Ippl.h>

#include "AbstractObjects/OpalData.h"
#include "Physics/Physics.h"
#include "Algorithms/PartBins.h"

#define RANLIBX
#include "Distribution/ranlib.h"
#include "Distribution/LaserProfile.h"
#include "disttest.h"
#include "myPartBunch.h"

class Distribution {

public:

    /// Exemplar constructor.
    Distribution(DistrTypeT actDist, size_t Np, double dEBins, int nBins, Vector_t distCutOff, 
		 bool hasLaserProfile,
		 string laserProfileFn,
		 string laserImage,
		 double intensityCut,
		 double tpulsefwhm, double trise, double tfall, Vector_t m,
		 double ekin,
		 Vektor<double,7>  corr,
		 Vector_t r,
		 Vector_t p,
		 double avrgt,
		 double avrgpt);

    ~Distribution();

    /// Return the embedded CLASSIC PartData.
    const PartData &getReference() const;

    void create(myPartBunch &p, size_t Np);

    Inform &print(Inform &os) const;

    bool isBinned() { return pbin_m != NULL; }
    PartBins *getPBin() { return pbin_m; }
    double getTBin() { return tBin_m;}

    pair<Vector_t, Vector_t> sample(int binNumber);
    pair<Vector_t, Vector_t> Distribution::sampleBinom();
    pair<Vector_t, Vector_t> Distribution::sampleGauss();

private:
    void sampleGauss(myPartBunch &beam, size_t Np);

    void createTimeBins(const int Np);

    void createBinom(myPartBunch &beam, size_t particles);


    inline Vector_t calcThermalEmittance (double pTot); 


    double eVtoBetaGamma(const double &valueIneV, const double &mass) {
        double tmp = 1. + valueIneV / mass;
        return sqrt(tmp * tmp - 1.);
    }

    double betaGammatoeV(const double &valueInbega, const double &mass) {
        return (sqrt(valueInbega * valueInbega + 1.) - 1.) * mass ;
    }

    PartBins *pbin_m;     // The structure for particle binning
	bool doEmission_m;
    double nBins_m;

    DistrTypeT distrType_m; 

    Vektor<double,7> corr_m;
    Vector_t sigx_m;
    Vector_t sigp_m;

    Vector_t binc_m;

    double gauss_offset_m[2];

    Vector_t distCutOff_m;


    /** 
	Temporary variables for binomial distribution
    */

    Vector_t M_m;
    Vector_t PM_m;
    Vector_t COSCHI_m;
    Vector_t SINCHI_m;
    Vector_t CHI_m;
    Vector_t AMI_m;
    Vector_t L_m;
    Vector_t PL_m;

    Vector_t emit_m;
    Vector_t alpha_m; 
    Vector_t beta_m; 
    Vector_t gamma_m;

    double pt_m;

    double avrgpt_m;
    double avrgt_m;


  
    double tEmission_m;
    double tPulseLengthFWHM_m;
    double tRise_m;
    double tFall_m;
    double sigmaRise_m;
    double sigmaFall_m;
    double cutoff_m;
    double tBin_m;

    double workf_m;         // eV
    double siglaser_m;      // m
    double elaser_m;        // eV
    double fe_m;            // Fermi energy eV
    double ag_m;            // Acceleration gradient eV/m
    double ekin_m;          // eV

    double schottky_m;      // eV
    double ptot_m;          // beta gamma

    RANLIB_class *rGen_m;

    gsl_histogram *h_m;

    /**
       Data for Laser Profile read in
    */
    string laserProfileFn_m;
    string laserImage_m;
    double intensityCut_m;

    LaserProfile *laserProfile_m;
};

inline Inform &operator<<(Inform &os, const Distribution &d) {
    return d.print(os);
}


#endif

