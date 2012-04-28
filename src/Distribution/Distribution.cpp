// ------------------------------------------------------------------------
// $RCSfile: Distribution.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.3.4.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Distribution
//   The class for the OPAL Distribution command.
//
// ------------------------------------------------------------------------
#include "Distribution/Distribution.h"
#include "AbstractObjects/Expressions.h"
#include "Attributes/Attributes.h"
#include "Utilities/Options.h"
#include "halton1d_sequence.hh"
#include "AbstractObjects/OpalData.h"
#include "Algorithms/PartBunch.h"
#include "Algorithms/PartBins.h"
#include "Algorithms/bet/EnvelopeBunch.h"
#include "Structure/BoundaryGeometry.h"
#include "Algorithms/PartBinsCyc.h"
#include "BasicActions/Option.h"
#include "Distribution/LaserProfile.h"

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_erf.h>

#include <cmath>
#include <cfloat>
#include <iomanip>   // Needed for I/O manipulators
#include <iostream>  // Neeeded for stream I/O
#include <string>
#include <vector>

using namespace std;
using namespace Expressions;
using namespace Physics;
using namespace Attributes;

extern Inform *gmsg;

//
// Class Distribution
// ------------------------------------------------------------------------

// The attributes of class Distribution.

namespace {
    enum {
        // DESCRIPTION OF THE DISTRIBUTION:
        DISTRIBUTION,
        FNAME,
        LASERPROFFN,
        IMAGENAME,
        INTENSITYCUT,
        XMULT,
        YMULT,
        TMULT,
        ZMULT,
        PXMULT,
        PYMULT,
        PTMULT,
        PZMULT,
        BETAX,
        BETAY,
        ALPHAX,
        ALPHAY,
        MX,
        MY,
        MT,
        MZ,
        DX,
        DDX,
        DY,
        DDY,
        R51,
        R52,
        R61,
        R62,
        PT,
        PZ,
        T,
        Z,
        SIGMAX,
        SIGMAY,
        SIGMAT,
        SIGMAZ,
        TRANSVCUTOFF,
        SIGMAPX,
        SIGMAPY,
        SIGMAPT,
        SIGMAPZ,
        TRISE,
        TFALL,
        CUTOFF,
        TPULSEFWHM,
        FTOSCAMPLITUDE,
        FTOSCPERIODS,
        WEIGHT,
        CORRX,
        CORRY,
        CORRT,
        CORRZ,
        OFFSETX,
        OFFSETY,
        TEMISSION,
        NBIN,
        SBIN,
        DEBIN,
        ELASER,
        SIGLASER,
        W,
        FE,
        AG,
        EKIN,
        NPDARKCUR,
        EINITHR,
        INWARDMARGIN,
        FNA,
        FNB,
        FNY,
        FNVYZERO,
        FNVYSECOND,
        FNPHIW,
        FNBETA,
        FNFIELDTHR,
        FNMAXEMI,
        SECONDARYFLAG,
        NEMISSIONMODE,
        VSEYZERO,// sey_0 in Vaughan's model
        VEZERO,// energy related to sey_0 in Vaughan's model
        VSEYMAX,// sey max in Vaughan's model
        VEMAX,// Emax in Vaughan's model
        VKENERGY,// fitting parameter denotes the roughness of surface for impact energy in Vaughan's model
        VKTHETA,// fitting parameter denotes the roughness of surface for impact angle in Vaughan's model
        VVTHERMAL,// thermal velocity of Maxwellian distribution of secondaries in Vaughan's model
        VW,
        SURFMATERIAL, // Add material type, currently 0 for copper and 1 for stainless steel.
        SIZE
    };
}

/**
 * Constructor
 *
 */
Distribution::Distribution():
    Definition(SIZE, "DISTRIBUTION", "The DISTRIBUTION statement defines data for the 6D particle distr."),
    distrTypeT_m(NODIST),
    rn_m(NULL),
    R_m(NULL),
    qrng_m(NULL),
    distributionTable_m(NULL) {
    itsAttr[DISTRIBUTION] = makeString("DISTRIBUTION", "Distribution type: GAUSS, BINOMIAL, FROMFILE,"
                                       "GUNGAUSSFLATTOPTH, ASTRAFLATTOPTH, SURFACEEMISSION, SURFACERANDCREATE", "GAUSS");

    itsAttr[FNAME] = makeString("FNAME", "File for reading in 6D particle coordinates");

    itsAttr[LASERPROFFN] = makeString("LASERPROFFN", "File for read in a measured laser profile (x,y)", "");
    itsAttr[IMAGENAME] = makeString("IMAGENAME", "Name of the image");
    itsAttr[INTENSITYCUT] = makeReal("INTENSITYCUT", "For background substraction, in % of max intensity", 0.0);


    itsAttr[XMULT] = makeReal("XMULT", "Multiplier for X", 1.0);
    itsAttr[YMULT] = makeReal("YMULT", "Multiplier for Y", 1.0);
    itsAttr[TMULT] = makeReal("TMULT", "Multiplier for T", 1.0);
    itsAttr[ZMULT] = makeReal("TMULT", "Multiplier for T", -99.0);
    itsAttr[TRANSVCUTOFF] = makeReal("TRANSVCUTOFF", "Transverse cut-off in units of sigma", 3.0);

    itsAttr[PXMULT] = makeReal("PXMULT", "Multiplier for PX", 1.0);
    itsAttr[PYMULT] = makeReal("PYMULT", "Multiplier for PY", 1.0);
    itsAttr[PTMULT] = makeReal("PTMULT", "Multiplier for PT", 1.0);
    itsAttr[PZMULT] = makeReal("PZMULT", "Multiplier for PZ", -99.0);

    itsAttr[ALPHAX] = makeReal("ALPHAX", "Courant Synder parameter", 1.0);
    itsAttr[ALPHAY] = makeReal("ALPHAY", "Courant Synder parameter", 1.0);

    itsAttr[BETAX] = makeReal("BETAX", "Courant Synder parameter", -1.0);
    itsAttr[BETAY] = makeReal("BETAY", "Courant Synder parameter", 1.0);

    itsAttr[MX]    = makeReal("MX", "Defines the distribution in x, 0+eps .. inf", 1.0);
    itsAttr[MY]    = makeReal("MY", "Defines the distribution in y, 0+eps .. inf", 1.0);
    itsAttr[MT]    = makeReal("MT", "Defines the distribution in t, 0+eps .. inf", 1.0);
    itsAttr[MZ]    = makeReal("MZ", "Defines the distribution in z, 0+eps .. inf", -99.0);

    itsAttr[DX]    = makeReal("DX", "Dispersion in x (R16 in Transport notation)", 0.0);
    itsAttr[DDX]   = makeReal("DDX", "First derivative of Dx", 0.0);

    itsAttr[DY]    = makeReal("DY", "DY", 0.0);
    itsAttr[DDY]   = makeReal("DDY", "DDY", 0.0);

    itsAttr[R51]    = makeReal("R51", "R51", 0.0);
    itsAttr[R52]   = makeReal("R52", "R52", 0.0);

    itsAttr[R61]    = makeReal("R61", "R61", 0.0);
    itsAttr[R62]   = makeReal("R62", "R62", 0.0);

    itsAttr[PT] = makeReal("PT", "average longitudinal momentum", 0.0);
    itsAttr[PZ] = makeReal("PZ", "average longitudinal momentum", -99.0);
    itsAttr[T] = makeReal("T", "average longitudinal position", 0.0);
    itsAttr[Z] = makeReal("Z", "average longitudinal position", -99.0);

    itsAttr[SIGMAX] = makeReal("SIGMAX", "SIGMAx (m)", 1.0e-2);
    itsAttr[SIGMAY] = makeReal("SIGMAY", "SIGMAy (m)", 1.0e-2);
    itsAttr[SIGMAT] = makeReal("SIGMAT", "SIGMAt (m)", 1.0e-2);
    itsAttr[SIGMAZ] = makeReal("SIGMAZ", "SIGMAz (m)", -99.0);

    itsAttr[SIGMAPX] = makeReal("SIGMAPX", "SIGMApx", 0.0);
    itsAttr[SIGMAPY] = makeReal("SIGMAPY", "SIGMApy", 0.0);
    itsAttr[SIGMAPT] = makeReal("SIGMAPT", "SIGMApt", 0.0);
    itsAttr[SIGMAPZ] = makeReal("SIGMAPZ", "SIGMApz", -99.0);

    itsAttr[CORRX] = makeReal("CORRX", "CORRx", -0.5);
    itsAttr[CORRY] = makeReal("CORRY", "CORRy", 0.5);
    itsAttr[CORRT] = makeReal("CORRT", "CORRt", 0.0);
    itsAttr[CORRZ] = makeReal("CORRZ", "CORRz", -99.0);

    itsAttr[OFFSETX] = makeReal("OFFSETX", "OFFSETx", 0.0);
    itsAttr[OFFSETY] = makeReal("OFFSETY", "OFFSETy", 0.0);

    itsAttr[TEMISSION] = makeReal("TEMISSION", "Time in seconds in which we have emission",  0.0);
    itsAttr[NBIN]      = makeReal("NBIN", "In case of emission how many energy bins should we use", 0.0);
    itsAttr[SBIN]      = makeReal("SBIN", "In case of emission how many sample bins an energy bin should use", 100.0);
    itsAttr[DEBIN]     = makeReal("DEBIN", "Energy band for a bin in keV, defines the rebinning", 1000000.0);

    itsAttr[TPULSEFWHM]  = makeReal("TPULSEFWHM", "Pulse FWHM (s)", 0.0);
    itsAttr[TRISE]       = makeReal("TRISE", "Rise time for GUNGAUSSFLATTOP distribution type (s)", 0.0);
    itsAttr[TFALL]       = makeReal("TFALL", "Fall time for GUNGAUSSFLATTOP distribution type (s)", 0.0);
    itsAttr[CUTOFF]      = makeReal("CUTOFF", "Cutoff for GUNGAUSSFLATTOP distribution type in sigmas", 3.0);

    itsAttr[FTOSCAMPLITUDE] = makeReal("FTOSCAMPLITUDE", "Amplitude of oscillations superimposed on flat top portion of GUNGAUSSFLATTOPTH distribtuion (in percent of flat top amplitude)", 0.0);
    itsAttr[FTOSCPERIODS]    = makeReal("FTOSCPERIODS", "Number of oscillations superimposed on flat top portion of GUNGAUSSFLATTOPTH distribution", 0.0);

    itsAttr[WEIGHT] = makeReal("WEIGHT", "Weight of this distribution when used in a distribution list.", 1.0);

    itsAttr[ELASER] = makeReal("ELASER", "Laser energy (eV)", 0.0);
    itsAttr[SIGLASER] = makeReal("SIGLASER", "Sigma of (uniform) laser spot size (m)", 0.0);
    itsAttr[W] = makeReal("W", "Workfunction of material (eV)", 0.0);
    itsAttr[FE] = makeReal("FE", "Fermi energy (eV)", 0.0);
    itsAttr[AG] = makeReal("AG", "Acceleration Gradient (MV/m)", 0.0);

    itsAttr[EKIN] = makeReal("EKIN", "Ekin used in ASTRA (eV)", -1.0);

    itsAttr[NPDARKCUR] = makeReal("NPDARKCUR", "Number of dark current particles", 1000.0);
    itsAttr[INWARDMARGIN] = makeReal("INWARDMARGIN", "Inward margin of initialized dark current particle positions", 0.001);
    itsAttr[EINITHR] = makeReal("EINITHR", "E field threshold (MV), only in position r with E(r)>EINITHR that particles will be initialized", 0.0);
    itsAttr[FNA] = makeReal("FNA", "Empirical constant A for Fowler-Nordheim emission model", 1.54e-6);
    itsAttr[FNB] = makeReal("FNB", "Empirical constant B for Fowler-Nordheim emission model", 6.83e9);
    itsAttr[FNY] = makeReal("FNY", "Constant for image charge effect parameter y(E) in Fowler-Nordheim emission model", 3.795e-5);
    itsAttr[FNVYZERO] = makeReal("FNVYZERO", "Zero order constant for v(y) function in Fowler-Nordheim emission model", 0.9632);
    itsAttr[FNVYSECOND] = makeReal("FNVYSECOND", "Second order constant for v(y) function in Fowler-Nordheim emission model", 1.065);
    itsAttr[FNPHIW] = makeReal("FNPHIW", "Work function of gun surface material (eV)", 4.65);
    itsAttr[FNBETA] = makeReal("FNBETA", "Field enhancement factor for Fowler-Nordheim emission", 50.0);
    itsAttr[FNFIELDTHR] = makeReal("FNFIELDTHR", "Field threshold for Fowler-Nordheim emission (MV/m)", 30.0);
    itsAttr[FNMAXEMI] = makeReal("FNMAXEMI", "Maximum number of electrons emitted from a single triangle for Fowler-Nordheim emission", 20.0);
    itsAttr[SECONDARYFLAG] = makeReal("SECONDARYFLAG", "Select the secondary model type(0:no secondary emission; 1:Furman-Pivi; 2 or larger: Vaughan's model", 0);
    itsAttr[NEMISSIONMODE] = makeBool("NEMISSIONMODE", "Secondary emission mode type(true: emit n true secondaries; false: emit one particle with n times charge", true);
    itsAttr[VSEYZERO] = makeReal("VSEYZERO", "Sey_0 in Vaughan's model", 0.5);
    itsAttr[VEZERO] = makeReal("VEZERO", "Energy related to sey_0 in Vaughan's model in eV", 12.5);
    itsAttr[VSEYMAX] = makeReal("VSEYMAX", "Sey max in Vaughan's model", 2.22);
    itsAttr[VEMAX] = makeReal("VEMAX", "Emax in Vaughan's model in eV", 165);
    itsAttr[VKENERGY] = makeReal("VKENERGY", "Fitting parameter denotes the roughness of surface for impact energy in Vaughan's model", 1.0);
    itsAttr[VKTHETA] = makeReal("VKTHETA", "Fitting parameter denotes the roughness of surface for impact angle in Vaughan's model", 1.0);
    itsAttr[VVTHERMAL] = makeReal("VVTHERMAL", "Thermal velocity of Maxwellian distribution of secondaries in Vaughan's model", 7.268929821 * 1e5); // electrons 1.5eV default
    itsAttr[VW] = makeReal("VW", "VW denote the velocity scalar for Parallel plate benchmark", 1.0);
    itsAttr[SURFMATERIAL] = makeReal("SURFMATERIAL", "Material type number of the cavity surface for Furman-Pivi's model, 0 for copper, 1 for stainless steel", 0);

    // Set up default beam.
    Distribution *defDistribution = clone("UNNAMED_Distribution");
    defDistribution->builtin = true;

    try {
        defDistribution->update();
        OpalData::getInstance()->define(defDistribution);
    } catch(...) {
        delete defDistribution;
    }
    pbin_m = NULL;
    lp_m = NULL;

    darkCurrentParts_m = (size_t) Attributes::getReal(itsAttr[NPDARKCUR]);
    darkInwardMargin_m = Attributes::getReal(itsAttr[INWARDMARGIN]);
    eInitThreshold_m = Attributes::getReal(itsAttr[EINITHR]);

    workFunction_m = Attributes::getReal(itsAttr[FNPHIW]);
    fieldEnhancement_m = Attributes::getReal(itsAttr[FNBETA]);
    maxFN_m = (size_t) Attributes::getReal(itsAttr[FNMAXEMI]);
    fieldThrFN_m = Attributes::getReal(itsAttr[FNFIELDTHR]);
    paraFNA_m = Attributes::getReal(itsAttr[FNA]);
    paraFNB_m = Attributes::getReal(itsAttr[FNB]);
    paraFNY_m = Attributes::getReal(itsAttr[FNY]);
    paraFNVYZe_m = Attributes::getReal(itsAttr[FNVYZERO]);
    paraFNVYSe_m = Attributes::getReal(itsAttr[FNVYSECOND]);

    secondaryFlag_m = Attributes::getReal(itsAttr[SECONDARYFLAG]);
    ppVw_m = Attributes::getReal(itsAttr[VW]);
    vVThermal_m = Attributes::getReal(itsAttr[VVTHERMAL]);
    tEmission_m = -1.0;
    qrng_m = NULL;
}
/**
 *
 *
 * @param name
 * @param parent
 */
Distribution::Distribution(const string &name, Distribution *parent):
    Definition(name, parent),
    reference(parent->reference),
    pbin_m(NULL),
    distrTypeT_m(NODIST),
    rn_m(parent->rn_m),
    R_m(parent->R_m),
    qrng_m(parent->qrng_m),
    tEmission_m(parent->tEmission_m),
    distributionTable_m(parent->distributionTable_m),
    lp_m(NULL),
    darkCurrentParts_m(parent->darkCurrentParts_m),
    darkInwardMargin_m(parent->darkInwardMargin_m),
    eInitThreshold_m(parent->eInitThreshold_m),
    workFunction_m(parent->workFunction_m),
    fieldEnhancement_m(parent->fieldEnhancement_m),
    fieldThrFN_m(parent->fieldThrFN_m),
    maxFN_m(parent->maxFN_m),
    paraFNA_m(parent-> paraFNA_m),
    paraFNB_m(parent-> paraFNB_m),
    paraFNY_m(parent-> paraFNY_m),
    paraFNVYSe_m(parent-> paraFNVYSe_m),
    paraFNVYZe_m(parent-> paraFNVYZe_m),
    secondaryFlag_m(parent->secondaryFlag_m) {

    darkCurrentParts_m = (size_t) Attributes::getReal(itsAttr[NPDARKCUR]);
    darkInwardMargin_m = Attributes::getReal(itsAttr[INWARDMARGIN]);
    eInitThreshold_m = Attributes::getReal(itsAttr[EINITHR]);

    workFunction_m = Attributes::getReal(itsAttr[FNPHIW]);
    fieldEnhancement_m = Attributes::getReal(itsAttr[FNBETA]);
    maxFN_m = (size_t) Attributes::getReal(itsAttr[FNMAXEMI]);
    fieldThrFN_m = Attributes::getReal(itsAttr[FNFIELDTHR]);
    paraFNA_m = Attributes::getReal(itsAttr[FNA]);
    paraFNB_m = Attributes::getReal(itsAttr[FNB]);
    paraFNY_m = Attributes::getReal(itsAttr[FNY]);
    paraFNVYZe_m = Attributes::getReal(itsAttr[FNVYZERO]);
    paraFNVYSe_m = Attributes::getReal(itsAttr[FNVYSECOND]);

    secondaryFlag_m = Attributes::getReal(itsAttr[SECONDARYFLAG]);
    ppVw_m = Attributes::getReal(itsAttr[VW]);
    vVThermal_m = Attributes::getReal(itsAttr[VVTHERMAL]);
}

/**
 * Destructor
 *
 */
Distribution::~Distribution() {
    if(pbin_m) {
        delete pbin_m;
        pbin_m = NULL;
    }

    if((Ippl::getNodes() == 1) && (os_m.is_open()))
        os_m.close();

    if(lp_m) {
        delete lp_m;
        lp_m = NULL;
    }

    if(distributionTable_m)
        delete[] distributionTable_m;

    if(rn_m) {
        gsl_rng_free(rn_m);
        gsl_qrng_free(R_m);
    }
    if(qrng_m)
        gsl_qrng_free(qrng_m);
}

/**
 * At the moment only write the header into the file dist.dat
 * PartBunch will then append (very uggly)
 * @param 
 * @param 
 * @param 
 */
void Distribution::writeToFile() {

    if(Ippl::getNodes() == 1) {
        if (os_m.is_open()) {
            ;
        }
        else {
            *gmsg << " Write distribution to file dist.dat" << endl;
            string file("data/dist.dat");
            os_m.open(file.c_str());
            if(os_m.bad()) {
                *gmsg << "Unable to open output file " <<  file << endl;
            }
            os_m << "# x y ti px py pz Ekin= " << ekin_m << " [eV] " << endl;
            os_m.close();
        }
    }
}

/**
 * This is the main entrypoint, called from PartBunch::setDistribution
 * The envelope trackes is doing it differently by calling createSliceBunch
 * @param beam
 * @param Np
 * @param scan
 */
void Distribution::setup(PartBunch &beam, size_t Np, bool scan) {

    scan_m = scan;
    nBins_m = static_cast<int>(fabs(Attributes::getReal(itsAttr[NBIN])));

    bool isBinned = (nBins_m > 0);

    if(isBinned) {
        if(pbin_m)
            delete pbin_m;
        pbin_m = new PartBins(static_cast<int>(fabs(Attributes::getReal(itsAttr[NBIN]))),
                              static_cast<int>(fabs(Attributes::getReal(itsAttr[SBIN]))));
    } else
        pbin_m = NULL;

    if(scan_m) {
        beam.destroy(beam.getLocalNum(), 0);
        beam.update();
        INFOMSG("In scan mode: delete all particles in the bunch" << endl;);
    }

    laserProfileFn_m = Attributes::getString(itsAttr[LASERPROFFN]);

    if(!(laserProfileFn_m == string(""))) {
        laserImage_m  = Attributes::getString(itsAttr[IMAGENAME]);
        intensityCut_m = Attributes::getReal(itsAttr[INTENSITYCUT]);
        lp_m = new LaserProfile(laserProfileFn_m, laserImage_m, intensityCut_m);
    }

    beam.setTEmission(Attributes::getReal(itsAttr[TEMISSION]));
    beam.setNumBunch(1);

    distT_m = Attributes::getString(itsAttr[DISTRIBUTION]);
    if(distT_m == "GAUSS")
        distrTypeT_m = GAUSS;
    else if(distT_m == "GUNGAUSSFLATTOPTH")
        distrTypeT_m = GUNGAUSSFLATTOPTH;
    else if(distT_m == "ASTRAFLATTOPTH")
        distrTypeT_m = ASTRAFLATTOPTH;
    else if(distT_m == "FROMFILE")
        distrTypeT_m = FROMFILE;
    else if(distT_m == "BINOMIAL")
        distrTypeT_m = BINOMIAL;
    else if(distT_m == "SURFACEEMISSION")
        distrTypeT_m = SURFACEEMISSION;
    else if(distT_m == "SURFACERANDCREATE")
        distrTypeT_m = SURFACERANDCREATE;

    /*
      Setup some common data:
    */
    if (distrTypeT_m == GUNGAUSSFLATTOPTH || distrTypeT_m == ASTRAFLATTOPTH) {
        corr_m[0] = Attributes::getReal(itsAttr[CORRX]);
        corr_m[1] = Attributes::getReal(itsAttr[CORRY]);
        corr_m[2] = Attributes::getReal(itsAttr[CORRT]);

        nBins_m = static_cast<int>(fabs(Attributes::getReal(itsAttr[NBIN])));
        sBins_m = static_cast<int>(fabs(Attributes::getReal(itsAttr[SBIN])));
        transvCutOff_m = Attributes::getReal(itsAttr[TRANSVCUTOFF]);
    
        sigx_m = Vector_t(Attributes::getReal(itsAttr[SIGMAX]),
                          Attributes::getReal(itsAttr[SIGMAY]),
                          Attributes::getReal(itsAttr[SIGMAT]));
        
        sigp_m = Vector_t(eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]), beam.getM()),
                          eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]), beam.getM()),
                          eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]), beam.getM()));

        tPulseLengthFWHM_m = Attributes::getReal(itsAttr[TPULSEFWHM]);
        cutoff_m = Attributes::getReal(itsAttr[CUTOFF]);
        tRise_m = Attributes::getReal(itsAttr[TRISE]);
        tFall_m = Attributes::getReal(itsAttr[TFALL]);
        double tratio = sqrt(2.0 * log(10.0)) - sqrt(2.0 * log(10.0 / 9.0));
        sigmaRise_m = tRise_m / tratio;
        sigmaFall_m = tFall_m / tratio;

        double dEBins = Attributes::getReal(itsAttr[DEBIN]);

        pbin_m->setRebinEnergy(dEBins);

        /*
          prepare quantities for thermal emittance calculation
        */
        workf_m = 0.0;         // eV
        siglaser_m = 0.0;      // m
        elaser_m = 0.0;        // eV
        fe_m = 0.0;            // Fermi energy eV
        ag_m = 0.0;            // Acceleration gradient eV/m
        ekin_m = 0.0;          // eV
        phimax_m = 0.0;        // rad
        schottky_m = 0.0;      // eV
        ptot_m = 0.0;          // beta gamma

        ekin_m = Attributes::getReal(itsAttr[EKIN]);
        ptot_m = eVtoBetaGamma(ekin_m, beam.getM());
    }

    if (distrTypeT_m == BINOMIAL || distrTypeT_m == GAUSS) {

        corr_m[0] = Attributes::getReal(itsAttr[CORRX]);
        corr_m[1] = Attributes::getReal(itsAttr[CORRY]);
        corr_m[2] = Attributes::getReal(itsAttr[CORRT]);
        corr_m[3] = Attributes::getReal(itsAttr[R61]);
        corr_m[4] = Attributes::getReal(itsAttr[R62]);
        corr_m[5] = Attributes::getReal(itsAttr[R51]);
        corr_m[6] = Attributes::getReal(itsAttr[R52]);
        gauss_offset_m[0] = Attributes::getReal(itsAttr[OFFSETX]);
        gauss_offset_m[1] = Attributes::getReal(itsAttr[OFFSETY]);

        sigx_m = Vector_t(Attributes::getReal(itsAttr[SIGMAX]),
                          Attributes::getReal(itsAttr[SIGMAY]),
                          Attributes::getReal(itsAttr[SIGMAT]));
        
        sigp_m = Vector_t(eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]), beam.getM()),
                          eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]), beam.getM()),
                          eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]), beam.getM()));

        binc_m = Vector_t(Attributes::getReal(itsAttr[MX]),
                          Attributes::getReal(itsAttr[MY]),
                          Attributes::getReal(itsAttr[MT]));
    }

    switch(distrTypeT_m) {
        case SURFACERANDCREATE: {
            darkCurrentParts_m = (size_t) Attributes::getReal(itsAttr[NPDARKCUR]);
            darkInwardMargin_m = Attributes::getReal(itsAttr[INWARDMARGIN]);
            //ppVw_m = Attributes::getReal(itsAttr[VW]);
            //vVThermal_m = Attributes::getReal(itsAttr[VVTHERMAL]);
        }
        break;
        case SURFACEEMISSION: {
            darkCurrentParts_m = (size_t) Attributes::getReal(itsAttr[NPDARKCUR]);
            darkInwardMargin_m = Attributes::getReal(itsAttr[INWARDMARGIN]);
            eInitThreshold_m = Attributes::getReal(itsAttr[EINITHR]);

            workFunction_m = Attributes::getReal(itsAttr[FNPHIW]);
            fieldEnhancement_m = Attributes::getReal(itsAttr[FNBETA]);
            maxFN_m = (size_t) Attributes::getReal(itsAttr[FNMAXEMI]);
            fieldThrFN_m = Attributes::getReal(itsAttr[FNFIELDTHR]);
            paraFNA_m = Attributes::getReal(itsAttr[FNA]);
            paraFNB_m = Attributes::getReal(itsAttr[FNB]);
            paraFNY_m = Attributes::getReal(itsAttr[FNY]);
            paraFNVYZe_m = Attributes::getReal(itsAttr[FNVYZERO]);
            paraFNVYSe_m = Attributes::getReal(itsAttr[FNVYSECOND]);

            secondaryFlag_m = Attributes::getReal(itsAttr[SECONDARYFLAG]);

        }
        break;
        case ASTRAFLATTOPTH: {

            if(Options::rngtype != string("RANDOM")) {
                INFOMSG("RNGTYPE= " << Options::rngtype << endl);
                if(Options::rngtype == string("HALTON"))
                    qrng_m = gsl_qrng_alloc(gsl_qrng_halton, 2);
                else if(Options::rngtype == string("SOBOL"))
                    qrng_m = gsl_qrng_alloc(gsl_qrng_sobol, 2);
                else if(Options::rngtype == string("NIEDERREITER"))
                    qrng_m = gsl_qrng_alloc(gsl_qrng_niederreiter_2, 2);
                else {
                    INFOMSG("RNGTYPE= " << Options::rngtype << " not known, using HALTON" << endl);
                    qrng_m = gsl_qrng_alloc(gsl_qrng_halton, 2);
                }
            }

            rGen_m = new RANLIB_class(265314159, 4);

            gsl_rng_env_setup();
            rn_m = gsl_rng_alloc(gsl_rng_1dhalton);
            R_m = gsl_qrng_alloc(gsl_qrng_halton, 2);
            int binTotal = sBins_m * nBins_m;

            h_m = gsl_histogram_alloc(binTotal);
            distributionTable_m = new double[binTotal + 1];

            double a = tPulseLengthFWHM_m / 2.;
            double sig = tRise_m / 2.;
            double inv_erf08 = 0.906193802436823; // erfinv(0.8)
            double sqr2 = sqrt(2.);
            double t = a - sqr2 * sig * inv_erf08;
            double tmps = sig;
            double tmpt = t;
            for(int i = 0; i < 10; ++ i) {
                sig = (t + tRise_m - a) / (sqr2 * inv_erf08);
                t = a - 0.5 * sqr2 * (sig + tmps) * inv_erf08;
                sig = (0.5 * (t + tmpt) + tRise_m - a) / (sqr2 * inv_erf08);
                tmps = sig;
                tmpt = t;
            }
            tEmission_m = tPulseLengthFWHM_m + 10. * sig;
            tBin_m = tEmission_m / nBins_m;

            double lo = -tBin_m / 2.0 * nBins_m;
            double hi = tBin_m / 2.0 * nBins_m;
            double dx = tBin_m / sBins_m;
            double x = lo;
            double tot = 0;
            double weight = 2.0;
            gsl_histogram_set_ranges_uniform(h_m, lo, hi);

            // sample the function that describes the histogram of the requested distribution
            for(int i = 0; i < binTotal + 1; ++ i, x += dx, weight = 6. - weight) {
                distributionTable_m[i] = gsl_sf_erf((x + a) / (sqrt(2) * sig)) - gsl_sf_erf((x - a) / (sqrt(2) * sig));
                tot += distributionTable_m[i] * weight;
            }
            tot -= distributionTable_m[binTotal] * (5. - weight);
            tot -= distributionTable_m[0];

            for(int k = 0; k < nBins_m; ++ k) {
                gsl_ran_discrete_t *table = gsl_ran_discrete_preproc(sBins_m, &(distributionTable_m[sBins_m * k]));
                double loc_fraction = -distributionTable_m[sBins_m * k] / tot;

                weight = 2.0;
                for(int i = sBins_m * k; i < sBins_m * (k + 1) + 1; ++ i, weight = 6. - weight) {
                    loc_fraction += distributionTable_m[i] * weight / tot;
                }
                loc_fraction -= distributionTable_m[sBins_m * (k + 1)] * (5. - weight) / tot;
                int bin_size = static_cast<int>(floor(loc_fraction * Np + 0.5)); //number of particles in bin!

                for(int i = 0; i < bin_size; i++) {
                    double xx[2];
                    gsl_qrng_get(R_m, xx);
                    gsl_histogram_increment(h_m, (hi * (xx[1] + 
                    static_cast<int>(gsl_ran_discrete(rn_m, table)) - binTotal / 2 + k * sBins_m) / (binTotal / 2)));
                }
                gsl_ran_discrete_free(table);
            }
            pbin_m->setHistogram(h_m);

            // ASTRA mode
            phimax_m = Physics::pi / 2.0;
            *gmsg << " -- B I N N I N G in T -----------------------------------------" << endl;
            *gmsg << " ---------------------I N P U T --------------------------------" << endl;
            *gmsg << " ASTRA FLAT TOP &  THERMAL EMITTANCE in ASTRA MODE" << endl;
            *gmsg << " Kinetic energy (thermal emittance) = " << ekin_m << " [eV]  " << endl;
            *gmsg << " Phi max = " << phimax_m * 180 / Physics::pi << " [deg]  " << endl;
            *gmsg << " tBin = " << tBin_m << " [sec]  nBins = " << nBins_m << " tEmission =  " << tEmission_m << " [sec] " << endl;

            writeToFile();
        }

        break;
        case GUNGAUSSFLATTOPTH: {

            tEmission_m = tPulseLengthFWHM_m + (cutoff_m - sqrt(2.0 * log(2.0))) * (sigmaRise_m + sigmaFall_m);
            tBin_m = tEmission_m / nBins_m;

            rGen_m = new RANLIB_class(265314159, 4);

            //FIXME: hack
            gsl_rng_env_setup();
            rn_m = gsl_rng_alloc(gsl_rng_1dhalton);
            R_m = gsl_qrng_alloc(gsl_qrng_halton, 2);
            int binTotal = sBins_m * nBins_m; // number of sampling bins

            h_m = gsl_histogram_alloc(binTotal);
            createTimeBins(Np);
            pbin_m->setHistogram(h_m);

            distributionTable_m = new double[binTotal];
            for(int i = 0; i < binTotal; i++)
                distributionTable_m[i] = gsl_histogram_get(h_m, i);

            // ASTRA mode
            phimax_m = Physics::pi / 2.0;
            *gmsg << " -- B I N N I N G in T -----------------------------------------" << endl;
            *gmsg << " ---------------------I N P U T --------------------------------" << endl;
            *gmsg << " GUNGAUSS FLAT TOP &  THERMAL EMITTANCE" << endl;
            *gmsg << " Kinetic energy (thermal emittance) = " << ekin_m << " [eV]  " << endl;
            *gmsg << " Phi max = " << phimax_m * 180 / Physics::pi << " [deg]  " << endl;
            *gmsg << " tBin = " << tBin_m << " [sec]  nBins = " << nBins_m << " tEmission =  " << tEmission_m << " [sec] " << endl;

            writeToFile();
        }
        break;

        case BINOMIAL: {

            for(int j = 0; j < 3; j++) {
                double chi = asin(corr_m[j]);
                emit_m[j] = sigx_m[j] * sigp_m[j] * cos(chi);
            }

            for(int j = 0; j < 3; j++) {
                beta_m[j]  = sigx_m[j] * sigx_m[j] / emit_m[j];
                gamma_m[j] = sigp_m[j] * sigp_m[j] / emit_m[j];
                alpha_m[j] = -corr_m[j] * sqrt(beta_m[j] * abs(gamma_m[j]));
            }
            createBinom(emit_m, alpha_m, beta_m, gamma_m, binc_m, beam, Np, isBinned);
        }
        break;

        case GAUSS: {

            avrgpt_m = eVtoBetaGamma(Attributes::getReal(itsAttr[PT]), beam.getM());
            avrgt_m  = Attributes::getReal(itsAttr[T]);

            /*
              give up the portability w.r.t. the rangen
              and hope to be more scalable
            */
            rGen_m = new RANLIB_class((Ippl::myNode() + 1) * 265314159, 4);

            IpplTimings::startTimer(beam.distrCreate_m);

            if(Np > 1E8) {
                int k = 10;
                Np = (size_t)Np / k;
                //Np = (size_t)Np/Ippl::getNodes()/k;
                *gmsg << "Sampl= " << Np *Ippl::getNodes() << " x " << k << " Total= " << k *Np *Ippl::getNodes() <<  endl;
                for(int kk = 0; kk < k; kk++) {
                    sampleGauss(beam, kk * Np);
                    beam.boundp();
                    *gmsg << "Sampl Gauss k= " << kk << " N= " << beam.getTotalNum() << endl;
                }
            } else {
                //Np = (size_t) Np / Ippl::getNodes();
                sampleGauss(beam, Np);
            }
            *gmsg << "Sample Gauss done ..." << endl;

            IpplTimings::stopTimer(beam.distrCreate_m);
        }
        break;
        case FROMFILE: {
            *gmsg << "\n-------------------------------------------------------------" << endl;
            *gmsg << "     READ ININITAL DISTRIBUTION FROM FILE    " << endl;
            *gmsg << "     BE AWARE OF THE FACT THAT ONLY NODE 0 IS READING IN " << endl;
            *gmsg << "-------------------------------------------------------------\n" << endl;

            if(isBinned) {
                *gmsg << "     DISTRIBUTION will be binned using " << nBins_m << " energy bins " << endl;
                const string fn;
                binnDistributionFromFile(beam, fn);

            } else {
                std::ofstream os;
                if(Ippl::getNodes() == 1) {
                    *gmsg << " Write distribution to file dist.dat" << endl;
                    string file("data/dist.dat");
                    os.open(file.c_str());
                    if(os.bad()) {
                        *gmsg << "Unable to open output file " <<  file << endl;
                    }
                    os << "# x px y py z pz " << endl;
                }

                if(Ippl::myNode() == 0) {
                    const string filename = Attributes::getString(itsAttr[FNAME]);
                    double x0, px0, y0, py0, psi0, del0;

                    std::ifstream fs;
                    fs.open(filename.c_str());

                    if(fs.fail()) {
                        throw OpalException("Distribution::Create()",
                                            "Open file operation failed, please check if \""
                                            + filename +  "\" really exists.");
                    }

                    fs >> Np;
                    if(Np <= 0) {
                        throw OpalException("Distribution::Create()",
                                            " The particle number should be bigger than zero! Please check the first line of file \""
                                            + filename +  "\".");
                    }

                    for(unsigned int i = 0; i < Np; i++) {
                        if(!fs.eof()) {
                            beam.create(1);
                            fs >> x0 >> px0 >> y0 >> py0 >> psi0 >> del0;
                            beam.R[i] = Vector_t(x0, y0, psi0);
                            beam.P[i] = Vector_t(px0, py0, del0);
                            beam.Bin[i] = 0; // not initialized
                            beam.Q[i] = beam.getChargePerParticle();
                            beam.PType[i] = 0;
                            if(Ippl::getNodes() == 1) {
                                os <<  beam.R[i](0) << "\t " <<  beam.P[i](0)    << "\t "
                                   <<  beam.R[i](1) << "\t " <<  beam.P[i](1)    << "\t "
                                   <<  beam.R[i](2) << "\t " <<  beam.P[i](2)    << "\t "
                                   << endl;
                            }
                        } else {
                            throw OpalException("Distribution::Create()",
                                                "End of file reached before all particles imported, please check file \""
                                                + filename +  "\".");
                            return;
                        }
                    }
                    fs.close();
                    os.close();
                }
            }
        }
        break;
        default:
            INFOMSG("Distribution unknown" << endl;);
    }

    /*
      In the case of a binned distribution (gun)
      we have to do the boundp after emission.
    */

    if(isBinned)
        beam.setPBins(pbin_m);
    else
        beam.boundp();
    beam.LastSection = 0;

}

/*
 * This adds other distributions to the main distribution. This is currently only
 * defined to work with distribution type "GUNGAUSSFLATTOPTH". All others will return
 * a false value.
 *
 * Also note that only the time structure of the added distributions are used. The
 * transverse profile is defined only by the main distribution (that is being added to).
 *
 */
bool Distribution::addDistributions(PartBunch &beam, vector<Distribution *> distributions, size_t numberOfParticles) {
    /// Check main distribution type.
    switch(distrTypeT_m) {

        case GUNGAUSSFLATTOPTH: {

            /// Find weight of each distribution. Also check that we have distributions to add.
            double totalWeight = Attributes::getReal(itsAttr[WEIGHT]);
            unsigned int numberOfDistToAdd = 0;
            vector<double> relativeWeight;
            relativeWeight.push_back(totalWeight);
            for(vector<Distribution *>::const_iterator distIterator = distributions.begin(); distIterator != distributions.end(); distIterator++) {

                if(distT_m != Attributes::getString(*((*distIterator)->findAttribute("DISTRIBUTION")))) {
                    *gmsg << " --- Mismatched Distribution types in Distribution List ---" << endl
                          << " Type: " << Attributes::getString(*((*distIterator)->findAttribute("DISTRIBUTION"))) << endl
                          << " Distribution will not be used. " << endl
                          << " ----------------------------------------------------------" << endl;
                    relativeWeight.push_back(0.0);
                } else {
                    relativeWeight.push_back(Attributes::getReal(*((*distIterator)->findAttribute("WEIGHT"))));
                    totalWeight += Attributes::getReal(*((*distIterator)->findAttribute("WEIGHT")));
                    numberOfDistToAdd++;
                }
            }

            if(numberOfDistToAdd == 0)
                return false;
            else {

                for(vector<double>::iterator weightIterator = relativeWeight.begin(); weightIterator != relativeWeight.end(); weightIterator++)
                    *weightIterator /= totalWeight;

                /*
                 * Find emission bounds of each distribution, including delay defined by itsAttr[T]".
                 * The total emission time will be defined by these bounds.
                 *
                 * We assume that particles start to emit at t = 0 i the absence of a distribution time
                 * shift.
                 */
                double maxT = 0.0;
                double minT = 0.0;
                for(unsigned int distIterator = 0; distIterator < relativeWeight.size(); distIterator++) {

                    double deltaT = 0.0;
                    if(distIterator == 0) {

                        deltaT = Attributes::getReal(itsAttr[T]);
                        maxT = tEmission_m + deltaT;
                        minT = deltaT;

                        // FIXME: floating point comparison
                    } else if(relativeWeight.at(distIterator) != 0.0) {

                        // // Find emission time without time shift.
                        // double pulseLengthFWHM = Attributes::getReal(*(distributions.at(distIterator - 1)->findAttribute("TPULSEFWHM")));
                        // double cutOff = Attributes::getReal(*(distributions.at(distIterator - 1)->findAttribute("CUTOFF")));
                        // double riseTime = Attributes::getReal(*(distributions.at(distIterator - 1)->findAttribute("TRISE")));
                        // double fallTime = Attributes::getReal(*(distributions.at(distIterator - 1)->findAttribute("TFALL")));

                        // double timeRatio = sqrt(2.0 * log(10.0)) - sqrt(2.0 * log(10.0 / 9.0));
                        // double sigmaRiseTime = riseTime / timeRatio;
                        // double sigmaFallTime = fallTime / timeRatio;

                        // double emissionTime = pulseLengthFWHM + (cutOff - sqrt(2.0 * log(2.0))) * (sigmaRiseTime + sigmaFallTime);

                        // Find max. and min. time.
                        deltaT = Attributes::getReal(*(distributions.at(distIterator - 1)->findAttribute("T")));
                        if(maxT < tEmission_m + deltaT) maxT = tEmission_m + deltaT;
                        if(minT > deltaT) minT = deltaT;
                    }
                }
                tEmission_m = maxT - minT;
                tBin_m = tEmission_m / nBins_m;

                /*
                 * Now we reset the main histogram and reallocate particles to the different time bins.
                 */
                gsl_histogram_reset(h_m);
                gsl_histogram_set_ranges_uniform(h_m, 0, tEmission_m);

                gsl_rng_env_setup();
                gsl_rng *ranNumberGen = gsl_rng_alloc(gsl_rng_default);

                unsigned int particleCount = 0;
                double deltaT = 0.0;

                for(int distIterator = relativeWeight.size() - 1; distIterator >= 0; distIterator--) {

                    // Calculate number of particles for this distribution.
                    unsigned int nParticles = 0;
                    if(distIterator != 0) {
                        nParticles = static_cast<unsigned int>(numberOfParticles * relativeWeight.at(distIterator));
                        particleCount += nParticles;
                    } else
                        nParticles = numberOfParticles - particleCount;

                    // Add particles to time histogram.
                    if(nParticles > 0) {
                        
                        double sigmaRiseTime = 0.0;
                        double sigmaFallTime = 0.0;
                        double timeFlat = 0.0;
                        double cutOff = 0.0;
                        
                        if(distIterator == 0) {
                            sigmaRiseTime = sigmaRise_m;
                            sigmaFallTime = sigmaFall_m;
                            timeFlat = tPulseLengthFWHM_m - sqrt(2.0 * log(2.0)) * (sigmaRise_m + sigmaFall_m);
                            cutOff = cutoff_m;
                                deltaT = Attributes::getReal(itsAttr[T]);
                        } else {
                            double riseTime = Attributes::getReal(*(distributions.at(distIterator - 1)->findAttribute("TRISE")));
                            double fallTime = Attributes::getReal(*(distributions.at(distIterator - 1)->findAttribute("TFALL")));
                            
                            cutOff = Attributes::getReal(*(distributions.at(distIterator - 1)->findAttribute("CUTOFF")));
                            deltaT = Attributes::getReal(*(distributions.at(distIterator - 1)->findAttribute("T")));
                            
                            double timeRatio = sqrt(2.0 * log(10.0)) - sqrt(2.0 * log(10.0 / 9.0));
                            sigmaRiseTime = riseTime / timeRatio;
                            sigmaFallTime = fallTime / timeRatio;
                            
                            timeFlat = Attributes::getReal(*(distributions.at(distIterator - 1)->findAttribute("TPULSEFWHM")))
                                - sqrt(2.0 * log(2.0)) * (sigmaRiseTime + sigmaFallTime);
                        }
                        if(timeFlat < 0.0) timeFlat = 0.0;
                        
                        const double totalArea = timeFlat + 0.5 * sqrt(2.0 * Physics::pi) * (sigmaRiseTime + sigmaFallTime);
                        
                        unsigned int numPartInRise = nParticles * 0.5 * gsl_sf_erf(cutOff / sqrt(2.0))
                            * sqrt(2.0 * Physics::pi) * sigmaRiseTime / totalArea;
                        unsigned int numPartInFall = nParticles * 0.5 * gsl_sf_erf(cutOff / sqrt(2.0))
                            * sqrt(2.0 * Physics::pi) * sigmaFallTime / totalArea;
                        unsigned int numPartInFlat = nParticles - numPartInRise - numPartInFall;
                        
                        if(timeFlat == 0.0) {
                            numPartInRise += numPartInFlat / 2;
                            numPartInFall = nParticles - numPartInRise;
                            numPartInFlat = 0;
                        }
                        
                        for(unsigned int partIterator = 0; partIterator < numPartInRise; partIterator++) {
                            double tRandom = gsl_ran_gaussian_tail(ranNumberGen, 0, sigmaRiseTime);
                            while(tRandom > cutOff * sigmaRiseTime)
                                tRandom = gsl_ran_gaussian_tail(ranNumberGen, 0, sigmaRiseTime);
                            gsl_histogram_increment(h_m, -tRandom + cutOff * sigmaRiseTime + deltaT - minT);
                        }
                        
                        for(unsigned int partIterator = 0; partIterator < numPartInFall; partIterator++) {
                            double tRandom = gsl_ran_gaussian_tail(ranNumberGen, 0, sigmaFallTime);
                            while(tRandom > cutOff * sigmaFallTime)
                                tRandom = gsl_ran_gaussian_tail(ranNumberGen, 0, sigmaFallTime);
                            gsl_histogram_increment(h_m, tRandom + cutOff * sigmaRiseTime + timeFlat + deltaT - minT);
                        }
                        
                        for(unsigned int partIterator = 0; partIterator < numPartInFlat; partIterator++) {
                            double tRandom = 0.0;
                            gsl_qrng_get(R_m, &tRandom);
                            tRandom *= timeFlat;
                            gsl_histogram_increment(h_m, tRandom + cutOff * sigmaRiseTime + deltaT - minT);
                        }
                    }
                }

                // Free random number generator.
                gsl_rng_free(ranNumberGen);
            }

            /*
             * Do a boundp if we don't have energy bins.
             */
            if(!(nBins_m > 0))
                beam.boundp();

            /*
             * Write time histogram to file.
             */
            if(Ippl::myNode() == 0) {
                FILE *fp;
                fp = fopen("data/hist.dat", "w");
                gsl_histogram_fprintf(fp, h_m, "%g", "%g");
                fclose(fp);
            }

            return true;
        }
        break;

        default:
            return false;
    }
}

/**
 * This is the generator for a Gaussian distribution
 *
 * @param beam
 * @param Np
 */
void Distribution::sampleGauss(PartBunch &beam, size_t Np) {
    int pc = 0;
    size_t count = 0;

    for(size_t i = beam.getTotalNum(); i < Np; i++) {
        double x, y;      // generate independent Gaussians, then correlate and finaly scale
        x  = rGen_m->gauss(0.0, 1.0);
        y  = rGen_m->gauss(0.0, 1.0);
        double xx = x;
        double yy = y;
        double px0  = x * corr_m[0] + y * sqrt(1.0 - corr_m[0] * corr_m[0]);
        double x0   = x * sigx_m[0] + gauss_offset_m[0];
        px0 *= sigp_m[0]; 

        x  = rGen_m->gauss(0.0, 1.0);
        y  = rGen_m->gauss(0.0, 1.0);
        double py0  = x * corr_m[1] + y * sqrt(1.0 - corr_m[1] * corr_m[1]);
        double y0   =  x * sigx_m[1] + gauss_offset_m[1];
        py0 *= sigp_m[1];

        double del0;
        double psi0;
        x  = rGen_m->gauss(0.0, 1.0);
        y  = rGen_m->gauss(0.0, 1.0);

        double tempa = 1.0 - corr_m[0] * corr_m[0];
        const double l32 = (corr_m[6] - corr_m[0] * corr_m[5]) / sqrt(abs(tempa)) * tempa / abs(tempa);
        double tempb = 1 - corr_m[5] * corr_m[5] - l32 * l32;
        const double l33 = sqrt(abs(tempb)) * tempb / abs(tempb);
        psi0 = xx * corr_m[5] + yy * l32 + x * l33;
        const double l42 = (corr_m[4] - corr_m[0] * corr_m[3]) / sqrt(abs(tempa)) * tempa / abs(tempa);
        const double l43 = (corr_m[2] - corr_m[5] * corr_m[3] - l42 * l32) / l33;
        double tempc = 1 - corr_m[3] * corr_m[3] - l42 * l42 - l43 * l43;
        const double l44 = sqrt(abs(tempc)) * tempc / abs(tempc);

        del0 = xx * corr_m[3] + yy * l42 + x * l43 + y * l44;
        psi0 = avrgt_m + psi0 * sigx_m[2];
        del0 = avrgpt_m + sigp_m[2] * del0;
        if(pc == Ippl::myNode()) {
            beam.create(1);
            beam.R[count] = Vector_t(x0, y0, psi0);
            beam.P[count] = Vector_t(px0, py0, del0);
            beam.Bin[count] = 0; // not initialized
            beam.Q[count] = beam.getChargePerParticle();
            beam.PType[count] = 0;
            beam.TriID[count] = 0;
            count++;
        }
        pc++;
        if(pc == Ippl::getNodes())
            pc = 0;
    }
}

/**
 *
 *
 * @param dt
 *
 * @return
 */
pair<Vector_t, Vector_t> Distribution::sample(double dt, int binNumber) {

    Vector_t r(0.0);
    Vector_t p(0.0);

    double x0, y0;

    switch(distrTypeT_m) {

        case ASTRAFLATTOPTH:
        case GUNGAUSSFLATTOPTH: {
            double x, y;
            double xy = 6;

            if(lp_m != NULL) {
                lp_m->GetXY(&x, &y);
                x = 2 * x - 1.0;
                y = 2 * y - 1.0;
            } else if(qrng_m != NULL) {
                while(xy > 1) {
                    double v0[2];
                    gsl_qrng_get(qrng_m, v0);
                    x = -1.0 + (2.0 * v0[0]);
                    y = -1.0 + (2.0 * v0[1]);
                    xy = sqrt(x * x + y * y);
                }
            } else {
                while(xy > 1) {
                    x  = rGen_m->uniform(-1.0, 1.0);
                    y  = rGen_m->uniform(-1.0, 1.0);
                    xy = sqrt(x * x + y * y);
                }
            }

            x0   =  x * sigx_m[0];
            y0   =  y * sigx_m[1];

            /*
              Now calculate the thermal emittances
            */

            const double phi   = 2.0 * acos(sqrt(rGen_m->uniform(0.0, 1.0)));
            const double theta = 2.0 * Physics::pi * rGen_m->uniform(0.0, 1.0);

            const double px0   = ptot_m * sin(phi) * cos(theta);
            const double py0   = ptot_m * sin(phi) * sin(theta);
            const double del0  = ptot_m * abs(cos(phi));

            p = Vector_t(px0, py0, del0);


            gsl_ran_discrete_t *table = gsl_ran_discrete_preproc(sBins_m, &(distributionTable_m[sBins_m * binNumber]));
            double xr[2] = {0.0, 0.0};
            gsl_qrng_get(R_m, xr);
            double s0 = dt * (xr[1] + static_cast<int>(gsl_ran_discrete(rn_m, table))) / sBins_m;
            r = Vector_t(x0, y0, s0);

            gsl_ran_discrete_free(table);
            break;
        }

        default:
            INFOMSG("Distribution unknown" << endl;);
    }
    return pair<Vector_t, Vector_t>(r, p);
}

pair<Vector_t, Vector_t> Distribution::sampleNEW(double dt, int binNumber) {

    Vector_t r(0.0);
    Vector_t p(0.0);

    double x0, y0;

    switch(distrTypeT_m) {

        case ASTRAFLATTOPTH:
        case GUNGAUSSFLATTOPTH: {
            double x, y;
            double xy = 6;

            if(lp_m != NULL) {
                lp_m->GetXY(&x, &y);
                x = 2 * x - 1.0;
                y = 2 * y - 1.0;
            } else if(qrng_m != NULL) {
                while(xy > 1) {
                    double v0[2];
                    gsl_qrng_get(qrng_m, v0);
                    x = -1.0 + (2.0 * v0[0]);
                    y = -1.0 + (2.0 * v0[1]);
                    xy = sqrt(x * x + y * y);
                }
            } else {
                while(xy > 1) {
                    x  = rGen_m->uniform(-1.0, 1.0);
                    y  = rGen_m->uniform(-1.0, 1.0);
                    xy = sqrt(x * x + y * y);
                }
            }

            x0   =  x * sigx_m[0];
            y0   =  y * sigx_m[1];

            /*
              Now calculate the thermal emittances
            */

            const double phi   = 2.0 * acos(sqrt(rGen_m->uniform(0.0, 1.0)));
            const double theta = 2.0 * Physics::pi * rGen_m->uniform(0.0, 1.0);

            const double px0   = ptot_m * sin(phi) * cos(theta);
            const double py0   = ptot_m * sin(phi) * sin(theta);
            const double del0  = ptot_m * abs(cos(phi));

            p = Vector_t(px0, py0, del0);

            //           gsl_ran_discrete_t *table = gsl_ran_discrete_preproc(1, &(distributionTable_m[binNumber]));
            double xr;
            gsl_qrng_get(R_m, &xr);

            double s0 = dt * xr;

            r = Vector_t(x0, y0, s0);

            //            gsl_ran_discrete_free(table);
            break;
        }

        default:
            INFOMSG("Distribution unknown" << endl;);
    }
    return pair<Vector_t, Vector_t>(r, p);
}

/**
 * This method fills the gsl histogram (h_m) with a binned
 * Gauss-Flattop-Distribution as specified in the manual.
 *
 * @param Np  number of total particles to generate
 */
void Distribution::createTimeBins(const int Np) {

    if(Options::rngtype != string("RANDOM")) {
        INFOMSG("RNGTYPE= " << Options::rngtype << endl);
        if(Options::rngtype == string("HALTON"))
            qrng_m = gsl_qrng_alloc(gsl_qrng_halton, 2);
        else if(Options::rngtype == string("SOBOL"))
            qrng_m = gsl_qrng_alloc(gsl_qrng_sobol, 2);
        else if(Options::rngtype == string("NIEDERREITER"))
            qrng_m = gsl_qrng_alloc(gsl_qrng_niederreiter_2, 2);
        else {
            INFOMSG("RNGTYPE= " << Options::rngtype << " not known, using HALTON" << endl);
            qrng_m = gsl_qrng_alloc(gsl_qrng_halton, 2);
        }
    }
    
    gsl_histogram_set_ranges_uniform(h_m, 0, tEmission_m);
    gsl_rng_env_setup();
    gsl_rng *r = gsl_rng_alloc(gsl_rng_default);
    const double sq2pi = sqrt(2.0 * Physics::pi);
    double tFlat = tPulseLengthFWHM_m - sqrt(2.0 * log(2.0)) * (sigmaRise_m + sigmaFall_m);
    if(tFlat < 0.0) tFlat = 0.0;
    const double totA = tFlat + 0.5 * sq2pi * (sigmaRise_m + sigmaFall_m);
    int nrise = Np * 0.5 * gsl_sf_erf(cutoff_m / sqrt(2.0)) * sq2pi * sigmaRise_m / totA;
    int nfall = Np * 0.5 * gsl_sf_erf(cutoff_m / sqrt(2.0)) * sq2pi * sigmaFall_m / totA;
    int nflat = Np - nrise - nfall;
    
    if(tFlat == 0.0) {
        nrise += nflat / 2;
        nfall = Np - nrise;
        nflat = 0;
    }
    
    // Rise: [0, c\sigma_R]
    for(int i = 0; i < nrise; i++) {
        double r1 = gsl_ran_gaussian_tail(r, 0, sigmaRise_m);
        while(r1 > cutoff_m * sigmaRise_m)
            r1 = gsl_ran_gaussian_tail(r, 0, sigmaRise_m);
        gsl_histogram_increment(h_m, -r1 + cutoff_m * sigmaRise_m);
    }
    // Fall: [c\sigma_R + tFlat, c\sigma_R + tFlat + c\sigma_F]
    for(int i = 0; i < nfall; i++) {
        double r1 = gsl_ran_gaussian_tail(r, 0, sigmaFall_m);
        while(r1 > cutoff_m * sigmaFall_m)
            r1 = gsl_ran_gaussian_tail(r, 0, sigmaFall_m);
        gsl_histogram_increment(h_m, r1 + cutoff_m * sigmaRise_m + tFlat);
    }
    // Flattop: [c\sigma_R, c\sigma_R + tFlat]
    //
    // The flat top can also have sinusoidal modulations.
    
    gsl_qrng *Qrng = gsl_qrng_alloc(gsl_qrng_halton, 1);
    
    gsl_qrng *R2_m = gsl_qrng_alloc(gsl_qrng_halton, 2);
    
    // Get modulation parameters.
    double modulationAmplitude = Attributes::getReal(itsAttr[FTOSCAMPLITUDE]) / 100.0;
    double numberOfModulationPeriods = fabs(Attributes::getReal(itsAttr[FTOSCPERIODS]));
    double modulationPeriod = 0.0;
    if(numberOfModulationPeriods != 0) modulationPeriod = tFlat / numberOfModulationPeriods;
    
    // Sample flat top.
    for(int i = 0; i < nflat; i++) {
        double r1 = 0.0;
        double r2 = 0.0;
        double rn[2] = {0.0, 0.0};
        if(modulationAmplitude == 0.0 || numberOfModulationPeriods == 0) {
            // r1 = gsl_ran_flat(r, 0, tFlat);
            gsl_qrng_get(Qrng, &r1);
            
            r1 *= tFlat;
            
        } else {
            bool accept = false;
            while(!accept) {
                gsl_qrng_get(R2_m, rn);
                // r1 = gsl_ran_flat(r, 0, tFlat);
                /// r2 = gsl_ran_flat(r, 0, 1.0);
                    r1 = rn[0] * tFlat;
                    r2 = rn[1];
                    double function = (1.0 + modulationAmplitude * sin(Physics::two_pi * r1 / modulationPeriod))
                        / (1.0 + fabs(modulationAmplitude));
                    if(r2 <= function) accept = true;
            }
        }
        gsl_histogram_increment(h_m, r1 + cutoff_m * sigmaRise_m);
    }
    if(Ippl::myNode() == 0) {
        FILE *fp;
        fp = fopen("data/hist.dat", "w");
        gsl_histogram_fprintf(fp, h_m, "%g", "%g");
        fclose(fp);
    }
    gsl_qrng_free(Qrng);
    gsl_rng_free(r);
}

/**
 *
 *
 * @param p
 */
void Distribution::createSlicedBunch(int sl, double charge, double gamma, double mass, double current, double center, double Bz0, EnvelopeBunch *p) {
    double beamWidth = 0.0;
    double beamEnergy = 0.0;
    //int sl = (int) Attributes::getReal(itsAttr[NBIN]);
    *gmsg << "About to create a sliced bunch with " << sl << " slices" << endl;
    *gmsg << "mass = " << mass << " gamma = " << gamma << endl;
    IpplTimings::startTimer(p->distrCreate_m);

    distT_m = Attributes::getString(itsAttr[DISTRIBUTION]);
    if(distT_m == "GAUSS")
        distrTypeT_m = GAUSS;
    else if(distT_m == "GUNGAUSSFLATTOPTH")
        distrTypeT_m = GUNGAUSSFLATTOPTH;
    else if(distT_m == "FROMFILE")
        distrTypeT_m = FROMFILE;
    else if(distT_m == "BINOMIAL")
        distrTypeT_m = BINOMIAL;

    corr_m[0] = Attributes::getReal(itsAttr[CORRX]);
    corr_m[1] = Attributes::getReal(itsAttr[CORRY]);
    corr_m[2] = Attributes::getReal(itsAttr[CORRT]);
    corr_m[3] = Attributes::getReal(itsAttr[R61]);
    corr_m[4] = Attributes::getReal(itsAttr[R62]);
    corr_m[5] = Attributes::getReal(itsAttr[R51]);
    corr_m[6] = Attributes::getReal(itsAttr[R52]);
    
    sigx_m = Vector_t(Attributes::getReal(itsAttr[SIGMAX]),
                      Attributes::getReal(itsAttr[SIGMAY]),
                      Attributes::getReal(itsAttr[SIGMAT]));
    
    switch(distrTypeT_m) {

    case GUNGAUSSFLATTOPTH: {
        
        nBins_m = static_cast<int>(fabs(Attributes::getReal(itsAttr[NBIN])));
        
        tPulseLengthFWHM_m = Attributes::getReal(itsAttr[TPULSEFWHM]);
        cutoff_m = Attributes::getReal(itsAttr[CUTOFF]);
        tRise_m = Attributes::getReal(itsAttr[TRISE]);
        tFall_m = Attributes::getReal(itsAttr[TFALL]);
        double tratio = sqrt(2.0 * log(10.0)) - sqrt(2.0 * log(10.0 / 9.0));
        sigmaRise_m = tRise_m / tratio;
        sigmaFall_m = tFall_m / tratio;
        tEmission_m = tPulseLengthFWHM_m + (cutoff_m - sqrt(2.0 * log(2.0))) * (sigmaRise_m + sigmaFall_m);

        ekin_m = Attributes::getReal(itsAttr[EKIN]);

        // EnvelopeTracker expects [eV]
        beamEnergy = ekin_m;
        beamWidth = tEmission_m * Physics::c * sqrt(1.0 - (1.0 / pow(gamma, 2)));
        
        break;
    }
    case GAUSS: {
        gauss_offset_m[0] = Attributes::getReal(itsAttr[OFFSETX]);
        gauss_offset_m[1] = Attributes::getReal(itsAttr[OFFSETY]);        
        avrgt_m  = Attributes::getReal(itsAttr[T]);

        beamEnergy = (gamma * mass - mass) * 1e9;         //FIXME: why 1e9??
        beamWidth = sigx_m[2];
        tEmission_m = 0.0;
        break;
    }
    default:
        ;
    }

    center = -1 * beamWidth / 2.0;
    *gmsg << "x = " << sigx_m[0] << " y = " << sigx_m[1] << endl;
    double frac = 0.9;

    // execute initialization command
    p->initialize(sl, charge, beamEnergy, beamWidth, tEmission_m, frac, current, center, sigx_m[0], sigx_m[1], 0, 0, Bz0, nBins_m);
    IpplTimings::stopTimer(p->distrCreate_m);
}

/**
 * restart and envelope tracker run
 *
 * @param beam
 * @param Np
 * @param restartStep
 */
void Distribution::doRestartEnvelope(EnvelopeBunch &beam, size_t Np, int restartStep) {
    h5_file_t *H5file;
    h5_int64_t rc;
    string fn;

    IpplTimings::startTimer(beam.distrReload_m);

    if(OpalData::getInstance()->hasRestartFile()) {
        fn = OpalData::getInstance()->getRestartFileName();
        *gmsg << "Restart from a specified file:" << fn << endl;

    } else {
        //        beam.setTEmission(Attributes::getReal(itsAttr[TEMISSION]));
        fn = OpalData::getInstance()->getInputBasename() + string(".h5");
    }

#ifdef PARALLEL_IO
    H5file = H5OpenFile(fn.c_str(), H5_O_RDONLY, Ippl::getComm());
#else
    H5file = H5PartOpenFile(fn.c_str(), H5_O_RDONLY, 0);
#endif

    if(!H5file) {
        ERRORMSG("could not open file '" << fn << "';  exiting!" << endl);
        exit(0);
    }

    if(restartStep == -1) {
        restartStep = H5GetNumSteps(H5file) - 1;
        OpalData::getInstance()->setRestartStep(restartStep);
    } else {
        if(restartStep != H5GetNumSteps(H5file) - 1 && !OpalData::getInstance()->hasRestartFile()) {
            ERRORMSG("can't append to the file '" << fn << "' exiting!" << endl);
            exit(0);
        }
    }

    rc = H5SetStep(H5file, restartStep);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    int N = (int)H5PartGetNumParticles(H5file);

    h5_int64_t totalSteps = H5GetNumSteps(H5file);
    *gmsg << "total number of slices = " << N << " total steps " << totalSteps << endl;

    beam.distributeSlices(N);
    beam.createBunch();
    long long starti = beam.mySliceStartOffset();
    long long endi = beam.mySliceEndOffset();

    rc = H5PartSetView(H5file, starti, endi);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    N = (int)H5PartGetNumParticles(H5file);
    assert(N >= 0 && (unsigned int) N != beam.numMySlices());

    double actualT;
    rc = H5ReadStepAttribFloat64(H5file, "TIME", &actualT);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    beam.setT(actualT);
    double dPhiGlobal;
    rc = H5ReadFileAttribFloat64(H5file, "dPhiGlobal", &dPhiGlobal);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    OpalData::getInstance()->setGlobalPhaseShift(dPhiGlobal);

    h5_int64_t ltstep;
    rc = H5ReadStepAttribInt64(H5file, "LocalTrackStep", &ltstep);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    beam.setLocalTrackStep((long long)ltstep);

    h5_int64_t gtstep;
    rc = H5ReadStepAttribInt64(H5file, "GlobalTrackStep", &gtstep);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    beam.setGlobalTrackStep((long long)gtstep);

    std::unique_ptr<char[]> varray(new char[(N)*sizeof(double)]);
    double *farray = reinterpret_cast<double *>(varray.get());
    h5_int64_t *larray = reinterpret_cast<h5_int64_t *>(varray.get());

    rc = H5PartReadDataFloat64(H5file, "x", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n) {
        beam.setX(n, farray[n]);
    }
    rc = H5PartReadDataFloat64(H5file, "y", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n)
        beam.setY(n, farray[n]);

    rc = H5PartReadDataFloat64(H5file, "z", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n)
        beam.setZ(n, farray[n]);

    rc = H5PartReadDataFloat64(H5file, "px", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n)
        beam.setPx(n, farray[n]);

    rc = H5PartReadDataFloat64(H5file, "py", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n)
        beam.setPy(n, farray[n]);

    rc = H5PartReadDataFloat64(H5file, "beta", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n)
        beam.setBeta(n, farray[n]);

    rc = H5PartReadDataFloat64(H5file, "X0", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n)
        beam.setX0(n, farray[n]);

    rc = H5PartReadDataFloat64(H5file, "pX0", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n)
        beam.setPx0(n, farray[n]);

    rc = H5PartReadDataFloat64(H5file, "Y0", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n)
        beam.setY0(n, farray[n]);

    rc = H5PartReadDataFloat64(H5file, "pY0", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n)
        beam.setPy0(n, farray[n]);

    rc = H5PartReadDataInt64(H5file, "lastsection", larray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n)
        beam.LastSection[n] = (short) larray[n];

    Ippl::Comm->barrier();
    rc = H5CloseFile(H5file);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    beam.setCharge(beam.getChargePerParticle());
    IpplTimings::stopTimer(beam.distrReload_m);
}

/**
 *
 *
 * @param object
 *
 * @return
 */
bool Distribution::canReplaceBy(Object *object) {
    // Can replace only by another Distribution.
    return dynamic_cast<Distribution *>(object) != 0;
}

/**
 *
 *
 * @param name
 *
 * @return
 */
Distribution *Distribution::clone(const string &name) {
    return new Distribution(name, this);
}

/**
 *
 *
 */
void Distribution::execute() {
    update();
}


/**
 *
 *
 * @param emit
 * @param alpha
 * @param beta
 * @param gamma
 * @param bincoef
 * @param beam
 * @param particles
 * @param isBinned
 */
void Distribution::createBinom(Vector_t emit, Vector_t alpha, Vector_t beta, Vector_t gamma,
                               Vector_t bincoef, PartBunch &beam, size_t particles,
                               bool isBinned) {
    const double two_pi = 2.0 * 4.0 * atan(1.0);

    Vector_t M;
    Vector_t PM;
    Vector_t COSCHI;
    Vector_t SINCHI;
    Vector_t CHI;
    Vector_t AMI;
    Vector_t L;
    Vector_t PL;

    size_t pc = 0;
    size_t count = 0;

    INFOMSG("BINCOEF= " << bincoef << endl);

    for(int i = 0; i < 3; i++) {
        gamma[i] *= 4.0;
        beta[i] *= 4.0;
        M[i]       =  sqrt(emit[i] * beta[i]);
        PM[i]      =  sqrt(emit[i] * gamma[i]);
        COSCHI[i]  =  sqrt(1.0 / (1.0 + alpha[i] * alpha[i]));
        SINCHI[i]  = -alpha[i] * COSCHI[i];
        CHI[i]     =  atan2(SINCHI[i], COSCHI[i]);
        AMI[i]     =  1.0 / bincoef[i];
        L[i]       =  sqrt((bincoef[i] + 1.0) / 2.0) * M[i];
        PL[i]      =  sqrt((bincoef[i] + 1.0) / 2.0) * PM[i];
    }

    Vector_t x;
    // now copy this over to the bunch
    Vector_t p;
    double betagamma_part, pos_part = 0.0;
    if(isBinned) {
        double ekin = Attributes::getReal(itsAttr[PT]);
        double mass = beam.getM();
        double gamma_part = 1. + ekin / mass;
        betagamma_part = sqrt(ekin * ekin / (mass * mass) + 2.*ekin / mass);
        pbin_m->setGamma(gamma_part);
        *gmsg << "* Gamma = " << gamma_part << "; Beta = " << betagamma_part / gamma_part  << endl;
    } else {
        betagamma_part = eVtoBetaGamma(Attributes::getReal(itsAttr[PT]), beam.getM());
        pos_part = Attributes::getReal(itsAttr[T]);
    }

    for(size_t n = 0; n < particles; ++n) {
        double S1, S2;
        double A, AL, U, V;
        //  for(int i = 0; i < 3; i++) {
        S1 = IpplRandom();
        S2 = IpplRandom();
        if(bincoef[0] <= 10000) {
            A = sqrt(1.0 - pow(S1, AMI[0]));
            AL = two_pi * S2;
            U = A * cos(AL);
            V = A * sin(AL);
            double Ucp = U;
            double Vcp = V;
            x[0] = L[0] * U;
            p[0] = PL[0] * (U * SINCHI[0] + V * COSCHI[0]);

            S1 = IpplRandom();
            S2 = IpplRandom();
            A = sqrt(1.0 - pow(S1, AMI[1]));
            AL = two_pi * S2;
            U = A * cos(AL);
            V = A * sin(AL);
            x[1] = L[1] * U;
            p[1] = PL[1] * (U * SINCHI[1] + V * COSCHI[1]);

            S1 = IpplRandom();
            S2 = IpplRandom();
            A = sqrt(1.0 - pow(S1, AMI[2]));
            AL = two_pi * S2;
            U = A * cos(AL);
            V = A * sin(AL);

            double tempa = 1.0 - corr_m[0] * corr_m[0];
            const double l32 = (corr_m[6] - corr_m[0] * corr_m[5]) / sqrt(abs(tempa)) * tempa / abs(tempa);
            double tempb = 1 - corr_m[5] * corr_m[5] - l32 * l32;
            const double l33 = sqrt(abs(tempb)) * tempb / abs(tempb);
            x[2] = Ucp * corr_m[5] + Vcp * l32 + U * l33;
            const double l42 = (corr_m[4] - corr_m[0] * corr_m[3]) / sqrt(abs(tempa)) * tempa / abs(tempa);
            const double l43 = (corr_m[2] - corr_m[5] * corr_m[3] - l42 * l32) / l33;
            double tempc = 1 - corr_m[3] * corr_m[3] - l42 * l42 - l43 * l43;
            const double l44 = sqrt(abs(tempc)) * tempc / abs(tempc);

            p[2] = Ucp * corr_m[3] + Vcp * l42 + U * l43 + V * l44;
            x[2]  *= L[2];
            p[2] *= PL[2];

        } else {
            A = sqrt(2.0) / 2.0 * sqrt(-log(S1));
            AL = two_pi * S2;
            U = A * cos(AL);
            V = A * sin(AL);
            double Ucp = U;
            double Vcp = V;
            x[0] = M[0] * U;
            p[0] = PM[0] * (U * SINCHI[0] + V * COSCHI[0]);


            S1 = IpplRandom();
            S2 = IpplRandom();
            A = sqrt(2.0) / 2.0 * sqrt(-log(S1));
            AL = two_pi * S2;
            U = A * cos(AL);
            V = A * sin(AL);
            x[1] = M[1] * U;
            p[1] = PM[1] * (U * SINCHI[1] + V * COSCHI[1]);

            S1 = IpplRandom();
            S2 = IpplRandom();
            A = sqrt(2.0) / 2.0 * sqrt(-log(S1));
            AL = two_pi * S2;
            U = A * cos(AL);
            V = A * sin(AL);

            double tempa = 1.0 - corr_m[0] * corr_m[0];
            const double l32 = (corr_m[6] - corr_m[0] * corr_m[5]) / sqrt(abs(tempa)) * tempa / abs(tempa);
            double tempb = 1 - corr_m[5] * corr_m[5] - l32 * l32;
            const double l33 = sqrt(abs(tempb)) * tempb / abs(tempb);
            x[2] = Ucp * corr_m[5] + Vcp * l32 + U * l33;
            const double l42 = (corr_m[4] - corr_m[0] * corr_m[3]) / sqrt(abs(tempa)) * tempa / abs(tempa);
            const double l43 = (corr_m[2] - corr_m[5] * corr_m[3] - l42 * l32) / l33;
            double tempc = 1 - corr_m[3] * corr_m[3] - l42 * l42 - l43 * l43;
            const double l44 = sqrt(abs(tempc)) * tempc / abs(tempc);

            p[2] = Ucp * corr_m[3] + Vcp * l42 + U * l43 + V * l44;
            x[2]  *= M[2];
            p[2] *= PM[2];

        }
        p[2] += betagamma_part;
        x[2] += pos_part;


        //       }
        if(pc == (size_t) Ippl::myNode()) {
            if(isBinned) {
                vector<double> tmp;
                tmp.push_back(x[0]);
                tmp.push_back(x[1]);
                tmp.push_back(x[2]);
                tmp.push_back(p[0]);
                tmp.push_back(p[1]);
                tmp.push_back(p[2]);
                tmp.push_back(0);
                pbin_m->fill(tmp);
            } else {
                beam.create(1);
                beam.R[count] = x;
                beam.P[count] = p;
                beam.Q[count] = beam.getChargePerParticle();
                beam.PType[count] = 0;
                beam.TriID[count] = 0;
                count++;
            }
        }

        if(!isBinned) {
            pc++;
            if(pc == (size_t) Ippl::getNodes())
                pc = 0;
        }
    }
    if(isBinned) {
        pbin_m->sortArray();
        // now copy this over to the bunch
        // so that we can emmit the particles
        beam.setPBins(pbin_m);
    }
}

/**
 *
 *
 * @param p   particle bunch
 * @param bg  boundary geometry
 */

void  Distribution::create(PartBunch &beam, BoundaryGeometry &bg) {

    int pc = 0;
    size_t N_mean = static_cast<size_t>(floor(bg.getN() / Ippl::getNodes()));
    size_t N_extra = static_cast<size_t>(bg.getN() - N_mean * Ippl::getNodes());
    if(Ippl::myNode() == 0)
        N_mean += N_extra;
    size_t count = 0;
    size_t lowMark = beam.getLocalNum();
    if(bg.getN() != 0) {

        for(size_t i = 0; i < bg.getN(); i++) {
            if(pc == Ippl::myNode()) {
                if(count < N_mean) {
                    beam.create(1);
                    if(pc != 0)
                        beam.R[lowMark + count] = bg.getCooridinate(Ippl::myNode() * N_mean + count + N_extra);
                    else
                        beam.R[lowMark + count] = bg.getCooridinate(count);
                    beam.P[lowMark + count] = Vector_t(0.0);
                    beam.Bin[lowMark + count] = 0;
                    beam.PType[lowMark + count] = 1; // create darkcurrent particles.
                    beam.TriID[lowMark + count] = 0;
                    beam.Q[lowMark + count] = beam.getChargePerParticle();
                    beam.LastSection[lowMark + count] = 0;
                    beam.Ef[lowMark + count] = Vector_t(0.0);
                    beam.Bf[lowMark + count] = Vector_t(0.0);
                    beam.dt[lowMark + count] = beam.getdT();
                    count ++;

                }
            }
            pc++;
            if(pc == Ippl::getNodes())
                pc = 0;

        }

    }
    bg.clearCooridinateArray();
    beam.boundp();
    *gmsg << &beam << endl;

}


void  Distribution::createPriPart(PartBunch *beam, BoundaryGeometry &bg) {

    if(Options::ppdebug) {  // This is Parallel Plate Benchmark.
        int pc = 0;
        size_t lowMark = beam->getLocalNum();
        double vw = this->getVw();
        double vt = this->getvVThermal();
        double f_max = vw / vt * exp(-0.5);
        double test_a = vt / vw;
        double test_asq = test_a * test_a;
        size_t count = 0;
        size_t N_mean = static_cast<size_t>(floor(bg.getN() / Ippl::getNodes()));
        size_t N_extra = static_cast<size_t>(bg.getN() - N_mean * Ippl::getNodes());
        if(Ippl::myNode() == 0)
            N_mean += N_extra;
        if(bg.getN() != 0) {
            for(size_t i = 0; i < bg.getN(); i++) {
                if(pc == Ippl::myNode()) {
                    if(count < N_mean) {
                        /*==============Parallel Plate Benchmark=====================================*/
                        double test_s = 1;
                        double f_x = 0;
                        double test_x = 0;
                        while(test_s > f_x) {
                            test_s = IpplRandom();
                            test_s *= f_max;
                            test_x = IpplRandom();
                            test_x *= 10 * test_a; //range for normalized emission speed(0,10*test_a);
                            f_x = test_x / test_asq * exp(-test_x * test_x / 2 / test_asq);
                        }
                        double v_emi = test_x * vw;

                        double betaemit = v_emi / Physics::c;
                        double betagamma = betaemit / sqrt(1 - betaemit * betaemit);
                        /*============================================================================ */
                        beam->create(1);
                        if(pc != 0) {
                            beam->R[lowMark + count] = bg.getCooridinate(Ippl::myNode() * N_mean + count + N_extra);
                            beam->P[lowMark + count] = betagamma * bg.getMomenta(Ippl::myNode() * N_mean + count);
                        } else {
                            beam->R[lowMark + count] = bg.getCooridinate(count);
                            beam->P[lowMark + count] = betagamma * bg.getMomenta(count);
                        }
                        beam->Bin[lowMark + count] = 0;
                        beam->PType[lowMark + count] = 0; // create primary particle bunch;
                        beam->TriID[lowMark + count] = 0;
                        beam->Q[lowMark + count] = beam->getChargePerParticle();
                        beam->LastSection[lowMark + count] = 0;
                        beam->Ef[lowMark + count] = Vector_t(0.0);
                        beam->Bf[lowMark + count] = Vector_t(0.0);
                        beam->dt[lowMark + count] = beam->getdT();
                        count ++;
                    }
                }
                pc++;
                if(pc == Ippl::getNodes())
                    pc = 0;
            }
            bg.clearCooridinateArray();
            bg.clearMomentaArray();
            beam->boundp();

        }

        *gmsg << *beam << endl;

    } else {// Normal procedure to create primary particles

        int pc = 0;
        size_t lowMark = beam->getLocalNum();
        size_t count = 0;
        size_t N_mean = static_cast<size_t>(floor(bg.getN() / Ippl::getNodes()));
        size_t N_extra = static_cast<size_t>(bg.getN() - N_mean * Ippl::getNodes());

        if(Ippl::myNode() == 0)
            N_mean += N_extra;
        if(bg.getN() != 0) {
            for(size_t i = 0; i < bg.getN(); i++) {

                if(pc == Ippl::myNode()) {
                    if(count < N_mean) {
                        beam->create(1);
                        if(pc != 0)
                            beam->R[lowMark + count] = bg.getCooridinate(Ippl::myNode() * N_mean + count + N_extra); // node 0 will emit the particle with coordinate ID from 0 to N_mean+N_extra, so other nodes should shift to node_number*N_mean+N_extra
                        else
                            beam->R[lowMark + count] = bg.getCooridinate(count); // for node0 the particle number N_mean =  N_mean + N_extra
                        beam->P[lowMark + count] = Vector_t(0.0);
                        beam->Bin[lowMark + count] = 0;
                        beam->PType[lowMark + count] = 0; // create primary particle bunch.
                        beam->TriID[lowMark + count] = 0;
                        beam->Q[lowMark + count] = beam->getChargePerParticle();
                        beam->LastSection[lowMark + count] = 0;
                        beam->Ef[lowMark + count] = Vector_t(0.0);
                        beam->Bf[lowMark + count] = Vector_t(0.0);
                        beam->dt[lowMark + count] = beam->getdT();
                        count++;

                    }

                }
                pc++;
                if(pc == Ippl::getNodes())
                    pc = 0;

            }

        }
        bg.clearCooridinateArray();
        beam->boundp();//fixme if bg.getN()==0?
    }
    *gmsg << *beam << endl;
}


/**
 *
 *
 * @param beam
 * @param Np
 * @param scan
 */
void Distribution::create(PartBunch &beam, size_t Np, bool scan) {

    if(beam.getTotalNum() != 0) {
        scan_m = scan;
        create(beam, beam.getLocalNum());
    } else {
        scan_m = false; // the first time we have to create particles

        // in case we create a symetric distribution
        // we only need to specify half of the particles
        if(Options::cZero)
            create(beam, Np / 2);
        else
            create(beam, Np);
    }
}

/**
 *
 *
 * @param beam
 * @param Np
 */
void Distribution::create(PartBunch &beam, size_t Np) {

    Inform msg("Distribution::create ");

    int ebins = (int) Attributes::getReal(itsAttr[NBIN]);
    bool isBinned = (ebins > 0);

    if(isBinned) {
        if(pbin_m)
            delete pbin_m;
        pbin_m = new PartBins(static_cast<int>(fabs(Attributes::getReal(itsAttr[NBIN]))), static_cast<int>(fabs(Attributes::getReal(itsAttr[SBIN]))));
        if(scan_m) {
            beam.destroy(beam.getLocalNum(), 0);
            beam.update();
            INFOMSG("In scan mode: deleted BIN structure and all particles in the bunch" << endl;);
        }
    } else {
        pbin_m = NULL;
    }

    beam.setTEmission(Attributes::getReal(itsAttr[TEMISSION]));
    beam.setNumBunch(1);
    const string disttype = Attributes::getString(itsAttr[DISTRIBUTION]);
    if(disttype == "GUNGAUSSFLATTOPTH") {
        binnDistribution(beam, Np, disttype);
    } else if(disttype == "BINOMIAL") {
        corr_m[0] = Attributes::getReal(itsAttr[CORRX]);
        corr_m[1] = Attributes::getReal(itsAttr[CORRY]);
        corr_m[2] = Attributes::getReal(itsAttr[CORRT]);
        corr_m[3] = Attributes::getReal(itsAttr[R61]);
        corr_m[4] = Attributes::getReal(itsAttr[R62]);
        corr_m[5] = Attributes::getReal(itsAttr[R51]);
        corr_m[6] = Attributes::getReal(itsAttr[R52]);

        Vector_t corr(Attributes::getReal(itsAttr[CORRX]),
                      Attributes::getReal(itsAttr[CORRY]),
                      Attributes::getReal(itsAttr[CORRT]));

        Vector_t sigX(Attributes::getReal(itsAttr[SIGMAX]),
                      Attributes::getReal(itsAttr[SIGMAY]),
                      Attributes::getReal(itsAttr[SIGMAT]));

        Vector_t sigPX(eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]), beam.getM()),
                       eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]), beam.getM()),
                       eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]), beam.getM()));

        Vector_t bincoef(Attributes::getReal(itsAttr[MX]),
                         Attributes::getReal(itsAttr[MY]),
                         Attributes::getReal(itsAttr[MT]));

        Vector_t emit;
        Vector_t alpha;
        Vector_t beta;
        Vector_t gamma;

        for(int j = 0; j < 3; j++) {
            double chi = asin(corr[j]);
            emit[j] = sigX[j] * sigPX[j] * cos(chi);
        }
        for(int j = 0; j < 3; j++) {
            beta[j]  = sigX[j] * sigX[j] / emit[j];
            gamma[j] = sigPX[j] * sigPX[j] / emit[j];
            alpha[j] = -corr[j] * sqrt(beta[j] * abs(gamma[j]));
        }
        msg << "About to create Binomial distribution -1 " << endl;
        createBinom(emit, alpha, beta, gamma, bincoef, beam, Np, isBinned);
    } else if(disttype == "GAUSS") {
        double corr[7];
        corr[0] = Attributes::getReal(itsAttr[CORRX]);
        corr[1] = Attributes::getReal(itsAttr[CORRY]);
        corr[2] = Attributes::getReal(itsAttr[CORRT]);
        corr[3] = Attributes::getReal(itsAttr[R61]);
        corr[4] = Attributes::getReal(itsAttr[R62]);
        corr[5] = Attributes::getReal(itsAttr[R51]);
        corr[6] = Attributes::getReal(itsAttr[R52]);
        double Hs2a = Attributes::getReal(itsAttr[SIGMAX]);
        double Hs2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]), beam.getM()); //in eV

        double Vs2a = Attributes::getReal(itsAttr[SIGMAY]);
        double Vs2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]), beam.getM()); //in eV

        double Ls2a = Attributes::getReal(itsAttr[SIGMAT]);
        double Ls2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]), beam.getM()); //in eV

        double avrgpt = eVtoBetaGamma(Attributes::getReal(itsAttr[PT]), beam.getM());
        double avrgt = Attributes::getReal(itsAttr[T]);

        RANLIB_class *rGen = new RANLIB_class(265314159, 4);

        unsigned int pc = 0;
        unsigned int count = 0;

        if(beam.isZPeriodic())
            INFOMSG("Distribution: uniform in z (periodic BC) z = " << -0.5 * beam.getGaBeLa() << " ... " << 0.5 * beam.getGaBeLa() << endl);
        for(unsigned int i = 0; i < Np; i++) {
            double x, y;      // generate independent Gaussians, then correlate and finaly scale

            x  = rGen->gauss(0.0, 1.0);
            y  = rGen->gauss(0.0, 1.0);
            double xx = x;
            double yy = y;
            double px0  = x * corr[0] + y * sqrt(1.0 - corr[0] * corr[0]);
            double x0   =  x * Hs2a;
            px0 *= Hs2b;

            x  = rGen->gauss(0.0, 1.0);
            y  = rGen->gauss(0.0, 1.0);
            double py0  = x * corr[1] + y * sqrt(1.0 - corr[1] * corr[1]);
            double y0   =  x * Vs2a;
            py0 *= Vs2b;

            double del0;
            double psi0;

            if(beam.isZPeriodic()) {
                /*
                  create uniform distribution
                */
                //del0 = (IpplRandom()-0.5)*p[2];
                //  psi0 = (IpplRandom()*beam.getGaBeLa()) - (0.5*beam.getGaBeLa());

            } else {
                x  = rGen->gauss(0.0, 1.0);
                y  = rGen->gauss(0.0, 1.0);
                double const tempa = 1.0 - corr[0] * corr[0];
                double const l32 = (corr[6] - corr[0] * corr[5]) / sqrt(abs(tempa)) * tempa / abs(tempa);
                double const tempb = 1 - corr[5] * corr[5] - l32 * l32;
                double const l33 = sqrt(abs(tempb)) * tempb / abs(tempb);
                psi0 = xx * corr[5] + yy * l32 + x * l33;
                double const l42 = (corr[4] - corr[0] * corr[3]) / sqrt(abs(tempa)) * tempa / abs(tempa);
                double const l43 = (corr[2] - corr[5] * corr[3] - l42 * l32) / l33;
                double const tempc = 1 - corr[3] * corr[3] - l42 * l42 - l43 * l43;
                double const l44 = sqrt(abs(tempc)) * tempc / abs(tempc);
                del0 = xx * corr[3] + yy * l42 + x * l43 + y * l44;
                psi0 = avrgt + psi0 * Ls2a;
                del0 = avrgpt + Ls2b * del0;
            }
            if(pc == (size_t) Ippl::myNode()) {
                if(!scan_m)
                    beam.create(1);
                beam.R[count] = Vector_t(x0, y0, psi0);
                beam.P[count] = Vector_t(px0, py0, del0);
                beam.Bin[count] = 0; // not initialized
                beam.PType[count] = 0;
                beam.TriID[count] = 0;
                count++;
            }
            pc++;
            if(pc == (size_t) Ippl::getNodes())
                pc = 0;
        }
    } else if(disttype == "FROMFILE") {
        *gmsg << "\n-------------------------------------------------------------" << endl;
        *gmsg << "     READ ININITAL DISTRIBUTION FROM FILE    " << endl;
        *gmsg << "     BE AWARE OF THE FACT THAT ONLY NODE 0 IS READING IN " << endl;
        *gmsg << "-------------------------------------------------------------\n" << endl;

        if(isBinned) {
            *gmsg << "     DISTRIBUTION will be binned using " << ebins << " energy bins " << endl;
            const string fn;
            binnDistributionFromFile(beam, fn);
        } else {
            double x0, px0, y0, py0, psi0, del0;
            if(Ippl::myNode() == 0) {
                const string filename = Attributes::getString(itsAttr[FNAME]);
                std::ifstream fs;
                fs.open(filename.c_str());

                if(fs.fail()) {
                    throw OpalException("Distribution::Create()",
                                        "Open file operation failed, please check if \""
                                        + filename +  "\" really exists.");
                }

                fs >> Np;
                if(Np <= 0) {
                    throw OpalException("Distribution::Create()",
                                        " The particle number should be bigger than zero! Please check the first line of file \""
                                        + filename +  "\".");
                }

                for(unsigned int i = 0; i < Np; i++) {
                    if(!fs.eof()) {
                        beam.create(1);
                        fs >> x0 >> px0 >> y0 >> py0 >> psi0 >> del0;
                        beam.R[i] = Vector_t(x0, y0, psi0);
                        beam.P[i] = Vector_t(px0, py0, del0);
                        beam.Bin[i] = 0; // not initialized
                        beam.Q[i] = beam.getChargePerParticle();
                    } else {
                        throw OpalException("Distribution::Create()",
                                            "End of file reached before all particles imported, please check file \""
                                            + filename +  "\".");
                        return;
                    }
                }
                fs.close();
                tEmission_m = 0.0;
            }
        }
        /*
          In the case of a binned distribution (gun)
          we have to do the boundp after emission.
        */
    }

    if(!(isBinned)) {
        beam.boundp();
        beam.LastSection = 0;
    }
}

double Distribution::getTEmission() {
    if(tEmission_m > 0.0) {
        return tEmission_m;
    }

    distT_m = Attributes::getString(itsAttr[DISTRIBUTION]);
    if(distT_m == "GAUSS")
        distrTypeT_m = GAUSS;
    else if(distT_m == "GUNGAUSSFLATTOPTH")
        distrTypeT_m = GUNGAUSSFLATTOPTH;
    else if(distT_m == "FROMFILE")
        distrTypeT_m = FROMFILE;
    else if(distT_m == "BINOMIAL")
        distrTypeT_m = BINOMIAL;

    tPulseLengthFWHM_m = Attributes::getReal(itsAttr[TPULSEFWHM]);
    cutoff_m = Attributes::getReal(itsAttr[CUTOFF]);
    tRise_m = Attributes::getReal(itsAttr[TRISE]);
    tFall_m = Attributes::getReal(itsAttr[TFALL]);
    double tratio = sqrt(2.0 * log(10.0)) - sqrt(2.0 * log(10.0 / 9.0));
    sigmaRise_m = tRise_m / tratio;
    sigmaFall_m = tFall_m / tratio;

    switch(distrTypeT_m) {
        case ASTRAFLATTOPTH: {
            double a = tPulseLengthFWHM_m / 2;
            double sig = tRise_m / 2;
            double inv_erf08 = 0.906193802436823; // erfinv(0.8)
            double sqr2 = sqrt(2.);
            double t = a - sqr2 * sig * inv_erf08;
            double tmps = sig;
            double tmpt = t;
            for(int i = 0; i < 10; ++ i) {
                sig = (t + tRise_m - a) / (sqr2 * inv_erf08);
                t = a - 0.5 * sqr2 * (sig + tmps) * inv_erf08;
                sig = (0.5 * (t + tmpt) + tRise_m - a) / (sqr2 * inv_erf08);
                tmps = sig;
                tmpt = t;
            }
            tEmission_m = tPulseLengthFWHM_m + 10 * sig;
            break;
        }
        case GUNGAUSSFLATTOPTH: {
            tEmission_m = tPulseLengthFWHM_m + (cutoff_m - sqrt(2.0 * log(2.0))) * (sigmaRise_m + sigmaFall_m);
            break;
        }
        default:
            tEmission_m = 0.0;
    }
    return tEmission_m;
}

/**
 *
 *
 * @param beam
 * @param Np
 * @param distType
 */
void Distribution::binnDistribution(PartBunch &beam, size_t Np, string distType) {
    const double &two_pi = Physics::two_pi;
    unsigned int pc = 0;

    double dEBins = Attributes::getReal(itsAttr[DEBIN]);
    pbin_m->setRebinEnergy(dEBins);

    double eneg = Attributes::getReal(itsAttr[PT]);
    double gamma = 1. + (eneg / beam.getM());
    double beta  = sqrt(1. - (1. / (gamma * gamma)));
    double bega  = beta * gamma;

    pbin_m->setGamma(gamma); //has to be done by all processors!

    double corr[3] = {Attributes::getReal(itsAttr[CORRX]),
                      Attributes::getReal(itsAttr[CORRY]),
                      Attributes::getReal(itsAttr[CORRT])
                     };

    double transvCutOff = Attributes::getReal(itsAttr[TRANSVCUTOFF]);
    double Hs2a = Attributes::getReal(itsAttr[SIGMAX]);

    double Vs2a = Attributes::getReal(itsAttr[SIGMAY]);

    double Ls2a = Attributes::getReal(itsAttr[SIGMAT]);
    double Ls2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]), beam.getM());

    double nBins = Attributes::getReal(itsAttr[NBIN]);
    tPulseLengthFWHM_m = Attributes::getReal(itsAttr[TPULSEFWHM]);
    cutoff_m = Attributes::getReal(itsAttr[CUTOFF]);
    tRise_m = Attributes::getReal(itsAttr[TRISE]);
    tFall_m = Attributes::getReal(itsAttr[TFALL]);
    sigmaRise_m = tRise_m / 1.6869;
    sigmaFall_m = tFall_m / 1.6869;

    tEmission_m = tPulseLengthFWHM_m + (cutoff_m - 1.77441) / 1.6869 * (tRise_m + tFall_m);
    tBin_m = tEmission_m / nBins;
    double riseTime = tRise_m;
    double fallTime = tFall_m;

    gsl_rng_env_setup();
    const double sq2pi = sqrt(2.0 * Physics::pi);
    const double flatTop = tPulseLengthFWHM_m - sqrt(2.0 * log(2.0)) * (sigmaRise_m + sigmaFall_m);
    const double totA = flatTop + 0.5 * sq2pi * (sigmaRise_m + sigmaFall_m);
    const unsigned int particlesFront = Np * 0.5 * gsl_sf_erf(cutoff_m / sqrt(2.0)) * sq2pi * sigmaRise_m / totA;
    const unsigned int particlesTail  = Np * 0.5 * gsl_sf_erf(cutoff_m / sqrt(2.0)) * sq2pi * sigmaFall_m / totA;

    RANLIB_class *rGen = new RANLIB_class(265314159, 4);

    /*
      prepare quantities for thermal emittance calculation
    */
    bool   astraMode = false;
    double workf = 0.0;         // eV
    double siglaser = 0.0;      // m
    double elaser = 0.0;        // eV
    double fe = 0.0;            // Fermi energy eV
    double ag = 0.0;            // Acceleration gradient eV/m
    double ekin = 0.0;          // eV
    double phimax = 0.0;        // rad
    double schottky = 0.0;      // eV
    double ptot = 0.0;          // beta gamma
    std::ofstream os;

    if(distType == "GUNGAUSSFLATTOPTH") {

        ekin = Attributes::getReal(itsAttr[EKIN]);
        workf = Attributes::getReal(itsAttr[W]);
        elaser = Attributes::getReal(itsAttr[ELASER]);

        if((astraMode = (ekin > 0.0))) {
            // ASTRA mode
            phimax = Physics::pi / 2.0;
            *gmsg << " ---------------------------------------------------------------" << endl;
            *gmsg << " ---------------------I N P U T --------------------------------" << endl;
            *gmsg << " GUNGAUSS FLAT TOP &  THERMAL EMITTANCE in ASTRA MODE" << endl;
            *gmsg << " Kinetic energy = " << ekin << " [eV]  " << endl;
            *gmsg << " Phi max = " << phimax * 180 / Physics::pi << " [deg]  " << endl;
            *gmsg << " ---------------------------------------------------------------" << endl;

            ptot = eVtoBetaGamma(ekin, beam.getM());

        } else {
            workf = Attributes::getReal(itsAttr[W]);
            siglaser = Attributes::getReal(itsAttr[SIGLASER]);
            elaser = Attributes::getReal(itsAttr[ELASER]);
            fe = Attributes::getReal(itsAttr[FE]);
            ag = Attributes::getReal(itsAttr[AG]) * 1000000;

            ekin = fe + elaser;

            schottky = sqrt(Physics::q_e * ag / 4 / Physics::pi / Physics::epsilon_0);

            workf -= schottky;

            phimax = acos(sqrt((fe + workf) / ekin));

            ptot = eVtoBetaGamma((workf + elaser - fe - schottky), beam.getM()),

            *gmsg << " ---------------------------------------------------------------" << endl;
            *gmsg << " ---------------------I N P U T --------------------------------" << endl;
            *gmsg << " GUNGAUSS FLAT TOP &  THERMAL EMITTANCE " << endl;
            *gmsg << " Laser energy = " << elaser << " [eV]  rrms= " << siglaser << endl;
            *gmsg << " Work function = " << Attributes::getReal(itsAttr[W]) << " [eV]  " << endl;
            *gmsg << " Fermi energy = " << fe << " [eV]  " << endl;
            *gmsg << " Shottky energy = " << schottky << " [eV]  " << endl;
            *gmsg << " ---------------------------------------------------------------" << endl;
            *gmsg << " -------D E R I V E D  Q U A N T I T I E S ---------------------" << endl;
            *gmsg << " Kinetic energy = " << ekin << " [eV]  " << endl;
            *gmsg << " Effective Work function = " << workf << " [eV]  " << endl;
            *gmsg << " Phi max = " << phimax * 180 / Physics::pi << " [deg]  " << endl;
            *gmsg << " Ptot = " << workf + elaser - fe - schottky << " [eV]  " << endl;

            *gmsg << " ---------------------------------------------------------------" << endl;
        }
    }

    if(Ippl::getNodes() == 1) {
        *gmsg << " Write distribution to file dist.dat" << endl;
        string file("data/dist.dat");
        os.open(file.c_str());
        if(os.bad()) {
            *gmsg << "Unable to open output file " <<  file << endl;
        }
        os << "# x px y py z pz phi theta Ethx [eV] Ethy [eV] Ethz [eV]" << endl;
    }

    for(unsigned int i = 0; i < Np; i++) {
        double x = 0.0, y = 0.0;      // generate independent Gaussians, then correlate and finally scale
        double u1, u2;
        double xy = 6;

        while(xy > transvCutOff || std::isnan(xy)) {
            u1 = rGen->uniform(0.0, 1.0);
            u2 = rGen->uniform(0.0, 1.0);
            x = sqrt(-2. * log(1. - u1)) * cos(two_pi * u2);
            y = sqrt(-2. * log(1. - u1)) * sin(two_pi * u2);
            xy = sqrt(x * x + y * y);
        }

        double x0   =  x * Hs2a;
        double y0   =  y * Vs2a;

        double del0 = 0.0;
        double psi0 = 0.0;

        if(distType == "GUNGAUSSFLATTOPTH") {
            if(i < particlesFront) {
                // Fill rise.
                psi0 = -1.0;
                while(psi0 < 0.0)
                    psi0 = rGen->gauss(0.0, 1.0);
                psi0 = psi0 * riseTime + flatTop;
            } else if(i >= particlesFront && i < particlesFront + particlesTail) {
                // Fill fall.
                psi0 = 1.0;
                while(psi0 > 0.0)
                    psi0 = rGen->gauss(0.0, 1.0);
                psi0 = psi0 * fallTime;
            } else {
                // Fill flat top.
                psi0 = rGen->uniform(0.0, 1.0);
                psi0 *= flatTop;
            }
        }

        if(distType != "GUNGAUSSFLATTOPTH") {
            del0  = x * corr[2] + y * sqrt(1.0 - corr[2] * corr[2]);
            psi0  = x * Ls2a;
        }

        /*
          This is the longitudinal momenta in units of beta gamma.
        */

        del0  = bega + Ls2b * del0;

        double px0 = 0.0;
        double py0 = 0.0;
        double phi = 0.0;
        double theta = 0.0;

        if(distType == "GUNGAUSSFLATTOPTH") {

            /*
              Now calculate the thermal emittances
            */

            phi   = 2.0 * acos(sqrt(rGen->uniform(0.0, 1.0)));
            theta = 2.0 * Physics::pi * rGen->uniform(0.0, 1.0);

            px0  = ptot * sin(phi) * cos(theta);
            py0  = ptot * sin(phi) * sin(theta);
            del0 = bega + (ptot * abs(cos(phi)));
        }

        if(pc == (size_t) Ippl::myNode()) {
            vector<double> tmp;
            tmp.push_back(x0);
            tmp.push_back(y0);
            tmp.push_back(psi0);
            tmp.push_back(px0);
            tmp.push_back(py0);
            tmp.push_back(del0);
            tmp.push_back(0);
            pbin_m->fill(tmp);
        }
        pc++;
        if(pc == (size_t) Ippl::getNodes())
            pc = 0;

        if(Ippl::getNodes() == 1) {
            os << x0 << "\t " << px0    << "\t "
               << y0 << "\t " << py0    << "\t "
               << psi0 << "\t " << del0 << "\t "
               << phi * 180. / Physics::pi << "\t " << theta * 180. / Physics::pi << "\t "
               << betaGammatoeV(px0, beam.getM())  << "\t "
               << betaGammatoeV(py0, beam.getM())  << "\t "
               << betaGammatoeV(del0, beam.getM()) << "\t " << endl;
        }

    }
    if(Ippl::getNodes() == 1)
        os.close();

    pbin_m->sortArray();
    // now copy this over to the bunch
    // so that we can emit the particles
    beam.setPBins(pbin_m);
}


/**
 *
 *
 * @param beam
 * @param fn
 */
void Distribution::binnDistributionFromFile(PartBunch &beam, const string fn) {
    unsigned int pc = 0;
    size_t Np;
    double x0, y0, psi0, px0, py0, del0;

    std::ofstream os;
    std::ifstream fs;

    double dEBins = Attributes::getReal(itsAttr[DEBIN]);
    pbin_m->setRebinEnergy(dEBins);

    double eneg = Attributes::getReal(itsAttr[PT]);
    double gamma = 1. + (eneg / beam.getM());

    pbin_m->setGamma(gamma); //has to be done by all processors!

    if(Ippl::getNodes() == 1) {
        *gmsg << " Write distribution to file dist.dat" << endl;
        string file("data/dist.dat");
        os.open(file.c_str());
        if(os.bad()) {
            *gmsg << "Unable to open output file " <<  file << endl;
        }
        os << "# x px y py z pz phi theta Ethx [eV] Ethy [eV] Ethz [eV]" << endl;
    }

    fs.open(fn.c_str());

    if(fs.fail()) {
        throw OpalException("Distribution::Create()",
                            "Open file operation failed, please check if \""
                            + fn +  "\" really exists.");
    }

    fs >> Np;

    if(Np <= 0) {
        throw OpalException("Distribution::Create()",
                            " The particle number should be bigger than zero! Please check the first line of file \""
                            + fn +  "\".");
    }

    for(unsigned int i = 0; i < Np; i++) {
        if(!fs.eof()) {
            fs >> x0 >> px0 >> y0 >> py0 >> psi0 >> del0;
            if(pc == (size_t) Ippl::myNode()) {
                vector<double> tmp;
                tmp.push_back(x0);
                tmp.push_back(y0);
                tmp.push_back(psi0);
                tmp.push_back(px0);
                tmp.push_back(py0);
                tmp.push_back(del0);
                tmp.push_back(0);
                pbin_m->fill(tmp);
            }
            pc++;
            if(pc == (size_t) Ippl::getNodes())
                pc = 0;

            if(Ippl::getNodes() == 1) {
                os << x0 << "\t " << px0    << "\t "
                   << y0 << "\t " << py0    << "\t "
                   << psi0 << "\t " << del0 << endl;
            }
        }
    }

    if(Ippl::getNodes() == 1)
        os.close();
    fs.close();

    pbin_m->sortArray();
    // now copy this over to the bunch
    // so that we can emit the particles
    beam.setPBins(pbin_m);
}

/**
 *
 *
 * @param beam
 * @param Np
 * @param restartStep
 */
void Distribution::doRestart(PartBunch &beam, size_t Np, int restartStep) {
    h5_file_t *H5file;
    h5_int64_t rc;
    string fn;

    IpplTimings::startTimer(beam.distrReload_m);

    //        beam.setTEmission(Attributes::getReal(itsAttr[TEMISSION]));
    fn = OpalData::getInstance()->getInputBasename() + string(".h5");

#ifdef PARALLEL_IO
    H5file = H5OpenFile(fn.c_str(), H5_O_RDONLY, Ippl::getComm());
#else
    H5file = H5OpenFile(fn.c_str(), H5_O_RDONLY, 0);
#endif

    if(!H5file) {
        ERRORMSG("could not open file '" << fn << "';  exiting!" << endl);
        exit(0);
    }

    if(restartStep == -1) {
        restartStep = H5GetNumSteps(H5file) - 1 ;
        OpalData::getInstance()->setRestartStep(restartStep);
    } else {
        if(restartStep != H5GetNumSteps(H5file) - 1 && !OpalData::getInstance()->hasRestartFile()) {
            ERRORMSG("can't append to the file '" << fn << "' exiting!" << endl);
            exit(0);
        }
    }

    rc = H5SetStep(H5file, restartStep);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    int N = (int)H5PartGetNumParticles(H5file);

    //    h5_int64_t totalSteps = H5GetNumSteps(H5file);

    int numberOfParticlesPerNode = (int) floor((double) N / Ippl::getNodes());
    long long starti = Ippl::myNode() * numberOfParticlesPerNode;
    long long endi = starti + numberOfParticlesPerNode - 1;

    // In case we miss some particles we add them at the end on the last core
    if(Ippl::myNode() == (Ippl::getNodes() - 1)) {
        if(Ippl::getNodes()*numberOfParticlesPerNode < N)
            endi += (N - (Ippl::getNodes() * numberOfParticlesPerNode));
    }

    rc = H5PartSetView(H5file, starti, endi);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    N = (int)H5PartGetNumParticles(H5file);
    assert(N >= 0);

    double actualT;
    rc = H5ReadStepAttribFloat64(H5file, "TIME", &actualT);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    beam.setT(actualT);

    double dPhiGlobal;
    rc = H5ReadFileAttribFloat64(H5file, "dPhiGlobal", &dPhiGlobal);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    OpalData::getInstance()->setGlobalPhaseShift(dPhiGlobal);

    h5_int64_t ltstep;
    rc = H5ReadStepAttribInt64(H5file, "LocalTrackStep", &ltstep);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    beam.setLocalTrackStep((long long)ltstep);

    h5_int64_t gtstep;
    rc = H5ReadStepAttribInt64(H5file, "GlobalTrackStep", &gtstep);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    beam.setGlobalTrackStep((long long)gtstep);

    std::unique_ptr<char[]> varray(new char[(N)*sizeof(double)]);
    h5_float64_t *farray = reinterpret_cast<h5_float64_t *>(varray.get());
    h5_int64_t *larray = reinterpret_cast<h5_int64_t *>(varray.get());

    beam.create(N);

    rc = H5PartReadDataFloat64(H5file, "x", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n) {
        beam.R[n](0) = farray[n];
        beam.Bin[n] = 0; // not initialized
    }
    rc = H5PartReadDataFloat64(H5file, "y", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n)
        beam.R[n](1) = farray[n];

    rc = H5PartReadDataFloat64(H5file, "z", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n)
        beam.R[n](2) = farray[n];

    rc = H5PartReadDataFloat64(H5file, "px", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n)
        beam.P[n](0) = farray[n];

    rc = H5PartReadDataFloat64(H5file, "py", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n)
        beam.P[n](1) = farray[n];

    rc = H5PartReadDataFloat64(H5file, "pz", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n)
        beam.P[n](2) = farray[n];

    rc = H5PartReadDataFloat64(H5file, "q", farray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n)
        beam.Q[n] = farray[n];

    rc = H5PartReadDataInt64(H5file, "lastsection", larray);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    for(unsigned long int n = 0; n < (unsigned int) N; ++n)
        beam.LastSection[n] = (short) larray[n];

    Ippl::Comm->barrier();
    rc = H5CloseFile(H5file);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    beam.boundp();

    IpplTimings::stopTimer(beam.distrReload_m);

    *gmsg << "Total number of particles in the h5 file = " << N << " NPerBunch= " << beam.getTotalNum()
          << " Global step " << gtstep << " Local step " << ltstep << endl
          << " restart step= " << restartStep << " time of restart = " << actualT
          << " phishift= " << OpalData::getInstance()->getGlobalPhaseShift() << endl;
}



/**
 *
 *
 * @param beam
 * @param Np
 * @param restartStep
 * @param specifiedNumBunch
 */
void Distribution::doRestart_cycl(PartBunch &beam, size_t Np, int restartStep, const int specifiedNumBunch) {
    h5_int64_t rc;
    IpplTimings::startTimer(beam.distrReload_m);
    *gmsg << "---------------- Start reading hdf5 file----------------" << endl;
    h5_file_t *H5file;

    string fn = OpalData::getInstance()->getInputBasename() + string(".h5");

#ifdef PARALLEL_IO
    H5file = H5OpenFile(fn.c_str(), H5_O_RDONLY, Ippl::getComm());
#else
    H5file = H5OpenFile(fn.c_str(), H5_O_RDONLY, 0);
#endif

    if(!H5file) {
        ERRORMSG("File open failed:  exiting!" << endl);
        exit(0);
    }

    if(restartStep == -1) {
        restartStep = H5GetNumSteps(H5file) - 1;
        OpalData::getInstance()->setRestartStep(restartStep);
    } else {
        if(restartStep != H5GetNumSteps(H5file) - 1 && !OpalData::getInstance()->hasRestartFile()) {
            ERRORMSG("can't append to the file '" << fn << "' exiting!" << endl);
            exit(0);
        }
    }

    *gmsg << "Restart from hdf5 format file " << fn << ", read phase space data of DumpStep " << (int)restartStep << endl;

    rc = H5SetStep(H5file, restartStep);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    const int globalN = (int)H5PartGetNumParticles(H5file);

    *gmsg << "total number of particles = " << globalN << endl;

    int numberOfParticlesPerNode = (int) floor((double) globalN / Ippl::getNodes());
    long long starti = Ippl::myNode() * numberOfParticlesPerNode;
    long long endi = 0;

    if(Ippl::myNode() == Ippl::getNodes() - 1)
        endi = -1;
    else
        endi = starti + numberOfParticlesPerNode;

    rc = H5PartSetView(H5file, starti, endi);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    const int localN = (int)H5PartGetNumParticles(H5file);
    assert(localN >= 0);

    double actualT;
    rc = H5ReadStepAttribFloat64(H5file, "TIME", &actualT);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    beam.setT(actualT);

    h5_int64_t ltstep;
    rc = H5ReadStepAttribInt64(H5file, "LocalTrackStep", &ltstep);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    beam.setLocalTrackStep((long long)ltstep);

    h5_int64_t gtstep;
    rc = H5ReadStepAttribInt64(H5file, "GlobalTrackStep", &gtstep);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    beam.setGlobalTrackStep((long long)gtstep);

    char opalFlavour[128];
    rc = H5ReadStepAttribString(H5file, "OPAL_flavour", opalFlavour);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);

    std::unique_ptr<char[]> varray(new char[(localN)*sizeof(double)]);
    double *farray = reinterpret_cast<double *>(varray.get());
    h5_int64_t *larray = reinterpret_cast<h5_int64_t *>(varray.get());

    beam.create(localN);

    if (strcmp(opalFlavour,"opal-t") == 0) {
        *gmsg<<"Restart from hdf5 file generated by OPAL-t"<<endl;
    
        rc = H5PartReadDataFloat64(H5file, "x", farray);
        if(rc != H5_SUCCESS)
          ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(unsigned long int n = 0; n < (unsigned int) localN; ++n)
          beam.R[n](0) = -farray[n];

        rc = H5PartReadDataFloat64(H5file, "y", farray);
        if(rc != H5_SUCCESS)
          ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(unsigned long int n = 0; n < (unsigned int) localN; ++n)
          beam.R[n](2) = farray[n];

        rc = H5PartReadDataFloat64(H5file, "z", farray);
        if(rc != H5_SUCCESS)
          ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(unsigned long int n = 0; n < (unsigned int) localN; ++n)
          beam.R[n](1) = farray[n];

        rc = H5PartReadDataFloat64(H5file, "px", farray);
        if(rc != H5_SUCCESS)
          ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(unsigned long int n = 0; n < (unsigned int) localN; ++n)
          beam.P[n](0) = -farray[n];

        rc = H5PartReadDataFloat64(H5file, "py", farray);
        if(rc != H5_SUCCESS)
          ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(unsigned long int n = 0; n < (unsigned int) localN; ++n)
          beam.P[n](2) = farray[n];

        rc = H5PartReadDataFloat64(H5file, "pz", farray);
        if(rc != H5_SUCCESS)
          ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(unsigned long int n = 0; n < (unsigned int) localN; ++n)
          beam.P[n](1) = farray[n];
        
        rc = H5PartReadDataInt64(H5file, "id", larray);
        if(rc != H5_SUCCESS)
          ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(unsigned long int n = 0; n < (unsigned int) localN; ++n)
          beam.ID[n] = larray[n];

    }else{
        *gmsg<<"Restart from hdf5 file generated by OPAL-cycl"<<endl;
   
        h5_int64_t SteptoLastInj;
        rc = H5ReadStepAttribInt64(H5file, "SteptoLastInj", &SteptoLastInj);
        if(rc != H5_SUCCESS)
          ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        beam.setSteptoLastInj((int)SteptoLastInj);
        *gmsg << "Tracking Step since last bunch injection is " << SteptoLastInj << endl;

        h5_int64_t numBunch;
        rc = H5ReadStepAttribInt64(H5file, "NumBunch", &numBunch);
        if(rc != H5_SUCCESS)
          ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        beam.setNumBunch((int)numBunch);
        *gmsg << numBunch << " Bunches(bins) exist in this file" << endl;

        rc = H5PartReadDataFloat64(H5file, "x", farray);
        if(rc != H5_SUCCESS)
          ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(unsigned long int n = 0; n < (unsigned int) localN; ++n)
          beam.R[n](0) = farray[n];

        rc = H5PartReadDataFloat64(H5file, "y", farray);
        if(rc != H5_SUCCESS)
          ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(unsigned long int n = 0; n < (unsigned int) localN; ++n)
          beam.R[n](1) = farray[n];

        rc = H5PartReadDataFloat64(H5file, "z", farray);
        if(rc != H5_SUCCESS)
          ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(unsigned long int n = 0; n < (unsigned int) localN; ++n)
          beam.R[n](2) = farray[n];

        rc = H5PartReadDataFloat64(H5file, "px", farray);
        if(rc != H5_SUCCESS)
          ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(unsigned long int n = 0; n < (unsigned int) localN; ++n)
          beam.P[n](0) = farray[n];

        rc = H5PartReadDataFloat64(H5file, "py", farray);
        if(rc != H5_SUCCESS)
          ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(unsigned long int n = 0; n < (unsigned int) localN; ++n)
          beam.P[n](1) = farray[n];

        rc = H5PartReadDataFloat64(H5file, "pz", farray);
        if(rc != H5_SUCCESS)
          ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(unsigned long int n = 0; n < (unsigned int) localN; ++n)
          beam.P[n](2) = farray[n];

        rc = H5PartReadDataInt64(H5file, "id", larray);
        if(rc != H5_SUCCESS)
          ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
        for(unsigned long int n = 0; n < (unsigned int) localN; ++n)
          beam.ID[n] = larray[n];
        
        // only for multi-bunch mode
        if(specifiedNumBunch > 1) {
          // the allowed maximal bin number is set to 1000
          beam.setPBins(new PartBinsCyc(1000, numBunch));
        
        }

    }

    Ippl::Comm->barrier();
    rc = H5CloseFile(H5file);
    if(rc != H5_SUCCESS)
        ERRORMSG("H5 rc= " << rc << " in " << __FILE__ << " @ line " << __LINE__ << endl);
    beam.boundp();
    beam.Q = beam.getChargePerParticle();

    if (strcmp(opalFlavour,"opal-t") == 0) {
        Vector_t meanR(0.0, 0.0, 0.0);
        Vector_t meanP(0.0, 0.0, 0.0);
        unsigned long int newLocalN = beam.getLocalNum();
        for(unsigned int i = 0; i < newLocalN; ++i) {
          for(int d = 0; d < 3; ++d) {
              meanR(d) += beam.R[i](d);
              meanP(d) += beam.P[i](d);
          }
        }
        reduce(meanR, meanR, OpAddAssign());
        meanR /= Vector_t(globalN);
        reduce(meanP, meanP, OpAddAssign());
        meanP /= Vector_t(globalN);
        *gmsg <<"Rmean = " <<meanR<<"[m], Pmean="<<meanP<<endl;
        
        for(unsigned int i = 0; i < newLocalN; ++i) {
          beam.R[i] -= meanR;
          beam.P[i] -= meanP;
        }
    }

    *gmsg << "----------------Finish reading hdf5 file----------------" << endl;
    IpplTimings::stopTimer(beam.distrReload_m);
}


/**
 *
 *
 * @param name
 *
 * @return
 */
Distribution *Distribution::find(const string &name) {
    Distribution *dist = dynamic_cast<Distribution *>(OpalData::getInstance()->find(name));

    if(dist == 0) {
        throw OpalException("Distribution::find()", "Distribution \"" + name + "\" not found.");
    }

    return dist;
}

/*
  double Distribution::getET() const
  {
  return Attributes::getReal(itsAttr[ET]);
  }

  void Distribution::setET(double value)
  {
  Attributes::setReal(itsAttr[ET], value);
  }

*/

/**
 *
 *
 *
 * @return
 */
const PartData &Distribution::getReference() const {
    // Cast away const, to allow logically constant Distribution to update.
    const_cast<Distribution *>(this)->update();
    return reference;
}

/**
 *
 *
 */
void Distribution::update() {

}

double Distribution::getEkin() const {return Attributes::getReal(itsAttr[EKIN]);}
double Distribution::getWorkFunctionRf() const {return Attributes::getReal(itsAttr[W]);}
double Distribution::getLaserEnergy() const {return Attributes::getReal(itsAttr[ELASER]);}

size_t Distribution::getNumberOfDarkCurrentParticles() { return (size_t) Attributes::getReal(itsAttr[NPDARKCUR]);}
double Distribution::getDarkCurrentParticlesInwardMargin() { return Attributes::getReal(itsAttr[INWARDMARGIN]);}
double Distribution::getEInitThreshold() { return Attributes::getReal(itsAttr[EINITHR]);}

double Distribution::getWorkFunction() { return Attributes::getReal(itsAttr[FNPHIW]); }
double Distribution::getFieldEnhancement() { return Attributes::getReal(itsAttr[FNBETA]); }
size_t Distribution::getMaxFNemissionPartPerTri() { return (size_t) Attributes::getReal(itsAttr[FNMAXEMI]);}
double Distribution::getFieldFNThreshold() { return Attributes::getReal(itsAttr[FNFIELDTHR]);}
double Distribution::getFNParameterA() { return Attributes::getReal(itsAttr[FNA]);}
double Distribution::getFNParameterB() { return Attributes::getReal(itsAttr[FNB]);}
double Distribution::getFNParameterY() { return Attributes::getReal(itsAttr[FNY]);}
double Distribution::getFNParameterVYZero() { return Attributes::getReal(itsAttr[FNVYZERO]);}
double Distribution::getFNParameterVYSecond() { return Attributes::getReal(itsAttr[FNVYSECOND]);}

int    Distribution::getSecondaryEmissionFlag() { return Attributes::getReal(itsAttr[SECONDARYFLAG]);}
bool   Distribution::getEmissionMode() { return Attributes::getBool(itsAttr[NEMISSIONMODE]);}
string Distribution::getTypeofDistribution() { return (string) Attributes::getString(itsAttr[DISTRIBUTION]);}

double Distribution::getvSeyZero() {return Attributes::getReal(itsAttr[VSEYZERO]);}// return sey_0 in Vaughan's model
double Distribution::getvEZero() {return Attributes::getReal(itsAttr[VEZERO]);}// return the energy related to sey_0 in Vaughan's model
double Distribution::getvSeyMax() {return Attributes::getReal(itsAttr[VSEYMAX]);}// return sey max in Vaughan's model
double Distribution::getvEmax() {return Attributes::getReal(itsAttr[VEMAX]);}// return Emax in Vaughan's model
double Distribution::getvKenergy() {return Attributes::getReal(itsAttr[VKENERGY]);}// return fitting parameter denotes the roughness of surface for impact energy in Vaughan's model
double Distribution::getvKtheta() {return Attributes::getReal(itsAttr[VKTHETA]);}// return fitting parameter denotes the roughness of surface for impact angle in Vaughan's model
double Distribution::getvVThermal() {return Attributes::getReal(itsAttr[VVTHERMAL]);}// thermal velocity of Maxwellian distribution of secondaries in Vaughan's model
double Distribution::getVw() {return Attributes::getReal(itsAttr[VW]);}// velocity scalar for parallel plate benchmark;

int Distribution::getSurfMaterial() {return (int)Attributes::getReal(itsAttr[SURFMATERIAL]);}// Surface material number for Furman-Pivi's Model;

/**
 *
 *
 * @param os
 *
 * @return
 */
Inform &Distribution::printInfo(Inform &os) const {
    double transvCutOff = Attributes::getReal(itsAttr[TRANSVCUTOFF]);
    os << "* ************* D I S T R I B U T I O N ******************************************** " << endl;
    if(!OpalData::getInstance()->inRestartRun()) {
        os << "* Distribution:\t" << getOpalName() << endl;
        if(!(Attributes::getString(itsAttr[LASERPROFFN]) == string(""))) {
            os << "* Distribution type:\t" << Attributes::getString(itsAttr[DISTRIBUTION]) << endl;
            os << "* Laser profile: " << Attributes::getString(itsAttr[LASERPROFFN])
               << " Image: " << Attributes::getString(itsAttr[IMAGENAME])
               << " Intensity cut: " << Attributes::getReal(itsAttr[INTENSITYCUT]) << endl;
        } else {
            os << "* Distribution type:\t" << Attributes::getString(itsAttr[DISTRIBUTION]) << endl;
        }

        if(Attributes::getString(itsAttr[DISTRIBUTION]) == "SURFACEEMISSION") {
            os << "* Number of electrons for surface emission  " <<     Attributes::getReal(itsAttr[NPDARKCUR]) << endl;
            os << "* Initialized electrons inward margin for surface emission  " << Attributes::getReal(itsAttr[INWARDMARGIN]) << endl;
            os << "* E field threshold (MV), only in position r with E(r)>EINITHR that particles will be initialized   " << Attributes::getReal(itsAttr[EINITHR]) << endl;
            os << "* Field enhancement for surface emission  " << Attributes::getReal(itsAttr[FNBETA]) << endl;
            os << "* Maximum number of electrons emitted from a single triangle for Fowler-Nordheim emission  " << Attributes::getReal(itsAttr[FNMAXEMI]) << endl;
            os << "* Field Threshold for Fowler-Nordheim emission (MV/m)  " <<  Attributes::getReal(itsAttr[FNFIELDTHR]) << endl;
            os << "* Empirical constant A for Fowler-Nordheim emission model  " <<  Attributes::getReal(itsAttr[FNA]) << endl;
            os << "* Empirical constant B for Fowler-Nordheim emission model  " <<  Attributes::getReal(itsAttr[FNB]) << endl;
            os << "* Constant for image charge effect parameter y(E) in Fowler-Nordheim emission model  " <<  Attributes::getReal(itsAttr[FNY]) << endl;
            os << "* Zero order constant for image charge effect parameter v(y) in Fowler-Nordheim emission model  " <<  Attributes::getReal(itsAttr[FNVYZERO]) << endl;
            os << "* Second order constant for image charge effect parameter v(y) in Fowler-Nordheim emission model  " <<  Attributes::getReal(itsAttr[FNVYSECOND]) << endl;
            os << "* Select secondary model type(0:no secondary emission; 1:Furman-Pivi; 2 or larger: Vaughan's model  " <<  Attributes::getReal(itsAttr[SECONDARYFLAG]) << endl;
            os << "* Secondary emission mode type(true: emit n true secondaries; false: emit one particle with n times charge  " << Attributes::getBool(itsAttr[NEMISSIONMODE]) << endl;
            os << "* Sey_0 in Vaughan's model " <<  Attributes::getReal(itsAttr[VSEYZERO]) << endl;
            os << "* Energy related to sey_0 in Vaughan's model in eV  " <<  Attributes::getReal(itsAttr[VEZERO]) << endl;
            os << "* Sey max in Vaughan's model  " <<  Attributes::getReal(itsAttr[VSEYMAX]) << endl;
            os << "* Emax in Vaughan's model in eV  " <<  Attributes::getReal(itsAttr[VEMAX]) << endl;
            os << "* Fitting parameter denotes the roughness of surface for impact energy in Vaughan's model  " <<  Attributes::getReal(itsAttr[VKENERGY]) << endl;
            os << "* Fitting parameter denotes the roughness of surface for impact angle in Vaughan's model  " <<  Attributes::getReal(itsAttr[VKTHETA]) << endl;
            os << "* Thermal velocity of Maxwellian distribution of secondaries in Vaughan's model  " <<  Attributes::getReal(itsAttr[VVTHERMAL]) << endl;
            os << "* VW denote the velocity scalar for Parallel plate benchmark  " <<  Attributes::getReal(itsAttr[VW]) << endl;
            os << "* Material type number of the cavity surface for Furman-Pivi's model, 0 for copper, 1 for stainless steel  " <<  Attributes::getReal(itsAttr[SURFMATERIAL]) << endl;

        } else if(Attributes::getString(itsAttr[DISTRIBUTION]) == "SURFACERANDCREATE") {
            os << "* Number of electrons initialized on the surface as primaries  " <<     Attributes::getReal(itsAttr[NPDARKCUR]) << endl;
            os << "* Initialized electrons inward margin for surface emission  " << Attributes::getReal(itsAttr[INWARDMARGIN]) << endl;
            os << "* E field threshold (MV), only in position r with E(r)>EINITHR that particles will be initialized   " << Attributes::getReal(itsAttr[EINITHR]) << endl;
        } else if(Attributes::getString(itsAttr[DISTRIBUTION]) != "FROMFILE") {


            os << "* sigmax=\t" << Attributes::getReal(itsAttr[SIGMAX]) << " [m];\t"
               << "sigmay=\t" << Attributes::getReal(itsAttr[SIGMAY]) << " [m];\t"
               << "sigmat=\t" << Attributes::getReal(itsAttr[SIGMAT]) << " [m];" << endl;
            os << "* transverse cut  " << transvCutOff << " [sigma] " << endl;
            os << "* sigmapx=\t" << Attributes::getReal(itsAttr[SIGMAPX]) << " [eV];\t"
               << "sigmapy=\t" << Attributes::getReal(itsAttr[SIGMAPY]) << " [eV];\t"
               << "pt +- sigmapt=\t" << Attributes::getReal(itsAttr[PT]) << "+-" << Attributes::getReal(itsAttr[SIGMAPT]) << " [eV]" << endl;

            os << "* corr x-px=\t" << Attributes::getReal(itsAttr[CORRX]);
            os << "\t corr y-py=\t" << Attributes::getReal(itsAttr[CORRY]);
            os << "\t corr t-pt=\t" << Attributes::getReal(itsAttr[CORRT]) << endl;
            if(Attributes::getReal(itsAttr[TEMISSION]) > 0.0) {
                os << "* -------- G U N -------------------------------------------------------------------" << endl;
                os << "* temission [s] \t" << Attributes::getReal(itsAttr[TEMISSION]);
                os << "\t n bins \t" << (int) Attributes::getReal(itsAttr[NBIN]) << endl;
                os << "\t dEbin [keV] \t" << Attributes::getReal(itsAttr[DEBIN]) << endl;
                os << "* ----------------------------------------------------------------------------------" << endl;
            }
        }
    } else
        os << "* Distribution from restart file" << endl;
    os << "* ********************************************************************************** " << endl;
    return os;
}

