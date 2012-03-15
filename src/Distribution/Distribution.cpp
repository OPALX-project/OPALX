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
#include "Distribution/LaserProfile.h"
#include "Distribution/Distribution.h"
#include "BasicActions/Option.h"
#include "AbstractObjects/OpalData.h"
#include <fstream>
#include <string>
#include <cmath>
#include <cfloat>
#include <vector>
#include <iostream>  // Neeeded for stream I/O
#include <fstream>   // Needed for file I/O
#include <iomanip>   // Needed for I/O manipulators

#include "Algorithms/bet/SLPartBunch.h"
#include "Algorithms/PartBins.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

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
        PXMULT,
        PYMULT,
        PTMULT,
        BETAX,
        BETAY,
        ALPHAX,
        ALPHAY,
        MX,
        MY,
        MT,
        DX,
        DDX,
        DY,
        DDY,
        R51,
        R52,
        R61,
        R62,
        PT,
        T,
        SIGMAX,
        SIGMAY,
        SIGMAT,
        TRANSVCUTOFF,
        SIGMAPX,
        SIGMAPY,
        SIGMAPT,
        SIGMARISE,
        SIGMAFALL,
        FLATTOPTIME,
        CUTOFFRISE,
        CUTOFFFALL,
        CORRX,
        CORRY,
        CORRT,
		OFFSETX,
		OFFSETY,
        TEMISSION,
        NBIN,
        DEBIN,
        ELASER,
        SIGLASER,
        W,
        FE,
        AG,
		EKIN,
        SIZE
    };
}

/** 
 * Construtor
 * 
 */
Distribution::Distribution():
    Definition(SIZE,"DISTRIBUTION","The DISTRIBUTION statement defines data for the 6D particle distr."),
    distrTypeT_m(NODIST)
{
    itsAttr[DISTRIBUTION] = makeString("DISTRIBUTION", "Distribution type: GAUSS, BINOMIAL, ROTSYMBINOMIAL, FROMFILE,"
                                       "GUNGAUSS, GUNGAUSS3D, GUNUNIFORM, GUNGAUSSFLATTOP, GUNGAUSSFLATTOPTH, UNIFORMXYZ ", "GAUSS");

    itsAttr[FNAME] = makeString("FNAME", "File for read in 6D particle coordinates");

    itsAttr[LASERPROFFN] = makeString("LASERPROFFN", "File for read in a measured laser profile (x,y)","");
    itsAttr[IMAGENAME] = makeString("IMAGENAME", "Name of the image");
    itsAttr[INTENSITYCUT] = makeReal("INTENSITYCUT", "For background substraction, in % of max intensity",0.0);


    itsAttr[XMULT] = makeReal("XMULT", "Multiplier for X",1.0);
    itsAttr[YMULT] = makeReal("YMULT", "Multiplier for Y",1.0);
    itsAttr[TMULT] = makeReal("TMULT", "Multiplier for T",1.0);
    itsAttr[TRANSVCUTOFF] = makeReal("TRANSVCUTOFF", "Transverse cut-off in units of sigma",3.0);

    itsAttr[PXMULT]= makeReal("PXMULT", "Multiplier for PX",1.0);
    itsAttr[PYMULT]= makeReal("PYMULT", "Multiplier for PY",1.0);
    itsAttr[PTMULT]= makeReal("PTMULT", "Multiplier for PT",1.0);

    itsAttr[ALPHAX]= makeReal("ALPHAX", "Courant Synder parameter",1.0);
    itsAttr[ALPHAY]= makeReal("ALPHAY", "Courant Synder parameter",1.0);

    itsAttr[BETAX] = makeReal("BETAX", "Courant Synder parameter",-1.0);
    itsAttr[BETAY] = makeReal("BETAY", "Courant Synder parameter",1.0);

    itsAttr[MX]    = makeReal("MX", "Defines the distribution in x, 0+eps .. inf",1.0);
    itsAttr[MY]    = makeReal("MY", "Defines the distribution in y, 0+eps .. inf",1.0);
    itsAttr[MT]    = makeReal("MT", "Defines the distribution in t, 0+eps .. inf",1.0);

    itsAttr[DX]    = makeReal("DX", "Dispersion in x (R16 in Transport notation",0.0);
    itsAttr[DDX]   = makeReal("DDX", "First derivative of Dx",0.0);


    itsAttr[DY]    = makeReal("DY", "DY",0.0);
    itsAttr[DDY]   = makeReal("DDY", "DDY",0.0);

    itsAttr[R51]    = makeReal("R51", "R51",0.0);
    itsAttr[R52]   = makeReal("R52", "R52",0.0);

    itsAttr[R61]    = makeReal("R61", "R61",0.0);
    itsAttr[R62]   = makeReal("R62", "R62",0.0);

    itsAttr[PT] = makeReal("PT", "average longitudinal momentum",0.0);
    itsAttr[T] = makeReal("T","average longitudinal position", 0.0);

    itsAttr[SIGMAX] = makeReal("SIGMAX", "SIGMAx (m)",1.0e-2);
    itsAttr[SIGMAY] = makeReal("SIGMAY", "SIGMAy (m)",1.0e-2);
    itsAttr[SIGMAT] = makeReal("SIGMAT", "SIGMAt (m)",1.0e-2);

    itsAttr[SIGMAPX] = makeReal("SIGMAPX", "SIGMApx",0.0);
    itsAttr[SIGMAPY] = makeReal("SIGMAPY", "SIGMApy",0.0);
    itsAttr[SIGMAPT] = makeReal("SIGMAPT", "SIGMApt",0.0);

    itsAttr[CORRX] = makeReal("CORRX", "CORRx",-0.5);
    itsAttr[CORRY] = makeReal("CORRY", "CORRy", 0.5);
    itsAttr[CORRT] = makeReal("CORRT", "CORRt", 0.0);

    itsAttr[OFFSETX] = makeReal("OFFSETX", "OFFSETx", 0.0);
    itsAttr[OFFSETY] = makeReal("OFFSETY", "OFFSETy", 0.0);

    itsAttr[TEMISSION] = makeReal("TEMISSION", "Time in seconds in which we have emission",  0.0);
    itsAttr[NBIN]      = makeReal("NBIN", "In case of emission how many bins should we use", 0.0);
    itsAttr[DEBIN]     = makeReal("DEBIN", "Energy band for a bin in keV, defines the rebinning", 1000000.0);

    itsAttr[SIGMARISE]   = makeReal("SIGMARISE", "Sigma rise for GUNGAUSSFLATTOP distribution type (m)", 1.0e-2);
    itsAttr[SIGMAFALL]   = makeReal("SIGMAFALL", "Sigma fall for GUNGAUSSFLATTOP distribution type (m)", 1.0e-2);
    itsAttr[FLATTOPTIME] = makeReal("FLATTOPTIME", "Flat top time for GUNGAUSSFLATTOP distribution type (m)", 1.0e-2);
    itsAttr[CUTOFFRISE]  = makeReal("CUTOFFRISE", "Rise cutoff for GUNGAUSSFLATTOP distribution type (m)", 1.0e-2);
    itsAttr[CUTOFFFALL]  = makeReal("CUTOFFFALL", "Fall cutoff for GUNGAUSSFLATTOP distribution type (m)", 1.0e-2);

    itsAttr[ELASER] = makeReal("ELASER", "Laser energy (eV)", 0.0);
    itsAttr[SIGLASER] = makeReal("SIGLASER", "Sigma of (uniform) laser spot size (m)", 0.0);
    itsAttr[W] = makeReal("W", "Workfunction of material (eV)", 0.0);
    itsAttr[FE] = makeReal("FE", "Fermi energy (eV)", 0.0);
    itsAttr[AG] = makeReal("AG", "Acceleration Gradient (MV/m)", 0.0);

    itsAttr[EKIN] = makeReal("EKIN", "Ekin used in ASTRA (eV)", -1.0);

    // Set up default beam.
    Distribution *defDistribution = clone("UNNAMED_Distribution");
    defDistribution->builtin = true;

    try {
        defDistribution->update();
        OPAL.define(defDistribution);
    } catch (...) {
        delete defDistribution;
    }
    pbin_m = NULL;
    lp_m = NULL;

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
    distrTypeT_m(NODIST)
{
    pbin_m = NULL;
    lp_m = NULL;
}

/** 
 * Destructor
 * 
 */
Distribution::~Distribution()
{
    if (pbin_m) {
        delete pbin_m;
        pbin_m = NULL;
    }

    if ((Ippl::getNodes()==1) && (os_m.is_open()))
		os_m.close();
	
    if(lp_m) {
		delete lp_m;
		lp_m = NULL;
    }
}

/** 
 * This is the main entrypoint! 
 * 
 * @param beam 
 * @param Np 
 * @param scan 
 */
void Distribution::setup(PartBunch &beam, size_t Np, bool scan) {

    scan_m = scan;
    int nBins_m = (int) Attributes::getReal(itsAttr[NBIN]);

    bool isBinned = (nBins_m > 0);

    if (isBinned) {
        if (pbin_m)
            delete pbin_m;
        pbin_m = new PartBins((int) Attributes::getReal(itsAttr[NBIN]));
    }
    else 
        pbin_m = NULL;

    if (scan_m) {
		beam.destroy(beam.getLocalNum(),0);
		beam.update();
		INFOMSG("In scan mode: delete all particles in the bunch" << endl;);
    }

    laserProfileFn_m=Attributes::getString(itsAttr[LASERPROFFN]);
    
    if (!(laserProfileFn_m==string(""))) {
		laserImage_m  =Attributes::getString(itsAttr[IMAGENAME]);
		intensityCut_m=Attributes::getReal(itsAttr[INTENSITYCUT]);    
		lp_m = new LaserProfile(laserProfileFn_m,laserImage_m,intensityCut_m);
    }

    beam.setTEmission(Attributes::getReal(itsAttr[TEMISSION]));
    beam.setNumBunch(1);

    distT_m = Attributes::getString(itsAttr[DISTRIBUTION]);
    if (distT_m == "GAUSS")
	distrTypeT_m = GAUSS;
    else if (distT_m == "GUNGAUSSFLATTOPTH")
	distrTypeT_m = GUNGAUSSFLATTOPTH;
    else if (distT_m == "FROMFILE")
	distrTypeT_m = FROMFILE;
    else if (distT_m == "UNIFORMXYZ")
	distrTypeT_m = UNIFORMXYZ;
    else if (distT_m == "BINOMIAL")
	distrTypeT_m = BINOMIAL;


    // Create a an initial beam bunch that is:
    // "GUNGAUSS": uniform in space transversely and with a Gaussian ("GUNGAUS") longitudinal profile
    // "GUNUNIFORM": uniform in space transversely and longitudinally.
    // "GUNGAUSS3D": Gaussian transversely and longitudinally.
    // "GUNGAUSSFLATTOP": uniform in space transversely, a Gaussian rise and fall time longitudinally with
    //                    a uniform flattop between.
    // "GUNGAUSSFLATTOPTH": uniform in space transversely, a Gaussian rise and fall time longitudinally with
    //                    a uniform flattop between, and a transvers thermal emittance

    switch (distrTypeT_m) {
	
    case GUNGAUSSFLATTOPTH:
		{
			const double &two_pi = Physics::two_pi;
	
			double dEBins = Attributes::getReal(itsAttr[DEBIN]);

			pbin_m->setRebinEnergy(dEBins);

			corr_m[0] = Attributes::getReal(itsAttr[CORRX]);
			corr_m[1] = Attributes::getReal(itsAttr[CORRY]);
			corr_m[2] = Attributes::getReal(itsAttr[CORRT]);

			nBins_m = Attributes::getReal(itsAttr[NBIN]);
			transvCutOff_m = Attributes::getReal(itsAttr[TRANSVCUTOFF]);

			Hs2a_m = Attributes::getReal(itsAttr[SIGMAX]);
			Hs2b_m = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]),beam.getM());
	
			Vs2a_m = Attributes::getReal(itsAttr[SIGMAY]);
			Vs2b_m = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]),beam.getM());

			Ls2a_m = Attributes::getReal(itsAttr[SIGMAT]);
			Ls2b_m = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]),beam.getM());

			tEmission_m = Attributes::getReal(itsAttr[TEMISSION]);
			sigmaRise_m = Attributes::getReal(itsAttr[SIGMARISE]);
			sigmaFall_m = Attributes::getReal(itsAttr[SIGMAFALL]);
			flatTop_m = Attributes::getReal(itsAttr[FLATTOPTIME]);
			cutoffRise_m = Attributes::getReal(itsAttr[CUTOFFRISE]);
			cutoffFall_m = Attributes::getReal(itsAttr[CUTOFFFALL]);

			//sigma_m = 3.0;
			tBin_m = tEmission_m/nBins_m;

			rGen_m = new RANLIB_class(265314159,4);

			h_m = gsl_histogram_alloc (nBins_m);
			//gsl_histogram_set_ranges_uniform (h_m, 0, tEmission_m);
			//createTimeBins(Np, sigma_m, nBins_m, riseTime_m, fallTime_m, tEmission_m, h_m);
			createTimeBins(Np, sigmaRise_m, sigmaFall_m, nBins_m, cutoffRise_m, cutoffFall_m, flatTop_m, h_m); 
			pbin_m->setHistogram(h_m);

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
			ptot_m = eVtoBetaGamma(ekin_m,beam.getM());
    
			// ASTRA mode
			phimax_m = Physics::pi/2.0;
			*gmsg << " -- B I N N I N G in T -----------------------------------------" << endl;
			*gmsg << " ---------------------I N P U T --------------------------------" << endl;
			*gmsg << " GUNGAUSS FLAT TOP &  THERMAL EMITTANCE in ASTRA MODE" << endl;
			*gmsg << " Kinetic energy (thermal emittance) = " << ekin_m << " [eV]  " << endl;
			*gmsg << " Phi max = " << phimax_m*180/Physics::pi << " [deg]  " << endl;
			*gmsg << " tBin = " << tBin_m << " [sec]  nBins = " << nBins_m << " tEmission =  " << tEmission_m << " [sec] " << endl;
    
			if (Ippl::getNodes() == 1) {
				*gmsg << " Write distribution to file dist.dat" << endl;
				string file("dist.dat");
				os_m.open(file.c_str());
				if (os_m.bad()) {
					*gmsg << "Unable to open output file " <<  file << endl;
				}
				os_m << "# x y ti px py pz Ekin= " << ekin_m << " [eV] " << endl;
			}
		}
		break;
      
    case BINOMIAL: 
		{

		    corr_m = Vector_t(Attributes::getReal(itsAttr[CORRX]),
				      Attributes::getReal(itsAttr[CORRY]),
				      Attributes::getReal(itsAttr[CORRT]));
		    
		    sigx_m = Vector_t(Attributes::getReal(itsAttr[SIGMAX]),
				      Attributes::getReal(itsAttr[SIGMAY]),
				      Attributes::getReal(itsAttr[SIGMAT]));
	
		    sigp_m = Vector_t(eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]),beam.getM()),
				      eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]),beam.getM()),
				      eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]),beam.getM()));
		    
		    binc_m = Vector_t(Attributes::getReal(itsAttr[MX]),
				      Attributes::getReal(itsAttr[MY]),
				      Attributes::getReal(itsAttr[MT]));
		    
		    
		    for(int j=0; j<3; j++) {
			double chi = asin(corr_m[j]);
			emit_m[j] = sigx_m[j]*sigp_m[j]*cos(chi);
		    }
		    for(int j=0; j<3; j++) {
			beta_m[j]  = sigx_m[j]*sigx_m[j]/emit_m[j];
			gamma_m[j] = sigp_m[j]*sigp_m[j]/emit_m[j];
			alpha_m[j] = -corr_m[j]*sqrt(beta_m[j]*abs(gamma_m[j]));
		    }
		    createBinom(emit_m, alpha_m, beta_m, gamma_m, binc_m, beam, Np, isBinned);
		}
		break;
      
    case GAUSS:
		{
			gauss_corr_m[0] = Attributes::getReal(itsAttr[CORRX]);
			gauss_corr_m[1] = Attributes::getReal(itsAttr[CORRY]);
			gauss_corr_m[2] = Attributes::getReal(itsAttr[CORRT]);
			gauss_corr_m[3] = Attributes::getReal(itsAttr[R61]);
			gauss_corr_m[4] = Attributes::getReal(itsAttr[R62]);
			gauss_corr_m[5] = Attributes::getReal(itsAttr[R51]);
			gauss_corr_m[6] = Attributes::getReal(itsAttr[R52]);
			gauss_offset_m[0]= Attributes::getReal(itsAttr[OFFSETX]);
			gauss_offset_m[1]= Attributes::getReal(itsAttr[OFFSETY]);

			Hs2a_m = Attributes::getReal(itsAttr[SIGMAX]);
			Hs2b_m = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]),beam.getM()); //in eV

			Vs2a_m = Attributes::getReal(itsAttr[SIGMAY]);
			Vs2b_m = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]),beam.getM()); //in eV

			Ls2a_m = Attributes::getReal(itsAttr[SIGMAT]);
			Ls2b_m = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]),beam.getM()); //in eV
	
			avrgpt_m = eVtoBetaGamma(Attributes::getReal(itsAttr[PT]),beam.getM());
			avrgt_m  = Attributes::getReal(itsAttr[T]);

			/* 
			   give up the portability w.r.t. the rangen
			   and hope to be more scalable
			*/
			rGen_m = new RANLIB_class((Ippl::myNode()+1)*265314159,4);

			IpplTimings::startTimer(beam.distrCreate_m);	
			
			if (Np > 1E8) {
				int k = 10;
				Np = (size_t)Np/Ippl::getNodes()/k;
				*gmsg << "Sampl= " << Np*Ippl::getNodes() << " x " << k << " Total= " << k*Np*Ippl::getNodes() <<  endl;
				for (int kk=0; kk<k; kk++) {
					sampleGauss(beam,Np);
					beam.boundp();
					*gmsg << "Sampl Gauss k= " << kk << " N= " << beam.getTotalNum() << endl;
				}
			}
			else {
				Np = (size_t) Np / Ippl::getNodes();
				sampleGauss(beam,Np);
			}
			*gmsg << "Sample Gauss done ..." << endl;
	
			IpplTimings::stopTimer(beam.distrCreate_m);
		}
		break;
    case UNIFORMXYZ: 
		{
			corr_m[0] = Attributes::getReal(itsAttr[CORRX]);
			corr_m[1] = Attributes::getReal(itsAttr[CORRY]);
			corr_m[2] = Attributes::getReal(itsAttr[CORRT]);
	
			Hs2a_m = Attributes::getReal(itsAttr[SIGMAX]);
			Hs2b_m = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]),beam.getM()); //in eV

			Vs2a_m = Attributes::getReal(itsAttr[SIGMAY]);
			Vs2b_m = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]),beam.getM()); //in eV

			Ls2a_m = Attributes::getReal(itsAttr[SIGMAT]);
			Ls2b_m = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]),beam.getM()); //in eV   

			nBins_m = Attributes::getReal(itsAttr[NBIN]);

			rGen_m = new RANLIB_class(265314159,4);
		}
		break;
    case FROMFILE:
		{
			*gmsg <<"---------------------------------------------"<<endl;
			*gmsg <<"     READ ININITAL DISTRIBUTION FROM FILE    "<<endl;
			*gmsg <<"     BE AWARE OF THE FACT THAT ONLY NODE 0 IS READING IN "<<endl;

			if (isBinned) {
				*gmsg <<"     DISTRIBUTION will be binned using " << nBins_m << " energy bins " << endl;
				const string fn;
				binnDistributionFromFile(beam, fn);
				/*

				std::ofstream os;
				if (Ippl::getNodes() == 1) {
				*gmsg << " Write distribution to file dist.dat" << endl;
				string file("dist.dat");
				os.open(file.c_str());
				if (os.bad()) {
				*gmsg << "Unable to open output file " <<  file << endl;
				}
				os << "# x px y py z pz " << endl;
				}
				os.close();	
				*/
			}
			else {

				if (Ippl::myNode() == 0) {
					const string filename=Attributes::getString(itsAttr[FNAME]);
					double x0,px0,y0,py0,psi0,del0;
		
					std::ifstream fs;
					fs.open(filename.c_str());

					if(fs.fail()){
						throw OpalException("Distribution::Create()",
											"Open file operation failed, please check if \""
											+ filename +  "\" really exists.");
					}
	    
					fs >> Np;
					if ( Np <= 0 ){
						throw OpalException("Distribution::Create()",
											" The particle number should be bigger than zero! Please check the first line of file \""
											+ filename +  "\".");
					}

					for(unsigned int i=0; i<Np; i++) {
						if ( !fs.eof() ){
							// create 1 particle
							beam.create(1);
							fs>> x0 >> px0 >> y0 >> py0 >> psi0 >> del0;
							beam.R[i] = Vector_t(x0,y0,psi0);
							beam.P[i] = Vector_t(px0,py0,del0);
							beam.Bin[i] = 0; // not initialized
							beam.Q[i] = beam.getChargePerParticle();
		    
							//if (Ippl::getNodes() == 1) {
							//   os << x0 << "\t " << px0    << "\t "
							//      << y0 << "\t " << py0    << "\t "
							//      << psi0 << "\t " << del0 << "\t " << endl;
							//}
						}
						else {
							throw OpalException("Distribution::Create()",
												"End of file reached before all particles imported, please check file \""
												+ filename +  "\".");
							return;
						}
					}
					fs.close();
	
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

    if (isBinned) 
        beam.setPBins(pbin_m);
    else
		beam.boundp();
    beam.LastSection=0;
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

    for(size_t i=beam.getLocalNum(); i<Np; i++) {
		double x,y;       // generate independent Gaussians, then correlate and finaly scale
		x  = rGen_m->gauss(0.0,1.0);
		y  = rGen_m->gauss(0.0,1.0);
		double xx =x;
		double yy =y;
		double px0  = x*gauss_corr_m[0] + y*sqrt(1.0 - gauss_corr_m[0]*gauss_corr_m[0]);
		double x0   = x*Hs2a_m+gauss_offset_m[0];
		px0 *= Hs2b_m;

		x  = rGen_m->gauss(0.0,1.0);
		y  = rGen_m->gauss(0.0,1.0);
		double py0  = x*gauss_corr_m[1] + y*sqrt(1.0 - gauss_corr_m[1]*gauss_corr_m[1]);
		double y0   =  x*Vs2a_m+gauss_offset_m[1];
		py0 *= Vs2b_m;
	
		double del0;
		double psi0;
	
		x  = rGen_m->gauss(0.0,1.0);
		y  = rGen_m->gauss(0.0,1.0);
		const double l32=(gauss_corr_m[6]-gauss_corr_m[0]*gauss_corr_m[5])/sqrt(1.0 - gauss_corr_m[0]*gauss_corr_m[0]);
		const double l33=sqrt(1-gauss_corr_m[5]*gauss_corr_m[5]-l32*l32);
		psi0 = xx*gauss_corr_m[5]+yy*l32+x*l33;
		const double l42=(gauss_corr_m[4]-gauss_corr_m[0]*gauss_corr_m[3])/sqrt(1.0 - gauss_corr_m[0]*gauss_corr_m[0]);
		const double l43=(gauss_corr_m[2]-gauss_corr_m[5]*gauss_corr_m[3]-l42*l32)/l33;
		const double l44=sqrt(1-gauss_corr_m[3]*gauss_corr_m[3]-l42*l42-l43*l43);
		del0 = xx*gauss_corr_m[3]+yy*l42+x*l43+y*l44;
		psi0 = avrgt_m + psi0*Ls2a_m;
		del0 = avrgpt_m + Ls2b_m*del0;
		if (pc == Ippl::myNode()) {
			beam.create(1);
			beam.R[count] = Vector_t(x0,y0,psi0);
			beam.P[count] = Vector_t(px0,py0,del0);
			beam.Bin[count] = 0; // not initialized
			beam.Q[count] = beam.getChargePerParticle();
			count++;
		}
		pc++;
		if (pc == Ippl::getNodes())
			pc = 0;        
    }
	INFOMSG("Sample Gauss done" << endl);
}

/** 
 * 
 * 
 * @param dt 
 * 
 * @return 
 */
pair<Vector_t,Vector_t> Distribution::sample(double dt) {

    Vector_t r(0.0);
    Vector_t p(0.0);
    

    double phi   = 0.0;
    double theta = 0.0;

    double px0   = 0.0;
    double py0   = 0.0;
    double del0  = 0.0;

    double x,y;
    double x0,y0;
    double xy = 6;

    switch (distrTypeT_m) {
	
    case GUNGAUSSFLATTOPTH:
		{
			double x,y;
			double xy = 6;
	
			if (lp_m != NULL) {
				lp_m->GetXY(&x,&y);
				x = 2*x - 1.0;
				y = 2*y - 1.0;
			}
			else {
				while (xy > 1){
					x  = rGen_m->uniform(-1.0,1.0);
					y  = rGen_m->uniform(-1.0,1.0);
					xy = sqrt(x*x + y*y);
				}
			}
	
			x0   =  x*Hs2a_m;
			y0   =  y*Vs2a_m;

			/*
			  Now calculate the thermal emittances
			*/
	    
			const double phi   = 2.0 * acos(sqrt(rGen_m->uniform(0.0, 1.0)));
			const double theta = 2.0 * Physics::pi * rGen_m->uniform(0.0, 1.0);

			const double px0   = ptot_m * sin(phi) * cos(theta);
			const double py0   = ptot_m * sin(phi) * sin(theta);
			const double del0  = ptot_m * abs(cos(phi));
	    
            r = Vector_t(x0,y0,rGen_m->uniform(0.0,1.0)*dt);
			p = Vector_t(px0,py0,del0);
	
			//if (Ippl::getNodes() == 1)
				//os_m << r(0) << "\t " 
					 //<< r(1) << "\t " 
					 //<< r(2) << "\t " 
					 //<< p(0) << "\t " 
					 //<< p(1) << "\t " 
					 //<< p(2) << "\t " 
					 //<< phi*180./Physics::pi << "\t " << theta*180./Physics::pi << endl;
		}
		break;
      
    default:
		INFOMSG("Distribution unknown" << endl;);
    }	
    return pair<Vector_t,Vector_t>(r,p);
}

/** 
 * 
 * 
 * @param Np 
 * @param sigma_rise 
 * @param sigma_fall 
 * @param bins 
 * @param cutoff_rise 
 * @param cutoff_fall 
 * @param tflat 
 * @param h 
 */
void Distribution::createTimeBins(const int Np, const double sigma_rise, const double sigma_fall, const unsigned int bins, 
								  const double cutoff_rise, const double cutoff_fall, const double tflat, gsl_histogram * h) 
{

    int i = 0;

    const double trise = (1+cutoff_rise)*sigma_rise;
    const double tfall = (1+cutoff_fall)*sigma_fall;
    gsl_histogram_set_ranges_uniform (h, 0, trise+tfall+tflat);

    gsl_rng_env_setup ();
    gsl_rng *r = gsl_rng_alloc (gsl_rng_default);

    const double pi = Physics::pi;
    //does not take into account cutoff
    const double totA = tflat + 0.5 * sqrt(2.0 * pi) * (sigma_rise + sigma_fall);
    //const int nrise = Np * 0.5 * sqrt(2.0 * pi) * sigma_rise / totA;
    //const int nfall = Np * 0.5 * sqrt(2.0 * pi) * sigma_fall / totA;
    //const int nflat = Np - nrise - nfall;

    //taking into account cutoff
    const int nrise = Np * (gsl_cdf_ugaussian_P(1+cutoff_rise) - 0.5) * sqrt(2.0*pi) * sigma_rise / totA;
    const int nfall = Np * (gsl_cdf_ugaussian_P(1+cutoff_fall) - 0.5) * sqrt(2.0*pi) * sigma_fall / totA;
    const int nflat = Np - nrise - nfall;

    //printf("nflat = %i\n", nflat);
    //printf("nrise = %i\n", nrise);
    //printf("nfall = %i\n", nfall);

    for (i=0; i<nrise; i++) {
        const double r1 = gsl_ran_gaussian_tail(r, 0, sigma_rise);
        gsl_histogram_increment(h,-r1+trise);
    }
    for (i=0; i<nfall; i++) {
        const double r1 = gsl_ran_gaussian_tail(r, 0, sigma_fall);
        gsl_histogram_increment(h,r1+tfall+tflat);
    }
    for (i=0; i<nflat; i++) {
        const double r1 = gsl_ran_flat(r, trise, trise+tflat);
        gsl_histogram_increment(h, r1);
    }  
    gsl_rng_free (r);

}

//void Distribution::createTimeBins(const int Np, const double sigma, const unsigned int bins, const double trise, const double tfall, 
//const double temis, gsl_histogram * h) 
//{
  
//const double tflat = temis - trise - tfall;
//int i = 0;

//gsl_rng_env_setup ();
//gsl_rng *r = gsl_rng_alloc (gsl_rng_default);
  
//int nback;
//const int nfront = nback = (int) Np * (2*gsl_cdf_ugaussian_P (1) - 1.0) / 2 / sigma; // (P-Q);
//const int nflat  = Np - nfront - nback;

//// INFORM("Np= " << Np << " nback=nfront= " << nback << " Sigma= " << sigma << " trise= " << trise << " tflat= " << tflat << endl;);

//for (i=0;i<nback;i++) {
//// front & back
//const double r1 = gsl_ran_gaussian_tail(r, 0, trise/sigma);
//gsl_histogram_increment (h, r1+tflat+tfall);
//gsl_histogram_increment (h,-r1+trise);
    
//}
//for (i=0;i<nflat;i++) {
//// flat top part
//const double r1 = gsl_ran_flat (r, trise, trise+tflat);
//gsl_histogram_increment (h, r1);
//}    
//gsl_rng_free (r);
  
//}


#ifdef HAVE_ENVELOPE_SOLVER
/** 
 * 
 * 
 * @param p 
 */
void Distribution::createSlicedBunch(double charge, double gamma, double mass, double current, double center, double Bz0, SLPartBunch* p)
{
    double beamWidth = 0.0;
    double beamEnergy = 0.0;
    int sl = (int) Attributes::getReal(itsAttr[NBIN]);
    *gmsg << "About to create a sliced bunch with " << sl << " slices" << endl;
    *gmsg << "mass = " << mass << " gamma = " << gamma << endl;
    
    distT_m = Attributes::getString(itsAttr[DISTRIBUTION]);
    if (distT_m == "GAUSS")
	distrTypeT_m = GAUSS;
    else if (distT_m == "GUNGAUSSFLATTOPTH")
	distrTypeT_m = GUNGAUSSFLATTOPTH;
    else if (distT_m == "FROMFILE")
	distrTypeT_m = FROMFILE;
    else if (distT_m == "UNIFORMXYZ")
	distrTypeT_m = UNIFORMXYZ;
    else if (distT_m == "BINOMIAL")
	distrTypeT_m = BINOMIAL;
    
    switch (distrTypeT_m) {
	
    case GUNGAUSSFLATTOPTH:
        {
			corr_m[0] = Attributes::getReal(itsAttr[CORRX]);
			corr_m[1] = Attributes::getReal(itsAttr[CORRY]);
			corr_m[2] = Attributes::getReal(itsAttr[CORRT]);

			nBins_m = Attributes::getReal(itsAttr[NBIN]);

			Hs2a_m = Attributes::getReal(itsAttr[SIGMAX]);
			Vs2a_m = Attributes::getReal(itsAttr[SIGMAY]);
			Ls2a_m = Attributes::getReal(itsAttr[SIGMAT]);

            tEmission_m = Attributes::getReal(itsAttr[TEMISSION]);
			ekin_m = Attributes::getReal(itsAttr[EKIN]);
            
            // EnvelopeTracker expects [eV]
            beamEnergy = ekin_m;
            beamWidth = tEmission_m*Physics::c*sqrt(1.0 - (1.0/pow(gamma,2)));

            break;
        }
    case GAUSS:
		{
			gauss_corr_m[0] = Attributes::getReal(itsAttr[CORRX]);
			gauss_corr_m[1] = Attributes::getReal(itsAttr[CORRY]);
			gauss_corr_m[2] = Attributes::getReal(itsAttr[CORRT]);
			gauss_corr_m[3] = Attributes::getReal(itsAttr[R61]);
			gauss_corr_m[4] = Attributes::getReal(itsAttr[R62]);
			gauss_corr_m[5] = Attributes::getReal(itsAttr[R51]);
			gauss_corr_m[6] = Attributes::getReal(itsAttr[R52]);
			gauss_offset_m[0]= Attributes::getReal(itsAttr[OFFSETX]);
			gauss_offset_m[1]= Attributes::getReal(itsAttr[OFFSETY]);

			Hs2a_m = Attributes::getReal(itsAttr[SIGMAX]);
			Vs2a_m = Attributes::getReal(itsAttr[SIGMAY]);
			Ls2a_m = Attributes::getReal(itsAttr[SIGMAT]);
	
			avrgt_m  = Attributes::getReal(itsAttr[T]);

            //FIXME: why 1e9??
            beamEnergy = (gamma*mass-mass)*1e9;
            beamWidth = Ls2a_m;

            break;
        }
    default:
    }
    
    center = -1*beamWidth/2.0;
    *gmsg << "x = " << Hs2a_m << " y = " << Vs2a_m << endl;
    *gmsg << "beamWidth = " << beamWidth << " beamCenter = " << center << " beamEnergy = " << beamEnergy << endl;
    //FIXME: what is this?
    double frac = 0.0;

    // get starting energy
    //double energy = gamma*mass-mass;
    //double width = sigmaT; // z width of the beam 
    //double Betenerg = energy *1e9;

    // execute initialization command 
    p->iniBetBunch(sl, charge, beamEnergy, beamWidth, frac, current, center, Hs2a_m, Vs2a_m, 0, 0, Bz0);
}
#endif

/** 
 * 
 * 
 * @param object 
 * 
 * @return 
 */
bool Distribution::canReplaceBy(Object *object)
{
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
Distribution *Distribution::clone(const string &name)
{
    return new Distribution(name, this);
}

/** 
 * 
 * 
 */
void Distribution::execute()
{
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
                               bool isBinned)
{
    const double two_pi = 2.0 * 4.0*atan(1.0);

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


    for (int i=0;i<3;i++) {
        gamma[i] *= 4.0;
        beta[i] *= 4.0;
        M[i]       =  sqrt(emit[i]*beta[i]);
        PM[i]      =  sqrt(emit[i]*gamma[i]);
        COSCHI[i]  =  sqrt(1.0/(1.0+alpha[i]*alpha[i]));
        SINCHI[i]  = -alpha[i]*COSCHI[i];
        CHI[i]     =  atan2(SINCHI[i],COSCHI[i]);
        AMI[i]     =  1.0/bincoef[i];
        L[i]       =  sqrt((bincoef[i]+1.0)/2.0)*M[i];
        PL[i]      =  sqrt((bincoef[i]+1.0)/2.0)*PM[i];
    }

    Vector_t x;
    // now copy this over to the bunch
    Vector_t p;
    double betagamma_part;
    if (isBinned)
        {
            double ekin = Attributes::getReal(itsAttr[PT]);
            double mass = beam.getM();
            double gamma_part = 1. + ekin/mass;
            betagamma_part = sqrt(ekin*ekin/(mass*mass) + 2.*ekin/mass );
            pbin_m->setGamma(gamma_part);
            *gmsg << "* Gamma = " << gamma_part << "; Beta = " << betagamma_part / gamma_part  << endl;
        }

    for (size_t n=0; n < particles; ++n) {
        double S1,S2;
        double A,AL,U,V;
        for (int i=0;i<3;i++) {
            S1=IpplRandom();
            S2=IpplRandom();
            if (bincoef[i] <= 10000)
                {
                    A=sqrt(1.0 - pow(S1,AMI[i]));
                    AL=two_pi*S2;
                    U=A*cos(AL);
                    V=A*sin(AL);
                    x[i] = L[i]*U;
                    p[i] = PL[i]*(U*SINCHI[i]+V*COSCHI[i]);
                }
            else
                {
                    A=sqrt(2.0)/2.0 * sqrt(-log(S1));
                    AL=two_pi*S2;
                    U=A*cos(AL);
                    V=A*sin(AL);
                    x[i]=M[i]*U;
                    p[i]=PM[i]*(U*SINCHI[i]+V*COSCHI[i]);
                }
            p[2] += betagamma_part;
        }
        if (pc == Ippl::myNode()) {
            if (isBinned) {
                vector<double> tmp;
                tmp.push_back(x[0]);
                tmp.push_back(x[1]);
                tmp.push_back(x[2]);
                tmp.push_back(p[0]);
                tmp.push_back(p[1]);
                tmp.push_back(p[2]);
                tmp.push_back(0);
                pbin_m->fill(tmp);
            }
            else {
                beam.create(1);
                beam.R[count] = x;
                beam.P[count] = p;
		beam.Q[count] = beam.getChargePerParticle();
                count++;
            }
        }
	
        if (!isBinned)
            {
                pc++;
                if (pc == Ippl::getNodes())
                    pc = 0;
            }
    }
    if (isBinned) {
        pbin_m->sortArray();
        // now copy this over to the bunch
        // so that we can emmit the particles
        beam.setPBins(pbin_m);
    }
}

/** 
 * 
 * 
 * @param emit 
 * @param alpha 
 * @param beta 
 * @param gamma 
 * @param beam 
 * @param particles 
 * @param isBinned 
 */
void Distribution::createUniformTUniformL(Vector_t emit, Vector_t alpha, Vector_t beta, Vector_t gamma,
                                          PartBunch &beam, size_t particles, bool isBinned)
{
    const double &two_pi = Physics::two_pi;

    gsl_rng_env_setup ();
    gsl_rng *r = gsl_rng_alloc (gsl_rng_default);

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

    for (int i=0;i<3;i++) {
        gamma[i]  *= 4.0;
        beta[i]   *= 4.0;
        M[i]       =  sqrt(emit[i]*beta[i]);
        PM[i]      =  sqrt(emit[i]*gamma[i]);
        COSCHI[i]  =  sqrt(1.0/(1.0+alpha[i]*alpha[i]));
        SINCHI[i]  = -alpha[i]*COSCHI[i];
        CHI[i]     =  atan2(SINCHI[i],COSCHI[i]);
    }

    Vector_t x;
    Vector_t p;
    double betagamma_part;
    if (isBinned)
        {
            double ekin = Attributes::getReal(itsAttr[PT]);
            double mass = beam.getM();
            double gamma_part = 1. + ekin/mass;
            betagamma_part = sqrt(ekin*ekin/(mass*mass) + 2.*ekin/mass );
            pbin_m->setGamma(gamma_part);
            *gmsg << "* Gamma = " << gamma_part << "; Beta = " << betagamma_part / gamma_part  << endl;
        }

    for (size_t n=0; n < particles; ++n) {
        double S1,S2,S3,S4,S5,S6;
        double A1,AL1,AC1,U1,V1,A2,AL2,AC2,U2,V2;
        double AB = 2.0;
        while (AB > 1.0)
            {
                S1=IpplRandom();
                S2=IpplRandom();
                S3=IpplRandom();
                S4=IpplRandom();
                S5=IpplRandom();
                S6=IpplRandom();
                AL1 = two_pi*S1;
                AL2 = two_pi*S2;
                A1=sqrt(1.0 - S3);
                A2=sqrt(1.0 - S4);
                AC1 = cos(AL1);
                U1=A1*AC1;
                AC2 = cos(AL2);
                U2=A2*AC2;
                AB = sqrt(U1*U1 + U2*U2);
            }

        V1=A1*sin(AL1);
        V2=A2*sin(AL2);
        x[0] = M[0]*U1;
        x[1] = M[1]*U2;
        p[0] = PM[0]*(U1*SINCHI[0]+V1*COSCHI[0]);
        p[1] = PM[1]*(U2*SINCHI[1]+V2*COSCHI[1]);


        S1=IpplRandom();
        S2=IpplRandom();

        A1=sqrt(1.0 - S5);
        AL1 = two_pi * S6;
        U1=A1*cos(AL1);
        // now copy this over to the bunch
        V1=A1*sin(AL1);
        x[2] = M[2] * U1;
        p[2] = betagamma_part + PM[2] * (U1 * SINCHI[2] + V1 * COSCHI[2]);


        if (pc == Ippl::myNode()) {
            if (isBinned) {
                vector<double> tmp;
                tmp.push_back(x[0]);
                tmp.push_back(x[1]);
                tmp.push_back(x[2]);
                tmp.push_back(p[0]);
                tmp.push_back(p[1]);
                tmp.push_back(p[2]);
                tmp.push_back(0);
                pbin_m->fill(tmp);
            }
            else {
                beam.create(1);
                beam.R[count] = x;
                beam.P[count] = p;
				beam.Q[count] = beam.getChargePerParticle();
                count++;
            }
        }

        if (!isBinned)
            {
                pc++;
                if (pc == Ippl::getNodes())
                    pc = 0;
            }
    }
    if (isBinned) {
        pbin_m->sortArray();
        // now copy this over to the bunch
        // so that we can emmit the particles
        beam.setPBins(pbin_m);
    }
}

/** 
 * 
 * 
 * @param beam 
 * @param Np 
 * @param scan 
 */
void Distribution::create(PartBunch &beam, size_t Np, bool scan) {

    if (beam.getTotalNum() != 0) {
        scan_m = scan;
	create(beam, beam.getLocalNum());
    }
    else {
	scan_m = false; // the first time we have to create particles
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

    if (isBinned) {
        if (pbin_m)
            delete pbin_m;
        pbin_m = new PartBins((int) Attributes::getReal(itsAttr[NBIN]));
        if (scan_m) {
	    beam.destroy(beam.getLocalNum(),0);
	    beam.update();
	    INFOMSG("In scan mode: deleted BIN structure and all particles in the bunch" << endl;);
	}
    }
    else {
        pbin_m = NULL;
    }
    
    beam.setTEmission(Attributes::getReal(itsAttr[TEMISSION]));
    beam.setNumBunch(1);
    const string disttype = Attributes::getString(itsAttr[DISTRIBUTION]);
    if (disttype == "GUNGAUSS" || disttype == "GUNUNIFORM" || disttype == "GUNGAUSS3D" || disttype == "GUNGAUSSFLATTOP" || disttype == "GUNGAUSSFLATTOPTH" || disttype == "GUNGAUSSFLATTOPTH-T")
        // Create a an initial beam bunch that is:
        // "GUNGAUSS": uniform in space transversely and with a Gaussian ("GUNGAUS") longitudinal profile
        // "GUNUNIFORM": uniform in space transversely and longitudinally.
        // "GUNGAUSS3D": Gaussian transversely and longitudinally.
        // "GUNGAUSSFLATTOP": uniform in space transversely, a Gaussian rise and fall time longitudinally with
        //                    a uniform flattop between.
        // "GUNGAUSSFLATTOPTH": uniform in space transversely, a Gaussian rise and fall time longitudinally with
        //                    a uniform flattop between, and a transvers thermal emittance
	
        {
	    if (disttype == "GUNGAUSSFLATTOPTH-T")
		binnDistributionT(beam, Np, disttype);
	    else
		binnDistributionZ(beam, Np, disttype);
        }
    else if (disttype == "BINOMIAL")
        {
	    
            Vector_t corr(Attributes::getReal(itsAttr[CORRX]),
                          Attributes::getReal(itsAttr[CORRY]),
                          Attributes::getReal(itsAttr[CORRT]));

            Vector_t sigX(Attributes::getReal(itsAttr[SIGMAX]),
                          Attributes::getReal(itsAttr[SIGMAY]),
                          Attributes::getReal(itsAttr[SIGMAT]));

            Vector_t sigPX(eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]),beam.getM()),
                           eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]),beam.getM()),
                           eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]),beam.getM()));

            Vector_t bincoef(Attributes::getReal(itsAttr[MX]),
                             Attributes::getReal(itsAttr[MY]),
                             Attributes::getReal(itsAttr[MT]));

            Vector_t emit;
            Vector_t alpha;
            Vector_t beta;
            Vector_t gamma;

            for(int j=0; j<3; j++) {
                double chi = asin(corr[j]);
                emit[j] = sigX[j]*sigPX[j]*cos(chi);
            }
            for(int j=0; j<3; j++) {
                beta[j]  = sigX[j]*sigX[j]/emit[j];
                gamma[j] = sigPX[j]*sigPX[j]/emit[j];
                alpha[j] = -corr[j]*sqrt(beta[j]*abs(gamma[j]));
            }
	    msg << "About to create Binomial distribution -1 " << endl;	    
            createBinom(emit, alpha, beta, gamma, bincoef, beam, Np, isBinned);
        }
    else if (disttype == "UNITUNIL")
        {

            Vector_t corr(Attributes::getReal(itsAttr[CORRX]),
                          Attributes::getReal(itsAttr[CORRY]),
                          Attributes::getReal(itsAttr[CORRT]));

            Vector_t sigX(Attributes::getReal(itsAttr[SIGMAX]),
                          Attributes::getReal(itsAttr[SIGMAY]),
                          Attributes::getReal(itsAttr[SIGMAT]));

            Vector_t sigPX(eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]),beam.getM()),
                           eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]),beam.getM()),
                           eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]),beam.getM()));

            Vector_t emit;
            Vector_t alpha;
            Vector_t beta;
            Vector_t gamma;

            for(int j=0; j<3; j++) {
                double chi = asin(corr[j]);
                emit[j] = sigX[j]*sigPX[j]*cos(chi);
            }
            for(int j=0; j<3; j++) {
                beta[j]  = sigX[j]*sigX[j]/emit[j];
                gamma[j] = sigPX[j]*sigPX[j]/emit[j];
                alpha[j] = -corr[j]*sqrt(beta[j]*abs(gamma[j]));
            }

            createUniformTUniformL(emit, alpha, beta, gamma, beam, Np, isBinned);
        }
    else if (disttype == "GAUSS")
        {
            double corr[7];
            corr[0] = Attributes::getReal(itsAttr[CORRX]);
            corr[1] = Attributes::getReal(itsAttr[CORRY]);
            corr[2] = Attributes::getReal(itsAttr[CORRT]);
            corr[3] = Attributes::getReal(itsAttr[R61]);
            corr[4] = Attributes::getReal(itsAttr[R62]);
            corr[5] = Attributes::getReal(itsAttr[R51]);
            corr[6] = Attributes::getReal(itsAttr[R52]);
            double Hs2a = Attributes::getReal(itsAttr[SIGMAX]);
            double Hs2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]),beam.getM()); //in eV

            double Vs2a = Attributes::getReal(itsAttr[SIGMAY]);
            double Vs2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]),beam.getM()); //in eV

            double Ls2a = Attributes::getReal(itsAttr[SIGMAT]);
            double Ls2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]),beam.getM()); //in eV

            double avrgpt = eVtoBetaGamma(Attributes::getReal(itsAttr[PT]),beam.getM());
            double avrgt = Attributes::getReal(itsAttr[T]);

            RANLIB_class *rGen = new RANLIB_class(265314159,4);

            unsigned int pc = 0;
            unsigned int count = 0;

            if (beam.isZPeriodic())
                INFOMSG("Distribution: uniform in z (periodic BC) z = " << -0.5*beam.getGaBeLa() << " ... " << 0.5*beam.getGaBeLa() << endl);
            for(int i=0; i<Np; i++) {
                double x,y;       // generate independent Gaussians, then correlate and finaly scale

                x  = rGen->gauss(0.0,1.0);
                y  = rGen->gauss(0.0,1.0);
                double xx =x;
                double yy =y;
                double px0  = x*corr[0] + y*sqrt(1.0 - corr[0]*corr[0]);
                double x0   =  x*Hs2a;
                px0 *= Hs2b;

                x  = rGen->gauss(0.0,1.0);
                y  = rGen->gauss(0.0,1.0);
                double py0  = x*corr[1] + y*sqrt(1.0 - corr[1]*corr[1]);
                double y0   =  x*Vs2a;
                py0 *= Vs2b;

                double del0;
                double psi0;

                if (beam.isZPeriodic()) {
                    /*
                      create uniform distribution
                    */
                    //del0 = (IpplRandom()-0.5)*p[2];
                    //	psi0 = (IpplRandom()*beam.getGaBeLa()) - (0.5*beam.getGaBeLa());

                } else {
                    x  = rGen->gauss(0.0,1.0);
                    y  = rGen->gauss(0.0,1.0);
                    double l32=(corr[6]-corr[0]*corr[5])/sqrt(1.0 - corr[0]*corr[0]);
                    double l33=sqrt(1-corr[5]*corr[5]-l32*l32);
                    psi0 = xx*corr[5]+yy*l32+x*l33;
                    double l42=(corr[4]-corr[0]*corr[3])/sqrt(1.0 - corr[0]*corr[0]);
                    double l43=(corr[2]-corr[5]*corr[3]-l42*l32)/l33;
                    double l44=sqrt(1-corr[3]*corr[3]-l42*l42-l43*l43);
                    del0 = xx*corr[3]+yy*l42+x*l43+y*l44;
                    //      	  del0 = xx*corr[3]+yy*0.268109486023689+y*0.693221684242576;
                    psi0 = avrgt + psi0*Ls2a;
                    del0 = avrgpt + Ls2b*del0;
                }
                if (pc == Ippl::myNode()) {
					if (!scan_m)
						beam.create(1);
                    beam.R[count] = Vector_t(x0,y0,psi0);
                    beam.P[count] = Vector_t(px0,py0,del0);
                    beam.Bin[count] = 0; // not initialized
                    //           dist.precision(8);
                    //           dist << x0 << "\t"
                    //                << y0 << "\t"
                    //                << psi0 << "\t"
                    //                << px0 << "\t"
                    //                << py0 << "\t";
                    //           dist.precision(15);
                    //           dist << del0 << endl;
                    count++;
                }
                pc++;
                if (pc == Ippl::getNodes())
                    pc = 0;
            }
        }
    else if (disttype == "UNIFORMXYZ")
        {

            double corr[3];

            corr[0] = Attributes::getReal(itsAttr[CORRX]);
            corr[1] = Attributes::getReal(itsAttr[CORRY]);
            corr[2] = Attributes::getReal(itsAttr[CORRT]);
	    
			double Hs2a = Attributes::getReal(itsAttr[SIGMAX]);
            double Hs2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]),beam.getM()); //in eV

            double Vs2a = Attributes::getReal(itsAttr[SIGMAY]);
            double Vs2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]),beam.getM()); //in eV

            double Ls2a = Attributes::getReal(itsAttr[SIGMAT]);
            double Ls2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]),beam.getM()); //in eV

            double nBins = Attributes::getReal(itsAttr[NBIN]);

            RANLIB_class *rGen = new RANLIB_class(265314159,4);

            unsigned int pc = 0;
            unsigned int count = 0;

            double gamma = 1 + (Ls2b / beam.getM());
            double beta  = sqrt(1- (1/(gamma*gamma)));
            double bega  = beta*gamma;

            for(int i=0; i<Np; i++) {

                double x,y,z,px,py,pz;       // generate independent Gaussians, then correlate and finaly scale
                double R =3;

                while(R >1){

                    x = rGen->uniform(-1.0,1.0);
                    y = rGen->uniform(-1.0,1.0);
                    z = rGen->uniform(-1.0,1.0);

                    px = rGen->uniform(-1.0,1.0);
                    py = rGen->uniform(-1.0,1.0);
                    pz = rGen->uniform(-1.0,1.0);

                    R = sqrt(x*x + y*y + z*z);

                    // or can generate uniform distribution in 6D phase space by
                    // using following line instead of above one.
                    // R = sqrt(x*x + y*y + z*z + px*px + py*py +pz*pz);
                }

                px  = x*corr[0] + px*sqrt(1.0 - corr[0]*corr[0]);
                x   = x*Hs2a;
                px *= Hs2b;

                py  = y*corr[1] + py*sqrt(1.0 - corr[1]*corr[1]);
                y   = y*Vs2a;
                py *= Vs2b;

                pz  = z*corr[2] + pz*sqrt(1.0 - corr[2]*corr[2]);
                z   = z*Ls2a;
                pz *= Ls2b;

                if (pc == Ippl::myNode()) {
                    beam.create(1);
                    beam.R[count] = Vector_t(x,y,z);
                    beam.P[count] = Vector_t(px,py,pz );
                    beam.Bin[count] = 0; // not initialized
                    count++;
                }
                pc++;
                if (pc == Ippl::getNodes())
                    pc = 0;
            }

        }
    else if (disttype == "FROMFILE") {


		/*

	    std::ofstream os;
	    if (Ippl::getNodes() == 1) {
		string file("dist.dat");
		os.open(file.c_str());
		if (os.bad()) {
		*gmsg << "Unable to open output file " <<  file << endl;
		}
		os << "# x px y py z pz " << endl;
	    }


		if (Ippl::getNodes() == 1) {
		os << x0 << "\t " << px0    << "\t "
		<< y0 << "\t " << py0    << "\t "
		<< psi0 << "\t " << del0 << "\t " << endl;
		}

		os.close();


		*/

		*gmsg <<"---------------------------------------------"<<endl;
		*gmsg <<"     READ ININITAL DISTRIBUTION FROM FILE    "<<endl;
		*gmsg <<"     BE AWARE OF THE FACT THAT ONLY NODE 0 IS READING IN "<<endl;
		if (isBinned) {
			*gmsg <<"     DISTRIBUTION will be binned using " << ebins << " energy bins " << endl;
			const string fn;
            binnDistributionFromFile(beam, fn);
		}
		else {
			double x0,px0,y0,py0,psi0,del0;
			if (Ippl::myNode() == 0) {
				const string filename=Attributes::getString(itsAttr[FNAME]);		
				std::ifstream fs;
				fs.open(filename.c_str());

				if(fs.fail()){
					throw OpalException("Distribution::Create()",
										"Open file operation failed, please check if \""
										+ filename +  "\" really exists.");
				}
		
				fs >> Np;
				if ( Np <= 0 ){
					throw OpalException("Distribution::Create()",
										" The particle number should be bigger than zero! Please check the first line of file \""
										+ filename +  "\".");
				}
		
				for(unsigned int i=0; i<Np; i++) {
					if ( !fs.eof() ){
						beam.create(1);
						fs >> x0 >> px0 >> y0 >> py0 >> psi0 >> del0;
						beam.R[i] = Vector_t(x0,y0,psi0);
						beam.P[i] = Vector_t(px0,py0,del0);
						beam.Bin[i] = 0; // not initialized
						beam.Q[i] = beam.getChargePerParticle();
					}
					else {
						throw OpalException("Distribution::Create()",
											"End of file reached before all particles imported, please check file \""
											+ filename +  "\".");
						return;
					}
				}
				fs.close();
			}
		}
		/*
		  In the case of a binned distribution (gun)
		  we have to do the boundp after emission.
		*/
    }

    if (!(isBinned)) {
		beam.boundp();
		beam.LastSection=0;
    }
}


/** 
 * 
 * 
 * @param beam 
 * @param Np 
 * @param distType 
 */
void Distribution::binnDistributionT(PartBunch &beam, size_t Np, string distType)
{
    const double &two_pi = Physics::two_pi;
    unsigned int pc = 0;

    double dEBins = Attributes::getReal(itsAttr[DEBIN]);
    pbin_m->setRebinEnergy(dEBins);

    double corr[3] = {Attributes::getReal(itsAttr[CORRX]),
					  Attributes::getReal(itsAttr[CORRY]),
					  Attributes::getReal(itsAttr[CORRT])};

    double nBins = Attributes::getReal(itsAttr[NBIN]);
    double transvCutOff = Attributes::getReal(itsAttr[TRANSVCUTOFF]);

    double Hs2a = Attributes::getReal(itsAttr[SIGMAX]);
    double Hs2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]),beam.getM());

    double Vs2a = Attributes::getReal(itsAttr[SIGMAY]);
    double Vs2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]),beam.getM());

    double Ls2a = Attributes::getReal(itsAttr[SIGMAT]);
    double Ls2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]),beam.getM());

    double tEmission = Attributes::getReal(itsAttr[TEMISSION]);
    double sigmaRise = Attributes::getReal(itsAttr[SIGMARISE]);
    double sigmaFall = Attributes::getReal(itsAttr[SIGMAFALL]);
    double flatTop = Attributes::getReal(itsAttr[FLATTOPTIME]);
    double cutoffRise = Attributes::getReal(itsAttr[CUTOFFRISE]);
    double cutoffFall = Attributes::getReal(itsAttr[CUTOFFFALL]);
    //double sigma = 3.0;

    double tBin = tEmission/nBins;
    tBin_m = tBin;

    RANLIB_class *rGen = new RANLIB_class(265314159,4);

    gsl_histogram * h = gsl_histogram_alloc (nBins);
    //gsl_histogram_set_ranges_uniform (h, 0, tEmission);

    //createTimeBins(Np, sigma, nBins, riseTime, fallTime, tEmission, h);
    createTimeBins(Np, sigmaRise, sigmaFall, nBins, cutoffRise, cutoffFall, flatTop, h); 

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

    ekin = Attributes::getReal(itsAttr[EKIN]);
    ptot = eVtoBetaGamma(ekin,beam.getM());
    
    // ASTRA mode
    phimax = Physics::pi/2.0;
    *gmsg << " -- B I N N I N G in T -----------------------------------------" << endl;
    *gmsg << " ---------------------I N P U T --------------------------------" << endl;
    *gmsg << " GUNGAUSS FLAT TOP &  THERMAL EMITTANCE in ASTRA MODE" << endl;
    *gmsg << " Kinetic energy = " << ekin << " [eV]  " << endl;
    *gmsg << " Phi max = " << phimax*180/Physics::pi << " [deg]  " << endl;
    *gmsg << " tBin = " << tBin << " [sec]  nBins = " << nBins << " tEmission =  " << tEmission << " [sec] " << endl;
    
    if (Ippl::getNodes() == 1) {
		*gmsg << " Write distribution to file dist.dat" << endl;
		string file("dist.dat");
		os.open(file.c_str());
		if (os.bad()) {
			*gmsg << "Unable to open output file " <<  file << endl;
		}
		os << "# x y ti px py pz phi theta Ekin= " << ekin << " [eV] " << endl;
    }

    for (int b=0; b<gsl_histogram_bins(h); b++) {
		/*
		  now many particles are in bin-number b?
	  
		*/
		*gmsg << "Fill bin " << b << " with n " << gsl_histogram_get(h,b) << " particles " 
			  << " myNode " << Ippl::myNode() << " getNodes " << Ippl::getNodes() << endl;
		pc = 0;
		for(int i=0; i<gsl_histogram_get(h,b); i++) {
	  
			double x,y;       // generate independent Gaussians, then correlate and finally scale
			double u1, u2;
			double xy = 6;


			while (xy > 1){
				x  = rGen->uniform(-1.0,1.0);
				y  = rGen->uniform(-1.0,1.0);
				xy = sqrt(x*x + y*y);
			}
	    
			double x0   =  x*Hs2a;
			double y0   =  y*Vs2a;

			/*
			  Now calculate the thermal emittances
			*/
	    
			const double phi   = 2.0 * acos(sqrt(rGen->uniform(0.0, 1.0)));
			const double theta = 2.0 * Physics::pi * rGen->uniform(0.0, 1.0);
			const double bega = 0.0;
			const double px0  = ptot * sin(phi) * cos(theta);
			const double py0  = ptot * sin(phi) * sin(theta);
			const double del0 = bega + (ptot * abs(cos(phi)));
	    
			if (pc == Ippl::myNode()) {
				vector<double> tmp;
				tmp.push_back(x0);
				tmp.push_back(y0);
				tmp.push_back(0.0);   
				tmp.push_back(px0);
				tmp.push_back(py0);
				tmp.push_back(del0);
				tmp.push_back((double)b);
				pbin_m->fill(tmp);
			}
	    
			pc++;
			if (pc == Ippl::getNodes())
				pc=0;

			if (Ippl::getNodes() == 1) {
				os << x0 << "\t " << px0    << "\t "
				   << y0 << "\t " << py0    << "\t "
				   << b << "\t " << del0 << "\t "
				   << phi*180./Physics::pi << "\t " << theta*180./Physics::pi << "\t "
				   << betaGammatoeV(px0, beam.getM())  << "\t "
				   << betaGammatoeV(py0, beam.getM())  << "\t "
				   << betaGammatoeV(del0, beam.getM()) << "\t " << endl;
			}

		}
    }
    if (Ippl::getNodes() == 1)
		os.close();
    
    pbin_m->setHistogram(h);
    pbin_m->sortArrayT();
    
    // now copy this over to the bunch
    // so that we can emit the particles
    beam.setPBins(pbin_m);

    *gmsg << " ---------------------------------------------------------------" << endl;
    *gmsg << " ----------- T - B I N N I N G  Done ---------------------------" << endl;		
    *gmsg << " ---------------------------------------------------------------" << endl;		


}

/** 
 * 
 * 
 * @param beam 
 * @param Np 
 * @param distType 
 */
void Distribution::binnDistributionZ(PartBunch &beam, size_t Np, string distType)
{
    const double &two_pi = Physics::two_pi;
    unsigned int pc = 0;

    double dEBins = Attributes::getReal(itsAttr[DEBIN]);
    pbin_m->setRebinEnergy(dEBins);

    double eneg= Attributes::getReal(itsAttr[PT]);
    double gamma = 1. + (eneg / beam.getM());
    double beta  = sqrt(1. - (1./(gamma*gamma)));
    double bega  = beta*gamma;

    pbin_m->setGamma(gamma); //has to be done by all processors!

    double corr[3] = {Attributes::getReal(itsAttr[CORRX]),
					  Attributes::getReal(itsAttr[CORRY]),
					  Attributes::getReal(itsAttr[CORRT])};

    double transvCutOff = Attributes::getReal(itsAttr[TRANSVCUTOFF]);
    double Hs2a = Attributes::getReal(itsAttr[SIGMAX]);
    double Hs2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPX]),beam.getM());

    double Vs2a = Attributes::getReal(itsAttr[SIGMAY]);
    double Vs2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPY]),beam.getM());

    double Ls2a = Attributes::getReal(itsAttr[SIGMAT]);
    double Ls2b = eVtoBetaGamma(Attributes::getReal(itsAttr[SIGMAPT]),beam.getM());

    double riseTime = Attributes::getReal(itsAttr[SIGMARISE]);
    double fallTime = Attributes::getReal(itsAttr[SIGMAFALL]);
    double flatTop = Attributes::getReal(itsAttr[FLATTOPTIME]);

    int particlesFront = Np * sqrt(2.0 * Physics::pi) * (riseTime / 2.0)
		/ (flatTop + sqrt(2.0 * Physics::pi) * (riseTime + fallTime) / 2.0);
    int particlesTail = Np * sqrt(2.0 * Physics::pi) * (fallTime / 2.0)
		/ (flatTop + sqrt(2.0 * Physics::pi) * (riseTime + fallTime) / 2.0);

    double nBins = Attributes::getReal(itsAttr[NBIN]);

    RANLIB_class *rGen = new RANLIB_class(265314159,4);

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

    if (distType == "GUNGAUSSFLATTOPTH") {

		ekin = Attributes::getReal(itsAttr[EKIN]);
		if ((astraMode=(ekin > 0.0))) {
			// ASTRA mode
			phimax = Physics::pi/2.0;
			*gmsg << " ---------------------------------------------------------------" << endl;
			*gmsg << " ---------------------I N P U T --------------------------------" << endl;
			*gmsg << " GUNGAUSS FLAT TOP &  THERMAL EMITTANCE in ASTRA MODE" << endl;
			*gmsg << " Kinetic energy = " << ekin << " [eV]  " << endl;
			*gmsg << " Phi max = " << phimax*180/Physics::pi << " [deg]  " << endl;
			*gmsg << " ---------------------------------------------------------------" << endl;		

			ptot = eVtoBetaGamma(ekin,beam.getM());

		}
		else {
			workf = Attributes::getReal(itsAttr[W]);
			siglaser = Attributes::getReal(itsAttr[SIGLASER]);
			elaser = Attributes::getReal(itsAttr[ELASER]);
			fe = Attributes::getReal(itsAttr[FE]);
			ag = Attributes::getReal(itsAttr[AG]) * 1000000;

			ekin = fe + elaser;

			schottky = sqrt(Physics::q_e*ag/4/Physics::pi/Physics::epsilon_0);

			workf -= schottky;

			phimax = acos(sqrt((fe+workf)/ekin));

			ptot = eVtoBetaGamma((workf+elaser-fe-schottky),beam.getM()),

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
			*gmsg << " Phi max = " << phimax*180/Physics::pi << " [deg]  " << endl;
			*gmsg << " Ptot = " << workf+elaser-fe-schottky << " [eV]  " << endl;
		
			*gmsg << " ---------------------------------------------------------------" << endl;
		}
    }

    if (Ippl::getNodes() == 1) {
		*gmsg << " Write distribution to file dist.dat" << endl;
		string file("dist.dat");
		os.open(file.c_str());
		if (os.bad()) {
			*gmsg << "Unable to open output file " <<  file << endl;
		}
		os << "# x px y py z pz phi theta Ethx [eV] Ethy [eV] Ethz [eV]" << endl;
    }


    for(int i=0; i<Np; i++) {
		double x,y;       // generate independent Gaussians, then correlate and finally scale
		double u1, u2;
		double xy = 6;

		if (distType != "GUNGAUSS3D") {
			while (xy > 1){
				x  = rGen->uniform(-1.0,1.0);
				y  = rGen->uniform(-1.0,1.0);
				xy = sqrt(x*x + y*y);
			}
		}
	    
		else {
			while (xy > transvCutOff || isnan(xy)) {
				u1 = rGen->uniform(0.0,1.0);
				u2 = rGen->uniform(0.0,1.0);
				x = sqrt(-2. * log(1. - u1)) * cos(two_pi * u2);
				y = sqrt(-2. * log(1. - u1)) * sin(two_pi * u2);
				xy = sqrt(x*x + y*y);
			}
		}

		double x0   =  x*Hs2a;
		double y0   =  y*Vs2a;

		double del0 = 0.0;
		double psi0 = 0.0;

		if (distType == "GUNGAUSS" || distType == "GUNGAUSS3D") {
			x  = rGen->gauss(0.0,1.0);
			y  = rGen->gauss(0.0,1.0);
		} else if (distType == "GUNUNIFORM") {
			x = rGen->uniform(0.0, 1.0);
			y = rGen->uniform(0.0, 1.0);
		} else if (distType == "GUNGAUSSFLATTOP" || distType == "GUNGAUSSFLATTOPTH") {
			if (i < particlesFront) {
				// Fill rise.
				psi0 = -1.0;
				while (psi0 < 0.0)
					psi0 = rGen->gauss(0.0, 1.0);
				psi0 = psi0 * riseTime + flatTop;
			} else if (i >= particlesFront && i < particlesFront + particlesTail) {
				// Fill fall.
				psi0 = 1.0;
				while (psi0 > 0.0)
					psi0 = rGen->gauss(0.0, 1.0);
				psi0 = psi0 * fallTime;
			} else {
				// Fill flat top.
				psi0 = rGen->uniform(0.0, 1.0);
				psi0 *= flatTop;
			}
		}

		if (distType != "GUNGAUSSFLATTOP" && distType != "GUNGAUSSFLATTOPTH") {
			del0  = x*corr[2] + y*sqrt(1.0 - corr[2]*corr[2]);
			psi0  = x*Ls2a;
		}

		/*
		  This is the longitudinal momenta in units of beta gamma.
		*/

		del0  = bega + Ls2b*del0;

		double px0 = 0.0;
		double py0 = 0.0;
		double phi = 0.0;
		double theta = 0.0;

		if (distType == "GUNGAUSSFLATTOPTH") {

			/*
			  Now calculate the thermal emittances
			*/

			phi   = 2.0 * acos(sqrt(rGen->uniform(0.0, 1.0)));
			theta = 2.0 * Physics::pi * rGen->uniform(0.0, 1.0);

			px0  = ptot * sin(phi) * cos(theta);
			py0  = ptot * sin(phi) * sin(theta);
			del0 = bega + (ptot * abs(cos(phi)));
		}

		if (pc == Ippl::myNode()) {
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
		if (pc == Ippl::getNodes())
			pc=0;

		if (Ippl::getNodes() == 1) {
			os << x0 << "\t " << px0    << "\t "
			   << y0 << "\t " << py0    << "\t "
			   << psi0 << "\t " << del0 << "\t "
			   << phi*180./Physics::pi << "\t " << theta*180./Physics::pi << "\t "
			   << betaGammatoeV(px0, beam.getM())  << "\t "
			   << betaGammatoeV(py0, beam.getM())  << "\t "
			   << betaGammatoeV(del0, beam.getM()) << "\t " << endl;
		}

    }
    if (Ippl::getNodes() == 1)
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
void Distribution::binnDistributionFromFile(PartBunch &beam, const string fn)
{
    unsigned int pc = 0;
    size_t Np;
    double x0,y0,psi0,px0,py0,del0;

    std::ofstream os;
    std::ifstream fs;

    double dEBins = Attributes::getReal(itsAttr[DEBIN]);
    pbin_m->setRebinEnergy(dEBins);

    double eneg= Attributes::getReal(itsAttr[PT]);
    double gamma = 1. + (eneg / beam.getM());
    double beta  = sqrt(1. - (1./(gamma*gamma)));
    double bega  = beta*gamma;

    pbin_m->setGamma(gamma); //has to be done by all processors!

    double nBins = Attributes::getReal(itsAttr[NBIN]);

    if (Ippl::getNodes() == 1) {
		*gmsg << " Write distribution to file dist.dat" << endl;
		string file("dist.dat");
		os.open(file.c_str());
		if (os.bad()) {
			*gmsg << "Unable to open output file " <<  file << endl;
		}
		os << "# x px y py z pz phi theta Ethx [eV] Ethy [eV] Ethz [eV]" << endl;
    }

    fs.open(fn.c_str());

    if(fs.fail()){
		throw OpalException("Distribution::Create()",
							"Open file operation failed, please check if \""
							+ fn +  "\" really exists.");
    }
	    
    fs >> Np;
	
    if ( Np <= 0 ){
		throw OpalException("Distribution::Create()",
							" The particle number should be bigger than zero! Please check the first line of file \""
							+ fn +  "\".");
    }
	
    for(unsigned int i=0; i<Np; i++) {
		if ( !fs.eof() ){
			fs>> x0 >> px0 >> y0 >> py0 >> psi0 >> del0;
			if (pc == Ippl::myNode()) {
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
			if (pc == Ippl::getNodes())
				pc=0;
		
			if (Ippl::getNodes() == 1) {
				os << x0 << "\t " << px0    << "\t "
				   << y0 << "\t " << py0    << "\t "
				   << psi0 << "\t " << del0 << endl;
			}
		}
    }

    if (Ippl::getNodes() == 1)
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
void Distribution::doRestart(PartBunch &beam, size_t Np, size_t restartStep)
{
    H5PartFile *H5file;
    string fn;

    IpplTimings::startTimer(beam.distrReload_m);

    if (OPAL.hasRestartFile()){
		fn=OPAL.getRestartFileName();
		*gmsg<<"Restart from a specified file:"<<fn<<endl;

    }
    else{
		fn= OPAL.getInputFn();
		int pos=fn.find(string("."),0);
		fn.erase(pos,fn.size()-pos);
		//        beam.setTEmission(Attributes::getReal(itsAttr[TEMISSION]));
		fn += string(".h5");
    }

#ifdef PARALLEL_IO
    H5file=H5PartOpenFileParallel(fn.c_str(),H5PART_READ,MPI_COMM_WORLD);
#else
    H5file=H5PartOpenFile(fn.c_str(),H5PART_READ);
#endif

    if(!H5file) {
		ERRORMSG("could not open file '" << fn << "';  exiting!" << endl);
		exit(0);
    }
    if (restartStep == -1) {
		restartStep = H5PartGetNumSteps(H5file) - 1;
		OPAL.setRestartStep(restartStep);
    } else {
		if (restartStep != H5PartGetNumSteps(H5file) - 1 && !OPAL.hasRestartFile()) {
			ERRORMSG("can't append to the file '" << fn << "' exiting!" << endl);
			exit(0);
		}
    }


    *gmsg << "restart step = " << restartStep << endl;
    H5PartSetStep(H5file,restartStep);
    int N=(int)H5PartGetNumParticles(H5file);

    h5part_int64_t totalSteps = H5PartGetNumSteps(H5file);
    *gmsg << "total number of particles = " << N << " total steps " << totalSteps << endl;

    //TODO: do a more sophisticated distribution of particles?
    //my guess is that the end range index is EXCLUSIVE!
    int numberOfParticlesPerNode = (int) floor((double) N / Ippl::getNodes());
    long long starti = Ippl::myNode() * numberOfParticlesPerNode;
    long long endi = 0;
    // ensure that we don't miss any particle in the end
    if(Ippl::myNode() == Ippl::getNodes() - 1)
		endi = -1;
    else
		endi = starti + numberOfParticlesPerNode;

    // cout << "Node " << Ippl::myNode() << " gets particles " << starti << " - " << endi << endl;

    H5PartSetView(H5file,starti,endi);
    N = (int)H5PartGetNumParticles(H5file);

    // cout << "Node " << Ippl::myNode() << " has " << N << " particles" << endl;

    double actualT;
    H5PartReadStepAttrib(H5file,"TIME",&actualT);

    beam.setT(actualT);

    /*
     * eventually read more attributes a la IMPACT-T
     */

    void *varray = malloc(N*sizeof(double));
    double *farray = (double*)varray;

    beam.create(N);

    H5PartReadDataFloat64(H5file,"x",farray);
    for (unsigned long int n=0; n < N; ++n) {
		beam.R[n](0)=farray[n];
		beam.Bin[n] = 0; // not initialized
    }
    H5PartReadDataFloat64(H5file,"y",farray);
    for (unsigned long int n=0; n < N; ++n)
		beam.R[n](1)=farray[n];

    H5PartReadDataFloat64(H5file,"z",farray);
    for (unsigned long int n=0; n < N; ++n)
		beam.R[n](2)=farray[n];

    H5PartReadDataFloat64(H5file,"px",farray);
    for (unsigned long int n=0; n < N; ++n)
		beam.P[n](0)=farray[n];

    H5PartReadDataFloat64(H5file,"py",farray);
    for (unsigned long int n=0; n < N; ++n)
		beam.P[n](1)=farray[n];

    H5PartReadDataFloat64(H5file,"pz",farray);
    for (unsigned long int n=0; n < N; ++n)
		beam.P[n](2)=farray[n];

    if(farray)
		free(farray);

    Ippl::Comm->barrier();
    H5PartCloseFile(H5file);
    beam.boundp();
    beam.LastSection=0;
	beam.Q = beam.getChargePerParticle();
    IpplTimings::stopTimer(beam.distrReload_m);
}

/** 
 * 
 * 
 * @param beam 
 * @param Np 
 * @param restartStep 
 * @param specifiedNumBunch 
 */
void Distribution::doRestart_cycl(PartBunch &beam, size_t Np, size_t restartStep, const int specifiedNumBunch)
{
    IpplTimings::startTimer(beam.distrReload_m);
    *gmsg<<"---------------- Start reading hdf5 file----------------"<<endl;
    H5PartFile *H5file;

    string fn;
    if (OPAL.hasRestartFile()){

		fn=OPAL.getRestartFileName();
		*gmsg<<"Restart from a specified file:"<<fn<<endl;
    }
    else{
		fn= OPAL.getInputFn();
		int pos=fn.find(string("."),0);
		fn.erase(pos,fn.size()-pos);

		beam.setTEmission(Attributes::getReal(itsAttr[TEMISSION]));

		fn += string(".h5");
    }

    *gmsg<<"Restart from hdf5 format file "<<fn<<", read phase space data of DumpStep "<<restartStep<<endl;

#ifdef PARALLEL_IO
    H5file=H5PartOpenFileParallel(fn.c_str(),H5PART_READ,MPI_COMM_WORLD);
#else
    H5file=H5PartOpenFile(fn.c_str(),H5PART_READ);
#endif

    if(!H5file) {
		ERRORMSG("File open failed:  exiting!" << endl);
		exit(0);
    }

    H5PartSetStep(H5file,restartStep);
    const int globalN=(int)H5PartGetNumParticles(H5file);

    h5part_int64_t totalSteps = H5PartGetNumSteps(H5file);
    *gmsg << "total number of particles = " <<globalN <<endl;

    int numberOfParticlesPerNode = (int) floor((double) globalN / Ippl::getNodes());
    long long starti = Ippl::myNode() * numberOfParticlesPerNode;
    long long endi = 0;

    if(Ippl::myNode() == Ippl::getNodes() - 1)
		endi = -1;
    else
		endi = starti + numberOfParticlesPerNode;

    H5PartSetView(H5file,starti,endi);
    const int localN = (int)H5PartGetNumParticles(H5file);

    // debug
    // Inform *gmsgAll;
    // gmsgAll = new  Inform("Message",INFORM_ALL_NODES);
    // *gmsgAll<< "total number of particles on this node = " <<localN <<endl;
    // end debug

    double actualT;
    H5PartReadStepAttrib(H5file,"TIME",&actualT);

    beam.setT(actualT);

    double lpath;
    H5PartReadStepAttrib(H5file,"LPATH",&lpath);
    beam.setLPath(lpath);
    
    h5part_int64_t tstep;
    H5PartReadStepAttrib(H5file,"TrackStep",&tstep);
    beam.setTrackStep((long long)tstep);
    
    h5part_int64_t SteptoLastInj;
    H5PartReadStepAttrib(H5file,"SteptoLastInj",&SteptoLastInj);
    beam.setSteptoLastInj((int)SteptoLastInj);
    *gmsg<<"Tracking Step since last bunch injection is "<<SteptoLastInj<<endl;
    
    h5part_int64_t numBunch;
    H5PartReadStepAttrib(H5file,"NumBunch",&numBunch);
    beam.setNumBunch((int)numBunch);
    *gmsg<<"There are "<<numBunch<<" Bunches(bins) exist in this file"<<endl;

    double gammaBin[numBunch];

    H5PartReadStepAttrib(H5file,"GammaBin",&gammaBin);

    void *varray = malloc(localN*sizeof(double));
    double *farray = (double*)varray;

    h5part_int64_t *larray = (h5part_int64_t *)varray;

    beam.create(localN);

    H5PartReadDataFloat64(H5file,"x",farray);
    for (unsigned long int n=0; n < localN; ++n)
		beam.R[n](0)=farray[n];

    H5PartReadDataFloat64(H5file,"y",farray);
    for (unsigned long int n=0; n < localN; ++n)
		beam.R[n](1)=farray[n];

    H5PartReadDataFloat64(H5file,"z",farray);
    for (unsigned long int n=0; n < localN; ++n)
		beam.R[n](2)=farray[n];

    H5PartReadDataFloat64(H5file,"px",farray);
    for (unsigned long int n=0; n < localN; ++n)
		beam.P[n](0)=farray[n];

    H5PartReadDataFloat64(H5file,"py",farray);
    for (unsigned long int n=0; n < localN; ++n)
		beam.P[n](1)=farray[n];

    H5PartReadDataFloat64(H5file,"pz",farray);
    for (unsigned long int n=0; n < localN; ++n)
		beam.P[n](2)=farray[n];

    H5PartReadDataInt64(H5file,"id",larray);
    for (unsigned long int n=0; n < localN; ++n)
		beam.ID[n]=larray[n];

    // only for multi-bunch mode
    if ( specifiedNumBunch > 1) {
		/*
		  size_t partInBin[numBunch];
		  for(int ii=0; ii<numBunch; ii++) partInBin[ii] = 0 ;
	      
		  // assign bin index for each particle and
		  // calculate total particles number for each bin
		  for (unsigned long int n=0; n < localN; ++n)
		  {
		  double deltgamma[numBunch];
		  double gamma = sqrt(1.0 + dot(beam.P[n], beam.P[n]));
		  int index = 0;
	      
		  for(int ii=0; ii<numBunch; ii++)
		  deltgamma[ii] = abs( gammaBin[ii] - gamma );
	      
		  for(int ii=0; ii<numBunch; ii++)
		  if( *(deltgamma+index) > *(deltgamma+ii) )
		  index = ii;
	      
		  beam.Bin[n]=index;
		  partInBin[index]++;
		  }
	      
		  for(int ii=0; ii<numBunch; ii++)
		  reduce(partInBin[ii],partInBin[ii],OpAddAssign());
	      
		  // instantiate PartBins class for restart run
		  beam.setPBins( new PartBins(specifiedNumBunch, (int)numBunch, partInBin ));
		  }
		  // instantiate PartBins class for restart run
		  beam.setPBins( new PartBins(specifiedNumBunch, (int)numBunch, partInBin ));
		*/
	    
		// the allowed maximal bin number is set to 100
		beam.setPBins( new PartBins(100));
    }
	
    if(farray) free(farray);
	
    Ippl::Comm->barrier();
    H5PartCloseFile(H5file);
    beam.boundp();
	beam.Q = beam.getChargePerParticle();
    *gmsg<<"----------------Finish reading hdf5 file----------------"<<endl;
    IpplTimings::stopTimer(beam.distrReload_m);
}


/** 
 * 
 * 
 * @param name 
 * 
 * @return 
 */
Distribution *Distribution::find(const string &name)
{
    Distribution *beam = dynamic_cast<Distribution*>(OPAL.find(name));

    if (beam == 0) {
		throw OpalException("Distribution::find()", "Distribution \"" + name + "\" not found.");
    }

    return beam;
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
const PartData &Distribution::getReference() const
{
    // Cast away const, to allow logically constant Distribution to update.
    const_cast<Distribution *>(this)->update();
    return reference;
}

/** 
 * 
 * 
 */
void Distribution::update()
{

}

/** 
 * 
 * 
 * @param os 
 */
void Distribution::tfsDescriptors(std::ostream &os) const
{
    os << "@ Distribution     %s  " << getOpalName() << '\n' ;
}

/** 
 * 
 * 
 * @param os 
 * 
 * @return 
 */
Inform &Distribution::print(Inform &os) const {
    double transvCutOff = Attributes::getReal(itsAttr[TRANSVCUTOFF]);
    os << "* ************* D I S T R I B U T I O N ******************************************** " << endl;
    if (!OPAL.inRestartRun()) {
        os << "* Distribution:\t" << getOpalName() << endl;
        if (!(Attributes::getString(itsAttr[LASERPROFFN])==string(""))) {
            os << "* Distribution type:\t" << Attributes::getString(itsAttr[DISTRIBUTION]) << endl;
            os << "* Laser profile: " << Attributes::getString(itsAttr[LASERPROFFN]) 
               << " Image: " << Attributes::getString(itsAttr[IMAGENAME])
               << " Intensity cut: " << Attributes::getReal(itsAttr[INTENSITYCUT]) << endl;
        }
        else {
            os << "* Distribution type:\t" << Attributes::getString(itsAttr[DISTRIBUTION]) << endl;
        }
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
    else
        os << "* Distribution from restart file" << endl;
    /*
	  switch (distT_m) {

	  case string("GUNGAUSSFLATTOPTH-T"):
	  os << "Selected Distribution " << distT_m << endl;
	  break;
      default:
	  os << "Selected Distribution " << distT_m << " not known" << endl;
	  }
    */
    os << "* ********************************************************************************** " << endl;
}

