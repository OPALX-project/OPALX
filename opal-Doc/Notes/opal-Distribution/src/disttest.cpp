/* Description
 * 
 *  This application is a testbed for all OPAL distributions
 *
 *  \author Andreas Adelmann
 * 
 *  \date 2010 nov 14
 * 
 *  \warning none
 * 
 *  \attention none required
 * 
 *  \bug this code is under permanent development.
 * 
 *  \todo ??
 */
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

#include "Algorithms/PartBins.h"

#include "Distribution.h"
#include "myPartBunch.h"
#include "disttest.h"

using namespace std;

int main(int argc, char *argv[]) {
    Ippl(argc, argv);
    Inform msg ("DistTest ");

    Vector_t rmsx(1.0E-3,2.0E-3,3.0E-3);

    Vector_t rmsp(1.0E-4,2.0E-4,3.0E-4);
	
    Vektor<double,7> corr(0.0);

    // for the binomial distribution define m 
    Vector_t m(1E6,0.5,9999.0);


    size_t Np           = 1E5; 
    double dEBins       = 1; 
    int nBins           = 0; 
    Vector_t distCutOff = Vector_t(3.0,3.0,3.0); 

    double tpulsefwhm   = 10.0E-12; 
    double trise        =  0.5E-12;
    double tfall        =  0.5E-12;
    double ekin         = 0.65;
    double avrgt        = 0; 
    double avrgpt       = 0; 

    bool hasLaserProfile  = false;
    string laserProfileFn = string("");
    string laserImage     = string("");
    double intensityCut   = 0;

    DistrTypeT actDist;
    string distfn;

    actDist =  GAUSS;
    distfn  = string("gauss");


    Distribution *dist = new Distribution(actDist, Np, dEBins, nBins, distCutOff,
					  hasLaserProfile, laserProfileFn, laserImage, intensityCut,
					  tpulsefwhm, trise, tfall, m,
					  ekin,
					  corr,
					  rmsx,
					  rmsp,
					  avrgt,
					  avrgpt); 

    myPartBunch *bunch = new myPartBunch(distfn, ekin, dist);

    msg << *dist << endl;
    dist->create(*bunch, Np);
    delete bunch;
    delete dist;



    actDist =  BINOMIAL;
    distfn  = string("binomial");

    dist = new Distribution(actDist, Np, dEBins, nBins, distCutOff,
					  hasLaserProfile, laserProfileFn, laserImage, intensityCut,
					  tpulsefwhm, trise, tfall, m,
					  ekin,
					  corr,
					  rmsx,
					  rmsp,
					  avrgt,
					  avrgpt); 

    bunch = new myPartBunch(distfn, ekin, dist);

    msg << *dist << endl;
    dist->create(*bunch, Np);
    delete bunch;
    delete dist;

    actDist   = GUNGAUSSFLATTOPTH;
    distfn    = string("gungaussft");
    nBins = 20; 
    dist = new Distribution(actDist, Np, dEBins, nBins, distCutOff,
			    hasLaserProfile, laserProfileFn, laserImage, intensityCut,
			    tpulsefwhm, trise, tfall, m,
			    ekin,
			    corr,
			    rmsx,
			    rmsp,
			    avrgt,
			    avrgpt); 

    bunch = new myPartBunch(distfn, ekin, dist);
    
    msg << *dist << endl;
    dist->create(*bunch, Np);

    size_t n  = 0;
    double t = 0.0;
    double dt = dist->getTBin();

    for (int step=0; step<nBins; step++){
	n += bunch->emitParticles();
    }
    msg << "Emission done particles emitted " << n << endl;
    delete bunch;
    delete dist;
    return 0;
}
