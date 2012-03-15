#include "disttest.h"
#include "Distribution.h"
#include "myPartBunch.h"


myPartBunch::myPartBunch(string distname, double ekin, Distribution *d): 
    dist_m(d),
    dt_m(0.0) 
{
    R.clear();
    P.clear();
    Q.clear();
    Bin.clear();
    PType.clear();
    dt.clear();
    pbin_m = NULL;
    
    if(Ippl::getNodes() == 1) {
	string fn = distname + string(".dist");
	os_m.open(fn.c_str());
	if(os_m.bad()) {
	    ERRORMSG("Unable to open output file " << fn << endl);
	}
	os_m << "# x y ti px py pz Ekin= " << ekin << " [eV] " << endl;
    }
    
}

size_t myPartBunch::emitParticles() {
    
    Inform msg ("emitParticles() ");

    int          binNumber  = pbin_m->getBinToEmit();
    const size_t partsInBin = pbin_m->getLocalBinCount(binNumber);

    if(binNumber != -1) {
	size_t numberOfEmittedParticles = this->getLocalNum();
	size_t oldNumberOfEmittedParticles = numberOfEmittedParticles;
	/*
	  now generate paticles in the interval t ... t+dt
	*/
	for(size_t particleNumber = 0; particleNumber < partsInBin; particleNumber++) {
	    
	    // sample the specifyed distribution
	    pair<Vector_t,Vector_t> xp = dist_m->sample(binNumber);

	    Vector_t x = xp.first;
	    Vector_t p = xp.second;

	    this->dt.push_back(x(2));

	    double oneOverCDt = 1.0 ; // / (Physics::c * this->dt[numberOfEmittedParticles]);
	    double particleGamma = sqrt(1.0 + pow(p(0), 2.0) + pow(p(1), 2.0) + pow(p(2), 2.0));
                
	    R.push_back(Vector_t(      oneOverCDt * (x(0) + p(0) * this->dt[numberOfEmittedParticles] * Physics::c / (2.0 * particleGamma)),
	    			       oneOverCDt * (x(1) + p(1) * this->dt[numberOfEmittedParticles] * Physics::c / (2.0 * particleGamma)),
				       x(2)));



			//				       oneOverCDt * (0.0  + p(2) * this->dt[numberOfEmittedParticles] * Physics::c / (2.0 * particleGamma))));

	    P.push_back(p);
	    Bin.push_back(binNumber);
	    Q.push_back(getChargePerParticle());
	    PType.push_back(0);

	    numberOfEmittedParticles++;                        /// Iterate number of particles that have been emitted.
	    binemitted_m[binNumber]++;                         /// Iterate number of particles in this bin that have been emitted.
	    
	}

	pbin_m->setBinEmitted(binNumber);
	
	msg << "* Bin number: " << binNumber << " has emitted all " << numberOfEmittedParticles - oldNumberOfEmittedParticles << " particles " << endl;
	
	/// Return number of particles emitted.
	return numberOfEmittedParticles - oldNumberOfEmittedParticles;
    } else
	return 0;
}
