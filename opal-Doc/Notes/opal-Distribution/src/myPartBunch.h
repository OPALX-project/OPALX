#ifndef myPartBunch_HH
#define myPartBunch_HH

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

#include "Distribution.h"

class myPartBunch {

public:

    myPartBunch(string distname, double ekin, Distribution *d);
    
    ~myPartBunch() { 
	if((Ippl::getNodes() == 1) && (os_m.is_open())) {
	    for (size_t i=0; i<R.size(); i++) {
		os_m << R[i](0) <<"\t"
		     << R[i](1) <<"\t"
		     << R[i](2) <<"\t"
		     << P[i](0) <<"\t"
		     << P[i](1) <<"\t"
		     << P[i](2) << endl;
	    }
	    os_m.close();
	}
    }
    
    double getChargePerParticle() {return 1.0;}
    double getM() { return 1.0;}
    
    size_t getTotalNum() {return R.size();}
    size_t getLocalNum() {return R.size();}

    void   setPBins(PartBins *pbin) {
	pbin_m = pbin;
	//	bingamma_m = new double[pbin_m->getNBins()];
	binemitted_m = new size_t[pbin_m->getNBins()];
	for(int i = 0; i < pbin_m->getNBins(); i++)
	    binemitted_m[i] = 0;
	//	pbin_m->sortArray();

    }

    bool weHaveBins() const {
        if(pbin_m)
            return pbin_m->weHaveBins();
        else
            return false;
    }

    double getdT() { return dt_m;}
    double setdT(double dt) { dt_m = dt;}

    size_t emitParticles();

    vector<Vector_t> R;
    vector<Vector_t> P;
    vector<double> dt;
    vector<int> Bin;
    vector<double> Q;
    vector<int> PType;

private:
    PartBins *pbin_m;
    size_t *binemitted_m;

    std::ofstream os_m;
    double dt_m;
    Distribution *dist_m;    
};
#endif
