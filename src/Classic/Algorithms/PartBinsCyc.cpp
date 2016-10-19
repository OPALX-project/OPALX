#include "Algorithms/PartBinsCyc.h"
#include "Physics/Physics.h"
#include <cfloat>
#include <vector>
extern Inform *gmsg;

// constructer function for cyclotron
PartBinsCyc::PartBinsCyc(int specifiedNumBins, int bins, size_t  partInBin[])
    : PartBins::PartBins(specifiedNumBins, 0) {

    bins_m = specifiedNumBins;        // max bin number
    nemittedBins_m = bins;            // the bin number with particles
    binsEmitted_m = new bool[bins_m]; // flag whether the bin contain particles or not
    nBin_m = new size_t[bins_m];      // number of particles in the bins on the local node

    for(int i = 0; i < bins_m; i++) {
        nBin_m[i] = 0;
        binsEmitted_m[i] = false;
    }


    for(int i = 0; i < nemittedBins_m; i++) {
        nBin_m[i] = partInBin[i];

        *gmsg << "Read in: Bin=" << i << " Particles Num=" << getTotalNumPerBin(i) << endl;
        binsEmitted_m[i] = true;
    }

}

// constructer function for cyclotron for restart run.
PartBinsCyc::PartBinsCyc(int specifiedNumBins, int bins)
    : PartBins::PartBins(specifiedNumBins, 0) {

    bins_m = specifiedNumBins;        // max bin number
    nemittedBins_m = bins;            // the bin number with particles
    binsEmitted_m = new bool[bins_m]; // flag whether the bin contain particles or not
    nBin_m = new size_t[bins_m];      // number of particles in the bins on the local node

    for(int i = 0; i < bins_m; i++) {
        nBin_m[i] = 0;
        binsEmitted_m[i] = false;
    }

    for(int i = 0; i < nemittedBins_m; i++) {
      binsEmitted_m[i] = true;
    }
}

PartBinsCyc::~PartBinsCyc() {
    if(nBin_m) {
        delete nBin_m;
        delete xbinmax_m;
        delete xbinmin_m;
        delete binsEmitted_m;
    }
    tmppart_m.clear();
    isEmitted_m.clear();
    if(h_m)
        delete h_m;
}
