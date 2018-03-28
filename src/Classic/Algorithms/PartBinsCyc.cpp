#include "Algorithms/PartBinsCyc.h"
#include "Physics/Physics.h"
#include <cfloat>
#include <vector>
extern Inform *gmsg;

// constructer function for cyclotron
PartBinsCyc::PartBinsCyc(int specifiedNumBins, int bins, size_t  partInBin[])
    : PartBins(specifiedNumBins, 0) {

    bins_m = specifiedNumBins;        // max bin number
    nemittedBins_m = bins;            // the bin number with particles
    
    for(int i = 0; i < nemittedBins_m; i++) {
        nBin_m[i] = partInBin[i];

        *gmsg << "Read in: Bin=" << i << " Particles Num=" << getTotalNumPerBin(i) << endl;
        binsEmitted_m[i] = true;
    }
}

// constructer function for cyclotron for restart run.
PartBinsCyc::PartBinsCyc(int specifiedNumBins, int bins)
    : PartBins(specifiedNumBins, 0) {

    bins_m = specifiedNumBins;        // max bin number
    nemittedBins_m = bins;            // the bin number with particles
    
    for(int i = 0; i < nemittedBins_m; i++) {
      binsEmitted_m[i] = true;
    }
}
