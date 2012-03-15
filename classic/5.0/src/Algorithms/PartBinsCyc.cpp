#include <cfloat>
#include <vector>
#include "Algorithms/PartBinsCyc.h"
//#include "Algorithms/PartBins.h"
#include "Physics/Physics.h"
extern Inform *gmsg;

// constructer function for cyclotron
PartBinsCyc::PartBinsCyc(int specifiedNumBins, int bins, size_t  partInBin[])
    : PartBins::PartBins(bins) {
    bins_m = specifiedNumBins;
    nemittedBins_m = bins;
    binsEmitted_m = new bool[bins_m];
    nBin_m = new size_t[bins_m];
    nDelBin_m = new size_t[bins_m];


    //dummy arrays, not used
    xbinmin_m = new double;
    xbinmax_m = new double;


    for(int i = 0; i < bins_m; i++) {
        nDelBin_m[i] = 0;
        nBin_m[i] = 0;
        binsEmitted_m[i] = false;
    }

    for(int i = 0; i < nemittedBins_m; i++) {
        nBin_m[i] = partInBin[i];
        *gmsg << "Read in: Bin=" << i << " Particles Num=" << nBin_m[i] << endl;
        binsEmitted_m[i] = true;
    }

}


size_t PartBinsCyc::getTotalNum() {
    size_t s = 0;
    size_t sd = 0;
    size_t gs = 0;

    for(int i = 0; i < getLastemittedBin(); i++) {
        s  += nBin_m[i];
        sd += nDelBin_m[i];
    }
    gs = s - sd;
    reduce(gs, gs, OpAddAssign());
    return gs;
}

size_t PartBinsCyc::getTotalNumPerBin(int b) {
    size_t s = 0;
    s  = nBin_m[b];
    reduce(s, s, OpAddAssign());
    return s;
}

void PartBinsCyc::updateStatus(const int bunchCount, const size_t partInBin) {
    // array index of binsEmitted_m[] starts from 0
    // nemittedBins_m and bins_m index starts from 1
    binsEmitted_m[bunchCount-1] = true;
    size_t NpartInBin = partInBin;
    reduce(NpartInBin, NpartInBin, OpAddAssign());
    nBin_m[bunchCount-1] = NpartInBin;
    nemittedBins_m++;
}

void PartBinsCyc::updateDeletedPartsInBin(size_t countLost[]) {
    Inform m2all("updateDeletedPartsInBin ", INFORM_ALL_NODES);

    for(int ii = 0; ii < getLastemittedBin(); ii++) {
        bool flagNeedUpdate = false;
        flagNeedUpdate = (countLost[ii] > 0);

        reduce(flagNeedUpdate, flagNeedUpdate, OpOr());

        if(flagNeedUpdate) {
            reduce(countLost[ii], countLost[ii], OpAddAssign());
            nDelBin_m[ii] = countLost[ii];
        }
        m2all << "In Bin: " << ii << ", " << nDelBin_m[ii] << " particle(s) lost" << endl;
    }
}

void PartBinsCyc::updatePartInBin(size_t countLost[]) {

    Inform msg0("updatePartInBin ");
    //  for (int ii=0; ii < bins_m; ii++){
    for(int ii = 0; ii < nemittedBins_m; ii++) {
        msg0 << "In Bin: " << ii << ", " << nBin_m[ii] << " particles " << endl;
    }
    for(int ii = 0; ii < nemittedBins_m; ii++) {
        bool flagNeedUpdate = false;
        if(countLost[ii] > 0)
            flagNeedUpdate = true;
        reduce(&flagNeedUpdate, &flagNeedUpdate + 1, &flagNeedUpdate, OpBitwiseOrAssign());
        //          reduce(flagNeedUpdate, flagNeedUpdate, OpOr());
        if(flagNeedUpdate) {
            reduce(countLost[ii], countLost[ii], OpAddAssign());
            nBin_m[ii] -= countLost[ii];
            msg0 << "In Bin: " << ii << ", " << countLost[ii] << " particle(s) lost" << endl;
        }
    }
}

void PartBinsCyc::resetPartInBin(size_t newPartNum[]) {
    for(int ii = 0; ii < nemittedBins_m; ii++) {
        reduce(newPartNum[ii], newPartNum[ii], OpAddAssign());
        nBin_m[ii] = newPartNum[ii];
        INFOMSG("After reset Bin: " << ii << ", particle(s): " << newPartNum[ii] << endl);
    }
}


void PartBinsCyc::resetPartInBin2(size_t newPartNum[], int maxbinIndex) {
    reduce(maxbinIndex, maxbinIndex, OpMaxAssign());
    // total number of bins nemittedBins_m
    nemittedBins_m =  maxbinIndex + 1;

    for(int ii = 0; ii < nemittedBins_m; ii++) {
        nBin_m[ii] = newPartNum[ii]; // only count particles on the local node
        setBinEmitted(ii);  // set true for this bin
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


bool PartBinsCyc::getPart(size_t n, int bin, vector<double> &p) {

    if(tmppart_m[n][6] == bin) {
        p = tmppart_m[n];
        return true;
    } else
        return false;
}

void PartBinsCyc::calcHBins() {

    for(int n = 0; n < tmppart_m.size(); n++)
        tmppart_m[n][6] = getBin(tmppart_m[n][2]);
    calcExtrema();
}

size_t PartBinsCyc::getSum() {
    size_t s = 0;
    for(int n = 0; n < bins_m; n++)
        s += nBin_m[n];
    return s;
}

Inform &PartBinsCyc::print(Inform &os) {

    os << "-----------------------------------------" << endl;
    os << "     CREATE BINNED GAUSS DISTRIBUTION DONE        " << endl;

    os << "Bins= " << bins_m << " hBin= " << hBin_m << " Particle vector length " << tmppart_m.size() << endl;

    for(int i = 0; i < gsl_histogram_bins(h_m); i++)
        os << "Bin # " << i << " val " << gsl_histogram_get(h_m, i) << endl;

    if(getLastemittedBin() >= 0)
        os << "Last emitted bin is " << getLastemittedBin() << endl;
    else
        os << "No bin is emitted !" << endl;
    return os;
}

int PartBinsCyc::getBin(double x) {
    /**
       returns the index of the bin to which the particle with z = 'x' belongs.
       If getBin returns b < 0 || b >= bins_m, then is x out of range!
    */
    int b = (int) floor(fabs(xmax_m - x) / hBin_m);
    nBin_m[b]++;
    return b;
}
