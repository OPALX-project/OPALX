/** \file
 *  \brief     Defines a structure to hold energy bins and their
 *             associated data for multi-bunch tracking in cyclotrons
 *
 *
 *
 *  \author    Jianjun Yang
 *  \date      01. June 2010
 *
 *  \warning   None.
 *  \attention
 *  \bug
 *  \todo
 */

#ifndef OPAL_BinsCyc_HH
#define OPAL_BinsCyc_HH

#ifndef PartBinTest
#include "Ippl.h"
#else
#include "ranlib.h"
#define Inform ostream
#endif

#include "Algorithms/PartBins.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

using namespace std;

class PartBinsCyc: PartBins {

public:

    /** constructer function for cyclotron*/
    PartBinsCyc(int bunches, int bins, size_t  partInBin[]);
    ~PartBinsCyc();


    /** \brief How many deleted particles are on one bin */
    inline int getDelBinCont(int bin) {
        int a = nDelBin_m[bin];
        reduce(a, a, OpAddAssign());
        return a;
    }


    /** \brief Add a particle to the temporary container */
    void fill(vector<double> &p) {
        tmppart_m.push_back(p);
        isEmitted_m.push_back(false);
    }

    /** \brief  get the number of particles in the temporary particle structure used for binning */
    size_t getNp() {
        //    size_t a = tmppart_m.size();
        // reduce(a, a, OpAddAssign());
        // return a;
        return tmppart_m.size();
    }

    /** set particles number in given bin */
    void setPartNum(int bin, long long num) {nBin_m[bin] = num;}

    /** assume we emmit in monotinic increasing order */
    void setBinEmitted(int bin) {binsEmitted_m[bin] = true;}

    bool getBinEmitted(int bin) {return binsEmitted_m[bin];}

    /** \brief Is true if we still have particles to emit */
    bool doEmission() {return getNp() > 0;}

    bool isEmitted(int n, int bin) {
        return isEmitted_m[n]; //(isEmitted_m[n][0]==1) && (isEmitted_m[n][1] == bin);
    }

    void setEmitted(int n, int bin) {
        isEmitted_m[n] = true;
    }

    void updatePartPos(int n, int bin, double z) {
        tmppart_m[n][2] = z;
    }

    void updateExtramePos(int bin, double dz) {
        xbinmax_m[bin] += dz;
        xbinmin_m[bin] += dz;
    }

    /** assigns the proper position of particle n if it belongs to bin 'bin' */
    bool getPart(size_t n, int bin, vector<double> &p);


    /** assigns the particles to the bins */
    void calcHBins();
    size_t getSum();

    /** update global bin parameters after inject a new bunch */
    void updateStatus(int bunchCount, size_t nPartInBin);
    /** update particles number in bin after reset Bin ID of PartBunch  */
    void resetPartInBin(size_t newPartNum[]);
    /** update particles number in bin after reset Bin ID of PartBunch  */
    void resetPartInBin2(size_t newPartNum[], int binID);
    /** update particles number in bin after particle deletion */
    void updatePartInBin(size_t countLost[]);

    void updateDeletedPartsInBin(size_t countLost[]) ;

    void setGamma(double gamma) { gamma_m = gamma;}
    double getGamma() {return gamma_m;}

private:

    double gamma_m;
    /**
       returns the index of the bin to which the particle with z = 'x' belongs.
       If getBin returns b < 0 || b >= bins_m, then x is out of range!
    */
    int getBin(double x);

    int bins_m;


    /** extremal particle positions */
    double xmin_m;
    double xmax_m;

    /** extremal particle position within the bins */
    double *xbinmin_m;
    double *xbinmax_m;

    /** bin size */
    double hBin_m;

    /** holds the particles not yet in the bunch */
    vector< vector<double> > tmppart_m;
    vector< bool > isEmitted_m;
    /** holds information whether all particles of a bin are emitted */
    //  vector< bool > binsEmitted_m;
    bool *binsEmitted_m;

    /**
        Here comes the new stuff, t-binning
    */

public:

    Inform &print(Inform &os);

    /** Set energy [keV] to define a rebin condition */

    void setRebinEnergy(double e) { dERebin_m = e; }

    double getRebinEnergy() { return dERebin_m; }

    /** get the number of used bin */
    int getNBins() {return gsl_histogram_bins(h_m); }

    int getBinToEmit() {
        int save;
        if(nemittedBins_m < getNBins()) {
            save =  nemittedBins_m;
            nemittedBins_m++;
            return save;
        } else
            return -1;
    }

    /** the last emitted bin is always smaller or equal getNbins */
    int getLastemittedBin() {return nemittedBins_m; }

    /** set the actual emitted bib */
    void setActualemittedBin(int bin) {nemittedBins_m = bin; }

    /** \brief If the bunch object rebins we need to call resetBins() */
    void resetBins() {
        if(h_m)
            delete h_m;
        h_m = NULL;
    }

    bool weHaveBins() {
        return h_m != NULL;
    }

    inline void setHistogram(gsl_histogram *h) { h_m = h;}

    /** \brief How many particles are on one bin */
    inline size_t getGlobalBinCount(int bin) {
        size_t a = gsl_histogram_get(h_m, bin);
        reduce(a, a, OpAddAssign());
        return a;
    }

    /** \brief How many particles are on one bin */
    inline size_t getLocalBinCount(int bin) {
        return gsl_histogram_get(h_m, bin);
    }


    /** \brief How many particles are in all the bins */
    size_t getTotalNum();

    /** \brief How many particles are in the bin b */
    size_t getTotalNumPerBin(int b);






private:

    /** Defines energy threshold for rebining */
    double dERebin_m;

    /** number of emitted bins */
    int nemittedBins_m;

    /** number of particles in the bins, the sum of all the nodes */
    size_t *nBin_m;

    /** number of deleted particles in the bins */
    size_t *nDelBin_m;

    gsl_histogram *h_m;

};

#endif // OPAL_BinsCyc_HH
