/** \file
 *  \brief     Defines a structure to hold energy bins and their
 *             associated data
 *
 *
 *
 *  \author    Andreas Adelmann
 *  \date      xx. November 2007
 *
 *  \warning   None.
 *  \attention
 *  \bug       sure no bug in this code :=)
 *  \todo
 */

#ifndef OPAL_Bins_HH
#define OPAL_Bins_HH

#ifndef PartBinTest
#include "Algorithms/PBunchDefs.h"
#else
#include "ranlib.h"
#define Inform ostream
#endif

#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>


class PartBins {

    /**

    Definition of a bin:  low <= x < hi
    Low water mark is included, high water
    is excluded.

    */


public:


    PartBins(int bins, int sbins);

    virtual ~PartBins();


    /** \brief How many deleted particles are on one bin */
    inline int getDelBinCont(int bin) {
        int a = nDelBin_m[bin];
        reduce(a, a, OpAddAssign());
        return a;
    }


    /** \brief Add a particle to the temporary container */
    void fill(std::vector<double> &p) {
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

    bool isEmitted(int n, int /*bin*/) {
        return isEmitted_m[n]; //(isEmitted_m[n][0]==1) && (isEmitted_m[n][1] == bin);
    }

    void setEmitted(int n, int /*bin*/) {
        isEmitted_m[n] = true;
    }

    void updatePartPos(int n, int /*bin*/, double z) {
        tmppart_m[n][2] = z;
    }

    void updateExtramePos(int bin, double dz) {
        xbinmax_m[bin] += dz;
        xbinmin_m[bin] += dz;
    }

    /** assigns the proper position of particle n if it belongs to bin 'bin' */
    bool getPart(size_t n, int bin, std::vector<double> &p);

    /** sort the vector of particles such that positions of the particles decrease with increasing index.
        Then push the particles back by xmax_m + jifactor * bunch_length. In order that the method getBin(double x) works xmin_m has to be lowered a bit more.
    */
    void sortArray();

    /** assigns the particles to the bins */
    void calcHBins();
    size_t getSum();
    void calcGlobalExtrema();
    void calcExtrema();
    void getExtrema(double &min, double &max) {
        min = xmin_m;
        max = xmax_m;
    }

    /** update global bin parameters after inject a new bunch */
    void updateStatus(int bunchCount, size_t nPartInBin);
    /** update particles number in bin after reset Bin ID of PartBunch  */
    void resetPartInBin(size_t newPartNum[]);
    /** update local particles number in bin after reset Bin ID of PartBunch  */
    void resetPartInBin_cyc(size_t newPartNum[], int binID);
    /** update particles number in bin after particle deletion */
    void updatePartInBin(size_t countLost[]);
    /** update local particles number in bin after particle deletion */
    void updatePartInBin_cyc(size_t countLost[]);

    void updateDeletedPartsInBin(size_t countLost[]) ;

    void setGamma(double gamma) { gamma_m = gamma;}
    double getGamma() {return gamma_m;}

protected:

    double gamma_m;
    /**
       returns the index of the bin to which the particle with z = 'x' belongs.
       If getBin returns b < 0 || b >= bins_m, then x is out of range!
    */
    int getBin(double x);

    int bins_m;
    int sBins_m;

    /** extremal particle positions */
    double xmin_m;
    double xmax_m;

    /** extremal particle position within the bins */
    std::unique_ptr<double[]> xbinmin_m;
    std::unique_ptr<double[]> xbinmax_m;

    /** bin size */
    double hBin_m;

    /** holds the particles not yet in the bunch */
    std::vector< std::vector<double> > tmppart_m;
    std::vector< bool > isEmitted_m;
    /** holds information whether all particles of a bin are emitted */
    //  std::vector< bool > binsEmitted_m;
    std::unique_ptr<bool[]> binsEmitted_m;

    /**
        Here comes the new stuff, t-binning
    */

public:

    Inform &print(Inform &os);

    int getSBins() { return sBins_m; };

    /** get the number of used bin */
    virtual int getNBins() {return gsl_histogram_bins(h_m.get()) / sBins_m; }

    /** Get the total number of sampled bins */
    virtual int getNSBins() {return gsl_histogram_bins(h_m.get()); }

    int getBinToEmit() {
        int save;
        if(nemittedBins_m < getNBins()) {
            save =  nemittedBins_m;
            nemittedBins_m++;
            return save;
        } else
            return -1;
    }

    int getSBinToEmit() {
        int save;
        if(nemittedSBins_m < getNSBins()) {
            save = nemittedSBins_m;
            nemittedSBins_m++;
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
        h_m.reset(nullptr);
    }

    virtual bool weHaveBins() {
        return h_m != nullptr;
    }

    /** sort the vector of particles according to the bin number */
    void sortArrayT();

    inline void setHistogram(gsl_histogram *h) { h_m.reset(h);}

    /** \brief How many particles are on one bin */
    virtual inline size_t getGlobalBinCount(int bin) {
        size_t a = gsl_histogram_get(h_m.get(), bin);
        reduce(a, a, OpAddAssign());
        return a;
    }

    /** \brief How many particles are on one energy bin */
    inline size_t getLocalBinCount(int bin) {
        size_t ret = 0;
        for(int i = sBins_m * bin; i < sBins_m * (bin + 1); i++) {
            ret += gsl_histogram_get(h_m.get(), i);
        }
        return ret;
    }

    /** \brief How many particles are in one sampling bin */
    inline size_t getLocalSBinCount(int bin) { return gsl_histogram_get(h_m.get(), bin);}


    /** \brief How many particles are in all the bins */
    size_t getTotalNum();

    /** \brief How many particles are in the bin b */
    size_t getTotalNumPerBin(int b);


protected:

    /** number of emitted bins */
    int nemittedBins_m;

    /** Number of total sampled bins emitted */
    int nemittedSBins_m;

    /** number of particles in the bins, the sum of all the nodes */
    std::unique_ptr<size_t[]> nBin_m;

    /** number of deleted particles in the bins */
    std::unique_ptr<size_t[]> nDelBin_m;

    std::unique_ptr<gsl_histogram> h_m;

};

inline Inform &operator<<(Inform &os, PartBins &p) {
    return p.print(os);
}


class AscendingLocationSort: public std::binary_function< std::vector<double>, std::vector<double>, bool> {
public:
    AscendingLocationSort(int direction = 0): direction_m(direction)
    {;}

    bool operator()(const std::vector<double>& first_part, const std::vector<double>& second_part) {
        return first_part[direction_m] < second_part[direction_m];
    }
private:
    int direction_m;
};

class DescendingLocationSort: public std::binary_function< std::vector<double>, std::vector<double>, bool> {
public:
    DescendingLocationSort(int direction = 0): direction_m(direction)
    {;}

    bool operator()(const std::vector<double>& first_part, const std::vector<double>& second_part) {
        return first_part[direction_m] > second_part[direction_m];
    }
private:
    int direction_m;
};

#endif // OPAL_Bins_HH
