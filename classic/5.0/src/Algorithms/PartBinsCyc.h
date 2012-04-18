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
#else
#include "ranlib.h"
#define Inform ostream
#endif

#include "Algorithms/PartBins.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>

class PartBinsCyc: public PartBins {

private:

    double gamma_m;

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
    std::vector< std::vector<double> > tmppart_m;
    std::vector< bool > isEmitted_m;
    /** holds information whether all particles of a bin are emitted */
    //  std::vector< bool > binsEmitted_m;
    bool *binsEmitted_m;


public:


    /** constructer function for cyclotron*/
    PartBinsCyc(int bunches, int bins, size_t  partInBin[]);
    PartBinsCyc(int specifiedNumBins, int bins);
    ~PartBinsCyc();

    /** get the number of used bin */
    int getNBins() {return bins_m; }

    /** \brief How many particles are on one bin */
    inline size_t getGlobalBinCount(int bin) {
      size_t a = nBin_m[bin];
      reduce(a, a, OpAddAssign());
      return a;
    }

    bool weHaveBins() {
      return  ( nemittedBins_m > 0 );
    }

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
