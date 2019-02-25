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

#include <cassert>

class PartBinsCyc: public PartBins {

public:


    /** constructer function for cyclotron*/
    PartBinsCyc(int bunches, int bins, size_t  partInBin[]);
    PartBinsCyc(int specifiedNumBins, int bins);

    /** get the number of used bin */
    int getNBins() {return bins_m; }

    /** \brief How many particles are on one bin */
    inline size_t getGlobalBinCount(int bin) {
      size_t a = nBin_m[bin];
      reduce(a, a, OpAddAssign());
      return a;
    }

    /** \brief How many particles are on one energy bin */
    inline size_t getLocalBinCount(int bin) {
        assert(bin < bins_m);
        return nBin_m[bin];
    }

    bool weHaveBins() {
      return  ( nemittedBins_m > 0 );
    }
};

#endif // OPAL_BinsCyc_HH
