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
 #include "Ippl.h"
#else
 #include "ranlib.h"
 #define Inform ostream
#endif


using namespace std;


class PartBins {

  /**

    Definition of a bin:  low <= x < hi
    Low water mark is included, high water
    is excluded.

  */

public:


  PartBins(int bins);

  /** constructer function for cyclotron*/
  PartBins(int bunches, int bins, size_t  partInBin[]);
  ~PartBins();

  /** \brief How many particles are on one bin */
  inline size_t getBinCont(int bin) { 
    size_t a = nBin_m[bin];
    reduce(a, a, OpAddAssign());
    return a;
  }

  /** \brief How many deleted particles are on one bin */
  inline int getDelBinCont(int bin) { 
    int a = nDelBin_m[bin]; 
    reduce(a, a, OpAddAssign());
    return a;
  }


  /** \brief How many particles are in all the bins */
  size_t getTotalNum();

  /** \brief How many particles are in the bin b */
  size_t getTotalNumPerBin(int b);

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

  /** get the number of used bin */
  int getNBins() {return bins_m; }

  /** set particles number in given bin */
  void setPartNum(int bin, long long num) {nBin_m[bin] = num;}

  /** assume we emmit in monotinic increasing order */
  void setBinEmitted(int bin) {binsEmitted_m[bin] = true;}

  bool getBinEmitted(int bin) {return binsEmitted_m[bin];}

  /** the last emitted bin is always smaller or equal getNbins */
  int getLastemittedBin() {return nemittedBins_m; }

  /** set the actual emitted bib */
  void setActualemittedBin(int bin) {nemittedBins_m=bin; }

  /** \brief Is true if we still have particles to emit */
  bool doEmission() {return getNp() != 0;}

  /** \brief If the bunch object rebins we need to call resetBins() */
  void resetBins() {bins_m = 0;}

  bool weHaveBins() { return bins_m != 0;}

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

  /** sort the vector of particles such that positions of the particles decrease with increasing index.
      Then push the particles back by xmax_m + jifactor * bunch_length. In order that the method getBin(double x) works xmin_m has to be lowered a bit more.
    */
  void sortArray();

  /** assigns the particles to the bins */
  void calcHBins();
  size_t getSum();
  void calcGlobalExtrema();
  void calcExtrema();
  void getExtrema(double &min, double &max)
    {
      min = xmin_m;
      max = xmax_m;
    }

  Inform &print(Inform &os);

  /** Set energy [keV] to define a rebin condition */

  void setRebinEnergy( double e ) { dERebin_m = e; }
  double getRebinEnergy() { return dERebin_m; }

  /** update global bin parameters after inject a new bunch */
  void updateStatus(int bunchCount, size_t nPartInBin);
  /** update particles number in bin after reset Bin ID of PartBunch  */
  void resetPartInBin( size_t newPartNum[]);
  /** update particles number in bin after particle deletion */
  void updatePartInBin(size_t countLost[]);

  void updateDeletedPartsInBin(size_t countLost[]) ;

  void setGamma(double gamma) { gamma_m = gamma;}
  double getGamma() {return gamma_m;}

private:

  /** Defines energy threshold for rebining */
  double dERebin_m;

  double gamma_m;
  /**
     returns the index of the bin to which the particle with z = 'x' belongs.
     If getBin returns b < 0 || b >= bins_m, then x is out of range!
  */
  int getBin(double x);

  int bins_m;
  int nemittedBins_m;

  /** extremal particle positions */
  double xmin_m;
  double xmax_m;

  /** extremal particle position within the bins */
  double *xbinmin_m;
  double *xbinmax_m;

  /** bin size */
  double hBin_m;

  /** number of particles in the bins, the sum of all the nodes */
  size_t *nBin_m;

  /** number of deleted particles in the bins */
  size_t *nDelBin_m;

  /** holds the particles not yet in the bunch */
  vector< vector<double> > tmppart_m;
  vector< bool > isEmitted_m;
  /** holds information whether all particles of a bin are emitted */
  //  vector< bool > binsEmitted_m;
  bool *binsEmitted_m;
};

inline Inform &operator<<(Inform &os, PartBins &p)
{
  return p.print(os);
}


class AscendingLocationSort: public binary_function< vector<double>, vector<double>, bool>
{
public:
  AscendingLocationSort(int direction = 0):direction_m(direction)
    {;}

  bool operator()(const vector<double>& first_part, const vector<double>& second_part)
    {
      return first_part[direction_m] < second_part[direction_m];
    }
 private:
  int direction_m;
};

class DescendingLocationSort: public binary_function< vector<double>, vector<double>, bool>
{
public:
  DescendingLocationSort(int direction = 0):direction_m(direction)
    {;}

  bool operator()(const vector<double>& first_part, const vector<double>& second_part)
    {
      return first_part[direction_m] > second_part[direction_m];
    }
 private:
  int direction_m;
};

#endif // OPAL_Bins_HH
