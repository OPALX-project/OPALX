#ifndef __TUNE__
#define __TUNE__
/*****************************************************************************/
/*                                                                           */
/* TUNE's Class Header                                                       */
/* =====================                                                     */
/*                                                                           */
/*****************************************************************************/

using namespace std;

class TUNE_class
/*---------------------------------------------------------------------------*/
{
private:  

  // lomb stuff:
  double ofac, hifac;
  double Qmin,Qmax;

public:

  TUNE_class(); //constructor
  virtual ~TUNE_class(void);       //destructor

  int TUNE_class::LombAnalysis(double *x, double *y, int Ndat, int nhis);
  int TUNE_class::LombAnalysis(vector<double> &x, vector<double> &y, int nhis, double Norm);
   
  
};
/*---------------------------------------------------------------------------*/

#endif
