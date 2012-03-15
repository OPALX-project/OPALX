#ifndef GREENWAKEFUNCTION_HH
#define GREENWAKEFUNCTION_HH

#include "Algorithms/Filters.h"
#include "Solvers/WakeFunction.hh"
#include <vector>
#include <fftw3.h>
#include "Physics/Physics.h"
#include "cassert"

using Physics::pi;

typedef map<string,int> FilterOptions;

class SavitzkyGolayFilter;

class GreenWakeFunction:public WakeFunction
{
 public:

  class S_G_FilterOptions
  {
  public:
    S_G_FilterOptions(const int &np = 33, const int &nl = 16, const int &nr = 16, const int &m = 2):
      np_m(np),
      nl_m(nl),
      nr_m(nr),
      m_m(m)
    {  
      check();
    }

    void check()
    {
      if (nl_m < 0)
        {
          cerr << "FilterOptions: the number of points to the left has to be greater than or \n"
               << "               equal to zero; resetting nl to 16" << endl;
          nl_m = 16;
        }
      if (nr_m < 0)
        {
          cerr << "FilterOptions: the number of points to the right has to be greater than or \n"
               << "               equal to zero; resetting nr to 16" << endl;
          nr_m = 16;
        }
      if ( nl_m + nr_m < m_m)
        {
          cerr << "FilterOptions: the sum of the number of points to the left and to the  right \n"
               << "has to be greater than or equal to polynomial order; resetting nr and nl to (m+1)/2" << endl;
          nl_m = nr_m = (int)((m_m+1)/2);
          np_m = nl_m + nr_m + 1;
        }
        
      if (np_m < nl_m + nr_m + 1)
        {
          cerr << "FilterOptions: the sum of the number of points to the left and to the  right \n"
               << "has to be less than or equal to the total number of points minus 1; resetting np = nl + nr" << endl;
          np_m = nl_m + nr_m + 1;
        }
    }
    int np_m; //number of points
    int nl_m; //number of points to the left
    int nr_m; //number of points to the right
    int m_m;  //polynomial order 

  };

  ~GreenWakeFunction();
  GreenWakeFunction(PartData &reference, int NBIN, double Z0, double radius, double sigma, double c, double acMode, double tau, double direction, bool constLength, const S_G_FilterOptions &options = S_G_FilterOptions(33,16,16,4));
  
  pair<int,int> distrIndices(int vectLen);
  
  void apply(PartBunch &bunch);
  void setWake(double* wakefield, double N);
  void setWake(string wFName, int NBin, double spacing);
  virtual const string getType() const;
  
 private:
  class Wake {

	public:
	
	Wake(double s, double Z0, double a, double sigma, double c, double acMode, double tau, double direction) 
	: Z0_(Z0), a_(a), sigma_(sigma), c_(c), s_(s), acMode_(acMode), tau_(tau), direction_(direction)
	{}
	
	/**
	* @brief	Used to integrate the function 
	*
	* @param[in]	k parameter
	*
	* @return	the function value at position k	
	*/
	double operator()(double k) {
		
		complex <double> i(0, 1);
		complex <double>  Z(0,0);
		double signK;
		signK = (k>0 ? 1 : -1);
		
	    if (acMode_ ==1){	// AC 
			Z = ( Z0_/(2 * pi*a_))* 1.0/(sqrt( Z0_*abs(k)/2)* sqrt(sigma_ /(1.0-i*c_*k*tau_)) *(i+signK)/k -(i*k*a_)/2.0);
	    }
	    if(acMode_ ==2){	//DC
	    	Z = ( Z0_/(2 * pi*a_))* 1.0/(sqrt(sigma_* Z0_*abs(k)/2) *(i+signK)/k -(i*k*a_)/2.0);
	    }
	    if (direction_==1){	// Longitudinal	
			return real(Z) *cos(k*s_) *2.0 *c_ / pi;  
	    }
	    else{	//Transversal
	    	return real(Z)*c_/k *cos(k*s_) *2.0 *c_ / pi; 
	    }
	}
	
	private: 
	
	/// impedanz
	double Z0_;
	/// radius
	double a_;
	/// material constant
	double sigma_;
	/// speed of light 
	double c_;
	/// distanc from the particle
	double s_;
	/// conductivity either 1="AC" or 2="DC"
	double acMode_;
	/// material constant
	double tau_;
	/// direction either 1="Longitudinal" 2= "Transversal"
	double direction_;
	
  };
  
/**
* @brief	Simpson-Integration from the function f from a to b with N steps
*
*
* @param[in] 	f the function to integrate
* @param[in] 	a integrate from a
* @param[in] 	b integrate to b
* @param[in] 	N Number of integration points
* @return 	function value of the integration
*
*/
  template<class F> double simpson(F &f, double a, double b, unsigned int N) 
	{
	assert(b>a);
	assert(N>0);

	double    result=0;
	double      h=(b-a)/N;

	// boundary values
	result += ( f(a) + 4*f(a+h/2) + f(b) ) / 2.0;

	// values between boundaries
	for ( unsigned int i = 1; i <= N-1; ++i ) {
		result += f(a+i*h) + 2*f(a+(i+0.5)*h);
	}

	result *= h/3.0;

	return result;

	}
  /// Filter to smooth the line Density of the particle bunch
  SavitzkyGolayFilter smoother_m;
  /// save the line Density of the particle bunch
  vector<double> lineDensity_m;
  /// FFT of the zerro paded wakefield
  fftw_complex  *FftWField;   
  /// divize the paricle bunch in NBin slizes
  int NBin;
  /// impedanz
  double Z0;
  /// radius
  double radius;
  /// material constant
  double sigma;
  /// speed of light 
  double c;
  /// conductivity either 1="AC" or 2="DC"
  double acMode;
  /// material constant
  double tau;
  /// direction either 1="Longitudinal" 2= "Transversal"
  double direction;
	
  /// thrue if the length of the particle bunch is considered as constant
  bool constLength;
  
  void testApply(PartBunch &bunch);
  void compEnergy(const double K, const double charge, const int LengthWake, const fftw_complex* FftWField, const double* lambda, double* OutEnergy);
  void compEnergy(const double K, const double charge, const int LengthWake, const fftw_complex* FftWField, vector<double> lambda, double* OutEnergy);
  void CalcWakeFFT(double Z0, double radius, double sigma, double c, double acMode, double tau, double direction, double spacing, int Lbunch, fftw_complex *FftWField);
};
#endif //GREENWAKEFUNCTION_HH
