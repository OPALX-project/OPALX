#ifndef GREENWAKEFUNCTION_HH
#define GREENWAKEFUNCTION_HH

#include "Filters/Filter.h"
#include "Utilities/OpalFilter.h"
#include "Solvers/WakeFunction.hh"
#include <vector>
#include <fftw3.h>
#include "Physics/Physics.h"
#include "cassert"

using Physics::pi;

enum { TRANSVERSAL, LONGITUDINAL };
typedef map<string,int> FilterOptions;

class SavitzkyGolayFilter;

class GreenWakeFunction:public WakeFunction
{
public:
    ~GreenWakeFunction();
    //IFF: changed direction to int (was double)
    //IFF: changed acMode to int (was double)
    //IFF: WHY do we pass c?
    GreenWakeFunction(PartData &reference,  vector<OpalFilter*> filters, int NBIN, double Z0, double radius, double sigma, double c, int acMode, double tau, int direction, bool constLength, string fname);

    pair<int,int> distrIndices(int vectLen);

    void apply(PartBunch &bunch);
    void setWake(double* wakefield, double N);
    void setWakeFromFile(int NBin, double spacing);
    virtual const string getType() const;

private:
    class Wake {

    public:

        Wake(double s, double Z0, double a, double sigma, double c, int acMode, double tau, int direction) 
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

            //1 == AC
            //2 == DC
            switch(acMode_){
                case 1: Z = ( Z0_/(2 * pi*a_))* 1.0/(sqrt( Z0_*abs(k)/2)* sqrt(sigma_ /(1.0-i*c_*k*tau_)) *(i+signK)/k -(i*k*a_)/2.0); break;
                case 2: Z = ( Z0_/(2 * pi*a_))* 1.0/(sqrt(sigma_* Z0_*abs(k)/2) *(i+signK)/k -(i*k*a_)/2.0); break;
            }
            switch(direction_){
                case LONGITUDINAL:  return real(Z) *cos(k*s_) *2.0 *c_ / pi; break;
                case TRANSVERSAL:   return real(Z)*c_/k *cos(k*s_) *2.0 *c_ / pi; break;
            }
        }

    private: 

        /// impedance
        double Z0_;
        /// radius
        double a_;
        /// material constant
        double sigma_;
        /// speed of light 
        double c_;
        /// distance from the particle
        double s_;
        /// conductivity either 1="AC" or 2="DC"
        int acMode_;
        /// material constant
        double tau_;
        /// direction either 1="Longitudinal" 0= "Transversal"
        int direction_;

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
    /// save the line Density of the particle bunch
    vector<double> lineDensity_m;
    /// FFT of the zero padded wakefield
    fftw_complex  *FftWField;   
    /// divides the particle bunch in NBin slices
    int NBin;
    /// impedance
    double Z0;
    /// radius
    double radius;
    /// material constant
    double sigma;
    /// speed of light 
    double c;
    /// conductivity either 1="AC" or 2="DC"
    int acMode;
    /// material constant
    double tau;
    /// direction either 1="Longitudinal" 2= "Transversal"
    int direction;
    /// true if the length of the particle bunch is considered as constant
    bool constLength;
    /// filename of the wakefield
    string filename_m;

    vector<OpalFilter*> filters_m;

    void testApply(PartBunch &bunch);
    void compEnergy(const double K, const double charge, const int LengthWake, const fftw_complex* FftWField, const double* lambda, double* OutEnergy);
    void compEnergy(const double K, const double charge, const int LengthWake, const fftw_complex* FftWField, vector<double> lambda, double* OutEnergy);
    void CalcWakeFFT(double Z0, double radius, double sigma, double c, double acMode, double tau, double direction, double spacing, int Lbunch, fftw_complex *FftWField);
};
#endif //GREENWAKEFUNCTION_HH
