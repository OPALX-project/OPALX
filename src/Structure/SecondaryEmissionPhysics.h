#ifndef OPAL_SECONDARY_EMISSION_PHYSICS_HH
#define OPAL_SECONDARY_EMISSION_PHYSICS_HH

#include <hdf5.h>
#include "Ippl.h"
#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <set>
#include <cassert>
#include <cmath>
#include <limits>
#include <sys/stat.h>
#include "Physics/Physics.h"
//#define assert2(_Expression, _Msg) (void)( (!!(_Expression)) || (_wassert(_Msg, _CRT_WIDE(__FILE__), __LINE__), 0) )
extern Inform *gmsg;
namespace myeps {
    const double EPS = numeric_limits<double>::epsilon();
    const double FPMIN = numeric_limits<double>::min() / EPS;
}
using namespace myeps;
using namespace Physics;
class SecondaryEmissionPhysics

{
public:

    /**
     * Exemplar Constructor
     */

    SecondaryEmissionPhysics();

    /**
     * Destructor.
     *
     * Delete the previous defined member arrays
     */
    ~SecondaryEmissionPhysics();

    void nSec(const double &incEnergy,  const double &cosTheta, const int &matNumber, int &seNum, int &seType, const double &incQ, const Vector_t &TriNorm, const Vector_t &inteCoords, const Vector_t &localX, PartBunch *itsBunch, double &seyNum, const double &ppVw, const double &vVThermal, const bool nEmissionMode);
    void nSec(const double &incEnergy,  const double &cosTheta, int &seNum, int &seType, const double &incQ, const Vector_t &TriNorm, const Vector_t &inteCoords, const Vector_t &localX, PartBunch *itsBunch, double &seyNum, const double &ppVw,const double &vSeyZero, const double &vEzero, const double &vSeyMax, const double &vEmax, const double &vKenergy, const double &vKtheta, const double &vVThermal, const bool nEmissionMode);
    double deltae_m;
    double deltar_m;
    double deltats_m;
    private:

    /**
     * @param TPnSec_m is the timmer of secondary emission module.
     */

    IpplTimings::TimerRef TPnSec_m;

    /**
     * @param seAlpha_m is the emitted angular spectrum parameter, i.e., the 1st parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seAlpha_m;


    /*==============================================================================================
      parameters for backscattered electrons
      ==============================================================================================*/



    /**
     * @param sePScat_m is the 2nd parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double sePScat_m;
    /**
     * @param sePScatPeak_m is the 3rd parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double sePScatPeak_m;
    /**
     * @param seEScatPeak_m is the 4th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seEScatPeak_m;
    /**
     * @param seW_m is the 5th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seW_m;
    /**
     * @param seP_m is the 6th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seP_m;
    /**
     * @param seDelta_m is the 7th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seDelta_m;
    /**
     * @param seEOne_m is the 8th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seEOne_m;
    /**
     * @param seETwo_m is the 9th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seETwo_m;



    /*==============================================================================================
      parameters for rediffused electrons
      ==============================================================================================*/



    /**
     * @param sePRed_m is the 10th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double sePRed_m;
    /**
     * @param seERed_m is the 11th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seERed_m;
    /**
     * @param seR_m is the 12th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seR_m;
    /**
     * @param seQ_m is the 13th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seQ_m;
    /**
     * @param seROne_m is the 14th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seROne_m;
    /**
     * @param seRTwo_m is the 15th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seRTwo_m;


    /*==============================================================================================
      parameters for true secondary electrons
      ==============================================================================================*/

    /**
     * @param seYPeakTS_m is the 16th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seYPeakTS_m;
    /**
     * @param seEPeakTS_m is the 17th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seEPeakTS_m;
    /**
     * @param seSTS_m is the 18th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seSTS_m;
    /**
     * @param seTOneTS_m is the 19th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seTOneTS_m;
    /**
     * @param seTTwoTS_m is the 20th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seTTwoTS_m;
    /**
     * @param seTThreeTS_m is the 21st parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seTThreeTS_m;
    /**
     * @param seTFourTS_m is the 22nd parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seTFourTS_m;
    /**
     * @param seEPeakTot_m is the 23rd parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seEPeakTot_m;
    /**
     * @param seYPeakTot_m is the 24th parameter in the TABLE I of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
     */
    double seYPeakTot_m;

    /*==============================================================================================
      parameters in he TABLE II of Phys. Rev. ST Accel. Beams 5, 124404 (2002)
      ==============================================================================================*/
    double sePn_m[10];
    double seEpsn_m[10];


    void calcEmiNum(const double incEnergy, const double cosTheta, int &seNum, const double &vSeyZero, const double &vEzero, const double &vSeyMax, const double &vEmax, const double &vKenergy, const double &vKtheta, double &seyNum) {// For Vaughan's model.
       
        double vSEY = 0;

        if (incEnergy<vEzero) {
            vSEY = vSeyZero;
        } else {

            double theta = acos(cosTheta);
            double delta_max = vSeyMax*(1.0+vKtheta*theta*theta/Physics::two_pi);//here the symbols k are different with reference.
            double E_max = vEmax*(1.0+vKenergy*theta*theta/Physics::two_pi);
            assert(E_max-vEzero>0);
            double v = (incEnergy-vEzero)/(E_max-vEzero);

            if (v<=3.6) {

                if (v<1.0) {
                    vSEY = delta_max*pow(v*exp(1.0-v),0.56);
                }else {
                    vSEY = delta_max*pow(v*exp(1.0-v),0.25);
                }

            }else {
                vSEY = delta_max*1.125/pow(v,0.35);
            }

        }
        double L = exp(-vSEY);// poisson distribution: Knuth's algorithm.
        int k = 0;
        double p = 1.0;
        do {
            k++;
            double u = IpplRandom();
            p*=u;
        }
        while (p>L);
        seNum = k-1;
        seyNum = vSEY;
     
    }
    
    void calcEmiNum(const double incEnergy, const double cosTheta, const double *prob, int &seNum) {// For Furman-Pivi's model

        double prob_max = 0.0;
        // Acceptance-rejection methods to generate random number with specified distribution.

        for(int i = 0; i < 11; i++) {

            if(prob[i] > prob_max) {
                prob_max = prob[i];
            }

        }
        double pY  = 1.0;
        double pX = 0.0;
        while(pY > pX) {

            double rand1 = IpplRandom();
            //double rand1 = (*rand)(rng);
            seNum = (int)(rand1 * 11.0);//fix me
            pX = prob[seNum];
            double rand2 = IpplRandom();
            //double rand2 = (*rand)(rng);
            pY = prob_max * rand2;
        }

    }
  

    double calcDeltats(const double incEnergy, const double cosTheta) {

        double seypeak = seYPeakTS_m * (1 + seTOneTS_m * (1.0 - pow(cosTheta, seTTwoTS_m))); //formula III.E (48a)
        double seepeak = seEPeakTS_m * (1 + seTThreeTS_m * (1.0 - pow(cosTheta, seTFourTS_m))); //formula III.E (48b)
        double tmpx = incEnergy / seepeak;
        double tmpD = seSTS_m * tmpx / (seSTS_m - 1 + pow(tmpx, seSTS_m)); //formula III.D (32)
        double ret = seypeak * tmpD; //formula III.D (31)
        return ret;

    }

    double calcDeltar(const double incEnergy, const double cosTheta) {

        double tmp = pow(incEnergy / seERed_m, seR_m);
        double ret = sePRed_m * (1.0 - exp(-1 * tmp)); //formula III.D (28)
        ret = ret * (1.0 + seROne_m * (1.0 - pow(cosTheta, seRTwo_m))); //formula III.E (47b)
        return ret;

    }


    double calcDeltae(const double incEnergy, const double cosTheta) {

        double tmp = pow(fabs(incEnergy - seEScatPeak_m) / seW_m, seP_m) / seP_m;
        double ret = sePScat_m + (sePScatPeak_m - sePScat_m) * exp(-1 * tmp); //formula III.D (25)
        ret = ret * (1.0 + seEOne_m * (1.0 - pow(cosTheta, seETwo_m))); //formula III.E (47a)
        return ret;

    }


    double calcProb(const double incEnergy, const double cosTheta, double *prob) {

        deltae_m = calcDeltae(incEnergy, cosTheta);
        deltar_m = calcDeltar(incEnergy, cosTheta);
        deltats_m = calcDeltats(incEnergy, cosTheta);
	
        double tmp = 1.0 - deltae_m - deltar_m;
        double p = deltats_m / tmp / 10.0;
        double q = 1.0 - p;
        double b[11];
        b[0]  = 1.0;
        b[1]  = 10.0;
        b[2]  = 45.0;
        b[3]  = 120.0;
        b[4]  = 210.0;
        b[5]  = 252.0;
        b[6]  = 210.0;
        b[7]  = 120.0;
        b[8]  = 45.0;
        b[9]  = 10.0;
        b[10] = 1.0;
        for(int i = 0; i < 11; i++) {
            prob[i] = tmp * b[i] * pow(p, i) * pow(q, (10 - i));
        }
        prob[1] = prob[1] + deltae_m + deltar_m;

        /*==============================================*/
    
	return (deltae_m+deltar_m+deltats_m);
        //cout << "sum prob: " << sum << endl;
        /*==============================================*/
    }

    void setSeMaterial(int material_num) {

        if(material_num == 0) {
            seAlpha_m = 1.0;
            sePScat_m = 0.02;
            sePScatPeak_m = 0.496;
            seEScatPeak_m = 0;
            seW_m = 60.86;
            seP_m = 1.0;
            seDelta_m = 2.0;
            seEOne_m = 0.26;
            seETwo_m = 2;

            sePRed_m = 0.2;
            seERed_m = 0.041;
            seR_m = 0.104;
            seQ_m = 0.5;
            seROne_m = 0.26;
            seRTwo_m = 2;

            seYPeakTS_m = 1.8848;
            seEPeakTS_m = 276.8;
            seSTS_m = 1.54;
            seTOneTS_m = 0.66;
            seTTwoTS_m = 0.8;
            seTThreeTS_m = 0.7;
            seTFourTS_m = 1.0;
            seEPeakTot_m = 271;
            seYPeakTot_m = 2.1;
            sePn_m[0] = 2.5;
            sePn_m[1] = 3.3;
            sePn_m[2] = 2.5;
            sePn_m[3] = 2.5;
            sePn_m[4] = 2.8;
            sePn_m[5] = 1.3;
            sePn_m[6] = 1.5;
            sePn_m[7] = 1.5;
            sePn_m[8] = 1.5;
            sePn_m[9] = 1.5;
            seEpsn_m[0] = 1.5;
            seEpsn_m[1] = 1.75;
            seEpsn_m[2] = 1.0;
            seEpsn_m[3] = 3.75;
            seEpsn_m[4] = 8.5;
            seEpsn_m[5] = 11.5;
            seEpsn_m[6] = 2.5;
            seEpsn_m[7] = 3.0;
            seEpsn_m[8] = 2.5;
            seEpsn_m[9] = 3.0;
        }
        if(material_num == 1) {
            seAlpha_m = 1.0;
            sePScat_m = 0.07;
            sePScatPeak_m = 0.5;
            seEScatPeak_m = 0;
            seW_m = 100;
            seP_m = 0.9;
            seDelta_m = 1.9;
            seEOne_m = 0.26;
            seETwo_m = 2;

            sePRed_m = 0.74;
            seERed_m = 40;
            seR_m = 1;
            seQ_m = 0.4;
            seROne_m = 0.26;
            seRTwo_m = 2;

            seYPeakTS_m = 1.22;
            seEPeakTS_m = 310;
            seSTS_m = 1.813;
            seTOneTS_m = 0.66;
            seTTwoTS_m = 0.8;
            seTThreeTS_m = 0.7;
            seTFourTS_m = 1.0;
            seEPeakTot_m = 292;
            seYPeakTot_m = 2.05;

            sePn_m[0] = 1.6;
            sePn_m[1] = 2.0;
            sePn_m[2] = 1.8;
            sePn_m[3] = 4.7;
            sePn_m[4] = 1.8;
            sePn_m[5] = 2.4;
            sePn_m[6] = 1.8;
            sePn_m[7] = 1.8;
            sePn_m[8] = 2.3;
            sePn_m[9] = 1.8;
            seEpsn_m[0] = 3.9;
            seEpsn_m[1] = 6.2;
            seEpsn_m[2] = 13.0;
            seEpsn_m[3] = 8.8;
            seEpsn_m[4] = 6.25;
            seEpsn_m[5] = 2.25;
            seEpsn_m[6] = 9.2;
            seEpsn_m[7] = 5.3;
            seEpsn_m[8] = 17.8;
            seEpsn_m[9] = 10;
        }
    }



    /*==========================================
     *http://www.taygeta.com/random/gaussian.html
     *return a gaussian distributed random number
     *==========================================*/

    double gaussRand() {

        double x1;
        double x2;
        double w;
        do {
            x1 = 2.0 * IpplRandom() - 1.0;
            //x1 = 2.0 * (*rand)(rng) - 1.0;
            x2 = 2.0 * IpplRandom() - 1.0;
            //x2 = 2.0 * (*rand)(rng) - 1.0;
            w = x1 * x1 + x2 * x2;
        } while(w >= 1.0);

        w = sqrt((-2.0 * log(w)) / w);
        double ret = x1 * w;
        return ret;

    }

    double gammp(const double a, const double x) {
        // Returns the incomplete gamma function P .a; x/.
        static const int ASWITCH = 100; //When to switch to quadrature method.
        if(x < 0.0 || a <= 0.0)
            throw("bad args in gammp");
        if(x == 0.0)
            return 0.0;
        else if((int)a >= ASWITCH)
            return gammpapprox(a, x, 1); //Quadrature.
        else if(x < a + 1.0)
            return gser(a, x); //Use the series representation.
        else
            return 1.0 - gcf(a, x); //Use the continued fraction representation.
    }

    double gser(const double a, const double x) {
        // Returns the incomplete gamma function P .a; x/ evaluated by its series representation.
        // Also sets ln .a/ as gln. User should not call directly.
        double sum, del, ap;
        double gln = gammln(a);
        ap = a;
        del = sum = 1.0 / a;
        for(;;) {
            ++ap;
            del *= x / ap;
            sum += del;
            if(fabs(del) < fabs(sum)*EPS) {
                return sum * exp(-x + a * log(x) - gln);
            }
        }
    }
    double gcf(const double a, const double x) {
        // Returns the incomplete gamma function Q.a; x/ evaluated by its continued fraction representation. Also sets ln .a/ as gln. User should not call directly.
        int i;
        double an, b, c, d, del, h;
        double gln = gammln(a);
        b = x + 1.0 - a; // Set up for evaluating continued fraction
        c = 1.0 / FPMIN; // by modified Lentz¡¯s method (5.2)
        d = 1.0 / b; // with b0 D 0.
        h = d;
        for(i = 1;; i++) {
            //Iterate to convergence.
            an = -i * (i - a);
            b += 2.0;
            d = an * d + b;
            if(fabs(d) < FPMIN) d = FPMIN;
            c = b + an / c;
            if(fabs(c) < FPMIN) c = FPMIN;
            d = 1.0 / d;
            del = d * c;
            h *= del;
            if(fabs(del - 1.0) <= EPS) break;
        }
        return exp(-x + a * log(x) - gln) * h; //Put factors in front.
    }
    double gammpapprox(double a, double x, int psig) {
        // Incomplete gamma by quadrature. Returns P .a; x/ or Q.a; x/, when psig is 1 or 0,respectively. User should not call directly.

        const double y[18] = {0.0021695375159141994,
                              0.011413521097787704, 0.027972308950302116, 0.051727015600492421,
                              0.082502225484340941, 0.12007019910960293, 0.16415283300752470,
                              0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
                              0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
                              0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
                              0.87126389619061517, 0.95698180152629142
                             };
        const double w[18] = {0.0055657196642445571,
                              0.012915947284065419, 0.020181515297735382, 0.027298621498568734,
                              0.034213810770299537, 0.040875750923643261, 0.047235083490265582,
                              0.053244713977759692, 0.058860144245324798, 0.064039797355015485,
                              0.068745323835736408, 0.072941885005653087, 0.076598410645870640,
                              0.079687828912071670, 0.082187266704339706, 0.084078218979661945,
                              0.085346685739338721, 0.085983275670394821
                             };
        int j;
        double xu, t, sum, ans;
        double a1 = a - 1.0, lna1 = log(a1), sqrta1 = sqrt(a1);
        int ngau = 18;
        double gln = gammln(a);// Set how far to integrate into the tail:
        if(x > a1) xu = max_Gamma(a1 + 11.5 * sqrta1, x + 6.0 * sqrta1);
        else xu = max_Gamma(0., min_Gamma(a1 - 7.5 * sqrta1, x - 5.0 * sqrta1));
        sum = 0;
        for(j = 0; j < ngau; j++) {
            //Gauss-Legendre.
            t = x + (xu - x) * y[j];
            sum += w[j] * exp(-(t - a1) + a1 * (log(t) - lna1));
        }
        ans = sum * (xu - x) * exp(a1 * (lna1 - 1.) - gln);
        return (psig ? (ans > 0.0 ? 1.0 - ans : -ans) : (ans >= 0.0 ? ans : 1.0 + ans));
    }
    double gammln(const double xx) {
        // Returns the value ln%GÅ%@.xx/ for xx > 0.
        int j;
        double x, tmp, y, ser;
        static const double cof[14] = {57.1562356658629235, -59.5979603554754912,
                                       14.1360979747417471, -0.491913816097620199, .339946499848118887e-4,
                                       .465236289270485756e-4, -.983744753048795646e-4, .158088703224912494e-3,
                                       -.210264441724104883e-3, .217439618115212643e-3, -.164318106536763890e-3,
                                       .844182239838527433e-4, -.261908384015814087e-4, .368991826595316234e-5
                                      };
        if(xx <= 0) throw("bad arg in gammln");
        y = x = xx;
        tmp = x + 5.24218750000000000; // Rational 671/128.
        tmp = (x + 0.5) * log(tmp) - tmp;
        ser = 0.999999999999997092;
        for(j = 0; j < 14; j++) ser += cof[j] / ++y;
        return tmp + log(2.5066282746310005 * ser / x);
    }

    double invgammp(double p, double a) {
        //Returns x such that P .a; x/ D p for an argument p between 0 and 1.
        int j;
        double x, err, t, u, pp, lna1, afac, a1 = a - 1;
        const double EPS = 1.e-8; //Accuracy is the square of EPS.
        double gln = gammln(a);
        if(a <= 0.) throw("a must be pos in invgammap");
        if(p >= 1.) return max_Gamma(100., a + 100.*sqrt(a));
        if(p <= 0.) return 0.0;
        if(a > 1.) {
            //Initial guess based on reference [1].
            lna1 = log(a1);
            afac = exp(a1 * (lna1 - 1.) - gln);
            pp = (p < 0.5) ? p : 1. - p;
            t = sqrt(-2.*log(pp));
            x = (2.30753 + t * 0.27061) / (1. + t * (0.99229 + t * 0.04481)) - t;
            if(p < 0.5) x = -x;
            x = max_Gamma(1.e-3, a * pow(1. - 1. / (9.*a) - x / (3.*sqrt(a)), 3));
        } else {
            //Initial guess based on equations (6.2.8) and (6.2.9).
            t = 1.0 - a * (0.253 + a * 0.12);

            if(p < t) x = pow(p / t, 1. / a);
            else x = 1. - log(1. - (p - t) / (1. - t));
        }
        for(j = 0; j < 12; j++) {
            if(x <= 0.0) return 0.0;
            //x too small to compute accurately.
            err = gammp(a, x) - p;
            if(a > 1.) t = afac * exp(-(x - a1) + a1 * (log(x) - lna1));
            else t = exp(-x + a1 * log(x) - gln);
            u = err / t;
            x -= (t = u / (1. - 0.5 * min_Gamma(1., u * ((a - 1.) / x - 1))));
            //Halley¡¯s method.
            if(x <= 0.) x = 0.5 * (x + t);
            //Halve old value if x tries to go negative.
            if(fabs(t) < EPS * x) break;
        }
        return x;
    }
    double min_Gamma(double x, double y) {
        if(x > y)
            return y;
        else
            return x;
    }
    double max_Gamma(double x, double y) {
        if(x < y)
            return y;
        else
            return x;
    }


    double betai(const double x, const double a, const double b) {
        //cout<<"betai called"<<endl;
        //Returns incomplete beta function Ix.a; b/ for positive a and b, and x between 0 and 1.
        static const int SWITCH = 3000;
        double bt;
        /*=========================debug code=====================
          if (a==b) {
          cout<< " x in betai "<<x<<endl;
          }
          =========================debug code=====================*/
        if(a <= 0.0 || b <= 0.0) {
            //cout<<"betai 1"<<endl;
            throw("Bad a or b in routine betai");
        }

        if(x < 0.0 || x > 1.0) {
            //cout<<"betai 2"<<endl;
            throw("Bad x in routine betai");

        }
        if(x == 0.0 || x == 1.0) {
            //cout<<"betai 3"<<endl;
            return x;
        }
        if(a > SWITCH && b > SWITCH) {
            //cout<<"betai 4"<<endl;
            return betaiapprox(a, b, x);
        }
        bt = exp(gammln(a + b) - gammln(a) - gammln(b) + a * log(x) + b * log(1.0 - x));
        if(x < (a + 1.0) / (a + b + 2.0)) {
            //cout<<"betai 5"<<endl;
            return bt * betacf(a, b, x) / a;
        } else {
            //cout<<"betai 6"<<endl;
            /*=========================debug code=====================
              if (a==b) {
              cout<<" betacf(b, a, 1.0 - x) = "<< betacf(b, a, 1.0 - x)<<" a = "<<a<<" b = "<<b<<endl;
              }
              =========================================================*/
            return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
        }
    }
    double betacf(const double a, const double b, const double x) {
        //Evaluates continued fraction for incomplete beta function by modified Lentzâs method (5.2). User should not call directly.
        int m, m2;
        double aa, c, d, del, h, qab, qam, qap;
        qab = a + b; //These qab will be used in factors that occur in the coefficients (6.4.6).
        qap = a + 1.0;
        qam = a - 1.0;
        c = 1.0; //First step of Lentzâs method.

        d = 1.0 - qab * x / qap;
        /*========debug====================
          if (a == b) {
          cout<<"d = "<<d<<" qap = "<<qap<<" qab = "<<qab<<" x = "<<x<<endl;
          }
          ========debug====================*/
        if(fabs(d) < FPMIN) d = FPMIN;
        d = 1.0 / d;
        h = d;
        // cout<<"h = "<<h<<" FPMIN = "<<FPMIN<<endl;
        for(m = 1; m < 10000; m++) {
            m2 = 2 * m;
            aa = m * (b - m) * x / ((qam + m2) * (a + m2));
            d = 1.0 + aa * d; //One step (the even one) of the recurrence.
            if(fabs(d) < FPMIN)
                d = FPMIN;
            c = 1.0 + aa / c;
            if(fabs(c) < FPMIN)
                c = FPMIN;
            d = 1.0 / d;
            h *= d * c;
            aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
            d = 1.0 + aa * d; //Next step of the recurrence (the odd one).
            if(fabs(d) < FPMIN)
                d = FPMIN;
            c = 1.0 + aa / c;
            if(fabs(c) < FPMIN)
                c = FPMIN;
            d = 1.0 / d;
            del = d * c;
            h *= del;
            if(fabs(del - 1.0) <= EPS)
                break;
        }
        return h;
    }
    double betaiapprox(double a, double b, double x) {
        //Incomplete beta by quadrature. Returns Ix.a; b/. User should not call directly.
        const double y[18] = {0.0021695375159141994,
                              0.011413521097787704, 0.027972308950302116, 0.051727015600492421,
                              0.082502225484340941, 0.12007019910960293, 0.16415283300752470,
                              0.21442376986779355, 0.27051082840644336, 0.33199876341447887,
                              0.39843234186401943, 0.46931971407375483, 0.54413605556657973,
                              0.62232745288031077, 0.70331500465597174, 0.78649910768313447,
                              0.87126389619061517, 0.95698180152629142
                             };
        const double w[18] = {0.0055657196642445571,
                              0.012915947284065419, 0.020181515297735382, 0.027298621498568734,
                              0.034213810770299537, 0.040875750923643261, 0.047235083490265582,
                              0.053244713977759692, 0.058860144245324798, 0.064039797355015485,
                              0.068745323835736408, 0.072941885005653087, 0.076598410645870640,
                              0.079687828912071670, 0.082187266704339706, 0.084078218979661945,
                              0.085346685739338721, 0.085983275670394821
                             };
        int j;
        double xu, t, sum, ans;
        double a1 = a - 1.0, b1 = b - 1.0, mu = a / (a + b);
        double lnmu = log(mu), lnmuc = log(1. - mu);
        t = sqrt(a * b / (sqrt(a + b) * (a + b + 1.0)));
        if(x > a / (a + b)) { //Set how far to integrate into the tail:
            if(x >= 1.0)
                return 1.0;
            xu = min_Gamma(1., max_Gamma(mu + 10.*t, x + 5.0 * t));
        } else {

            if(x <= 0.0)
                return 0.0;
            xu = max_Gamma(0., min_Gamma(mu - 10.*t, x - 5.0 * t));
        }
        sum = 0;
        for(j = 0; j < 18; j++) { //Gauss-Legendre.
            t = x + (xu - x) * y[j];
            sum += w[j] * exp(a1 * (log(t) - lnmu) + b1 * (log(1 - t) - lnmuc));
        }
        ans = sum * (xu - x) * exp(a1 * lnmu - gammln(a) + b1 * lnmuc - gammln(b) + gammln(a + b));
        return ans > 0.0 ? 1.0 - ans : -ans;
    }
    double invbetai(double p, double a, double b) {
        // Inverse of incomplete beta function. Returns x such that Ix.a; b/ D p for argument p between 0 and 1.
        const double EPS = 1.e-8;
        double pp, t, u, err, x, al, h, w, afac, a1 = a - 1., b1 = b - 1.;
        int j;

        if(p <= 0.)
            return 0.;
        else if(p >= 1.)
            return 1.;
        else if(a >= 1. && b >= 1.) { // Set initial guess. See text.
            pp = (p < 0.5) ? p : 1. - p;
            t = sqrt(-2.*log(pp));
            x = (2.30753 + t * 0.27061) / (1. + t * (0.99229 + t * 0.04481)) - t;
            /*========debug====================
              if (a==b) {

              cout<<"x = "<<x<<" p = "<<p<<" w = "<<w<<" t = "<<t<<" al = "<<al<<" h = "<<h<<endl;

              }
              ========debug====================*/
            // if(p < 0.5)//origin code from numerical ricipes.
            if(x < 0.0)//fixed bug.
                x = -x;
            al = (sqrt(x) - 3.) / 6.;
            h = 2. / (1. / (2.*a - 1.) + 1. / (2.*b - 1.));
            w = (x * sqrt(al + h) / h) - (1. / (2.*b - 1) - 1. / (2.*a - 1.)) * (al + 5. / 6. - 2. / (3.*h));
            x = a / (a + b * exp(2.*w));

        } else {
            double lna = log(a / (a + b)), lnb = log(b / (a + b));
            t = exp(a * lna) / a;
            u = exp(b * lnb) / b;
            w = t + u;
            if(p < t / w)
                x = pow(a * w * p, 1. / a);
            else
                x = 1. - pow(b * w * (1. - p), 1. / b);
        }
        afac = -gammln(a) - gammln(b) + gammln(a + b);
        for(j = 0; j < 10; j++) {
            if(x == 0. || x == 1.) // a or b too small for accurate calculation.
                return x;
            err = betai(x, a, b) - p;
            /*=========================debug code=====================
              if (a==b) {
              cout<<"p = "<<p<<" err = "<<err<<" t = "<<t<<" u = "<<u<<" x = "<<x<<endl;
              }
              =========================================================*/
            t = exp(a1 * log(x) + b1 * log(1. - x) + afac);
            u = err / t; //Halley:
            x -= (t = u / (1. - 0.5 * min_Gamma(1., u * (a1 / x - b1 / (1. - x)))));
            /*=========================debug code=====================
              if (a==b) {
              cout<<"err = "<<err<<" t = "<<t<<" u = "<<u<<" x = "<<x<<endl;
              }
              =========================================================*/
            if(x <= 0.)
                x = 0.5 * (x + t); // Bisect if x tries to go neg or > 1.
            if(x >= 1.)
                x = 0.5 * (x + t + 1.);
            if(fabs(t) < EPS * x && j > 0)
                break;
        }
        return x;
    }

    void coordConverter (const Vector_t &TriNormal, Vector_t &x) {
        


    }



};





/*  ==========================================================================*/







#endif //OPAL_SECONDARY_EMISSION_PHYSICS_HH
