// ------------------------------------------------------------------------
//
// Class:BeamStrippingPhysics
// Defines the beam stripping physics models
//
// ------------------------------------------------------------------------
// Class category: ParticleMatterInteraction
// ------------------------------------------------------------------------
// $Date: 2018/11 $
// $Author: PedroCalvo$
//-------------------------------------------------------------------------

#include "Algorithms/PartBunchBase.h"
#include "Algorithms/PartData.h"
#include "Physics/Physics.h"
#include "Solvers/BeamStrippingPhysics.hh"
#include "Structure/LossDataSink.h"
#include "Utilities/LogicalError.h"
#include "Utilities/Options.h"
#include "Utilities/GeneralClassicException.h"
#include "Utilities/Util.h"
#include "Utilities/Timer.h"

#include "Ippl.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <math.h>
#include <boost/math/special_functions/chebyshev.hpp>

using namespace std;

using Physics::kB;
using Physics::q_e;
using Physics::m_e;
using Physics::c;
using Physics::h_bar;
using Physics::alpha;
using Physics::m_hm;
using Physics::m_p;
using Physics::pi;

namespace {
    struct InsideTester {
        virtual ~InsideTester() {}

        virtual bool checkHit(const Vector_t &R) = 0;
        virtual double getPressure(const Vector_t &R) = 0;
    };
    struct BeamStrippingInsideTester: public InsideTester {
        BeamStrippingInsideTester(ElementBase * el) {
            bstp_m = static_cast<BeamStripping*>(el);
        }
        virtual bool checkHit(const Vector_t &R) {
            return bstp_m->checkPoint(R(0), R(1), R(2));
        }
        double getPressure(const Vector_t &R) {
            return bstp_m->checkPressure(R(0)*0.001, R(1)*0.001);
        }
    private:
        BeamStripping *bstp_m;
    };
}


BeamStrippingPhysics::BeamStrippingPhysics(const std::string &name, ElementBase *element):
    ParticleMatterInteractionHandler(name, element),
    T_m(0.0),
    dT_m(0.0),
    mass_m(0.0),
    charge_m(0.0),
    m_h(0.0),
    pressure_m(0.0),
    rhoinit_m(0.0),
    rhofinal_m(0.0),
    NCS_a(0.0),
    NCS_b(0.0),
    NCS_c(0.0),
    NCS_total(0.0),
    bunchToMatStat_m(0),
    stoppedPartStat_m(0),
    rediffusedStat_m(0),
    locPartsInMat_m(0)
{
    bstp_m = dynamic_cast<BeamStripping *>(getElement()->removeWrappers());    
    bstpshape_m = element_ref_m->getType();
    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(element_ref_m->getName(), !Options::asciidump));

    const gsl_rng_type * T;
    gsl_rng_env_setup();
    struct timeval tv; // Seed generation based on time
    gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
    T = gsl_rng_default; // Generator setup
    r_m = gsl_rng_alloc (T);
    gsl_rng_set(r_m, mySeed);

    gas_m = bstp_m->getResidualGas();
}


BeamStrippingPhysics::~BeamStrippingPhysics() {
    lossDs_m->save();
    if(r_m)
        gsl_rng_free(r_m);
}

const std::string BeamStrippingPhysics::getType() const {
    return "BeamStrippingPhysics";
}

void BeamStrippingPhysics::apply(PartBunchBase<double, 3> *bunch,
                                 const std::pair<Vector_t, double> &boundingSphere,
                                 size_t numParticlesInSimulation) {

    Inform m ("BeamStrippingPhysics::apply ", INFORM_ALL_NODES);

    bunchToMatStat_m  = 0;
    rediffusedStat_m  = 0;
    stoppedPartStat_m = 0;
    locPartsInMat_m   = 0;
    bool onlyOneLoopOverParticles = ! (allParticleInMat_m);

    dT_m = bunch->getdT();

    double mass = bunch->getM()*1E-9;

    double uam = 0.9314940954;  // Unified atomic mass unit in GeV
    m_h = 1.00794 * uam;        // Hydrogen atom mass in GeV

    rhoinit_m = 1000 * cycl_m->getMinR();
    rhofinal_m = 1000 * cycl_m->getMaxR();

    do {
        if(abs(mass-m_hm)  < 1E-6 || 
           abs(mass-m_p)   < 1E-6 || 
           abs(mass-2*m_p) < 1E-6 || 
           abs(mass-m_h)   < 1E-6)
            doPhysics(bunch);
        else {
            Inform gmsgALL("OPAL ", INFORM_ALL_NODES);
            gmsgALL << getName() << ": Unsupported type of particle for residual gas interactions!"<< endl;
            gmsgALL << getName() << "-> Beam Stripping Physics not apply"<< endl;
        }
    } while (onlyOneLoopOverParticles == false);
}


void BeamStrippingPhysics::doPhysics(PartBunchBase<double, 3> *bunch) {
    /*
      Do physics if
      -- particle in material
      -- particle not dead (bunch->Bin[i] != -1)
      Delete particle i: bunch->Bin[i] != -1;
    */

    double Eng = 0.0;
    double gamma = 0.0;
    double beta = 0.0;
    double deltas = 0.0;
    double Lpath = bunch->getLPath();
    Vector_t extE = Vector_t(0.0, 0.0, 0.0);
    Vector_t extB = Vector_t(0.0, 0.0, 0.0); //kGauss
    double E = 0.0;
    double B = 0.0;

    bool pdead_GS = false;
    bool pdead_LS = false;
    bool stop = bstp_m->getStop();

    InsideTester *tester;
    tester = new BeamStrippingInsideTester(element_ref_m);

    for (size_t i = 0; i < bunch->getLocalNum(); ++i) {
        if ( (bunch->Bin[i] != -1) && (tester->checkHit(bunch->R[i])) && (Lpath != 0) ) {

            pressure_m = tester->getPressure(bunch->R[i]);
            mass_m = bunch->M[i];
            charge_m = bunch->Q[i];

            Eng = (sqrt(1.0  + dot(bunch->P[i], bunch->P[i])) - 1) * mass_m; //GeV

            gamma = (Eng + mass_m) / mass_m;
            beta = sqrt(1.0 - 1.0 / (gamma * gamma));
            deltas = dT_m * beta * c;

            CrossSection(bunch->R[i], Eng);
            pdead_GS = GasStripping(deltas);

            if(abs(mass_m-m_hm) < 1E-6 && charge_m == -q_e) {
                cycl_m->apply(bunch->R[i]*0.001, bunch->P[i], T_m, extE, extB);
                B = 0.1 * sqrt(extB[0]*extB[0] + extB[1]*extB[1] + extB[2]*extB[2]); //T
                E = gamma * beta * c * B;
                pdead_LS = LorentzStripping(gamma, E);
            }

            if (pdead_GS == true || pdead_LS == true) {
                lossDs_m->addParticle(bunch->R[i], bunch->P[i], bunch->ID[i], bunch->getT()*1e9);
                if (stop) {
                    bunch->Bin[i] = -1;
                    stoppedPartStat_m++;
                    *gmsg << level4 << "* Particle " << bunch->ID[i] << " is deleted by beam stripping" << endl;
                }
                else {
                    bunch->updateNumTotal();
                    *gmsg << level4 << "* Total number of particles after beam stripping = " << bunch->getTotalNum() << endl;
                    SecondaryParticles(bunch, i, pdead_LS);
                }
            }
        }
    }
    delete tester;
}


void BeamStrippingPhysics::CrossSection(const Vector_t &R, double &Eng){

    const double temperature = bstp_m->getTemperature();    // K

    if(pressure_m <= 0.0)
        throw LogicalError("BeamStrippingPhysics::CrossSection()", "Pressure must not be zero");
    if(temperature <= 0.0)
        throw LogicalError("BeamStrippingPhysics::CrossSection()", "Temperature must not be zero");

    Eng *=1E6; //keV
    double Emin = 0.0;
    double Emax = 0.0;
    double Eth = 0.0;
    double a1 = 0.0;
    double a2 = 0.0;
    double a3 = 0.0;
    double a4 = 0.0;
    double a5 = 0.0;
    double a6 = 0.0;
    double a7 = 0.0;
    double a8 = 0.0;
    double a9 = 0.0;
    double a10 = 0.0;
    double a11 = 0.0;
    double a12 = 0.0;
    
    if (gas_m == "H2") {

        totalmolecularDensity_m = 100 * pressure_m / (kB * q_e * temperature);

        double CS_a = 0.0;
        double CS_b = 0.0;
        double CS_c = 0.0;
        double CS_total = 0.0;

        if(abs(mass_m-m_hm) < 1E-6 && charge_m == -q_e) {
            // Single-electron detachment - Hydrogen Production
            Emin = CSCoefSingle_Hminus_Chebyshev[0];
            Emax = CSCoefSingle_Hminus_Chebyshev[1];
            for (int i = 0; i < 9; ++i)
                a_m[i] = {CSCoefSingle_Hminus_Chebyshev[i+2]};
            CS_a = CSChebyshevFitting(Eng*1E3, Emin, Emax);
            
            // Double-electron detachment - Proton Production
            Emin = CSCoefDouble_Hminus_Chebyshev[0];
            Emax = CSCoefDouble_Hminus_Chebyshev[1];
            for (int i = 0; i < 9; ++i)
                a_m[i] = {CSCoefDouble_Hminus_Chebyshev[i+2]};
            CS_b = CSChebyshevFitting(Eng*1E3, Emin, Emax);
        }
        else if(abs(mass_m-m_p) < 1E-6 && charge_m == q_e) {
            // Single-electron capture - Hydrogen Production
            Emin = CSCoefSingle_Hplus_Chebyshev[0];
            Emax = CSCoefSingle_Hplus_Chebyshev[1];
            for (int i = 0; i < 9; ++i)
                a_m[i] = {CSCoefSingle_Hplus_Chebyshev[i+2]};
            CS_a = CSChebyshevFitting(Eng*1E3, Emin, Emax);
        
            // Double-electron capture - Hminus Production
            Emin = CSCoefDouble_Hplus_Chebyshev[0];
            Emax = CSCoefDouble_Hplus_Chebyshev[1];
            for (int i = 0; i < 9; ++i)
                a_m[i] = {CSCoefDouble_Hplus_Chebyshev[i+2]};
            CS_b = CSChebyshevFitting(Eng*1E3, Emin, Emax);
        }
        else if(abs(mass_m-m_h) < 1E-6 && charge_m == 0.0) {
            // Single-electron detachment - Proton Production
            Eth = CSCoefProtonProduction_H_Tabata[0];
            a1 = CSCoefProtonProduction_H_Tabata[1];
            a2 = CSCoefProtonProduction_H_Tabata[2];
            a3 = CSCoefProtonProduction_H_Tabata[3];
            a4 = CSCoefProtonProduction_H_Tabata[4];
            a5 = CSCoefProtonProduction_H_Tabata[5];
            a6 = CSCoefProtonProduction_H_Tabata[6];
            CS_a = CSAnalyticFunctionTabata(Eng, Eth, a1, a2, 0.0, 0.0, a3, a4) +
                a5 * CSAnalyticFunctionTabata(Eng/a6, Eth/a6, a1, a2, 0.0, 0.0, a3, a4);
            
            // Single-electron capture - Hminus Production
            Eth = CSCoefHminusProduction_H_Tabata[0];
            a1 = CSCoefHminusProduction_H_Tabata[1];
            a2 = CSCoefHminusProduction_H_Tabata[2];
            a3 = CSCoefHminusProduction_H_Tabata[3];
            a4 = CSCoefHminusProduction_H_Tabata[4];
            a5 = CSCoefHminusProduction_H_Tabata[5];
            a6 = CSCoefHminusProduction_H_Tabata[6];
            a7 = CSCoefHminusProduction_H_Tabata[7];
            a8 = CSCoefHminusProduction_H_Tabata[8];
            a9 = CSCoefHminusProduction_H_Tabata[9];
            a10 = CSCoefHminusProduction_H_Tabata[10];
            a11 = CSCoefHminusProduction_H_Tabata[11];
            a12 = CSCoefHminusProduction_H_Tabata[12];
            CS_b = CSAnalyticFunctionTabata(Eng, Eth, a1, a2, a3, a4, a5, a6) +
                CSAnalyticFunctionTabata(Eng, Eth, a7, a8, a9, a10, a11, a12);
        }
        else if(abs(mass_m-2*m_p) < 1E-6 && charge_m == q_e) {
            // Proton production
            Eth = CSCoefProtonProduction_H2plus_Tabata[0];
            a1 = CSCoefProtonProduction_H2plus_Tabata[1];
            a2 = CSCoefProtonProduction_H2plus_Tabata[2];
            a3 = CSCoefProtonProduction_H2plus_Tabata[3];
            a4 = CSCoefProtonProduction_H2plus_Tabata[4];
            a5 = CSCoefProtonProduction_H2plus_Tabata[5];
            a6 = CSCoefProtonProduction_H2plus_Tabata[6];
            a7 = CSCoefProtonProduction_H2plus_Tabata[7];
            a8 = CSCoefProtonProduction_H2plus_Tabata[8];
            a9 = CSCoefProtonProduction_H2plus_Tabata[9];
            a10 = CSCoefProtonProduction_H2plus_Tabata[10];
            CS_a = CSAnalyticFunctionTabata(Eng, Eth, a1, a2, 0.0, 0.0, a3, a4) + 
                CSAnalyticFunctionTabata(Eng, Eth, a5, a6, 0.0, 0.0, a7, a8) +
                a9 * CSAnalyticFunctionTabata(Eng/a10, Eth/a10, a5, a6, 0.0, 0.0, a7, a8);
            
            // Hydrogen production
            Emin = CSCoefHydrogenProduction_H2plus_Chebyshev[0];
            Emax = CSCoefHydrogenProduction_H2plus_Chebyshev[1];
            for (int i = 0; i < 9; ++i)
                a_m[i] = {CSCoefHydrogenProduction_H2plus_Chebyshev[i+2]};
            CS_b = CSChebyshevFitting(Eng*1E3, Emin, Emax);
        
            // H3+, H
            Eth = CSCoefH3plusProduction_H2plus_Tabata[0];
            a1 = CSCoefH3plusProduction_H2plus_Tabata[1];
            a2 = CSCoefH3plusProduction_H2plus_Tabata[2];
            a3 = CSCoefH3plusProduction_H2plus_Tabata[3];
            a4 = CSCoefH3plusProduction_H2plus_Tabata[4];
            a5 = CSCoefH3plusProduction_H2plus_Tabata[5];
            a6 = CSCoefH3plusProduction_H2plus_Tabata[6];
            CS_c = CSAnalyticFunctionTabata(Eng, Eth, a1, a2, a3, a4, a5, a6);
        }
        CS_total = CS_a + CS_b + CS_c;
        //*gmsg << "* CS_total = " << CS_total << " cm2" << endl;
    
        NCS_a = CS_a * 1E-4 * totalmolecularDensity_m;
        NCS_b = CS_b * 1E-4 * totalmolecularDensity_m;
        NCS_c = CS_c * 1E-4 * totalmolecularDensity_m;
        NCS_total = CS_total * 1E-4 * totalmolecularDensity_m;
    }
    
    
    else if (gas_m == "AIR") {
        
        static const double fMolarFraction[3] = {
            // Nitrogen
            78.084/100,
            //Oxygen
            20.947/100,
            //Argon
            0.934/100
        };
        double CS_single[3];
        double CS_double[3];
        double CS_total[3];
        double NCS_single[3];
        double NCS_double[3];
        double NCS_all[3];
        double NCS_single_all = 0.0;
        double NCS_double_all = 0.0;
        double NCS_total_all = 0.0;

        for (int i = 0; i < 3; ++i) {

            molecularDensity[i] = 100 * pressure_m * fMolarFraction[i] / (kB * q_e * temperature);

            if(abs(mass_m-m_hm) < 1E-6 && charge_m == -q_e) {
                // Single-electron detachment - Hydrogen Production
                Eth = CSCoefSingle_Hminus[i][0];
                for (int j = 0; j < 8; ++j)
                    b_m[i][j] = CSCoefSingle_Hminus[i][j+1];
                CS_single[i] = CSAnalyticFunctionNakai(Eng, Eth, i);

                // Double-electron detachment - Proton Production
                Eth = CSCoefDouble_Hminus[i][0];
                for (int j = 0; j < 8; ++j)
                    b_m[i][j] = CSCoefDouble_Hminus[i][j+1];
                CS_double[i] = CSAnalyticFunctionNakai(Eng, Eth, i);
            }
            else if(abs(mass_m-m_p) < 1E-6 && charge_m == q_e) {
                // Single-electron capture - Hydrogen Production
                Eth = CSCoefSingle_Hplus[i][0];
                for (int j = 0; j < 8; ++j)
                    b_m[i][j] = CSCoefSingle_Hplus[i][j+1];
                CS_single[i] = CSAnalyticFunctionNakai(Eng, Eth, i) +
                    b_m[i][6] * CSAnalyticFunctionNakai(Eng/b_m[i][7], Eth/b_m[i][7], i);

                // Double-electron capture - Hminus Production
                Eth = CSCoefDouble_Hplus[i][0];
                for (int j = 0; j < 8; ++j)
                    b_m[i][j] = CSCoefDouble_Hplus[i][j+1];
                if(b_m[i][7] != 0) {
                    CS_double[i] = CSAnalyticFunctionNakai(Eng, Eth, i) +
                        b_m[i][6] * CSAnalyticFunctionNakai(Eng/b_m[i][7], Eth/b_m[i][7], i);
                }
                else {
                    CS_double[i] = CSAnalyticFunctionNakai(Eng, Eth, i);
                }
            }
            else if(abs(mass_m-m_h) < 1E-6 && charge_m == 0.0) {
                // Single-electron detachment - Proton Production
                Eth = CSCoefSingleLoss_H[i][0];
                for (int j = 0; j < 8; ++j)
                    b_m[i][j] = CSCoefSingleLoss_H[i][j+1];
                CS_single[i] = CSAnalyticFunctionNakai(Eng, Eth, i);

                // Single-electron capture - Hminus Production
                Eth = CSCoefSingleCapt_H[i][0];
                for (int j = 0; j < 8; ++j)
                    b_m[i][j] = CSCoefSingleCapt_H[i][j+1];
                if(b_m[i][7] != 0) {
                    CS_double[i] = CSAnalyticFunctionNakai(Eng, Eth, i) +
                        b_m[i][6] * CSAnalyticFunctionNakai(Eng/b_m[i][7], Eth/b_m[i][7], i);
                }
                else {
                    CS_double[i] = CSAnalyticFunctionNakai(Eng, Eth, i);
                }
            }
            else {
                CS_single[i] = {0.0};
                CS_double[i] = {0.0};
            }
            CS_total[i] = CS_single[i] + CS_double[i];
            //*gmsg << "* CS_total[" << i <<"]        = " << CS_total[i] << " cm2" << endl;

            NCS_single[i] = CS_single[i] * 1E-4 * molecularDensity[i];
            NCS_double[i] = CS_double[i] * 1E-4 * molecularDensity[i];
            NCS_all[i] = CS_total[i] * 1E-4 * molecularDensity[i];

            NCS_single_all += NCS_single[i];
            NCS_double_all += NCS_double[i];
            NCS_total_all += NCS_all[i];
        }
        NCS_a = NCS_single_all;
        NCS_b = NCS_double_all;
        NCS_total = NCS_total_all;
    }
    else {
        throw GeneralClassicException("BeamStrippingPhysics::CrossSection", "Residual gas not found ...");
    }
    //*gmsg << "* NCS_a      = " << NCS_a << endl;
    //*gmsg << "* NCS_b      = " << NCS_b << endl;
    //*gmsg << "* NCS_total       = " << NCS_total << endl;
    //*gmsg << "**************** " << endl;
}


double BeamStrippingPhysics::CSAnalyticFunctionNakai(double Eng, double Eth, int &i) {

    //const double eps2 = epsilon_0*epsilon_0;
    //const double h = 2*pi*h_bar*1E9*q_e; //Kg m2 s-1
    //const double E_inf = alpha2 * m_e / ( 4*pi*h_bar*c ); //m-1
    //const double E_R = m_h*1e6 * pow(q_e,4) / ( 8*c*c*eps2*pow(h,2) ); //keV

    const double E_ryd = alpha*alpha * m_e*1e6 * 0.5; //keV
    const double E_R = E_ryd * m_h / m_e; //keV
    const double sigma_0 = 1E-16;
    double E1 = 0.0;
    double f1 = 0.0;
    double sigma = 0.0;
    if(Eng > Eth) {
        E1 = (Eng-Eth);
        f1 = sigma_0 * b_m[i][0] * pow((E1/E_R),b_m[i][1]);
        if(b_m[i][2] != 0.0 && b_m[i][3] != 0.0)
            sigma = f1 / (1 + pow((E1/b_m[i][2]),(b_m[i][1]+b_m[i][3])) + pow((E1/b_m[i][4]),(b_m[i][1]+b_m[i][5])));
        else
            sigma = f1 / (1 + pow((E1/b_m[i][4]),(b_m[i][1]+b_m[i][5])));
    }
    return sigma;
}

double BeamStrippingPhysics::CSAnalyticFunctionTabata(double Eng, double Eth,
                                                double a1, double a2, double a3, double a4, double a5, double a6) {
    
    const double E_ryd = alpha*alpha * m_e*1e6 * 0.5; //keV
    const double sigma_0 = 1E-16;
    double E1 = 0.0;
    double f1 = 0.0;
    double sigma = 0.0;
    if(Eng > Eth) {
        E1 = (Eng-Eth);
        f1 = sigma_0 * a1 * pow((E1/E_ryd),a2);
        if(a3 != 0.0 && a4 != 0.0)
            sigma = f1 / (1 + pow((E1/a3),(a2+a4)) + pow((E1/a5),(a2+a6)));
        else
            sigma = f1 / (1 + pow((E1/a5),(a2+a6)));
    }
    return sigma;
}

double BeamStrippingPhysics::CSChebyshevFitting(double Eng, double Emin, double Emax) {
    
    double sum_aT = 0.0;
    double sigma = 0.0;
    double aT[8] = {0.0};
    if(Eng > Emin && Eng < Emax) {
        double X = ((log(Eng)-log(Emin)) - (log(Emax)-log(Eng))) / (log(Emax)- log(Emin));
        for (int i = 0; i < 8; ++i) {
            aT[i] = (a_m[i+1] * boost::math::chebyshev_t(i+1, X));
            sum_aT += aT[i];
        }
        sigma = exp(0.5*a_m[0] + sum_aT);
    }
    return sigma;
}


bool BeamStrippingPhysics::GasStripping(double &deltas) {

    double xi = gsl_rng_uniform(r_m);
    double fg = 1-exp(-NCS_total*deltas);

    return (fg >= xi);
}

bool BeamStrippingPhysics::LorentzStripping(double &gamma, double &E) {

    double tau = 0.0;
    double fL = 0.0;
    double xi = gsl_rng_uniform(r_m);

//      //Parametrization
//      const double A1 = 3.073E-6;
//      const double A2 = 4.414E9;
//      tau = (A1/E) * exp(A2/E);
//      fL = 1 - exp( - dT_m / (gamma * tau));

    //Theoretical
    const double eps0 = 0.75419 * q_e;
    const double hbar = h_bar*1E9*q_e;
    const double me = m_e*1E9*q_e/(c*c);
    const double a0 = hbar / (me * c * alpha);
    const double p = 0.0126;
    const double S0 = 0.783;
    const double a = 2.01407/a0;
    const double k0 = sqrt(2 * me * eps0)/hbar;
    const double N = (sqrt(2 * k0 * (k0+a) * (2*k0+a)))/a;
    double zT = eps0 / (q_e * E);
    tau = (4 * me * zT)/(S0 * N * N * hbar * (1+p)*(1+p) * (1-1/(2*k0*zT))) * exp(4*k0*zT/3);
    fL = 1 - exp( - dT_m / (gamma * tau));

    return (fL > xi);
}

void BeamStrippingPhysics::SecondaryParticles(PartBunchBase<double, 3> *bunch, size_t &i, bool pdead_LS) {

    double opyield_m = 1.0;
    size_t tempnum = bunch->getLocalNum();
    size_t count = 0;
    double r = gsl_rng_uniform(r_m);

    // change the mass_m and charge_m
    if(abs(mass_m-m_hm) < 1E-6 && charge_m == -q_e) {
        if (pdead_LS == true) {
            TransformToHydrogen(bunch, i);
        }
        else {
            if(r > NCS_b/NCS_total)
                TransformToHydrogen(bunch, i);
            else
                TransformToProton(bunch, i);
        }
    }

    else if(abs(mass_m-m_p) < 1E-6 && charge_m == q_e) {
        if(r > NCS_b/NCS_total)
            TransformToHydrogen(bunch, i);
        else
            TransformToHminus(bunch, i);
    }
    
    else if(abs(mass_m-m_h) < 1E-6 && charge_m == 0.0) {
        if(r > NCS_b/NCS_total)
            TransformToProton(bunch, i);
        else
            TransformToHminus(bunch, i);
    }

    else if(abs(mass_m-2*m_p) < 1E-6 && charge_m == q_e) {
        if(NCS_c>NCS_b && NCS_b>NCS_a){
            if(r > (NCS_a+NCS_b)/NCS_total)
                TransformToH3plus(bunch, i);
            else if(r > NCS_a/NCS_total && r < (NCS_a+NCS_b)/NCS_total)
                TransformToHydrogen(bunch, i);
            else if(r < NCS_a/NCS_total)
                TransformToProton(bunch, i);
        }
        else if(NCS_a>NCS_b && NCS_b>NCS_c) {
            if(r > (NCS_c+NCS_b)/NCS_total)
                TransformToProton(bunch, i);
            else if(r > NCS_c/NCS_total && r < (NCS_c+NCS_b)/NCS_total)
                TransformToHydrogen(bunch, i);
            else if(r < NCS_c/NCS_total)
                TransformToH3plus(bunch, i);
        }
        else if(NCS_a>NCS_b && NCS_c>NCS_a) {
            if(r > (NCS_a+NCS_b)/NCS_total)
                TransformToH3plus(bunch, i);
            else if(r > NCS_b/NCS_total && r < (NCS_a+NCS_b)/NCS_total)
                TransformToProton(bunch, i);
            else if(r < NCS_b/NCS_total)
                TransformToHydrogen(bunch, i);
        }
        else if(NCS_a>NCS_c && NCS_c>NCS_b) {
            if(r > (NCS_c+NCS_b)/NCS_total)
                TransformToProton(bunch, i);
            else if(r > NCS_b/NCS_total && r < (NCS_c+NCS_b)/NCS_total)
                TransformToH3plus(bunch, i);
            else if(r < NCS_b/NCS_total)
                TransformToHydrogen(bunch, i);
        }
        else if(NCS_b>NCS_c && NCS_b>NCS_a && NCS_a>NCS_c) {
            if(r > (NCS_c+NCS_a)/NCS_total)
                TransformToHydrogen(bunch, i);
            else if(r > NCS_c/NCS_total && r < (NCS_c+NCS_a)/NCS_total)
                TransformToProton(bunch, i);
            else if(r < NCS_c/NCS_total)
                TransformToH3plus(bunch, i);
        }
        else {
            if(r > (NCS_c+NCS_a)/NCS_total)
                TransformToHydrogen(bunch, i);
            else if(r > NCS_a/NCS_total && r < (NCS_c+NCS_a)/NCS_total)
                TransformToH3plus(bunch, i);
            else if(r < NCS_a/NCS_total)
                TransformToProton(bunch, i);
        }
    }

    bunch->PType[i] = ParticleType::SECONDARY;
    int j = 1;
    //create new particles
    while (j < opyield_m){
        bunch->create(1);
        bunch->R[tempnum+count] = bunch->R[i];
        bunch->P[tempnum+count] = bunch->P[i];
        bunch->Q[tempnum+count] = bunch->Q[i];
        bunch->M[tempnum+count] = bunch->M[i];
        bunch->PType[tempnum+count] = ParticleType::SECONDARY;
        count++;
        j++;
    }
    if(bunch->weHaveBins())
        bunch->Bin[bunch->getLocalNum()-1] = bunch->Bin[i];
}

void BeamStrippingPhysics::TransformToProton(PartBunchBase<double, 3> *bunch, size_t &i) {
    INFOMSG(level4 << "* Particle " << bunch->ID[i] << " is transformed to proton" << endl;);
    bunch->M[i] = m_p;
    bunch->Q[i] = q_e;
}
void BeamStrippingPhysics::TransformToHydrogen(PartBunchBase<double, 3> *bunch, size_t &i) {
    INFOMSG(level4 << "* Particle " << bunch->ID[i] << " is transformed to neutral hydrogen" << endl;);
    bunch->M[i] = m_h;
    bunch->Q[i] = 0.0;
}
void BeamStrippingPhysics::TransformToHminus(PartBunchBase<double, 3> *bunch, size_t &i) {
    INFOMSG(level4 << "* Particle " << bunch->ID[i] << " is transformed to negative hydrogen ion" << endl;);
    bunch->M[i] = m_hm;
    bunch->Q[i] = -q_e;
}
void BeamStrippingPhysics::TransformToH3plus(PartBunchBase<double, 3> *bunch, size_t &i) {
    INFOMSG(level4 << "* Particle " << bunch->ID[i] << " is transformed to H3+" << endl;);  //c
    bunch->M[i] = 3*m_p;
    bunch->Q[i] = q_e;
}

void BeamStrippingPhysics::print(Inform& msg) {
}

bool BeamStrippingPhysics::stillActive() {
    return locPartsInMat_m != 0;
}

bool BeamStrippingPhysics::stillAlive(PartBunchBase<double, 3> *bunch) {
    bool beamstrippingAlive = true;
    return beamstrippingAlive;
}


/*
    Cross sections parameters for interaction with air
    -- [1] -> Nitrogen
    -- [2] -> Oxygen
    -- [3] -> Argon
*/
const double BeamStrippingPhysics::CSCoefSingle_Hminus[3][9] = {
    {7.50E-04, 4.38E+02, 7.28E-01, 8.40E-01, 2.82E-01, 4.10E+01, 1.37E+00, 0.00E+00, 0.00E+00},
    {-2.00E-04, 3.45E+02, 4.80E-01, 5.30E-02, 8.40E-02, 1.00E+01, 9.67E-01, 0.00E+00, 0.00E+00},
    {1.70E-3, 2.47E+01, 3.36E-01, 7.00E+01, 5.00E-01, 9.90E+01, 7.80E-01, 0.00E+00, 0.00E+00}
};
const double BeamStrippingPhysics::CSCoefDouble_Hminus[3][9] = {
    {1.40E-02, 1.77E+00, 4.80E-01, 0.00E+00, 0.00E+00, 1.52E+02, 1.52E+00, 0.00E+00, 0.00E+00},
    {1.30E-02, 1.90E+00, 6.20E-01, 0.00E+00, 0.00E+00, 5.20E+01, 9.93E-01, 0.00E+00, 0.00E+00},
    {1.50E-02, 1.97E+00, 8.90E-01, 0.00E+00, 0.00E+00, 5.10E+01, 9.37E-01, 0.00E+00, 0.00E+00}
};
const double BeamStrippingPhysics::CSCoefSingle_Hplus[3][9] = {
    {2.00E-03, 1.93E+03, 1.64E+00, 1.98E+00, 6.69E-01, 2.19E+01, 4.15E+00, 3.23E-04, 1.00E+01},
    {-1.50E-03, 3.86E+05, 1.60E+00, 6.93E-02, 3.28E-01, 7.86E+00, 3.92E+00, 3.20E-04, 1.00E+01},
    {2.18E-03, 1.61E+04, 2.12E+00, 1.16E+00, 4.44E-01, 1.39E+01, 4.07E+00, 2.99E-04, 1.45E+01}
};
const double BeamStrippingPhysics::CSCoefDouble_Hplus[3][9] = {
    {2.90E-02, 2.90E-01, 1.50E+00, 1.39E+01, 1.65E+00, 3.94E+01, 5.79E+00, 0.00E+00, 0.00E+00},
    {2.20E-02, 2.64E-01, 1.50E+00, 1.40E+01, 1.70E+00, 4.00E+01, 5.80E+00, 0.00E+00, 0.00E+00},
    {2.90E-02, 1.40E-01, 1.87E+00, 2.40E+01, 1.60E+00, 3.37E+01, 4.93E+00, 4.50E-01, 2.00E-01}
};
const double BeamStrippingPhysics::CSCoefSingleLoss_H[3][9] = {
    {1.36E-02, 2.59E+03, 1.78E+00, 2.69E-01, -3.52E-01, 5.60E+00, 8.99E-01, 0.00E+00, 0.00E+00},
    {1.36E-02, 3.29E+03, 1.85E+00, 2.89E-01, -2.72E-01, 8.20E+00, 1.09E+00, 0.00E+00, 0.00E+00},
    {1.36E-02, 1.16E+12, 7.40E+00, 5.36E-01, -5.42E-01, 1.23E+00, 7.26E-01, 0.00E+00, 0.00E+00}
};
const double BeamStrippingPhysics::CSCoefSingleCapt_H[3][9] = {
    {1.48E-02, 1.15E+00, 1.18E+00, 1.15E+01, 1.19E+00, 3.88E+01, 3.38E+00, 1.00E-01, 2.82E-02},
    {1.13E-02, 3.26E+00, 1.02E+00, 2.99E+00, 4.50E-01, 1.90E+01, 3.42E+00, 0.00E+00, 0.00E+00},
    {1.50E-02, 1.24E+03, 3.38E+00, 2.46E+00, 5.20E-01, 7.00E+00, 2.56E+00, 9.10E-02, 1.95E-02}
};

// Cross sections parameters for interaction with hydrogen gas (H2)
const double BeamStrippingPhysics::CSCoefHminusProduction_H_Tabata[13] = {
    2.1E-02, 9.73E-03, 2.38E+00, 1.39E-02, -5.51E-01, 7.7E-02, 2.12E+00, 1.97E-06, 2.051E+00, 5.5E+00, 6.62E-01, 2.02E+01, 3.62E+00
};
const double BeamStrippingPhysics::CSCoefProtonProduction_H_Tabata[9] = {
    2.0E-02, 2.53E-04, 1.728E+00, 2.164E+00, 7.74E-01, 1.639E+00, 1.43E+01, 0, 0
};
const double BeamStrippingPhysics::CSCoefProtonProduction_H2plus_Tabata[11] = {
    5.0E-03, 6.34E+01, 1.78E+00, 1.38E-03, 4.06E-01, 1.63E-01, 3.27E-01, 1.554E+01, 3.903E+00, 1.735, 1.02E+01
};
const double BeamStrippingPhysics::CSCoefH3plusProduction_H2plus_Tabata[7] = {
    0.00E+00, 6.05E+00, -5.247E-01, 4.088E-03, 2.872E+00, 7.3E-03, 6.99E+00
};
const double BeamStrippingPhysics::CSCoefSingle_Hminus_Chebyshev[11] = {
    2.3, 1.7E+07, -73.1505813599, -1.7569509745, -2.0016760826, -0.1902804971, 0.0171353221, 0.1270833164, -0.1523126215, 0, 0
};
const double BeamStrippingPhysics::CSCoefDouble_Hminus_Chebyshev[11] = {
    1.00E+03, 1.00E+07, -79.0158996582, -2.1025185585, -1.2407282591, 0.174798578, 0.1062489152, -0.0004342734, -0.0465673618, 0, 0
};
const double BeamStrippingPhysics::CSCoefSingle_Hplus_Chebyshev[11] = {
    2.6E+00, 4.00E+06, -82.5164, -6.70755, -6.10977, -2.6281, 0.709759, 0.639033, 0.10298, 0.26124, -0.263817
};
const double BeamStrippingPhysics::CSCoefDouble_Hplus_Chebyshev[11] = {
    200, 1.00E+06, -95.8165, -7.17049, -7.48288, -1.93034, 0.761153, 0.556689, -0.0542859, -0.270184, -0.0147
};
const double BeamStrippingPhysics::CSCoefHydrogenProduction_H2plus_Chebyshev[11] = {
    2.00E+03, 1.00E+05, -70.670173645, -0.632612288, -0.6065212488, -0.0915143117, -0.0121710282, 0.0168179292, 0.0104796877, 0, 0
};

double BeamStrippingPhysics::a_m[9] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
double BeamStrippingPhysics::b_m[3][9] = {0.0};
