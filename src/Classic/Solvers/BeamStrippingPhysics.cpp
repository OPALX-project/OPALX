// Class:BeamStrippingPhysics
// Defines the beam stripping physics models
// ------------------------------------------------------------------------
// Class category:
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

using namespace std;

using Physics::kB;
using Physics::q_e;
using Physics::m_e;
using Physics::Avo;
using Physics::epsilon_0;
using Physics::c;
using Physics::h_bar;
using Physics::alpha;
using Physics::m_hm;
using Physics::m_p;
using Physics::z_p;

namespace {
    struct InsideTester {
        virtual ~InsideTester() {}

        virtual bool checkHit(const Vector_t &R) = 0;
    };
    struct BeamStrippingInsideTester: public InsideTester {
        BeamStrippingInsideTester(ElementBase * el) {
            bstp_m = static_cast<BeamStripping*>(el);
        }
        virtual bool checkHit(const Vector_t &R) {
            return bstp_m->checkPoint(R(0), R(1), R(2));
        }
    private:
        BeamStripping *bstp_m;
    };
}


BeamStrippingPhysics::BeamStrippingPhysics(const std::string &name, ElementBase *element, std::string &material):
    ParticleMatterInteractionHandler(name, element),
    material_m(material),
    T_m(0.0),
    dT_m(0.0),
    mass_m(0.0),
    charge_m(0.0),
    m_h(0.0),
    NbComponents(0.0),
    NCS_single_all(0.0),
    NCS_double_all(0.0),
    NCS_total_all(0.0),
    bunchToMatStat_m(0),
    stoppedPartStat_m(0),
    rediffusedStat_m(0),
    locPartsInMat_m(0)
{
    bstp_m = dynamic_cast<BeamStripping *>(getElement()->removeWrappers());
    NbComponents = 3;
    Material();
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
    m_h = 1.00794 * uam;                // Hydrogen atom mass in GeV

    do {
        if(mass-m_hm < 1E-6 || mass-m_p < 1E-6 || mass-m_h < 1E-6)
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
    Vector_t extB = Vector_t(0.0, 0.0, 0.0); //kG
    double E = 0.0;
    double B = 0.0;

    bool pdead_GS = false;
    bool pdead_LS = false;

    bool stop = bstp_m->getStop();

    InsideTester *tester;
    tester = new BeamStrippingInsideTester(element_ref_m);

    for (size_t i = 0; i < bunch->getLocalNum(); ++i) {
        if ( (bunch->Bin[i] != -1) && (tester->checkHit(bunch->R[i])) && (Lpath != 0) ) {

            mass_m = bunch->M[i];
            charge_m = bunch->Q[i];
//              *gmsg << "* bunch->M[" << bunch->ID[i] << "] = " << bunch->M[i]
//                                << "  bunch->Q[" << bunch->ID[i] << "] = " << bunch->Q[i] << endl;

            Eng = (sqrt(1.0  + dot(bunch->P[i], bunch->P[i])) - 1) * mass_m; //GeV
//              *gmsg << "* Energy =  " << Eng  << endl;

            gamma = (Eng + mass_m) / mass_m;
            beta = sqrt(1.0 - 1.0 / (gamma * gamma));
            deltas = dT_m * beta * c;
//              *gmsg << "* deltas =  " << deltas  << endl;

            CrossSection(Eng);

            pdead_GS = GasStripping(deltas);

            if(mass_m-m_hm < 1E-6 && charge_m == -q_e) {
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
                    INFOMSG(level2 << "* Particle " << bunch->ID[i] << " is deleted by beam stripping" << endl;);
                }
                else {
                    bunch->updateNumTotal();
                    INFOMSG(level2 << "* Total number of particles after beam stripping = " << bunch->getTotalNum() << endl;);
                    SecondaryParticles(bunch, i);
                }
            }
        }
    }
    delete tester;
}


void BeamStrippingPhysics::MolecularDensity(const double &pressure, const double &temperature, int &iComp) {
    if(pressure == 0.0)
        throw LogicalError("BeamStrippingPhysics::setPressure()", "Pressure must not be zero");
    if(temperature == 0.0)
        throw LogicalError("BeamStrippingPhysics::setTemperature()", "Temperature must not be zero");

//    double TotalmolecularDensity_m = 100 * pressure / (kB * q_e * temperature);
    molecularDensity[iComp] = 100 * pressure * fMolarFraction[iComp] / (kB * q_e * temperature);
//    *gmsg << "* Molecular Density     = " << iComp << "  " << molecularDensity[iComp]  << endl;
}

void BeamStrippingPhysics::CrossSection(double &Eng){

    Eng *=1E6; //keV

    double Eth = 0.0;
    double a1 = 0.0;
    double a2 = 0.0;
    double a3 = 0.0;
    double a4 = 0.0;
    double a5 = 0.0;
    double a6 = 0.0;
    double a7 = 0.0;
    double a8 = 0.0;

    double CS_single[3];
    double CS_double[3];
    double CS_total[3];

    double NCS_single[3];
    double NCS_double[3];
    double NCS_total[3];

    NCS_single_all = 0.0;
    NCS_double_all = 0.0;
    NCS_total_all = 0.0;

    for (int i = 0; i < NbComponents; ++i) {

        if(mass_m-m_hm < 1E-6 && charge_m == -q_e) {
            // Single-electron detachment
            Eth = CSCoefSingle_Hminus[i][0];
            a1 = CSCoefSingle_Hminus[i][1];
            a2 = CSCoefSingle_Hminus[i][2];
            a3 = CSCoefSingle_Hminus[i][3];
            a4 = CSCoefSingle_Hminus[i][4];
            a5 = CSCoefSingle_Hminus[i][5];
            a6 = CSCoefSingle_Hminus[i][6];
            CS_single[i] = CSAnalyticFunction(Eng, Eth, a1, a2, a3, a4, a5, a6);

            // Double-electron detachment
            Eth = CSCoefDouble_Hminus[i][0];
            a1 = CSCoefDouble_Hminus[i][1];
            a2 = CSCoefDouble_Hminus[i][2];
            a3 = CSCoefDouble_Hminus[i][3];
            a4 = CSCoefDouble_Hminus[i][4];
            a5 = CSCoefDouble_Hminus[i][5];
            a6 = CSCoefDouble_Hminus[i][6];
            CS_double[i] = CSAnalyticFunction(Eng, Eth, a1, a2, a3, a4, a5, a6);
        }
        else if(mass_m-m_p < 1E-6 && charge_m == q_e) {
            // Single-electron capture
            Eth = CSCoefSingle_Hplus[i][0];
            a1 = CSCoefSingle_Hplus[i][1];
            a2 = CSCoefSingle_Hplus[i][2];
            a3 = CSCoefSingle_Hplus[i][3];
            a4 = CSCoefSingle_Hplus[i][4];
            a5 = CSCoefSingle_Hplus[i][5];
            a6 = CSCoefSingle_Hplus[i][6];
            a7 = CSCoefSingle_Hplus[i][7];
            a8 = CSCoefSingle_Hplus[i][8];
            CS_single[i] = CSAnalyticFunction(Eng, Eth, a1, a2, a3, a4, a5, a6) +
                a7 * CSAnalyticFunction(Eng/a8, Eth, a1, a2, a3, a4, a5, a6);

            // Double-electron capture
            Eth = CSCoefDouble_Hplus[i][0];
            a1 = CSCoefDouble_Hplus[i][1];
            a2 = CSCoefDouble_Hplus[i][2];
            a3 = CSCoefDouble_Hplus[i][3];
            a4 = CSCoefDouble_Hplus[i][4];
            a5 = CSCoefDouble_Hplus[i][5];
            a6 = CSCoefDouble_Hplus[i][6];
            a7 = CSCoefDouble_Hplus[i][7];
            a8 = CSCoefDouble_Hplus[i][8];
            if(a8 != 0) {
                CS_double[i] = CSAnalyticFunction(Eng, Eth, a1, a2, a3, a4, a5, a6) +
                    a7 * CSAnalyticFunction(Eng, Eth/a8, a1, a2, a3, a4, a5, a6);
            }
            else {
                CS_double[i] = CSAnalyticFunction(Eng, Eth, a1, a2, a3, a4, a5, a6);
            }
        }
        else if(mass_m-m_h < 1E-6 && charge_m == 0.0) {
            // Single-electron detachment
            Eth = CSCoefSingleLoss_H[i][0];
            a1 = CSCoefSingleLoss_H[i][1];
            a2 = CSCoefSingleLoss_H[i][2];
            a3 = CSCoefSingleLoss_H[i][3];
            a4 = CSCoefSingleLoss_H[i][4];
            a5 = CSCoefSingleLoss_H[i][5];
            a6 = CSCoefSingleLoss_H[i][6];
            CS_single[i] = CSAnalyticFunction(Eng, Eth, a1, a2, a3, a4, a5, a6);

            // Single-electron capture
            Eth = CSCoefSingleCapt_H[i][0];
            a1 = CSCoefSingleCapt_H[i][1];
            a2 = CSCoefSingleCapt_H[i][2];
            a3 = CSCoefSingleCapt_H[i][3];
            a4 = CSCoefSingleCapt_H[i][4];
            a5 = CSCoefSingleCapt_H[i][5];
            a6 = CSCoefSingleCapt_H[i][6];
            a7 = CSCoefSingleCapt_H[i][7];
            a8 = CSCoefSingleCapt_H[i][8];
            if(a8 != 0) {
                CS_double[i] = CSAnalyticFunction(Eng, Eth, a1, a2, a3, a4, a5, a6) +
                    a7 * CSAnalyticFunction(Eng, Eth/a8, a1, a2, a3, a4, a5, a6);
            }
            else {
                CS_double[i] = CSAnalyticFunction(Eng, Eth, a1, a2, a3, a4, a5, a6);
            }
        }
//              else if(mass_m != m_hm || mass_m != m_p || mass_m != m_h) {
        else {
            CS_single[i] = {0.0};
            CS_double[i] = {0.0};
        }

        CS_total[i] = CS_single[i] + CS_double[i];

//              *gmsg << "* CS_single[" << i <<"]       = " << CS_single[i] << " cm2" << endl;
//              *gmsg << "* CS_double[" << i <<"]       = " << CS_double[i] << " cm2" << endl;
//              *gmsg << "* CS_total[" << i <<"]        = " << CS_total[i] << " cm2" << endl;

        CS_single[i] *= 1E-4;
        CS_double[i] *= 1E-4;
        CS_total[i] *= 1E-4;

        NCS_single[i] = CS_single[i] * molecularDensity[i];
        NCS_double[i] = CS_double[i] * molecularDensity[i];
        NCS_total[i] = CS_total[i] * molecularDensity[i];

//              *gmsg << "* NCS_single[" << i <<"]      = " << NCS_single[i] << endl;
//              *gmsg << "* NCS_double[" << i <<"]      = " << NCS_double[i] << endl;
//              *gmsg << "* NCS_total[" << i <<"]       = " << NCS_total[i]  << endl;

        NCS_single_all += NCS_single[i];
        NCS_double_all += NCS_double[i];
        NCS_total_all += NCS_total[i];
    }
//      *gmsg << "* NCS_single_all      = " << NCS_single_all << endl;
//      *gmsg << "* NCS_double_all      = " << NCS_double_all << endl;
//      *gmsg << "* NCS_total_all       = " << NCS_total_all << endl;
}

double BeamStrippingPhysics::CSAnalyticFunction(double Eng, double Eth,
                                                double a1, double a2, double a3, double a4, double a5, double a6) {

    const double eps2 = epsilon_0*epsilon_0;
    const double h = 2*pi*h_bar*1E9*q_e;
    const double c2 = c*c;
    const double E_R = m_h*1e6 * pow(q_e,4) / ( 8*c2*eps2*pow(h,2) ); //keV
    const double sigma_0 = 1E-16;
    double E1 = 0.0;
    double f1 = 0.0;
    double f = 0.0;
    if(Eng > Eth) {
        E1 = (Eng-Eth);
        f1 = sigma_0 * a1 * pow((E1/E_R),a2);
        if(a3 != 0.0 && a4 != 0.0)
            f = f1 / (1 + pow((E1/a3),(a2+a4)) + pow((E1/a5),(a2+a6)));
        else
            f = f1 / (1 + pow((E1/a5),(a2+a6)));
    }
    return f;
}


bool BeamStrippingPhysics::GasStripping(double &deltas) {

    double xi = gsl_rng_uniform(r_m);

    double fg = 1-exp(-NCS_total_all*deltas);

//      *gmsg << "* deltas = " << deltas << endl;
//      *gmsg << "* 1/lambda = " << NCS_total_all << endl;
//      *gmsg << "* fg = " << fg << "    Random number = " << xi << endl;

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

//      *gmsg << "* fL =    " << fL << endl;

    return (fL > xi);
}

void BeamStrippingPhysics::SecondaryParticles(PartBunchBase<double, 3> *bunch, size_t &i) {

    double opyield_m = 1.0;
    size_t tempnum = bunch->getLocalNum();
    size_t count = 0;

    double r = gsl_rng_uniform(r_m);

//      *gmsg << "* random number = " << r
//                << " NCS_single_all/NCS_total_all = " << NCS_single_all/NCS_total_all
//        << " NCS_double_all/NCS_total_all = " << NCS_double_all/NCS_total_all << endl;

    // change the mass_m and charge_m
    if(mass_m-m_hm < 1E-6 && charge_m == -q_e) {
        if(r > NCS_double_all/NCS_total_all) {
            INFOMSG(level2 << "* Particle " << bunch->ID[i] << " is transformed to neutral hydrogen" << endl;);
            bunch->M[i] = m_h;
            bunch->Q[i] = 0.0;
        }
        else {
            INFOMSG(level2 << "* Particle " << bunch->ID[i] << " is transformed to proton" << endl;);
            bunch->M[i] = m_p;
            bunch->Q[i] = q_e;
        }
    }
    else if(mass_m-m_p < 1E-6 && charge_m == q_e) {
        if(r > NCS_double_all/NCS_total_all) {
            INFOMSG(level2 << "* Particle " << bunch->ID[i] << " is transformed to neutral hydrogen" << endl;);
            bunch->M[i] = m_h;
            bunch->Q[i] = 0.0;
        }
        else {
            INFOMSG(level2 << "* Particle " << bunch->ID[i] << " is transformed to negative hydrogen ion" << endl;);
            bunch->M[i] = m_hm;
            bunch->Q[i] = -q_e;
        }
    }
    else if(mass_m-m_h < 1E-6 && charge_m == 0.0) {
        if(r > NCS_double_all/NCS_total_all) {
            INFOMSG(level2 << "* Particle " << bunch->ID[i] << " is transformed to proton" << endl;);
            bunch->M[i] = m_p;
            bunch->Q[i] = q_e;
        }
        else {
            INFOMSG(level2 << "* Particle " << bunch->ID[i] << " is transformed to negative hydrogen ion" << endl;);
            bunch->M[i] = m_hm;
            bunch->Q[i] = -q_e;
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


void BeamStrippingPhysics::print(Inform& msg) {
}

bool BeamStrippingPhysics::stillActive() {
    return locPartsInMat_m != 0;
}

bool BeamStrippingPhysics::stillAlive(PartBunchBase<double, 3> *bunch) {
    bool beamstrippingAlive = true;
    return beamstrippingAlive;
}


//  ------------------------------------------------------------------------
/// Vacuum for beam stripping
//  ------------------------------------------------------------------------
void  BeamStrippingPhysics::Material() {

    material_m == "VACUUM";

    const double pressure = bstp_m->getPressure();                              // mbar
    const double temperature = bstp_m->getTemperature();                // K

    for(int iComp = 0; iComp<NbComponents; ++iComp) {
        MolecularDensity(pressure, temperature, iComp);
    }
}

const double BeamStrippingPhysics::fMolarFraction[3] = {
    // Nitrogen
    78.084/100,
    //Oxygen
    20.947/100,
    //Argon
    0.934/100
};

const double BeamStrippingPhysics::CSCoefSingle_Hminus[3][7] = {
    // Nitrogen
    {7.50E-04, 4.38E+02, 7.28E-01, 8.40E-01, 2.82E-01, 4.10E+01, 1.37E+00},
    //Oxygen
    {-2.00E-04, 3.45E+02, 4.80E-01, 5.30E-02, 8.40E-02, 1.00E+01, 9.67E-01},
    //Argon
    {1.70E-3, 2.47E+01, 3.36E-01, 7.00E+01, 5.00E-01, 9.90E+01, 7.80E-01}
};

const double BeamStrippingPhysics::CSCoefDouble_Hminus[3][7] = {
    // Nitrogen
    {1.40E-02, 1.77E+00, 4.80E-01, 0.00E+00, 0.00E+00, 1.52E+02, 1.52E+00},
    //Oxygen
    {1.30E-02, 1.90E+00, 6.20E-01, 0.00E+00, 0.00E+00, 5.20E+01, 9.93E-01},
    //Argon
    {1.50E-02, 1.97E+00, 8.90E-01, 0.00E+00, 0.00E+00, 5.10E+01, 9.37E-01}
};

const double BeamStrippingPhysics::CSCoefSingle_Hplus[3][9] = {
    // Nitrogen
    {2.00E-03, 1.93E+03, 1.64E+00, 1.98E+00, 6.69E-01, 2.19E+01, 4.15E+00, 3.23E-04, 1.00E+01},
    //Oxygen
    {-1.50E-03, 3.86E+05, 1.60E+00, 6.93E-02, 3.28E-01, 7.86E+00, 3.92E+00, 3.20E-04, 1.00E+01},
    //Argon
    {2.18E-03, 1.61E+04, 2.12E+00, 1.16E+00, 4.44E-01, 1.39E+01, 4.07E+00, 2.99E-04, 1.45E+01}
};

const double BeamStrippingPhysics::CSCoefDouble_Hplus[3][9] = {
    // Nitrogen
    {2.90E-02, 2.90E-01, 1.50E+00, 1.39E+01, 1.65E+00, 3.94E+01, 5.79E+00, 0.00E+00, 0.00E+00},
    //Oxygen
    {2.20E-02, 2.64E-01, 1.50E+00, 1.40E+01, 1.70E+00, 4.00E+01, 5.80E+00, 0.00E+00, 0.00E+00},
    //Argon
    {2.90E-02, 1.40E-01, 1.87E+00, 2.40E+01, 1.60E+00, 3.37E+01, 4.93E+00, 4.50E-01, 2.00E-01}
};

const double BeamStrippingPhysics::CSCoefSingleLoss_H[3][7] = {
    // Nitrogen
    {1.36E-02, 2.59E+03, 1.78E+00, 2.69E-01, -3.52E-01, 5.60E+00, 8.99E-01},
    //Oxygen
    {1.36E-02, 3.29E+03, 1.85E+00, 2.89E-01, -2.72E-01, 8.20E+00, 1.09E+00},
    //Argon
    {1.36E-02, 1.16E+12, 7.40E+00, 5.36E-01, -5.42E-01, 1.23E+00, 7.26E-01}
};

const double BeamStrippingPhysics::CSCoefSingleCapt_H[3][9] = {
    // Nitrogen
    {1.48E-02, 1.15E+00, 1.18E+00, 1.15E+01, 1.19E+00, 3.88E+01, 3.38E+00, 1.00E-01, 2.82E-02},
    //Oxygen
    {1.13E-02, 3.26E+00, 1.02E+00, 2.99E+00, 4.50E-01, 1.90E+01, 3.42E+00, 0.00E+00, 0.00E+00},
    //Argon
    {1.50E-02, 1.24E+03, 3.38E+00, 2.46E+00, 5.20E-01, 7.00E+00, 2.56E+00, 9.10E-02, 1.95E-02}
};
