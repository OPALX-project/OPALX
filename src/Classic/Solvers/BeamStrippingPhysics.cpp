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
    FN_m = element_ref_m->getName();
    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(FN_m, !Options::asciidump));

	double uam = 0.9314940954;	// Unified atomic mass unit in GeV
	m_h = 1.00794 * uam;		// Hydrogen atom mass in GeV
}


BeamStrippingPhysics::~BeamStrippingPhysics() {
    lossDs_m->save();
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

    dT_m = bunch->getdT();  // own time, maybe larger than in the bunch object
    T_m  = bunch->getT(); 	// dt from bunch

    double mass = bunch->getM()*1E-9;

    do {
    	if(mass-m_hm < 1E-6 || mass-m_p < 1E-6 || mass-m_h < 1E-6)
    		doPhysics(bunch);
        else {
    		Inform gmsgALL("OPAL ", INFORM_ALL_NODES);
    		gmsgALL << getName() << ": Unsupported type of particle for residual gas interactions!"<< endl;
    		gmsgALL << getName() << "-> Beam Stripping Physics not apply"<< endl;
        }
//        lossDs_m->save();
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

    Vector_t extE = Vector_t(0.0, 0.0, 0.0);
    Vector_t extB = Vector_t(0.0, 0.0, 0.0); //kG
	double E = 0.0;
	double B = 0.0;

	bool pdead_GS = false;
	bool pdead_LS = false;

	bool stop = bstp_m->getStop();
	double r = 0.0;
//	r = RandomGenerator();

    InsideTester *tester;
    tester = new BeamStrippingInsideTester(element_ref_m);

    for (size_t i = 0; i < bunch->getLocalNum(); ++i) {
    	if ( (bunch->Bin[i] != -1) && (tester->checkHit(bunch->R[i])) ) {

    		mass_m = bunch->M[i];
    		charge_m = bunch->Q[i];
//    		*gmsg << "* bunch->M[" << bunch->ID[i] << "] = " << bunch->M[i]
//				  << "  bunch->Q[" << bunch->ID[i] << "] = " << bunch->Q[i] << endl;

    		cycl_m->apply(bunch->R[i]*0.001, bunch->P[i], T_m, extE, extB);
    		//			*gmsg << "* extB_m " << extB_m << endl;

    		Eng = (sqrt(1.0  + dot(bunch->P[i], bunch->P[i])) - 1) * mass_m;

    		CrossSection(Eng);

    		Eng *= 1E-9; // Eng GeV
    		gamma = (Eng + mass_m) / mass_m;
    		beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    		deltas = dT_m * beta * c;

    		B = 0.1 * sqrt(extB[0]*extB[0] + extB[1]*extB[1] + extB[2]*extB[2]); //T
    		E = gamma * beta * c * B;
    		//		    *gmsg << "* B =    " << B  << endl;
    		//		    *gmsg << "* E =    " << E  << endl;
    		//		    *gmsg << "* gamma =    " << gamma  << endl;
    		//		    *gmsg << "* beta =    " << beta  << endl;
    		//		    *gmsg << "* deltas =    " << deltas  << endl;

    		pdead_GS = GasStripping(deltas);

//    		if(mass_m-m_hm < 1E-6 && charge_m == -q_e)
    			pdead_LS = LorentzStripping(gamma, E);
    		//			*gmsg << "pdead_GS = " << pdead_GS << endl;
    		//			*gmsg << "pdead_LS = " << pdead_LS << endl;

    		if (pdead_GS == true || pdead_LS == true) {
//    			*gmsg << "* The particle " << bunch->ID[i] << " is stripped in remainder gas" << endl;
    			lossDs_m->addParticle(bunch->R[i], bunch->P[i], bunch->ID[i], T_m);
    			if (stop) {
    				bunch->Bin[i] = -1;
    				stoppedPartStat_m++;
    			}
    			else {
    				bunch->updateNumTotal();
    				*gmsg << "* Total number of particles after beam stripping = " << bunch->getTotalNum() << endl;
    				r = RandomGenerator();
    				SecondaryParticles(bunch, i, r);
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
//    *gmsg << "* Molecular Density	= " << iComp << "  " << molecularDensity[iComp]  << endl;
}

void BeamStrippingPhysics::CrossSection(double &Eng){

	// Analytic function
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
//		else if(mass_m != m_hm || mass_m != m_p || mass_m != m_h) {
		else {
			CS_single[i] = {0.0};
			CS_double[i] = {0.0};
		}

		CS_total[i] = CS_single[i] + CS_double[i];

//		*gmsg << "* CS_single[i] 	= " << CS_single[i] << " cm2" << endl;
//		*gmsg << "* CS_double[i] 	= " << CS_double[i] << " cm2" << endl;
//		*gmsg << "* CS_total[i] 	= " << CS_total[i] << " cm2" << endl;

	    CS_single[i] *= 1E-4;
	    CS_double[i] *= 1E-4;
	    CS_total[i] *= 1E-4;

	    NCS_single[i] = CS_single[i] * molecularDensity[i];
	    NCS_double[i] = CS_double[i] * molecularDensity[i];
	    NCS_total[i] = CS_total[i] * molecularDensity[i];

//		*gmsg << "* NCS_single[i] 	= " << NCS_single[i] << endl;
//		*gmsg << "* NCS_double[i] 	= " << NCS_double[i] << endl;
//		*gmsg << "* NCS_total[i] 	= " << NCS_total[i]  << endl;

	    NCS_single_all += NCS_single[i];
	    NCS_double_all += NCS_double[i];
	    NCS_total_all += NCS_total[i];
	}
//	*gmsg << "* NCS_single_all 	= " << NCS_single_all << endl;
//	*gmsg << "* NCS_double_all 	= " << NCS_double_all << endl;
//	*gmsg << "* NCS_total_all	= " << NCS_total_all << endl;
}

double BeamStrippingPhysics::CSAnalyticFunction(double Eng, double Eth,
		double a1, double a2, double a3, double a4, double a5, double a6) {

	double E_R = 25.00;
	double sigma_0 = 1E-16;
	double E1 = 0.0;
	double f1 = 0.0;
	double f = 0.0;

	if(Eng > Eth) {
		E1 = (Eng-Eth);
		f1 = sigma_0 * a1 * pow((E1/E_R),a2);
		if(a3 != 0.0 && a4 != 0.0){
			f = f1 / (1 + pow((E1/a3),(a2+a4)) + pow((E1/a5),(a2+a6)));
		}
		else {
			f = f1 / (1 + pow((E1/a5),(a2+a6)));
		}
	}
	return f;
}


bool BeamStrippingPhysics::GasStripping(double &deltas) {
	double r = RandomGenerator();
	double fg = 1 - exp(-NCS_total_all * deltas);
//	*gmsg << "* fg = " << fg << "    Random number = " << r << endl;
	return (fg > r);
}


bool BeamStrippingPhysics::LorentzStripping(double &gamma, double &E) {

	double tau = 0.0;
	double fL = 0.0;
	double r = RandomGenerator();

	//Parametrization 1
	const double A1 = 3.073E-6;
	const double A2 = 4.414E9;
	tau = (A1/E) * exp(A2/E);
	fL = 1 - exp( - dT_m / (gamma * tau));

//    	//Parametrization 2
//    	double fL2 = 0.0;
//    	const double aF = 2.653E-6;
//    	const double bF = 4.448E9;
//    	const double eta = 1.4901E-10;
//    	tau = ( aF / (E * (1-eta * E)) ) * exp(bF/E);
//    	fL = 1 - exp( - dT_m / (gamma * tau));
//
//    	//Parametrization 3 (theoretical)
//    	const double eps0 = 0.75419 * q_e;
//    	const double hbar = h_bar*1E9*q_e;
//    	const double me = m_e*1E9*q_e/(c*c);
//    	const double a0 = hbar / (me * c * alpha);
//    	const double p = 0.0126;
//    	const double S0 = 0.783;
//    	const double a = 2.01407/a0;
//    	const double k0 = sqrt(2 * me * eps0)/hbar;
//    	const double N = (sqrt(2 * k0 * (k0+a) * (2*k0+a)))/a;
//    	double zT = eps0 / (q_e * E);
//    	tau = (4 * me * zT)/(S0 * N * N * hbar * (1+p)*(1+p) * (1-1/(2*k0*zT))) * exp(4*k0*zT/3);
//    	fL = 1 - exp( - dT_m / (gamma * tau));

//	*gmsg << "* fL =    " << fL << endl;
    return (fL > r);
}

void BeamStrippingPhysics::SecondaryParticles(PartBunchBase<double, 3> *bunch, size_t &i, double &r) {

	double opyield_m = 1.0;
    size_t tempnum = bunch->getLocalNum();
	size_t count = 0;

//	*gmsg << "* random number = " << r
//		  << " NCS_single_all/NCS_total_all = " << NCS_single_all/NCS_total_all
//    	  << " NCS_double_all/NCS_total_all = " << NCS_double_all/NCS_total_all << endl;

	// change the mass_m and charge_m
	if(mass_m-m_hm < 1E-6 && charge_m == -q_e) {
		if(r > NCS_double_all/NCS_total_all) {
			*gmsg << "* Particle " << bunch->ID[i] << " is transformed to neutral hydrogen" << endl;
			bunch->M[i] = m_h;
			bunch->Q[i] = 0.0;
		}
		else {
			*gmsg << "* Particle " << bunch->ID[i] << " is transformed to proton" << endl;
			bunch->M[i] = m_p;
			bunch->Q[i] = q_e;
		}
	}
	else if(mass_m-m_p < 1E-6 && charge_m == q_e) {
		if(r > NCS_double_all/NCS_total_all) {
			*gmsg << "* Particle " << bunch->ID[i] << " is transformed to neutral hydrogen" << endl;
			bunch->M[i] = m_h;
			bunch->Q[i] = 0.0;
		}
		else {
			*gmsg << "* Particle " << bunch->ID[i] << " is transformed to negative hydrogen ion" << endl;
			bunch->M[i] = m_hm;
			bunch->Q[i] = -q_e;
		}
	}
	else if(mass_m-m_h < 1E-6 && charge_m == 0.0) {
		if(r > NCS_double_all/NCS_total_all) {
			*gmsg << "* Particle " << bunch->ID[i] << " is transformed to proton" << endl;
			bunch->M[i] = m_p;
			bunch->Q[i] = q_e;
		}
		else {
			*gmsg << "* Particle " << bunch->ID[i] << " is transformed to negative hydrogen ion" << endl;
			bunch->M[i] = m_hm;
			bunch->Q[i] = -q_e;
		}
	}
	bunch->PType[i] = ParticleType::NEWSECONDARY;

    int j = 1;
    //create new particles
    while (j < opyield_m){
      bunch->create(1);
      bunch->R[tempnum+count] = bunch->R[i];
      bunch->P[tempnum+count] = bunch->P[i];
      bunch->Q[tempnum+count] = bunch->Q[i];
      bunch->M[tempnum+count] = bunch->M[i];
      // once the particle is stripped, change PType from 0 to 1 as a flag so as to avoid repetitive stripping.
      bunch->PType[tempnum+count] = ParticleType::NEWSECONDARY;
      count++;
      j++;
    }

    if(bunch->weHaveBins())
        bunch->Bin[bunch->getLocalNum()-1] = bunch->Bin[i];
}

//void BeamStrippingPhysics::ResetMQ(PartBunchBase<double, 3> *bunch, double &r) {
//
//	if(mass_m-m_hm < 1E-6 && charge_m == -q_e) {
//		if(r > NCS_double_all/NCS_total_all) {
//			bunch->resetM(m_h * 1.0e9);
//			bunch->resetQ(0.0);
//		}
//		else {
//			bunch->resetM(m_p * 1.0e9);
//			bunch->resetQ(z_p);
//		}
//	}
//	else if(mass_m-m_p < 1E-6 && charge_m == q_e) {
//		if(r > NCS_double_all/NCS_total_all) {
//			bunch->resetM(m_h * 1.0e9);
//			bunch->resetQ(0.0);
//		}
//		else {
//			bunch->resetM(m_hm * 1.0e9);
//			bunch->resetQ(-z_p);
//		}
//	}
//	else if(mass_m-m_h < 1E-6 && charge_m == 0.0) {
//		if(r > NCS_double_all/NCS_total_all) {
//			bunch->resetM(m_p * 1.0e9);
//			bunch->resetQ(z_p);
//		}
//		else {
//			bunch->resetM(m_hm * 1.0e9);
//			bunch->resetQ(-z_p);
//		}
//	}
//}

double BeamStrippingPhysics::RandomGenerator() {
    // Random number function based on the GNU Scientific Library
    // Returns a random float between 0 and 1, exclusive: (0,1)
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    struct timeval tv; // Seed generation based on time
    gettimeofday(&tv,0);
    unsigned long mySeed = tv.tv_sec + tv.tv_usec;
    T = gsl_rng_default; // Generator setup
    r = gsl_rng_alloc (T);
    gsl_rng_set(r, mySeed);
    double u = gsl_rng_uniform(r); // Generate it!
    gsl_rng_free (r);
    return (double)u;
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
/// The material of the vacuum for beam stripping
//  ------------------------------------------------------------------------
void  BeamStrippingPhysics::Material() {

	material_m == "VACUUM";

    const double pressure = bstp_m->getPressure();				// mbar
	const double temperature = bstp_m->getTemperature();		// K

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

/*
const double BeamStrippingPhysics::fCrossSectionSingle[3][48] = {
		// Nitrogen
		{
				1.05E-15,
				1.31E-15,
				1.48E-15,
				1.59E-15,
				1.67E-15,
				1.74E-15,
				1.80E-15,
				1.85E-15,
				1.90E-15,
				2.08E-15,
				2.09E-15,
				2.05E-15,
				1.99E-15,
				1.96E-15,
				1.91E-15,
				1.86E-15,
				1.83E-15,
				1.77E-15,
				9.23E-16,
				6.87E-16,
				5.83E-16,
				5.42E-16,
				4.83E-16,
				3.50E-16,
				2.80E-16,
				2.40E-16,
				1.83E-17,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0
		},
		//Oxygen
		{
				5.990E-16,
				5.880E-16,
				6.350E-16,
				6.980E-16,
				6.970E-16,
				7.000E-16,
				7.220E-16,
				6.970E-16,
				7.430E-16,
				8.200E-16,
				8.220E-16,
				9.280E-16,
				9.580E-16,
				1.006E-15,
				1.225E-15,
				1.153E-15,
				1.320E-15,
				1.300E-15,
				1.260E-15,
				1.220E-15,
				1.200E-15,
				1.180E-15,
				1.170E-15,
				1.090E-15,
				9.750E-16,
				1.006E-15,
				9.885E-16,
				1.019E-15,
				1.033E-15,
				1.062E-15,
				1.078E-15,
				1.078E-15,
				9.890E-16,
				9.950E-16,
				1.028E-15,
				1.050E-15,
				1.040E-15,
				3.440E-16,
				2.760E-16,
				2.350E-16,
				1.910E-16,
				1.760E-16,
				1.260E-16,
				1.990E-17,
				0.0,
				0.0,
				0.0,
				0.0
		},
		//Argon
		{
				1.820E-16,
				2.170E-16,
				2.440E-16,
				2.670E-16,
				2.860E-16,
				3.100E-16,
				3.310E-16,
				3.520E-16,
				3.730E-16,
				3.940E-16,
				4.150E-16,
				4.360E-16,
				4.560E-16,
				4.760E-16,
				4.940E-16,
				5.110E-16,
				5.280E-16,
				5.430E-16,
				5.580E-16,
				5.750E-16,
				5.950E-16,
				6.190E-16,
				6.500E-16,
				7.380E-16,
				8.050E-16,
				8.820E-16,
				9.730E-16,
				1.079E-15,
				1.202E-15,
				1.343E-15,
				1.503E-15,
				1.684E-15,
				1.460E-15,
				1.330E-15,
				1.190E-15,
				9.700E-16,
				8.874E-16,
				8.150E-16,
				6.860E-16,
				5.680E-16,
				5.705E-16,
				5.005E-16,
				5.100E-16,
				4.290E-16,
				4.170E-16,
				3.630E-16,
				3.270E-16,
				5.560E-17
		}
//		//Hydrogen (molecular)
//		{
//			0.000E+00,
//			1.270E-17,
//			8.330E-17,
//			1.430E-16,
//			1.940E-16,
//			2.840E-16,
//			3.140E-16,
//			3.510E-16,
//			3.760E-16,
//			4.040E-16,
//			4.120E-16,
//			4.240E-16,
//			4.250E-16,
//			4.180E-16,
//			4.250E-16,
//			4.250E-16,
//			4.350E-16,
//			4.390E-16,
//			4.480E-16,
//			4.510E-16,
//			4.510E-16,
//			4.540E-16,
//			4.630E-16,
//			4.720E-16,
//			4.750E-16,
//			4.820E-16,
//			4.800E-16,
//			4.880E-16,
//			4.850E-16,
//			4.940E-16,
//			5.050E-16,
//			5.130E-16,
//			5.140E-16,
//			5.150E-16,
//			5.260E-16,
//			5.410E-16,
//			7.450E-16,
//			8.280E-16,
//			8.850E-16,
//			9.260E-16,
//			9.520E-16,
//			9.830E-16,
//			1.010E-15,
//			1.040E-15,
//			1.046E-15,
//			1.070E-15,
//			1.086E-15,
//			1.070E-15,
//			1.075E-15,
//			1.103E-15,
//			1.120E-15,
//			1.053E-15,
//			9.900E-16,
//			9.570E-16,
//			8.907E-16,
//			8.440E-16,
//			8.070E-16
//			7.690E-16,
//			7.677E-16,
//			7.203E-16,
//			6.600E-16,
//			5.960E-16,
//			5.577E-16,
//			4.650E-16,
//			4.230E-16,
//			3.850E-16,
//			3.310E-16,
//			2.840E-16,
//			2.430E-16,
//			2.340E-16,
//			1.690E-16,
//			1.345E-16,
//			1.170E-16,
//			8.650E-17,
//			6.820E-17,
//			6.760E-17,
//			5.840E-17,
//			5.500E-17,
//			4.810E-17,
//			4.400E-17,
//			3.860E-17,
//			3.440E-17,
//			3.000E-17,
//			2.490E-17,
//			2.070E-17,
//			1.570E-17,
//			5.770E-18,
//			3.430E-18,
//			2.070E-18,
//			1.600E-18,
//			1.290E-18,
//		}
};

const double BeamStrippingPhysics::fEnergyCSSingle[3][48] = {
		// Nitrogen
		{
				2.00E+02,
				3.00E+02,
				4.00E+02,
				5.00E+02,
				6.00E+02,
				7.00E+02,
				8.00E+02,
				9.00E+02,
				1.00E+03,
				2.00E+03,
				3.00E+03,
				4.00E+03,
				5.00E+03,
				6.00E+03,
				7.00E+03,
				8.00E+03,
				9.00E+03,
				1.00E+04,
				1.00E+05,
				2.00E+05,
				3.00E+05,
				4.00E+05,
				5.00E+05,
				9.00E+05,
				1.10E+06,
				1.30E+06,
				1.46E+07,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0
		},
		//Oxygen
		{
				7.390E+00,
				9.290E+00,
				1.440E+01,
				1.750E+01,
				2.220E+01,
				2.460E+01,
				2.950E+01,
				3.350E+01,
				3.810E+01,
				5.270E+01,
				6.900E+01,
				8.810E+01,
				1.070E+02,
				1.400E+02,
				1.950E+02,
				3.000E+02,
				4.000E+02,
				5.000E+02,
				6.000E+02,
				7.000E+02,
				8.000E+02,
				9.000E+02,
				1.000E+03,
				2.000E+03,
				3.000E+03,
				4.000E+03,
				5.000E+03,
				6.000E+03,
				7.000E+03,
				8.000E+03,
				9.000E+03,
				1.000E+04,
				1.200E+04,
				1.600E+04,
				2.000E+04,
				2.400E+04,
				3.000E+04,
				4.000E+05,
				6.000E+05,
				8.000E+05,
				1.000E+06,
				1.250E+06,
				1.500E+06,
				1.460E+07,
				0.0,
				0.0,
				0.0,
				0.0
		},
		//Argon
		{
				1.000E+00,
				1.334E+00,
				1.778E+00,
				2.370E+00,
				3.160E+00,
				4.220E+00,
				5.620E+00,
				7.500E+00,
				1.000E+01,
				1.335E+01,
				1.778E+01,
				2.370E+01,
				3.160E+01,
				4.220E+01,
				5.620E+01,
				7.500E+01,
				1.000E+02,
				1.334E+02,
				1.778E+02,
				2.370E+02,
				3.160E+02,
				4.220E+02,
				5.620E+02,
				1.000E+03,
				1.334E+03,
				1.778E+03,
				2.370E+03,
				3.160E+03,
				4.220E+03,
				5.620E+03,
				7.500E+03,
				1.000E+04,
				2.500E+04,
				3.000E+04,
				5.000E+04,
				7.500E+04,
				1.000E+05,
				1.500E+05,
				2.000E+05,
				3.000E+05,
				4.000E+05,
				5.000E+05,
				6.000E+05,
				8.000E+05,
				1.000E+06,
				1.500E+06,
				1.750E+06,
				1.460E+07
		}
//		//Hydrogen (molecular)
//		{
//			1.450E+00,
//			1.520E+00,
//			2.150E+00,
//			2.650E+00,
//			3.260E+00,
//			4.560E+00,
//			5.110E+00,
//			6.510E+00,
//			7.960E+00,
//			1.030E+01,
//			1.150E+01,
//			1.280E+01,
//			1.410E+01,
//			1.550E+01,
//			1.670E+01,
//			1.800E+01,
//			2.150E+01,
//			2.460E+01,
//			2.810E+01,
//			3.140E+01,
//			3.430E+01,
//			3.760E+01,
//			4.090E+01,
//			4.440E+01,
//			4.760E+01,
//			5.070E+01,
//			5.400E+01,
//			5.770E+01,
//			6.380E+01,
//			7.030E+01,
//			7.720E+01,
//			8.410E+01,
//			9.020E+01,
//			9.800E+01,
//			1.090E+02,
//			1.220E+02,
//			2.000E+02,
//			3.000E+02,
//			4.000E+02,
//			5.000E+02,
//			6.000E+02,
//			7.000E+02,
//			8.000E+02,
//			9.000E+02,
//			3.000E+03,
//			4.000E+03,
//			5.000E+03,
//			6.000E+03,
//			7.000E+03,
//			8.000E+03,
//			9.000E+03,
//			1.000E+04,
//			1.100E+04,
//			1.300E+04,
//			1.500E+04,
//			1.800E+04,
//			2.000E+04,
//			2.300E+04,
//			2.500E+04,
//			3.000E+04,
//			3.300E+04,
//			4.000E+04,
//			5.000E+04,
//			6.500E+04,
//			8.000E+04,
//			1.000E+05,
//			1.300E+05,
//			1.500E+05,
//			2.000E+05,
//			2.500E+05,
//			3.000E+05,
//			4.000E+05,
//			5.000E+05,
//			6.000E+05,
//			8.000E+05,
//			9.000E+05,
//			1.000E+06,
//			1.100E+06,
//			1.250E+06,
//			1.300E+06,
//			1.500E+06,
//			1.750E+06,
//			2.000E+06,
//			2.500E+06,
//			3.000E+06,
//			3.500E+06,
//			4.200E+06,
//			7.400E+06,
//			9.800E+06,
//			1.460E+07,
//			1.790E+07
//		}
};

const double BeamStrippingPhysics::fCrossSectionDouble[3][40] = {
		// Nitrogen
		{
				1.07E-16,
				1.10E-16,
				1.13E-16,
				1.43E-16,
				1.54E-16,
				1.74E-16,
				2.20E-16,
				2.60E-16,
				2.80E-16,
				3.00E-16,
				2.80E-16,
				2.00E-16,
				8.95E-17,
				6.70E-17,
				4.80E-17,
				3.45E-17,
				1.70E-17,
				1.40E-17,
				9.00E-18,
				5.45E-19,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0
		},
		//Oxygen
		{
				9.150E-17,
				9.570E-17,
				1.030E-16,
				1.090E-16,
				1.320E-16,
				1.400E-16,
				1.480E-16,
				1.640E-16,
				2.200E-16,
				1.760E-16,
				1.480E-16,
				1.160E-16,
				9.800E-17,
				9.200E-17,
				7.000E-17,
				5.600E-17,
				5.600E-17,
				7.000E-17,
				1.180E-18,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0
		},
		//Argon
		{
				8.000E-20,
				3.100E-19,
				5.400E-19,
				8.000E-19,
				1.110E-18,
				1.480E-18,
				1.940E-18,
				2.500E-18,
				3.200E-18,
				4.000E-18,
				5.100E-18,
				6.300E-18,
				7.900E-18,
				9.900E-18,
				1.230E-17,
				1.540E-17,
				1.910E-17,
				2.370E-17,
				2.950E-17,
				3.660E-17,
				4.550E-17,
				5.650E-17,
				6.813E-17,
				9.035E-17,
				1.150E-16,
				1.380E-16,
				1.490E-16,
				1.610E-16,
				1.680E-16,
				1.945E-16,
				2.370E-16,
				1.847E-16,
				1.700E-16,
				1.055E-16,
				8.975E-17,
				8.000E-17,
				1.020E-16,
				6.250E-17,
				3.480E-17,
				1.960E-18
		}
//		//Hydrogen 9molecular)
//		{
//			1.24E-17,
//			2.89E-17,
//			3.57E-17,
//			3.78E-17,
//			3.88E-17,
//			3.90E-17,
//			3.96E-17,
//			4.08E-17,
//			4.10E-17,
//			4.09E-17,
//			4.12E-17,
//			4.14E-17,
//			4.19E-17,
//			4.16E-17,
//			4.16E-17,
//			4.20E-17,
//			4.18E-17,
//			3.80E-17,
//			3.39E-17,
//			2.88E-17,
//			2.81E-17,
//			2.32E-17,
//			2.14E-17,
//			1.53E-17,
//			1.45E-17,
//			1.00E-17,
//			4.30E-18,
//			3.00E-18
//		}
};

const double BeamStrippingPhysics::fEnergyCSDouble[3][40] = {
		// Nitrogen
		{
				5.40E+03,
				1.00E+04,
				1.50E+04,
				2.50E+04,
				3.00E+04,
				4.00E+04,
				4.50E+04,
				5.00E+04,
				6.00E+04,
				7.00E+04,
				1.00E+05,
				1.50E+05,
				2.00E+05,
				3.00E+05,
				4.00E+05,
				5.00E+05,
				9.00E+05,
				1.10E+06,
				1.30E+06,
				1.46E+07,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0
		},
		//Oxygen
		{
				5.500E+03,
				1.000E+04,
				1.500E+04,
				2.000E+04,
				2.500E+04,
				3.000E+04,
				3.500E+04,
				4.000E+04,
				5.000E+04,
				6.000E+04,
				7.000E+04,
				1.100E+05,
				1.300E+05,
				1.500E+05,
				1.700E+05,
				1.900E+05,
				2.100E+05,
				2.300E+05,
				1.460E+07,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0,
				0.0
		},
		//Argon
		{
				1.778E+01,
				2.370E+01,
				3.160E+01,
				4.220E+01,
				5.620E+01,
				7.500E+01,
				1.000E+02,
				1.334E+02,
				1.778E+02,
				2.370E+02,
				3.160E+02,
				4.220E+02,
				5.620E+02,
				7.500E+02,
				1.000E+03,
				1.334E+03,
				1.778E+03,
				2.370E+03,
				3.160E+03,
				4.220E+03,
				5.620E+03,
				7.500E+03,
				1.000E+04,
				1.500E+04,
				2.000E+04,
				2.500E+04,
				3.000E+04,
				3.500E+04,
				4.000E+04,
				5.000E+04,
				7.750E+04,
				1.000E+05,
				1.100E+05,
				1.500E+05,
				2.000E+05,
				2.200E+05,
				3.000E+05,
				4.000E+05,
				5.000E+05,
				1.460E+07
		}
		//Hydrogen (molecular)
//		{
//			1.000E+03,
//			2.000E+03,
//			3.000E+03,
//			4.000E+03,
//			5.000E+03,
//			7.000E+03,
//			9.000E+03,
//			1.100E+04,
//			1.300E+04,
//			1.500E+04,
//			1.800E+04,
//			2.000E+04,
//			2.300E+04,
//			2.500E+04,
//			3.000E+04,
//			3.300E+04,
//			4.000E+04,
//			5.000E+04,
//			6.500E+04,
//			8.000E+04,
//			1.000E+05,
//			1.300E+05,
//			1.600E+05,
//			2.000E+05,
//			2.500E+05,
//			3.000E+05,
//			9.000E+05,
//			1.100E+06
//		}
};
*/
