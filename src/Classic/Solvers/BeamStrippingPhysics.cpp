// Class:BeamStrippingPhysics
// Defines the beam stripping physics models
// ------------------------------------------------------------------------
// Class category:
// ------------------------------------------------------------------------
// $Date: 2009/07/20 09:32:31 $
// $Author: Bi, Yang Stachel, Adelmann$
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

using namespace std;

using Physics::kB;
using Physics::q_e;
using Physics::Avo;
using Physics::c;

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
    T_m(0.0),
    dT_m(0.0),
    material_m(material),
    Z_m(0),
    A_m(0.0),
    A2_c(0.0),
    A3_c(0.0),
    A4_c(0.0),
    A5_c(0.0),
    rho_m(0.0),
    X0_m(0.0),
    I_m(0.0),
    n_m(0.0),
    bunchToMatStat_m(0),
    stoppedPartStat_m(0),
    rediffusedStat_m(0),
    locPartsInMat_m(0)
{

//    gsl_rng_env_setup();
//    rGen_m = gsl_rng_alloc(gsl_rng_default);
//    gsl_rng_set(rGen_m, Options::seed);

    Material();
    bstpshape_m = element_ref_m->getType();
    FN_m = element_ref_m->getName();

    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(FN_m, !Options::asciidump));
}


BeamStrippingPhysics::~BeamStrippingPhysics() {
    lossDs_m->save();
//    if (rGen_m)
//        gsl_rng_free(rGen_m);
}


void BeamStrippingPhysics::doPhysics(PartBunchBase<double, 3> *bunch) {
    /*
        Do physics if
        -- particle in material
        -- particle not dead (bunch->Bin[i] != -1)

        Delete particle i: bunch->Bin[i] != -1;
    */

    InsideTester *tester;
    tester = new BeamStrippingInsideTester(element_ref_m);

    BeamStripping *bstp = dynamic_cast<BeamStripping *>(getElement()->removeWrappers());
	double pressure = bstp->getPressure();				// mbar
	double temperature = bstp->getTemperature();		// K
    vector<double> sigma = bstp->getCrossSection(); 	// cm2
    vector<double> energycs = bstp->getEnergyCS();		// eV

    calcNumMolecules(pressure, temperature);

    mass = bunch->getM()*1E-9;

    double r = RandomGenerator();

	size_t tempnum = bunch->getLocalNum();
	for (unsigned int i = 0; i < tempnum; ++i) {
		if ( (bunch->Bin[i] != -1) && (tester->checkHit(bunch->R[i])) ) {
			double Eng = (sqrt(1.0  + dot(bunch->P[i], bunch->P[i])) - 1) * mass;

			CrossSection(Eng, sigma, energycs);

			FractionLost(Eng);

			bool pdead = GasStripping(r);

			if (pdead) {
				// The particle is stripped in the residual gas, set lable_m to -1
				bunch->Bin[i] = -1;
				stoppedPartStat_m++;
				lossDs_m->addParticle(bunch->R[i], bunch->P[i], bunch->ID[i]);
				*gmsg << "* The particle " << bunch->ID[i] << " is stripped in remainder gas" << endl;
			}
		}
	}
	delete tester;
}


void BeamStrippingPhysics::calcNumMolecules(double &pressure, double &temperature) {
	if(pressure == 0.0)
        throw LogicalError("BeamStrippingPhysics::setPressure()", "Pressure must not be zero");
	if(temperature == 0.0)
		throw LogicalError("BeamStrippingPhysics::setTemperature()", "Temperature must not be zero");

    NumMolecules_m = 100 * pressure / (kB * q_e * temperature);
//    *gmsg << "* NumberMolecules	= " << NumMolecules_m  << endl;
}

void BeamStrippingPhysics::CrossSection(double &Eng, vector<double> &sigma, vector<double> &energycs){

	double m;
	double n;

	int a = energycs.size();
	Eng *= 1E9; // Eng eV
//	*gmsg << "* Energy = " << Eng << " eV" << endl;

	if(Eng < energycs[0]) {
		m = (sigma[1]-sigma[0]) / (energycs[1]-energycs[0]);
		n = sigma[0] - m * energycs[0];
		CS = m * Eng + n;
	}
	else if(Eng > energycs[a-1]) {
		m = (sigma[a-1]-sigma[a-2]) / (energycs[a-1]-energycs[a-2]);
		n = sigma[a-1] - m * energycs[a-1];
		CS = m * Eng + n;
	}
	else if(Eng > energycs[0] && Eng < energycs[a-1]){
		for (int i=0; i<a; i++){
			if(Eng == energycs[i])
				CS = sigma[i];
			else if(Eng > energycs[i] && Eng < energycs[i+1]) {
				m = (sigma[i+1]-sigma[i]) / (energycs[i+1]-energycs[i]);
				n = sigma[i] - m * energycs[i];
				CS = m * Eng + n;
			}
		}
	}
	else {
		*gmsg << "**Cross section not calculated**" << endl;
	}
//	*gmsg << "* Cross section = " << CS << " cm2" << endl;
}

void BeamStrippingPhysics::FractionLost(double &Eng) {

	Eng *= 1E-9; // Eng GeV
    const double gamma = (Eng + mass) / mass;
    const double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    const double deltas = dT_m * beta * c;
//    *gmsg << "* deltas =    " << deltas << " m" << endl;

    CS *= 1E-4; // CS en m2
    fg = 1 - exp(-(CS * NumMolecules_m * deltas));
//    *gmsg << "* fg =    " << fg << endl;
}


bool BeamStrippingPhysics::GasStripping(double &r) {

//	double r = gsl_rng_uniform (rGen_m);
//    *gmsg << "random uniform number = " << r << endl;

    return (fg > r);
}

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

void BeamStrippingPhysics::apply(PartBunchBase<double, 3> *bunch,
                              const std::pair<Vector_t, double> &boundingSphere,
                              size_t numParticlesInSimulation) {

    Inform m ("BeamStrippingPhysics::apply ", INFORM_ALL_NODES);

    bunchToMatStat_m  = 0;
    rediffusedStat_m   = 0;
    stoppedPartStat_m = 0;
    locPartsInMat_m   = 0;

    bool onlyOneLoopOverParticles = ! (allParticleInMat_m);

    dT_m = bunch->getdT();
    T_m  = bunch->getT();

    do {
        doPhysics(bunch);
    } while (onlyOneLoopOverParticles == false);
}


const std::string BeamStrippingPhysics::getType() const {
    return "BeamStrippingPhysics";
}

/// The material of the vacuum for beam stripping
//  ------------------------------------------------------------------------
void  BeamStrippingPhysics::Material() {

    if (material_m == "VACUUM") {
        Z_m = 4.0;
        A_m = 9.012;
        rho_m = 1.848;

        X0_m = 65.19 / rho_m / 100;
        I_m = 10 * Z_m;
        n_m = rho_m / A_m * Avo;

        A2_c = 2.590;
        A3_c = 9.660e2;
        A4_c = 1.538e2;
        A5_c = 3.475e-2;
    }

    else if (material_m == "AIR") {
    	Z_m = 7;
    	A_m = 14;
        rho_m = 0.0012;

        X0_m = 37.99 / rho_m / 100;
        I_m = 10 * Z_m;
        n_m = rho_m / A_m * Avo;

        A2_c = 3.350;
        A3_c = 1.683e3;
        A4_c = 1.900e3;
        A5_c = 2.513e-2;
    }

    else {
        throw GeneralClassicException("BeamStrippingPhysics::Material", "Material not found ...");
    }
//    // mean exitation energy from Leo
//    if (Z_m < 13.0)
//        I_m = 12 * Z_m + 7.0;
//    else
//        I_m = 9.76 * Z_m + (Z_m * 58.8 * std::pow(Z_m, -1.19));
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
