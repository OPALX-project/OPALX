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
    bunchToMatStat_m(0),
    stoppedPartStat_m(0),
    rediffusedStat_m(0),
    locPartsInMat_m(0)
{

//    gsl_rng_env_setup();
//    rGen_m = gsl_rng_alloc(gsl_rng_default);
//    gsl_rng_set(rGen_m, Options::seed);

	NbComponents = 3;
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

    mass = bunch->getM()*1E-9;

    double r = RandomGenerator();

	size_t tempnum = bunch->getLocalNum();
	for (unsigned int i = 0; i < tempnum; ++i) {
		if ( (bunch->Bin[i] != -1) && (tester->checkHit(bunch->R[i])) ) {
			double Eng = (sqrt(1.0  + dot(bunch->P[i], bunch->P[i])) - 1) * mass;

			CrossSection(Eng);

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


void BeamStrippingPhysics::GasDensity(double &pressure, double &temperature, int &iComp) {
	if(pressure == 0.0)
        throw LogicalError("BeamStrippingPhysics::setPressure()", "Pressure must not be zero");
	if(temperature == 0.0)
		throw LogicalError("BeamStrippingPhysics::setTemperature()", "Temperature must not be zero");

//    double TotalgasDensity_m = 100 * pressure / (kB * q_e * temperature);
	gasDensity[iComp] = 100 * pressure * fMolarFraction[iComp] / (kB * q_e * temperature);
//    *gmsg << "* Gas Density	= " << gasDensity[iComp]  << endl;
}

void BeamStrippingPhysics::CrossSection(double &Eng){

	Eng *= 1E9; // Eng eV
//	*gmsg << "* Energy = " << Eng << " eV" << endl;

	 int size_single = (sizeof(fCrossSectionSingle[0])/sizeof(fCrossSectionSingle[0][0]));
	 int size_double = (sizeof(fCrossSectionDouble[0])/sizeof(fCrossSectionDouble[0][0]));

	 // Single-electron detachment
	 for (int i = 0; i < NbComponents; ++i) {
		 double CS_single[i+1];

		 double max_single = 0.0;
		 int max_index_single = 0;
		 double max_double = 0.0;
		 int max_index_double = 0;

		 double m = 0.0;
		 double n = 0.0;

		 for (int j = 0; j < size_single; ++j){
			 if (fEnergyCSSingle[i][j] > max_single){
				 max_single = fEnergyCSSingle[i][j];
				 max_index_single = j;
			 }
		 }
		 if(Eng < fEnergyCSSingle[i][0]) {
			 m = (fCrossSectionSingle[i][1]-fCrossSectionSingle[i][0]) /
					 (fEnergyCSSingle[i][1]-fEnergyCSSingle[i][0]);
			 n = fCrossSectionSingle[i][0] - m * fEnergyCSSingle[i][0];
			 CS_single[i] = m * Eng + n;
		 }
		 else if(Eng > fEnergyCSSingle[i][max_index_single]) {
			 m = (fCrossSectionSingle[i][max_index_single]-fCrossSectionSingle[i][max_index_single-1]) /
					 (fEnergyCSSingle[i][max_index_single]-fEnergyCSSingle[i][max_index_single-1]);
			 n = fCrossSectionSingle[i][max_index_single] - m * fEnergyCSSingle[i][max_index_single];
			 CS_single[i] = m * Eng + n;
		 }
		 else if(Eng > fEnergyCSSingle[i][0] && Eng < fEnergyCSSingle[i][max_index_single]){
			 for (int j=0; j<(max_index_single); ++j){
				 if(Eng == fEnergyCSSingle[i][j])
					 CS_single[i] = fEnergyCSSingle[i][j];
				 else if(Eng > fEnergyCSSingle[i][j] && Eng < fEnergyCSSingle[i][j+1]) {
					 m = (fCrossSectionSingle[i][j+1]-fCrossSectionSingle[i][j]) /
							 (fEnergyCSSingle[i][j+1]-fEnergyCSSingle[i][j]);
					 n = fCrossSectionSingle[i][j] - m * fEnergyCSSingle[i][j];
					 CS_single[i] = m * Eng + n;
				 }
			 }
		 }
//		 *gmsg << "* CS_single[i] 	= " << CS_single[i]   << endl;

	 // Double-electron detachment
		 double CS_double[i+1];
		 for (int j = 0; j < size_double; ++j){
			 if (fEnergyCSDouble[i][j] > max_double){
				 max_double = fEnergyCSDouble[i][j];
				 max_index_double = j;
			 }
		 }
		 if(Eng < fEnergyCSDouble[i][0]) {
			 m = (fCrossSectionDouble[i][1]-fCrossSectionDouble[i][0]) /
					 (fEnergyCSDouble[i][1]-fEnergyCSDouble[i][0]);
			 n = fCrossSectionDouble[i][0] - m * fEnergyCSDouble[i][0];
			 CS_double[i] = m * Eng + n;
		 }
		 else if(Eng > fEnergyCSDouble[i][max_index_double]) {
			 m = (fCrossSectionDouble[i][max_index_double]-fCrossSectionDouble[i][max_index_double-1]) /
					 (fEnergyCSDouble[i][max_index_double]-fEnergyCSDouble[i][max_index_double-1]);
			 n = fCrossSectionDouble[i][max_index_double] - m * fEnergyCSDouble[i][max_index_double];
			 CS_double[i] = m * Eng + n;
		 }
		 else if(Eng > fEnergyCSDouble[i][0] && Eng < fEnergyCSDouble[i][max_index_double]){
			 for (int j=0; j<(max_index_double); ++j){
				 if(Eng == fEnergyCSDouble[i][j])
					 CS_double[i] = fEnergyCSDouble[i][j];
				 else if(Eng > fEnergyCSDouble[i][j] && Eng < fEnergyCSDouble[i][j+1]) {
					 m = (fCrossSectionDouble[i][j+1]-fCrossSectionDouble[i][j]) /
							 (fEnergyCSDouble[i][j+1]-fEnergyCSDouble[i][j]);
					 n = fCrossSectionDouble[i][j] - m * fEnergyCSDouble[i][j];
					 CS_double[i] = m * Eng + n;
				 }
			 }
		 }
//		 *gmsg << "* CS_double[i] 	= " << CS_double[i]   << endl;

		 CS[i] = CS_single[i] + CS_double[i];
//		 *gmsg << "* CS[i] = " << CS[i] << " cm2" << endl;
	 }

}

void BeamStrippingPhysics::FractionLost(double &Eng) {

	Eng *= 1E-9; // Eng GeV
    const double gamma = (Eng + mass) / mass;
    const double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    const double deltas = dT_m * beta * c;
//    *gmsg << "* deltas =    " << deltas << " m" << endl;

    double NCS = 0.0;
    for (int i = 0; i < NbComponents; ++i) {
    	CS[i] *= 1E-4; // CS en m2
    	NCS += CS[i] * gasDensity[i];
    }
    fg = 1 - exp(-(NCS * deltas));
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

    BeamStripping *bstp = dynamic_cast<BeamStripping *>(getElement()->removeWrappers());
	double pressure = bstp->getPressure();				// mbar
	double temperature = bstp->getTemperature();		// K

	for(int iComp = 0; iComp<NbComponents; ++iComp) {
//		*gmsg << "fMolarFraction = " << fMolarFraction[iComp] << endl;
		GasDensity(pressure, temperature, iComp);
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
};
