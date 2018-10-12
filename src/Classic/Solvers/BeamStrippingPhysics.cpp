// Class:BeamStrippingPhysics
// Defines the beam stripping physics models
// ------------------------------------------------------------------------
// Class category:
// ------------------------------------------------------------------------
// $Date: 2009/07/20 09:32:31 $
// $Author: Bi, Yang Stachel, Adelmann$
//-------------------------------------------------------------------------

//#include "AbsBeamline/Degrader.h"
//#include "AbsBeamline/Drift.h"
//#include "AbsBeamline/SBend.h"
//#include "AbsBeamline/RBend.h"
//#include "AbsBeamline/Multipole.h"
#include "Algorithms/PartBunchBase.h"
#include "Physics/Physics.h"
#include "Solvers/BeamStrippingPhysics.hh"
#include "Structure/LossDataSink.h" // OPAL file
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

using Physics::pi;
using Physics::m_p;
using Physics::m_e;
using Physics::r_e;
using Physics::z_p;
using Physics::Avo;
using Physics::kB;
using Physics::q_e;

#ifdef OPAL_DKS
const int BeamStrippingPhysics::numpar = 13;
#endif


namespace {
    struct InsideTester {
        virtual ~InsideTester() {}

        virtual bool checkHit(const Vector_t &R, const Vector_t &P, double dt) = 0;
    };

    struct BeamStrippingInsideTester: public InsideTester {
        BeamStrippingInsideTester(ElementBase * el) {
            bstp_m = static_cast<BeamStripping*>(el);
        }
        virtual bool checkHit(const Vector_t &R, const Vector_t &P, double dt) {
            return bstp_m->checkPoint(R(0), R(1), R(2));
        }

    private:
        BeamStripping *bstp_m;
    };
}


BeamStrippingPhysics::BeamStrippingPhysics(const std::string &name, ElementBase *element, std::string &material):
    ParticleMatterInteractionHandler(name, element),
    // :FIXME: unused
    //allParticlesIn_m(false),
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
    locPartsInMat_m(0),
    Eavg_m(0.0),
    Emax_m(0.0),
    Emin_m(0.0)
//#ifdef OPAL_DKS
//    , curandInitSet(0)
//    , ierr(0)
//    , maxparticles(0)
//    , numparticles(0)
//    , numlocalparts(0)
//    , par_ptr(NULL)
//    , mem_ptr(NULL)
//#endif
{

    gsl_rng_env_setup();
    rGen_m = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rGen_m, Options::seed);

    Material();
    bstpshape_m = element_ref_m->getType();
    FN_m = element_ref_m->getName();

    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(FN_m, !Options::asciidump));

//#ifdef OPAL_DKS
//    if (IpplInfo::DKSEnabled) {
//        dksbase.setAPI("Cuda", 4);
//        dksbase.setDevice("-gpu", 4);
//        dksbase.initDevice();
//        curandInitSet = -1;
//    }
//
//    BeamStrippingInitTimer_m = IpplTimings::getTimer("BeamStrippingInit");
//#endif

    BeamStrippingApplyTimer_m = IpplTimings::getTimer("BeamStrippingApply");
    BeamStrippingLoopTimer_m = IpplTimings::getTimer("BeamStrippingLoop");
    BeamStrippingDestroyTimer_m = IpplTimings::getTimer("BeamStrippingDestroy");
}


BeamStrippingPhysics::~BeamStrippingPhysics() {
    locParts_m.clear();
    lossDs_m->save();
    if (rGen_m)
        gsl_rng_free(rGen_m);

#ifdef OPAL_DKS
    if (IpplInfo::DKSEnabled)
        clearCollimatorDKS();
#endif

}


void BeamStrippingPhysics::doPhysics(PartBunchBase<double, 3> *bunch) {
    /***
        Do physics if
        -- particle in material
        -- particle not dead (locParts_m[i].label != -1.0)

        Absorbed particle i: locParts_m[i].label = -1.0;

        Particle goes back to beam if
        -- not absorbed and out of material
    */

	*gmsg << "* BeamStrippingPhysics::doPhysics " << endl;

    InsideTester *tester;
    tester = new BeamStrippingInsideTester(element_ref_m);

    BeamStripping *bstp = dynamic_cast<BeamStripping *>(getElement()->removeWrappers());
	double pressure = bstp->getPressure();
	double temperature = bstp->getTemperature();

	calcNumMolecules(pressure, temperature);

    for (size_t i = 0; i < locParts_m.size(); ++i) {
        Vector_t &R = locParts_m[i].Rincol;
        Vector_t &P = locParts_m[i].Pincol;

        double Eng = (sqrt(1.0  + dot(P, P)) - 1) * m_p;
        CrossSection(Eng);
        if (locParts_m[i].label != -1) {
            if (tester->checkHit(R, P, dT_m)) {
//                bool pdead = EnergyLoss(Eng, dT_m);
//                if (!pdead) {
//                    double ptot = sqrt((m_p + Eng) * (m_p + Eng) - (m_p) * (m_p)) / m_p;
//                    P = ptot * P / sqrt(dot(P, P));
//                    /*
//                      Now scatter and transport particle in material.
//                      The checkInColl call just above will detect if the
//                      particle is rediffused from the material into vacuum.
//                    */
//                    // INFOMSG("final energy: " << (sqrt(1.0  + dot(P, P)) - 1) * m_p /1000 << " MeV" <<endl);
//                    CoulombScat(R, P, dT_m);
//                } else {
//                    // The particle is stopped in the material, set lable_m to -1
//                    locParts_m[i].label = -1.0;
//                    stoppedPartStat_m++;
//                    lossDs_m->addParticle(R, P, -locParts_m[i].IDincol);
//                }
            } else {
                /* The particle exits the material but is still in the loop of the substep,
                   Finish the timestep by letting the particle drift and after the last
                   substep call addBackToBunch
                */
                double gamma = (Eng + m_p) / m_p;
                double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
                if (bstpshape_m == ElementBase::BEAMSTRIPPING) {
                    R = R + dT_m * beta * Physics::c * P / sqrt(dot(P, P)) * 1000;
                } else {
                    R = R + dT_m * Physics::c * P / sqrt(1.0 + dot(P, P)) ;
                    addBackToBunch(bunch, i);
                }
            }
        }
    }

    delete tester;
}


double BeamStrippingPhysics::calcNumMolecules(double &pressure, double &temperature) {

	*gmsg << "* BeamStrippingPhysics::calcNumMolecules " << endl;

	if(pressure == 0.0)
        throw LogicalError("BeamStrippingPhysics::setPressure()", "Pressure must not be zero");
	if(temperature == 0.0)
		throw LogicalError("BeamStrippingPhysics::setTemperature()", "Temperature must not be zero");

    NumMolecules_m = 100 * pressure / (kB * q_e * temperature);
    *gmsg << "* NumberMolecules	= " << NumMolecules_m  << endl;
    return NumMolecules_m;
}

void BeamStrippingPhysics::CrossSection(double &Eng){

	*gmsg << "* BeamStrippingPhysics::CrossSection " << endl;

    BeamStripping *bstp = dynamic_cast<BeamStripping *>(getElement()->removeWrappers());
    vector<double> sigma = bstp->getCrossSection();
    vector<double> energycs = bstp->getEnergyCS();

	double CS;
	double m;
	double n;

	int a = energycs.size();
	*gmsg << "* Energy = " << Eng << endl;

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
			if(Eng == energycs[i]) {
				CS = sigma[i];
			}
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
	*gmsg << "* Cross section = " << CS << endl;
}


///// Energy Loss:  using the Bethe-Bloch equation.
///// Energy straggling: For relatively thick absorbers such that the number of collisions is large,
///// the energy loss distribution is shown to be Gaussian in form.
//// -------------------------------------------------------------------------
//bool BeamStrippingPhysics::EnergyLoss(double &Eng, const double &deltat) {
//    /// Eng GeV
//
//    double dEdx = 0.0;
//    const double gamma = (Eng + m_p) / m_p;
//    const double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
//    const double gamma2 = gamma * gamma;
//    const double beta2 = beta * beta;
//
//    const double deltas = deltat * beta * Physics::c;
//    const double deltasrho = deltas * 100 * rho_m;
//    const double K = 4.0 * pi * Avo * r_e * r_e * m_e * 1e7;
//    const double sigma_E = sqrt(K * m_e * rho_m * (Z_m / A_m) * deltas * 1e5);
//
//    if ((Eng > 0.00001) && (Eng < 0.0006)) {
//        const double Ts = (Eng * 1e6) / 1.0073; // 1.0073 is the proton mass divided by the atomic mass number. T is in KeV
//        const double epsilon_low = A2_c * pow(Ts, 0.45);
//        const double epsilon_high = (A3_c / Ts) * log(1 + (A4_c / Ts) + (A5_c * Ts));
//        const double epsilon = (epsilon_low * epsilon_high) / (epsilon_low + epsilon_high);
//        dEdx = -epsilon / (1e21 * (A_m / Avo)); // Stopping_power is in MeV INFOMSG("stopping power: " << dEdx << " MeV" << endl);
//        const double delta_Eave = deltasrho * dEdx;
//        const double delta_E = delta_Eave + gsl_ran_gaussian(rGen_m, sigma_E);
//        Eng = Eng + delta_E / 1e3;
//    }
//
//    if (Eng >= 0.0006) {
//        const double Tmax = (2.0 * m_e * 1e9 * beta2 * gamma2 /
//                             (1.0 + 2.0 * gamma * m_e / m_p + (m_e / m_p) * (m_e / m_p)));
//        dEdx = (-K * z_p * z_p * Z_m / (A_m * beta2) *
//                (1.0 / 2.0 * std::log(2 * m_e * 1e9 * beta2 * gamma2 * Tmax / I_m / I_m) - beta2));
//
//        // INFOMSG("stopping power_BB: " << dEdx << " MeV" << endl);
//        const double delta_Eave = deltasrho * dEdx;
//        double tmp = gsl_ran_gaussian(rGen_m, sigma_E);
//        const double delta_E = delta_Eave + tmp;
//        Eng = Eng + delta_E / 1e3;
//    }
//
//    // INFOMSG("final energy: " << Eng/1000 << " MeV" <<endl);
//    return ((Eng < 1e-4) || (dEdx > 0));
//}
//
//
void BeamStrippingPhysics::apply(PartBunchBase<double, 3> *bunch,
                              const std::pair<Vector_t, double> &boundingSphere,
                              size_t numParticlesInSimulation) {

//	IpplTimings::startTimer(BeamStrippingApplyTimer_m);

	*gmsg << "* BeamStrippingPhysics::apply " << endl;

    Inform m ("BeamStrippingPhysics::apply ", INFORM_ALL_NODES);
    /*
      Particles that have entered material are flagged as Bin[i] == -1.
      Fixme: should use PType

      Flagged particles are copied to a local structue within Beam Stripping Physics locParts_m.

      Particles in that structure will be pushed in the material and either come
      back to the bunch or will be fully stopped in the material. For the push in the
      material we use sub-timesteps.

      Newely entered particles will be copied to locParts_m at the end of apply.
    */

    Eavg_m = 0.0;
    Emax_m = 0.0;
    Emin_m = 0.0;

    bunchToMatStat_m  = 0;
    rediffusedStat_m   = 0;
    stoppedPartStat_m = 0;
    locPartsInMat_m   = 0;

    bool onlyOneLoopOverParticles = ! (allParticleInMat_m);

    dT_m = bunch->getdT();
    T_m  = bunch->getT();

    /*
      Because this is not propper set in the Component class when calling in the Constructor
    */

//#ifdef OPAL_DKS
//
//    if (collshape_m == DEGRADER && IpplInfo::DKSEnabled) {
//
//        //if firs call to apply setup needed accelerator resources
//        setupCollimatorDKS(bunch, numParticlesInSimulation);
//
//        int numaddback;
//        do {
//            IpplTimings::startTimer(DegraderLoopTimer_m);
//
//            //write particles to GPU if there are any to write
//            if (dksParts_m.size() > 0) {
//                //wrtie data from dksParts_m to the end of mem_ptr (offset = numparticles)
//                dksbase.writeDataAsync<PART_DKS>(mem_ptr, &dksParts_m[0],
//                                                 dksParts_m.size(), -1, numparticles);
//
//                //update number of particles on Device
//                numparticles += dksParts_m.size();
//
//                //free locParts_m vector
//                dksParts_m.erase(dksParts_m.begin(), dksParts_m.end());
//            }
//
//            //execute CollimatorPhysics kernels on GPU if any particles are there
//            if (numparticles > 0) {
//                dksbase.callCollimatorPhysics2(mem_ptr, par_ptr, numparticles);
//            }
//
//            //sort device particles and get number of particles comming back to bunch
//            numaddback = 0;
//            if (numparticles > 0) {
//                dksbase.callCollimatorPhysicsSort(mem_ptr, numparticles, numaddback);
//            }
//
//            //read particles from GPU if any are comming out of material
//            if (numaddback > 0) {
//
//                //resize dksParts_m to hold particles that need to go back to bunch
//                dksParts_m.resize(numaddback);
//
//                //read particles that need to be added back to bunch
//                //particles that need to be added back are at the end of Device array
//                dksbase.readData<PART_DKS>(mem_ptr, &dksParts_m[0], numaddback,
//                                           numparticles - numaddback);
//
//                //add particles back to the bunch
//                for (unsigned int i = 0; i < dksParts_m.size(); ++i) {
//                    if (dksParts_m[i].label == -2) {
//                        addBackToBunchDKS(bunch, i);
//                        rediffusedStat_m++;
//                    } else {
//                        stoppedPartStat_m++;
//                        lossDs_m->addParticle(dksParts_m[i].Rincol, dksParts_m[i].Pincol,
//                                              -locParts_m[dksParts_m[i].localID].IDincol);
//                    }
//                }
//
//                //erase particles that came from device from host array
//                dksParts_m.erase(dksParts_m.begin(), dksParts_m.end());
//
//                //update number of particles on Device
//                numparticles -= numaddback;
//            }
//
//            IpplTimings::stopTimer(DegraderLoopTimer_m);
//
//            if (onlyOneLoopOverParticles)
//	      copyFromBunchDKS(bunch, boundingSphere);
//
//            T_m += dT_m;
//
//            locPartsInMat_m = numparticles;
//            reduce(locPartsInMat_m, locPartsInMat_m, OpAddAssign());
//
//            int maxPerNode = bunch->getLocalNum();
//            reduce(maxPerNode, maxPerNode, OpMaxAssign());
//
//            //more than one loop only if all the particles are in this degrader
//            if (allParticleInMat_m) {
//                onlyOneLoopOverParticles = ( (unsigned)maxPerNode > bunch->getMinimumNumberOfParticlesPerCore() || locPartsInMat_m <= 0);
//            } else {
//                onlyOneLoopOverParticles = true;
//            }
//
//        } while (onlyOneLoopOverParticles == false);
//
//    } else {
//
//        do {
//            IpplTimings::startTimer(DegraderLoopTimer_m);
//
//            doPhysics(bunch);
//            /*
//              delete absorbed particles and particles that went to the bunch
//            */
//            deleteParticleFromLocalVector();
//            IpplTimings::stopTimer(DegraderLoopTimer_m);
//
//            /*
//              if we are not looping copy newly arrived particles
//            */
//            if (onlyOneLoopOverParticles)
//                copyFromBunch(bunch, boundingSphere);
//
//            T_m += dT_m;              // update local time
//
//            locPartsInMat_m = locParts_m.size();
//            reduce(locPartsInMat_m, locPartsInMat_m, OpAddAssign());
//
//
//            int maxPerNode = bunch->getLocalNum();
//            reduce(maxPerNode, maxPerNode, OpMaxAssign());
//            if (allParticleInMat_m) {
//                onlyOneLoopOverParticles = ( (unsigned)maxPerNode > bunch->getMinimumNumberOfParticlesPerCore() || locPartsInMat_m <= 0);
//            } else {
//                onlyOneLoopOverParticles = true;
//            }
//
//        } while (onlyOneLoopOverParticles == false);
//
//    }
//#else

    do {
//        IpplTimings::startTimer(BeamStrippingLoopTimer_m);
        doPhysics(bunch);
        /*
          delete absorbed particles and particles that went to the bunch
        */
        deleteParticleFromLocalVector();
//        IpplTimings::stopTimer(BeamStrippingLoopTimer_m);

        /*
          if we are not looping copy newly arrived particles
        */
        if (onlyOneLoopOverParticles)
            copyFromBunch(bunch, boundingSphere);

        T_m += dT_m;              // update local time

        locPartsInMat_m = locParts_m.size();
        reduce(locPartsInMat_m, locPartsInMat_m, OpAddAssign());

        int maxPerNode = bunch->getLocalNum();
        reduce(maxPerNode, maxPerNode, OpMaxAssign());
        if (allParticleInMat_m) {
            onlyOneLoopOverParticles = ( (unsigned)maxPerNode > bunch->getMinimumNumberOfParticlesPerCore() || locPartsInMat_m <= 0);
        } else {
            onlyOneLoopOverParticles = true;
        }

    } while (onlyOneLoopOverParticles == false);

//#endif

//    IpplTimings::stopTimer(BeamStrippingApplyTimer_m);
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
    // mean exitation energy from Leo
    if (Z_m < 13.0)
        I_m = 12 * Z_m + 7.0;
    else
        I_m = 9.76 * Z_m + (Z_m * 58.8 * std::pow(Z_m, -1.19));
}

//// Implement the rotation in 2 dimensions here
//// For details: see J. Beringer et al. (Particle Data Group), Phys. Rev. D 86, 010001 (2012), "Passage of particles through matter"
//void  BeamStrippingPhysics::Rot(double &px, double &pz, double &x, double &z, double xplane, double normP, double thetacou, double deltas, int coord) {
//    // Calculate the angle between the px and pz momenta to change from beam coordinate to lab coordinate
//    double Psixz = std::fmod(std::atan2(px,pz) + Physics::two_pi, Physics::two_pi);
//
//    // Apply the rotation about the random angle thetacou & change from beam
//    // coordinate system to the lab coordinate system using Psixz (2 dimensions)
//    if (coord == 1) {
//        x = x + deltas * px / normP + xplane * cos(Psixz);
//        z = z - xplane * sin(Psixz);
//    }
//    if (coord == 2) {
//        x = x + deltas * px / normP + xplane * cos(Psixz);
//        z = z - xplane * sin(Psixz) + deltas * pz / normP;
//    }
//
//    double pxz = sqrt(px * px + pz * pz);
//    px = pxz * cos(Psixz) * sin(thetacou) + pxz * sin(Psixz) * cos(thetacou);
//    pz = -pxz * sin(Psixz) * sin(thetacou) + pxz * cos(Psixz) * cos(thetacou);
//}
//
///// Coulomb Scattering: Including Multiple Coulomb Scattering and large angle Rutherford Scattering.
///// Using the distribution given in Classical Electrodynamics, by J. D. Jackson.
////--------------------------------------------------------------------------
//void  BeamStrippingPhysics::CoulombScat(Vector_t &R, Vector_t &P, const double &deltat) {
//
//    const double Eng = sqrt(dot(P, P) + 1.0) * m_p - m_p;
//    const double gamma = (Eng + m_p) / m_p;
//    const double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
//    const double normP = sqrt(dot(P, P));
//    const double deltas = deltat * beta * Physics::c;
//    const double theta0 = 13.6e6 / (beta * normP * m_p * 1e9) * z_p * sqrt(deltas / X0_m) * (1.0 + 0.038 * log(deltas / X0_m));
//
//    // x-direction: See Physical Review, "Multiple Scattering"
//    double z1 = gsl_ran_gaussian(rGen_m, 1.0);
//    double z2 = gsl_ran_gaussian(rGen_m, 1.0);
//    double thetacou = z2 * theta0;
//
//    while(fabs(thetacou) > 3.5 * theta0) {
//        z1 = gsl_ran_gaussian(rGen_m, 1.0);
//        z2 = gsl_ran_gaussian(rGen_m, 1.0);
//        thetacou = z2 * theta0;
//    }
//    double xplane = (z1 / sqrt(3.0) + z2) * deltas * theta0 / 2.0;
//    // Apply change in coordinates for multiple scattering but not for Rutherford
//    // scattering (take the deltas step only once each turn)
//    int coord = 1;
//    Rot(P(0), P(2), R(0), R(2), xplane, normP, thetacou, deltas, coord);
//
//
//    // Rutherford-scattering in x-direction
//    if (collshape_m == ElementBase::CCOLLIMATOR)
//        R = R * 1e-3;
//
//    double P2 = gsl_rng_uniform(rGen_m);
//    if (P2 < 0.0047) {
//        double P3 = gsl_rng_uniform(rGen_m);
//        double thetaru = 2.5 * sqrt(1 / P3) * sqrt(2.0) * theta0;
//        double P4 = gsl_rng_uniform(rGen_m);
//        if (P4 > 0.5)
//            thetaru = -thetaru;
//        coord = 0; // no change in coordinates but one in momenta-direction
//        Rot(P(0), P(2), R(0), R(2), xplane, normP, thetaru, deltas, coord);
//    }
//
//    // y-direction: See Physical Review, "Multiple Scattering"
//    z1 = gsl_ran_gaussian(rGen_m, 1.0);
//    z2 = gsl_ran_gaussian(rGen_m, 1.0);
//    thetacou = z2 * theta0;
//
//    while(fabs(thetacou) > 3.5 * theta0) {
//        z1 = gsl_ran_gaussian(rGen_m, 1.0);
//        z2 = gsl_ran_gaussian(rGen_m, 1.0);
//        thetacou = z2 * theta0;
//    }
//    double yplane = z1 * deltas * theta0 / sqrt(12.0) + z2 * deltas * theta0 / 2.0;
//    // Apply change in coordinates for multiple scattering but not for Rutherford
//    // scattering (take the deltas step only once each turn)
//    coord = 2;
//    Rot(P(1), P(2), R(1), R(2), yplane, normP, thetacou, deltas, coord);
//
//    // Rutherford-scattering in x-direction
//    if (collshape_m == ElementBase::CCOLLIMATOR)
//        R = R * 1e3;
//
//    P2 = gsl_rng_uniform(rGen_m);
//    if (P2 < 0.0047) {
//        double P3 = gsl_rng_uniform(rGen_m);
//        double thetaru = 2.5 * sqrt(1 / P3) * sqrt(2.0) * theta0;
//        double P4 = gsl_rng_uniform(rGen_m);
//        if (P4 > 0.5)
//            thetaru = -thetaru;
//        coord = 0; // no change in coordinates but one in momenta-direction
//        Rot(P(1), P(2), R(1), R(2), yplane, normP, thetaru, deltas, coord);
//    }
//}


void BeamStrippingPhysics::addBackToBunch(PartBunchBase<double, 3> *bunch, unsigned i) {

    bunch->createWithID(locParts_m[i].IDincol);
    /*
      Binincol is still <0, but now the particle is rediffused
      from the material and hence this is not a "lost" particle anymore
    */
    const unsigned int last = bunch->getLocalNum() - 1;
    bunch->Bin[last] = 1;
    bunch->R[last]   = locParts_m[i].Rincol;
    bunch->P[last]   = locParts_m[i].Pincol;
    bunch->Q[last]   = locParts_m[i].Qincol;
    bunch->Bf[last]  = locParts_m[i].Bfincol;
    bunch->Ef[last]  = locParts_m[i].Efincol;
    bunch->dt[last]  = locParts_m[i].DTincol;

    /*
      This particle is back to the bunch, by set
      ting the lable to -1.0
      the particle will be deleted.
    */
    locParts_m[i].label = -1.0;
    ++ rediffusedStat_m;
}


void BeamStrippingPhysics::copyFromBunch(PartBunchBase<double, 3> *bunch,
                                      const std::pair<Vector_t, double> &boundingSphere)
{
    const size_t nL = bunch->getLocalNum();
    if (nL == 0) return;

//    IpplTimings::startTimer(BeamStrippingDestroyTimer_m);
//    double zmax = boundingSphere.first(2) + boundingSphere.second;
//    double zmin = boundingSphere.first(2) - boundingSphere.second;
//    if (zmax < 0.0 || zmin > element_ref_m->getElementLength()) {
//        IpplTimings::stopTimer(BeamStrippingDestroyTimer_m);
//        return;
//    }

    InsideTester *tester;
    switch (bstpshape_m) {
    case ElementBase::BEAMSTRIPPING:
        tester = new BeamStrippingInsideTester(element_ref_m);
        break;
    default:
        throw OpalException("BeamStrippingPhysics::doPhysics","Unsupported element type");
    }

    size_t ne = 0;
    std::set<size_t> partsToDel;
    const unsigned int minNumOfParticlesPerCore = bunch->getMinimumNumberOfParticlesPerCore();
    for (size_t i = 0; i < nL; ++ i) {
        if ((bunch->Bin[i] == -1 || bunch->Bin[i] == 1) &&
            ((nL - ne) > minNumOfParticlesPerCore) &&
            tester->checkHit(bunch->R[i], bunch->P[i], dT_m))
        {
            PARTbstp x;
            x.localID      = i;
            x.DTincol      = bunch->dt[i];
            x.IDincol      = bunch->ID[i];
            x.Binincol     = bunch->Bin[i];
            x.Rincol       = bunch->R[i];
            x.Pincol       = bunch->P[i];
            x.Qincol       = bunch->Q[i];
            x.Bfincol      = bunch->Bf[i];
            x.Efincol      = bunch->Ef[i];
            x.label        = 0;            // allive in matter

            locParts_m.push_back(x);
            ne++;
            bunchToMatStat_m++;

            partsToDel.insert(i);
        }
    }

    for (auto it = partsToDel.begin(); it != partsToDel.end(); ++ it) {
        bunch->destroy(1, *it);
    }

    if (ne > 0) {
        bunch->performDestroy(true);
    }

    delete tester;
//    IpplTimings::stopTimer(BeamStrippingDestroyTimer_m);
}

void BeamStrippingPhysics::print(Inform& msg) {
    Inform::FmtFlags_t ff = msg.flags();

    // ToDo: need to move that to a statistics function
//#ifdef OPAL_DKS
//    if (strippingshape_m == DEGRADER && IpplInfo::DKSEnabled)
//        locPartsInMat_m = numparticles + dksParts_m.size();
//    else
//        locPartsInMat_m = locParts_m.size();
//#else
//    locPartsInMat_m = locParts_m.size();
//#endif
    reduce(locPartsInMat_m, locPartsInMat_m, OpAddAssign());
    reduce(bunchToMatStat_m, bunchToMatStat_m, OpAddAssign());
    reduce(rediffusedStat_m, rediffusedStat_m, OpAddAssign());
    reduce(stoppedPartStat_m, stoppedPartStat_m, OpAddAssign());

    if (locPartsInMat_m + bunchToMatStat_m + rediffusedStat_m + stoppedPartStat_m > 0) {
        OPALTimer::Timer time;
        msg << level2
            << "--- BeamStrippingPhysics - Name " << FN_m
            << " Material " << material_m << "\n"
            << "Particle Statistics @ " << time.time() << "\n"
            << std::setw(21) << "entered: " << Util::toStringWithThousandSep(bunchToMatStat_m) << "\n"
            << std::setw(21) << "rediffused: " << Util::toStringWithThousandSep(rediffusedStat_m) << "\n"
            << std::setw(21) << "stopped: " << Util::toStringWithThousandSep(stoppedPartStat_m) << "\n"
            << std::setw(21) << "total in material: " << Util::toStringWithThousandSep(locPartsInMat_m)
            << endl;
//         msg << "BeamStripping statistics: "
//             << " bunch to material " << bunchToMatStat_m << " rediffused " << rediffusedStat_m
//             << " stopped " << stoppedPartStat_m << endl;
    }
    msg.flags(ff);
}

bool BeamStrippingPhysics::stillActive() {
    return locPartsInMat_m != 0;
}

bool BeamStrippingPhysics::stillAlive(PartBunchBase<double, 3> *bunch) {
    bool degraderAlive = true;
    return degraderAlive;
}

namespace {
    bool myCompF(PARTbstp x, PARTbstp y) {
      return x.label > y.label;
    }
}

void BeamStrippingPhysics::deleteParticleFromLocalVector() {
    /*
      the particle to be deleted (label < 0) are all at the end of
      the vector.
    */
    sort(locParts_m.begin(), locParts_m.end(), myCompF);

    // find start of particles to delete
    std::vector<PARTbstp>::iterator inv = locParts_m.begin();

    for (; inv != locParts_m.end(); ++inv) {
        if ((*inv).label == -1)
            break;
    }
    locParts_m.erase(inv, locParts_m.end());
    locParts_m.resize(inv - locParts_m.begin());

    // update statistics
    if (locParts_m.size() > 0) {
        Eavg_m /= locParts_m.size();
        Emin_m /= locParts_m.size();
        Emax_m /= locParts_m.size();
    }
}

//#ifdef OPAL_DKS
//
//namespace {
//    bool myCompFDKS(PART_DKS x, PART_DKS y) {
//        return x.label > y.label;
//    }
//}
//
//void BeamStrippingPhysics::addBackToBunchDKS(PartBunchBase<double, 3> *bunch, unsigned i) {
//
//    int id = dksParts_m[i].localID;
//
//    bunch->createWithID(locParts_m[id].IDincol);
//
//    /*
//      Binincol is still <0, but now the particle is rediffused
//      from the material and hence this is not a "lost" particle anymore
//    */
//    unsigned int last = bunch->getLocalNum() - 1;
//    bunch->Bin[last] = 1;
//
//    bunch->R[last]   = dksParts_m[i].Rincol;
//    bunch->P[last]   = dksParts_m[i].Pincol;
//
//    bunch->Q[last]   = locParts_m[id].Qincol;
//    bunch->Bf[last]  = locParts_m[id].Bfincol;
//    bunch->Ef[last]  = locParts_m[id].Efincol;
//    bunch->dt[last]  = locParts_m[id].DTincol;
//
//    dksParts_m[i].label = -1.0;
//
//    ++ rediffusedStat_m;
//
//}
//
//void BeamStrippingPhysics::copyFromBunchDKS(PartBunchBase<double, 3> *bunch,
//					 const std::pair<Vector_t, double> &boundingSphere)
//{
//    const size_t nL = bunch->getLocalNum();
//    if (nL == 0) return;
//
//    IpplTimings::startTimer(DegraderDestroyTimer_m);
//    double zmax = boundingSphere.first(2) + boundingSphere.second;
//    double zmin = boundingSphere.first(2) - boundingSphere.second;
//    if (zmax < 0.0 || zmin > element_ref_m->getElementLength()) {
//        IpplTimings::stopTimer(DegraderDestroyTimer_m);
//        return;
//    }
//
//    InsideTester *tester;
//    switch (collshape_m) {
//    case ElementBase::DEGRADER:
//        tester = new DegraderInsideTester(element_ref_m);
//        break;
//    case ElementBase::CCOLLIMATOR:
//        tester = new CollimatorInsideTester(element_ref_m);
//        break;
//    case ElementBase::FLEXIBLECOLLIMATOR:
//        tester = new FlexCollimatorInsideTester(element_ref_m);
//        break;
//    default:
//        throw OpalException("CollimatorPhysics::doPhysics",
//                            "Unsupported element type");
//    }
//
//    size_t ne = 0;
//    const unsigned int minNumOfParticlesPerCore = bunch->getMinimumNumberOfParticlesPerCore();
//
//    for (unsigned int i = 0; i < nL; ++i) {
//        if ((bunch->Bin[i] == -1 || bunch->Bin[i] == 1) &&
//            ((nL - ne) > minNumOfParticlesPerCore) &&
//            tester->checkHit(bunch->R[i], bunch->P[i], dT_m);
//        {
//
//            PART x;
//            x.localID      = numlocalparts; //unique id for each particle
//            x.DTincol      = bunch->dt[i];
//            x.IDincol      = bunch->ID[i];
//            x.Binincol     = bunch->Bin[i];
//            x.Rincol       = bunch->R[i];
//            x.Pincol       = bunch->P[i];
//            x.Qincol       = bunch->Q[i];
//            x.Bfincol      = bunch->Bf[i];
//            x.Efincol      = bunch->Ef[i];
//            x.label        = 0;            // allive in matter
//
//            PART_DKS x_gpu;
//            x_gpu.label = x.label;
//            x_gpu.localID = x.localID;
//            x_gpu.Rincol = x.Rincol;
//            x_gpu.Pincol = x.Pincol;
//
//            locParts_m.push_back(x);
//            dksParts_m.push_back(x_gpu);
//
//            ne++;
//            bunchToMatStat_m++;
//            numlocalparts++;
//
//            //mark particle to be deleted from bunch as soon as it enters the material
//            bunch->destroy(1, i);
//        }
//    }
//    if (ne > 0) {
//      bunch->performDestroy(true);
//    }
//    IpplTimings::stopTimer(DegraderDestroyTimer_m);
//
//}
//
//void BeamStrippingPhysics::setupCollimatorDKS(PartBunchBase<double, 3> *bunch,
//        size_t numParticlesInSimulation)
//{
//
//    if (curandInitSet == -1) {
//
//        IpplTimings::startTimer(DegraderInitTimer_m);
//
//        //int size = bunch->getLocalNum() + 0.5 * bunch->getLocalNum();
//        //int size = bunch->getTotalNum() + 0.5 * bunch->getTotalNum();
//        int size = numParticlesInSimulation;
//
//        //allocate memory for parameters
//        par_ptr = dksbase.allocateMemory<double>(numpar, ierr);
//
//        //allocate memory for particles
//        mem_ptr = dksbase.allocateMemory<PART_DKS>((int)size, ierr);
//
//        maxparticles = (int)size;
//        numparticles = 0;
//        numlocalparts = 0;
//
//        //reserve space for locParts_m vector
//        locParts_m.reserve(size);
//
//        //init curand
//        dksbase.callInitRandoms(size, Options::seed);
//        curandInitSet = 1;
//
//        //create and transfer parameter array
//        Degrader *deg = static_cast<Degrader *>(element_ref_m);
//        double zBegin, zEnd;
//        deg->getDimensions(zBegin, zEnd);
//
//        double params[numpar] = {zBegin, zEnd, rho_m, Z_m,
//                                 A_m, A2_c, A3_c, A4_c, A5_c, X0_m, I_m, dT_m, 1e-4
//                                };
//        dksbase.writeDataAsync<double>(par_ptr, params, numpar);
//
//        IpplTimings::stopTimer(DegraderInitTimer_m);
//
//    }
//
//}
//
//void BeamStrippingPhysics::clearCollimatorDKS() {
//    if (curandInitSet == 1) {
//        dksbase.freeMemory<double>(par_ptr, numpar);
//        dksbase.freeMemory<PART_DKS>(mem_ptr, maxparticles);
//        curandInitSet = -1;
//    }
//
//}
//
//void BeamStrippingPhysics::deleteParticleFromLocalVectorDKS() {
//
//    /*
//      the particle to be deleted (label < 0) are all at the end of
//      the vector.
//    */
//    sort(dksParts_m.begin(), dksParts_m.end(), myCompFDKS);
//
//    // find start of particles to delete
//    std::vector<PART_DKS>::iterator inv = dksParts_m.begin() + stoppedPartStat_m + rediffusedStat_m;
//
//    /*
//      for (; inv != dksParts_m.end(); inv++) {
//      if ((*inv).label == -1)
//      break;
//      }
//    */
//
//    dksParts_m.erase(inv, dksParts_m.end());
//
//}
//
//#endif
