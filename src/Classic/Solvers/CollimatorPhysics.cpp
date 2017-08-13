// Class:CollimatorPhysics
// Defines the collimator physics models
// ------------------------------------------------------------------------
// Class category:
// ------------------------------------------------------------------------
// $Date: 2009/07/20 09:32:31 $
// $Author: Bi, Yang Stachel, Adelmann$
//-------------------------------------------------------------------------
#include "Solvers/CollimatorPhysics.hh"
#include "Physics/Physics.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/Collimator.h"
#include "AbsBeamline/Degrader.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/Multipole.h"
#include "Structure/LossDataSink.h" // OPAL file
#include "Utilities/Options.h"
#include "Utilities/GeneralClassicException.h"

#include "Ippl.h"

#include <iostream>
#include <fstream>
#include <algorithm>

using Physics::pi;
using Physics::m_p;
using Physics::m_e;
using Physics::h_bar;
using Physics::epsilon_0;
using Physics::r_e;
using Physics::z_p;
using Physics::Avo;

#ifdef OPAL_DKS
const int CollimatorPhysics::numpar = 13;
#endif

Vector_t ArbitraryRotation(Vector_t &W, Vector_t &Rorg, double Theta);

CollimatorPhysics::CollimatorPhysics(const std::string &name,
                                     ElementBase *element,
				     std::string &material,
                                     bool enableRutherfordScattering,
				     double lowEnergyThr):
    SurfacePhysicsHandler(name, element),
    T_m(0.0),
    dT_m(0.0),
    rGen_m(NULL),
    material_m(material),
    FN_m(""),
    collShapeStr_m(""),
    collShape_m(NOSHAPE),
    Z_m(0),
    A_m(0.0),
    A2_c(0.0),
    A3_c(0.0),
    A4_c(0.0),
    A5_c(0.0),
    rho_m(0.0),
    X0_m(0.0),
    I_m(0.0),
    enableRutherfordScattering_m(enableRutherfordScattering),
    lowEnergyThr_m(lowEnergyThr),
    bunchToMatStat_m(0),
    stoppedPartStat_m(0),
    redifusedStat_m(0),
    localPartsInMat_m(0),
    globalPartsInMat_m(0),
    Eavg_m(0.0),
    Emax_m(0.0),
    Emin_m(0.0),
    nextLocalID_m(0)
#ifdef OPAL_DKS
    , curandInitSet(0)
    , ierr(0)
    , maxparticles(0)
    , numparticles(0)
    , par_ptr(NULL)
    , mem_ptr(NULL)
#endif
{

    gsl_rng_env_setup();
    rGen_m = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(rGen_m, Options::seed);

    Material();

    if(dynamic_cast<Collimator *>(element_ref_m)) {
        Collimator *coll = dynamic_cast<Collimator *>(element_ref_m);
        FN_m = coll->getName();
        collShapeStr_m = coll->getCollimatorShape();
        collShape_m = (collShapeStr_m == "CCOLLIMATOR"? CYCLCOLLIMATOR: COLLIMATOR);
    } else if(dynamic_cast<Drift *>(element_ref_m)) {
        Drift *drf = dynamic_cast<Drift *>(element_ref_m);
        FN_m = drf->getName();
    } else if(dynamic_cast<SBend *>(element_ref_m)) {
        ERRORMSG("SBend Begin_m and End_m not defined");
    } else if(dynamic_cast<RBend *>(element_ref_m)) {
        ERRORMSG("RBend Begin_m and End_m not defined");
    } else if(dynamic_cast<Multipole *>(element_ref_m)) {
        Multipole *quad = dynamic_cast<Multipole *>(element_ref_m);
        FN_m = quad->getName();
    } else if(dynamic_cast<Degrader *>(element_ref_m)) {
        Degrader *deg = dynamic_cast<Degrader *>(element_ref_m);
        FN_m = deg->getName();
        collShapeStr_m = "DEGRADER";
        collShape_m = DEGRADER;
    }

    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(FN_m, !Options::asciidump));

#ifdef OPAL_DKS
    if (IpplInfo::DKSEnabled) {
        dksbase_m.setAPI("Cuda", 4);
        dksbase_m.setDevice("-gpu", 4);
        dksbase_m.initDevice();
        curandInitSet_m = -1;
    }

    DegraderInitTimer_m = IpplTimings::getTimer("DegraderInit");
#endif

    DegraderApplyTimer_m = IpplTimings::getTimer("DegraderApply");
    DegraderLoopTimer_m = IpplTimings::getTimer("DegraderLoop");
}

CollimatorPhysics::~CollimatorPhysics() {
    locParts_m.clear();
    lossDs_m->save();
    if (rGen_m)
        gsl_rng_free(rGen_m);

#ifdef OPAL_DKS
    if (IpplInfo::DKSEnabled)
        clearCollimatorDKS();
#endif

}

bool CollimatorPhysics::checkHit(Vector_t R, Vector_t P, double dt, Degrader *deg, Collimator *coll) {
    if(collShape_m == CYCLCOLLIMATOR) {
        return coll->checkPoint(R(0), R(1));
    } else if (collShape_m == DEGRADER) {
        return deg->isInMaterial(R);
    } else {
        return coll->isInColl(R, P, Physics::c * dt/sqrt(1.0  + dot(P, P)));
    }
}

void CollimatorPhysics::apply(PartBunch &bunch, size_t numParticlesInSimulation) {
    IpplTimings::startTimer(DegraderApplyTimer_m);

    Inform m ("CollimatorPhysics::apply ", INFORM_ALL_NODES);
    /*
      Particles that have entered material are flagged as Bin[i] == INMATERIAL.
      Fixme: should use PType

      Flagged particles are copied to a local structue within Collimator Physics locParts_m.

      Particles in that structure will be pushed in the material and either come
      back to the bunch or will be fully stopped in the material. For the push in the
      material we use sub-timesteps.

      Newely entered particles will be copied to locParts_m at the end of apply.
    */

    Eavg_m = 0.0;
    Emax_m = 0.0;
    Emin_m = 0.0;

    bunchToMatStat_m   = 0;
    redifusedStat_m    = 0;
    stoppedPartStat_m  = 0;
    localPartsInMat_m  = 0;
    globalPartsInMat_m = 0;

    dT_m = bunch.getdT();
    T_m  = bunch.getT();

#ifdef OPAL_DKS

    if (collShape_m == DEGRADER && IpplInfo::DKSEnabled) {

        applyDKS(bunch, numParticlesInSimulation);

    } else
#endif
    {
        bool onlyOneLoopOverParticles = ! (allParticleInMat_m);

        Degrader   *deg  = NULL;
        Collimator *coll = NULL;

        if(collShape_m == DEGRADER) {
            deg = dynamic_cast<Degrader *>(element_ref_m);
        }
        else {
            coll = dynamic_cast<Collimator *>(element_ref_m);
        }


        do{
            IpplTimings::startTimer(DegraderLoopTimer_m);

            doPhysics(bunch, deg, coll);
            /*
              delete absorbed particles and particles that went to the bunch
            */
            deleteParticleFromLocalVector();
            IpplTimings::stopTimer(DegraderLoopTimer_m);

            /*
              if we are not looping copy newly arrived particles
            */
            if (onlyOneLoopOverParticles)
                copyFromBunch(bunch);

            T_m += dT_m;              // update local time

            localPartsInMat_m = locParts_m.size();
            reduce(localPartsInMat_m, globalPartsInMat_m, OpAddAssign());

            int maxPerNode = bunch.getLocalNum();
            reduce(maxPerNode, maxPerNode, OpMaxAssign());
            if (allParticleInMat_m) {
                onlyOneLoopOverParticles = ( (unsigned)maxPerNode > bunch.getMinimumNumberOfParticlesPerCore() ||
                                             globalPartsInMat_m == 0);
            } else {
                onlyOneLoopOverParticles = true;
            }

        } while (!onlyOneLoopOverParticles);
    }

    IpplTimings::stopTimer(DegraderApplyTimer_m);
}

void CollimatorPhysics::doPhysics(PartBunch &bunch, Degrader *deg, Collimator *col) {
    /***
        Do physics if
        -- particle in material
        -- particle not dead (locParts_m[i].label != -1.0)

        Absorbed particle i: locParts_m[i].label = -1.0;

        Particle goes back to beam if
        -- not absorbed and out of material
    */
    static const double m2mm = 1000.0;
    const double mass = bunch.getM() * 1e-9;
    for(size_t i = 0; i < locParts_m.size(); ++ i) {
        Vector_t &R = locParts_m[i].Rincol;
        Vector_t &P = locParts_m[i].Pincol;
        double Eng = (sqrt(1.0  + dot(P, P)) - 1) * mass;

        if (locParts_m[i].label != -1) {
            if (checkHit(R, P, dT_m, deg, col)) {
                bool pdead = computeEnergyLoss(Eng, dT_m);
                if(!pdead) {
                    double ptot = sqrt((mass + Eng) * (mass + Eng) - mass * mass) / mass;
                    P = P * ptot / sqrt(dot(P, P));
                    /*
                      Now scatter and transport particle in material.
                      The checkInColl call just above will detect if the
                      particle is rediffused from the material into vacuum.
                    */
                    applyCoulombScat(R, P, dT_m);
                } else {
                    // The particle is stopped in the material, set lable_m to -1
                    locParts_m[i].label = -1.0;
                    ++ stoppedPartStat_m;
                    lossDs_m->addParticle(R, P, -locParts_m[i].IDincol);
                }
            } else {
                /* The particle exits the material but is still in the loop of the substep,
                   Finish the timestep by letting the particle drift and after the last
                   substep call addBackToBunch
                */
                if(collShape_m == CYCLCOLLIMATOR) {
                    double gamma = (Eng + mass) / mass;
                    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
                    R = R + dT_m * beta * Physics::c * P / sqrt(dot(P, P)) * m2mm;
                } else {
                    R = R + dT_m * Physics::c * P / sqrt(1.0 + dot(P, P)) ;
                }

                addBackToBunch(bunch, i);

                ++ redifusedStat_m;
            }
        }
    }
}

/// Energy Loss:  using the Bethe-Bloch equation.
/// Energy straggling: For relatively thick absorbers such that the number of collisions is large,
/// the energy loss distribution is shown to be Gaussian in form.
// -------------------------------------------------------------------------
bool CollimatorPhysics::computeEnergyLoss(double &Eng /* in GeV */, double &deltat) {

    double dEdx = 0.0;
    const double gamma = (Eng + m_p) / m_p;
    const double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    const double gamma2 = gamma * gamma;
    const double beta2 = beta * beta;

    const double deltas = deltat * beta * Physics::c;
    const double deltasrho = deltas * 100 * rho_m;
    static const double K = 4.0 * pi * Avo * r_e * r_e * m_e * 1e7;
    const double sigma_E = sqrt(K * m_e * rho_m * (Z_m / A_m) *  deltas * 1e5);

    if ((Eng > 0.00001) && (Eng < 0.0006)) {
        const double Ts = (Eng * 1e6) / 1.0073;  //  1.0073 is the mass of the proton in dalton. T is in keV
        const double epsilon_low = A2_c * pow(Ts, 0.45);
        const double epsilon_high = (A3_c / Ts) * log(1 + (A4_c / Ts) + (A5_c * Ts));
        const double epsilon = (epsilon_low * epsilon_high) / (epsilon_low + epsilon_high);
        dEdx = -epsilon  / (1e21 * A_m / Avo);

        const double delta_E = deltasrho * dEdx + gsl_ran_gaussian(rGen_m, sigma_E);
        Eng = Eng + delta_E / 1e3;
    }

    if (Eng >= 0.0006) {
        const double Tmax = (2.0 * m_e * 1e9 * beta2 * gamma2  /
                             (1.0 + 2.0 * gamma * m_e / m_p + (m_e / m_p) * (m_e / m_p)));
        dEdx = (-K * z_p * z_p * Z_m / (A_m * beta2)  *
                (0.5 * std::log(2 * m_e * 1e9 * beta2 * gamma2 * Tmax / I_m / I_m) - beta2));

        const double delta_E = deltasrho * dEdx + gsl_ran_gaussian(rGen_m, sigma_E);;
        Eng = Eng + delta_E / 1e3;
    }
    return ((Eng < lowEnergyThr_m) || (dEdx > 0));
}

// splitting the scattering into x and y plane as it is done now seems odd. Trying to
// rotate the momentum about e_z by a random angle, then scatter and finally rotate back
// but somehow this doesn't seem to work well yet; uncomment TRYNEWWAY to try this out.
//
// #define TRYNEWWAY
//
/// Coulomb Scattering: Including Multiple Coulomb Scattering and large angle Rutherford Scattering.
/// Using the distribution given in Classical Electrodynamics, by J. D. Jackson.
/// For details: see J. Beringer et al. (Particle Data Group),
///                  Phys. Rev. D 86, 010001 (2012),
///                  "Passage of particles through matter"
//--------------------------------------------------------------------------
void  CollimatorPhysics::applyCoulombScat(Vector_t &R, Vector_t &P, double &deltat) {
    static const double mm2m = 1e-3, m2mm = 1000.0;
    const double Eng = sqrt(dot(P, P) + 1.0) * m_p - m_p;
    const double gamma = (Eng + m_p) / m_p;
    const double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    const double betaGamma = sqrt(dot(P, P));
    const double deltas = deltat * beta * Physics::c;
    const double mass = m_p * 1e9; // in eV
    // (30.15) in [Beringer]
    const double theta0 = 13.6e6 / (beta * betaGamma * mass) * z_p * sqrt(deltas / X0_m) * (1.0 + 0.038 * log(deltas / X0_m));

    // Rutherford-scattering in x-direction
    if(collShape_m == CYCLCOLLIMATOR)
        R = R * mm2m;

    // x-direction: See Physical Review, "Multiple Scattering"
    double z1 = gsl_ran_gaussian(rGen_m,1.0);
    double z2 = gsl_ran_gaussian(rGen_m,1.0);
    double thetacou = z2 * theta0;

    while (std::abs(thetacou) > 3.5 * theta0) {
        z1 = gsl_ran_gaussian(rGen_m,1.0);
        z2 = gsl_ran_gaussian(rGen_m,1.0);
        thetacou = z2 * theta0;
    }

#ifdef TRYNEWWAY
    double xplane = 0.5 * deltas * theta0 * (z1 / sqrt(3.0) + z2);
    double phi = Physics::pi * gsl_rng_uniform(rGen_m);
    Vector_t X(0.0);
    Quaternion M(cos(phi), sin(phi) * Vector_t(0, 0, 1)); // rotation about (0, 0, 1) by 2 * phi
    Vector_t origP = P;
    P = M.rotate(P);
    computeScatteringEffect(P, X, xplane, betaGamma, thetacou, deltas, 0);
    P = M.conjugate().rotate(P);
    R += M.conjugate().rotate(X);

    // z1 = gsl_ran_gaussian(rGen_m,1.0);
    // z2 = gsl_ran_gaussian(rGen_m,1.0);
    // thetacou = z2 * theta0;

    // while(std::abs(thetacou) > 3.5 * theta0) {
    //     z1 = gsl_ran_gaussian(rGen_m,1.0);
    //     z2 = gsl_ran_gaussian(rGen_m,1.0);
    //     thetacou = z2 * theta0;
    // }

#else
    double xplane = 0.5 * deltas * theta0 * (z1 / sqrt(3.0) + z2);
    computeScatteringEffect(P, R, xplane, betaGamma, thetacou, deltas, 0);

    // y-direction: See Physical Review, "Multiple Scattering"
    z1 = gsl_ran_gaussian(rGen_m,1.0);
    z2 = gsl_ran_gaussian(rGen_m,1.0);
    thetacou = z2 * theta0;

    while(std::abs(thetacou) > 3.5 * theta0) {
        z1 = gsl_ran_gaussian(rGen_m,1.0);
        z2 = gsl_ran_gaussian(rGen_m,1.0);
        thetacou = z2 * theta0;
    }

    double yplane = 0.5 * deltas * theta0 * (z1 / sqrt(3.0) + z2);
    computeScatteringEffect(P, R, yplane, betaGamma, thetacou, deltas, 1);
#endif

    // Rutherford-scattering in x-direction
    if(collShape_m == CYCLCOLLIMATOR)
        R = R * m2mm;

    // Rutherford-scattering
    double P2 = gsl_rng_uniform(rGen_m);
    if ((P2 < 0.0047) && enableRutherfordScattering_m) {
        double P3 = gsl_rng_uniform(rGen_m);
        //double thetaru = 2.5 * sqrt(1 / P3) * sqrt(2.0) * theta0;
        double thetaru = 2.5 * sqrt(1 / P3) * 2.0 * theta0;
        double phiru = 2.0 * pi * gsl_rng_uniform(rGen_m);
        double th0 = atan2(sqrt(P(0)*P(0) + P(1)*P(1)), std::abs(P(2)));
        Vector_t W,X;
        X(0) = cos(phiru) * sin(thetaru);
        X(1) = sin(phiru) * sin(thetaru);
        X(2) = cos(thetaru);
        X *= sqrt(dot(P,P));
        W(0) = -P(1);
        W(1) = P(0);
        W(2) = 0.0;
        P = ArbitraryRotation(W, X, th0);
    }
}

// Implement the rotation in 2 dimensions here
void  CollimatorPhysics::computeScatteringEffect(Vector_t &P,
                                                 Vector_t &R,
                                                 double plane,
                                                 double betaGamma,
                                                 double theta,
                                                 double deltas,
                                                 unsigned char coord) {
#ifdef TRYNEWWAY
    // Calculate the angle between the px and pz momenta to change from beam coordinate to lab coordinate
    const double Psi = atan2(P(0), P(2));
    const double pxz = sqrt(P(0)*P(0) + P(2)*P(2));
    const double cosPsi = cos(Psi);
    const double sinPsi = sin(Psi);
    const double cosTheta = cos(theta);
    const double sinTheta = sin(theta);

    // Apply the rotation about the random angle thetacou & change from beam
    // coordinate system to the lab coordinate system using Psixz (2 dimensions)
    R(0) += deltas * P(0) / betaGamma + plane * cosPsi;
    R(2) += deltas * P(2) / betaGamma - plane * sinPsi;

    P(0) = pxz * (cosPsi * sinTheta + sinPsi * cosTheta);
    P(2) = pxz * (-sinPsi * sinTheta + cosPsi * cosTheta);
#else
    coord = coord % 3;
    double &px = P(coord);
    double &pz = P(2);
    double &x = R(coord);
    double &z = R(2);
    // Calculate the angle between the px and pz momenta to change from beam coordinate to lab coordinate
    const double Psi = atan2(px, pz);
    const double pxz = sqrt(px*px + pz*pz);
    const double cosPsi = cos(Psi);
    const double sinPsi = sin(Psi);
    const double cosTheta = cos(theta);
    const double sinTheta = sin(theta);

    // Apply the rotation about the random angle thetacou & change from beam
    // coordinate system to the lab coordinate system using Psixz (2 dimensions)
    x += deltas * px / betaGamma + plane * cosPsi;
    z -= plane * sinPsi;

    if (coord == 1) {
        z += deltas * pz / betaGamma;
    }

    px = pxz * (cosPsi * sinTheta + sinPsi * cosTheta);
    pz = pxz * (-sinPsi * sinTheta + cosPsi * cosTheta);
#endif
}

/// The material of the collimator
//  ------------------------------------------------------------------------
void  CollimatorPhysics::Material() {
    // mean exitation energy (I_m) from Leo

    if(material_m == "Berilium") {
        Z_m = 4.0;
        A_m = 9.012;
        rho_m = 1.848;

        X0_m = 65.19 / rho_m / 100;
        I_m = 12 * Z_m + 7.0;

        A2_c = 2.590;
        A3_c = 9.660e2;
        A4_c = 1.538e2;
        A5_c = 3.475e-2;

    }

    else if(material_m == "Graphite") {
        Z_m = 6.0;
        A_m = 12.0107;
        rho_m = 2.210;

        X0_m = 42.70 / rho_m / 100;
        I_m = 12 * Z_m + 7.0;

        A2_c = 2.601;
        A3_c = 1.701e3;
        A4_c = 1.279e3;
        A5_c = 1.638e-2;

    }

    else if(material_m == "GraphiteR6710") {
        Z_m = 6.0;
        A_m = 12.0107;
        rho_m = 1.88;

        X0_m = 42.70 / rho_m / 100;
        I_m = 12 * Z_m + 7.0;

        A2_c = 2.601;
        A3_c = 1.701e3;
        A4_c = 1.279e3;
        A5_c = 1.638e-2;

    }

    else if(material_m == "Molybdenum") {
        Z_m = 42.0;
        A_m = 95.94;
        rho_m = 10.22;

        X0_m = 9.8 / rho_m / 100;
        I_m = 9.76 * Z_m + (Z_m * 58.8 * std::pow(Z_m, -1.19));

        A2_c = 7.248;
        A3_c = 9.545e3;
        A4_c = 4.802e2;
        A5_c = 5.376e-3;
    }

    /*
      needs to be checked !

      Z from http://journals.aps.org/prb/pdf/10.1103/PhysRevB.40.8530
    */

    else if(material_m == "Mylar") {
        Z_m = 6.702;
        A_m = 12.88;
        rho_m = 1.4;

        X0_m = 39.95 / rho_m / 100;
        I_m = 12 * Z_m + 7.0;

	A2_c = 3.350;
	A3_c = 1683;
        A4_c = 1900;
        A5_c = 2.513e-02;
    }

    else if (material_m == "Aluminum") {
        Z_m = 13;
        A_m = 26.98;
        rho_m = 2.7;

        X0_m = 24.01 / rho_m / 100;
        I_m = 9.76 * Z_m + (Z_m * 58.8 * std::pow(Z_m, -1.19));

        A2_c = 4.739;
        A3_c = 2.766e3;
        A4_c = 1.645e2;
        A5_c = 2.023e-2;
    }

    else if (material_m == "Copper") {
        Z_m = 29;
        A_m = 63.54;
        rho_m = 8.96;

        X0_m = 12.86 / rho_m / 100;
        I_m = 9.76 * Z_m + (Z_m * 58.8 * std::pow(Z_m, -1.19));

        A2_c = 4.194;
        A3_c = 4.649e3;
        A4_c = 8.113e1;
        A5_c = 2.242e-2;
    }

    else if (material_m == "Titan") {
        Z_m = 22;
        A_m = 47.8;
        rho_m = 4.54;

        X0_m = 16.16 / rho_m / 100;
        I_m = 9.76 * Z_m + (Z_m * 58.8 * std::pow(Z_m, -1.19));

        A2_c = 5.489;
        A3_c = 5.260e3;
        A4_c = 6.511e2;
        A5_c = 8.930e-3;
    }

    else if (material_m == "AluminaAl2O3") {
        Z_m = 50;
        A_m = 101.96;
        rho_m = 3.97;

        X0_m = 27.94 / rho_m / 100;
        I_m = 9.76 * Z_m + (Z_m * 58.8 * std::pow(Z_m, -1.19));

        A2_c = 7.227;
        A3_c = 1.121e4;
        A4_c = 3.864e2;
        A5_c = 4.474e-3;
    }

    else if (material_m == "Air") {
        Z_m = 7;
        A_m = 14;
        rho_m = 0.0012;

        X0_m = 37.99 / rho_m / 100;
        I_m = 12 * Z_m + 7.0;

        A2_c = 3.350;
        A3_c = 1.683e3;
        A4_c = 1.900e3;
        A5_c = 2.513e-2;
    }

    else if (material_m == "Kapton") {
        Z_m = 6;
        A_m = 12;
        rho_m = 1.4;

        X0_m = 39.95 / rho_m / 100;
        I_m = 12 * Z_m + 7.0;

        A2_c = 2.601;
        A3_c = 1.701e3;
        A4_c = 1.279e3;
        A5_c = 1.638e-2;
    }

    else if (material_m == "Gold") {
        Z_m = 79;
        A_m = 197;
        rho_m = 19.3;

        X0_m = 6.46 / rho_m / 100;
        I_m = 9.76 * Z_m + (Z_m * 58.8 * std::pow(Z_m, -1.19));

        A2_c = 5.458;
        A3_c = 7.852e3;
        A4_c = 9.758e2;
        A5_c = 2.077e-2;
    }

    else if (material_m == "Water") {
        Z_m = 10;
        A_m = 18;
        rho_m = 1;

        X0_m = 36.08 / rho_m / 100;
        I_m = 12 * Z_m + 7.0;

        A2_c = 2.199;
        A3_c = 2.393e3;
        A4_c = 2.699e3;
        A5_c = 1.568e-2;
    }

    else {
        throw GeneralClassicException("CollimatorPhysics::Material",
                                      "Material '" + material_m + "' not found.");
    }
}

void CollimatorPhysics::addBackToBunch(PartBunch &bunch, unsigned i) {

    bunch.createWithID(locParts_m[i].IDincol);

    /*
      Binincol is still <0, but now the particle is rediffused
      from the material and hence this is not a "lost" particle anymore
    */
    const size_t idx = bunch.getLocalNum() - 1;
    bunch.Bin[idx] = 1;

    bunch.R[idx]           = locParts_m[i].Rincol;
    bunch.P[idx]           = locParts_m[i].Pincol;
    bunch.Q[idx]           = locParts_m[i].Qincol;
    bunch.LastSection[idx] = locParts_m[i].LastSecincol;
    bunch.Bf[idx]          = locParts_m[i].Bfincol;
    bunch.Ef[idx]          = locParts_m[i].Efincol;
    bunch.dt[idx]          = locParts_m[i].DTincol;

    /*
      This particle is back to the bunch, by set
      ting the lable to -1.0
      the particle will be deleted.
    */
    locParts_m[i].label = -1.0;
}

void CollimatorPhysics::copyFromBunch(PartBunch &bunch)
{
    Degrader   *deg  = NULL;
    Collimator *coll = NULL;

    if(collShape_m == DEGRADER)
        deg = dynamic_cast<Degrader *>(element_ref_m);
    else
        coll = dynamic_cast<Collimator *>(element_ref_m);

    const size_t nL = bunch.getLocalNum();
    size_t ne = 0;
    std::list<size_t> partsToDel;
    const unsigned int minNumOfParticlesPerCore = bunch.getMinimumNumberOfParticlesPerCore();
    for(unsigned int i = 0; i < nL; ++ i) {
        if ((bunch.Bin[i] == -1 ||
             bunch.Bin[i] == 1) &&
            ((nL-ne) > minNumOfParticlesPerCore) &&
	    checkHit(bunch.R[i], bunch.P[i], dT_m, deg, coll)) {

            PART x;
            x.localID      = nextLocalID_m;
            x.DTincol      = bunch.dt[i];
            x.IDincol      = bunch.ID[i];
            x.Binincol     = bunch.Bin[i];
            x.Rincol       = bunch.R[i];
            x.Pincol       = bunch.P[i];
            x.Qincol       = bunch.Q[i];
            x.LastSecincol = bunch.LastSection[i];
            x.Bfincol      = bunch.Bf[i];
            x.Efincol      = bunch.Ef[i];
            x.label        = 0;            // alive in matter

            locParts_m.push_back(x);
            ++ ne;
            ++ bunchToMatStat_m;
            ++ nextLocalID_m;

            //mark particle to be deleted from bunch as soon as it enters the material
            partsToDel.push_front(i);
        }
    }

    for (auto it = partsToDel.begin(); it != partsToDel.end(); ++ it) {
        bunch.destroy(1, *it);
    }
}


void CollimatorPhysics::print(Inform &msg){
    Inform::FmtFlags_t ff = msg.flags();

    // ToDo: need to move that to a statistics function
#ifdef OPAL_DKS
    if (collShape_m == DEGRADER && IpplInfo::DKSEnabled)
        localPartsInMat_m = numparticles_m + dksParts_m.size();
    else
        localPartsInMat_m = locParts_m.size();
#else
    localPartsInMat_m = locParts_m.size();
#endif

    reduce(localPartsInMat_m, globalPartsInMat_m, OpAddAssign());
    reduce(bunchToMatStat_m, bunchToMatStat_m, OpAddAssign());
    reduce(redifusedStat_m, redifusedStat_m, OpAddAssign());
    reduce(stoppedPartStat_m, stoppedPartStat_m, OpAddAssign());

    /*
      Degrader   *deg  = NULL;
      deg = dynamic_cast<Degrader *>(element_ref_m);
      double zBegin, zEnd;
      deg->getDimensions(zBegin, zEnd);
    */

    msg << std::scientific;
    msg << "Name " << FN_m
        << " material " << material_m << " particles in material " << globalPartsInMat_m << endl;
    msg << collShapeStr_m
        << " stats: bunch to material " << bunchToMatStat_m << " redifused " << redifusedStat_m
        << " stopped " << stoppedPartStat_m << endl;

    msg.flags(ff);
}

bool CollimatorPhysics::stillActive() { return globalPartsInMat_m != 0;}

bool CollimatorPhysics::stillAlive(PartBunch &bunch) {

    bool degraderAlive = true;

    //free GPU memory in case element is degrader, it is empty and bunch has moved past it
    if(collShape_m == DEGRADER && globalPartsInMat_m == 0) {
        Degrader   *deg  = NULL;
        deg = dynamic_cast<Degrader *>(element_ref_m);

        //get the size of the degrader
        double zBegin, zEnd;
        deg->getDimensions(zBegin, zEnd);

        //get the average Z position of the bunch
        Vector_t bunchOrigin = bunch.get_origin();

        //if bunch has moved past degrader free GPU memory
        if (bunchOrigin[2] > zBegin) {
            degraderAlive = false;
#ifdef OPAL_DKS
            if (IpplInfo::DKSEnabled)
                clearCollimatorDKS();
#endif
        }
    }

    return degraderAlive;

}


bool myCompF(PART x, PART y) {
    return x.label > y.label;
}

void CollimatorPhysics::deleteParticleFromLocalVector() {
    /*
      the particle to be deleted (label < 0) are all at the end of
      the vector.
    */
    sort(locParts_m.begin(),locParts_m.end(),myCompF);

    // find start of particles to delete
    std::vector<PART>::iterator inv = locParts_m.begin();

    for (; inv != locParts_m.end(); ++ inv) {
        if ((*inv).label == -1)
            break;
    }
    locParts_m.erase(inv,locParts_m.end());

    // update statistics
    if (locParts_m.size() > 0) {
        Eavg_m /= locParts_m.size();
        Emin_m /= locParts_m.size();
        Emax_m /= locParts_m.size();
    }
}

#ifdef OPAL_DKS

void CollimatorPhysics::applyDKS(PartBunch &bunch, size_t numParticlesInSimulation) {

    bool onlyOneLoopOverParticles = ! (allParticleInMat_m);

    if(collShape_m != DEGRADER) return;

    Degrader deg = dynamic_cast<Degrader *>(element_ref_m);

    //if firs call to apply setup needed accelerator resources
    setupCollimatorDKS(bunch, deg, numParticlesInSimulation);

    int numaddback;
    do {
        IpplTimings::startTimer(DegraderLoopTimer_m);

        //write particles to GPU if there are any to write
        if (dksParts_m.size() > 0) {
            //wrtie data from dksParts_m to the end of mem_mp (offset = numparticles)
            dksbase_m.writeDataAsync<PART_DKS>(mem_mp, &dksParts_m[0],
                                             dksParts_m.size(), -1, numparticles_m);

            //update number of particles on Device
            numparticles_m += dksParts_m.size();

            dksParts_m.erase(dksParts_m.begin(), dksParts_m.end());
        }

        //execute CollimatorPhysics kernels on GPU if any particles are there
        if (numparticles_m > 0) {
            dksbase_m.callCollimatorPhysics2(mem_mp, par_mp,
                                           numparticles_m, enableRutherfordScattering_m);
        }

        //sort device particles and get number of particles comming back to bunch
        numaddback = 0;
        if (numparticles_m > 0) {
            dksbase_m.callCollimatorPhysicsSort(mem_mp, numparticles_m, numaddback);
        }

        //read particles from GPU if any are comming out of material
        if (numaddback > 0) {

            //resize dksParts_m to hold particles that need to go back to bunch
            dksParts_m.resize(numaddback);

            //read particles that need to be added back to bunch
            //particles that need to be added back are at the end of Device array
            dksbase_m.readData<PART_DKS>(mem_mp, &dksParts_m[0], numaddback,
                                       numparticles_m - numaddback);

            //add particles back to the bunch
            for (unsigned int i = 0; i < dksParts_m.size(); ++ i) {
                if (dksParts_m[i].label == -2) {
                    addBackToBunchDKS(bunch, i);
                    ++ redifusedStat_m;
                } else {
                    ++ stoppedPartStat_m;
                    lossDs_m->addParticle(dksParts_m[i].Rincol, dksParts_m[i].Pincol,
                                          -locParts_m[dksParts_m[i].localID].IDincol);
                }
            }

            //erase particles that came from device from host array
            dksParts_m.erase(dksParts_m.begin(), dksParts_m.end());

            //update number of particles on Device
            numparticles_m -= numaddback;
        }

        IpplTimings::stopTimer(DegraderLoopTimer_m);

        if (onlyOneLoopOverParticles)
            copyFromBunchDKS(bunch);

        //bunch.boundp();

        T_m += dT_m;

        localPartsInMat_m = numparticles_m;
        reduce(localPartsInMat_m, globalPartsInMat_m, OpAddAssign());

        int maxPerNode = bunch.getLocalNum();
        reduce(maxPerNode, maxPerNode, OpMaxAssign());

        //more than one loop only if all the particles are in this degrader
        if (allParticleInMat_m) {
            onlyOneLoopOverParticles = ( (unsigned)maxPerNode > bunch.getMinimumNumberOfParticlesPerCore() ||
                                         globalPartsInMat_m == 0);
        } else {
            onlyOneLoopOverParticles = true;
        }

    } while (onlyOneLoopOverParticles == false);
}

bool myCompFDKS(PART_DKS x, PART_DKS y) {
    return x.label > y.label;
}

void CollimatorPhysics::addBackToBunchDKS(PartBunch &bunch, unsigned i) {

    int id = dksParts_m[i].localID;

    bunch.createWithID(locParts_m[id].IDincol);

    /*
      Binincol is still <0, but now the particle is rediffused
      from the material and hence this is not a "lost" particle anymore
    */
    bunch.Bin[bunch.getLocalNum()-1] = 1;

    bunch.R[bunch.getLocalNum()-1]           = dksParts_m[i].Rincol;
    bunch.P[bunch.getLocalNum()-1]           = dksParts_m[i].Pincol;

    bunch.Q[bunch.getLocalNum()-1]           = locParts_m[id].Qincol;
    bunch.LastSection[bunch.getLocalNum()-1] = locParts_m[id].LastSecincol;
    bunch.Bf[bunch.getLocalNum()-1]          = locParts_m[id].Bfincol;
    bunch.Ef[bunch.getLocalNum()-1]          = locParts_m[id].Efincol;
    bunch.dt[bunch.getLocalNum()-1]          = locParts_m[id].DTincol;

    dksParts_m[i].label = -1.0;

}

void CollimatorPhysics::copyFromBunchDKS(PartBunch &bunch)
{
    Degrader   *deg  = NULL;
    Collimator *coll = NULL;

    if(collShape_m == DEGRADER)
        deg = dynamic_cast<Degrader *>(element_ref_m);
    else
        coll = dynamic_cast<Collimator *>(element_ref_m);


    const size_t nL = bunch.getLocalNum();
    size_t ne = 0;
    std::list<size_t> partsToDel;
    const unsigned int minNumOfParticlesPerCore = bunch.getMinimumNumberOfParticlesPerCore();

    for(unsigned int i = 0; i < nL; ++ i) {
	if ((bunch.Bin[i]==-1 || bunch.Bin[i]==1) && ((nL-ne)>minNumOfParticlesPerCore)
	    && checkHit(bunch.R[i], bunch.P[i], dT_m, deg, coll))
            {

                PART x;
                x.localID      = nextLocalID_m; //unique id for each particle
                x.DTincol      = bunch.dt[i];
                x.IDincol      = bunch.ID[i];
                x.Binincol     = bunch.Bin[i];
                x.Rincol       = bunch.R[i];
                x.Pincol       = bunch.P[i];
                x.Qincol       = bunch.Q[i];
                x.LastSecincol = bunch.LastSection[i];
                x.Bfincol      = bunch.Bf[i];
                x.Efincol      = bunch.Ef[i];
                x.label        = 0;            // allive in matter

                PART_DKS x_gpu;
                x_gpu.label = x.label;
                x_gpu.localID = x.localID;
                x_gpu.Rincol = x.Rincol;
                x_gpu.Pincol = x.Pincol;

                locParts_m.push_back(x);
                dksParts_m.push_back(x_gpu);

                ++ ne;
                ++ bunchToMatStat_m;
                ++ nextLocalID_m;

                //mark particle to be deleted from bunch as soon as it enters the material
                partsToDel.push_front(i); // delete first those in the back of the container
            }
    }

    for (auto it = partsToDel.begin(); it != partsToDel.end(); ++ it) {
        bunch.destroy(1, *it);
    }
}

void CollimatorPhysics::setupCollimatorDKS(PartBunch &bunch, Degrader *deg,
					   size_t numParticlesInSimulation)
{

    if (curandInitSet_m == -1) {

        IpplTimings::startTimer(DegraderInitTimer_m);

        //int size = bunch.getLocalNum() + 0.5 * bunch.getLocalNum();
	//int size = bunch.getTotalNum() + 0.5 * bunch.getTotalNum();
	int size = numParticlesInSimulation;

        //allocate memory for parameters
        par_mp = dksbase_m.allocateMemory<double>(numpar, ierr_m);

        //allocate memory for particles
        mem_mp = dksbase_m.allocateMemory<PART_DKS>((int)size, ierr_m);

        maxparticles_m = (int)size;
        numparticles_m = 0;

        //reserve space for locParts_m vector
        locParts_m.reserve(size);

        //init curand
        dksbase_m.callInitRandoms(size, Options::seed);
        curandInitSet_m = 1;

        //create and transfer parameter array
        double zBegin, zEnd;
        deg->getDimensions(zBegin, zEnd);

        double params[numpar_ms] = {zBegin, deg->getZSize(), rho_m, Z_m,
                                    A_m, A2_c, A3_c, A4_c, A5_c, X0_m, I_m, dT_m, lowEnergyThr_m};
        dksbase_m.writeDataAsync<double>(par_mp, params, numpar_ms);

        IpplTimings::stopTimer(DegraderInitTimer_m);

    }

}

void CollimatorPhysics::clearCollimatorDKS() {
    if (curandInitSet_m == 1) {
        dksbase_m.freeMemory<double>(par_mp, numpar_ms);
        dksbase_m.freeMemory<PART_DKS>(mem_mp, maxparticles_m);
        curandInitSet_m = -1;
    }

}

void CollimatorPhysics::applyHost(PartBunch &bunch, Degrader *deg, Collimator *coll) {

    //loop trough particles in dksParts_m
    static const double m2mm = 1000.0;
    for (unsigned int i = 0; i < dksParts_m.size(); ++ i) {
        if(dksParts_m[i].label != -1) {
            Vector_t &R = dksParts_m[i].Rincol;
            Vector_t &P = dksParts_m[i].Pincol;
            double Eng = (sqrt(1.0  + dot(P, P)) - 1) * m_p;

            if(checkHit(R, P, dT_m, deg, coll)) {
                bool pdead = computeEnergyLoss(Eng, dT_m);

                if(!pdead) {

                    double ptot =  sqrt((m_p + Eng) * (m_p + Eng) - (m_p) * (m_p)) / m_p;
                    P = P * ptot / sqrt(dot(P, P));
                    /*
                      Now scatter and transport particle in material.
                      The checkInColl call just above will detect if the
                      particle is rediffused from the material into vacuum.
                    */

                    applyCoulombScat(R, P, dT_m);

                    dksParts_m[i].Rincol = R;
                    dksParts_m[i].Pincol = P;

                    calcStat(Eng);

                } else {
                    // The particle is stopped in the material, set lable_m to -1
                    dksParts_m[i].label = -1.0;
                    ++ stoppedPartStat_m;
                    lossDs_m->addParticle(R,P,-locParts_m[dksParts_m[i].localID].IDincol);
                }
            } else {
                /* The particle exits the material but is still in the loop of the substep,
                   Finish the timestep by letting the particle drift and after the last
                   substep call addBackToBunch
                */
                if(collShape_m == CYCLCOLLIMATOR) {
                    double gamma = (Eng + m_p) / m_p;
                    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
                    R = R + dT_m * beta * Physics::c * P / sqrt(dot(P, P)) * m2mm;
                } else {
                    R = R + dT_m * Physics::c * P / sqrt(1.0 + dot(P, P)) ;
                }

                addBackToBunchDKS(bunch, i);
                ++ redifusedStat_m;
            }
        }
    }

}

void CollimatorPhysics::deleteParticleFromLocalVectorDKS() {

    /*
      the particle to be deleted (label < 0) are all at the end of
      the vector.
    */
    sort(dksParts_m.begin(),dksParts_m.end(),myCompFDKS);

    // find start of particles to delete
    std::vector<PART_DKS>::iterator inv = dksParts_m.begin() + stoppedPartStat_m + redifusedStat_m;

    /*
      for (; inv != dksParts_m.end(); ++ inv) {
      if ((*inv).label == -1)
      break;
      }
    */

    dksParts_m.erase(inv,dksParts_m.end());

}

#endif

Vector_t ArbitraryRotation(Vector_t &W, Vector_t &Rorg, double Theta) {
    double C=cos(Theta);
    double S=sin(Theta);
    W = W / sqrt(dot(W,W));
    return Rorg * C + cross(W,Rorg) * S + W * dot(W,Rorg) * (1.0-C);
}