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


CollimatorPhysics::CollimatorPhysics(const std::string &name, ElementBase *element, std::string &material):
    SurfacePhysicsHandler(name, element),
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
    dT_m(0.0),
    T_m(0.0),
    allParticlesIn_m(false)
{


    gsl_rng_env_setup();
    rGen_m = gsl_rng_alloc(gsl_rng_default);

    if(dynamic_cast<Collimator *>(element_ref_m)) {
        Collimator *coll = dynamic_cast<Collimator *>(element_ref_m);
        FN_m = coll->getName();
        collshape_m = coll->getCollimatorShape();
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
        collshape_m = deg->getDegraderShape();
    }
    // statistics counters
    bunchToMatStat_m = 0;
    stoppedPartStat_m = 0;
    redifusedStat_m   = 0;

    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(FN_m, !Options::asciidump));

#ifdef OPAL_DKS
    dksbase.setAPI("Cuda", 4);
    dksbase.setDevice("-gpu", 4);
    dksbase.initDevice();
    curandInitSet = -1;
#endif

    DegraderApplyTimer_m = IpplTimings::getTimer("DegraderApply");
    DegraderLoopTimer_m = IpplTimings::getTimer("DegraderLoop");
    DegraderInitTimer_m = IpplTimings::getTimer("DegraderInit");
}

CollimatorPhysics::~CollimatorPhysics() {
    locParts_m.clear();
    lossDs_m->save();
    if (rGen_m)
        gsl_rng_free(rGen_m);

#ifdef OPAL_DKS
    clearCollimatorDKS();
#endif

}

bool CollimatorPhysics::checkHit(Vector_t R, Vector_t P, double dt, Degrader *deg, Collimator *coll) {
    bool hit = false;
    if(collshape_m == "CCollimator")
        hit = coll->checkPoint(R(0),R(1));
    else if (collshape_m == "DEGRADER") {
        hit = deg->isInMaterial(R(2));
    }
    else
        hit = coll->isInColl(R,P,Physics::c * dt/sqrt(1.0  + dot(P, P)));
    return hit;
}

void CollimatorPhysics::apply(PartBunch &bunch) {
    IpplTimings::startTimer(DegraderApplyTimer_m);

    Inform m ("CollimatorPhysics::apply ");
    /*
      Particles that have entered material are flagged as Bin[i] == -1.
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

    bunchToMatStat_m  = 0;
    redifusedStat_m   = 0;
    stoppedPartStat_m = 0;
    locPartsInMat_m   = 0;

    bool onlyOneLoopOverParticles = true;

    dT_m = bunch.getdT();
    T_m  = bunch.getT();

    /*
      Because this is not propper set in the Component class when calling in the Constructor
    */
    Degrader   *deg  = NULL;
    Collimator *coll = NULL;

    if(collshape_m == "DEGRADER") {
        deg = dynamic_cast<Degrader *>(element_ref_m);
    }
    else {
        coll = dynamic_cast<Collimator *>(element_ref_m);
    }

    if(allParticleInMat_m) {
        /*
          WE ARE HERE because on one node only 1 particles,
          not in the material i.e. in the bunch are left.

          Now integrate particles in the matierial
          until on every node more than 2 partciles - out
          of the material - exist.
        */
        m << "All Particles are in the material ... " << endl;
        onlyOneLoopOverParticles = false;
    }

#ifdef OPAL_DKS

    //if firs call to apply setup needed accelerator resources
    setupCollimatorDKS(bunch, deg);

    int numaddback;
    do {

        IpplTimings::startTimer(DegraderLoopTimer_m);

        //write particles to GPU if there are any to write
        if (dksParts_m.size() >= 0) {
            //wrtie data from dksParts_m to the end of mem_ptr (offset = numparticles)
            dksbase.writeDataAsync<PART_DKS>(mem_ptr, &dksParts_m[0],
                                             dksParts_m.size(), -1, numparticles);

            //update number of particles on Device
            numparticles += dksParts_m.size();

            //free locParts_m vector
            dksParts_m.erase(dksParts_m.begin(), dksParts_m.end());
        }

        //execute CollimatorPhysics kernels on GPU if any particles are there
        if (numparticles > 0) {
            dksbase.callCollimatorPhysics2(mem_ptr, par_ptr, numparticles);
        }

        //execute CollimatorPhysics loop and sort on particles that are on host
        if (dksParts_m.size() > 0) {
            applyHost(bunch, deg, coll);
            if (stoppedPartStat_m > 0 || redifusedStat_m > 0)
                deleteParticleFromLocalVectorDKS();
        }

        //sort device particles and get number of particles comming back to bunch
        numaddback = 0;
        if (numparticles > 0) {
            dksbase.callCollimatorPhysicsSort(mem_ptr, numparticles, numaddback);
        }

        //read particles from GPU if any are comming out of material
        if (numaddback > 0) {
            //get particles that are in material on host
            int parthost = dksParts_m.size();

            //resize dksParts_m to hold particles that need to go back to bunch
            dksParts_m.resize(parthost + numaddback);

            //read particles that need to be added back to bunch
            //particles that need to be added back are at the end of Device array
            dksbase.readData<PART_DKS>(mem_ptr, &dksParts_m[parthost], numaddback,
                                       numparticles - numaddback);

            //add particles back to the bunch
            for (unsigned int i = parthost; i < dksParts_m.size(); ++i) {
                if (dksParts_m[i].label == -2) {
                    addBackToBunchDKS(bunch, i);
                    redifusedStat_m++;
                } else {
                    stoppedPartStat_m++;
                    lossDs_m->addParticle(dksParts_m[i].Rincol, dksParts_m[i].Pincol,
                                          -locParts_m[dksParts_m[i].localID].IDincol);
                }
            }

            //erase particles that came from device from host array
            dksParts_m.erase(dksParts_m.begin() + parthost, dksParts_m.end());

            //update number of particles on Device
            numparticles -= numaddback;
        }

        IpplTimings::stopTimer(DegraderLoopTimer_m);

        bunch.boundp();

        if (onlyOneLoopOverParticles)
            copyFromBunchDKS(bunch);

        T_m += dT_m;

        if (!onlyOneLoopOverParticles) {
            bunch.gatherLoadBalanceStatistics();
            onlyOneLoopOverParticles = (bunch.getMinLocalNum() > 1);
        }
    } while (onlyOneLoopOverParticles == false);

    if (Ippl::myNode() == 0)
        nvtxRangePop();
#else

    do{
        IpplTimings::startTimer(DegraderLoopTimer_m);
        for(unsigned int i = 0; i < locParts_m.size(); ++i) {
            if(locParts_m[i].label != -1) {
                bool pdead = false;
                Vector_t &R = locParts_m[i].Rincol;
                Vector_t &P = locParts_m[i].Pincol;
                double Eng = (sqrt(1.0  + dot(P, P)) - 1) * m_p;

                if(checkHit(R,P,dT_m, deg, coll)) {
                    EnergyLoss(Eng, pdead, dT_m);

                    if(!pdead) {

                        double ptot =  sqrt((m_p + Eng) * (m_p + Eng) - (m_p) * (m_p)) / m_p;
                        P = P * ptot / sqrt(dot(P, P));
                        /*
                          Now scatter and transport particle in material.
                          The checkInColl call just above will detect if the
                          particle is rediffused from the material into vacuum.
                        */

                        CoulombScat(R, P, dT_m);

                        locParts_m[i].Rincol = R;
                        locParts_m[i].Pincol = P;

                        calcStat(Eng);

                    } else {
                        // The particle is stopped in the material, set lable_m to -1
                        locParts_m[i].label = -1.0;
                        stoppedPartStat_m++;
                        lossDs_m->addParticle(R,P,-locParts_m[i].IDincol);
                    }
                } else {
                    /* The particle exits the material but is still in the loop of the substep,
                       Finish the timestep by letting the particle drift and after the last
                       substep call addBackToBunch
                    */
                    double gamma = (Eng + m_p) / m_p;
                    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
                    if(collshape_m == "CCollimator") {
                        R = R + dT_m * beta * Physics::c * P / sqrt(dot(P, P)) * 1000;
                    } else {
                        locParts_m[i].Rincol = locParts_m[i].Rincol + dT_m * Physics::c * P / sqrt(1.0+dot(P, P)) ;
                        addBackToBunch(bunch, i);
                        redifusedStat_m++;
                    }
                }
            }
        }

        /*
          delete absorbed particles and particles that went to the bunch
        */
        deleteParticleFromLocalVector();

        IpplTimings::stopTimer(DegraderLoopTimer_m);

        /*
          because we have particles going back from material to the bunch
        */
        bunch.boundp();

        /*
          if we are not looping copy newly arrived particles
        */
        if (onlyOneLoopOverParticles)
            copyFromBunch(bunch);

        T_m += dT_m;              // update local time

        if (!onlyOneLoopOverParticles) {
            bunch.gatherLoadBalanceStatistics();
            onlyOneLoopOverParticles = (bunch.getMinLocalNum() > 1);
            //	m << "npl,npt= " << bunch.getLocalNum() << "," << bunch.getTotalNum() << " npm= " << locParts_m.size() << " redifused = " << redifusedStat_m << endl;
        }
    } while (onlyOneLoopOverParticles == false);

#endif

    IpplTimings::stopTimer(DegraderApplyTimer_m);
}


const std::string CollimatorPhysics::getType() const {
    return "CollimatorPhysics";
}

/// The material of the collimator
//  ------------------------------------------------------------------------
void  CollimatorPhysics::Material() {

    if(material_m == "Berilium") {
        Z_m = 4.0;
        A_m = 9.012;
        rho_m = 1.848;

        X0_m = 65.19 / rho_m / 100;
        I_m = 10 * Z_m;
        n_m = rho_m / A_m * Avo;

        A2_c = 2.590;
        A3_c = 9.660E2;
        A4_c = 1.538E2;
        A5_c = 3.475E-2;

    }

    if(material_m == "Graphite") {
        Z_m = 6.0;
        A_m = 12.0107;
        rho_m = 2.210;

        X0_m = 42.70 / rho_m / 100;
        I_m = 10 * Z_m;
        n_m = rho_m / A_m * Avo;

        A2_c = 2.601;
        A3_c = 1.701E3;
        A4_c = 1.279E3;
        A5_c = 1.638E-2;

    }

    if(material_m == "GraphiteR6710") {
        Z_m = 6.0;
        A_m = 12.0107;
        rho_m = 1.88;

        X0_m = 42.70 / rho_m / 100;
        I_m = 10 * Z_m;
        n_m = rho_m / A_m * Avo;

        A2_c = 2.601;
        A3_c = 1.701E3;
        A4_c = 1.279E3;
        A5_c = 1.638E-2;

    }

    if(material_m == "Molybdenum") {
        Z_m = 42.0;
        A_m = 95.94;
        rho_m = 10.22;

        X0_m = 9.8 / rho_m / 100;
        I_m = 10 * Z_m;
        n_m = rho_m / A_m * Avo;

        A2_c = 7.248;
        A3_c = 9.545E3;
        A4_c = 4.802E2;
        A5_c = 5.376E-3;
    }

    /*
      needs to be checked !

      Z from http://journals.aps.org/prb/pdf/10.1103/PhysRevB.40.8530
    */

    if(material_m == "Mylar") {
        Z_m = 6.702;
        A_m = 12.88;
        rho_m = 1.4;

        X0_m = 39.95 / rho_m / 100;
        I_m = 10 * Z_m;
        n_m = rho_m / A_m * Avo;
	A2_c = 3.350;
	A3_c = 1683;
        A4_c = 1900;
        A5_c = 2.513e-02;
    }


    //new materials
    if (material_m == "Aluminum") {
        Z_m = 13;
        A_m = 26.98;
        rho_m = 2.7;

        X0_m = 24.01 / rho_m / 100;
        I_m = 10 * Z_m;
        n_m = rho_m / A_m * Avo;

        A2_c = 4.739;
        A3_c = 2.766e3;
        A4_c = 1.645e2;
        A5_c = 2.023e-2;
    }

    if (material_m == "Copper") {
        Z_m = 29;
        A_m = 63.54;
        rho_m = 8.96;

        X0_m = 12.86 / rho_m / 100;
        I_m = 10 * Z_m;
        n_m = rho_m / A_m * Avo;

        A2_c = 4.194;
        A3_c = 4.649e3;
        A4_c = 8.113e1;
        A5_c = 2.242e-2;
    }

    if (material_m == "Titan") {
        Z_m = 22;
        A_m = 47.8;
        rho_m = 4.54;

        X0_m = 16.16 / rho_m / 100;
        I_m = 10 * Z_m;
        n_m = rho_m / A_m * Avo;

        A2_c = 5.489;
        A3_c = 5.260e3;
        A4_c = 6.511e2;
        A5_c = 8.930e-3;
    }

    if (material_m == "AluminaAl2O3") {
        Z_m = 50;
        A_m = 101.96;
        rho_m = 3.97;

        X0_m = 27.94 / rho_m / 100;
        I_m = 10 * Z_m;
        n_m = rho_m / A_m * Avo;

        A2_c = 7.227;
        A3_c = 1.121e4;
        A4_c = 3.864e2;
        A5_c = 4.474e-3;
    }

    if (material_m == "Air") {
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

    if (material_m == "Kapton") {
        Z_m = 6;
        A_m = 12;
        rho_m = 1.4;

        X0_m = 39.95 / rho_m / 100;
        I_m = 10 * Z_m;
        n_m = rho_m / A_m * Avo;

        A2_c = 2.601;
        A3_c = 1.701e3;
        A4_c = 1.279e3;
        A5_c = 1.638e-2;
    }

    if (material_m == "Gold") {
        Z_m = 79;
        A_m = 197;
        rho_m = 19.3;

        X0_m = 6.46 / rho_m / 100;
        I_m = 10 * Z_m;
        n_m = rho_m / A_m * Avo;

        A2_c = 5.458;
        A3_c = 7.852e3;
        A4_c = 9.758e2;
        A5_c = 2.077e-2;
    }

    if (material_m == "Water") {
        Z_m = 10;
        A_m = 18;
        rho_m = 1;

        X0_m = 36.08 / rho_m / 100;
        I_m = 10 * Z_m;
        n_m = rho_m / A_m * Avo;

        A2_c = 2.199;
        A3_c = 2.393e3;
        A4_c = 2.699e3;
        A5_c = 1.568e-2;
    }

    // mean exitation energy from Leo
    if (Z_m < 13.0)
        I_m = 12 * Z_m + 7.0;
    else
        I_m = 9.76 * Z_m + (Z_m * 58.8 * std::pow(Z_m, -1.19));
}

/// Energy Loss:  using the Bethe-Bloch equation.
/// Energy straggling: For relatively thick absorbers such that the number of collisions is large,
/// the energy loss distribution is shown to be Gaussian in form.
// -------------------------------------------------------------------------

void  CollimatorPhysics::EnergyLoss(double &Eng, bool &pdead, double &deltat) {
    /// Eng GeV

    Material();
    double dEdx = 0.0;
    const double gamma = (Eng + m_p) / m_p;
    const double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    const double gamma2 = gamma * gamma;
    const double beta2 = beta * beta;

    const double deltas = deltat * beta * Physics::c;
    const double deltasrho = deltas * 100 * rho_m;
    const double K = 4.0 * pi * Avo * r_e * r_e * m_e * 1E7;
    const double sigma_E = sqrt(K * m_e * rho_m * (Z_m/A_m)* deltas * 1E5);

    if ((Eng > 0.00001) && (Eng < 0.0006)) {
        const double Ts = (Eng*1E6)/1.0073; // 1.0073 is the proton mass divided by the atomic mass number. T is in KeV
        const double epsilon_low = A2_c*pow(Ts,0.45);
        const double epsilon_high = (A3_c/Ts)*log(1+(A4_c/Ts)+(A5_c*Ts));
        const double epsilon = (epsilon_low*epsilon_high)/(epsilon_low + epsilon_high);
        dEdx = - epsilon /(1E21*(A_m/Avo)); // Stopping_power is in MeV INFOMSG("stopping power: " << dEdx << " MeV" << endl);
        const double delta_Eave = deltasrho * dEdx;
        const double delta_E = delta_Eave + gsl_ran_gaussian(rGen_m,sigma_E);
        Eng = Eng + delta_E / 1E3;
    }

    if (Eng >= 0.0006) {
        const double Tmax = 2.0 * m_e * 1e9 * beta2 * gamma2 /
            (1.0 + 2.0 * gamma * m_e / m_p + (m_e / m_p) * (m_e / m_p));
        dEdx = -K * z_p * z_p * Z_m / (A_m * beta2) *
            (1.0 / 2.0 * std::log(2 * m_e * 1e9 * beta2 * gamma2 * Tmax / I_m / I_m) - beta2);

        // INFOMSG("stopping power_BB: " << dEdx << " MeV" << endl);
        const double delta_Eave = deltasrho * dEdx;
	double tmp = gsl_ran_gaussian(rGen_m,sigma_E);
        const double delta_E = delta_Eave + tmp;
        Eng = Eng+delta_E / 1E3;
    }

    //INFOMSG("final energy: " << Eng/1000 << " MeV" <<endl);
    pdead = ((Eng<1E-4) || (dEdx>0));
}

// Implement the rotation in 2 dimensions here
// For details: see J. Beringer et al. (Particle Data Group), Phys. Rev. D 86, 010001 (2012), "Passage of particles through matter"
void  CollimatorPhysics::Rot(double &px, double &pz, double &x, double &z, double xplane, double normP, double thetacou, double deltas, int coord) {
    double Psixz;
    double pxz;
    // Calculate the angle between the px and pz momenta to change from beam coordinate to lab coordinate
    if (px>=0 && pz>=0) Psixz = atan(px/pz);
    else if (px>0 && pz<0)
	Psixz = atan(px/pz) + Physics::pi;
    else if (px<0 && pz>0)
	Psixz = atan(px/pz) + 2*Physics::pi;
    else
	Psixz = atan(px/pz) + Physics::pi;

    // Apply the rotation about the random angle thetacou & change from beam
    // coordinate system to the lab coordinate system using Psixz (2 dimensions)
    pxz = sqrt(px*px + pz*pz);
    if(coord==1) {
    	x = x + deltas * px/normP + xplane*cos(Psixz);
    	z = z - xplane * sin(Psixz);
    }
    if(coord==2) {
	x = x + deltas * px/normP + xplane*cos(Psixz);
	z = z - xplane * sin(Psixz) + deltas * pz / normP;
    }
    px = pxz*cos(Psixz)*sin(thetacou) + pxz*sin(Psixz)*cos(thetacou);
    pz = -pxz*sin(Psixz)*sin(thetacou) + pxz*cos(Psixz)*cos(thetacou);
}

/// Coulomb Scattering: Including Multiple Coulomb Scattering and large angle Rutherford Scattering.
/// Using the distribution given in Classical Electrodynamics, by J. D. Jackson.
//--------------------------------------------------------------------------
void  CollimatorPhysics::CoulombScat(Vector_t &R, Vector_t &P, double &deltat) {
    Material();
    double Eng = sqrt(dot(P, P) + 1.0) * m_p - m_p;
    double gamma = (Eng + m_p) / m_p;
    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    double normP = sqrt(dot(P, P));
    double deltas = deltat * beta * Physics::c;
    double theta0 = 13.6e6 / (beta * sqrt(dot(P, P)) * m_p * 1e9) * z_p * sqrt(deltas / X0_m) * (1.0 + 0.038 * log(deltas / X0_m));

    // x-direction: See Physical Review, "Multiple Scattering"
    double z1 = gsl_ran_gaussian(rGen_m,1.0);
    double z2 = gsl_ran_gaussian(rGen_m,1.0);
    double thetacou = z2 * theta0;

    while(fabs(thetacou) > 3.5 * theta0) {
        z1 = gsl_ran_gaussian(rGen_m,1.0);
        z2 = gsl_ran_gaussian(rGen_m,1.0);
        thetacou = z2 * theta0;
    }
    double xplane = z1 * deltas * theta0 / sqrt(12.0) + z2 * deltas * theta0 / 2.0;
    int coord = 1; // Apply change in coordinates for multiple scattering but not for Rutherford scattering (take the deltas step only once each turn)
    Rot(P(0),P(2),R(0),R(2), xplane, normP, thetacou, deltas, coord);


    // Rutherford-scattering in x-direction
    if(collshape_m == "CCollimator")
        R = R * 1000.0;

    double P2 = gsl_rng_uniform(rGen_m);
    if(P2 < 0.0047) {
        double P3 = gsl_rng_uniform(rGen_m);
        double thetaru = 2.5 * sqrt(1 / P3) * sqrt(2.0) * theta0;
        double P4 = gsl_rng_uniform(rGen_m);
        if(P4 > 0.5)
            thetaru = -thetaru;
	coord = 0; // no change in coordinates but one in momenta-direction
	Rot(P(0),P(2),R(0),R(2), xplane, normP, thetaru, deltas, coord);
    }

    // y-direction: See Physical Review, "Multiple Scattering"
    z1 = gsl_ran_gaussian(rGen_m,1.0);
    z2 = gsl_ran_gaussian(rGen_m,1.0);
    thetacou = z2 * theta0;

    while(fabs(thetacou) > 3.5 * theta0) {
        z1 = gsl_ran_gaussian(rGen_m,1.0);
        z2 = gsl_ran_gaussian(rGen_m,1.0);
        thetacou = z2 * theta0;
    }
    double yplane = z1 * deltas * theta0 / sqrt(12.0) + z2 * deltas * theta0 / 2.0;
    coord = 2; // Apply change in coordinates for multiple scattering but not for Rutherford scattering (take the deltas step only once each turn)
    Rot(P(1),P(2),R(1),R(2), yplane, normP, thetacou, deltas, coord);

    // Rutherford-scattering in x-direction
    if(collshape_m == "CCollimator")
        R = R * 1000.0;

    P2 = gsl_rng_uniform(rGen_m);
    if(P2 < 0.0047) {
        double P3 = gsl_rng_uniform(rGen_m);
        double thetaru = 2.5 * sqrt(1 / P3) * sqrt(2.0) * theta0;
        double P4 = gsl_rng_uniform(rGen_m);
        if(P4 > 0.5)
            thetaru = -thetaru;
	coord = 0; // no change in coordinates but one in momenta-direction
	Rot(P(1),P(2),R(1),R(2), yplane, normP, thetaru, deltas, coord);
    }
}

void CollimatorPhysics::addBackToBunch(PartBunch &bunch, unsigned i) {

    bunch.createWithID(locParts_m[i].IDincol);

    /*
      Binincol is still <0, but now the particle is rediffused
      from the material and hence this is not a "lost" particle anymore
    */
    bunch.Bin[bunch.getLocalNum()-1] = -1*locParts_m[i].Binincol;

    bunch.R[bunch.getLocalNum()-1]           = locParts_m[i].Rincol;
    bunch.P[bunch.getLocalNum()-1]           = locParts_m[i].Pincol;
    bunch.Q[bunch.getLocalNum()-1]           = locParts_m[i].Qincol;
    bunch.LastSection[bunch.getLocalNum()-1] = locParts_m[i].LastSecincol;
    bunch.Bf[bunch.getLocalNum()-1]          = locParts_m[i].Bfincol;
    bunch.Ef[bunch.getLocalNum()-1]          = locParts_m[i].Efincol;
    bunch.dt[bunch.getLocalNum()-1]          = locParts_m[i].DTincol;

    /*
      This particle is back to the bunch, by set
      ting the lable to -1.0
      the particle will be deleted.
    */
    locParts_m[i].label = -1.0;
}

void CollimatorPhysics::copyFromBunch(PartBunch &bunch)
{
    for(unsigned int i = 0; i < bunch.getLocalNum(); ++i) {
        if (bunch.Bin[i]<0) {
            PART x;
            x.localID      = i;
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

            locParts_m.push_back(x);

            bunchToMatStat_m++;
        }
    }
}


void CollimatorPhysics::print(Inform &msg){
    Inform::FmtFlags_t ff = msg.flags();

    // ToDo: need to move that to a statistics function
#ifdef OPAL_DKS
    locPartsInMat_m = numparticles + dksParts_m.size();
#else
    locPartsInMat_m = locParts_m.size();
#endif
    reduce(locPartsInMat_m,locPartsInMat_m, OpAddAssign());
    reduce(bunchToMatStat_m,bunchToMatStat_m, OpAddAssign());
    reduce(redifusedStat_m,redifusedStat_m, OpAddAssign());
    reduce(stoppedPartStat_m,stoppedPartStat_m, OpAddAssign());

    msg << std::scientific;
    msg << "--- CollimatorPhysics - Type is " << collshape_m << " Name " << FN_m
        << " Material " << material_m << " Particles in material " << locPartsInMat_m << endl;
    msg << "Coll/Deg statistics: "
        << " bunch to material " << bunchToMatStat_m << " redifused " << redifusedStat_m
        << " stopped " << stoppedPartStat_m << endl;

    msg.flags(ff);
}

bool CollimatorPhysics::stillActive() { return locPartsInMat_m != 0;}


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

    for (; inv != locParts_m.end(); inv++) {
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
    bunch.Bin[bunch.getLocalNum()-1] = -1*locParts_m[id].Binincol;

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
    for(unsigned int i = 0; i < bunch.getLocalNum(); ++i) {
        if (bunch.Bin[i]<0) {
            PART x;
            x.localID      = numlocalparts; //unique id for each particle
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

            bunchToMatStat_m++;
            numlocalparts++;
        }
    }
}

void CollimatorPhysics::setupCollimatorDKS(PartBunch &bunch, Degrader *deg) {

    if (curandInitSet == -1) {

        IpplTimings::startTimer(DegraderInitTimer_m);

        int size = bunch.getLocalNum() + 0.05 * bunch.getLocalNum();
        //allocate memory for parameters
        par_ptr = dksbase.allocateMemory<double>(numpar, ierr);

        //allocate memory for particles
        mem_ptr = dksbase.allocateMemory<PART_DKS>(size, ierr);

        maxparticles = size;
        numparticles = 0;
        numlocalparts = 0;

        //reserve space for locParts_m vector
        locParts_m.reserve(size);

        //init curand
        dksbase.callInitRandoms(size);
        curandInitSet = 1;

        //create and transfer parameter array
        Material();
        double zBegin, zEnd;
        deg->getDimensions(zBegin, zEnd);

        double params[numpar] = {zBegin, deg->getZSize(), rho_m, Z_m,
                                 A_m, A2_c, A3_c, A4_c, A5_c, X0_m, I_m, dT_m};
        dksbase.writeDataAsync<double>(par_ptr, params, numpar);

        IpplTimings::stopTimer(DegraderInitTimer_m);

    }

}

void CollimatorPhysics::clearCollimatorDKS() {

    if (curandInitSet == 1) {
        dksbase.freeMemory<double>(par_ptr, numpar);
        dksbase.freeMemory<PART_DKS>(mem_ptr, maxparticles);
        curandInitSet = -1;
    }

}

void CollimatorPhysics::applyHost(PartBunch &bunch, Degrader *deg, Collimator *coll) {

    //loop trough particles in dksParts_m
    for (unsigned int i = 0; i < dksParts_m.size(); ++i) {
        if(dksParts_m[i].label != -1) {
            bool pdead = false;
            Vector_t &R = dksParts_m[i].Rincol;
            Vector_t &P = dksParts_m[i].Pincol;
            double Eng = (sqrt(1.0  + dot(P, P)) - 1) * m_p;

            if(checkHit(R,P,dT_m, deg, coll)) {
                EnergyLoss(Eng, pdead, dT_m);

                if(!pdead) {

                    double ptot =  sqrt((m_p + Eng) * (m_p + Eng) - (m_p) * (m_p)) / m_p;
                    P = P * ptot / sqrt(dot(P, P));
                    /*
                      Now scatter and transport particle in material.
                      The checkInColl call just above will detect if the
                      particle is rediffused from the material into vacuum.
                    */

                    CoulombScat(R, P, dT_m);

                    dksParts_m[i].Rincol = R;
                    dksParts_m[i].Pincol = P;

                    calcStat(Eng);

                } else {
                    // The particle is stopped in the material, set lable_m to -1
                    dksParts_m[i].label = -1.0;
                    stoppedPartStat_m++;
                    lossDs_m->addParticle(R,P,-locParts_m[dksParts_m[i].localID].IDincol);
                }
            } else {
                /* The particle exits the material but is still in the loop of the substep,
                   Finish the timestep by letting the particle drift and after the last
                   substep call addBackToBunch
                */
                double gamma = (Eng + m_p) / m_p;
                double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
                if(collshape_m == "CCollimator") {
                    R = R + dT_m * beta * Physics::c * P / sqrt(dot(P, P)) * 1000;
                } else {
                    dksParts_m[i].Rincol = dksParts_m[i].Rincol + dT_m * Physics::c * P / sqrt(1.0+dot(P, P)) ;
                    addBackToBunchDKS(bunch, i);
                    redifusedStat_m++;
                }
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
      for (; inv != dksParts_m.end(); inv++) {
      if ((*inv).label == -1)
      break;
      }
    */

    dksParts_m.erase(inv,dksParts_m.end());

}

#endif


