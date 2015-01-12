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
    allParticlesIn_m(false) {
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
}

CollimatorPhysics::~CollimatorPhysics() {
  locParts_m.clear();
  lossDs_m->save();
  if (rGen_m)
    gsl_rng_free(rGen_m);
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
    do{
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
	because we have particles going back from material to the bunch
      */
      bunch.boundp();

      /*
	delete absorbed particles and particles that went to the bunch
      */
      deleteParticleFromLocalVector();
      
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
}

const std::string CollimatorPhysics::getType() const {
    return "CollimatorPhysics";
}

/// The material of the collimator
//  ------------------------------------------------------------------------
void  CollimatorPhysics::Material() {

    if(material_m == "Cu") {
        Z_m = 29.0;
        A_m = 63.546;
        rho_m = 8.96;

        X0_m = 12.86 / rho_m / 100;
        I_m = 10. * Z_m;
        n_m = rho_m / A_m * Avo;

        A2_c = 4.194;
        A3_c = 4.649E3;
        A4_c = 8.113E1;
        A5_c = 2.42E-2;
    }

    if(material_m == "Be") {
        Z_m = 4.0;
        A_m = 9.012;
        rho_m = 1.848;

        X0_m = 65.19 / rho_m / 100;
        I_m = 10 * Z_m;
        n_m = rho_m / A_m * Avo;

        A2_c = 2.590;
        A3_c = 9.660E2;
        A4_c = 1.538E2;
        A5_c =3.475E-2;

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

    if(material_m == "Mo") {
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

    if(material_m == "My") {
        Z_m = 4.54;
        A_m = 12.88;       
        rho_m = 1.4;

        X0_m = 9.8 / rho_m / 100;
        I_m = 10 * Z_m;
        n_m = rho_m / A_m * Avo;

        A2_c = 1.0;
        A3_c = 1.0;
        A4_c = 1.0;
        A5_c = 1.0;
    }

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
        const double delta_E = delta_Eave + gsl_ran_gaussian(rGen_m,sigma_E);
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
  locPartsInMat_m = locParts_m.size();
  reduce(locPartsInMat_m,locPartsInMat_m, OpAddAssign());
  reduce(bunchToMatStat_m,bunchToMatStat_m, OpAddAssign());
  reduce(redifusedStat_m,redifusedStat_m, OpAddAssign());
  reduce(stoppedPartStat_m,stoppedPartStat_m, OpAddAssign());

  msg << std::scientific;
  msg << "--- CollimatorPhysics - Type is " << collshape_m
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
