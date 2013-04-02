// Class:CollimatorPhysics
// Defines the collimator physics models
// ------------------------------------------------------------------------
// Class category:
// ------------------------------------------------------------------------
// $Date: 2009/07/20 09:32:31 $
// $Author: Bi, Yang $
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
#include "Structure/LossDataSink.h"
#include "Distribution/ranlib.h"

#include <iostream>
#include <fstream>
#include <algorithm>

//using Physics::c;
using Physics::pi;
using Physics::m_p;
using Physics::m_e;
using Physics::h_bar;
using Physics::epsilon_0;
using Physics::r_e;
using Physics::z_p;
using Physics::Avo;


CollimatorPhysics::CollimatorPhysics(const string &name, ElementBase *element, const double &major, const double &minor,  string &material):
    SurfacePhysicsHandler(name, element),
    a_m(major),
    b_m(minor),
    xp_m(0.0),
    yp_m(0.0),
    xstart_m(0.0),
    xend_m(0.0),
    ystart_m(0.0),
    yend_m(0.0),
    zstart_m(0.0),
    zend_m(0.0),
    width_m(0.0),
    Begin_m(0.0),
    End_m(0.0),
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
    N_m(0),
    incoll_m(false),
    time_m(0.0) {
    rGen_m = new RANLIB_class(115314159, 4);
    // initialize DataSink with H5Part output enabled
    bool doH5 = false;
    lossDs_m = new LossDataSink(FN_m, doH5);

    if(dynamic_cast<Collimator *>(element_ref_m)) {
        Collimator *coll = dynamic_cast<Collimator *>(element_ref_m);
        FN_m = coll->getName();
        collshape_m = coll->getCollimatorShape();
        xp_m = coll->getXpos();
        yp_m = coll->getYpos();
        xstart_m = coll->getXStart();
        ystart_m = coll->getYStart();
        zstart_m = coll->getZStart();
        xend_m = coll->getXEnd();
        yend_m = coll->getYEnd();
        zend_m = coll->getZEnd();
        width_m = coll->getWidth();
        setCColimatorGeom();    
    } else if(dynamic_cast<Drift *>(element_ref_m)) {
        Drift *drf = dynamic_cast<Drift *>(element_ref_m);
        drf->getDimensions(Begin_m, End_m);
        FN_m = drf->getName();
    } else if(dynamic_cast<SBend *>(element_ref_m)) {
        SBend *sben = dynamic_cast<SBend *>(element_ref_m);
        Begin_m = sben->getStartElement();
        End_m = Begin_m + sben->getElementLength();
        FN_m = sben->getName();
    } else if(dynamic_cast<RBend *>(element_ref_m)) {
        RBend *rben = dynamic_cast<RBend *>(element_ref_m);
        Begin_m = rben->getStartElement();
        End_m = Begin_m + rben->getElementLength();
        FN_m = rben->getName();
    } else if(dynamic_cast<Multipole *>(element_ref_m)) {
        Multipole *quad = dynamic_cast<Multipole *>(element_ref_m);
        quad->getDimensions(Begin_m, End_m);
        FN_m = quad->getName();
    } else if(dynamic_cast<Degrader *>(element_ref_m)) {
        Degrader *deg = dynamic_cast<Degrader *>(element_ref_m);
        FN_m = deg->getName();
        collshape_m = deg->getDegraderShape();
	width_m = deg->getZSize();
    }
}

CollimatorPhysics::~CollimatorPhysics() {
 
    locParts_m.clear();
    if(lossDs_m)
        delete lossDs_m;
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

    bunchToMatStat_m = 0;
    stoppedPartStat_m = 0;
    redifusedStat_m   = 0;

    Eavg_m = 0.0;
    Emax_m = 0.0;
    Emin_m = 0.0;
    
    /// the deltat in collimator is 1/n of the main tracking timestep.


    dT_m = bunch.getdT();

    /*    
    while(dT_m > 1.01e-12)
        dT_m = dT_m / 10;
    N_m = bunch.getdT() / dT_m;

    */

    N_m = 1;

    time_m = bunch.getT();
   
    /*
      Because this is not propper set in the Component class when calling in the Constructor
    */

    Degrader *deg = NULL;
    Collimator *coll = NULL;

    if(collshape_m == "DEGRADER") {
        deg = dynamic_cast<Degrader *>(element_ref_m);
        deg->getDimensions(Begin_m, End_m);
        End_m = Begin_m + width_m;
    }
    else {
        coll = dynamic_cast<Collimator *>(element_ref_m);
        coll->getDimensions(Begin_m, End_m);
    }

    for(int ii = 0; ii < N_m; ++ii) {
        for(unsigned int i = 0; i < locParts_m.size(); ++i) {
            if(locParts_m[i].label != -1) {
                bool pdead = false;
                Vector_t &R = locParts_m[i].Rincol;
                Vector_t &P = locParts_m[i].Pincol;
                double Eng = (sqrt(1.0  + dot(P, P)) - 1) * m_p;

		if(collshape_m == "CCollimator") 
                    incoll_m = checkInColl(R);    // fixme this needs to go to the Collimator
                else if (collshape_m == "DEGRADER") {
                    incoll_m = deg->isInMaterial(R(2));
                }
		else
                    incoll_m = coll->isInColl(R,P,Physics::c * bunch.getdT()/sqrt(1.0  + dot(P, P)));

		if(incoll_m) {
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

                        Eavg_m += Eng;
                        if (Emin_m > Eng)
                            Emin_m = Eng;
                        if (Emax_m < Eng)
                            Emax_m = Eng;

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
                        //                        m << "2 Out of mat but still in the loop" << endl;
                        //locParts_m[i].Rincol = locParts_m[i].Rincol + dT_m * beta * Physics::c * P / sqrt(dot(P, P)) ;
                        locParts_m[i].Rincol = locParts_m[i].Rincol + dT_m * Physics::c * P / sqrt( 1+ dot(P, P)) ;

                        if (ii == N_m-1) {
			    addBackToBunch(bunch, i);
                            redifusedStat_m++;
                        }   
                    }
                }
            }
        }
    }


    
    /* 
       add new (lost particles) to local data structure
    */

    copyFromBunch(bunch); 
 
    /*
      delete absorbed particles
    */
    deleteParticleFromLocalVector();

}

const string CollimatorPhysics::getType() const {
    return "CollimatorPhysics";
}

/// The material of the collimator
//  ------------------------------------------------------------------------
void  CollimatorPhysics::Material() {

    if(material_m == "Cu") {
        Z_m = 29;
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
        Z_m = 4;
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
        Z_m = 6;
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

    if(material_m == "Mo") {
        Z_m = 42;
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

}

/// Energy Loss:  using the Bethe-Bloch equation.
/// Energy straggling: For relatively thick absorbers such that the number of collisions is large,
/// the energy loss distribution is shown to be Gaussian in form.
// -------------------------------------------------------------------------

void  CollimatorPhysics::EnergyLoss(double &Eng, bool &pdead, double &deltat) {
    /// Eng GeV

    Material();
    double dEdx = 0.0;
    double gamma = (Eng + m_p) / m_p;
    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    double gamma2 = gamma * gamma;
    double beta2 = beta * beta;
    
    double deltas = deltat * beta * Physics::c;
    double deltasrho = deltas * 100 * rho_m; 
    double K = 4.0 * pi * Avo * r_e * r_e * m_e * 1E7;
    double sigma_E = sqrt(K * m_e * rho_m * (Z_m/A_m))* deltas * 1E5;
    
    
     
    if ((Eng > 0.00001) && (Eng < 0.0006)) { 
        double Ts = (Eng*1E6)/1.0073; // 1.0073 is the proton mass divided by the atomic mass number. T is in KeV
        double epsilon_low = A2_c*pow(Ts,0.45);		
        double epsilon_high = (A3_c/Ts)*log(1+(A4_c/Ts)+(A5_c*Ts));
        double epsilon = (epsilon_low*epsilon_high)/(epsilon_low + epsilon_high);
        dEdx = - epsilon /(1E21*(A_m/Avo)); // Stopping_power is in MeV
        // INFOMSG("stopping power: " << dEdx << " MeV" << endl);
        double delta_Eave = deltasrho * dEdx;
        double delta_E = rGen_m->gauss(delta_Eave, sigma_E);
        Eng = Eng + delta_E / 1E3;
    }
    
    if (Eng >= 0.0006) {
        double Tmax = 2.0 * m_e * 1e9 * beta2 * gamma2 / 
                      (1.0 + 2.0 * gamma * m_e / m_p + (m_e / m_p) * (m_e / m_p));
               dEdx = -K * z_p * z_p * Z_m / (A_m * beta2) * 
                      (1.0 / 2.0 * std::log(2 * m_e * 1e9 * beta2 * gamma2 * Tmax / I_m / I_m) - beta2);

	       // INFOMSG("stopping power_BB: " << dEdx << " MeV" << endl);
        double delta_Eave = deltasrho * dEdx;
        double delta_E = rGen_m->gauss(delta_Eave, sigma_E);
        Eng = Eng+delta_E / 1E3;
    }	

    //INFOMSG("final energy: " << Eng/1000 << " MeV" <<endl);

    pdead = ((Eng<1E-4) || (dEdx>0));


}



///From beam coordinate to lab coordinate
// --------------------------------------------------------------------------
void  CollimatorPhysics::Rot(Vector_t &P, Vector_t Prot, double normP) {
    Vector_t A0, A1, A2;
    if(fabs(P(0)) < 1e-8) {
        A0 = Vector_t(1.0, 0.0, 0.0);
        A1 = Vector_t(0.0, P(2), -P(1));
        A2 = Vector_t(0.0, P(1), P(2));
    } else {
        A0 = Vector_t(-(P(1) * P(1) + P(2) * P(2)) / P(0), P(1), P(2));
        A1 = Vector_t(0.0, P(2), -P(1));
        A2 = Vector_t(P(0), P(1), P(2));
    }
    A0 /= sqrt(dot(A0, A0));
    A1 /= sqrt(dot(A1, A1));
    A2 /= sqrt(dot(A2, A2));

    P(0) = A0(0) * Prot(0) + A1(0) * Prot(1) + A2(0) * Prot(2);
    P(1) = A0(1) * Prot(0) + A1(1) * Prot(1) + A2(1) * Prot(2);
    P(2) = A0(2) * Prot(0) + A1(2) * Prot(1) + A2(2) * Prot(2);

    P = P * normP;
}

void  CollimatorPhysics::CoulombScat()
{}
void  CollimatorPhysics::EnergyLoss(double &Eng, bool &pdead)
{}

///Coulomb Scattering: Including Multiple Coulomb Scattering and large angle Rutherford Scattering.
///Using the distribution given in Classical Electrodynamics, by J. D. Jackson.
//--------------------------------------------------------------------------
void  CollimatorPhysics::CoulombScat(Vector_t &R, Vector_t &P, double &deltat) {
    Material();
    double Eng = sqrt(dot(P, P) + 1.0) * m_p - m_p;
    double gamma = (Eng + m_p) / m_p;
    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    double normP = sqrt(dot(P, P));

    double deltas = deltat * beta * Physics::c;

    double theta0 = 13.6e6 / (beta * sqrt(dot(P, P)) * m_p * 1e9) * z_p * sqrt(deltas / X0_m) * (1.0 + 0.038 * log(deltas / X0_m));

    double z1 = rGen_m->gauss(0.0, 1.0);
    double z2 = rGen_m->gauss(0.0, 1.0);
    double thetacou = z2 * theta0;
   
    while(fabs(thetacou) > 3.5 * theta0) {
        z1 = rGen_m->gauss(0.0, 1.0);
        z2 = rGen_m->gauss(0.0, 1.0);
        thetacou = z2 * theta0;
    }
    double poscou = z1 * deltas * theta0 / sqrt(12.0) + z2 * deltas * theta0 / 2.0;

    double z3 = rGen_m->uniform(0.0, 1.0);
    double phicou = z3 * 2 * pi;

    Vector_t Pcou(sin(thetacou)*cos(phicou), sin(thetacou)*sin(phicou), cos(thetacou));
    Vector_t Pcoutmp(cos(phicou), sin(phicou), 0.0);

    Vector_t tmpP = P;

    Rot(P, Pcou, normP);

    Rot(tmpP, Pcoutmp, normP);

    //    R = R + P * Physics::c / gamma * deltat + tmpP / beta / gamma * poscou;


    R = R + P * Physics::c * deltat /  sqrt(1.0  + dot(P, P)) + tmpP / beta / gamma * poscou / (Physics::c * deltat); // old



    if(collshape_m == "CCollimator") 
        R = R * 1000.0;

    double P2 = rGen_m->uniform(0.0, 1.0);
    if(P2 < 0.0047) {
        double P3 = rGen_m->uniform(0.0, 1.0);
        double thetaru = 2.5 * sqrt(1 / P3) * sqrt(2.0) * theta0;
        double P4 = rGen_m->uniform(0.0, 1.0);
        if(P4 > 0.5) 
            thetaru = -thetaru;
        double phiru = z3 * 2 * pi;
        Vector_t Pru(sin(thetaru)*cos(phiru), sin(thetaru)*sin(phiru), cos(thetaru));
        Rot(P, Pru, normP);
    }
}

void CollimatorPhysics::setCColimatorGeom() {

    double slope;
    if (xend_m == xstart_m)
      slope = 1.0e12;
    else
      slope = (yend_m - ystart_m) / (xend_m - xstart_m);

    double coeff2 = sqrt(1 + slope * slope);
    double coeff1 = slope / coeff2;
    double halfdist = width_m / 2.0;
    geom_m[0].x = xstart_m - halfdist * coeff1;
    geom_m[0].y = ystart_m + halfdist / coeff2;

    geom_m[1].x = xstart_m + halfdist * coeff1;
    geom_m[1].y = ystart_m - halfdist / coeff2;

    geom_m[2].x = xend_m + halfdist * coeff1;
    geom_m[2].y = yend_m - halfdist  / coeff2;

    geom_m[3].x = xend_m - halfdist * coeff1;
    geom_m[3].y = yend_m + halfdist / coeff2;

    geom_m[4].x = geom_m[0].x;
    geom_m[4].y = geom_m[0].y;

    if (zstart_m > zend_m){
      double tempz = 0.0;
      tempz = zstart_m;
      zstart_m = zend_m;
      zend_m = tempz;
    }
}


int CollimatorPhysics::checkPoint(const double &x, const double &y) {
    int    cn = 0;

    for(int i = 0; i < 4; i++) {
        if(((geom_m[i].y <= y) && (geom_m[i+1].y > y))
           || ((geom_m[i].y > y) && (geom_m[i+1].y <= y))) {

            float vt = (float)(y - geom_m[i].y) / (geom_m[i+1].y - geom_m[i].y);
            if(x < geom_m[i].x + vt * (geom_m[i+1].x - geom_m[i].x))
                ++cn;
        }
    }
    return (cn & 1);  // 0 if even (out), and 1 if odd (in)
}

bool  CollimatorPhysics::checkInColl(Vector_t R)
{            
  // this is for CCollimator FixMe: need to go to Collimator class 
  if(R(2) < zend_m && R(2) > zstart_m ) 
    return (checkPoint(R(0), R(1)) == 1 );
  else
    return false;
   
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
      This particle is back to the bunch, by setting the lable to -1.0
      the particle will be deleted.
    */
    locParts_m[i].label = -1.0;
}

void CollimatorPhysics::copyFromBunch(PartBunch &bunch)
{
    if ((min(bunch.Bin) < 0) && (bunch.getTotalNum() > 10)){   // we have particles entering mater 
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
                x.label        = 0; // allive in matter

                locParts_m.push_back(x);

                bunchToMatStat_m++;
            }
        }
    }    
}

void CollimatorPhysics::print(Inform &msg){
    Inform::FmtFlags_t ff = msg.flags();
    msg << std::scientific;

    /*
    if(dynamic_cast<Degrader *>(element_ref_m)) {
        Degrader *deg = dynamic_cast<Degrader *>(element_ref_m);
        width_m = deg->getZStart();
    }
    */
    // msg << "\n--- CollimatorPhysics - Type is " << collshape_m << " ------------------------------------------\n" << endl;
    //    msg << "StartElement= " << std::setw(8) << std::setprecision(3) << Begin_m  
    //	<< " (m) EndElement= " << std::setw(8) << std::setprecision(3) << Begin_m + width_m << endl;
    
    // msg << "Material " << material_m
    //   << " a= " << a_m << " (m) b= " << b_m << " (m)" << endl;
    //msg << "dTm= " << std::setw(8) << std::setprecision(3) << dT_m << " sub-timesteps " << N_m << endl;
    
    msg << "Coll/Deg statistics:  t= " << time_m << " (s) total particles in material " << locParts_m.size() 
        << " new hits " << bunchToMatStat_m << " redifused " << redifusedStat_m 
        << " stopped " << stoppedPartStat_m 
        << " Eavg= " << Eavg_m*1E3 << " (MeV)" << endl;
    
    // msg << "\n--- CollimatorPhysics -------------------------------------------------\n" << endl;
    msg.flags(ff);
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
