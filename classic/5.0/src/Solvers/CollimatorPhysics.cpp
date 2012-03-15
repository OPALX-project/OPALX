//Class:CollimatorPhysics
//  Defines the collimator physics models
// ------------------------------------------------------------------------
// Class category: 
// ------------------------------------------------------------------------
// $Date: 2009/07/20 09:32:31 $
// $Author: bi $
//-------------------------------------------------------------------------
#include "Solvers/CollimatorPhysics.hh"
#include "Physics/Physics.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/Collimator.h"
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/Multipole.h"
#include "Structure/LossDataSink.h"
#include <iostream>
#include <fstream>
using Physics::c;
using Physics::pi;
using Physics::m_p;
using Physics::m_e;
using Physics::h_bar;
using Physics::epsilon_0;
using Physics::r_e;
using Physics::z_p;
using Physics::Avo;


CollimatorPhysics::CollimatorPhysics(const string &name, ElementBase* element, const double& major, const double& minor,  string &material):
    SurfacePhysicsHandler(name, element),
    a_m(major),
    b_m(minor),
    x0_m(0.0),
    y0_m(0.0),
    Z_m(0),
    A_m(0.0),
    rho_m(0.0),
    X0_m(0.0),
    I_m(0.0),
    n_m(0.0),
    rGen_m(NULL),
    index_m(0),
    Begin_m(0.0),
    End_m(0.0),
    material_m(material)
{ 
    if (!rGen_m)
        rGen_m = new RANLIB_class(115314159,4);
}

CollimatorPhysics::~CollimatorPhysics()
{
    label_m.clear();
    Rincol_m.clear();
    Pincol_m.clear();
    IDincol_m.clear();
    Binincol_m.clear();
    DTincol_m.clear();
    Qincol_m.clear();
    LastSecincol_m.clear();
    Bfincol_m.clear();
    Efincol_m.clear();
    time_m.clear();
    steps_m.clear();
}


void CollimatorPhysics::apply(PartBunch &bunch)
{
    double scalefactor=Physics::c * bunch.getdT();
    double deltatbunch=bunch.getdT();
    double deltat=deltatbunch;
    double rr,rr1,rr2;

    if (End_m - Begin_m < 1e-8) {
        if (dynamic_cast<Collimator*>(element_ref_m)) {
            Collimator* coll = dynamic_cast<Collimator*>(element_ref_m);
            coll->getDimensions(Begin_m, End_m);
            FN_m=coll->getName();
            collshape_m=coll->getCollimatorShape();
	    x0_m=coll->getXpos();
	    y0_m=coll->getYpos();

        } else if (dynamic_cast<Drift*>(element_ref_m)) {
            Drift* drf = dynamic_cast<Drift*>(element_ref_m);
            drf->getDimensions(Begin_m, End_m);
            FN_m=drf->getName();
        } else if (dynamic_cast<SBend*>(element_ref_m)) {
            SBend* sben = dynamic_cast<SBend*>(element_ref_m);
            Begin_m=sben->getStartElement();
            End_m=Begin_m+sben->getElementLength();
            FN_m=sben->getName();
        } else if (dynamic_cast<RBend*>(element_ref_m)) {
            RBend* rben = dynamic_cast<RBend*>(element_ref_m);
            Begin_m=rben->getStartElement();
            End_m=Begin_m+rben->getElementLength();
            FN_m=rben->getName();
        } else if (dynamic_cast<Multipole*>(element_ref_m)) {
            Multipole* quad = dynamic_cast<Multipole*>(element_ref_m);
            quad->getDimensions(Begin_m, End_m);
            FN_m=quad->getName();
        }
    }

    ///the deltat in collimator is 1/n of the main tracking timestep.
    while (deltat>1.01e-12)
        deltat=deltat/10;
    int n=deltatbunch/deltat;
    /// check if the particle in partColArray is still within the collimator,
    /// if it is, track one step in collimator, if not, creat a new particle in the Bunch.
    for (int i=0;i<index_m;++i) {
        if (!label_m[i]){
            bool pdead = false; 
            Vector_t &R = Rincol_m[i];
            Vector_t &P = Pincol_m[i];
            if (collshape_m == "Slit"){


                if(a_m>0){
                    rr1=-R(0)*scalefactor/std::fabs(a_m);
                    rr2=R(0)*scalefactor/std::fabs(b_m);	
                    rr=0.0;
                    if(rr1>1.0||rr2>1.0) rr=2.0;
                }
                else{
                    rr1=-R(1)*scalefactor/std::fabs(a_m);
                    rr2=R(1)*scalefactor/std::fabs(b_m);	
                    rr=0.0;
                    if(rr1>1.0||rr2>1.0) rr=2.0;
                }
            }
            else if(collshape_m == "RCollimator"){
                rr1=std::fabs(R(0)*scalefactor/a_m);
                rr2=std::fabs(R(1)*scalefactor/b_m);
                rr=0.0;
                if(rr1>1.0||rr2>1.0) rr=2.0;

            }
	    else if(collshape_m == "Wire"){

	      rr1=std::fabs((R(0)*scalefactor-x0_m)/a_m);
	      //	      rr2=std::fabs((R(1)*scalefactor-y0_m)/b_m);
	      rr=0.0;
	      if (rr1<1.0){
		rr=2.0;
	      }

	    }
            else{
                rr = std::sqrt(std::pow(R(0)*scalefactor/a_m,2) + std::pow(R(1)*scalefactor/b_m,2));

            }
	
            if (R(2)*scalefactor>Begin_m&&R(2)*scalefactor<End_m&&rr>1.0) {
                steps_m[i]++;     
                double Eng=(sqrt(1.0  + dot(P,P))-1)*m_p;
                EnergyLoss(Eng, pdead,deltat);
                if (!pdead) {
                    double ptot=sqrt((m_p+Eng)*(m_p+Eng)-(m_p)*(m_p))/m_p;	
                    P(0)=P(0)*ptot/sqrt(dot(P,P));
                    P(1)=P(1)*ptot/sqrt(dot(P,P));
                    P(2)=P(2)*ptot/sqrt(dot(P,P));     
                    CoulombScat(R,P,deltat,scalefactor);	
                }
                else {
                    // dead 
                    //The lable of this particle is 1.0, so it is not in the partColArray any more
                    bunch.lossDs_m->addParticle(R*scalefactor,P,-IDincol_m[i]);  
                    label_m[i]=1.0;
                }
            }
            else{
                //create a new particle
                bunch.createWithID(IDincol_m[i]);
                bunch.Bin[bunch.getLocalNum()-1] =Binincol_m[i];
                bunch.R[bunch.getLocalNum()-1] = Rincol_m[i];
                bunch.P[bunch.getLocalNum()-1] = Pincol_m[i];
                bunch.Q[bunch.getLocalNum()-1] =  Qincol_m[i];
                bunch.LastSection[bunch.getLocalNum()-1] =LastSecincol_m[i];
                bunch.Bf[bunch.getLocalNum()-1] =Bfincol_m[i];
                bunch.Ef[bunch.getLocalNum()-1] = Efincol_m[i];    
                bunch.dt[bunch.getLocalNum()-1] = DTincol_m[i];    
                //remove this particle from the partColArray 	  
                label_m[i]=1.0;

            }
	
        }
      
    }    
    
    ///Check if the partilce in Bunch enters the Collimator, if it does, tracking one step in Collimator.
    Vector_t rmin,rmax;
    bunch.get_bounds(rmin,rmax);
    Vector_t rrms= bunch.get_rrms();
    Vector_t rmean=bunch.get_rmean();
     
    if (rmax(2)*scalefactor>Begin_m&&rmin(2)*scalefactor<End_m) {
        for (int i = 0; i < bunch.getLocalNum(); ++i) {   
            bool pdead = false;     
            Vector_t R = bunch.R[i];
            Vector_t P = bunch.P[i];
            if (collshape_m == "Slit"){
                if(a_m>0){
                    rr1=-R(0)*scalefactor/std::fabs(a_m);
                    rr2=R(0)*scalefactor/std::fabs(b_m);	
                    rr=0.0;
                    if(rr1>1.0||rr2>1.0) rr=2.0;
                }
                else{
                    rr1=-R(1)*scalefactor/std::fabs(a_m);
                    rr2=R(1)*scalefactor/std::fabs(b_m);	
                    rr=0.0;
                    if(rr1>1.0||rr2>1.0) rr=2.0;
                }
            }
            else if (collshape_m == "RCollimator"){
                rr1=std::fabs(R(0)*scalefactor/a_m);
                rr2=std::fabs(R(1)*scalefactor/b_m);
                rr=0.0;
                if(rr1>1.0||rr2>1.0) rr=2.0;
            }
      else if(collshape_m == "Wire"){

	rr1=std::fabs((R(0)*scalefactor-x0_m)/a_m);
	rr=0.0;
	if (rr1<1.0){
	  rr=2.0;
	}

      }
            else{
                rr = std::sqrt(std::pow(R(0)*scalefactor/a_m,2) + std::pow(R(1)*scalefactor/b_m,2));
            }
            if (R(2)*scalefactor>Begin_m&&R(2)*scalefactor<End_m&&rr>1.0) {
                //particle enters the collimator

                bunch.lossDs_m->addParticle(R*scalefactor,P,bunch.ID[i]);
	   
                double Eng=(sqrt(1.0  + dot(P,P))-1)*m_p;      
                EnergyLoss(Eng, pdead,deltat);

                if (!pdead) {
                    double ptot=sqrt((m_p+Eng)*(m_p+Eng)-(m_p)*(m_p))/m_p;
	             
                    P(0)=P(0)*ptot/sqrt(dot(P,P));
                    P(1)=P(1)*ptot/sqrt(dot(P,P));
                    P(2)=P(2)*ptot/sqrt(dot(P,P));
	    
                    CoulombScat(R,P,deltat,scalefactor);
                    //record the ID and R,P Bin ...  in partColArray
                    DTincol_m.push_back(bunch.dt[i]);
                    IDincol_m.push_back(bunch.ID[i]);
                    Binincol_m.push_back(bunch.Bin[i]);
                    Rincol_m.push_back(R);
                    Pincol_m.push_back(P);
                    Qincol_m.push_back(bunch.Q[i]);
                    LastSecincol_m.push_back(bunch.LastSection[i]);
                    Bfincol_m.push_back(bunch.Bf[i]);
                    Efincol_m.push_back(bunch.Ef[i]);
	    
                    steps_m.push_back(1);
                    label_m.push_back(0);
	    
                }
                else {
                    // dead 
                    long temp1=bunch.ID[i];
                    temp1=-temp1;
                    bunch.lossDs_m->addParticle(R*scalefactor,P,temp1);
	    
                }
                bunch.Bin[i] = -1;
 
            }
            if (std::fabs(R(2)-rmean(2))>7*rrms(2)&&rrms(2)>0)   bunch.Bin[i] = -1;
      
        }
    }


    ///index_m:The number of particle in partColArray.
    index_m= IDincol_m.size();

    ///Since the stepsize in Collimator is 1/n of main tracking timestep, loop over another n-1 step.
    for (int ii=0;ii<n-1;++ii){
        for (int i=0;i<index_m;++i) {
            if (!label_m[i]){
                bool pdead = false; 
                Vector_t &R = Rincol_m[i];
                Vector_t &P = Pincol_m[i];
                if (collshape_m == "Slit"){
                    if(a_m>0){
                        rr1=-R(0)*scalefactor/std::fabs(a_m);
                        rr2=R(0)*scalefactor/std::fabs(b_m);	
                        rr=0.0;
                        if(rr1>1.0||rr2>1.0) rr=2.0;
                    }
                    else{
                        rr1=-R(1)*scalefactor/std::fabs(a_m);
                        rr2=R(1)*scalefactor/std::fabs(b_m);	
                        rr=0.0;
                        if(rr1>1.0||rr2>1.0) rr=2.0;
                    }
                }
                else if (collshape_m == "RCollimator"){
                    rr1=std::fabs(R(0)*scalefactor/a_m);
                    rr2=std::fabs(R(1)*scalefactor/b_m);
                    rr=0.0;
                    if(rr1>1.0||rr2>1.0) rr=2.0;

                }
      else if(collshape_m == "Wire"){

	rr1=std::fabs((R(0)*scalefactor-x0_m)/a_m);
	rr=0.0;
	if (rr1<1.0){
	  rr=2.0;
	}

      }
                else{
                    rr = std::sqrt(std::pow(R(0)*scalefactor/a_m,2) + std::pow(R(1)*scalefactor/b_m,2));
                }
                double Eng=(sqrt(1.0  + dot(P,P))-1)*m_p;
                if (R(2)*scalefactor>Begin_m&&R(2)*scalefactor<End_m&&rr>1.0) {
                    steps_m[i]++;     

                    EnergyLoss(Eng, pdead,deltat);
                    if (!pdead) {
                        double ptot=sqrt((m_p+Eng)*(m_p+Eng)-(m_p)*(m_p))/m_p;	
                        P(0)=P(0)*ptot/sqrt(dot(P,P));
                        P(1)=P(1)*ptot/sqrt(dot(P,P));
                        P(2)=P(2)*ptot/sqrt(dot(P,P));     
                        CoulombScat(R,P,deltat,scalefactor);	
                    }
                    else {
                        //  dead 
                        //The lable of this particle is 1.0, so it is not in the partColArray any more.
                        bunch.lossDs_m->addParticle(R*scalefactor,P,-IDincol_m[i]);  
                        label_m[i]=1.0;
                    }
                }
                else{
                    ///The particle exits the Collimator, but still in the loop of the substep, 
                    ///so before it goes back to the Bunch, tracking it as in a drift.
                    double gamma=(Eng+m_p)/m_p;
                    double beta=sqrt(1.0-1.0/(gamma*gamma));
                    R(0)=R(0)+deltat*beta*c*P(0)/sqrt(dot(P,P))/scalefactor;
                    R(1)=R(1)+deltat*beta*c*P(1)/sqrt(dot(P,P))/scalefactor;
                    R(2)=R(2)+deltat*beta*c*P(2)/sqrt(dot(P,P))/scalefactor;

	  
                }
	
            }
      
        }
    }

    bunch.lossDs_m->save(FN_m);  


}
   

const string CollimatorPhysics::getType() const
{
    return "CollimatorPhysics";
}

/// The material of the collimator
//  ------------------------------------------------------------------------
void  CollimatorPhysics::Material()
{


    if (material_m == "Cu")
        {
            Z_m=29;
            A_m=63.546;
            rho_m=8.96;

            X0_m=12.86/rho_m/100;
            I_m=10*Z_m;
            n_m=rho_m/A_m*Avo;
        }

    if (material_m == "Be")
        {
            Z_m=4;
            A_m=9.012;
            rho_m=1.848;

            X0_m=65.19/rho_m/100;
            I_m=10*Z_m;
            n_m=rho_m/A_m*Avo;
        }

    if (material_m == "Graphite")
        {
            Z_m=6;
            A_m=12.0107;
            rho_m=2.210;

            X0_m=42.70/rho_m/100;
            I_m=10*Z_m;
            n_m=rho_m/A_m*Avo;
        }

    if (material_m == "Mo")
        {
            Z_m=42;
            A_m=95.94;
            rho_m=10.22;

            X0_m=9.8/rho_m/100;
            I_m=10*Z_m;
            n_m=rho_m/A_m*Avo;
        }

}

/// Energy Loss:  using the Bethe-Bloch equation.
/// Energy straggling: For relatively thick absorbers such that the number of collisions is large,
/// the energy loss distribution is shown to be Gaussian in form.
// -------------------------------------------------------------------------
void  CollimatorPhysics::EnergyLoss(double &Eng,bool &pdead, double &deltat)
{
 
    Material();
    double gamma=(Eng+m_p)/m_p;
    double beta=sqrt(1.0-1.0/(gamma*gamma));
    double gamma2=gamma*gamma;
    double beta2=beta*beta;
  
    double deltas=deltat*beta*c;
    double deltasrho=deltas*100*rho_m;
  
    double K=4.0*pi*Avo*r_e*r_e*m_e*1E7;
    double Tmax=2.0*m_e*1e9*beta2*gamma2/(1.0+2.0*gamma*m_e/m_p+(m_e/m_p)*(m_e/m_p));
    double dEdx=-K*z_p*z_p*Z_m/(A_m*beta2)*(1.0/2.0*log(2*m_e*1e9*beta2*gamma2*Tmax/I_m/I_m)-beta2);
    double delta_Eave=deltasrho*dEdx;
    double sigma_E=sqrt(K*m_e*rho_m*Z_m/A_m*deltas*1e5);
    double delta_E=rGen_m->gauss(delta_Eave,sigma_E);

    Eng=Eng+delta_E/1e3;
    if (Eng<1e-4)  pdead=true;

}
///From beam coordinate to lab coordinate
// --------------------------------------------------------------------------
double  CollimatorPhysics::Rot(double &p1, double &p2, double &scatang)
{
    double rotang;
    if (p1 >= 0&&p2>=0) rotang=atan(p1/p2);
    else if (p1>0&&p2<0) rotang=atan(p1/p2)+pi;
    else if (p1<0&&p2>0) rotang=atan(p1/p2)+2.0*pi;
    else rotang=atan(p1/p2)+pi;
    double p12=sqrt(p1*p1+p2*p2);
    double p1par=p12*sin(scatang);
    double p2par=p12*cos(scatang);
    p1=p1par*cos(rotang)+p2par*sin(rotang);
    p2=-p1par*sin(rotang)+p2par*cos(rotang);
    return rotang;
}

void  CollimatorPhysics::CoulombScat()
{}
void  CollimatorPhysics::EnergyLoss(double &Eng, bool &pdead)
{}

///Coulomb Scattering: Including Multiple Coulomb Scattering and large angle Rutherford Scattering.
///Using the distribution given in Classical Electrodynamics, by J. D. Jackson.
//--------------------------------------------------------------------------
void  CollimatorPhysics::CoulombScat(Vector_t &R,Vector_t &P, double &deltat,double scalefactor)
{


    Material();

    double Eng=sqrt(dot(P,P)+1)*m_p-m_p;
    double gamma=(Eng+m_p)/m_p;
    double beta=sqrt(1.0-1.0/(gamma*gamma));

    double deltas=deltat*beta*c;
  
    double theta0=13.6e6/(beta*sqrt( dot(P,P) )*m_p*1e9 )*z_p*sqrt(deltas/X0_m)*(1.0+0.038*log(deltas/X0_m));

    //x direction------------------------------
    double z1=rGen_m->gauss(0.0,1.0);
    double z2=rGen_m->gauss(0.0,1.0);
    double thetaplanex=z2*theta0;
    while (fabs(thetaplanex)>3.5*theta0)
        {
            z1=rGen_m->gauss(0.0,1.0);
            z2=rGen_m->gauss(0.0,1.0);
            thetaplanex=z2*theta0;
        }
    double xplane=z1*deltas*theta0/sqrt(12.0)+z2*deltas*theta0/2.0;
    double rotmcx=Rot(P(0), P(2), thetaplanex);

    R(0)=R(0)+P(0)/sqrt( dot(P,P) )*deltas/scalefactor+xplane*cos(rotmcx)/scalefactor;
    R(2)=R(2)-xplane*sin(rotmcx)/scalefactor;


    double P2=rGen_m->uniform(0.0,1.0);
    if (P2<0.0047)
        {
            double P3=rGen_m->uniform(0.0,1.0);
            double thetarux=2.5*sqrt(1/P3)*sqrt(2.0)*theta0;  
            double P4=rGen_m->uniform(0.0,1.0);
            if (P4>0.5) thetarux=-thetarux;
            double rotrux=Rot(P(0), P(2), thetarux);
        }

    //y direction------------------------------
    double z3=rGen_m->gauss(0.0,1.0);
    double z4=rGen_m->gauss(0.0,1.0); 
    double thetaplaney=z4*theta0;
    while (fabs(thetaplaney)>3.5*theta0)
        {
            z3=rGen_m->gauss(0.0,1.0);
            z4=rGen_m->gauss(0.0,1.0);
            thetaplaney=z4*theta0;
        }
    double yplane=z3*deltas*theta0/sqrt(12.0)+z4*deltas*theta0/2.0;
    double rotmcy=Rot(P(1), P(2), thetaplaney);

    R(1)=R(1)+P(1)/sqrt( dot(P,P) )*deltas/scalefactor+yplane*cos(rotmcy)/scalefactor;
    R(2)=R(2)+P(2)/sqrt( dot(P,P) )*deltas/scalefactor-yplane*sin(rotmcy)/scalefactor;


    double P5=rGen_m->uniform(0.0,1.0);
    if (P5<0.0047)
        {
            double P6=rGen_m->uniform(0.0,1.0);
            double thetaruy=2.5*sqrt(1/P6)*sqrt(2.0)*theta0;  
            double P7=rGen_m->uniform(0.0,1.0);
            if (P7>0.5) thetaruy=-thetaruy;
            double rotruy=Rot(P(1), P(2), thetaruy);
        }

}


