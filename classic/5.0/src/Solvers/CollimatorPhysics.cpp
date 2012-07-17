//Class:CollimatorPhysics
//  Defines the collimator physics models
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
#include "AbsBeamline/Drift.h"
#include "AbsBeamline/SBend.h"
#include "AbsBeamline/RBend.h"
#include "AbsBeamline/Multipole.h"
#include "Structure/LossDataSink.h"
#include "Distribution/ranlib.h"

#include <iostream>
#include <fstream>
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
    incoll_m(false),
    index_m(0) {
    rGen_m = new RANLIB_class(115314159, 4);
    // initialize DataSink with H5Part output enabled
    bool doH5 = false;
    lossDs_m = new LossDataSink(1000000, doH5);
    lossDs_m->openH5(FN_m);
    if(dynamic_cast<Collimator *>(element_ref_m)) {
        Collimator *coll = dynamic_cast<Collimator *>(element_ref_m);
        coll->getDimensions(Begin_m, End_m);
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
    }

}

CollimatorPhysics::~CollimatorPhysics() {
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
    if(lossDs_m)
        delete lossDs_m;
}

void CollimatorPhysics::apply(PartBunch &bunch) {
    double scalefactor;
    // int flagcoll; !! carefull not initialized !!
    // Vector_t tmploss;
    double deltatbunch = bunch.getdT();
    double deltat = deltatbunch;
    double rr = 0.0, rr1 = 0.0, rr2 = 0.0;

    if(collshape_m == "CCollimator") {
        scalefactor = 1;
    } else {
        scalefactor = Physics::c * bunch.getdT();
    }

    ///the deltat in collimator is 1/n of the main tracking timestep.
    while(deltat > 1.01e-12)
        deltat = deltat / 10;
    int n = deltatbunch / deltat;
    /// check if the particle in partColArray is still within the collimator,
    /// if it is, track one step in collimator, if not, creat a new particle in the Bunch.
    for(int i = 0; i < index_m; ++i) {

        if(!label_m[i]) {
            bool pdead = false;
            Vector_t &R = Rincol_m[i];
            Vector_t &P = Pincol_m[i];

            if(collshape_m == "CCollimator") {
                if(R(2) < zend_m && R(2) > zstart_m ) {
                    if(checkPoint(R(0), R(1)) == 1 )
                        incoll_m = true;
                }
                else
                    incoll_m = false;
            } else {

                if(collshape_m == "Slit") {
                    if(a_m > 0) {
                        rr1 = -R(0) * scalefactor / std::fabs(a_m);
                        rr2 = R(0) * scalefactor / std::fabs(b_m);
                        rr = 0.0;
                        if(rr1 > 1.0 || rr2 > 1.0) rr = 2.0;

                    } else {
                        rr1 = -R(1) * scalefactor / std::fabs(a_m);
                        rr2 = R(1) * scalefactor / std::fabs(b_m);
                        rr = 0.0;
                        if(rr1 > 1.0 || rr2 > 1.0) rr = 2.0;

                    }
                } else if(collshape_m == "RCollimator") {
                    rr1 = std::fabs(R(0) * scalefactor / a_m);
                    rr2 = std::fabs(R(1) * scalefactor / b_m);
                    rr = 0.0;
                    if(rr1 > 1.0 || rr2 > 1.0) rr = 2.0;

                } else if(collshape_m == "Wire") {

                    rr1 = std::fabs((R(0) * scalefactor - xp_m) / a_m);
                    //        rr2=std::fabs((R(1)*scalefactor-yp_m)/b_m);
                    rr = 0.0;
                    if(rr1 < 1.0) {
                        rr = 2.0;
                    }
                } else if(collshape_m == "ECollimator") {
                    rr = std::sqrt(std::pow(R(0) * scalefactor / a_m, 2) + std::pow(R(1) * scalefactor / b_m, 2));
                }
                if(R(2)*scalefactor > Begin_m && R(2)*scalefactor < End_m && rr > 1.0) incoll_m = true;
            }


            if(incoll_m) {
                incoll_m = false;
                steps_m[i]++;
                double Eng = (sqrt(1.0  + dot(P, P)) - 1) * m_p;

                EnergyLoss(Eng, pdead, deltat);

                if(!pdead) {
                    double ptot = sqrt((m_p + Eng) * (m_p + Eng) - (m_p) * (m_p)) / m_p;
                    P(0) = P(0) * ptot / sqrt(dot(P, P));
                    P(1) = P(1) * ptot / sqrt(dot(P, P));
                    P(2) = P(2) * ptot / sqrt(dot(P, P));

                    CoulombScat(R, P, deltat, scalefactor);

                } else {
                    // dead
                    //The lable of this particle is 1.0, so it is not in the partColArray any more
                    // tmploss = Vector_t(P(0), P(1), flagcoll);
                    //      lossDs_m->addParticle(R*scalefactor,tmploss,-IDincol_m[i]);
                    label_m[i] = 1.0;
                }
            } else {

                //create a new particle
                bunch.createWithID(IDincol_m[i]);
                bunch.Bin[bunch.getLocalNum()-1] = Binincol_m[i];
                bunch.R[bunch.getLocalNum()-1] = Rincol_m[i];
                bunch.P[bunch.getLocalNum()-1] = Pincol_m[i];
                bunch.Q[bunch.getLocalNum()-1] =  Qincol_m[i];
                bunch.LastSection[bunch.getLocalNum()-1] = LastSecincol_m[i];
                bunch.Bf[bunch.getLocalNum()-1] = Bfincol_m[i];
                bunch.Ef[bunch.getLocalNum()-1] = Efincol_m[i];
                bunch.dt[bunch.getLocalNum()-1] = DTincol_m[i];
                //remove this particle from the partColArray
                label_m[i] = 1.0;

            }

        }

    }

    ///Check if the partilce in Bunch enters the Collimator, if it does, tracking one step in Collimator.
    Vector_t rmin, rmax;
    bunch.get_bounds(rmin, rmax);
    Vector_t rrms = bunch.get_rrms();
    Vector_t rmean = bunch.get_rmean();

    if(rmax(2)*scalefactor > Begin_m && rmin(2)*scalefactor < End_m) {
        for(unsigned int i = 0; i < bunch.getLocalNum(); ++i) {
            bool pdead = false;
            Vector_t R = bunch.R[i];
            Vector_t P = bunch.P[i];

            if(collshape_m == "CCollimator") {
                if(R(2) < zend_m && R(2) > zstart_m ) {
                    if(checkPoint(R(0), R(1)) == 1 )
                        incoll_m = true;
                }
                else
                    incoll_m = false;
            } else {
 
                if(collshape_m == "Slit") {
                    if(a_m > 0) {
                        rr1 = -R(0) * scalefactor / std::fabs(a_m);
                        rr2 = R(0) * scalefactor / std::fabs(b_m);
                        rr = 0.0;
                        if(rr1 > 1.0 || rr2 > 1.0) rr = 2.0;
                    } else {
                        rr1 = -R(1) * scalefactor / std::fabs(a_m);
                        rr2 = R(1) * scalefactor / std::fabs(b_m);
                        rr = 0.0;
                        if(rr1 > 1.0 || rr2 > 1.0) rr = 2.0;
                    }
                } else if(collshape_m == "RCollimator") {
                    rr1 = std::fabs(R(0) * scalefactor / a_m);
                    rr2 = std::fabs(R(1) * scalefactor / b_m);
                    rr = 0.0;
                    if(rr1 > 1.0 || rr2 > 1.0) rr = 2.0;
                } else if(collshape_m == "Wire") {

                    rr1 = std::fabs((R(0) * scalefactor - xp_m) / a_m);
                    rr = 0.0;
                    if(rr1 < 1.0) {
                        rr = 2.0;
                    }

                } else if(collshape_m == "ECollimator") {
                    rr = std::sqrt(std::pow(R(0) * scalefactor / a_m, 2) + std::pow(R(1) * scalefactor / b_m, 2));
                }
                if(R(2)*scalefactor > Begin_m && R(2)*scalefactor < End_m && rr > 1.0) incoll_m = true;
            }
            if(incoll_m) {

                //particle enters the collimator

                incoll_m = false;
                // tmploss = Vector_t(P(0), P(1), flagcoll);

                lossDs_m->addParticle(R * scalefactor, P, bunch.ID[i]);

                double Eng = (sqrt(1.0  + dot(P, P)) - 1) * m_p;

                EnergyLoss(Eng, pdead, deltat);

                if(!pdead) {
                    double ptot = sqrt((m_p + Eng) * (m_p + Eng) - (m_p) * (m_p)) / m_p;

                    P(0) = P(0) * ptot / sqrt(dot(P, P));
                    P(1) = P(1) * ptot / sqrt(dot(P, P));
                    P(2) = P(2) * ptot / sqrt(dot(P, P));

                    CoulombScat(R, P, deltat, scalefactor);

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

                } else {
                    // dead
                    long temp1 = bunch.ID[i];
                    temp1 = -temp1;
                    // tmploss = Vector_t(P(0), P(1), flagcoll);
                    //   lossDs_m->addParticle(R*scalefactor,tmploss,temp1);

                }
                bunch.Bin[i] = -1;

            }
            if(collshape_m != "CCollimator") {
                if(std::fabs(R(2) - rmean(2)) > 7 * rrms(2) && rrms(2) > 0)   bunch.Bin[i] = -1;
            }

        }
    }


    ///index_m:The number of particle in partColArray.
    index_m = IDincol_m.size();

    ///Since the stepsize in Collimator is 1/n of main tracking timestep, loop over another n-1 step.
    for(int ii = 0; ii < n - 1; ++ii) {
        for(int i = 0; i < index_m; ++i) {
            if(!label_m[i]) {
                bool pdead = false;
                Vector_t &R = Rincol_m[i];
                Vector_t &P = Pincol_m[i];

                if(collshape_m == "CCollimator") {
                    if(R(2) < zend_m && R(2) > zstart_m ) {
                        if(checkPoint(R(0), R(1)) == 1 )
                            incoll_m = true;
                    }
                    else
                        incoll_m = false;
                } else {
                  
                    if(collshape_m == "Slit") {
                        if(a_m > 0) {
                            rr1 = -R(0) * scalefactor / std::fabs(a_m);
                            rr2 = R(0) * scalefactor / std::fabs(b_m);
                            rr = 0.0;
                            if(rr1 > 1.0 || rr2 > 1.0) rr = 2.0;
                        } else {
                            rr1 = -R(1) * scalefactor / std::fabs(a_m);
                            rr2 = R(1) * scalefactor / std::fabs(b_m);
                            rr = 0.0;
                            if(rr1 > 1.0 || rr2 > 1.0) rr = 2.0;
                        }
                    } else if(collshape_m == "RCollimator") {
                        rr1 = std::fabs(R(0) * scalefactor / a_m);
                        rr2 = std::fabs(R(1) * scalefactor / b_m);
                        rr = 0.0;
                        if(rr1 > 1.0 || rr2 > 1.0) rr = 2.0;

                    } else if(collshape_m == "Wire") {

                        rr1 = std::fabs((R(0) * scalefactor - xp_m) / a_m);
                        rr = 0.0;
                        if(rr1 < 1.0) {
                            rr = 2.0;
                        }

                    } else if(collshape_m == "ECollimator") {
                        rr = std::sqrt(std::pow(R(0) * scalefactor / a_m, 2) + std::pow(R(1) * scalefactor / b_m, 2));
                    }

                    if(R(2)*scalefactor > Begin_m && R(2)*scalefactor < End_m && rr > 1.0) incoll_m = true;
                }
                double Eng = (sqrt(1.0  + dot(P, P)) - 1) * m_p;
                if(incoll_m) {
                    incoll_m = false;
                    steps_m[i]++;

                    EnergyLoss(Eng, pdead, deltat);

                    if(!pdead) {
                        double ptot = sqrt((m_p + Eng) * (m_p + Eng) - (m_p) * (m_p)) / m_p;
                        P(0) = P(0) * ptot / sqrt(dot(P, P));
                        P(1) = P(1) * ptot / sqrt(dot(P, P));
                        P(2) = P(2) * ptot / sqrt(dot(P, P));

                        CoulombScat(R, P, deltat, scalefactor);

                    } else {
                        //  dead
                        //The lable of this particle is 1.0, so it is not in the partColArray any more.
                        // tmploss = Vector_t(P(0), P(1), flagcoll);
                        //    lossDs_m->addParticle(R*scalefactor,tmploss,-IDincol_m[i]);
                        label_m[i] = 1.0;
                    }
                } else {
                    ///The particle exits the Collimator, but still in the loop of the substep,
                    ///so before it goes back to the Bunch, tracking it as in a drift.
                    double gamma = (Eng + m_p) / m_p;
                    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
                    if(collshape_m == "CCollimator") {
                        R(0) = R(0) + deltat * beta * Physics::c * P(0) / sqrt(dot(P, P)) / scalefactor * 1000;
                        R(1) = R(1) + deltat * beta * Physics::c * P(1) / sqrt(dot(P, P)) / scalefactor * 1000;
                        R(2) = R(2) + deltat * beta * Physics::c * P(2) / sqrt(dot(P, P)) / scalefactor * 1000;
                    } else {
                        R(0) = R(0) + deltat * beta * Physics::c * P(0) / sqrt(dot(P, P)) / scalefactor;
                        R(1) = R(1) + deltat * beta * Physics::c * P(1) / sqrt(dot(P, P)) / scalefactor;
                        R(2) = R(2) + deltat * beta * Physics::c * P(2) / sqrt(dot(P, P)) / scalefactor;
                    }

                }

            }

        }
    }
    lossDs_m->save(FN_m);

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
        I_m = 10 * Z_m;
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

    Material();
    double gamma = (Eng + m_p) / m_p;
    double beta = sqrt(1.0 - 1.0 / (gamma * gamma));
    double gamma2 = gamma * gamma;
    double beta2 = beta * beta;

    double deltas = deltat * beta * Physics::c;
    double deltasrho = deltas * 100 * rho_m; 
    double K = 4.0 * pi * Avo * r_e * r_e * m_e * 1E7;
    double sigma_E = sqrt(K * m_e * rho_m * (Z_m/A_m))* deltas * 1E5;
    if (Eng > 0.00001&& Eng < 0.0006) //GeV
    { 
	double Ts = (Eng*1E6)/1.0073; // the 1.0073 comes from the proton mass divided by the atomic mass number. T is in KeV
	double epsilon_low = A2_c*pow(Ts,0.45);		
	double epsilon_high = (A3_c/Ts)*log(1+(A4_c/Ts)+(A5_c*Ts));
	double epsilon = (epsilon_low*epsilon_high)/(epsilon_low + epsilon_high);
	double stopping_power = - epsilon /(1E21*(A_m/Avo)); // Stopping_power is in MeV
	INFOMSG("stopping power: " << stopping_power << " MeV" <<endl);
	double delta_Eave = (deltasrho * stopping_power);
	double delta_E1 = rGen_m->gauss(delta_Eave, sigma_E);
	double delta_E = delta_E1/1000.; // Delta_E is in GeV
	double Eng = Eng + delta_E;
    }
    
    if (Eng >= 0.0006)
    {
	double Tmax = 2.0 * m_e * 1e9 * beta2 * gamma2 / (1.0 + 2.0 * gamma * m_e / m_p + (m_e / m_p) * (m_e / m_p));
	double dEdx = -K * z_p * z_p * Z_m / (A_m * beta2) * (1.0 / 2.0 * log(2 * m_e * 1e9 * beta2 * gamma2 * Tmax / I_m / I_m) - beta2);
	INFOMSG("stopping power_BB: " << dEdx << " MeV");
	double delta_Eave = deltasrho * dEdx;
	double delta_E = rGen_m->gauss(delta_Eave, sigma_E);
    Eng = Eng+delta_E / 1E3;
    }	
    else
      INFOMSG("final energy: " << Eng/1000 << " MeV" <<endl);
	 pdead = true;

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
void  CollimatorPhysics::CoulombScat(Vector_t &R, Vector_t &P, double &deltat, double scalefactor) {

    Material();
    double Eng = sqrt(dot(P, P) + 1) * m_p - m_p;
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
    if(collshape_m == "CCollimator") {
        R = R + (P * Physics::c / gamma * deltat / scalefactor + tmpP / beta / gamma * poscou / scalefactor) * 1000;
    } else {
        R = R + P * Physics::c / gamma * deltat / scalefactor + tmpP / beta / gamma * poscou / scalefactor;
    }
    double P2 = rGen_m->uniform(0.0, 1.0);
    if(P2 < 0.0047) {
        double P3 = rGen_m->uniform(0.0, 1.0);
        double thetaru = 2.5 * sqrt(1 / P3) * sqrt(2.0) * theta0;
        double P4 = rGen_m->uniform(0.0, 1.0);
        if(P4 > 0.5) thetaru = -thetaru;
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
