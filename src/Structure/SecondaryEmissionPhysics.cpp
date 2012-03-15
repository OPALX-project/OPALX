//#include <fstream>
#include <vector>
//#include <iostream>
#include <cmath>
#include <algorithm>
#include "hdf5.h"
#include "Distribution/Distribution.h"
#include <sys/stat.h>
#include "Ippl.h"
#include "Structure/SecondaryEmissionPhysics.h"
#include "BasicActions/Option.h"
typedef Vektor<double, 3> Vector_t;

using namespace myeps;
using namespace Physics;

SecondaryEmissionPhysics::SecondaryEmissionPhysics() {

    TPnSec_m = IpplTimings::getTimer("Secondary emission");
}
/**
 * Destructor.
 *
 * Delete the previous defined member arrays
 */
SecondaryEmissionPhysics::~SecondaryEmissionPhysics() {

}

void SecondaryEmissionPhysics::nSec(const double &incEnergy,  const double &cosTheta, const int &matNumber, int &seNum, int &seType, const double &incQ, const Vector_t &TriNorm, const Vector_t &inteCoords, const Vector_t &localX, PartBunch *itsBunch, double &seyNum, const double &ppVw, const double &vVThermal, const bool nEmissionMode) {

    IpplTimings::startTimer(TPnSec_m);
    double prob[11] = {0};
    std::vector<Vector_t> se_P;
    setSeMaterial(matNumber);//set material based parameters
    seyNum=calcProb(incEnergy, cosTheta, prob);//calculate probability
    calcEmiNum(incEnergy, cosTheta, prob, seNum);//calculate emitted number
    PAssert(seNum<11);//debug
    PAssert(seNum>=2);//debug
    double Eemit[10];
    double emiTheta[10];
    double emiPhi[10];
    Vector_t interCoords_l = inteCoords;
    Vector_t TriNorm_l = TriNorm;
    double incQ_l = incQ;
   
    
    /*===========================Definations for benchmark===================================*/
    double vw=ppVw; //1.6*1e-19*1200/9.10938188*1e-31/(2*3.1415926*2.0*1e8)/0.03;//benchmark
    double vt=vVThermal;//7.268929821*1e5;//1.5eV//benchmark
    double f_max=vw/vt*exp(-0.5);//benchmark
    double test_a=vt/vw;//benchmark
    double test_asq=test_a*test_a;//benchmark
    /*---------------------------------------------------------------------------------------*/
    if( Options::ppdebug ) {

    } else {
           
        if(seNum != 0) {
            
            for(int i = 0; i < seNum; i++) {
                
                double tmp1 = IpplRandom();
                double tmp2 = IpplRandom();
                double temp = 1.0 / (1.0 + seAlpha_m);
                emiTheta[i] = acos(pow(tmp1, temp));
                emiPhi[i] = Physics::two_pi * tmp2;
                
            }
        }
    }

  
    if(seNum == 0) {
        
        // The incident particle will be marked for deletion
	if (!nEmissionMode) {
	    if( Options::ppdebug ) {
		/*=========(Velocity with Maxwellian Distribution For Parallel Plate Benchmark)========*/
		double test_s=1;
		double f_x=0;
		double test_x=0;
		while (test_s>f_x) {
		    test_s=IpplRandom();
		    test_s*=f_max;
		    test_x=IpplRandom();
		    test_x*=10*test_a;// range for normalized emission speed(0,10*test_a);
		    f_x=test_x/test_asq*exp(-test_x*test_x/2/test_asq);
		}
		double v_emi=test_x*vw;
		Eemit[0]=(1.0/sqrt(1.0-v_emi*v_emi/9.0/1e16)-1)*Physics::m_e*1.0e9;
		// cout<<"Single Eemit[0]: "<<Eemit[0]<<endl;
		/*---------------------End Of Maxwellian Distribution(For Benchmark)--------------------*/
	    } else {

		// For absorption case in constant simulation particle mode, just use true secondary emission with 1 particle energy distribution for Furman Pivi model.
		double u2 = IpplRandom();
                
                double temp = incEnergy / seEpsn_m[0];
                double p0 = gammp(sePn_m[0], temp);
                temp = p0 * u2;
                Eemit[0] = seEpsn_m[0] * invgammp(temp, sePn_m[0]) ;
	

	    }

	}


    } else if(seNum == 1) {

        if( Options::ppdebug ) {
            /*=========(Velocity with Maxwellian Distribution For Parallel Plate Benchmark)========*/
            double test_s=1;
            double f_x=0;
            double test_x=0;
            while (test_s>f_x) {
                test_s=IpplRandom();
                test_s*=f_max;
                test_x=IpplRandom();
                test_x*=10*test_a;//range for normalized emission speed(0,10*test_a);
                f_x=test_x/test_asq*exp(-test_x*test_x/2/test_asq);
            }
            double v_emi=test_x*vw;
            Eemit[0]=(1.0/sqrt(1.0-v_emi*v_emi/9.0/1e16)-1)*Physics::m_e*1.0e9;
          //cout<<"Single Eemit[0]: "<<Eemit[0]<<endl;
            /*---------------------End Of Maxwellian Distribution(For Benchmark)--------------------*/
        } else {
            
            //double delta_e = calcDeltae(incEnergy, cosTheta);
            //double delta_r = calcDeltar(incEnergy, cosTheta);
            double tmp = prob[1] + deltae_m + deltar_m; 
            //double ae = delta_e / tmp;
	    double ae = deltae_m / tmp;
            //double ar = delta_r / tmp;
	    double ar = deltar_m / tmp;
            double a_re = ae + ar;
            double urand = IpplRandom();
            
            
            if(urand < ae) {
                int t_count = 0;
                do {
                    Eemit[0] = incEnergy - seDelta_m * fabs(gaussRand()) ;
                    t_count++;
                } while(Eemit[0] < 0&&t_count<200);
                if(Eemit[0]<0)// if the above do - while loops over 200 times, the loop will break out, and Eemit will be its mean value, i.e., incident energy.
                    Eemit[0]=incEnergy;
                seType = 0;
                
            } else if(urand >= ae && urand < a_re) {
                double u1 = IpplRandom();
                
                double powArg = 1.0 / (1.0 + seQ_m);
                Eemit[0] = incEnergy * pow(u1, powArg);
                seType = 1;
                
                
            } else {
                
                double u2 = IpplRandom();
                
                double temp = incEnergy / seEpsn_m[0];
                double p0 = gammp(sePn_m[0], temp);
                temp = p0 * u2;
                Eemit[0] = seEpsn_m[0] * invgammp(temp, sePn_m[0]) ;
                seType = 2;
             
            }
        }
       
    } else {
        seType = 2;
       
        if( Options::ppdebug ) {
            /*==========(Velocity with Maxwellian Distribution For Parallel Plate Benchmark)========*/
	    if (!nEmissionMode) {
		/*double Eemit_mean = 0.0;
		  for(int i = 0; i < seNum; i++) {
		  Eemit_mean += Eemit[i];
		  }*/
		//Eemit[0] = Eemit_mean/seNum;
		double test_s=1.0;
		double f_x=0;
		double test_x=0;
		while (test_s>f_x) {
		    test_s=IpplRandom();
		    test_s*=f_max;
		    test_x=IpplRandom();
		    test_x*=10*test_a;//range for normalized emission speed(0,10*test_a);
		    f_x=test_x/test_asq*exp(-test_x*test_x/2/test_asq);
		}
		double v_emi=test_x*vw;
		Eemit[0]=(1.0/sqrt(1.0-v_emi*v_emi/9.0/1e16)-1)*Physics::m_e*1.0e9;
		
		
	    } else {
		for(int i = 0; i < seNum; i++) {
		    double test_s=1.0;
		    double f_x=0;
		    double test_x=0;
		    while (test_s>f_x) {
			test_s=IpplRandom();
			test_s*=f_max;
			test_x=IpplRandom();
			test_x*=10*test_a;//range for normalized emission speed(0,10*test_a);
			f_x=test_x/test_asq*exp(-test_x*test_x/2/test_asq);
		    }
		    double v_emi=test_x*vw;
		    Eemit[i]=(1.0/sqrt(1.0-v_emi*v_emi/9.0/1e16)-1)*Physics::m_e*1.0e9;
		    
		}
	    
	
	    }
            /*---------------------End Of Maxwellian Distribution(For Benchmark)--------------------*/
        } else {
            /*=======================3D velocity for Furman-Pivi's Model====================*/
            double x0 = incEnergy / seEpsn_m[seNum-1];
            double parg = seNum * sePn_m[seNum-1];
            double p0 =  gammp(parg, x0);
            double sin2theta_n[seNum];
            double cos2theta_n[seNum];
            double rand_y =  IpplRandom();
            double invarg = rand_y * p0;
            double y2 = invgammp(invarg, parg);
            double y2_n[seNum];
            double multisin = 1.0;
            double Eemisum = 0.0;
            
	    if (!nEmissionMode) {// only emit 1 particle

		for(int i = 0; i < seNum - 1; i++) {
		    double mu = sePn_m[seNum-1] * (seNum  - i);
		    double nu =  sePn_m[seNum-1];
		    double rand_n = IpplRandom();
		    sin2theta_n[i] = invbetai(rand_n, mu, nu);
		    cos2theta_n[i] = 1.0 - sin2theta_n[i];
		    if(i != 0) {
			
			multisin *=  sin2theta_n[i-1];
		    }
		    y2_n[i] = y2 *  multisin * cos2theta_n[i];
		    Eemit[i] = seEpsn_m[seNum-1] * y2_n[i];
		    Eemisum += Eemit[i];
		    
		    
		}
            
		Eemit[seNum-1] = Eemit[seNum-2] / cos2theta_n[seNum-2] * sin2theta_n[seNum-2];
		Eemisum += Eemit[seNum-1];
	   
		Eemit[0] = Eemisum/seNum;
	    
	    } else {// emit seNum particles
		for(int i = 0; i < seNum - 1; i++) {
		    double mu = sePn_m[seNum-1] * (seNum  - i);
		    double nu =  sePn_m[seNum-1];
		    double rand_n = IpplRandom();
		    sin2theta_n[i] = invbetai(rand_n, mu, nu);
		    cos2theta_n[i] = 1.0 - sin2theta_n[i];
		    if(i != 0) {
			
			multisin *=  sin2theta_n[i-1];
		    }
		    y2_n[i] = y2 *  multisin * cos2theta_n[i];
		    Eemit[i] = seEpsn_m[seNum-1] * y2_n[i];
		    		    
		    
		}
		
		Eemit[seNum-1] = Eemit[seNum-2] / cos2theta_n[seNum-2] * sin2theta_n[seNum-2];
		
        /*=====================================================================================================*/
	    }
	}
               
    }
        
    Vector_t z_unit = TriNorm_l;
    Vector_t x_unit = localX - interCoords_l;
    double tmp = sqrt(dot(x_unit,x_unit));
    x_unit /= tmp;
    Vector_t y_unit = cross(z_unit,x_unit);

    size_t lowMark = itsBunch->getLocalNum();
    if (!nEmissionMode) {// emit only 1 particle with larger charge instead of emit seNum secondaries.
	//if (seNum>0) {// we dont delete particles even when seNum==0.
	double gamma_const = Eemit[0] / Physics::m_e/1.0e9 + 1.0;
	double beta_const = sqrt(1.0 - 1.0 / pow(gamma_const, 2.0));
	double P_emitted = gamma_const * beta_const;
	
	Vector_t P_local;
	Vector_t P_global = (0.0);
	
	if( Options::ppdebug ) {
	    
	    P_global = P_emitted*TriNorm_l;// 1D for parallel plate benchmark
	    
	} else {
	    /*==================3D for Furman-Pivi's Model====================*/
	    P_local[2] = P_emitted * cos(emiTheta[0]);
	    P_local[1] = P_emitted * sin(emiTheta[0]) * sin(emiPhi[0]);
	    P_local[0] = P_emitted * sin(emiTheta[0]) * cos(emiPhi[0]);
	    P_global = P_local[0]*x_unit + P_local[1]*y_unit +  P_local[2]*z_unit;//Pivi's model
	    /*================================================================*/
	    
	    
	}
	itsBunch->create(1);
	itsBunch->R[lowMark] = interCoords_l;
	    
	itsBunch->P[lowMark] = P_global;
	itsBunch->Bin[lowMark] = 0;
	itsBunch->PType[lowMark] = 3;// 3 to denote the newly generated secondaries
	itsBunch->TriID[lowMark] = 0;
	//itsBunch->Q[lowMark] = incQ_l*seNum;// charge of simulation particle will be sum of secondaies
	itsBunch->Q[lowMark] = incQ_l*seyNum;// charge of simulation particle will be multiplied by SEY.
	itsBunch->LastSection[lowMark] = 0;
	itsBunch->Ef[lowMark] = Vector_t(0.0);
	itsBunch->Bf[lowMark] = Vector_t(0.0);
	itsBunch->dt[lowMark] = itsBunch->getdT();
	    //}
	
	
	
    } else {
	for(size_t i = 0; i < (size_t) seNum; i++) {
	    
	    double gamma_const = Eemit[i] / Physics::m_e/1.0e9 + 1.0;
	    double beta_const = sqrt(1.0 - 1.0 / pow(gamma_const, 2.0));
	    double P_emitted = gamma_const * beta_const;
	    
	    Vector_t P_local;
	    Vector_t P_global = (0.0);
	    
	    if( Options::ppdebug ) {
		
		P_global = P_emitted*TriNorm_l;// 1D for parallel plate benchmark
		
	    } else {
		/*==================3D for Furman-Pivi's Model====================*/
		P_local[2] = P_emitted * cos(emiTheta[i]);
		P_local[1] = P_emitted * sin(emiTheta[i]) * sin(emiPhi[i]);
		P_local[0] = P_emitted * sin(emiTheta[i]) * cos(emiPhi[i]);
		P_global = P_local[0]*x_unit + P_local[1]*y_unit +  P_local[2]*z_unit;//Pivi's model
		/*================================================================*/
           
		
	    }
	    itsBunch->create(1);
	    itsBunch->R[lowMark+i] = interCoords_l;
	    
	    itsBunch->P[lowMark+i] = P_global;
	    itsBunch->Bin[lowMark+i] = 0;
	    itsBunch->PType[lowMark+i] = 3;//3 to denote the newly generated secondaries
	    itsBunch->TriID[lowMark+i] = 0;
	    itsBunch->Q[lowMark+i] = incQ_l;
	    itsBunch->LastSection[lowMark+i] = 0;
	    itsBunch->Ef[lowMark+i] = Vector_t(0.0);
	    itsBunch->Bf[lowMark+i] = Vector_t(0.0);
	    itsBunch->dt[lowMark+i] = itsBunch->getdT();
	    
	}
    }
        
    IpplTimings::stopTimer(TPnSec_m);
}


//Vaughan's secondary emission model.
void SecondaryEmissionPhysics::nSec(const double &incEnergy,  const double &cosTheta, int &seNum, int &seType, const double &incQ, const Vector_t &TriNorm, const Vector_t &inteCoords, const Vector_t &localX, PartBunch *itsBunch, double &seyNum, const double &ppVw,const double &vSeyZero, const double &vEzero, const double &vSeyMax, const double &vEmax, const double &vKenergy, const double &vKtheta, const double &vVThermal, const bool nEmissionMode) {


    IpplTimings::startTimer(TPnSec_m);
       
    std::vector<Vector_t> se_P;
    calcEmiNum(incEnergy, cosTheta, seNum, vSeyZero, vEzero, vSeyMax, vEmax, vKenergy, vKtheta, seyNum);//calculate emitted number & SEY factor
    double Eemit[seNum];
    double emiTheta[seNum];
    double emiPhi[seNum];
    Vector_t interCoords_l = inteCoords;
    Vector_t TriNorm_l = TriNorm;//simpler model
    double vw=ppVw; //1.6*1e-19*1200/9.10938188*1e-31/(2*3.1415926*2.0*1e8)/0.03;//benchmark
    double vt=vVThermal;//7.268929821*1e5 1.5eV//benchmark
    double f_max;//benchmark
    double test_a;//benchmark
    double test_asq;//benchmark
    if( Options::ppdebug ) {
	f_max=vw/vt*exp(-0.5);// velocity a Maxwellian distribution. See Anza et al.,Phys. Plasmas 17, 062110 (2010)

	test_a=vt/vw;
	test_asq=test_a*test_a;
    } else {// Energy Maxwell-Boltzmann distribution f(E)=E/a^2*exp(-E/a), a= kinetic energy w.r.t thermal speed vt. See www.nuc.berkeley.edu/dept/Courses/NE-255/minicourse.pdf p.20
	test_a = Physics::m_e*(1.0/sqrt(1-vt*vt/Physics::c/Physics::c)-1.0)*1.0e9;// m_e GeV change it to eV
	test_asq=test_a*test_a;
	f_max= 1.0/test_a*exp(-1.0);

    }
    
    double incQ_l = incQ;
    if( Options::ppdebug ) {
        // 1D emission angle along the surface triangle normal

    } else {
	if (!nEmissionMode) {

	    double tmp1 = IpplRandom();
	    double tmp2 = IpplRandom();
	    double seAlpha = 1.0;
	    double temp = 1.0 / (1.0 + seAlpha);// pow(cosine(theta),seAlpha) distribution. Here seAlpha=1.0
	    emiTheta[0] = acos(pow(tmp1, temp));
	    emiPhi[0] = Physics::two_pi * tmp2;


	} else {

	    for(int i = 0; i < seNum; i++) {
		    
		double tmp1 = IpplRandom();
		double tmp2 = IpplRandom();
		double seAlpha = 1.0;
		double temp = 1.0 / (1.0 + seAlpha);// pow(cosine(theta),seAlpha) distribution. Here seAlpha=1.0
		emiTheta[i] = acos(pow(tmp1, temp));
		emiPhi[i] = Physics::two_pi * tmp2;
		
	    }
	    
	}
    }
    
    if(seNum == 0) {
	if (!nEmissionMode) {
	    if( Options::ppdebug ) {
		double test_s=1.0;
		double f_x=0;
		double test_x=0;
		while (test_s>f_x) {
		    test_s=IpplRandom();
		    test_s*=f_max;
		    test_x=IpplRandom();
		    test_x*=10*test_a;// truncation range for emission speed(0,10*test_a);
		    f_x=test_x/test_asq*exp(-test_x*test_x/2/test_asq);
		}
		double v_emi=test_x*vw;
		Eemit[0]=(1.0/sqrt(1.0-v_emi*v_emi/Physics::c/Physics::c)-1)*Physics::m_e*1.0e9;
            } else {
		double test_s=1.0;
		double f_x=0;
		double test_x=0;
		while (test_s>f_x) {// rejection & acceptance method
		    test_s=IpplRandom();
		    test_s*=f_max;
		    test_x=IpplRandom();
		    test_x*=10*test_a;// truncation range for emission energy(0,10*test_a);
		    f_x=test_x/test_asq*exp(-test_x/test_a);// thermal distribution for energy
		}
	
		Eemit[0]=test_x;


	    }

	}
        // else The incident particle will be marked for deletion


    } else {
        seType = 2;
        if (!nEmissionMode) {
	    if( Options::ppdebug ) {
		double test_s=1.0;
		double f_x=0;
		double test_x=0;
		while (test_s>f_x) {
		    test_s=IpplRandom();
		    test_s*=f_max;
		    test_x=IpplRandom();
		    test_x*=10*test_a;// truncation range for emission speed(0,10*test_a);
		    f_x=test_x/test_asq*exp(-test_x*test_x/2/test_asq);
		}
		double v_emi=test_x*vw;
		Eemit[0]=(1.0/sqrt(1.0-v_emi*v_emi/Physics::c/Physics::c)-1)*Physics::m_e*1.0e9;
            } else {
		double test_s=1.0;
		double f_x=0;
		double test_x=0;
		while (test_s>f_x) {
		    test_s=IpplRandom();
		    test_s*=f_max;
		    test_x=IpplRandom();
		    test_x*=10*test_a;// truncation range for emission energy(0,10*test_a);
		    f_x=test_x/test_asq*exp(-test_x/test_a);// thermal distribution for energy
		}
	
		Eemit[0]=test_x;

	    }

	} else {
	    if( Options::ppdebug ) {
		/*=======================Maxwellian Distribution====================*/
		// the velocity distribution for Vaughan's model:fu=u*vw*vw/vt/vt*exp(-u*u*vw*vw/vt/vt) valid for benchmarking; 
		
		for(int i = 0; i < seNum; i++) {
		    double test_s=1.0;
		    double f_x=0;
		    double test_x=0;
		    while (test_s>f_x) {
			test_s=IpplRandom();
			test_s*=f_max;
			test_x=IpplRandom();
			test_x*=10*test_a;//range for normalized emission speed(0,10*test_a);
			f_x=test_x/test_asq*exp(-test_x*test_x/2/test_asq);
		    }
		    double v_emi=test_x*vw;
		    Eemit[i]=(1.0/sqrt(1.0-v_emi*v_emi/9.0/1e16)-1)*Physics::m_e*1.0e9;
		    
		}
		
		/*---------------------End Of Maxwellian Distribution--------------------*/
	    } else {// energy thermal distribution for Vaughan model
		for(int i = 0; i < seNum; i++) {
		    double test_s=1.0;
		    double f_x=0;
		    double test_x=0;
		    while (test_s>f_x) {
			test_s=IpplRandom();
			test_s*=f_max;
			test_x=IpplRandom();
			test_x*=10*test_a;// truncation range for emission energy(0,10*test_a);
			f_x=test_x/test_asq*exp(-test_x/test_a);// thermal distribution for energy
		    }
		    
		    Eemit[i]=test_x;
		}
	    }
	}
    }
        
    Vector_t z_unit = TriNorm;
    Vector_t x_unit = localX - interCoords_l;
    double tmp = sqrt(dot(x_unit,x_unit));
    x_unit /= tmp;
    Vector_t y_unit = cross(z_unit,x_unit);

    size_t lowMark = itsBunch->getLocalNum();
    if (!nEmissionMode) {
	double gamma_const = Eemit[0] / Physics::m_e/1.0e9 + 1.0;
	double beta_const = sqrt(1.0 - 1.0 / pow(gamma_const, 2.0));
	double P_emitted = gamma_const * beta_const;
	
	
	Vector_t P_local;
	Vector_t P_global = (0.0);
	if( Options::ppdebug ) {
	    
	    P_global = P_emitted*TriNorm_l;//1D for parallel plate benchmark
		
	} else {
	    /*==================3D Vaughan' Model==========================================*/
	    P_local[2] = P_emitted * cos(emiTheta[0]);
	    P_local[1] = P_emitted * sin(emiTheta[0]) * sin(emiPhi[0]);
	    P_local[0] = P_emitted * sin(emiTheta[0]) * cos(emiPhi[0]);
	    P_global = P_local[0]*x_unit + P_local[1]*y_unit +  P_local[2]*z_unit;// the same cosine distribution as Pivi's model
	    /*=============================================================================*/
	}
	    
	itsBunch->create(1);
	itsBunch->R[lowMark] = interCoords_l;
	
	itsBunch->P[lowMark] = P_global;
	itsBunch->Bin[lowMark] = 0;
	itsBunch->PType[lowMark] = 3;//3 to denote the newly generated secondaries
	itsBunch->TriID[lowMark] = 0;
	itsBunch->Q[lowMark] = incQ_l*seyNum;
	itsBunch->LastSection[lowMark] = 0;// fixme: what about last section !=0 ? 
	itsBunch->Ef[lowMark] = Vector_t(0.0);
	itsBunch->Bf[lowMark] = Vector_t(0.0);
	itsBunch->dt[lowMark] = itsBunch->getdT();

    } else {
	for(size_t i = 0; i < (size_t) seNum; i++) {
	    
	    double gamma_const = Eemit[i] / Physics::m_e/1.0e9 + 1.0;
	    double beta_const = sqrt(1.0 - 1.0 / pow(gamma_const, 2.0));
	    double P_emitted = gamma_const * beta_const;
	    
	    
	    Vector_t P_local;
	    Vector_t P_global = (0.0);
	    if( Options::ppdebug ) {
		
		P_global = P_emitted*TriNorm_l;//1D for parallel plate benchmark
		
	    } else {
		/*==================3D Vaughan' Model==========================================*/
		P_local[2] = P_emitted * cos(emiTheta[i]);
		P_local[1] = P_emitted * sin(emiTheta[i]) * sin(emiPhi[i]);
		P_local[0] = P_emitted * sin(emiTheta[i]) * cos(emiPhi[i]);
		P_global = P_local[0]*x_unit + P_local[1]*y_unit +  P_local[2]*z_unit;//the same cosine distribution as Pivi's model
		/*=============================================================================*/
	    }
	    
	    itsBunch->create(1);
	    itsBunch->R[lowMark+i] = interCoords_l;
	    
	    itsBunch->P[lowMark+i] = P_global;
	    itsBunch->Bin[lowMark+i] = 0;
	    itsBunch->PType[lowMark+i] = 3;//3 to denote the newly generated secondaries
	    itsBunch->TriID[lowMark+i] = 0;
	    itsBunch->Q[lowMark+i] = incQ_l;
	    itsBunch->LastSection[lowMark+i] = 0;
	    itsBunch->Ef[lowMark+i] = Vector_t(0.0);
	    itsBunch->Bf[lowMark+i] = Vector_t(0.0);
	    itsBunch->dt[lowMark+i] = itsBunch->getdT();
	    
            
	}
    }
            
    IpplTimings::stopTimer(TPnSec_m);

}
