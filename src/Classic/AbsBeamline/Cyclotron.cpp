// ------------------------------------------------------------------------
// $RCSfile: Cyclotron.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: Cyclotron
//   Defines the abstract interface for a sector bend magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2007/08/01 $
// $Author: Yang, Adelmann $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Cyclotron.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Physics/Physics.h"
#include "Structure/LossDataSink.h"
#include "Utilities/Options.h"
#include "Fields/Fieldmap.h"
#include "Utilities/GeneralClassicException.h"
#include <fstream>

#define CHECK_CYC_FSCANF_EOF(arg) if(arg == EOF)\
throw GeneralClassicException("Cyclotron::getFieldFromFile",\
                              "fscanf returned EOF at " #arg);

extern Inform *gmsg;

using Physics::pi;
using namespace std;

// Class Cyclotron
// ------------------------------------------------------------------------

Cyclotron::Cyclotron():
    Component() {
}


Cyclotron::Cyclotron(const Cyclotron &right):
    Component(right),
    fmapfn_m(right.fmapfn_m),
    rffrequ_m(right.rffrequ_m),
    rfphi_m(right.rfphi_m),
    escale_m(right.escale_m),
    superpose_m(right.superpose_m),
    symmetry_m(right.symmetry_m),
    rinit_m(right.rinit_m),
    prinit_m(right.prinit_m),
    phiinit_m(right.phiinit_m),
    zinit_m(right.zinit_m),
    pzinit_m(right.pzinit_m),
    spiral_flag_m(right.spiral_flag_m),
    type_m(right.type_m),
    harm_m(right.harm_m),
    bscale_m(right.bscale_m),
    tcr1_m(right.tcr1_m),
    tcr2_m(right.tcr2_m),
    mbtc_m(right.mbtc_m),
    slptc_m(right.slptc_m),
    tcr1V_m(right.tcr1V_m),
    tcr2V_m(right.tcr2V_m),
    mbtcV_m(right.mbtcV_m),
    slptcV_m(right.slptcV_m),
    minr_m(right.minr_m),
    maxr_m(right.maxr_m),
    minz_m(right.minz_m),
    maxz_m(right.maxz_m),
    RFfilename_m(right.RFfilename_m),
    RFFCoeff_fn_m(right.RFFCoeff_fn_m),
    RFVCoeff_fn_m(right.RFVCoeff_fn_m),
    fmLowE_m(right.fmLowE_m),
    fmHighE_m(right.fmHighE_m) {
}


Cyclotron::Cyclotron(const std::string &name):
    Component(name) {
}


Cyclotron::~Cyclotron() {
    //~ if(BP.rarr) delete[] BP.rarr;
}


void Cyclotron::applyTrimCoil(double r, double z, double slope, double tcr1, double tcr2, double magnitude, double *br, double *bz) {
  /// updated bz and br with trim coil contributions
  if(tcr1 != 0.0 && tcr2 != 0.0) {
    const double Amax1 = 1;
    const double Amax2 = 3;
    const double Amin = -2;
    const double x01 = 4;
    const double x02 = 8;
    const double h1 = 0.03;
    const double h2 = 0.2;
    double part1;
    double part2;

    if (r < ((tcr2+tcr1)/2)) {
      part1 = pow(10.0, (r / slope - tcr1 / slope - x01) * h1);
      part2 = pow(10.0, (x02 - r / slope + tcr1 / slope) * h2);
    }
    else {
      part1 = pow(10.0, (tcr2 / slope - r / slope - x01) * h1);
      part2 = pow(10.0, (x02 - tcr2 / slope + r / slope) * h2);
    }
    double part3 = -(Amax1 - Amin) * h1 * log(10) / slope / (1 + part1) / (1 + part1) * part1;
    double part4 = (Amax2 - Amin) * h2 * log(10) / slope / (1 + part2) / (1 + part2) * part2;
    double dr = magnitude / 2.78 * (part3 + part4);
    double btr = magnitude / 2.78 * (Amin + (Amax1 - Amin) / (1 + part1) + (Amax2 - Amin) / (1 + part2) - 1.0);

    if(r < tcr2)       {
      *bz -= btr;
      *br -= dr * z;
    }
  }
}


void Cyclotron::accept(BeamlineVisitor &visitor) const {
    visitor.visitCyclotron(*this);
}

void Cyclotron::setRinit(double rinit) {
    rinit_m = rinit;
}

double Cyclotron::getRinit() const {
    return rinit_m;
}

void Cyclotron::setPRinit(double prinit) {
    prinit_m = prinit;
}

double Cyclotron::getPRinit() const {
    return prinit_m;
}

void Cyclotron::setPHIinit(double phiinit) {
    phiinit_m = phiinit;
}

double Cyclotron::getPHIinit() const {
    return phiinit_m;
}

void Cyclotron::setZinit(double zinit){
    zinit_m = zinit;
}

double Cyclotron::getZinit() const {
    return zinit_m;
}

void Cyclotron::setPZinit(double pzinit){
    pzinit_m = pzinit;
}

double Cyclotron::getPZinit() const {
    return pzinit_m;
}

void Cyclotron::setSpiralFlag(bool spiral_flag) {
    spiral_flag_m = spiral_flag;
}

bool Cyclotron::getSpiralFlag() const {
    return spiral_flag_m;
}

void Cyclotron::setFieldMapFN(std::string f) {
    fmapfn_m = f;
}

string Cyclotron::getFieldMapFN() const {
    return fmapfn_m;
}

void Cyclotron::setRfFieldMapFN(vector<string> f) {
    RFfilename_m = f;
}

void Cyclotron::setRFFCoeffFN(vector<string> f) {
    RFFCoeff_fn_m = f;
}

void Cyclotron::setRFVCoeffFN(vector<string> f) {
    RFVCoeff_fn_m = f;
}

void Cyclotron::setRfPhi(vector<double> f) {
    rfphi_m = f;
}

void Cyclotron::setRfFrequ(vector<double> f) {
    rffrequ_m = f;
}

double Cyclotron::getRfFrequ() const {
  return rffrequ_m[0];
}

void Cyclotron::setSuperpose(std::vector<bool> flag) {
  superpose_m = flag;
}

//bool Cyclotron::getSuperpose() const {
//    return superpose_m;
//}

void Cyclotron::setSymmetry(double s) {
    symmetry_m = s;
}

double Cyclotron::getSymmetry() const {
    return symmetry_m;
}


void Cyclotron::setType(std::string t) {
    type_m = t;
}

const std::string &Cyclotron::getCyclotronType() const {
    return type_m;
}

ElementBase::ElementType Cyclotron::getType() const {
    return CYCLOTRON;
}

void Cyclotron::setCyclHarm(double h) {
    harm_m = h;
}

void Cyclotron::setBScale(double s) {
    bscale_m = s;
}

double Cyclotron::getBScale() const {
    return bscale_m;
}

void Cyclotron::setEScale(vector<double> s) {
    escale_m = s;
}

unsigned int Cyclotron::getNumberOfTrimcoils() const {
  return tcr1V_m.size();
}

double Cyclotron::getCyclHarm() const {
    return harm_m;
}

double Cyclotron::getRmin() const {
    return BP.rmin;
}


double Cyclotron::getRmax() const {
    return BP.rmin + (Bfield.nrad - 1) * BP.delr;
}

// This function aims at obtaining magentic field at any given location R by interpolation.
// arguments t is useless here.
// arguments E is set to zero.


void Cyclotron::setTCr1(double tcr1) {
    tcr1_m = tcr1;
}

double Cyclotron::getTCr1() const {
    return tcr1_m;
}

void Cyclotron::setTCr2(double tcr2) {
    tcr2_m = tcr2;
}

double Cyclotron::getTCr2() const {
    return tcr2_m;
}
void Cyclotron::setMBtc(double mbtc) {
    mbtc_m = mbtc;
}

double Cyclotron::getMBtc() const {
    return mbtc_m;
}

void Cyclotron::setSLPtc(double slptc) {
    slptc_m = slptc;
}

double Cyclotron::getSLPtc() const {
    return slptc_m;
}


void Cyclotron::setMinR(double r) {
    minr_m = r;
}

void Cyclotron::setMaxR(double r) {
    maxr_m = r;
}
double Cyclotron::getMinR() const {
    return minr_m;
}

double Cyclotron::getMaxR() const {
    return maxr_m;
}

void  Cyclotron::setMinZ(double z) {
    minz_m = z;
}
double Cyclotron::getMinZ() const {
    return minz_m;
}
void Cyclotron::setMaxZ(double z) {
    maxz_m = z;
}
double Cyclotron::getMaxZ() const {
    return maxz_m;
}

void Cyclotron::setTCr1V(vector<double> tcr1) {
    tcr1V_m = tcr1;
}

void Cyclotron::setTCr2V(vector<double> tcr2) {
    tcr2V_m = tcr2;
}

void Cyclotron::setMBtcV(vector<double> mbtc) {
    mbtcV_m = mbtc;
}

void Cyclotron::setSLPtcV(vector<double> slptc) {
    slptcV_m = slptc;
}

void Cyclotron::setFMLowE(double e) { fmLowE_m = e;}
double Cyclotron::getFMLowE() const { return fmLowE_m;}

void Cyclotron::setFMHighE(double e) { fmHighE_m = e;}
double Cyclotron::getFMHighE() const { return fmHighE_m;}


bool Cyclotron::apply(const size_t &id, const double &t, Vector_t &E, Vector_t &B) {

  bool flagNeedUpdate = false;

  const double rpos = sqrt(RefPartBunch_m->R[id](0) * RefPartBunch_m->R[id](0)
                           + RefPartBunch_m->R[id](1) * RefPartBunch_m->R[id](1));
  const double zpos = RefPartBunch_m->R[id](2);

  if (zpos > maxz_m || zpos < minz_m || rpos > maxr_m || rpos < minr_m){
      flagNeedUpdate = true;
      Inform gmsgALL("OPAL ", INFORM_ALL_NODES);
      gmsgALL << getName() << ": particle "<< id <<" out of the global aperture of cyclotron!"<< endl;

  } else{

      flagNeedUpdate = apply(RefPartBunch_m->R[id], RefPartBunch_m->P[id], t, E, B);
      if(flagNeedUpdate){
          Inform gmsgALL("OPAL ", INFORM_ALL_NODES);
          gmsgALL << getName() << ": particle "<< id <<" out of the field map boundary!"<< endl;
      }
  }

  if (flagNeedUpdate) {
      lossDs_m->addParticle(RefPartBunch_m->R[id], RefPartBunch_m->P[id],id);
      RefPartBunch_m->Bin[id] = -1;
  }

  return flagNeedUpdate;
}

bool Cyclotron::apply(const Vector_t &R, const Vector_t &P, const double &t, Vector_t &E, Vector_t &B) {

    const double rad = sqrt(R[0] * R[0] + R[1] * R[1]);
    const double xir = (rad - BP.rmin) / BP.delr;

    // ir : the mumber of path whoes radius is less then the 4 points of cell which surrond particle.
    const int    ir = (int)xir;

    // wr1 : the relative distance to the inner path radius
    const double wr1 = xir - (double)ir;
    // wr2 : the relative distance to the outer path radius
    const double wr2 = 1.0 - wr1;

    const double tempv = atan(R[1] / R[0]);
    double tet = tempv, tet_map, xit;

    /* if((R[0] > 0) && (R[1] >= 0)) tet = tempv;
       else*/
    if((R[0] < 0) && (R[1] >= 0)) tet = pi + tempv;
    else if((R[0] < 0) && (R[1] <= 0)) tet = pi + tempv;
    else if((R[0] > 0) && (R[1] <= 0)) tet = 2.0 * pi + tempv;
    else if((R[0] == 0) && (R[1] > 0)) tet = pi / 2.0;
    else if((R[0] == 0) && (R[1] < 0)) tet = 1.5 * pi;

    double tet_rad = tet;

    // the actual angle of particle
    tet = tet / pi * 180.0;

    // Necessary for gap phase output -DW
    if (0 <= tet && tet <= 45) waiting_for_gap = 1;

    // the corresponding angle on the field map
    // Note: this does not work if the start point of field map does not equal zero.
    tet_map = fmod(tet, 360.0 / symmetry_m);

    xit = tet_map / BP.dtet;

    int it = (int) xit;

    //    *gmsg << R << " tet_map= " << tet_map << " ir= " << ir << " it= " << it << " bf= " << Bfield.bfld[idx(ir,it)] << endl;

    const double wt1 = xit - (double)it;
    const double wt2 = 1.0 - wt1;

    // it : the number of point on the inner path whoes angle is less then the particle' corresponding angle.
    // include zero degree point
    it = it + 1;

    int r1t1, r2t1, r1t2, r2t2;
    int ntetS = Bfield.ntet + 1;

    // r1t1 : the index of the "min angle, min radius" point in the 2D field array.
    // considering  the array start with index of zero, minus 1.

    if(myBFieldType_m != FFAGBF) {
        /*
          For FFAG this does not work
        */
        r1t1 = it + ntetS * ir - 1;
        r1t2 = r1t1 + 1;
        r2t1 = r1t1 + ntetS;
        r2t2 = r2t1 + 1 ;

    } else {
        /*
          With this we have B-field AND this is far more
          intuitive for me ....
        */
        r1t1 = idx(ir, it);
        r2t1 = idx(ir + 1, it);
        r1t2 = idx(ir, it + 1);
        r2t2 = idx(ir + 1, it + 1);
    }

    double bzf = 0.0, bz = 0.0 /*, bzcub = 0.0*/;
    double brf = 0.0, br = 0.0 /*, brcub = 0.0*/;
    double btf = 0.0, bt = 0.0 /*, btcub = 0.0*/;

    if((it >= 0) && (ir >= 0) && (it < Bfield.ntetS) && (ir < Bfield.nrad)) {

        /* Bz */
        bzf = (Bfield.bfld[r1t1] * wr2 * wt2 + Bfield.bfld[r2t1] * wr1 * wt2 +
               Bfield.bfld[r1t2] * wr2 * wt1 + Bfield.bfld[r2t2] * wr1 * wt1);

        // bzcub = (Bfield.f2[r1t1] * wr2 * wt2 +
        //          Bfield.f2[r2t1] * wr1 * wt2 +
        //          Bfield.f2[r1t2] * wr2 * wt1 +
        //          Bfield.f2[r2t2] * wr1 * wt1) * pow(R[2], 2.0);

        // bz = -( bzf - bzcub );
        bz = - bzf ;


        /* Br */
        brf = (Bfield.dbr[r1t1] * wr2 * wt2 +
               Bfield.dbr[r2t1] * wr1 * wt2 +
               Bfield.dbr[r1t2] * wr2 * wt1 +
               Bfield.dbr[r2t2] * wr1 * wt1) * R[2];


        // brcub = (Bfield.f3[r1t1] * wr2 * wt2 +
        //          Bfield.f3[r2t1] * wr1 * wt2 +
        //          Bfield.f3[r1t2] * wr2 * wt1 +
        //          Bfield.f3[r2t2] * wr1 * wt1) * pow(R[2], 3.0);

        // br = -( brf - brcub );
        br = - brf;


        /* Btheta */
        btf = (Bfield.dbt[r1t1] * wr2 * wt2 +
               Bfield.dbt[r2t1] * wr1 * wt2 +
               Bfield.dbt[r1t2] * wr2 * wt1 +
               Bfield.dbt[r2t2] * wr1 * wt1) / rad * R[2];


        // btcub = (Bfield.g3[r1t1] * wr2 * wt2 +
        //          Bfield.g3[r2t1] * wr1 * wt2 +
        //          Bfield.g3[r1t2] * wr2 * wt1 +
        //          Bfield.g3[r2t2] * wr1 * wt1) / rad * pow(R[2], 3.0);

        // bt = -( btf - btcub );
        bt = - btf;


        /* Br Btheta -> Bx By */

	if (slptcV_m.size() != 0) {
	  for (unsigned int i=0; i<slptcV_m.size(); i++)
	    applyTrimCoil(rad, R[2], slptcV_m[i], tcr1V_m[i], tcr2V_m[i], mbtcV_m[i], &br, &bz);
	}
	else
	  applyTrimCoil(rad, R[2], slptc_m, tcr1_m, tcr2_m, mbtc_m, &br, &bz);

        B[0] = br * cos(tet_rad) - bt * sin(tet_rad);
        B[1] = br * sin(tet_rad) + bt * cos(tet_rad);
        B[2] = bz;

	//*gmsg << "R = " << rad << ", Theta = " << tet << ", B = (" << B[0] << "/" << B[1] << "/" << B[2] << ")" << endl;

    } else {
      return true;
    }
    if(myBFieldType_m == BANDRF) {
      //The RF field is suppose to be sampled on a cartesian grid
        vector<Fieldmap *>::const_iterator fi  = RFfields_m.begin();
        vector<double>::const_iterator rffi    = rffrequ_m.begin();
        vector<double>::const_iterator rfphii  = rfphi_m.begin();
        vector<double>::const_iterator escali  = escale_m.begin();
	vector<bool>::const_iterator superposei = superpose_m.begin();
        double xBegin(0), xEnd(0), yBegin(0), yEnd(0), zBegin(0), zEnd(0);
        int fcount = 0;
        for(; fi != RFfields_m.end(); ++fi, ++rffi, ++rfphii, ++escali, ++superposei) {
            (*fi)->getFieldDimensions(xBegin, xEnd, yBegin, yEnd, zBegin, zEnd);
	    bool SuperPose = *superposei;
            if (fcount > 0 && !SuperPose) {
	      //INFOMSG ("Field maps taken : " << fcount << "Superpose false" << endl);
	      break;
            }
            if (R(0) >= xBegin && R(0) <= xEnd && R(1) >= yBegin && R(1) <= yEnd && R(2) >= zBegin && R(2) <= zEnd) {
                Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);
                if(!(*fi)->getFieldstrength(R, tmpE, tmpB)) {
		  ++fcount;
                  double phase = 2.0 * pi * (1E-3 * (*rffi)) * t + *rfphii;

		  double ebscale = *escali;
                  E += ebscale * cos(phase) * tmpE;
                  B -= ebscale * sin(phase) * tmpB;

//                INFOMSG("Field " << fcount << " BANDRF E= " << tmpE << " R= " << R << " phase " << phase << endl);
                }
            }
	}
    } else if(myBFieldType_m == SYNCHRO) {
        //The RF field is suppose to be sampled on a cartesian grid
        vector<Fieldmap *>::const_iterator fi  = RFfields_m.begin();
        vector<double>::const_iterator rffi    = rffrequ_m.begin();
        vector< vector<double> >::const_iterator rffci = rffc_m.begin();
        vector< vector<double> >::const_iterator rfvci = rfvc_m.begin();
        vector<double>::const_iterator rfphii  = rfphi_m.begin();
        vector<double>::const_iterator escali  = escale_m.begin();
	vector<bool>::const_iterator superposei = superpose_m.begin();
        double xBegin(0), xEnd(0), yBegin(0), yEnd(0), zBegin(0), zEnd(0);
        int fcount = 0;
	int ecount = 0;
	double frequency, ebscale;

        for(; fi != RFfields_m.end(); ++fi, ++rffi, ++rfphii, ++escali, ++superposei, ++rffci, ++rfvci) {
            (*fi)->getFieldDimensions(xBegin, xEnd, yBegin, yEnd, zBegin, zEnd);
	    bool SuperPose = *superposei;
            if (fcount > 0 && !SuperPose) {
	      //INFOMSG ("Field maps taken : " << fcount << "Superpose false" << endl);
	      break;
            }
            if (R(0) >= xBegin && R(0) <= xEnd && R(1) >= yBegin && R(1) <= yEnd && R(2) >= zBegin && R(2) <= zEnd) {
                Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);
                if(!(*fi)->getFieldstrength(R, tmpE, tmpB)) {

		    ++fcount;

		    frequency = (*rffi);   // frequency in MHz
		    ebscale = (*escali);   // E and B field scaling

		    ecount = 0;
                    for(vector<double>::const_iterator fcoefi = (*rffci).begin(); fcoefi != (*rffci).end(); ++fcoefi) {

			++ecount;
			frequency = frequency + (*fcoefi) * pow((t * 1e-9), ecount); // Add frequency ramp (in MHz/s^n)

		    }

		    ecount = 0;
                    for(vector<double>::const_iterator vcoefi = (*rfvci).begin(); vcoefi != (*rfvci).end(); ++vcoefi) {

			++ecount;
			ebscale = ebscale + (*vcoefi) * pow((t * 1e-9), ecount); // Add frequency ramp (in MHz/s^n)

		    }

                    double phase = 2.0 * pi * 1.0E-3 * frequency * t + (*rfphii);  // f in [MHz], t in [ns]

                    E += ebscale * cos(phase) * tmpE;
                    B -= ebscale * sin(phase) * tmpB;

		    // Some phase output -DW

		    if (tet >= 90.0 && waiting_for_gap == 1) {

                        double phase_print = 180.0 * phase / pi;
			phase_print = fmod(phase_print, 360) - 360.0;

			*gmsg << endl << "Gap 1 phase = " << phase_print << " Deg" << endl;
			*gmsg << "Gap 1 E-Field = (" << E[0] << "/" << E[1] << "/" << E[2] << ")" << endl;
			*gmsg << "Gap 1 B-Field = (" << B[0] << "/" << B[1] << "/" << B[2] << ")" << endl;
			*gmsg << "RF Frequency = " << frequency << " MHz" << endl;

			waiting_for_gap = 2;
		    }
		    else if (tet >= 270.0 && waiting_for_gap == 2) {

		        double phase_print = 180.0 * phase / pi;
			phase_print = fmod(phase_print, 360) - 360.0;

			*gmsg << endl << "Gap 2 phase = " << phase_print << " Deg" << endl;
			*gmsg << "Gap 2 E-Field = (" << E[0] << "/" << E[1] << "/" << E[2] << ")" << endl;
			*gmsg << "Gap 2 B-Field = (" << B[0] << "/" << B[1] << "/" << B[2] << ")" << endl;
			*gmsg << "RF Frequency = " << frequency << " MHz" << endl;
			waiting_for_gap = 0;
                    }

//                  INFOMSG("Field " << fcount << " BANDRF E= " << tmpE << " R= " << R << " phase " << phase << endl);
                }
            }
    	}
    }
    return false;
}

void Cyclotron::finalise() {

    online_m = false;
    lossDs_m->save();
    *gmsg << "* Finalize cyclotron" << endl;

}

bool Cyclotron::bends() const {
    return true;
}

// calculate derivatives with 5-point lagrange's formula.
double Cyclotron::gutdf5d(double *f, double dx, const int kor, const int krl, const int lpr)

{
    double C[5][5][3], FAC[3];
    double result;
    int j;
    /* CALCULATE DERIVATIVES WITH 5-POINT LAGRANGE FORMULA
     * PARAMETERS:
     * F  STARTADDRESS FOR THE 5 SUPPORT POINTS
     * DX STEPWIDTH FOR ARGUMENT
     * KOR        ORDER OF DERIVATIVE (KOR=1,2,3).
     * KRL        NUMBER OF SUPPORT POINT, WHERE THE DERIVATIVE IS TO BE CALCULATED
     *  (USUALLY 3, USE FOR BOUNDARY 1 ,2, RESP. 4, 5)
     * LPR        DISTANCE OF THE 5 STORAGE POSITIONS (=1 IF THEY ARE NEIGHBORS OR LENGTH
     * OF COLUMNLENGTH OF A MATRIX, IF THE SUPPORT POINTS ARE ON A LINE).
     * ATTENTION! THE INDICES ARE NOW IN C-FORMAT AND NOT IN FORTRAN-FORMAT.*/

    /* COEFFICIENTS FOR THE 1ST DERIVATIVE: */
    C[0][0][0] = -50.0;
    C[1][0][0] = 96.0;
    C[2][0][0] = -72.0;
    C[3][0][0] = 32.0;
    C[4][0][0] = -6.0;
    C[0][1][0] = -6.0;
    C[1][1][0] = -20.0;
    C[2][1][0] = 36.0;
    C[3][1][0] = -12.0;
    C[4][1][0] =  2.0;
    C[0][2][0] =  2.0;
    C[1][2][0] = -16.0;
    C[2][2][0] =  0.0;
    C[3][2][0] = 16.0;
    C[4][2][0] = -2.0;
    C[0][3][0] = -2.0;
    C[1][3][0] = 12.0;
    C[2][3][0] = -36.0;
    C[3][3][0] = 20.0;
    C[4][3][0] =  6.0;
    C[0][4][0] =  6.0;
    C[1][4][0] = -32.0;
    C[2][4][0] = 72.0;
    C[3][4][0] = -96.0;
    C[4][4][0] = 50.0;

    /* COEFFICIENTS FOR THE 2ND DERIVATIVE: */
    C[0][0][1] = 35.0;
    C[1][0][1] = -104;
    C[2][0][1] = 114.0;
    C[3][0][1] = -56.0;
    C[4][0][1] = 11.0;
    C[0][1][1] = 11.0;
    C[1][1][1] = -20.0;
    C[2][1][1] =  6.0;
    C[3][1][1] =  4.0;
    C[4][1][1] = -1.0;
    C[0][2][1] = -1.0;
    C[1][2][1] = 16.0;
    C[2][2][1] = -30.0;
    C[3][2][1] = 16.0;
    C[4][2][1] = -1.0;
    C[0][3][1] = -1.0;
    C[1][3][1] =  4.0;
    C[2][3][1] =  6.0;
    C[3][3][1] = -20.0;
    C[4][3][1] = 11.0;
    C[0][4][1] = 11.0;
    C[1][4][1] = -56.0;
    C[2][4][1] = 114.0;
    C[3][4][1] = -104;
    C[4][4][1] = 35.0;


    /* COEFFICIENTS FOR THE 3RD DERIVATIVE: */
    C[0][0][2] = -10.0;
    C[1][0][2] = 36.0;
    C[2][0][2] = -48.0;
    C[3][0][2] = 28.0;
    C[4][0][2] = -6.0;
    C[0][1][2] = -6.0;
    C[1][1][2] = 20.0;
    C[2][1][2] = -24.0;
    C[3][1][2] = 12.0;
    C[4][1][2] = -2.0;
    C[0][2][2] = -2.0;
    C[1][2][2] =  4.0;
    C[2][2][2] =  0.0;
    C[3][2][2] = -4.0;
    C[4][2][2] =  2.0;
    C[0][3][2] =  2.0;
    C[1][3][2] = -12.0;
    C[2][3][2] = 24.0;
    C[3][3][2] = -20.0;
    C[4][3][2] =  6.0;
    C[0][4][2] =  6.0;
    C[1][4][2] = -28.0;
    C[2][4][2] = 48.0;
    C[3][4][2] = -36.0;
    C[4][4][2] = 10.0;

    /* FACTOR: */
    FAC[0] = 24.0;
    FAC[1] = 12.0;
    FAC[2] = 4.0;

    result = 0.0;
    for(j = 0; j < 5; j++) {
        result += C[j][krl][kor] * *(f + j * lpr);
    }

    return result / (FAC[kor] * pow(dx, (kor + 1)));
}


// evaulate other derivative of magnetic field.
void Cyclotron::getdiffs() {

    //~ if(Bfield.dbr) delete[] Bfield.dbr;
    //~ Bfield.dbr   = new double[Bfield.ntot];
    //~ if(Bfield.dbrr) delete[] Bfield.dbrr;
    //~ Bfield.dbrr  = new double[Bfield.ntot];
    //~ if(Bfield.dbrrr) delete[] Bfield.dbrrr;
    //~ Bfield.dbrrr = new double[Bfield.ntot];
//~
    //~ if(Bfield.dbrt) delete[] Bfield.dbrt;
    //~ Bfield.dbrt  = new double[Bfield.ntot];
    //~ if(Bfield.dbrrt) delete[] Bfield.dbrrt;
    //~ Bfield.dbrrt = new double[Bfield.ntot];
    //~ if(Bfield.dbrtt) delete[] Bfield.dbrtt;
    //~ Bfield.dbrtt = new double[Bfield.ntot];
//~
    //~ if(Bfield.f2) delete[] Bfield.f2;
    //~ Bfield.f2    = new double[Bfield.ntot];
    //~ if(Bfield.f3) delete[] Bfield.f3;
    //~ Bfield.f3    = new double[Bfield.ntot];
    //~ if(Bfield.g3) delete[] Bfield.g3;
    //~ Bfield.g3    = new double[Bfield.ntot];
    Bfield.dbr.resize(Bfield.ntot);
    Bfield.dbrr.resize(Bfield.ntot);
    Bfield.dbrrr.resize(Bfield.ntot);

    Bfield.dbrt.resize(Bfield.ntot);
    Bfield.dbrrt.resize(Bfield.ntot);
    Bfield.dbrtt.resize(Bfield.ntot);

    Bfield.f2.resize(Bfield.ntot);
    Bfield.f3.resize(Bfield.ntot);
    Bfield.g3.resize(Bfield.ntot);

    for(int i = 0; i < Bfield.nrad; i++) {

        for(int k = 0; k < Bfield.ntet; k++) {

            double dtheta = pi / 180.0 * BP.dtet;

            int kEdge;

            kEdge = max(k - 2, 0);
            kEdge = min(kEdge, Bfield.ntet - 5);

            int dkFromEdge = k - kEdge;
            int index = idx(i, k);
            int indexkEdge = idx(i, kEdge);


            Bfield.dbt[index]    = gutdf5d(&Bfield.bfld[indexkEdge], dtheta, 0, dkFromEdge, 1);
            Bfield.dbtt[index]   = gutdf5d(&Bfield.bfld[indexkEdge], dtheta, 1, dkFromEdge, 1);
            Bfield.dbttt[index]  = gutdf5d(&Bfield.bfld[indexkEdge], dtheta, 2, dkFromEdge, 1);
        }
    }



    for(int k = 0; k < Bfield.ntet; k++) {
        // inner loop varies R
        for(int i = 0; i < Bfield.nrad; i++) {
            double rac = BP.rarr[i];
            // define iredg, the reference index for radial interpolation
            // standard: i-2 minimal: 0 (not negative!)  maximal: nrad-4
            int iredg = max(i - 2, 0);
            iredg = min(iredg, Bfield.nrad - 5);
            int irtak = i - iredg;
            int index = idx(i, k);
            int indexredg = idx(iredg, k);


            Bfield.dbr[index]    = gutdf5d(&Bfield.bfld[indexredg], BP.delr, 0, irtak, Bfield.ntetS);
            Bfield.dbrr[index]   = gutdf5d(&Bfield.bfld[indexredg], BP.delr, 1, irtak, Bfield.ntetS);
            Bfield.dbrrr[index]  = gutdf5d(&Bfield.bfld[indexredg], BP.delr, 2, irtak, Bfield.ntetS);

            Bfield.dbrt[index]   = gutdf5d(&Bfield.dbt[indexredg], BP.delr, 0, irtak, Bfield.ntetS);
            Bfield.dbrrt[index]  = gutdf5d(&Bfield.dbt[indexredg], BP.delr, 1, irtak, Bfield.ntetS);
            Bfield.dbrtt[index]  = gutdf5d(&Bfield.dbtt[indexredg], BP.delr, 0, irtak, Bfield.ntetS);

            // fehlt noch!! f2,f3,g3,
            Bfield.f2[index] = (Bfield.dbrr[index]
                                + Bfield.dbr[index] / rac
                                + Bfield.dbtt[index] / rac / rac) / 2.0;

            Bfield.f3[index] = (Bfield.dbrrr[index]
                                + Bfield.dbrr[index] / rac
                                + (Bfield.dbrtt[index] - Bfield.dbr[index]) / rac / rac
                                - 2.0 * Bfield.dbtt[index] / rac / rac / rac) / 6.0;

            Bfield.g3[index] = (Bfield.dbrrt[index]
                                + Bfield.dbrt[index] / rac
                                + Bfield.dbttt[index] / rac / rac) / 6.0;
        } // Radius Loop
    } // Azimuth loop

    // copy 1st azimuth to last + 1 to always yield an interval
    for(int i = 0; i < Bfield.nrad; i++) {
        int iend = idx(i, Bfield.ntet);
        int istart = idx(i, 0);

        Bfield.bfld[iend]   = Bfield.bfld[istart];
        Bfield.dbt[iend]    = Bfield.dbt[istart];
        Bfield.dbtt[iend]   = Bfield.dbtt[istart];
        Bfield.dbttt[iend]  = Bfield.dbttt[istart];

        Bfield.dbr[iend]    = Bfield.dbr[istart];
        Bfield.dbrr[iend]   = Bfield.dbrr[istart];
        Bfield.dbrrr[iend]  = Bfield.dbrrr[istart];

        Bfield.dbrt[iend]   = Bfield.dbrt[istart];
        Bfield.dbrtt[iend]  = Bfield.dbrtt[istart];
        Bfield.dbrrt[iend]  = Bfield.dbrrt[istart];

        Bfield.f2[iend]     = Bfield.f2[istart];
        Bfield.f3[iend]     = Bfield.f3[istart];
        Bfield.g3[iend]     = Bfield.g3[istart];

    }

    /* debug

    for(int i = 0; i< Bfield.nrad; i++){
      for(int j = 0; j< Bfield.ntetS; j++){
    int index = idx(i,j);
    double x = (BP.rmin+i*BP.delr) * sin(j*BP.dtet*pi/180.0);
    double y = (BP.rmin+i*BP.delr) * cos(j*BP.dtet*pi/180.0);
    *gmsg<<"x= "<<x<<" y= "<<y<<" B= "<<Bfield.bfld[index]<<endl;
      }
    }
    */
}

// read field map from external file.
void Cyclotron::getFieldFromFile(const double &scaleFactor) {

    FILE *f = NULL;
    int lpar;
    char fout[100];
    double dtmp;

    *gmsg << "* ----------------------------------------------" << endl;
    *gmsg << "*             READ IN RING FIELD MAP            " << endl;
    *gmsg << "*      (The first data block is useless)        " << endl;
    *gmsg << "* ----------------------------------------------" << endl;

    BP.Bfact = scaleFactor;

    if((f = fopen(fmapfn_m.c_str(), "r")) == NULL) {
        throw GeneralClassicException("Cyclotron::getFieldFromField",
                                      "failed to open file '" + fmapfn_m + "', please check if it exists");
    }

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.rmin));
    *gmsg << "* Minimal radius of measured field map: " << BP.rmin << " [mm]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.delr));
    //if the value is nagtive, the actual value is its reciprocal.
    if(BP.delr < 0.0) BP.delr = 1.0 / (-BP.delr);
    *gmsg << "* Stepsize in radial direction: " << BP.delr << " [mm]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.tetmin));
    *gmsg << "* Minimal angle of measured field map: " << BP.tetmin << " [deg.]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.dtet));
    //if the value is nagtive, the actual value is its reciprocal.
    if(BP.dtet < 0.0) BP.dtet = 1.0 / (-BP.dtet);
    *gmsg << "* Stepsize in azimuth direction: " << BP.dtet << " [deg.]" << endl;

    for(int i = 0; i < 13; i++)
        CHECK_CYC_FSCANF_EOF(fscanf(f, "%s", fout));

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%d", &Bfield.nrad));
    *gmsg << "* Index in radial direction: " << Bfield.nrad << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%d", &Bfield.ntet));
    *gmsg << "* Index in azimuthal direction: " << Bfield.ntet << endl;

    Bfield.ntetS = Bfield.ntet + 1;
    *gmsg << "* Accordingly, total grid point along azimuth:  " << Bfield.ntetS << endl;

    for(int i = 0; i < 5; i++) {
        CHECK_CYC_FSCANF_EOF(fscanf(f, "%s", fout));
    }
    CHECK_CYC_FSCANF_EOF(fscanf(f, "%d", &lpar));
    // msg<< "READ"<<lpar<<" DATA ENTRIES"<<endl;

    for(int i = 0; i < 4; i++) {
        CHECK_CYC_FSCANF_EOF(fscanf(f, "%s", fout));
    }

    for(int i = 0; i < lpar; i++) {
        CHECK_CYC_FSCANF_EOF(fscanf(f, "%16lE", &dtmp));
    }
    for(int i = 0; i < 6; i++) {
        CHECK_CYC_FSCANF_EOF(fscanf(f, "%s", fout));
    }
    //*gmsg << "* READ FILE DESCRIPTION..." <<endl;
    for(int i = 0; i < 10000; i++) {
        CHECK_CYC_FSCANF_EOF(fscanf(f, "%s", fout));
        if(strcmp(fout, "LREC=") == 0)break;
    }

    for(int i = 0; i < 5; i++) {
        CHECK_CYC_FSCANF_EOF(fscanf(f, "%s", fout));
    }
    Bfield.ntot = idx(Bfield.nrad - 1, Bfield.ntet) + 1;
    //jjyang
    *gmsg << "* Total stored grid point number ( ntetS * nrad ) : " << Bfield.ntot << endl;

    Bfield.bfld.resize(Bfield.ntot);
    Bfield.dbt.resize(Bfield.ntot);
    Bfield.dbtt.resize(Bfield.ntot);
    Bfield.dbttt.resize(Bfield.ntot);

    *gmsg << "* Read-in loop one block per radius" << endl;
    *gmsg << "* Rescaling of the fields with factor: " << BP.Bfact << endl;
    for(int i = 0; i < Bfield.nrad; i++) {

        if(i > 0) {
            for(int dummy = 0; dummy < 6; dummy++) {
                CHECK_CYC_FSCANF_EOF(fscanf(f, "%s", fout)); // INFO-LINE
            }
        }
        for(int k = 0; k < Bfield.ntet; k++) {
            CHECK_CYC_FSCANF_EOF(fscanf(f, "%16lE", &(Bfield.bfld[idx(i, k)])));
            Bfield.bfld[idx(i, k)] *= BP.Bfact;
        }
        for(int k = 0; k < Bfield.ntet; k++) {
            CHECK_CYC_FSCANF_EOF(fscanf(f, "%16lE", &(Bfield.dbt[idx(i, k)])));
            Bfield.dbt[idx(i, k)] *= BP.Bfact;
        }
        for(int k = 0; k < Bfield.ntet; k++) {
            CHECK_CYC_FSCANF_EOF(fscanf(f, "%16lE", &(Bfield.dbtt[idx(i, k)])));
            Bfield.dbtt[idx(i, k)] *= BP.Bfact;
        }
        for(int k = 0; k < Bfield.ntet; k++) {
            CHECK_CYC_FSCANF_EOF(fscanf(f, "%16lE", &(Bfield.dbttt[idx(i, k)])));
            Bfield.dbttt[idx(i, k)] *= BP.Bfact;
        }
    }
    fclose(f);


    *gmsg << "* Field Map read successfully!" << endl << endl;
}



// Calculates Radiae of initial grid.
// dimensions in [mm]!
void Cyclotron::initR(double rmin, double dr, int nrad) {
    BP.rarr.resize(nrad);
    for(int i = 0; i < nrad; i++) {
        BP.rarr[i] = rmin + i * dr;
    }
    BP.delr = dr;
}

void Cyclotron::initialise(PartBunch *bunch, double &startField, double &endField) {
    RefPartBunch_m = bunch;
    online_m = true;
}

void Cyclotron::initialise(PartBunch *bunch, const int &fieldflag, const double &scaleFactor) {
    RefPartBunch_m = bunch;
    lossDs_m = std::unique_ptr<LossDataSink>(new LossDataSink(getName(), !Options::asciidump));

    //    PSIBF, AVFEQBF, ANSYSBF, FFAGBF
    // for your own format field, you should add your own getFieldFromFile() function by yourself.

    if(fieldflag == 1) {
        //*gmsg<<"Read field data from PSI format field map file."<<endl;
        myBFieldType_m = PSIBF;
        getFieldFromFile(scaleFactor);

    } else if(fieldflag == 2) {
        // *gmsg<<"Read data from 450MeV Carbon cyclotron field file"<<endl
        myBFieldType_m = CARBONBF;
        getFieldFromFile_Carbon(scaleFactor);

    } else if(fieldflag == 3) {
        // *gmsg<<"Read data from 100MeV H- cyclotron CYCIAE-100 field file"<<endl;
        myBFieldType_m = ANSYSBF;
        getFieldFromFile_CYCIAE(scaleFactor);

    } else if(fieldflag == 4) {
        *gmsg << "* Read AVFEQ data (Riken) use bfield scale factor bs = " << getBScale() << endl;
        myBFieldType_m = AVFEQBF;
        getFieldFromFile_AVFEQ(scaleFactor);

    } else if(fieldflag == 5) {
        *gmsg << "* Read FFAG data MSU/FNAL " << getBScale() << endl;
        myBFieldType_m = FFAGBF;
        getFieldFromFile_FFAG(scaleFactor);

    } else if(fieldflag == 6) {
        *gmsg << "* Read both median plane B field map and 3D E field map of RF cavity for compact cyclotron" << getBScale() << endl;
        myBFieldType_m = BANDRF;
        getFieldFromFile_BandRF(scaleFactor);

    } else if(fieldflag == 7) {
        *gmsg << "* Read midplane B-field, 3D RF fieldmaps, and text files with RF frequency/Voltage coefficients for Synchrocyclotron. (Midplane scaling = " << getBScale() << ")" << endl;
        myBFieldType_m = SYNCHRO;
        getFieldFromFile_Synchrocyclotron(scaleFactor);

    } else
        ERRORMSG("* The field reading function of this TYPE of CYCLOTRON has not implemented yet!" << endl);

    // calculate the radii of initial grid.
    initR(BP.rmin, BP.delr, Bfield.nrad);

    // calculate the remaining derivatives
    getdiffs();

}


void Cyclotron::getFieldFromFile_FFAG(const double &scaleFactor) {

    /*
      Field is read in from ascci file (COSY output) in the oder:
      R(m) theta(Deg) x(m) y(m) Bz(T).

      Theta is the fast varing variable

      2.0000   0.0  2.0000  0.0000      0.0000000000000000
      2.0000   1.0  1.9997  0.0349      0.0000000000000000
      2.0000   2.0  1.9988  0.0698      0.0000000000000000
      2.0000   3.0  1.9973  0.1047      0.0000000000000000

      ......
      <blank line>

      2.1000   0.0  2.1000  0.0000      0.0000000000000000
      2.1000   1.0  2.0997  0.0367      0.0000000000000000
      2.1000   2.0  2.0987  0.0733      0.0000000000000000
    */



    vector<double> rv;
    vector<double> thv;
    vector<double> xv;
    vector<double> yv;
    vector<double> bzv;
    vector<double>::iterator vit;

    *gmsg << "* ----------------------------------------------" << endl;
    *gmsg << "*    READ IN FFAG FIELD MAP     " << endl;
    *gmsg << "* ----------------------------------------------" << endl;

    BP.Bfact = -10.0; // T->kG and H- for the current FNAL FFAG

    ifstream file_to_read(fmapfn_m.c_str());
    const int max_num_of_char_in_a_line = 128;
    const int num_of_header_lines = 1;

    // STEP2: SKIP ALL THE HEADER LINES
    for(int i = 0; i < num_of_header_lines; ++i)
        file_to_read.ignore(max_num_of_char_in_a_line, '\n');

    while(!file_to_read.eof()) {
        double r, th, x, y, bz;
        file_to_read >> r >> th >> x >> y >> bz;
        if((int)th != 360) {
            rv.push_back(r * 1000.0);
            thv.push_back(th);
            xv.push_back(x * 1000.0);
            yv.push_back(y * 1000.0);
            bzv.push_back(bz);
        }
    }

    double maxtheta = 360.0;
    BP.dtet = thv[1] - thv[0];
    BP.rmin = *(rv.begin());
    double rmax = rv.back();

    // find out dR
    for(vit = rv.begin(); *vit <= BP.rmin; vit++) {}
    BP.delr = *vit - BP.rmin;

    BP.tetmin = thv[0];

    Bfield.ntet = (int)((maxtheta - thv[0]) / BP.dtet);
    Bfield.nrad  = (int)(rmax - BP.rmin) / BP.delr + 1;
    Bfield.ntetS  = Bfield.ntet + 1;
    *gmsg << "* Minimal radius of measured field map: " << BP.rmin << " [mm]" << endl;
    *gmsg << "* Maximal radius of measured field map: " << rmax << " [mm]" << endl;
    *gmsg << "* Stepsize in radial direction: " << BP.delr << " [mm]" << endl;
    *gmsg << "* Minimal angle of measured field map: " << BP.tetmin << " [deg.]" << endl;
    *gmsg << "* Maximal angle of measured field map: " << maxtheta << " [deg.]" << endl;

    //if the value is negtive, the actual value is its reciprocal.
    if(BP.dtet < 0.0) BP.dtet = 1.0 / (-BP.dtet);
    *gmsg << "* Stepsize in azimuth direction: " << BP.dtet << " [deg.]" << endl;
    *gmsg << "* Total grid point along azimuth:  " << Bfield.ntetS << endl;
    *gmsg << "* Total grid point along radius: " << Bfield.nrad << endl;

    Bfield.ntot = Bfield.ntetS * Bfield.nrad;
    *gmsg << "* Total stored grid point number ( ntetS * nrad ) : " << Bfield.ntot << endl;

    Bfield.bfld.resize(Bfield.ntot);
    Bfield.dbt.resize(Bfield.ntot);
    Bfield.dbtt.resize(Bfield.ntot);
    Bfield.dbttt.resize(Bfield.ntot);

    *gmsg << "* Rescaling of the fields with factor: " << BP.Bfact << endl;

    int count = 0;
    if((Ippl::getNodes()) == 1 && Options::info) {
        fstream fp;
        fp.open("data/gnu.out", ios::out);

        for(int r = 0; r < Bfield.nrad; r++) {
            for(int k = 0; k < Bfield.ntet; k++) {
                Bfield.bfld[idx(r, k)] = bzv[count] * BP.Bfact;
                fp << BP.rmin + (r * BP.delr) << " \t " << k*(BP.tetmin + BP.dtet) << " \t " << Bfield.bfld[idx(r, k)] << endl;
                count++;
            }
        }
        fp.close();
    }

    *gmsg << "* Field Map read successfully nelem= " << count << endl << endl;
}

void Cyclotron::getFieldFromFile_AVFEQ(const double &scaleFactor) {

    FILE *f = NULL;
    *gmsg << "* ----------------------------------------------" << endl;
    *gmsg << "*    READ IN AVFEQ CYCLOTRON FIELD MAP     " << endl;
    *gmsg << "* ----------------------------------------------" << endl;

    /*  From Hiroki-san
        The first line tells r minimum (500mm),
                             r maximum(4150mm),
                             r step(50mm),
                             theta minimum(0deg),
                             theta maximum(90deg)
                             theta step(0.5deg).

        From the next line data repeat the block for a given r which the first line of the block tells.
        Each block consists of the data Bz from theta minimum (0deg) to theta maximum(90deg) with theta step(0.5deg).
    */

    BP.Bfact = scaleFactor / 1000.;

    if((f = fopen(fmapfn_m.c_str(), "r")) == NULL) {
        throw GeneralClassicException("Cyclotron::getFieldFromFile_AVFEQ",
                                      "failed to open file '" + fmapfn_m + "', please check if it exists");
    }

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.rmin));
    *gmsg << "* Minimal radius of measured field map: " << BP.rmin << " [mm]" << endl;

    double rmax;
    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &rmax));
    *gmsg << "* Maximal radius of measured field map: " << rmax << " [mm]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.delr));
    *gmsg << "* Stepsize in radial direction: " << BP.delr << " [mm]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.tetmin));
    *gmsg << "* Minimal angle of measured field map: " << BP.tetmin << " [deg.]" << endl;

    double tetmax;
    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &tetmax));
    *gmsg << "* Maximal angle of measured field map: " << tetmax << " [deg.]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.dtet));
    //if the value is nagtive, the actual value is its reciprocal.

    if(BP.dtet < 0.0) BP.dtet = 1.0 / (-BP.dtet);
    *gmsg << "* Stepsize in azimuth direction: " << BP.dtet << " [deg.]" << endl;

    Bfield.ntetS = (int)((tetmax - BP.tetmin) / BP.dtet + 1);
    *gmsg << "* Total grid point along azimuth:  " << Bfield.ntetS << endl;

    Bfield.nrad = (int)(rmax - BP.rmin) / BP.delr;


    int ntotidx = idx(Bfield.nrad, Bfield.ntetS) + 1;

    Bfield.ntot = Bfield.ntetS * Bfield.nrad;
    *gmsg << "* Total stored grid point number ( ntetS * nrad ) : " << Bfield.ntot << " ntot-idx= " << ntotidx << endl;

    Bfield.bfld.resize(Bfield.ntot);
    Bfield.dbt.resize(Bfield.ntot);
    Bfield.dbtt.resize(Bfield.ntot);
    Bfield.dbttt.resize(Bfield.ntot);

    *gmsg << "* rescaling of the fields with factor: " << BP.Bfact << endl;

    fstream fp;
    if((Ippl::getNodes()) == 1 && Options::info)
      fp.open("data/gnu.out", ios::out);

    double tmp;
    int count = 0;

    for(int r = 0; r < Bfield.nrad; r++) {
      CHECK_CYC_FSCANF_EOF(fscanf(f, "%16lE", &tmp));   // over read
      for(int k = 0; k < Bfield.ntetS; k++) {
	CHECK_CYC_FSCANF_EOF(fscanf(f, "%16lE", &(Bfield.bfld[idx(r, k)])));
	Bfield.bfld[idx(r, k)] *= BP.Bfact;
	if((Ippl::getNodes()) == 1 && Options::info)
	  fp << BP.rmin + (r * BP.delr) << " \t " << k*(BP.tetmin + BP.dtet) << " \t " << Bfield.bfld[idx(r, k)] << " idx= " << idx(r, k)  << endl;
	count++;
      }
    }
    if((Ippl::getNodes()) == 1 && Options::info)
      fp.close();
    fclose(f);
    *gmsg << "* Field Map read successfully nelem= " << count << endl << endl;
}


// read field map from external file.
void Cyclotron::getFieldFromFile_Carbon(const double &scaleFactor) {

    FILE *f = NULL;
    *gmsg << "* ----------------------------------------------" << endl;
    *gmsg << "*      READ IN CARBON CYCLOTRON FIELD MAP       " << endl;
    *gmsg << "* ----------------------------------------------" << endl;

    BP.Bfact = scaleFactor;

    if((f = fopen(fmapfn_m.c_str(), "r")) == NULL) {
        throw GeneralClassicException("Cyclotron::getFieldFromFile_Carbon",
                                      "failed to open file '" + fmapfn_m + "', please check if it exists");
    }

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.rmin));
    *gmsg << "* Minimal radius of measured field map: " << BP.rmin << " [mm]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.delr));
    //if the value is negative, the actual value is its reciprocal.
    if(BP.delr < 0.0) BP.delr = 1.0 / (-BP.delr);
    *gmsg << "* Stepsize in radial direction: " << BP.delr << " [mm]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.tetmin));
    *gmsg << "* Minimal angle of measured field map: " << BP.tetmin << " [deg]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.dtet));
    //if the value is negative, the actual value is its reciprocal.
    if(BP.dtet < 0.0) BP.dtet = 1.0 / (-BP.dtet);
    *gmsg << "* Stepsize in azimuthal direction: " << BP.dtet << " [deg]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%d", &Bfield.ntet));
    *gmsg << "* Grid points along azimuth (ntet): " << Bfield.ntet << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%d", &Bfield.nrad));
    *gmsg << "* Grid points along radius (nrad): " << Bfield.nrad << endl;

//    Bfield.ntetS = Bfield.ntet;
    Bfield.ntetS = Bfield.ntet + 1;
    //*gmsg << "* Accordingly, total grid point along azimuth:  " << Bfield.ntetS << endl;

    //Bfield.ntot = idx(Bfield.nrad - 1, Bfield.ntet) + 1;
    Bfield.ntot = Bfield.nrad * Bfield.ntetS;

    *gmsg << "* Adding a guard cell along azimuth" << endl;
    *gmsg << "* Total stored grid point number ((ntet+1) * nrad) : " << Bfield.ntot << endl;
    Bfield.bfld.resize(Bfield.ntot);
    Bfield.dbt.resize(Bfield.ntot);
    Bfield.dbtt.resize(Bfield.ntot);
    Bfield.dbttt.resize(Bfield.ntot);

    *gmsg << "* rescaling of the fields with factor: " << BP.Bfact << endl;

    for(int i = 0; i < Bfield.nrad; i++) {
        for(int k = 0; k < Bfield.ntet; k++) {
            CHECK_CYC_FSCANF_EOF(fscanf(f, "%16lE", &(Bfield.bfld[idx(i, k)])));
            Bfield.bfld[idx(i, k)] *= BP.Bfact;
        }
    }

    if((Ippl::getNodes()) == 1 && Options::info) {
        fstream fp1, fp2;
        fp1.open("data/gnu.out", ios::out);
        fp2.open("data/eb.out", ios::out);
        for(int i = 0; i < Bfield.nrad; i++) {
            for(int k = 0; k < Bfield.ntet; k++) {
                fp1 << BP.rmin + (i * BP.delr) << " \t " << k * (BP.tetmin + BP.dtet) << " \t " << Bfield.bfld[idx(i, k)] << endl;
              
                Vector_t tmpR = Vector_t (BP.rmin + (i * BP.delr), 0.0, k * (BP.tetmin + BP.dtet));
                Vector_t tmpE(0.0, 0.0, 0.0), tmpB(0.0, 0.0, 0.0);
                tmpR /= 1000.0; // -> mm to m
                vector<Fieldmap *>::const_iterator fi  = RFfields_m.begin();
                vector<double>::const_iterator rffi    = rffrequ_m.begin();
                vector<double>::const_iterator rfphii  = rfphi_m.begin();
                vector<double>::const_iterator escali  = escale_m.begin();
                for(; fi != RFfields_m.end(); ++fi, ++rffi, ++rfphii, ++escali) {
                    Vector_t E(0.0, 0.0, 0.0), B(0.0, 0.0, 0.0);
                    if(!(*fi)->getFieldstrength(tmpR, tmpE, tmpB)) {
                        tmpE += E;
                        tmpB -= B;
                    }
                }
                fp2 << tmpR  <<  " \t E= " << tmpE << "\t B= " << tmpB << endl;
            }
        }
        fp1.close();
        fp2.close();
    }
    
    fclose(f);
    
    *gmsg << "* Field Maps read successfully!" << endl << endl;
}


// read field map from external file.
void Cyclotron::getFieldFromFile_CYCIAE(const double &scaleFactor) {

    FILE *f = NULL;
    char fout[100];
    int dtmp;

    *gmsg << "* ----------------------------------------------" << endl;
    *gmsg << "*    READ IN CYCIAE-100 CYCLOTRON FIELD MAP     " << endl;
    *gmsg << "* ----------------------------------------------" << endl;

    BP.Bfact = scaleFactor;

    if((f = fopen(fmapfn_m.c_str(), "r")) == NULL) {
        throw GeneralClassicException("Cyclotron::getFieldFromFile_CYCIAE",
                                      "failed to open file '" + fmapfn_m + "', please check if it exists");
    }

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.rmin));
    *gmsg << "* Minimal radius of measured field map: " << BP.rmin << " [mm]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.delr));
    *gmsg << "* Stepsize in radial direction: " << BP.delr << " [mm]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.tetmin));
    *gmsg << "* Minimal angle of measured field map: " << BP.tetmin << " [deg.]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &BP.dtet));
    //if the value is nagtive, the actual value is its reciprocal.
    if(BP.dtet < 0.0) BP.dtet = 1.0 / (-BP.dtet);
    *gmsg << "* Stepsize in azimuth direction: " << BP.dtet << " [deg.]" << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%d", &Bfield.ntet));
    *gmsg << "* Index in azimuthal direction: " << Bfield.ntet << endl;

    CHECK_CYC_FSCANF_EOF(fscanf(f, "%d", &Bfield.nrad));
    *gmsg << "* Index in radial direction: " << Bfield.nrad << endl;

    Bfield.ntetS = Bfield.ntet + 1;
    *gmsg << "* Accordingly, total grid point along azimuth:  " << Bfield.ntetS << endl;

    Bfield.ntot = idx(Bfield.nrad - 1, Bfield.ntet) + 1;

    *gmsg << "* Total stored grid point number ( ntetS * nrad ) : " << Bfield.ntot << endl;
    Bfield.bfld.resize(Bfield.ntot);
    Bfield.dbt.resize(Bfield.ntot);
    Bfield.dbtt.resize(Bfield.ntot);
    Bfield.dbttt.resize(Bfield.ntot);

    *gmsg << "* rescaling of the fields with factor: " << BP.Bfact << endl;

    int nHalfPoints = Bfield.ntet / 2.0 + 1;

    for(int i = 0; i < Bfield.nrad; i++) {

        for(int ii = 0; ii < 13; ii++)
			CHECK_CYC_FSCANF_EOF(fscanf(f, "%s", fout));

        for(int k = 0; k < nHalfPoints; k++) {
            CHECK_CYC_FSCANF_EOF(fscanf(f, "%d", &dtmp));
            CHECK_CYC_FSCANF_EOF(fscanf(f, "%d", &dtmp));
            CHECK_CYC_FSCANF_EOF(fscanf(f, "%d", &dtmp));
            CHECK_CYC_FSCANF_EOF(fscanf(f, "%lf", &(Bfield.bfld[idx(i, k)])));
            Bfield.bfld[idx(i, k)] = Bfield.bfld[idx(i, k)] * (-10.0); //  T --> kGs, minus for minus hydrongen
        }

        for(int k = nHalfPoints; k < Bfield.ntet; k++) {
            Bfield.bfld[idx(i, k)] = Bfield.bfld[idx(i, Bfield.ntet-k)];
        }
    }
    //  for(int i=0; i < 300; i++) msg <<"i="<<i<<", Bfield = "<< Bfield.bfld[i]<<endl;

    fclose(f);

    *gmsg << "* Field Map read successfully!" << endl << endl;
}

void Cyclotron::getFieldFromFile_BandRF(const double &scaleFactor) {

    // read 3D E&B field data file
    vector<string>::const_iterator fm    = RFfilename_m.begin();
    // loop over all field maps and superpose fields
    vector<double>::const_iterator rffi    = rffrequ_m.begin();
    vector<double>::const_iterator rfphii  = rfphi_m.begin();
    vector<double>::const_iterator escali  = escale_m.begin();
    int fcount = 0;

    for(; fm != RFfilename_m.end(); ++fm, ++rffi, ++rfphii, ++escali, ++fcount) {
        Fieldmap *f = Fieldmap::getFieldmap(*fm, false);
        if(f == NULL) {
            throw GeneralClassicException("Cyclotron::getFieldFromFile_BandRF",
                                          "failed to open file '" + *fm + "', please check if it exists");
        }
        f->readMap();
	//	if (IPPL::Comm->getOutputLevel() != 0)
	//  f->getInfo(gmsg);
        RFfields_m.push_back(f);
    }
    // read CARBON type B field
    getFieldFromFile_Carbon(scaleFactor);
}

void Cyclotron::getFieldFromFile_Synchrocyclotron(const double &scaleFactor) {

    // read 3D E&B field data file
    vector<string>::const_iterator fm    = RFfilename_m.begin();
    vector<string>::const_iterator rffcfni = RFFCoeff_fn_m.begin();
    vector<string>::const_iterator rfvcfni = RFVCoeff_fn_m.begin();
    // loop over all field maps and superpose fields
    vector<double>::const_iterator rffi    = rffrequ_m.begin();
    vector<double>::const_iterator rfphii  = rfphi_m.begin();
    vector<double>::const_iterator escali  = escale_m.begin();
    int fcount = 0;
    FILE *rffcf = NULL;
    FILE *rfvcf = NULL;

    *gmsg << endl;
    *gmsg << "* ------------------------------------------------------------" << endl;
    *gmsg << "*      READ IN 3D RF Fields and Frequency Coefficients        " << endl;
    *gmsg << "* ------------------------------------------------------------" << endl;

    for(; fm != RFfilename_m.end(); ++fm, ++rffi, ++rfphii, ++escali, ++rffcfni, ++rfvcfni, ++fcount) {
        Fieldmap *f = Fieldmap::getFieldmap(*fm, false);
        if(f == NULL) {
            throw GeneralClassicException("Cyclotron::getFieldFromFile_Synchrocyclotron",
                                          "failed to open file '" + *fm + "', please check if it exists");
        }
        f->readMap();
	// if (IPPL::Comm->getOutputLevel() != 0) f->getInfo(gmsg);
        RFfields_m.push_back(f);

        // Read RF Frequency Coefficients from file
	*gmsg << "RF Frequency Coefficient Filename: " << (*rffcfni) << endl;

	rffcf = fopen((*rffcfni).c_str(), "r");

        if(rffcf == NULL) {
            throw GeneralClassicException("Cyclotron::getFieldFromFile_Synchrocyclotron",
                                          "failed to open file '" + *rffcfni + "', please check if it exists");
        }

	vector<double> fcoeff;

	int nc; //Number of coefficients
	double value;

        CHECK_CYC_FSCANF_EOF(fscanf(rffcf, "%d", &nc));
        *gmsg << "* Number of coefficients in file: " << nc << endl;
        for(int k = 0; k < nc; k++) {
	    CHECK_CYC_FSCANF_EOF(fscanf(rffcf, "%16lE", &value));
	    fcoeff.push_back(value);
            //*gmsg << "* Coefficient " << k << ": " << value << endl;
        }
	rffc_m.push_back(fcoeff);

	fclose(rffcf);

        // Read RF Voltage Coefficients from file
	*gmsg << "RF Voltage Coefficient Filename: " << (*rfvcfni) << endl;

	rfvcf = fopen((*rfvcfni).c_str(), "r");

        if(rfvcf == NULL) {
            throw GeneralClassicException("Cyclotron::getFieldFromFile_Synchrocyclotron",
                                          "failed to open file '" + *rfvcfni + "', please check if it exists");
        }

	vector<double> vcoeff;

        CHECK_CYC_FSCANF_EOF(fscanf(rfvcf, "%d", &nc));
        *gmsg << "* Number of coefficients in file: " << nc << endl;
        for(int k = 0; k < nc; k++) {
	    CHECK_CYC_FSCANF_EOF(fscanf(rfvcf, "%16lE", &value));
	    vcoeff.push_back(value);
            //*gmsg << "* Coefficient " << k << ": " << value << endl;
        }
	rfvc_m.push_back(vcoeff);

	fclose(rfvcf);
    }

    // read CARBON type B field for mid-plane field
    getFieldFromFile_Carbon(scaleFactor);
}

void Cyclotron::getDimensions(double &zBegin, double &zEnd) const
{ }

#undef CHECK_CYC_FSCANF_EOF