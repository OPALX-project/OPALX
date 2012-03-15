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
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Cyclotron.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Physics/Physics.h"
#include "assert.h"

extern Inform *gmsg;

using Physics::pi;

// Class Cyclotron
// ------------------------------------------------------------------------

Cyclotron::Cyclotron():
    Component()
{}


Cyclotron::Cyclotron(const Cyclotron &right):
    Component(right),
    rinit_m(right.rinit_m),
    prinit_m(right.prinit_m),
    phiinit_m(right.phiinit_m),
    fmapfn_m(right.fmapfn_m),
    rffrequ_m(right.rffrequ_m),
    symmetry_m(right.symmetry_m),
    type_m(right.type_m),
    harm_m(right.harm_m),
    bscale_m(right.bscale_m),
    tcr1_m(right.tcr1_m),
    tcr2_m(right.tcr2_m),
    mbtc_m(right.mbtc_m),
    slptc_m(right.slptc_m)
{}



Cyclotron::Cyclotron(const string &name):
    Component(name)
{}


Cyclotron::~Cyclotron() {
    delete[] Bfield.bfld;
    delete[] Bfield.dbt;
    delete[] Bfield.dbtt;
    delete[] Bfield.dbttt;

    delete[] BP.rarr;

    delete[] Bfield.dbr;
    delete[] Bfield.dbrr;
    delete[] Bfield.dbrrr;

    delete[] Bfield.dbrt;
    delete[] Bfield.dbrrt;
    delete[] Bfield.dbrtt;

    delete[] Bfield.f2;
    delete[] Bfield.f3;
    delete[] Bfield.g3;

}


void Cyclotron::accept(BeamlineVisitor &visitor) const {
    visitor.visitCyclotron(*this);
}

/**
 *
 *
 * @param rinit initial set the radius of the beam (m)
 */

void Cyclotron::setRinit(double rinit) {
    rinit_m = rinit;
}
/**
 *
 *
 *
 * @return get the initial radius of the beam (m)
 */

double Cyclotron::getRinit() {
    return rinit_m;
}


void Cyclotron::setPRinit(double prinit) {
    prinit_m = prinit;
}

double Cyclotron::getPRinit() {
    return prinit_m;
}

void Cyclotron::setPHIinit(double phiinit) {
    phiinit_m = phiinit;
}

double Cyclotron::getPHIinit() {
    return phiinit_m;
}

void Cyclotron::setFieldMapFN(string f) {
    fmapfn_m = f;
}

string Cyclotron::getFieldMapFN() {
    return fmapfn_m;
}

void Cyclotron::setRfFrequ(double f) {
    rffrequ_m = f;
}

double Cyclotron::getRfFrequ() {
    return rffrequ_m;
}

void Cyclotron::setSymmetry(double s) {
    symmetry_m = s;
}

double Cyclotron::getSymmetry() {
    return symmetry_m;
}


void Cyclotron::setType(string t) {
    type_m = t;
}

const string &Cyclotron::getType() const {
    return type_m;
}

void Cyclotron::setCyclHarm(double h) {
    harm_m = h;
}

void Cyclotron::setBScale(double s) {
    bscale_m = s;
}

double Cyclotron::getBScale() {
    return bscale_m;
}


double Cyclotron::getCyclHarm() {
    return harm_m;
}

double Cyclotron::getRmin() {
    return BP.rmin;
}


double Cyclotron::getRmax() {
    return BP.rmin + (Bfield.nrad - 1) * BP.delr;
}

// This function aims at obtaining magentic field at any given location R by interpolation.
// arguments t is useless here.
// arguments E is set to zero.


void Cyclotron::setTCr1(double tcr1) {
    tcr1_m = tcr1;
}

double Cyclotron::getTCr1() {
    return tcr1_m;
}

void Cyclotron::setTCr2(double tcr2) {
    tcr2_m = tcr2;
}

double Cyclotron::getTCr2() {
    return tcr2_m;
}
void Cyclotron::setMBtc(double mbtc) {
    mbtc_m = mbtc;
}

double Cyclotron::getMBtc() {
    return mbtc_m;
}

void Cyclotron::setSLPtc(double slptc) {
    slptc_m = slptc;
}

double Cyclotron::getSLPtc() {
    return slptc_m;
}

bool Cyclotron::apply(const int &i, const double &t, double E[], double B[]) {
    Vector_t Ev(0, 0, 0), Bv(0, 0, 0);
    if(apply(RefPartBunch_m->R[i], Vector_t(0.0), t, Ev, Bv)) return true;

    E[0] = Ev(0);
    E[1] = Ev(1);
    E[2] = Ev(2);
    B[0] = Bv(0);
    B[1] = Bv(1);
    B[2] = Bv(2);

    return false;
}

bool Cyclotron::apply(const int &i, const double &t, Vector_t &E, Vector_t &B) {
    return apply(RefPartBunch_m->R[i], Vector_t(0.0), t, E, B);
}

bool Cyclotron::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    // ada

    const double rad = sqrt(R[0] * R[0] + R[1] * R[1]);
    const double xir = (rad - BP.rmin) / BP.delr;

    // ir : the mumber of path whoes radius is less then the 4 points of cell which surrond particle.
    // that is Rmin = 1900, dR = 20 , r = 1911, then ir = 0
    const int    ir = (int)xir;

    // wr1 : the relative distance to the inner path radius
    const double wr1 = xir - (double)ir;
    // wr2 : the relative distance to the outer path radius
    const double wr2 = 1.0 - wr1;

    double tet, tet_map, xit;
    int it;

    if((R[0] > 0) && (R[1] >= 0)) tet = atan(R[1] / R[0]);
    else if((R[0] < 0) && (R[1] >= 0)) tet = pi + atan(R[1] / R[0]);
    else if((R[0] < 0) && (R[1] <= 0)) tet = pi + atan(R[1] / R[0]);
    else if((R[0] > 0) && (R[1] <= 0)) tet = 2.0 * pi + atan(R[1] / R[0]);
    else if((R[0] == 0) && (R[1] > 0)) tet = pi / 2.0;
    else if((R[0] == 0) && (R[1] < 0)) tet = 3.0 / 2.0 * pi;

    double tet_rad = tet;

    // the actual angle of particle
    tet = tet / pi * 180.0;

    // the corresponding angle on the field map
    // Note: this does not work if the start point of field map does not equal zero.
    tet_map = fmod(tet, 360.0 / symmetry_m);

    xit = tet_map / BP.dtet;

    it = (int) xit;

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
          With t his we have B-field AND this is far more
          intuitive for me ....
        */
        r1t1 = idx(ir, it);
        r2t1 = idx(ir + 1, it);
        r1t2 = idx(ir, it + 1);
        r2t2 = idx(ir + 1, it + 1);
    }


    /* debug
     *gmsg << "x,y,z  = ( "<< R[0] << " , "<< R[1] << " , "<< R[2] << ")[mm] **********, "<<endl;
     *gmsg << "r2t2, r2t1, r1t2, r1t1 = (" << r2t2 <<" , "<<r2t1 <<" , "<< r1t2<<" , "<<r1t1 <<" ) " <<endl;
     *gmsg << "wr1, wt1= (" << wr1 <<" , "<<wt1 <<" ) " <<endl;
     */


    double bzcub = 0.0, bzf = 0.0, bz = 0.0;
    double brcub = 0.0, brf = 0.0, br = 0.0;
    double btcub = 0.0, btf = 0.0, bt = 0.0;

    if((it >= 0) && (ir >= 0) && (it < Bfield.ntetS) && (ir < Bfield.nrad)) {

        /* Bz */
        bzf = (Bfield.bfld[r1t1] * wr2 * wt2 + Bfield.bfld[r2t1] * wr1 * wt2 +
               Bfield.bfld[r1t2] * wr2 * wt1 + Bfield.bfld[r2t2] * wr1 * wt1);

        bzcub = (Bfield.f2[r1t1] * wr2 * wt2 +
                 Bfield.f2[r2t1] * wr1 * wt2 +
                 Bfield.f2[r1t2] * wr2 * wt1 +
                 Bfield.f2[r2t2] * wr1 * wt1) * pow(R[2], 2.0);

        // bz = -( bzf - bzcub );
        bz = - bzf ;


        /* Br */
        brf = (Bfield.dbr[r1t1] * wr2 * wt2 +
               Bfield.dbr[r2t1] * wr1 * wt2 +
               Bfield.dbr[r1t2] * wr2 * wt1 +
               Bfield.dbr[r2t2] * wr1 * wt1) * R[2];


        brcub = (Bfield.f3[r1t1] * wr2 * wt2 +
                 Bfield.f3[r2t1] * wr1 * wt2 +
                 Bfield.f3[r1t2] * wr2 * wt1 +
                 Bfield.f3[r2t2] * wr1 * wt1) * pow(R[2], 3.0);

        // br = -( brf - brcub );
        br = - brf;


        /* Btheta */
        btf = (Bfield.dbt[r1t1] * wr2 * wt2 +
               Bfield.dbt[r2t1] * wr1 * wt2 +
               Bfield.dbt[r1t2] * wr2 * wt1 +
               Bfield.dbt[r2t2] * wr1 * wt1) / rad * R[2];


        btcub = (Bfield.g3[r1t1] * wr2 * wt2 +
                 Bfield.g3[r2t1] * wr1 * wt2 +
                 Bfield.g3[r1t2] * wr2 * wt1 +
                 Bfield.g3[r2t2] * wr1 * wt1) / rad * pow(R[2], 3.0);

        // bt = -( btf - btcub );
        bt = - btf;

        /* Br Btheta -> Bx By */

        double partR = sqrt(R[0] * R[0] + R[1] * R[1]);
        double Amax1 = 1;
        double Amax2 = 3;
        double Amin = -2;
        double x01 = 4;
        double x02 = 8;
        double h1 = 0.03;
        double h2 = 0.2;
        double ftc = slptc_m;
        double part1 = pow(10.0, (partR / ftc - tcr1_m / ftc - x01) * h1);
        double part2 = pow(10.0, (x02 - partR / ftc + tcr1_m / ftc) * h2);
        double part3 = -(Amax1 - Amin) * h1 * log(10) / ftc / (1 + part1) / (1 + part1) * part1;
        double part4 = (Amax2 - Amin) * h2 * log(10) / ftc / (1 + part2) / (1 + part2) * part2;
        double dr = mbtc_m / 2.78 * (part3 + part4);
        double btr = mbtc_m / 2.78 * (Amin + (Amax1 - Amin) / (1 + part1) + (Amax2 - Amin) / (1 + part2) - 1.0);

        // trim coil
        double tca1 = 4;
        double tca2 = 22;
        // bz=bz*(1-btr/20);
        //for (int kk=0;kk<8;++kk) {
        //if  (tet>tca1+kk&&tet<tca2+kk) {

        if(partR < tcr2_m)       {
            bz -= btr;
            br -= dr * R[2];
        }



        //}
        //}

        B[0] = br * cos(tet_rad) - bt * sin(tet_rad);
        B[1] = br * sin(tet_rad) + bt * cos(tet_rad);
        B[2] = bz;

        /* Test for homo field.
           B(0) =0.0;
           B(1) =0.0;
           B(2) = -5.00000;// 5000Gauss

           *gmsg << "r2t2, r2t1, r1t2, r1t1 = (" << r2t2 <<" , "<<r2t1 <<" , "<< r1t2<<" , "<<r1t1 <<" ) b= " << B << endl;
           *gmsg << "wr1, wt1, wr2, wt2= (" << wr1 <<" , "<<wt1 <<" , "<< wr2 <<" , "<< wt2 << " ) " <<endl;
        */

    }

    /* error output, out of field! */
    else {
      ERRORMSG("Error!" << getName() << ".getFieldstrength(): out of boudaries (z,r) = (" << R[0] << "," << R[1] << "," << R[2] << ")" << endl
	       << " rad=" << rad << " [mm], theta=" << tet << " [deg], it=" << it << ", ir=" << ir << endl);
        return true;

    }
    return false;
}

void Cyclotron::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    RefPartBunch_m = bunch;
    online_m = true;
}

void Cyclotron::finalise() {
    online_m = false;
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

    assert(Bfield.dbr   = new double[Bfield.ntot]);
    assert(Bfield.dbrr  = new double[Bfield.ntot]);
    assert(Bfield.dbrrr = new double[Bfield.ntot]);

    assert(Bfield.dbrt  = new double[Bfield.ntot]);
    assert(Bfield.dbrrt = new double[Bfield.ntot]);
    assert(Bfield.dbrtt = new double[Bfield.ntot]);

    assert(Bfield.f2    = new double[Bfield.ntot]);
    assert(Bfield.f3    = new double[Bfield.ntot]);
    assert(Bfield.g3    = new double[Bfield.ntot]);

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
    /*
      debug

    for(int i = 0; i< Bfield.nrad; i++){
      for(int j = 0; j< Bfield.ntetS; j++){
	int index = idx(i,j);
	double x = i*BP.delr * sin(j);
	double y = i*BP.delr * cos(j);
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
    *gmsg << "*      (The first data block is useless)       " << endl;
    *gmsg << "* ----------------------------------------------" << endl;

    //  Bfield.filename = "s03av.nar"
    BP.Bfact = scaleFactor;

    if((f = fopen(fmapfn_m.c_str(), "r")) == NULL) {
        ERRORMSG("Error in Cyclotron::getFieldFromFile()!" << endl);
        ERRORMSG(" Cannot open file, please check if it really exists." << endl);
        exit(1);
    }

    assert(fscanf(f, "%lf", &BP.rmin));
    *gmsg << "* Minimal radius of measured field map: " << BP.rmin << " [mm]" << endl;

    assert(fscanf(f, "%lf", &BP.delr));
    *gmsg << "* Stepsize in radial direction: " << BP.delr << " [mm]" << endl;

    assert(fscanf(f, "%lf", &BP.tetmin));
    *gmsg << "* Minimal angle of measured field map: " << BP.tetmin << " [deg.]" << endl;

    assert(fscanf(f, "%lf", &BP.dtet));
    //if the value is nagtive, the actual value is its reciprocal.
    if(BP.dtet < 0.0) BP.dtet = 1.0 / (-BP.dtet);
    *gmsg << "* Stepsize in azimuth direction: " << BP.dtet << " [deg.]" << endl;

    for(int i = 0; i < 13; i++)assert(fscanf(f, "%s", fout));

    assert(fscanf(f, "%d", &Bfield.nrad));
    *gmsg << "* Index in radial direction: " << Bfield.nrad << endl;

    assert(fscanf(f, "%d", &Bfield.ntet));
    *gmsg << "* Index in azimuthal direction: " << Bfield.ntet << endl;

    Bfield.ntetS = Bfield.ntet + 1;
    *gmsg << "* Accordingly, total grid point along azimuth:  " << Bfield.ntetS << endl;

    for(int i = 0; i < 5; i++) {
        assert(fscanf(f, "%s", fout));
    }
    assert(fscanf(f, "%d", &lpar));
    // msg<< "READ"<<lpar<<" DATA ENTRIES"<<endl;

    for(int i = 0; i < 4; i++) {
        assert(fscanf(f, "%s", fout));
    }

    for(int i = 0; i < lpar; i++) {
        assert(fscanf(f, "%16lE", &dtmp));
    }
    for(int i = 0; i < 6; i++) {
        assert(fscanf(f, "%s", fout));
    }
    //*gmsg << "* READ FILE DESCRIPTION..." <<endl;
    for(int i = 0; i < 10000; i++) {
        assert(fscanf(f, "%s", fout));
        if(strcmp(fout, "LREC=") == 0)break;
    }

    for(int i = 0; i < 5; i++) {
        assert(fscanf(f, "%s", fout));
    }
    Bfield.ntot = idx(Bfield.nrad - 1, Bfield.ntet) + 1;
    //jjyang
    *gmsg << "* Total stored grid point number ( ntetS * nrad ) : " << Bfield.ntot << endl;
    assert(Bfield.bfld  = new double[Bfield.ntot]);
    assert(Bfield.dbt   = new double[Bfield.ntot]);
    assert(Bfield.dbtt  = new double[Bfield.ntot]);
    assert(Bfield.dbttt = new double[Bfield.ntot]);

    *gmsg << "* read -in loop one block per radius" << endl;
    *gmsg << "* rescaling of the fields with factor: " << BP.Bfact << endl;
    for(int i = 0; i < Bfield.nrad; i++) {

        if(i > 0) {
            for(int dummy = 0; dummy < 6; dummy++) {
                assert(fscanf(f, "%s", fout)); // INFO-LINE
            }
        }
        for(int k = 0; k < Bfield.ntet; k++) {
            assert(fscanf(f, "%16lE", &(Bfield.bfld[idx(i, k)])));
            Bfield.bfld[idx(i, k)] *= BP.Bfact;
        }
        for(int k = 0; k < Bfield.ntet; k++) {
            assert(fscanf(f, "%16lE", &(Bfield.dbt[idx(i, k)])));
            Bfield.dbt[idx(i, k)] *= BP.Bfact;
        }
        for(int k = 0; k < Bfield.ntet; k++) {
            assert(fscanf(f, "%16lE", &(Bfield.dbtt[idx(i, k)])));
            Bfield.dbtt[idx(i, k)] *= BP.Bfact;
        }
        for(int k = 0; k < Bfield.ntet; k++) {
            assert(fscanf(f, "%16lE", &(Bfield.dbttt[idx(i, k)])));
            Bfield.dbttt[idx(i, k)] *= BP.Bfact;
        }
    }
    fclose(f);


    *gmsg << "* Field Map read successfully!" << endl << endl;
}

// Calculates Radiae of initial grid.
// dimensions in [mm]!
void Cyclotron::initR(double rmin, double dr, int nrad) {
    assert(BP.rarr = new double[nrad]);
    for(int i = 0; i < nrad; i++) {
        BP.rarr[i] = rmin + i * dr;
    }
    BP.delr = dr;
}


void Cyclotron::initialise(PartBunch *bunch, const int &fieldflag, const double &scaleFactor) {
    RefPartBunch_m = bunch;
    //    PSIBF, AVFEQBF, ANSYSBF, FFAGBF
    // for your own format field, you should add your own getFieldFromFile() function by yourself.

    if(fieldflag == 1) {
        //*gmsg<<"Read field data from PSI format field map file."<<endl;
        myBFieldType_m = PSIBF;
        getFieldFromFile(scaleFactor);

    } else if(fieldflag == 2) {
        // *gmsg<<"Read data from 450MeV Carbon cyclotron field file"<<endl
        myBFieldType_m = PSIBF;
        getFieldFromFile_Carbon(scaleFactor);

    } else if(fieldflag == 3) {
        // *gmsg<<"Read data from 100MeV H- cyclotron CYCIAE-100 field file"<<endl;
        myBFieldType_m = ANSYSBF;
        getFieldFromFile_CYCIAE(scaleFactor);

    } else if(fieldflag == 4) {
        *gmsg << "* Read AVFEQ data (Riken) use bfield scale factor bs= " << getBScale() << endl;
        myBFieldType_m = AVFEQBF;
        getFieldFromFile_AVFEQ(scaleFactor);

    } else if(fieldflag == 5) {
        *gmsg << "Read FFAG data MSU/FNAL " << getBScale() << endl;
        myBFieldType_m = FFAGBF;
        getFieldFromFile_FFAG(scaleFactor);

    } else
        ERRORMSG("The field reading function of this TYPE of CYCLOTRON has not implemented yet.!" << endl);

    // calculate the radii of initial grid.
    initR(BP.rmin, BP.delr, Bfield.nrad);

    // calculate the remaining derivatives
    getdiffs();

    //*gmsg<<"----------------------------------------------"<<endl;

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

    assert(Bfield.bfld  = new double[Bfield.ntot]);
    assert(Bfield.dbt   = new double[Bfield.ntot]);
    assert(Bfield.dbtt  = new double[Bfield.ntot]);
    assert(Bfield.dbttt = new double[Bfield.ntot]);

    *gmsg << "* rescaling of the fields with factor: " << BP.Bfact << endl;

    fstream fp("gnu.out", ios::out);

    int count = 0;

    for(int r = 0; r < Bfield.nrad; r++) {
        for(int k = 0; k < Bfield.ntet; k++) {
            Bfield.bfld[idx(r, k)] = bzv[count] * BP.Bfact;
            fp << BP.rmin + (r * BP.delr) << " \t " << k*(BP.tetmin + BP.dtet) << " \t " << Bfield.bfld[idx(r, k)] << endl;
            count++;
        }
    }
    fp.close();
    *gmsg << "* Field Map read successfully nelem= " << count << endl << endl;
}

void Cyclotron::getFieldFromFile_AVFEQ(const double &scaleFactor) {

    FILE *f = NULL;
    int lpar;
    char fout[100];
    double dtmp;

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
        ERRORMSG("Error in Cyclotron::getFieldFromFile_AVFEQ()!" << endl);
        ERRORMSG(" Cannot open file, please check if it really exists." << endl);
        exit(1);
    }

    assert(fscanf(f, "%lf", &BP.rmin));
    *gmsg << "* Minimal radius of measured field map: " << BP.rmin << " [mm]" << endl;

    double rmax;
    assert(fscanf(f, "%lf", &rmax));
    *gmsg << "* Maximal radius of measured field map: " << rmax << " [mm]" << endl;

    assert(fscanf(f, "%lf", &BP.delr));
    *gmsg << "* Stepsize in radial direction: " << BP.delr << " [mm]" << endl;

    assert(fscanf(f, "%lf", &BP.tetmin));
    *gmsg << "* Minimal angle of measured field map: " << BP.tetmin << " [deg.]" << endl;

    double tetmax;
    assert(fscanf(f, "%lf", &tetmax));
    *gmsg << "* Maximal angle of measured field map: " << tetmax << " [deg.]" << endl;

    assert(fscanf(f, "%lf", &BP.dtet));
    //if the value is nagtive, the actual value is its reciprocal.

    if(BP.dtet < 0.0) BP.dtet = 1.0 / (-BP.dtet);
    *gmsg << "* Stepsize in azimuth direction: " << BP.dtet << " [deg.]" << endl;

    Bfield.ntetS = (int)((tetmax - BP.tetmin) / BP.dtet + 1);
    *gmsg << "* Total grid point along azimuth:  " << Bfield.ntetS << endl;

    Bfield.nrad = (int)(rmax - BP.rmin) / BP.delr;


    int ntotidx = idx(Bfield.nrad, Bfield.ntetS) + 1;

    Bfield.ntot = Bfield.ntetS * Bfield.nrad;
    *gmsg << "* Total stored grid point number ( ntetS * nrad ) : " << Bfield.ntot << " ntot-idx= " << ntotidx << endl;

    assert(Bfield.bfld  = new double[Bfield.ntot]);
    assert(Bfield.dbt   = new double[Bfield.ntot]);
    assert(Bfield.dbtt  = new double[Bfield.ntot]);
    assert(Bfield.dbttt = new double[Bfield.ntot]);

    *gmsg << "* rescaling of the fields with factor: " << BP.Bfact << endl;

    fstream fp("gnu.out", ios::out);

    double tmp;
    int count = 0;

    for(int r = 0; r < Bfield.nrad; r++) {
        assert(fscanf(f, "%16lE", &tmp));   // over read
        for(int k = 0; k < Bfield.ntetS; k++) {
            assert(fscanf(f, "%16lE", &(Bfield.bfld[idx(r, k)])));
            Bfield.bfld[idx(r, k)] *= BP.Bfact;

            fp << BP.rmin + (r * BP.delr) << " \t " << k*(BP.tetmin + BP.dtet) << " \t " << Bfield.bfld[idx(r, k)] << " idx= " << idx(r, k)  << endl;
            count++;
        }
    }
    fp.close();
    fclose(f);
    *gmsg << "* Field Map read successfully nelem= " << count << endl << endl;
}


// read field map from external file.
void Cyclotron::getFieldFromFile_Carbon(const double &scaleFactor) {

    FILE *f = NULL;
    int lpar;
    char fout[100];
    double dtmp;

    *gmsg << "* ----------------------------------------------" << endl;
    *gmsg << "*      READ IN CARBON CYCLOTRON FIELD MAP       " << endl;
    *gmsg << "* ----------------------------------------------" << endl;

    BP.Bfact = scaleFactor;

    if((f = fopen(fmapfn_m.c_str(), "r")) == NULL) {
        ERRORMSG("* Error in Cyclotron::getFieldFromFile_Carbon()!" << endl);
        ERRORMSG(" Cannot open file, please check if it really exists." << endl);
        exit(1);
    }

    assert(fscanf(f, "%lf", &BP.rmin));
    *gmsg << "* Minimal radius of measured field map: " << BP.rmin << " [mm]" << endl;

    assert(fscanf(f, "%lf", &BP.delr));
    *gmsg << "* Stepsize in radial direction: " << BP.delr << " [mm]" << endl;

    assert(fscanf(f, "%lf", &BP.tetmin));
    *gmsg << "* Minimal angle of measured field map: " << BP.tetmin << " [deg.]" << endl;

    assert(fscanf(f, "%lf", &BP.dtet));
    //if the value is nagtive, the actual value is its reciprocal.
    if(BP.dtet < 0.0) BP.dtet = 1.0 / (-BP.dtet);
    *gmsg << "* Stepsize in azimuth direction: " << BP.dtet << " [deg.]" << endl;

    assert(fscanf(f, "%d", &Bfield.ntet));
    *gmsg << "* Index in azimuthal direction: " << Bfield.ntet << endl;

    assert(fscanf(f, "%d", &Bfield.nrad));
    *gmsg << "* Index in radial direction: " << Bfield.nrad << endl;

    Bfield.ntetS = Bfield.ntet + 1;
    *gmsg << "* Accordingly, total grid point along azimuth:  " << Bfield.ntetS << endl;

    Bfield.ntot = idx(Bfield.nrad - 1, Bfield.ntet) + 1;

    *gmsg << "* Total stored grid point number ( ntetS * nrad ) : " << Bfield.ntot << endl;
    assert(Bfield.bfld  = new double[Bfield.ntot]);
    assert(Bfield.dbt   = new double[Bfield.ntot]);
    assert(Bfield.dbtt  = new double[Bfield.ntot]);
    assert(Bfield.dbttt = new double[Bfield.ntot]);

    *gmsg << "* rescaling of the fields with factor: " << BP.Bfact << endl;

    fstream fp("gnu.out", ios::out);

    for(int i = 0; i < Bfield.nrad; i++) {
        for(int k = 0; k < Bfield.ntet; k++) {
            assert(fscanf(f, "%16lE", &(Bfield.bfld[idx(i, k)])));
            Bfield.bfld[idx(i, k)] *= BP.Bfact;

            fp << BP.rmin + (i * BP.delr) << " \t " << k*(BP.tetmin + BP.dtet) << " \t " << Bfield.bfld[idx(i, k)] << endl;
        }
    }
    fclose(f);
    fp.close();
    *gmsg << "* Field Map read successfully!" << endl << endl;
}


// read field map from external file.
void Cyclotron::getFieldFromFile_CYCIAE(const double &scaleFactor) {

    FILE *f = NULL;
    int lpar;
    char fout[100];
    double dtmp;

    *gmsg << "* ----------------------------------------------" << endl;
    *gmsg << "*    READ IN CYCIAE-100 CYCLOTRON FIELD MAP     " << endl;
    *gmsg << "* ----------------------------------------------" << endl;

    BP.Bfact = scaleFactor;

    if((f = fopen(fmapfn_m.c_str(), "r")) == NULL) {
        ERRORMSG("* Error in Cyclotron::getFieldFromFile_Carbon()!" << endl);
        ERRORMSG(" Cannot open file, please check if it really exists." << endl);
        exit(1);
    }

    assert(fscanf(f, "%lf", &BP.rmin));
    *gmsg << "* Minimal radius of measured field map: " << BP.rmin << " [mm]" << endl;

    assert(fscanf(f, "%lf", &BP.delr));
    *gmsg << "* Stepsize in radial direction: " << BP.delr << " [mm]" << endl;

    assert(fscanf(f, "%lf", &BP.tetmin));
    *gmsg << "* Minimal angle of measured field map: " << BP.tetmin << " [deg.]" << endl;

    assert(fscanf(f, "%lf", &BP.dtet));
    //if the value is nagtive, the actual value is its reciprocal.
    if(BP.dtet < 0.0) BP.dtet = 1.0 / (-BP.dtet);
    *gmsg << "* Stepsize in azimuth direction: " << BP.dtet << " [deg.]" << endl;

    assert(fscanf(f, "%d", &Bfield.ntet));
    *gmsg << "* Index in azimuthal direction: " << Bfield.ntet << endl;

    assert(fscanf(f, "%d", &Bfield.nrad));
    *gmsg << "* Index in radial direction: " << Bfield.nrad << endl;

    Bfield.ntetS = Bfield.ntet + 1;
    *gmsg << "* Accordingly, total grid point along azimuth:  " << Bfield.ntetS << endl;

    Bfield.ntot = idx(Bfield.nrad - 1, Bfield.ntet) + 1;

    *gmsg << "* Total stored grid point number ( ntetS * nrad ) : " << Bfield.ntot << endl;
    assert(Bfield.bfld  = new double[Bfield.ntot]);
    assert(Bfield.dbt   = new double[Bfield.ntot]);
    assert(Bfield.dbtt  = new double[Bfield.ntot]);
    assert(Bfield.dbttt = new double[Bfield.ntot]);

    *gmsg << "* rescaling of the fields with factor: " << BP.Bfact << endl;

    int nHalfPoints = Bfield.ntet / 2.0 + 1;

    for(int i = 0; i < Bfield.nrad; i++) {

        for(int ii = 0; ii < 13; ii++)assert(fscanf(f, "%s", fout));
        for(int k = 0; k < nHalfPoints; k++) {
            assert(fscanf(f, "%d", &dtmp));
            assert(fscanf(f, "%d", &dtmp));
            assert(fscanf(f, "%d", &dtmp));
            assert(fscanf(f, "%lf", &(Bfield.bfld[idx(i, k)])));
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

void Cyclotron::getDimensions(double &zBegin, double &zEnd) const
{ }

