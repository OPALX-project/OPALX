#include "CycMagneticField.h"

extern Inform *gmsg;

#define CHECK_CYC_FSCANF_EOF(arg) if(arg == EOF)\
throw GeneralClassicException("Cyclotron::getFieldFromFile",\
                              "fscanf returned EOF at " #arg);

CycMagneticField::CycMagneticField(const std::string fmapfn,
                                   const double& symmetry) :
    fmapfn_m(fmapfn),
    symmetry_m(symmetry)
{ }

void CycMagneticField::read(const double& scaleFactor) {
    
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

    //Bfield.ntetS = Bfield.ntet;
    Bfield.ntetS = Bfield.ntet + 1;
    *gmsg << "* Accordingly, total grid point along azimuth:  " << Bfield.ntetS << endl;

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
    
    std::fstream fp1;
    fp1.open("gnu.out", std::ios::out);
    for (int i = 0; i < Bfield.nrad; ++i) {
        for (int k = 0; k < Bfield.ntet; ++k) {
            fp1 << BP.rmin + (i * BP.delr) << " \t " << k * (BP.tetmin + BP.dtet)
                << " \t " << Bfield.bfld[idx(i, k)] << std::endl;
        }
    }
    fp1.close();
    
    fclose(f);

    *gmsg << "* Field Maps read successfully!" << endl << endl;
    
    // calculate the radii of initial grid.
    initR_m(BP.rmin, BP.delr, Bfield.nrad);

    // calculate the remaining derivatives
    getdiffs_m();
}


void CycMagneticField::interpolate(double& bint,
                                   double& brint,
                                   double& btint,
                                   const double& rad,
                                   const double& theta
                                  )
{
    // x horizontal
    // y longitudinal
    // z is vertical
    
//     std::cout << "rad = " << rad << " theta = " << theta << std::endl;
    const double xir = (rad * 1000 - BP.rmin /** 0.001*/) / (BP.delr /** 0.001*/);

    // ir : the mumber of path whoes radius is less then the 4 points of cell which surrond particle.
    const int    ir = (int)xir;

    // wr1 : the relative distance to the inner path radius
    const double wr1 = xir - (double)ir;
    // wr2 : the relative distance to the outer path radius
    const double wr2 = 1.0 - wr1;
    

    double tet_rad = theta;

    // the actual angle of particle
//     tet_rad = theta / Physics::pi * 180.0;
    
    // the corresponding angle on the field map
    // Note: this does not work if the start point of field map does not equal zero.
    double tet_map = fmod(tet_rad, 360.0 / symmetry_m);

    double xit = tet_map / BP.dtet;

    int it = (int) xit;

//        *gmsg << " tet_map= " << tet_map << " ir= " << ir << " it= " << it << " bf= " << Bfield.bfld[idx(ir,it)] << endl;

    const double wt1 = xit - (double)it;
    const double wt2 = 1.0 - wt1;

    // it : the number of point on the inner path whoes angle is less then the particle' corresponding angle.
    // include zero degree point
    it = it + 1;

    int r1t1, r2t1, r1t2, r2t2;
//     int ntetS = Bfield.ntet + 1;

    // r1t1 : the index of the "min angle, min radius" point in the 2D field array.
    // considering  the array start with index of zero, minus 1.

//     if(myBFieldType_m != FFAGBF) {
        /*
          For FFAG this does not work
        */
//         r1t1 = it + ntetS * ir - 1;
//         r1t2 = r1t1 + 1;
//         r2t1 = r1t1 + ntetS;
//         r2t2 = r2t1 + 1 ;

//     } else {
//         /*
//           With this we have B-field AND this is far more
//           intuitive for me ....
//         */
        r1t1 = idx(ir, it);
        r2t1 = idx(ir + 1, it);
        r1t2 = idx(ir, it + 1);
        r2t2 = idx(ir + 1, it + 1);
//     }

    bint = 0.0;
    brint = 0.0;
    btint = 0.0;

    if((it >= 0) && (ir >= 0) && (it < Bfield.ntetS) && (ir < Bfield.nrad)) {
        
        // dB_{z}/dr
        brint = (Bfield.dbr[r1t1] * wr2 * wt2 +
                Bfield.dbr[r2t1] * wr1 * wt2 +
                Bfield.dbr[r1t2] * wr2 * wt1 +
                Bfield.dbr[r2t2] * wr1 * wt1) * 1000.0;
        
        // dB_{z}/dtheta
        btint = Bfield.dbt[r1t1] * wr2 * wt2 +
                Bfield.dbt[r2t1] * wr1 * wt2 +
                Bfield.dbt[r1t2] * wr2 * wt1 +
                Bfield.dbt[r2t2] * wr1 * wt1;
        
        // B_{z}
        bint = Bfield.bfld[r1t1] * wr2 * wt2 +
               Bfield.bfld[r2t1] * wr1 * wt2 +
               Bfield.bfld[r1t2] * wr2 * wt1 +
               Bfield.bfld[r2t2] * wr1 * wt1;
    }
}


void CycMagneticField::initR_m(double rmin, double dr, int nrad) {
    BP.rarr.resize(nrad);
    for(int i = 0; i < nrad; i++) {
        BP.rarr[i] = rmin + i * dr;
    }
    BP.delr = dr;
}

// calculate derivatives with 5-point lagrange's formula.
double CycMagneticField::gutdf5d_m(double *f, double dx, const int kor, const int krl, const int lpr) {
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



void CycMagneticField::getdiffs_m() {

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

            double dtheta = Physics::pi / 180.0 * BP.dtet;

            int kEdge;

            kEdge = std::max(k - 2, 0);
            kEdge = std::min(kEdge, Bfield.ntet - 5);

            int dkFromEdge = k - kEdge;
            int index = idx(i, k);
            int indexkEdge = idx(i, kEdge);


            Bfield.dbt[index]    = gutdf5d_m(&Bfield.bfld[indexkEdge], dtheta, 0, dkFromEdge, 1);
            Bfield.dbtt[index]   = gutdf5d_m(&Bfield.bfld[indexkEdge], dtheta, 1, dkFromEdge, 1);
            Bfield.dbttt[index]  = gutdf5d_m(&Bfield.bfld[indexkEdge], dtheta, 2, dkFromEdge, 1);
        }
    }



    for(int k = 0; k < Bfield.ntet; k++) {
        // inner loop varies R
        for(int i = 0; i < Bfield.nrad; i++) {
            double rac = BP.rarr[i];
            // define iredg, the reference index for radial interpolation
            // standard: i-2 minimal: 0 (not negative!)  maximal: nrad-4
            int iredg = std::max(i - 2, 0);
            iredg = std::min(iredg, Bfield.nrad - 5);
            int irtak = i - iredg;
            int index = idx(i, k);
            int indexredg = idx(iredg, k);


            Bfield.dbr[index]    = gutdf5d_m(&Bfield.bfld[indexredg], BP.delr, 0, irtak, Bfield.ntetS);
            Bfield.dbrr[index]   = gutdf5d_m(&Bfield.bfld[indexredg], BP.delr, 1, irtak, Bfield.ntetS);
            Bfield.dbrrr[index]  = gutdf5d_m(&Bfield.bfld[indexredg], BP.delr, 2, irtak, Bfield.ntetS);

            Bfield.dbrt[index]   = gutdf5d_m(&Bfield.dbt[indexredg], BP.delr, 0, irtak, Bfield.ntetS);
            Bfield.dbrrt[index]  = gutdf5d_m(&Bfield.dbt[indexredg], BP.delr, 1, irtak, Bfield.ntetS);
            Bfield.dbrtt[index]  = gutdf5d_m(&Bfield.dbtt[indexredg], BP.delr, 0, irtak, Bfield.ntetS);

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
}