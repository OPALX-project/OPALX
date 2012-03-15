#include <fstream>
#include <ios>

#include "Fields/FM1DMagnetoStatic_fast.hh"
#include "Fields/Fieldmap.icc"
#include "Physics/Physics.h"
#include "gsl/gsl_fft_real.h"

using namespace std;
using Physics::mu_0;
using Physics::c;
using Physics::two_pi;
using Physics::pi;

FM1DMagnetoStatic_fast::FM1DMagnetoStatic_fast(string aFilename)
    : Fieldmap(aFilename),
      FieldstrengthBz_m(NULL),
      FieldstrengthBr_m(NULL) {
    Inform msg("FM1DMS ");
    ifstream file;
    int tmpInt;
    string tmpString;
    double tmpDouble;

    Type = T1DMagnetoStatic;

    file.open(Filename_m.c_str());
    if(file.good()) {
        bool parsing_passed = interpreteLine<string, int>(file, tmpString, tmpInt);
        parsing_passed = parsing_passed &&
                         interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
        parsing_passed = parsing_passed &&
                         interpreteLine<double, double, int>(file, rbegin_m, rend_m, num_gridpr_m);
        for(int i = 0; (i <= num_gridpz_m) && parsing_passed; ++ i) {
            parsing_passed = parsing_passed &&
                             interpreteLine<double>(file, tmpDouble);
        }

        parsing_passed = parsing_passed &&
                         interpreteEOF(file);

        file.close();
        lines_read_m = 0;

        if(!parsing_passed) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
        } else {
            // conversion cm to m
            rbegin_m /= 100.;
            rend_m /= 100.;
            zbegin_m /= 100.;
            zend_m /= 100.;

            hr_m = (rend_m - rbegin_m) / num_gridpr_m;
            hz_m = (zend_m - zbegin_m) / num_gridpz_m;

            // num spacings -> num grid points
            num_gridpr_m++;
            num_gridpz_m++;
        }
    } else {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m = -1e-3;
    }
}


FM1DMagnetoStatic_fast::~FM1DMagnetoStatic_fast() {
    if(FieldstrengthBz_m != NULL) {
        delete[] FieldstrengthBz_m;
        delete[] FieldstrengthBr_m;
    }
}

void FM1DMagnetoStatic_fast::readMap() {
    if(FieldstrengthBz_m == NULL) {
        Inform msg("FM1DMS ");
        ifstream in;

        string tmpString;

        int tmpInt;
        int accuracy;

        double tmpDouble;
        double ez, ezp, ezpp, ezppp, kz;
        double somefactor, somefactor_base;
        double R;
        double Bz_max = 0.0;

        double *RealValues = new double[2*num_gridpz_m];

        gsl_fft_real_wavetable *real = gsl_fft_real_wavetable_alloc(2 * num_gridpz_m);
        gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(2 * num_gridpz_m);

        FieldstrengthBz_m = new double[num_gridpz_m * num_gridpr_m];
        FieldstrengthBr_m = new double[num_gridpz_m * num_gridpr_m];

        // read in field map and parse it
        in.open(Filename_m.c_str());
        interpreteLine<string, int>(in, tmpString, accuracy);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);

        for(int i = 0; i < num_gridpz_m; i++) {
            interpreteLine<double>(in, RealValues[num_gridpz_m + i]);
            if(fabs(RealValues[num_gridpz_m + i]) > Bz_max) {
                Bz_max = fabs(RealValues[num_gridpz_m + i]);
            }
            RealValues[num_gridpz_m - 1 - i] = RealValues[num_gridpz_m + i];
        }
        in.close();

        // apply FFT
        gsl_fft_real_transform(RealValues, 1, 2 * num_gridpz_m, real, work);

        // normalize the fourier coeficients such that [H] = A/m after interpolation
        RealValues[0] /= (Bz_max * 2 * num_gridpz_m);
        for(int i = 1; i < 2 * accuracy - 1; i++) {
            RealValues[i] /= (Bz_max * num_gridpz_m);
        }

        // do fourier interpolation for on-axis value
        int n = 1;
        for(int i = 0; i < num_gridpz_m; i++) {
            kz = (pi * i) / num_gridpz_m + pi;
            ez = RealValues[0];
            ezp = 0.0;
            ezpp = 0.0;
            ezppp = 0.0;
            for(int l = 1; l < accuracy ; l++, n += 2) {
                somefactor_base = pi * l / (num_gridpz_m * hz_m);
                somefactor = 1.0;
                ez    += (RealValues[n] * cos(kz * l) - RealValues[n+1] * sin(kz * l));
                somefactor *= somefactor_base;
                ezp   += somefactor * (-RealValues[n] * sin(kz * l) - RealValues[n+1] * cos(kz * l));
                somefactor *= somefactor_base;
                ezpp  += somefactor * (-RealValues[n] * cos(kz * l) + RealValues[n+1] * sin(kz * l));
                somefactor *= somefactor_base;
                ezppp += somefactor * (RealValues[n] * sin(kz * l) + RealValues[n+1] * cos(kz * l));
            }
            // expand off-axis
            R = 0.0;
            for(int j = 0; j < num_gridpr_m; j++) {
                FieldstrengthBz_m[i + j * num_gridpz_m] =  ez - ezpp * R * R / 4.;
                FieldstrengthBr_m[i + j * num_gridpz_m] = -ezp / 2. + ezppp * R * R / 16.;
                R = hr_m * (j + 1);
            }
        }

        gsl_fft_real_workspace_free(work);
        gsl_fft_real_wavetable_free(real);

        delete[] RealValues;

        msg << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << "\n"
            << endl;
    }
}

void FM1DMagnetoStatic_fast::freeMap() {
    if(FieldstrengthBz_m != NULL) {
        Inform msg("FM1DD ");

        delete[] FieldstrengthBz_m;
        delete[] FieldstrengthBr_m;

        msg << typeset_msg("freed fieldmap '" + Filename_m  + "'", "info") << "\n"
            << endl;

    }
}

bool FM1DMagnetoStatic_fast::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const {
    // do bi-linear interpolation
    const double RR = sqrt(R(0) * R(0) + R(1) * R(1));

    const int indexr = (int)floor(RR / hr_m);
    const double leverr = (RR / hr_m) - indexr;

    const int indexz = (int)floor((R(2)) / hz_m);
    const double leverz = (R(2) / hz_m) - indexz;

    if((indexz < 0) || (indexz + 2 > num_gridpz_m) || (indexr < 0) || (indexr + 2 > num_gridpr_m)) {
        return true;
    }

    const int index1 = indexz + indexr * num_gridpz_m;
    const int index2 = index1 + num_gridpz_m;
    const double BfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBr_m[index1]
                           + leverz         * (1.0 - leverr) * FieldstrengthBr_m[index1 + 1]
                           + (1.0 - leverz) * leverr         * FieldstrengthBr_m[index2]
                           + leverz         * leverr         * FieldstrengthBr_m[index2 + 1];

    const double BfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBz_m[index1]
                           + leverz         * (1.0 - leverr) * FieldstrengthBz_m[index1 + 1]
                           + (1.0 - leverz) * leverr         * FieldstrengthBz_m[index2]
                           + leverz         * leverr         * FieldstrengthBz_m[index2 + 1];

    B(0) += BfieldR * R(0);
    B(1) += BfieldR * R(1);
    B(2) += BfieldZ;
    return false;
}

bool FM1DMagnetoStatic_fast::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const
{ }


void FM1DMagnetoStatic_fast::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = rbegin_m;
    rEnd = rend_m;
}

void FM1DMagnetoStatic_fast::swap()
{ }

void FM1DMagnetoStatic_fast::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (1D magnetostatic fast); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double FM1DMagnetoStatic_fast::getFrequency() const {
    return 0.0;
}

void FM1DMagnetoStatic_fast::setFrequency(double freq)
{ }
