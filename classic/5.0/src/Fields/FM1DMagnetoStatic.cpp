#include <fstream>
#include <ios>

#include "Fields/FM1DMagnetoStatic.hh"
#include "Fields/Fieldmap.icc"
#include "Physics/Physics.h"
#include "gsl/gsl_fft_real.h"

using namespace std;
using Physics::two_pi;

FM1DMagnetoStatic::FM1DMagnetoStatic(string aFilename)
    : Fieldmap(aFilename),
      FourCoefs_m(NULL) {
    int tmpInt;
    string tmpString;
    double tmpDouble;
    ifstream file;

    Type = T1DMagnetoStatic;

    // open field map, parse it and disable element on error
    file.open(Filename_m.c_str());
    if(file.good()) {
        bool parsing_passed = interpreteLine<string, int>(file, tmpString, accuracy_m);
        parsing_passed = parsing_passed &&
                         interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
        parsing_passed = parsing_passed &&
                         interpreteLine<double, double, int>(file, rbegin_m, rend_m, tmpInt);
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

            // num spacings -> num grid points
            num_gridpz_m++;
        }
        length_m = 2.0 * num_gridpz_m * (zend_m - zbegin_m) / (num_gridpz_m - 1);
    } else {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m = -1e-3;
    }
}

FM1DMagnetoStatic::~FM1DMagnetoStatic() {
    if(FourCoefs_m != NULL) {
        delete[] FourCoefs_m;
    }
}

void FM1DMagnetoStatic::readMap() {
    if(FourCoefs_m == NULL) {
        // declare variables and allocate memory
        ifstream in;

        string tmpString;

        int tmpInt;

        double tmpDouble;
        double Bz_max = 0.0;

        double *RealValues = new double[2*num_gridpz_m];

        gsl_fft_real_wavetable *real = gsl_fft_real_wavetable_alloc(2 * num_gridpz_m);
        gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(2 * num_gridpz_m);

        FourCoefs_m = new double[2*accuracy_m - 1];

        // read in field map and parse it
        in.open(Filename_m.c_str());
        interpreteLine<string, int>(in, tmpString, accuracy_m);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);

        for(int i = 0; i < num_gridpz_m; i++) {
            interpreteLine<double>(in, RealValues[num_gridpz_m + i]);
            if(fabs(RealValues[num_gridpz_m + i]) > Bz_max) {
                Bz_max = fabs(RealValues[num_gridpz_m + i]);
            }
            // mirror the field map to make sure that it is periodic
            RealValues[num_gridpz_m - 1 - i] = RealValues[num_gridpz_m + i];
        }
        in.close();

        // apply FFT
        gsl_fft_real_transform(RealValues, 1, 2 * num_gridpz_m, real, work);

        // normalize the fourier coeficients such that [H] = A/m after interpolation
        FourCoefs_m[0] = RealValues[0] / (Bz_max * 2. * num_gridpz_m);
        for(int i = 1; i < 2 * accuracy_m - 1; i++) {
            FourCoefs_m[i] = RealValues[i] / (Bz_max * num_gridpz_m);
        }

        gsl_fft_real_workspace_free(work);
        gsl_fft_real_wavetable_free(real);

        delete[] RealValues;

        INFOMSG(typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << endl);
    }
}

void FM1DMagnetoStatic::freeMap() {
    if(FourCoefs_m != NULL) {

        delete[] FourCoefs_m;

        INFOMSG(typeset_msg("freed fieldmap '" + Filename_m  + "'", "info") << endl);

    }
}

bool FM1DMagnetoStatic::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const {
    // do fourier interpolation for on-axis value
    const double RR2 = R(0) * R(0) + R(1) * R(1);

    const double kz = two_pi * R(2) / length_m + Physics::pi;

    double ez = FourCoefs_m[0];
    double ezp = 0.0;
    double ezpp = 0.0;
    double ezppp = 0.0;
    double somefactor_base, somefactor;
    double coskzl;
    double sinkzl;

    int n = 1;
    for(int l = 1; l < accuracy_m ; l++, n += 2) {
        somefactor_base =  two_pi / length_m * l;
        somefactor = 1.0;
        coskzl = cos(kz * l);
        sinkzl = sin(kz * l);
        ez    += (FourCoefs_m[n] * coskzl - FourCoefs_m[n+1] * sinkzl);
        somefactor *= somefactor_base;
        ezp   += somefactor * (-FourCoefs_m[n] * sinkzl - FourCoefs_m[n+1] * coskzl);
        somefactor *= somefactor_base;
        ezpp  += somefactor * (-FourCoefs_m[n] * coskzl + FourCoefs_m[n+1] * sinkzl);
        somefactor *= somefactor_base;
        ezppp += somefactor * (FourCoefs_m[n] * sinkzl + FourCoefs_m[n+1] * coskzl);
    }

    // expand to off-axis
    const double BfieldR = -ezp / 2. + ezppp / 16. * RR2;

    B(0) += BfieldR * R(0);
    B(1) += BfieldR * R(1);
    B(2) += ez - ezpp * RR2 / 4.;
    return false;
}

bool FM1DMagnetoStatic::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const {

    const double kz = two_pi * R(2) / length_m + Physics::pi;
    double ezp = 0.0;
    int n = 1;
    for(int l = 1; l < accuracy_m ; l++, n += 2)
        ezp += two_pi / length_m * l * (-FourCoefs_m[n] * sin(kz * l) - FourCoefs_m[n+1] * cos(kz * l));

    B(2) += ezp;
    return false;
}

void FM1DMagnetoStatic::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = rbegin_m;
    rEnd = rend_m;
}
void FM1DMagnetoStatic::getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const {}
void FM1DMagnetoStatic::swap()
{ }

void FM1DMagnetoStatic::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (1D magnetostatic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double FM1DMagnetoStatic::getFrequency() const {
    return 0.0;
}

void FM1DMagnetoStatic::setFrequency(double freq)
{ }
