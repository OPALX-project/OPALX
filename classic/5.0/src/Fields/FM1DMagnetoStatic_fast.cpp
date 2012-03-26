#include <fstream>
#include <ios>

#include "Fields/FM1DMagnetoStatic_fast.hh"
#include "Fields/Fieldmap.icc"
#include "Physics/Physics.h"
#include "gsl/gsl_fft_real.h"

using namespace std;
using Physics::two_pi;

FM1DMagnetoStatic_fast::FM1DMagnetoStatic_fast(string aFilename):
    Fieldmap(aFilename) {
    int tmpInt;
    string tmpString;
    double tmpDouble;
    ifstream file;

    onAxisField_m = NULL;

    Type = T1DMagnetoStatic;

    // open field map, parse it and disable element on error
    file.open(Filename_m.c_str());
    if(file.good()) {
        bool parsing_passed = interpreteLine<string, int>(file, tmpString, tmpInt);
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
        hz_m = (zend_m - zbegin_m) / (num_gridpz_m - 1);
        length_m = 2.0 * num_gridpz_m * hz_m;
    } else {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m = -1e-3;
    }
}

FM1DMagnetoStatic_fast::~FM1DMagnetoStatic_fast() {
    if(onAxisField_m != NULL) {
        for(int i = 0; i < 4; ++i) {
            gsl_spline_free(onAxisInterpolants_m[i]);
            gsl_interp_accel_free(onAxisAccel_m[i]);
        }
        delete[] onAxisField_m;
    }
}

void FM1DMagnetoStatic_fast::readMap() {
    if(onAxisField_m == NULL) {
        // declare variables and allocate memory
        ifstream in;

        string tmpString;

        int tmpInt;
        int accuracy;

        double tmpDouble;
        double Bz_max = 0.0;
        double z = 0.0;
        double interior_derivative, base;
        double coskzl, sinkzl;

        double *RealValues = new double[2*num_gridpz_m];
        onAxisField_m = new double[num_gridpz_m];
        double *higherDerivatives[3];
        higherDerivatives[0] = new double[num_gridpz_m];
        higherDerivatives[1] = new double[num_gridpz_m];
        higherDerivatives[2] = new double[num_gridpz_m];
        double *zvals = new double[num_gridpz_m];

        gsl_fft_real_wavetable *real = gsl_fft_real_wavetable_alloc(2 * num_gridpz_m);
        gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(2 * num_gridpz_m);

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
            // mirror the field map to make sure that it is periodic
            RealValues[num_gridpz_m - 1 - i] = RealValues[num_gridpz_m + i];
            onAxisField_m[i] = RealValues[num_gridpz_m + i];
        }
        in.close();

        // apply FFT
        gsl_fft_real_transform(RealValues, 1, 2 * num_gridpz_m, real, work);

        // normalize the fourier coeficients such that [H] = A/m after interpolation
        RealValues[0] /= (Bz_max * 2. * num_gridpz_m);
        onAxisField_m[0] /= Bz_max;
        for(int i = 1; i < accuracy; i++) {
            RealValues[i] /= (Bz_max * num_gridpz_m);
            onAxisField_m[i] /= Bz_max;
        }
        for(int i = accuracy; i < 2 * accuracy - 1; i++) {
            RealValues[i] /= (Bz_max * num_gridpz_m);
        }

        for(int i = 0; i < num_gridpz_m; ++ i, z += hz_m) {
            const double kz = two_pi * (z / length_m + 0.5);
            higherDerivatives[0][i] = 0.0;
            higherDerivatives[1][i] = 0.0;
            higherDerivatives[2][i] = 0.0;
            int n = 1;
            zvals[i] = z;
            for(int l = 1; l < accuracy; ++l, n += 2) {
                base = two_pi / length_m * l;
                interior_derivative = base;
                coskzl = cos(kz * l);
                sinkzl = sin(kz * l);
                
                higherDerivatives[0][i] += interior_derivative * (-RealValues[n] * sinkzl - RealValues[n+1] * coskzl);
                interior_derivative *= base;
                higherDerivatives[1][i] += interior_derivative * (-RealValues[n] * coskzl + RealValues[n+1] * sinkzl);
                interior_derivative *= base;
                higherDerivatives[2][i] += interior_derivative * (RealValues[n] * sinkzl + RealValues[n+1] * coskzl);
            }
        }

        gsl_fft_real_workspace_free(work);
        gsl_fft_real_wavetable_free(real);

        delete[] RealValues;

        onAxisInterpolants_m[0] = gsl_spline_alloc(gsl_interp_cspline, num_gridpz_m);
        gsl_spline_init(onAxisInterpolants_m[0], zvals, onAxisField_m, num_gridpz_m);
        onAxisAccel_m[0] = gsl_interp_accel_alloc();
        for(int i = 1; i < 4; ++i) {
            onAxisInterpolants_m[i] = gsl_spline_alloc(gsl_interp_cspline, num_gridpz_m);
            gsl_spline_init(onAxisInterpolants_m[i], zvals, higherDerivatives[i-1], num_gridpz_m);
            onAxisAccel_m[i] = gsl_interp_accel_alloc();
        }

        delete[] higherDerivatives[0];
        delete[] higherDerivatives[1];
        delete[] higherDerivatives[2];
        delete[] zvals;

        INFOMSG(typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << endl);
    }
}

void FM1DMagnetoStatic_fast::freeMap() {
    if(onAxisField_m != NULL) {
        for(int i = 0; i < 4; ++i) {
            gsl_spline_free(onAxisInterpolants_m[i]);
            gsl_interp_accel_free(onAxisAccel_m[i]);
        }

        delete[] onAxisField_m;

        INFOMSG(typeset_msg("freed fieldmap '" + Filename_m  + "'", "info") << endl);
    }
}

bool FM1DMagnetoStatic_fast::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const {
    // do fourier interpolation for on-axis value
    const double RR2 = R(0) * R(0) + R(1) * R(1);

    double bz = gsl_spline_eval(onAxisInterpolants_m[0], R(2), onAxisAccel_m[0]);
    double bzp = gsl_spline_eval(onAxisInterpolants_m[1], R(2), onAxisAccel_m[1]);
    double bzpp = gsl_spline_eval(onAxisInterpolants_m[2], R(2), onAxisAccel_m[2]);
    double bzppp = gsl_spline_eval(onAxisInterpolants_m[3], R(2), onAxisAccel_m[3]);

    // expand to off-axis
    const double BfieldR = -bzp / 2. + bzppp / 16. * RR2;

    B(0) += BfieldR * R(0);
    B(1) += BfieldR * R(1);
    B(2) += bz - bzpp * RR2 / 4.;
    return false;
}

bool FM1DMagnetoStatic_fast::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const {
    double bzp = gsl_spline_eval(onAxisInterpolants_m[1], R(2), onAxisAccel_m[1]);

    B(2) += bzp;
    return false;
}

void FM1DMagnetoStatic_fast::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = rbegin_m;
    rEnd = rend_m;
}

void FM1DMagnetoStatic_fast::getFieldDimensions(double &xIni, double &xFinal, double &yIni, double &yFinal, double &zIni, double &zFinal) const {}

void FM1DMagnetoStatic_fast::swap()
{ }

void FM1DMagnetoStatic_fast::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (1D magnetostatic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double FM1DMagnetoStatic_fast::getFrequency() const {
    return 0.0;
}

void FM1DMagnetoStatic_fast::setFrequency(double freq)
{ }
