#include <fstream>
#include <ios>


#include "Fields/Astra1DMagnetoStatic_fast.hh"
#include "Fields/Fieldmap.icc"
#include "Physics/Physics.h"
#include "gsl/gsl_fft_real.h"

using namespace std;
using Physics::mu_0;
using Physics::c;
using Physics::two_pi;

Astra1DMagnetoStatic_fast::Astra1DMagnetoStatic_fast(string aFilename):
    Fieldmap(aFilename){
    ifstream file;
    int tmpInt;
    int skippedValues = 0;
    string tmpString;
    double tmpDouble;
    double tmpDouble2;

    onAxisField_m = NULL;

    Type = TAstraMagnetoStatic;

    // open field map, parse it and disable element on error
    file.open(Filename_m.c_str());
    if(file.good()) {
        bool parsing_passed = interpreteLine<string, int>(file, tmpString, tmpInt);
        parsing_passed = parsing_passed &&
                         interpreteLine<double, double>(file, zbegin_m, tmpDouble);

        tmpDouble2 = zbegin_m;
        while(!file.eof() && parsing_passed) {
            parsing_passed = interpreteLine<double, double>(file, zend_m, tmpDouble, false);
            if(zend_m - tmpDouble2 > 1e-10) {
                tmpDouble2 = zend_m;
            } else if(parsing_passed) {
                ++ skippedValues;
            }
        }

        file.close();
        num_gridpz_m = lines_read_m - 2 - skippedValues;
        lines_read_m = 0;

        if(!parsing_passed && !file.eof()) {
            disableFieldmapWarning();
            zend_m = zbegin_m - 1e-3;
        }
        hz_m = (zend_m - zbegin_m) / (num_gridpz_m - 1);
        length_m = 2.0 * num_gridpz_m * hz_m;
    } else {
        noFieldmapWarning();
        zbegin_m = 0.0;
        zend_m = -1e-3;
    }
}

Astra1DMagnetoStatic_fast::~Astra1DMagnetoStatic_fast() {
    if(onAxisField_m != NULL) {
        for(int i = 0; i < 4; ++i) {
            gsl_spline_free(onAxisInterpolants_m[i]);
            gsl_interp_accel_free(onAxisAccel_m[i]);
        }
        delete[] onAxisField_m;
    }
}

void Astra1DMagnetoStatic_fast::readMap() {
    if(onAxisField_m == NULL) {
        // declare variables and allocate memory
        ifstream in;

        bool parsing_passed = true;
        int accuracy;
        int nsp;
        int ii;

        string tmpString;

        double Bz_max = 0.0;
        double *higherDerivatives[3];
        double *zvals[2];
        double tmpDouble = zbegin_m - hz_m;
        // double z = 0.0;
        double interior_derivative, base;
        double coskzl, sinkzl;

        double *RealValues = new double[2*num_gridpz_m];
        onAxisField_m = new double[num_gridpz_m];
        zvals[0] = new double[num_gridpz_m];

        // read in and parse field map
        in.open(Filename_m.c_str());
        interpreteLine<string, int>(in, tmpString, accuracy);

        for(nsp = 0; nsp < num_gridpz_m && parsing_passed; /* skip increment of nsp here */) {
            parsing_passed = interpreteLine<double, double>(in, zvals[0][nsp], RealValues[nsp]);
            // the sequence of z-position should be strictly increasing
            // drop sampling points that don't comply to this
            if(zvals[0][nsp] - tmpDouble > 1e-10) {
                if(fabs(RealValues[nsp]) > Bz_max) {
                    Bz_max = fabs(RealValues[nsp]);
                }
                tmpDouble = zvals[0][nsp];
                ++ nsp; // increment nsp only if sampling point is accepted
            }
        }
        in.close();
        num_gridpz_m = nsp;
        hz_m = (zend_m - zbegin_m) / (num_gridpz_m - 1);

        higherDerivatives[0] = new double[num_gridpz_m];
        higherDerivatives[1] = new double[num_gridpz_m];
        higherDerivatives[2] = new double[num_gridpz_m];
        zvals[1] = new double[num_gridpz_m];

        for(int i = 0; i < 4; ++i) {
            onAxisInterpolants_m[i] = gsl_spline_alloc(gsl_interp_cspline, num_gridpz_m);
            onAxisAccel_m[i] = gsl_interp_accel_alloc();
        }

        gsl_fft_real_wavetable *real = gsl_fft_real_wavetable_alloc(2 * num_gridpz_m);
        gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(2 * num_gridpz_m);

        ii = num_gridpz_m - 1;
        for(int i = 0; i < num_gridpz_m; ++ i, -- ii) {
            onAxisField_m[i] = RealValues[i] / Bz_max;
            zvals[0][ii] -= zvals[0][0];
        }
        gsl_spline_init(onAxisInterpolants_m[0], zvals[0], onAxisField_m, num_gridpz_m);

        // get equidistant sampling from the, possibly, non-equidistant sampling
        // using cubic spline for this
        ii = num_gridpz_m;
        double z = 0.0;

        for(int i = 0; i < num_gridpz_m - 1; ++ i, ++ ii, z += hz_m) {
            zvals[1][i] = z;
            RealValues[ii] = gsl_spline_eval(onAxisInterpolants_m[0], z, onAxisAccel_m[0]);
        }
        zvals[1][num_gridpz_m - 1] = zvals[0][num_gridpz_m - 1];
        RealValues[ii ++] = onAxisField_m[num_gridpz_m - 1];
        // prepend mirror sampling points such that field values are periodic for sure
        // ii == 2*num_gridpz_m at the moment
        -- ii; 
        for(int i = 0; i < num_gridpz_m; ++ i, -- ii) {
            RealValues[i] = RealValues[ii];
        }

        gsl_fft_real_transform(RealValues, 1, 2 * num_gridpz_m, real, work);

        // normalize to Bz_max = 1 A/m
        RealValues[0] /= (2. * num_gridpz_m);
        for(int i = 1; i < 2 * accuracy - 1; i++) {
            RealValues[i] /= num_gridpz_m;
        }

        for(int i = 0; i < num_gridpz_m; ++ i) {
            const double kz = two_pi * (zvals[1][i] / length_m + 0.5);
            higherDerivatives[0][i] = 0.0;
            higherDerivatives[1][i] = 0.0;
            higherDerivatives[2][i] = 0.0;
            int n = 1;
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

        for(int i = 1; i < 4; ++i) {
            gsl_spline_init(onAxisInterpolants_m[i], zvals[1], higherDerivatives[i-1], num_gridpz_m);
        }

        delete[] RealValues;
        delete[] higherDerivatives[0];
        delete[] higherDerivatives[1];
        delete[] higherDerivatives[2];
        delete[] zvals[0];
        delete[] zvals[1];

        INFOMSG(typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << endl);
    }
}

void Astra1DMagnetoStatic_fast::freeMap() {
    if(onAxisField_m != NULL) {
        for(int i = 0; i < 4; ++i) {
            gsl_spline_free(onAxisInterpolants_m[i]);
            gsl_interp_accel_free(onAxisAccel_m[i]);
        }
        delete[] onAxisField_m;

        INFOMSG(typeset_msg("freed fieldmap '" + Filename_m  + "'", "info") << endl);
    }
}

bool Astra1DMagnetoStatic_fast::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const {
    // do fourier interpolation in z-direction
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

bool Astra1DMagnetoStatic_fast::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const {
    return false;
}

void Astra1DMagnetoStatic_fast::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
}

void Astra1DMagnetoStatic_fast::swap()
{ }

void Astra1DMagnetoStatic_fast::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (1D magnetostatic); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double Astra1DMagnetoStatic_fast::getFrequency() const {
    return 0.0;
}

void Astra1DMagnetoStatic_fast::setFrequency(double freq)
{ }
