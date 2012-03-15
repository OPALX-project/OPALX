#include <fstream>
#include <ios>

#include "Fields/FM1DDynamic_fast.hh"
#include "Fields/Fieldmap.icc"
#include "Physics/Physics.h"
#include "gsl/gsl_fft_real.h"

using namespace std;
using Physics::mu_0;
using Physics::c;
using Physics::two_pi;

FM1DDynamic_fast::FM1DDynamic_fast(string aFilename)
    : Fieldmap(aFilename),
      FieldstrengthEz_m(NULL),
      FieldstrengthEr_m(NULL),
      FieldstrengthBt_m(NULL) {
    Inform msg("FM1DD ");
    ifstream file;
    Type = T1DDynamic;

    // open field map, parse it and disable element on error
    file.open(Filename_m.c_str());
    if(file.good()) {
        int tmpInt;
        string tmpString;
        double tmpDouble;

        bool parsing_passed = interpreteLine<string, int>(file, tmpString, tmpInt);
        parsing_passed = parsing_passed &&
                         interpreteLine<double, double, int>(file, zbegin_m, zend_m, num_gridpz_m);
        parsing_passed = parsing_passed &&
                         interpreteLine<double>(file, frequency_m);
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
            // conversion MHz to Hz and frequency to angular frequency
            frequency_m *= Physics::two_pi * 1e6;

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
        num_gridpz_m = num_gridpr_m = 0;
        zbegin_m = 0.0;
        zend_m = -1e-3;
    }
}

FM1DDynamic_fast::~FM1DDynamic_fast() {
    if(FieldstrengthEz_m != NULL) {
        delete[] FieldstrengthEz_m;
        delete[] FieldstrengthEr_m;
        delete[] FieldstrengthBt_m;
    }
}

void FM1DDynamic_fast::readMap() {
    if(FieldstrengthEz_m == NULL) {
        // declare variables and allocate memory
        Inform msg("FM1DD ");
        ifstream in;

        string tmpString;

        int tmpInt;
        int accuracy;

        double tmpDouble;
        double kz, R;
        double xlrep = frequency_m / c;
        double ez, ezp, ezpp, ezppp, f, fp;
        double somefactor, somefactor_base;
        double Ez_max = 0.0;

        double *RealValues = new double[2*num_gridpz_m];

        gsl_fft_real_wavetable *real = gsl_fft_real_wavetable_alloc(2 * num_gridpz_m);
        gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc(2 * num_gridpz_m);

        FieldstrengthEz_m = new double[num_gridpz_m * num_gridpr_m];
        FieldstrengthEr_m = new double[num_gridpz_m * num_gridpr_m];
        FieldstrengthBt_m = new double[num_gridpz_m * num_gridpr_m];

        // read in field map and parse it
        in.open(Filename_m.c_str());
        interpreteLine<string, int>(in, tmpString, accuracy);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);
        interpreteLine<double>(in, tmpDouble);
        interpreteLine<double, double, int>(in, tmpDouble, tmpDouble, tmpInt);

        for(int i = 0; i < num_gridpz_m; ++ i) {
            interpreteLine<double>(in, RealValues[num_gridpz_m + i]);
            if(fabs(RealValues[num_gridpz_m + i]) > Ez_max) {
                Ez_max = fabs(RealValues[num_gridpz_m + i]);
            }
            // mirror the field map to make sure that it is periodic
            RealValues[num_gridpz_m - 1 - i] = RealValues[num_gridpz_m + i];
        }
        in.close();

        // apply FFT
        gsl_fft_real_transform(RealValues, 1, 2 * num_gridpz_m, real, work);

        // normalize the fourier coeficients such that [E] = V/m after interpolation
        // the 1.e6 stems from conversion MV/m to V/m
        RealValues[0] *= 1.e6 / (Ez_max * 2 * num_gridpz_m);
        for(int i = 1; i < 2 * accuracy - 1; i++) {
            RealValues[i] *=  1.e6 / (Ez_max * num_gridpz_m);
        }

        // interpolate the on-axis field
        for(int i = 0; i < num_gridpz_m; i++) {
            kz = Physics::pi / num_gridpz_m * i + Physics::pi; // additional phase of pi due to mirroring
            ez = RealValues[0];
            ezp = 0.0;
            ezpp = 0.0;
            ezppp = 0.0;
            int n = 1;
            for(int l = 1; l < accuracy ; l++, n += 2) {
                somefactor_base = l * Physics::pi / (num_gridpz_m * hz_m);
                somefactor = 1.0;
                const double coskzl = cos(kz * l);
                const double sinkzl = sin(kz * l);
                ez    += (RealValues[n] * coskzl - RealValues[n+1] * sinkzl);
                somefactor *= somefactor_base;
                ezp   += somefactor * (-RealValues[n] * sinkzl - RealValues[n+1] * coskzl);
                somefactor *= somefactor_base;
                ezpp  += somefactor * (-RealValues[n] * coskzl + RealValues[n+1] * sinkzl);
                somefactor *= somefactor_base;
                ezppp += somefactor * (RealValues[n] * sinkzl + RealValues[n+1] * coskzl);
            }

            f  = -(ezpp  + ez *  xlrep * xlrep) / 16.;
            fp = -(ezppp + ezp * xlrep * xlrep) / 16.;

            // expand the field at off-axis sampling points
            for(int j = 0; j < num_gridpr_m; j++) {
                R = j * hr_m;
                FieldstrengthEz_m[i + j * num_gridpz_m] =   ez + 4 * f * R * R;
                FieldstrengthEr_m[i + j * num_gridpz_m] = -(ezp / 2. + fp * R * R);
                FieldstrengthBt_m[i + j * num_gridpz_m] = (ez / 2. + f * R * R) * xlrep / c;
            }
        }

        // free the memory
        gsl_fft_real_workspace_free(work);
        gsl_fft_real_wavetable_free(real);

        delete[] RealValues;


        msg << typeset_msg("read in fieldmap '" + Filename_m + "'", "info") << "\n"
            << endl;
    }
}

void FM1DDynamic_fast::freeMap() {
    if(FieldstrengthEz_m != NULL) {
        Inform msg("FM1DD ");
        delete[] FieldstrengthEz_m;
        delete[] FieldstrengthEr_m;
        delete[] FieldstrengthBt_m;

        msg << typeset_msg("freed fieldmap '" + Filename_m  + "'", "info") << "\n"
            << endl;
    }
}


bool FM1DDynamic_fast::getFieldstrength(const Vector_t &R, Vector_t &E, Vector_t &B) const {
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
    const double EfieldR = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEr_m[index1]
                           + leverz         * (1.0 - leverr) * FieldstrengthEr_m[index1 + 1]
                           + (1.0 - leverz) * leverr         * FieldstrengthEr_m[index2]
                           + leverz         * leverr         * FieldstrengthEr_m[index2 + 1];

    const double EfieldZ = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthEz_m[index1]
                           + leverz         * (1.0 - leverr) * FieldstrengthEz_m[index1 + 1]
                           + (1.0 - leverz) * leverr         * FieldstrengthEz_m[index2]
                           + leverz         * leverr         * FieldstrengthEz_m[index2 + 1];

    const double BfieldT = (1.0 - leverz) * (1.0 - leverr) * FieldstrengthBt_m[index1]
                           + leverz         * (1.0 - leverr) * FieldstrengthBt_m[index1 + 1]
                           + (1.0 - leverz) * leverr         * FieldstrengthBt_m[index2]
                           + leverz         * leverr         * FieldstrengthBt_m[index2 + 1];

    E(0) += EfieldR * R(0);
    E(1) += EfieldR * R(1);
    E(2) += EfieldZ;
    B(0) += -BfieldT * R(1);
    B(1) +=  BfieldT * R(0);

    return false;
}

bool FM1DDynamic_fast::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const
{ }


void FM1DDynamic_fast::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const {
    zBegin = zbegin_m;
    zEnd = zend_m;
    rBegin = rbegin_m;
    rEnd = rend_m;
}

void FM1DDynamic_fast::swap()
{ }

void FM1DDynamic_fast::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (1D dynamic fast); zini= " << zbegin_m << " m; zfinal= " << zend_m << " m;" << endl;
}

double FM1DDynamic_fast::getFrequency() const {
    return frequency_m;
}

void FM1DDynamic_fast::setFrequency(double freq) {
    frequency_m = freq;
}
