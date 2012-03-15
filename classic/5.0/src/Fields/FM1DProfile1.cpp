#include <fstream>
#include <ios>

#include "Fields/FM1DProfile1.hh"
#include "Fields/Fieldmap.icc"
#include "Physics/Physics.h"

using namespace std;
using Physics::mu_0;
using Physics::c;
using Physics::two_pi;

FM1DProfile1::FM1DProfile1(string aFilename)
    : Fieldmap(aFilename),
      EngeCoefs_entry_m(NULL),
      EngeCoefs_exit_m(NULL),
      exit_slope_m(0.0) {
    Inform msg("FM1DP ");
    int tmpInt;
    string tmpString;
    double tmpDouble;

    Type = T1DProfile1;
    ifstream file(Filename_m.c_str());

    if(file.good()) {
        bool parsing_passed =                               \
                interpreteLine<string, int, int, double>(file,
                        tmpString,
                        polynomialOrder_entry_m,
                        polynomialOrder_exit_m,
                        gapHeight_m);
        parsing_passed = parsing_passed &&
                         interpreteLine<double, double, double, int>(file,
                                 zbegin_entry_m,
                                 polynomialOrigin_entry_m,
                                 zend_entry_m,
                                 tmpInt);
        parsing_passed = parsing_passed &&
                         interpreteLine<double, double, double, int>(file,
                                 zbegin_exit_m,
                                 polynomialOrigin_exit_m,
                                 zend_exit_m,
                                 tmpInt);
        for(int i = 0; (i < polynomialOrder_entry_m + polynomialOrder_exit_m + 2) && parsing_passed; ++ i) {
            parsing_passed = parsing_passed &&
                             interpreteLine<double>(file, tmpDouble);
        }

        parsing_passed = parsing_passed &&
                         interpreteEOF(file);

        file.close();
        lines_read_m = 0;

        if(!parsing_passed) {
            disableFieldmapWarning();
            zend_exit_m = zbegin_entry_m - 1e-3;
            zend_entry_m = zbegin_entry_m - 1e-3;
            zbegin_exit_m = zbegin_entry_m - 1e-3;
        } else {
            // conversion cm to m
            zbegin_entry_m /= 100.;
            zend_entry_m /= 100.;
            polynomialOrigin_entry_m /= 100.;
            zbegin_exit_m /= 100.;
            zend_exit_m /= 100.;
            polynomialOrigin_exit_m /= 100.;
            gapHeight_m /= 100.0;
        }
        length_m = zend_exit_m - zbegin_entry_m;
    } else {
        Inform msg("FM1DMS ");
        noFieldmapWarning();
        zbegin_entry_m = 0.0;
        zend_exit_m = zbegin_entry_m - 1e-3;
        zend_entry_m = zbegin_entry_m - 1e-3;
        zbegin_exit_m = zbegin_entry_m - 1e-3;
    }
}

FM1DProfile1::~FM1DProfile1() {
    if(EngeCoefs_entry_m != NULL) {
        delete[] EngeCoefs_entry_m;
        delete[] EngeCoefs_exit_m;
    }
}

void FM1DProfile1::readMap() {
    if(EngeCoefs_entry_m == NULL) {
        double tolerance = 1e-8;

        Inform msg("FM1DP1 ");
        ifstream in(Filename_m.c_str());

        int tmpInt;
        string tmpString;
        double tmpDouble;

        interpreteLine<string, int, int, double>(in, tmpString, tmpInt, tmpInt, tmpDouble);
        interpreteLine<double, double, double, int>(in, tmpDouble, tmpDouble, tmpDouble, tmpInt);
        interpreteLine<double, double, double, int>(in, tmpDouble, tmpDouble, tmpDouble, tmpInt);

        EngeCoefs_entry_m = new double[polynomialOrder_entry_m + 1];
        EngeCoefs_exit_m = new double[polynomialOrder_exit_m + 1];

        for(int i = 0; i < polynomialOrder_entry_m + 1; ++i) {
            interpreteLine<double>(in,  EngeCoefs_entry_m[i]);
        }
        for(int i = 0; i < polynomialOrder_exit_m + 1; ++i) {
            interpreteLine<double>(in, EngeCoefs_exit_m[i]);
        }
        in.close();

        msg << typeset_msg("read in fieldmap '" + Filename_m  + "'", "info") << "\n"
            << endl;

    }
}

void FM1DProfile1::freeMap() {
    if(EngeCoefs_entry_m != NULL) {
        Inform msg("FM1DMS ");

        delete[] EngeCoefs_entry_m;
        delete[] EngeCoefs_exit_m;

        msg << typeset_msg("freed fieldmap '" + Filename_m  + "'", "info") << "\n"
            << endl;

    }
}

bool FM1DProfile1::getFieldstrength(const Vector_t &R, Vector_t &strength, Vector_t &info) const {
    info = Vector_t(0.0);
    const Vector_t tmpR(R(0), R(1), R(2) + zbegin_entry_m);

    if(tmpR(2) >= zend_entry_m && tmpR(2) <= exit_slope_m * tmpR(0) + zbegin_exit_m) {
        strength = Vector_t(1.0, 0.0, 0.0);
    } else {
        double S, dSdz, d2Sdz2 = 0.0;
        double expS, f, dfdz, d2fdz2;
        double z;
        double *EngeCoefs;
        int polynomialOrder;
        info(0) = 1.0;
        if(tmpR(2) >= zbegin_entry_m && tmpR(2) < zend_entry_m) {
            z = -(tmpR(2) - polynomialOrigin_entry_m) / gapHeight_m;
            EngeCoefs = EngeCoefs_entry_m;
            polynomialOrder = polynomialOrder_entry_m;
        } else if(tmpR(2) > exit_slope_m * tmpR(0) + zbegin_exit_m && tmpR(2) <= exit_slope_m * tmpR(0) + zend_exit_m) {
            z = (tmpR(2) - exit_slope_m * tmpR(0) - polynomialOrigin_exit_m) / sqrt(exit_slope_m * exit_slope_m + 1) / gapHeight_m;
            EngeCoefs = EngeCoefs_exit_m;
            polynomialOrder = polynomialOrder_exit_m;
            info(1) = 1.0;
        } else {
            return true;
        }

        S = EngeCoefs[polynomialOrder] * z;
        S += EngeCoefs[polynomialOrder - 1];
        dSdz = polynomialOrder * EngeCoefs[polynomialOrder];

        for(int i = polynomialOrder - 2; i >= 0; i--) {
            S = S * z + EngeCoefs[i];
            dSdz = dSdz * z + (i + 1) * EngeCoefs[i+1];
            d2Sdz2 = d2Sdz2 * z + (i + 2) * (i + 1) * EngeCoefs[i+2];
        }
        expS = exp(S);
        f = 1.0 / (1.0 + expS);
        if(f > 1.e-30) {
            dfdz = - f * ((f * expS) * dSdz); // first derivative of f
            d2fdz2 = ((-d2Sdz2 - dSdz * dSdz * (1. - 2. * (expS * f))) * (f * expS) * f) / (gapHeight_m * gapHeight_m);  // second derivative of f

            strength(0) = f;
            strength(1) = dfdz / gapHeight_m;
            strength(2) = d2fdz2;
        } else {
            strength = Vector_t(0.0);
        }

    }
    info(2) = exit_slope_m;
    return true;

}

bool FM1DProfile1::getFieldstrength_fdiff(const Vector_t &R, Vector_t &E, Vector_t &B, const DiffDirection &dir) const {

}

void FM1DProfile1::getFieldDimensions(double &zBegin, double &zEnd, double &rBegin, double &rEnd) const {
    zBegin = zbegin_entry_m;
    zEnd = zend_exit_m;
}

void FM1DProfile1::swap()
{}

void FM1DProfile1::getInfo(Inform *msg) {
    (*msg) << Filename_m << " (1D Profile type 1); zini= " << zbegin_entry_m << " m; zfinal= " << zend_exit_m << " m;" << endl;
}

double FM1DProfile1::getFrequency() const {
    return 0.0;
}

void FM1DProfile1::setFrequency(double freq)
{}

void FM1DProfile1::setExitFaceSlope(const double &m) {
    exit_slope_m = m;
}
