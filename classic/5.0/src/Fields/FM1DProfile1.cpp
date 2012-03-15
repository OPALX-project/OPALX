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
      exit_slope_m(0.0),
      xExit_m(0.0),
      zExit_m(0.0),
      cosEntranceRotation_m(1.0),
      sinEntranceRotation_m(0.0),
      cosExitRotation_m(1.0),
      sinExitRotation_m(0.0) {

    // Change inform message.
    Inform msg("FM1DP ");

    // Read file header information. Set field type.
    Type = T1DProfile1;

    // Read file. Check if we are using the default field profile.
    if(Filename_m == "1DPROFILE1-DEFAULT") {

        // Use default field profile coefficients. This
        // is a bend that is 100.0 cm long with a 2.0 cm
        // gap. Both the length and the gap can be modified
        // in the bend declaration.
        polynomialOrder_entry_m = 5;
        polynomialOrder_exit_m = 5;
        gapHeight_m = 2.0 / 100.0;

        zbegin_entry_m = -10.0 / 100.0;
        polynomialOrigin_entry_m = 0.0;
        zend_entry_m = 10.0 / 100.0;

        zbegin_exit_m = 90.0 / 100.0;
        polynomialOrigin_exit_m = 100.0 / 100.0;
        zend_exit_m = 110.0 / 100.0;


    } else {

        ifstream file(Filename_m.c_str());

        if(file.good()) {

            // Declares some dummy variables.
            int tmpInt = 0;
            string tmpString = "";
            double tmpDouble = 0.0;

            // Parse field map file header. Check that it at
            // least has the right number of values on each line.
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

            if(!parsing_passed) {

                // Disable map if parsing detects errors.
                disableFieldmapWarning();
                zend_exit_m = zbegin_entry_m - 1e-3;
                zend_entry_m = zbegin_entry_m - 1e-3;
                zbegin_exit_m = zbegin_entry_m - 1e-3;

            } else {

                // Convert from cm to m.
                zbegin_entry_m /= 100.;
                zend_entry_m /= 100.;
                polynomialOrigin_entry_m /= 100.;
                zbegin_exit_m /= 100.;
                zend_exit_m /= 100.;
                polynomialOrigin_exit_m /= 100.;
                gapHeight_m /= 100.0;

            }

            // Compute map length.
            length_m = zend_exit_m - zbegin_entry_m;

        } else {

            // No field map file or field map file is somehow flawed.
            Inform msg("FM1DMS ");
            noFieldmapWarning();
            zbegin_entry_m = 0.0;
            zend_exit_m = zbegin_entry_m - 1e-3;
            zend_entry_m = zbegin_entry_m - 1e-3;
            zbegin_exit_m = zbegin_entry_m - 1e-3;
        }
    }

    // Set origin of exit.
    xExit_m = 0.0;
    zExit_m = polynomialOrigin_exit_m;
}

FM1DProfile1::~FM1DProfile1() {
    if(EngeCoefs_entry_m != NULL) {
        delete[] EngeCoefs_entry_m;
        delete[] EngeCoefs_exit_m;
    }
}

void FM1DProfile1::readMap() {
    if(EngeCoefs_entry_m == NULL) {

        Inform msg("FM1DP1 ");

        // Check if using default profile.
        if(Filename_m == "1DPROFILE1-DEFAULT") {

            EngeCoefs_entry_m = new double[polynomialOrder_entry_m + 1];
            EngeCoefs_exit_m = new double[polynomialOrder_exit_m + 1];

            EngeCoefs_entry_m[0] = 0.478959;
            EngeCoefs_exit_m[0] = 0.478959;

            EngeCoefs_entry_m[1] = 1.911289;
            EngeCoefs_exit_m[1] = 1.911289;

            EngeCoefs_entry_m[2] = -1.185953;
            EngeCoefs_exit_m[2] = -1.185953;

            EngeCoefs_entry_m[3] = 1.630554;
            EngeCoefs_exit_m[3] = 1.630554;

            EngeCoefs_entry_m[4] = -1.082657;
            EngeCoefs_exit_m[4] = -1.082657;

            EngeCoefs_entry_m[5] = 0.318111;
            EngeCoefs_exit_m[5] = 0.31811;

        } else {

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

    // Find coordinates in the entrance frame.
    Vector_t REntrance(R(0), 0.0, 0.0);
    REntrance(0) = R(0) * cosEntranceRotation_m + (R(2) + zbegin_entry_m - polynomialOrigin_entry_m) * sinEntranceRotation_m;
    REntrance(2) = -R(0) * sinEntranceRotation_m + (R(2) + zbegin_entry_m - polynomialOrigin_entry_m) * cosEntranceRotation_m + polynomialOrigin_entry_m;

    // Find coordinates in the exit frame centered on (xExit_m, zExit_m)
    // and rotated so that z is perpendicular to exit face.
    Vector_t RExit(0.0, R(1), 0.0);

    RExit(0) = (R(0) - xExit_m) * cosExitRotation_m + (R(2) + zbegin_entry_m - zExit_m) * sinExitRotation_m;
    RExit(2) = -(R(0) - xExit_m) * sinExitRotation_m + (R(2) + zbegin_entry_m - zExit_m) * cosExitRotation_m + polynomialOrigin_exit_m;

    if(REntrance(2) >= zend_entry_m && RExit(2) <= zbegin_exit_m) {
        strength = Vector_t(1.0, 0.0, 0.0);
    } else {
        double S, dSdz, d2Sdz2 = 0.0;
        double expS, f, dfdz, d2fdz2;
        double z;
        double *EngeCoefs;
        int polynomialOrder;
        info(0) = 1.0;
        if(REntrance(2) >= zbegin_entry_m && REntrance(2) < zend_entry_m) {
            z = -(REntrance(2) - polynomialOrigin_entry_m) / gapHeight_m;
            EngeCoefs = EngeCoefs_entry_m;
            polynomialOrder = polynomialOrder_entry_m;
        } else if(RExit(2) > zbegin_exit_m && RExit(2) <= zend_exit_m) {
            z = (RExit(2) - polynomialOrigin_exit_m) / gapHeight_m;
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
            dSdz = dSdz * z + (i + 1) * EngeCoefs[i + 1];
            d2Sdz2 = d2Sdz2 * z + (i + 2) * (i + 1) * EngeCoefs[i + 2];
        }
        expS = exp(S);
        f = 1.0 / (1.0 + expS);
        if(f > 1.e-30) {
            // First derivative of Enge function, f.
            dfdz = - f * ((f * expS) * dSdz);

            // Second derivative of Enge functioin, f.
            d2fdz2 = ((-d2Sdz2 - dSdz * dSdz * (1. - 2. * (expS * f))) * (f * expS) * f) / (gapHeight_m * gapHeight_m);

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
    return false;
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

void FM1DProfile1::setEdgeConstants(const double &bendAngle, const double &entranceAngle, const double &exitAngle) {

    double deltaZ = polynomialOrigin_exit_m - polynomialOrigin_entry_m;
    zExit_m = polynomialOrigin_entry_m + deltaZ * cos(bendAngle / 2.0);
    xExit_m = -deltaZ * sin(bendAngle / 2.0);

    cosEntranceRotation_m = cos(entranceAngle);
    sinEntranceRotation_m = sin(entranceAngle);

    cosExitRotation_m = cos(bendAngle - exitAngle);
    sinExitRotation_m = sin(bendAngle - exitAngle);
}

void FM1DProfile1::setFieldGap(const double &gap) {
    gapHeight_m = gap;
}

void FM1DProfile1::setFieldLength(const double &length) {
    double deltaExitEntrance = zbegin_exit_m - polynomialOrigin_exit_m;
    double deltaExitExit = zend_exit_m - polynomialOrigin_exit_m;
    polynomialOrigin_exit_m = polynomialOrigin_entry_m + length;
    zbegin_exit_m = polynomialOrigin_exit_m + deltaExitEntrance;
    zend_exit_m = polynomialOrigin_exit_m + deltaExitExit;
}

bool FM1DProfile1::adjustFringeFields() {

    bool adjustedPositions = false;

    // Adjust zbegin_entry_m. If the field strength is greater than 1.0e-4, move fringe entry out
    // by 1 mm steps.
    Vector_t R(0.0);
    Vector_t strength(0.0);
    Vector_t info(0.0);

    getFieldstrength(R, strength, info);
    while(strength(0) > 1.0e-4) {
        zbegin_entry_m -= 0.001;
        getFieldstrength(R, strength, info);
        adjustedPositions = true;
    }

    // Adjust zend_entry_m. Must move in slightly (we'll do 1 mm) to make
    // sure we are in the fringe field region. Adjust fringe region exit
    // by 1 mm steps until field strength is 0.9999 or we are at the midpoint
    // of the magnet.
    R(2) = -zbegin_entry_m + zend_entry_m - 0.001;
    getFieldstrength(R, strength, info);
    while(strength(0) < 0.9999 && zend_entry_m < (polynomialOrigin_entry_m + polynomialOrigin_exit_m) / 2.0) {
        zend_entry_m += 0.001;
        R(2) = -zbegin_entry_m + zend_entry_m - 0.001;
        getFieldstrength(R, strength, info);
        adjustedPositions = true;
    }
    if(zend_entry_m > (polynomialOrigin_entry_m + polynomialOrigin_exit_m) / 2.0)
        zend_entry_m = (polynomialOrigin_entry_m + polynomialOrigin_exit_m) / 2.0;

    // Adjust zbegin_exit_m. Similar to zend_entry_m adjustment.
    R(2) = -zbegin_entry_m + zbegin_exit_m + 0.001;
    getFieldstrength(R, strength, info);
    while(strength(0) < 0.9999 && zbegin_exit_m > (polynomialOrigin_entry_m + polynomialOrigin_exit_m) / 2.0) {
        zbegin_exit_m -= 0.001;
        R(2) = -zbegin_entry_m + zbegin_exit_m + 0.001;
        getFieldstrength(R, strength, info);
        adjustedPositions = true;
    }
    if(zbegin_exit_m < (polynomialOrigin_entry_m + polynomialOrigin_exit_m) / 2.0)
        zbegin_exit_m = (polynomialOrigin_entry_m + polynomialOrigin_exit_m) / 2.0;

    // Adjust zend_exit_m. Similar to zbegin_entry_m adjustment.
    R(2) = -zbegin_entry_m + zend_exit_m - 0.001;
    getFieldstrength(R, strength, info);
    while(strength(0) > 1.0e-4) {
        zend_exit_m += 0.001;
        R(2) = -zbegin_entry_m + zend_exit_m - 0.001;
        cout << R(2) << " " << strength(0) << endl;
        getFieldstrength(R, strength, info);
        adjustedPositions = true;
    }

    return adjustedPositions;
}

