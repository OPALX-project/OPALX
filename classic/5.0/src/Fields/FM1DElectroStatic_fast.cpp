#include <fstream>
#include <ios>

#include "Fields/FM1DElectroStatic_fast.hh"
#include "Fields/Fieldmap.icc"
#include "Physics/Physics.h"

#include "gsl/gsl_fft_real.h"

FM1DElectroStatic_fast::FM1DElectroStatic_fast(std::string aFilename):
    Fieldmap(aFilename) {

    Type = T1DElectroStatic;
    onAxisField_m = NULL;

    std::ifstream fieldFile(Filename_m.c_str());
    if(fieldFile.good()) {

        bool parsingPassed = readFileHeader(fieldFile);
        parsingPassed = checkFileData(fieldFile, parsingPassed);
        fieldFile.close();

        if(!parsingPassed) {
            disableFieldmapWarning();
            zEnd_m = zBegin_m - 1.0e-3;
        } else
            convertHeaderData();

        deltaZ_m = (zEnd_m - zBegin_m) / (numberOfGridPoints_m - 1);
        length_m = 2.0 * numberOfGridPoints_m * deltaZ_m;

    } else {
        noFieldmapWarning();
        zBegin_m = 0.0;
        zEnd_m = -1.0e-3;
    }
}

FM1DElectroStatic_fast::~FM1DElectroStatic_fast() {
    if(onAxisField_m != NULL) {
        delete [] onAxisField_m;
        gsl_spline_free(onAxisFieldInterpolants_m);
        gsl_spline_free(onAxisFieldPInterpolants_m);
        gsl_spline_free(onAxisFieldPPInterpolants_m);
        gsl_spline_free(onAxisFieldPPPInterpolants_m);
        gsl_interp_accel_free(onAxisFieldAccel_m);
        gsl_interp_accel_free(onAxisFieldPAccel_m);
        gsl_interp_accel_free(onAxisFieldPPAccel_m);
        gsl_interp_accel_free(onAxisFieldPPPAccel_m);
    }
}

void FM1DElectroStatic_fast::readMap() {
    if(onAxisField_m == NULL) {

        std::ifstream fieldFile(Filename_m.c_str());
        int accuracy = stripFileHeader(fieldFile);

        onAxisField_m = new double[numberOfGridPoints_m];
        double maxEz = readFileData(fieldFile, onAxisField_m);
        fieldFile.close();

        std::vector<double> fourierCoefs
        = computeFourierCoefficients(accuracy,
                                     onAxisField_m);
        normalizeField(maxEz, fourierCoefs);

        double *onAxisFieldP = new double[numberOfGridPoints_m];
        double *onAxisFieldPP = new double[numberOfGridPoints_m];
        double *onAxisFieldPPP = new double[numberOfGridPoints_m];
        computeFieldDerivatives(accuracy, fourierCoefs, onAxisFieldP,
                                onAxisFieldPP, onAxisFieldPPP);
        computeInterpolationVectors(onAxisFieldP, onAxisFieldPP,
                                    onAxisFieldPPP);

        delete [] onAxisFieldP;
        delete [] onAxisFieldPP;
        delete [] onAxisFieldPPP;

        INFOMSG(typeset_msg("read in fieldmap '" + Filename_m  + "'", "info")
                << endl);
    }
}

void FM1DElectroStatic_fast::freeMap() {
    if(onAxisField_m != NULL) {
        delete [] onAxisField_m;
        gsl_spline_free(onAxisFieldInterpolants_m);
        gsl_spline_free(onAxisFieldPInterpolants_m);
        gsl_spline_free(onAxisFieldPPInterpolants_m);
        gsl_spline_free(onAxisFieldPPPInterpolants_m);
        gsl_interp_accel_free(onAxisFieldAccel_m);
        gsl_interp_accel_free(onAxisFieldPAccel_m);
        gsl_interp_accel_free(onAxisFieldPPAccel_m);
        gsl_interp_accel_free(onAxisFieldPPPAccel_m);

        INFOMSG(typeset_msg("freed fieldmap '" + Filename_m  + "'", "info")
                << endl);
    }
}

bool FM1DElectroStatic_fast::getFieldstrength(const Vector_t &R, Vector_t &E,
        Vector_t &B) const {

    std::vector<double> fieldComponents;
    computeFieldOnAxis(R(2), fieldComponents);
    computeFieldOffAxis(R, E, B, fieldComponents);

    return false;

}

bool FM1DElectroStatic_fast::getFieldDerivative(const Vector_t &R,
        Vector_t &E,
        Vector_t &B,
        const DiffDirection &dir) const {

    E(2) += gsl_spline_eval(onAxisFieldPInterpolants_m, R(2),
                            onAxisFieldPAccel_m);

    return false;

}

void FM1DElectroStatic_fast::getFieldDimensions(double &zBegin, double &zEnd,
        double &rBegin, double &rEnd) const {
    zBegin = zBegin_m;
    zEnd = zEnd_m;
    rBegin = rBegin_m;
    rEnd = rEnd_m;
}
void FM1DElectroStatic_fast::getFieldDimensions(double &xIni, double &xFinal,
        double &yIni, double &yFinal,
        double &zIni, double &zFinal) const {}

void FM1DElectroStatic_fast::swap()
{ }

void FM1DElectroStatic_fast::getInfo(Inform *msg) {
    (*msg) << Filename_m
           << " (1D electrotostatic); zini= "
           << zBegin_m << " m; zfinal= "
           << zEnd_m << " m;" << endl;
}

double FM1DElectroStatic_fast::getFrequency() const {
    return 0.0;
}

void FM1DElectroStatic_fast::setFrequency(double freq)
{ }

bool FM1DElectroStatic_fast::checkFileData(std::ifstream &fieldFile,
        bool parsingPassed) {

    double tempDouble;
    for(int dataIndex = 0; dataIndex <= numberOfGridPoints_m; dataIndex++)
        parsingPassed = parsingPassed
                        && interpreteLine<double>(fieldFile, tempDouble);

    return parsingPassed && interpreteEOF(fieldFile);

}

void FM1DElectroStatic_fast::computeFieldDerivatives(int accuracy,
        std::vector<double> fourierCoefs,
        double onAxisFieldP[],
        double onAxisFieldPP[],
        double onAxisFieldPPP[]) {

    for(int zStepIndex = 0; zStepIndex < numberOfGridPoints_m; zStepIndex++) {

        double z = deltaZ_m * zStepIndex;
        double kz = Physics::two_pi * z / length_m + Physics::pi;
        onAxisFieldP[zStepIndex] = 0.0;
        onAxisFieldPP[zStepIndex] = 0.0;
        onAxisFieldPPP[zStepIndex] = 0.0;

        int coefIndex = 1;
        for(int n = 1; n < accuracy; n++) {

            double kn = n * Physics::two_pi / length_m;
            double coskzn = cos(kz * n);
            double sinkzn = sin(kz * n);

            onAxisFieldP[zStepIndex] += kn * (-fourierCoefs.at(coefIndex) * sinkzn
                                              - fourierCoefs.at(coefIndex + 1) * coskzn);

            double derivCoeff = pow(kn, 2.0);
            onAxisFieldPP[zStepIndex] += derivCoeff * (-fourierCoefs.at(coefIndex) * coskzn
                                         + fourierCoefs.at(coefIndex + 1) * sinkzn);
            derivCoeff *= kn;
            onAxisFieldPPP[zStepIndex] += derivCoeff * (fourierCoefs.at(coefIndex) * sinkzn
                                          + fourierCoefs.at(coefIndex + 1) * coskzn);

            coefIndex += 2;
        }
    }
}

void FM1DElectroStatic_fast::computeFieldOffAxis(const Vector_t &R, Vector_t &E, Vector_t &B,
        std::vector<double> fieldComponents) const {

    double radiusSq = pow(R(0), 2.0) + pow(R(1), 2.0);
    double transverseEFactor = -fieldComponents.at(1) / 2.0
                               + radiusSq * fieldComponents.at(3) / 16.0;

    E(0) += R(0) * transverseEFactor;
    E(1) += R(1) * transverseEFactor;
    E(2) += fieldComponents.at(0) - fieldComponents.at(2) * radiusSq / 4.0;

}

void FM1DElectroStatic_fast::computeFieldOnAxis(double z,
        std::vector<double> &fieldComponents)
const {

    fieldComponents.push_back(gsl_spline_eval(onAxisFieldInterpolants_m,
                              z, onAxisFieldAccel_m));
    fieldComponents.push_back(gsl_spline_eval(onAxisFieldPInterpolants_m,
                              z, onAxisFieldPAccel_m));
    fieldComponents.push_back(gsl_spline_eval(onAxisFieldPPInterpolants_m,
                              z, onAxisFieldPPAccel_m));
    fieldComponents.push_back(gsl_spline_eval(onAxisFieldPPPInterpolants_m,
                              z, onAxisFieldPPPAccel_m));
}

std::vector<double> FM1DElectroStatic_fast::computeFourierCoefficients(int accuracy,
        double fieldData[]) {

    gsl_fft_real_wavetable *waveTable = gsl_fft_real_wavetable_alloc(2
                                        * numberOfGridPoints_m - 1);
    gsl_fft_real_workspace *workSpace = gsl_fft_real_workspace_alloc(2
                                        * numberOfGridPoints_m - 1);

    // Reflect field about minimum z value to ensure that it is periodic.
    double *fieldDataReflected = new double[2 * numberOfGridPoints_m - 1];
    for(int dataIndex = 0; dataIndex < numberOfGridPoints_m; dataIndex++) {
        fieldDataReflected[numberOfGridPoints_m - 1 + dataIndex]
        = fieldData[dataIndex];
        if(dataIndex != 0)
            fieldDataReflected[numberOfGridPoints_m - 1 - dataIndex]
            = fieldData[dataIndex];
    }

    gsl_fft_real_transform(fieldDataReflected, 1,
                           2 * numberOfGridPoints_m - 1,
                           waveTable, workSpace);

    std::vector<double> fourierCoefs;
    fourierCoefs.push_back(fieldDataReflected[0]
                           / (2.0 * numberOfGridPoints_m - 1));
    for(int coefIndex = 1; coefIndex < 2 * accuracy - 1; coefIndex++)
        fourierCoefs.push_back(2.0 * fieldDataReflected[coefIndex]
                               / (2.0 * numberOfGridPoints_m - 1));

    delete [] fieldDataReflected;
    gsl_fft_real_workspace_free(workSpace);
    gsl_fft_real_wavetable_free(waveTable);

    return fourierCoefs;

}

void FM1DElectroStatic_fast::computeInterpolationVectors(double onAxisFieldP[],
        double onAxisFieldPP[],
        double onAxisFieldPPP[]) {

    onAxisFieldInterpolants_m = gsl_spline_alloc(gsl_interp_cspline,
                                numberOfGridPoints_m);
    onAxisFieldPInterpolants_m = gsl_spline_alloc(gsl_interp_cspline,
                                 numberOfGridPoints_m);
    onAxisFieldPPInterpolants_m = gsl_spline_alloc(gsl_interp_cspline,
                                  numberOfGridPoints_m);
    onAxisFieldPPPInterpolants_m = gsl_spline_alloc(gsl_interp_cspline,
                                   numberOfGridPoints_m);

    double *z = new double[numberOfGridPoints_m];
    for(int zStepIndex = 0; zStepIndex < numberOfGridPoints_m; zStepIndex++)
        z[zStepIndex] = deltaZ_m * zStepIndex;
    gsl_spline_init(onAxisFieldInterpolants_m, z,
                    onAxisField_m, numberOfGridPoints_m);
    gsl_spline_init(onAxisFieldPInterpolants_m, z,
                    onAxisFieldP, numberOfGridPoints_m);
    gsl_spline_init(onAxisFieldPPInterpolants_m, z,
                    onAxisFieldPP, numberOfGridPoints_m);
    gsl_spline_init(onAxisFieldPPPInterpolants_m, z,
                    onAxisFieldPPP, numberOfGridPoints_m);

    onAxisFieldAccel_m = gsl_interp_accel_alloc();
    onAxisFieldPAccel_m = gsl_interp_accel_alloc();
    onAxisFieldPPAccel_m = gsl_interp_accel_alloc();
    onAxisFieldPPPAccel_m = gsl_interp_accel_alloc();

    delete [] z;
}

void FM1DElectroStatic_fast::convertHeaderData() {

    // Convert to m.
    rBegin_m /= 100.0;
    rEnd_m /= 100.0;
    zBegin_m /= 100.0;
    zEnd_m /= 100.0;

    // Convert number of grid spacings to number of grid points.
    numberOfGridPoints_m++;
}

void FM1DElectroStatic_fast::normalizeField(double maxEz, std::vector<double> &fourierCoefs) {

    for(int dataIndex = 0; dataIndex < numberOfGridPoints_m; dataIndex++)
        onAxisField_m[dataIndex] /= maxEz;

    for(std::vector<double>::iterator fourierIterator = fourierCoefs.begin();
        fourierIterator < fourierCoefs.end(); fourierIterator++)
        *fourierIterator /= maxEz;

}

double FM1DElectroStatic_fast::readFileData(std::ifstream &fieldFile,
        double fieldData[]) {

    double maxEz = 0.0;
    for(int dataIndex = 0; dataIndex < numberOfGridPoints_m; dataIndex++) {
        interpreteLine<double>(fieldFile, fieldData[dataIndex]);
        if(std::abs(fieldData[dataIndex]) > maxEz)
            maxEz = std::abs(fieldData[dataIndex]);
    }

    return maxEz;
}

bool FM1DElectroStatic_fast::readFileHeader(std::ifstream &fieldFile) {

    std::string tempString;
    int tempInt;

    bool parsingPassed = interpreteLine<std::string, int>(fieldFile, tempString,
                         tempInt);
    parsingPassed = parsingPassed &&
                    interpreteLine<double, double, int>(fieldFile,
                            zBegin_m,
                            zEnd_m,
                            numberOfGridPoints_m);
    parsingPassed = parsingPassed &&
                    interpreteLine<double, double, int>(fieldFile, rBegin_m,
                            rEnd_m, tempInt);
    return parsingPassed;
}

int FM1DElectroStatic_fast::stripFileHeader(std::ifstream &fieldFile) {

    std::string tempString;
    int tempInt;
    int accuracy;
    double tempDouble;

    interpreteLine<std::string, int>(fieldFile, tempString, accuracy);
    interpreteLine<double, double, int>(fieldFile, tempDouble,
                                        tempDouble, tempInt);
    interpreteLine<double, double, int>(fieldFile, tempDouble,
                                        tempDouble, tempInt);

    if(accuracy > numberOfGridPoints_m)
        accuracy = numberOfGridPoints_m;

    return accuracy;
}
