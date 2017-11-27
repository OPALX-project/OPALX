// ------------------------------------------------------------------------
// $RCSfile: Bend.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: Bend
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

#include "AbsBeamline/Bend.h"
#include "Algorithms/PartPusher.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Utilities/Options.h"
#include "Fields/Fieldmap.h"
#include "AbstractObjects/OpalData.h"
#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class Bend
// ------------------------------------------------------------------------

Bend::Bend():
    Component(),
    bX_m(0.0),
    bY_m(0.0),
    designEnergy_m(0.0),
    messageHeader_m(" * "),
    pusher_m(),
    fileName_m(""),
    fieldmap_m(NULL),
    fast_m(false),
    angle_m(0.0),
    aperture_m(0.0),
    designRadius_m(0.0),
    fieldAmplitude_m(0.0),
    angleGreaterThanPi_m(false),
    entranceAngle_m(0.0),
    exitAngle_m(0.0),
    fieldIndex_m(0.0),
    elementEdge_m(0.0),
    startField_m(0.0),
    endField_m(0.0),
    reinitialize_m(false),
    recalcRefTraj_m(false),
    length_m(0.0),
    gap_m(0.0),
    refTrajMapSize_m(0),
    refTrajMapStepSize_m(0.0),
    entranceParameter1_m(0.0),
    entranceParameter2_m(0.0),
    entranceParameter3_m(0.0),
    exitParameter1_m(0.0),
    exitParameter2_m(0.0),
    exitParameter3_m(0.0),
    xOriginEngeEntry_m(0.0),
    zOriginEngeEntry_m(0.0),
    deltaBeginEntry_m(0.0),
    deltaEndEntry_m(0.0),
    polyOrderEntry_m(0),
    xExit_m(0.0),
    zExit_m(0.0),
    xOriginEngeExit_m(0.0),
    zOriginEngeExit_m(0.0),
    deltaBeginExit_m(0.0),
    deltaEndExit_m(0.0),
    polyOrderExit_m(0),
    cosEntranceAngle_m(1.0),
    sinEntranceAngle_m(0.0),
    exitEdgeAngle_m(0.0),
    cosExitAngle_m(1.0),
    sinExitAngle_m(0.0) {

    setElType(isDipole);

}

Bend::Bend(const Bend &right):
    Component(right),
    bX_m(right.bX_m),
    bY_m(right.bY_m),
    designEnergy_m(right.designEnergy_m),
    messageHeader_m(" * "),
    pusher_m(right.pusher_m),
    fileName_m(right.fileName_m),
    fieldmap_m(right.fieldmap_m),
    fast_m(right.fast_m),
    angle_m(right.angle_m),
    aperture_m(right.aperture_m),
    designRadius_m(right.designRadius_m),
    fieldAmplitude_m(right.fieldAmplitude_m),
    angleGreaterThanPi_m(right.angleGreaterThanPi_m),
    entranceAngle_m(right.entranceAngle_m),
    exitAngle_m(right.exitAngle_m),
    fieldIndex_m(right.fieldIndex_m),
    elementEdge_m(right.elementEdge_m),
    startField_m(right.startField_m),
    endField_m(right.endField_m),
    reinitialize_m(right.reinitialize_m),
    recalcRefTraj_m(right.recalcRefTraj_m),
    length_m(right.length_m),
    gap_m(right.gap_m),
    refTrajMapX_m(right.refTrajMapX_m),
    refTrajMapY_m(right.refTrajMapY_m),
    refTrajMapZ_m(right.refTrajMapZ_m),
    refTrajMapSize_m(right.refTrajMapSize_m),
    refTrajMapStepSize_m(right.refTrajMapStepSize_m),
    entranceParameter1_m(right.entranceParameter1_m),
    entranceParameter2_m(right.entranceParameter2_m),
    entranceParameter3_m(right.entranceParameter3_m),
    exitParameter1_m(right.exitParameter1_m),
    exitParameter2_m(right.exitParameter2_m),
    exitParameter3_m(right.exitParameter3_m),
    xOriginEngeEntry_m(right.xOriginEngeEntry_m),
    zOriginEngeEntry_m(right.zOriginEngeEntry_m),
    deltaBeginEntry_m(right.deltaBeginEntry_m),
    deltaEndEntry_m(right.deltaEndEntry_m),
    polyOrderEntry_m(right.polyOrderEntry_m),
    xExit_m(right.xExit_m),
    zExit_m(right.zExit_m),
    xOriginEngeExit_m(right.xOriginEngeExit_m),
    zOriginEngeExit_m(right.zOriginEngeExit_m),
    deltaBeginExit_m(right.deltaBeginExit_m),
    deltaEndExit_m(right.deltaEndExit_m),
    polyOrderExit_m(right.polyOrderExit_m),
    cosEntranceAngle_m(right.cosEntranceAngle_m),
    sinEntranceAngle_m(right.sinEntranceAngle_m),
    exitEdgeAngle_m(right.exitEdgeAngle_m),
    cosExitAngle_m(right.cosExitAngle_m),
    sinExitAngle_m(right.sinExitAngle_m) {

    setElType(isDipole);

}

Bend::Bend(const std::string &name):
    Component(name),
    bX_m(0.0),
    bY_m(0.0),
    designEnergy_m(0.0),
    messageHeader_m(" * "),
    pusher_m(),
    fileName_m(""),
    fieldmap_m(NULL),
    fast_m(false),
    angle_m(0.0),
    aperture_m(0.0),
    designRadius_m(0.0),
    fieldAmplitude_m(0.0),
    angleGreaterThanPi_m(false),
    entranceAngle_m(0.0),
    exitAngle_m(0.0),
    fieldIndex_m(0.0),
    elementEdge_m(0.0),
    startField_m(0.0),
    endField_m(0.0),
    reinitialize_m(false),
    recalcRefTraj_m(false),
    length_m(0.0),
    gap_m(0.0),
    refTrajMapSize_m(0),
    refTrajMapStepSize_m(0.0),
    entranceParameter1_m(0.0),
    entranceParameter2_m(0.0),
    entranceParameter3_m(0.0),
    exitParameter1_m(0.0),
    exitParameter2_m(0.0),
    exitParameter3_m(0.0),
    xOriginEngeEntry_m(0.0),
    zOriginEngeEntry_m(0.0),
    deltaBeginEntry_m(0.0),
    deltaEndEntry_m(0.0),
    polyOrderEntry_m(0),
    xExit_m(0.0),
    zExit_m(0.0),
    xOriginEngeExit_m(0.0),
    zOriginEngeExit_m(0.0),
    deltaBeginExit_m(0.0),
    deltaEndExit_m(0.0),
    polyOrderExit_m(0),
    cosEntranceAngle_m(1.0),
    sinEntranceAngle_m(0.0),
    exitEdgeAngle_m(0.0),
    cosExitAngle_m(1.0),
    sinExitAngle_m(0.0) {

    setElType(isDipole);

}

Bend::~Bend() {
}

/*
 * OPAL-T Methods.
 * ===============
 */

/*
 *  This function merely repackages the field arrays as type Vector_t and calls
 *  the equivalent method but with the Vector_t data types.
 */
bool Bend::apply(const size_t &i, const double &t, double E[], double B[]) {

    Vector_t Ev(0.0, 0.0, 0.0);
    Vector_t Bv(0.0, 0.0, 0.0);
    if(apply(RefPartBunch_m->R[i], RefPartBunch_m->get_rmean(), t, Ev, Bv))
        return true;

    E[0] = Ev(0);
    E[1] = Ev(1);
    E[2] = Ev(2);
    B[0] = Bv(0);
    B[1] = Bv(1);
    B[2] = Bv(2);

    return false;
}

bool Bend::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {

    if(designRadius_m > 0.0) {

        // Shift position to magnet frame.
        Vector_t X = RefPartBunch_m->X[i];
        X(2) += startField_m - elementEdge_m;

        /*
         * Add in transverse bend displacements. (ds is already
         * accounted for.)
         */
        X(0) -= dx_m;
        X(1) -= dy_m;

        // Get field from field map.
        Vector_t eField(0.0, 0.0, 0.0);
        Vector_t bField(0.0, 0.0, 0.0);
        CalculateMapField(X, eField, bField);
        bField *= fieldAmplitude_m;

        B(0) += bField(0);
        B(1) += bField(1);
        B(2) += bField(2);

    }

    return false;
}

bool Bend::apply(const Vector_t &R,
                  const Vector_t &centroid,
                  const double &t,
                  Vector_t &E,
                  Vector_t &B) {

    if(designRadius_m > 0.0) {

        int index = static_cast<int>
            (std::floor((R(2) - startField_m) / refTrajMapStepSize_m));

        if(index > 0 && index + 1 < refTrajMapSize_m) {

            // Find indices for position in pre-computed central trajectory map.
            double lever = (R(2) - startField_m) / refTrajMapStepSize_m - index;
            double x = (1.0 - lever) * refTrajMapX_m.at(index)
                + lever * refTrajMapX_m.at(index + 1);
            double y = (1.0 - lever) * refTrajMapY_m.at(index)
                + lever * refTrajMapY_m.at(index + 1);
            double z = (1.0 - lever) * refTrajMapZ_m.at(index)
                + lever * refTrajMapZ_m.at(index + 1);

            // Adjust position relative to pre-computed central trajectory map.
            Vector_t X(0.0, 0.0, 0.0);
            X(0) = R(0) + x;
            X(1) = R(1) + y;
            X(2) = z;

            Vector_t tempE(0.0, 0.0, 0.0);
            Vector_t tempB(0.0, 0.0, 0.0);
            Vector_t XInBendFrame = RotateToBendFrame(X);

            /*
             * Add in transverse bend displacements. (ds is already
             * accounted for.)
             */
            XInBendFrame(0) -= dx_m;
            XInBendFrame(1) -= dy_m;

            CalculateMapField(XInBendFrame, tempE, tempB);
            tempB = fieldAmplitude_m * RotateOutOfBendFrame(tempB);

            B(0) += tempB(0);
            B(1) += tempB(1);
            B(2) += tempB(2);

        }
    }

    return false;

}

bool Bend::bends() const {
    return true;
}

void Bend::finalise() {
    online_m = false;
}

void Bend::goOnline(const double &) {

    // Check if we need to reinitialize the bend field amplitude.
    if(reinitialize_m) {
        reinitialize_m = Reinitialize();
        recalcRefTraj_m = false;
    }

    /*
     * Always recalculate the reference trajectory on first call even
     * if we do not reinitialize the bend. The reference trajectory
     * has to be calculated at the same energy as the actual beam or
     * we do not get accurate values for the magnetic field in the output
     * file.
     */
    if(recalcRefTraj_m) {
        double angleX = 0.0;
        double angleY = 0.0;
        CalculateRefTrajectory(angleX, angleY);
        recalcRefTraj_m = false;
    }

    online_m = true;
}

void Bend::getDimensions(double &sBegin, double &sEnd) const {
    sBegin = startField_m;
    sEnd = endField_m;
}

void Bend::initialise(PartBunch *bunch,
                       double &startField,
                       double &endField,
                       const double &scaleFactor) {

    Inform msg(messageHeader_m.c_str(), *gmsg);

    if(InitializeFieldMap(msg)) {
        retrieveDesignEnergy(startField);

        SetupPusher(bunch);
        ReadFieldMap(msg);
        SetupBendGeometry(msg, startField, endField);

        double bendAngleX = 0.0;
        double bendAngleY = 0.0;
        CalculateRefTrajectory(bendAngleX, bendAngleY);
        recalcRefTraj_m = true;
        Print(msg, bendAngleX, bendAngleY);

        // Pass start and end of field to calling function.
        startField = startField_m;
        endField = endField_m;

    } else {
        ERRORMSG("There is something wrong with your field map \""
                 << fileName_m
                 << "\"");
        endField = startField - 1.0e-3;
    }
}

double Bend::GetBendAngle() const {
    return angle_m;
}

double Bend::GetBendRadius() const {
    return designRadius_m;
}

double Bend::GetEffectiveCenter() const {
    return elementEdge_m + designRadius_m * angle_m / 2.0;
}

double Bend::GetEffectiveLength() const {
    return designRadius_m * angle_m;
}

std::string Bend::GetFieldMapFN() const {
    return fileName_m;
}

double Bend::GetStartElement() const {
    return elementEdge_m;
}

void Bend::SetAngleGreaterThanPiFlag(bool angleGreaterThanPi) {
    angleGreaterThanPi_m = angleGreaterThanPi;
}

void Bend::SetAperture(double aperture) {
    aperture_m = std::abs(aperture);
}

void Bend::SetBendAngle(double angle) {
    angle_m = angle;
}

void Bend::SetBeta(double beta) {
    Orientation_m(1) = beta;
}

void Bend::SetDesignEnergy(double energy) {
    designEnergy_m = std::abs(energy);
}

void Bend::SetEntranceAngle(double entranceAngle) {
    entranceAngle_m = entranceAngle;
}

void Bend::setExitAngle(double exitAngle) {
    exitAngle_m = exitAngle;
}

void Bend::SetFieldAmplitude(double k0, double k0s) {
    bY_m = k0;
    bX_m = k0s;
}

void Bend::SetFieldMapFN(std::string fileName) {
    fileName_m = fileName;
}

void Bend::SetFullGap(double gap) {
    gap_m = std::abs(gap);
}

void Bend::SetK1(double k1) {
    fieldIndex_m = k1;
}

void Bend::SetLength(double length) {
    length_m = std::abs(length);
}

void Bend::SetRotationAboutZ(double rotation) {
    Orientation_m(2) = rotation;
}

void Bend::AdjustFringeFields(double ratio) {

    double delta = std::abs(entranceParameter1_m - entranceParameter2_m);
    entranceParameter1_m = entranceParameter2_m - delta * ratio;

    delta = std::abs(entranceParameter2_m - entranceParameter3_m);
    entranceParameter3_m = entranceParameter2_m + delta * ratio;

    delta = std::abs(exitParameter1_m - exitParameter2_m);
    exitParameter1_m = exitParameter2_m - delta * ratio;

    delta = std::abs(exitParameter2_m - exitParameter3_m);
    exitParameter3_m = exitParameter2_m + delta * ratio;

}

double Bend::CalculateBendAngle() {

    const double mass = RefPartBunch_m->getM();
    const double gamma = designEnergy_m / mass + 1.0;
    const double betaGamma = sqrt(pow(gamma, 2) - 1.0);
    const double beta = betaGamma / gamma;
    const double deltaT = RefPartBunch_m->getdT();

    // Integrate through field for initial angle.
    Vector_t X(0.0, 0.0, startField_m - elementEdge_m);
    Vector_t P(0.0, 0.0, betaGamma);
    double deltaS = 0.0;
    double bendLength = endField_m - startField_m;

    while(deltaS < bendLength) {

        X /= Vector_t(Physics::c * deltaT);
        pusher_m.push(X, P, deltaT);
        X *= Vector_t(Physics::c * deltaT);

        Vector_t eField(0.0, 0.0, 0.0);
        Vector_t bField(0.0, 0.0, 0.0);
        CalculateMapField(X, eField, bField);
        bField = fieldAmplitude_m * bField;

        X /= Vector_t(Physics::c * deltaT);
        pusher_m.kick(X, P, eField, bField, deltaT);

        pusher_m.push(X, P, deltaT);
        X *= Vector_t(Physics::c * deltaT);

        deltaS += deltaT * beta * Physics::c;

    }

    double angle =  -atan2(P(0), P(2));

    return angle;

}

void Bend::CalcCentralField(const Vector_t &R,
                             double deltaX,
                             double angle,
                             Vector_t &B) {

    double nOverRho = fieldIndex_m / designRadius_m;
    double expFactor = exp(-nOverRho * deltaX);
    double bxBzFactor = expFactor * nOverRho * R(1);

    B(0) = -bxBzFactor * cos(angle);
    B(1) = expFactor * (1.0 - pow(nOverRho * R(1), 2) / 2.0);
    B(2) = -bxBzFactor * sin(angle);

}

void Bend::CalcEngeFunction(double zNormalized,
                             const std::vector<double> &engeCoeff,
                             int polyOrder,
                             double &engeFunc,
                             double &engeFuncDeriv,
                             double &engeFuncSecDerivNorm) {

    double expSum = 0.0;
    double expSumDeriv = 0.0;
    double expSumSecDeriv = 0.0;

    if(polyOrder >= 2) {

        expSum = engeCoeff.at(0)
            + engeCoeff.at(1) * zNormalized;
        expSumDeriv = engeCoeff.at(1);

        for(int index = 2; index <= polyOrder; index++) {
            expSum += engeCoeff.at(index) * pow(zNormalized, index);
            expSumDeriv += index * engeCoeff.at(index)
                * pow(zNormalized, index - 1);
            expSumSecDeriv += index * (index - 1) * engeCoeff.at(index)
                * pow(zNormalized, index - 2);
        }

    } else if(polyOrder == 1) {

        expSum = engeCoeff.at(0)
            + engeCoeff.at(1) * zNormalized;
        expSumDeriv = engeCoeff.at(1);

    } else
        expSum = engeCoeff.at(0);

    double engeExp = exp(expSum);
    engeFunc = 1.0 / (1.0 + engeExp);

    if(!std::isnan(engeFunc)) {

        expSumDeriv /= gap_m;
        expSumSecDeriv /= pow(gap_m, 2);
        double engeExpDeriv = expSumDeriv * engeExp;
        double engeExpSecDeriv = (expSumSecDeriv + pow(expSumDeriv, 2))
            * engeExp;
        double engeFuncSq = pow(engeFunc, 2);

        engeFuncDeriv = -engeExpDeriv * engeFuncSq;
        if (std::isnan(engeFuncDeriv) || std::isinf(engeFuncDeriv))
            engeFuncDeriv = 0.0;

        engeFuncSecDerivNorm = -engeExpSecDeriv * engeFunc
            + 2.0 * pow(engeExpDeriv, 2)
            * engeFuncSq;
        if (std::isnan(engeFuncSecDerivNorm) || std::isinf(engeFuncSecDerivNorm))
            engeFuncSecDerivNorm = 0.0;

    } else {
        engeFunc = 0.0;
        engeFuncDeriv = 0.0;
        engeFuncSecDerivNorm = 0.0;

    }
}

void Bend::CalcEntranceFringeField(const Vector_t &REntrance,
                                    double deltaX,
                                    Vector_t &B) {

    double zNormalized = -REntrance(2) / gap_m;
    double engeFunc = 0.0;
    double engeFuncDeriv = 0.0;
    double engeFuncSecDerivNorm = 0.0;

    CalcEngeFunction(zNormalized,
                     engeCoeffsEntry_m,
                     polyOrderEntry_m,
                     engeFunc,
                     engeFuncDeriv,
                     engeFuncSecDerivNorm);

    double nOverRho = fieldIndex_m / designRadius_m;
    double expFactor = exp(-nOverRho * deltaX);
    double trigFactor = pow(nOverRho, 2) + engeFuncSecDerivNorm;

    double bXEntrance = -engeFunc * nOverRho * expFactor* REntrance(1);
    double bYEntrance = expFactor * engeFunc
        * (1.0  - trigFactor * pow(REntrance(1), 2) / 2.0);
    double bZEntrance = -expFactor * engeFuncDeriv * REntrance(1);

    B(0) = bXEntrance * cosEntranceAngle_m - bZEntrance * sinEntranceAngle_m;
    B(1) = bYEntrance;
    B(2) = bXEntrance * sinEntranceAngle_m + bZEntrance * cosEntranceAngle_m;
}

void Bend::CalcExitFringeField(const Vector_t &RExit, double deltaX, Vector_t &B) {

    double zNormalized = RExit(2) / gap_m;
    double engeFunc = 0.0;
    double engeFuncDeriv = 0.0;
    double engeFuncSecDerivNorm = 0.0;
    CalcEngeFunction(zNormalized,
                     engeCoeffsExit_m,
                     polyOrderExit_m,
                     engeFunc,
                     engeFuncDeriv,
                     engeFuncSecDerivNorm);

    double nOverRho = fieldIndex_m / designRadius_m;
    double expFactor = exp(-nOverRho * deltaX);
    double trigFactor = pow(nOverRho, 2) + engeFuncSecDerivNorm;

    double bXExit = -engeFunc * nOverRho * expFactor* RExit(1);
    double bYExit = expFactor * engeFunc
        * (1.0 - trigFactor * pow(RExit(1), 2) / 2.0);
    double bZExit = expFactor * engeFuncDeriv * RExit(1);

    B(0) = bXExit * cosExitAngle_m - bZExit * sinExitAngle_m;
    B(1) = bYExit;
    B(2) = bXExit * sinExitAngle_m + bZExit * cosExitAngle_m;
}

void Bend::CalculateMapField(const Vector_t &R, Vector_t &E, Vector_t &B) {

    E = Vector_t(0.0);
    B = Vector_t(0.0);

    //    Vector_t REntrance(0.0, 0.0, 0.0);
    //    Vector_t RExit(0.0, 0.0, 0.0);
    //    if (IsPositionInEntranceField(R, REntrance)) {
    //        CalcEntranceFringeField(REntrance, 0.0, B);
    //    } else if (IsPositionInExitField(R, RExit)) {
    //        CalcExitFringeField(RExit, 0.0, B);
    //    } else {
    //        CalcCentralField(R, 0.0, 0.0, B);
    //    }

    double deltaXEntrance = 0.0;
    double deltaXExit = 0.0;
    bool inEntranceRegion = InMagnetEntranceRegion(R, deltaXEntrance);
    bool inExitRegion = InMagnetExitRegion(R, deltaXExit);

    if(!inEntranceRegion && !inExitRegion) {

        double deltaX = 0.0;
        double angle = 0.0;
        if(InMagnetCentralRegion(R, deltaX, angle)) {
            Vector_t REntrance(0.0, 0.0, 0.0);
            Vector_t RExit(0.0, 0.0, 0.0);
            if(IsPositionInEntranceField(R, REntrance))
                CalcEntranceFringeField(REntrance, deltaX, B);
            else if(IsPositionInExitField(R, RExit))
                CalcExitFringeField(RExit, deltaX, B);
            else
                CalcCentralField(R, deltaX, angle, B);

        }

    } else if(inEntranceRegion && !inExitRegion) {

        Vector_t REntrance(0.0, 0.0, 0.0);
        if(IsPositionInEntranceField(R, REntrance)) {
            CalcEntranceFringeField(REntrance, deltaXEntrance, B);
        } else if(REntrance(2) > 0.0)
            CalcCentralField(R, deltaXEntrance, 0.0, B);

    } else if(!inEntranceRegion && inExitRegion) {

        Vector_t RExit(0.0, 0.0, 0.0);
        if(IsPositionInExitField(R, RExit)) {
            CalcExitFringeField(RExit, deltaXExit, B);
        } else if(RExit(2) < 0.0)
            CalcCentralField(R, deltaXExit, angle_m, B);

    } else if(inEntranceRegion && inExitRegion) {

        /*
         * This is an unusual condition and should only happen with
         * a sector magnet that bends more than 180 degrees. Here, we
         * have the possibility that the particle sees both the
         * entrance and exit fringe fields.
         */
        Vector_t BEntrance(0.0, 0.0, 0.0);
        Vector_t REntrance(0.0, 0.0, 0.0);
        if(IsPositionInEntranceField(R, REntrance))
            CalcEntranceFringeField(REntrance, deltaXEntrance, BEntrance);

        Vector_t BExit(0.0, 0.0, 0.0);
        Vector_t RExit(0.0, 0.0, 0.0);
        if(IsPositionInExitField(R, RExit))
            CalcExitFringeField(RExit, deltaXExit, BExit);

        B(0) = BEntrance(0) + BExit(0);
        B(1) = BEntrance(1) + BExit(1);
        B(2) = BEntrance(2) + BExit(2);

    }
}

void Bend::CalculateRefTrajectory(double &angleX, double &angleY) {

    const double mass = RefPartBunch_m->getM();
    const double gamma = designEnergy_m / mass + 1.;
    const double betaGamma = sqrt(gamma * gamma - 1.);
    const double dt = RefPartBunch_m->getdT();

    Vector_t X(0.0, 0.0, startField_m - elementEdge_m);
    Vector_t P(0.0, 0.0, betaGamma);

    if(!refTrajMapX_m.empty())
        refTrajMapX_m.clear();
    if(!refTrajMapY_m.empty())
        refTrajMapY_m.clear();
    if(!refTrajMapZ_m.empty())
        refTrajMapZ_m.clear();

    refTrajMapX_m.push_back(X(0));
    refTrajMapY_m.push_back(X(1));
    refTrajMapZ_m.push_back(X(2));

    refTrajMapStepSize_m = betaGamma / gamma * Physics::c * dt;
    double deltaS = 0.0;
    double bendLength = endField_m - startField_m;

    while(deltaS < bendLength) {

        X /= Vector_t(Physics::c * dt);
        pusher_m.push(X, P, dt);
        X *= Vector_t(Physics::c * dt);

        Vector_t eField(0.0, 0.0, 0.0);
        Vector_t bField(0.0, 0.0, 0.0);
        Vector_t XInBendFrame = RotateToBendFrame(X);

        /*
         * Add in transverse bend displacements. (ds is already
         * accounted for.)
         */
        XInBendFrame(0) -= dx_m;
        XInBendFrame(1) -= dy_m;

        CalculateMapField(XInBendFrame, eField, bField);
        bField = fieldAmplitude_m * RotateOutOfBendFrame(bField);

        X /= Vector_t(Physics::c * dt);
        pusher_m.kick(X, P, eField, bField, dt);

        pusher_m.push(X, P, dt);
        X *= Vector_t(Physics::c * dt);

        refTrajMapX_m.push_back(X(0));
        refTrajMapY_m.push_back(X(1));
        refTrajMapZ_m.push_back(X(2));

        deltaS += refTrajMapStepSize_m;

    }

    refTrajMapSize_m = refTrajMapX_m.size();

    if(std::abs(Orientation_m(2)) == Physics::pi / 2.0
       || Orientation_m(2) == 3.0 * Physics::pi / 2.0)
        angleX = 0.0;
    else
        angleX = -atan2(P(0), P(2));

    if(Orientation_m(2) == 0.0
       || Orientation_m(2) == Physics::pi)
        angleY = 0.0;
    else
        angleY = atan2(P(1), P(2));


    std::ofstream traj(std::string("data/") + getName() + std::string("_traj.dat"));
    traj.precision(8);

    double E1 = getEntranceAngle();
    double cosE1 = cos(0.5 * E1), sinE1 = sin(0.5 * E1);
    double zRotation = Orientation_m(2);
    Quaternion toStandard = (Quaternion(cosE1, sinE1 * Vector_t(0, 1, 0)) *
                             Quaternion(cos(0.5 * zRotation), sin(0.5 * zRotation) * Vector_t(0, 0, -1)));

    X = Vector_t(0.0, 0.0, startField_m - elementEdge_m);
    P = Vector_t(0.0, 0.0, betaGamma);

    deltaS = 0.0;

    while(deltaS < bendLength) {

        X /= Vector_t(Physics::c * dt);
        pusher_m.push(X, P, dt);
        X *= Vector_t(Physics::c * dt);

        Vector_t eField(0.0, 0.0, 0.0);
        Vector_t bField(0.0, 0.0, 0.0);
        Vector_t XInBendFrame = RotateToBendFrame(X);

        /*
         * Add in transverse bend displacements. (ds is already
         * accounted for.)
         */
        CalculateMapField(XInBendFrame, eField, bField);
        bField = fieldAmplitude_m * RotateOutOfBendFrame(bField);

        X /= Vector_t(Physics::c * dt);
        pusher_m.kick(X, P, eField, bField, dt);
        Vector_t R = toStandard.rotate(X * Physics::c * dt);
        traj << std::setw(16) << R(0)
             << std::setw(16) << R(2)
             << std::setw(16) << bField(1)
             << std::endl;


        pusher_m.push(X, P, dt);
        X *= Vector_t(Physics::c * dt);

        deltaS += refTrajMapStepSize_m;

    }

}

double Bend::EstimateFieldAdjustmentStep(double actualBendAngle,
                                          double mass,
                                          double betaGamma) {

    double amplitude1 = fieldAmplitude_m;
    double bendAngle1 = actualBendAngle;

    // Estimate field adjustment step.
    double effectiveLength = angle_m * designRadius_m;
    double fieldStep = (angle_m - bendAngle1) * betaGamma * mass / (2.0 * effectiveLength * Physics::c);
    if(pow(fieldAmplitude_m * effectiveLength * Physics::c / (betaGamma * mass), 2) < 1.0)
        fieldStep = (angle_m - bendAngle1) * betaGamma * mass / (2.0 * effectiveLength * Physics::c)
            * std::sqrt(1.0 - pow(fieldAmplitude_m * effectiveLength * Physics::c / (betaGamma * mass), 2));

    fieldStep *= amplitude1 / std::abs(amplitude1);

    return fieldStep;

}

void Bend::FindBendEffectiveLength(double startField, double endField) {

    /*
     * Use an iterative procedure to set the width of the
     * default field map for the defined field amplitude
     * and bend angle.
     */
    SetEngeOriginDelta(0.0);
    SetFieldCalcParam(false);
    SetFieldBoundaries(startField, endField);

    double actualBendAngle = CalculateBendAngle();
    double error = std::abs(actualBendAngle - angle_m);

    if(error > 1.0e-6) {

        double deltaStep = 0.0;
        if(std::abs(actualBendAngle) < std::abs(angle_m))
            deltaStep = -gap_m / 2.0;
        else
            deltaStep = gap_m / 2.0;

        double delta1 = 0.0;
        double bendAngle1 = actualBendAngle;

        double delta2 = deltaStep;
        SetEngeOriginDelta(delta2);
        SetFieldCalcParam(false);
        SetFieldBoundaries(startField, endField);
        double bendAngle2 = CalculateBendAngle();

        if(std::abs(bendAngle1) > std::abs(angle_m)) {
            while(std::abs(bendAngle2) > std::abs(angle_m)) {
                delta2 += deltaStep;
                SetEngeOriginDelta(delta2);
                SetFieldCalcParam(false);
                SetFieldBoundaries(startField, endField);
                bendAngle2 = CalculateBendAngle();
            }
        } else {
            while(std::abs(bendAngle2) < std::abs(angle_m)) {
                delta2 += deltaStep;
                SetEngeOriginDelta(delta2);
                SetFieldCalcParam(false);
                SetFieldBoundaries(startField, endField);
                bendAngle2 = CalculateBendAngle();
            }
        }

        // Now we should have the proper field map width bracketed.
        unsigned int iterations = 1;
        double delta = 0.0;
        error = std::abs(actualBendAngle - angle_m);
        while(error > 1.0e-6 && iterations < 100) {

            delta = (delta1 + delta2) / 2.0;
            SetEngeOriginDelta(delta);
            SetFieldCalcParam(false);
            SetFieldBoundaries(startField, endField);
            double newBendAngle = CalculateBendAngle();

            error = std::abs(newBendAngle - angle_m);

            if(error > 1.0e-6) {

                if(bendAngle1 - angle_m < 0.0) {

                    if(newBendAngle - angle_m < 0.0) {
                        bendAngle1 = newBendAngle;
                        delta1 = delta;
                    } else {
                        bendAngle2 = newBendAngle;
                        delta2 = delta;
                    }

                } else {

                    if(newBendAngle - angle_m < 0.0) {
                        bendAngle2 = newBendAngle;
                        delta2 = delta;
                    } else {
                        bendAngle1 = newBendAngle;
                        delta1 = delta;
                    }
                }
            }
            iterations++;
        }
    }
}

void Bend::FindBendStrength(double mass,
                             double gamma,
                             double betaGamma,
                             double charge) {

    /*
     * Use an iterative procedure to set the magnet field amplitude
     * for the defined bend angle.
     */
    double actualBendAngle = CalculateBendAngle();
    double fieldStep = EstimateFieldAdjustmentStep(actualBendAngle,
                                                   mass,
                                                   betaGamma);
    double amplitude1 = fieldAmplitude_m;
    double bendAngle1 = actualBendAngle;

    double amplitude2 = fieldAmplitude_m + fieldStep;
    fieldAmplitude_m = amplitude2;
    double bendAngle2 = CalculateBendAngle();

    if(std::abs(bendAngle1) > std::abs(angle_m)) {
        while(std::abs(bendAngle2) > std::abs(angle_m)) {
            amplitude2 += fieldStep;
            fieldAmplitude_m = amplitude2;
            bendAngle2 = CalculateBendAngle();
        }
    } else {
        while(std::abs(bendAngle2) < std::abs(angle_m)) {
            amplitude2 += fieldStep;
            fieldAmplitude_m = amplitude2;
            bendAngle2 = CalculateBendAngle();
        }
    }

    // Now we should have the proper field amplitude bracketed.
    unsigned int iterations = 1;
    double error = std::abs(actualBendAngle - angle_m);
    while(error > 1.0e-6 && iterations < 100) {

        fieldAmplitude_m = (amplitude1 + amplitude2) / 2.0;
        double newBendAngle = CalculateBendAngle();

        error = std::abs(newBendAngle - angle_m);

        if(error > 1.0e-6) {

            if(bendAngle1 - angle_m < 0.0) {

                if(newBendAngle - angle_m < 0.0) {
                    bendAngle1 = newBendAngle;
                    amplitude1 = fieldAmplitude_m;
                } else {
                    bendAngle2 = newBendAngle;
                    amplitude2 = fieldAmplitude_m;
                }

            } else {

                if(newBendAngle - angle_m < 0.0) {
                    bendAngle2 = newBendAngle;
                    amplitude2 = fieldAmplitude_m;
                } else {
                    bendAngle1 = newBendAngle;
                    amplitude1 = fieldAmplitude_m;
                }
            }
        }
        iterations++;
    }
}

bool Bend::FindIdealBendParameters(double chordLength) {

    double refMass = RefPartBunch_m->getM();
    double refGamma = designEnergy_m / refMass + 1.0;
    double refBetaGamma = sqrt(pow(refGamma, 2) - 1.0);
    double refCharge = RefPartBunch_m->getQ();

    if(angle_m != 0.0) {

        if(angle_m < 0.0) {
            // Negative angle is a positive bend rotated 180 degrees.
            angle_m = std::abs(angle_m);
            fieldIndex_m *= -1.0;
            Orientation_m(2) += Physics::pi;
        }
        designRadius_m = chordLength / (2.0 * std::sin(angle_m / 2.0));
        fieldAmplitude_m = (refCharge / std::abs(refCharge))
            * refBetaGamma * refMass
            / (Physics::c * designRadius_m);
        return true;

    } else if(bX_m == 0.0) {

        // Negative angle is a positive bend rotated 180 degrees.
        if((refCharge > 0.0 && bY_m < 0.0)
           || (refCharge < 0.0 && bY_m > 0.0)) {
            fieldIndex_m *= -1.0;
            Orientation_m(2) += Physics::pi;
        }

        fieldAmplitude_m = refCharge * std::abs(bY_m / refCharge);
        designRadius_m = std::abs(refBetaGamma * refMass / (Physics::c * fieldAmplitude_m));
        double bendAngle = 2.0 * std::asin(chordLength / (2.0 * designRadius_m));

        if(angleGreaterThanPi_m)
            bendAngle = 2.0 * Physics::pi - bendAngle;

        angle_m = bendAngle;

        return false;

    } else {

        Orientation_m(2) += atan2(bX_m, bY_m);
        if(refCharge < 0.0) {
            fieldIndex_m *= -1.0;
            Orientation_m(2) -= Physics::pi;
        }

        fieldAmplitude_m = refCharge
            * std::abs(sqrt(pow(bY_m, 2) + pow(bX_m, 2))
                       / refCharge);
        designRadius_m = std::abs(refBetaGamma * refMass / (Physics::c * fieldAmplitude_m));
        double bendAngle = 2.0 * std::asin(chordLength / (2.0 * designRadius_m));

        if(angleGreaterThanPi_m)
            bendAngle = 2.0 * Physics::pi - bendAngle;

        angle_m = bendAngle;

        return false;
    }
}

void Bend::FindReferenceExitOrigin(double &x, double &z) {

    /*
     * Find x,z coordinates of reference trajectory as it passes exit edge
     * of the bend magnet. This assumes an entrance position of (x,z) = (0,0).
     */
    if(angle_m <= Physics::pi / 2.0) {
        x = - designRadius_m * (1.0 - std::cos(angle_m));
        z = designRadius_m * std::sin(angle_m);
    } else if(angle_m <= Physics::pi) {
        x = -designRadius_m * (1.0 + std::sin(angle_m - Physics::pi / 2.0));
        z = designRadius_m * std::cos(angle_m - Physics::pi / 2.0);
    } else if(angle_m <= 3.0 * Physics::pi / 2.0) {
        x = -designRadius_m * (2.0 - std::cos(angle_m - Physics::pi));
        z = -designRadius_m * std::sin(angle_m - Physics::pi);
    } else {
        x = -designRadius_m * (1.0 - std::cos(angle_m - 3.0 * Physics::pi / 2.0));
        z = -designRadius_m * std::sin(angle_m - 3.0 * Physics::pi / 2.0);
    }
}

bool Bend::InitializeFieldMap(Inform &msg) {

    fieldmap_m = Fieldmap::getFieldmap(fileName_m, fast_m);

    if(fieldmap_m != NULL) {
        if(fileName_m != "1DPROFILE1-DEFAULT")
            return true;
        else
            return SetupDefaultFieldMap(msg);

    } else
        return false;

}

bool Bend::InMagnetCentralRegion(const Vector_t &R, double &deltaX, double &angle) {

    deltaX = sqrt(pow(R(2), 2) + pow(R(0) + designRadius_m, 2)) - designRadius_m;
    if(std::abs(deltaX) <= aperture_m / 2.0) {

        angle = atan2(R(2), R(0) + designRadius_m);
        return true;

    } else
        return false;

}

bool Bend::InMagnetEntranceRegion(const Vector_t &R, double &deltaX) {

    if(std::abs(R(0)) <= aperture_m / 2.0) {

        Vector_t RTransformed(0.0, R(1), 0.0);
        RTransformed(0) = (R(0) - xOriginEngeEntry_m) * cosEntranceAngle_m
            + (R(2) - zOriginEngeEntry_m) * sinEntranceAngle_m;
        RTransformed(2) = -(R(0) - xOriginEngeEntry_m) * sinEntranceAngle_m
            + (R(2) - zOriginEngeEntry_m) * cosEntranceAngle_m;

        if(RTransformed(2) <= 0.0) {
            deltaX = R(0);
            return true;
        } else
            return false;

    } else
        return false;

}

bool Bend::InMagnetExitRegion(const Vector_t &R, double &deltaX) {

    Vector_t RTransformed(0.0, R(1), 0.0);
    RTransformed(0) = (R(0) - xExit_m) * cosExitAngle_m
        + (R(2) - zExit_m) * sinExitAngle_m;
    RTransformed(2) = -(R(0) - xExit_m) * sinExitAngle_m
        + (R(2) - zExit_m) * cosExitAngle_m;

    if(RTransformed(2) >= 0.0) {

        deltaX = (R(0) - xExit_m) * cos(angle_m)
            + (R(2) - zExit_m) * sin(angle_m);
        if(std::abs(deltaX) <= aperture_m / 2.0)
            return true;
        else
            return false;

    } else
        return false;
}

bool Bend::IsPositionInEntranceField(const Vector_t &R, Vector_t &REntrance) {
    if (polyOrderEntry_m < 0)
        return false;

    REntrance(1) = R(1);

    REntrance(0) = (R(0) - xOriginEngeEntry_m) * cosEntranceAngle_m
        + (R(2) - zOriginEngeEntry_m) * sinEntranceAngle_m;
    REntrance(2) = -(R(0) - xOriginEngeEntry_m) * sinEntranceAngle_m
        + (R(2) - zOriginEngeEntry_m) * cosEntranceAngle_m;

    if(REntrance(2) >= -deltaBeginEntry_m && REntrance(2) <= deltaEndEntry_m)
        return true;
    else
        return false;
}

bool Bend::IsPositionInExitField(const Vector_t &R, Vector_t &RExit) {
    if (polyOrderExit_m < 0)
        return false;

    RExit(1) = R(1);

    RExit(0) = (R(0) - xOriginEngeExit_m) * cosExitAngle_m
        + (R(2) - zOriginEngeExit_m) * sinExitAngle_m;
    RExit(2) = -(R(0) - xOriginEngeExit_m) * sinExitAngle_m
        + (R(2) - zOriginEngeExit_m) * cosExitAngle_m;

    if(RExit(2) >= -deltaBeginExit_m && RExit(2) <= deltaEndExit_m)
        return true;
    else
        return false;

}

void Bend::Print(Inform &msg, double bendAngleX, double bendAngleY) {
    msg << level2 << "\n"
        << "Start of field map:      "
        << startField_m
        << " m (in s coordinates)"
        << "\n";
    msg << "End of field map:        "
        << endField_m
        << " m (in s coordinates)"
        << "\n";
    msg << "Entrance edge of magnet: "
        << elementEdge_m
        << " m (in s coordinates)"
        << "\n";
    msg << "\n"
        << "Reference Trajectory Properties"
        << "\n"
        << "==============================="
        << "\n\n";
    msg << "Bend angle magnitude:    "
        << angle_m
        << " rad ("
        << angle_m * 180.0 / Physics::pi
        << " degrees)"
        << "\n";
    msg << "Entrance edge angle:     "
        << entranceAngle_m
        << " rad ("
        << entranceAngle_m * 180.0 / Physics::pi
        << " degrees)"
        << "\n";
    msg << "Exit edge angle:         "
        << exitAngle_m
        << " rad ("
        << exitAngle_m * 180.0 / Physics::pi
        << " degrees)"
        << "\n";
    msg << "Bend design radius:      "
        << designRadius_m
        << " m"
        << "\n";
    msg << "Bend design energy:      "
        << designEnergy_m
        << " eV"
        << "\n";
    msg << "\n"
        << "Bend Field and Rotation Properties"
        << "\n"
        << "=================================="
        << "\n\n";
    msg << "Field amplitude:         "
        << fieldAmplitude_m
        << " T"
        << "\n";
    msg << "Field index:  "
        << fieldIndex_m
        << "\n";
    msg << "Rotation about x axis:   "
        << Orientation_m(1)
        << " rad ("
        << Orientation_m(1) * 180.0 / Physics::pi
        << " degrees)"
        << "\n";
    msg << "Rotation about y axis:   "
        << Orientation_m(0)
        << " rad ("
        << Orientation_m(0) * 180.0 / Physics::pi
        << " degrees)"
        << "\n";
    msg << "Rotation about z axis:   "
        << Orientation_m(2)
        << " rad ("
        << Orientation_m(2) * 180.0 / Physics::pi
        << " degrees)"
        << "\n";
    msg << "\n"
        << "Reference Trajectory Properties Through Bend Magnet with Fringe Fields"
        << "\n"
        << "======================================================================"
        << "\n\n";
    msg << "Reference particle is bent: "
        << bendAngleX
        << " rad ("
        << bendAngleX * 180.0 / Physics::pi
        << " degrees) in x plane"
        << "\n";
    msg << "Reference particle is bent: "
        << bendAngleY
        << " rad ("
        << bendAngleY * 180.0 / Physics::pi
        << " degrees) in y plane"
        << "\n" << endl;
}

void Bend::ReadFieldMap(Inform &msg) {
    msg << level2 << getName() << " using file ";
    fieldmap_m->getInfo(&msg);

    Fieldmap::readMap(fileName_m);
    fieldmap_m->Get1DProfile1EntranceParam(entranceParameter1_m,
                                           entranceParameter2_m,
                                           entranceParameter3_m);
    fieldmap_m->Get1DProfile1ExitParam(exitParameter1_m,
                                       exitParameter2_m,
                                       exitParameter3_m);
    SetGapFromFieldMap();
    fieldmap_m->Get1DProfile1EngeCoeffs(engeCoeffsEntry_m,
                                        engeCoeffsExit_m);
    polyOrderEntry_m = engeCoeffsEntry_m.size() - 1;
    polyOrderExit_m = engeCoeffsExit_m.size() - 1;

}

bool Bend::Reinitialize() {

    designEnergy_m = RefPartBunch_m->get_meanEnergy() * 1.0e6;
    SetBendStrength();
    double bendAngleX = 0.0;
    double bendAngleY = 0.0;
    CalculateRefTrajectory(bendAngleX, bendAngleY);

    INFOMSG(level2 << "Bend design energy changed to: "
            << designEnergy_m * 1.0e-6
            << " MeV"
            << endl);
    Print(*Ippl::Info, bendAngleX, bendAngleY);

    return false;
}

Vector_t Bend::RotateOutOfBendFrame(const Vector_t &X) {

    /*
     * Rotate vector out of the bend's local coordinate system back to
     * the lab frame.
     *
     * 1) Rotate about the x axis by angle negative Orientation_m(1).
     * 2) Rotate about the y axis by angle Orientation_m(0).
     * 3) Rotate about the z axis by angle Orientation_m(3).
     */

    double sina = sin(Orientation_m(0));
    double cosa = cos(Orientation_m(0));
    double sinb = sin(Orientation_m(1));
    double cosb = cos(Orientation_m(1));
    double sinc = sin(Orientation_m(2));
    double cosc = cos(Orientation_m(2));

    Vector_t temp(0.0, 0.0, 0.0);

    temp(0) = (cosa * cosc) *                       X(0)
        + (-sina * sinb * cosc - cosb * sinc) * X(1)
        + (sina * cosb * cosc - sinb * sinc)  * X(2);
    temp(1) = (cosa * sinc) *                       X(0)
        + (-sina * sinb * sinc + cosb * cosc) * X(1)
        + (sina * cosb * sinc + sinb * cosc)  * X(2);
    temp(2) =   -sina *                               X(0)
        + (-cosa * sinb) *                      X(1)
        + (cosa * cosb) *                       X(2);

    return temp;

}

Vector_t Bend::RotateToBendFrame(const Vector_t &X) {

    /*
     * Rotate vector to the bend's local coordinate system.
     *
     * 1) Rotate about the z axis by angle negative Orientation_m(2).
     * 2) Rotate about the y axis by angle negative Orientation_m(0).
     * 3) Rotate about the x axis by angle Orientation_m(1).
     */

    double sina = sin(Orientation_m(0));
    double cosa = cos(Orientation_m(0));
    double sinb = sin(Orientation_m(1));
    double cosb = cos(Orientation_m(1));
    double sinc = sin(Orientation_m(2));
    double cosc = cos(Orientation_m(2));

    Vector_t temp(0.0, 0.0, 0.0);

    temp(0) = (cosa * cosc) * X(0)
        + (cosa * sinc) * X(1)
        -  sina *         X(2);
    temp(1) = (-cosb * sinc - sina * sinb * cosc) * X(0)
        + (cosb * cosc - sina * sinb * sinc)  * X(1)
        -  cosa * sinb *                        X(2);
    temp(2) = (-sinb * sinc + sina * cosb * cosc) * X(0)
        + (sinb * cosc + sina * cosb * sinc)  * X(1)
        + cosa * cosb *                         X(2);

    return temp;
}

void Bend::SetBendEffectiveLength(double startField, double endField) {

    // Find initial angle.
    double actualBendAngle = CalculateBendAngle();

    // Adjust field map to match bend angle.
    double error = std::abs(actualBendAngle - angle_m);
    if(error > 1.0e-6)
        FindBendEffectiveLength(startField, endField);

}

void Bend::SetBendStrength() {

    // Estimate bend field magnitude.
    double mass = RefPartBunch_m->getM();
    double gamma = designEnergy_m / mass + 1.0;
    double betaGamma = sqrt(pow(gamma, 2) - 1.0);
    double charge = RefPartBunch_m->getQ();

    fieldAmplitude_m = (charge / std::abs(charge)) * betaGamma * mass
        / (Physics::c * designRadius_m);

    // Find initial angle.
    double actualBendAngle = CalculateBendAngle();

    // Search for angle if initial guess is not good enough.
    double error = std::abs(actualBendAngle - angle_m);
    if(error > 1.0e-6)
        FindBendStrength(mass, gamma, betaGamma, charge);

}

void Bend::SetEngeOriginDelta(double delta) {
    /*
     * This function is used to shift the perpendicular distance of the
     * entrance and exit Enge function origins with respect to the entrance
     * and exit points in the magnet. A positive delta shifts them towards
     * the center of the magnet.
     */
    entranceParameter1_m = delta - std::abs(entranceParameter1_m
                                            - entranceParameter2_m);
    entranceParameter3_m = delta + std::abs(entranceParameter2_m
                                            - entranceParameter3_m);
    entranceParameter2_m = delta;

    exitParameter1_m = -delta - std::abs(exitParameter1_m - exitParameter2_m);
    exitParameter3_m = -delta + std::abs(exitParameter2_m - exitParameter3_m);
    exitParameter2_m = -delta;
}

void Bend::SetFieldCalcParam(bool chordLengthFromMap) {

    cosEntranceAngle_m = cos(entranceAngle_m);
    sinEntranceAngle_m = sin(entranceAngle_m);

    deltaBeginEntry_m = std::abs(entranceParameter1_m - entranceParameter2_m);
    deltaEndEntry_m = std::abs(entranceParameter2_m - entranceParameter3_m);

    exitEdgeAngle_m = angle_m - exitAngle_m;
    cosExitAngle_m = cos(exitEdgeAngle_m);
    sinExitAngle_m = sin(exitEdgeAngle_m);

    deltaBeginExit_m = std::abs(exitParameter1_m - exitParameter2_m);
    deltaEndExit_m = std::abs(exitParameter2_m - exitParameter3_m);

    if(chordLengthFromMap) {
        /*
         * The magnet chord length is taken from this field map. In this case,
         * we assume that the origin points for the entrance and exit Enge
         * functions correspond to the physical edges of the bend magnet. These
         * are the points where the ideal reference particle trajectory intercepts
         * the entrance and exit edges.
         */
        xOriginEngeEntry_m = 0.0;
        zOriginEngeEntry_m = 0.0;
        xOriginEngeExit_m = xExit_m;
        zOriginEngeExit_m = zExit_m;

    } else {
        /*
         * In this case, the chord length of the magnet is set in the input file
         * and is not taken from the field map file. This allows us to set the
         * origin points for the entrance and exit Enge functions some perpendicular
         * distance away from the physical edges of the magnet.
         */
        xOriginEngeEntry_m = -entranceParameter2_m * sinEntranceAngle_m;
        zOriginEngeEntry_m = entranceParameter2_m * cosEntranceAngle_m;

        xOriginEngeExit_m = xExit_m - exitParameter2_m * sinExitAngle_m;
        zOriginEngeExit_m = zExit_m + exitParameter2_m * cosExitAngle_m;

    }
}

void Bend::SetGapFromFieldMap() {

    if(gap_m <= 0.0)
        gap_m = fieldmap_m->GetFieldGap();
    else if(gap_m != fieldmap_m->GetFieldGap())
        AdjustFringeFields(gap_m / fieldmap_m->GetFieldGap());

}

bool Bend::SetupBendGeometry(Inform &msg, double &startField, double &endField) {

    double chordLength = 0.0;
    bool chordLengthFromMap = false;
    if(!FindChordLength(msg, chordLength, chordLengthFromMap))
        return false;

    if(TreatAsDrift(msg, chordLength)) {
        startField_m = startField;
        endField_m = startField + chordLength;
        return true;
    }

    reinitialize_m = FindIdealBendParameters(chordLength);
    if (getType() == RBEND) {
        SetEntranceAngle(getEntranceAngle());
    }
    FindReferenceExitOrigin(xExit_m, zExit_m);

    /*
     * Set field map geometry.
     */
    if(aperture_m <= 0.0)
        aperture_m = designRadius_m / 2.0;
    SetFieldCalcParam(chordLengthFromMap);

    /*
     * If we are using the default field map, then the origins for the
     * Enge functions will be shifted so that we get the desired bend
     * angle for the given bend strength. (We match the effective length
     * of our field map to the ideal bend length.)
     *
     * If we are not using the default field map, we assume it cannot
     * change so we either leave everything alone (if the user defines
     * the bend strength) or we adjust the bend field to get the right
     * angle.
     */
    elementEdge_m = startField + ds_m;
    SetFieldBoundaries(startField, endField);

    if(fileName_m != "1DPROFILE1-DEFAULT") {
        if(reinitialize_m)
            SetBendStrength();
    } else {
        SetBendEffectiveLength(startField, endField);
    }

    startField = startField_m;
    endField = endField_m;
    return true;

}

bool Bend::SetupDefaultFieldMap(Inform &msg) {

    if(length_m <= 0.0) {
        ERRORMSG("If using \"1DPROFILE1-DEFAULT\" field map you must set the "
                 "chord length of the bend using the L attribute in the OPAL "
                 "input file."
                 << endl);
        return false;
    } else {
        msg << level2;
        fieldmap_m->getInfo(&msg);
        return true;
    }

}

void Bend::SetFieldBoundaries(double startField, double endField) {

    startField_m = startField - deltaBeginEntry_m / cos(entranceAngle_m) + ds_m;
    endField_m = startField + angle_m * designRadius_m
        + deltaEndExit_m / cos(exitAngle_m)
        + ds_m;

}

void Bend::SetupPusher(PartBunch *bunch) {

    RefPartBunch_m = bunch;
    pusher_m.initialise(bunch->getReference());

}

bool Bend::TreatAsDrift(Inform &msg, double chordLength) {
    if(designEnergy_m <= 0.0) {
        WARNMSG("Warning: bend design energy is zero. Treating as drift."
                << endl);
        designRadius_m = 0.0;
        return true;
    } else if(angle_m == 0.0) {

        double refMass = RefPartBunch_m->getM();
        double refGamma = designEnergy_m / refMass + 1.0;
        double refBetaGamma = sqrt(pow(refGamma, 2) - 1.0);

        double amplitude = std::abs(fieldAmplitude_m);
        double radius = std::abs(refBetaGamma * refMass / (Physics::c * amplitude));
        double sinArgument = chordLength / (2.0 * radius);

        if(std::abs(sinArgument) > 1.0) {
            WARNMSG("Warning: bend strength and defined reference trajectory "
                    << "chord length are not consistent. Treating bend as drift."
                    << endl);
            designRadius_m = 0.0;
            return true;
        } else
            return false;

    } else if(angle_m == 0.0 &&
              pow(bX_m, 2) + pow(bY_m, 2) == 0.0) {

        WARNMSG("Warning bend angle/strength is zero. Treating as drift."
                << endl);
        designRadius_m = 0.0;
        return true;

    } else
        return false;
}

void Bend::retrieveDesignEnergy(double startField) {
    energyEvolution_t::iterator it = OpalData::getInstance()->getFirstEnergyData();
    energyEvolution_t::iterator end = OpalData::getInstance()->getLastEnergyData();
    for (; it != end; ++ it) {
        if ((*it).first > startField) break;
        designEnergy_m = (*it).second;
    }
}
