// ------------------------------------------------------------------------
// $RCSfile: SBend.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Definitions for class: SBend
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

#include "Algorithms/PartPusher.h"
#include "AbsBeamline/SBend.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.hh"
#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class SBend
// ------------------------------------------------------------------------

SBend::SBend():
    Component(),
    pusher_m(),
    field_orientation_m(1.0, 0.0),
    filename_m(""),
    fieldmap_m(NULL),
    fast_m(false),
    reinitialize_m(false),
    startField_m(0.0),
    endField_m(0.0),
    length_m(0.0),
    gap_m(0.0),
    ElementEdge_m(0.0),
    startElement_m(0.0),
    amplitude_m(0.0),
    gradient_m(0.0),
    angle_m(0.0),
    design_energy_m(0.0),
    alpha_m(0.0),
    exitAngle_m(0.0),
    sin_face_alpha_m(0.0),
    cos_face_alpha_m(1.0),
    tan_face_alpha_m(0.0),
    sin_face_beta_m(0.0),
    cos_face_beta_m(1.0),
    tan_face_beta_m(0.0),
    map_m(NULL),
    map_size_m(0),
    map_step_size_m(0.0),
    R_m(0.0),
    effectiveLength_m(0.0),
    effectiveCenter_m(0.0),
    effectiveStart_m(0.0) {
    setElType(isDipole);
}

SBend::SBend(const SBend &right):
    Component(right),
    pusher_m(right.pusher_m),
    field_orientation_m(right.field_orientation_m),
    filename_m(right.filename_m),
    fieldmap_m(right.fieldmap_m),
    fast_m(right.fast_m),
    reinitialize_m(right.reinitialize_m),
    startField_m(right.startField_m),
    endField_m(right.endField_m),
    length_m(right.length_m),
    gap_m(right.gap_m),
    ElementEdge_m(right.ElementEdge_m),
    startElement_m(right.startElement_m),
    amplitude_m(right.amplitude_m),
    gradient_m(right.gradient_m),
    angle_m(right.angle_m),
    design_energy_m(right.design_energy_m),
    alpha_m(right.alpha_m),
    exitAngle_m(right.exitAngle_m),
    sin_face_alpha_m(right.sin_face_alpha_m),
    cos_face_alpha_m(right.cos_face_alpha_m),
    tan_face_alpha_m(right.tan_face_alpha_m),
    sin_face_beta_m(right.sin_face_beta_m),
    cos_face_beta_m(right.cos_face_beta_m),
    tan_face_beta_m(right.tan_face_beta_m),
    map_m(NULL),
    map_size_m(right.map_size_m),
    map_step_size_m(right.map_step_size_m),
    R_m(right.R_m),
    effectiveLength_m(right.effectiveLength_m),
    effectiveCenter_m(right.effectiveCenter_m),
    effectiveStart_m(right.effectiveStart_m) {
    setElType(isDipole);
    if(map_size_m > 0) {
        map_m = new double[3 * (map_size_m + 1)];
        for(int i = 0; i < 3 * (map_size_m + 1); ++i)
            map_m[i] = right.map_m[i];
    }
}


SBend::SBend(const std::string &name):
    Component(name),
    pusher_m(),
    field_orientation_m(1.0, 0.0),
    filename_m(""),
    fieldmap_m(NULL),
    fast_m(false),
    reinitialize_m(false),
    startField_m(0.0),
    endField_m(0.0),
    length_m(0.0),
    gap_m(0.0),
    ElementEdge_m(0.0),
    startElement_m(0.0),
    amplitude_m(0.0),
    gradient_m(0.0),
    angle_m(0.0),
    design_energy_m(0.0),
    alpha_m(0.0),
    exitAngle_m(0.0),
    sin_face_alpha_m(0.0),
    cos_face_alpha_m(1.0),
    tan_face_alpha_m(0.0),
    sin_face_beta_m(0.0),
    cos_face_beta_m(1.0),
    tan_face_beta_m(0.0),
    map_m(NULL),
    map_size_m(0),
    map_step_size_m(0.0),
    R_m(0.0),
    effectiveLength_m(0.0),
    effectiveCenter_m(0.0),
    effectiveStart_m(0.0) {
    setElType(isDipole);
}

SBend::~SBend() {
    if(map_m) {
        delete[] map_m;
        map_m = NULL;
    }
}

void SBend::accept(BeamlineVisitor &visitor) const {
    visitor.visitSBend(*this);
}

/*
 * OPAL-MAP methods
 * ================
 */
double SBend::getNormalComponent(int n) const {
    return getField().getNormalComponent(n);
}

double SBend::getSkewComponent(int n) const {
    return getField().getSkewComponent(n);
}

void SBend::setNormalComponent(int n, double v) {
    getField().setNormalComponent(n, v);
}

void SBend::setSkewComponent(int n, double v) {
    getField().setSkewComponent(n, v);
}

/*
 * OPAL-T Methods.
 * ===============
 */

/*
 *  This function merely repackages the field arrays as type Vector_t and calls
 *  the equivalent method but with the Vector_t data types.
 */

bool SBend::apply(const size_t &i, const double &t, double E[], double B[]) {
    Vector_t Ev(0, 0, 0), Bv(0, 0, 0);
    if(apply(RefPartBunch_m->R[i], RefPartBunch_m->get_rmean(), t, Ev, Bv)) return true;

    E[0] = Ev(0);
    E[1] = Ev(1);
    E[2] = Ev(2);
    B[0] = Bv(0);
    B[1] = Bv(1);
    B[2] = Bv(2);

    return false;
}

bool SBend::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {

    // Check if we need to reinitialize the bend field amplitude.
    if(reinitialize_m)
        reinitialize_m = reinitialize();

    // Get field from field map.
    const Vector_t &X = RefPartBunch_m->X[i];
    double bX = 0.0;
    double bY = 0.0;
    double bZ = 0.0;
    calculateMapField(X, bX, bY, bZ);

    B(0) += bX;
    B(1) += bY;
    B(2) += bZ;

    //    Vector_t strength(0.0);
    //    Vector_t info(0.0);
    //    fieldmap_m->getFieldstrength(X, strength, info);
    //
    //    double dd, dx, dz, rho;
    //    const double &k34 = info(2);
    //
    //    if(k34 > 0) {
    //        dx = X(0) + R_m * cos_face_alpha_m;
    //        dz = X(2) - R_m * sin_face_alpha_m;
    //        rho = sqrt(dx * dx + dz * dz);
    //        dd = 1.0 - rho / R_m;
    //    } else {
    //        dx = -X(0) + R_m * cos_face_alpha_m;
    //        dz = -X(2) + R_m * sin_face_alpha_m;
    //        rho = sqrt(dx * dx + dz * dz);
    //        dd = -1.0 + rho / R_m;
    //    }
    //
    //    if(info(0) > 0.99) {
    //
    //        B(1) += amplitude_m * (strength(0) - strength(2) / 2.0 * pow(X(1), 2.0)) * (1.0 - gradient_m * dd);
    //        double bX = amplitude_m * (strength(0) - strength(2) / 2.0 * pow(X(1), 2.0)) * gradient_m * X(1) / (rho * R_m);
    //
    //        if(info(1) > 0.99) {
    //            B(0) -= amplitude_m * strength(1) * X(1) * sin(angle_m - alpha_m - exitAngle_m) + bX * dx;
    //            B(2) += amplitude_m * strength(1) * X(1) * cos(angle_m - alpha_m - exitAngle_m) + bX * dz;
    //        } else {
    //            B(0) += bX * dx;
    //            B(2) -= amplitude_m * strength(1) * X(1) + bX * dz;
    //        }
    //
    //    } else if(fabs(info(0)) < 0.01) {
    //        B(0) += amplitude_m * gradient_m * X(1) * dx / (rho * R_m);
    //        B(1) += amplitude_m * (1.0 - gradient_m * dd);
    //        B(2) += amplitude_m * gradient_m * X(1) * dz / (rho * R_m);
    //    }
    //    if(info(0) > 0.99) {
    //        B(1) +=  amplitude_m * (strength(0) - strength(2) / 2. * X(1) * X(1)) * (1. - gradient_m * dd);
    //        double Bx = amplitude_m * (strength(0) - strength(2) / 2. * X(1) * X(1)) * gradient_m * X(1) / (rho * R_m);
    //        double Bz = amplitude_m * strength(1) * X(1);
    //        if(info(1) > 0.99) {
    //            B(0) += - Bz * k34 / sqrt(1. + k34 * k34) + Bx * dx;
    //            B(2) += Bz / sqrt(1. + k34 * k34) + Bx * dz;
    //        } else {
    //            B(0) += Bx * dx;
    //            B(2) += -Bz + Bx * dz;
    //        }
    //    } else if(info(0) < 0.01) {
    //        B(0) += amplitude_m * gradient_m * X(1) * dx / (rho * R_m);
    //        B(1) += amplitude_m * (1. - gradient_m * dd);
    //        B(2) += amplitude_m * gradient_m * X(1) * dz / (rho * R_m);
    //    }

    return false;
}

bool SBend::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {

    int index = (int)floor((R(2) - startField_m) / map_step_size_m);
    if(index > 0 && index + 1 < map_size_m) {

        //        // Find indices for position in pre-computed central trajectory map.
        //        double lever = (R(2) - startField_m) / map_step_size_m - index;
        //        double z = (1. - lever) * map_m[3 * index + 2] + lever * map_m[3 * index + 5];
        //        double x = (1. - lever) * map_m[3 * index] + lever * map_m[3 * index + 3];
        //        double y = (1. - lever) * map_m[3 * index + 1] + lever * map_m[3 * index + 4];

        // Rotate x and y to the the bend's local coordinate system.
        //
        // 1) Rotate about the z axis by angle negative ori(2).
        // 2) Rotate about the y axis by angle negative ori(0).
        // 3) Rotate about the x axis by angle ori(1).

        const double sina = sin(Orientation_m(0));
        const double cosa = cos(Orientation_m(0));
        const double sinb = sin(Orientation_m(1));
        const double cosb = cos(Orientation_m(1));
        const double sinc = sin(Orientation_m(2));
        const double cosc = cos(Orientation_m(2));

        Vector_t X(0.0);
        //        X(0) = (cosa * cosc) * (x + R(0)) + (cosa * sinc) * (y + R(1)) - sina *        R(2);
        //        X(1) = (-cosb * sinc - sina * sinb * cosc) * (x + R(0)) + (cosb * cosc - sina * sinb * sinc) * (y + R(1)) - cosa * sinb * R(2);
        //        X(2) = z;

        X(0) = (cosa * cosc) *                       R(0) + (cosa * sinc) *                      R(1) - sina *        R(2);
        X(1) = (-cosb * sinc - sina * sinb * cosc) * R(0) + (cosb * cosc - sina * sinb * sinc) * R(1) - cosa * sinb * R(2);
        X(2) = (-sinb * sinc + sina * cosb * cosc) * R(0) + (sinb * cosc + sina * cosb * sinc) * R(1) + cosa * cosb * R(2);

        // Find indices for position in pre-computed central trajectory map.
        double lever = (R(2) - startField_m) / map_step_size_m - index;
        double z = (1. - lever) * map_m[3 * index + 2] + lever * map_m[3 * index + 5];
        double x = (1. - lever) * map_m[3 * index] + lever * map_m[3 * index + 3];
        double y = (1. - lever) * map_m[3 * index + 1] + lever * map_m[3 * index + 4];

        // Adjust position relative to pre-computed central trajectory map.
        X(0) += x;
        X(1) += y;
        X(2) = z;

        Vector_t strength(0.0), info(0.0);

        fieldmap_m->getFieldstrength(X, strength, info);
        Vector_t tempB(0.0);

        if(info(0) > 0.99) {
            tempB(1) += amplitude_m * (strength(0) - strength(2) / 2.0 * pow(X(1), 2.0));

            if(info(1) > 0.99) {
                tempB(0) -= amplitude_m * strength(1) * X(1) * sin(exitAngle_m);
                tempB(2) += amplitude_m * strength(1) * X(1) * cos(exitAngle_m);
            } else {
                tempB(2) -= amplitude_m * strength(1) * X(1);
            }
        } else if(fabs(info(0)) < 0.01)
            tempB(1) += amplitude_m;

        // double dd, dx, dz, rho;
        //        const double &k34 = info(2);
        //
        //        if(k34 > 0) {
        //            dx = X(0) + R_m * cos_face_alpha_m;
        //            dz = X(2) - R_m * sin_face_alpha_m;
        //            rho = sqrt(dx * dx + dz * dz);
        //            dd = 1.0 - rho / R_m;
        //        } else {
        //            dx = -X(0) + R_m * cos_face_alpha_m;
        //            dz = -X(2) + R_m * sin_face_alpha_m;
        //            rho = sqrt(dx * dx + dz * dz);
        //            dd = -1.0 + rho / R_m;
        //        }
        //
        //        if(info(0) > 0.99) {
        //            tempB(1) =  amplitude_m * (strength(0) - strength(2) / 2. * X(1) * X(1)) * (1. - gradient_m * dd);
        //            double Bx = amplitude_m * (strength(0) - strength(2) / 2. * X(1) * X(1)) * gradient_m * X(1) / (rho * R_m);
        //            double Bz = amplitude_m * strength(1) * X(1);
        //            if(info(1) > 0.99) {
        //                tempB(0) = - Bz * k34 / sqrt(1. + k34 * k34) + Bx * dx;
        //                tempB(2) = Bz / sqrt(1. + k34 * k34) + Bx * dz;
        //            } else {
        //                tempB(0) = Bx * dx;
        //                tempB(2) = -Bz + Bx * dz;
        //            }
        //        } else if(info(0) < 0.01) {
        //            tempB(0) = amplitude_m * gradient_m * X(1) * dx / (rho * R_m);
        //            tempB(1) = amplitude_m * (1. - gradient_m * dd);
        //            tempB(2) = amplitude_m * gradient_m * X(1) * dz / (rho * R_m);
        //        }

        // Rotate field out of the bend's local coordinate system and back to lab frame.
        //
        // 1) Rotate about the x axis by angle ori(1).
        // 2) Rotate about the y axis by angle ori(0).
        // 3) Rotate about the z axis by angle negative ori(3).

        B(0) +=  cosa * cosc * tempB(0) + (-sina * sinb * cosc - cosb * sinc) * tempB(1) + (sina * cosb * cosc - sinb * sinc) * tempB(2);
        B(1) +=  cosa * sinc * tempB(0) + (-sina * sinb * sinc + cosb * cosc) * tempB(1) + (sina * cosb * sinc + sinb * cosc) * tempB(2);
        B(2) += -sina *        tempB(0) + (-cosa * sinb) * tempB(1) + (cosa * cosb) * tempB(2);
    }


    return false;

}

/// Does initial setup of the bend.
void SBend::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {

    Inform msg("SBend ");

    double zBegin = 0.0;
    double zEnd = 0.0;
    double rBegin = 0.0;
    double rEnd = 0.0;

    startElement_m = startField;

    RefPartBunch_m = bunch;
    pusher_m.initialise(bunch->getReference());

    // Initialize field map.
    fieldmap_m = Fieldmap::getFieldmap(filename_m, fast_m);
    fieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);

    if((fieldmap_m != NULL) && (zEnd > zBegin)) {

        // Read in field map.
        msg << getName() << " using file ";
        fieldmap_m->getInfo(&msg);
        Fieldmap::readMap(filename_m);

        // Check that the design energy is greater then zero.
        if(design_energy_m <= 0.0) {
            msg << "The bend must have a design energy greater than zero set in the input file." << endl;
            return;
        }

        // If using default field map, set length and gap.
        if(filename_m == "1DPROFILE1-DEFAULT") {
            if(gap_m <= 0.0 || length_m <= 0.0) {
                msg << "If using \"1DPROFILE1-DEFAULT\" field map you must set the full magnet gap, GAP, "
                    << " and magnet length, L, in the OPAL input file."
                    << endl;
                return;
            } else {
                fieldmap_m->setFieldGap(gap_m);
                fieldmap_m->setFieldLength(length_m);
                fieldmap_m->setEdgeConstants(0.0, 0.0, 0.0);
                fieldmap_m->adjustFringeFields();
                msg << "Adjusted fringe field parameters." << endl;
                fieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);
                fieldmap_m->getInfo(&msg);

                // Primitive check to make sure field map is ok.
                if(zEnd < zBegin)
                    return;
            }
        }

        length_m = zEnd - zBegin;

        /*
         * If the bend angle is specified, find the proper field amplitude.
         * Otherwise, calculate the bend angle for the given field strength.
         */
        if(angle_m != 0.0) {
            // Negative bend angle is just a positive bend rotated 180 degrees.
            if(angle_m < 0.0) {
                angle_m *= -1.0;
                Orientation_m(2) += Physics::pi;
            }
            setBendStrength();
            reinitialize_m = true;
        } else if(amplitude_m != 0.0) {
            angle_m = calculateBendAngle(length_m, true);

            // Make negative bend angle a positive bend rotated 180 degrees.
            if(angle_m < 0.0) {
                amplitude_m *= -1.0;
                Orientation_m(2) += Physics::pi;
                angle_m = calculateBendAngle(length_m, true);
            }
            reinitialize_m = false;
        }

        // Calculate the reference particle trajectory map. Make sure field map edge constants
        // are reset.
        double bendAngle = calculateRefTrajectory(zBegin);

        // Set effective length and effective center.
        effectiveLength_m = R_m * std::abs(bendAngle);
        calculateEffectiveCenter();
        effectiveStart_m = startField_m + effectiveCenter_m - effectiveLength_m / 2.0;

        // Pass start and end of field to calling function.
        startField = startField_m;
        endField = endField_m;

        // Output bend information.
        msg << "Start:                " << startField_m << " m (in floor coordinates)" << endl;
        msg << "End:                  " << endField_m << " m (in floor coordinates)" << endl;
        msg << "Bend angle magnitude: " << bendAngle * 180.0 / Physics::pi << " degrees" << endl;
        msg << "Field amplitude:      " << amplitude_m << " T" << endl;
        msg << "Bend radius:          " << R_m << " m" << endl;
        msg << "Effective length:     " << effectiveLength_m << " m (in s coordinates)" << endl;
        msg << "Effective center:     " << effectiveCenter_m << " m (in s coordinates with respect to bend field map start position)" << endl;
        msg << "Effective start:      " << effectiveStart_m << " m (in floor coordinates)" << endl;

    } else {
        endField = startField - 1e-3;
    }
}

void SBend::finalise() {
    online_m = false;
}

bool SBend::bends() const
{ return true; }

std::string SBend::getFieldMapFN() const {
    return filename_m;
}

const std::string &SBend::getType() const {
    static const std::string type("SBend");
    return type;
}

void SBend::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = endField_m;
}

double SBend::getBendAngle() const {
    return angle_m;
}

double SBend::getEffectiveLength() const {
    return effectiveLength_m;
}

double SBend::getEffectiveCenter() const {
    return effectiveCenter_m;
}

double SBend::getR() const {
    return R_m;
}

double SBend::getStartElement() const {
    return startElement_m;
}

void SBend::setAmplitudem(double vPeak) {
    amplitude_m = vPeak;
}

void SBend::setBendAngle(const double &angle) {
    angle_m = angle * Physics::pi / 180.0;
}

void SBend::setAlpha(const double &alpha) {
    alpha_m = alpha * Physics::pi / 180.0;

    /*
     * We set the orientation of the element to negative
     * of this angle in order to preserve the standard
     * definition of the entrance edge angle of a dipole.
     */
    Orientation_m(0) = -alpha * Physics::pi / 180.0;
    sin_face_alpha_m = sin(Orientation_m(0));
    cos_face_alpha_m = cos(Orientation_m(0));
    tan_face_alpha_m = tan(Orientation_m(0));
}

void SBend::setExitAngle(const double &exitAngle) {
    exitAngle_m = exitAngle * Physics::pi / 180.0;
}

void SBend::setBeta(const double &beta) {
    Orientation_m(1) = beta * Physics::pi / 180.0;
    sin_face_beta_m = sin(Orientation_m(1));
    cos_face_beta_m = cos(Orientation_m(1));
    tan_face_beta_m = tan(Orientation_m(1));
}

void SBend::setDesignEnergy(const double &energy) {
    design_energy_m = energy;
}

void SBend::setK1(const double &k1) {
    gradient_m = k1;
}

void SBend::setFieldMapFN(std::string fmapfn) {
    filename_m = fmapfn;
}

void SBend::setFullGap(const double &gap) {
    gap_m = gap;
}

void SBend::setLength(const double &length) {
    length_m = length;
}

void SBend::setLongitudinalRotation(const double &rotation) {
    Orientation_m(2) = rotation * Physics::pi / 180.0;
}

void SBend::setLongitudinalRotation(const double &k0, const double &k0s) {

    Orientation_m(2) = atan2(k0s, k0);
    amplitude_m = sqrt(pow(k0, 2.0) + pow(k0s, 2.0));

    if(fabs(k0s) > 1.e-8) {
        if(amplitude_m > 1.e-8) {
            field_orientation_m(0) = k0 / amplitude_m;
            field_orientation_m(1) = k0s / amplitude_m;
        }
    } else {
        field_orientation_m(0) = 1.0;
        field_orientation_m(1) = 0.0;
    }
}

double SBend::calculateBendAngle(double bendLength, bool modifyField) {

    // Calculate current bend angle given fixed field parameters.
    const double mass = RefPartBunch_m->getM();
    const double gamma = design_energy_m / mass + 1.0;
    const double betaGamma = sqrt(pow(gamma, 2.0) - 1.0);
    const double deltaT = RefPartBunch_m->getdT();

    // Integrate through field for initial angle.
    Vector_t X(0.0, 0.0, 0.0);
    Vector_t P(-betaGamma * sin_face_alpha_m, 0.0, betaGamma * cos_face_alpha_m);
    Vector_t strength(0.0, 0.0, 0.0);
    Vector_t bField(0.0, 0.0, 0.0);
    Vector_t temp(0.0, 0.0, 0.0);

    while(P(2) > 0.0 && X(2) < bendLength) {

        strength = Vector_t(0.0);
        X /= Vector_t(Physics::c * deltaT);
        pusher_m.push(X, P, deltaT);
        X *= Vector_t(Physics::c * deltaT);

        fieldmap_m->getFieldstrength(X, strength, temp);
        bField(1) = amplitude_m * strength(0);
        temp = Vector_t(0.0);

        X /= Vector_t(Physics::c * deltaT);
        pusher_m.kick(X, P, temp, bField, deltaT);

        pusher_m.push(X, P, deltaT);
        X *= Vector_t(Physics::c * deltaT);

    }

    double angle =  -atan2(P(0), P(2)) - Orientation_m(0);

    // If given permission, modify field map parameters to iteratively converge on a solution.
    double error = 1.0;
    while(error > 1.0e-6 && modifyField) {

        fieldmap_m->setEdgeConstants(angle, alpha_m, exitAngle_m);

        X = Vector_t(0.0);
        P = Vector_t(-betaGamma * sin_face_alpha_m, 0.0, betaGamma * cos_face_alpha_m);
        temp = Vector_t(0.0);

        while(P(2) > 0.0 && X(2) < bendLength) {

            strength = Vector_t(0.0);
            X /= Vector_t(Physics::c * deltaT);
            pusher_m.push(X, P, deltaT);
            X *= Vector_t(Physics::c * deltaT);

            fieldmap_m->getFieldstrength(X, strength, temp);
            bField(1) = amplitude_m * strength(0);
            temp = Vector_t(0.0);

            X /= Vector_t(Physics::c * deltaT);
            pusher_m.kick(X, P, temp, bField, deltaT);

            pusher_m.push(X, P, deltaT);
            X *= Vector_t(Physics::c * deltaT);

        }

        double newAngle =  -atan2(P(0), P(2)) - Orientation_m(0);

        error = fabs(newAngle - angle);
        angle = newAngle;
    }

    return angle;
}

// Transverse distance of particle from reference trajectory.
void SBend::calculateDistFromRef(Vector_t X, double &deltaX, double &angle) {
    double x = X(0) * cos(alpha_m) - X(2) * sin(alpha_m) + R_m;
    double z = X(0) * sin(alpha_m) + X(2) * cos(alpha_m) - effectiveStart_m;

    if(z <= 0.0) {
        deltaX = x - R_m;
        angle = 0.0;
    } else {
        deltaX = sqrt(pow(z, 2.0) + pow(x, 2.0)) - R_m;
        angle = atan2(z, x);
        if(std::abs(angle) > std::abs(angle_m))
            angle = angle_m;
    }
}

void SBend::calculateEffectiveLength() {

    double zBegin = 0.0;
    double zEnd = 0.0;
    double rBegin = 0.0;
    double rEnd = 0.0;
    fieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);

    // Uses Simpson's rule to integrate field. Make step size about 1 mm.

    // This must be odd.
    unsigned int numberOfIntSteps = 2 * static_cast<unsigned int>(floor((zEnd - zBegin) * 1000.0 / 2.0)) + 1;

    double deltaZ = (zEnd - zBegin) / numberOfIntSteps;
    effectiveLength_m = 0.0;

    for(unsigned int integralIndex = 1; integralIndex <= (numberOfIntSteps - 1) / 2; integralIndex++) {

        Vector_t strength(0.0);
        Vector_t info(0.0);
        Vector_t X(0.0);
        X(2) = (2 * integralIndex - 1) * deltaZ;
        fieldmap_m->getFieldstrength(X, strength, info);
        double field1 = strength(0);

        X(2) = 2 * integralIndex * deltaZ;
        fieldmap_m->getFieldstrength(X, strength, info);
        double field2 = strength(0);

        X(2) = (2 * integralIndex + 1) * deltaZ;
        fieldmap_m->getFieldstrength(X, strength, info);
        double field3 = strength(0);

        effectiveLength_m += deltaZ * (field1 + 4.0 * field2 + field3) / 3.0;
    }
}

void SBend::calculateEffectiveCenter() {

    double zBegin = 0.0;
    double zEnd = 0.0;
    double rBegin = 0.0;
    double rEnd = 0.0;
    fieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);

    // Initial guess for effective center.
    double effectiveCenter = fabs(R_m * angle_m / 2.0) - zBegin;

    // Find initial angle.
    double actualBendAngle = calculateBendAngle(effectiveCenter, false);

    // Adjust effective center to get a bend angle 0.5 times the full bend angle.
    int iterations = 1;
    double lengthAdjustment = effectiveCenter / 10.0;

    if(fabs(actualBendAngle) > fabs(angle_m / 2.0))
        lengthAdjustment *= -1.0;

    bool lastGreater = true;
    if(fabs(actualBendAngle) < fabs(angle_m / 2.0))
        lastGreater = false;

    while(fabs(actualBendAngle - angle_m / 2.0) > 1.0e-8 && iterations <= 100) {

        actualBendAngle = calculateBendAngle(effectiveCenter, false);
        iterations++;

        if((!lastGreater && fabs(actualBendAngle) > fabs(angle_m / 2.0)) || (lastGreater && fabs(actualBendAngle) < fabs(angle_m / 2.0)))
            lengthAdjustment /= -10.0;

        if(fabs(actualBendAngle) > fabs(angle_m / 2.0)) lastGreater = true;
        else lastGreater = false;

        effectiveCenter += lengthAdjustment;

    }
    effectiveCenter_m = effectiveCenter - R_m * sin(angle_m / 2.0) + R_m * angle_m / 2.0;
}

void SBend::calculateMapField(Vector_t X, double &bX, double &bY, double &bZ) {

    // Get field from field map.
    Vector_t strength(0.0);
    Vector_t info(0.0);
    fieldmap_m->getFieldstrength(X, strength, info);

    //    double dd, dx, dz, rho;
    //    const double &k34 = info(2);
    //
    //    if(k34 > 0) {
    //        dx = X(0) + R_m * cos_face_alpha_m;
    //        dz = X(2) - R_m * sin_face_alpha_m;
    //        rho = sqrt(dx * dx + dz * dz);
    //        dd = 1.0 - rho / R_m;
    //    } else {
    //        dx = -X(0) + R_m * cos_face_alpha_m;
    //        dz = -X(2) + R_m * sin_face_alpha_m;
    //        rho = sqrt(dx * dx + dz * dz);
    //        dd = -1.0 + rho / R_m;
    //    }
    // Calculations for field index.
    double deltaX = 0.0;
    double angle = 0.0;
    if(gradient_m != 0.0)
        calculateDistFromRef(X, deltaX, angle);

    // Calculate field components based on where we are in the field map.
    if(info(0) > 0.99) {
        // We are in fringe field region.
        bY = amplitude_m * (strength(0) - strength(2) / 2.0 * pow(X(1), 2.0)) * (1.0 - gradient_m * deltaX / R_m);
        double bXTemp = -amplitude_m * (strength(0) - strength(2) / 2.0 * pow(X(1), 2.0)) * gradient_m * X(1) / R_m;
        double bZTemp = amplitude_m * strength(1) * X(1);

        if(info(1) > 0.99) {
            // Exit fringe field.
            bX = -bZTemp * sin(angle_m - alpha_m - exitAngle_m) + bXTemp * cos(angle + alpha_m);
            bZ = bZTemp * cos(angle_m - alpha_m - exitAngle_m) + bXTemp * sin(angle + alpha_m);
        } else {
            // Entrance fringe field.
            bX = bXTemp * cos(angle + alpha_m);
            bZ = -bZTemp + bXTemp * sin(angle + alpha_m);
        }

    } else if(fabs(info(0)) < 0.01) {
        // Central portion of magnet.
        bY = amplitude_m * (1.0 - gradient_m * deltaX / R_m);
        bX = -amplitude_m * gradient_m * X(1) * cos(angle - alpha_m) / R_m;
        bZ = -amplitude_m * gradient_m * X(1) * sin(angle - alpha_m) / R_m;
    }
}

double SBend::calculateRefTrajectory(const double zBegin) {

    // Calculate the reference trajectory map.
    const double mass = RefPartBunch_m->getM();
    const double gamma = design_energy_m / mass + 1.;
    const double betaGamma = sqrt(gamma * gamma - 1.);
    const double dt = RefPartBunch_m->getdT();

    int j = 0;

    Vector_t tmp(0.0);
    Vector_t Bfield(0.0);
    Vector_t strength(0.0);
    Vector_t X(0.0);
    Vector_t P(-betaGamma * sin_face_alpha_m, 0.0, betaGamma * cos_face_alpha_m); // TODO: make it 3D

    bool entryFringePassed = false;
    double pathLengthEntryFringe = 0.0;  // in S coordinates. This value is different from zBegin due to the curvature!

    if(map_m != NULL) delete map_m;

    map_step_size_m = betaGamma / gamma * Physics::c * dt;
    map_size_m = (int)floor(length_m / 2. * Physics::pi / map_step_size_m);
    map_m = new double[3 * (map_size_m + 1)];
    map_m[0] = map_m[1] = map_m[2] = 0.0;

    while(map_m[3 * j + 2] < length_m && j < map_size_m) {
        strength = Vector_t(0.0);
        X /= Vector_t(Physics::c * dt);
        pusher_m.push(X, P, dt);
        X *= Vector_t(Physics::c * dt);

        fieldmap_m->getFieldstrength(X, strength, tmp);
        if(X(2) >= fabs(zBegin) && !entryFringePassed) {
            /*
             * This is the point where we pass ELEMEDGEnot the end of the entry
             * fringe field as the name suggests.
             */
            entryFringePassed = true;
            pathLengthEntryFringe = j * map_step_size_m;
        }
        Bfield(1) = amplitude_m * strength(0);
        tmp = Vector_t(0.0);
        X /= Vector_t(Physics::c * dt);
        pusher_m.kick(X, P, tmp, Bfield, dt);
        pusher_m.push(X, P, dt);
        X *= Vector_t(Physics::c * dt);

        ++ j;
        map_m[3 * j] = X(0);
        map_m[3 * j + 1] = X(1);
        map_m[3 * j + 2] = X(2);
    }

    map_size_m = j;
    double angle = -atan2(P(0), P(2)) - Orientation_m(0);

    // Set the radius in the center of the bend.
    R_m = fabs(betaGamma * mass / (Physics::c * amplitude_m));

    // Adjust start and end of field map.
    startField_m = startElement_m - pathLengthEntryFringe;
    endField_m = startField_m + map_step_size_m * j;

    return angle;
}

bool SBend::reinitialize() {

    if(design_energy_m != RefPartBunch_m->get_meanEnergy() * 1.0e6) {
        design_energy_m = RefPartBunch_m->get_meanEnergy() * 1.0e6;

        setBendStrength();

        double zBegin = 0.0;
        double zEnd = 0.0;
        double rBegin = 0.0;
        double rEnd = 0.0;
        fieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);
        calculateRefTrajectory(zBegin);

        Inform msg("SBend ");
        msg << "Bend design energy changed to: " << design_energy_m * 1.0e-6 << " MeV" << endl;
        msg << "Field amplitude:               " << amplitude_m << " T" << endl;
    }

    return false;

}

void SBend::setBendStrength() {
    // This routine uses an iterative procedure to set the bend strength
    // so that the bend angle is the one we want.

    // Estimate bend field magnitude.
    const double mass = RefPartBunch_m->getM();
    const double gamma = design_energy_m / mass + 1.0;
    const double betaGamma = sqrt(pow(gamma, 2.0) - 1.0);
    const double charge = RefPartBunch_m->getQ();

    fieldmap_m->setEdgeConstants(0.0, 0.0, 0.0);
    calculateEffectiveLength();
    double radius = effectiveLength_m / (2.0 * sin(angle_m / 2.0));

    // We want a positive bend angle to bend in the negative x direction independent of charge.
    amplitude_m = (charge / fabs(charge)) * betaGamma * mass / (Physics::c * radius);

    // Find initial angle.
    double zBegin = 0.0;
    double zEnd = 0.0;
    double rBegin = 0.0;
    double rEnd = 0.0;
    fieldmap_m->setEdgeConstants(angle_m, alpha_m, exitAngle_m);
    fieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);
    double actualBendAngle = calculateBendAngle(zEnd - zBegin, false);

    // Search for angle if initial guess is not good enough.
    double error = fabs(actualBendAngle - angle_m);
    if(error > 1.0e-6) {

        double amplitude1 = amplitude_m;
        double bendAngle1 = actualBendAngle;

        // Estimate field adjustment step.
        double fieldStep = (angle_m - bendAngle1) * betaGamma * mass / (2.0 * effectiveLength_m * Physics::c);
        if(pow(amplitude_m * effectiveLength_m * Physics::c / (betaGamma * mass), 2.0) < 1.0)
            fieldStep = (angle_m - bendAngle1) * betaGamma * mass / (2.0 * effectiveLength_m * Physics::c)
                        * std::sqrt(1.0 - pow(amplitude_m * effectiveLength_m * Physics::c / (betaGamma * mass), 2.0));

        fieldStep *= amplitude1 / std::abs(amplitude1);
        double amplitude2 = amplitude_m + fieldStep;
        amplitude_m = amplitude2;
        double bendAngle2 = calculateBendAngle(zEnd - zBegin, false);

        if(fabs(bendAngle1) > fabs(angle_m)) {
            while(fabs(bendAngle2) > fabs(angle_m)) {
                amplitude2 += fieldStep;
                amplitude_m = amplitude2;
                bendAngle2 = calculateBendAngle(zEnd - zBegin, false);
            }
        } else {
            while(fabs(bendAngle2) < fabs(angle_m)) {
                amplitude2 += fieldStep;
                amplitude_m = amplitude2;
                bendAngle2 = calculateBendAngle(zEnd - zBegin, false);
            }
        }

        // Now we should have the proper field amplitude bracketed.
        unsigned int iterations = 1;
        while(error > 1.0e-6 && iterations < 100) {

            amplitude_m = (amplitude1 + amplitude2) / 2.0;
            double newBendAngle = calculateBendAngle(zEnd - zBegin, false);

            error = fabs(newBendAngle - angle_m);

            if(error > 1.0e-6) {

                if(bendAngle1 - angle_m < 0.0) {

                    if(newBendAngle - angle_m < 0.0) {
                        bendAngle1 = newBendAngle;
                        amplitude1 = amplitude_m;
                    } else {
                        bendAngle2 = newBendAngle;
                        amplitude2 = amplitude_m;
                    }

                } else {

                    if(newBendAngle - angle_m < 0.0) {
                        bendAngle2 = newBendAngle;
                        amplitude2 = amplitude_m;
                    } else {
                        bendAngle1 = newBendAngle;
                        amplitude1 = amplitude_m;
                    }
                }
            }
            iterations++;
        }
    }
}
