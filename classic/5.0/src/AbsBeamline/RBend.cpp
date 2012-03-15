// ------------------------------------------------------------------------
// $RCSfile: RBend.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: RBend
//   Defines the abstract interface for a rectangular bend magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/RBend.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.hh"
#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class RBend
// ------------------------------------------------------------------------

int RBend::RBend_counter_m = 0;

RBend::RBend():
    Component(),
    filename_m(""),
    fieldmap_m(NULL),
    amplitude_m(0.0),
    field_orientation_m(1.0, 0.0),
    startField_m(0.0),
    endField_m(0.0),
    length_m(0.0),
    fast_m(false),
    sin_face_alpha_m(0.0),
    cos_face_alpha_m(1.0),
    tan_face_alpha_m(0.0),
    sin_face_beta_m(0.0),
    cos_face_beta_m(1.0),
    tan_face_beta_m(0.0),
    design_energy_m(0.0),
    angle_m(0.0),
    map_m(NULL),
    map_size_m(0),
    map_step_size_m(0.0),
    pusher_m(),
    startElement_m(0.0),
    R_m(0.0),
    effectiveLength_m(0.0),
    effectiveCenter_m(0.0) {
    setElType(isDipole);
}


RBend::RBend(const RBend &right):
    Component(right),
    filename_m(right.filename_m),
    fieldmap_m(right.fieldmap_m),
    amplitude_m(right.amplitude_m),
    field_orientation_m(right.field_orientation_m),
    startField_m(right.startField_m),
    endField_m(right.endField_m),
    length_m(right.length_m),
    fast_m(right.fast_m),
    sin_face_alpha_m(right.sin_face_alpha_m),
    cos_face_alpha_m(right.cos_face_alpha_m),
    tan_face_alpha_m(right.tan_face_alpha_m),
    sin_face_beta_m(right.sin_face_beta_m),
    cos_face_beta_m(right.cos_face_beta_m),
    tan_face_beta_m(right.tan_face_beta_m),
    design_energy_m(right.design_energy_m),
    angle_m(right.angle_m),
    map_size_m(right.map_size_m),
    map_step_size_m(right.map_step_size_m),
    pusher_m(right.pusher_m),
    startElement_m(right.startElement_m),
    R_m(right.R_m),
    effectiveLength_m(right.effectiveLength_m),
    effectiveCenter_m(right.effectiveCenter_m) {
    setElType(isDipole);
    if(map_size_m > 0) {
        map_m = new double[map_size_m + 1];
        for(int i = 0; i < map_size_m + 1; ++i)
            map_m[i] = right.map_m[i];
    } else {
        map_m = NULL;
    }
}


RBend::RBend(const string &name):
    Component(name),
    filename_m(""),
    fieldmap_m(NULL),
    amplitude_m(0.0),
    field_orientation_m(1.0, 0.0),
    startField_m(0.0),
    endField_m(0.0),
    length_m(0.0),
    fast_m(false),
    sin_face_alpha_m(0.0),
    cos_face_alpha_m(1.0),
    tan_face_alpha_m(0.0),
    sin_face_beta_m(0.0),
    cos_face_beta_m(1.0),
    tan_face_beta_m(0.0),
    design_energy_m(0.0),
    angle_m(0.0),
    map_m(NULL),
    map_size_m(0),
    map_step_size_m(0.0),
    pusher_m(),
    startElement_m(0.0),
    R_m(0.0),
    effectiveLength_m(0.0),
    effectiveCenter_m(0.0) {
    setElType(isDipole);
}


RBend::~RBend() {
    if(map_m)
        delete[] map_m;
}


void RBend::accept(BeamlineVisitor &visitor) const {
    visitor.visitRBend(*this);
}


double RBend::getNormalComponent(int n) const {
    return getField().getNormalComponent(n);
}


double RBend::getSkewComponent(int n) const {
    return getField().getSkewComponent(n);
}


void RBend::setNormalComponent(int n, double v) {
    getField().setNormalComponent(n, v);
}


void RBend::setSkewComponent(int n, double v) {
    getField().setSkewComponent(n, v);
}

//Implemented BET functions for dipole
//ff
//transverse kick
void RBend::addKR(int i, double t, Vector_t &K) {
    Inform msg("RBend::addK()");

    Vector_t tmpE(0.0, 0.0, 0.0);
    Vector_t tmpB(0.0, 0.0, 0.0);
    Vector_t tmpE_diff(0.0, 0.0, 0.0);
    Vector_t tmpB_diff(0.0, 0.0, 0.0);
    double pz = RefPartBunch_m->getZ(i) - startField_m - ds_m;
    const Vector_t tmpA(RefPartBunch_m->getX(i) - dx_m, RefPartBunch_m->getY(i) - dy_m, pz);

    DiffDirection zdir(DZ);
    myFieldmap->getFieldstrength(tmpA, tmpE, tmpB);
    myFieldmap->getFieldstrength_fdiff(tmpA, tmpE_diff, tmpB_diff, zdir);

    double g = RefPartBunch_m->getGamma(i);

    if(fabs(scale_m * tmpB_diff(2)) > 0.1) {
        double cf = Physics::q_e * tmpB(2) / (g * Physics::EMASS);
        K += Vector_t(-pow(cf * scale_m * tmpB(0), 2) / 3.0, -pow(cf * scale_m * tmpB(1), 2) / 3.0, 0.0);
    }
}

void RBend::addKT(int i, double t, Vector_t &K) {
    Inform msg("RBend::addK()");

    Vector_t tmpE(0.0, 0.0, 0.0);
    Vector_t tmpB(0.0, 0.0, 0.0);
    double pz = RefPartBunch_m->getZ(i) - startField_m - ds_m;
    const Vector_t tmpA(RefPartBunch_m->getX(i) - dx_m, RefPartBunch_m->getY(i) - dy_m, pz);
    myFieldmap->getFieldstrength(tmpA, tmpE, tmpB);

    double b = RefPartBunch_m->getBeta(i);
    double g = 1 / sqrt(1 - b * b);

    double cf = -Physics::q_e * Physics::c * b * tmpB(2) * scale_m / (g * Physics::EMASS);
    Vector_t temp(cf * tmpB(1), cf * tmpB(0), 0.0);

    //FIXME: K += ??
}

bool RBend::apply(const int &i, const double &t, double E[], double B[]) {
    Vector_t Ev(0, 0, 0), Bv(0, 0, 0);
    if(apply(i, t, Ev, Bv)) return true;

    E[0] = Ev(0);
    E[1] = Ev(1);
    E[2] = Ev(2);
    B[0] = Bv(0);
    B[1] = Bv(1);
    B[2] = Bv(2);

    return false;
}

bool RBend::apply(const int &i, const double &t, Vector_t &E, Vector_t &B) {

    const Vector_t &X = RefPartBunch_m->X[i];
    Vector_t strength(0.0), info(0.0);

    fieldmap_m->getFieldstrength(X, strength, info);

    if(info(0) > 0.99) {

        B(1) += amplitude_m * (strength(0) - strength(2) / 2. * X(1) * X(1));
        if(info(1) > 0.99) {
            B(2) += amplitude_m * strength(1) * X(1);
        } else {
            B(2) -= amplitude_m * strength(1) * X(1);
        }
    } else if(fabs(info(0)) < 0.01)
        B(1) += amplitude_m;

    return false;
}

bool RBend::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {

    int index = (int)floor((R(2) - startField_m) / map_step_size_m);
    if(index > 0 && index + 1 < map_size_m) {

        // Find indices for z position in pre-computed central trajectory map.
        double lever = (R(2) - startField_m) / map_step_size_m - index;
        double z = (1. - lever) * map_m[index] + lever * map_m[index + 1];

        // Rotate x and y to the the bend's local coordinate system.
        //
        // 1) Rotate about the z axis by angle negative Orientation_m(2).
        // 2) Rotate about the y axis by angle negative Orientation_m(0).
        // 3) Rotate about the x axis by angle Orientation_m(1).

        const double sina = sin(Orientation_m(0));
        const double cosa = cos(Orientation_m(0));
        const double sinb = sin(Orientation_m(1));
        const double cosb = cos(Orientation_m(1));
        const double sinc = sin(Orientation_m(2));
        const double cosc = cos(Orientation_m(2));

        Vector_t X(0.0);
        X(0) = (cosa * cosc) * R(0) + (cosa * sinc) * R(1) - sina *        R(2);
        X(1) = (-cosb * sinc - sina * sinb * cosc) * R(0) + (cosb * cosc - sina * sinb * sinc) * R(1) - cosa * sinb * R(2);
        X(2) = z;

        Vector_t strength(0.0), info(0.0);

        fieldmap_m->getFieldstrength(X, strength, info);
        Vector_t tempB(0.0);

        if(info(0) > 0.99) {

            tempB(1) = amplitude_m * (strength(0) - strength(2) / 2. * X(1) * X(1));
            if(info(1) > 0.99) {
                tempB(2) = amplitude_m * strength(1) * X(1);
            } else {
                tempB(2) = -amplitude_m * strength(1) * X(1);
            }
        } else if(fabs(info(0)) < 0.01)
            tempB(1) = amplitude_m;

        // Rotate field out of the bend's local coordinate system and back to lab frame.
        //
        // 1) Rotate about the x axis by angle Orientation_m(1).
        // 2) Rotate about the y axis by angle Orientation_m(0).
        // 3) Rotate about the z axis by angle negative Orientation_(2).

        B(0) +=  cosa * cosc * tempB(0) + (-sina * sinb * cosc - cosb * sinc) * tempB(1) + (sina * cosb * cosc - sinb * sinc) * tempB(2);
        B(1) +=  cosa * sinc * tempB(0) + (-sina * sinb * sinc + cosb * cosc) * tempB(1) + (sina * cosb * sinc + sinb * cosc) * tempB(2);
        B(2) += -sina *        tempB(0) + (-cosa * sinb) * tempB(1) + (cosa * cosb) * tempB(2);
    }
    return false;
}

void RBend::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    using Physics::c;

    Inform msg("RBend ");

    double zBegin = 0.0;
    double zEnd = 0.0;
    double rBegin = 0.0;
    double rEnd = 0.0;
    startElement_m = startField;

    RefPartBunch_m = bunch;
    pusher_m.initialise(bunch->getReference());

    fieldmap_m = Fieldmap::getFieldmap(filename_m, fast_m);
    fieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);

    if((fieldmap_m != NULL) && (zEnd > zBegin)) {
        const double mass = bunch->getM();
        const double gamma = design_energy_m / mass + 1.;
        const double betagamma = sqrt(gamma * gamma - 1.);
        const double dt = bunch->getdT();
        int j = 0;
        Vector_t tmp(0.0), Bfield(0.0), strength(0.0);
        Vector_t X;
        Vector_t P(-betagamma * sin_face_alpha_m, 0.0, betagamma * cos_face_alpha_m); // TODO: make it 3D
        bool EntryFringe_passed = false;
        double PathLengthEntryFringe = 0.0;  // in S coordinates. This value is different from zBegin due to the curvature!

        msg << getName() << " using file ";
        fieldmap_m->getInfo(&msg);
        Fieldmap::readMap(filename_m);

        // Find effective length and effective center.
        calculateEffectiveLength();
        calculateEffectiveCenter();

        X = Vector_t(0.0, 0.0, 0.0);

        length_m = zEnd - zBegin;

        map_step_size_m = betagamma / gamma * c * dt;

        map_size_m = (int)floor(length_m / 2. * Physics::pi / map_step_size_m);
        map_m = new double[map_size_m + 1];
        map_m[0] = 0.0;

        while(map_m[j] < length_m && j < map_size_m) {

            strength = Vector_t(0.0);
            X /= Vector_t(c * dt);
            pusher_m.push(X, P, dt);
            X *= Vector_t(c * dt);
            fieldmap_m->getFieldstrength(X, strength, tmp);
            if(X(2) >= fabs(zBegin) && !EntryFringe_passed) {
                EntryFringe_passed = true;                       // this is the point where we enter we pass ELEMEDGE
                PathLengthEntryFringe = j * map_step_size_m;     // not the end of the entry fringe field as the name
                // suggests
            }
            Bfield(1) = amplitude_m * strength(0);
            tmp = Vector_t(0.0);
            X /= Vector_t(c * dt);
            pusher_m.kick(X, P, tmp, Bfield, dt);
            pusher_m.push(X, P, dt);
            X *= Vector_t(c * dt);

            map_m[++j] = X(2);

            angle_m = atan2(P(0), P(2));
        }
        map_size_m = j;
        angle_m += Orientation_m(0);

        startField -= PathLengthEntryFringe;
        endField = startField + map_step_size_m * j;

        startField_m = startField;
        endField_m = endField;

        msg << "Start:      " << startField_m << " End: " << endField_m << endl;
        msg << "Bend angle: " << angle_m * 180.0 / Physics::pi << " degrees" << endl;

        R_m = fabs(betagamma * mass / (c * amplitude_m));

    } else {
        endField = startField - 1e-3;
    }

}

void RBend::finalise() {
    online_m = false;
}

bool RBend::bends() const {
    return true;
}

void RBend::setAmplitudem(double vPeak) {
    amplitude_m = vPeak;
}

void RBend::setFieldMapFN(string fmapfn) {
    filename_m = fmapfn;
}

string RBend::getFieldMapFN() const {
    return filename_m;
}

void RBend::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = endField_m;
}

double RBend::getEffectiveLength() const {
    return effectiveLength_m;
}

double RBend::getEffectiveCenter() const {
    return effectiveCenter_m;
}

double RBend::getBendAngle() const {
    return angle_m;
}

double RBend::getStartElement() const {
    return startElement_m;
}

double RBend::getR() const {
    return R_m;
}


const string &RBend::getType() const {
    static const string type("RBend");
    return type;
}

//void RBend::calculateEffectiveLength() {
//
//    double zBegin = 0.0;
//    double zEnd = 0.0;
//    double rBegin = 0.0;
//    double rEnd = 0.0;
//    fieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);
//
//    // Uses Simpson's rule to integrate field. Make step size about 1 mm.
//
//    // This must be odd.
//    unsigned int numberOfIntSteps = 2 * static_cast<unsigned int>(floor((zEnd - zBegin) * 1000.0 / 2.0)) + 1;
//
//    double deltaZ = (zBegin - zEnd) / numberOfIntSteps;
//    effectiveLength_m = 0.0;
//    for(unsigned int integralIndex = 1; integralIndex <= (numberOfIntSteps - 1) / 2; integralIndex++) {
//
//        Vector_t strength(0.0);
//        Vector_t info(0.0);
//        Vector_t X(0.0);
//        X(2) = (2 * integralIndex - 1) * deltaZ;
//        fieldmap_m->getFieldstrength(X, strength, info);
//        double field1 = strength(0);
//
//        X(2) = 2 * integralIndex * deltaZ;
//        fieldmap_m->getFieldstrength(X, strength, info);
//        double field2 = strength(0);
//
//        X(2) = (2 * integralIndex + 1) * deltaZ;
//        fieldmap_m->getFieldstrength(X, strength, info);
//        double field3 = strength(0);
//
//        effectiveLength_m += deltaZ * (field1 + 4.0 * field2 + field3) / 3.0;
//    }
//}

void RBend::calculateEffectiveLength() {

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

void RBend::calculateEffectiveCenter() {

    double zBegin = 0.0;
    double zEnd = 0.0;
    double rBegin = 0.0;
    double rEnd = 0.0;
    fieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);

    // Uses Simpson's rule to integrate field. Make step size 1 mm.

    // This must be odd.
    unsigned int numberOfIntSteps = 2 * static_cast<unsigned int>(floor((zEnd - zBegin) * 1000.0 / 2.0)) + 1;

    double deltaZ = 0.001;
    double integral = 0.0;

    double zBracket1 = 0.0;
    double integralOld = 0.0;

    unsigned int integralIndex = 1;
    while(effectiveCenter_m == 0.0 && integralIndex <= (numberOfIntSteps - 1) / 2) {

        Vector_t strength(0.0);
        Vector_t info(0.0);
        Vector_t X(0.0);
        X(2) = (2 * integralIndex - 1) * deltaZ;
        zBracket1 = X(2);
        fieldmap_m->getFieldstrength(X, strength, info);
        double field1 = strength(0);

        X(2) = 2 * integralIndex * deltaZ;
        fieldmap_m->getFieldstrength(X, strength, info);
        double field2 = strength(0);

        X(2) = (2 * integralIndex + 1) * deltaZ;
        fieldmap_m->getFieldstrength(X, strength, info);
        double field3 = strength(0);

        integralOld = integral;
        integral += deltaZ * (field1 + 4.0 * field2 + field3) / 3.0;

        if(integral >= effectiveLength_m / 2.0)
            effectiveCenter_m = zBracket1 + deltaZ * 2.0 * (effectiveLength_m / 2.0 - integralOld) / (integral - integralOld);

        integralIndex++;
    }
}

