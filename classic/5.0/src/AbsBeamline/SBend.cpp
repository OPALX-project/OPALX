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
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.hh"
#include <iostream>
#include <fstream>

extern Inform *gmsg;

// Class SBend
// ------------------------------------------------------------------------

SBend::SBend():
    Component(),
    fast_m(false),
    filename_m(""),
    amplitude_m(0.0),
    field_orientation_m(1.0, 0.0),
    startField_m(0.0),
    length_m(0.0),
    fieldmap_m(NULL),
    sin_face_alpha_m(0.0),
    cos_face_alpha_m(1.0),
    tan_face_alpha_m(0.0),
    sin_face_beta_m(0.0),
    cos_face_beta_m(1.0),
    tan_face_beta_m(0.0),
    map_step_size_m(0.0),
    map_size_m(0),
    map_m(NULL),
    pusher_m(),
    design_energy_m(0.0),
    startElement_m(0.0) {
    setElType(isDipole);
}


SBend::SBend(const SBend &right):
    Component(right),
    fast_m(right.fast_m),
    filename_m(right.filename_m),
    amplitude_m(right.amplitude_m),
    field_orientation_m(right.field_orientation_m),
    startField_m(right.startField_m),
    length_m(right.length_m),
    fieldmap_m(right.fieldmap_m),
    gradient_m(right.gradient_m),
    sin_face_alpha_m(right.sin_face_alpha_m),
    cos_face_alpha_m(right.cos_face_alpha_m),
    tan_face_alpha_m(right.tan_face_alpha_m),
    sin_face_beta_m(right.sin_face_beta_m),
    cos_face_beta_m(right.cos_face_beta_m),
    tan_face_beta_m(right.tan_face_beta_m),
    map_step_size_m(right.map_step_size_m),
    map_size_m(right.map_size_m),
    pusher_m(right.pusher_m),
    design_energy_m(right.design_energy_m),
    startElement_m(right.startElement_m) {
    setElType(isDipole);
    if(map_size_m > 0) {
        if(map_m)
            delete[] map_m;

        map_m = new double[3*(map_size_m + 1)];
        for(int i = 0; i < 3 * (map_size_m + 1); ++i)
            map_m[i] = right.map_m[i];
    } else {
        map_m = NULL;
    }
}


SBend::SBend(const string &name):
    Component(name),
    fast_m(false),
    filename_m(""),
    amplitude_m(0.0),
    field_orientation_m(1.0, 0.0),
    startField_m(0.0),
    length_m(0.0),
    fieldmap_m(NULL),
    gradient_m(0.0),
    sin_face_alpha_m(0.0),
    cos_face_alpha_m(1.0),
    tan_face_alpha_m(0.0),
    sin_face_beta_m(0.0),
    cos_face_beta_m(1.0),
    tan_face_beta_m(0.0),
    map_step_size_m(0.0),
    map_size_m(0),
    map_m(NULL),
    pusher_m(),
    design_energy_m(0.0),
    startElement_m(0.0) {
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

bool SBend::apply(const int &i, const double &t, double E[], double B[]) {
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

bool SBend::apply(const int &i, const double &t, Vector_t &E, Vector_t &B) {
    double dd, dx, dz, rho;
    const Vector_t &X = RefPartBunch_m->X[i];
    Vector_t strength(0.0), info(0.0);
    fieldmap_m->getFieldstrength(X, strength, info);
    const double &k34 = info(2);

    if(k34 > 0) {
        dx = X(0) + R_m * cos_face_alpha_m;
        dz = X(2) - R_m * sin_face_alpha_m;
        rho = sqrt(dx * dx + dz * dz);
        dd = 1.0 - rho / R_m;
    } else {
        dx = -X(0) + R_m * cos_face_alpha_m;
        dz = -X(2) + R_m * sin_face_alpha_m;
        rho = sqrt(dx * dx + dz * dz);
        dd = -1.0 + rho / R_m;
    }

    if(info(0) > 0.99) {
        B(1) +=  amplitude_m * (strength(0) - strength(2) / 2. * X(1) * X(1)) * (1. - gradient_m * dd);
        double Bx = amplitude_m * (strength(0) - strength(2) / 2. * X(1) * X(1)) * gradient_m * X(1) / (rho * R_m);
        double Bz = amplitude_m * strength(1) * X(1);
        if(info(1) > 0.99) {
            B(0) += - Bz * k34 / sqrt(1. + k34 * k34) + Bx * dx;
            B(2) += Bz / sqrt(1. + k34 * k34) + Bx * dz;
        } else {
            B(0) += Bx * dx;
            B(2) += -Bz + Bx * dz;
        }
    } else if(info(0) < 0.01) {
        B(0) += amplitude_m * gradient_m * X(1) * dx / (rho * R_m);
        B(1) += amplitude_m * (1. - gradient_m * dd);
        B(2) += amplitude_m * gradient_m * X(1) * dz / (rho * R_m);
    }
    return false;
}

bool SBend::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    int index = (int)floor((R(2) - startField_m) / map_step_size_m);
    if(index > 0 && index + 1 < map_size_m) {
        double dd, dx, dz, rho;
        double lever = (R(2) - startField_m) / map_step_size_m - index;
        double z = (1. - lever) * map_m[3*index + 2] + lever * map_m[3*index + 5];
        double x = (1. - lever) * map_m[3*index] + lever * map_m[3*index + 3];
        double y = (1. - lever) * map_m[3*index + 1] + lever * map_m[3*index + 4];

        Vector_t X = Vector_t(x + R(0), y + R(1), z);
        Vector_t strength(0.0), info(0.0);
        fieldmap_m->getFieldstrength(X, strength, info);
        const double &k34 = info(2);

        if(k34 > 0) {
            dx = X(0) + R_m * cos_face_alpha_m;
            dz = X(2) - R_m * sin_face_alpha_m;
            rho = sqrt(dx * dx + dz * dz);
            dd = 1.0 - rho / R_m;
        } else {
            dx = -X(0) + R_m * cos_face_alpha_m;
            dz = -X(2) + R_m * sin_face_alpha_m;
            rho = sqrt(dx * dx + dz * dz);
            dd = -1.0 + rho / R_m;
        }

        if(info(0) > 0.99) {
            B(1) +=  amplitude_m * (strength(0) - strength(2) / 2. * X(1) * X(1)) * (1. - gradient_m * dd);
            double Bx = amplitude_m * (strength(0) - strength(2) / 2. * X(1) * X(1)) * gradient_m * X(1) / (rho * R_m);
            double Bz = amplitude_m * strength(1) * X(1);
            if(info(1) > 0.99) {
                B(0) += - Bz * k34 / sqrt(1. + k34 * k34) + Bx * dx;
                B(2) += Bz / sqrt(1. + k34 * k34) + Bx * dz;
            } else {
                B(0) += Bx * dx;
                B(2) += -Bz + Bx * dz;
            }
        } else if(info(0) < 0.01) {
            B(0) += amplitude_m * gradient_m * X(1) * dx / (rho * R_m);
            B(1) += amplitude_m * (1. - gradient_m * dd);
            B(2) += amplitude_m * gradient_m * X(1) * dz / (rho * R_m);
        }
    }


    return false;

}

void SBend::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    using Physics::c;
    Inform msg("SBend ");

    double tmpDouble;
    double alpha, beta;
    double tolerance = 1e-10;
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
        const double charge = bunch->getQ();
        int j = 0;
        Vector_t tmp(0.0), Bfield(0.0), strength(0.0);
        Vector_t X;
        Vector_t P(-betagamma * sin_face_alpha_m, 0.0, betagamma * cos_face_alpha_m); // TODO: make it 3D
        bool EntryFringe_passed = false;
        double PathLengthEntryFringe = 0.0;  // in S coordinates. This value is different from zBegin due to the curvature!


        msg << getName() << " using file "; fieldmap_m->getInfo(&msg);
        Fieldmap::readMap(filename_m);
        fieldmap_m->setExitFaceSlope(exit_face_slope_m);

        X = Vector_t(0.0, 0.0, 0.0);

        length_m =  zEnd - zBegin;
        if(length_m < 0.0) {
            // there is probably something wrong with the fieldmap
            return;
        }

        map_step_size_m = betagamma / gamma * c * dt;

        map_size_m = (int)floor(length_m / 2. * Physics::pi / map_step_size_m);
        map_m = new double[3*(map_size_m + 1)];
        map_m[0] = map_m[1] = map_m[2] = 0.0;

        while(map_m[3*j + 2] < length_m + X(0) * exit_face_slope_m && j < map_size_m) {
            strength = Vector_t(0.0);
            X /= Vector_t(c * dt);
            pusher_m.push(X, P, dt);
            X *= Vector_t(c * dt);

            fieldmap_m->getFieldstrength(X, strength, tmp);
            if(X(2) >= fabs(zBegin) && !EntryFringe_passed) {
                EntryFringe_passed = true;
                PathLengthEntryFringe = j * map_step_size_m;
            }
            Bfield(1) = amplitude_m * strength(0);
            tmp = Vector_t(0.0);
            X /= Vector_t(c * dt);
            pusher_m.kick(X, P, tmp, Bfield, dt);
            pusher_m.push(X, P, dt);
            X *= Vector_t(c * dt);

            ++ j;
            map_m[3*j] = X(0);
            map_m[3*j + 1] = X(1);
            map_m[3*j + 2] = X(2);
            //          *gmsg << "poix \t"<<X(0)<<"\t posz \t" <<X(2)<<"\t poi \t"<<j* map_step_size_m<<"\t field \t"<<Bfield(1)<< endl;
        }
        map_size_m = j;
        startField -= PathLengthEntryFringe;
        endField = startField + map_step_size_m * j;

        startField_m = startField;
        endField_m = endField;
        R_m = fabs(betagamma * mass / (c * amplitude_m));

    } else {
        endField = startField - 1e-3;
    }
}

void SBend::finalise() {
    online_m = false;
}

bool SBend::bends() const
{ return true; }


void SBend::setAmplitudem(double vPeak) {
    amplitude_m = vPeak;
}

void SBend::setFieldMapFN(string fmapfn) {
    filename_m = fmapfn;
}

string SBend::getFieldMapFN() const {
    return filename_m;
}

void SBend::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = endField_m;
}

double SBend::getStartElement() const {
    return startElement_m;
}

double SBend::getR() const {
    return R_m;
}


const string &SBend::getType() const {
    static const string type("SBend");
    return type;
}

