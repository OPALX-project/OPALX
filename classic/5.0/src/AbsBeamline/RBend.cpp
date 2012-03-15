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
    map_m(NULL),
    map_size_m(0),
    map_step_size_m(0.0),
    pusher_m(),
    startElement_m(0.0),
    R_m(0.0)
{
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
    map_size_m(right.map_size_m),
    map_step_size_m(right.map_step_size_m),
    pusher_m(right.pusher_m),
    startElement_m(right.startElement_m),
    R_m(right.R_m)
{
    setElType(isDipole);
    if (map_size_m > 0) {
        map_m = new double[map_size_m + 1];
        for (int i = 0; i < map_size_m + 1; ++i)
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
    map_m(NULL),
    map_size_m(0),
    map_step_size_m(0.0),
    pusher_m(),
    startElement_m(0.0),
    R_m(0.0)
{
    setElType(isDipole);
}


RBend::~RBend()
{
    if (map_m)
        delete[] map_m;
}


void RBend::accept(BeamlineVisitor &visitor) const
{
    visitor.visitRBend(*this);
}


double RBend::getNormalComponent(int n) const
{
    return getField().getNormalComponent(n);
}


double RBend::getSkewComponent(int n) const
{
    return getField().getSkewComponent(n);
}


void RBend::setNormalComponent(int n, double v)
{
    getField().setNormalComponent(n, v);
}


void RBend::setSkewComponent(int n, double v)
{
    getField().setSkewComponent(n, v);
}

//Implemented BET functions for dipole
//ff
//transverse kick
void RBend::addKR(int i, double t, Vector_t &K) 
{
    Inform msg("RBend::addK()");

    Vector_t tmpE(0.0,0.0,0.0); 
    Vector_t tmpB(0.0,0.0,0.0);
    Vector_t tmpE_diff(0.0,0.0,0.0); 
    Vector_t tmpB_diff(0.0,0.0,0.0);
    double pz = RefPartBunch_m->getZ(i) - startField_m - ds_m;
    const Vector_t tmpA(RefPartBunch_m->getX(i) - dx_m, RefPartBunch_m->getY(i) - dy_m, pz);

    bool out_of_bounds;
    DiffDirection zdir(DZ);
    out_of_bounds = myFieldmap->getFieldstrength(tmpA,tmpE,tmpB);
    out_of_bounds = myFieldmap->getFieldstrength_fdiff(tmpA,tmpE_diff,tmpB_diff,zdir);

    double g = RefPartBunch_m->getGamma(i);

    if(fabs(scale_m*tmpB_diff(2)) > 0.1) {
        double dbdz = scale_m * tmpB_diff(2);
        double cf = Physics::q_e*tmpB(2)/(g*Physics::m_bet);
        K += Vector_t(-pow(cf*scale_m*tmpB(0),2)/3.0,-pow(cf*scale_m*tmpB(1),2)/3.0,0.0); 
    }
}

void RBend::addKT(int i, double t, Vector_t &K) 
{
    Inform msg("RBend::addK()");

    Vector_t tmpE(0.0,0.0,0.0); 
    Vector_t tmpB(0.0,0.0,0.0);
    double pz = RefPartBunch_m->getZ(i) - startField_m - ds_m;
    const Vector_t tmpA(RefPartBunch_m->getX(i) - dx_m, RefPartBunch_m->getY(i) - dy_m, pz);
    bool out_of_bounds = myFieldmap->getFieldstrength(tmpA,tmpE,tmpB);

    double b = RefPartBunch_m->getBeta(i);
    double g = 1/sqrt(1-b*b);

    double cf = -Physics::q_e*Physics::c*b*tmpB(2)*scale_m/(g*Physics::m_bet); 
    Vector_t temp(cf*tmpB(1),cf*tmpB(0),0.0);

    //FIXME: K += ??
}

bool RBend::apply(const int &i, const double &t, double E[], double B[])
{
    Vector_t Ev(0,0,0), Bv(0,0,0);
    if (apply(i,t,Ev,Bv)) return true;
      
    E[0] = Ev(0); E[1] = Ev(1); E[2] = Ev(2);
    B[0] = Bv(0); B[1] = Bv(1); B[2] = Bv(2);
      
    return false;
}

bool RBend::apply(const int &i, const double &t, Vector_t &E, Vector_t &B)
{
    Vector_t temp1(RefPartBunch_m->getX(i),RefPartBunch_m->getY(i),RefPartBunch_m->getZ(i));
    Vector_t temp2(RefPartBunch_m->getPx(i),RefPartBunch_m->getPy(i),RefPartBunch_m->getPz(i));

    const Vector_t &R(temp1);
    const Vector_t &P(temp2);
    double fraction = 1.0;

    const Vector_t &X = RefPartBunch_m->X[i];
    Vector_t strength(0.0), info(0.0);

    fieldmap_m->getFieldstrength(X, strength, info);

    if (info(0) > 0.99) {

        double y_tilde = field_orientation_m(1) * X(0) + field_orientation_m(0) * X(1);
        B(0) += field_orientation_m(1) * amplitude_m * (strength(0) - strength(2)/2. * y_tilde * y_tilde);
        B(1) += field_orientation_m(0) * amplitude_m * (strength(0) - strength(2)/2. * y_tilde * y_tilde);
        if (info(1) > 0.99) {
            B(2) += amplitude_m * strength(1) * y_tilde;
        } else {
            B(2) -= amplitude_m * strength(1) * y_tilde;
        }
    } else if (fabs(info(0)) < 0.01) {
        B(0) += field_orientation_m(1) * amplitude_m;
        B(1) += field_orientation_m(0) * amplitude_m;
    }

    return false;
}
  
bool RBend::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B)
{
    int index = (int)floor((R(2) - startField_m) / map_step_size_m);
    if (index > 0 && index + 1 < map_size_m) {
        double lever = (R(2) - startField_m) / map_step_size_m - index;
        double z = (1. - lever) * map_m[index] + lever * map_m[index + 1];
        
        Vector_t X = Vector_t(R(0), R(1), z);
        Vector_t strength(0.0), info(0.0);
        fieldmap_m->getFieldstrength(X, strength, info);
        
        if (info(0) > 0.99) {
            double y_tilde = field_orientation_m(1) * X(0) + field_orientation_m(0) * X(1);
            B(0) += field_orientation_m(1) * amplitude_m * (strength(0) - strength(2)/2. * y_tilde * y_tilde);
            B(1) += field_orientation_m(0) * amplitude_m * (strength(0) - strength(2)/2. * y_tilde * y_tilde);
            if (info(1) > 0.99) {
                B(2) +=  amplitude_m * strength(1) * y_tilde;
            } else {
                B(2) -=  amplitude_m * strength(1) * y_tilde;
            }
        } else if (fabs(info(0)) < 0.01) {
            B(0) += field_orientation_m(1) * amplitude_m;
            B(1) += field_orientation_m(0) * amplitude_m;
        }
    }
    return false;
}

void RBend::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor)
{
    using Physics::c;

    Inform msg("RBend ");

    double tmpDouble;
    double alpha, beta;
    double tolerance = 1e-10;
    double zBegin = 0.0;
    double zEnd = 0.0;
    double rBegin = 0.0;
    double rEnd = 0.0;
    startElement_m= startField;

    RefPartBunch_m = bunch;
    pusher_m.initialise(bunch->getReference());

    fieldmap_m = Fieldmap::getFieldmap(filename_m, fast_m);
    fieldmap_m->getFieldDimensions(zBegin, zEnd, rBegin, rEnd);
    if ((fieldmap_m != NULL) && (zEnd > zBegin)) {
        msg << "RBend " << getName() << " using file ";
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

        fieldmap_m->getInfo(&msg);
        Fieldmap::readMap(filename_m);

        X = Vector_t(0.0, 0.0, 0.0);

        length_m = zEnd - zBegin;

        map_step_size_m = betagamma / gamma * c * dt;

        map_size_m = (int)floor(length_m / 2. * Physics::pi / map_step_size_m);
        map_m = new double[map_size_m + 1];
        map_m[0] = 0.0;
        while (map_m[j] < length_m && j < map_size_m) {
            strength = Vector_t(0.0);
            X /= Vector_t(c * dt);
            pusher_m.push(X, P, dt);
            X *= Vector_t(c * dt);
            fieldmap_m->getFieldstrength(X, strength, tmp);
            if (X(2) >= fabs(zBegin) && !EntryFringe_passed) {
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

void RBend::finalise()
{
    online_m = false;
}

bool RBend::bends() const
{
    return true;
}

void RBend::setAmplitudem(double vPeak)
{
    amplitude_m = vPeak;
}

void RBend::setFieldMapFN(string fmapfn)
{
    filename_m = fmapfn;
}

string RBend::getFieldMapFN() const
{
    return filename_m;
}

void RBend::getDimensions(double &zBegin, double &zEnd) const
{
    zBegin = startField_m;
    zEnd = endField_m;
}

double RBend::getStartElement() const
{
    return startElement_m;
}

double RBend::getR() const
{
    return R_m;
}


const string& RBend::getType() const
{
    static const string type("RBend");
    return type;
}

