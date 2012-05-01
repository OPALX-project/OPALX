// ------------------------------------------------------------------------
// $RCSfile: Multipole.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Multipole
//   Defines the abstract interface for a Multipole magnet.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Multipole.h"
#include "Algorithms/PartBunch.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "Fields/Fieldmap.hh"
#include "Physics/Physics.h"

extern Inform *gmsg;

// Class Multipole
// ------------------------------------------------------------------------

Multipole::Multipole():
    Component(),
    myFieldmap_m(NULL) {
    setElType(isMultipole);
}


Multipole::Multipole(const Multipole &right):
    Component(right),
    NormalComponents(right.NormalComponents),
    SkewComponents(right.SkewComponents),
    max_SkewComponent_m(right.max_SkewComponent_m),
    max_NormalComponent_m(right.max_NormalComponent_m),
    myFieldmap_m(right.myFieldmap_m) {
    setElType(isMultipole);
}


Multipole::Multipole(const string &name):
    Component(name),
    NormalComponents(4, 0.0),
    SkewComponents(4, 0.0),
    max_SkewComponent_m(0),
    max_NormalComponent_m(0),
    myFieldmap_m(NULL) {
    setElType(isMultipole);
}


Multipole::~Multipole()
{}


void Multipole::accept(BeamlineVisitor &visitor) const {
    visitor.visitMultipole(*this);
}


double Multipole::getNormalComponent(int n) const {
    return getField().getNormalComponent(n);
}


double Multipole::getSkewComponent(int n) const {
    return getField().getSkewComponent(n);
}


void Multipole::setNormalComponent(int n, double v) {
    //   getField().setNormalComponent(n, v);
    NormalComponents[n-2] = v; //starting from the quad (= 2)
    if(n - 1 > max_NormalComponent_m)
        max_NormalComponent_m = n - 1;
}


void Multipole::setSkewComponent(int n, double v) {
    //   getField().setSkewComponent(n, v);
    SkewComponents[n-2] = v;  //starting from the quad (= 2)
    if(n - 1 > max_SkewComponent_m)
        max_SkewComponent_m = n - 1;
}

//ff
// radial focussing term
void Multipole::addKR(int i, double t, Vector_t &K) {
    Inform msg("Multipole::addK()");

    double b = RefPartBunch_m->getBeta(i);
    double g = RefPartBunch_m->getGamma(i); //1 / sqrt(1 - b * b);

    // calculate the average of all normal components, to obtain the gradient
    double l = NormalComponents.size();
    double temp_n = 0;
    for(int j = 0; j < l; j++)
        temp_n += NormalComponents.at(j);

    double Grad = temp_n / l;
    double k = -Physics::q_e * b * Physics::c * Grad / (g * Physics::EMASS);
    //FIXME: factor? k *= 5?

    //FIXME: sign?
    Vector_t temp(k, -k, 0.0);

    K += temp;
}

//ff
//transverse kick
void Multipole::addKT(int i, double t, Vector_t &K) {
    Inform msg("Multipole::addK()");

    Vector_t tmpE(0.0, 0.0, 0.0);
    Vector_t tmpB(0.0, 0.0, 0.0);
    Vector_t tmpE_diff(0.0, 0.0, 0.0);
    Vector_t tmpB_diff(0.0, 0.0, 0.0);

    double b = RefPartBunch_m->getBeta(i);
    double g = RefPartBunch_m->getGamma(i);

    // calculate the average of all normal components, to obtain the gradient

    double l = NormalComponents.size();
    double temp_n = 0;
    for(int j = 0; j < l; j++)
        temp_n += NormalComponents.at(j);

    double G = temp_n / l;
    double cf = -Physics::q_e * b * Physics::c * G / (g * Physics::EMASS);
    double dx = RefPartBunch_m->getX0(i) - dx_m;
    double dy = RefPartBunch_m->getY0(i) - dy_m;

    K += Vector_t(cf * dx, -cf * dy, 0.0);
}

bool Multipole::apply(const size_t &i, const double &t, double E[], double B[]) {
    Vector_t Ev(0, 0, 0), Bv(0, 0, 0);
    Vector_t Rt(RefPartBunch_m->getX(i), RefPartBunch_m->getY(i), RefPartBunch_m->getZ(i));
    if(apply(Rt, Vector_t(0.0), t, Ev, Bv)) return true;

    E[0] = Ev(0);
    E[1] = Ev(1);
    E[2] = Ev(2);
    B[0] = Bv(0);
    B[1] = Bv(1);
    B[2] = Bv(2);

    return false;
}

bool Multipole::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    Vector_t temp(RefPartBunch_m->getX(i), RefPartBunch_m->getY(i), RefPartBunch_m->getZ(i));

    const Vector_t &R(temp);
    //const Vector_t &R = RefPartBunch_m->R[i];

    Vector_t FieldFactor(1.0, 0.0, 0.0);
    if(myFieldmap_m) {
        if(myFieldmap_m->getType() == T3DMagnetoStatic) {
            return false;
        } else if(myFieldmap_m->getType() == T1DProfile1 || myFieldmap_m->getType() == T1DProfile2) {
            Vector_t tmpE(0, 0, 0);
            myFieldmap_m->getFieldstrength(R, tmpE, FieldFactor);
        }
    }

    if(R(2) > startField_m && R(2) <= endField_m) {
        if(max_NormalComponent_m > 0) {
            B(0) += NormalComponents[0] * (FieldFactor(0) * R(1) - FieldFactor(2) * R(1) * R(1) * R(1) / 6.);
            B(1) += NormalComponents[0] * (FieldFactor(0) * R(0) - FieldFactor(2) * R(0) * R(0) / 2.);
            B(2) += NormalComponents[0] * FieldFactor(1) * R(0) * R(1);

            if(max_NormalComponent_m > 1) {
                const double R02 = R(0) * R(0);
                const double R12 = R(1) * R(1);
                B(0) += NormalComponents[1] * R(0) * R(1);
                B(1) += NormalComponents[1] * (R02 - R12) / 2.;

                if(max_NormalComponent_m > 2) {
                    B(0) += NormalComponents[2] * (3. * R02 * R(1) - R12 * R(1)) / 6.;
                    B(1) += NormalComponents[2] * (R02 * R(0) - 3. * R(0) * R12) / 6.;

                    if(max_NormalComponent_m > 3) {
                        B(0) += NormalComponents[3] * (R02 * R(0) * R(1) - R(0) * R12) / 6.;
                        B(1) += NormalComponents[3] * (R02 * R02 - 6. * R02 * R12 + R12 * R12) / 24.;

                        if(max_NormalComponent_m > 4) {
                            Inform msg("Multipole ");
                            msg << Fieldmap::typeset_msg("HIGHER MULTIPOLES THAN DECAPOLE NOT IMPLEMENTED!", "warning") << "\n"
                                << endl;
                        }
                    }
                }
            }
        }

        if(max_SkewComponent_m > 0) {
            B(0) += -SkewComponents[0] * R(0);
            B(1) += SkewComponents[0] * R(1);

            if(max_SkewComponent_m > 1) {
                const double R02 = R(0) * R(0);
                const double R12 = R(1) * R(1);

                B(0) += -SkewComponents[1] * (R02 - R12) / 2.;
                B(1) += SkewComponents[1] * R(0) * R(1);

                if(max_SkewComponent_m > 2) {
                    B(0) += -SkewComponents[2] * (R02 * R(0) - 3. * R(0) * R12) / 6.;
                    B(1) += -SkewComponents[2] * (3. * R02 * R(1) - R12 * R(1)) / 6.;

                    if(max_SkewComponent_m > 3) {
                        B(0) += -SkewComponents[3] * (R02 * R02 - 6. * R02 * R12 + R12 * R12) / 24.;
                        B(1) += SkewComponents[3] * (R02 * R(0) * R(1) - R(0) * R12 * R(1)) / 6.;

                        if(max_SkewComponent_m > 4) {
                            Inform msg("Multipole ");
                            msg << Fieldmap::typeset_msg("HIGHER MULTIPOLES THAN DECAPOLE NOT IMPLEMENTED!", "warning") << "\n"
                                << endl;
                        }
                    }
                }
            }
        }
    }

    return false;
}

bool Multipole::apply(const Vector_t &R, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {
    if(R(2) > startField_m && R(2) <= endField_m) {
        if(max_NormalComponent_m > 0) {
            B(0) += NormalComponents[0] * R(1);
            B(1) += NormalComponents[0] * R(0);

            if(max_NormalComponent_m > 1) {
                const double R02 = R(0) * R(0);
                const double R12 = R(1) * R(1);
                B(0) += NormalComponents[1] * R(0) * R(1);
                B(1) += NormalComponents[1] * (R02 - R12) / 2.;

                if(max_NormalComponent_m > 2) {
                    B(0) += NormalComponents[2] * (3. * R02 * R(1) - R12 * R(1)) / 6.;
                    B(1) += NormalComponents[2] * (R02 * R(0) - 3. * R(0) * R12) / 6.;

                    if(max_NormalComponent_m > 3) {
                        B(0) += NormalComponents[3] * (R02 * R(0) * R(1) - R(0) * R12) / 6.;
                        B(1) += NormalComponents[3] * (R02 * R02 - 6. * R02 * R12 + R12 * R12) / 24.;

                        if(max_NormalComponent_m > 4) {
                            Inform msg("Multipole ");
                            msg << Fieldmap::typeset_msg("HIGHER MULTIPOLES THAN DECAPOLE NOT IMPLEMENTED!", "warning") << "\n"
                                << endl;
                        }
                    }
                }
            }
        }

        if(max_SkewComponent_m > 0) {
            B(0) += -SkewComponents[0] * R(0);
            B(1) += SkewComponents[0] * R(1);

            if(max_SkewComponent_m > 1) {
                const double R02 = R(0) * R(0);
                const double R12 = R(1) * R(1);

                B(0) += -SkewComponents[1] * (R02 - R12) / 2.;
                B(1) += SkewComponents[1] * R(0) * R(1);

                if(max_SkewComponent_m > 2) {
                    B(0) += -SkewComponents[2] * (R02 * R(0) - 3. * R(0) * R12) / 6.;
                    B(1) += -SkewComponents[2] * (3. * R02 * R(1) - R12 * R(1)) / 6.;

                    if(max_SkewComponent_m > 3) {
                        B(0) += -SkewComponents[3] * (R02 * R02 - 6. * R02 * R12 + R12 * R12) / 24.;
                        B(1) += SkewComponents[3] * (R02 * R(0) * R(1) - R(0) * R12 * R(1)) / 6.;

                        if(max_SkewComponent_m > 4) {
                            Inform msg("Multipole ");
                            msg << Fieldmap::typeset_msg("HIGHER MULTIPOLES THAN DECAPOLE NOT IMPLEMENTED!", "warning") << "\n"
                                << endl;
                        }
                    }
                }
            }
        }
    }

    return false;
}
void Multipole::initialise(PartBunch *bunch, double &startField, double &endField, const double &scaleFactor) {
    RefPartBunch_m = bunch;
    endField = startField + getElementLength();
    startField_m = startField;
    endField_m = endField;
    online_m = true;
}


void Multipole::finalise() {
    online_m = false;
}

bool Multipole::bends() const {
    return false;
}


void Multipole::getDimensions(double &zBegin, double &zEnd) const {
    zBegin = startField_m;
    zEnd = endField_m;
}


const string &Multipole::getType() const {
    static const string type("Multipole");
    return type;
}

