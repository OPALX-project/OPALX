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
    NormalComponents(1, 0.0),
    SkewComponents(1, 0.0),
    max_SkewComponent_m(1),
    max_NormalComponent_m(1),
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


Multipole::Multipole(const std::string &name):
    Component(name),
    NormalComponents(1, 0.0),
    SkewComponents(1, 0.0),
    max_SkewComponent_m(1),
    max_NormalComponent_m(1),
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
    PAssert(n > 1);

    if(n - 1 >  max_NormalComponent_m) {
        max_NormalComponent_m = n - 1;
        NormalComponents.resize(max_NormalComponent_m, 0.0);
    }
    NormalComponents[n - 2] = v; //starting from the quad (= 2)
}

void Multipole::setSkewComponent(int n, double v) {
    //   getField().setSkewComponent(n, v);
    PAssert(n > 1);

    if(n - 1 > max_SkewComponent_m) {
        max_SkewComponent_m = n - 1;
        SkewComponents.resize(max_SkewComponent_m, 0.0);
    }
    SkewComponents[n - 2] = v;  //starting from the quad (= 2)
}

double Multipole::EngeFunc(double z) {
  const double a1 = 0.296417;
  const double a2 = 4.533;
  const double a3 = -2.27;
  const double a4 = 1.06;
  const double a5 = -0.03;
  const double a6 = 0.02;					\
  const double DD = 0.99;

  const double y = z - 1.0;
  return 1.0/(1 + std::exp(a1 + a2*(y/DD) + a3*std::pow(y/DD,2) + a4*std::pow(y/DD,3) + a5*std::pow(y/DD,4) + a6*std::pow(y/DD,5)));
}

double Multipole::EngeFact(double z) {
  // Normalize
  const double lFringe = std::abs(endField_m-startField_m);
  const double zn = (z - startField_m) / lFringe;

  const double scale = 1.0; // 1.4289638937448055; // to make integrated field strenght like the hard edge approximation

  if(zn > 0.5) {
    return scale*EngeFunc((zn-0.5)*10 - 3.) ;
  }
  else
    return scale*EngeFunc((0.5-zn)*10 - 3.) ;
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

    const Vector_t Rt(RefPartBunch_m->getX(i) - dx_m, RefPartBunch_m->getY(i) - dy_m , RefPartBunch_m->getZ(i) - ds_m);
    // before misalignment    Vector_t Rt(RefPartBunch_m->getX(i), RefPartBunch_m->getY(i), RefPartBunch_m->getZ(i));

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

    const Vector_t temp(RefPartBunch_m->getX(i) - dx_m, RefPartBunch_m->getY(i) - dy_m , RefPartBunch_m->getZ(i) - ds_m);
    // before misalignment Vector_t temp(RefPartBunch_m->getX(i), RefPartBunch_m->getY(i), RefPartBunch_m->getZ(i));

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

      FieldFactor(0) = 1.0; // FixMe: need to define a way to incooperate fringe fields EngeFact(R(2));

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

bool Multipole::apply(const Vector_t &R0, const Vector_t &centroid, const double &t, Vector_t &E, Vector_t &B) {

   const Vector_t R = R0 - Vector_t(dx_m,dy_m,ds_m);

   const Vector_t tmpR(R(0) - dx_m, R(1) - dy_m, R(2) - ds_m);

    if(tmpR(2) > startField_m && tmpR(2) <= endField_m) {
        if(max_NormalComponent_m > 0) {
            B(0) += NormalComponents[0] * tmpR(1);
            B(1) += NormalComponents[0] * tmpR(0);

            if(max_NormalComponent_m > 1) {
                const double R02 = tmpR(0) * tmpR(0);
                const double R12 = tmpR(1) * tmpR(1);
                B(0) += NormalComponents[1] * tmpR(0) * tmpR(1);
                B(1) += NormalComponents[1] * (R02 - R12) / 2.;

                if(max_NormalComponent_m > 2) {
                    B(0) += NormalComponents[2] * (3. * R02 * tmpR(1) - R12 * tmpR(1)) / 6.;
                    B(1) += NormalComponents[2] * (R02 * tmpR(0) - 3. * tmpR(0) * R12) / 6.;

                    if(max_NormalComponent_m > 3) {
                        B(0) += NormalComponents[3] * (R02 * tmpR(0) * tmpR(1) - tmpR(0) * R12) / 6.;
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

            B(0) += -SkewComponents[0] * tmpR(0);
            B(1) += SkewComponents[0] * tmpR(1);

            if(max_SkewComponent_m > 1) {
                const double R02 = tmpR(0) * tmpR(0);
                const double R12 = tmpR(1) * tmpR(1);

                B(0) += -SkewComponents[1] * (R02 - R12) / 2.;
                B(1) += SkewComponents[1] * tmpR(0) * tmpR(1);

                if(max_SkewComponent_m > 2) {
                    B(0) += -SkewComponents[2] * (R02 * tmpR(0) - 3. * tmpR(0) * R12) / 6.;
                    B(1) += -SkewComponents[2] * (3. * R02 * tmpR(1) - R12 * tmpR(1)) / 6.;

                    if(max_SkewComponent_m > 3) {
                        B(0) += -SkewComponents[3] * (R02 * R02 - 6. * R02 * R12 + R12 * R12) / 24.;
                        B(1) += SkewComponents[3] * (R02 * tmpR(0) * tmpR(1) - tmpR(0) * R12 * tmpR(1)) / 6.;

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


ElementBase::ElementType Multipole::getType() const {
    return MULTIPOLE;
}
