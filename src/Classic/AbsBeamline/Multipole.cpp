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
#include "Fields/Fieldmap.h"
#include "Physics/Physics.h"

extern Inform *gmsg;

namespace{
    enum {
        DIPOLE,
        QUADRUPOLE,
        SEXTUPOLE,
        OCTUPOLE,
        DECAPOLE
    };
}
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
    PAssert(n >= 1);

    if(n >  max_NormalComponent_m) {
        max_NormalComponent_m = n;
        NormalComponents.resize(max_NormalComponent_m, 0.0);
    }
    switch(n-1) {
    case SEXTUPOLE:
        NormalComponents[n - 1] = v / 2;
        break;
    case OCTUPOLE:
    case DECAPOLE:
        NormalComponents[n - 1] = v / 24;
        break;
    default:
        NormalComponents[n - 1] = v;
    }
}

void Multipole::setSkewComponent(int n, double v) {
    //   getField().setSkewComponent(n, v);
    PAssert(n >= 1);

    if(n > max_SkewComponent_m) {
        max_SkewComponent_m = n;
        SkewComponents.resize(max_SkewComponent_m, 0.0);
    }
    switch(n-1) {
    case SEXTUPOLE:
        SkewComponents[n - 1] = v / 2;
        break;
    case OCTUPOLE:
        SkewComponents[n - 1] = v / 6;
        break;
    case DECAPOLE:
        SkewComponents[n - 1] = v / 24;
        break;
    default:
        SkewComponents[n - 1] = v;
    }
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

void Multipole::computeField(Vector_t R, const double &t, Vector_t &E, Vector_t &B) {

    if(R(2) > startField_m && R(2) <= endField_m) {

        {
            std::vector<Vector_t> Rn(max_NormalComponent_m + 1);
            std::vector<double> fact(max_NormalComponent_m + 1);
            Rn[0] = Vector_t(1.0);
            fact[0] = 1;
            for (int i = 0; i < max_NormalComponent_m; ++ i) {
                switch(i) {
                case DIPOLE:
                    B(1) += NormalComponents[i];
                    break;

                case QUADRUPOLE:
                    B(0) += NormalComponents[i] * R(1);
                    B(1) += NormalComponents[i] * R(0);
                    break;

                case SEXTUPOLE:
                    B(0) += 2 * NormalComponents[i] * R(0) * R(1);
                    B(1) += NormalComponents[i] * (Rn[2](0) - Rn[2](1));
                    break;

                case OCTUPOLE:
                    B(0) += NormalComponents[i] * (3 * Rn[2](0) * Rn[1](1) - Rn[3](1));
                    B(1) += NormalComponents[i] * (Rn[3](0) - 3 * Rn[1](0) * Rn[2](1));
                    break;

                case DECAPOLE:
                    B(0) += 4 * NormalComponents[i] * (Rn[3](0) * Rn[1](1) - Rn[1](0) * Rn[3](1));
                    B(1) += NormalComponents[i] * (Rn[4](0) - 6 * Rn[2](0) * Rn[2](1) + Rn[4](1));
                    break;

                default:
                    {
                        double powMinusOne = 1;
                        double Bx = 0.0, By = 0.0;
                        for (int j = 1; j <= (i + 1) / 2; ++ j) {
                            Bx += powMinusOne * NormalComponents[i] * (Rn[i - 2 * j + 1](0) * fact[i - 2 * j + 1] *
                                                                       Rn[2 * j - 1](1) * fact[2 * j - 1]);
                            By += powMinusOne * NormalComponents[i] * (Rn[i - 2 * j + 2](0) * fact[i - 2 * j + 2] *
                                                                       Rn[2 * j - 2](1) * fact[2 * j - 2]);
                            powMinusOne *= -1;
                        }

                        if ((i + 1) / 2 == i / 2) {
                            int j = (i + 2) / 2;
                            By += powMinusOne * NormalComponents[i] * (Rn[i - 2 * j + 2](0) * fact[i - 2 * j + 2] *
                                                                       Rn[2 * j - 2](1) * fact[2 * j - 2]);
                        }
                        B(0) += Bx;
                        B(1) += By;
                    }
                }

                Rn[i + 1](0) = Rn[i](0) * R(0);
                Rn[i + 1](1) = Rn[i](1) * R(1);
                fact[i + 1] = fact[i] / (i + 1);
            }
        }

        {

            std::vector<Vector_t> Rn(max_SkewComponent_m + 1);
            std::vector<double> fact(max_SkewComponent_m + 1);
            Rn[0] = Vector_t(1.0);
            fact[0] = 1;
            for (int i = 0; i < max_SkewComponent_m; ++ i) {
                switch(i) {
                case DIPOLE:
                    B(0) -= SkewComponents[i];
                    break;

                case QUADRUPOLE:
                    B(0) -= SkewComponents[i] * R(0);
                    B(1) += SkewComponents[i] * R(1);
                    break;

                case SEXTUPOLE:
                    B(0) -= SkewComponents[i] * (Rn[2](0) - Rn[2](1));
                    B(1) += 2 * SkewComponents[i] * R(0) * R(1);
                    break;

                case OCTUPOLE:
                    B(0) -= SkewComponents[i] * (Rn[3](0) - 3 * Rn[1](0) * Rn[2](1));
                    B(1) += SkewComponents[i] * (3 * Rn[2](0) * Rn[1](1) - Rn[3](1));
                    break;

                case DECAPOLE:
                    B(0) -= SkewComponents[i] * (Rn[4](0) - 6 * Rn[2](0) * Rn[2](1) + Rn[4](1));
                    B(1) += 4 * SkewComponents[i] * (Rn[3](0) * Rn[1](1) - Rn[1](0) * Rn[3](1));
                    break;

                default:
                    {
                        double powMinusOne = 1;
                        double Bx = 0, By = 0;
                        for (int j = 1; j <= (i + 1) / 2; ++ j) {
                            Bx -= powMinusOne * SkewComponents[i] * (Rn[i - 2 * j + 2](0) * fact[i - 2 * j + 2] *
                                                                     Rn[2 * j - 2](1) * fact[2 * j - 2]);
                            By += powMinusOne * SkewComponents[i] * (Rn[i - 2 * j + 1](0) * fact[i - 2 * j + 1] *
                                                                     Rn[2 * j - 1](1) * fact[2 * j - 1]);
                            powMinusOne *= -1;
                        }

                        if ((i + 1) / 2 == i / 2) {
                            int j = (i + 2) / 2;
                            Bx -= powMinusOne * SkewComponents[i] * (Rn[i - 2 * j + 2](0) * fact[i - 2 * j + 2] *
                                                                     Rn[2 * j - 2](1) * fact[2 * j - 2]);
                        }

                        B(0) += Bx;
                        B(1) += By;
                    }
                }

                Rn[i + 1](0) = Rn[i](0) * R(0);
                Rn[i + 1](1) = Rn[i](1) * R(1);
                fact[i + 1] = fact[i] / (i + 1);
            }
        }
    }

}

bool Multipole::apply(const size_t &i, const double &t, double E[], double B[]) {
    Vector_t Ef(0.0), Bf(0.0);
    Vector_t R(RefPartBunch_m->getX(i) - dx_m, RefPartBunch_m->getY(i) - dy_m, RefPartBunch_m->getZ(i) - ds_m);
    computeField(R, t, Ef, Bf);
    for (unsigned int d = 0; d < 3; ++ d) {
        E[d] += Ef(d);
        B[d] += Bf(d);
    }

    return false;
}
bool Multipole::apply(const size_t &i, const double &t, Vector_t &E, Vector_t &B) {
    Vector_t R(RefPartBunch_m->getX(i) - dx_m, RefPartBunch_m->getY(i) - dy_m, RefPartBunch_m->getZ(i) - ds_m);
    computeField(R, t, E, B);

    return false;
}

bool Multipole::apply(const Vector_t &R0, const Vector_t &, const double &t, Vector_t &E, Vector_t &B) {
    Vector_t R = R0 - Vector_t(dx_m, dy_m, ds_m);
    computeField(R, t, E, B);

    return false;
}
void Multipole::initialise(PartBunch *bunch, double &startField, double &endField, const double&) {
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