#include "AbsBeamline/ConstantEFieldCavity.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "PartBunch/PartBunch.h"

ConstantEFieldCavity::ConstantEFieldCavity() : ConstantEFieldCavity("") {}

ConstantEFieldCavity::ConstantEFieldCavity(const ConstantEFieldCavity& right)
    : Component(right),
      Ex_m(right.Ex_m),
      Ey_m(right.Ey_m),
      Ez_m(right.Ez_m),
      startField_m(right.startField_m) {}

ConstantEFieldCavity::ConstantEFieldCavity(const std::string& name)
    : Component(name),
      Ex_m(0.0),
      Ey_m(0.0),
      Ez_m(0.0),
      startField_m(0.0) {}

ConstantEFieldCavity::~ConstantEFieldCavity() {}

void ConstantEFieldCavity::accept(BeamlineVisitor& visitor) const {
    visitor.visitConstantEFieldCavity(*this);
}

void ConstantEFieldCavity::initialise(PartBunch_t* bunch, double& startField, double& endField) {
    RefPartBunch_m = bunch;
    startField_m   = startField;
    endField       = startField + getElementLength();
}

void ConstantEFieldCavity::finalise() {}

bool ConstantEFieldCavity::bends() const {
    return false;
}

ElementType ConstantEFieldCavity::getType() const {
    return ElementType::CONSTANTEFIELDCAVITY;
}

void ConstantEFieldCavity::getDimensions(double& zBegin, double& zEnd) const {
    zBegin = startField_m;
    zEnd   = startField_m + getElementLength();
}

double ConstantEFieldCavity::getEx() const {
    return Ex_m;
}

double ConstantEFieldCavity::getEy() const {
    return Ey_m;
}

double ConstantEFieldCavity::getEz() const {
    return Ez_m;
}

void ConstantEFieldCavity::setEx(double ex) {
    Ex_m = ex;
}

void ConstantEFieldCavity::setEy(double ey) {
    Ey_m = ey;
}

void ConstantEFieldCavity::setEz(double ez) {
    Ez_m = ez;
}

bool ConstantEFieldCavity::apply() {
    std::shared_ptr<ParticleContainer_t> pc = RefPartBunch_m->getParticleContainer();
    auto Rview                              = pc->R.getView();
    auto Eview                              = pc->E.getView();
    const size_t nLocal                     = pc->getLocalNum();

    const double elemLength = getElementLength();
    const double Ex         = Ex_m;
    const double Ey         = Ey_m;
    const double Ez         = Ez_m;

    Kokkos::parallel_for(
        "ConstantEFieldCavity::apply()", Kokkos::RangePolicy<>(0, nLocal),
        KOKKOS_LAMBDA(const size_t i) {
            if (Rview(i)[2] > 0.0 && Rview(i)[2] <= elemLength) {
                Eview(i)[0] += Ex;
                Eview(i)[1] += Ey;
                Eview(i)[2] += Ez;
            }
        });
    return false;
}

bool ConstantEFieldCavity::apply(const size_t& i, const double& /*t*/,
                                 Vector_t<double, 3>& E, Vector_t<double, 3>& /*B*/) {
    std::shared_ptr<ParticleContainer_t> pc = RefPartBunch_m->getParticleContainer();
    auto Rview                              = pc->R.getView();
    const Vector_t<double, 3> R             = Rview(i);

    if (R(2) < 0.0 || R(2) > getElementLength())
        return false;
    if (!isInsideTransverse(R))
        return getFlagDeleteOnTransverseExit();

    E(0) += Ex_m;
    E(1) += Ey_m;
    E(2) += Ez_m;
    return false;
}

bool ConstantEFieldCavity::apply(const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/,
                                 const double& /*t*/, Vector_t<double, 3>& E,
                                 Vector_t<double, 3>& /*B*/) {
    if (R(2) < 0.0 || R(2) > getElementLength())
        return false;
    if (!isInsideTransverse(R))
        return getFlagDeleteOnTransverseExit();

    E(0) += Ex_m;
    E(1) += Ey_m;
    E(2) += Ez_m;
    return false;
}

bool ConstantEFieldCavity::applyToReferenceParticle(const Vector_t<double, 3>& R,
                                                    const Vector_t<double, 3>& /*P*/,
                                                    const double& /*t*/,
                                                    Vector_t<double, 3>& E,
                                                    Vector_t<double, 3>& /*B*/) {
    if (R(2) < 0.0 || R(2) > getElementLength())
        return false;
    if (!isInsideTransverse(R))
        return true;

    E(0) += Ex_m;
    E(1) += Ey_m;
    E(2) += Ez_m;
    return false;
}

