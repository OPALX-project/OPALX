#include "AbsBeamline/ConstantEz.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "PartBunch/PartBunch.h"

ConstantEz::ConstantEz() : ConstantEz("") {}

ConstantEz::ConstantEz(const ConstantEz& right)
    : Component(right),
      Ez_m(right.Ez_m),
      startField_m(right.startField_m) {}

ConstantEz::ConstantEz(const std::string& name)
    : Component(name),
      Ez_m(0.0),
      startField_m(0.0) {}

ConstantEz::~ConstantEz() {}

void ConstantEz::accept(BeamlineVisitor& visitor) const {
    visitor.visitConstantEz(*this);
}

void ConstantEz::initialise(PartBunch_t* bunch, double& startField, double& endField) {
    RefPartBunch_m = bunch;
    startField_m   = startField;
    endField       = startField + getElementLength();
}

void ConstantEz::finalise() {}

bool ConstantEz::bends() const {
    return false;
}

ElementType ConstantEz::getType() const {
    return ElementType::CONSTANTEZ;
}

void ConstantEz::getDimensions(double& zBegin, double& zEnd) const {
    zBegin = startField_m;
    zEnd   = startField_m + getElementLength();
}

double ConstantEz::getEz() const {
    return Ez_m;
}

void ConstantEz::setEz(double ez) {
    Ez_m = ez;
}

bool ConstantEz::apply() {
    std::shared_ptr<ParticleContainer_t> pc = RefPartBunch_m->getParticleContainer();
    auto Rview                               = pc->R.getView();
    auto Eview                               = pc->E.getView();

    const double elemLength = getElementLength();
    const double Ez         = Ez_m;

    Kokkos::parallel_for(
        "ConstantEz::apply()", ippl::getRangePolicy(Rview),
        KOKKOS_LAMBDA(const int i) {
            if (Rview(i)(2) > 0.0 && Rview(i)(2) <= elemLength) {
                Eview(i)(2) += Ez;
            }
        });
    return false;
}

bool ConstantEz::apply(const size_t& i, const double& /*t*/,
                       Vector_t<double, 3>& E, Vector_t<double, 3>& /*B*/) {
    std::shared_ptr<ParticleContainer_t> pc = RefPartBunch_m->getParticleContainer();
    auto Rview                               = pc->R.getView();
    const Vector_t<double, 3> R             = Rview(i);

    if (R(2) < 0.0 || R(2) > getElementLength())
        return false;
    if (!isInsideTransverse(R))
        return getFlagDeleteOnTransverseExit();

    E(2) += Ez_m;
    return false;
}

bool ConstantEz::apply(const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/,
                       const double& /*t*/, Vector_t<double, 3>& E,
                       Vector_t<double, 3>& /*B*/) {
    if (R(2) < 0.0 || R(2) > getElementLength())
        return false;
    if (!isInsideTransverse(R))
        return getFlagDeleteOnTransverseExit();

    E(2) += Ez_m;
    return false;
}

bool ConstantEz::applyToReferenceParticle(const Vector_t<double, 3>& R,
                                          const Vector_t<double, 3>& /*P*/,
                                          const double& /*t*/,
                                          Vector_t<double, 3>& E,
                                          Vector_t<double, 3>& /*B*/) {
    if (R(2) < 0.0 || R(2) > getElementLength())
        return false;
    if (!isInsideTransverse(R))
        return true;

    E(2) += Ez_m;
    return false;
}
