#include "AbsBeamline/ConstantEFieldCavity.h"
#include "AbsBeamline/BeamlineVisitor.h"
#include "PartBunch/PartBunch.h"

ConstantEFieldCavity::ConstantEFieldCavity() : ConstantEFieldCavity("") {}

ConstantEFieldCavity::ConstantEFieldCavity(const ConstantEFieldCavity& right)
    : ElementBase(right), Ex_m(right.Ex_m), Ey_m(right.Ey_m), Ez_m(right.Ez_m) {}

ConstantEFieldCavity::ConstantEFieldCavity(const std::string& name)
    : ElementBase(name), Ex_m(0.0), Ey_m(0.0), Ez_m(0.0) {}

ConstantEFieldCavity::~ConstantEFieldCavity() {}

void ConstantEFieldCavity::accept(BeamlineVisitor& visitor) const {
    visitor.visitConstantEFieldCavity(*this);
}

void ConstantEFieldCavity::initialise(PartBunch_t* bunch) { RefPartBunch_m = bunch; }

void ConstantEFieldCavity::finalise() {}

ElementType ConstantEFieldCavity::getType() const { return ElementType::CONSTANTEFIELDCAVITY; }

void ConstantEFieldCavity::getFieldExtent(double& zBegin, double& zEnd) const {
    // Local-chart field-support interval.
    zBegin = 0.0;
    zEnd   = getGeometry().getElementLength();
}

double ConstantEFieldCavity::getEx() const { return Ex_m; }

double ConstantEFieldCavity::getEy() const { return Ey_m; }

double ConstantEFieldCavity::getEz() const { return Ez_m; }

void ConstantEFieldCavity::setEx(double ex) { Ex_m = ex; }

void ConstantEFieldCavity::setEy(double ey) { Ey_m = ey; }

void ConstantEFieldCavity::setEz(double ez) { Ez_m = ez; }

bool ConstantEFieldCavity::apply(const std::shared_ptr<ParticleContainer_t>& pc) {
    auto Rview          = pc->R.getView();
    auto Eview          = pc->E.getView();
    const size_t nLocal = pc->getLocalNum();

    const double elemLength = getGeometry().getElementLength();
    const double Ex         = Ex_m;
    const double Ey         = Ey_m;
    const double Ez         = Ez_m;

    Kokkos::parallel_for(
            "ConstantEFieldCavity::apply()", nLocal, KOKKOS_LAMBDA(const size_t i) {
                if (Rview(i)[2] > 0.0 && Rview(i)[2] <= elemLength) {
                    Eview(i)[0] += Ex;
                    Eview(i)[1] += Ey;
                    Eview(i)[2] += Ez;
                }
            });
    return false;
}

bool ConstantEFieldCavity::apply(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& /*t*/,
        Vector_t<double, 3>& E, Vector_t<double, 3>& /*B*/) {
    if (R(2) < 0.0 || R(2) > getGeometry().getElementLength()) return false;
    if (!isInsideTransverse(R)) return getFlagDeleteOnTransverseExit();

    E(0) += Ex_m;
    E(1) += Ey_m;
    E(2) += Ez_m;
    return false;
}

bool ConstantEFieldCavity::applyToReferenceParticle(
        const Vector_t<double, 3>& R, const Vector_t<double, 3>& /*P*/, const double& /*t*/,
        Vector_t<double, 3>& E, Vector_t<double, 3>& /*B*/) {
    if (R(2) < 0.0 || R(2) > getGeometry().getElementLength()) return false;
    if (!isInsideTransverse(R)) return true;

    E(0) += Ex_m;
    E(1) += Ey_m;
    E(2) += Ez_m;
    return false;
}
