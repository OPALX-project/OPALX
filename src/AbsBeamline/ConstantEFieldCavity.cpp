#include "AbsBeamline/ConstantEFieldCavity.h"
#include "PartBunch/PartBunch.h"

ConstantEFieldCavity::ConstantEFieldCavity() : ConstantEFieldCavity("") {}

ConstantEFieldCavity::ConstantEFieldCavity(const ConstantEFieldCavity& right)
    : ElementBase(right), Ex_m(right.Ex_m), Ey_m(right.Ey_m), Ez_m(right.Ez_m) {}

ConstantEFieldCavity::ConstantEFieldCavity(const std::string& name)
    : ElementBase(name), Ex_m(0.0), Ey_m(0.0), Ez_m(0.0) {}

ConstantEFieldCavity::~ConstantEFieldCavity() {}

void ConstantEFieldCavity::initialise(PartBunch_t* bunch) { RefPartBunch_m = bunch; }

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

namespace {
    /// The field of the cavity, for the kernel and the host: (Ex, Ey, Ez) inside the field
    /// window and the aperture, nothing elsewhere.
    KOKKOS_INLINE_FUNCTION void addConstantField(
            const ElementRegions& regions, const Vector_t<double, 3>& R, const double Ex,
            const double Ey, const double Ez, Vector_t<double, 3>& E) {
        if (!regions.fieldAppliesAt(R)) return;
        E(0) += Ex;
        E(1) += Ey;
        E(2) += Ez;
    }
}  // namespace

void ConstantEFieldCavity::apply(
        const std::shared_ptr<ParticleContainer_t>& pc, const StepTime& /*time*/) {
    auto Rview          = pc->R.getView();
    auto Eview          = pc->E.getView();
    const size_t nLocal = pc->getLocalNum();

    const ElementRegions regions = getRegions();
    const double Ex              = Ex_m;
    const double Ey              = Ey_m;
    const double Ez              = Ez_m;

    Kokkos::parallel_for(
            "ConstantEFieldCavity::apply()", nLocal, KOKKOS_LAMBDA(const size_t i) {
                addConstantField(regions, Rview(i), Ex, Ey, Ez, Eview(i));
            });
}

void ConstantEFieldCavity::applyOnHost(
        const Vector_t<double, 3>& R, double /*t*/, Vector_t<double, 3>& E,
        Vector_t<double, 3>& /*B*/) {
    addConstantField(getRegions(), R, Ex_m, Ey_m, Ez_m, E);
}
