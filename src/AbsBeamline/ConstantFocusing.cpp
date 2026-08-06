#include "AbsBeamline/ConstantFocusing.h"

#include "AbsBeamline/BeamlineVisitor.h"
#include "PartBunch/PartBunch.h"
#include "Utilities/OpalException.h"

#include <cmath>
#include <functional>

ConstantFocusing::ConstantFocusing() : ConstantFocusing("") {}

ConstantFocusing::ConstantFocusing(const ConstantFocusing& right)
    : ElementBase(right),
      strength_m(right.strength_m),
      radius_m(right.radius_m),
      gradient_m(right.gradient_m),
      gradientInitialized_m(right.gradientInitialized_m) {}

ConstantFocusing::ConstantFocusing(const std::string& name)
    : ElementBase(name),
      strength_m(1.0),
      radius_m(1.0),
      gradient_m(0.0),
      gradientInitialized_m(false) {}

ConstantFocusing::~ConstantFocusing() {}

void ConstantFocusing::accept(BeamlineVisitor& visitor) const {
    visitor.visitConstantFocusing(*this);
}

void ConstantFocusing::initialise(PartBunch_t* bunch) {
    if (radius_m <= 0.0) {
        throw OpalException(
                "ConstantFocusing::initialise", "CONSTANTFOCUSING requires RADIUS > 0.");
    }
    RefPartBunch_m        = bunch;
    gradient_m            = 0.0;
    gradientInitialized_m = false;
}

void ConstantFocusing::finalise() {}

ElementType ConstantFocusing::getType() const { return ElementType::CONSTANTFOCUSING; }

void ConstantFocusing::getFieldExtent(double& zBegin, double& zEnd) const {
    zBegin = 0.0;
    zEnd   = getGeometry().getElementLength();
}

double ConstantFocusing::getStrength() const { return strength_m; }

double ConstantFocusing::getRadius() const { return radius_m; }

double ConstantFocusing::getGradient() const { return gradient_m; }

void ConstantFocusing::setStrength(double strength) {
    strength_m            = strength;
    gradientInitialized_m = false;
}

void ConstantFocusing::setRadius(double radius) {
    radius_m              = radius;
    gradientInitialized_m = false;
}

void ConstantFocusing::apply(const std::shared_ptr<ParticleContainer_t>& pc) {
    const size_t nTotal = pc->getTotalNum();
    const size_t nLocal = pc->getLocalNum();
    if (nTotal == 0 || strength_m == 0.0) {
        return;
    }

    const auto R = pc->R.getView();
    const auto E = pc->E.getView();

    if (!gradientInitialized_m) {
        Vector_t<double, 3> avgAbsE(0.0);
        Kokkos::parallel_reduce(
                "ConstantFocusing::initialField", nLocal,
                KOKKOS_LAMBDA(const size_t i, Vector_t<double, 3>& sum) {
                    sum[0] += Kokkos::abs(E(i)[0]);
                    sum[1] += Kokkos::abs(E(i)[1]);
                    sum[2] += Kokkos::abs(E(i)[2]);
                },
                avgAbsE);
        Kokkos::fence();
        ippl::Comm->allreduce(avgAbsE[0], 3, std::plus<double>());
        avgAbsE /= static_cast<double>(nTotal);

        gradient_m            = strength_m * std::sqrt(dot(avgAbsE, avgAbsE)) / radius_m;
        gradientInitialized_m = true;
    }

    Vector_t<double, 3> centroid(0.0);
    Kokkos::parallel_reduce(
            "ConstantFocusing::centroid", nLocal,
            KOKKOS_LAMBDA(const size_t i, Vector_t<double, 3>& sum) { sum += R(i); }, centroid);
    Kokkos::fence();
    ippl::Comm->allreduce(centroid[0], 3, std::plus<double>());
    centroid /= static_cast<double>(nTotal);

    const double gradient = gradient_m;
    Kokkos::parallel_for(
            "ConstantFocusing::apply", nLocal,
            KOKKOS_LAMBDA(const size_t i) { E(i) += gradient * (R(i) - centroid); });
}

void ConstantFocusing::apply(
        const Vector_t<double, 3>& /*R*/, const Vector_t<double, 3>& /*P*/, const double& /*t*/,
        Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& /*B*/) {}

bool ConstantFocusing::applyToReferenceParticle(
        const Vector_t<double, 3>& /*R*/, const Vector_t<double, 3>& /*P*/, const double& /*t*/,
        Vector_t<double, 3>& /*E*/, Vector_t<double, 3>& /*B*/) {
    return false;
}
