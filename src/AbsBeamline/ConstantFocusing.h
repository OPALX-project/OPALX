#ifndef OPALX_ConstantFocusing_HH
#define OPALX_ConstantFocusing_HH

#include "AbsBeamline/ElementBase.h"

/**
 * Linear three-dimensional focusing used by disorder-induced-heating studies.
 *
 * On its first particle application, the element measures the component-wise mean absolute
 * space-charge field and fixes the gradient to
 *
 *   strength * |<|E_sc|>| / radius.
 *
 * The resulting field is applied about the current bunch centroid. This is the moving-reference
 * equivalent of the origin-centred force used by the IPPL P3MHeating example.
 */
class ConstantFocusing : public ElementBase {
public:
    explicit ConstantFocusing(const std::string& name);
    ConstantFocusing();
    ConstantFocusing(const ConstantFocusing& right);
    ~ConstantFocusing() override;

    void accept(BeamlineVisitor& visitor) const override;

    void initialise(PartBunch_t* bunch) override;
    void finalise() override;
    ElementType getType() const override;
    void getFieldExtent(double& zBegin, double& zEnd) const override;

    void apply(const std::shared_ptr<ParticleContainer_t>& pc) override;
    void apply(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;
    bool applyToReferenceParticle(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    double getStrength() const;
    double getRadius() const;
    double getGradient() const;

    void setStrength(double strength);
    void setRadius(double radius);

protected:
    double strength_m;
    double radius_m;
    double gradient_m;
    bool gradientInitialized_m;

private:
    void operator=(const ConstantFocusing&);
};

#endif  // OPALX_ConstantFocusing_HH
