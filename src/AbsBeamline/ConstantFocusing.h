#ifndef OPALX_ConstantFocusing_HH
#define OPALX_ConstantFocusing_HH

#include "AbsBeamline/ElementBase.h"

/**
 * @brief Linear three-dimensional focusing for disorder-induced-heating studies.
 *
 * The first bunch application fixes the gradient to
 * $g=\mathrm{STRENGTH}\,\lVert\langle|\mathbf{E}_{sc}|\rangle\rVert_2/\mathrm{RADIUS}$.
 * The scale and moving centroid are reduced over the distributed bunch, while the resulting field
 * $g(\mathbf{R}-\langle\mathbf{R}\rangle)$ is applied only to particles with element-local
 * $z\in[0,L)$.
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
