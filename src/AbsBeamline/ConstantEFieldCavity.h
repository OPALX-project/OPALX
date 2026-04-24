#ifndef OPALX_ConstantEFieldCavity_HH
#define OPALX_ConstantEFieldCavity_HH

#include "AbsBeamline/Component.h"

/**
 * @class ConstantEFieldCavity
 * @brief Component applying a constant accelerating electric field (Ex,Ey,Ez).
 *
 * Homogeneous electrostatic field E = (Ex, Ey, Ez) over the element length.
 * GPU-compatible apply() via Kokkos parallel_for.
 */
class ConstantEFieldCavity : public Component {
public:
    explicit ConstantEFieldCavity(const std::string& name);
    ConstantEFieldCavity();
    ConstantEFieldCavity(const ConstantEFieldCavity& right);
    virtual ~ConstantEFieldCavity();

    virtual void accept(BeamlineVisitor& visitor) const override;

    virtual void initialise(PartBunch_t* bunch, double& startField, double& endField) override;
    virtual void finalise() override;
    virtual bool bends() const override;
    virtual ElementType getType() const override;
    virtual void getDimensions(double& zBegin, double& zEnd) const override;

    virtual bool apply(const std::shared_ptr<ParticleContainer_t>& pc) override;
    virtual bool apply(const size_t& i, const double& t,
                       Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;
    virtual bool apply(const Vector_t<double, 3>& R, const Vector_t<double, 3>& P,
                       const double& t, Vector_t<double, 3>& E,
                       Vector_t<double, 3>& B) override;
    virtual bool applyToReferenceParticle(const Vector_t<double, 3>& R,
                                          const Vector_t<double, 3>& P,
                                          const double& t,
                                          Vector_t<double, 3>& E,
                                          Vector_t<double, 3>& B) override;

    double getEx() const;
    double getEy() const;
    double getEz() const;

    virtual void setEx(double ex);
    virtual void setEy(double ey);
    virtual void setEz(double ez);

protected:
    double Ex_m;
    double Ey_m;
    double Ez_m;
    double startField_m;

private:
    void operator=(const ConstantEFieldCavity&);
};

#endif  // OPALX_ConstantEFieldCavity_HH

