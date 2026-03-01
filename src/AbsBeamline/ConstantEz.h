#ifndef CLASSIC_ConstantEz_HH
#define CLASSIC_ConstantEz_HH

#include "AbsBeamline/Component.h"

/**
 * @class ConstantEz
 * @brief Component applying a constant accelerating electric field in z.
 *
 * Homogeneous electrostatic field E = (0, 0, Ez) over the element length.
 * GPU-compatible apply() via Kokkos parallel_for.
 */
class ConstantEz : public Component {
public:
    explicit ConstantEz(const std::string& name);
    ConstantEz();
    ConstantEz(const ConstantEz& right);
    virtual ~ConstantEz();

    virtual void accept(BeamlineVisitor& visitor) const override;

    virtual void initialise(PartBunch_t* bunch, double& startField, double& endField) override;
    virtual void finalise() override;
    virtual bool bends() const override;
    virtual ElementType getType() const override;
    virtual void getDimensions(double& zBegin, double& zEnd) const override;

    virtual bool apply() override;
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

    double getEz() const;
    virtual void setEz(double ez);

protected:
    double Ez_m;
    double startField_m;

private:
    void operator=(const ConstantEz&);
};

#endif  // CLASSIC_ConstantEz_HH
