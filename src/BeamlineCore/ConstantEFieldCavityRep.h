#ifndef OPALX_ConstantEFieldCavityRep_HH
#define OPALX_ConstantEFieldCavityRep_HH

#include "AbsBeamline/ConstantEFieldCavity.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/ConstEzField.h"

class ConstantEFieldCavityRep : public ConstantEFieldCavity {
public:
    explicit ConstantEFieldCavityRep(const std::string& name);
    ConstantEFieldCavityRep();
    ConstantEFieldCavityRep(const ConstantEFieldCavityRep&);
    virtual ~ConstantEFieldCavityRep();

    ElementBase* clone() const override;

    Channel* getChannel(const std::string& aKey, bool = false) override;

    ConstEzField& getField() override;
    const ConstEzField& getField() const override;

    StraightGeometry& getGeometry() override;
    const StraightGeometry& getGeometry() const override;

    void setElementLength(double length) override;
    void setEz(double ez) override;

private:
    void operator=(const ConstantEFieldCavityRep&);

    StraightGeometry geometry;
    ConstEzField field;
};

#endif  // OPALX_ConstantEFieldCavityRep_HH
