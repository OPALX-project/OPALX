#ifndef CLASSIC_ConstantEFieldCavityRep_HH
#define CLASSIC_ConstantEFieldCavityRep_HH

#include "AbsBeamline/ConstantEFieldCavity.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/ConstEzField.h"

class ConstantEFieldCavityRep : public ConstantEFieldCavity {
public:
    explicit ConstantEFieldCavityRep(const std::string& name);
    ConstantEFieldCavityRep();
    ConstantEFieldCavityRep(const ConstantEFieldCavityRep&);
    virtual ~ConstantEFieldCavityRep();

    virtual ElementBase* clone() const;

    virtual Channel* getChannel(const std::string& aKey, bool = false);

    virtual ConstEzField& getField();
    virtual const ConstEzField& getField() const;

    virtual StraightGeometry& getGeometry();
    virtual const StraightGeometry& getGeometry() const;

    void setElementLength(double length);
    void setEz(double ez) override;

private:
    void operator=(const ConstantEFieldCavityRep&);

    StraightGeometry geometry;
    ConstEzField field;
};

#endif  // CLASSIC_ConstantEFieldCavityRep_HH

