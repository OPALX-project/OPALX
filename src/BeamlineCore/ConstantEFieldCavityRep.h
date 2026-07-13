#ifndef OPALX_ConstantEFieldCavityRep_HH
#define OPALX_ConstantEFieldCavityRep_HH

#include "AbsBeamline/ConstantEFieldCavity.h"
#include "BeamlineGeometry/Geometry.h"

class ConstantEFieldCavityRep : public ConstantEFieldCavity {
public:
    explicit ConstantEFieldCavityRep(const std::string& name);
    ConstantEFieldCavityRep();
    ConstantEFieldCavityRep(const ConstantEFieldCavityRep&);
    virtual ~ConstantEFieldCavityRep();

    ElementBase* clone() const override;

    Channel* getChannel(const std::string& aKey, bool = false) override;

    Geometry& getGeometry() override;
    const Geometry& getGeometry() const override;

private:
    void operator=(const ConstantEFieldCavityRep&);

    Geometry geometry;
};

#endif  // OPALX_ConstantEFieldCavityRep_HH
