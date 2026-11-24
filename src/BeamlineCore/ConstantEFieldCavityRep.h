#ifndef OPALX_ConstantEFieldCavityRep_HH
#define OPALX_ConstantEFieldCavityRep_HH

#include "AbsBeamline/ConstantEFieldCavity.h"

class ConstantEFieldCavityRep : public ConstantEFieldCavity {
public:
    explicit ConstantEFieldCavityRep(const std::string& name);
    ConstantEFieldCavityRep();
    ConstantEFieldCavityRep(const ConstantEFieldCavityRep&);
    virtual ~ConstantEFieldCavityRep();

    ElementBase* clone() const override;

private:
    void operator=(const ConstantEFieldCavityRep&);
};

#endif  // OPALX_ConstantEFieldCavityRep_HH
