#ifndef OPAL_OpalConstantEFieldCavity_HH
#define OPAL_OpalConstantEFieldCavity_HH

#include "Elements/OpalElement.h"

class OpalConstantEFieldCavity : public OpalElement {
public:
    enum { EX = COMMON, EY, EZ, SIZE };

    OpalConstantEFieldCavity();
    virtual ~OpalConstantEFieldCavity();

    virtual OpalConstantEFieldCavity* clone(const std::string& name);

    virtual void update();

private:
    OpalConstantEFieldCavity(const OpalConstantEFieldCavity&);
    void operator=(const OpalConstantEFieldCavity&);

    OpalConstantEFieldCavity(const std::string& name, OpalConstantEFieldCavity* parent);
};

#endif  // OPAL_OpalConstantEFieldCavity_HH
