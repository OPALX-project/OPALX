#ifndef OPAL_OpalConstantEFieldCavity_HH
#define OPAL_OpalConstantEFieldCavity_HH

#include <memory>

#include "Elements/OpalElement.h"

class OpalConstantEFieldCavity : public OpalElement {
public:
    enum { EX = COMMON, EY, EZ, SIZE };

    OpalConstantEFieldCavity();
    virtual ~OpalConstantEFieldCavity();

    virtual OpalConstantEFieldCavity* clone(const std::string& name);

    /// Build the runtime element of one occurrence from its attribute values.
    static std::shared_ptr<ElementBase> build(const ElementValues& values);

private:
    /// Set the fields of @p element from @p values.
    static void fill(ElementBase& element, const ElementValues& values);

    OpalConstantEFieldCavity(const OpalConstantEFieldCavity&);
    void operator=(const OpalConstantEFieldCavity&);

    OpalConstantEFieldCavity(const std::string& name, OpalConstantEFieldCavity* parent);
};

#endif  // OPAL_OpalConstantEFieldCavity_HH
