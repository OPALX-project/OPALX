#ifndef OPAL_OpalConstantFocusing_HH
#define OPAL_OpalConstantFocusing_HH

#include <memory>

#include "Elements/OpalElement.h"

class OpalConstantFocusing : public OpalElement {
public:
    enum { STRENGTH = COMMON, RADIUS, SIZE };

    OpalConstantFocusing();
    ~OpalConstantFocusing() override;

    OpalConstantFocusing* clone(const std::string& name) override;
    void update() override;

    /// Build the runtime element of one occurrence from its attribute values.
    static std::shared_ptr<ElementBase> build(const ElementValues& values);

private:
    /// Set the fields of @p element from @p values.
    static void fill(ElementBase& element, const ElementValues& values);

    OpalConstantFocusing(const OpalConstantFocusing&);
    void operator=(const OpalConstantFocusing&);

    OpalConstantFocusing(const std::string& name, OpalConstantFocusing* parent);
};

#endif  // OPAL_OpalConstantFocusing_HH
