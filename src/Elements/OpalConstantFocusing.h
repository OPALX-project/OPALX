#ifndef OPAL_OpalConstantFocusing_HH
#define OPAL_OpalConstantFocusing_HH

#include "Elements/OpalElement.h"

class OpalConstantFocusing : public OpalElement {
public:
    enum { STRENGTH = COMMON, RADIUS, SIZE };

    OpalConstantFocusing();
    ~OpalConstantFocusing() override;

    OpalConstantFocusing* clone(const std::string& name) override;
    void update() override;

private:
    OpalConstantFocusing(const OpalConstantFocusing&);
    void operator=(const OpalConstantFocusing&);

    OpalConstantFocusing(const std::string& name, OpalConstantFocusing* parent);
};

#endif  // OPAL_OpalConstantFocusing_HH
