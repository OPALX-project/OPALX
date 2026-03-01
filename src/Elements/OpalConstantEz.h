#ifndef OPAL_OpalConstantEz_HH
#define OPAL_OpalConstantEz_HH

#include "Elements/OpalElement.h"

class OpalConstantEz : public OpalElement {
public:
    enum {
        EZ = COMMON,
        SIZE
    };

    OpalConstantEz();
    virtual ~OpalConstantEz();

    virtual OpalConstantEz* clone(const std::string& name);

    virtual void update();

private:
    OpalConstantEz(const OpalConstantEz&);
    void operator=(const OpalConstantEz&);

    OpalConstantEz(const std::string& name, OpalConstantEz* parent);
};

#endif  // OPAL_OpalConstantEz_HH
