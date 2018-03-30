#ifndef OPAL_TRIM_COIL_H
#define OPAL_TRIM_COIL_H

#include <vector>
#include <memory>
#include "AbstractObjects/Definition.h"
#include "TrimCoils/TrimCoil.h"

// Class OpalTrimCoil
// ------------------------------------------------------------------------
/// The TRIMCOIL definition.
//  A TRIMCOIL definition is used to define a trim coil which can be applied
//  to a Cyclotron

class OpalTrimCoil: public Definition {

public:

    /// Exemplar constructor.
    OpalTrimCoil();

    virtual ~OpalTrimCoil();

    /// Test if replacement is allowed.
    //  Can replace only by another OpalTrimCoil
    virtual bool canReplaceBy(Object *object);

    /// Make clone.
    virtual OpalTrimCoil *clone(const std::string &name);

    /// Check the OpalTrimCoil data.
    virtual void execute();

    /// Find named FILTER.
    static OpalTrimCoil *find(const std::string &name);

    /// Update the OpalTrimCoil data.
    virtual void update();

    Inform& print(Inform& os) const;

    void initOpalTrimCoil();

    TrimCoil *trimcoil_m;
private:

    // Not implemented.
    OpalTrimCoil(const OpalTrimCoil &);
    void operator=(const OpalTrimCoil &);

    // Clone constructor.
    OpalTrimCoil(const std::string &name, OpalTrimCoil *parent);

};

inline Inform &operator<<(Inform &os, const OpalTrimCoil &b) {
    return b.print(os);
}

#endif // OPAL_TRIM_COIL_H
