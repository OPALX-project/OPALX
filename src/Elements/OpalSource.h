#ifndef OPAL_OPALSOURCE_HH
#define OPAL_OPALSOURCE_HH

#include "Elements/OpalElement.h"


// Class OpalSource
// ------------------------------------------------------------------------
/// The SOURCE element.

class OpalSource: public OpalElement {

public:

    /// The attributes of class OpalSource.
    enum {
        DISTRIBUTION = COMMON,  // The longitudinal magnetic field.
        SIZE
    };

    /// Exemplar constructor.
    OpalSource();

    virtual ~OpalSource();

    /// Make clone.
    virtual OpalSource *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC solenoid.
    virtual void update();

private:

    // Not implemented.
    OpalSource(const OpalSource &);
    void operator=(const OpalSource &);

    // Clone constructor.
    OpalSource(const std::string &name, OpalSource *parent);
};

#endif // OPAL_OPALSOURCE_HH
