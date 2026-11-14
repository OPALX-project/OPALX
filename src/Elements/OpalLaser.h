#ifndef OPAL_OpalLaser_HH
#define OPAL_OpalLaser_HH

#include <memory>

#include "Elements/OpalElement.h"

class OpalLaser : public OpalElement {
public:
    enum { WAVELENGTH = COMMON, PULSEENERGY, PULSELENGTH, WAISTX, WAISTY, DIR, STOKES, SIZE };

    OpalLaser();
    ~OpalLaser() override;

    /// Build the runtime element of one occurrence from its attribute values.
    static std::shared_ptr<ElementBase> build(const ElementValues& values);

private:
    /// Set the fields of @p element from @p values.
    static void fill(ElementBase& element, const ElementValues& values);

    OpalLaser(const OpalLaser&);
    void operator=(const OpalLaser&);
};

#endif  // OPAL_OpalLaser_HH
