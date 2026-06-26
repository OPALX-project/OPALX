#ifndef OPAL_OpalLaser_HH
#define OPAL_OpalLaser_HH

#include "Elements/OpalElement.h"

class OpalLaser : public OpalElement {
public:
    enum { WAVELENGTH = COMMON, PULSEENERGY, PULSELENGTH, WAISTX, WAISTY, DIR, STOKES, SIZE };

    OpalLaser();
    ~OpalLaser() override;

    OpalLaser* clone(const std::string& name) override;
    void update() override;

private:
    OpalLaser(const OpalLaser&);
    void operator=(const OpalLaser&);
    OpalLaser(const std::string& name, OpalLaser* parent);
};

#endif  // OPAL_OpalLaser_HH
