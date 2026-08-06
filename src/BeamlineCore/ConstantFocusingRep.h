#ifndef OPALX_ConstantFocusingRep_HH
#define OPALX_ConstantFocusingRep_HH

#include "AbsBeamline/ConstantFocusing.h"
#include "BeamlineGeometry/Geometry.h"

class ConstantFocusingRep : public ConstantFocusing {
public:
    explicit ConstantFocusingRep(const std::string& name);
    ConstantFocusingRep();
    ConstantFocusingRep(const ConstantFocusingRep& right);
    ~ConstantFocusingRep() override;

    ElementBase* clone() const override;
    Channel* getChannel(const std::string& key, bool create = false) override;

    Geometry& getGeometry() override;
    const Geometry& getGeometry() const override;

private:
    void operator=(const ConstantFocusingRep&);

    Geometry geometry;
};

#endif  // OPALX_ConstantFocusingRep_HH
