#ifndef OPALX_RBendRep_HH
#define OPALX_RBendRep_HH

#include "AbsBeamline/RBend.h"

/**
 * @class RBendRep
 * @brief Concrete OPALX representation of an analytic rectangular bend.
 *
 * The representation owns the rectangular bend geometry and the local analytic
 * multipole field expansion used by the first OPALX-native `RBEND` port.
 */
class RBendRep : public RBend {
public:
    RBendRep();
    explicit RBendRep(const std::string& name);
    RBendRep(const RBendRep&);
    ~RBendRep() override;

    ElementBase* clone() const override;
    Channel* getChannel(const std::string& aKey, bool create = false) override;

    BMultipoleField& getField() override;
    const BMultipoleField& getField() const override;

    RBendGeometry& getGeometry() override;
    const RBendGeometry& getGeometry() const override;

    void setField(const BMultipoleField& field);

private:
    RBendGeometry geometry_m;
    BMultipoleField field_m;
};

#endif  // OPALX_RBendRep_HH
