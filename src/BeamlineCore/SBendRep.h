#ifndef OPALX_SBendRep_HH
#define OPALX_SBendRep_HH

#include "AbsBeamline/SBend.h"

/**
 * @class SBendRep
 * @brief Concrete OPALX representation of an analytic sector bend.
 */
class SBendRep : public SBend {
public:
    SBendRep();
    explicit SBendRep(const std::string& name);
    SBendRep(const SBendRep&);
    ~SBendRep() override;

    ElementBase* clone() const override;
};

#endif  // OPALX_SBendRep_HH
