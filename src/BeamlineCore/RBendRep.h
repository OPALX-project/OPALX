#ifndef OPALX_RBendRep_HH
#define OPALX_RBendRep_HH

#include "AbsBeamline/RBend.h"

/**
 * @class RBendRep
 * @brief Concrete OPALX representation of an analytic rectangular bend.
 */
class RBendRep : public RBend {
public:
    RBendRep();
    explicit RBendRep(const std::string& name);
    RBendRep(const RBendRep&);
    ~RBendRep() override;
};

#endif  // OPALX_RBendRep_HH
