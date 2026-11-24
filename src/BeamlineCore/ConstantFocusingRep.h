#ifndef OPALX_ConstantFocusingRep_HH
#define OPALX_ConstantFocusingRep_HH

#include "AbsBeamline/ConstantFocusing.h"

class ConstantFocusingRep : public ConstantFocusing {
public:
    explicit ConstantFocusingRep(const std::string& name);
    ConstantFocusingRep();
    ConstantFocusingRep(const ConstantFocusingRep& right);
    ~ConstantFocusingRep() override;

    ElementBase* clone() const override;

private:
    void operator=(const ConstantFocusingRep&);
};

#endif  // OPALX_ConstantFocusingRep_HH
