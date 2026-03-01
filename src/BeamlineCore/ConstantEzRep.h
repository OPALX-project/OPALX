#ifndef CLASSIC_ConstantEzRep_HH
#define CLASSIC_ConstantEzRep_HH

#include "AbsBeamline/ConstantEz.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/ConstEzField.h"

class ConstantEzRep : public ConstantEz {
public:
    explicit ConstantEzRep(const std::string& name);
    ConstantEzRep();
    ConstantEzRep(const ConstantEzRep&);
    virtual ~ConstantEzRep();

    virtual ElementBase* clone() const;

    virtual Channel* getChannel(const std::string& aKey, bool = false);

    virtual ConstEzField& getField();
    virtual const ConstEzField& getField() const;

    virtual StraightGeometry& getGeometry();
    virtual const StraightGeometry& getGeometry() const;

    void setElementLength(double length);
    void setEz(double ez) override;

private:
    void operator=(const ConstantEzRep&);

    StraightGeometry geometry;
    ConstEzField field;
};

#endif  // CLASSIC_ConstantEzRep_HH
