#ifndef CLASSIC_Ip_HH
#define CLASSIC_Ip_HH

// ------------------------------------------------------------------------
// $RCSfile: Ip.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Ip
//   Defines the abstract interface for a drift space.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Component.h"

// Class Ip
// ------------------------------------------------------------------------
/// Interface for drift space.
//  Class Ip defines the abstract interface for a drift space.

class Ip : public Component {
public:
    /// Constructor with given name.
    explicit Ip(const std::string& name);

    Ip();
    Ip(const Ip& right);
    virtual ~Ip();

    /// Apply visitor to Ip.
    virtual void accept(BeamlineVisitor&) const override;

    virtual void initialise(PartBunch_t* bunch, double& startField, double& endField) override;

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual ElementType getType() const override;

    virtual void getDimensions(double& zBegin, double& zEnd) const override;

    // set number of slices for map tracking
    void setNSlices(const std::size_t& nSlices);  // Philippe was here

    // set number of slices for map tracking
    std::size_t getNSlices() const;  // Philippe was here

    virtual int getRequiredNumberOfTimeSteps() const override;

private:
    double startField_m;
    std::size_t nSlices_m;

    // Not implemented.
    void operator=(const Ip&);
};

inline int Ip::getRequiredNumberOfTimeSteps() const {
    return 1;
}

#endif  // CLASSIC_Ip_HH
