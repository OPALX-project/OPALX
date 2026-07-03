#ifndef CLASSIC_BeamBeam_HH
#define CLASSIC_BeamBeam_HH

// ------------------------------------------------------------------------
// $RCSfile: Ip.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamBeam
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

// Class BeamBeam
// ------------------------------------------------------------------------
/// Interface for a beam-beam interaction element.
//  Class BeamBeam defines the abstract interface for a beam-beam element.

class BeamBeam : public Component {
public:
    /// Constructor with given name.
    explicit BeamBeam(const std::string& name);

    BeamBeam();
    BeamBeam(const BeamBeam& right);
    virtual ~BeamBeam();

    /// Apply visitor to BeamBeam.
    virtual void accept(BeamlineVisitor&) const override;

    virtual void initialise(PartBunch_t* bunch, double& startField, double& endField) override;

    virtual void finalise() override;

    virtual bool bends() const override;

    virtual ElementType getType() const override;

    virtual void getFieldExtend(double& zBegin, double& zEnd) const override;

    // set number of slices for map tracking
    void setNSlices(const std::size_t& nSlices);  // Philippe was here

    // set number of slices for map tracking
    std::size_t getNSlices() const;  // Philippe was here

    virtual int getRequiredNumberOfTimeSteps() const override;

private:
    double startField_m;
    std::size_t nSlices_m;

    // Not implemented.
    void operator=(const BeamBeam&);
};

inline int BeamBeam::getRequiredNumberOfTimeSteps() const { return 1; }

#endif  // CLASSIC_BeamBeam_HH
