#ifndef OPALX_Marker_HH
#define OPALX_Marker_HH

// ------------------------------------------------------------------------
// $RCSfile: Marker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Marker
//   Defines a marker element.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/ElementBase.h"

// Class Marker
// ------------------------------------------------------------------------
/// Interface for a marker.
//  Class Marker defines a marker element.

class Marker : public FramelessElement {
public:
    /// Constructor with given name.
    explicit Marker(const std::string& name);

    Marker();
    Marker(const Marker&);
    virtual ~Marker();

    virtual void initialise(PartBunch_t* bunch) override;

    virtual ElementType getType() const override;
    TrackedAs trackedAs() const override { return TrackedAs::Ignored; }

    virtual void getFieldExtent(double& zBegin, double& zEnd) const override;

    virtual int getRequiredNumberOfTimeSteps() const override;

private:
    // Not implemented.
    void operator=(const Marker&);
};

inline int Marker::getRequiredNumberOfTimeSteps() const { return 1; }

#endif  // OPALX_Marker_HH