#ifndef OPALX_MarkerRep_HH
#define OPALX_MarkerRep_HH

// ------------------------------------------------------------------------
// $RCSfile: MarkerRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MarkerRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/Marker.h"

// Class MarkerRep
// ------------------------------------------------------------------------
/// Representation for a marker element.

class MarkerRep : public Marker {
public:
    /// Constructor with given name.
    explicit MarkerRep(const std::string& name);

    MarkerRep();
    MarkerRep(const MarkerRep&);
    virtual ~MarkerRep();

private:
    // Not implemented.
    void operator=(const MarkerRep&);
};

#endif  // OPALX_MarkerRep_HH
