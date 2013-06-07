#ifndef CLASSIC_XMonitorRep_HH
#define CLASSIC_XMonitorRep_HH

// ------------------------------------------------------------------------
// $RCSfile: XMonitorRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: XMonitorRep
//
// ------------------------------------------------------------------------
// Class category: BeamlineCore
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:34 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "BeamlineCore/MonitorRep.h"


// Class XMonitorRep
// ------------------------------------------------------------------------
/// Representation for a orbit position monitor.
//  Acts on the horizontal plane.

class XMonitorRep: public MonitorRep {

public:

    /// Constructor with given name.
    explicit XMonitorRep(const string &name);

    XMonitorRep();
    XMonitorRep(const XMonitorRep &);
    virtual ~XMonitorRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns NULL.
    virtual Channel *getChannel(const string &aKey, bool = false);

    /// Construct an image.
    //  Return the image of the element, containing the name and type string
    //  of the element, and a copy of the user-defined attributes.
    virtual ElementImage *getImage() const;

    /// Get plane.
    //  Return the x-plane for this class.
    virtual Plane getPlane() const;

private:

    // Not implemented.
    void operator=(const XMonitorRep &);
};

#endif // CLASSIC_XMonitorRep_HH
