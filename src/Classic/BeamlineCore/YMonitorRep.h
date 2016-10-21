#ifndef CLASSIC_YMonitorRep_HH
#define CLASSIC_YMonitorRep_HH

// ------------------------------------------------------------------------
// $RCSfile: YMonitorRep.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: YMonitorRep
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


// Class YMonitorRep
// ------------------------------------------------------------------------
/// Representation for a orbit position monitor.
//  Acts on the vertical plane.

class YMonitorRep: public MonitorRep {

public:

    /// Constructor with given name.
    explicit YMonitorRep(const std::string &name);

    YMonitorRep();
    YMonitorRep(const YMonitorRep &);
    virtual ~YMonitorRep();

    /// Return clone.
    //  Return an identical deep copy of the element.
    virtual ElementBase *clone() const;

    /// Construct a read/write channel.
    //  This method constructs a Channel permitting read/write access to
    //  the attribute [b]aKey[/b] and returns it.
    //  If the attribute does not exist, it returns NULL.
    virtual Channel *getChannel(const std::string &aKey, bool = false);

    /// Construct an image.
    //  Return the image of the element, containing the name and type string
    //  of the element, and a copy of the user-defined attributes.
    virtual ElementImage *getImage() const;

    /// Get plane.
    //  Return y-plane for this class.
    virtual Plane getPlane() const;

private:

    // Not implemented.
    void operator=(const YMonitorRep &);
};

#endif // CLASSIC_YMonitorRep_HH
