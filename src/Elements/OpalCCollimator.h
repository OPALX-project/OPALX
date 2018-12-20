#ifndef OPAL_OpalCCollimator_HH
#define OPAL_OpalCCollimator_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalSlit.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalCCollimator
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann, Jianjun Yang $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"

class ParticleMatterInteraction;
// Class OpalCCollimator
// ------------------------------------------------------------------------
/// The CCOLLIMATOR element.

class OpalCCollimator: public OpalElement {

public:

    /// The attributes of class OpalCCollimator.
    enum {
        XSTART = COMMON,  // Start of x coordinate
        XEND,             // End of x coordinate
        YSTART,           // Start of y coordinate
        YEND,             // End of y coordinate
        ZSTART,           // Top boundary
        ZEND,             // Bottom boundary
        WIDTH,            // The width of collimator
        OUTFN,            // Output filename
        SIZE
    };

    /// Exemplar constructor.
    OpalCCollimator();

    virtual ~OpalCCollimator();

    /// Make clone.
    virtual OpalCCollimator *clone(const std::string &name);

    /// Fill in all registered attributes.
    virtual void fillRegisteredAttributes(const ElementBase &, ValueFlag);

    /// Update the embedded CLASSIC collimator.
    virtual void update();

private:

    // Not implemented.
    OpalCCollimator(const OpalCCollimator &);
    void operator=(const OpalCCollimator &);

    // Clone constructor.
    OpalCCollimator(const std::string &name, OpalCCollimator *parent);
    ParticleMatterInteraction *parmatint_m;
};

#endif // OPAL_OpalCCollimator_HH
