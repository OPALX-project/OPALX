#ifndef OPAL_OpalDrift_HH
#define OPAL_OpalDrift_HH

// ------------------------------------------------------------------------
// $RCSfile: OpalDrift.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalDrift
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:39 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Elements/OpalElement.h"
#include "Structure/BoundaryGeometry.h"

// Class OpalDrift
// ------------------------------------------------------------------------
/// The DRIFT element.

class OpalWake;
class SurfacePhysics;

class OpalDrift: public OpalElement {

public:

    enum {
         GEOMETRY = COMMON,       // geometry of boundary, one more enum member besides the common ones in OpalElement.
	 SIZE

    };
    /// Exemplar constructor.
    OpalDrift();

    virtual ~OpalDrift();

    /// Make clone.
    virtual OpalDrift *clone(const string &name);

    /// Test for drift.
    //  Return true.
    virtual bool isDrift() const;

    /// Update the embedded CLASSIC drift.
    virtual void update();

private:

    // Not implemented.
    OpalDrift(const OpalDrift &);
    void operator=(const OpalDrift &);

    // Clone constructor.
    OpalDrift(const string &name, OpalDrift *parent);

    OpalWake *owk_m;
    SurfacePhysics *sphys_m;
    BoundaryGeometry *obgeo_m;
};

#endif // OPAL_OpalDrift_HH
