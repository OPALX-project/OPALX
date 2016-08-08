#ifndef CLASSIC_SimpleBeamline_HH
#define CLASSIC_SimpleBeamline_HH

// ------------------------------------------------------------------------
// $RCSfile: SimpleBeamline.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Type definition: SimpleBeamline
//
// ------------------------------------------------------------------------
// Class category: Beamlines
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:35 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Beamlines/ElmPtr.h"
#include "Beamlines/TBeamline.h"


// Typedef SimpleBeamline
// ------------------------------------------------------------------------
/// A beam line representation.
//  An example of a beamline sequence constructed from ElmPtr
//  objects.  Such a sequence could be used to build a simple design model
//  of a beamline.

typedef TBeamline<ElmPtr> SimpleBeamline;

#endif // CLASSIC_SimpleBeamline_HH
