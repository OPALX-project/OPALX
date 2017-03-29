// ------------------------------------------------------------------------
// $RCSfile: OpalParticle.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: OpalParticle
//   A OpalParticle represents the phase space coordinates of a particle.
//   It can be propagated through a beamline.
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/OpalParticle.h"


// Class OpalParticle
// ------------------------------------------------------------------------

OpalParticle::OpalParticle()
{}


OpalParticle::OpalParticle
(double x, double px, double y, double py, double t, double pt)

{
    phase[X]  = x;
    phase[PX] = px;
    phase[Y]  = y;
    phase[PY] = py;
    phase[T]  = t;
    phase[PT] = pt;
}
