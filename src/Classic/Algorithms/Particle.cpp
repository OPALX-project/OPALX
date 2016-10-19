// ------------------------------------------------------------------------
// $RCSfile: Particle.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Particle
//   A Particle represents the phase space coordinates of a particle.
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

#include "Algorithms/Particle.h"


// Class Particle
// ------------------------------------------------------------------------

Particle::Particle()
{}


Particle::Particle
(double x, double px, double y, double py, double t, double pt)

{
    phase[X]  = x;
    phase[PX] = px;
    phase[Y]  = y;
    phase[PY] = py;
    phase[T]  = t;
    phase[PT] = pt;
}
