// ------------------------------------------------------------------------
// $RCSfile: AlignRemover.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AlignRemover
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:40 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Errors/AlignRemover.h"
#include "AbsBeamline/AlignWrapper.h"
#include "BeamlineGeometry/Euclid3D.h"


// Class AlignRemover
// ------------------------------------------------------------------------

AlignRemover::AlignRemover()
{}


AlignRemover::~AlignRemover()
{}


void AlignRemover::misalignment(const AlignWrapper &wrap, int)
{
  wrap.offset() = Euclid3D();
}
