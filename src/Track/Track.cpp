// ------------------------------------------------------------------------
// $RCSfile: Track.cpp,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Struct: Track
//   This structure holds all data for tracking.
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:46 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Track/Track.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"
// Class Track
// ------------------------------------------------------------------------

Track *Track::block = 0;


/**

Track is asking the dictionary if already a
particle bunch was allocated. If that is the
case Track is using the already allocated bunch,
otherwise a new bunch is allocated in the dictionary.
*/


Track::Track(BeamSequence *u, const PartData &ref, double dt, int maxtsteps, int stepsperturn, double zStop, int timeintegrator, int nslices):
    reference(ref),
    use(u),
    dT(dt),
    zstop(zStop),
    maxTSteps(maxtsteps),
    stepsPerTurn(stepsperturn),
    timeIntegrator(timeintegrator),
    parser() {
    if(nslices > 0) {
        if(!OPAL.hasSLBunchAllocated())
            OPAL.setSLPartBunch(new EnvelopeBunch(&ref));

        slbunch = OPAL.getSLPartBunch();
    } else {
        if(!OPAL.hasBunchAllocated())
            OPAL.setPartBunch(new PartBunch(&ref));

        bunch = OPAL.getPartBunch();
    }
}


Track::~Track()
{}
