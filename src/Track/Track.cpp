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
// #include "Algorithms/PartBunchBase.h"
#include "Algorithms/PartBunch.h" //FIXME
#ifdef ENABLE_AMR
    #include "Algorithms/AmrPartBunch.h"
#endif
#include "Algorithms/bet/EnvelopeBunch.h"
#include "AbstractObjects/OpalData.h"
#include "Utilities/Options.h"
// Class Track
// ------------------------------------------------------------------------

Track *Track::block = 0;
std::stack<Track*> Track::stashedTrack;

/**

Track is asking the dictionary if already a
particle bunch was allocated. If that is the
case Track is using the already allocated bunch,
otherwise a new bunch is allocated in the dictionary.
*/


Track::Track(BeamSequence *u, const PartData &ref, const std::vector<double> & dt,
             const std::vector<unsigned long long> & maxtsteps, int stepsperturn,
             double zStart, const std::vector<double> & zStop, int timeintegrator,
             int nslices, double t0, double dtScInit, double deltaTau):
    bunch(nullptr),
    slbunch(nullptr),
    reference(ref),
    use(u),
    parser(),
    dT(dt),
    dtScInit(dtScInit),
    deltaTau(deltaTau),
    t0_m(t0),
    localTimeSteps(maxtsteps),
    stepsPerTurn(stepsperturn),
    zstart(zStart),
    zstop(zStop),
    timeIntegrator(timeintegrator) {
    if(nslices > 0) {
        if(!OpalData::getInstance()->hasSLBunchAllocated())
            OpalData::getInstance()->setSLPartBunch(new EnvelopeBunch(&ref));

        if(!OpalData::getInstance()->hasBunchAllocated()) {           // we need this for Autophasing
#ifdef ENABLE_AMR
            if ( Options::amr )
                OpalData::getInstance()->setPartBunch(new AmrPartBunch(&ref));
            else
#endif
                OpalData::getInstance()->setPartBunch(new PartBunch(&ref));
        }

        bunch = OpalData::getInstance()->getPartBunch();
        slbunch = OpalData::getInstance()->getSLPartBunch();
    } else {
        if(!OpalData::getInstance()->hasBunchAllocated()) {
#ifdef ENABLE_AMR
            if ( Options::amr )
                OpalData::getInstance()->setPartBunch(new AmrPartBunch(&ref));
            else
#endif
                OpalData::getInstance()->setPartBunch(new PartBunch(&ref));
        }

        bunch = OpalData::getInstance()->getPartBunch();
    }

}


Track::~Track()
{}

void Track::stash() {
    PAssert_EQ(stashedTrack.size(), 0);

    stashedTrack.push(block);
    block = 0;
}

Track* Track::pop() {
    delete block;
    block = stashedTrack.top();
    stashedTrack.pop();

    return block;
}
