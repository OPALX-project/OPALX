#ifndef CLASSIC_Integrator_HH
#define CLASSIC_Integrator_HH

// ------------------------------------------------------------------------
// $RCSfile: Integrator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: Integrator
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:31 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "AbsBeamline/ElementBase.h"
#include "MemoryManagement/Pointer.h"

class PartBunch;
class PartData;

namespace classic {
class Particle;

}

template <class T, int N> class FVps;


// Class Integrator
// ------------------------------------------------------------------------
/// Base class for special integrators.
//  Integrator is a pure abstract class.  It forms the base class for all
//  special propagators through an element or a beam line.  If present, it
//  overrides a complete Component or Beamline object.

class Integrator: public ElementBase {

public:

    explicit Integrator(ElementBase *);
    Integrator(const Integrator &rhs);
    virtual ~Integrator();

    /// Return the embedded element.
    inline ElementBase *getElement() const;

    /// Set sharable flag.
    //  The whole structure depending on [b]this[/b] is marked as sharable.
    //  After this call a [b]copyStructure()[/b] call reuses the element.
    virtual void makeSharable();

    /// Track a particle.
    //  The first argument describes the particle's phase space position,
    //  the second argument describes the particle's momentum and mass,
    //  [b]revBeam[/b] true, means that the beam runs backwards, and
    //  [b]revTrack[/b] true, means that we track against the beam.
    virtual void trackParticle(classic::Particle &, const PartData &,
                               bool revBeam, bool revTrack) const = 0;

    /// Track a particle bunch.
    //  The first argument describes the particles in the bunch,
    //  the second argument describes the reference momentum and mass,
    //  [b]revBeam[/b] true, means that the beam runs backwards, and
    //  [b]revTrack[/b] true, means that we track against the beam.
    virtual void trackBunch(PartBunch &, const PartData &,
                            bool revBeam, bool revTrack) const = 0;

    /// Track a map.
    //  The first argument describes the map to be tracked,
    //  the second argument describes the reference momentum and mass,
    //  [b]revBeam[/b] true, means that the beam runs backwards, and
    //  [b]revTrack[/b] true, means that we track against the beam.
    virtual void trackMap(FVps<double, 6> &, const PartData &,
                          bool revBeam, bool revTrack) const = 0;

protected:

    /// Pointer to the replaced element.
    Pointer<ElementBase> itsElement;

private:

    // Not implemented.
    Integrator();
    void operator=(const Integrator &);
};


// Implementation.
// ------------------------------------------------------------------------

inline ElementBase *Integrator::getElement() const {
    return &*itsElement;
}

#endif // CLASSIC_Integrator_HH
