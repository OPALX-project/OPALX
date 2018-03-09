#ifndef CLASSIC_MapIntegrator_HH
#define CLASSIC_MapIntegrator_HH

// ------------------------------------------------------------------------
// $RCSfile: MapIntegrator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MapIntegrator
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:32:33 $
// $Author: fci $
//
// ------------------------------------------------------------------------

#include "Algorithms/TrackIntegrator.h"

template <class T, unsigned Dim>
class PartBunchBase;
class PartData;
class OpalParticle;

template <class T, int N> class FVps;


// Class MapIntegrator
// ------------------------------------------------------------------------
/// Integrate a map.
//  An abstract base class for all integrators capable of tracking a
//  transfer map through a beam element.
//  Implements some default behaviour for such integrators.


class MapIntegrator: public TrackIntegrator {

public:

    explicit MapIntegrator(ElementBase *);
    MapIntegrator(const MapIntegrator &right);
    virtual ~MapIntegrator();


    /// Apply visitor.
    virtual void accept(BeamlineVisitor &visitor) const;

    /// Make a clone.
    virtual MapIntegrator *clone() const = 0;

    /// Track a particle.
    //  The particle is stored in [b]part[/b].
    //  The particle reference data are taken from [b]data[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    virtual void trackParticle(OpalParticle &part, const PartData &,
                               bool revBeam, bool revTrack) const;

    /// Track a particle bunch.
    //  The bunch is stored in [b]bunch[/b].
    //  The particle reference data are taken from [b]data[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    virtual void trackBunch(PartBunchBase<double, 3> *, const PartData &,
                            bool revBeam, bool revTrack) const;

    /// Track a map.
    //  The map is stored in [b]map[/b].
    //  The particle reference data are taken from [b]data[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    virtual void trackMap(FVps<double, 6> &, const PartData &,
                          bool revBeam, bool revTrack) const;

private:

    // Not implemented.
    MapIntegrator();
    void operator=(const MapIntegrator &);

    // Return the map for this MapIntegrator.
    virtual void getMap(FVps<double, 6> &, const PartData &,
                        bool revBeam, bool revTrack) const = 0;
};

#endif // CLASSIC_MapIntegrator_HH
