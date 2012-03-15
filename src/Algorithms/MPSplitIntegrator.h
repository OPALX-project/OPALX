#ifndef CLASSIC_MPSplitIntegrator_HH
#define CLASSIC_MPSplitIntegrator_HH

// ------------------------------------------------------------------------
// $RCSfile: MPSplitIntegrator.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: MPSplitIntegrator
//
// ------------------------------------------------------------------------
//
// $Date: 2000/03/27 09:33:36 $
// $Author: Andreas Adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/MapIntegrator.h"
#include <vector>

class BeamlineVisitor;
class BMultipoleField;
class Multipole;
class PartBunch;
class PartData;
class Particle;

template <class T, int N> class FVps;


// Class MPSplitIntegrator
// ------------------------------------------------------------------------
/// Integrator replacing each multipole by a set of thin lenses.
// The integrator may be used by algorithms to replace the normal
// procedure for traversing an element.  Phase space coordinates numbering:
// [tab 3 b]
// [row]number [&]name        [&]unit  [/row]
// [row]0      [&]x           [&]metres [/row]
// [row]1      [&]p_x/p_r     [&]1      [/row]
// [row]2      [&]y           [&]metres [/row]
// [row]3      [&]p_y/p_r     [&]1      [/row]
// [row]4      [&]v*delta_t   [&]metres [/row]
// [row]5      [&]delta_p/p_r [&]1      [/row]
// [/tab][p]
// Where $p_r$ is the constant reference momentum defining the reference
// frame velocity, $m$ is the rest mass of the particles, and $v$ is the
// instantaneous velocity of the particle.
// [p]
// Other units used:
// [tab 2 b]
// [row]quantity             [&]unit           [/row]
// [row]reference momentum   [&]electron-volts [/row]
// [row]velocity             [&]metres/second  [/row]
// [row]accelerating voltage [&]volts          [/row]
// [row]separator voltage    [&]volts          [/row]
// [row]frequencies          [&]hertz          [/row]
// [row]phase lags           [&]$2*pi$         [/row]
// [/tab][p]
// A MPSplitIntegrator performs integration through an element using two
// thin lenses of force 1/2, one placed at 1/6 and the other at 5/6 of
// the length respectively.

class MPSplitIntegrator: public MapIntegrator {

public:

    /// Constructor.
    // Attach this integrator to the given Multipole, using [b]slices[/b]
    // subdivisions.
    MPSplitIntegrator(Multipole *, int slices);

    MPSplitIntegrator(const MPSplitIntegrator &);
    virtual ~MPSplitIntegrator();

    /// Make clone.
    virtual MPSplitIntegrator *clone() const;

    /// Get geometry.
    //  Return the element geometry.
    //  Version for non-constant object.
    virtual Geometry &getGeometry();

    /// Get geometry.
    //  Return the element geometry
    //  Version for constant object.
    virtual const Geometry &getGeometry() const;

    /// Get element type string.
    virtual const string &getType() const;

    /// Get map from MPSplitIntegrator.
    // The map is returned in [b]map[/b], the other values are the same
    // as in the calling mapper.
    virtual void getMap(FVps<double, 6> &map, const PartData &data,
                        bool revBeam, bool revTrack) const;

    /// Track map through MPSplitIntegrator.
    // The map tracked is [b]map[/b], the other values are the same
    // as in the calling mapper.
    virtual void trackMap(FVps<double, 6> &map, const PartData &data,
                          bool revBeam, bool revTrack) const;

    /// Track particle through MPSplitIntegrator.
    // The particle tracked is [b]part[/b], the other values are the same
    // as in the calling mapper.
    virtual void trackParticle(Particle &part, const PartData &data,
                               bool revBeam, bool revTrack) const;

    /// Track particle bunch through MPSplitIntegrator.
    // The bunch tracked is [b]buch[/b], the other values are the same
    // as in the calling mapper.
    virtual void trackBunch(PartBunch &bunch, const PartData &data,
                            bool revBeam, bool revTrack) const;

    /// Return slice positions.
    // Build a vector [b]v[/b] containing the longitudinal positions
    // of the thin lens slices.
    void getSlices(std::vector<double> &v) const;

private:

    // Not implemented.
    MPSplitIntegrator();
    void operator=(const MPSplitIntegrator &);

    // Track a map through a particular element.
    void applyDrift(FVps<double, 6> &map, double, const PartData &) const;
    void applyMultipole(FVps<double, 6> &map,
                        const BMultipoleField &field, double factor) const;

    // Track a particle through a particular element.
    void applyDrift(Particle &, double, const PartData &) const;
    void applyMultipole(Particle &,
                        const BMultipoleField &field, double factor) const;

    // The embedded multipole.
    Multipole *itsMultipole;

    // The number of slices.
    int itsSlices;
};

#endif // CLASSIC_MPSplitIntegrator_HH
