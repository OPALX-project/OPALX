#ifndef CLASSIC_AbstractMapper_HH
#define CLASSIC_AbstractMapper_HH

// ------------------------------------------------------------------------
// $RCSfile: AbstractMapper.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AbstractMapper
//
// ------------------------------------------------------------------------
// Class category: Algorithms
// ------------------------------------------------------------------------
//
// $Date: 2004/11/12 18:57:53 $
// $Author: adelmann $
//
// ------------------------------------------------------------------------

#include "Algorithms/DefaultVisitor.h"
#include "Algorithms/PartData.h"

class BMultipoleField;
class Euclid3D;
template <class T, int N> class LinearMap;
template <class T, int N> class FTps;
template <class T, int N> class FVps;


// Class AbstractMapper
// ------------------------------------------------------------------------
/// Build transfer map.
//  An abstract visitor class implementing the default behaviour for all
//  visitors capable of tracking a transfer map through a beam line.
//  It implements access to the accumulated map, and keeps track of the
//  beam reference data.
//  This class redefines all visitXXX() methods for elements as pure
//  to force their implementation in derived classes.

class AbstractMapper: public DefaultVisitor {

public:

    // Particle coordinate numbers.
    enum { X, PX, Y, PY, T, PT };

    /// Constructor.
    //  The beam line to be tracked is [b]bl[/b].
    //  The particle reference data are taken from [b]data[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    AbstractMapper(const Beamline &bl, const PartData &data,
                   bool revBeam, bool revTrack);

    virtual ~AbstractMapper();


    /// Return the linear part of the accumulated map.
    virtual void getMap(LinearMap<double, 6> &) const = 0;

    /// Return the full map accumulated so far.
    virtual void getMap(FVps<double, 6> &) const = 0;

    /// Reset the linear part of the accumulated map for restart.
    virtual void setMap(const LinearMap<double, 6> &) = 0;

    /// Reset the full map for restart.
    virtual void setMap(const FVps<double, 6> &) = 0;


    /// Apply the algorithm to a beam-beam.
    virtual void visitBeamBeam(const BeamBeam &) = 0;

    /// Apply the algorithm to a collimator.
    virtual void visitCollimator(const Collimator &) = 0;

    /// Apply the algorithm to an arbitrary component.
    virtual void visitComponent(const Component &) = 0;

    /// Apply the algorithm to a corrector.
    virtual void visitCorrector(const Corrector &) = 0;

    /// Apply the algorithm to a diagnostic.
    virtual void visitDiagnostic(const Diagnostic &) = 0;

    /// Apply the algorithm to a drift.
    virtual void visitDrift(const Drift &) = 0;

    /// Apply the algorithm to a Lambertson.
    virtual void visitLambertson(const Lambertson &) = 0;

    /// Apply the algorithm to a marker.
    virtual void visitMarker(const Marker &) = 0;

    /// Apply the algorithm to a monitor.
    virtual void visitMonitor(const Monitor &) = 0;

    /// Apply the algorithm to a multipole.
    virtual void visitMultipole(const Multipole &) = 0;

    /// Apply the algorithm to a patch.
    virtual void visitPatch(const Patch &pat) = 0;

    /// Apply the algorithm to a probe.
    virtual void visitProbe(const Probe &prob) = 0;

    /// Apply the algorithm to a rectangular bend.
    virtual void visitRBend(const RBend &) = 0;

    /// Apply the algorithm to a RF cavity.
    virtual void visitRFCavity(const RFCavity &) = 0;

    /// Apply the algorithm to a RF quadrupole.
    virtual void visitRFQuadrupole(const RFQuadrupole &) = 0;

    /// Apply the algorithm to a sector bend.
    virtual void visitSBend(const SBend &) = 0;

    /// Apply the algorithm to a separator.
    virtual void visitSeparator(const Separator &) = 0;

    /// Apply the algorithm to a septum.
    virtual void visitSeptum(const Septum &) = 0;

    /// Apply the algorithm to a solenoid.
    virtual void visitSolenoid(const Solenoid &) = 0;

protected:

    /// Construct the vector potential for a Multipole.
    FTps<double, 6> buildMultipoleVectorPotential(const BMultipoleField &);

    /// Construct the vector potential for an SBend.
    FTps<double, 6> buildSBendVectorPotential(const BMultipoleField &, double h);

    /// The reference information.
    const PartData itsReference;

private:

    // Not implemented.
    AbstractMapper();
    AbstractMapper(const AbstractMapper &);
    void operator=(const AbstractMapper &);
};

#endif // CLASSIC_AbstractMapper_HH
