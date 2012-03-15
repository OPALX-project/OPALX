#ifndef CLASSIC_AbstractTracker_HH
#define CLASSIC_AbstractTracker_HH

// ------------------------------------------------------------------------
// $RCSfile: AbstractTracker.h,v $
// ------------------------------------------------------------------------
// $Revision: 1.1.1.1.2.1 $
// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: AbstractTracker
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
#include "Algorithms/PartBunch.h"
#include "Algorithms/PartData.h"
#include "Algorithms/Particle.h"

class PartBunch;


// Class AbstractTracker
// ------------------------------------------------------------------------
/// Track particles or bunches.
//  An abstract base class for all visitors capable of tracking particles
//  through a beam element.
//  This class redefines all visitXXX() methods for elements as pure
//  to force their implementation in derived classes.

class AbstractTracker: public DefaultVisitor {

public:

    // Particle coordinate numbers.
    enum { X, PX, Y, PY, T, PT };

    /// Constructor.
    //  The beam line to be tracked is [b]bl[/b].
    //  The particle reference data are taken from [b]data[/b].
    //  If [b]revBeam[/b] is true, the beam runs from s = C to s = 0.
    //  If [b]revTrack[/b] is true, we track against the beam.
    AbstractTracker(const Beamline &, const PartData &,
                    bool backBeam, bool backTrack);

    virtual ~AbstractTracker();


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

    /// Apply the algorithm to a ParallelPlate.
    virtual void visitParallelPlate(const ParallelPlate &) = 0;

    /// Apply the algorithm to a CyclotronValley.
    virtual void visitCyclotronValley(const CyclotronValley &) = 0;

    /// Apply the algorithm to an align wrapper..
    virtual void visitAlignWrapper(const AlignWrapper &) = 0;

protected:

    /// The reference information.
    const PartData itsReference;

private:

    // Not implemented.
    AbstractTracker();
    AbstractTracker(const AbstractTracker &);
    void operator=(const AbstractTracker &);
};

#endif // CLASSIC_AbstractTracker_HH
