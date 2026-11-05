#ifndef OPALX_BeamBeam_HH
#define OPALX_BeamBeam_HH

// ------------------------------------------------------------------------
// Copyright: see Copyright.readme
// ------------------------------------------------------------------------
//
// Class: BeamBeam
//   Defines the interface for a placed collective beam-beam element.
//
// ------------------------------------------------------------------------
// Class category: AbsBeamline
// ------------------------------------------------------------------------
//

#include "AbsBeamline/ElementBase.h"

/**
 * @brief Placed beam-beam element; creates collective behavior for each tracking run.
 *
 * Geometry and input attributes belong to the element, while BeamBeamInteraction
 * owns mutable state. The interaction point is the element midpoint. The current
 * model supports one interaction in a single-pass line, with a physical primary,
 * its optional same-charge mirrored primary, and passive witness containers.
 * The fixed numerical field window is defined by BEAMBEAM::fieldWindowLength;
 * it is distinct from the element's geometric length. See BeamBeamInteraction
 * and BEAMBEAM::Config for field and rigid-primary policies.
 */
class BeamBeam : public ElementBase {
public:
    /// Constructor with given name.
    /// @param name Placed element name.
    explicit BeamBeam(const std::string& name);

    /// Construct an unnamed element.
    BeamBeam();
    /// Copy placement/attributes, not per-run interaction state.
    /// @param right Element whose configuration is copied.
    BeamBeam(const BeamBeam& right);
    virtual ~BeamBeam();

    /// Apply visitor to BeamBeam.
    /// @param visitor Beamline visitor dispatched through visitBeamBeam().
    virtual void accept(BeamlineVisitor& visitor) const override;

    /// Borrow the tracking bunch without taking ownership.
    /// @param bunch Bunch whose lifetime encloses element use.
    virtual void initialise(PartBunch_t* bunch) override;

    /// No element-owned finalization; runtime state belongs to the interaction.
    virtual void finalise() override;

    virtual ElementType getType() const override;

    /// Create the per-run collective behavior for this placed element.
    /// @return Caller-owned runtime object borrowing this element, which must outlive it.
    std::unique_ptr<ElementInteraction> createInteraction() const override;

    virtual void getFieldExtent(double& zBegin, double& zEnd) const override;

    /// Store legacy map-slice metadata; collective time tracking does not use it.
    /// @param nSlices Map-slice metadata value.
    void setNSlices(const std::size_t& nSlices);

    /// @return Legacy map-slice metadata (not a collective-solver discretization).
    std::size_t getNSlices() const;

    /// @return One required tracking step; the tracker controls the physical time step.
    virtual int getRequiredNumberOfTimeSteps() const override;

private:
    std::size_t nSlices_m;

    // Not implemented.
    void operator=(const BeamBeam&);
};

inline int BeamBeam::getRequiredNumberOfTimeSteps() const { return 1; }

#endif  // OPALX_BeamBeam_HH
