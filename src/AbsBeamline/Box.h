//
// Class Box
//   Interface for a rectangular block that absorbs every particle inside it.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL. If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPALX_Box_HH
#define OPALX_Box_HH

#include "AbsBeamline/ElementBase.h"

/**
 * @class Box
 * @brief Field-free solid block that removes every particle inside it.
 *
 * The block occupies |x| < halfWidth, |y| < halfHeight, 0 <= z < L in the
 * element-local frame. The boundary counts as outside. It is the opposite of
 * a Collimator, which removes particles outside its aperture: a Box models a
 * piece of material such as a collimator jaw, which a particle can pass around.
 * Offset and rotation come from the usual element placement.
 *
 * @note The block size is kept apart from the common aperture (aperture_m).
 *       ElementBase::isInside only registers an element along the reference
 *       path if the reference particle is inside the aperture, and a block that
 *       sits beside the reference path must still be registered. So the
 *       aperture stays at the wide-open default.
 */
class Box : public ElementBase {
public:
    /// Constructor with given name.
    explicit Box(const std::string& name);

    Box();
    Box(const Box& right);
    virtual ~Box();

    /// Apply visitor to Box.
    virtual void accept(BeamlineVisitor&) const override;

    virtual void initialise(PartBunch_t* bunch) override;

    virtual void finalise() override;

    virtual ElementType getType() const override;

    virtual void getFieldExtent(double& zBegin, double& zEnd) const override;

    virtual int getRequiredNumberOfTimeSteps() const override;

    /**
     * @brief Mark particles inside the block in the container's InvalidMask.
     *
     * The tracker calls this for every element near the bunch, with the bunch
     * already moved into the element-local frame. Despite the inherited name,
     * a Box marks the particles inside its body, not outside an aperture.
     *
     * @param pc Particle container in the element-local frame.
     * @return Number of particles newly marked on this rank.
     */
    virtual size_t markOutsideAperture(const std::shared_ptr<ParticleContainer_t>& pc) override;

    /**
     * @brief Set the transverse size of the block.
     * @param halfWidth  Half of the full width along local x [m].
     * @param halfHeight Half of the full height along local y [m].
     */
    void setHalfSize(double halfWidth, double halfHeight);

    double getHalfWidth() const;

    double getHalfHeight() const;

private:
    // Not implemented.
    void operator=(const Box&);

    /// Half of the block width along local x [m].
    double halfWidth_m;

    /// Half of the block height along local y [m].
    double halfHeight_m;
};

inline int Box::getRequiredNumberOfTimeSteps() const { return 1; }

inline void Box::setHalfSize(double halfWidth, double halfHeight) {
    halfWidth_m  = halfWidth;
    halfHeight_m = halfHeight;
}

inline double Box::getHalfWidth() const { return halfWidth_m; }

inline double Box::getHalfHeight() const { return halfHeight_m; }

#endif  // OPALX_Box_HH
