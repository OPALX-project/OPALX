//
// Class Collimator
//   Interface for a transverse-aperture collimator.
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
#ifndef OPALX_Collimator_HH
#define OPALX_Collimator_HH

#include "AbsBeamline/ElementBase.h"

/**
 * @class Collimator
 * @brief Field-free element that removes particles hitting its aperture.
 *
 * A collimator carries no field. Everything it does is inherited: the tracker
 * calls ElementBase::markOutsideAperture() once per time step for every element
 * near the bunch, which marks particles inside the body window [0, L) but
 * outside the aperture in the container's InvalidMask. Those marks are flushed
 * by the single end-of-step deletion in ParallelTracker::execute(), so absorbed
 * particles never appear in any output.
 *
 * @note The aperture shape comes from the common APERTURE attribute stored in
 *       ElementBase. OpalCollimator requires both APERTURE and L > 0 to be
 *       written in the inputfile, since a collimator without either does nothing.
 */
class Collimator : public ElementBase {
public:
    /// Constructor with given name.
    explicit Collimator(const std::string& name);

    Collimator();
    Collimator(const Collimator& right);
    virtual ~Collimator();

    /// Apply visitor to Collimator.
    virtual void accept(BeamlineVisitor&) const override;

    virtual void initialise(PartBunch_t* bunch) override;

    virtual void finalise() override;

    virtual ElementType getType() const override;

    virtual void getFieldExtent(double& zBegin, double& zEnd) const override;

    virtual int getRequiredNumberOfTimeSteps() const override;

private:
    // Not implemented.
    void operator=(const Collimator&);
};

inline int Collimator::getRequiredNumberOfTimeSteps() const { return 1; }

#endif  // OPALX_Collimator_HH
