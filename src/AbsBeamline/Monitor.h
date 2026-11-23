//
// Class Monitor
//   Defines the abstract interface for a beam position monitor.
//
// Copyright (c) 2026, Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// This file is part of OPAL.
//
// OPAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// You should have received a copy of the GNU General Public License
// along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
//
#ifndef OPALX_Monitor_HH
#define OPALX_Monitor_HH

#include "AbsBeamline/ElementBase.h"
#include "PartBunch/PartBunch.h"
#include "Structure/LossDataSink.h"

#include <cmath>
#include <memory>
#include <string>

class Monitor : public ElementBase {
public:
    /// Constructor with given name.
    explicit Monitor(const std::string& name);

    Monitor();
    Monitor(const Monitor&);
    ~Monitor() override;

    virtual void apply(
            const std::shared_ptr<ParticleContainer_t>& pc, const StepTime& time) override;

    virtual void apply(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    virtual bool applyToReferenceParticle(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    /// Takes ASCIIDUMP, COMPUTEPERCENTILES, the restart flag and the open mode of the loss file.
    void setRunValues(
            const opalx::run::OptionValues& options,
            const opalx::run::RunContext& context) override;

    virtual void initialise(PartBunch_t* bunch) override;

    virtual void goOffline() override;

    virtual ElementType getType() const override;

    virtual void getFieldExtent(double& zBegin, double& zEnd) const override;

    void setCollectionType(CollectionType type);

    /**
     * @brief Whether a reference particle that crosses the monitor is recorded.
     * @note The tracker turns this off for the reference start and the reference pass; only
     *       the step loop records.
     */
    void setRecordReference(bool record) { recordReference_m = record; }

    virtual int getRequiredNumberOfTimeSteps() const override;

    virtual bool isInside(const Vector_t<double, 3>& r) const override;

    /// The isInside() window |z| <= L/2 around the monitor plane joined with the body [0, L).
    BoundingBox getLocalSupportBox(double reach) const override;

private:
    void driftToCorrectPositionAndSave(const Vector_t<double, 3>& R, const Vector_t<double, 3>& P);

    Monitor& operator=(const Monitor&) = delete;

    std::string filename_m; /**< The name of the output file */
    CollectionType type_m;
    unsigned int numPassages_m;
    bool recordReference_m = true;  ///< See setRecordReference().

    // Run values, see setRunValues(); the defaults are those of a fresh run without OPTION.
    bool asciiDump_m          = false;
    bool computePercentiles_m = false;
    bool restart_m            = false;
    bool appendOutput_m       = false;  ///< Append to an existing loss file (restart).

    std::unique_ptr<LossDataSink> lossDs_m;

    static const double halfLength_s;
};

inline void Monitor::setCollectionType(CollectionType type) { type_m = type; }

inline int Monitor::getRequiredNumberOfTimeSteps() const { return 1; }

inline bool Monitor::isInside(const Vector_t<double, 3>& r) const {
    const double length = getGeometry().getElementLength();
    return std::abs(r(2)) <= 0.5 * length && ApertureHelper::isInsideAperture(r, aperture_m);
}

#endif  // OPALX_Monitor_HH
