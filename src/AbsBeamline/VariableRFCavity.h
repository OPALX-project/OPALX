//
// Class VariableRFCavity
//   Defines the abstract interface for a RF Cavity
//   with Time Dependent Parameters.
//
// Copyright (c) 2014 - 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
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
#ifndef ABSBEAMLINE_VARIABLERFCAVITY_HH
#define ABSBEAMLINE_VARIABLERFCAVITY_HH

#include "AbsBeamline/Component.h"
#include "Algorithms/AbstractTimeDependence.h"
#include "BeamlineGeometry/StraightGeometry.h"
#include "Fields/EMField.h"
#include "Physics/Physics.h"

/** @class VariableRFCavity
 *
 *  Generates a field like
 *      E = a(t) * sin{2*pi * integral(f(tau),0,t) + q(t)}
 *      B = 0
 *  where E0, B0 are user defined fields, a(t), f(t), q(t) are time
 *  dependent amplitude, frequency, phase respectively; it is assumed that these
 *  quantities vary sufficiently slowly that Maxwell is satisfied.
 *
 *  The time dependent quantities are
 */
class VariableRFCavity : public Component {
public:
    /// Constructor with given name.
    explicit VariableRFCavity(const std::string& name);
    /** Copy Constructor; performs deepcopy on time-dependence models */
    VariableRFCavity(const VariableRFCavity&);
    /** Default constructor */
    VariableRFCavity();
    /** Assignment operator; performs deepcopy on time-dependence models*/
    VariableRFCavity& operator=(const VariableRFCavity&);
    /** Destructor does nothing
     *
     * The shared_ptrs will self-destruct when reference count goes to 0
     */
    ~VariableRFCavity() override = default;

    /** Apply visitor to RFCavity.
     *
     *  The RF cavity finds the "time dependence" models by doing a string
     *  lookup against a list held by AbstractTimeDependence at accept time.
     */
    void accept(BeamlineVisitor&) const override;

    /** Inheritable deepcopy method */
    ElementBase* clone() const override;

    /** Calculate the field for all particles */
    bool apply(const std::shared_ptr<ParticleContainer_t>& pc) override;

    /** Calculate the field at the position of the i^th particle
     *
     *  @param i indexes the particle whose field we need
     *  @param t the time at which the field is calculated
     *  @param E return value with electric field strength
     *  @param B return value with magnetic field strength
     *
     *  @returns True if particle is outside the boundaries; else False
     */
    bool apply(const size_t& i, const double& t, Vector_t<double, 3>& E, Vector_t<double, 3>& B)
            override;

    /** Calculate the field at a given position
     *
     *  @param R the position at which the field is calculated
     *  @param P the momentum (not used)
     *  @param t the time at which the field is calculated
     *  @param E return value; filled with electric field strength
     *  @param B return value; filled with magnetic field strength
     *
     *  @returns True if particle is outside the boundaries; else False
     */
    bool apply(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    /** Calculate the field at a given position. This is identical to "apply".
     *
     *  @param R the position at which the field is calculated
     *  @param P the momentum (not used)
     *  @param t the time at which the field is calculated
     *  @param E return value; filled with electric field strength
     *  @param B return value; filled with magnetic field strength
     *
     *  @returns True if particle is outside the boundaries; else False
     */
    bool applyToReferenceParticle(
            const Vector_t<double, 3>& R, const Vector_t<double, 3>& P, const double& t,
            Vector_t<double, 3>& E, Vector_t<double, 3>& B) override;

    /** Initialise ready for tracking
     *
     *  Just sets RefPartBunch_m
     */
    void initialise(PartBunch_t* bunch, double& startField, double& endField) override;

    /** Finalise following tracking
     *
     *  Just sets RefPartBunch_m to nullptr
     */
    void finalise() override;

    /** @returns false (cavity does not bend the trajectory) */
    bool bends() const override { return false; }

    /** Return the longitudinal field-support extent.
     *
     *  For a VariableRFCavity the field support coincides with the nominal
     *  body extent, i.e. the interval
     *  latexmath:[z \in [0, L)] in the local element frame.
     */
    void getFieldExtend(double& zBegin, double& zEnd) const override {
        zBegin = 0.0;
        zEnd   = getElementLength();
    }

    /** Get the amplitude at a given time
     *
     *  @param time: the time at which the amplitude is calculated
     *
     *  @returns the RF field gradient.
     */
    virtual double getAmplitude(const double time) const { return amplitudeTD_m->getValue(time); }

    /** Get the frequency at a given time
     *
     *  @param time: the time at which the frequency is calculated
     *
     *  @returns the RF cavity frequency.
     */
    virtual double getFrequency(const double time) const { return frequencyTD_m->getValue(time); }

    /** Get the phase at a given time
     *
     *  @param time: the time at which the phase is calculated
     *
     *  @returns the RF cavity phase.
     */
    virtual double getPhase(const double time) const { return phaseTD_m->getValue(time); }

    /** @returns the full height of the cavity */
    virtual double getHeight() const { return halfHeight_m * 2; }
    /** @returns the full width of the cavity */
    virtual double getWidth() const { return halfWidth_m * 2; }
    /** @returns the length of the cavity */
    virtual double getLength() const { return length_m; }
    /** Set the full height of the cavity */
    virtual void setHeight(const double fullHeight) { halfHeight_m = fullHeight / 2; }
    /** Set the full width of the cavity */
    virtual void setWidth(const double fullWidth) { halfWidth_m = fullWidth / 2; }
    /** Set the length of the cavity */
    virtual void setLength(double length);

    /** @returns shared_ptr to the amplitude (field gradient) time dependence */
    virtual std::shared_ptr<AbstractTimeDependence> getAmplitudeModel() const;
    /** @returns shared_ptr to the phase time dependence */
    virtual std::shared_ptr<AbstractTimeDependence> getPhaseModel() const;
    /** @returns shared_ptr to the frequency */
    virtual std::shared_ptr<AbstractTimeDependence> getFrequencyModel() const;

    /** Set the amplitude (field gradient) time dependence */
    virtual void setAmplitudeModel(std::shared_ptr<AbstractTimeDependence> amplitude_td);
    /** Set the phase time dependence */
    virtual void setPhaseModel(std::shared_ptr<AbstractTimeDependence> phase_td);
    /** Set the frequency time dependence */
    virtual void setFrequencyModel(std::shared_ptr<AbstractTimeDependence> frequency_td);

    /** Set the amplitude time dependence name
     *
     *  The name is used to find the amplitude model at accept time
     */
    virtual void setAmplitudeName(const std::string& amplitude) { amplitudeName_m = amplitude; }

    /** Set the phase time dependence name
     *
     *  The name is used to find the phase model at accept time
     */
    virtual void setPhaseName(const std::string& phase) { phaseName_m = phase; }

    /** Set the frequency time dependence name
     *
     *  The name is used to find the frequency model at accept time
     */
    virtual void setFrequencyName(const std::string& frequency) { frequencyName_m = frequency; }

    /** Set the cavity geometry */
    StraightGeometry& getGeometry() override;
    /** @returns the cavity geometry */
    const StraightGeometry& getGeometry() const override;

    /** Lookup the time dependencies and update.
     *
     *  Throws if the named time dependencies can't be found. Also throws if the
     *  width or height is < 1 nm
     */
    void initialiseTimeDependencies() const;

    /// Not implemented
    EMField& getField() override;
    /// Not implemented
    const EMField& getField() const override;

protected:
    void initNull();
    mutable std::shared_ptr<AbstractTimeDependence> phaseTD_m;
    mutable std::shared_ptr<AbstractTimeDependence> amplitudeTD_m;
    mutable std::shared_ptr<AbstractTimeDependence> frequencyTD_m;
    std::string phaseName_m;
    std::string amplitudeName_m;
    std::string frequencyName_m;
    double halfWidth_m;
    double halfHeight_m;
    double length_m;

    /// The cavity's geometry.
    StraightGeometry geometry;

    /* The host/device compute function */
    static KOKKOS_INLINE_FUNCTION bool computeField(
            const Vector_t<double, 3>& R, Vector_t<double, 3>& E, double E0, double integralF,
            double phi, double halfWidth, double halfHeight, double length);
};

KOKKOS_INLINE_FUNCTION bool VariableRFCavity::computeField(
        const Vector_t<double, 3>& R, Vector_t<double, 3>& E, const double E0,
        const double integralF, const double phi, const double halfWidth, const double halfHeight,
        const double length) {
    if (R[2] < 0.0 || R[2] >= length) {
        return false;
    }
    E[2] += E0 * Kokkos::sin(Physics::two_pi * integralF + phi);
    const bool outsideAperture = Kokkos::abs(R[0]) > halfWidth || Kokkos::abs(R[1]) > halfHeight;
    return outsideAperture;
}

#endif  // ABSBEAMLINE_VARIABLERFCAVITY_HH
