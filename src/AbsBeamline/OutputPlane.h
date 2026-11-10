// Copyright (c) 2023, Chris Rogers
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

#ifndef CLASSIC_OutputPlane_HH
#define CLASSIC_OutputPlane_HH

#include <list>
#include <string>
#include "AbsBeamline/PluginElement.h"

/** Class OutputPlane
 *  ------------------------------------------------------------------------
 *  Class OutputPlane defines an output plane, which enables output of tracking
 *  data at an arbitrary spatial plane. This is a much simplified version of
 *  the original OutputPlance class from old Opal.
 *
 */

class OutputPlane : public PluginElement {
public:
    using VectorView_t  = Kokkos::View<ippl::Vector<double, 3>*>;
    using BooleanView_t = Kokkos::View<bool*>;
    using DoubleView_t  = Kokkos::View<double*>;
    using Vector3D_t    = Vector_t<double, 3>;

    /** Constructor with given name. */
    explicit OutputPlane(const std::string& name);

    /** Default constructor (initialise to empty */
    OutputPlane();

    /** Copy constructor */
    OutputPlane(const OutputPlane&);

    /** Destructor */
    ~OutputPlane() override;

    /** Inheritable copy operation */
    ElementBase* clone() const override;

    /** Apply visitor to output plane. */
    void accept(BeamlineVisitor&) const override;

    /** @return The element type (ElementType::OUTPUTPLANE). */
    ElementType getType() const override;
    /** Return the cell geometry */
    Geometry& getGeometry() override { return geometry_m; }
    /** Return the cell geometry */
    const Geometry& getGeometry() const override { return geometry_m; }

    /** Return the longitudinal field-support extent.
     *
     *  The body interval latexmath:[[0, L]] plus the reach of the tanh fringe past each
     *  pole face, measured along the design arc. A hard edge (zero fringe lengths) gives
     *  the plain body extent.
     */
    void getFieldExtent(double& zBegin, double& zEnd) const override;

    /** Get the normal to the plane */
    Vector3D_t getNormal() const { return normal_m; }
    /** Set the normal to the plane */
    void setNormal(Vector3D_t normal) { normal_m = normal; }

    /** Get the centre of the plane */
    Vector3D_t getCentre() const { return centre_m; }
    /** Set the centre of the plane */
    void setCentre(const Vector3D_t centre) { centre_m = centre; }

    /** Get the radial extent from the plane centre */
    double getRadialExtent() const { return radialExtent_m; }
    /** Set the radial extent from the centre */
    void setRadialExtent(const double r) { radialExtent_m = r; }

    /** Record particles that have crossed the plane in this most recent step. */
    static void observe(PartBunch_t* bunch);
    static void store();

    // Types
    class StashedContainer {
    public:
        explicit StashedContainer(ParticleContainer_t* pc);
        void stash(ParticleContainer_t* pc);
        VectorView_t& r() { return r_m; }
        VectorView_t& p() { return p_m; }
        BooleanView_t& crossed() { return crossed_m; }
        DoubleView_t& t() { return t_m; }

    private:
        VectorView_t r_m;
        VectorView_t p_m;
        BooleanView_t crossed_m;
        DoubleView_t t_m;
    };

    /** Accessors for test cases */
    const std::vector<std::unique_ptr<StashedContainer>>& stash() const { return stash_m; }

private:
    // Overrides of PluginElement
    void doInitialise(PartBunch_t* bunch) override;
    void doFinalise() override;
    void doGoOffline() override;
    bool doCheck(PartBunch_t* bunch, int turnnumber, double t, double tstep) override;

protected:
    // Members
    Vector3D_t normal_m;           // normal to the output plane
    Vector3D_t centre_m;           // centre of the output plane
    double radialExtent_m = -1.0;  // maximum radial extent of the plane (circular)
    Geometry geometry_m;
    std::vector<std::unique_ptr<StashedContainer>> stash_m;
    static std::list<OutputPlane*> outputPlanes_m;

    // Helpers
    void stashParticles(const PartBunch_t* bunch);
    void detectCrossings(const PartBunch_t* bunch, double t, double tStep) const;
    void writeParticles(int turnNumber) const;
    KOKKOS_FUNCTION static void doDetectCrossing(
            size_t i, VectorView_t r1, VectorView_t r2, VectorView_t p, VectorView_t crossing,
            VectorView_t crossingP, BooleanView_t crossed, BooleanView_t invalid,
            DoubleView_t crossingT, double t, double tStep, Vector3D_t centre, Vector3D_t normal,
            double radialExtent, Vector3D_t refR);
};

#endif  // CLASSIC_OutputPlane_HH
