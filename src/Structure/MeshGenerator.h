//
// Mesh Generator
//
// Copyright (c) 2008-2020
// Paul Scherrer Institut, Villigen PSI, Switzerland
// All rights reserved.
//
// OPAL is licensed under GNU GPL version 3.
//

#ifndef MESHGENERATOR_H_
#define MESHGENERATOR_H_

#include <string>
#include <utility>
#include <vector>
#include "OPALTypes.h"

class ElementBase;

class MeshData {
public:
    std::vector<Vector_t<double, 3>> vertices_m;
    std::vector<Vector_t<unsigned int, 3>> triangles_m;
    std::vector<std::pair<Vector_t<double, 3>, Vector_t<double, 3>>> decorations_m;
    int type_m;
};

class MeshGenerator {
public:
    MeshGenerator();

    /**
     * @brief Extract a transverse support size suitable for placement/export meshes.
     *
     * The returned radii are in the element-local transverse frame. This is used
     * both for direct meshing and for deriving a representative drift radius from
     * the first non-drift element in the lattice.
     *
     * @param element Beamline element to inspect.
     * @param minor Output minor transverse radius.
     * @param major Output major transverse radius.
     * @return true if a finite support size is available, false otherwise.
     */
    static bool getTransverseSupport(const ElementBase& element, double& minor, double& major);

    /**
     * @brief Set the drift support radius used when meshing drift spaces.
     *
     * Drifts do not carry their own physical support envelope. For visualisation
     * we derive a representative cylindrical radius from the first non-drift
     * element in the beamline and reuse it for all drifts.
     *
     * @param minor Minor transverse radius.
     * @param major Major transverse radius.
     */
    void setDriftReference(double minor, double major);

    void add(const ElementBase& element);

    void write(const std::string& fname);

private:
    enum MeshType {
        OTHER = 0,
        DIPOLE,
        QUADRUPOLE,
        SEXTUPOLE,
        OCTUPOLE,
        SOLENOID,
        RFCAVITY,
        TRAVELINGWAVE,
        DRIFT
    };

    static MeshData getCylinder(
            double length, double minor, double major, double formFactor,
            const unsigned int numSegments = 36);

    /**
     * @brief Build a hollow tube aligned with the local z-axis.
     *
     * The tube is represented by triangulated inner and outer walls without end
     * caps. This is used for drift sections to visualise the beam pipe rather
     * than a solid support body.
     *
     * @param length Tube length along the local z-axis.
     * @param innerMinor Inner semi-axis in the local y-direction.
     * @param innerMajor Inner semi-axis in the local x-direction.
     * @param outerMinor Outer semi-axis in the local y-direction.
     * @param outerMajor Outer semi-axis in the local x-direction.
     * @param numSegments Number of azimuthal segments.
     * @return Hollow triangulated tube mesh.
     */
    static MeshData getTube(
            double length, double innerMinor, double innerMajor, double outerMinor,
            double outerMajor, const unsigned int numSegments = 36);

    /**
     * @brief Build a curved bend body by sweeping an elliptical ring along a centerline.
     *
     * The generated surface follows the sampled bend design path and keeps a
     * constant transverse support ellipse. This is the first visualization mesh
     * for analytic bends and therefore models the bend body as a smooth curved
     * tube without end caps.
     *
     * @param centerline Sampled bend centerline in the local body chart.
     * @param minor Outer semi-axis in the local y-direction.
     * @param major Outer semi-axis in the local x-direction.
     * @param numSegments Number of azimuthal segments per ring.
     * @return Triangulated curved bend body mesh.
     */
    static MeshData getBendTube(
            const std::vector<Vector_t<double, 3>>& centerline, double minor, double major,
            const unsigned int numSegments = 24);

    /**
     * @brief Build a quadrupole-like body from four longitudinal pole blocks.
     *
     * The pole blocks are scaled from the element aperture and extruded along
     * the full element length so that the resulting shape remains meaningful
     * when the quadrupole length changes.
     *
     * @param length Body length along the local z-axis.
     * @param minor Aperture semi-axis in the local y-direction.
     * @param major Aperture semi-axis in the local x-direction.
     * @param formFactor Exit scaling factor used for conic apertures.
     * @return Triangulated quadrupole body mesh.
     */
    static MeshData getQuadrupole(double length, double minor, double major, double formFactor);

    /**
     * @brief Build a rounded prism from a superellipse cross-section.
     *
     * This is used for visual bodies that should read as machined pole pieces
     * rather than sharp-edged boxes. The cross-section is
     * @f[
     *   \left|\frac{x}{a}\right|^n + \left|\frac{y}{b}\right|^n = 1
     * @f]
     * with @f$n > 2@f$ so that the shape has flat-ish faces and rounded edges.
     *
     * @param length Body length along the local z-axis.
     * @param halfWidth Cross-section half-width in the local x-direction.
     * @param halfHeight Cross-section half-height in the local y-direction.
     * @param formFactor Exit scaling factor used for conic apertures.
     * @param exponent Superellipse exponent.
     * @param numSegments Number of azimuthal segments.
     * @return Triangulated rounded prism mesh.
     */
    static MeshData getSuperellipsePrism(
            double length, double halfWidth, double halfHeight, double formFactor,
            double exponent = 4.0, unsigned int numSegments = 28);

    /**
     * @brief Build a pointed quadrupole pole tip with a rounded nose.
     *
     * The tip is an extruded 2D profile that narrows toward the beam aperture.
     * A short rounded nose is added by collapsing the final profile width to a
     * finite tip radius instead of a sharp point.
     *
     * @param length Body length along the local z-axis.
     * @param baseHalfWidth Half-width of the pole root.
     * @param baseHalfHeight Half-height of the pole root.
     * @param tipReach Longitudinal reach of the inward-facing nose.
     * @param tipRadius Half-width of the rounded nose.
     * @param formFactor Exit scaling factor used for conic apertures.
     * @return Triangulated pole-tip prism mesh.
     */
    static MeshData getPoleTipPrism(
            double length, double baseHalfWidth, double baseHalfHeight, double tipReach,
            double tipRadius, double formFactor);

    /**
     * @brief Build a solenoid body as a hollow tube with short end collars.
     *
     * The body is scaled from the transverse support envelope and extruded over
     * the full element length. Short end collars make the solenoid visually
     * distinct from both drifts and plain cylindrical supports.
     *
     * @param length Body length along the local z-axis.
     * @param minor Outer semi-axis in the local y-direction.
     * @param major Outer semi-axis in the local x-direction.
     * @return Triangulated solenoid body mesh.
     */
    static MeshData getSolenoid(double length, double minor, double major);

    /**
     * @brief Build a standing-wave cavity body from a sequence of bulged cells.
     *
     * The body is assembled from hollow tube segments with alternating neck and
     * cell radii. This produces a programmer-facing standing-wave cavity shape
     * while remaining purely geometric and length-scaled.
     *
     * @param length Body length along the local z-axis.
     * @param minor Outer semi-axis in the local y-direction.
     * @param major Outer semi-axis in the local x-direction.
     * @return Triangulated standing-wave cavity mesh.
     */
    static MeshData getRFCavity(double length, double minor, double major);

    /**
     * @brief Build a hollow tube with a smoothly varying outer profile.
     *
     * The inner bore stays constant while the outer support follows a sequence
     * of longitudinal stations @f$(z_i, s_i)@f$ with linearly interpolated
     * scale factors @f$s_i@f$.
     *
     * @param length Body length along the local z-axis.
     * @param innerMinor Inner semi-axis in the local y-direction.
     * @param innerMajor Inner semi-axis in the local x-direction.
     * @param profile Ordered list of normalized longitudinal stations and outer-radius scales.
     * @param numSegments Number of azimuthal segments per ring.
     * @return Triangulated profiled tube mesh.
     */
    static MeshData getProfiledTube(
            double length, double innerMinor, double innerMajor,
            const std::vector<std::pair<double, double>>& profile, unsigned int numSegments = 36);

    /**
     * @brief Build a traveling-wave structure with repeated shallow corrugations.
     *
     * Compared to the standing-wave cavity, the modulation is more uniform and
     * elongated to suggest a periodic accelerating structure.
     *
     * @param length Body length along the local z-axis.
     * @param minor Outer semi-axis in the local y-direction.
     * @param major Outer semi-axis in the local x-direction.
     * @return Triangulated traveling-wave cavity mesh.
     */
    static MeshData getTravelingWave(double length, double minor, double major);

    static MeshData getBox(double length, double width, double height, double formFactor);

    std::vector<MeshData> elements_m;
    bool hasDriftReference_m;
    double driftMinor_m;
    double driftMajor_m;
};

#endif
