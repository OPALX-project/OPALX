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
class ScalingFFAMagnet;

/** Display-only geometry in the element frame, transformed once when added.
 * Empty triangles with nonempty vertices represent a consecutive polyline in
 * VTK exports (used for zero-field design arcs).
 */
class MeshData {
public:
    std::vector<Vector_t<double, 3>> vertices_m;
    std::vector<Vector_t<unsigned int, 3>> triangles_m;
    std::vector<std::pair<Vector_t<double, 3>, Vector_t<double, 3>>> decorations_m;
    int type_m{0};
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
     * OPALX's 1e6-m unbounded-aperture sentinel is not a physical support and
     * is rejected.
     *
     * @return true if a positive, bounded support size is available, false otherwise.
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

    /** Build display geometry without changing the element or tracking aperture.
     * SBENDs without a finite aperture use the representative support size, or
     * 0.05 m transverse half-sizes when no representative support is available.
     */
    void add(const ElementBase& element);

    /** Add a grey detector plane at global laboratory X-Z endpoints [m].
     * PROBE is vertically unbounded for tracking. The 1 m display height only
     * makes its orientation visible; it imposes no detector acceptance limit.
     * These vertices are already global and must not undergo element placement.
     */
    void addProbe(const Vector_t<double, 3>& start, const Vector_t<double, 3>& end);

    void write(const std::string& fname);

private:
    friend class MeshGeneratorTest;
    enum MeshType {
        OTHER = 0,
        DIPOLE,
        QUADRUPOLE,
        SEXTUPOLE,
        OCTUPOLE,
        SOLENOID,
        RFCAVITY,
        TRAVELINGWAVE,
        DRIFT,
        SCALINGFFA,
        PROBE
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
     * @brief Build a solid sector-bend body following the reference arc.
     *
     * Cross sections remain normal to the local design-orbit tangent. The
     * horizontal and vertical half sizes come from the element aperture.
     *
     * @param arcLength Reference-orbit arc length.
     * @param curvature Signed reference curvature.
     * @param horizontalHalfSize Horizontal half size in the local radial direction.
     * @param verticalHalfSize Vertical half size.
     * @return Triangulated curved dipole body.
     */
    static MeshData getSBend(
            double arcLength, double curvature, double horizontalHalfSize, double verticalHalfSize);

    /** Nominal ScalingFFA display volume, not a mechanical yoke or fringe envelope.
     * Uses the resolved field chart and nominal length; spiral faces are sampled
     * radially. The nominal spiral interval is centred on the resolved phiStart;
     * its width is arcLength/abs(R0). Height and radial limits are the configured
     * support dimensions, not inferred mechanical dimensions. Both signs of R0
     * use the inverse of the field-support chart. Sampling is bounded and purely
     * for display; it does not set tracking accuracy or an iso-field threshold.
     * Zero-field geometry arcs produce only a design-orbit polyline.
     * No aperture, field or tracking configuration is changed.
     */
    static MeshData getScalingFFA(const ScalingFFAMagnet& magnet);

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
