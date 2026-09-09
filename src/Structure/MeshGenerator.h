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
#include "Algorithms/Matrix.h"
#include "OPALTypes.h"

class ElementBase;

class MeshData {
public:
    std::vector<Vector_t<double, 3>> vertices_m;
    std::vector<Vector_t<unsigned int, 3>> triangles_m;
    std::vector<std::pair<Vector_t<double, 3>, Vector_t<double, 3>>> decorations_m;
    int type_m;
};

/**
 * @brief Serializable description of one beamline element for external tools.
 *
 * Carries everything a downstream consumer needs to build geometry on demand:
 * the element's global pose (local-to-global rotation and the local origin in
 * global coordinates), its length along the local z-axis, and its transverse
 * extent in the local x/y semi-axes.
 */
struct ElementDescription {
    std::string name;
    std::string type;
    double length;
    Vector_t<double, 3> position;
    matrix3x3_t rotation;
    bool hasTransverse;
    std::string transverseKind;
    double semiX;
    double semiY;
    int multipoleOrder;
    std::vector<Vector_t<double, 3>> designPath;
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

    /**
     * @brief Record a lightweight description of an element for external consumers.
     *
     * Unlike add(), describe() never builds geometry: it records the element's
     * name, type, length, global pose, and transverse extent so that a separate
     * tool (e.g. a python script) can generate meshes on demand. Every element of
     * the beamline is recorded, including those that add() skips or would mesh
     * at a fallback radius. Transverse extents of unset ("infinite") apertures
     * are recorded as-is; the consumer is expected to clamp them.
     */
    void describe(const ElementBase& element);

    void write(const std::string& fname);

    /**
     * @brief Write the recorded element descriptions as JSON.
     *
     * The file <fname>_ElementPositions.json is placed in the auxiliary output
     * directory next to the legacy python mesh output. It also reports the
     * representative drift radius (reference_transverse) that consumers should
     * use to clamp elements without a physical transverse extent.
     */
    void writeDescriptions(const std::string& fname) const;

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
    std::vector<ElementDescription> descriptions_m;
    bool hasDriftReference_m;
    double driftMinor_m;
    double driftMajor_m;
};

#endif
