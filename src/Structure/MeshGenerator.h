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

    static MeshData getBox(double length, double width, double height, double formFactor);

    std::vector<MeshData> elements_m;
    bool hasDriftReference_m;
    double driftMinor_m;
    double driftMajor_m;
};

#endif
